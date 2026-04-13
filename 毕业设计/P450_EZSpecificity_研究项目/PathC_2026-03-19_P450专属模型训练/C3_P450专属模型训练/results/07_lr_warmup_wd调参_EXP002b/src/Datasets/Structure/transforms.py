import sys

root_dir = "/projects/bbto/suyufeng/enzyme_specificity"
sys.path.append(f"{root_dir}/src")

import torch
import torch.nn.functional as F
import numpy as np

from rdkit.Chem.rdchem import BondType

from Datasets.Structure.utils import StructureComplexData
from Datasets.Structure.protein_ligand import ATOM_FEATS
from torch_geometric.nn import knn_graph
import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()})


class FeaturizeProteinAtom(object):

    def __init__(self):
        super().__init__()
        self.atomic_numbers = torch.LongTensor([1, 6, 7, 8, 16, 34, 26])    # H, C, N, O, S, Se, Fe
        self.max_num_aa = 22

    @property
    def feature_dim(self):
        return self.atomic_numbers.size(0) + self.max_num_aa + 2

    def __call__(self, data: StructureComplexData):
        element = data.protein_element.view(-1, 1) == self.atomic_numbers.view(1, -1)   # (N_atoms, N_elements)
        amino_acid = F.one_hot(data.protein_atom_to_aa_type.long(), num_classes=self.max_num_aa)
        is_backbone = data.protein_is_backbone.view(-1, 1).long()
        is_hetero = data.protein_is_hetero.view(-1, 1).long()
        x = torch.cat([element, amino_acid, is_backbone, is_hetero], dim=-1)
        data.protein_atom_feature = x
        return data

class FeaturizeLigandAtom(object):
    
    def __init__(self):
        super().__init__()
        self.atomic_numbers = torch.LongTensor([1, 6, 7, 8, 9, 15, 16, 17, 26, 35, 53])  # H, C, N, O, F, P, S, Cl, Fe, Br, I
        # self.n_degree = torch.LongTensor([0, 1, 2, 3, 4, 5])  # 0 - 5
        # self.n_num_hs = 6  # 0 - 5

    @property
    def num_properties(self):
        return sum(ATOM_FEATS.values())

    @property
    def feature_dim(self):
        return self.atomic_numbers.size(0) + self.num_properties

    def __call__(self, data: StructureComplexData):
        element = data.ligand_element.view(-1, 1) == self.atomic_numbers.view(1, -1)   # (N_atoms, N_elements)
        # convert some features to one-hot vectors
        atom_feature = []
        for i, (k, v) in enumerate(ATOM_FEATS.items()):
            feat = data.ligand_atom_feature[:, i:i+1]
            if v > 1:
                feat = (feat == torch.LongTensor(list(range(v))).view(1, -1))
            else:
                if k == 'AtomicNumber':
                    feat = feat / 100.
            atom_feature.append(feat)

        # data.ligand_atom_feature: atomic_number, aromatic, degree, num_hs
        # try:
        #     # degree = F.one_hot(data.ligand_atom_feature[:, 2], num_classes=self.n_degree)
        #     degree = data.ligand_atom_feature[:, 2:3] == self.n_degree.view(1, -1)
        #     num_hs = F.one_hot(data.ligand_atom_feature[:, 3], num_classes=self.n_num_hs)
        # except:
        #     print(torch.max(data.ligand_atom_feature[:, 2]))
        #     print(torch.max(data.ligand_atom_feature[:, 3]))
        #     raise ValueError
        atom_feature = torch.cat(atom_feature, dim=-1)
        data.ligand_atom_feature_full = torch.cat([element, atom_feature], dim=-1)
        return data


from rdkit.Chem.rdchem import BondType
from Models.common import GaussianSmearing

class EdgeConnection(object):
    def __init__(self, dist_noise=True, cutoff=12, num_r_gaussian=24, k=48):
        super(EdgeConnection, self).__init__()
        self.dist_noise = dist_noise
        self.k = k
        self.num_r_gaussian = num_r_gaussian
        self.cutoff = cutoff
        self.distance_expansion = GaussianSmearing(stop=cutoff, num_gaussians=num_r_gaussian)

    @property
    def feature_dim(self):
        return 6 + 1 + self.num_r_gaussian
    
    def _build_knn_graph(self, data):
        pos = torch.cat([data.ligand_pos, data.protein_pos], dim=0)
        dict = {(x, y): 1 for x, y in zip(data.ligand_bond_index[0, :], data.ligand_bond_index[1, :])}
        knn_index = knn_graph(pos, k=self.k, flow='target_to_source')
        result = []
        for edge in knn_index.permute(1, 0):
            if (edge[0], edge[1]) not in dict:
                result.append(edge[:, None])
        if len(result) == 0:
            return torch.zeros((2, 0), dtype=torch.long)
        return torch.cat(result, dim=1)
    
    def _get_dist(self, edge_index, pos):
        dst, src = edge_index
        dist = torch.norm(pos[dst] - pos[src], p=2, dim=-1)
        if self.dist_noise:
            dist = dist + torch.from_numpy(np.random.laplace(0.001994, 0.031939, (dist.size(0),))).float()
        # mask = dist < self.cutoff
        # dist = dist[mask]
        # edge_index = edge_index[:, mask]
        return edge_index, self.distance_expansion(dist)
    
    def __call__(self, data):
        N_substrate = data.ligand_pos.size(0)
        N_protein = data.protein_pos.size(0)

        if data.str_tag == 'ligand':
            knn_index = torch.zeros((2, 0), dtype=torch.long)
        else:
            knn_index = self._build_knn_graph(data)

        # dist matrix
        pos = torch.cat([data.ligand_pos, data.protein_pos], dim=0)
        knn_index, dist_feat_knn = self._get_dist(knn_index, pos)
        data.ligand_bond_index, dist_feat_real =self._get_dist(data.ligand_bond_index, pos)
        dist_feat = torch.cat([dist_feat_knn, dist_feat_real], dim=0)
        # bond type
        data.ligand_bond_type[data.ligand_bond_type==12] = 5
        ligand_bond_feature_real = F.one_hot((data.ligand_bond_type - 1).long(), num_classes=6)
        ligand_bond_feature_knn = F.one_hot(torch.ones((knn_index.shape[1]),dtype=torch.long) * 5, num_classes=6)
        bond_feat = torch.cat([ligand_bond_feature_knn, ligand_bond_feature_real], dim=0)
        # cross bond
        cross_bond = ((knn_index[0, :] < N_substrate).logical_and(knn_index[1, :] >= N_substrate)).logical_or((knn_index[0, :] >= N_substrate).logical_and(knn_index[1, :] < N_substrate))
        cross_bond = torch.cat([cross_bond, torch.zeros([data.ligand_bond_index.size(1)], dtype=bool)], dim=0)[:, None].int()
        
        data.complex_edge_attr = torch.cat([bond_feat, dist_feat, cross_bond], dim=-1)  # [E, F]
        data.protein_x = torch.cat([torch.zeros((N_substrate, data.protein_atom_feature.shape[1])), data.protein_atom_feature], dim=0)
        data.ligand_x = torch.cat([data.ligand_atom_feature_full, torch.zeros((N_protein, data.ligand_atom_feature_full.shape[1]))], dim=0)
        data.ligand_mask = torch.cat([torch.ones((N_substrate,)), torch.zeros((N_protein,))], dim=0)
        data.protein_mask = torch.cat([torch.zeros((N_substrate,)), torch.ones((N_protein,))], dim=0)
        data.complex_edge_index = torch.cat([data.ligand_bond_index, knn_index], dim=-1)
        data.ligand_index = torch.cat([data.ligand_index, torch.ones((N_protein, )).long() * 280], dim=-1)

        # print(data)
        data.ligand_atom_feature_full = None
        data.protein_atom_feature = None
        data.protein_pos = None
        data.ligand_pos = None
        data.ligand_bond_index = None
        data.ligand_bond_type = None
                
        data.protein_molecule_name = None
        data.protein_filename = None
        # data.y = None
        
        return data
    

if __name__ == "__main__":
    N_protein = 10
    N_substrate = 20

    protein_pos = torch.rand([N_protein, 3])
    ligand_pos = torch.rand([N_substrate, 3])

    protein_atom_feature = torch.randint(0, 2, [N_protein, 10])
    ligand_atom_feature_full = torch.randint(0, 2, [N_substrate, 10])

    ligand_bond_index = torch.randint(0, N_substrate, [2, 10])
    ligand_bond_type = torch.randint(1, len(BondType.names.values()), [10])
    ligand_index = torch.randint(0, 280, [N_substrate])

    data = StructureComplexData(protein_pos=protein_pos, ligand_pos=ligand_pos, ligand_bond_index=ligand_bond_index, ligand_bond_type=ligand_bond_type, protein_atom_feature=protein_atom_feature, ligand_atom_feature_full=ligand_atom_feature_full, ligand_index=ligand_index)
    data = EdgeConnection()(data)

    print(data)