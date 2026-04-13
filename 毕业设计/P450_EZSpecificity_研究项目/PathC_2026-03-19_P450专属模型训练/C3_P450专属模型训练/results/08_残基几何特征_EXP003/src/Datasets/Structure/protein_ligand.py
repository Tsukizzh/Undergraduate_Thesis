import os
import numpy as np
from rdkit import Chem
from rdkit.Chem.rdchem import BondType, HybridizationType
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import torch
import torch.nn.functional as F
from torch_scatter import scatter

ATOM_FAMILIES = ['Acceptor', 'Donor', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe', 'NegIonizable', 'PosIonizable', 'ZnBinder']
ATOM_FAMILIES_ID = {s: i for i, s in enumerate(ATOM_FAMILIES)}
ATOM_FEATS = {'AtomicNumber': 1, 'Aromatic': 1, 'Degree': 6, 'NumHs': 6, 'Hybridization': len(HybridizationType.values)}
BOND_TYPES = {t: i for i, t in enumerate(BondType.names.values())}
BOND_NAMES = {i: t for i, t in enumerate(BondType.names.keys())}


def _calc_dihedral(p0, p1, p2, p3):
    """Compute dihedral angle in radians from 4 points. Returns None if any is missing."""
    if any(p is None for p in (p0, p1, p2, p3)):
        return None
    p0 = np.asarray(p0, dtype=np.float32)
    p1 = np.asarray(p1, dtype=np.float32)
    p2 = np.asarray(p2, dtype=np.float32)
    p3 = np.asarray(p3, dtype=np.float32)
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1_norm = np.linalg.norm(b1)
    if b1_norm < 1e-6:
        return None
    b1 = b1 / b1_norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    v_norm = np.linalg.norm(v)
    w_norm = np.linalg.norm(w)
    if v_norm < 1e-6 or w_norm < 1e-6:
        return None
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.arctan2(y, x))



class PDBProtein(object):

    AA_NAME_SYM = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
        'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
        'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'UNK': 'X'
    }

    AA_NAME_NUMBER = {
        k: i for i, (k, _) in enumerate(AA_NAME_SYM.items())
    }
    HEM_AA_TYPE = len(AA_NAME_NUMBER)      # = 21
    AA_NAME_NUMBER['HEM'] = HEM_AA_TYPE
    AA_NAME_NUMBER['HEC'] = HEM_AA_TYPE    # heme C variant
    HEM_RES_NAMES = {'HEM', 'HEC'}

    BACKBONE_NAMES = ["CA", "C", "N", "O"]

    RESIDUE_ANGLE_DIM = 6
    CHI1_ATOMS = {
        'ARG': ("N", "CA", "CB", "CG"),
        'ASN': ("N", "CA", "CB", "CG"),
        'ASP': ("N", "CA", "CB", "CG"),
        'CYS': ("N", "CA", "CB", "SG"),
        'GLN': ("N", "CA", "CB", "CG"),
        'GLU': ("N", "CA", "CB", "CG"),
        'HIS': ("N", "CA", "CB", "CG"),
        'ILE': ("N", "CA", "CB", "CG1"),
        'LEU': ("N", "CA", "CB", "CG"),
        'LYS': ("N", "CA", "CB", "CG"),
        'MET': ("N", "CA", "CB", "CG"),
        'PHE': ("N", "CA", "CB", "CG"),
        'PRO': ("N", "CA", "CB", "CG"),
        'SER': ("N", "CA", "CB", "OG"),
        'THR': ("N", "CA", "CB", "OG1"),
        'TRP': ("N", "CA", "CB", "CG"),
        'TYR': ("N", "CA", "CB", "CG"),
        'VAL': ("N", "CA", "CB", "CG1"),
    }

    def __init__(self, data, mode='auto'):
        super().__init__()
        self.fn = data
        if (data[-4:].lower() == '.pdb' and mode == 'auto') or mode == 'path':
            with open(data, 'r') as f:
                self.block = f.read()
        else:
            self.block = data

        self.ptable = Chem.GetPeriodicTable()

        # Molecule properties
        self.title = None
        # Atom properties
        self.atoms = []
        self.element = []
        self.atomic_weight = []
        self.pos = []
        self.atom_name = []
        self.is_backbone = []
        self.is_hetero = []
        self.atom_to_aa_type = []
        # Residue properties
        self.residues = []
        self.amino_acid = []
        self.center_of_mass = []
        self.pos_CA = []
        self.pos_C = []
        self.pos_N = []
        self.pos_O = []

        self._parse()

    def _enum_formatted_atom_lines(self):
        for line in self.block.splitlines():
            record_type = line[0:6].strip()
            if record_type in ('ATOM', 'HETATM'):
                res_name = line[17:20].strip()
                if record_type == 'HETATM' and res_name not in self.HEM_RES_NAMES:
                    continue
                element_symb = line[76:78].strip().capitalize()
                if len(element_symb) == 0:
                    element_symb = line[13:14]
                yield {
                    'line': line,
                    'type': record_type,
                    'atom_id': int(line[6:11]),
                    'atom_name': line[12:16].strip(),
                    'res_name': res_name,
                    'chain': line[21:22].strip(),
                    'res_id': int(line[22:26]),
                    'res_insert_id': line[26:27].strip(),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'occupancy': float(line[54:60]),
                    'segment': line[72:76].strip(),
                    'element_symb': element_symb,
                    'charge': line[78:80].strip(),
                }
            elif line[0:6].strip() == 'HEADER':
                yield {
                    'type': 'HEADER',
                    'value': line[10:].strip()
                }
            elif line[0:6].strip() == 'ENDMDL':
                break   # Some PDBs have more than 1 model.

    def _parse(self):
        # Process atoms
        residues_tmp = {}
        for atom in self._enum_formatted_atom_lines():
            if atom['type'] == 'HEADER':
                self.title = atom['value'].lower()
                continue
            self.atoms.append(atom)
            atomic_number = self.ptable.GetAtomicNumber(atom['element_symb'])
            next_ptr = len(self.element)
            self.element.append(atomic_number)
            self.atomic_weight.append(self.ptable.GetAtomicWeight(atomic_number))
            self.pos.append(np.array([atom['x'], atom['y'], atom['z']], dtype=np.float32))
            self.atom_name.append(atom['atom_name'])
            is_hetero = atom['type'] == 'HETATM'
            self.is_backbone.append((atom['atom_name'] in self.BACKBONE_NAMES) and not is_hetero)
            self.is_hetero.append(is_hetero)
            if is_hetero:
                self.atom_to_aa_type.append(self.HEM_AA_TYPE)
            elif atom['res_name'] not in self.AA_NAME_NUMBER:
                self.atom_to_aa_type.append(self.AA_NAME_NUMBER['UNK'])
            else:
                self.atom_to_aa_type.append(self.AA_NAME_NUMBER[atom['res_name']])

            chain_res_id = '%s_%s_%d_%s' % (atom['chain'], atom['segment'], atom['res_id'], atom['res_insert_id'])
            if chain_res_id not in residues_tmp:
                residues_tmp[chain_res_id] = {
                    'name': atom['res_name'],
                    'atoms': [next_ptr],
                    'chain': atom['chain'],
                    'segment': atom['segment'],
                    'is_hetero': is_hetero,
                    'atom_pos': {atom['atom_name']: self.pos[next_ptr]},
                }
            else:
                assert residues_tmp[chain_res_id]['name'] == atom['res_name']
                assert residues_tmp[chain_res_id]['chain'] == atom['chain']
                residues_tmp[chain_res_id]['atoms'].append(next_ptr)
                residues_tmp[chain_res_id]['atom_pos'].setdefault(atom['atom_name'], self.pos[next_ptr])

        # Process residues
        self.residues = [r for _, r in residues_tmp.items()]
        for residue in self.residues:
            sum_pos = np.zeros([3], dtype=np.float32)
            sum_mass = 0.0
            for atom_idx in residue['atoms']:
                sum_pos += self.pos[atom_idx] * self.atomic_weight[atom_idx]
                sum_mass += self.atomic_weight[atom_idx]
                if self.atom_name[atom_idx] in self.BACKBONE_NAMES:
                    residue['pos_%s' % self.atom_name[atom_idx]] = self.pos[atom_idx]
            residue['center_of_mass'] = sum_pos / sum_mass

        self._annotate_residue_angle_features()
        
        # Process backbone atoms of residues
        for residue in self.residues:
            if residue['name'] not in self.AA_NAME_NUMBER:
                self.amino_acid.append(self.AA_NAME_NUMBER['UNK'])
            else:
                self.amino_acid.append(self.AA_NAME_NUMBER[residue['name']])
            self.center_of_mass.append(residue['center_of_mass'])
            for name in self.BACKBONE_NAMES:
                pos_key = 'pos_%s' % name   # pos_CA, pos_C, pos_N, pos_O
                if pos_key in residue:
                    getattr(self, pos_key).append(residue[pos_key])
                else:
                    getattr(self, pos_key).append(residue['center_of_mass'])


    def _has_peptide_link(self, left_residue, right_residue):
        left_c = left_residue.get('pos_C')
        right_n = right_residue.get('pos_N')
        if left_c is None or right_n is None:
            return False
        return np.linalg.norm(left_c - right_n) < 2.0

    def _annotate_residue_angle_features(self):
        for residue in self.residues:
            residue['residue_angle_feature'] = np.zeros(self.RESIDUE_ANGLE_DIM, dtype=np.float32)

        chain_groups = {}
        for residue in self.residues:
            if residue.get('is_hetero', False):
                continue
            key = (residue['chain'], residue['segment'])
            chain_groups.setdefault(key, []).append(residue)

        for residues in chain_groups.values():
            for i, residue in enumerate(residues):
                feat = np.zeros(self.RESIDUE_ANGLE_DIM, dtype=np.float32)
                prev_r = residues[i - 1] if i > 0 else None
                next_r = residues[i + 1] if i + 1 < len(residues) else None

                phi = None
                if prev_r is not None and self._has_peptide_link(prev_r, residue):
                    phi = _calc_dihedral(
                        prev_r.get('pos_C'), residue.get('pos_N'),
                        residue.get('pos_CA'), residue.get('pos_C'))

                psi = None
                if next_r is not None and self._has_peptide_link(residue, next_r):
                    psi = _calc_dihedral(
                        residue.get('pos_N'), residue.get('pos_CA'),
                        residue.get('pos_C'), next_r.get('pos_N'))

                chi1 = None
                chi1_atoms = self.CHI1_ATOMS.get(residue['name'])
                if chi1_atoms is not None:
                    apos = residue.get('atom_pos', {})
                    chi1 = _calc_dihedral(*[apos.get(n) for n in chi1_atoms])

                if phi is not None:
                    feat[0], feat[1] = np.sin(phi), np.cos(phi)
                if psi is not None:
                    feat[2], feat[3] = np.sin(psi), np.cos(psi)
                if chi1 is not None:
                    feat[4], feat[5] = np.sin(chi1), np.cos(chi1)

                residue['residue_angle_feature'] = feat

    def to_dict_atom(self):
        residue_angle_feature = np.zeros((len(self.atoms), self.RESIDUE_ANGLE_DIM), dtype=np.float32)
        for residue in self.residues:
            residue_angle_feature[residue['atoms']] = residue.get(
                'residue_angle_feature',
                np.zeros(self.RESIDUE_ANGLE_DIM, dtype=np.float32),
            )
        return {
            'element': np.array(self.element, dtype=int),
            'molecule_name': self.title,
            'pos': np.array(self.pos, dtype=np.float32),
            'is_backbone': np.array(self.is_backbone, dtype=bool),
            'is_hetero': np.array(self.is_hetero, dtype=bool),
            'atom_name': self.atom_name,
            'atom_to_aa_type': np.array(self.atom_to_aa_type, dtype=int),
            'residue_angle_feature': residue_angle_feature,
        }

    def to_dict_residue(self):
        return {
            'amino_acid': np.array(self.amino_acid, dtype=int),
            'center_of_mass': np.array(self.center_of_mass, dtype=np.float32),
            'pos_CA': np.array(self.pos_CA, dtype=np.float32),
            'pos_C': np.array(self.pos_C, dtype=np.float32),
            'pos_N': np.array(self.pos_N, dtype=np.float32),
            'pos_O': np.array(self.pos_O, dtype=np.float32),
        }

    def query_residues_radius(self, center, radius, criterion='center_of_mass'):
        center = np.array(center).reshape(3)
        selected = []
        for residue in self.residues:
            distance = np.linalg.norm(residue[criterion] - center, ord=2)
            print(residue[criterion], distance)
            if distance < radius:
                selected.append(residue)
        return selected

    def query_residues_ligand(self, ligand, radius, criterion='min_atom_distance'):
        selected = []
        sel_idx = set()
        # The time-complexity is O(mn).
        for center in ligand['pos']:
            for i, residue in enumerate(self.residues):
                if criterion == 'min_atom_distance':
                    atom_pos = np.array(
                        [self.pos[atom_idx] for atom_idx in residue['atoms']],
                        dtype=np.float32,
                    )
                    distance = np.linalg.norm(atom_pos - center[None, :], axis=1).min()
                else:
                    distance = np.linalg.norm(residue[criterion] - center, ord=2)
                if distance < radius and i not in sel_idx:
                    selected.append(residue)
                    sel_idx.add(i)
        return selected

    def residues_to_pdb_block(self, residues, name='POCKET'):
        block = "HEADER    %s\n" % name
        block += "COMPND    %s\n" % name
        for residue in residues:
            for atom_idx in residue['atoms']:
                block += self.atoms[atom_idx]['line'] + "\n"
        block += "END\n"
        return block

def get_zero_protein_feature():
    return {
        'element': np.array([], dtype=int),
        'molecule_name': 'None',
        'pos': np.array([], dtype=np.float32),
        'is_backbone': np.array([], dtype=bool),
        'is_hetero': np.array([], dtype=bool),
        'atom_name': [],
        'atom_to_aa_type': np.array([], dtype=int),
        'residue_angle_feature': np.zeros((0, PDBProtein.RESIDUE_ANGLE_DIM), dtype=np.float32),
    }

def get_ligand_atom_features(rdmol):
    num_atoms = rdmol.GetNumAtoms()
    atomic_number = []
    aromatic = []
    # sp, sp2, sp3 = [], [], []
    hybrid = []
    degree = []
    for atom_idx in range(num_atoms):
        atom = rdmol.GetAtomWithIdx(atom_idx)
        atomic_number.append(atom.GetAtomicNum())
        aromatic.append(1 if atom.GetIsAromatic() else 0)
        hybridization = atom.GetHybridization()
        HYBRID_TYPES = {t: i for i, t in enumerate(HybridizationType.names.values())}
        hybrid.append(HYBRID_TYPES[hybridization])
        # sp.append(1 if hybridization == HybridizationType.SP else 0)
        # sp2.append(1 if hybridization == HybridizationType.SP2 else 0)
        # sp3.append(1 if hybridization == HybridizationType.SP3 else 0)
        degree.append(atom.GetDegree())
    node_type = torch.tensor(atomic_number, dtype=torch.long)

    row, col = [], []
    for bond in rdmol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
    row = torch.tensor(row, dtype=torch.long)
    col = torch.tensor(col, dtype=torch.long)
    hs = (node_type == 1).to(torch.float)
    num_hs = scatter(hs[row], col, dim_size=num_atoms).numpy()
    # need to change ATOM_FEATS accordingly
    feat_mat = np.array([atomic_number, aromatic, degree, num_hs, hybrid], dtype=int).transpose()
    return feat_mat


# used for fixing some errors in sdf file
def parse_sdf_file_text(path):
    with open(path, 'r') as f:
        sdf = f.read()

    sdf = sdf.splitlines()
    num_atoms, num_bonds = map(int, [sdf[3][0:3], sdf[3][3:6]])
    ptable = Chem.GetPeriodicTable()

    element, pos = [], []
    accum_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    accum_mass = 0.0
    for atom_line in map(lambda x:x.split(), sdf[4:4+num_atoms]):
        x, y, z = map(float, atom_line[:3])
        symb = atom_line[3]
        atomic_number = ptable.GetAtomicNumber(symb.capitalize())
        element.append(atomic_number)
        pos.append([x, y, z])

        atomic_weight = ptable.GetAtomicWeight(atomic_number)
        accum_pos += np.array([x, y, z]) * atomic_weight
        accum_mass += atomic_weight

    center_of_mass = np.array(accum_pos / accum_mass, dtype=np.float32)

    element = np.array(element, dtype=int)
    pos = np.array(pos, dtype=np.float32)
    BOND_TYPES = {t: i for i, t in enumerate(BondType.names.values())}
    bond_type_map = {
        1: BOND_TYPES[BondType.SINGLE],
        2: BOND_TYPES[BondType.DOUBLE],
        3: BOND_TYPES[BondType.TRIPLE],
        4: BOND_TYPES[BondType.AROMATIC],
        8: BOND_TYPES[BondType.UNSPECIFIED]
    }
    row, col, edge_type = [], [], []
    for bond_line in sdf[4+num_atoms:4+num_atoms+num_bonds]:
        start, end = int(bond_line[0:3])-1, int(bond_line[3:6])-1
        row += [start, end]
        col += [end, start]
        edge_type += 2 * [bond_type_map[int(bond_line[6:9])]]

    edge_index = np.array([row, col], dtype=int)
    edge_type = np.array(edge_type, dtype=int)

    perm = (edge_index[0] * num_atoms + edge_index[1]).argsort()
    edge_index = edge_index[:, perm]
    edge_type = edge_type[perm]

    data = {
        'element': element,
        'pos': pos,
        'bond_index': edge_index,
        'bond_type': edge_type,
        'center_of_mass': center_of_mass
    }
    return data


# used for preparing the dataset
def parse_sdf_file_mol(path, mol=None, heavy_only=True):
    if mol is None:
        mol = next(iter(Chem.SDMolSupplier(path, removeHs=heavy_only, sanitize=False)))
    feat_mat = get_ligand_atom_features(mol)

    # fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    # factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    # rdmol = next(iter(Chem.SDMolSupplier(path, removeHs=heavy_only)))
    # rd_num_atoms = rdmol.GetNumAtoms()
    # feat_mat = np.zeros([rd_num_atoms, len(ATOM_FAMILIES)], dtype=int)
    # for feat in factory.GetFeaturesForMol(rdmol):
    #     feat_mat[feat.GetAtomIds(), ATOM_FAMILIES_ID[feat.GetFamily()]] = 1

    ptable = Chem.GetPeriodicTable()

    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    pos = mol.GetConformer().GetPositions()

    element = []
    indexs = []
    accum_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    accum_mass = 0.0
    for atom_idx in range(num_atoms):
        atom = mol.GetAtomWithIdx(atom_idx)
        atomic_number = atom.GetAtomicNum()
        element.append(atomic_number)
        indexs.append(atom.GetAtomMapNum())
        x, y, z = pos[atom_idx]
        atomic_weight = ptable.GetAtomicWeight(atomic_number)
        accum_pos += np.array([x, y, z]) * atomic_weight
        accum_mass += atomic_weight
    center_of_mass = np.array(accum_pos / accum_mass, dtype=np.float32)
    element = np.array(element, dtype=int)
    pos = np.array(pos, dtype=np.float32)
    indexs = np.array(indexs, dtype=int)
    
    row, col, edge_type = [], [], []
    
    BOND_TYPES = {}
    BOND_TYPES[BondType.SINGLE] = 1
    BOND_TYPES[BondType.DOUBLE] = 2
    BOND_TYPES[BondType.TRIPLE] = 3
    BOND_TYPES[BondType.AROMATIC] = 4
    BOND_TYPES[BondType.UNSPECIFIED] = 5

    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        edge_type += 2 * [BOND_TYPES[bond.GetBondType()]]
    edge_index = np.array([row, col], dtype=int)
    edge_type = np.array(edge_type, dtype=int)
    perm = (edge_index[0] * num_atoms + edge_index[1]).argsort()
    edge_index = edge_index[:, perm]
    edge_type = edge_type[perm]

    data = {
        'element': element,
        'pos': pos,
        'bond_index': edge_index,
        'bond_type': edge_type,
        'center_of_mass': center_of_mass,
        'atom_feature': feat_mat,
        'index': indexs
    }
    return data