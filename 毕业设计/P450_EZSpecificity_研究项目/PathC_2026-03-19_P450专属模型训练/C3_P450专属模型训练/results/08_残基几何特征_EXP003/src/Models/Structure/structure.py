import pytorch_lightning as pl
import torch
import torch.nn as nn
from torch_scatter import scatter

from Models.Structure.gnn import GNN
from Models.Structure.egnn import EGNN
import Datasets.Structure.transforms as trans

def get_encoder(config, edge_dim):
    if config.name == 'gnn':
        net = GNN(
            config=config,
            edge_dim=edge_dim
        )
    elif config.name == 'egnn':
        net = EGNN(
            config=config,
            edge_dim=edge_dim
        )
    else:
        raise ValueError(config.name)
    return net

class Graph(nn.Module):
    def __init__(self, config, trans_config):
        super().__init__()

        # Build the network
        self.config = config
        self.hidden_dim = config.hidden_dim

        protein_featurizer = trans.FeaturizeProteinAtom()
        ligand_featurizer = trans.FeaturizeLigandAtom()
        add_knn_edge_featurizer = trans.EdgeConnection(dist_noise=trans_config.dist_noise, cutoff=trans_config.cutoff, num_r_gaussian=trans_config.num_r_gaussian, k=trans_config.k)

        self.protein_atom_emb = nn.Linear(protein_featurizer.feature_dim, self.hidden_dim)
        self.ligand_atom_emb = nn.Linear(ligand_featurizer.feature_dim, self.hidden_dim)

        if config.mode == 'complex' or config.mode == 'both':
            self.encoder_complex = get_encoder(config, edge_dim = add_knn_edge_featurizer.feature_dim)
        
        if config.mode == 'sub' or config.mode == 'both':
            self.encoder_sub = get_encoder(config, edge_dim = add_knn_edge_featurizer.feature_dim)

    def forward(self, G):
        h_protein = self.protein_atom_emb(G.protein_x)
        h_ligand = self.ligand_atom_emb(G.ligand_x)
        
        h = h_protein * G.protein_mask[:, None] + h_ligand * G.ligand_mask[:, None]
        batch = G.ligand_index_batch

        if self.config.mode == 'both' or self.config.mode == 'complex':
            h_complex = self.encoder_complex(
                x = h,
                edge_index = G.complex_edge_index,
                edge_attr = G.complex_edge_attr,
                batch = batch
            )  # (N_p+N_l, H)

        if self.config.mode == 'both' or self.config.mode == 'sub':
            h_sub = self.encoder_sub(
                x = h,
                edge_index = G.complex_edge_index,
                edge_attr = G.complex_edge_attr,
                batch = batch
            )
        
        if self.config.mode == 'both':
            h = h_complex * G.mask_use_complex[:, None] + h_sub * G.mask_use_ligand[:, None]
        elif self.config.mode == 'complex':
            h = h_complex
        elif self.config.mode == 'sub':
            h = h_sub
       
        # Aggregate messages
        pre_atom_output = (h, batch, G.ligand_index)

        output = scatter(h, index=batch, dim=0, reduce='mean')  # (N, F)
        # output_ligand = scatter(h * G.ligand_mask[:, None], index=batch, dim=0, reduce='sum') / (scatter(G.ligand_mask, index=batch, dim=0, reduce='sum') + 1e-07)[:, None]  # (N, F)
        # output_protein = scatter(h * G.protein_mask[:, None], index=batch, dim=0, reduce='sum') / (scatter(G.protein_mask, index=batch, dim=0, reduce='sum') + 1e-07)[:, None]  # (N, F)
        
        # return output, output_ligand, output_protein, pre_atom_output
        return output, pre_atom_output