from torch_geometric.nn import GATConv
from torch_geometric.nn.conv import PDNConv
import torch_geometric.nn as gnn
import torch.nn as nn

from Models.common import NONLINEARITIES


class GNN(nn.Module):   
    def __init__(self, config, edge_dim):
        super().__init__()
        self.convs = nn.ModuleList()
        self.use_graph_norm = config.use_graph_norm
        self.edge_dim = edge_dim

        if self.use_graph_norm:
            self.batch_norms = nn.ModuleList()
            for i in range(config.num_layers):
                self.batch_norms.append(gnn.norm.GraphNorm(config.hidden_dim))

        for i in range(config.num_layers):
            self.convs.append(self._get_layer(config))

        self.activation = NONLINEARITIES[config.act_fn]

    def _get_layer(self, config):
        if config.layer_name == 'GAT':
            return GATConv(config.hidden_dim, config.hidden_dim, heads = config.n_head, dropout = config.dropout, edge_dim = self.edge_dim)
        elif config.layer_name == 'PDN':
            return PDNConv(config.hidden_dim, config.hidden_dim, hidden_channels=config.hidden_dim, dropout = config.dropout, edge_dim = self.edge_dim)

    def forward(self, x, edge_index, edge_attr, batch):

        # Run everything through the graph convolutions and activations
        for i in range(len(self.convs)):
            x = self.convs[i](x = x, edge_index = edge_index, edge_attr = edge_attr)

            if self.use_graph_norm:
                x = self.batch_norms[i](x, batch=batch)

            x = self.activation(x)
        
        # pool with respect to batches
        return x