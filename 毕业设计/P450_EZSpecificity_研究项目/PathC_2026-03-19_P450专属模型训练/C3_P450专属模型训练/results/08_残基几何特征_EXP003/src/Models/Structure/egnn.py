from torch import nn
import torch
from torch_scatter import scatter_sum, scatter_mean

from Models.common import MLP

class EnBaseLayer(nn.Module):
    def __init__(self, hidden_dim, edge_feat_dim, act_fn='relu', norm='None', attention=False, residual=True):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.edge_feat_dim = edge_feat_dim
        self.act_fn = act_fn
        self.norm = norm
        self.attention = attention
        self.residual = residual

        self.edge_mlp = MLP(2 * hidden_dim + edge_feat_dim, hidden_dim, hidden_dim,
                            num_layer=2, norm=norm, act_fn=act_fn, act_last=True)
        if self.attention:
            self.edge_inf = nn.Sequential(nn.Linear(hidden_dim, 1), nn.Sigmoid())

        self.node_mlp = MLP(2 * hidden_dim, hidden_dim, hidden_dim, num_layer=2, norm=norm, act_fn=act_fn)

    def forward(self, h, edge_index, edge_attr):
        """Forward pass of the linear layer
        Args:
            h: dict of node-features
            x: coordinates
            G: minibatch of (homo)graphs
        Returns:
            tensor with new features [B, n_points, n_features_out]
        """
        edge_index = edge_index
        src, dst = edge_index
        if self.edge_feat_dim > 0:
            edge_feat = edge_attr  # shape: [#edges_in_batch, #bond_types]
        else:
            edge_feat = None

        hi, hj = h[dst], h[src]

        if edge_feat is None:
            mij = self.edge_mlp(torch.cat([hi, hj], -1))
        else:
            mij = self.edge_mlp(torch.cat([edge_feat, hi, hj], -1))

        if self.attention:
            eij = mij * self.edge_inf(mij)
        else:
            eij = mij

        mi = scatter_sum(eij, dst, dim=0, dim_size=h.shape[0])

        h = h + self.node_mlp(torch.cat([mi, h], -1))

        return h
    
class EGNN(nn.Module):
    def __init__(self, config, edge_dim):
        super().__init__()
        # Build the network
        self.num_layers = config.num_layers
        self.hidden_dim = config.hidden_dim
        self.edge_feat_dim = edge_dim
        self.act_fn = config.act_fn
        self.norm = config.norm
        self.attention = config.attention
        self.residual = config.residual
        
        self.net = self._build_network()

    def _build_network(self):
        # Equivariant layers
        layers = []
        for l_idx in range(self.num_layers):
            layer = EnBaseLayer(self.hidden_dim, self.edge_feat_dim, act_fn=self.act_fn, norm=self.norm, attention=self.attention, residual=self.residual)
            layers.append(layer)
        return nn.ModuleList(layers)

    def forward(self, x, edge_index, edge_attr, batch):
        edge_index = edge_index
        h = x
        for interaction in self.net:
            h = h + interaction(h, edge_index, edge_attr)
        return h