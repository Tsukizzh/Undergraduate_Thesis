# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

import torch
import numpy as np


def pad_1d_unsqueeze(x, padlen):
    x = x + 1 # pad id = 0
    xlen = x.size(0)
    if xlen < padlen:
        new_x = x.new_zeros([padlen], dtype=x.dtype)
        new_x[:xlen] = x
        x = new_x
    return x.unsqueeze(0)


def pad_2d_unsqueeze(x, padlen):
    x = x + 1 # pad id = 0
    xlen, xdim = x.size()
    if xlen < padlen:
        new_x = x.new_zeros([padlen, xdim], dtype=x.dtype)
        new_x[:xlen, :] = x
        x = new_x
    return x.unsqueeze(0)


def pad_attn_bias_unsqueeze(x, padlen):
    xlen = x.size(0)
    if xlen < padlen:
        new_x = x.new_zeros([padlen, padlen], dtype=x.dtype).fill_(float('-inf'))
        new_x[:xlen, :xlen] = x
        new_x[xlen:, :xlen] = 0
        x = new_x
    return x.unsqueeze(0)


def pad_edge_type_unsqueeze(x, padlen):
    xlen = x.size(0)
    if xlen < padlen:
        new_x = x.new_zeros([padlen, padlen, x.size(-1)], dtype=x.dtype)
        new_x[:xlen, :xlen, :] = x
        x = new_x
    return x.unsqueeze(0)


def pad_rel_pos_unsqueeze(x, padlen):
    x = x + 1
    xlen = x.size(0)
    if xlen < padlen:
        new_x = x.new_zeros([padlen, padlen], dtype=x.dtype)
        new_x[:xlen, :xlen] = x
        x = new_x
    return x.unsqueeze(0)


def pad_rel_pos_3d_unsqueeze(x, padlen):
    xlen = x.size(0)
    if xlen < padlen:
        new_x = x.new_zeros([padlen, padlen], dtype=x.dtype)
        new_x[:xlen, :xlen] = x
        x = new_x
    return x.unsqueeze(0)


def pad_3d_unsqueeze(x, padlen1, padlen2, padlen3):
    x = x + 1
    xlen1, xlen2, xlen3, xlen4 = x.size()
    if xlen1 < padlen1 or xlen2 < padlen2 or xlen3 < padlen3:
        new_x = x.new_zeros([padlen1, padlen2, padlen3, xlen4], dtype=x.dtype)
        new_x[:xlen1, :xlen2, :xlen3, :] = x
        x = new_x
    return x.unsqueeze(0)


class Batch():
    def __init__(self, idx, attn_bias, attn_edge_type, graph_dist, geo_dist, in_degree, out_degree,
                 protein_x, ligand_x, edge_input, y, num_protein_nodes, num_ligand_nodes):
        super(Batch, self).__init__()
        self.idx = idx
        self.in_degree, self.out_degree = in_degree, out_degree
        self.protein_x, self.ligand_x, self.y = protein_x, ligand_x, y
        self.attn_bias, self.attn_edge_type, self.rel_pos = attn_bias, attn_edge_type, graph_dist
        self.edge_input = edge_input
        self.geo_dist = geo_dist
        self.graph_dist = graph_dist
        self.num_protein_nodes = num_protein_nodes
        self.num_ligand_nodes = num_ligand_nodes

    def to(self, device):
        self.idx = self.idx.to(device)
        self.in_degree, self.out_degree = self.in_degree.to(device), self.out_degree.to(device)
        self.x, self.y = self.x.to(device), self.y.to(device)
        self.protein_x = self.protein_x.to(device)
        self.ligand_x = self.ligand_x.to(device)
        self.attn_bias, self.attn_edge_type, self.rel_pos = self.attn_bias.to(device), self.attn_edge_type.to(device), self.rel_pos.to(device)
        self.edge_input = self.edge_input.to(device)
        self.geo_dist = self.geo_dist.to(device)
        self.graph_dist = self.graph_dist.to(device)
        return self

    def __len__(self):
        return self.in_degree.size(0)


def collator(items, max_node=1024, multi_hop_max_dist=20, rel_pos_max=20):
    items = [item for item in items if item is not None and item.protein_x.size(0) + item.ligand_x.size(0) <= max_node]
    items_ = [(item.id, item.attn_bias, item.attn_edge_type, item.rel_pos, item.in_degree, item.out_degree,
               item.protein_x, item.ligand_x, item.edge_input[:, :, :multi_hop_max_dist, :], item.y) for item in items]
    idxs, attn_biases, attn_edge_types, rel_poses, in_degrees, out_degrees, p_xs, l_xs, edge_inputs, ys = zip(*items_)

    items_ = [(item.geo_dist,) for item in items]
    all_rel_pos_3d_1s, = zip(*items_)

    for idx, _ in enumerate(attn_biases):
        attn_biases[idx][1:, 1:][rel_poses[idx] >= rel_pos_max] = float('-inf')
    max_dist = max(i.size(-2) for i in edge_inputs)
    y = torch.cat(ys)
    # max_node_num = max(i.size(0) for i in xs)
    # x = torch.cat([pad_2d_unsqueeze(i, max_node_num) for i in xs])
    max_p_node_num = max(i.size(0) for i in p_xs)
    max_l_node_num = max(i.size(0) for i in l_xs)
    protein_x = torch.cat([pad_2d_unsqueeze(i, max_p_node_num) for i in p_xs])
    ligand_x = torch.cat([pad_2d_unsqueeze(i, max_l_node_num) for i in l_xs])
    max_node_num = max_p_node_num + max_l_node_num

    edge_input = torch.cat([pad_3d_unsqueeze(i, max_node_num, max_node_num, max_dist) for i in edge_inputs])
    attn_bias = torch.cat([pad_attn_bias_unsqueeze(i, max_node_num + 1) for i in attn_biases])
    attn_edge_type = torch.cat([pad_edge_type_unsqueeze(i, max_node_num) for i in attn_edge_types])
    rel_pos = torch.cat([pad_rel_pos_unsqueeze(i, max_node_num) for i in rel_poses])
    all_rel_pos_3d_1 = torch.cat([pad_rel_pos_3d_unsqueeze(i, max_node_num) for i in all_rel_pos_3d_1s])
    in_degree = torch.cat([pad_1d_unsqueeze(i, max_node_num) for i in in_degrees])
    out_degree = torch.cat([pad_1d_unsqueeze(i, max_node_num) for i in out_degrees])
    return Batch(
        idx=torch.LongTensor(idxs),
        attn_bias=attn_bias,
        attn_edge_type=attn_edge_type,
        graph_dist=rel_pos,
        geo_dist=all_rel_pos_3d_1,
        in_degree=in_degree,
        out_degree=out_degree,
        protein_x=protein_x,
        ligand_x=ligand_x,
        edge_input=edge_input,
        y=y,
        num_protein_nodes=[i.size(0) for i in p_xs],
        num_ligand_nodes=[i.size(0) for i in l_xs]
    )