"""
GVP (Geometric Vector Perceptron) encoder — ported from EnzymeCAGE for EXP005.

Source: github.com/drorlab/gvp-pytorch (via EnzymeCAGE) + enzymecage/gvp/__init__.py
       + enzymecage/model.py:GVP_embedding

Port scope (training-time only, no data-pipeline helpers):
  - Tuple helpers: tuple_sum, tuple_cat, tuple_index, _norm_no_nan, _split, _merge
  - Primitives: GVP, _VDropout, Dropout, LayerNorm
  - Conv stack: GVPConv, GVPConvLayer
  - Top wrapper: GVP_embedding

Deliberately NOT ported:
  - torch_cluster.knn_graph (we use offline cache from `knn_graph_custom`)
  - autoregressive branch (not needed for EXP005)
  - atom3d-specific classes

Bug fix vs upstream: `GVP_embedding.forward` used `if seq:` which raises
RuntimeError on multi-element tensors in modern PyTorch; changed to
`if seq is not None:`.
"""
import functools

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing
from torch_scatter import scatter_add


# ---------------------------------------------------------------------------
# Tuple (s, V) helpers
# ---------------------------------------------------------------------------
def tuple_sum(*args):
    """Sum any number of (s, V) tuples elementwise."""
    return tuple(map(sum, zip(*args)))


def tuple_cat(*args, dim=-1):
    """Concatenate any number of (s, V) tuples along `dim`.
    `dim=-1` on the scalar tensor is mapped to `dim=-2` on the vector tensor.
    """
    dim %= len(args[0][0].shape)
    s_args, v_args = list(zip(*args))
    return torch.cat(s_args, dim=dim), torch.cat(v_args, dim=dim)


def tuple_index(x, idx):
    """Index into a (s, V) tuple along the first dim."""
    return x[0][idx], x[1][idx]


def _norm_no_nan(x, axis=-1, keepdims=False, eps=1e-8, sqrt=True):
    """L2 norm clamped from below to avoid NaN gradients."""
    out = torch.clamp(torch.sum(torch.square(x), axis, keepdims), min=eps)
    return torch.sqrt(out) if sqrt else out


def _split(x, nv):
    """Split merged [s | V] tensor back into (s, V)."""
    v = torch.reshape(x[..., -3 * nv:], x.shape[:-1] + (nv, 3))
    s = x[..., :-3 * nv]
    return s, v


def _merge(s, v):
    """Merge (s, V) into single tensor along last dim."""
    v = torch.reshape(v, v.shape[:-2] + (3 * v.shape[-2],))
    return torch.cat([s, v], -1)


# ---------------------------------------------------------------------------
# GVP primitive
# ---------------------------------------------------------------------------
class GVP(nn.Module):
    """Geometric Vector Perceptron: (s_in, V_in) -> (s_out, V_out)."""

    def __init__(self, in_dims, out_dims, h_dim=None,
                 activations=(F.relu, torch.sigmoid), vector_gate=False):
        super().__init__()
        self.si, self.vi = in_dims
        self.so, self.vo = out_dims
        self.vector_gate = vector_gate
        if self.vi:
            self.h_dim = h_dim or max(self.vi, self.vo)
            self.wh = nn.Linear(self.vi, self.h_dim, bias=False)
            self.ws = nn.Linear(self.h_dim + self.si, self.so)
            if self.vo:
                self.wv = nn.Linear(self.h_dim, self.vo, bias=False)
                if self.vector_gate:
                    self.wsv = nn.Linear(self.so, self.vo)
        else:
            self.ws = nn.Linear(self.si, self.so)

        self.scalar_act, self.vector_act = activations
        # device-pinning dummy parameter (used only when creating zero vectors)
        self.dummy_param = nn.Parameter(torch.empty(0))

    def forward(self, x):
        if self.vi:
            s, v = x
            v = torch.transpose(v, -1, -2)
            vh = self.wh(v)
            vn = _norm_no_nan(vh, axis=-2)
            s = self.ws(torch.cat([s, vn], -1))
            if self.vo:
                v = self.wv(vh)
                v = torch.transpose(v, -1, -2)
                if self.vector_gate:
                    if self.vector_act:
                        gate = self.wsv(self.vector_act(s))
                    else:
                        gate = self.wsv(s)
                    v = v * torch.sigmoid(gate).unsqueeze(-1)
                elif self.vector_act:
                    v = v * self.vector_act(_norm_no_nan(v, axis=-1, keepdims=True))
        else:
            s = self.ws(x)
            if self.vo:
                v = torch.zeros(s.shape[0], self.vo, 3, device=self.dummy_param.device)
        if self.scalar_act:
            s = self.scalar_act(s)
        return (s, v) if self.vo else s


# ---------------------------------------------------------------------------
# Dropout / LayerNorm for (s, V)
# ---------------------------------------------------------------------------
class _VDropout(nn.Module):
    """Per-channel dropout for vector features: each vector channel is
    dropped as a whole (all 3 components together)."""

    def __init__(self, drop_rate):
        super().__init__()
        self.drop_rate = drop_rate
        self.dummy_param = nn.Parameter(torch.empty(0))

    def forward(self, x):
        if not self.training:
            return x
        device = self.dummy_param.device
        mask = torch.bernoulli(
            (1 - self.drop_rate) * torch.ones(x.shape[:-1], device=device)
        ).unsqueeze(-1)
        return mask * x / (1 - self.drop_rate)


class Dropout(nn.Module):
    """Joint dropout for (s, V) tuples or raw scalar tensors."""

    def __init__(self, drop_rate):
        super().__init__()
        self.sdropout = nn.Dropout(drop_rate)
        self.vdropout = _VDropout(drop_rate)

    def forward(self, x):
        if isinstance(x, torch.Tensor):
            return self.sdropout(x)
        s, v = x
        return self.sdropout(s), self.vdropout(v)


class LayerNorm(nn.Module):
    """LayerNorm for (s, V) tuples. Scalar part gets standard LayerNorm;
    vector part gets division by its RMS magnitude (SO(3)-equivariant)."""

    def __init__(self, dims):
        super().__init__()
        self.s, self.v = dims
        self.scalar_norm = nn.LayerNorm(self.s)

    def forward(self, x):
        if not self.v:
            return self.scalar_norm(x)
        s, v = x
        vn = _norm_no_nan(v, axis=-1, keepdims=True, sqrt=False)
        vn = torch.sqrt(torch.mean(vn, dim=-2, keepdim=True))
        return self.scalar_norm(s), v / vn


# ---------------------------------------------------------------------------
# GVP message passing
# ---------------------------------------------------------------------------
class GVPConv(MessagePassing):
    """Graph convolution (message passing) with GVP message function."""

    def __init__(self, in_dims, out_dims, edge_dims,
                 n_layers=3, module_list=None, aggr="mean",
                 activations=(F.relu, torch.sigmoid), vector_gate=False):
        super().__init__(aggr=aggr)
        self.si, self.vi = in_dims
        self.so, self.vo = out_dims
        self.se, self.ve = edge_dims

        GVP_ = functools.partial(GVP, activations=activations, vector_gate=vector_gate)

        module_list = module_list or []
        if not module_list:
            if n_layers == 1:
                module_list.append(
                    GVP_((2 * self.si + self.se, 2 * self.vi + self.ve),
                         (self.so, self.vo), activations=(None, None))
                )
            else:
                module_list.append(
                    GVP_((2 * self.si + self.se, 2 * self.vi + self.ve), out_dims)
                )
                for _ in range(n_layers - 2):
                    module_list.append(GVP_(out_dims, out_dims))
                module_list.append(GVP_(out_dims, out_dims, activations=(None, None)))
        self.message_func = nn.Sequential(*module_list)

    def forward(self, x, edge_index, edge_attr):
        x_s, x_v = x
        message = self.propagate(
            edge_index,
            s=x_s,
            v=x_v.reshape(x_v.shape[0], 3 * x_v.shape[1]),
            edge_attr=edge_attr,
        )
        return _split(message, self.vo)

    def message(self, s_i, v_i, s_j, v_j, edge_attr):
        v_j = v_j.view(v_j.shape[0], v_j.shape[1] // 3, 3)
        v_i = v_i.view(v_i.shape[0], v_i.shape[1] // 3, 3)
        message = tuple_cat((s_j, v_j), edge_attr, (s_i, v_i))
        message = self.message_func(message)
        return _merge(*message)


class GVPConvLayer(nn.Module):
    """Full GVP conv block: GVPConv + residual + LayerNorm + feedforward."""

    def __init__(self, node_dims, edge_dims,
                 n_message=3, n_feedforward=2, drop_rate=0.1,
                 autoregressive=False,
                 activations=(F.relu, torch.sigmoid), vector_gate=False):
        super().__init__()
        self.conv = GVPConv(
            node_dims, node_dims, edge_dims, n_message,
            aggr="add" if autoregressive else "mean",
            activations=activations, vector_gate=vector_gate,
        )
        GVP_ = functools.partial(GVP, activations=activations, vector_gate=vector_gate)
        self.norm = nn.ModuleList([LayerNorm(node_dims) for _ in range(2)])
        self.dropout = nn.ModuleList([Dropout(drop_rate) for _ in range(2)])

        ff_func = []
        if n_feedforward == 1:
            ff_func.append(GVP_(node_dims, node_dims, activations=(None, None)))
        else:
            hid_dims = 4 * node_dims[0], 2 * node_dims[1]
            ff_func.append(GVP_(node_dims, hid_dims))
            for _ in range(n_feedforward - 2):
                ff_func.append(GVP_(hid_dims, hid_dims))
            ff_func.append(GVP_(hid_dims, node_dims, activations=(None, None)))
        self.ff_func = nn.Sequential(*ff_func)

    def forward(self, x, edge_index, edge_attr,
                autoregressive_x=None, node_mask=None):
        """EXP005 only uses the non-autoregressive path; the autoregressive
        branch is preserved to keep the port faithful but will never trigger."""
        if autoregressive_x is not None:
            src, dst = edge_index
            mask = src < dst
            edge_index_forward = edge_index[:, mask]
            edge_index_backward = edge_index[:, ~mask]
            edge_attr_forward = tuple_index(edge_attr, mask)
            edge_attr_backward = tuple_index(edge_attr, ~mask)
            dh = tuple_sum(
                self.conv(x, edge_index_forward, edge_attr_forward),
                self.conv(autoregressive_x, edge_index_backward, edge_attr_backward),
            )
            count = scatter_add(
                torch.ones_like(dst), dst, dim_size=dh[0].size(0)
            ).clamp(min=1).unsqueeze(-1)
            dh = dh[0] / count, dh[1] / count.unsqueeze(-1)
        else:
            dh = self.conv(x, edge_index, edge_attr)

        if node_mask is not None:
            x_ = x
            x, dh = tuple_index(x, node_mask), tuple_index(dh, node_mask)

        x = self.norm[0](tuple_sum(x, self.dropout[0](dh)))
        dh = self.ff_func(x)
        x = self.norm[1](tuple_sum(x, self.dropout[1](dh)))

        if node_mask is not None:
            x_[0][node_mask], x_[1][node_mask] = x[0], x[1]
            x = x_
        return x


# ---------------------------------------------------------------------------
# Top-level wrapper
# ---------------------------------------------------------------------------
class GVP_embedding(nn.Module):
    """Residue-graph encoder. Input: precomputed node/edge features +
    optional integer seq. Output: per-node scalar embedding [N, 2 * n_scalar].

    For EXP005 we use node_h_dim=(128, 16), so the output is [N, 256].
    """

    def __init__(self, node_in_dim, node_h_dim,
                 edge_in_dim, edge_h_dim,
                 seq_in=False, num_layers=3, drop_rate=0.1):
        super().__init__()

        if seq_in:
            self.W_s = nn.Embedding(20, 20)
            node_in_dim = (node_in_dim[0] + 20, node_in_dim[1])
        else:
            self.W_s = None

        self.W_v = nn.Sequential(
            LayerNorm(node_in_dim),
            GVP(node_in_dim, node_h_dim, activations=(None, None)),
        )
        self.W_e = nn.Sequential(
            LayerNorm(edge_in_dim),
            GVP(edge_in_dim, edge_h_dim, activations=(None, None)),
        )

        self.layers = nn.ModuleList(
            GVPConvLayer(node_h_dim, edge_h_dim, drop_rate=drop_rate)
            for _ in range(num_layers)
        )

        ns, _ = node_h_dim
        self.W_out = nn.Sequential(
            LayerNorm(node_h_dim),
            GVP(node_h_dim, (2 * ns, 0)),
        )

    def forward(self, h_V, edge_index, h_E, seq=None):
        """
        :param h_V: tuple (s, V) of node embeddings
                    s: [N, node_in_scalar] (float), V: [N, node_in_vector, 3] (float)
        :param edge_index: [2, E] long (or int — will be cast to long)
        :param h_E: tuple (s, V) of edge embeddings
                    s: [E, edge_in_scalar] (float), V: [E, edge_in_vector, 3] (float)
        :param seq: optional tensor [N] of AA indices (0..19) used when
                    `seq_in=True`; will be cast to long for nn.Embedding.
                    BUG FIX vs upstream: upstream used `if seq:` which fails
                    on multi-element tensors.

        DTYPE POLICY: the GVP cache stores node_s/edge_s as fp16 and
        node_v/edge_v as fp32, but Linear/LayerNorm layers here are fp32.
        All four geometric inputs are cast to float32 at entry. `seq` and
        `edge_index` are cast to long.
        """
        # Normalize dtypes at the boundary (cache -> model).
        h_V = (h_V[0].float(), h_V[1].float())
        h_E = (h_E[0].float(), h_E[1].float())
        if edge_index.dtype != torch.long:
            edge_index = edge_index.long()

        if seq is not None and self.W_s is not None:
            seq_emb = self.W_s(seq.long())
            h_V = (torch.cat([h_V[0], seq_emb], dim=-1), h_V[1])
        h_V = self.W_v(h_V)
        h_E = self.W_e(h_E)
        for layer in self.layers:
            h_V = layer(h_V, edge_index, h_E)
        out = self.W_out(h_V)
        return out  # [N, 2 * ns] pure scalar tensor
