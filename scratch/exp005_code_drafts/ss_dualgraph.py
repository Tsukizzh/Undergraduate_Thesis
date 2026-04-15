"""
EXP005 dualgraph 2+ model: SS + residue-level GVP channel.

Two fusion paths (option 2+):
  1. h_res -> x_pro injection  (deep fusion, before cross-attention)
     For every pocket residue with pocket_residue_idx >= 0 and gvp_valid=1,
     add `alpha * h_res_proj(h_res)` into `x_pro` at the corresponding
     (sample, uniprot_pos) slot. Duplicates (e.g. homodimer chains mapping
     to the same UniProt position) are averaged first via scatter_mean.
  2. g_res bypass             (shallow bypass, 8th vector in final concat)
     `g_res = scatter_mean(h_res, gvp_batch)` per sample, gated by gvp_valid,
     projected to hidden_dim, appended as the 8th item of `embeddings`.

Baseline preservation at step 0:
  - `h_res_proj[-1]` (Linear 256->hidden_dim) is weight-and-bias zero-init
    so the injection contributes exactly zero regardless of `alpha`.
  - `specificity_header.mlp.net[0]`'s new 128-dim block (for g_res) is
    weight-zero-init so the bypass contributes exactly zero at step 0.
  - Net effect: SSDualgraph at step 0 computes the identical logits as the
    baseline SS, while gradients flow into the GVP branch immediately and
    let the model start using it when it becomes informative.

Expected batch fields (from PtCacheDualgraphDataset / DualgraphData):
  # baseline fields (unchanged):
  G.embedding, G.grover, G.grover_mean, G.morgan,
  G.enzyme_padding_mask, G.reaction_padding_mask,
  G.tag, G.str_tag, G.label, G.sample_weight, plus the atom graph fields.
  # new GVP fields:
  G.gvp_node_s            [sum_N, 6]        float32 (cast inside gvp encoder)
  G.gvp_node_v            [sum_N, 3, 3]     float32
  G.gvp_edge_index        [2, sum_E]        long (offset by node count)
  G.gvp_edge_s            [sum_E, 32]       float32
  G.gvp_edge_v            [sum_E, 1, 3]     float32
  G.gvp_aa_type           [sum_N]           long
  G.gvp_pocket_residue_idx[sum_N]           long, 0..seq_len-1 (or -1 sentinel)
  G.gvp_node_s_batch      [sum_N]           long, sample id in [0, B)
                                            (auto-created by DataLoader when
                                             `follow_batch=['gvp_node_s']`)
  G.gvp_valid             [B]               bool, False for placeholder samples
"""
import torch
import torch.nn as nn
from torch_scatter import scatter_mean

from Models.ss import SS
from Models.Structure.gvp import GVP_embedding

# Max enzyme sequence length, matches SS.forward view size.
MAX_SEQ_LEN = 1450


class SSDualgraph(SS):
    def __init__(self, config):
        super().__init__(config)
        gvp_cfg = config.model.gvp

        # ---- GVP encoder ----
        node_h_dim = (
            int(gvp_cfg.get("node_h_scalar", 128)),
            int(gvp_cfg.get("node_h_vector", 16)),
        )
        edge_h_dim = (
            int(gvp_cfg.get("edge_h_scalar", 32)),
            int(gvp_cfg.get("edge_h_vector", 1)),
        )
        self.gvp_encoder = GVP_embedding(
            node_in_dim=(6, 3),        # 6 dihedral cos/sin + 3 orientation/sidechain vectors
            node_h_dim=node_h_dim,     # (128, 16)
            edge_in_dim=(32, 1),       # 32 (RBF 16 + positional 16) + 1 CA-CA unit vector
            edge_h_dim=edge_h_dim,     # (32, 1)
            seq_in=bool(gvp_cfg.get("seq_in", True)),
            num_layers=int(gvp_cfg.get("num_layers", 3)),
            drop_rate=float(gvp_cfg.get("drop_rate", 0.1)),
        )

        # GVP_embedding output is [N, 2 * node_h_dim[0]] = [N, 256]
        gvp_out_dim = 2 * node_h_dim[0]

        # ---- h_res projection (for injection into x_pro) ----
        # LayerNorm + Linear; the Linear is zero-initialized so injection
        # starts at exactly zero regardless of the alpha gate.
        self.h_res_proj = nn.Sequential(
            nn.LayerNorm(gvp_out_dim),
            nn.Linear(gvp_out_dim, self.hidden_dim),
        )
        nn.init.zeros_(self.h_res_proj[-1].weight)
        nn.init.zeros_(self.h_res_proj[-1].bias)

        # Learnable alpha gate on injection magnitude (starts ~0.10).
        # Kept separately even with zero-init projection so the gate becomes
        # meaningful once the projection learns non-zero weights.
        self.h_res_gate_logit = nn.Parameter(
            torch.tensor(float(gvp_cfg.get("alpha_logit_init", -2.2)))
        )

        # ---- g_res projection (for 8th embedding vector) ----
        self.g_res_proj = nn.Linear(gvp_out_dim, self.hidden_dim)

        # ---- Expand specificity_header input 7H -> 8H ----
        # The new 8th block corresponds to g_res. Zero-init the new column
        # slice so the bypass contributes exactly zero at step 0.
        old_first = self.specificity_header.mlp.net[0]
        assert isinstance(old_first, nn.Linear), \
            f"expected specificity_header.mlp.net[0] to be Linear, got {type(old_first)}"
        old_in = old_first.in_features
        new_in = old_in + self.hidden_dim
        new_first = nn.Linear(new_in, old_first.out_features)
        with torch.no_grad():
            new_first.weight[:, :old_in].copy_(old_first.weight)
            new_first.weight[:, old_in:].zero_()
            new_first.bias.copy_(old_first.bias)
        self.specificity_header.mlp.net[0] = new_first
        self._gvp_head_offset = old_in  # cached for debugging

    # -------------------------------------------------------------------
    # Forward (reimplemented from SS.forward with GVP insertions)
    # -------------------------------------------------------------------
    def forward(self, G):
        # ---------- Baseline step 1: protein MLP + reshape ----------
        x_pro = self.protein_mlp(G.embedding)
        x_pro = x_pro.view(-1, MAX_SEQ_LEN, self.hidden_dim)
        B = x_pro.shape[0]

        # ---------- EXP005 STEP A: GVP encoder ----------
        # h_res: [sum_N, 256] (pure scalar output)
        h_res = self.gvp_encoder(
            (G.gvp_node_s, G.gvp_node_v),
            G.gvp_edge_index,
            (G.gvp_edge_s, G.gvp_edge_v),
            seq=G.gvp_aa_type,
        )

        # ---------- EXP005 STEP B: h_res -> x_pro injection ----------
        # Only inject nodes where (a) pocket_residue_idx >= 0 AND
        # (b) the sample is gvp_valid. Aggregate duplicate (sample, seq_idx)
        # writes via scatter_mean so homodimer multiplicity does not inflate
        # the signal.
        # NOTE: PyG auto-creates `gvp_node_s_batch` because the DataLoader is
        # set up with `follow_batch=['gvp_node_s']`.
        gvp_batch = G.gvp_node_s_batch.long()
        prix = G.gvp_pocket_residue_idx.long()
        valid_sample = G.gvp_valid.to(torch.bool)  # [B] bool
        valid_node = (prix >= 0) & valid_sample[gvp_batch]

        if valid_node.any():
            h_res_proj_valid = self.h_res_proj(h_res[valid_node])  # [M, H]
            sample_idx = gvp_batch[valid_node]
            seq_idx = prix[valid_node]
            # Build a unique sample*MAX_SEQ_LEN + seq_idx key for duplicates
            flat_key = sample_idx * MAX_SEQ_LEN + seq_idx
            unique_keys, inverse = flat_key.unique(sorted=True, return_inverse=True)
            h_agg = scatter_mean(
                h_res_proj_valid, inverse, dim=0, dim_size=unique_keys.numel()
            )
            sid_unique = torch.div(unique_keys, MAX_SEQ_LEN, rounding_mode="floor")
            rix_unique = unique_keys.remainder(MAX_SEQ_LEN)
            alpha = torch.sigmoid(self.h_res_gate_logit)
            # Explicit dtype cast avoids AMP ambiguity when x_pro is fp16
            # autocast and h_res_proj stays fp32.
            delta = (alpha * h_agg).to(x_pro.dtype)
            # Advanced indexing into [B, 1450, H]
            x_pro[sid_unique, rix_unique] = x_pro[sid_unique, rix_unique] + delta
        # else: all nodes are placeholders, injection is a no-op

        # ---------- Baseline step 2: atom features + EGNN ----------
        atom_features = []
        if self.config.model.use_gnn:
            str_mean, x_str = self._get_graph_output(G, x_pro.shape[0])
            atom_features.append(x_str)
        if "grover" in self.config.data.atom_features:
            grover = G.grover.view(-1, 280, int(self.config.data.atom_features[1]))
            grover = self.atom_feature_mlp_dict["grover"](grover)
            atom_features.append(grover)
        x_reaction = self.atom_feature_aggregator(atom_features)

        # ---------- Baseline step 3: cross-attention ----------
        if self.use_attention:
            weighted_sum_pro, _ = self.enzyme_attention(
                x_pro, x_reaction, x_reaction,
                need_weights=True,
                key_padding_mask=G.reaction_padding_mask,
            )
            weighted_sum_reaction, _ = self.reaction_attention(
                x_reaction, x_pro, x_pro,
                need_weights=True,
                key_padding_mask=G.enzyme_padding_mask,
            )

        # ---------- Baseline step 4: masked pooling ----------
        reaction_mask = (~G.reaction_padding_mask).unsqueeze(-1).float()
        enzyme_mask = (~G.enzyme_padding_mask).unsqueeze(-1).float()
        x_reaction = (x_reaction * reaction_mask).sum(dim=1) / reaction_mask.sum(dim=(1, 2)).unsqueeze(-1)
        x_pro = (x_pro * enzyme_mask).sum(dim=1) / enzyme_mask.sum(dim=(1, 2)).unsqueeze(-1)
        if self.use_attention:
            weighted_sum_reaction = (weighted_sum_reaction * reaction_mask).sum(dim=1) / reaction_mask.sum(dim=(1, 2)).unsqueeze(-1)
            weighted_sum_pro = (weighted_sum_pro * enzyme_mask).sum(dim=1) / enzyme_mask.sum(dim=(1, 2)).unsqueeze(-1)

        # ---------- Baseline step 5: embedding list ----------
        embeddings = [x_pro, x_reaction]
        if self.use_attention:
            embeddings.extend([weighted_sum_pro, weighted_sum_reaction])
        if self.config.model.use_gnn:
            embeddings.extend([str_mean])
        if "grover_mean" in self.config.data.features:
            embeddings.append(self.feature_mlp_dict["grover_mean"](G.grover_mean))
        if "morgan" in self.config.data.features:
            embeddings.append(self.feature_mlp_dict["morgan"](G.morgan))

        # ---------- EXP005 STEP C: g_res bypass ----------
        # Per-sample mean of h_res, gated by gvp_valid, projected to H.
        g_res = scatter_mean(h_res, gvp_batch, dim=0, dim_size=B)
        g_res = g_res * G.gvp_valid.to(g_res.dtype).unsqueeze(-1)
        g_res_out = self.g_res_proj(g_res)
        embeddings.append(g_res_out)  # 8th vector

        # ---------- Baseline step 6: specificity header ----------
        specificity_output = self.specificity_header(embeddings)
        return specificity_output, [G.tag, G.str_tag]
