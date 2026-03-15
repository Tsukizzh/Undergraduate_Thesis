"""
Step 10: PyG Dataset reading from .pt shard cache (ezspec_pt_v1/).

Codex-reviewed (3 rounds). Fixes applied:
- Shard field names match build_pt_cache.py exactly
- Global enzyme/substrate IDs (dataset_id * 10M + local_idx)
- ligand_x feature order matches FeaturizeLigandAtom exactly
- Per-entity tensor cache (enzyme/substrate decoded+padded)
- Pre-built vocab tensors (no per-call allocation)
- Reuses single GaussianSmearing instance
- Worker-safe via __getstate__
"""
from __future__ import annotations

import os
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

# Ensure src/ is importable
_SRC_CANDIDATES = [
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "src"),
    os.path.join(os.path.dirname(__file__), "..", "..", "src"),
    os.path.join(os.path.dirname(__file__), "..", "src"),
]
for _p in _SRC_CANDIDATES:
    if os.path.isdir(_p) and os.path.abspath(_p) not in sys.path:
        sys.path.insert(0, os.path.abspath(_p))
        break

from Datasets.data_representer import StructureSequenceData
from Datasets.Structure.protein_ligand import ATOM_FEATS
from Models.common import GaussianSmearing

# ---------------------------------------------------------------------------
# Constants — derived from source code, not hardcoded
# ---------------------------------------------------------------------------

PROTEIN_ELEMENTS = [1, 6, 7, 8, 16, 34]
PROTEIN_NUM_AA = 21
PROTEIN_FEAT_DIM = len(PROTEIN_ELEMENTS) + PROTEIN_NUM_AA + 1  # 28

LIGAND_ELEMENTS = [1, 6, 7, 8, 9, 15, 16, 17, 26, 35, 53]
LIGAND_DEGREE_CLASSES = int(ATOM_FEATS["Degree"])
LIGAND_NUMHS_CLASSES = int(ATOM_FEATS["NumHs"])
LIGAND_HYBRID_CLASSES = int(ATOM_FEATS["Hybridization"])
LIGAND_FEAT_DIM = len(LIGAND_ELEMENTS) + sum(int(v) for v in ATOM_FEATS.values())

BOND_NUM_CLASSES = 6
KNN_BOND_CLASS = 5
GLOBAL_ID_STRIDE = 10_000_000

DATASET_TAGS = ["0", "1", "2", "3", "4", "5", "6"]
STR_TAG_VOCAB = ["complex", "ligand"]

MAX_ENZYME_LEN = 1450
MAX_SUBSTRATE_LEN = 280
LIGAND_INDEX_SENTINEL = 280

# Pre-built vocab tensors (avoid per-call allocation)
PROTEIN_ELEMENTS_T = torch.tensor(PROTEIN_ELEMENTS, dtype=torch.long)
LIGAND_ELEMENTS_T = torch.tensor(LIGAND_ELEMENTS, dtype=torch.long)


# ---------------------------------------------------------------------------
# Caches
# ---------------------------------------------------------------------------

class ShardCache:
    """LRU cache for loaded .pt shards."""
    def __init__(self, max_size: int = 8):
        self._cache: OrderedDict[str, dict] = OrderedDict()
        self._max_size = max_size

    def get(self, path: str) -> dict:
        if path in self._cache:
            self._cache.move_to_end(path)
            return self._cache[path]
        data = torch.load(path, map_location="cpu", weights_only=False)
        self._cache[path] = data
        if len(self._cache) > self._max_size:
            self._cache.popitem(last=False)
        return data

    def clear(self):
        self._cache.clear()


class EntityTensorCache:
    """LRU cache for decoded+padded enzyme/substrate tensors."""
    def __init__(self, max_size: int = 4096):
        self._cache: OrderedDict[int, object] = OrderedDict()
        self._max_size = max_size

    def get(self, key: int):
        if key not in self._cache:
            return None
        self._cache.move_to_end(key)
        return self._cache[key]

    def put(self, key: int, value) -> None:
        self._cache[key] = value
        self._cache.move_to_end(key)
        if len(self._cache) > self._max_size:
            self._cache.popitem(last=False)


# ---------------------------------------------------------------------------
# Vectorized feature reconstruction
# ---------------------------------------------------------------------------

def _one_hot_from_values(values: torch.Tensor, vocab_t: torch.Tensor) -> torch.Tensor:
    """Vectorized one-hot using pre-built vocab tensor."""
    return (values.unsqueeze(-1) == vocab_t.unsqueeze(0)).float()


def _one_hot_clamp(values: torch.Tensor, num_classes: int) -> torch.Tensor:
    return F.one_hot(values.clamp(0, num_classes - 1).long(), num_classes).float()


def rebuild_protein_x(
    protein_element: torch.Tensor,
    protein_aa_type: torch.Tensor,
    protein_is_backbone: torch.Tensor,
    n_lig: int,
) -> torch.Tensor:
    elem_oh = _one_hot_from_values(protein_element.long(), PROTEIN_ELEMENTS_T)
    aa_oh = _one_hot_clamp(protein_aa_type, PROTEIN_NUM_AA)
    bb = protein_is_backbone.float().unsqueeze(-1)
    prot_feat = torch.cat([elem_oh, aa_oh, bb], dim=-1)
    return torch.cat([torch.zeros(n_lig, PROTEIN_FEAT_DIM), prot_feat], dim=0)


def rebuild_ligand_x(
    ligand_element: torch.Tensor,
    ligand_atom_aux: torch.Tensor,
    n_prot: int,
) -> torch.Tensor:
    """Exact order from FeaturizeLigandAtom:
    [element_oh(11), atomic_number/100(1), aromatic(1), degree_oh(6), numhs_oh(6), hybrid_oh(H)]
    """
    elem_oh = _one_hot_from_values(ligand_element.long(), LIGAND_ELEMENTS_T)
    atomic_scaled = ligand_element.float().unsqueeze(-1) / 100.0
    aromatic = ligand_atom_aux[:, 0].float().unsqueeze(-1)
    degree_oh = _one_hot_clamp(ligand_atom_aux[:, 1], LIGAND_DEGREE_CLASSES)
    numhs_oh = _one_hot_clamp(ligand_atom_aux[:, 2], LIGAND_NUMHS_CLASSES)
    hybrid_oh = _one_hot_clamp(ligand_atom_aux[:, 3], LIGAND_HYBRID_CLASSES)
    lig_feat = torch.cat(
        [elem_oh, atomic_scaled, aromatic, degree_oh, numhs_oh, hybrid_oh],
        dim=-1,
    )
    return torch.cat([lig_feat, torch.zeros(n_prot, LIGAND_FEAT_DIM)], dim=0)


def rebuild_edge_features(
    ligand_pos: torch.Tensor,
    protein_pos: torch.Tensor,
    bond_index: torch.Tensor,
    bond_type: torch.Tensor,
    knn_edge_index: torch.Tensor,
    n_lig: int,
    edge_mode: str,
    dist_noise: bool,
    smear: GaussianSmearing,
) -> tuple[torch.Tensor, torch.Tensor]:
    pos = torch.cat([ligand_pos, protein_pos], dim=0)
    num_r_gaussian = int(smear.offset.numel())

    def _dist(ei):
        if ei.numel() == 0:
            return torch.zeros(0, dtype=torch.float32)
        return torch.norm(pos[ei[1]] - pos[ei[0]], p=2, dim=-1)

    dist_knn = _dist(knn_edge_index)
    dist_bond = _dist(bond_index)

    if dist_noise:
        if dist_knn.numel() > 0:
            dist_knn = dist_knn + torch.from_numpy(
                np.random.laplace(0.001994, 0.031939, size=(dist_knn.shape[0],))
            ).float()
        if dist_bond.numel() > 0:
            dist_bond = dist_bond + torch.from_numpy(
                np.random.laplace(0.001994, 0.031939, size=(dist_bond.shape[0],))
            ).float()

    empty_dist = torch.zeros(0, num_r_gaussian)
    dist_feat_knn = smear(dist_knn) if dist_knn.numel() > 0 else empty_dist
    dist_feat_bond = smear(dist_bond) if dist_bond.numel() > 0 else empty_dist

    bt = bond_type.clone().long()
    bt[bt == 12] = 5
    bond_oh = F.one_hot((bt - 1).clamp(0, BOND_NUM_CLASSES - 1), BOND_NUM_CLASSES).float()
    knn_oh = F.one_hot(
        torch.full((knn_edge_index.shape[1],), KNN_BOND_CLASS, dtype=torch.long),
        BOND_NUM_CLASSES,
    ).float()

    if knn_edge_index.numel() > 0:
        cross_knn = (
            ((knn_edge_index[0] < n_lig) & (knn_edge_index[1] >= n_lig)) |
            ((knn_edge_index[0] >= n_lig) & (knn_edge_index[1] < n_lig))
        ).int().unsqueeze(-1)
    else:
        cross_knn = torch.zeros(0, 1, dtype=torch.int)
    cross_bond = torch.zeros(bond_index.shape[1], 1, dtype=torch.int)

    bond_attr = torch.cat([bond_oh, dist_feat_bond, cross_bond], dim=-1)
    knn_attr = torch.cat([knn_oh, dist_feat_knn, cross_knn], dim=-1)

    complex_edge_index = torch.cat([bond_index, knn_edge_index], dim=-1)
    if edge_mode == "legacy_bug":
        complex_edge_attr = torch.cat([knn_attr, bond_attr], dim=0)
    elif edge_mode == "fixed":
        complex_edge_attr = torch.cat([bond_attr, knn_attr], dim=0)
    else:
        raise ValueError(f"Unknown edge_mode: {edge_mode}")

    return complex_edge_index.long(), complex_edge_attr.float()


# ---------------------------------------------------------------------------
# Main Dataset
# ---------------------------------------------------------------------------

class PtCacheDataset(torch.utils.data.Dataset):
    """PyG-compatible dataset reading from .pt shard cache."""

    def __init__(
        self,
        cache_dir: str | Path,
        split: str,
        edge_mode: str = "fixed",
        dist_noise: bool = False,
        cutoff: float = 10.0,
        num_r_gaussian: int = 32,
        max_enzyme_len: int = MAX_ENZYME_LEN,
        max_substrate_len: int = MAX_SUBSTRATE_LEN,
        graph_cache_size: int = 4,
        enzyme_cache_size: int = 8,
        substrate_cache_size: int = 8,
    ):
        super().__init__()
        self.cache_dir = Path(cache_dir)
        self.split = split
        self.edge_mode = edge_mode
        self.dist_noise = dist_noise
        self.max_enzyme_len = max_enzyme_len
        self.max_substrate_len = max_substrate_len

        # Split index (always in memory, small)
        self._index = torch.load(
            self.cache_dir / split / "index.pt", map_location="cpu", weights_only=False
        )
        self._n_samples = len(self._index["sample_ids"])

        # Enzyme & substrate indices
        enz_idx = torch.load(
            self.cache_dir / "enzymes" / "index.pt", map_location="cpu", weights_only=False
        )
        sub_idx = torch.load(
            self.cache_dir / "substrates" / "index.pt", map_location="cpu", weights_only=False
        )

        # Build lookup dicts: global_id → (shard_id, row_id)
        self._enz_lookup = {
            int(enz_idx["enzyme_ids"][i]): (int(enz_idx["shard_ids"][i]), int(enz_idx["row_ids"][i]))
            for i in range(len(enz_idx["enzyme_ids"]))
        }
        self._sub_lookup = {
            int(sub_idx["substrate_ids"][i]): (int(sub_idx["shard_ids"][i]), int(sub_idx["row_ids"][i]))
            for i in range(len(sub_idx["substrate_ids"]))
        }

        # Caches (per-worker, cleared on spawn)
        self._graph_cache = ShardCache(graph_cache_size)
        self._enz_cache = ShardCache(enzyme_cache_size)
        self._sub_cache = ShardCache(substrate_cache_size)
        self._enzyme_tensor_cache = EntityTensorCache(4096)
        self._substrate_tensor_cache = EntityTensorCache(4096)

        # Reusable GaussianSmearing (stateless)
        self._smear = GaussianSmearing(stop=cutoff, num_gaussians=num_r_gaussian)

    def __len__(self) -> int:
        return self._n_samples

    def __getstate__(self):
        """Clear caches before pickling (Windows spawn mode)."""
        state = self.__dict__.copy()
        state["_graph_cache"] = ShardCache(state["_graph_cache"]._max_size)
        state["_enz_cache"] = ShardCache(state["_enz_cache"]._max_size)
        state["_sub_cache"] = ShardCache(state["_sub_cache"]._max_size)
        state["_enzyme_tensor_cache"] = EntityTensorCache(4096)
        state["_substrate_tensor_cache"] = EntityTensorCache(4096)
        return state

    def _load_graph_sample(self, idx: int) -> dict:
        shard_id = int(self._index["graph_shards"][idx])
        row_id = int(self._index["graph_rows"][idx])
        shard_path = str(self.cache_dir / self.split / f"graph_{shard_id:04d}.pt")
        shard = self._graph_cache.get(shard_path)

        def _b(name):
            p = shard[name]
            return int(p[row_id]), int(p[row_id + 1])

        s_l, e_l = _b("ligand_ptr")
        s_p, e_p = _b("protein_ptr")
        s_b, e_b = _b("bond_ptr")
        s_k, e_k = _b("knn_ptr")

        ds_id = int(shard["dataset_ids"][row_id])
        local_enz = int(shard["enzyme_ids"][row_id])
        local_sub = int(shard["substrate_ids"][row_id])

        return {
            "ligand_pos": shard["ligand_pos"][s_l:e_l].float(),
            "protein_pos": shard["protein_pos"][s_p:e_p].float(),
            "ligand_element": shard["ligand_element"][s_l:e_l].long(),
            "protein_element": shard["protein_element"][s_p:e_p].long(),
            "protein_aa_type": shard["protein_aa_type"][s_p:e_p].long(),
            "protein_is_backbone": shard["protein_is_backbone"][s_p:e_p].long(),
            "ligand_atom_aux": shard["ligand_atom_aux"][s_l:e_l].long(),
            "ligand_index_raw": shard["ligand_index_raw"][s_l:e_l].long(),
            "bond_index": shard["bond_index"][:, s_b:e_b].long(),
            "bond_type": shard["bond_type"][s_b:e_b].long(),
            "knn_edge_index": shard["knn_edge_index"][:, s_k:e_k].long(),
            "enzyme_global_id": local_enz,  # already global in shards
            "substrate_global_id": local_sub,  # already global in shards
            "dataset_id": ds_id,
            "label": int(shard["labels"][row_id]),
            "str_tag_code": int(shard["str_tag_codes"][row_id]),
            "sample_weight": float(shard["sample_weights"][row_id]),
        }

    def _load_enzyme_embedding(self, enzyme_global_id: int) -> tuple[torch.Tensor, torch.Tensor]:
        cached = self._enzyme_tensor_cache.get(enzyme_global_id)
        if cached is not None:
            return cached

        if enzyme_global_id not in self._enz_lookup:
            out = (
                torch.zeros(self.max_enzyme_len, 1280, dtype=torch.float32),
                torch.ones(1, self.max_enzyme_len, dtype=torch.bool),
            )
            self._enzyme_tensor_cache.put(enzyme_global_id, out)
            return out

        shard_id, row_id = self._enz_lookup[enzyme_global_id]
        shard_path = str(self.cache_dir / "enzymes" / f"esm_{shard_id:04d}.pt")
        shard = self._enz_cache.get(shard_path)

        offsets = shard["offsets"]
        start, end = int(offsets[row_id]), int(offsets[row_id + 1])
        raw_emb = shard["embedding"][start:end].float()

        length = min(raw_emb.shape[0], self.max_enzyme_len)
        emb = torch.zeros(self.max_enzyme_len, 1280, dtype=torch.float32)
        emb[:length] = raw_emb[:length]

        mask = torch.zeros(1, self.max_enzyme_len, dtype=torch.bool)
        mask[0, length:] = True

        out = (emb, mask)
        self._enzyme_tensor_cache.put(enzyme_global_id, out)
        return out

    def _load_substrate_features(self, substrate_global_id: int) -> dict:
        cached = self._substrate_tensor_cache.get(substrate_global_id)
        if cached is not None:
            return cached

        if substrate_global_id not in self._sub_lookup:
            out = {
                "grover": torch.zeros(self.max_substrate_len, 2400, dtype=torch.float32),
                "reaction_padding_mask": torch.ones(1, self.max_substrate_len, dtype=torch.bool),
                "grover_mean": torch.zeros(1, 4885, dtype=torch.float32),
                "morgan": torch.zeros(1, 1024, dtype=torch.float32),
            }
            self._substrate_tensor_cache.put(substrate_global_id, out)
            return out

        shard_id, row_id = self._sub_lookup[substrate_global_id]
        shard_path = str(self.cache_dir / "substrates" / f"grover_{shard_id:04d}.pt")
        shard = self._sub_cache.get(shard_path)

        atom_offsets = shard["atom_offsets"]
        a_start, a_end = int(atom_offsets[row_id]), int(atom_offsets[row_id + 1])
        raw_grover = shard["grover_atom"][a_start:a_end].float()

        atom_len = min(raw_grover.shape[0], self.max_substrate_len)
        grover = torch.zeros(self.max_substrate_len, 2400, dtype=torch.float32)
        grover[:atom_len] = raw_grover[:atom_len]

        r_mask = torch.zeros(1, self.max_substrate_len, dtype=torch.bool)
        r_mask[0, atom_len:] = True

        out = {
            "grover": grover,
            "reaction_padding_mask": r_mask,
            "grover_mean": shard["grover_mean"][row_id:row_id + 1].float(),
            "morgan": shard["morgan"][row_id:row_id + 1].float(),
        }
        self._substrate_tensor_cache.put(substrate_global_id, out)
        return out

    def __getitem__(self, idx: int) -> StructureSequenceData:
        g = self._load_graph_sample(idx)

        n_lig = g["ligand_pos"].shape[0]
        n_prot = g["protein_pos"].shape[0]

        protein_x = rebuild_protein_x(
            g["protein_element"], g["protein_aa_type"], g["protein_is_backbone"], n_lig
        )
        ligand_x = rebuild_ligand_x(g["ligand_element"], g["ligand_atom_aux"], n_prot)

        ligand_mask = torch.cat([torch.ones(n_lig), torch.zeros(n_prot)])
        protein_mask = torch.cat([torch.zeros(n_lig), torch.ones(n_prot)])

        ligand_index = torch.cat([
            g["ligand_index_raw"],
            torch.full((n_prot,), LIGAND_INDEX_SENTINEL, dtype=torch.long),
        ])

        complex_edge_index, complex_edge_attr = rebuild_edge_features(
            ligand_pos=g["ligand_pos"],
            protein_pos=g["protein_pos"],
            bond_index=g["bond_index"],
            bond_type=g["bond_type"],
            knn_edge_index=g["knn_edge_index"],
            n_lig=n_lig,
            edge_mode=self.edge_mode,
            dist_noise=self.dist_noise,
            smear=self._smear,
        )

        embedding, enzyme_padding_mask = self._load_enzyme_embedding(g["enzyme_global_id"])
        sub = self._load_substrate_features(g["substrate_global_id"])

        return StructureSequenceData(
            embedding=embedding,
            enzyme_padding_mask=enzyme_padding_mask,
            grover=sub["grover"],
            reaction_padding_mask=sub["reaction_padding_mask"],
            grover_mean=sub["grover_mean"],
            morgan=sub["morgan"],
            protein_x=protein_x,
            ligand_x=ligand_x,
            protein_mask=protein_mask,
            ligand_mask=ligand_mask,
            ligand_index=ligand_index,
            complex_edge_index=complex_edge_index,
            complex_edge_attr=complex_edge_attr,
            tag=DATASET_TAGS[g["dataset_id"]],
            str_tag=STR_TAG_VOCAB[g["str_tag_code"]],
            label=torch.tensor(g["label"], dtype=torch.float64),
            sample_weight=torch.tensor(g["sample_weight"], dtype=torch.float32),
        )

    @property
    def valid_idx(self) -> list[int]:
        """Original dataframe row IDs for prediction alignment."""
        return self._index["sample_ids"].tolist()


# ---------------------------------------------------------------------------
# Convenience factory
# ---------------------------------------------------------------------------

def build_pt_dataloaders(
    cache_dir: str | Path,
    batch_size: int = 14,
    edge_mode: str = "fixed",
    dist_noise_train: bool = True,
    cutoff: float = 10.0,
    num_r_gaussian: int = 32,
    num_workers: int = 4,
    pin_memory: bool = True,
    prefetch_factor: int = 4,
) -> dict[str, torch.utils.data.DataLoader]:
    from torch_geometric.loader import DataLoader

    loaders = {}
    for split in ["train", "val", "test"]:
        split_dir = Path(cache_dir) / split
        if not split_dir.exists():
            continue

        dataset = PtCacheDataset(
            cache_dir=cache_dir,
            split=split,
            edge_mode=edge_mode,
            dist_noise=(dist_noise_train if split == "train" else False),
            cutoff=cutoff,
            num_r_gaussian=num_r_gaussian,
        )

        kw = dict(
            batch_size=batch_size,
            shuffle=(split == "train"),
            num_workers=num_workers,
            pin_memory=pin_memory,
            follow_batch=["ligand_index"],
            persistent_workers=(num_workers > 0),
        )
        if num_workers > 0:
            kw["prefetch_factor"] = prefetch_factor

        loaders[split] = DataLoader(dataset, **kw)

    return loaders
