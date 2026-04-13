"""
Structure cache utilities for EZSpecificity training optimization.

Eliminates the CPU bottleneck (k-NN graph construction, ~50-70% of per-sample time)
by pre-computing and caching structure features in LMDB. Sequence features remain live.

Also FIXES a bug in the original EdgeConnection (transforms.py:130-147) where
complex_edge_attr and complex_edge_index had mismatched edge ordering.

Components:
    BuildStructureCacheData  - offline transform: raw structure → cache payload
    RebuildComplexEdgeAttr   - runtime transform: cache payload → model-ready tensors
    CachedStructureSequenceDataset - Dataset combining live sequence + cached structure
    ordered_intersection     - deterministic intersection preserving left order
"""
from __future__ import annotations

import os
import pickle
import random
import sys
from typing import Any, Collection, Literal, Mapping, Optional, Sequence

# Ensure src/ is importable
_SRC_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..", "src"))
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

import lmdb
import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset
from torch_geometric.data import Data
from torch_geometric.nn import knn_graph

from Datasets.data_representer import Reaction, StructureSequenceData
from Datasets.Structure.protein_ligand import ATOM_FEATS as _SRC_ATOM_FEATS
from Models.common import GaussianSmearing

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
CACHE_SCHEMA_VERSION = 1
LMDB_MAP_SIZE = 256 * 1024 ** 3  # 256 GB
MAX_SUBSTRATE_LENGTH = 280

# Mirrors src/Datasets/Structure/transforms.py:23
PROTEIN_ATOMIC_NUMBERS = torch.LongTensor([1, 6, 7, 8, 16, 34])  # H,C,N,O,S,Se
MAX_NUM_AA = 21

# Mirrors src/Datasets/Structure/transforms.py:42
LIGAND_ATOMIC_NUMBERS = torch.LongTensor([1, 6, 7, 8, 9, 15, 16, 17, 26, 35, 53])

# Use the authoritative ATOM_FEATS from source to stay in sync with RDKit version
ATOM_FEATS = dict(_SRC_ATOM_FEATS)  # {'AtomicNumber':1, 'Aromatic':1, 'Degree':6, 'NumHs':6, 'Hybridization':N}

PARTITION_FIELDS = (
    "real_edge_index", "knn_edge_index",
    "real_dist", "knn_dist",
    "real_bond_onehot", "knn_bond_onehot",
    "real_cross", "knn_cross",
)


def ordered_intersection(left: Sequence[int], right: Collection[int]) -> list[int]:
    """Return items from *left* that appear in *right*, preserving *left* order."""
    right_set = right if isinstance(right, set) else set(right)
    return [x for x in left if x in right_set]


# ---------------------------------------------------------------------------
# LMDB helpers
# ---------------------------------------------------------------------------
def open_readonly_lmdb(path: str) -> lmdb.Environment:
    return lmdb.open(
        path, map_size=LMDB_MAP_SIZE,
        create=False, subdir=False, readonly=True,
        lock=False, readahead=False, meminit=False,
    )


def _scan_keys(path: Optional[str]) -> set[bytes]:
    if path is None or not os.path.exists(path):
        return set()
    env = open_readonly_lmdb(path)
    try:
        with env.begin(write=False) as txn:
            return set(txn.cursor().iternext(values=False))
    finally:
        env.close()


def _normalize_paths(raw) -> list[Optional[str]]:
    if raw is None:
        return []
    if isinstance(raw, (str, os.PathLike)):
        p = os.fspath(raw)
        return [None if p.strip().lower() in ("", "none") else p]
    result = []
    for item in raw:
        if item is None:
            result.append(None)
        else:
            p = str(item).strip()
            result.append(None if p.lower() in ("", "none") else p)
    return result


# ---------------------------------------------------------------------------
# BuildStructureCacheData — offline transform
# ---------------------------------------------------------------------------
class BuildStructureCacheData:
    """Build raw cache payload from one StructureComplexData sample.

    Replicates FeaturizeProteinAtom + FeaturizeLigandAtom + EdgeConnection
    but stores split real/knn partitions without noise or smearing.
    """

    def __init__(self, k: int = 48, max_substrate_length: int = MAX_SUBSTRATE_LENGTH):
        self.k = k
        self.max_substrate_length = max_substrate_length
        assert max_substrate_length == MAX_SUBSTRATE_LENGTH, (
            f"Model hardcodes 280 in multiple places; got {max_substrate_length}"
        )

    # -- Featurization (exact replicas of transforms.py) --------------------

    @staticmethod
    def _featurize_protein(data: Data) -> torch.Tensor:
        element = data.protein_element.view(-1, 1) == PROTEIN_ATOMIC_NUMBERS.view(1, -1)
        amino_acid = F.one_hot(data.protein_atom_to_aa_type.long(), num_classes=MAX_NUM_AA)
        is_backbone = data.protein_is_backbone.view(-1, 1).long()
        return torch.cat([element, amino_acid, is_backbone], dim=-1)

    @staticmethod
    def _featurize_ligand(data: Data) -> torch.Tensor:
        element = data.ligand_element.view(-1, 1) == LIGAND_ATOMIC_NUMBERS.view(1, -1)
        parts = []
        for i, (key, size) in enumerate(ATOM_FEATS.items()):
            feat = data.ligand_atom_feature[:, i:i + 1]
            if size > 1:
                feat = (feat == torch.LongTensor(list(range(size))).view(1, -1))
            elif key == "AtomicNumber":
                feat = feat / 100.0
            parts.append(feat)
        atom_feat = torch.cat(parts, dim=-1)
        return torch.cat([element, atom_feat], dim=-1)

    # -- k-NN graph (exact replica of transforms.py:97-107) -----------------

    def _build_knn_graph(self, data: Data) -> torch.Tensor:
        if getattr(data, "str_tag", "complex") == "ligand":
            return torch.zeros((2, 0), dtype=torch.long)
        pos = torch.cat([data.ligand_pos, data.protein_pos], dim=0)
        real_set = {
            (int(a), int(b))
            for a, b in zip(data.ligand_bond_index[0], data.ligand_bond_index[1])
        }
        knn_idx = knn_graph(pos, k=self.k, flow="target_to_source")
        kept = [e[:, None] for e in knn_idx.T if (int(e[0]), int(e[1])) not in real_set]
        if not kept:
            return torch.zeros((2, 0), dtype=torch.long)
        return torch.cat(kept, dim=1).long()

    # -- Raw distances (no noise, no smearing) ------------------------------

    @staticmethod
    def _raw_dist(edge_index: torch.Tensor, pos: torch.Tensor) -> torch.Tensor:
        if edge_index.numel() == 0:
            return torch.zeros(0, dtype=torch.float32)
        dst, src = edge_index
        return torch.norm(pos[dst] - pos[src], p=2, dim=-1).float()

    # -- Main entry ---------------------------------------------------------

    def __call__(self, data: Data) -> dict[str, Any]:
        prot_feat = self._featurize_protein(data)
        lig_feat = self._featurize_ligand(data)

        n_sub = data.ligand_pos.size(0)
        n_prot = data.protein_pos.size(0)
        pos = torch.cat([data.ligand_pos, data.protein_pos], dim=0)

        real_ei = data.ligand_bond_index.long()
        knn_ei = self._build_knn_graph(data).long()

        real_dist = self._raw_dist(real_ei, pos)
        knn_dist = self._raw_dist(knn_ei, pos)

        # Bond type one-hot (remap 12→5, then one_hot of bond_type-1)
        bt = data.ligand_bond_type.clone().long()
        bt[bt == 12] = 5
        real_bond_oh = F.one_hot((bt - 1).long(), num_classes=6).float()
        knn_bond_oh = F.one_hot(
            torch.full((knn_ei.size(1),), 5, dtype=torch.long), num_classes=6
        ).float()

        # Cross-bond flags
        knn_cross = (
            ((knn_ei[0] < n_sub) & (knn_ei[1] >= n_sub))
            | ((knn_ei[0] >= n_sub) & (knn_ei[1] < n_sub))
        ).view(-1, 1).float()
        real_cross = torch.zeros(real_ei.size(1), 1, dtype=torch.float32)

        # Node tensors (same as transforms.py:143-148)
        protein_x = torch.cat([
            torch.zeros(n_sub, prot_feat.shape[1]),
            prot_feat.float(),
        ], dim=0)
        ligand_x = torch.cat([
            lig_feat.float(),
            torch.zeros(n_prot, lig_feat.shape[1]),
        ], dim=0)
        ligand_mask = torch.cat([torch.ones(n_sub), torch.zeros(n_prot)])
        protein_mask = torch.cat([torch.zeros(n_sub), torch.ones(n_prot)])
        ligand_index = torch.cat([
            data.ligand_index.long(),
            torch.full((n_prot,), self.max_substrate_length, dtype=torch.long),
        ])

        return {
            "str_tag": str(getattr(data, "str_tag", "complex")),
            "sample_weight": getattr(data, "sample_weight", None),
            "mask_use_complex": getattr(data, "mask_use_complex", None),
            "mask_use_ligand": getattr(data, "mask_use_ligand", None),
            "y": getattr(data, "y", None),
            "protein_x": protein_x,
            "ligand_x": ligand_x,
            "ligand_mask": ligand_mask,
            "protein_mask": protein_mask,
            "ligand_index": ligand_index,
            "real_edge_index": real_ei,
            "knn_edge_index": knn_ei,
            "real_dist": real_dist,
            "knn_dist": knn_dist,
            "real_bond_onehot": real_bond_oh,
            "knn_bond_onehot": knn_bond_oh,
            "real_cross": real_cross,
            "knn_cross": knn_cross,
        }


# ---------------------------------------------------------------------------
# RebuildComplexEdgeAttr — runtime transform
# ---------------------------------------------------------------------------
class RebuildComplexEdgeAttr:
    """Rebuild complex_edge_index and complex_edge_attr from cached partitions.

    mode='fixed':      attr aligned with index [real, knn] — corrects original bug
    mode='legacy_bug': attr is [knn, real] while index is [real, knn] — original behavior
    """

    def __init__(
        self,
        mode: Literal["fixed", "legacy_bug"] = "fixed",
        *,
        is_train: bool,
        dist_noise: bool,
        cutoff: float,
        num_r_gaussian: int,
        clear_partitions: bool = True,
    ):
        assert mode in ("fixed", "legacy_bug"), f"Unknown mode: {mode}"
        self.mode = mode
        self.apply_noise = is_train and dist_noise
        self.clear_partitions = clear_partitions
        self.smearing = GaussianSmearing(stop=cutoff, num_gaussians=num_r_gaussian)

    def _add_noise(self, dist: torch.Tensor) -> torch.Tensor:
        if not self.apply_noise or dist.numel() == 0:
            return dist
        noise = torch.from_numpy(
            np.random.laplace(0.001994, 0.031939, (dist.size(0),))
        ).float()
        return dist + noise

    def __call__(self, data: Data) -> Data:
        real_ei = data.real_edge_index.long()
        knn_ei = data.knn_edge_index.long()

        real_dist = self._add_noise(data.real_dist.float())
        knn_dist = self._add_noise(data.knn_dist.float())

        real_attr = torch.cat([
            data.real_bond_onehot.float(),
            self.smearing(real_dist),
            data.real_cross.float(),
        ], dim=-1)
        knn_attr = torch.cat([
            data.knn_bond_onehot.float(),
            self.smearing(knn_dist),
            data.knn_cross.float(),
        ], dim=-1)

        # Index is always [real, knn]
        data.complex_edge_index = torch.cat([real_ei, knn_ei], dim=-1)

        # Attr order depends on mode
        if self.mode == "fixed":
            data.complex_edge_attr = torch.cat([real_attr, knn_attr], dim=0)
        else:  # legacy_bug: original mismatch
            data.complex_edge_attr = torch.cat([knn_attr, real_attr], dim=0)

        if self.clear_partitions:
            for key in PARTITION_FIELDS:
                setattr(data, key, None)

        return data


# ---------------------------------------------------------------------------
# CachedStructureSequenceDataset
# ---------------------------------------------------------------------------
class _CachedData(Data):
    """Thin wrapper so cached dicts behave as PyG Data objects."""
    pass


def _dict_to_data(payload: dict[str, Any]) -> _CachedData:
    """Convert a cache dict to a Data-like object, coercing types."""
    kw = {}
    for key, val in payload.items():
        if isinstance(val, np.ndarray):
            val = torch.from_numpy(val)
            if val.dtype == torch.float64:
                val = val.float()
        elif isinstance(val, np.generic):
            val = torch.as_tensor(val)
        elif key in ("sample_weight", "y") and val is not None and not isinstance(val, torch.Tensor):
            val = torch.tensor(val, dtype=torch.float32)
        kw[key] = val
    return _CachedData(**kw)


def _close_reaction_lmdb(reaction_db: Reaction) -> None:
    """Close LMDB handles on a Reaction object, enabling pickle for spawn workers."""
    for attr in ("grover_dbs", "enzyme_dbs", "reaction_dbs"):
        dbs = getattr(reaction_db, attr, None)
        if dbs is not None:
            for db in dbs:
                if hasattr(db, "close"):
                    db.close()
            setattr(reaction_db, attr, None)


class CachedStructureSequenceDataset(Dataset):
    """Sequence-live, structure-cached dataset with Windows-safe lazy LMDB opens.

    Sequence features read live from Reaction LMDBs (enzyme/reaction/grover/morgan).
    Structure features read from pre-built cache LMDBs (built by build_structure_cache.py).
    Merged via StructureSequenceData.from_sequence_structure, then runtime transform applied.
    """

    def __init__(
        self,
        df,
        config,
        complex_cache_paths,
        *,
        structure_transform: Optional[RebuildComplexEdgeAttr] = None,
        is_train: bool = False,
        ligand_cache_paths=None,
    ):
        super().__init__()
        self.df = df
        self.config = config
        self.is_train = is_train
        self.structure_transform = structure_transform
        self.full_data = bool(getattr(config.data, "full_data", False))

        # Cache paths
        self.complex_cache_paths = _normalize_paths(complex_cache_paths)
        assert self.complex_cache_paths, "At least one complex cache path required"
        lc = _normalize_paths(ligand_cache_paths)
        while len(lc) < len(self.complex_cache_paths):
            lc.append(None)
        self.ligand_cache_paths = lc

        # Build sequence DB and get valid idx (opens LMDB temporarily)
        self.sequence_db = Reaction(config=config, df=df)
        self.sequence_valid_idx = list(self.sequence_db.valid_idx)
        # Close handles so the main process can pickle this to workers
        _close_reaction_lmdb(self.sequence_db)

        # High quality IDs
        self._hq_dicts = self._load_hq_dicts()

        # Scan cache keys and build structure valid idx
        self._complex_keys = [_scan_keys(p) for p in self.complex_cache_paths]
        self._ligand_keys = [_scan_keys(p) for p in self.ligand_cache_paths]
        self.idx_valid_dict: dict[int, tuple[bool, bool]] = {}
        structure_valid_idx = self._build_structure_valid_idx()

        # Final valid idx = ordered intersection (fixes set() nondeterminism)
        self.valid_idx = ordered_intersection(self.sequence_valid_idx, set(structure_valid_idx))
        print(f"[CachedDataset] seq_valid={len(self.sequence_valid_idx)}, "
              f"str_valid={len(structure_valid_idx)}, final={len(self.valid_idx)}")

        # LMDB envs — opened lazily in workers
        self._complex_envs: Optional[list] = None
        self._ligand_envs: Optional[list] = None

    def __len__(self) -> int:
        return len(self.valid_idx)

    def __getstate__(self) -> dict:
        """Strip LMDB handles before Windows spawn pickling."""
        state = self.__dict__.copy()
        state["_complex_envs"] = None
        state["_ligand_envs"] = None
        # Sequence DB will be recreated lazily in workers
        if state.get("sequence_db") is not None:
            _close_reaction_lmdb(state["sequence_db"])
        return state

    def _ensure_connections(self) -> None:
        """Lazily open LMDB envs in the current worker process."""
        if self._complex_envs is None:
            self._complex_envs = [
                open_readonly_lmdb(p) if p and os.path.exists(p) else None
                for p in self.complex_cache_paths
            ]
        if self._ligand_envs is None:
            self._ligand_envs = [
                open_readonly_lmdb(p) if p and os.path.exists(p) else None
                for p in self.ligand_cache_paths
            ]

    def _ensure_sequence_db(self) -> None:
        """Reconnect sequence LMDB in the current worker process."""
        if self.sequence_db is None:
            self.sequence_db = Reaction(config=self.config, df=self.df)
            # Set our pre-computed valid_idx to avoid re-scanning
            self.sequence_db.valid_idx = self.sequence_valid_idx

    # -- High quality ID loading -------------------------------------------

    def _load_hq_dicts(self) -> list[Optional[dict[int, bool]]]:
        paths = _normalize_paths(getattr(self.config.data, "high_quality_id_path", None))
        while len(paths) < len(self.complex_cache_paths):
            paths.append(None)
        result = []
        for p in paths[:len(self.complex_cache_paths)]:
            if p is None or not os.path.exists(p):
                result.append(None)
                continue
            with open(p, "r", encoding="utf-8") as f:
                result.append({int(line.strip()): True for line in f if line.strip()})
        return result

    def _check_hq(self, dataset_id: int, dock_index: int) -> bool:
        if dataset_id >= len(self._hq_dicts):
            return True
        hq = self._hq_dicts[dataset_id]
        return hq is None or dock_index in hq

    # -- Structure valid idx -----------------------------------------------

    def _build_structure_valid_idx(self) -> list[int]:
        valid = []
        for idx in range(self.df.index.shape[0]):
            dock_idx = int(self.df["Dock Index"].values[idx])
            ds_id = int(self.df["dataset_id"].values[idx])
            sub_idx = int(self.df["Substrate Index"].values[idx])

            has_complex = (
                ds_id < len(self._complex_keys)
                and str(dock_idx).encode() in self._complex_keys[ds_id]
            )
            if not self.full_data:
                self.idx_valid_dict[idx] = (has_complex, False)
                if has_complex and self._check_hq(ds_id, dock_idx):
                    valid.append(idx)
            else:
                has_ligand = (
                    ds_id < len(self._ligand_keys)
                    and str(sub_idx).encode() in self._ligand_keys[ds_id]
                )
                self.idx_valid_dict[idx] = (has_complex, has_ligand)
                if has_complex or has_ligand:
                    valid.append(idx)
        return valid

    # -- Branch selection (mirrors structure.py:317) -----------------------

    def _select_branch(self, idx: int) -> Literal["complex", "ligand"]:
        dock_idx = int(self.df["Dock Index"].values[idx])
        ds_id = int(self.df["dataset_id"].values[idx])

        no_structure = False
        try:
            no_structure = bool(self.config.data.no_structure)
        except Exception:
            pass
        rand_val = random.random() if no_structure else 1.0
        fake_ratio = float(getattr(self.config.data, "fake_sequence_ratio", 0.0)) if self.is_train else 0.0

        fc, fl = self.idx_valid_dict[idx]
        use_complex = (
            (self._check_hq(ds_id, dock_idx) and rand_val > fake_ratio) or (not fl)
        ) and fc
        return "complex" if use_complex else "ligand"

    # -- Load cached structure payload -------------------------------------

    def _load_cached(self, env: Optional[lmdb.Environment], key: bytes) -> _CachedData:
        assert env is not None, f"Cache env unavailable for key {key!r}"
        with env.begin(write=False) as txn:
            raw = txn.get(key)
        assert raw is not None, f"Key not found in cache: {key!r}"
        return _dict_to_data(pickle.loads(raw))

    # -- Main getitem ------------------------------------------------------

    def getitem_with_real_idx(self, idx: int) -> StructureSequenceData:
        self._ensure_sequence_db()
        self._ensure_connections()

        # Sequence data (live from LMDB)
        sequence_data = self.sequence_db.getitem_with_real_idx(idx)

        # Structure data (from cache)
        branch = self._select_branch(idx)
        ds_id = int(self.df["dataset_id"].values[idx])
        if branch == "complex":
            dock_idx = int(self.df["Dock Index"].values[idx])
            structure_data = self._load_cached(self._complex_envs[ds_id], str(dock_idx).encode())
        else:
            sub_idx = int(self.df["Substrate Index"].values[idx])
            structure_data = self._load_cached(self._ligand_envs[ds_id], str(sub_idx).encode())

        # Merge
        data = StructureSequenceData.from_sequence_structure(sequence_data, structure_data)

        # Runtime edge rebuild
        if self.structure_transform is not None:
            data = self.structure_transform(data)

        # Post-merge checks (matches data_representer.py:244-248)
        if getattr(data, "y", None) is not None and getattr(data, "label", None) is not None:
            assert torch.as_tensor(data.y).item() == torch.as_tensor(data.label).item()
        data.active_site = None
        data.y = None
        data.id = idx
        return data

    def __getitem__(self, item_idx: int) -> StructureSequenceData:
        return self.getitem_with_real_idx(self.valid_idx[item_idx])
