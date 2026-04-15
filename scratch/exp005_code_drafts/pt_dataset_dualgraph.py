"""
EXP005 dualgraph dataset wrapper.

Extends the baseline `PtCacheDataset` to also load per-dock GVP feature files
from `gvp_cache/samples/NNN/gvp_DDDDDD.pt`. Samples that are in the
`gvp_invalid_docks.pt` manifest (64 docks with empty pocket or N=1) receive a
canonical placeholder GVP payload with `gvp_valid=False`, so the model can
gate them at forward time without any training-loop special-casing.

Directory layout (the EXP005 cache path IS the dualgraph overlay):

    <overlay_root>/
    ├── random/                    <--- cache_dir passed to Dataset
    │   ├── enzymes    -> symlink  (baseline ESM flatbin)
    │   ├── substrates -> symlink
    │   ├── train/
    │   │   ├── index.pt   -> symlink (baseline index)
    │   │   ├── samples    -> symlink
    │   │   └── dock_sidecar.pt    (new)
    │   ├── val/  (same)
    │   └── test/ (same)
    ├── enzyme_resid_map.pt
    └── gvp_cache/
        ├── samples/
        ├── manifest.pt
        └── gvp_invalid_docks.pt

So relative to the dataset's `cache_dir` (= overlay_root/random):
    overlay_root   = Path(cache_dir).parent
    gvp_cache_dir  = overlay_root / "gvp_cache"
    sidecar_path   = Path(cache_dir) / split / "dock_sidecar.pt"

Design (codex rounds 2-5):
  - Subclass `PtCacheDataset.__getitem__` to attach GVP fields.
  - Use `DualgraphData(StructureSequenceData)` subclass with __inc__/__cat_dim__.
  - DataLoader uses `follow_batch=['gvp_node_s']` so PyG creates
    `batch.gvp_node_s_batch` automatically.
  - Placeholder payload: 1 node, 0 edges, pocket_residue_idx=-1, gvp_valid=0.
  - `gvp_valid` stored as bool (cat dim 0; per-sample shape [1] → batch [B]).
  - `preload=True` loads all GVP files into RAM at init for hot training.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Optional

import torch

from Datasets.data_representer import StructureSequenceData
from pt_dataset import PtCacheDataset  # baseline


# ---------------------------------------------------------------------------
# Data subclass with dualgraph batching rules
# ---------------------------------------------------------------------------
class DualgraphData(StructureSequenceData):
    """PyG Data subclass with batching rules for the GVP residue subgraph.

    Baseline fields retain their existing `__inc__` / `__cat_dim__`
    semantics via `super().__inc__`. New GVP fields:

      gvp_node_s       [N, 6]        cat dim 0, no offset
      gvp_node_v       [N, 3, 3]     cat dim 0, no offset
      gvp_edge_index   [2, E]        cat dim 1, offset by gvp_node_s.size(0)
      gvp_edge_s       [E, 32]       cat dim 0, no offset
      gvp_edge_v       [E, 1, 3]     cat dim 0, no offset
      gvp_aa_type      [N]           cat dim 0, no offset
      gvp_pocket_residue_idx [N]     cat dim 0, NO offset (absolute index into
                                      padded x_pro sequence, or -1 sentinel)
      gvp_valid        [1]           cat dim 0, no offset (per-sample bool flag)
    """

    def __inc__(self, key, value, *args, **kwargs):
        if key == "gvp_edge_index":
            return self["gvp_node_s"].size(0)
        if key in {
            "gvp_node_s", "gvp_node_v",
            "gvp_edge_s", "gvp_edge_v",
            "gvp_aa_type", "gvp_pocket_residue_idx",
            "gvp_valid",
        }:
            return 0
        return super().__inc__(key, value, *args, **kwargs)

    def __cat_dim__(self, key, value, *args, **kwargs):
        if key == "gvp_edge_index":
            return 1
        if key in {
            "gvp_node_s", "gvp_node_v",
            "gvp_edge_s", "gvp_edge_v",
            "gvp_aa_type", "gvp_pocket_residue_idx",
            "gvp_valid",
        }:
            return 0
        return super().__cat_dim__(key, value, *args, **kwargs)


# ---------------------------------------------------------------------------
# Placeholder + real GVP payload loaders
# ---------------------------------------------------------------------------
_PLACEHOLDER_TEMPLATE = None  # lazy, per-process

def _placeholder_gvp() -> dict:
    """1-node-0-edge canonical placeholder with gvp_valid=False."""
    global _PLACEHOLDER_TEMPLATE
    if _PLACEHOLDER_TEMPLATE is None:
        _PLACEHOLDER_TEMPLATE = {
            "gvp_node_s":          torch.zeros(1, 6, dtype=torch.float32),
            "gvp_node_v":          torch.zeros(1, 3, 3, dtype=torch.float32),
            "gvp_edge_index":      torch.zeros(2, 0, dtype=torch.long),
            "gvp_edge_s":          torch.zeros(0, 32, dtype=torch.float32),
            "gvp_edge_v":          torch.zeros(0, 1, 3, dtype=torch.float32),
            "gvp_aa_type":         torch.zeros(1, dtype=torch.long),
            "gvp_pocket_residue_idx": torch.full((1,), -1, dtype=torch.long),
            "gvp_valid":           torch.zeros(1, dtype=torch.bool),
        }
    # Return fresh tensors so downstream modification doesn't poison the template
    return {k: v.clone() for k, v in _PLACEHOLDER_TEMPLATE.items()}


def _load_real_gvp(path: Path) -> dict:
    """Load a real gvp_{dock}.pt file and coerce into canonical schema."""
    raw = torch.load(str(path), map_location="cpu", weights_only=False)
    return {
        "gvp_node_s":          raw["node_s"].to(torch.float32),
        "gvp_node_v":          raw["node_v"].to(torch.float32),
        "gvp_edge_index":      raw["edge_index"].to(torch.long),
        "gvp_edge_s":          raw["edge_s"].to(torch.float32),
        "gvp_edge_v":          raw["edge_v"].to(torch.float32),
        "gvp_aa_type":         raw["aa_type"].to(torch.long),
        "gvp_pocket_residue_idx": raw["pocket_residue_idx"].to(torch.long),
        "gvp_valid":           torch.ones(1, dtype=torch.bool),
    }


# ---------------------------------------------------------------------------
# Dataset wrapper
# ---------------------------------------------------------------------------
class PtCacheDualgraphDataset(PtCacheDataset):
    """PtCacheDataset + per-dock GVP feature cache + invalid manifest fallback.

    Additional args:
      gvp_cache_dir:         Path to gvp_cache/. If None, defaults to
                             `Path(cache_dir).parent / 'gvp_cache'` which
                             assumes cache_dir is an overlay-root/random
                             directory.
      gvp_invalid_manifest:  Path to gvp_invalid_docks.pt. If None, defaults
                             to `gvp_cache_dir / 'gvp_invalid_docks.pt'`.
      sidecar_path:          Path to dock_sidecar.pt for this split. If None,
                             defaults to `Path(cache_dir) / split / 'dock_sidecar.pt'`.
      preload_gvp:           If True, load all GVP .pt into RAM at init
                             (recommended; total size ~350 MB).

    Path contract: `cache_dir` points at an overlay root/random directory
    containing dock_sidecar.pt per split and symlinking back to the baseline
    index.pt / samples / enzymes / substrates. The default locations derive
    from `cache_dir` assuming this layout. For a non-standard layout, pass
    gvp_cache_dir / gvp_invalid_manifest / sidecar_path explicitly.
    """

    def __init__(
        self,
        cache_dir,
        split,
        gvp_cache_dir: Optional[str] = None,
        gvp_invalid_manifest: Optional[str] = None,
        sidecar_path: Optional[str] = None,
        preload_gvp: bool = True,
        **kwargs,
    ):
        super().__init__(cache_dir=cache_dir, split=split, **kwargs)

        overlay_root = Path(cache_dir).parent  # e.g. .../pt_cache_dualgraph_allfix_unified

        # ---- GVP cache dir ----
        if gvp_cache_dir is None:
            gvp_cache_dir = overlay_root / "gvp_cache"
        self._gvp_cache_dir = Path(gvp_cache_dir)
        self._gvp_samples_dir = self._gvp_cache_dir / "samples"
        if not self._gvp_samples_dir.exists():
            raise FileNotFoundError(
                f"GVP samples dir missing: {self._gvp_samples_dir}"
            )

        # ---- Invalid manifest ----
        if gvp_invalid_manifest is None:
            gvp_invalid_manifest = self._gvp_cache_dir / "gvp_invalid_docks.pt"
        self._gvp_invalid_manifest_path = Path(gvp_invalid_manifest)
        if self._gvp_invalid_manifest_path.exists():
            invalid = torch.load(
                str(self._gvp_invalid_manifest_path),
                map_location="cpu",
                weights_only=False,
            )
            self._invalid_docks: set[int] = set(int(k) for k in invalid.keys())
        else:
            self._invalid_docks = set()
            print(
                f"[PtCacheDualgraphDataset] WARNING: no invalid manifest at "
                f"{self._gvp_invalid_manifest_path}, assuming every dock is valid"
            )

        # ---- Dock sidecar ----
        if sidecar_path is None:
            sidecar_path = Path(cache_dir) / split / "dock_sidecar.pt"
        self._sidecar_path = Path(sidecar_path)
        if not self._sidecar_path.exists():
            raise FileNotFoundError(f"sidecar missing: {self._sidecar_path}")
        sidecar = torch.load(
            str(self._sidecar_path), map_location="cpu", weights_only=False
        )
        # Safety: sidecar sample_ids must match the baseline index.pt sample_ids.
        base_sids = self._index["sample_ids"]
        side_sids = sidecar["sample_ids"]
        if not torch.equal(torch.as_tensor(side_sids), torch.as_tensor(base_sids)):
            raise RuntimeError(
                f"dock_sidecar.pt sample_ids diverge from index.pt sample_ids "
                f"for split={split}. This should never happen — the sidecar "
                f"was built to be row-aligned. Rebuild the overlay."
            )
        self._dock_indices = sidecar["dock_indices"]
        if len(self._dock_indices) != self._n_samples:
            raise RuntimeError(
                f"sidecar/index length mismatch: "
                f"sidecar={len(self._dock_indices)} vs index={self._n_samples}"
            )

        # ---- Optional preload of GVP files into RAM ----
        self._gvp_preloaded: Optional[dict[int, dict]] = None
        if preload_gvp:
            self._preload_gvp_files()

    def _preload_gvp_files(self) -> None:
        """Load all GVP .pt files for this split into RAM. ~350MB total
        for 44026 successful files, cheap compared to baseline preload."""
        t0 = time.time()
        preloaded: dict[int, dict] = {}
        hit = miss = 0
        for k in range(self._n_samples):
            dock = int(self._dock_indices[k])
            if dock in self._invalid_docks:
                continue  # no real file
            shard = dock // 1000
            path = self._gvp_samples_dir / f"{shard:03d}" / f"gvp_{dock:06d}.pt"
            if path.exists():
                preloaded[dock] = _load_real_gvp(path)
                hit += 1
            else:
                miss += 1
        self._gvp_preloaded = preloaded
        print(
            f"[PtCacheDualgraphDataset {self.split}] GVP preload: "
            f"{hit} files into RAM, {miss} missing, {len(self._invalid_docks)} in invalid set, "
            f"{time.time()-t0:.1f}s"
        )

    def _load_gvp_for_idx(self, idx: int) -> dict:
        dock = int(self._dock_indices[idx])
        if dock in self._invalid_docks:
            return _placeholder_gvp()
        if self._gvp_preloaded is not None and dock in self._gvp_preloaded:
            # Shallow clone — tensors are immutable during training
            return dict(self._gvp_preloaded[dock])
        shard = dock // 1000
        path = self._gvp_samples_dir / f"{shard:03d}" / f"gvp_{dock:06d}.pt"
        if not path.exists():
            print(
                f"[PtCacheDualgraphDataset] WARN: missing gvp file for dock={dock},"
                f" falling back to placeholder"
            )
            return _placeholder_gvp()
        return _load_real_gvp(path)

    def __getitem__(self, idx: int) -> DualgraphData:
        base: StructureSequenceData = super().__getitem__(idx)
        # Promote to DualgraphData (same underlying tensors, new batching rules)
        out = DualgraphData(**base.to_dict())
        # Preserve baseline's optional None attrs if present (mirrors the
        # legacy StructureSequence.__getitem__ contract even though
        # PtCacheDataset does not emit them today — defensive).
        out.active_site = getattr(base, "active_site", None)
        out.y = getattr(base, "y", None)
        # Attach GVP fields
        for k, v in self._load_gvp_for_idx(idx).items():
            out[k] = v
        return out

    def __getstate__(self):
        """Inherit baseline drop-file-handles behavior. Preloaded GVP dict
        survives the pickle (workers inherit it via fork on Linux)."""
        return super().__getstate__()


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------
def build_pt_dualgraph_dataloaders(
    cache_dir,
    gvp_cache_dir: Optional[str] = None,
    gvp_invalid_manifest: Optional[str] = None,
    batch_size: int = 14,
    edge_mode: str = "fixed",
    dist_noise_train: bool = True,
    cutoff: float = 10.0,
    num_r_gaussian: int = 32,
    num_workers: int = 4,
    pin_memory: bool = True,
    prefetch_factor: int = 4,
    preload: bool = False,
    preload_gvp: bool = True,
):
    from torch_geometric.loader import DataLoader

    loaders = {}
    for split in ["train", "val", "test"]:
        split_dir = Path(cache_dir) / split
        if not split_dir.exists():
            continue

        dataset = PtCacheDualgraphDataset(
            cache_dir=cache_dir,
            split=split,
            gvp_cache_dir=gvp_cache_dir,
            gvp_invalid_manifest=gvp_invalid_manifest,
            edge_mode=edge_mode,
            dist_noise=(dist_noise_train if split == "train" else False),
            cutoff=cutoff,
            num_r_gaussian=num_r_gaussian,
            preload=preload,
            preload_gvp=preload_gvp,
        )

        kw = dict(
            batch_size=batch_size,
            shuffle=(split == "train"),
            num_workers=num_workers,
            pin_memory=pin_memory,
            # PyG will auto-create:
            #   batch.ligand_index_batch  (baseline)
            #   batch.gvp_node_s_batch    (new)
            follow_batch=["ligand_index", "gvp_node_s"],
            persistent_workers=(num_workers > 0),
        )
        if num_workers > 0:
            kw["prefetch_factor"] = prefetch_factor

        loaders[split] = DataLoader(dataset, **kw)

    return loaders
