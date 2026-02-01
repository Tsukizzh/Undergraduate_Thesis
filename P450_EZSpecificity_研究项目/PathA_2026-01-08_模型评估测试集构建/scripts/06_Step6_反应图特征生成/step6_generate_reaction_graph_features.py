"""
Step 6: Reaction Graph Feature Generation (substrate-only)

Input:
  - Substrates.csv with columns:
      Substrate_Index, Substrate_SMILES

Output:
  - reaction_features.lmdb
      key   = str(Substrate_Index).encode()
      value = pickle.dumps(dict) where dict contains:
          element, edge_index, edge_type, atom_feature
          reaction_attention_label (all zeros)
          num_nodes

Note: This script uses a custom parse_smile implementation that doesn't require torch_scatter.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import pickle
import platform
import sys
from dataclasses import dataclass
from multiprocessing import Pool, freeze_support
from pathlib import Path
from typing import Any, Iterable

import lmdb
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.rdchem import BondType, HybridizationType
from tqdm import tqdm


def _find_repo_root(start: Path) -> Path:
    for p in [start] + list(start.parents):
        if (p / "src" / "Datasets").exists():
            return p
    raise FileNotFoundError("Could not locate repo root.")


def _ensure_src_on_syspath() -> Path:
    script_dir = Path(__file__).parent.resolve()
    repo_root = _find_repo_root(script_dir)
    src_path = str(repo_root / "src")
    if src_path not in sys.path:
        sys.path.insert(0, src_path)
    return repo_root


def _resolve_input_csv(cli_path: str | None, project_dir: Path) -> Path:
    if cli_path:
        p = Path(cli_path)
        if not p.exists():
            raise FileNotFoundError(f"--input-csv not found: {p}")
        return p

    step4_default = project_dir / "data" / "04_Step4_格式修正后数据" / "Substrates.csv"
    if step4_default.exists():
        return step4_default

    cwd_candidate = Path("Substrates.csv")
    if cwd_candidate.exists():
        return cwd_candidate

    data_dir = project_dir / "data"
    candidates = list(data_dir.rglob("Substrates.csv")) if data_dir.exists() else []
    if not candidates:
        raise FileNotFoundError("Could not find Substrates.csv.")
    candidates.sort(key=lambda p: (0 if "04_Step4" in str(p) else 1, len(str(p))))
    return candidates[0]


def _default_output_dir(project_dir: Path) -> Path:
    return project_dir / "data" / "06_Step6_反应图特征"


def _default_output_lmdb(project_dir: Path) -> Path:
    return _default_output_dir(project_dir) / "reaction_features.lmdb"


def _default_report_csv(project_dir: Path) -> Path:
    return _default_output_dir(project_dir) / "reaction_features_report.csv"


def _default_session_log(project_dir: Path) -> Path:
    return project_dir / "sessions" / "06_Step6_反应图特征生成" / "session_log.md"


# ============================================================
# Custom parse_smile implementation (no torch_scatter required)
# ============================================================

def get_ligand_atom_features_no_scatter(rdmol):
    """
    Extract atom features without using torch_scatter.
    Returns a (num_atoms, 7) feature matrix.
    """
    num_atoms = rdmol.GetNumAtoms()
    atomic_number = []
    aromatic = []
    sp, sp2, sp3 = [], [], []
    degree = []
    index_to_newindex_dict = {}

    for index, atom in enumerate(rdmol.GetAtoms()):
        atom_idx = atom.GetIdx()
        index_to_newindex_dict[atom_idx] = index
        atom = rdmol.GetAtomWithIdx(atom_idx)
        atomic_number.append(atom.GetAtomicNum())
        aromatic.append(1 if atom.GetIsAromatic() else 0)
        hybridization = atom.GetHybridization()
        sp.append(1 if hybridization == HybridizationType.SP else 0)
        sp2.append(1 if hybridization == HybridizationType.SP2 else 0)
        sp3.append(1 if hybridization == HybridizationType.SP3 else 0)
        degree.append(atom.GetDegree())

    # Build edge lists
    row, col = [], []
    for bond in rdmol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        start = index_to_newindex_dict[start]
        end = index_to_newindex_dict[end]
        row.append(start)
        col.append(end)
        row.append(end)
        col.append(start)

    # Compute num_hs (number of adjacent hydrogen atoms) without torch_scatter
    # This replaces: num_hs = scatter(hs[row], col, dim_size=num_atoms)
    num_hs = np.zeros(num_atoms, dtype=np.float32)
    for i in range(len(row)):
        if atomic_number[row[i]] == 1:  # If source atom is hydrogen
            num_hs[col[i]] += 1

    feat_mat = np.array([atomic_number, aromatic, degree, num_hs, sp, sp2, sp3], dtype=np.int64).transpose()
    return feat_mat


def parse_smile_no_scatter(string: str):
    """
    Parse SMILES string and extract molecular graph features.
    This is a self-contained version that doesn't require torch_scatter.
    """
    mol = Chem.MolFromSmiles(string)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {string}")

    # Canonicalize
    tmp_mol = Chem.MolFromSmiles(string)
    Chem.SanitizeMol(tmp_mol)
    for atom in tmp_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    canonical_mol = Chem.MolFromSmiles(Chem.MolToSmiles(tmp_mol))

    # Handle _smilesAtomOutputOrder - may be empty for simple molecules
    output_order_prop = tmp_mol.GetProp("_smilesAtomOutputOrder")
    # Remove brackets and split, filter empty strings
    output_order_str = output_order_prop[1:-2] if len(output_order_prop) > 2 else ""
    if output_order_str:
        mapping_idx = list(map(int, [x for x in output_order_str.split(",") if x.strip()]))
    else:
        # Fallback: identity mapping
        mapping_idx = list(range(mol.GetNumAtoms()))

    for atom in canonical_mol.GetAtoms():
        idx = atom.GetIdx()
        if idx < len(mapping_idx):
            atom.SetAtomMapNum(mol.GetAtomWithIdx(mapping_idx[idx]).GetAtomMapNum())
    mol = canonical_mol

    # Get atom features
    feat_mat = get_ligand_atom_features_no_scatter(mol)

    num_atoms = mol.GetNumAtoms()

    # Extract element array
    element = []
    for index in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(index)
        element.append(atom.GetAtomicNum())
    element = np.array(element, dtype=np.int64)

    # Extract edge info
    row, col, edge_type_list = [], [], []
    BOND_TYPES = {t: i for i, t in enumerate(BondType.names.values())}
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        edge_type_list += 2 * [BOND_TYPES[bond.GetBondType()]]

    edge_index = np.array([row, col], dtype=np.int64)
    edge_type = np.array(edge_type_list, dtype=np.int64)

    # Sort edges
    if edge_index.shape[1] > 0:
        perm = (edge_index[0] * num_atoms + edge_index[1]).argsort()
        edge_index = edge_index[:, perm]
        edge_type = edge_type[perm]

    data = {
        'element': element,
        'edge_index': edge_index,
        'edge_type': edge_type,
        'atom_feature': feat_mat
    }
    return data


# ============================================================
# Feature generation logic
# ============================================================

@dataclass(frozen=True)
class FeaturizeTask:
    substrate_index: int
    smiles: str


@dataclass(frozen=True)
class FeaturizeResult:
    substrate_index: int
    smiles: str
    status: str
    reason: str
    num_atoms: int | None
    data: dict[str, Any] | None


def _rdkit_precheck(smiles: str, max_atoms: int) -> tuple[Chem.Mol | None, str | None]:
    if not isinstance(smiles, str):
        return None, "SMILES is not a string"
    smiles = smiles.strip()
    if not smiles:
        return None, "Empty SMILES"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES (RDKit MolFromSmiles returned None)"

    mol_no_h = Chem.rdmolops.RemoveAllHs(mol)
    if mol_no_h.GetNumAtoms() == 0:
        return None, "Molecule has no non-H atoms after RemoveAllHs"

    if mol.GetNumAtoms() > max_atoms:
        return None, f"Too many atoms ({mol.GetNumAtoms()} > {max_atoms})"

    return mol, None


def _featurize_one(args: tuple[FeaturizeTask, int]) -> FeaturizeResult:
    task, max_atoms = args
    smiles = task.smiles

    mol, err = _rdkit_precheck(smiles, max_atoms=max_atoms)
    if err is not None:
        return FeaturizeResult(
            substrate_index=task.substrate_index,
            smiles=smiles,
            status="skipped",
            reason=err,
            num_atoms=None if mol is None else int(mol.GetNumAtoms()),
            data=None,
        )

    try:
        substrate_data = parse_smile_no_scatter(smiles)
        n = int(substrate_data["element"].shape[0])
        data: dict[str, Any] = {
            "element": substrate_data["element"],
            "edge_index": substrate_data["edge_index"],
            "edge_type": substrate_data["edge_type"],
            "atom_feature": substrate_data["atom_feature"],
            "reaction_attention_label": np.zeros((n,), dtype=np.float32),
            "num_nodes": np.array(n, dtype=np.int64),
        }
        return FeaturizeResult(
            substrate_index=task.substrate_index,
            smiles=smiles,
            status="ok",
            reason="OK",
            num_atoms=n,
            data=data,
        )
    except Exception as e:
        return FeaturizeResult(
            substrate_index=task.substrate_index,
            smiles=smiles,
            status="failed",
            reason=f"parse_smile failed: {type(e).__name__}: {e}",
            num_atoms=int(mol.GetNumAtoms()) if mol is not None else None,
            data=None,
        )


def _iter_results(tasks: list[FeaturizeTask], max_atoms: int, num_workers: int) -> Iterable[FeaturizeResult]:
    if num_workers <= 1:
        for t in tasks:
            yield _featurize_one((t, max_atoms))
        return

    chunksize = 16
    with Pool(processes=num_workers) as pool:
        for res in pool.imap_unordered(_featurize_one, ((t, max_atoms) for t in tasks), chunksize=chunksize):
            yield res


def _load_existing_keys(env: lmdb.Environment) -> set[bytes]:
    with env.begin(write=False) as txn:
        return set(txn.cursor().iternext(values=False))


def _validate_columns(df: pd.DataFrame) -> None:
    required = {"Substrate_Index", "Substrate_SMILES"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Substrates.csv missing required columns: {sorted(missing)}")


def _validate_unique_indices(df: pd.DataFrame) -> None:
    try:
        idx_int = df["Substrate_Index"].astype("int64")
    except Exception as e:
        raise ValueError(f"Substrate_Index column cannot be converted to int: {e}")
    dup = idx_int[idx_int.duplicated()].unique()
    if len(dup) > 0:
        raise ValueError(f"Duplicate Substrate_Index values found: {dup.tolist()}")


def _open_lmdb(path: Path, map_size_gb: float) -> lmdb.Environment:
    path.parent.mkdir(parents=True, exist_ok=True)
    map_size = int(map_size_gb * 1024 * 1024 * 1024)
    return lmdb.open(str(path), map_size=map_size, subdir=False, create=True, readonly=False)


def _validate_lmdb_contents(output_lmdb: Path, expected_keys: set[bytes]) -> tuple[int, list[str]]:
    env = lmdb.open(
        str(output_lmdb),
        subdir=False,
        create=False,
        readonly=True,
        lock=False,
        readahead=False,
        meminit=False,
    )
    try:
        present = _load_existing_keys(env)
        missing = sorted(k.decode() for k in (expected_keys - present))
        return len(present & expected_keys), missing
    finally:
        env.close()


def _write_session_log(
    session_log_path: Path,
    input_csv: Path,
    output_lmdb: Path,
    report_csv: Path,
    args: argparse.Namespace,
    stats: dict[str, int],
    report_df: pd.DataFrame,
) -> None:
    session_log_path.parent.mkdir(parents=True, exist_ok=True)
    now = _dt.datetime.now().isoformat(timespec="seconds")
    reason_counts = (
        report_df[report_df["status"] != "ok"]
        .groupby(["status", "reason"])
        .size()
        .reset_index(name="count")
        .sort_values(["status", "count"], ascending=[True, False])
    )

    lines: list[str] = []
    lines.append(f"\n## Run {now}\n")
    lines.append(f"- input_csv: {input_csv}\n")
    lines.append(f"- output_lmdb: {output_lmdb}\n")
    lines.append(f"- report_csv: {report_csv}\n")
    lines.append(f"- args: {vars(args)}\n")
    lines.append(f"- python: {sys.version.split()[0]} ({platform.platform()})\n")
    lines.append(
        f"- stats: total={stats['total']} ok={stats['ok']} skipped={stats['skipped']} failed={stats['failed']} skipped_exists={stats['skipped_exists']}\n"
    )
    if len(reason_counts) > 0:
        lines.append("\nTop non-ok reasons:\n")
        for _, r in reason_counts.head(20).iterrows():
            lines.append(f"- {r['status']}: {r['count']} - {r['reason']}\n")

    with open(session_log_path, "a", encoding="utf-8") as f:
        f.writelines(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Step 6: generate substrate graph features into LMDB.")
    parser.add_argument("--input-csv", default=None, help="Path to Substrates.csv.")
    parser.add_argument("--output-lmdb", default=None, help="Path to reaction_features.lmdb.")
    parser.add_argument("--report-csv", default=None, help="Path to write a CSV report.")
    parser.add_argument("--session-log", default=None, help="Path to session_log.md.")
    parser.add_argument("--no-session-log", action="store_true", help="Disable writing session log.")
    parser.add_argument("--max-atoms", type=int, default=280, help="Skip molecules with > max-atoms atoms.")
    parser.add_argument("--resume", action="store_true", help="Skip keys that already exist.")
    parser.add_argument("--num-workers", type=int, default=1, help="Featurization workers.")
    parser.add_argument("--map-size-gb", type=float, default=4.0, help="LMDB map size in GB.")
    parser.add_argument("--validate", action="store_true", help="Validate keys after generation.")
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")

    script_dir = Path(__file__).parent.resolve()
    project_dir = script_dir.parent.parent

    input_csv = _resolve_input_csv(args.input_csv, project_dir=project_dir)
    output_lmdb = Path(args.output_lmdb) if args.output_lmdb else _default_output_lmdb(project_dir)
    report_csv = Path(args.report_csv) if args.report_csv else _default_report_csv(project_dir)
    session_log = Path(args.session_log) if args.session_log else _default_session_log(project_dir)

    print("=" * 60)
    print("Step 6: Reaction Graph Feature Generation")
    print("=" * 60)
    print(f"Input CSV:    {input_csv}")
    print(f"Output LMDB:  {output_lmdb}")
    print(f"Report CSV:   {report_csv}")
    print(f"Max atoms:    {args.max_atoms}")
    print(f"Num workers:  {args.num_workers}")
    print("=" * 60)

    df = pd.read_csv(input_csv)
    _validate_columns(df)
    _validate_unique_indices(df)

    tasks: list[FeaturizeTask] = []
    for _, row in df.iterrows():
        tasks.append(FeaturizeTask(substrate_index=int(row["Substrate_Index"]), smiles=row["Substrate_SMILES"]))

    print(f"Total substrates to process: {len(tasks)}")

    env = _open_lmdb(output_lmdb, map_size_gb=args.map_size_gb)
    try:
        existing_keys: set[bytes] = set()
        if args.resume and output_lmdb.exists():
            existing_keys = _load_existing_keys(env)

        report_rows: list[dict[str, Any]] = []
        stats = {"total": len(tasks), "ok": 0, "skipped": 0, "failed": 0, "skipped_exists": 0}

        pending: list[FeaturizeTask] = []
        for t in tasks:
            k = str(t.substrate_index).encode()
            if args.resume and k in existing_keys:
                stats["skipped_exists"] += 1
                report_rows.append({
                    "Substrate_Index": t.substrate_index,
                    "Substrate_SMILES": t.smiles,
                    "status": "skipped",
                    "reason": "Key already exists in LMDB (resume)",
                    "num_atoms": None,
                })
            else:
                pending.append(t)

        for res in tqdm(_iter_results(pending, max_atoms=args.max_atoms, num_workers=max(1, args.num_workers)),
                        total=len(pending), desc="Featurizing"):
            if res.status == "ok" and res.data is not None:
                with env.begin(write=True) as txn:
                    txn.put(key=str(res.substrate_index).encode(),
                            value=pickle.dumps(res.data, protocol=pickle.HIGHEST_PROTOCOL))
                stats["ok"] += 1
            elif res.status == "skipped":
                stats["skipped"] += 1
            else:
                stats["failed"] += 1

            report_rows.append({
                "Substrate_Index": res.substrate_index,
                "Substrate_SMILES": res.smiles,
                "status": res.status,
                "reason": res.reason,
                "num_atoms": res.num_atoms,
            })

        report_csv.parent.mkdir(parents=True, exist_ok=True)
        report_df = pd.DataFrame(report_rows).sort_values(["Substrate_Index"])
        report_df.to_csv(report_csv, index=False, encoding="utf-8")

        print("\n" + "=" * 60)
        print("Step 6 complete")
        print(f"Input CSV:    {input_csv}")
        print(f"Output LMDB:  {output_lmdb}")
        print(f"Report CSV:   {report_csv}")
        print(
            f"Stats: total={stats['total']} ok={stats['ok']} skipped={stats['skipped']} "
            f"failed={stats['failed']} skipped_exists={stats['skipped_exists']}"
        )
        print("=" * 60)

        if not args.no_session_log:
            _write_session_log(session_log, input_csv, output_lmdb, report_csv, args, stats, report_df)
            print(f"Session log:  {session_log}")

    finally:
        env.close()

    if args.validate:
        expected = {str(int(x)).encode() for x in df["Substrate_Index"].astype("int64").tolist()}
        present_cnt, missing = _validate_lmdb_contents(output_lmdb, expected_keys=expected)
        if missing:
            print(f"\nValidation: missing {len(missing)} keys (showing up to 20): {missing[:20]}")
            return 2
        print(f"\nValidation: {present_cnt} / {len(expected)} keys present.")

    return 0


if __name__ == "__main__":
    freeze_support()
    raise SystemExit(main())
