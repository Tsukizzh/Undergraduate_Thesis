#!/usr/bin/env python3
"""Regenerate structure_features_heme.lmdb with parallel parsing.

Uses the PATCHED PDBProtein from EXP002a (reads HEM HETATM, feature_dim=31).
Existing pocket PDBs already contain HEM lines — no re-extraction needed.

LMDB key = str(row_idx) where row_idx is the row index in data.csv.
Resume-safe: skips keys already in LMDB.

Codex-reviewed (session 019d4552).
"""
import argparse
import os
import pickle
import sys
import time
import traceback
from multiprocessing import Pool

import lmdb
import pandas as pd
import torch


# Globals set once per worker by init_worker().
PDBProtein = None
parse_sdf_file_mol = None
StructureComplexData = None
torchify_dict = None


def expand(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


def init_worker(exp_src: str) -> None:
    global PDBProtein, parse_sdf_file_mol, StructureComplexData, torchify_dict

    exp_src = expand(exp_src)
    if exp_src not in sys.path:
        sys.path.insert(0, exp_src)

    from Datasets.Structure.protein_ligand import PDBProtein as _PDBProtein
    from Datasets.Structure.protein_ligand import parse_sdf_file_mol as _parse_sdf_file_mol
    from Datasets.Structure.utils import StructureComplexData as _StructureComplexData
    from Datasets.Structure.utils import torchify_dict as _torchify_dict

    PDBProtein = _PDBProtein
    parse_sdf_file_mol = _parse_sdf_file_mol
    StructureComplexData = _StructureComplexData
    torchify_dict = _torchify_dict


def worker_process(item):
    row_idx, enzyme, reaction, label, pocket_path, ligand_path = item
    key = str(row_idx)

    if not os.path.exists(pocket_path):
        return {"status": "missing_pocket", "key": key, "row_idx": row_idx}

    if not os.path.exists(ligand_path):
        return {"status": "missing_ligand", "key": key, "row_idx": row_idx}

    try:
        pocket_dict = PDBProtein(pocket_path).to_dict_atom()
        ligand_dict = parse_sdf_file_mol(ligand_path, heavy_only=True)

        data = StructureComplexData.from_protein_ligand_dicts(
            protein_dict=torchify_dict(pocket_dict),
            ligand_dict=torchify_dict(ligand_dict),
        )

        if data.protein_pos.size(0) == 0:
            return {"status": "parse_error", "key": key, "row_idx": row_idx,
                    "error": "protein_pos is empty"}

        data.protein_filename = enzyme
        data.ligand_filename = reaction
        data.y = torch.tensor(float(label))

        payload = pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL)
        return {"status": "ok", "key": key, "row_idx": row_idx, "payload": payload}
    except Exception:
        return {"status": "parse_error", "key": key, "row_idx": row_idx,
                "error": traceback.format_exc()}


def load_existing_keys(db) -> set:
    existing = set()
    with db.begin(write=False) as txn:
        with txn.cursor() as cursor:
            for key, _ in cursor:
                existing.add(key.decode("utf-8"))
    return existing


def flush_batch(db, batch):
    if not batch:
        return 0
    with db.begin(write=True) as txn:
        for key, payload in batch:
            txn.put(key.encode("utf-8"), payload)
    return len(batch)


def build_tasks(df, pocket_dir, ligand_dir, existing_keys):
    tasks = []
    skipped_existing = 0

    for row_idx in range(len(df)):
        key = str(row_idx)
        if key in existing_keys:
            skipped_existing += 1
            continue

        enzyme = int(df.iloc[row_idx]["enzyme"])
        reaction = int(df.iloc[row_idx]["reaction"])
        label = int(df.iloc[row_idx]["label"])

        pocket_path = os.path.join(pocket_dir, f"{row_idx}.pdb")
        ligand_path = os.path.join(ligand_dir, f"{row_idx}.sdf")

        tasks.append((row_idx, enzyme, reaction, label, pocket_path, ligand_path))

    return tasks, skipped_existing


def parse_args():
    p = argparse.ArgumentParser(
        description="Regenerate structure_features_heme.lmdb with patched PDBProtein (Fe/HEM support)."
    )
    p.add_argument("--exp-src",
        default="~/rivermind-data/EZSpecificity/PathC/P450/experiments/EXP002a_fe_heme_10A/src")
    p.add_argument("--csv-path",
        default="~/rivermind-data/EZSpecificity/PathC/P450/data/data.csv")
    p.add_argument("--pocket-dir",
        default="~/rivermind-data/EZSpecificity/PathC/P450/data/structure/str_tmp_data/pocket")
    p.add_argument("--ligand-dir",
        default="~/rivermind-data/EZSpecificity/PathC/P450/data/structure/str_tmp_data/ligand")
    p.add_argument("--output-lmdb",
        default="~/rivermind-data/EZSpecificity/PathC/P450/data/structure/structure_features_heme.lmdb")
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--batch-size", type=int, default=256,
        help="Buffer N results before one LMDB write transaction")
    p.add_argument("--chunksize", type=int, default=32)
    p.add_argument("--progress-every", type=int, default=2000)
    p.add_argument("--map-size-gb", type=int, default=32)
    return p.parse_args()


def main():
    args = parse_args()

    exp_src = expand(args.exp_src)
    csv_path = expand(args.csv_path)
    pocket_dir = expand(args.pocket_dir)
    ligand_dir = expand(args.ligand_dir)
    output_lmdb = expand(args.output_lmdb)

    # Verify patched src
    sys.path.insert(0, exp_src)
    init_worker(exp_src)
    from Datasets.Structure.transforms import FeaturizeProteinAtom
    fdim = FeaturizeProteinAtom().feature_dim
    assert fdim == 31, f"Expected feature_dim=31 (patched), got {fdim}. Wrong src?"
    print(f"Verified: FeaturizeProteinAtom.feature_dim = {fdim} ✓")

    df = pd.read_csv(csv_path, encoding="utf-8-sig")

    db = lmdb.open(
        output_lmdb,
        map_size=args.map_size_gb * (1024 ** 3),
        create=True, subdir=False, readonly=False, lock=True,
        readahead=False, meminit=False, max_readers=2048,
    )

    start_time = time.time()
    existing_keys = load_existing_keys(db)
    tasks, skipped_existing = build_tasks(df, pocket_dir, ligand_dir, existing_keys)

    print(f"Total rows: {len(df)}")
    print(f"Already in LMDB: {skipped_existing}")
    print(f"To process: {len(tasks)}")
    print(f"Workers: {args.workers}, batch_size: {args.batch_size}, chunksize: {args.chunksize}")

    if len(tasks) == 0:
        print("Nothing to do.")
        db.close()
        return

    written = 0
    missing_pocket = 0
    missing_ligand = 0
    parse_error = 0
    processed = 0
    batch = []

    with Pool(processes=args.workers, initializer=init_worker, initargs=(exp_src,)) as pool:
        for result in pool.imap_unordered(worker_process, tasks, chunksize=args.chunksize):
            processed += 1
            status = result["status"]

            if status == "ok":
                batch.append((result["key"], result["payload"]))
                if len(batch) >= args.batch_size:
                    written += flush_batch(db, batch)
                    batch.clear()
            elif status == "missing_pocket":
                missing_pocket += 1
            elif status == "missing_ligand":
                missing_ligand += 1
            else:
                parse_error += 1

            if processed % args.progress_every == 0:
                if batch:
                    written += flush_batch(db, batch)
                    batch.clear()
                elapsed = time.time() - start_time
                rate = processed / elapsed if elapsed > 0 else 0
                print(
                    f"[{processed}/{len(tasks)}] written={written} "
                    f"missing_pocket={missing_pocket} missing_ligand={missing_ligand} "
                    f"parse_error={parse_error} rate={rate:.1f}/s elapsed={elapsed:.0f}s",
                    flush=True,
                )

    if batch:
        written += flush_batch(db, batch)
        batch.clear()

    db.sync()
    db.close()

    elapsed = time.time() - start_time
    print(f"\n{'='*50}")
    print(f"DONE")
    print(f"Total rows:      {len(df)}")
    print(f"Already existed: {skipped_existing}")
    print(f"Processed:       {processed}")
    print(f"Written:         {written}")
    print(f"Missing pocket:  {missing_pocket}")
    print(f"Missing ligand:  {missing_ligand}")
    print(f"Parse errors:    {parse_error}")
    print(f"Elapsed:         {elapsed:.0f}s")


if __name__ == "__main__":
    main()
