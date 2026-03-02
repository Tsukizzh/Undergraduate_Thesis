"""Generate random negative enzyme-substrate pairs (fix substrate, random enzyme)."""
from __future__ import annotations

import csv
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

TARGET_NEG_RATIO = 9.36


@dataclass
class NegativeSamplingResult:
    output_path: Path
    success: bool
    error: str = ""
    positive_count: int = 0
    negative_count: int = 0
    achieved_ratio: float = 0.0


def _norm_key(key: str) -> str:
    """Normalize CSV column name: strip, lowercase, spaces → underscores."""
    return key.strip().lower().replace(" ", "_")


def _load_positive_pairs(data_csv: Path) -> list[tuple[int, int]]:
    """Load (enzyme_index, substrate_index) pairs where Label == 1."""
    pairs: list[tuple[int, int]] = []
    with data_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            normed = {_norm_key(k): v.strip() for k, v in row.items()}
            if int(float(normed["label"])) != 1:
                continue
            pairs.append((int(normed["enzyme_index"]), int(normed["substrate_index"])))
    return pairs


def _load_enzyme_pdb_map(enzymes_csv: Path) -> dict[int, str]:
    """Map Enzyme_Index → PDB_ID from Enzymes.csv."""
    mapping: dict[int, str] = {}
    with enzymes_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            normed = {_norm_key(k): v.strip() for k, v in row.items()}
            mapping[int(normed["enzyme_index"])] = normed.get("pdb_id", "")
    return mapping


def _sample_negatives(
    positives: list[tuple[int, int]],
    all_enzyme_indices: list[int],
    ratio: float,
    seed: int,
) -> list[tuple[int, int]]:
    """Single-direction sampling: for each substrate, pick random 'fake' enzymes."""
    rng = random.Random(seed)

    true_enzymes: dict[int, set[int]] = defaultdict(set)
    pos_per_sub: Counter[int] = Counter()
    for enz, sub in positives:
        true_enzymes[sub].add(enz)
        pos_per_sub[sub] += 1

    selected: set[tuple[int, int]] = set()
    for sub, n_pos in pos_per_sub.items():
        candidates = [e for e in all_enzyme_indices if e not in true_enzymes[sub]]
        if not candidates:
            continue
        n_neg = min(int(round(n_pos * ratio)), len(candidates))
        picked = rng.sample(candidates, n_neg) if n_neg < len(candidates) else candidates
        selected.update((e, sub) for e in picked)

    # Adjust to global target ratio
    desired = min(int(round(len(positives) * ratio)),
                  sum(len(all_enzyme_indices) - len(true_enzymes[s])
                      for s in pos_per_sub))
    if len(selected) < desired:
        universe = [
            (e, s) for s in pos_per_sub for e in all_enzyme_indices
            if e not in true_enzymes[s] and (e, s) not in selected
        ]
        selected.update(rng.sample(universe, min(desired - len(selected), len(universe))))
    elif len(selected) > desired:
        selected = set(rng.sample(list(selected), desired))

    return sorted(selected, key=lambda x: (x[1], x[0]))


def generate_negative_pairs(
    data_csv: Path,
    enzymes_csv: Path,
    output_csv: Path,
    ratio: float = TARGET_NEG_RATIO,
    random_seed: int = 2026,
    dock_index_start: int = 1000,
) -> NegativeSamplingResult:
    """Generate random negative pairs and save to CSV.

    Args:
        data_csv: B6 data.csv with columns [Dock Index, Enzyme Index, Substrate Index, Label].
        enzymes_csv: B6 Enzymes.csv with columns [Enzyme_Index, ..., PDB_ID, ...].
        output_csv: Output path for negative_pairs.csv.
        ratio: Target negative/positive ratio (default 9.36 matching ESIBank).
        random_seed: For reproducibility.
        dock_index_start: Starting Dock_Index for negatives (avoid collision with 0-515).
    """
    result = NegativeSamplingResult(output_path=output_csv, success=False)

    if not data_csv.exists():
        result.error = f"data.csv not found: {data_csv}"
        return result
    if not enzymes_csv.exists():
        result.error = f"Enzymes.csv not found: {enzymes_csv}"
        return result

    try:
        positives = _load_positive_pairs(data_csv)
        enzyme_pdb = _load_enzyme_pdb_map(enzymes_csv)
        all_enzymes = sorted(enzyme_pdb.keys())

        if not positives:
            result.error = "No positive pairs (Label==1) found"
            return result
        if not all_enzymes:
            result.error = "No enzymes found in Enzymes.csv"
            return result

        negatives = _sample_negatives(positives, all_enzymes, ratio, random_seed)

        output_csv.parent.mkdir(parents=True, exist_ok=True)
        with output_csv.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Dock_Index", "PDB_ID", "Enzyme_Index", "Substrate_Index", "Label"])
            for i, (enz, sub) in enumerate(negatives):
                w.writerow([dock_index_start + i, enzyme_pdb.get(enz, ""), enz, sub, 0])

        result.positive_count = len(positives)
        result.negative_count = len(negatives)
        result.achieved_ratio = result.negative_count / result.positive_count if positives else 0
        result.success = True
        return result

    except Exception as exc:
        result.error = f"Negative sampling failed: {exc}"
        return result
