"""
Phase 6 Stage A: Build pair_registry_master.csv
- Maps 52,254 enzyme-substrate pairs to structure readiness
- Assigns each pair to batch_1 / batch_2 / blocked
- Creates enzyme_registry, substrate_registry, and pair_registry_master

Codex-reviewed: 2026-03-25
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

BASE = Path(__file__).resolve().parent.parent.parent
DATA = BASE / "data"
P450 = DATA / "P450_Family"
STRUCTS = DATA / "structures" / "manifests"
OUT = DATA / "registries"
OUT.mkdir(parents=True, exist_ok=True)

BATCH1_STATUSES = {"downloaded_with_heme", "heme_transplanted", "pcpd_predicted", "experimental_pdb"}
# All 219 enzymes (131 needs_colabfold + 54 heme_transplant_failed + 26 no_template + 8 alphafill_not_found)
# are uploaded to ColabFold and currently being predicted. They are ALL batch_2, not blocked.
BATCH2_STATUSES = {"needs_colabfold", "heme_transplant_failed", "no_template", "alphafill_not_found"}
BLOCKED_STATUSES = set()  # no permanently blocked enzymes for now; reassess after ColabFold finishes


def build_enzyme_registry():
    """Build enzyme_registry.csv: one row per enzyme_index (0-1621) with structure status."""
    print("\n=== Building enzyme_registry ===")
    enzymes = pd.read_csv(P450 / "Enzymes.csv")
    manifest = pd.read_csv(STRUCTS / "receptor_manifest.csv")

    # Validate manifest uniqueness
    assert manifest["global_enzyme_id"].is_unique, "global_enzyme_id not unique in manifest!"
    assert manifest["canonical_uniprot_id"].dropna().is_unique, "canonical_uniprot_id not unique!"

    enzymes["enzyme_index"] = enzymes.index
    enzymes_id_col = "uniprots"

    # Step 1: join on canonical_uniprot_id (real UniProt IDs)
    manifest_cols = ["global_enzyme_id", "canonical_uniprot_id", "status",
                     "has_heme", "quality_tier", "existing_pdb_path",
                     "alphafill_cif_path", "heme_transplant_pdb_path",
                     "alphafold_pdb_path", "heme_fe_x", "heme_fe_y", "heme_fe_z"]
    manifest_sub = manifest[manifest_cols].copy()

    # First join: uniprots -> canonical_uniprot_id
    real_uniprot_mask = ~enzymes[enzymes_id_col].str.startswith("ENZ_G", na=False)
    enz_real = enzymes[real_uniprot_mask].copy()
    enz_ensg = enzymes[~real_uniprot_mask].copy()

    merged1 = enz_real.merge(manifest_sub, left_on=enzymes_id_col,
                             right_on="canonical_uniprot_id", how="left")
    matched1 = merged1["status"].notna().sum()
    unmatched1 = merged1["status"].isna().sum()
    print(f"  Step 1 (UniProt join): {matched1} matched, {unmatched1} unmatched")

    # Second join: uniprots (ENZ_G*) -> global_enzyme_id
    merged2 = enz_ensg.merge(manifest_sub, left_on=enzymes_id_col,
                             right_on="global_enzyme_id", how="left")
    matched2 = merged2["status"].notna().sum()
    unmatched2 = merged2["status"].isna().sum()
    print(f"  Step 2 (ENZ_G join): {matched2} matched, {unmatched2} unmatched")

    # Combine
    registry = pd.concat([merged1, merged2], ignore_index=True).sort_values("enzyme_index").reset_index(drop=True)

    # Validate
    assert len(registry) == 1622, f"Expected 1622 rows, got {len(registry)}"
    assert registry["enzyme_index"].is_unique, "Duplicate enzyme_index!"
    assert (registry["enzyme_index"] == np.arange(1622)).all(), "enzyme_index not 0-1621!"
    assert registry["status"].notna().all(), f"{registry['status'].isna().sum()} enzymes have no manifest match!"

    # Assign enzyme_batch
    def assign_batch(status):
        if status in BATCH1_STATUSES:
            return "batch_1"
        elif status in BATCH2_STATUSES:
            return "batch_2"
        elif status in BLOCKED_STATUSES:
            return "blocked"
        else:
            return "unknown"

    registry["enzyme_batch"] = registry["status"].apply(assign_batch)
    assert (registry["enzyme_batch"] != "unknown").all(), f"Unknown statuses: {registry[registry['enzyme_batch']=='unknown']['status'].unique()}"

    # Select output columns
    out_cols = ["enzyme_index", "uniprots", "global_enzyme_id", "canonical_uniprot_id",
                "status", "enzyme_batch", "has_heme", "quality_tier",
                "existing_pdb_path", "alphafill_cif_path", "heme_transplant_pdb_path",
                "alphafold_pdb_path", "heme_fe_x", "heme_fe_y", "heme_fe_z"]
    registry = registry[[c for c in out_cols if c in registry.columns]]
    registry.rename(columns={"status": "receptor_status", "uniprots": "enzyme_input_id"}, inplace=True)

    registry.to_csv(OUT / "enzyme_registry.csv", index=False)
    print(f"  Saved: {OUT / 'enzyme_registry.csv'}")
    print(f"  Batch 1: {(registry['enzyme_batch']=='batch_1').sum()}")
    print(f"  Batch 2: {(registry['enzyme_batch']=='batch_2').sum()}")
    print(f"  Blocked: {(registry['enzyme_batch']=='blocked').sum()}")
    return registry


def build_substrate_registry():
    """Build substrate_registry.csv: one row per substrate_index (0-2124) with 3D status."""
    print("\n=== Building substrate_registry ===")
    substrates = pd.read_csv(P450 / "Substrates.csv")
    lig_manifest = pd.read_csv(STRUCTS / "ligand_manifest.csv")

    substrates["substrate_index"] = substrates.index

    # Drop NaN compound_id rows from manifest
    lig_clean = lig_manifest.dropna(subset=["global_compound_id"]).copy()
    print(f"  ligand_manifest: {len(lig_manifest)} rows, after dropping NaN: {len(lig_clean)}")

    # Join by SMILES (normalized)
    # Substrates.csv has 'Substrate_SMILES', ligand_manifest has 'canonical_smiles'
    merged = substrates.merge(
        lig_clean[["global_compound_id", "canonical_smiles", "status", "sdf_path",
                    "num_atoms", "final_energy", "force_field_used"]],
        left_on="Substrate_SMILES", right_on="canonical_smiles", how="left"
    )

    matched = merged["status"].notna().sum()
    unmatched = merged["status"].isna().sum()
    print(f"  SMILES join: {matched} matched, {unmatched} unmatched")

    # For unmatched, check if it's the wildcard substrates
    if unmatched > 0:
        unmatched_smiles = merged[merged["status"].isna()]["Substrate_SMILES"].tolist()
        print(f"  Unmatched SMILES samples: {unmatched_smiles[:5]}")
        # Mark as blocked
        merged.loc[merged["status"].isna(), "status"] = "no_3d"

    # Validate
    assert len(merged) == 2125, f"Expected 2125 rows, got {len(merged)}"
    assert merged["substrate_index"].is_unique, "Duplicate substrate_index!"

    # Assign ligand_batch
    merged["ligand_status"] = merged["status"].apply(lambda s: "ok" if s == "ok" else "failed")

    out_cols = ["substrate_index", "Substrate_SMILES", "global_compound_id",
                "status", "ligand_status", "sdf_path", "num_atoms"]
    merged = merged[[c for c in out_cols if c in merged.columns]]
    merged.rename(columns={"status": "ligand_3d_status"}, inplace=True)

    merged.to_csv(OUT / "substrate_registry.csv", index=False)
    print(f"  Saved: {OUT / 'substrate_registry.csv'}")
    print(f"  OK: {(merged['ligand_status']=='ok').sum()}")
    print(f"  Failed: {(merged['ligand_status']=='failed').sum()}")
    return merged


def build_pair_registry(enz_reg, sub_reg):
    """Build pair_registry_master.csv from data.csv + enzyme/substrate registries."""
    print("\n=== Building pair_registry_master ===")
    data = pd.read_csv(P450 / "random_split" / "data.csv")

    # Freeze dock_index as explicit column
    data["dock_index"] = np.arange(len(data))
    data.rename(columns={"enzyme": "enzyme_index", "reaction": "substrate_index"}, inplace=True)

    # Validate input
    assert len(data) == 52254, f"Expected 52254 rows, got {len(data)}"
    assert data["dock_index"].is_unique
    assert data["enzyme_index"].between(0, 1621).all()
    assert data["substrate_index"].between(0, 2124).all()
    assert set(data["label"].unique()) == {0, 1}

    # Join enzyme registry
    enz_join_cols = ["enzyme_index", "enzyme_input_id", "global_enzyme_id",
                     "receptor_status", "enzyme_batch",
                     "heme_fe_x", "heme_fe_y", "heme_fe_z"]
    enz_join = enz_reg[[c for c in enz_join_cols if c in enz_reg.columns]]
    pairs = data.merge(enz_join, on="enzyme_index", how="left")
    assert len(pairs) == 52254, f"Row count changed after enzyme join: {len(pairs)}"

    # Join substrate registry
    sub_join_cols = ["substrate_index", "Substrate_SMILES", "ligand_3d_status", "ligand_status"]
    sub_join = sub_reg[[c for c in sub_join_cols if c in sub_reg.columns]]
    pairs = pairs.merge(sub_join, on="substrate_index", how="left")
    assert len(pairs) == 52254, f"Row count changed after substrate join: {len(pairs)}"

    # Assign final batch
    def assign_pair_batch(row):
        if row.get("ligand_status") == "failed":
            return "blocked"
        eb = row.get("enzyme_batch", "unknown")
        if eb == "blocked":
            return "blocked"
        elif eb == "batch_2":
            return "batch_2"
        elif eb == "batch_1":
            return "batch_1"
        return "unknown"

    def assign_blocked_reason(row):
        reasons = []
        if row.get("ligand_status") == "failed":
            reasons.append(f"substrate_3d_failed({row.get('ligand_3d_status','')})")
        if row.get("enzyme_batch") == "blocked":
            reasons.append(f"enzyme_structure_failed({row.get('receptor_status','')})")
        return "; ".join(reasons) if reasons else ""

    pairs["batch"] = pairs.apply(assign_pair_batch, axis=1)
    pairs["blocked_reason"] = pairs.apply(assign_blocked_reason, axis=1)
    pairs["pair_status"] = "pending"

    # Validate final
    assert len(pairs) == 52254
    assert pairs["dock_index"].is_unique
    assert (pairs["dock_index"] == np.arange(52254)).all()
    assert pairs["batch"].notna().all()
    assert (pairs["batch"] != "unknown").all(), f"Unknown batches: {(pairs['batch']=='unknown').sum()}"

    batch_counts = pairs["batch"].value_counts()
    print(f"  Batch distribution:")
    for b, c in batch_counts.items():
        print(f"    {b}: {c} pairs")
    assert batch_counts.sum() == 52254

    # Save master
    pairs.to_csv(OUT / "pair_registry_master.csv", index=False)
    print(f"  Saved: {OUT / 'pair_registry_master.csv'}")

    # Save splits
    for batch_name in ["batch_1", "batch_2", "blocked"]:
        subset = pairs[pairs["batch"] == batch_name]
        fname = f"pair_registry_{batch_name}.csv"
        subset.to_csv(OUT / fname, index=False)
        print(f"  Saved: {OUT / fname} ({len(subset)} rows)")

    return pairs


def main():
    print("=" * 60)
    print("Phase 6 Stage A: Build Pair Registry")
    print("=" * 60)

    enz_reg = build_enzyme_registry()
    sub_reg = build_substrate_registry()
    pairs = build_pair_registry(enz_reg, sub_reg)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Enzyme registry:    {len(enz_reg)} enzymes")
    print(f"Substrate registry: {len(sub_reg)} substrates")
    print(f"Pair registry:      {len(pairs)} pairs")
    print(f"  Batch 1:          {(pairs['batch']=='batch_1').sum()}")
    print(f"  Batch 2:          {(pairs['batch']=='batch_2').sum()}")
    print(f"  Blocked:          {(pairs['batch']=='blocked').sum()}")
    print(f"\nOutput directory: {OUT}")


if __name__ == "__main__":
    main()
