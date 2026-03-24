#!/usr/bin/env python3
"""
Phase 6: Download AlphaFold raw structures + batch heme transplantation.

For 540 enzymes that have UniProt but failed AlphaFill:
1. Download raw AlphaFold PDB (no heme)
2. Find best heme template from 951 existing heme-containing structures
3. Structural alignment (BioPython Superimposer) → copy heme → save merged PDB
4. Update receptor_manifest.csv

Design (Codex review):
- Template selection by sequence identity (k-mer based fast screening)
- BioPython Superimposer for Cα alignment
- QC: RMSD < 3Å = good, 3-5Å = warn, >5Å = flag
- ThreadPoolExecutor for downloads, sequential for transplantation
- Resume support, incremental manifest update
"""
from __future__ import annotations

import argparse, csv, hashlib, os, re, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import requests

SCRIPT = Path(__file__).resolve()
PROJECT = SCRIPT.parents[2]
COMBINED = PROJECT / "data" / "combined"
STRUCTURES = PROJECT / "data" / "structures"
MANIFEST_PATH = STRUCTURES / "manifests" / "receptor_manifest.csv"
ALPHAFOLD_PDB_DIR = STRUCTURES / "alphafold" / "pdb"
HEME_TRANSPLANT_DIR = STRUCTURES / "heme_transplant" / "pdb"

ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v6.pdb"
USER_AGENT = "P450-AlphaFold-HemeTransplant/1.0"
DOWNLOAD_INTERVAL = 0.5  # seconds between downloads
DOWNLOAD_TIMEOUT = 60

VALID_UNIPROT = re.compile(r"^[A-Z][A-Z0-9]{4,9}$")
HEME_NAMES = {"HEM", "HEC", "HEA", "HEO"}
GOOD_RMSD = 3.0
WARN_RMSD = 5.0

def log(msg): print(msg, file=sys.stderr, flush=True)
def norm(v):
    t = str(v).strip() if v is not None else ""
    return "" if t.lower() in {"","na","n/a","none","null"} else t

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
def load_manifest():
    rows = []
    with MANIFEST_PATH.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    return rows

EXTRA_FIELDS = [
    "alphafold_pdb_path", "heme_transplant_pdb_path",
    "heme_fe_x", "heme_fe_y", "heme_fe_z",
]

def save_manifest(rows):
    if not rows: return
    # Collect ALL fields from all rows (handles new fields added during processing)
    all_fields = []
    seen = set()
    for r in rows:
        for k in r:
            if k not in seen:
                all_fields.append(k)
                seen.add(k)
    # Ensure extra fields are included
    for f in EXTRA_FIELDS:
        if f not in seen:
            all_fields.append(f)
            seen.add(f)
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST_PATH.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=all_fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in all_fields})

def load_enzyme_sequences():
    """Load global_enzyme_id → sequence mapping."""
    seqs = {}
    with (COMBINED / "global_enzymes.csv").open("r", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            seqs[r["global_enzyme_id"]] = norm(r.get("canonical_sequence",""))
    return seqs

# ---------------------------------------------------------------------------
# Template library: build from existing heme-containing structures
# ---------------------------------------------------------------------------
def build_template_library(manifest_rows, enzyme_sequences):
    """Build template library from enzymes that already have heme."""
    templates = []
    heme_statuses = {"experimental_pdb", "pcpd_predicted", "downloaded_with_heme"}
    for r in manifest_rows:
        if r.get("status") not in heme_statuses: continue
        if r.get("has_heme") != "1": continue
        gid = r["global_enzyme_id"]
        seq = enzyme_sequences.get(gid, "")
        if not seq: continue

        # Determine structure file path
        structure_path = ""
        if r.get("alphafill_cif_path"):
            structure_path = str(PROJECT / r["alphafill_cif_path"])
        elif r.get("existing_pdb_path"):
            # These are relative paths like "s1_rcsb/5NCB" or "pcpd/CYP2A6.pdb"
            # Need to resolve based on source
            pass  # Will handle in find_best_template

        templates.append({
            "global_enzyme_id": gid,
            "uniprot": r.get("canonical_uniprot_id",""),
            "sequence": seq,
            "status": r["status"],
            "alphafill_cif": r.get("alphafill_cif_path",""),
            "existing_pdb": r.get("existing_pdb_path",""),
        })
    return templates

def sequence_identity(seq1, seq2):
    """Quick sequence identity estimate using k-mer overlap (k=3)."""
    if not seq1 or not seq2: return 0.0
    k = 3
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1)-k+1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2)-k+1))
    if not kmers1 or not kmers2: return 0.0
    return len(kmers1 & kmers2) / max(len(kmers1), len(kmers2))

def find_best_template(target_seq, templates, top_n=5):
    """Find the best heme template by sequence similarity."""
    scored = []
    for t in templates:
        sim = sequence_identity(target_seq, t["sequence"])
        scored.append((sim, t))
    scored.sort(key=lambda x: -x[0])
    return scored[:top_n]

# ---------------------------------------------------------------------------
# Download AlphaFold PDB
# ---------------------------------------------------------------------------
def download_alphafold_pdb(uniprot_id, sess):
    """Download raw AlphaFold PDB. Returns (path, status)."""
    url = ALPHAFOLD_URL.format(uniprot=uniprot_id)
    out_path = ALPHAFOLD_PDB_DIR / f"AF-{uniprot_id}-F1.pdb"

    if out_path.exists() and out_path.stat().st_size > 0:
        return str(out_path), "already_downloaded"

    try:
        resp = sess.get(url, timeout=DOWNLOAD_TIMEOUT)
        if resp.status_code == 404:
            return "", "alphafold_not_found"
        if resp.status_code == 429:
            return "", "rate_limited"
        resp.raise_for_status()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(resp.content)
        return str(out_path), "downloaded"
    except Exception as e:
        return "", f"error:{str(e)[:100]}"

# ---------------------------------------------------------------------------
# Heme transplantation using BioPython
# ---------------------------------------------------------------------------
def transplant_heme(target_pdb_path, template_structure_path, output_path):
    """
    Align target to template, copy heme residues, save merged PDB.
    Returns (rmsd, n_matched_ca, fe_coords, status).
    """
    from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Superimposer

    # Parse target
    parser_pdb = PDBParser(QUIET=True)
    target_struct = parser_pdb.get_structure("target", target_pdb_path)

    # Parse template (could be CIF or PDB)
    if template_structure_path.endswith(".cif"):
        parser_cif = MMCIFParser(QUIET=True)
        template_struct = parser_cif.get_structure("template", template_structure_path)
    else:
        template_struct = parser_pdb.get_structure("template", template_structure_path)

    # Get Cα atoms from both
    target_model = target_struct[0]
    template_model = template_struct[0]

    def get_ca_atoms(model):
        cas = {}
        for chain in model:
            for res in chain:
                if res.id[0] != " ": continue  # skip hetero
                if "CA" in res:
                    cas[res.id[1]] = res["CA"]
        return cas

    target_cas = get_ca_atoms(target_model)
    template_cas = get_ca_atoms(template_model)

    # Find common residue positions
    common = sorted(set(target_cas.keys()) & set(template_cas.keys()))
    if len(common) < 50:
        return None, len(common), None, "too_few_common_residues"

    # Align
    target_atoms = [target_cas[i] for i in common]
    template_atoms = [template_cas[i] for i in common]

    sup = Superimposer()
    sup.set_atoms(target_atoms, template_atoms)
    sup.apply(template_model.get_atoms())
    rmsd = sup.rms

    # Find heme residues in template
    heme_residues = []
    for chain in template_model:
        for res in chain:
            if res.get_resname().strip() in HEME_NAMES:
                heme_residues.append(res)

    if not heme_residues:
        return rmsd, len(common), None, "no_heme_in_template"

    # Find Fe coordinates
    fe_coords = None
    for res in heme_residues:
        for atom in res:
            if atom.get_name().strip().upper() == "FE":
                fe_coords = tuple(atom.get_vector().get_array())
                break

    # Copy heme to target: add to first chain
    target_chain = list(target_model.get_chains())[0]
    max_resid = max(r.id[1] for r in target_chain if r.id[0] == " ")

    import copy
    for i, heme_res in enumerate(heme_residues):
        new_res = copy.deepcopy(heme_res)
        new_id = (" ", max_resid + 100 + i, " ")
        new_res.id = new_id
        # Detach from old parent
        new_res.detach_parent()
        target_chain.add(new_res)

    # Save
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(target_struct)
    io.save(str(output_path))

    qc = "good" if rmsd < GOOD_RMSD else ("warn" if rmsd < WARN_RMSD else "bad")
    return rmsd, len(common), fe_coords, qc

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-workers", type=int, default=4)
    args = parser.parse_args()

    log("[setup] Loading data...")
    manifest = load_manifest()
    enzyme_seqs = load_enzyme_sequences()

    # Identify targets: alphafill_not_found or no_heme, with valid UniProt
    targets = []
    for i, r in enumerate(manifest):
        if r.get("status") not in ("alphafill_not_found", "no_heme"): continue
        uid = norm(r.get("canonical_uniprot_id",""))
        if not uid or not VALID_UNIPROT.match(uid): continue
        # Skip if already has heme_transplant
        if r.get("heme_transplant_pdb_path", "").strip(): continue
        targets.append((i, r))

    log(f"[setup] {len(targets)} enzymes need AlphaFold download + heme transplant")

    # Build template library
    templates = build_template_library(manifest, enzyme_seqs)
    log(f"[templates] {len(templates)} heme-containing structures available as templates")

    if not templates:
        log("ERROR: No heme templates available!")
        return

    # Create dirs
    ALPHAFOLD_PDB_DIR.mkdir(parents=True, exist_ok=True)
    HEME_TRANSPLANT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Batch download AlphaFold PDBs
    log(f"\n[download] Downloading {len(targets)} AlphaFold structures ({args.download_workers} workers)...")
    sess = requests.Session()
    sess.headers["User-Agent"] = USER_AGENT

    download_results = {}  # manifest_idx → (path, status)
    dl_ok = dl_fail = 0

    with ThreadPoolExecutor(max_workers=args.download_workers) as executor:
        futures = {}
        for idx, r in targets:
            uid = norm(r["canonical_uniprot_id"])
            futures[executor.submit(download_alphafold_pdb, uid, sess)] = (idx, uid)
            time.sleep(DOWNLOAD_INTERVAL / args.download_workers)

        for future in as_completed(futures):
            idx, uid = futures[future]
            try:
                path, status = future.result(timeout=120)
            except Exception as e:
                path, status = "", f"error:{e}"

            download_results[idx] = (path, status)
            if "downloaded" in status or "already" in status:
                dl_ok += 1
            else:
                dl_fail += 1

            if (dl_ok + dl_fail) % 50 == 0:
                log(f"  [download] {dl_ok+dl_fail}/{len(targets)}: {dl_ok} ok, {dl_fail} fail")

    log(f"[download] Complete: {dl_ok} ok, {dl_fail} fail")

    # Step 2: Transplant heme for successfully downloaded structures
    log(f"\n[transplant] Starting heme transplantation...")
    tp_ok = tp_fail = 0

    for idx, r in targets:
        gid = r["global_enzyme_id"]
        dl_path, dl_status = download_results.get(idx, ("", "not_downloaded"))

        if not dl_path or "error" in dl_status or dl_status == "alphafold_not_found":
            # Update manifest with download failure
            manifest[idx]["alphafold_pdb_path"] = ""
            manifest[idx]["notes"] = f"AlphaFold download: {dl_status}"
            tp_fail += 1
            continue

        # Store download path
        manifest[idx]["alphafold_pdb_path"] = str(Path(dl_path).relative_to(PROJECT))

        # Find best template
        target_seq = enzyme_seqs.get(gid, "")
        best_templates = find_best_template(target_seq, templates)

        if not best_templates:
            manifest[idx]["notes"] = "No suitable template found"
            tp_fail += 1
            continue

        # Try templates in order until one works
        success = False
        for sim, tmpl in best_templates:
            # Resolve template structure path
            tmpl_path = ""
            if tmpl["alphafill_cif"]:
                p = PROJECT / tmpl["alphafill_cif"]
                if p.exists(): tmpl_path = str(p)
            if not tmpl_path and tmpl["existing_pdb"]:
                # Try PCPD path
                pcpd_dir = PROJECT / "downloads" / "PCPD" / "PDB"
                ep = tmpl["existing_pdb"]
                if "pcpd" in ep.lower():
                    fname = ep.split("/")[-1]
                    if not fname.endswith(".pdb"): fname += ".pdb"
                    p = pcpd_dir / fname
                    if p.exists(): tmpl_path = str(p)

            if not tmpl_path:
                continue

            output_path = HEME_TRANSPLANT_DIR / f"{gid}_heme.pdb"

            try:
                rmsd, n_ca, fe_coords, qc = transplant_heme(dl_path, tmpl_path, str(output_path))
            except Exception as e:
                continue

            if rmsd is None:
                continue

            # Success!
            manifest[idx]["status"] = "heme_transplanted"
            manifest[idx]["has_heme"] = "1"
            manifest[idx]["structure_source"] = "alphafold_heme_transplant"
            manifest[idx]["heme_transplant_pdb_path"] = str(output_path.relative_to(PROJECT))
            manifest[idx]["template_pdb"] = tmpl["global_enzyme_id"]
            manifest[idx]["sequence_identity"] = f"{sim:.3f}"
            manifest[idx]["global_rmsd"] = f"{rmsd:.2f}"
            manifest[idx]["notes"] = f"qc={qc}, n_ca={n_ca}, template={tmpl['uniprot']}"
            if fe_coords:
                manifest[idx]["heme_fe_x"] = f"{fe_coords[0]:.3f}"
                manifest[idx]["heme_fe_y"] = f"{fe_coords[1]:.3f}"
                manifest[idx]["heme_fe_z"] = f"{fe_coords[2]:.3f}"
            success = True
            tp_ok += 1
            break

        if not success:
            manifest[idx]["notes"] = f"All {len(best_templates)} templates failed"
            tp_fail += 1

        if (tp_ok + tp_fail) % 20 == 0:
            save_manifest(manifest)
            log(f"  [transplant] {tp_ok+tp_fail}/{len(targets)}: {tp_ok} ok, {tp_fail} fail (checkpoint)")

    # Final save
    save_manifest(manifest)

    log(f"\n{'='*60}")
    log(f"[DONE] Download: {dl_ok} ok, {dl_fail} fail")
    log(f"[DONE] Transplant: {tp_ok} ok, {tp_fail} fail")
    log(f"[DONE] Manifest updated: {MANIFEST_PATH}")
    log(f"{'='*60}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted"); raise SystemExit(130)
    except Exception as e:
        log(f"ERROR: {e}"); import traceback; traceback.print_exc(); raise SystemExit(1)
