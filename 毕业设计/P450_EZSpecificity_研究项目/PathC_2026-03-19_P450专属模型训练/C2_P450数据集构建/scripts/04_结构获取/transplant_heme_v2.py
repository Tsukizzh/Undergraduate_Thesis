#!/usr/bin/env python3
"""
Phase 6 v2: Heme transplantation with sequence-alignment-driven Cα pairing.

v1 bugs fixed:
1. Residue-number matching → sequence alignment matching (BioPython PairwiseAligner)
2. k-mer template selection → real pairwise sequence identity
3. Experimental PDB templates now included (downloaded from RCSB on demand)

Quality tiers:
  Tier 1: RMSD < 2Å, identity > 50%  → high confidence
  Tier 2: RMSD < 3Å, identity > 30%  → good
  Tier 3: RMSD < 5Å, identity > 25%  → acceptable
  Tier 4: RMSD >= 5Å                  → reject, try next template

Fallback: AlphaFill API upload for enzymes where all local templates fail.
"""
from __future__ import annotations

import argparse, copy, csv, json, os, re, sys, tempfile, time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import requests
from Bio import Align
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa

SCRIPT = Path(__file__).resolve()
PROJECT = SCRIPT.parents[2]
COMBINED = PROJECT / "data" / "combined"
STRUCTURES = PROJECT / "data" / "structures"
DOWNLOADS = PROJECT / "downloads"

MANIFEST_PATH = STRUCTURES / "manifests" / "receptor_manifest.csv"
ALPHAFOLD_PDB_DIR = STRUCTURES / "alphafold" / "pdb"
HEME_DIR = STRUCTURES / "heme_transplant" / "pdb"
ALPHAFILL_CIF_DIR = STRUCTURES / "alphafill" / "cif"
PCPD_PDB_DIR = DOWNLOADS / "PCPD" / "PDB"
RCSB_CACHE = DOWNLOADS / "RCSB" / "PDB"

ALPHAFILL_API = "https://alphafill.eu/v1/aff"
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
USER_AGENT = "P450-Heme-v2/1.0"
TIMEOUT = 120

VALID_UNIPROT = re.compile(r"^[A-Z][A-Z0-9]{4,9}$")
HEME_NAMES = {"HEM", "HEC", "HEA", "HEB", "HEO"}
MIN_IDENTITY = 0.25
MIN_ALIGNED = 50

def log(msg): print(msg, file=sys.stderr, flush=True)
def norm(v):
    t = str(v).strip() if v is not None else ""
    return "" if t.lower() in {"","na","n/a","none","null"} else t

# ---------------------------------------------------------------------------
# Sequence alignment
# ---------------------------------------------------------------------------
def make_aligner():
    a = Align.PairwiseAligner(mode="global")
    a.match_score = 2.0; a.mismatch_score = 0.0
    a.open_gap_score = -10.0; a.extend_gap_score = -0.5
    return a

def align_sequences(seq_a, seq_b, aligner):
    """Returns (identity, aligned_pairs_indices, n_aligned)."""
    aln = aligner.align(seq_a, seq_b)[0]
    pairs = []
    t_blocks, q_blocks = aln.aligned
    for tb, qb in zip(t_blocks, q_blocks):
        for off in range(min(int(tb[1]-tb[0]), int(qb[1]-qb[0]))):
            pairs.append((int(tb[0])+off, int(qb[0])+off))
    matches = sum(1 for i, j in pairs if seq_a[i] == seq_b[j])
    identity = matches / len(pairs) if pairs else 0.0
    return identity, pairs, len(pairs)

# ---------------------------------------------------------------------------
# Structure parsing helpers
# ---------------------------------------------------------------------------
def parse_structure(path):
    if str(path).endswith(".cif"):
        return MMCIFParser(QUIET=True).get_structure("s", str(path))
    return PDBParser(QUIET=True).get_structure("s", str(path))

def extract_chain_ca(model):
    """Returns list of (chain_id, sequence_str, list_of_residues_with_CA)."""
    chains = []
    for chain in model:
        seq, residues = [], []
        for res in chain:
            if res.id[0] != " ": continue
            rn = res.get_resname().strip().upper()
            aa = protein_letters_3to1_extended.get(rn)
            if not aa: continue
            if "CA" not in res: continue
            seq.append(aa); residues.append(res)
        if residues:
            chains.append((chain.id, "".join(seq), residues))
    return chains

def find_heme(model):
    """Returns list of heme residues across all chains."""
    hemes = []
    for chain in model:
        for res in chain:
            if res.get_resname().strip().upper() in HEME_NAMES:
                hemes.append(res)
    return hemes

def get_fe_coords(residue):
    for atom in residue:
        if atom.get_name().strip().upper().startswith("FE"):
            return tuple(float(x) for x in atom.coord)
    return None

# ---------------------------------------------------------------------------
# Template library
# ---------------------------------------------------------------------------
def build_template_library(manifest, enzyme_seqs):
    """Build list of templates from all heme-containing enzymes."""
    templates = []
    for r in manifest:
        if r.get("has_heme") != "1": continue
        if r.get("status") not in ("experimental_pdb","pcpd_predicted","downloaded_with_heme"): continue
        gid = r["global_enzyme_id"]
        seq = enzyme_seqs.get(gid, "").upper()
        if not seq or len(seq) < 100: continue

        # Resolve structure file path
        path = ""
        source = ""

        # 1. AlphaFill CIF
        af_cif = norm(r.get("alphafill_cif_path",""))
        if af_cif:
            p = PROJECT / af_cif
            if p.exists(): path, source = str(p), "alphafill_cif"

        # 2. PCPD PDB
        if not path:
            ep = norm(r.get("existing_pdb_path",""))
            if ep and "pcpd" in ep.lower():
                fname = Path(ep).name
                if not fname.endswith(".pdb"): fname += ".pdb"
                p = PCPD_PDB_DIR / fname
                if p.exists(): path, source = str(p), "pcpd_pdb"

        # 3. S1 experimental PDB (need PDB ID to download from RCSB)
        if not path:
            ep = norm(r.get("existing_pdb_path",""))
            if ep and "s1_rcsb" in ep.lower():
                pdb_id = Path(ep).stem.upper()
                if len(pdb_id) == 4:  # valid PDB ID
                    source = "experimental_pdb"
                    # Will download on demand
                    path = f"RCSB:{pdb_id}"

        if not path: continue
        templates.append({
            "gid": gid, "uid": norm(r.get("canonical_uniprot_id","")),
            "seq": seq, "path": path, "source": source,
        })
    return templates

def ensure_rcsb_pdb(pdb_id, sess):
    """Download PDB from RCSB if not cached."""
    RCSB_CACHE.mkdir(parents=True, exist_ok=True)
    out = RCSB_CACHE / f"{pdb_id}.pdb"
    if out.exists() and out.stat().st_size > 0: return out
    resp = sess.get(RCSB_URL.format(pdb_id=pdb_id), timeout=TIMEOUT)
    resp.raise_for_status()
    out.write_bytes(resp.content)
    return out

def resolve_template_path(tmpl, sess):
    """Resolve template structure file path, downloading RCSB if needed."""
    path = tmpl["path"]
    if path.startswith("RCSB:"):
        pdb_id = path[5:]
        return str(ensure_rcsb_pdb(pdb_id, sess))
    return path

# ---------------------------------------------------------------------------
# Core: transplant heme using sequence-alignment-driven Cα pairing
# ---------------------------------------------------------------------------
def transplant_heme(target_pdb, template_path, target_seq, template_seq, aligner):
    """
    1. Sequence-align target and template
    2. Extract Cα atoms from aligned residue pairs
    3. Superimpose template onto target
    4. Copy heme from aligned template to target
    5. Return (rmsd, n_ca, fe_coords, quality_tier) or raise on failure
    """
    target_struct = parse_structure(Path(target_pdb))
    template_struct = parse_structure(Path(template_path))

    target_model = target_struct[0]
    template_model = template_struct[0]

    # Get chains with CA atoms
    target_chains = extract_chain_ca(target_model)
    template_chains = extract_chain_ca(template_model)
    if not target_chains or not template_chains:
        raise ValueError("No protein chains found")

    # Find best chain pair by alignment
    best = None
    for t_cid, t_seq, t_res in target_chains:
        for q_cid, q_seq, q_res in template_chains:
            identity, pairs, n_aln = align_sequences(t_seq, q_seq, aligner)
            if n_aln < MIN_ALIGNED: continue
            score = (identity, n_aln)
            if best is None or score > best["score"]:
                best = {"score": score, "identity": identity, "pairs": pairs,
                        "t_cid": t_cid, "t_res": t_res, "q_cid": q_cid, "q_res": q_res}

    if best is None:
        raise ValueError(f"No chain pair with >= {MIN_ALIGNED} aligned Cα")

    # Build Cα atom lists from aligned pairs
    fixed_atoms, moving_atoms = [], []
    for t_idx, q_idx in best["pairs"]:
        if t_idx >= len(best["t_res"]) or q_idx >= len(best["q_res"]): continue
        t_res = best["t_res"][t_idx]
        q_res = best["q_res"][q_idx]
        try:
            t_ca = t_res["CA"]
            q_ca = q_res["CA"]
            fixed_atoms.append(t_ca)
            moving_atoms.append(q_ca)
        except KeyError:
            continue

    if len(fixed_atoms) < MIN_ALIGNED:
        raise ValueError(f"Only {len(fixed_atoms)} Cα pairs after filtering")

    # Superimpose
    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    sup.apply(template_model.get_atoms())
    rmsd = float(sup.rms)

    # Find heme in template
    hemes = find_heme(template_model)
    if not hemes:
        raise ValueError("No heme in template after alignment")

    # Copy heme to target
    target_chain = target_model[best["t_cid"]]
    used_ids = {r.id[1] for r in target_chain}
    next_id = max(used_ids, default=0) + 100

    fe_coords = None
    for heme_res in hemes:
        new_res = copy.deepcopy(heme_res)
        new_res.detach_parent()
        while next_id in used_ids: next_id += 1
        het = f"H_{new_res.get_resname().strip()}"
        new_res.id = (het, next_id, " ")
        target_chain.add(new_res)
        used_ids.add(next_id)
        next_id += 1
        if not fe_coords:
            fe_coords = get_fe_coords(new_res)

    if not fe_coords:
        raise ValueError("Transplanted heme has no Fe atom")

    # Save
    output_path = HEME_DIR / f"{Path(target_pdb).stem}_heme.pdb"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(target_struct)
    io.save(str(output_path))

    # Quality tier
    identity = best["identity"]
    if rmsd < 2.0 and identity > 0.50: tier = "tier1"
    elif rmsd < 3.0 and identity > 0.30: tier = "tier2"
    elif rmsd < 5.0 and identity > 0.25: tier = "tier3"
    else: tier = "tier4"

    return {
        "rmsd": rmsd, "n_ca": len(fixed_atoms), "fe": fe_coords,
        "tier": tier, "identity": identity, "output": str(output_path),
        "template_chain": best["q_cid"],
    }

# ---------------------------------------------------------------------------
# Process one enzyme
# ---------------------------------------------------------------------------
def process_one(gid, uid, target_pdb, target_seq, templates, use_fallback):
    """Process one enzyme: try templates in order, fallback to AlphaFill API."""
    aligner = make_aligner()
    sess = requests.Session()
    sess.headers["User-Agent"] = USER_AGENT

    # Rank templates by sequence identity
    ranked = []
    for tmpl in templates:
        identity, pairs, n_aln = align_sequences(target_seq, tmpl["seq"], aligner)
        if identity >= MIN_IDENTITY and n_aln >= MIN_ALIGNED:
            ranked.append((identity, tmpl))
    ranked.sort(key=lambda x: -x[0])

    if not ranked:
        return {"gid": gid, "status": "no_template", "has_heme": "0",
                "notes": f"No template with identity >= {MIN_IDENTITY}"}

    # Try templates in order
    best_fail = None
    for identity, tmpl in ranked[:10]:  # try top 10
        try:
            tmpl_path = resolve_template_path(tmpl, sess)
            result = transplant_heme(target_pdb, tmpl_path, target_seq, tmpl["seq"], aligner)

            if result["tier"] != "tier4":
                return {
                    "gid": gid, "status": "heme_transplanted", "has_heme": "1",
                    "structure_source": "alphafold_heme_transplant_v2",
                    "sequence_identity": f"{result['identity']:.4f}",
                    "global_rmsd": f"{result['rmsd']:.2f}",
                    "quality_tier": result["tier"],
                    "aligned_ca_count": str(result["n_ca"]),
                    "template_pdb": tmpl["gid"],
                    "template_source": tmpl["source"],
                    "heme_transplant_pdb_path": result["output"],
                    "heme_fe_x": f"{result['fe'][0]:.3f}",
                    "heme_fe_y": f"{result['fe'][1]:.3f}",
                    "heme_fe_z": f"{result['fe'][2]:.3f}",
                    "notes": f"{result['tier']}, rmsd={result['rmsd']:.2f}, identity={result['identity']:.3f}, n_ca={result['n_ca']}",
                }
            # tier4 — try next template
            if best_fail is None or result["rmsd"] < float(best_fail.get("global_rmsd", "999")):
                best_fail = {
                    "gid": gid, "status": "tier4_best",
                    "global_rmsd": f"{result['rmsd']:.2f}",
                    "sequence_identity": f"{result['identity']:.4f}",
                    "template_pdb": tmpl["gid"],
                }
        except Exception as e:
            continue

    # All local templates failed → try AlphaFill API
    if use_fallback:
        try:
            with open(target_pdb, "rb") as f:
                resp = sess.post(ALPHAFILL_API, files={"structure": (Path(target_pdb).name, f)}, timeout=TIMEOUT)
            if resp.status_code == 200:
                data = resp.json()
                entry_id = data.get("id", "")
                if entry_id:
                    # Wait and download
                    time.sleep(5)
                    cif_resp = sess.get(f"{ALPHAFILL_API}/{entry_id}", timeout=TIMEOUT)
                    if cif_resp.status_code == 200:
                        cif_dir = STRUCTURES / "alphafill_api" / "cif"
                        cif_dir.mkdir(parents=True, exist_ok=True)
                        cif_path = cif_dir / f"{entry_id}.cif"
                        cif_path.write_bytes(cif_resp.content)
                        # Check for heme
                        struct = parse_structure(cif_path)
                        hemes = find_heme(struct[0])
                        if hemes:
                            fe = None
                            for h in hemes:
                                fe = get_fe_coords(h)
                                if fe: break
                            if fe:
                                return {
                                    "gid": gid, "status": "downloaded_with_heme",
                                    "has_heme": "1", "structure_source": "alphafill_api",
                                    "alphafill_cif_path": str(cif_path),
                                    "heme_fe_x": f"{fe[0]:.3f}", "heme_fe_y": f"{fe[1]:.3f}", "heme_fe_z": f"{fe[2]:.3f}",
                                    "notes": f"AlphaFill API fallback, entry={entry_id}",
                                }
        except Exception as e:
            pass

    # Return best failed attempt info
    notes = f"All {min(len(ranked),10)} templates gave tier4"
    if best_fail:
        notes += f"; best rmsd={best_fail.get('global_rmsd','?')}"
    return {"gid": gid, "status": "heme_transplant_failed", "has_heme": "0", "notes": notes}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    parser.add_argument("--no-fallback", action="store_true")
    parser.add_argument("--limit", type=int, default=0, help="Process only first N targets (for testing)")
    args = parser.parse_args()

    log("[setup] Loading data...")
    # Load manifest
    with MANIFEST_PATH.open("r", encoding="utf-8-sig", newline="") as f:
        manifest = list(csv.DictReader(f))
    # Load sequences
    with (COMBINED / "global_enzymes.csv").open("r", encoding="utf-8-sig") as f:
        enzyme_seqs = {r["global_enzyme_id"]: norm(r.get("canonical_sequence","")).upper() for r in csv.DictReader(f)}

    # Build template library
    templates = build_template_library(manifest, enzyme_seqs)
    log(f"[setup] {len(templates)} templates (AlphaFill CIF + PCPD PDB + S1 RCSB)")

    # Identify targets: all enzymes that need heme transplant
    targets = []
    for i, r in enumerate(manifest):
        if r.get("status") not in ("alphafill_not_found", "no_heme", "heme_transplanted", "tier4_best", "heme_transplant_failed"):
            continue
        uid = norm(r.get("canonical_uniprot_id",""))
        if not uid or not VALID_UNIPROT.match(uid): continue
        # Must have AlphaFold PDB on disk
        af_path = ALPHAFOLD_PDB_DIR / f"AF-{uid}-F1.pdb"
        if not af_path.exists(): continue
        seq = enzyme_seqs.get(r["global_enzyme_id"], "")
        if not seq: continue
        targets.append((i, r["global_enzyme_id"], uid, str(af_path), seq))

    if args.limit > 0:
        targets = targets[:args.limit]

    log(f"[setup] {len(targets)} targets to process, {args.workers} workers")

    if not targets:
        log("[done] Nothing to process")
        return

    # Process
    ok = fail = 0
    results = {}  # manifest_idx → result dict

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for idx, gid, uid, af_path, seq in targets:
            fut = executor.submit(process_one, gid, uid, af_path, seq, templates, not args.no_fallback)
            futures[fut] = idx

        for fut in as_completed(futures):
            idx = futures[fut]
            try:
                result = fut.result(timeout=300)
            except Exception as e:
                result = {"gid": manifest[idx]["global_enzyme_id"], "status": "error", "has_heme": "0", "notes": str(e)[:200]}

            results[idx] = result
            if result.get("has_heme") == "1": ok += 1
            else: fail += 1

            if (ok + fail) % 20 == 0:
                log(f"[progress] {ok+fail}/{len(targets)}: {ok} ok, {fail} fail")

    # Update manifest
    all_fields = set()
    for r in manifest: all_fields.update(r.keys())
    for _, result in results.items():
        for k in result:
            if k != "gid": all_fields.add(k)

    for idx, result in results.items():
        for k, v in result.items():
            if k != "gid":
                manifest[idx][k] = str(v)
        manifest[idx]["transplant_version"] = "seqalign_v2"

    # Save
    fields = list(manifest[0].keys())
    for f in sorted(all_fields - set(fields)):
        fields.append(f)

    with MANIFEST_PATH.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in manifest:
            w.writerow({k: r.get(k, "") for k in fields})

    log(f"\n{'='*60}")
    log(f"[DONE] {ok+fail}/{len(targets)} processed: {ok} ok, {fail} fail")
    log(f"[DONE] Manifest updated: {MANIFEST_PATH}")
    log(f"{'='*60}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted"); raise SystemExit(130)
    except Exception as e:
        log(f"ERROR: {e}"); import traceback; traceback.print_exc(); raise SystemExit(1)
