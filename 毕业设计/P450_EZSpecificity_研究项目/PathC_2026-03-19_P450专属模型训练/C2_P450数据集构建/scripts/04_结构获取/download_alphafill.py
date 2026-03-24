#!/usr/bin/env python3
"""
Phase 6 Step 2: Batch download AlphaFill structures for P450 enzymes.

For each enzyme with UniProt ID:
  1. Download JSON metadata from AlphaFill API → check if HEM was transplanted
  2. If HEM found → download the CIF structure file (contains protein + heme)
  3. Record everything in receptor_manifest.csv

Also handles:
  - 349 existing PDB enzymes (S1 RCSB + PCPD) → mark in manifest
  - 108 no-UniProt enzymes → mark as needs_colabfold

Resume support: skips enzymes already processed in manifest.
"""
from __future__ import annotations

import csv, json, os, sys, time
from collections import defaultdict
from pathlib import Path
from typing import Any

import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
SCRIPT = Path(__file__).resolve()
PROJECT = SCRIPT.parents[2]
COMBINED = PROJECT / "data" / "combined"
STRUCTURES = PROJECT / "data" / "structures"

ALPHAFILL_BASE = "https://alphafill.eu/v1/aff"
REQUEST_INTERVAL = 1.5  # seconds between requests
TIMEOUT = 60

MANIFEST_PATH = STRUCTURES / "manifests" / "receptor_manifest.csv"
AF_JSON_DIR = STRUCTURES / "alphafill" / "json"
AF_CIF_DIR = STRUCTURES / "alphafill" / "cif"

MANIFEST_FIELDS = [
    "global_enzyme_id", "canonical_uniprot_id", "source_category",
    "structure_source", "status", "has_heme",
    "alphafill_id", "template_pdb", "sequence_identity", "global_rmsd",
    "local_rmsd", "clash_score",
    "alphafill_json_path", "alphafill_cif_path",
    "existing_pdb_path", "needs_colabfold", "notes",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def log(msg): print(msg, file=sys.stderr, flush=True)
def norm(v):
    t = str(v).strip() if v is not None else ""
    return "" if t.lower() in {"","na","n/a","none","null"} else t

def read_csv(path):
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as f:
        return list(csv.DictReader(f))

def write_manifest(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        w.writeheader()
        for r in rows: w.writerow(r)

# ---------------------------------------------------------------------------
# Classify enzymes: existing PDB vs AlphaFill candidate vs no source
# ---------------------------------------------------------------------------
def classify_enzymes():
    """Returns dict: global_enzyme_id → {category, uid, existing_path, ...}"""
    enzymes = read_csv(COMBINED / "global_enzymes.csv")
    xref = read_csv(COMBINED / "enzyme_xref.csv")

    # S1 RCSB enzymes with pdb_id
    s1_src = PROJECT / "data" / "sources" / "Source_RCSB" / "enzymes.csv"
    s1_enz = {r["enzyme_id"]: r.get("pdb_id","").strip() for r in read_csv(s1_src)}

    # PCPD source enzymes: enzyme_id → p450_symbol
    pcpd_src = PROJECT / "data" / "sources" / "Source_PCPD" / "enzymes.csv"
    pcpd_eid_sym = {r["enzyme_id"]: r.get("p450_symbol","").strip() for r in read_csv(pcpd_src)}

    # PCPD PDB files on disk
    pcpd_pdb_dir = PROJECT / "downloads" / "PCPD" / "PDB"
    pcpd_on_disk = set()
    if pcpd_pdb_dir.exists():
        pcpd_on_disk = {f.replace(".pdb","").upper() for f in os.listdir(pcpd_pdb_dir) if f.endswith(".pdb")}

    # Build gid → source info from xref
    gid_src = defaultdict(lambda: {"s1_eid": "", "pcpd_eid": ""})
    for x in xref:
        gid = x["global_enzyme_id"]
        if x["source"] == "S1_RCSB": gid_src[gid]["s1_eid"] = x["source_enzyme_id"]
        if x["source"] == "S9_PCPD": gid_src[gid]["pcpd_eid"] = x["source_enzyme_id"]

    result = {}
    for e in enzymes:
        gid = e["global_enzyme_id"]
        uid = norm(e.get("canonical_uniprot_id",""))
        info = gid_src[gid]

        # Check S1 PDB
        s1_pdb = s1_enz.get(info["s1_eid"], "") if info["s1_eid"] else ""

        # Check PCPD PDB
        pcpd_sym = pcpd_eid_sym.get(info["pcpd_eid"], "") if info["pcpd_eid"] else ""
        pcpd_pdb = pcpd_sym if pcpd_sym and pcpd_sym.upper() in pcpd_on_disk else ""

        if s1_pdb:
            result[gid] = {"category": "experimental_pdb", "uid": uid,
                           "existing_path": f"s1_rcsb/{s1_pdb}", "pdb_id": s1_pdb}
        elif pcpd_pdb:
            result[gid] = {"category": "pcpd_predicted", "uid": uid,
                           "existing_path": f"pcpd/{pcpd_pdb}.pdb", "pdb_id": pcpd_pdb}
        elif uid:
            result[gid] = {"category": "alphafill_candidate", "uid": uid,
                           "existing_path": "", "pdb_id": ""}
        else:
            result[gid] = {"category": "no_source", "uid": "",
                           "existing_path": "", "pdb_id": ""}
    return result

# ---------------------------------------------------------------------------
# AlphaFill download
# ---------------------------------------------------------------------------
def find_heme_in_json(data):
    """Search AlphaFill JSON for HEM transplant. Returns best HEM hit or None."""
    best = None
    hits = data.get("hits") or []
    for hit in hits:
        transplants = hit.get("transplants") or []
        for t in transplants:
            if not isinstance(t, dict): continue
            aid = norm(t.get("analogue_id",""))
            cid = norm(t.get("compound_id",""))
            if aid == "HEM" or cid == "HEM":
                alignment = hit.get("alignment") or {}
                candidate = {
                    "template_pdb": norm(hit.get("pdb_id","")),
                    "identity": alignment.get("identity", ""),
                    "global_rmsd": hit.get("global_rmsd", ""),
                    "local_rmsd": t.get("local_rmsd", ""),
                    "clash_score": ((t.get("clash") or {}).get("score", "")),
                }
                cur_id = candidate["identity"] or 0
                best_id = best["identity"] if best else 0
                if best is None or cur_id > best_id:
                    best = candidate
    return best

def download_alphafill(sess, uid):
    """Download AlphaFill JSON + CIF for a UniProt ID. Returns manifest row fields."""
    afill_id = f"{uid}-F1"
    json_path = AF_JSON_DIR / f"{afill_id}.json"
    cif_path = AF_CIF_DIR / f"{afill_id}.cif"
    result = {"alphafill_id": afill_id, "alphafill_json_path": "", "alphafill_cif_path": ""}

    # Step 1: Download JSON metadata
    json_url = f"{ALPHAFILL_BASE}/{afill_id}/json"
    try:
        resp = sess.get(json_url, timeout=TIMEOUT)
    except requests.RequestException as e:
        return {**result, "status": "error", "has_heme": "", "notes": str(e)}

    if resp.status_code == 404:
        return {**result, "status": "alphafill_not_found", "has_heme": "0",
                "notes": "Entry not in AlphaFill DB"}
    if resp.status_code == 429:
        return {**result, "status": "rate_limited", "has_heme": "", "notes": "Try again later"}
    if resp.status_code != 200:
        return {**result, "status": f"http_{resp.status_code}", "has_heme": "",
                "notes": resp.text[:200]}

    # Save JSON
    data = resp.json()
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    result["alphafill_json_path"] = str(json_path.relative_to(PROJECT))

    # Step 2: Check for HEM
    heme = find_heme_in_json(data)
    if not heme:
        return {**result, "status": "no_heme", "has_heme": "0",
                "template_pdb": "", "sequence_identity": "", "global_rmsd": "",
                "local_rmsd": "", "clash_score": "",
                "notes": f"AlphaFill has {len(data.get('hits') or [])} hits but no HEM"}

    # Step 3: Download CIF structure
    cif_url = f"{ALPHAFILL_BASE}/{afill_id}"
    try:
        cif_resp = sess.get(cif_url, timeout=TIMEOUT)
    except requests.RequestException as e:
        return {**result, "status": "cif_error", "has_heme": "1",
                "template_pdb": heme["template_pdb"],
                "notes": f"JSON ok but CIF download failed: {e}"}

    if cif_resp.status_code != 200:
        return {**result, "status": f"cif_http_{cif_resp.status_code}", "has_heme": "1",
                "template_pdb": heme["template_pdb"],
                "notes": "JSON ok but CIF download failed"}

    cif_path.parent.mkdir(parents=True, exist_ok=True)
    cif_path.write_bytes(cif_resp.content)
    result["alphafill_cif_path"] = str(cif_path.relative_to(PROJECT))

    return {
        **result,
        "status": "downloaded_with_heme",
        "has_heme": "1",
        "template_pdb": heme["template_pdb"],
        "sequence_identity": heme["identity"],
        "global_rmsd": heme["global_rmsd"],
        "local_rmsd": heme["local_rmsd"],
        "clash_score": heme["clash_score"],
        "notes": "",
    }

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log("[setup] classifying 1,622 enzymes...")
    classifications = classify_enzymes()

    # Load existing manifest for resume
    existing = {}
    if MANIFEST_PATH.exists():
        for r in read_csv(MANIFEST_PATH):
            gid = r.get("global_enzyme_id","")
            st = r.get("status","")
            if st in ("downloaded_with_heme","no_heme","alphafill_not_found",
                       "experimental_pdb","pcpd_predicted","needs_colabfold"):
                existing[gid] = r

    sess = requests.Session()
    sess.headers["User-Agent"] = "P450-AlphaFill-Batch/1.0"

    results = []
    cats = defaultdict(int)
    total = len(classifications)
    af_to_download = sum(1 for v in classifications.values() if v["category"]=="alphafill_candidate" and v["uid"])

    log(f"[setup] total={total}, existing_pdb=349, alphafill_candidates={af_to_download}, no_source=108")
    log(f"[setup] already processed (resume): {len(existing)}")
    log(f"[setup] AlphaFill downloads needed: {af_to_download - sum(1 for gid,v in classifications.items() if v['category']=='alphafill_candidate' and gid in existing)}")

    idx = 0
    for gid in sorted(classifications):
        idx += 1
        info = classifications[gid]
        cat = info["category"]

        # Resume: skip already processed
        if gid in existing:
            results.append(existing[gid])
            cats[existing[gid].get("status","")] += 1
            continue

        row = {f: "" for f in MANIFEST_FIELDS}
        row["global_enzyme_id"] = gid
        row["canonical_uniprot_id"] = info["uid"]
        row["source_category"] = cat

        if cat in ("experimental_pdb", "pcpd_predicted"):
            row["structure_source"] = cat
            row["status"] = cat
            row["has_heme"] = "1"
            row["existing_pdb_path"] = info["existing_path"]
            row["needs_colabfold"] = "0"
            cats[cat] += 1

        elif cat == "alphafill_candidate":
            log(f"[{idx}/{total}] {gid} {info['uid']} → AlphaFill...")
            af_result = download_alphafill(sess, info["uid"])
            row["structure_source"] = "alphafill"
            row.update(af_result)
            row["needs_colabfold"] = "0"
            cats[af_result["status"]] += 1
            log(f"  status={af_result['status']} heme={af_result.get('has_heme','')}")
            time.sleep(REQUEST_INTERVAL)

        elif cat == "no_source":
            row["structure_source"] = ""
            row["status"] = "needs_colabfold"
            row["has_heme"] = "0"
            row["needs_colabfold"] = "1"
            row["notes"] = "No UniProt, no existing PDB"
            cats["needs_colabfold"] += 1

        results.append(row)

        # Periodic save (every 50 enzymes)
        if idx % 50 == 0:
            write_manifest(results, MANIFEST_PATH)
            log(f"  [checkpoint] saved manifest ({idx}/{total})")

    # Final save
    write_manifest(results, MANIFEST_PATH)

    # Summary
    log(f"\n{'='*60}")
    log(f"[DONE] receptor_manifest.csv: {len(results)} enzymes")
    for status, count in sorted(cats.items()):
        log(f"  {status}: {count}")
    log(f"{'='*60}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted"); raise SystemExit(130)
    except Exception as e:
        log(f"ERROR: {e}"); import traceback; traceback.print_exc(); raise SystemExit(1)
