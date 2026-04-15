"""
EXP005 full run: build enzyme_resid_map.pt for all 1622 enzymes.

Algorithm per enzyme (proven on 5 smoke cases):
  1. Read one complex PDB (any dock belonging to this enzyme)
  2. Multi-chain global alignment to UniProt sequence
  3. For experimental_pdb: cross-check with SIFTS
  4. Classify tier: gold / trusted / partial / suspect / skip
  5. Store per-residue (chain, resid, icode) -> uniprot_pos mapping

Runs on server with 20-worker process pool (leaves 5 cores for SIFTS HTTP).
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv, json, os, sys, time, traceback, urllib.request, urllib.error
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import torch

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
ENZ_CSV = BASE / "data/Enzymes.csv"
COMPLEX_DIR = BASE / "data/structure/complex"
SPLIT_CSVS = [
    BASE / "data/splits/random/training_datas_0_pt.csv",
    BASE / "data/splits/random/val_datas_0_pt.csv",
    BASE / "data/splits/random/testing_datas_0_pt.csv",
]
# Use dock_sidecar.pt from pt_cache_allfix_unified (authoritative: only docks
# that actually end up in training)
PT_CACHE = BASE / "data/pt_cache_allfix_unified/random"
SIDECAR_PATHS = [
    Path("/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_dualgraph_allfix_unified/random/train/dock_sidecar.pt"),
    Path("/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_dualgraph_allfix_unified/random/val/dock_sidecar.pt"),
    Path("/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_dualgraph_allfix_unified/random/test/dock_sidecar.pt"),
]
META_PATH = Path("/tmp/enzyme_meta.json")
OVERLAY_DIR = BASE / "data/pt_cache_dualgraph_allfix_unified"
SIFTS_CACHE = OVERLAY_DIR / "sifts_cache"
OUT_PT = OVERLAY_DIR / "enzyme_resid_map.pt"
OUT_REPORT = OVERLAY_DIR / "enzyme_alignment_report.md"

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M",
}

N_WORKERS = 20  # leave ~5 cores for SIFTS HTTP phase later

# -------------------------------------------------------------------
# Parsing
# -------------------------------------------------------------------
def parse_complex_receptor_residues(pdb_path):
    """Return list of (chain, resid, icode, resname, aa1). One entry per unique
    (chain, resid, icode). Stops at COMPND (ligand separator)."""
    seen = {}
    order = []
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "COMPND":
                break
            if rec != "ATOM":
                continue
            resname = line[17:20].strip()
            if resname not in THREE_TO_ONE:
                continue
            chain = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26]
            key = (chain, resseq, icode)
            if key not in seen:
                seen[key] = resname
                order.append(key)
    result = []
    for key in order:
        chain, resseq, icode = key
        resname = seen[key]
        result.append((chain, resseq, icode, resname, THREE_TO_ONE[resname]))
    return result


# -------------------------------------------------------------------
# Alignment (per-chain, multi-chain aware)
# -------------------------------------------------------------------
def build_aligner():
    from Bio import Align
    a = Align.PairwiseAligner(mode="global")
    a.match_score = 2.0
    a.mismatch_score = 0.0
    a.open_gap_score = -10.0
    a.extend_gap_score = -0.5
    return a


def align_residues(residues, uniprot_seq, aligner,
                   min_chain_identity=0.70, min_chain_coverage=0.50):
    """Multi-chain alignment. Merges all chains that pass the per-chain gate."""
    by_chain = defaultdict(list)
    for entry in residues:
        by_chain[entry[0]].append(entry)

    combined = {}
    per_chain = {}
    used = []
    rejected = []
    n_struct_total = sum(len(v) for v in by_chain.values())
    n_mapped_total = 0
    n_match_total = 0

    for chain_id, chain_res in by_chain.items():
        n_struct = len(chain_res)
        if n_struct < 10:
            rejected.append([chain_id, "too_few_residues"])
            continue
        struct_seq = "".join(e[4] for e in chain_res)
        try:
            aln = aligner.align(struct_seq, uniprot_seq)[0]
        except Exception as exc:
            rejected.append([chain_id, f"align_err:{type(exc).__name__}"])
            continue

        s_blocks, u_blocks = aln.aligned
        mapped_this = {}
        n_match = 0
        for (s_beg, s_end), (u_beg, u_end) in zip(s_blocks, u_blocks):
            for off in range(min(int(s_end - s_beg), int(u_end - u_beg))):
                s_idx = int(s_beg) + off
                u_idx = int(u_beg) + off
                chain, resid, icode, _, aa1 = chain_res[s_idx]
                mapped_this[(chain, resid, icode)] = u_idx
                if uniprot_seq[u_idx] == aa1:
                    n_match += 1

        n_mapped = len(mapped_this)
        identity = (n_match / n_mapped) if n_mapped else 0.0
        coverage = (n_mapped / n_struct) if n_struct else 0.0
        per_chain[chain_id] = {
            "n_struct": n_struct, "n_mapped": n_mapped, "n_match": n_match,
            "identity": round(identity, 4), "coverage": round(coverage, 4),
        }
        if identity >= min_chain_identity and coverage >= min_chain_coverage:
            used.append(chain_id)
            combined.update(mapped_this)
            n_mapped_total += n_mapped
            n_match_total += n_match
        else:
            rejected.append([chain_id, f"low_id_cov(id={identity:.3f},cov={coverage:.3f})"])

    if not used:
        return None

    return {
        "chains_used": used,
        "per_chain": per_chain,
        "rejected_chains": rejected,
        "n_struct_total": n_struct_total,
        "n_mapped_total": n_mapped_total,
        "n_match_total": n_match_total,
        "identity": (n_match_total / n_mapped_total) if n_mapped_total else 0.0,
        "coverage": (n_mapped_total / n_struct_total) if n_struct_total else 0.0,
        "resid_map": combined,
    }


# -------------------------------------------------------------------
# Worker entry
# -------------------------------------------------------------------
_aligner = None
_enz_rows = None

def _worker_init():
    global _aligner, _enz_rows
    _aligner = build_aligner()
    with open(ENZ_CSV, encoding="utf-8-sig") as f:
        _enz_rows = list(csv.DictReader(f))


def _worker_process(job):
    enz_id, dock_index = job
    try:
        row = _enz_rows[enz_id]
        seq = row["Protein sequence"].strip()
        cpath = COMPLEX_DIR / f"{dock_index}.pdb"
        if not cpath.exists():
            return {"enz_id": enz_id, "status": "no_complex_pdb", "dock": dock_index}
        residues = parse_complex_receptor_residues(cpath)
        if not residues:
            return {"enz_id": enz_id, "status": "no_residues", "dock": dock_index}
        ali = align_residues(residues, seq, _aligner)
        if ali is None:
            return {"enz_id": enz_id, "status": "align_failed", "dock": dock_index}
        return {
            "enz_id": enz_id,
            "status": "ok",
            "dock": dock_index,
            "seq_len": len(seq),
            "ali": ali,
        }
    except Exception as exc:
        return {
            "enz_id": enz_id,
            "status": "exception",
            "dock": dock_index,
            "error": f"{type(exc).__name__}: {exc}",
            "trace": traceback.format_exc()[:500],
        }


# -------------------------------------------------------------------
# SIFTS cross-check
# -------------------------------------------------------------------
SIFTS_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb_id}"

def fetch_sifts(pdb_id, retries=3):
    pdb_id = pdb_id.lower()
    SIFTS_CACHE.mkdir(parents=True, exist_ok=True)
    cache_path = SIFTS_CACHE / f"{pdb_id}.json"
    if cache_path.exists():
        try:
            return json.loads(cache_path.read_text())
        except Exception:
            pass
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                SIFTS_URL.format(pdb_id=pdb_id),
                headers={"User-Agent": "EZSpec-EXP005/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode())
            cache_path.write_text(json.dumps(data))
            return data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                cache_path.write_text(json.dumps(None))  # negative cache
                return None
            last_err = e
        except Exception as e:
            last_err = e
        time.sleep(1.5 * (attempt + 1))
    print(f"  [SIFTS] FAIL {pdb_id}: {last_err}", flush=True)
    return None


def sifts_build_map(residues, sifts_data, target_uniprot, pdb_id, seq):
    """Returns dict (chain, resid, icode) -> uniprot_pos_0idx + stats."""
    pdb_id = pdb_id.lower()
    if not sifts_data:
        return None, "no_data"
    entry = sifts_data.get(pdb_id, {})
    upmap = entry.get("UniProt", {})
    if target_uniprot not in upmap:
        return None, f"uniprot_{target_uniprot}_not_in_sifts"
    mappings = upmap[target_uniprot].get("mappings", [])

    by_chain = defaultdict(list)
    for c, rs, ic, rn, a1 in residues:
        by_chain[c].append((rs, ic, rn, a1))

    result = {}
    for chain_id, ch_res in by_chain.items():
        segs = [m for m in mappings if m.get("chain_id") == chain_id]
        if not segs:
            continue
        for m in segs:
            unp_start = int(m["unp_start"])
            author_start_d = m.get("start", {})
            author_start = author_start_d.get("author_residue_number")
            if author_start is None:
                continue
            unp_end = int(m["unp_end"])
            seg_len = unp_end - unp_start + 1
            delta = unp_start - author_start  # 1-indexed
            seg_author_min = author_start
            seg_author_max = author_start + seg_len - 1
            for resid, icode, rn, a1 in ch_res:
                if seg_author_min <= resid <= seg_author_max:
                    u_pos = resid + delta - 1  # 0-indexed
                    if 0 <= u_pos < len(seq):
                        result[(chain_id, resid, icode)] = u_pos
    return (result if result else None), "ok" if result else "no_match_in_chains"


# -------------------------------------------------------------------
# Tier classification
# -------------------------------------------------------------------
def classify_tier(ali, sifts_agree=None):
    """
    gold:    SIFTS 100% agreement OR (identity >= 0.99 AND coverage >= 0.95)
    trusted: identity >= 0.90 AND coverage >= 0.90
    partial: identity >= 0.75 AND coverage >= 0.75
    suspect: otherwise
    """
    if sifts_agree is not None:
        agree, total = sifts_agree
        if total > 0 and agree / total >= 0.99:
            return "gold"
    iden = ali["identity"]
    cov = ali["coverage"]
    if iden >= 0.99 and cov >= 0.95:
        return "gold"
    if iden >= 0.90 and cov >= 0.90:
        return "trusted"
    if iden >= 0.75 and cov >= 0.75:
        return "partial"
    return "suspect"


# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------
def main():
    print("=" * 70)
    print("EXP005: Full enzyme_resid_map build")
    print("=" * 70)
    t_start = time.time()

    # Load enzyme meta
    meta = json.loads(META_PATH.read_text())
    meta_int = {int(k): v for k, v in meta.items()}
    print(f"loaded enzyme meta for {len(meta_int)} enzymes", flush=True)

    # Build ENZYME -> list of docks from dock_sidecar.pt (authoritative pt_cache source)
    enzyme_to_docks = defaultdict(list)
    for sp in SIDECAR_PATHS:
        if not sp.exists():
            print(f"[WARN] sidecar missing: {sp}", flush=True)
            continue
        sc = torch.load(str(sp), weights_only=False)
        dock_idx = sc["dock_indices"].tolist() if hasattr(sc["dock_indices"], "tolist") else list(sc["dock_indices"])
        # Also need enzyme_ids - loaded via sample_ids + CSV
        # Actually easier: read base index.pt which has enzyme_ids aligned by k
        base_idx_path = Path(str(sp).replace("dock_sidecar.pt", "index.pt"))
        base_idx = torch.load(str(base_idx_path), weights_only=False)
        enz_ids = base_idx["enzyme_ids"].tolist() if hasattr(base_idx["enzyme_ids"], "tolist") else list(base_idx["enzyme_ids"])
        assert len(dock_idx) == len(enz_ids), f"len mismatch in {sp}"
        for e, d in zip(enz_ids, dock_idx):
            enzyme_to_docks[int(e)].append(int(d))
    print(f"built enzyme_to_docks from pt_cache sidecars: {len(enzyme_to_docks)} enzymes", flush=True)

    # Pick FIRST EXISTING complex PDB for each enzyme
    enzyme_to_dock = {}
    missing_all = []
    for e, docks in enzyme_to_docks.items():
        for d in docks:
            if (COMPLEX_DIR / f"{d}.pdb").exists():
                enzyme_to_dock[e] = d
                break
        else:
            missing_all.append(e)
    print(f"enzymes with valid dock PDB: {len(enzyme_to_dock)}", flush=True)
    print(f"enzymes with NO complex PDB at all: {len(missing_all)}", flush=True)

    # Build job list
    jobs = []
    no_dock = []
    for enz_id in range(len(meta_int)):
        if enz_id in enzyme_to_dock:
            jobs.append((enz_id, enzyme_to_dock[enz_id]))
        else:
            no_dock.append(enz_id)
    print(f"jobs: {len(jobs)}, no_dock/no_pdb: {len(no_dock)}", flush=True)

    # Load Enzymes.csv here for SIFTS verify
    with open(ENZ_CSV, encoding="utf-8-sig") as f:
        enz_rows = list(csv.DictReader(f))

    # ---- Phase 1: Alignment (parallel) ----
    print(f"\n[Phase 1] Alignment with {N_WORKERS} workers...", flush=True)
    results = {}
    t_ali = time.time()
    with ProcessPoolExecutor(max_workers=N_WORKERS, initializer=_worker_init) as ex:
        futs = [ex.submit(_worker_process, j) for j in jobs]
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            results[r["enz_id"]] = r
            done += 1
            if done % 200 == 0:
                print(f"  alignment {done}/{len(jobs)}  ({time.time()-t_ali:.1f}s)", flush=True)
    print(f"  alignment done in {time.time()-t_ali:.1f}s", flush=True)

    # ---- Phase 2: SIFTS cross-check for experimental_pdb ----
    print(f"\n[Phase 2] SIFTS cross-check for experimental_pdb...", flush=True)
    t_sifts = time.time()
    sifts_stats = {"queried": 0, "hit": 0, "miss": 0, "agree_100": 0, "disagree": 0}
    for enz_id, r in results.items():
        if r["status"] != "ok":
            continue
        m = meta_int[enz_id]
        if m["source"] != "experimental_pdb" or not m["pdb_id"]:
            continue
        sifts_stats["queried"] += 1
        try:
            sifts_data = fetch_sifts(m["pdb_id"])
        except Exception as exc:
            print(f"  SIFTS fetch error for enz={enz_id} pdb={m['pdb_id']}: {exc}", flush=True)
            continue
        if sifts_data is None:
            sifts_stats["miss"] += 1
            r["sifts_status"] = "no_data"
            continue
        sifts_stats["hit"] += 1
        seq = enz_rows[enz_id]["Protein sequence"].strip()
        cpath = COMPLEX_DIR / f"{r['dock']}.pdb"
        residues = parse_complex_receptor_residues(cpath)
        sm_map, sm_notes = sifts_build_map(residues, sifts_data, m["uniprot"], m["pdb_id"], seq)
        if not sm_map:
            r["sifts_status"] = f"no_mapping:{sm_notes}"
            continue
        # Compare against alignment map
        ali_map = r["ali"]["resid_map"]
        common = set(sm_map) & set(ali_map)
        if not common:
            r["sifts_status"] = "no_overlap"
            continue
        agree = sum(1 for k in common if sm_map[k] == ali_map[k])
        r["sifts_agree"] = [agree, len(common)]
        r["sifts_status"] = "ok"
        if agree == len(common):
            sifts_stats["agree_100"] += 1
        else:
            sifts_stats["disagree"] += 1
            print(f"  SIFTS DISAGREE enz={enz_id} pdb={m['pdb_id']}: {agree}/{len(common)}", flush=True)
    print(f"  SIFTS done in {time.time()-t_sifts:.1f}s: {sifts_stats}", flush=True)

    # ---- Phase 3: Tier classification + output ----
    print(f"\n[Phase 3] Classifying tiers & saving output...", flush=True)
    final = {}
    tier_count = defaultdict(int)
    source_tier = defaultdict(lambda: defaultdict(int))
    suspect_list = []

    for enz_id, r in results.items():
        m = meta_int.get(enz_id, {"source": "unknown", "uniprot": "", "pdb_id": None})
        if r["status"] != "ok":
            tier_count[r["status"]] += 1
            source_tier[m["source"]][r["status"]] += 1
            continue
        ali = r["ali"]
        sifts_agree = r.get("sifts_agree")
        tier = classify_tier(ali, sifts_agree)
        tier_count[tier] += 1
        source_tier[m["source"]][tier] += 1
        if tier == "suspect":
            suspect_list.append({
                "enz_id": enz_id,
                "source": m["source"],
                "uniprot": m["uniprot"],
                "identity": round(ali["identity"], 4),
                "coverage": round(ali["coverage"], 4),
                "n_mapped": ali["n_mapped_total"],
                "n_struct": ali["n_struct_total"],
            })

        final[enz_id] = {
            "uniprot": m["uniprot"],
            "source": m["source"],
            "pdb_id": m["pdb_id"],
            "seq_len": r["seq_len"],
            "dock": r["dock"],
            "chains_used": ali["chains_used"],
            "per_chain": ali["per_chain"],
            "rejected_chains": ali["rejected_chains"],
            "n_struct_total": ali["n_struct_total"],
            "n_mapped_total": ali["n_mapped_total"],
            "n_match_total": ali["n_match_total"],
            "identity": ali["identity"],
            "coverage": ali["coverage"],
            "tier": tier,
            "resid_map": {f"{c}|{r_}|{ic}": u for (c, r_, ic), u in ali["resid_map"].items()},
            "sifts_agree": sifts_agree,
            "sifts_status": r.get("sifts_status", "not_queried"),
        }

    OUT_PT.parent.mkdir(parents=True, exist_ok=True)
    torch.save(final, str(OUT_PT))
    print(f"  saved: {OUT_PT}  ({OUT_PT.stat().st_size} bytes)", flush=True)

    # ---- Phase 4: Report ----
    print(f"\n[Phase 4] Writing report...", flush=True)
    lines = []
    lines.append("# EXP005 enzyme_resid_map alignment report")
    lines.append("")
    lines.append(f"- generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"- total enzymes: {len(meta_int)}")
    lines.append(f"- built: {len(final)}")
    lines.append(f"- no_dock (not in any split): {len(no_dock)}")
    lines.append(f"- total time: {time.time()-t_start:.1f}s")
    lines.append("")
    lines.append("## Tier distribution")
    lines.append("")
    lines.append("| tier | count |")
    lines.append("|---|---|")
    for t in ["gold", "trusted", "partial", "suspect"]:
        lines.append(f"| {t} | {tier_count.get(t, 0)} |")
    lines.append(f"| (align_failed / no_complex_pdb / no_residues / exception) | "
                 f"{sum(v for k, v in tier_count.items() if k not in ('gold','trusted','partial','suspect'))} |")
    lines.append("")
    lines.append("## By structure source")
    lines.append("")
    lines.append("| source | gold | trusted | partial | suspect | failed |")
    lines.append("|---|---|---|---|---|---|")
    for src in sorted(source_tier.keys()):
        row = source_tier[src]
        fail = sum(v for k, v in row.items() if k not in ("gold", "trusted", "partial", "suspect"))
        lines.append(f"| {src} | {row.get('gold',0)} | {row.get('trusted',0)} | "
                     f"{row.get('partial',0)} | {row.get('suspect',0)} | {fail} |")
    lines.append("")
    lines.append("## SIFTS cross-check (experimental_pdb only)")
    lines.append("")
    for k, v in sifts_stats.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    if suspect_list:
        lines.append("## Suspect enzymes (need manual review)")
        lines.append("")
        lines.append("| enz_id | source | uniprot | identity | coverage | mapped/struct |")
        lines.append("|---|---|---|---|---|---|")
        for s in suspect_list[:50]:
            lines.append(f"| {s['enz_id']} | {s['source']} | {s['uniprot']} | "
                         f"{s['identity']:.4f} | {s['coverage']:.4f} | "
                         f"{s['n_mapped']}/{s['n_struct']} |")
        if len(suspect_list) > 50:
            lines.append(f"")
            lines.append(f"(+ {len(suspect_list)-50} more)")

    report_text = "\n".join(lines)
    OUT_REPORT.write_text(report_text)
    print(f"  saved: {OUT_REPORT}", flush=True)

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"total built:    {len(final)}")
    print(f"no_dock:        {len(no_dock)}")
    print(f"total time:     {time.time()-t_start:.1f}s")
    print(f"tier breakdown: {dict(tier_count)}")
    print(f"SIFTS stats:    {sifts_stats}")
    print(f"suspect count:  {len(suspect_list)}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
'''


def main():
    import subprocess, sys
    print("[local driver] running full resid_map build on server...")
    result = subprocess.run(
        ["ssh", "autodl-4x5090-bj",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=3600,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write("\n[STDERR]\n")
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
