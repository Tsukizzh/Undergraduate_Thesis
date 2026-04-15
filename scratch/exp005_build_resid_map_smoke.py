"""
EXP005 Step 3 (fix) — smoke test on 5 known cases.

Goal: for each of 5 enzymes (3/93, 8444/659, 24/838, 985/54, 30789/952),
build `resid -> uniprot_pos` via:
  1. pairwise global alignment of structure residues to Enzymes.csv UniProt seq
  2. for experimental_pdb: cross-check with SIFTS REST API

Run entirely on server (has Enzymes.csv + complex PDBs + BioPython + network).
Local driver streams REMOTE_SCRIPT via ssh stdin.

If 5/5 pass with tier=gold/trusted and 100% aa_match on pocket residues,
proceed to full 1622-enzyme version.
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv, json, gzip, sys, urllib.request, urllib.error, time
from pathlib import Path
from collections import defaultdict
from Bio import Align

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
ENZ_CSV = BASE / "data/Enzymes.csv"
COMPLEX_DIR = BASE / "data/structure/complex"
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"
SPLIT_CSVS = [
    BASE / "data/splits/random/training_datas_0_pt.csv",
    BASE / "data/splits/random/val_datas_0_pt.csv",
    BASE / "data/splits/random/testing_datas_0_pt.csv",
]

# Test cases: dock_index -> (enzyme_id, label, source_hint, pdb_id_if_experimental)
TEST_CASES = {
    3:     (93,  "PASS", "alphafill",        None),
    8444:  (659, "PASS", "pcpd_predicted",   None),
    24:    (838, "FAIL", "experimental_pdb", "1JPZ"),
    985:   (54,  "FAIL", "pcpd_predicted",   None),  # pcpd with weird numbering
    30789: (952, "FAIL", "alphafill",        None),  # alphafill with delta=+68
}

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M",  # selenomethionine -> methionine (cleanup normalizes)
}

# ---------- Load CSV-level data ----------
with open(ENZ_CSV, encoding="utf-8-sig") as f:
    ENZ_ROWS = list(csv.DictReader(f))

DOCK_TO_ENZYME = {}
for p in SPLIT_CSVS:
    with open(p, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            DOCK_TO_ENZYME[int(row["Dock Index"])] = int(row["Enzyme Index"])


# ---------- Structure parsing ----------
def parse_complex_receptor_residues(pdb_path):
    """
    Return list of (chain, resid, icode, aa_3letter, aa_1letter) for the
    receptor portion of a complex PDB. Only ATOM records, only standard/MSE
    residues, one entry per unique (chain, resid, icode, resname).
    Stops at COMPND boundary (ligand separator).
    """
    seen = {}
    order = []
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "COMPND":
                break  # ligand section starts
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


# ---------- Alignment ----------
def make_aligner():
    a = Align.PairwiseAligner(mode="global")
    a.match_score = 2.0
    a.mismatch_score = 0.0
    a.open_gap_score = -10.0
    a.extend_gap_score = -0.5
    return a


def align_and_build_map(struct_residues, uniprot_seq, aligner,
                        min_chain_identity=0.80, min_chain_coverage=0.50):
    """
    struct_residues: list of (chain, resid, icode, resname, aa1)
    Align EACH chain separately. Accept all chains that pass
    (identity >= min_chain_identity AND coverage >= min_chain_coverage).
    Merge their per-residue maps into a single combined map.

    Returns:
      {
        "chains_used": [chain_id, ...],
        "per_chain": { chain_id: {n_struct, n_mapped, n_match, identity, coverage} },
        "n_struct_total": int,
        "n_mapped_total": int,
        "n_match_total": int,
        "identity": float (weighted by mapped count),
        "coverage": float (overall),
        "resid_map": { (chain, resid, icode): uniprot_pos_0indexed },
        "rejected_chains": [chain_id, ...],
      }
    """
    by_chain = defaultdict(list)
    for entry in struct_residues:
        by_chain[entry[0]].append(entry)

    combined_map = {}
    per_chain = {}
    chains_used = []
    rejected = []
    n_struct_total = 0
    n_mapped_total = 0
    n_match_total = 0

    for chain_id, residues in by_chain.items():
        n_struct = len(residues)
        n_struct_total += n_struct
        if n_struct < 10:
            rejected.append((chain_id, "too_few_residues"))
            continue
        struct_seq = "".join(e[4] for e in residues)
        try:
            aln = aligner.align(struct_seq, uniprot_seq)[0]
        except Exception as exc:
            rejected.append((chain_id, f"align_err:{exc}"))
            continue

        s_blocks, u_blocks = aln.aligned
        mapped_this = {}
        n_match = 0
        for (s_beg, s_end), (u_beg, u_end) in zip(s_blocks, u_blocks):
            for off in range(min(int(s_end - s_beg), int(u_end - u_beg))):
                s_idx = int(s_beg) + off
                u_idx = int(u_beg) + off
                chain, resid, icode, _, aa1 = residues[s_idx]
                mapped_this[(chain, resid, icode)] = u_idx
                if uniprot_seq[u_idx] == aa1:
                    n_match += 1

        n_mapped = len(mapped_this)
        identity = (n_match / n_mapped) if n_mapped else 0.0
        coverage = (n_mapped / n_struct) if n_struct else 0.0
        per_chain[chain_id] = {
            "n_struct": n_struct,
            "n_mapped": n_mapped,
            "n_match": n_match,
            "identity": identity,
            "coverage": coverage,
        }
        if identity >= min_chain_identity and coverage >= min_chain_coverage:
            chains_used.append(chain_id)
            combined_map.update(mapped_this)
            n_mapped_total += n_mapped
            n_match_total += n_match
        else:
            rejected.append((chain_id, f"low_id_cov(id={identity:.3f},cov={coverage:.3f})"))

    if not chains_used:
        return None

    overall_identity = (n_match_total / n_mapped_total) if n_mapped_total else 0.0
    overall_coverage = (n_mapped_total / n_struct_total) if n_struct_total else 0.0
    return {
        "chains_used": chains_used,
        "per_chain": per_chain,
        "n_struct_total": n_struct_total,
        "n_mapped_total": n_mapped_total,
        "n_match_total": n_match_total,
        "identity": overall_identity,
        "coverage": overall_coverage,
        "resid_map": combined_map,
        "rejected_chains": rejected,
    }


# ---------- SIFTS ----------
SIFTS_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb_id}"

def fetch_sifts(pdb_id, cache_dir, retries=3):
    pdb_id = pdb_id.lower()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{pdb_id}.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(SIFTS_URL.format(pdb_id=pdb_id),
                                         headers={"User-Agent": "EZSpec-EXP005/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode())
            cache_path.write_text(json.dumps(data))
            return data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last_err = e
        except Exception as e:
            last_err = e
        time.sleep(1.5)
    raise RuntimeError(f"SIFTS fetch {pdb_id} failed: {last_err}")


def sifts_build_map(struct_residues, sifts_data, target_uniprot, pdb_id):
    """
    Build {(chain, resid, icode): uniprot_pos_0indexed} from SIFTS segments,
    applied to the residues we have in the structure.
    Returns (resid_map, chain_id_used, n_mapped, n_struct_in_chain, n_aa_match, notes).
    """
    pdb_id = pdb_id.lower()
    entry = sifts_data.get(pdb_id, {})
    upmap = entry.get("UniProt", {})
    if target_uniprot not in upmap:
        return None, None, 0, 0, 0, f"uniprot_{target_uniprot}_not_in_sifts"
    mappings = upmap[target_uniprot].get("mappings", [])

    # Group struct residues by chain for convenience
    by_chain = defaultdict(list)
    for c, rs, ic, rn, a1 in struct_residues:
        by_chain[c].append((rs, ic, rn, a1))

    # For each chain in our structure, try each SIFTS segment of the SAME chain
    # Use the first chain that yields the most matches.
    best_map = {}
    best_chain = None
    best_n_struct = 0
    best_n_match = 0
    for chain_id, residues in by_chain.items():
        # SIFTS may describe chain_id as-is (author chain)
        segs = [m for m in mappings if m.get("chain_id") == chain_id]
        if not segs:
            continue
        mapped = {}
        for m in segs:
            unp_start = int(m["unp_start"])
            start = m.get("start", {})
            author_start = start.get("author_residue_number")
            if author_start is None:
                continue
            # residue_number is the internal SIFTS counter (sequential modeled pos)
            # we need (unp_end - unp_start) span
            unp_end = int(m["unp_end"])
            seg_len = unp_end - unp_start + 1
            # Map: author residue (author_start + k) -> uniprot (unp_start + k)
            # for k in [0, seg_len), BUT the structure may have non-contiguous author
            # numbering inside. We'll apply delta to the RANGE and rely on alignment
            # agreement to catch any missed residues.
            delta = unp_start - author_start  # unp_pos (1-indexed) = author_resid + delta
            seg_author_min = author_start
            seg_author_max = author_start + seg_len - 1
            for resid, icode, rn, a1 in residues:
                if seg_author_min <= resid <= seg_author_max:
                    unp_pos_1idx = resid + delta
                    unp_pos_0idx = unp_pos_1idx - 1  # convert to 0-indexed
                    if 0 <= unp_pos_0idx < len(ENZ_ROWS[0]["Protein sequence"]):  # dummy len check replaced later
                        mapped[(chain_id, resid, icode)] = unp_pos_0idx
        if len(mapped) > len(best_map):
            best_map = mapped
            best_chain = chain_id
            best_n_struct = len(residues)
    # Final aa match count done outside using uniprot seq
    return best_map, best_chain, len(best_map), best_n_struct, 0, "ok"


# ---------- Main smoke test ----------
def run():
    aligner = make_aligner()
    sifts_cache = BASE / "data/pt_cache_dualgraph_allfix_unified/sifts_cache"

    print("=" * 80)
    print("EXP005 resid-map smoke test: 5 cases")
    print("=" * 80)

    results = []
    for dock, (enz_id, label, src_hint, pdb_id) in TEST_CASES.items():
        print(f"\n--- dock={dock}  enzyme={enz_id}  [{label}]  source_hint={src_hint} ---")
        row = ENZ_ROWS[enz_id]
        uniprot_id = row["uniprots"].strip()
        seq = row["Protein sequence"].strip()
        print(f"  uniprot={uniprot_id}, seq_len={len(seq)}, seq[0:20]={seq[:20]}")

        cpath = COMPLEX_DIR / f"{dock}.pdb"
        if not cpath.exists():
            print(f"  complex/{dock}.pdb NOT FOUND")
            continue

        struct_res = parse_complex_receptor_residues(cpath)
        print(f"  struct residues: {len(struct_res)}")
        if not struct_res:
            print("  EMPTY, skipping")
            continue

        # Alignment-based mapping (multi-chain aware)
        ali = align_and_build_map(struct_res, seq, aligner)
        if ali is None:
            print("  alignment FAILED")
            continue
        chain_summary = []
        for c, v in ali["per_chain"].items():
            chain_summary.append(
                "%s:m=%d/%d,id=%.3f,cov=%.3f" % (c, v["n_mapped"], v["n_struct"], v["identity"], v["coverage"])
            )
        print(f"  [alignment] chains_used={ali['chains_used']}")
        print(f"  [per_chain] {chain_summary}")
        print(f"  [total]    mapped={ali['n_mapped_total']}/{ali['n_struct_total']}, "
              f"identity={ali['identity']:.4f}, coverage={ali['coverage']:.4f}")
        if ali["rejected_chains"]:
            print(f"  rejected: {ali['rejected_chains']}")

        # SIFTS cross-check if experimental_pdb
        sifts_agree = None
        if src_hint == "experimental_pdb" and pdb_id is not None:
            try:
                sifts_data = fetch_sifts(pdb_id, sifts_cache)
                if sifts_data is None:
                    print(f"  [SIFTS] no entry for {pdb_id}")
                else:
                    sm_map, sm_chain, sm_n, sm_struct, _, sm_notes = sifts_build_map(
                        struct_res, sifts_data, uniprot_id, pdb_id)
                    # Recompute aa match for SIFTS map using the real seq
                    aa_ok = 0
                    for (c, r, ic), up in sm_map.items():
                        rn = next((e[3] for e in struct_res if (e[0], e[1], e[2]) == (c, r, ic)), None)
                        if rn is None:
                            continue
                        if 0 <= up < len(seq) and seq[up] == THREE_TO_ONE.get(rn, "?"):
                            aa_ok += 1
                    sm_ident = aa_ok / sm_n if sm_n else 0.0
                    print(f"  [SIFTS]     chain={sm_chain}, mapped={sm_n}/{sm_struct}, "
                          f"aa_match={aa_ok}/{sm_n} ({100*sm_ident:.1f}%), notes={sm_notes}")
                    # Compare with alignment
                    if sm_map and ali["resid_map"]:
                        common = set(sm_map.keys()) & set(ali["resid_map"].keys())
                        agree = sum(1 for k in common if sm_map[k] == ali["resid_map"][k])
                        sifts_agree = (agree, len(common))
                        print(f"  [compare] SIFTS vs alignment on {len(common)} common keys: "
                              f"{agree} agree ({100*agree/len(common) if common else 0:.1f}%)")
            except Exception as e:
                print(f"  [SIFTS] ERROR: {e}")

        # Verify using POCKET residues (these are what GVP actually uses)
        ppath = POCKET_DIR / f"{dock}.pdb"
        pocket_residues = parse_complex_receptor_residues(ppath)  # same parser works
        if pocket_residues:
            hits = 0
            miss_examples = []
            for c, r, ic, rn, a1 in pocket_residues:
                up = ali["resid_map"].get((c, r, ic))
                if up is None:
                    continue
                if 0 <= up < len(seq) and seq[up] == a1:
                    hits += 1
                else:
                    if len(miss_examples) < 3:
                        miss_examples.append((c, r, ic, rn, seq[up] if 0 <= up < len(seq) else "?"))
            total_mapped = sum(1 for e in pocket_residues
                               if (e[0], e[1], e[2]) in ali["resid_map"])
            total = len(pocket_residues)
            print(f"  [pocket check] pocket N={total}, mapped={total_mapped}, "
                  f"aa_match={hits}/{total_mapped} "
                  f"({100*hits/total_mapped if total_mapped else 0:.1f}%)")
            if miss_examples:
                print(f"  pocket miss examples: {miss_examples}")

        results.append({
            "dock": dock,
            "enzyme": enz_id,
            "identity": ali["identity"],
            "coverage": ali["coverage"],
            "sifts_agree": sifts_agree,
        })

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for r in results:
        print(f"  dock={r['dock']:6d} enz={r['enzyme']:4d}: "
              f"identity={r['identity']:.4f} coverage={r['coverage']:.4f} "
              f"sifts_agree={r['sifts_agree']}")


if __name__ == "__main__":
    run()
'''


def main():
    import subprocess, sys
    print("[local driver] running smoke test on server...")
    result = subprocess.run(
        ["ssh", "autodl-4x5090-bj",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=600,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write("\n[STDERR]\n")
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
