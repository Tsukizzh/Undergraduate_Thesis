#!/usr/bin/env python3
"""
Uni-Dock batch docking: 2 GPU workers, one receptor at a time, all ligands via --gpu_batch.
Results in results/<enz_XXXXXX>/. Resume by checking existing output files.
"""

import os, sys, csv, time, logging, subprocess
from pathlib import Path
from collections import defaultdict
from threading import Thread
from queue import Queue

UNIDOCK_BIN = "/opt/conda/envs/unidock/bin/unidock"
SEARCH_MODE = "fast"  # fast: exhaustiveness=128, max_step=20 (~0.10s/ligand)
NUM_MODES = 1
SCORING = "vina"
BOX_SIZE = 22.5

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()],
)
log = logging.getLogger()


def build_jobs(base_dir):
    reg_dir = base_dir / "registries"

    enz_reg = {}
    with open(reg_dir / "enzyme_registry.csv") as f:
        for row in csv.DictReader(f):
            enz_reg[int(row["enzyme_index"])] = row

    jobs = defaultdict(lambda: {"ligands": []})
    skipped = 0
    with open(reg_dir / "pair_registry_batch_1.csv") as f:
        for row in csv.DictReader(f):
            eidx = int(row["enzyme_index"])
            sidx = int(row["substrate_index"])

            rec = base_dir / "receptors_pdbqt" / f"ENZ_G{eidx+1:06d}.pdbqt"
            lig = base_dir / "ligands_pdbqt" / f"CMP_G{sidx+1:06d}.pdbqt"
            if not rec.exists() or not lig.exists():
                skipped += 1
                continue

            enz = enz_reg.get(eidx, {})
            fe_x, fe_y, fe_z = enz.get("heme_fe_x", ""), enz.get("heme_fe_y", ""), enz.get("heme_fe_z", "")
            if not fe_x or not fe_y or not fe_z:
                skipped += 1
                continue

            jobs[eidx]["receptor"] = str(rec)
            jobs[eidx]["center"] = (float(fe_x), float(fe_y), float(fe_z))
            jobs[eidx]["ligands"].append(str(lig))

    total = sum(len(j["ligands"]) for j in jobs.values())
    log.info(f"Built {len(jobs)} receptor jobs, {total} pairs, {skipped} skipped")
    return dict(jobs)


def run_one_receptor(eidx, job, gpu_id, results_dir):
    out_dir = results_dir / f"enz_{eidx:06d}"
    out_dir.mkdir(parents=True, exist_ok=True)

    all_ligands = job["ligands"]
    cx, cy, cz = job["center"]

    # skip ligands that already have output
    todo = []
    for lig in all_ligands:
        out_file = out_dir / f"{Path(lig).stem}_out.pdbqt"
        if not out_file.exists():
            todo.append(lig)

    if not todo:
        return len(all_ligands), 0

    # write ligand list to temp file for --ligand_index
    ligand_index = out_dir / "_ligand_index.txt"
    ligand_index.write_text("\n".join(todo) + "\n")

    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    cmd = [
        UNIDOCK_BIN,
        "--receptor", job["receptor"],
        "--ligand_index", str(ligand_index),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(BOX_SIZE), "--size_y", str(BOX_SIZE), "--size_z", str(BOX_SIZE),
        "--scoring", SCORING,
        "--search_mode", SEARCH_MODE,
        "--num_modes", str(NUM_MODES),
        "--dir", str(out_dir),
    ]

    ok_count = 0
    fail_count = 0
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
        try:
            stdout, stderr = proc.communicate(timeout=1800)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            log.warning(f"GPU{gpu_id} enz={eidx} timeout, killed")
            fail_count = len(todo)
            return fail_count, fail_count

        if proc.returncode == 0:
            for lig in todo:
                if (out_dir / f"{Path(lig).stem}_out.pdbqt").exists():
                    ok_count += 1
                else:
                    fail_count += 1
        else:
            err = (stderr.decode(errors="replace") or stdout.decode(errors="replace") or "unknown")[:200]
            log.warning(f"GPU{gpu_id} enz={eidx} failed: {err}")
            fail_count = len(todo)
    except Exception as e:
        log.warning(f"GPU{gpu_id} enz={eidx} error: {e}")
        fail_count = len(todo)
    finally:
        ligand_index.unlink(missing_ok=True)

    already_done = len(all_ligands) - len(todo)
    return ok_count + already_done, fail_count


def worker(gpu_id, job_queue, results_dir, stats):
    total_ok = 0
    total_fail = 0
    total_receptors = 0

    while True:
        try:
            eidx, job = job_queue.get(timeout=3)
        except Exception:
            break
        try:
            ok, fail = run_one_receptor(eidx, job, gpu_id, results_dir)
            total_ok += ok
            total_fail += fail
            total_receptors += 1
            if total_receptors % 10 == 0:
                log.info(f"GPU{gpu_id}: {total_receptors} receptors done, {total_ok} ok, {total_fail} fail")
        except Exception as e:
            log.error(f"GPU{gpu_id} enz={eidx} crashed: {e}")
        finally:
            job_queue.task_done()

    stats[gpu_id] = {"ok": total_ok, "fail": total_fail, "receptors": total_receptors}
    log.info(f"GPU{gpu_id} finished: {total_receptors} receptors, {total_ok} ok, {total_fail} fail")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", required=True)
    parser.add_argument("--shutdown", action="store_true")
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    results_dir = base_dir / "results"
    results_dir.mkdir(exist_ok=True)

    log.info(f"search_mode={SEARCH_MODE}, box={BOX_SIZE}, scoring={SCORING}, num_modes={NUM_MODES}")

    jobs = build_jobs(base_dir)
    if not jobs:
        log.error("No jobs found")
        return

    # sort by ligand count descending (big jobs first for better load balance)
    sorted_jobs = sorted(jobs.items(), key=lambda x: len(x[1]["ligands"]), reverse=True)

    q = Queue()
    for eidx, job in sorted_jobs:
        q.put((eidx, job))

    stats = {}
    t0 = time.time()
    threads = []
    for gpu_id in range(2):
        t = Thread(target=worker, args=(gpu_id, q, results_dir, stats), daemon=True)
        t.start()
        threads.append(t)
        log.info(f"Started worker on GPU {gpu_id}")

    q.join()
    for t in threads:
        t.join(timeout=10)

    elapsed = time.time() - t0
    total_ok = sum(s["ok"] for s in stats.values())
    total_fail = sum(s["fail"] for s in stats.values())
    log.info(f"DONE. {total_ok} ok, {total_fail} fail, {elapsed:.0f}s ({elapsed/60:.1f}min)")

    if args.shutdown:
        log.info("Shutting down in 60s...")
        subprocess.run(["shutdown", "-h", "+1"], capture_output=True)


if __name__ == "__main__":
    main()
