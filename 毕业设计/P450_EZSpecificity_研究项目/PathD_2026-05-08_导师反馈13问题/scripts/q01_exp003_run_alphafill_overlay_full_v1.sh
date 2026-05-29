#!/usr/bin/env bash
set -euo pipefail

cd /root/autodl-tmp/EZSpecificity/PathD/P450/scripts

BASE=/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch
BASE_CACHE=$BASE/exp001_esibank_aonly_baseline/pt_cache/full_legacy_cuda_20260526_161909_per_sample_graph_v1
OUT_CACHE=$BASE/exp003_full_p450_fe_heme_overlay/pt_cache/full_exp003_alphafill_heme_fe_v1
TRAINABLE=$BASE/exp003_full_p450_fe_heme_overlay/manifests/exp003_p450389_trainable_split_sample_manifest_20260527_214331.csv
OFF_MAN=$BASE/exp003_full_p450_fe_heme_overlay/official_esibank_p450_pdb_20260528_package/exp003_p450_official_esibank_pdb_manifest_20260528.csv
OFF_PDB=$BASE/exp003_full_p450_fe_heme_overlay/official_esibank_p450_pdb_20260528_package/pdb
AF_SUM=$BASE/exp003_full_p450_fe_heme_overlay/alphafill_p450_216_20260528/manifests/alphafill_hem_candidate_summary_20260528.csv
AF_CIF=$BASE/exp003_full_p450_fe_heme_overlay/alphafill_p450_216_20260528/cif
REPORT=$BASE/exp003_full_p450_fe_heme_overlay/overlay_reports/alphafill_full_v1

mkdir -p "$REPORT"

for i in $(seq 0 7); do
  mkdir -p "$REPORT/shard_$i"
  /root/miniconda3/bin/python q01_exp003_make_alphafill_heme_overlay_samples_20260528.py \
    --base-cache "$BASE_CACHE" \
    --out-cache "$OUT_CACHE" \
    --trainable-manifest "$TRAINABLE" \
    --official-pdb-manifest "$OFF_MAN" \
    --official-pdb-dir "$OFF_PDB" \
    --alphafill-summary "$AF_SUM" \
    --alphafill-cif-dir "$AF_CIF" \
    --report-dir "$REPORT/shard_$i" \
    --shard-count 8 \
    --shard-index "$i" \
    --no-manifest-update \
    --progress-every 200 \
    > "$REPORT/shard_$i/run.log" 2>&1 &
  echo $! > "$REPORT/shard_$i/pid.txt"
done

wait
