#!/usr/bin/env bash
set -euo pipefail

Q10_DIR="/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528"
Q10_EXP="/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528"
EXP008="/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q02_sequence_similarity_split/EXP008_random_gdtable_b80_retry_20260525"
PY="/root/miniconda3/bin/python"
CF_PY="/root/autodl-tmp/envs/colabfold_q10/bin/python"

echo "[Q10] Step 1: collect ColabFold manifest"
"${CF_PY}" "${Q10_EXP}/scripts/q10_collect_colabfold_manifest_20260528.py" \
  --q10-dir "${Q10_DIR}" \
  --link-selected

SUMMARY="${Q10_DIR}/structures/manifests/q10_colabfold_structure_summary_msa_m1r1_v1.json"
"${PY}" -c "import json,sys; p='${SUMMARY}'; d=json.load(open(p,encoding='utf-8')); print('[Q10] ColabFold summary:', d); assert d.get('total')==204, d; assert d.get('complete')==204, d"

echo "[Q10] Step 2: build HEM/Fe receptor PDBs"
"${PY}" "${Q10_EXP}/scripts/q10_build_heme_fe_receptors_20260528.py" \
  --q10-dir "${Q10_DIR}" \
  --alphafill-summary "/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/alphafill_p450_216_20260528/manifests/alphafill_hem_candidate_summary_20260528.csv" \
  --alphafill-cif-dir "/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/alphafill_p450_216_20260528/cif" \
  --alignment-mode ce \
  --max-rmsd 10

RECEPTOR_SUMMARY="${Q10_DIR}/structures/manifests/q10_heme_fe_receptor_summary_template_v1.json"
"${PY}" -c "import json,sys; p='${RECEPTOR_SUMMARY}'; d=json.load(open(p,encoding='utf-8')); print('[Q10] HEM/Fe receptor summary:', d); c=d.get('counts',{}); assert c.get('target_complete')==204, d; assert c.get('failed')==0, d"

echo "[Q10] Step 3: prepare Uni-Dock boxes from HEM/Fe receptors"
"${PY}" "${Q10_EXP}/scripts/q10_prepare_unidock_boxes_20260528.py" \
  --q10-dir "${Q10_DIR}" \
  --allow-non-good-fe-cys

BOX_SUMMARY="${Q10_DIR}/docking/manifests/q10_unidock_box_summary_msa_m1r1_v1.json"
"${PY}" -c "import json,sys; p='${BOX_SUMMARY}'; d=json.load(open(p,encoding='utf-8')); print('[Q10] box summary:', d); assert d.get('complete_receptors')==204, d"

echo "[Q10] Step 4: run Uni-Dock for all HEM/Fe receptors"
"${PY}" "${Q10_EXP}/scripts/q10_run_unidock_batch_20260528.py" \
  --q10-dir "${Q10_DIR}" \
  --gpu-ids "0,1" \
  --rerun-existing \
  --fresh-manifest

DOCK_SUMMARY="${Q10_DIR}/docking/manifests/q10_unidock_results_summary_msa_m1r1_v1.json"
"${PY}" -c "import json,sys; p='${DOCK_SUMMARY}'; d=json.load(open(p,encoding='utf-8')); print('[Q10] Uni-Dock summary:', d); total=d.get('total_result_rows'); statuses=d.get('by_status',{}); complete=sum(v for k,v in statuses.items() if str(k).startswith('complete')); assert total==204, d; assert complete==204, d"

echo "[Q10] Step 5: build model PT cache"
"${PY}" "${Q10_EXP}/scripts/q10_build_model_pt_cache_20260528.py" \
  --q10-dir "${Q10_DIR}"

CACHE_DIR="${Q10_DIR}/model_cache/q10_unidock_pt_msa_m1r1_v1"
CACHE_SUMMARY="${CACHE_DIR}/manifests/q10_model_cache_summary_test_msa_m1r1_v1.json"
"${PY}" -c "import json,sys; p='${CACHE_SUMMARY}'; d=json.load(open(p,encoding='utf-8')); print('[Q10] PT cache summary:', d); assert d.get('samples_written')==204, d; assert d.get('failures')==0, d"

echo "[Q10] Step 6: validate cache with EXP008 dataset"
"${PY}" "${Q10_EXP}/scripts/q10_validate_model_cache_20260528.py" \
  --cache-dir "${CACHE_DIR}" \
  --exp008-scripts "${EXP008}/scripts" \
  --split test \
  --batch-size 2

echo "[Q10] Step 7: infer scores with Q2 random best checkpoint"
"${PY}" "${Q10_EXP}/scripts/q10_predict_with_q2_random_ckpt_20260528.py" \
  --cache-dir "${CACHE_DIR}" \
  --checkpoint "${EXP008}/results/checkpoints/pt-Q2_EXP008_random_gdtable_b80_retry_full_20260525_174150-ep79-auc0.9316.ckpt" \
  --config "${EXP008}/configs/config.yml" \
  --exp008-scripts "${EXP008}/scripts" \
  --output-dir "${Q10_DIR}/results/model_scores" \
  --split test \
  --batch-size 16 \
  --device cuda:0

SCORE_SUMMARY="${Q10_DIR}/results/model_scores/q10_model_scores_summary_test_msa_m1r1_v1.json"
"${PY}" -c "import json,sys; p='${SCORE_SUMMARY}'; d=json.load(open(p,encoding='utf-8')); print('[Q10] Score summary:', d); assert d.get('n_samples')==204, d; assert d.get('load_state_dict_missing_count')==0, d; assert d.get('load_state_dict_unexpected_count')==0, d"

echo "[Q10] Step 8: audit candidate-id consistency across outputs"
"${PY}" "${Q10_EXP}/scripts/q10_audit_full_outputs_20260528.py" \
  --q10-dir "${Q10_DIR}" \
  --split test

echo "[Q10] Full pipeline finished."
