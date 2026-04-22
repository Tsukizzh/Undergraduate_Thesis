"""
Preflight: verify paper checkpoint can be loaded strict=True into SS(config).
Local CPU only. No GPU, no dataset, no server.

Exits 0 on clean load (possibly after one explicit prefix normalization).
Exits nonzero on any remaining key/shape mismatch.
"""
from __future__ import annotations
import os, sys, json, time, traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CKPT = ROOT / "saved_model/model/run_0/models/best-checkpoint.ckpt"
CFG = ROOT / "saved_model/model/run_0/complete-full-random-all-0-complex.yml"
REPORT_OUT = ROOT / "scratch/ckpt_preflight_report.json"

# Add local src/ to path
sys.path.insert(0, str(ROOT / "src"))


def log(msg):
    print(msg, flush=True)


def main():
    report = {}

    log("=" * 70)
    log("Paper Checkpoint Preflight (strict load check)")
    log("=" * 70)

    # 1. Library versions
    import torch
    import pytorch_lightning as pl
    log(f"Python         : {sys.version.split()[0]}")
    log(f"torch          : {torch.__version__}")
    log(f"pytorch_lightning: {pl.__version__}")
    report["versions"] = {
        "python": sys.version.split()[0],
        "torch": torch.__version__,
        "pytorch_lightning": pl.__version__,
    }

    # 2. Load config
    log(f"\n[1/5] Loading config: {CFG.relative_to(ROOT)}")
    from utils import load_config  # EasyDict loader from paper's src
    cfg = load_config(str(CFG))
    log(f"  top-level keys: {list(cfg.keys())}")
    report["config_top_keys"] = list(cfg.keys())

    # 3. Load checkpoint (just the state_dict, don't instantiate yet)
    log(f"\n[2/5] Loading ckpt: {CKPT.relative_to(ROOT)}")
    t0 = time.time()
    try:
        ckpt = torch.load(str(CKPT), map_location="cpu", weights_only=False)
    except Exception as e:
        log(f"  ckpt load FAILED: {e}")
        report["ckpt_load_error"] = str(e)
        with open(REPORT_OUT, "w") as f:
            json.dump(report, f, indent=2, default=str)
        return 2
    log(f"  loaded in {time.time()-t0:.1f}s")
    log(f"  ckpt top-level keys: {list(ckpt.keys()) if isinstance(ckpt, dict) else type(ckpt).__name__}")
    report["ckpt_top_keys"] = list(ckpt.keys()) if isinstance(ckpt, dict) else [type(ckpt).__name__]

    if isinstance(ckpt, dict) and "state_dict" in ckpt:
        ckpt_sd = ckpt["state_dict"]
        log(f"  using ckpt['state_dict'] ({len(ckpt_sd)} keys)")
        report["ckpt_epoch"] = ckpt.get("epoch")
        report["ckpt_global_step"] = ckpt.get("global_step")
        report["ckpt_pl_version"] = ckpt.get("pytorch-lightning_version")
    else:
        ckpt_sd = ckpt
        log(f"  raw state_dict ({len(ckpt_sd)} keys)")
    report["ckpt_state_dict_keys"] = len(ckpt_sd)

    # Sample a few keys
    ckpt_keys_sample = list(ckpt_sd.keys())[:10]
    log(f"  sample keys: {ckpt_keys_sample}")
    report["ckpt_sample_keys"] = ckpt_keys_sample

    # 4. Instantiate SS model
    log(f"\n[3/5] Instantiating SS(config)...")
    t0 = time.time()
    try:
        from Models.ss import SS
        model = SS(cfg)
    except Exception as e:
        log(f"  SS instantiation FAILED: {e}")
        traceback.print_exc()
        report["ss_instantiation_error"] = str(e)
        with open(REPORT_OUT, "w") as f:
            json.dump(report, f, indent=2, default=str)
        return 3
    log(f"  built in {time.time()-t0:.1f}s")
    model_sd = model.state_dict()
    log(f"  model state_dict keys: {len(model_sd)}")
    model_keys_sample = list(model_sd.keys())[:10]
    log(f"  sample keys: {model_keys_sample}")
    report["model_state_dict_keys"] = len(model_sd)
    report["model_sample_keys"] = model_keys_sample

    # 5. Key diff
    log(f"\n[4/5] Key diff analysis")
    model_keys = set(model_sd.keys())
    ckpt_keys = set(ckpt_sd.keys())
    missing = sorted(model_keys - ckpt_keys)
    unexpected = sorted(ckpt_keys - model_keys)
    matched = model_keys & ckpt_keys
    log(f"  matched:    {len(matched)}")
    log(f"  missing:    {len(missing)}  (in model, not in ckpt)")
    log(f"  unexpected: {len(unexpected)}  (in ckpt, not in model)")
    report["matched_keys"] = len(matched)
    report["missing_keys_count"] = len(missing)
    report["unexpected_keys_count"] = len(unexpected)
    report["missing_keys_first20"] = missing[:20]
    report["unexpected_keys_first20"] = unexpected[:20]

    # 6. Candidate prefix rewrite detection
    normalized = False
    if missing and unexpected and len(missing) == len(unexpected):
        # Check if all unexpected keys share a common prefix that, if stripped, gives missing keys
        prefixes = {k.split(".", 1)[0] for k in unexpected}
        if len(prefixes) == 1:
            p = prefixes.pop() + "."
            stripped = {k[len(p):] if k.startswith(p) else k for k in unexpected}
            if stripped == set(missing):
                log(f"  ** candidate prefix normalization: strip '{p}' from ckpt keys **")
                log(f"  ckpt keys all have prefix '{p}', model keys don't — deterministic rewrite")
                # Apply normalization
                ckpt_sd = {(k[len(p):] if k.startswith(p) else k): v for k, v in ckpt_sd.items()}
                normalized = True
                report["prefix_normalization"] = {"stripped": p}

    # Shape check on matched keys
    log(f"\n[5/5] Shape check on matched keys")
    shape_mismatches = []
    for k in sorted(model_keys & set(ckpt_sd.keys())):
        m_shape = tuple(model_sd[k].shape)
        c_shape = tuple(ckpt_sd[k].shape)
        if m_shape != c_shape:
            shape_mismatches.append({"key": k, "model": m_shape, "ckpt": c_shape})
    log(f"  matched shape mismatches: {len(shape_mismatches)}")
    for sm in shape_mismatches[:5]:
        log(f"    {sm['key']}: model={sm['model']} ckpt={sm['ckpt']}")
    report["shape_mismatches_count"] = len(shape_mismatches)
    report["shape_mismatches_first5"] = shape_mismatches[:5]

    # Strict load attempt
    log(f"\n[final] strict=True load attempt")
    try:
        result = model.load_state_dict(ckpt_sd, strict=True)
        log(f"  [OK] strict=True PASSED")
        report["strict_load"] = "PASSED"
        clean = True
    except RuntimeError as e:
        log(f"  [FAIL] strict=True FAILED: {str(e)[:300]}")
        report["strict_load_error"] = str(e)[:500]
        # Fallback info via strict=False
        r = model.load_state_dict(ckpt_sd, strict=False)
        log(f"  strict=False: missing={len(r.missing_keys)}, unexpected={len(r.unexpected_keys)}")
        report["strict_false_missing"] = len(r.missing_keys)
        report["strict_false_unexpected"] = len(r.unexpected_keys)
        clean = False

    # Summary
    log("\n" + "=" * 70)
    if clean and not shape_mismatches:
        if normalized:
            log("RESULT: PASS (after explicit prefix normalization)")
            report["result"] = "PASS_NORMALIZED"
            exit_code = 0
        else:
            log("RESULT: PASS (clean strict load)")
            report["result"] = "PASS_CLEAN"
            exit_code = 0
    else:
        log("RESULT: FAIL")
        report["result"] = "FAIL"
        exit_code = 1
    log("=" * 70)

    with open(REPORT_OUT, "w") as f:
        json.dump(report, f, indent=2, default=str)
    log(f"\nReport: {REPORT_OUT}")
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
