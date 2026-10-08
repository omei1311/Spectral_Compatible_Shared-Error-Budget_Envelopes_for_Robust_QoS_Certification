# -*- coding: utf-8 -*-
"""Freeze the final dataset: code hashes, result-directory inventory and the
pinned final table numbers for V6.1 submission.  Policy after freeze:
error-correction only; no seed/parameter changes, no result optimization."""
import hashlib
import json
import os
from datetime import datetime

import numpy as np

CODE_FILES = [
    "core.py", "experiments_a.py", "experiments_b.py", "experiments_c.py",
    "experiments_v41.py", "experiments_v6.py", "optimization_v41.py",
    "sdp_baseline.py", "sdp_baseline_v2.py", "audit_and_stress.py",
    "run_experiments.py", "run_v41.py", "run_optimization.py",
    "run_sdp.py", "run_sdp_v2.py", "v61_data_audit.py",
    "paper/build_v41_docx.py", "paper/build_v5_docx.py",
    "paper/build_v6_docx.py", "paper/build_v6_1_docx.py",
]
RESULT_DIRS = ["results", "results_review", "results_v41",
               "results_v41_optimization", "results_v41_sdp",
               "results_v41_sdp_corrected", "results_v6", "results_v61"]
PAPER_FILES = ["paper/paper_v6_1.md",
               "paper/Sharp_Robust_QoS_Certification_from_Six-Scalar_"
               "Spectral_Summaries_V6.1_EN.docx"]

FINAL_TABLE_NUMBERS = {
    "Table 3 (sensitivity, 1080 units)": {
        "feasible": 472, "infeasible": 608, "D": 326, "old": 380,
        "Q_G": 390, "L_quad": 408, "added": 10},
    "Table 4 (validation, 540 units)": {
        "feasible": 159, "D": 111, "old": 125, "Q_G": 126, "L_quad": 138},
    "Table 6 (heterogeneous links, 540 units)": {
        "feasible": 211, "old": 156, "Q_G": 163, "L_quad": 186,
        "joint": 190, "flips": 7},
    "Table 7 (joint/staged)": {
        "exploratory": "feas 472 | Q_G 390 | Lq 408 | joint 424",
        "validation": "feas 159 | 126 | 138 | 139",
        "heterogeneous": "feas 211 | 163 | 186 | 190",
        "confirmation": "final accepted 202 all pipelines; calls 1424->394",
        "speedup": "2.12x (25/30) / 2.01x (11/15) vs tolerance-stopped ref"},
    "Table 8 (controlled scaling)": {
        "bound_ms_per_user": {"4": 0.041, "8": 0.050, "16": 0.044,
                              "32": 0.042},
        "ref_ms_per_user": {"4": 0.356, "8": 0.368, "16": 0.381,
                            "32": 0.435},
        "pipeline_s": {"4": [0.0192, 0.0193, 0.0254],
                       "8": [0.0378, 0.0165, 0.0236],
                       "16": [0.0563, 0.0175, 0.0161],
                       "32": [0.0859, 0.0310, 0.0310]},
        "order": ["TR_pure", "L_quad+ref", "Lq->QG+ref"]},
    "Table 9 (SDP, 2160 problems)": {
        "L_minus_qstar_max": "-1.2e-13", "U_deficit_max": "1.3e-12",
        "decisions": "1247 feasible, 913 infeasible, 0 unresolved",
        "witness_U_minus_qstar": "median 6.7e-6, max 1.6e-4",
        "max_abs_t_minus_qstar": "5.6e-8"},
}


def sha(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    os.makedirs("docs", exist_ok=True)
    man = dict(
        frozen_at=datetime.now().isoformat(timespec="seconds"),
        policy="FINAL DATASET for the V6.1 submission. After this freeze, "
               "table numbers may only be corrected for proven errors; no "
               "seed, parameter or protocol changes; no re-optimization of "
               "results.",
        code_hashes={f: sha(f) for f in CODE_FILES if os.path.exists(f)},
        paper_hashes={f: sha(f) for f in PAPER_FILES if os.path.exists(f)},
        result_dirs={},
        final_table_numbers=FINAL_TABLE_NUMBERS,
        experiment_seeds=dict(
            A1_direction=20261006, easy=[11, 22, 33, 44, 55],
            sensitivity=[11, 22, 33, 44, 55, 66, 77, 88, 99, 110],
            validation=[211, 212, 213, 214, 215],
            confirmation=[221, 222, 223, 224, 225],
            heterogeneous=[301, 302, 303, 304, 305],
            scaling=[311, 312, 313],
            v41_construction_battery=202610061),
        parameters=dict(
            gammas=[0.5, 1.0, 2.0], rhos=[0.06, 0.15, 0.30],
            easy_rhos=[0.01, 0.03, 0.06], noise=0.05,
            interval_J=256, ref_mu_rel_tol=1e-12,
            tau="1e-8*max(1,|c_q|,2r||g||,r^2||A||_2)"))
    for d in RESULT_DIRS:
        if os.path.isdir(d):
            files = []
            for root, _, names in os.walk(d):
                for n in sorted(names):
                    p = os.path.join(root, n)
                    files.append((os.path.relpath(p, d), sha(p)))
            man["result_dirs"][d] = dict(n_files=len(files),
                                         files=files)
    out = "docs/final_freeze_manifest.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(man, f, indent=1)
    # compact md summary
    with open("docs/final_freeze_manifest.md", "w", encoding="utf-8") as f:
        f.write("# Final dataset freeze (V6.1 submission)\n\n")
        f.write(f"Frozen: {man['frozen_at']}\n\nPolicy: {man['policy']}\n\n")
        f.write("## Pinned final table numbers\n\n")
        for k, v in FINAL_TABLE_NUMBERS.items():
            f.write(f"- **{k}**: {json.dumps(v, ensure_ascii=False)}\n")
        f.write("\n## Result directories (content-hashed)\n\n")
        for d, info in man["result_dirs"].items():
            f.write(f"- `{d}/`: {info['n_files']} files\n")
        f.write("\nFull hashes: final_freeze_manifest.json "
                "(code, paper, every result file).\n")
    print("freeze manifest written:", out)
    print("code files hashed:", len(man["code_hashes"]),
          "| paper:", len(man["paper_hashes"]),
          "| result files:",
          sum(i["n_files"] for i in man["result_dirs"].values()))


if __name__ == "__main__":
    main()
