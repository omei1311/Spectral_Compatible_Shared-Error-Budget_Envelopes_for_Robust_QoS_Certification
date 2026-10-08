# Final dataset freeze (V6.1 submission)

Frozen: 2026-10-07T15:13:46

Policy: FINAL DATASET for the V6.1 submission. After this freeze, table numbers may only be corrected for proven errors; no seed, parameter or protocol changes; no re-optimization of results.

## Pinned final table numbers

- **Table 3 (sensitivity, 1080 units)**: {"feasible": 472, "infeasible": 608, "D": 326, "old": 380, "Q_G": 390, "L_quad": 408, "added": 10}
- **Table 4 (validation, 540 units)**: {"feasible": 159, "D": 111, "old": 125, "Q_G": 126, "L_quad": 138}
- **Table 6 (heterogeneous links, 540 units)**: {"feasible": 211, "old": 156, "Q_G": 163, "L_quad": 186, "joint": 190, "flips": 7}
- **Table 7 (joint/staged)**: {"exploratory": "feas 472 | Q_G 390 | Lq 408 | joint 424", "validation": "feas 159 | 126 | 138 | 139", "heterogeneous": "feas 211 | 163 | 186 | 190", "confirmation": "final accepted 202 all pipelines; calls 1424->394", "speedup": "2.12x (25/30) / 2.01x (11/15) vs tolerance-stopped ref"}
- **Table 8 (controlled scaling)**: {"bound_ms_per_user": {"4": 0.041, "8": 0.05, "16": 0.044, "32": 0.042}, "ref_ms_per_user": {"4": 0.356, "8": 0.368, "16": 0.381, "32": 0.435}, "pipeline_s": {"4": [0.0192, 0.0193, 0.0254], "8": [0.0378, 0.0165, 0.0236], "16": [0.0563, 0.0175, 0.0161], "32": [0.0859, 0.031, 0.031]}, "order": ["TR_pure", "L_quad+ref", "Lq->QG+ref"]}
- **Table 9 (SDP, 2160 problems)**: {"L_minus_qstar_max": "-1.2e-13", "U_deficit_max": "1.3e-12", "decisions": "1247 feasible, 913 infeasible, 0 unresolved", "witness_U_minus_qstar": "median 6.7e-6, max 1.6e-4", "max_abs_t_minus_qstar": "5.6e-8"}

## Result directories (content-hashed)

- `results/`: 9 files
- `results_review/`: 31 files
- `results_v41/`: 27 files
- `results_v41_optimization/`: 18 files
- `results_v41_sdp/`: 10 files
- `results_v41_sdp_corrected/`: 10 files
- `results_v6/`: 5 files
- `results_v61/`: 3 files

Full hashes: final_freeze_manifest.json (code, paper, every result file).
