# Spectral-Compatible Shared-Error-Budget Envelopes for Robust QoS Certification — lightweight simulation

**Current manuscript**: `paper/Sharp_Robust_QoS_Certification_from_Six-Scalar_Spectral_Summaries_V6.5_EN.docx`
(paper/paper_v6_5.md + paper/build_v6_5_docx.py; docs/V6.5_修改说明.md
records the V6.4→V6.5 repositioning: partial-information extremal framing,
scope subsection with restricted-class monotonicity, bounded-nominal-norm
bracketed experiments → results_v65/, L_quad richer-interface framing, and
the strict Appendix-A novelty-defense table). Earlier versions
(v3 … V6.4) are preserved in `paper/`.

## Directory map

| location | contents |
|---|---|
| `paper/` | all manuscript versions (v3 draft, v4 Chinese, V4.1 EN, V5 EN) with their Markdown sources and builders (`build_v41_docx.py`, `build_v5_docx.py`, `docx_io.py`, `make_revision.py`, `paper_content_*.py`, `render_v41_word.ps1`) |
| `docs/` | revision notes (修改说明_v4 / _v4.1 / _v5), review and audit reports (review_report.md, v4_review.md, 实验审查报告_v41.md, 优化报告_v41.md, SDP基线审查报告.md / _修正版.md, findings.md), paper-ready English material (english_*.md), and the V6 preparation set (V6_novelty_audit_and_plan.md, V6_replacement_text.md) |
| `logs/` | run logs and text extractions (draft_text.txt, revised_check.txt, v5_check.txt, run_*_log.txt) |
| root | experiment code only: core.py, experiments_*.py, run_*.py entry points, audit_and_stress.py, optimization_v41.py, sdp_baseline*.py, verify_*.py, review_v4_checks.py, utils.py, figures.py |
| `results/`, `results_review/`, `results_v41/`, `results_v41_optimization/`, `results_v41_sdp/`, `results_v41_sdp_corrected/` | immutable experiment outputs, one directory per round, never overwritten |

Reproducible NumPy/SciPy/Matplotlib simulation of direction-envelope
correctness, six-scalar information-class construction, and radius-safety
certification.  No CVXPY/SDP, no GPU, no parallelism.

## Run

```powershell
python run_experiments.py
```

Additional audit and fixed-grid sensitivity experiment: `python audit_and_stress.py`.
It preserves `results/`, writes to `results_review/`, and checks final pipeline
decisions against an early-stopping exact-reference baseline. See `review_report.md`
for corrections to the original findings, including the separate boundary Q_partial.

**S-procedure-SDP baseline (v2, corrected)**: `python run_sdp_v2.py`
(`--quick` = tests + smallcheck) writes into `results_v41_sdp_corrected/`.
Fixes the v1 decision-logic error (negative lower bounds never reject;
infeasibility requires a feasible-witness upper bound U < −τ recovered
from the corrected multiplier, projected into the ball and re-evaluated).
Includes the error-exposing test battery, the pre-fixed smallcheck, the
2160-user batch with two-sided bound assertions, a v1→v2 decision-basis
comparison, and pipelines with explicitly named SDP+TR-completion
variants. See `SDP基线审查报告_修正版.md` and
`english_sdp_baseline_v2.md`. The v1 entry `python run_sdp.py`
(results_v41_sdp/) is retained as history with correction headers on its
reports.

**S-procedure-SDP baseline (v1, superseded)**: `python run_sdp.py` writes
into `results_v41_sdp/` — superseded by v2 above; its decision logic and
timing claims must not be cited.

**Optimization round**: `python run_optimization.py` (or `--quick`) writes
into `results_v41_optimization/`: joint-certification complementarity
analysis (max of Q_G and L_quad), staged screening cascades with fair paired
timing (all references tolerance-stopped), adaptive interval refinement
(64→256→1024→4096) vs fixed resolutions, a frozen strategy
(`frozen_strategy_config.json`, hashes + unused confirmation seeds 221–225),
and small/medium/large scale configurations. See `优化报告_v41.md` (Chinese)
and `english_optimization_v41.md` (paper-ready). Exits nonzero on any failed
check.

**V4.1 unified suite**: `python run_v41.py` (or `--quick` smoke mode) writes
everything into `results_v41/`: V4.1 Theorem-2/3 construction checks
(`theory_constructions.csv`, `boundary_instances.csv`, `boundary_family.csv`),
degenerate/endpoints battery, spectral scan diagnostics, the three candidate
pools (easy baseline + exploration grid + **pre-registered validation seeds
211–215** with `validation_config.json` written before the run), the
J = 256/1024/4096 resolution ablation with loss decomposition, overlap
accounting, and paired timing including the validated tolerance-stopped
trust-region reference. The suite exits nonzero on any failed check. See
`实验审查报告_v41.md` (Chinese review) and `english_sections_v41.md`
(paper-ready paragraphs and tables). Original `results/` and `results_review/`
are never overwritten.

Runs experiment A (checks), B (spectral scan), C (RIS candidate pool),
writes all CSV/JSON to `results/`, figures (PDF+PNG) to `figures/`, and
prints a post-run verification checklist.  All seeds and parameters are
fixed; nothing is tuned against results.  Individual experiments can be run
via `python experiments_a.py` / `experiments_b.py` / `experiments_c.py`.

## Files

| file | content |
|---|---|
| `core.py` | spectral summary (six scalars via Hermitian eigendecomposition), envelope `R`, bounds `D`/`Q`/`Q_G`/`L_quad`, certified interval lower bounds, exact sphere-constrained trust-region reference (PSD interior, singular PSD, negative curvature, hard case), A2/A3 constructions |
| `experiments_a.py` | A1 direction extrema, A2 three-branch construction, A3 certification-flip instance + deterministic flip, A4 boundary family |
| `experiments_b.py` | B: 21-point beta0 scan with per-point explicit matrices |
| `experiments_c.py` | C: 2 BS x 4 users, 16-element 2-bit RIS, 12 candidates x 5 seeds x 3 radii (720 user-level ball problems, dim 8) |
| `figures.py` | the two figures |
| `results/` | CSV/JSON outputs (see below) |
| `findings.md` | actual results, anomalies, supported/unsupported conclusions |

Outputs: `direction_checks.csv`, `construction_checks.csv`,
`boundary_checks.csv`, `a3_flip.json`, `spectral_scan.csv`,
`pool_user_details.csv` (720 rows), `pool_summary.csv`, `timings.csv`,
`config.json`; figures `fig1_direction_construction`, `fig2_spectral_scan`.

## Metric definitions (experiment C)

- user-level status from the exact reference: `clearly_feasible`
  (q\* > tau), `clearly_infeasible` (q\* < −tau), `ambiguous`
  (|q\*| <= tau), `reference_invalid` (KKT checks failed).
- a method *directly accepts* a candidate-radius unit iff its lower bound
  exceeds tau for **all four** users; otherwise the unit is *unknown* for
  that method.
- *coverage* = accepted AND reference-clearly-feasible units / units whose
  reference status is clearly feasible (NA if the denominator is 0).
- *misaccept* = accepted while at least one user is clearly infeasible.
  Any occurrence suspends performance conclusions.
- *fallback*: serial per-user exact reference with early stop at the first
  clearly-infeasible user; counted separately from full-reference
  evaluation.

## Numerical-safety limitations

- Double precision only; **not** a machine-checked formal certificate and
  no directed-rounding claim is made.
- Certification uses only interval lower bounds (default J = 256 cells on
  theta in [0, pi/2]; escalation to 1024/4096 if a gap check requires);
  grid/golden-section minima are labelled *approximations* and are never
  used for acceptance.
- Engineering threshold tau = 1e-8 · max(1, |c|, 2r‖g‖, r²‖A‖₂); accept
  only if LB > tau; near-zero outcomes are `ambiguous`.
- Experiment C uses dimensionless normalized links (no path-loss or
  hardware model); r_u = rho·‖h_u‖ is a specified effective-channel ball,
  with no claim of calibrated statistical coverage.  The 720 user-level
  problems are correlated (60 candidates x 3 radii x 4 users), not 720
  independent experiments.
- Timing: warm-up then 3 repetitions with `perf_counter`; medians with
  ranges; no small-sample significance tests.  End-to-end times are
  measured directly, never summed from parts.  Reference solves reuse one
  eigendecomposition per instance across radii (fair caching).
