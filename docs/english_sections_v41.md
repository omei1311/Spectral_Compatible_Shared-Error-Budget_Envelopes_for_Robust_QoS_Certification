# Paper-ready English sections and tables for V4.1 (from results_v41/)

All numbers below are produced by `python run_v41.py` (results_v41/, full
mode, manifest failures empty). Insert or replace in the indicated places.

## Replacement for the multiplier-search sentence in Section 6.2

> The reference solves the Hermitian trust region problem by eigenvalue
> decomposition and a multiplier search, handling PSD interior solutions,
> singular cases and the indefinite hard case. The global optimality
> conditions are stationarity, ball feasibility, $A+\mu I\succeq 0$,
> $\mu\geq 0$ and complementary slackness [2,3]. Numerical checks inspect the
> full stationarity residual, shifted minimum eigenvalue, ball violation and
> complementarity. A 106 case audit reports a maximum primal–dual gap of
> $4.97\times10^{-14}$. The fixed 200 iteration multiplier bisection used in
> the original timing over-resolves the multiplier: stopping the bisection
> once the bracket width falls below $10^{-12}$ of its magnitude changes no
> computed value by more than $5.7\times10^{-10}$, preserves every validity
> check, and uses a median of 40 iterations. A $10^{-10}$ stop already shifts
> some values by up to $2.8\times10^{-8}$ and is not used. Timing sections
> therefore report both the original 200 iteration reference and the
> validated tolerance stopped reference.

## New subsection for Section 7 (after the exploratory grid): pre-registered validation batch

> **Pre-registered validation batch.** Before computing any result on new
> seeds, we wrote a fixed configuration file recording the unused seeds
> 211–215, the grid $\gamma\in\{0.5,1,2\}$, $\rho\in\{0.06,0.15,0.30\}$, the
> candidate generation rules, the acceptance rule (interval lower bound
> exceeding the threshold for all four users at $N_\theta=256$), the
> reference validity criteria, the timing protocol (one warm-up plus three
> interleaved repetitions, per user eigendecomposition caching, certified
> user skipping, infeasibility early stopping, candidate generation outside
> the timer) and SHA-256 hashes of the four code files. The batch adds 540
> candidate–radius–threshold units and 2160 user ball evaluations with zero
> overlap with the previous pools; the union of the three pools contains 6960
> distinct user problems. All references are valid and no false acceptance
> or bound violation occurs. Of the 540 units, 159 are clearly feasible; the
> implemented bounds directly accept D 111, previous 125, new 126 and
> $L_{\rm quad}$ 138 units. The new bound adds one unit over the previous
> bound, raising feasible coverage from 78.62% to 79.24% (+0.63 percentage
> points). Combined with the exploratory grid (+10 units, +2.12 points), the
> certification gain is positive but small and varies between batches; both
> batches are finite correlated samples and no distributional claim is made.
> Reference calls fall from 511 (previous bound pipeline) to 496 (new bound
> pipeline) against 1319 for the full reference baseline.

Validation table (per cell; 60 units each):

| γ | ρ | feasible | D | previous | new | L_quad |
|---|---|---:|---:|---:|---:|---:|
| 0.5 | 0.06 | 54 | 53 | 53 | 53 | 54 |
| 0.5 | 0.15 | 44 | 28 | 36 | 37 | 40 |
| 0.5 | 0.30 | 13 | 0 | 4 | 4 | 5 |
| 1 | 0.06 | 29 | 22 | 23 | 23 | 28 |
| 1 | 0.15 | 15 | 5 | 6 | 6 | 8 |
| 1 | 0.30 | 0 | 0 | 0 | 0 | 0 |
| 2 | 0.06 | 4 | 3 | 3 | 3 | 3 |
| 2 | 0.15 | 0 | 0 | 0 | 0 | 0 |
| 2 | 0.30 | 0 | 0 | 0 | 0 | 0 |
| total | | 159 | 111 | 125 | 126 | 138 |

## Replacement/addition for Section 7.5 (resolution and timing)

> **Resolution and loss decomposition.** Table R separates interval
> resolution effects. In the exploratory grid the previous envelope gains one
> additional accept at 4096 cells while the new envelope is unchanged, so the
> ten additional accepts are tied to the implemented 256 cell setting. In the
> validation batch the counts are unchanged at every resolution. Decomposing
> the gap on the 2160 validation rows, the numerical interval loss (one
> dimensional upper estimate minus the 256 cell certified bound) has median
> $1.22\times10^{-3}$ and maximum $2.53\times10^{-2}$, whereas the
> theoretical envelope gap (reference minus upper estimate) has median
> $6.08\times10^{-2}$ and maximum $2.68$: certification outcomes are limited
> by the theoretical envelope rather than by the interval discretization,
> which supports the 256 cell main comparison.

Table R — direct accepts by interval resolution:

| cohort | cells | previous | new | added | feasible |
|---|---:|---:|---:|---:|---:|
| exploration | 256 | 380 | 390 | 10 | 472 |
| exploration | 1024 | 380 | 390 | 10 | 472 |
| exploration | 4096 | 381 | 390 | 9 | 472 |
| validation | 256 | 125 | 126 | 1 | 159 |
| validation | 1024 | 125 | 126 | 1 | 159 |
| validation | 4096 | 125 | 126 | 1 | 159 |

> **End to end timing.** Timing measures the same final candidate decision
> task for all pipelines (summary construction, all radii, caching, validity
> checks, early stopping; candidate generation excluded; one warm-up plus
> three interleaved repetitions per seed and threshold; medians over five
> validation seeds and three thresholds). Table T reports the original 200
> iteration reference, the validated tolerance stopped reference, and the
> screening pipelines with both fallback variants. With the tolerance
> stopped reference the full baseline drops from 0.150 s to 0.041 s per
> seed and threshold, so the apparent advantage of screening pipelines over
> the exact baseline narrows substantially: the best configuration is the
> new bound with tolerance stopped fallback (0.032 s), a saving of about 20%
> relative to the tolerance stopped pure reference, while the screening
> pipeline with the original 200 iteration fallback (0.075 s) is slower than
> the tolerance stopped pure reference. Across seed–threshold pairs the new
> pipeline is faster than the previous pipeline in 6 of 15 validation pairs
> and in 11–13 of 30 exploratory pairs across two full runs (paired counts
> vary with machine timing noise), and reduces total reference calls by
> 3.41% (850→821) and 2.9% (511→496) respectively. The efficiency
> conclusions are therefore limited to the tested pool, scale and reference
> implementation; in particular, part of the previously reported saving is
> attributable to multiplier search over-iteration in the reference rather
> than to screening.

Table T — same decision task, medians over validation seeds (seconds;
parentheses: total reference calls over the batch):

| pipeline | 200-iteration reference | tolerance reference |
|---|---:|---:|
| pure reference | 0.1496 (1319) | 0.0411 (1319) |
| new bound + fallback | 0.0750 (496) | 0.0323 (496) |

## Suggested abstract/conclusion wording adjustments

- Replace "The resulting reduction in reference calls does not establish a
  stable runtime advantage over the previous envelope" with: "The reduction
  in reference calls (3.41% exploratory, 2.9% pre-registered validation)
  does not establish a stable runtime advantage over the previous envelope,
  and a validated tolerance stop of the reference multiplier search removes
  most of the remaining gap between screening and the exact baseline."
- Add to the abstract's numerical sentence: "a pre-registered validation
  batch on five unused seeds reproduces the direction of the certification
  gain at a smaller magnitude (+0.63 versus +2.12 percentage points)."

## Reproducibility addendum

> The unified entry point `python run_v41.py` regenerates theory checks,
> spectral scan, all three candidate pools, the resolution ablation, overlap
> accounting and timing into `results_v41/`; `python run_v41.py --quick`
> provides a small-scale smoke mode. Original `results/` and `results_review/`
> are preserved. The suite exits nonzero if any registered check fails.
