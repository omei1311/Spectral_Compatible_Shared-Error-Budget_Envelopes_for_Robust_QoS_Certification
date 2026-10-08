# Paper-ready English sections: joint certification, staged screening and cost (results_v41_optimization/)

All numbers come from `python run_optimization.py` (full mode; manifest
failures empty). Every exact reference in this round uses the validated
tolerance-stopped multiplier search (value deviation at most
5.7×10⁻¹⁰ from the 200-iteration version, all validity checks preserved).

## Insert after the sensitivity/validation subsections (Section 7)

> **Joint certification and staged screening.** The quadratic bound
> $L_{\rm quad}$ uses full quadratic geometry, while the new envelope bound
> uses only the six summaries; neither dominates. Taking the per-user
> maximum $L_{\rm joint}=\max\{Q_{G,\rm LB},L_{\rm quad}\}$ and accepting a
> candidate only when every user passes the common threshold yields a
> strictly larger direct-acceptance set: 424 of 472 feasible units (89.83%)
> in the exploratory grid versus 408 (86.44%) for $L_{\rm quad}$ alone and
> 390 (82.63%) for the envelope alone, and 139 of 159 (87.42%) versus 138
> and 126 in the pre-registered validation batch. The gain over
> $L_{\rm quad}$ comes entirely from units certified only by the six scalar
> envelope (16 and 1 units respectively; 113 and 52 user-level problems are
> certified only by the envelope), so it is a composition benefit and is not
> attributed to the envelope alone. No unit is accepted by the joint rule
> while both individual rules fail. The separable bound is dominated by the
> envelope at every user problem in all pools and is dropped from the flow.
>
> We then time staged screening pipelines that all complete the same final
> decision task: bounds are evaluated only on users not yet certified, a
> failed bound is never a rejection, the exact reference (tolerance stopped,
> validity checked, invalid references left unresolved) ends a candidate at
> the first clearly infeasible user, and each user's eigendecomposition is
> cached across radii and stages. With one warm-up and three interleaved
> repetitions per seed and threshold, the ordering $L_{\rm quad}\rightarrow
> Q_G\rightarrow$ reference — selected by a fixed rule on the exploratory
> cohort and frozen, with code hashes and unused seeds recorded before any
> new computation — attains a median paired speedup of 2.12 over the pure
> tolerance-stopped reference (25 of 30 pairs; ratio of medians 2.10) on
> the exploratory grid and 2.01 (11 of 15; ratio of medians 2.02) on five
> fresh confirmation seeds (221–225), while reducing total reference calls
> by 75% (2946→731) and 72% (1424→394) respectively. Screening with
> $L_{\rm quad}$ alone is the fastest single-bound pipeline (2.65× and
> 2.52×); the envelope stage contributes the additional certified units at
> about 7% further reference-call savings, not speed. All pipelines produce
> identical final decisions in every run.

Table J — direct certification and cost (units; reference calls):

| batch | feasible | $Q_G$ | $L_{\rm quad}$ | joint | selected pipeline speedup | calls (pure→selected) |
|---|---:|---:|---:|---:|---:|---:|
| exploratory (1080 units) | 472 | 390 | 408 | 424 | 2.12× (25/30) | 2946→731 |
| validation (540 units) | 159 | 126 | 138 | 139 | — | — |
| confirmation (540 units, fresh seeds) | 202* | 471 calls | 422 calls | 394 calls | 2.01× (11/15) | 1424→394 |

*202 is the common final accepted set of every pipeline (identical decisions);
the direct-acceptance split is the reference-call column.

> **Adaptive interval refinement.** Starting at 64 cells and refining
> through 256, 1024 and 4096 only while the certified lower bound remains
> below the threshold and a rounding-guarded feasible-point value of the
> one dimensional objective stays above it, reproduces the acceptance of
> fixed 4096 cells exactly while evaluating 3.9× fewer interval cells
> (176 768 versus 686 592 at fixed 256 over the exploratory grid). At this
> problem scale the additional Python-level bookkeeping outweighs the
> arithmetic saving, so the adaptive variant is slower in wall clock
> (1.45× versus 1.69× speedup over the pure reference) and is not part of
> the selected pipeline; it remains useful where interval arithmetic
> dominates, for example at much finer resolutions or in compiled
> implementations.

> **Scale behaviour.** With the candidate generation rules held fixed (per
> base station total power 1, per stream power 1/users per cell, identical
> normalization), three configurations — small (4 stacked dimensions, 4
> users), medium (8, 4) and large (16, 6) — show that the screening
> advantage over the tolerance-stopped pure reference grows with dimension:
> 1.1× at small, 4.1× at medium and 3.5× at large (median per seed; 144
> units per configuration; three seeds). Reference cost grows with the
> problem dimension while bound evaluation is nearly flat. Per-stream power
> differs between the large configuration (1/3) and the others (1/2), so
> acceptance rates are not comparable across configurations; only within-
> configuration pipeline differences are compared. All pipelines again
> produce identical final decisions.

Table S — scale configurations (γ=1.0, ρ=0.15; medians per seed):

| config | dim | users | pure ref (s) | $L_{\rm quad}$+ref (s) | $Q_G$+ref (s) | calls pure/$Q_G$/$L_{\rm quad}$ |
|---|---:|---:|---:|---:|---:|---|
| S | 4 | 4 | 0.0071 | 0.0066 | 0.0098 | 57/38/36 |
| M | 8 | 4 | 0.0167 | 0.0041 | 0.0102 | 104/36/25 |
| L | 16 | 6 | 0.0229 | 0.0065 | 0.0148 | 148/38/25 |

## Suggested wording for the discussion/limitations

> The optimization round adds engineering composition benefits (the joint
> maximum of the two bounds and the $L_{\rm quad}$-first staging) and
> implementation benefits (tolerance-stopped reference, shared caching).
> These are distinct from the theoretical contributions: the six scalar
> envelope remains the component that certifies units beyond the quadratic
> bound at equal information, while the quadratic bound, which uses more
> information, remains the cheaper first stage. No complexity order claim
> is made; all timings are paired, interleaved and limited to the tested
> pools, scales and Python implementation.
