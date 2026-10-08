# Paper-ready English sections: corrected S-procedure-SDP baseline (v2)

From `python run_sdp_v2.py` (results_v41_sdp_corrected/, exit code 0, all
registered checks passed). Replaces the v1 text in
english_sdp_baseline.md (kept with a correction header). CLARABEL via
CVXPY 1.7.5; identical shared error ball, thresholds and caching.

## Decision rule paragraph (replaces the v1 safety paragraph)

> The S-procedure SDP provides a certified lower bound $L$ on $q^*(r)$
> after residual correction, but a negative lower bound proves nothing
> about infeasibility. Infeasibility is certified by an upper bound: from
> the corrected multiplier $\lambda'$ we recover the Lagrangian minimizer
> $(A+\lambda'I)e=-g$ (least squares if singular), project it onto the
> ball if $\|e\|>r$ (required in 956 of 2160 batch cases), and re-evaluate
> the true objective $U=q(e_w)$. Since $e_w$ is ball-feasible by
> construction, $U$ is a valid upper bound of $q^*$ and $U<-\tau$
> certifies infeasibility; no optimality of $e_w$ is claimed. All
> remaining cases return *unresolved* and are never treated as
> rejections; the raw SDP value $t$ is never used as an upper bound, and
> an "optimal" solver status alone certifies nothing. The residual
> correction itself is an implementation detail — the shift
> $(t,\lambda)\mapsto(t-\eta(1+r^2),\lambda+\eta)$ adds
> $\mathrm{blockdiag}(\eta I,\eta)$ to the LMI (Weyl), with an empirical
> rounding allowance and a numerical re-check of positive semidefiniteness
> (achieved strictly in 2160/2160 cases; $\eta$ at most
> $6.5\times10^{-9}$); it is not a formal error bound and no
> interval-arithmetic certificate is claimed. An explicitly named
> "SDP plus trust-region completion" pipeline exists for users left
> unresolved by the SDP, with those calls counted separately; it was not
> triggered in this batch.

## Consistency table (validation batch, 2160 user problems)

| quantity | value |
|---|---:|
| $L \le q^*$ (lower-bound safety), max excess | −1.2e-13 (never above) |
| $U \ge q^*$ (upper-bound validity), max deficit | 1.3e-12 (fp noise) |
| bound ordering $L \le U$, violations | 0 |
| false acceptance / false rejection | 0 / 0 |
| witnesses inside the ball | 2160/2160 (956 projected) |
| decisions | 1247 feasible, 913 infeasible, 0 unresolved |
| witness quality $U-q^*$ over the 913 infeasible | median 6.7e-6, max 1.6e-4 |
| max $\|t-q^*\|$ | 5.6e-8 |
| per solve: build / solver / wall | 0.33 / 1.6 / 5.7 ms |

> A targeted test battery covers the failure mode of the earlier
> implementation: a valid negative lower bound on a truly feasible problem
> returns unresolved instead of a rejection; a feasible witness with
> $U<-\tau$ certifies infeasibility; straddling bounds return unresolved;
> out-of-ball witnesses are projected and re-evaluated; non-finite,
> non-Hermitian and negative-radius inputs are rejected; and decisions are
> invariant under overall problem rescaling by $10^{-3}$–$10^{3}$ (the
> threshold is scale-relative). In the validation batch, the 913 instances
> that the earlier implementation classified as infeasible from a negative
> lower bound are all re-confirmed by witness upper bounds — the same
> outcomes, now with a valid justification.

## Timing table (same final decision task; this run's numbers)

| pipeline | median s per (seed, γ) | SDP calls | TR completion calls |
|---|---:|---:|---:|
| trust-region reference (tolerance stopped) | 0.042 | — | — |
| $L_{\rm quad}\rightarrow Q_G\rightarrow$ SDP | 0.268 | 435 | 0 |
| old $Q\rightarrow$ SDP | 0.284 | 511 | 0 |
| $Q_G\rightarrow$ SDP | 0.302 | 496 | 0 |
| pure SDP (unresolved allowed) | 0.719 | 1319 | 0 |
| SDP + TR completion | 0.730 | 1319 | 0 |

> Screening with SDP fallback saves about 63% of the time and 67% of the
> solver calls of pure SDP certification on this batch. No unresolved
> units occurred, so the SDP-only pipelines happened to be complete
> decision procedures here; where unresolved units occur, only the
> explicitly named completion pipeline gives an unconditional end-to-end
> comparison. The trust-region reference remains faster than every SDP
> pipeline in all paired comparisons; the slowness of the generic conic
> path (CVXPY canonicalization plus interior-point iterations) is not a
> statement about full-information methods — the trust-region reference is
> itself full-information and is the fastest method overall.

## Attribution sentence (corrected)

> The single-sphere S-procedure is a classical technique, used here as an
> adaptation representative of robust RIS design [Zheng et al., IEEE Trans.
> Commun. 73(8):6306–6321, 2025]; the residual shift and the
> witness-recovery procedure are implementation details of this project
> that make interior-point output usable as acceptance certificates;
> agreement between the two independent exact methods is numerical
> evidence, not an analytical proof.
