# Paper-ready English sections: S-procedure-SDP baseline (results_v41_sdp/)

> **CORRECTION (2026-10-06).** The v1 decision logic classified a
> certified lower bound $t_{\rm cert}<-\tau$ as infeasible, which is
> invalid: a negative lower bound proves nothing about $q^*<0$. The
> corrected sections are in **english_sdp_baseline_v2.md** (results in
> results_v41_sdp_corrected/): infeasibility is now certified only by a
> feasible-witness upper bound $U<-\tau$. All 913 affected instances were
> re-confirmed by witnesses. Do not cite the v1 status-agreement or
> timing claims ("70% saving") from this file.

Numbers from `python run_sdp.py` (full mode; exit code 0, all registered
checks passed). CLARABEL via CVXPY 1.7.5; no warm start; identical shared
error ball, thresholds and caching across all pipelines.

## Baseline description (for the reference list / method section)

> As a full-information baseline independent of our eigendecomposition
> reference, we adapt the standard single-sphere S-procedure LMI to the
> fixed-candidate certification problem: for each user,
> maximize $t$ subject to $\lambda\geq0$ and
> $[[A+\lambda I, g],[g^H, c-t-\lambda r^2]]\succeq0$, whose optimum
> equals $q^*(r)$ by the S-lemma; $A=bb^H-\gamma CC^H$, $g=A\hat h$,
> $c=\hat h^HA\hat h-\gamma\sigma^2$ use the same shared error ball as all
> other methods (no separate signal/interference spheres). The LMI is
> solved on the real embedding of the Hermitian data with CLARABEL. This
> is a classical technique — representative of its use in RIS robust
> design [Zheng et al., IEEE Trans. Commun. 73(8):6306–6321, 2025; DOI
> 10.1109/TCOMM.2024.3516491] — adapted to our certification interface;
> it is not a reproduction of that paper's energy-efficiency algorithm and
> is not claimed as a new method.

## Numerical safety paragraph (Section on implementation)

> A solver status of "optimal" is not a safety certificate: in the
> validation batch, 933 of 2160 optimal solutions have a slightly negative
> minimum eigenvalue of the LMI (at worst $-6.5\times10^{-9}$). We
> therefore convert each solution into a dual-feasible certified lower
> bound by the shift
> $(t,\lambda)\mapsto(t-\eta(1+r^2),\,\lambda+\eta)$ with
> $\eta=\max(0,-\lambda_{\min}(M))+K(2n{+}1)\|M\|_2 2^{-52}$, $K=10$,
> which adds $\mathrm{blockdiag}(\eta I,\eta)$ to the LMI matrix so that
> Weyl's inequality restores positive semidefiniteness up to eigensolver
> rounding; the corrected minimum eigenvalue is re-checked numerically
> (one retry with doubled $\eta$ if needed). Corrections are negligible
> (median $3.6\times10^{-13}$, maximum $8.9\times10^{-9}$). The corrected
> bound is used for all acceptance decisions; the raw optimal value $t$ is
> a numerical estimate only, and no interval-arithmetic formal certificate
> is claimed. Across all 2160 user problems the corrected bound never
> exceeded the trust-region reference (maximum excess
> $-1.2\times10^{-13}$), the maximum absolute deviation of $t$ from the
> reference was $5.6\times10^{-8}$, and per-user feasibility statuses
> agreed in 100% of cases — two independent exact methods (eigenvalue
> multiplier search and conic SDP) thus cross-validate each other and the
> certification pipeline.

Table SDP-1 — validation batch (2160 user problems, 540 candidate units):

| quantity | value |
|---|---:|
| solved / valid / failed / timeout / unresolved | 2160 / 2160 / 0 / 0 / 0 |
| max abs and rel deviation of $t$ from $q^*$ | 5.6e-8 |
| corrected bound minus $q^*$, maximum | −1.2e-13 (never above) |
| status agreement with trust-region reference | 2160/2160 |
| per solve: model build / solver / wall | 0.7 / 2.7 / 10.7 ms |
| LMI $\lambda_{\min}<0$ at solver output | 933/2160 (correction applied) |

## Timing paragraph (Section on cost; same final decision task)

> Timing compares five pipelines on the same candidates, radii, thresholds
> and final decision target, with identical caching, certified-user
> skipping, infeasibility early stopping, one warm-up and three interleaved
> repetitions per seed and threshold, and no warm start anywhere. The
> tolerance-stopped trust-region reference is the fastest pure method
> (median 0.088 s per seed and threshold); the generic conic SDP pure
> certification is about 12.5× slower (1.18 s), with the gap dominated by
> CVXPY canonicalization and interior-point iterations rather than problem
> information — the trust-region reference is itself a full-information
> method. Screening restores most of the SDP cost (the
> $L_{\rm quad}\rightarrow Q_G\rightarrow$ SDP pipeline runs at 0.358 s,
> a 70% saving over pure SDP with 435 instead of 1319 solver calls, and
> the cheap-bound-first ordering remains best), but every SDP-based
> pipeline remains 4–13× slower than the trust-region reference, which is
> faster in all 15 paired comparisons. All pipelines reproduce the
> reference decisions on all 540 units. The S-procedure-SDP baseline is
> therefore valuable as an independent cross-check of correctness, not as
> a runtime alternative in this problem structure.

Table SDP-2 — same decision task, medians over 15 (seed, γ) pairs:

| pipeline | median s | solver calls | decision agreement vs TR |
|---|---:|---:|---:|
| trust-region reference (tolerance stopped) | 0.088 | 1319 | — |
| $L_{\rm quad}\rightarrow Q_G\rightarrow$ SDP | 0.358 | 435 | 540/540 |
| old $Q\rightarrow$ SDP | 0.460 | 511 | 540/540 |
| $Q_G\rightarrow$ SDP | 0.478 | 496 | 540/540 |
| pure S-procedure SDP | 1.180 | 1319 | 540/540 |

## Attribution sentence (for the discussion)

> The single-sphere S-procedure LMI is a classical technique adapted from
> the robust RIS literature [Zheng et al., 2025]; the residual-corrected
> dual-feasible lower bound construction (Weyl shift with rounding
> allowance and numerical re-check) is our numerical-safety contribution,
> and it is what allows an interior-point SDP output to be used as an
> acceptance certificate at all.
