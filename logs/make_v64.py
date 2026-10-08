# -*- coding: utf-8 -*-
"""Patch paper_v6_3.md -> paper/paper_v6_4.md (V6.4 changes)."""
src = open("paper/paper_v6_3.md", encoding="utf-8").read()
n = [0]


def rep(old, new, cnt=1):
    assert src.count(old) == cnt, (src.count(old), old[:70])
    globals()['src'] = src.replace(old, new)
    n[0] += 1


# 1. 'interior family' -> endpoint instance (wording fix)
rep(
    "divergence may "
    "occur in particular parameter families (Section 5 exhibits one), "
    "but is not asserted for every nondegenerate summary vector, and "
    "finite attainment is proved below for at least one interior "
    "family.",
    "divergence may "
    "occur in particular parameter families (Section 5 exhibits one), "
    "but is not asserted for every nondegenerate summary vector, and "
    "finite attainment is proved below for at least one nondegenerate "
    "instance (at the upper spectral endpoint).")

# 2. interface study paragraph rewrite (strict clients, batch-measured)
rep(
    "A focused interface study asks why a user would retain the "
    "six-scalar summary rather than adding $κ$ or keeping the full "
    "geometry, in the cached-descriptor, repeated-query scenario "
    "(seeds 401–402, preregistered configuration). Per user problem, "
    "the six-scalar descriptor occupies 48 bytes and answers a radius "
    "query in about 47 µs (median), the $6+κ$ descriptor 56 bytes and "
    "about 3.5 µs, and the cached full-geometry branch (data plus one "
    "eigendecomposition) about 1.7 KB and 340–440 µs per query. The "
    "storage difference between six and seven scalars is negligible "
    "(8 bytes); the practical difference is between any fixed-size "
    "summary and the full geometry. Over the 240 units of the study, "
    "direct certification coverage is 119/240 for $Q_G$, 123/240 for "
    "$L_quad$ and 124/240 for the joint rule (denominator: all "
    "candidate–radius units), so the summary interface answers most "
    "queries without fallback, and the exact-decision task requires "
    "fallback for the remainder in all three branches. The study does "
    "not show a runtime advantage of six scalars over $6+κ$ — the "
    "quadratic certificate is both cheaper per query and slightly "
    "wider in coverage; the six-scalar envelope contributes one to "
    "two additional certified units in this pool.",
    "A focused interface study (seeds 401–402, preregistered "
    "configuration) implements strict query-side clients — the "
    "six-scalar client receives only $(p,a,I_0,α,β,β_0)$ plus the "
    "external radius and parameters, the seven-scalar client adds "
    "$κ$, and a compact exact client receives only the eigendecomposed "
    "spectrum $(λ,w,c)$ of $(A,g)$ — and measures actual batch totals "
    "for 1, 10 and 100 repeated radius queries (longer batches cycle "
    "the same deterministic grid). Numerical payloads are 48 and 56 "
    "float64 bytes for the summary clients, 144 bytes for the compact "
    "exact cache ($λ$, $w$, $c$ plus scale, actual array nbytes), and "
    "about 1.7 KB for the full-geometry cache. Amortized per query, "
    "the seven-scalar quadratic certificate is the cheapest branch "
    "(well under 1 µs), the six-scalar envelope client costs about "
    "0.9 µs, and the compact exact cache about 25 µs — so the strict "
    "seven-scalar client is faster than the six-scalar one, and the "
    "compact exact cache shrinks the storage gap between summaries "
    "and exact certificates to a factor of about three. Direct "
    "certification coverage over the 240 candidate–radius units is "
    "119 for $Q_G$, 123 for $L_quad$ and 124 for the joint rule "
    "(91.5%, 94.6% and 95.4% of the 130 reference-feasible units; the "
    "joint rule adds exactly one unit over $L_quad$ in this pool). "
    "With summary screening and exact fallback, 217 compact-cache "
    "queries complete all 240 exact decisions with no unresolved "
    "units. The study shows no runtime or storage advantage of the "
    "six-scalar interface over the seven-scalar one; the envelope's "
    "practical contribution here is coverage supplementation.")

# 3. Discussion: storage claim contraction
rep(
    "Empirically, the certificate gain over the separate-gain envelope "
    "is positive in the sensitivity, validation and heterogeneous-link "
    "batches and zero in the saturated easy batch; no stable runtime "
    "advantage over the separate-gain pipeline is observed, and "
    "measured speedups are tied to the tested implementation.",
    "Empirically, the certificate gain over the separate-gain envelope "
    "is positive in the sensitivity, validation and heterogeneous-link "
    "batches and zero in the saturated easy batch; no runtime "
    "advantage of the six-scalar interface over the seven-scalar one "
    "was observed in the strict-client study, and a compact exact "
    "cache of the spectrum reduces the storage gap to within a small "
    "factor, so the engineering case for the summary interface rests "
    "on fixed-size descriptors rather than on speed or large storage "
    "savings.")

rep("Version 6.3 • English revision • 7 October 2026",
    "Version 6.4 • English revision • 7 October 2026")

open("paper/paper_v6_4.md", "w", encoding="utf-8", newline="").write(src)
print(f"paper_v6_4.md written ({n[0]} patches)")
