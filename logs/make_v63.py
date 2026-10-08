# -*- coding: utf-8 -*-
"""Patch paper_v6_2.md -> paper/paper_v6_3.md (V6.3 changes)."""
import re

src = open("paper/paper_v6_2.md", encoding="utf-8").read()
n = [0]


def rep(old, new, cnt=1):
    assert src.count(old) == cnt, (src.count(old), old[:70])
    src_ = src
    globals()['src'] = src.replace(old, new)
    n[0] += 1


# ---- 1. Theorem 2 attainment: remove blanket divergence claim ----------
rep(
    "For the remaining nondegenerate cases ($α>0$, $r>0$), the "
    "constructed $ε$-family approaches the infimum along nominal "
    "channels whose norms grow without bound as $ε→0$; this is a "
    "property of that family, not an impossibility claim — whether "
    "some other instance of finite nominal-channel norm attains the "
    "infimum is left open.",
    "For the remaining nondegenerate cases ($α>0$, $r>0$), the "
    "construction does not furnish a uniform bound on the "
    "nominal-channel norm over the information class; divergence may "
    "occur in particular parameter families (Section 5 exhibits one), "
    "but is not asserted for every nondegenerate summary vector, and "
    "finite attainment is proved below for at least one interior "
    "family. Whether every summary vector admits a finite-norm "
    "attaining instance is open.")

# ---- 2. bounded endpoint remark + attainment remark after Thm 2 proof --
rep(
    "**Zero alpha and zero radius.**",
    "**Remark (bounded construction at the upper spectral endpoint).** "
    "At $β_0=√{α^2+β^2}$ we have $δ=β$ and $ν=α$, so the projection "
    "gives $k=α$ for every $(x,y)$: $d=αe_1$ and $ℓ=δe_1$ are "
    "collinear, $t_ε=1$, and $u_ε=−(αx+(1−ε)δy)e_1$ points along $e_1$ "
    "for every $ε$. The right-hand side of (20) then lies in "
    "$\\mathrm{span}\\{e_1\\}$, and the minimum-norm solution has "
    "$ε$-independent norm; the construction family therefore keeps "
    "$∥ĥ_ε∥$ bounded. Numerically, for $p=2$, $a=I_0=α=β=γ=1$, "
    "$β_0=√2$, $σ^2=0.05$, $r=0.2$ (minimizer $x^*≈0.18332$, "
    "$y^*≈0.07995$), $∥ĥ_ε∥=√5$ for every tested $ε$, while "
    "$q_ε(e_ε)−Q_G=O(ε)$ as guaranteed by (22). At this endpoint the "
    "infimum is in fact attained finitely: with "
    "$C^H=[[1,1,0,0,0],[0,0,0,1,0],[0,0,0,0,1]]$, $b=√2e_1$ and "
    "$ĥ=(1,−2,0,0,0)^T$, the six summaries are exactly "
    "$(2,1,1,1,1,√2)$ and the error $e=(−x^*,−y^*,0,0,0)^T$ satisfies "
    "$C^He=−(x^*+y^*)e_1$ parallel to $C^Hĥ=−e_1$, so "
    "$q(e)=2(1−x^*)^2−(1+x^*+y^*)^2−0.05=Q_G=q^*$ (verified to "
    "$7×10^{−16}$). We do not claim finite attainment for the whole "
    "upper-endpoint family.\n\n"
    "**Zero alpha and zero radius.**")

# ---- 3. Lemma 1 space-confusion fix ------------------------------------
rep(
    "Necessity: (3) is classical, and the two amplitude "
    "conditions are forced. If $β=0$, then $C^HP=0$, so the range of "
    "$C^H$ lies in span{v} and $C^Hx=(C^Hv)(v^Hx)$ for every $x$; hence ",
    "Necessity: (3) is classical, and the two amplitude "
    "conditions are forced. Throughout, $v∈ℂ^n$ and $C∈ℂ^{n×m}$, so "
    "$C^H:ℂ^n→ℂ^m$ maps input space to output space. If $β=0$, then "
    "$C^HP=0$, so $C^H=(C^Hv)v^H$ (equivalently, the range of $C$ "
    "lies in span{v} in $ℂ^n$) and $C^Hx=(C^Hv)(v^Hx)$ for every "
    "$x∈ℂ^n$; hence ")

rep(
    "Fix $n≥5$, $m≥3$, $p>0$ and a summary vector "
    "$s=(p,a,I_0,α,β,β_0)$ with $a,I_0≥0$.",
    "Fix $n≥5$, $m≥3$, $p>0$ and a summary vector "
    "$s=(p,a,I_0,α,β,β_0)$ with $a,I_0,α,β,β_0≥0$.")

# ---- 4. Lemma 2 norm/eigenvalue statement fix ---------------------------
rep(
    "and $∥C∥_2=∥[d,ℓ,ℓ_2]∥_2$ equals the "
    "largest eigenvalue of the block diagonal Gram "
    "$[[α^2,δk],[δk,∥ℓ∥^2]]⊕[β^2]$, which is "
    "$max{β_0^2,β^2}=β_0^2$.",
    "and $∥C∥_2^2=λ_{max}(CC^H)$ equals the "
    "largest eigenvalue of the block diagonal Gram "
    "$[[α^2,δk],[δk,∥ℓ∥^2]]⊕[β^2]$ — the other eigenvalue on the "
    "$(d,ℓ)$ block is $α^2−k^2≤α^2≤β_0^2$ since $k∈[ν,α]$ — so "
    "$∥C∥_2^2=max{β_0^2,β^2}=β_0^2$, i.e. $∥C∥_2=β_0$.")

rep(
    "A direct "
    "computation gives $∥d∥=α$, $d^Hℓ=δk$, "
    "$∥ℓ∥^2=β_0^2−k^2$, and the Gram matrix of $[d,ℓ]$ has "
    "$β_0^2$ as an eigenvalue — its characteristic polynomial "
    "$(α^2−λ)(β_0^2−k^2−λ)−δ^2k^2$ vanishes at $λ=β_0^2$ because "
    "$α^2+δ^2=β_0^2$ — so $∥[d,ℓ]∥_2=β_0$ and",
    "A direct "
    "computation gives $∥d∥=α$, $d^Hℓ=δk$, "
    "$∥ℓ∥^2=β_0^2−k^2$, and the Gram matrix $[d,ℓ]^H[d,ℓ]$ has "
    "$β_0^2$ as an eigenvalue — its characteristic polynomial "
    "$(α^2−λ)(β_0^2−k^2−λ)−δ^2k^2$ vanishes at $λ=β_0^2$ because "
    "$α^2+δ^2=β_0^2$, and the other eigenvalue $α^2−k^2$ does not "
    "exceed it since $k∈[ν,α]$ — so $∥[d,ℓ]∥_2^2=β_0^2$ and")

# ---- 5. Lemma 2 beta=0 rectangular construction ------------------------
rep(
    "If $β=0$ (so $β_0=α$), take "
    "$C^H=αvv^H$, for which $∥C^H(xv+yu)∥=αx=R(x,y)$; the zero "
    "operator is trivial. Zero padding reaches any $n≥3,m≥3$. □",
    "If $β=0$ (so $β_0=α$), take $C^H=αzv^H$ with any unit "
    "$z∈ℂ^m$: then $∥C^Hv∥=α$, $C^HP=0$ and $∥C∥_2=α$, and "
    "$∥C^H(xv+yu)∥=αx=R(x,y)$ for every $m,n$; the zero operator "
    "is trivial. Zero padding reaches any $n≥3,m≥3$. □")

# ---- 6. Lipschitz / supremum reference now to Lemma 2 -------------------
rep(
    "the same expression yields the Lipschitz control used in Section 6.",
    "the same expression, together with the attainment construction of "
    "Lemma 2, yields the Lipschitz control used in Section 6: the "
    "envelope is a supremum of realized linear output norms, each "
    "Lipschitz with constants $(α,β)$.")

# ---- 7. interface experiment subsection in 7.5 --------------------------
rep(
    "A controlled scaling study isolates dimension from load:",
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
    "two additional certified units in this pool.\n\n"
    "A controlled scaling study isolates dimension from load:")

# ---- 8. SDP timing wording (numbers unchanged, semantics refined) -------
rep(
    "On the same batch, the generic SDP "
    "pipeline is 17.3 times slower than the tolerance-stopped "
    "reference as a ratio of medians (median per-pair ratio 15.6; the "
    "reference is faster in 15 of 15 pairs), and the staged "
    "$L_quad→Q_G→$SDP pipeline is 6.5 times slower (median per-pair "
    "ratio 7.2); these are single-round figures from the corrected SDP "
    "study.",
    "On the same batch, the generic SDP "
    "pipeline is 17.3 times slower than the tolerance-stopped "
    "reference as a ratio of medians (median per-pair slowdown 15.7, "
    "range 13.2–18.8; the reference is faster in 15 of 15 pairs), and "
    "the staged $L_quad→Q_G→$SDP pipeline is 6.5 times slower (median "
    "per-pair slowdown 7.2, range 1.5–12.4); these are single-round "
    "figures from the corrected SDP study, with all ratios defined as "
    "pipeline time divided by reference time.")

# ---- 9. Discussion: divergence + kappa tradeoff concentration -----------
rep(
    "The sharpness statement is intentionally strong and narrow: it "
    "holds for the unrestricted six-scalar information class at fixed "
    "$r,γ,σ^2$ and dimensions, whose extremal constructions may "
    "require large nominal channel norms that a physical RIS channel "
    "may not admit.",
    "The sharpness statement is intentionally strong and narrow: it "
    "holds for the unrestricted six-scalar information class at fixed "
    "$r,γ,σ^2$ and dimensions; the approximating construction gives "
    "no uniform nominal-channel-norm bound over the class, divergence "
    "occurs in particular families, and bounded or attaining "
    "constructions exist in others (Remark after Theorem 2).")

rep("Version 6.2 • English revision • 7 October 2026",
    "Version 6.3 • English revision • 7 October 2026")

open("paper/paper_v6_3.md", "w", encoding="utf-8", newline="").write(src)
print(f"paper_v6_3.md written ({n[0]} patches)")
