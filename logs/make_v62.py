# -*- coding: utf-8 -*-
"""Patch paper_v6_1.md -> paper/paper_v6_2.md with the V6.2 changes."""
import re

src = open("paper/paper_v6_1.md", encoding="utf-8").read()
n_patch = 0


def rep(old, new, count=1):
    global src, n_patch
    assert src.count(old) == count, (src.count(old), old[:80])
    src = src.replace(old, new)
    n_patch += 1


# ---- 1. Lemma 1 proof: full A-E branch construction --------------------
rep(
    "**Proof.** Necessity: if $β=0$, then $C^HP=0$, so $C^Hx=(C^Hv)(v^Hx)$ "
    "for every $x$, hence $I_0=∥C^Hv∥|v^Hĥ|=αa$. If $β_0=α>0$, then "
    "$CC^Hv=α^2v$ (norm attainment at $v$), so "
    "$I_0^2=∥C^Hĥ∥^2=α^2a^2+∥(I−d̂d̂^H)C^Hĥ∥^2≥α^2a^2$ with $d̂=C^Hv/α$. "
    "Sufficiency: we exhibit instances with $v=e_1$ and $b=√p e_1$. "
    "(a) If $α=0$, then $β_0=β$; for $β>0$ take $C^H=[0,β]$ and "
    "$ĥ=(a,I_0/β)$, and for $β=0$ (so $I_0=0$ by (i)) take $C=0$ and "
    "$ĥ=(a,0)$. (b) If $α>0$ and $β_0>α$, set $δ=√{β_0^2−α^2}>0$, so "
    "$β≥δ>0$ by (3), and take $C^H=mat{α&δ&0;0&0&β}$, for which "
    "$∥C^Hv∥=α$, $∥C∥_2=β_0$ and $∥C^HP∥_2=max{δ,β}=β$ by the "
    "two-by-two determinant computation of Section 3; the nominal channel "
    "is $ĥ=(a,0,√{I_0^2−α^2a^2}/β)$ if $I_0≥αa$, and "
    "$ĥ=(a,(I_0−αa)/δ,0)$ if $I_0<αa$, and in both cases $∥C^Hĥ∥=I_0$ "
    "and $|v^Hĥ|=a$. (c) On the boundary $β_0=α>0$ with $I_0≥αa$ "
    "(condition (ii)), take $C^H=diag{α,β}$ and $ĥ=(a,J_⊥/β)$ with "
    "$J_⊥=√{I_0^2−α^2a^2}$ if $β>0$; if $β=0$, then (i) forces "
    "$I_0=αa$, and $C^H=diag{α,0}$ with $ĥ=(a,0)$ realizes the "
    "summaries. Zero padding of inputs and outputs reaches any fixed "
    "$n≥5,m≥3$; the padded error coordinates enter neither the desired "
    "nor the interference term, so the infimum in (1) is unchanged. □",
    "**Proof.** Necessity: (3) is classical, and the two amplitude "
    "conditions are forced. If $β=0$, then $C^HP=0$, so the range of "
    "$C^H$ lies in span{v} and $C^Hx=(C^Hv)(v^Hx)$ for every $x$; hence "
    "$I_0=∥C^Hv∥|v^Hĥ|=αa$, and in particular $α=β=β_0=0$ forces "
    "$I_0=0$. If $β_0=α>0$, norm attainment at $v$ gives $CC^Hv=α^2v$, "
    "and splitting $ĥ$ along $v$ and $v^⊥$ yields "
    "$I_0^2=α^2a^2+∥(I−d̂d̂^H)C^Hĥ∥^2≥α^2a^2$ with $d̂=C^Hv/α$. "
    "Sufficiency: all instances below use $v=e_1$ and $b=√p e_1$, "
    "followed by zero padding to the fixed $(n,m)$; the padded error "
    "coordinates enter neither the desired nor the interference term, so "
    "the infimum in (1) is unchanged. The branches are mutually exclusive "
    "and exhaustive: **(A)** zero operator ($α=β=β_0=0$, hence $I_0=0$ "
    "by (i)): $C=0$, $ĥ=ae_1$. **(B)** $β=0<α$ (so $β_0=α$ and "
    "$I_0=αa$): $C^H=[α,0,…]$, $ĥ=ae_1$. **(C)** $α=0<β$ (so "
    "$β_0=β$): $C^H=[0,β]$, $ĥ=(a,I_0/β)^T$. **(D)** $0<α<β_0$: with "
    "$δ=√{β_0^2−α^2}>0$ and $β≥δ$ by (3), take "
    "$C^H=mat{α&δ&0;0&0&β}$ and $ĥ=(a,(I_0−αa)/δ,0)^T$, for which "
    "$C^Hĥ=(I_0,0)^T$, $∥C^Hv∥=α$, $∥C^HP∥_2=max{δ,β}=β$ and "
    "$∥C∥_2=max{√{α^2+δ^2},β}=β_0$. **(E)** $β_0=α>0$, $β>0$ (so "
    "$I_0≥αa$): $C^H=diag{α,β}$, $ĥ=(a,J_⊥/β)^T$ with "
    "$J_⊥=√{I_0^2−α^2a^2}$, giving $I_0^2=α^2a^2+J_⊥^2$. Every branch "
    "realizes the six summaries exactly; (D) covers all $I_0≥0$ with a "
    "single formula. □")

# ---- 2. drop the unproven minimality claim in 2.3 -----------------------
rep(
    "The six-scalar interface is studied because it is the smallest "
    "summary set that already exhibits nontrivial spectral compatibility, "
    "and because the answer is closed-form.",
    "The six-scalar interface is studied because, together with the "
    "nominal amplitudes required by the QoS envelope, it retains the "
    "three spectral quantities of one operator, so that joint spectral "
    "compatibility and the strongest uniform guarantee can be analyzed "
    "in closed form at a fixed summary size.")

# ---- 3. Lemma 2 (direction attainment) replacing the 'two readings' ----
rep(
    "Two readings of (6) should be distinguished. For a fixed operator, "
    "$R$ is a safe upper bound on the interference output, possibly "
    "strict. Over the information class, $R$ is also the exact extremum: "
    "the constructions of Section 4 realize (6) with equality, so no "
    "smaller uniform bound is possible from the three summaries. The "
    "envelope is nondecreasing in either nonnegative input amplitude, as "
    "follows from (5); the same expression yields the Lipschitz control "
    "used in Section 6.",
    "For a fixed operator, $R$ is a safe upper bound on the "
    "interference output, possibly strict. The next lemma shows that "
    "over the class of operators with the same three spectral summaries "
    "the bound is exact, so no smaller uniform interference bound "
    "follows from $(α,β,β_0)$ alone. This is a statement about the "
    "direction problem only; the six-scalar QoS sharpness of Section 4 "
    "concerns the smaller class that also fixes $a$ and $I_0$, where "
    "the nominal channel enters. The envelope is nondecreasing in "
    "either nonnegative input amplitude, as follows from (5); the same "
    "expression yields the Lipschitz control used in Section 6.\n\n"
    "**Lemma 2 (attainment over the three-spectral class).** Fix "
    "compatible $(α,β,β_0)$ and $x,y≥0$. Let $𝒟$ be the class of "
    "operators $C$ on $n≥3$ inputs and $m≥3$ outputs with "
    "$∥C^Hv∥=α$, $∥C^HP∥_2=β$, $∥C∥_2=β_0$, let $u⊥v$ be unit, and "
    "let $e=xv+yu$. Then $max_{C∈𝒟}∥C^He∥=R(x,y)$.\n\n"
    "**Proof.** Theorem 1 gives $∥C^He∥≤R(x,y)$ for every $C∈𝒟$. For "
    "attainment, suppose first $α>0$. Let $k$ be the projection of "
    "$δx/y$ onto $[ν,α]$ (with $k=α$ if $y=0$) and set "
    "$d=(α,0,0)^T$, $ℓ=(δk/α,β_0√{1−k^2/α^2},0)^T$. A direct "
    "computation gives $∥d∥=α$, $d^Hℓ=δk$, $∥ℓ∥^2=β_0^2−k^2$, and "
    "the Gram matrix of $[d,ℓ]$ has $β_0^2$ as an eigenvalue — its "
    "characteristic polynomial $(α^2−λ)(β_0^2−k^2−λ)−δ^2k^2$ vanishes "
    "at $λ=β_0^2$ because $α^2+δ^2=β_0^2$ — so $∥[d,ℓ]∥_2=β_0$ and "
    "$∥xd+yℓ∥^2=α^2x^2+y^2(β_0^2−k^2)+2xyδk=R^2(x,y)$ by the branch "
    "evaluation in the proof of Theorem 1. Let $u'$ be a unit vector "
    "orthogonal to both $v$ and $u$, and set "
    "$C^H=[d,ℓ,ℓ_2][v,u,u']^H$ with $ℓ_2=βe_3$; this is allowed since "
    "$∥ℓ∥=√{β_0^2−k^2}≤β$ (as $k≥ν$) and $e_3$ is orthogonal to $d$ "
    "and $ℓ$. Then $C^Hv=d$ and $C^Hu=ℓ$, so $∥C^Hv∥=α$; "
    "$C^HP$ acts as $ℓ(u^H\\cdot)+ℓ_2(u'^H\\cdot)$ on $\\{v\\}^⊥$, so "
    "$∥C^HP∥_2=max{∥ℓ∥,β}=β$; and $∥C∥_2=∥[d,ℓ,ℓ_2]∥_2$ equals the "
    "largest eigenvalue of the block diagonal Gram "
    "$[[α^2,δk],[δk,∥ℓ∥^2]]⊕[β^2]$, which is $max{β_0^2,β^2}=β_0^2$. "
    "Finally $C^He=xd+yℓ$ has norm $R(x,y)$. If $α=0$ (so $β_0=β$), "
    "take $d=0$, $ℓ=(β,0,0)^T$, $ℓ_2=(0,β,0)^T$ in the same "
    "construction: $∥C^He∥=βy=R(x,y)$. If $β=0$ (so $β_0=α$), take "
    "$C^H=αvv^H$, for which $∥C^H(xv+yu)∥=αx=R(x,y)$; the zero "
    "operator is trivial. Zero padding reaches any $n≥3,m≥3$. □")

# ---- 4. Theorem 2 attainment enumeration -------------------------------
rep(
    "The theorem asserts three distinct facts: (a) $Q_G$ is safe on the "
    "whole class (Proposition 1); (b) the class infimum of the true "
    "worst-case margin equals $Q_G$; and (c) an explicit construction "
    "sequence approaches the infimum. Attainment of the infimum by a "
    "single instance is not claimed in general: at $r=0$ every instance "
    "attains it trivially, while for other summaries the constructed "
    "family approaches the infimum only along nominal channels of "
    "unbounded norm, and whether some other instance attains it "
    "finitely is left open.",
    "The theorem asserts three distinct facts: (a) $Q_G$ is safe on the "
    "whole class (Proposition 1); (b) the class infimum of the true "
    "worst-case margin equals $Q_G$; and (c) an explicit construction "
    "sequence approaches the infimum. Finite attainment is proved in "
    "three cases: $r=0$ (every instance attains it); $α=0$ with $β>0$ "
    "(the explicit construction below attains (9) exactly); and the "
    "boundary class of Theorem 3 (attained by an explicit instance). "
    "For the remaining nondegenerate cases ($α>0$, $r>0$), the "
    "constructed $ε$-family approaches the infimum along nominal "
    "channels whose norms grow without bound as $ε→0$; this is a "
    "property of that family, not an impossibility claim — whether "
    "some other instance of finite nominal-channel norm attains the "
    "infimum is left open.")

# ---- 5. SDP timing numbers, corrected round only ------------------------
rep(
    "The results do not establish complexity advantages or superiority "
    "over an optimized solver. The generic SDP pipeline on the same "
    "batch is a further 4–13× slower than the tolerance-stopped "
    "reference; its value is cross validation, and the slowness of the "
    "generic conic path is not a statement about full-information "
    "methods.",
    "The results do not establish complexity advantages or superiority "
    "over an optimized solver. On the same batch, the generic SDP "
    "pipeline is 17.3 times slower than the tolerance-stopped "
    "reference as a ratio of medians (median per-pair ratio 15.6; the "
    "reference is faster in 15 of 15 pairs), and the staged "
    "$L_quad→Q_G→$SDP pipeline is 6.5 times slower (median per-pair "
    "ratio 7.2); these are single-round figures from the corrected SDP "
    "study. The SDP pipeline's value is cross validation, and the "
    "slowness of the generic conic path is not a statement about "
    "full-information methods.")

rep(
    "With the tolerance stop, the pure reference runs at 0.041 s per "
    "seed and threshold on the validation batch; the best tested "
    "configuration (new bound with tolerance-stopped fallback, 0.032 "
    "s) is about 20% below it,",
    "With the tolerance stop, the pure reference runs at 0.041 s per "
    "seed and threshold on the validation batch; among the two "
    "pipelines of Table 11, the new bound with tolerance-stopped "
    "fallback (0.032 s) is about 20% below it,")

rep("Version 6.1 • English revision • 7 October 2026",
    "Version 6.2 • English revision • 7 October 2026")

open("paper/paper_v6_2.md", "w", encoding="utf-8", newline="").write(src)
print(f"paper_v6_2.md written ({n_patch} patches)")
