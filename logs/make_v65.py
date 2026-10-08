# -*- coding: utf-8 -*-
"""Patch paper_v6_4.md -> paper/paper_v6_5.md (V6.5 repositioning)."""
src = open("paper/paper_v6_4.md", encoding="utf-8").read()
n = [0]


def rep(old, new, cnt=1):
    assert src.count(old) == cnt, (src.count(old), old[:60])
    globals()['src'] = src.replace(old, new)
    n[0] += 1


# ---- 1. Abstract rewrite (theory-first ordering) ------------------------
i0 = src.index("## Abstract")
i1 = src.index("Keywords:")
rep(src[i0:i1], """## Abstract
Certifying that a fixed beamforming and reconfigurable intelligent
surface (RIS) configuration satisfies quality-of-service constraints
under a norm-bounded channel error is, given the complete channel and
interference geometry, a trust-region problem that exact methods already
solve; this paper therefore does not address full-information
optimization. Instead it asks a partial-information question: when the
certification interface retains only six scalar summaries — the
desired-stream power, the nominal desired and interference amplitudes,
and three gains of the interference operator — what is the strongest
worst-case QoS guarantee that holds uniformly over every instance
consistent with those summaries? The three gains originate from one
operator and obey spectral compatibility that a separate-gain envelope
ignores. Exploiting this compatibility yields a continuous three-branch
directional interference extremum in closed form. The induced QoS
certificate is safe, dominates the separate-gain envelope pointwise
under the same six scalars, and — the main result — equals the infimum
of the true worst-case margin over the entire nondegenerate information
class, so no uniformly better certificate exists at this information
budget. At the degenerate boundary where the desired direction attains
the total norm, a distinct phase-coupled envelope is derived and
attained explicitly. Experiments verify the analytical formulas,
exhibit strict but modest and batch-dependent empirical gains over the
separate-gain envelope, and quantify how norm-restricted information
classes can certify strictly more.

""")

# ---- 2. Intro: sharpen the research question ---------------------------
rep(
    "If only six scalar summaries are retained — the desired-stream "
    "power, the nominal desired and interference amplitudes, and three "
    "gains of the interference operator — a precise question arises: "
    "what is the strongest worst-case QoS guarantee that is uniform "
    "over every instance consistent with these summaries?",
    "This paper fixes the information budget deliberately: six scalar "
    "summaries — the desired-stream power, the nominal desired and "
    "interference amplitudes, and three gains of the interference "
    "operator. The research question is then exact: what is the "
    "strongest worst-case QoS guarantee that can be certified "
    "uniformly over all full channel/interference instances "
    "consistent with exactly these summaries? The paper does not "
    "study full-information robust beamforming optimization, optimal "
    "RIS phase design, the optimal compression dimension, or whether "
    "six scalars are the best descriptor; it studies what the "
    "spectral compatibility induced by the common operator origin of "
    "the three gains contributes to the uniform certificate at this "
    "fixed budget.")

# ---- 3. Contributions: 3 points, theory-first --------------------------
i0 = src.index("This work makes three contributions.")
i1 = src.index("## 2 Fixed configuration model")
rep(src[i0:i1], """This work makes three contributions. First, a closed-form
spectral-compatible directional extremum: from the joint constraints
imposed by the desired-direction gain, the orthogonal-subspace gain
and the total spectral norm of one interference operator, an explicit
continuous three-branch maximum of the interference output under a
shared error budget is derived and attained exactly over the
three-spectral operator class. Second, information-class optimality:
the induced shared-budget QoS envelope is safe, dominates the envelope
that uses the three gains separately under the same six scalars, and
equals the information-class infimum away from the degenerate
boundary, with a separately derived and explicitly attained
phase-coupled envelope at that boundary. Third, numerical verification
of the analytical claims — formula correctness, construction
attainment, boundary behaviour — together with an honest evaluation of
the conservatism and complementarity of the six-scalar interface,
including bounded-nominal-norm experiments that quantify what
additional physical information would buy.

""")

# ---- 4. Sec 2.3: three levels + sufficiency sentence -------------------
rep(
    "Sections 3–5 characterize this quantity from the retained "
    "spectral information:",
    "Three levels should be kept apart: the full physical instance "
    "$𝒳$; the information map $ℐ$ that compresses it to $s$; and the "
    "information class $𝒞(s)$. The quantity $B^*(s)$ is not the true "
    "robust margin of any particular physical instance; it is the best "
    "uniform certificate that any method observing only the retained "
    "summaries can guarantee. Information-class sharpness is "
    "therefore a statement of information sufficiency, not "
    "instance-wise exactness: $Q_G$ may be conservative for a "
    "particular fixed operator while remaining optimal uniformly over "
    "the entire information class. Sections 3–5 characterize this "
    "quantity from the retained spectral information:")

# ---- 5. New subsection: scope of sharpness + restricted class ----------
rep(
    "## 6 Numerical certification and reference interface",
    """### 5.2 Scope of information-class sharpness and norm-restricted
classes

Theorem 2 characterizes $B^*(s)$ for the *unrestricted* information
class $𝒞(s)$, in which the nominal-channel norm is unconstrained. For
any additional physically meaningful restriction the admissible class
shrinks, and the optimal certificate can only improve. The most
natural restriction retains the nominal-channel norm bound
$∥ĥ∥≤H_{max}$:

$𝒞(s,H_{max})={𝒳∈𝒞(s):∥ĥ∥≤H_{max}}$,   $B_{H}^*(s)=min_{𝒳∈𝒞(s,H_{max})} q^*(𝒳)$.

**Remark (monotonicity).** Since $𝒞(s,H_{max})⊆𝒞(s)$, immediately

$B_{H}^*(s)≥B^*(s)=Q_G(s)$  whenever  $H_{max}≥\\min_{𝒳∈𝒞(s)}∥ĥ∥$,

and $B_{H}^*$ is nonincreasing in $H_{max}$. $Q_G$ remains *safe* on
$𝒞(s,H_{max})$, but it need not remain *sharp*: sharpness of the
unrestricted class does not transfer to restricted classes, and the
complete envelope of the $H$-branch is open. The divergence observed
in the families of Section 5 does not invalidate Theorem 2; it "
"reveals that the unrestricted class admits nominal-channel "
"geometries that a norm constraint may exclude, so the boundary "
"discontinuity is information-class dependent. Divergence occurs in "
"particular families; bounded and finitely attaining constructions "
"exist in others; whether every summary vector admits a finite-norm "
"attaining instance remains open. Section 7.6 quantifies empirically "
"how much a norm restriction can tighten the class optimum.

## 6 Numerical certification and reference interface""")

# ---- 6. New experiment subsection 7.6 (before 7.7 timing) --------------
rep(
    "### 7.7 Timing and implementation sensitivity",
    """### 7.6 Effect of a bounded nominal-channel norm

Two numerical studies estimate how much the unrestricted sharp
envelope loses when the class is restricted by $∥ĥ∥≤H_{max}$. The "
"estimates are bracketed: $Q_G(s)=B^*(s)$ is a rigorous lower bound "
"for $B_{H}^*(s)$, and a multi-start feasible search over class "
"instances (random unitary conjugations fixing $v$, exact matching of "
"$a$ and $I_0$ on random orthogonal directions, exact trust-region "
"evaluation) provides an *upper* estimate; no global optimality of "
"the restricted search is claimed.

*Boundary-approach family (Section 5).* For the family with "
"$δ∈\\{0.1,0.03,0.01\\}$ and $H_{max}∈\\{2,…,100\\}$, the upper "
"estimate stays near $0.59$ — close to the boundary-class value — "
"whenever the divergent family member (norms $20.0$, $66.7$, "
"$200.0$) is excluded, and drops to $Q_G≈−0.21$ exactly when "
"$H_{max}$ admits it. The restricted gap upper estimate reaches "
"$≈0.8$: for these summaries, the unrestricted sharp value is driven "
"by large-norm geometry, and a norm bound that excludes it can "
"tighten the class optimum back toward the boundary value.

*Pool instances (validation seeds 211–213, three users each, "
"$H_{max}=c∥ĥ∥$, $c∈\\{1,…,5\\}$, 120 starts).* The search found no "
"class member below the actual instance even at $c=5$, so the gap "
"upper estimate equals the actual-instance envelope looseness "
"(median $0.14$, 90th percentile $0.33$, maximum $0.54$ of 45 "
"cases); none was near zero. Whether the restricted infimum for "
"these summaries genuinely exceeds $Q_G$ remains unresolved — the "
"upper estimate is not tight — so the pool study neither proves nor "
"excludes a substantial restricted-class improvement.

Together, the two studies show that the unrestricted sharpness of "
"Theorem 2 is tight in a class-dependent way: norm restrictions can "
"certify strictly more for boundary-approaching summaries, while for "
"generic pool summaries the question is open. A complete "
"$H$-restricted envelope is left as future work.

### 7.7 Timing and implementation sensitivity""")

# ---- 7. L_quad reframing sentence --------------------------------------
rep(
    "The bounds therefore have no universal dominance relation; their "
    "maximum is safe and yields the highest direct-certification "
    "coverage among tested bounds.",
    "The bounds therefore have no universal dominance relation. This "
    "does not contradict the sharpness of $Q_G$: $L_quad$ observes "
    "an additional scalar $κ$ and operates on a strictly richer "
    "information interface, whereas Theorem 2 says that no "
    "certificate using *only* the six retained summaries can "
    "uniformly improve $Q_G$ — it says nothing about certificates "
    "with more information. The joint rule $max{Q_G,L_quad}$ is "
    "accordingly reported as evidence that the two information "
    "branches carry complementary information, not as a new "
    "algorithm; it is safe and yields the highest "
    "direct-certification coverage among the tested bounds.")

# ---- 8. Discussion: four limitations, replace current text -------------
i0 = src.index("## 9 Discussion")
i1 = src.index("## 10 Conclusion")
rep(src[i0:i1], """## 9 Discussion
Four limitations bound the claims. (i) The six-scalar interface is
fixed by the problem formulation; the paper does not prove that six
scalars are Pareto-optimal among compressed representations, and the
interface study shows no runtime or storage advantage over a "
"seven-scalar branch. (ii) Information-class sharpness holds for the "
"unrestricted class; additional physical restrictions — including a "
"bounded nominal-channel norm — can only tighten the class optimum, "
"and the boundary-approach experiments show the tightening can be "
"large for some summaries. (iii) The empirical advantage over the "
"separate-gain envelope is strict but modest and dataset dependent "
"(zero in a saturated batch; +1 to +10 units elsewhere). (iv) The "
"floating-point implementation is numerically guarded but is not a "
"formally verified interval-arithmetic certificate. The paper's "
"retained significance is that it completely characterizes what can "
"and cannot be guaranteed from this particular information map: the "
"envelope, its sharpness, and the degenerate boundary are all "
"explicit, so any future richer-interface certificate can be "
"measured against a known optimum rather than an open gap.

""")

# ---- 9. Conclusion rewrite ----------------------------------------------
i0 = src.index("## 10 Conclusion")
i1 = src.index("## Appendix A")
rep(src[i0:i1], """## 10 Conclusion
This paper studies robust QoS certification at a fixed information
budget of six scalar summaries. Because the three spectral gains
originate from one interference operator, their compatibility admits
an explicit three-branch directional extremum; lifting it through the
shared error budget yields a QoS certificate that is safe, pointwise
dominates the separate-gain envelope under the same six scalars, and
— the central result — equals the information-class infimum away from
the degenerate boundary, with an explicitly attained phase-coupled
envelope at that boundary. The class-optimal value is a statement of
information sufficiency, not instance-wise exactness, and norm-"
"restricted classes can certify strictly more, as the bounded-norm "
"experiments quantify. A natural next step is to characterize the "
"optimal certificate under richer or physically restricted "
"information classes, including classes that retain a bound on the "
"nominal-channel norm.

""")

rep("Version 6.4 • English revision • 7 October 2026",
    "Version 6.5 • English revision • 8 October 2026")

open("paper/paper_v6_5.md", "w", encoding="utf-8", newline="").write(src)
print(f"paper_v6_5.md written ({n[0]} patches)")
