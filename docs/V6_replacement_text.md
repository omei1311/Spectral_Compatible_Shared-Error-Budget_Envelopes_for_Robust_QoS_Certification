# V6 Replacement Text (Part E) — IEEE style, ready to paste

All novelty statements use the bounded form; no "first/novel" claims.

**DECISION APPLIED (2026-10-07, per author)**: Title = "Sharp Robust QoS
Certification from Six-Scalar Spectral Summaries" (RIS-suffixed variant
retained as fallback only). §2.3 rewritten as "Six-Scalar Information Model
and Certification Hierarchy" with the information map as the problem
definition. §2.4 restructured (math → wireless meaning → why the old bound
is loose) with over-strong wordings removed. Terminology conventions and
the first-three-pages spine added below. Candidate-title analysis kept for
the record; the fused title supersedes it.

## Terminology conventions (apply throughout V6)

Prefer: "six-scalar information class", "six-scalar information budget",
"partial-information certificate", "certificate based on the retained
summaries", "information map". The phrase "six-scalar information model"
is acceptable when first introduced. Avoid "six-scalar channel model" —
the full system still exists; only the certificate's observation is
restricted.

## First-three-pages spine (assembly checklist)

By page 3 the reader must have traversed: full-information QoS
verification is solved → we deliberately restrict the information
available to the certificate → the six scalars define an information
class → the mathematically relevant target is inf{q^*(X): I(X)=s} →
existing separate-gain bounds ignore operator compatibility → we solve
the compatible directional extremum → we prove the resulting QoS bound
equals the information-class infimum (nondegenerate case), with the
boundary class treated separately.

## Title (decided)

Sharp Robust QoS Certification from Six-Scalar Spectral Summaries

Fallback if RIS retrieval keywords are required:
Sharp Robust QoS Certification from Six-Scalar Spectral Summaries for
Fixed RIS Configurations

Rationale (author's): the four terms map to the four selling points —
"Sharp" to Theorem 2 sharpness and Theorem 3 attainment (not mere
tightness); "Robust QoS Certification" to fixed-configuration
verification rather than beamforming optimization; "Six-Scalar" promotes
the information budget to a title-level object; "Spectral Summaries" (not
"Shared Error Budget", which is a technique) names the
information-class formulation. Rejected: "Optimal Six-Scalar
Certification..." (invites a global-optimality misreading that the
abstract would have to fight); the pure "Information-Limited ... for
Fixed Beamforming and RIS Configurations" (too safe; hides six-scalar and
sharpness).

## Earlier title candidates (superseded; kept for the record)

A. Information-limited certification emphasis
1. "Sharp Robust QoS Certification under Spectral Summary Information"
   — math-forward, "sharp" is earned (Thm 2), no application context;
     risk: reviewers expect operator theory.
2. "Information-Limited Robust QoS Certification with Spectral-Compatible
   Shared Error Budgets"
   — keeps the original technical term as subtitle; longest, precise.

B. Six-scalar emphasis
3. "Six-Scalar Spectral Envelopes for Robust QoS Certification of Fixed
   RIS Configurations"
   — immediately signals the information model + application; slight risk
     of reading as "yet another RIS heuristic".
4. "Optimal Six-Scalar Certification of Worst-Case QoS under Bounded CSI
   Error"
   — "optimal" is defensible only with the sharpness qualifier nearby.

C. RIS application context retained
5. "Information-Limited Robust QoS Certification for Fixed Beamforming
   and RIS Configurations"
   — clearest scope statement; least mathematical signal.
6. "Robust QoS Certification from Scalar Spectral Summaries: Envelopes,
   Sharpness, and a Degenerate Boundary"
   — signals all three theory pillars; no wireless keywords (title-level
     discoverability cost).

Recommendation rationale (not a decision): candidates 2/5 pair the
information-limit framing with the retained technical vocabulary; 6 is the
most accurate for a theory-heavy venue; 3 for a wireless venue.

## Abstract (201 words)

Certifying that a fixed beamforming and reconfigurable intelligent surface
(RIS) configuration satisfies quality-of-service constraints under
norm-bounded channel-state error is, with full geometry, a trust region
problem solvable exactly. This paper asks what can be guaranteed when the
certification interface retains only six scalar summaries of the nominal
channel and the interference operator: desired and nominal amplitudes,
desired-direction gain, orthogonal-subspace gain, and total spectral norm.
Because these gains originate from one operator, they satisfy spectral
compatibility that previous shared-budget envelopes ignore. We solve the
resulting directional interference extremum in closed form, obtaining a
continuous three-branch envelope and a safe worst-case QoS bound that
dominates the previous envelope under the same information. For the
nondegenerate information class we prove the bound equals the infimum of
the true worst-case margin over all instances sharing the six summaries;
at the degenerate boundary where the desired direction attains the total
norm, a distinct phase-coupled envelope is derived and attained
explicitly. Numerical checks reproduce all analytical envelopes, two
independent exact methods agree on every validation problem, and the
certification gain over the previous six-scalar bound is strict but
modest and batch-dependent; the certificate is complementary to a cheap
full-geometry quadratic bound in a staged screening interface.

## 1 Introduction (5 paragraphs)

[P1 — setting] After beamforming vectors and discrete RIS phases have been
fixed for a candidate configuration, deployment, reuse, and screening
still require verifying that every user's QoS constraint survives
channel-state uncertainty. With the complete channel and interference
geometry available, this verification is a ball-constrained quadratic
problem that trust-region methods and the S-procedure solve exactly
[2],[3],[4]; robust transmission design with imperfect CSI, including
bounded-error RIS settings, is correspondingly mature [1],[4]. The
existence of exact full-information solvers is therefore not the issue
this paper addresses.

[P2 — information limit] In candidate screening over many configurations,
in distributed base-station and RIS coordination, and in control-plane
storage of configuration certificates, maintaining or exchanging the full
quadratic geometry of every candidate may be undesirable; a compact,
fixed-size summary is attractive, and evaluating one certificate at many
uncertainty radii is common. If only six scalar summaries are retained —
the desired-stream power, the nominal desired and interference amplitudes,
and three gains of the interference operator — a natural question arises:
what is the strongest worst-case QoS guarantee that is uniform over every
instance consistent with these summaries? Envelopes that use the three
gains separately, as is natural in the bounded-error robust beamforming
literature, ignore the fact that all three originate from one operator
and therefore cannot attain their individual worst cases simultaneously.

[P3 — method] We answer the question in three steps. First, we solve the
directional interference extremum under the joint spectral constraints:
given the desired-direction gain α, the orthogonal-subspace gain β, and
the total spectral norm β0, we evaluate the largest interference output
that a parallel error amplitude x and an orthogonal budget y can jointly
produce. The extremum is a continuous three-branch closed form obtained
from a two-by-two positive semidefinite completion, a specialization of
classical contraction-completion and refined Cauchy–Schwarz (Buzano-type)
arguments [5],[6],[7]. Second, we embed the extremum into the worst-case
QoS margin, obtaining a shared-error-budget lower envelope. Third, we ask
whether any safer certificate using the same six scalars can be strictly
better.

[P4 — results] The envelope is safe and dominates the previous
shared-budget envelope pointwise. Over the nondegenerate information
class (β0 > α) the envelope equals the infimum of the true worst-case
margin over all instances with the same six summaries — the sharpness
statement, in the sense of the radius of information from optimal
recovery [8],[9]. At the degenerate boundary β0 = α, the desired
direction already attains the operator norm; singular-vector geometry
forces an additional orthogonality, hence a phase coupling between the
desired signal and the parallel interference, and the optimal envelope of
the boundary class is different and is attained by an explicit diagonal
construction. The two envelopes can differ sharply for the same boundary
summaries (0.59 versus −0.21 in a constructed example); the discontinuity
is a property of the information classes, not of any physical channel.

[P5 — numerics] The numerical study has five roles: verifying the
analytical implementation to machine precision; exhibiting constructed
instances where the new envelope certifies feasibility that the previous
envelope cannot; validating on preregistered, previously unused seeds;
showing that the six-scalar certificate and a cheap full-geometry
quadratic bound are complementary (no universal dominance, safe to
combine by maximization); and cross-validating every certification
statement with an independent exact S-procedure SDP pipeline. The
experiments support implementation correctness and complementary use;
they do not prove the theorems, and the empirical acceptance gain over
the previous six-scalar envelope, while positive in every batch, is
modest and batch-dependent.

## Contributions paragraph (3 points)

This work makes three contributions. First, we derive the
spectral-compatible directional envelope: from the joint constraints
imposed by the desired-direction gain, the orthogonal-subspace gain, and
the total spectral norm of one interference operator, we obtain an
explicit continuous three-branch extremum for the worst interference
output under a shared error budget. Second, we establish
information-class optimality: the induced shared-budget QoS envelope is
safe, dominates the envelope that uses the three gains separately, equals
the information-class infimum away from the degenerate boundary, and is
complemented by a separate, explicitly attained phase-coupled envelope at
the boundary β0 = α. Third, we specify and validate a certification
interface: exact trust-region and S-procedure references confirm the
implementation on every validation problem, the envelope produces
certifications the previous six-scalar envelope cannot, and, having no
universal dominance relation with a full-geometry quadratic bound, it
serves as a complementary stage between cheap bounds and exact fallback.

## 2.3 Six-Scalar Information Model and Certification Hierarchy

The six retained summaries define an information map, not a physical
channel model. Let ℐ(𝒳) denote the six-scalar summary of a full
configuration instance 𝒳=(b,ĥ,C) at a user, and let 𝒞(s)={𝒳:ℐ(𝒳)=s} be
the class of all instances indistinguishable to a certificate that
observes only s. A six-scalar certificate is any function B(s) of the
retained summaries that must remain valid across the entire class,

    B(s) ≤ q*(𝒳)   for all 𝒳 ∈ 𝒞(s),

and the strongest possible uniform certificate under this information
budget is therefore

    B*(s) = inf_{𝒳 ∈ 𝒞(s)} q*(𝒳).

Sections 3–5 characterize this quantity from the retained spectral
information: Section 3 solves the directional extremum forced by the
operator origin of the summaries; Section 4 proves that the resulting
envelope equals B*(s) on the nondegenerate class and characterizes the
separate optimum of the boundary class β0=α.

This formulation places the paper in an explicit certification
hierarchy: increasing retained information shrinks the admissible
information class and can only tighten the strongest uniform
certificate, so certificate strength trades against information
retained. At the coarse end, separable norm bounds use the full error
radius twice (Level 0). The shared-budget envelope that uses the three
interference gains separately (Level 1) is the natural envelope of the
bounded-error robust beamforming literature. Level 2 — a shared budget
with joint spectral compatibility — is the object of this paper. Level 3
retains one additional scalar, the nominal channel norm, for which we
give a necessary compatibility condition but no complete envelope.
Level 4 is the full quadratic geometry, where cheap bounds such as the
quadratic certificate L_quad and exact trust-region or S-procedure
solvers operate. The hierarchy also explains a recurring empirical
observation: a full-geometry bound sometimes certifies more than the
six-scalar envelope not because it is a better algorithm but because it
retains more information. Levels 2 and 4 are complementary, and our
experiments combine them.

Such a fixed-size interface can arise when candidate configurations are
screened repeatedly, when configuration descriptors are cached or
exchanged across network components, or when one wishes to evaluate a
stored candidate at multiple uncertainty radii without retaining the
complete quadratic geometry. These scenarios motivate the abstraction
but are not assumptions required by the theory.

## 2.4 Physical interpretation of the retained summaries

The retained quantities separate the desired direction from its
orthogonal complement and summarize how one interference operator acts
on these two subspaces. In wireless terms: p is the desired-stream
transmit power; a is the nominal desired-signal amplitude along the
desired beam; I0 is the nominal interference amplitude, whose square
contributes to the interference power. The three spectral quantities
describe how the interference operator maps channel error to
interference: α is its response to error parallel to the desired beam;
β is the largest amplification any error component orthogonal to the
beam can suffer; β0 is the largest amplification any error at all can
suffer.

The three spectral quantities should not be viewed as independent
uncertainty parameters. They are coupled observations of the same
interference operator. The desired-direction gain describes the
operator response to error parallel to the desired beam, the orthogonal
gain bounds amplification over the complementary subspace, and the total
spectral norm limits both simultaneously. Spectral compatibility is
therefore a structural consistency condition: a large response in one
direction restricts what the same operator can realize in orthogonal
directions.

A bound that uses these gains separately may combine individually valid
worst cases that cannot be realized simultaneously by one operator and
one shared error vector. The proposed directional extremum removes
exactly this incompatibility while retaining the same six scalar
summaries. Its improvement over the previous envelope therefore comes
from using more of the structure already contained in the retained
information, rather than from introducing additional channel
information. This incompatible combination of separate worst cases is
the specific source of conservatism removed by the spectral-compatible
envelope.
## Boundary intuition paragraph (insert before Theorem 3)

The degenerate case β0 = α is not a technicality but a change of regime.
Whenever the beam direction attains the total spectral norm, v is a
singular direction of the interference operator, and the equality case
of the singular-vector characterization forces the interference generated
by any orthogonal error component to be orthogonal to the interference
generated along the beam. The parallel interference output is then a
fixed complex multiple of the desired signal output — they share one
phase — whereas in the nondegenerate regime the two phases can be chosen
independently, which is exactly what the completion construction of
Theorem 2 exploits. Consequently the boundary information class carries
an additional phase compatibility that the six scalars do not display,
and its optimal envelope must be derived separately. In particular, when
nondegenerate families approach the boundary, the class over which the
infimum is taken changes: the construction families drift to unbounded
nominal channel norm while the boundary class contains bounded channels,
so the optimum jumps (0.59 versus −0.21 in Section 5's example) — an
information-class discontinuity, not a physical one.

## Related Work section (with table)

Robust transmission under imperfect CSI splits along the information
available. Full-geometry worst-case design under norm-bounded errors is
classical: S-procedure and semidefinite reformulations deliver exact or
conservative constraints for MISO downlink, relay, secure, and RIS-aided
systems, with [1] an early IRS-aided framework and [4] a current
multi-functional RIS design under bounded cascaded CSI error; exact
solution of the ball-constrained quadratic rests on the trust-region and
S-lemma literature [2],[3]. In this line, decomposing the channel error
into components parallel and orthogonal to the desired channel and
bounding their effects separately is a recurring device; the envelope we
call "previous" is the shared-budget bound obtained when the three
interference gains are used separately, implicit in this literature.
A second line compresses information: probabilistic SINR constraints use
distributional summaries with conservative approximations. Our question
differs from both: the uncertainty is deterministic and bounded, and the
compressed object is a fixed set of spectral summaries of one operator.

The mathematics we need is classical but dispersed. The cross-term bound
in our extremum is a Buzano-type refinement of Cauchy–Schwarz [5],[6];
the completion step — extending a contraction with prescribed action on
a vector — belongs to the Parrott/operator-ball and chordal PSD
completion theorems [7],[10]; and interpreting an infimum over an
information class as the worst-case radius of a certificate is the
radius-of-information viewpoint of optimal recovery [8],[9]. To the best
of our knowledge, we have not found an existing result that evaluates
the directional extremum in closed form under exactly the three retained
spectral summaries and a shared error budget, nor an existing treatment
of the β0 = α boundary's phase-coupled envelope; we do not claim these
ingredients individually as new.

[Insert comparison table from Part B, translated to English columns.]

## Discussion (replacement, consolidated; each limitation stated once)

The sharpness statement is intentionally strong and narrow: it holds for
the unrestricted information class, where extremal constructions may
require large nominal channel norms that a physical RIS channel may not
admit. Restricting the class — for example by retaining the channel norm
H, for which (28) gives a necessary compatibility condition — can only
shrink the admissible set, and the six-scalar envelope need not remain
sharp; the complete seven-scalar envelope is open. On the numerical
side, the analytical inequalities are exact, while the implementation
uses floating-point arithmetic with conservative numerical safeguards
rather than formally verified directed rounding; the experiments support
correctness and safety on all tested instances but do not constitute a
proof. Empirically, the certificate gain over the previous six-scalar
envelope is positive in every batch yet small and batch-dependent
(+2.12 and +0.63 percentage points of feasible coverage), and no stable
runtime advantage over the previous envelope is observed; measured
speedups of screening pipelines are tied to the tested implementation,
and part of the apparent saving in early measurements traced to
multiplier-search over-iteration in the reference solver. A formula-level
comparison with operator interpolation, contraction completion, and
information-based complexity [7],[8],[9] has been carried out to the
extent above; the positioning claims are limited accordingly.

## Conclusion (three layers)

Mathematically, joint spectral compatibility of three gains of one
interference operator yields an explicit three-branch directional
extremum, a safe shared-budget QoS envelope dominating the
separately-used-gain envelope, equality with the information-class
infimum away from the degenerate boundary, and a separately derived,
explicitly attained phase-coupled envelope at that boundary. Practically,
the envelope is a fixed-size certificate whose empirical gain over the
previous six-scalar bound is strict but modest, and which is
complementary to a cheap full-geometry quadratic bound — the two have no
universal dominance relation and combine safely by maximization, giving
the highest direct-certification coverage among the tested bounds. In
scope, the certificate does not replace exact full-information solvers,
does not optimize candidate generation, and its sharpness concerns the
stated unrestricted information class; natural extensions include
additional retained information such as the channel norm, physically
restricted information classes, validated numerical certificates, and
lifecycle management of stored configuration certificates.
