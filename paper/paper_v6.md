# Sharp Robust QoS Certification from Six-Scalar Spectral Summaries
Fixed beamforming and RIS configurations under a partial-information model
Version 6.0 • English revision • 7 October 2026

## Abstract
Certifying that a fixed beamforming and reconfigurable intelligent surface (RIS) configuration satisfies quality-of-service constraints under norm-bounded channel-state error is, with full geometry, a trust region problem solvable exactly. This paper asks what can be guaranteed when the certification interface retains only six scalar summaries of the nominal channel and the interference operator: desired and nominal amplitudes, desired-direction gain, orthogonal-subspace gain, and total spectral norm. Because these gains originate from one operator, they satisfy spectral compatibility that previous shared-budget envelopes ignore. We solve the resulting directional interference extremum in closed form, obtaining a continuous three-branch envelope and a safe worst-case QoS bound that dominates the previous envelope under the same information. For the nondegenerate information class we prove the bound equals the infimum of the true worst-case margin over all instances sharing the six summaries; at the degenerate boundary where the desired direction attains the total norm, a distinct phase-coupled envelope is derived and attained explicitly. Numerical checks reproduce all analytical envelopes, two independent exact methods agree on every validation problem, and the certification gain over the previous six-scalar bound is strict but modest and batch-dependent; the certificate is complementary to a cheap full-geometry quadratic bound in a staged screening interface.

Keywords: robust QoS certification; partial information; information class; spectral compatibility; RIS; fixed configuration.

## 1 Introduction
After beamforming vectors and discrete RIS phases have been fixed for a candidate configuration, deployment, reuse, and screening still require verifying that every user's QoS constraint survives channel-state uncertainty. With the complete channel and interference geometry available, this verification is a ball-constrained quadratic problem that trust-region methods and the S-procedure solve exactly [2],[3],[4]; robust transmission design with imperfect CSI, including bounded-error RIS settings, is correspondingly mature [1],[4]. The existence of exact full-information solvers is therefore not the issue this paper addresses.

In candidate screening over many configurations, when configuration descriptors are cached or exchanged across network components, or when a stored candidate is evaluated at several uncertainty radii, a compact fixed-size summary can be preferable to carrying the full quadratic geometry. If only six scalar summaries are retained — the desired-stream power, the nominal desired and interference amplitudes, and three gains of the interference operator — a precise question arises: what is the strongest worst-case QoS guarantee that is uniform over every instance consistent with these summaries? Envelopes that use the three gains separately, as is natural in the bounded-error robust beamforming literature [9], ignore the fact that all three originate from one operator and therefore cannot attain their individual worst cases simultaneously.

We answer the question in three steps. First, we solve the directional interference extremum under the joint spectral constraints: given the desired-direction gain α, the orthogonal-subspace gain β, and the total spectral norm β_0, we evaluate the largest interference output that a parallel error amplitude x and an orthogonal budget y can jointly produce. The extremum is a continuous three-branch closed form obtained from a two-by-two positive semidefinite completion, a specialization of classical contraction-completion [7] and refined Cauchy–Schwarz (Buzano-type) arguments [5], with the chordal completion characterization [6] in the background. Second, we embed the extremum into the worst-case QoS margin, obtaining a shared-error-budget lower envelope. Third, we ask whether any safer certificate using the same six scalars can be strictly better, and we answer in the negative away from a degenerate spectral boundary.

The envelope is safe and dominates the previous shared-budget envelope pointwise. Over the nondegenerate information class (β_0>α) the envelope equals the infimum of the true worst-case margin over all instances with the same six summaries — the sharpness statement, in the sense of the radius of information from optimal recovery [8]. At the degenerate boundary β_0=α, the desired direction already attains the operator norm; singular-vector geometry forces an additional orthogonality, hence a phase coupling between the desired signal and the parallel interference, and the optimal envelope of the boundary class is different and is attained by an explicit diagonal construction. The two envelopes can differ sharply for the same boundary summaries (0.59 versus −0.21 in a constructed example); the discontinuity is a property of the information classes, not of any physical channel.

The numerical study has five roles: verifying the analytical implementation to machine precision; exhibiting constructed instances where the new envelope certifies feasibility that the previous envelope cannot; validating on preregistered, previously unused seeds, including a pool with distance-dependent path loss and Rician components; showing that the six-scalar certificate and a cheap full-geometry quadratic bound are complementary (no universal dominance, safe to combine by maximization); and cross-validating every certification statement with an independent exact S-procedure SDP pipeline. The experiments support implementation correctness and complementary use; they do not prove the theorems, and the empirical acceptance gain over the previous six-scalar envelope, while positive in every batch, is modest and batch-dependent.

This work makes three contributions. First, we derive the spectral-compatible directional envelope: from the joint constraints imposed by the desired-direction gain, the orthogonal-subspace gain, and the total spectral norm of one interference operator, we obtain an explicit continuous three-branch extremum for the worst interference output under a shared error budget. Second, we establish information-class optimality: the induced shared-budget QoS envelope is safe, dominates the envelope that uses the three gains separately, equals the information-class infimum away from the degenerate boundary, and is complemented by a separate, explicitly attained phase-coupled envelope at the boundary β_0=α. Third, we specify and validate a certification interface: exact trust-region and S-procedure references confirm the implementation on every validation problem, the envelope produces certifications the previous six-scalar envelope cannot, and, having no universal dominance relation with a full-geometry quadratic bound, it serves as a complementary stage between cheap bounds and exact fallback.

## 2 Fixed configuration model, information class, and related work
### 2.1 Robust QoS margin
For a user, let $b∈ℂ^n$ be the nonzero desired beam, $C∈ℂ^{n×m}$ the matrix of interfering beams, and $ĥ∈ℂ^n$ the nominal stacked channel. The actual channel is $ĥ+e$, with $∥e∥≤r$. The SINR requirement is $γ>0$, and $σ^2≥0$ denotes noise power, not noise amplitude. We suppress user indices and define

EQ (1): q^*(r)=min_{∥e∥≤r} { |b^H(ĥ+e)|^2−γ∥C^H(ĥ+e)∥^2−γσ^2 }.

All users satisfying $q_u^*(r_u)≥0$ makes the fixed candidate robustly feasible. A negative lower bound alone does not establish infeasibility.

### 2.2 Six scalar summaries
Set $v=b/∥b∥$ and $P=I−vv^H$. The retained summaries are

EQ (2a): p=∥b∥^2,   a=|v^Hĥ|,   I_0=∥C^Hĥ∥.
EQ (2b): α=∥C^Hv∥,   β=∥C^HP∥_2,   β_0=∥C∥_2.

Their necessary spectral compatibility condition is

EQ (3): max{α,β}≤β_0≤√{α^2+β^2}.

### 2.3 Six-scalar information model and certification hierarchy
The six retained summaries define an information map, not a physical channel model. Let $ℐ(𝒳)$ denote the six-scalar summary of a full configuration instance $𝒳=(b,ĥ,C)$ at a user, and let

EQ (2c): 𝒞(s)={𝒳:ℐ(𝒳)=s}

be the class of all instances indistinguishable to a certificate that observes only $s$. A six-scalar certificate is any function $B(s)$ of the retained summaries that must remain valid across the entire class, $B(s)≤q^*(𝒳)$ for all $𝒳∈𝒞(s)$, and the strongest possible uniform certificate under this information budget is therefore

EQ (2d): B^*(s)=min_{𝒳∈𝒞(s)} q^*(𝒳).

Sections 3–5 characterize this quantity from the retained spectral information: Section 3 solves the directional extremum forced by the operator origin of the summaries; Section 4 proves that the resulting envelope equals $B^*(s)$ on the nondegenerate class and characterizes the separate optimum of the boundary class $β_0=α$.

This formulation places the paper in an explicit certification hierarchy: increasing retained information shrinks the admissible information class and can only tighten the strongest uniform certificate, so certificate strength trades against information retained. At the coarse end, separable norm bounds use the full error radius twice (Level 0). The shared-budget envelope that uses the three interference gains separately (Level 1) is the natural envelope of the bounded-error robust beamforming literature. Level 2 — a shared budget with joint spectral compatibility — is the object of this paper. Level 3 retains one additional scalar, the nominal channel norm, for which we give a necessary compatibility condition but no complete envelope. Level 4 is the full quadratic geometry, where cheap bounds such as the quadratic certificate $L_quad$ and exact trust-region or S-procedure solvers operate. The hierarchy also explains a recurring empirical observation: a full-geometry bound sometimes certifies more than the six-scalar envelope not because it is a better algorithm but because it retains more information. Levels 2 and 4 are complementary, and our experiments combine them.

Such a fixed-size interface can arise when candidate configurations are screened repeatedly, when configuration descriptors are cached or exchanged across network components, or when one wishes to evaluate a stored candidate at multiple uncertainty radii without retaining the complete quadratic geometry. These scenarios motivate the abstraction but are not assumptions required by the theory.

### 2.4 Physical interpretation of the retained summaries
The retained quantities separate the desired direction from its orthogonal complement and summarize how one interference operator acts on these two subspaces. In wireless terms: $p$ is the desired-stream transmit power; $a$ is the nominal desired-signal amplitude along the desired beam; $I_0$ is the nominal interference amplitude, whose square contributes to the interference power. The three spectral quantities describe how the interference operator maps channel error to interference: $α$ is its response to error parallel to the desired beam; $β$ is the largest amplification any error component orthogonal to the beam can suffer; $β_0$ is the largest amplification any error at all can suffer.

The three spectral quantities should not be viewed as independent uncertainty parameters. They are coupled observations of the same interference operator. Spectral compatibility is therefore a structural consistency condition: a large response in one direction restricts what the same operator can realize in orthogonal directions.

A bound that uses these gains separately may combine individually valid worst cases that cannot be realized simultaneously by one operator and one shared error vector. The proposed directional extremum removes exactly this incompatibility while retaining the same six scalar summaries. Its improvement over the previous envelope therefore comes from using more of the structure already contained in the retained information, rather than from introducing additional channel information.

### 2.5 Related work
Robust transmission under imperfect CSI splits along the information available. Full-geometry worst-case design under norm-bounded errors is classical: S-procedure and semidefinite reformulations deliver exact or conservative constraints for MISO downlink, relay, secure, and RIS-aided systems, with [1] an early IRS-aided framework and [4] a current multi-functional RIS design under bounded cascaded CSI error; exact solution of the ball-constrained quadratic rests on the trust-region and S-lemma literature [2],[3]. In this line, decomposing the channel error into components parallel and orthogonal to the desired channel and bounding their effects separately is a recurring device [9]; the envelope we call "previous" is the shared-budget bound obtained when the three interference gains are used separately, implicit in this literature. A second line compresses information probabilistically, using distributional summaries with conservative approximations. Our question differs from both: the uncertainty is deterministic and bounded, and the compressed object is a fixed set of spectral summaries of one operator.

TABLE relatedwork

The mathematics we need is classical but dispersed. The cross-term bound in our extremum is a Buzano-type refinement of Cauchy–Schwarz [5]; the completion step — extending a contraction with prescribed action on a vector — belongs to the Parrott/operator-ball theorems [7], with the chordal positive-semidefinite completion characterization [6] behind the two-by-two case; and interpreting an infimum over an information class as the worst-case radius of a certificate is the radius-of-information viewpoint of optimal recovery [8]. To the best of our knowledge, we have not found an existing result that evaluates the directional extremum in closed form under exactly the three retained spectral summaries and a shared error budget, nor an existing treatment of the β_0=α boundary's phase-coupled envelope; we do not claim these ingredients individually as new.

## 3 Spectral compatible directional envelope
Write $e=tv+w$, where $v^Hw=0$, $x=|t|$ and $∥w∥≤y$. Put $d=C^Hv$. For an orthogonal unit or subunit input $u$, put $ℓ=C^Hu$ and $s=∥ℓ∥≤β$. Since the map with columns $v,u$ is a contraction, its compressed output Gram matrix is bounded by $β_0^2I$. Define $δ=√{β_0^2−α^2}$ and $ν=√{β_0^2−β^2}$. Cauchy–Schwarz and the determinant condition give

EQ (4): |d^Hℓ|≤min{αs, δ√{β_0^2−s^2}}.

Consequently, define the directional envelope by the scalar extremum

EQ (5): R^2(x,y)=max_{0≤s≤β}{α^2x^2+s^2y^2+2xy min{αs,δ√{β_0^2−s^2}}}.

**Theorem 1 — Directional envelope.** Under (3), the extremum in (5) is given by

EQ (6): R(x,y)=cases{√{α^2x^2+β^2y^2+2δνxy}&δx≤νy;β_0√{x^2+y^2}&νy<δx<αy;αx+δy&δx≥αy}.

For every admissible error, $∥C^He∥≤R(x,y)$, and

EQ (7): R(x,y)≤min{β_0√{x^2+y^2}, αx+βy}.

**Proof.** Compatibility implies $δ≤β$ and $ν≤α$. The two cross term bounds in (4) exchange order at $s=δ$. On $[0,δ]$, the objective is $(αx+sy)^2$, whose maximum is at $s=δ$. For $y>0$, on $[δ,β]$, change variables to $k=√{β_0^2−s^2}∈[ν,α]$. The objective becomes

EQ (8): α^2x^2+β_0^2y^2−y^2k^2+2xyδk.

Its maximizing argument is $δx/y$ projected onto $[ν,α]$. The lower endpoint, interior and upper endpoint give the three branches in (6). Values agree on both branch boundaries. If $y=0$, direct substitution into (5) gives $R=αx$, so division by $y$ is unnecessary. If $α=0$, compatibility forces $β_0=β$ and $R=βy$. If $β=0$, it forces $β_0=α$ and $R=αx$. When $β_0^2=α^2+β^2$, the middle region disappears and $R=αx+βy$. The all zero operator is included by these limits. Inequality (7) follows from the norm cap and the linear cross term cap in (5). Finally, the compression and (4) bound every actual output, proving safety. □

The envelope is nondecreasing in either nonnegative input amplitude. This follows directly from (5), because each objective being maximized has nonnegative coefficients. The same expression, or the constructive characterization below, yields the Lipschitz control used in Section 6.

## 4 QoS envelopes and sharpness
### 4.1 Safe shared budget bound
For $x=|v^He|$, the reverse triangle inequality gives $|b^H(ĥ+e)|≥√p|a−x|$, including $x>a$. The orthogonal budget is at most $√{r^2−x^2}$. Define

EQ (9): Q_G(r)=min_{0≤x≤r}{p(a−x)^2−γ[I_0+R(x,√{r^2−x^2})]^2−γσ^2}.

The previous shared budget bound $Q(r)$ replaces $R$ in (9) with

EQ (10): R_old(x,y)=min{β_0r, αx+βy},   y=√{r^2−x^2}.

**Proposition 1 — Safety and dominance.** For every compatible fixed instance,

EQ (11): Q(r)≤Q_G(r)≤q^*(r).

**Proof.** For any ball feasible error, square the nonnegative lower bound $√p|a−x|$. Apply the interference triangle inequality and the directional envelope with the remaining orthogonal budget. This gives a pointwise lower bound by the objective in (9). Minimizing proves $Q_G≤q^*$. Inequality (7), with $x^2+y^2=r^2$, makes the objective in (9) no smaller than the old objective, proving $Q≤Q_G$. □

Thus, a positive safe numerical lower bound can accept a candidate. A nonpositive value returns unknown and requires a reference decision if a final classification is needed.

### 4.2 Nondegenerate information class sharpness
**Theorem 2 — Sharpness away from the degenerate boundary.** Fix $n≥5$, $m≥3$, $β>0$, $β_0>α≥0$, $p>0$, $γ>0$, $r≥0$, $a,I_0≥0$ and compatible spectral summaries. Over the corresponding nonempty information class,

EQ (12): min_{(b,ĥ,C)∈S_6} q^*(r)=Q_G(r).

The dimensions are sufficient for the construction; no minimal dimension claim is made. The infimum need not be attained by an instance of finite nominal channel norm.

**Proof for positive alpha and radius.** First take $α>0$ and $r>0$. Continuity gives a minimizing $x∈[0,r]$ in (9); put $y=√{r^2−x^2}$. If $y>0$, let $k$ be the projection of $δx/y$ onto $[ν,α]$. If $y=0$, choose $k=α$. In three output coordinates define

EQ (13): d=mat{α;0;0},   ℓ=mat{δk/α;β_0√{1−k^2/α^2};0}.

These vectors satisfy $∥d∥=α$, $∥ℓ∥≤β$, $d^Hℓ=δk$ and $∥xd+yℓ∥=R(x,y)$. For $Z=[d,ℓ]$, the two by two matrix $Z^HZ$ has largest eigenvalue $β_0^2$: the determinant of $β_0^2I−Z^HZ$ is zero, and its diagonal entries are nonnegative. Therefore,

EQ (14): dd^H+ℓℓ^H≼β_0^2I_3.

Choose $0<ε<1$ and put $ℓ_ε=(1−ε)ℓ$. In output space, rank one scaling gives the valid PSD ordering

EQ (15): dd^H+ℓ_εℓ_ε^H≼dd^H+ℓℓ^H≼β_0^2I_3.

No ordering between the original and contracted two by two Gram matrices is required. In fact, their difference generally has a negative determinant when the cross term is nonzero.

Let $d_2,ℓ_{ε,2}$ be the first two coordinates, and define a planar interpolant

EQ (16): S(t)=d_2d_2^H+(1−t)ℓ_{ε,2}ℓ_{ε,2}^H+tβ^2I_2.

We claim that $λ_max(S(0))<β_0^2$. Otherwise a unit vector attaining $β_0^2$ for $S(0)$ also attains equality in (14). Since the removed PSD term is $(2ε−ε^2)ℓ_2ℓ_2^H$, that vector is orthogonal to $ℓ_2$. Its eigenvalue for $S(0)$ must then come entirely from $d_2d_2^H$, which is at most $α^2<β_0^2$, a contradiction. At the other endpoint, $λ_max(S(1))=α^2+β^2≥β_0^2$. Continuity therefore gives a $t_ε∈(0,1]$ for which $λ_max(S(t_ε))=β_0^2$. An explicit root or uniqueness is not needed. Using the planar eigenvalue is essential when $β=β_0$: an auxiliary third direction could otherwise make the full maximum eigenvalue constant over an initial interval of $t$.

Complete the orthogonal operator through

EQ (17): T_ε=diag{(1−t_ε)ℓ_{ε,2}ℓ_{ε,2}^H+t_εβ^2I_2, β^2}.

The planar block is bounded by $β^2I_2$, and the third direction equals $β^2$. Hence $∥T_ε∥_2=β^2$ and $∥dd^H+T_ε∥_2=β_0^2$. Moreover,

EQ (18): M_ε=T_ε−ℓ_εℓ_ε^H=diag{t_ε(β^2I_2−ℓ_{ε,2}ℓ_{ε,2}^H),β^2}≻0.

Positive definiteness follows from $t_ε>0$, $β>0$ and $∥ℓ_ε∥<β$. Let $K_ε=M_ε^{1/2}$, which is invertible, and set

EQ (19): B_ε=[ℓ_ε,K_ε],   C_ε^H=[d,B_ε].

Here $B_ε$ has three rows and four columns and full row rank, while $C_ε^H$ has three rows and five columns. With $v=e_1$ and $b=√p e_1$, we obtain $B_εB_ε^H=T_ε$ and the required exact values of $α,β,β_0$.

Set $u_ε=−xd−yℓ_ε$. Its norm is nonzero: $δ>0$, the chosen $ℓ$ is nonzero, its first coordinate is nonnegative, and $x,y$ are not both zero. Define

EQ (20): z_ε=I_0 frac{u_ε}{∥u_ε∥},   B_εĥ_⊥=z_ε−ad,   ĥ_ε=(a,ĥ_⊥)^T.

Full row rank ensures that this linear equation has a solution. Thus the nominal desired amplitude is exactly $a$ and the interference amplitude is exactly $I_0$, including $I_0=0$. For $e_ε=(−x,−y,0,0,0)^T$, the ball norm equals $r$ and $C_ε^He_ε=u_ε$. Nominal and perturbation interference are aligned. Consequently,

EQ (21): q_ε(e_ε)=p(a−x)^2−γ(I_0+∥u_ε∥)^2−γσ^2.

Let $R=R(x,y)$ and $R_ε=∥u_ε∥$. Triangle inequalities imply $0≤R−R_ε≤εy∥ℓ∥≤εyβ$, because the constructed vectors have nonnegative real coordinates. Therefore,

EQ (22): 0≤q_ε(e_ε)−Q_G(r)≤γ(2I_0+2R)εyβ.

Since $Q_G≤q_ε^*≤q_ε(e_ε)$, the compatible sequence proves that the information class infimum is at most $Q_G$. Proposition 1 proves the reverse inequality.

**Zero alpha and zero radius.** If $α=0$, compatibility gives $β_0=β$. Use $C^H=[0,β]$, $ĥ=(a,I_0/β)^T$ and $b=√p e_1$. For $e=(−x,y)^T$, the signal equals $p(a−x)^2$ and interference equals $(I_0+βy)^2$. Since $R=βy$, choosing a minimizing $x$ attains (9) exactly. This handles $y=0$ without normalizing a zero output. If $r=0$, every compatible instance has the same margin $pa^2−γI_0^2−γσ^2$, equal to $Q_G(0)$. Finally, append zero input coordinates and zero output rows to reach any fixed $n≥5,m≥3$; these operations preserve all six summaries and the evaluated error norm. The proof covers the upper spectral endpoint and $β=β_0$ without excluding them. □

### 4.3 Degenerate boundary and explicit attainment
Suppose $β_0=α>0$. Whenever the beam direction attains the total spectral norm, $v$ is a singular direction of the interference operator, and the equality case of the singular-vector characterization forces the interference generated by any orthogonal error component to be orthogonal to the interference generated along the beam. The parallel interference output is then a fixed complex multiple of the desired signal output — they share one phase — whereas in the nondegenerate regime the two phases can be chosen independently, which is exactly what the completion construction of Theorem 2 exploits. Consequently the boundary information class carries an additional phase compatibility that the six scalars do not display, and its optimal envelope must be derived separately. Put $s=v^Hĥ$, $d̂=d/α$ and

EQ (23): J_⊥=√{I_0^2−α^2a^2},   c_∥=p−γα^2.

Boundary compatibility requires $I_0≥αa$ and $0≤β≤α$. If $β=0$, it additionally requires $I_0=αa$, so $J_⊥=0$. The margin decomposes exactly as

EQ (24): q(e)=c_∥|s+t|^2−γ∥(I−d̂d̂^H)C^H(ĥ+w)∥^2−γσ^2.

**Theorem 3 — Boundary sharpness.** On a nonempty boundary information class with $β_0=α>0$ and $β>0$, for fixed dimensions $n,m≥2$, the minimum information class margin is

EQ (25): Q_∂(r)=min_{0≤x≤r}{c_∥F_{c_∥}(a,x)−γ[J_⊥+β√{r^2−x^2}]^2−γσ^2},
EQ (26): F_{c_∥}(a,x)=cases{(a−x)^2&c_∥≥0;(a+x)^2&c_∥<0}.

In particular, $min_{S_6}q^*(r)=Q_∂(r)$ and this value is attained. For $β=0$, the same formula with $J_⊥=0$ reduces to the parallel problem and is attainable in every compatible dimension.

**Proof of safety.** In (24), if $c_∥≥0$, the reverse triangle inequality gives $c_∥|s+t|^2≥c_∥(a−x)^2$. If $c_∥<0$, use $|s+t|≤a+x$ and reverse the inequality upon multiplication. The nominal orthogonal interference has norm $J_⊥$; its perturbation has norm at most $βy$, with $y≤√{r^2−x^2}$. Applying the triangle inequality and minimizing yields $Q_∂≤q^*$.

**Explicit construction.** For $β>0$, take

EQ (27): C^H=mat{α&0;0&β},   ĥ=mat{a;J_⊥/β},   b=mat{√p;0}.

The six summaries are exactly $p,a,I_0,α,β,α$, since $β≤α$ and $α^2a^2+J_⊥^2=I_0^2$. For a minimizing $x$ in (25), put $y=√{r^2−x^2}$ and choose $e=(−x,y)^T$ for $c_∥≥0$, or $e=(x,y)^T$ for $c_∥<0$. The parallel term in (24) then equals $c_∥F_{c_∥}(a,x)$, and the orthogonal interference norm equals $J_⊥+βy$. The error is feasible and its objective equals $Q_∂$. Safety gives equality with $q^*$ for this instance. Zero padding establishes attainment in larger dimensions. If $β=0$, take $C^H=[α,0,…]$, $ĥ=ae_1$ and $b=√p e_1$; the same phase choice attains the parallel formula. At $r=0$, no limiting argument is needed. □

The general safe bound $Q_G$ remains valid on the boundary, while $Q_∂$ incorporates additional phase compatibility forced by that boundary. Comparing them for the same boundary summaries is legitimate. Applying $Q_∂$ to an instance with $β_0>α$ is not justified.

## 5 Boundary behavior and additional channel norm information
Take $p=2$, $a=I_0=α=γ=1$, $β=0.7$, $σ^2=0.05$ and $r=0.2$. The boundary construction (27) has $J_⊥=0$, $c_∥=1$ and $ĥ=(1,0)^T$. Its one dimensional objective is $(1−x)^2−0.49(0.04−x^2)−0.05$, which decreases on $[0,0.2]$. Thus $Q_∂=q^*=0.59$, attained at $x=0.2$. The general safe bound for these boundary summaries is $Q_G=−0.21$.

For comparison, consider a nondegenerate family with $C^H=mat{1&δ&0;0&0&0.7}$, $ĥ=(1,−2/δ,0)^T$ and $b=(√2,0,0)^T$. It preserves the other five summaries while $β_0=√{1+δ^2}$ approaches $α$ from above. For $δ=0.1,0.03,0.01$, the reference and $Q_G$ agree to numerical accuracy at −0.211028, −0.210093 and −0.210010. The nominal channel norms are 20.02, 66.67 and 200.00. The objective at $e=(−0.2,0,0)^T$ is always −0.21. The nominal channels do not approach a finite fixed channel: when nondegenerate families approach the boundary, the class over which the infimum is taken changes, so the optimum jumps — an information-class discontinuity, not a physical one.

TABLE boundary

If the channel norm $H=∥ĥ∥$ is additionally retained, set $z=C^Hĥ$. A necessary phase compatibility condition is

EQ (28): |d^Hz−α^2s|≤α√{β_0^2−α^2}√{H^2−a^2}.

Indeed, $(CC^H−α^2I)v$ is orthogonal to $v$ and its squared norm is at most $α^2(β_0^2−α^2)$. Cauchy–Schwarz with the orthogonal nominal channel component proves (28). This condition excludes some constructions requiring an unbounded channel norm; the complete seven-scalar envelope is left open.

## 6 Numerical certification and reference interface
### 6.1 Interval lower bounds
A sampled minimum is an upper estimate of a continuous minimum and must not be used for acceptance. We set $x=r sinθ$, $y=r cosθ$ and partition $[0,π/2]$ into $N_θ$ cells. For a cell midpoint $θ_m$ and width $Δθ$, use

EQ (29): R_U=min{β_0r, R(r sinθ_m,r cosθ_m)+r(α+β)Δθ/2}.

The Lipschitz bound $|ΔR|≤α|Δx|+β|Δy|$ follows by representing the direction envelope as the supremum of the compatible linear output norms; the completion in Theorem 2 approaches the scalar extremum for each direction. Degenerate spectra have explicit compatible diagonal or collinear constructions. Each fixed norm is Lipschitz with these constants, so their supremum is as well.

For the signal, use $p$ times the squared distance from $a$ to $[r sinθ_lo,r sinθ_hi]$. Subtract $γ(I_0+R_U)^2+γσ^2$ and take the minimum cell lower bound. This produces $Q_{G,LB}≤Q_G$ in real arithmetic. The old envelope uses the same partition and cap. Main comparisons use $N_θ=256$; refinements to 1024 and 4096 are reported separately. The analytical inequality is exact; the reported implementation uses floating-point arithmetic with conservative numerical safeguards rather than formally verified directed rounding.

### 6.2 Engineering threshold and reference solver
For reference evaluation define $A=bb^H−γCC^H$, $g=Aĥ$ and $c_q=ĥ^H Aĥ−γσ^2$. The margin is $c_q+e^HAe+2Re(e^Hg)$. The common engineering threshold is

EQ (30): τ=10^{−8}max{1,|c_q|,2r∥g∥,r^2∥A∥_2}.

The theoretical envelope depends only on six summaries; this experimental threshold and the reference validity checks additionally use full quadratic data. A method directly accepts a candidate only if every user's implemented lower bound exceeds its threshold. Otherwise it returns unknown. Reference statuses are clearly feasible for $q^*>τ$, clearly infeasible for $q^*<−τ$, ambiguous near zero, and invalid if numerical validity checks fail. Invalid references remain explicitly counted and cannot establish a rejection or a safety conclusion.

The reference solves the Hermitian trust region problem by eigenvalue decomposition and a multiplier search, handling PSD interior solutions, singular cases and the indefinite hard case. The global optimality conditions are stationarity, ball feasibility, $A+μI≽0$, $μ≥0$ and complementary slackness [2],[3]. Numerical checks inspect the full stationarity residual, shifted minimum eigenvalue, ball violation and complementarity. A 106 case audit reports a maximum primal–dual gap of $4.97×10^{−14}$. The fixed 200 iteration multiplier bisection used in early timing over resolves the multiplier: stopping the bisection once the bracket width falls below $10^{−12}$ of its magnitude changes no computed value by more than $5.7×10^{−10}$, preserves every validity check, and uses a median of 40 iterations. A $10^{−10}$ stop already shifts some values by up to $2.8×10^{−8}$ and is not used. Timing sections therefore report both the original 200 iteration reference and the validated tolerance stopped reference, and screening savings are quoted against both.

Fallback skips users already certified by the lower bound and stops at the first clearly infeasible user. Each user's eigenvalue factorization is reused across radii. Screening and reference calculations do not read the truth pass cache during timing. These implementation choices affect measured runtime and are not theoretical contributions.

### 6.3 S procedure SDP baseline
As a full information baseline independent of our eigendecomposition reference, we adapt the standard single sphere S procedure LMI to the fixed candidate certification problem: for each user, maximize $t$ subject to $λ≥0$ and

EQ (31): mat{A+λI&g;g^H&c_q−t−λr^2}≽0,

whose optimum equals $q^*(r)$ by the S lemma; the LMI is solved on the real embedding of the Hermitian data with an interior point solver [2],[3],[4]. This is a classical technique, representative of its use in robust RIS design [4], adapted to our certification interface; it is not a reproduction of that literature's system level algorithms and is not claimed as a new method.

The SDP provides a lower bound, but a negative lower bound proves nothing about infeasibility, and a solver status of optimal is not a certificate. Our decision rule is two sided. The lower bound $L$ is the residual corrected dual feasible value: the shift $(t,λ)→(t−η(1+r^2),λ+η)$ adds blockdiag$(ηI,η)$ to the LMI, so Weyl's inequality restores positive semidefiniteness up to eigensolver rounding, with an empirical rounding allowance and a numerical recheck; this is an implementation detail, not a formal error bound. Infeasibility is certified only by an upper bound: from the corrected multiplier we recover the Lagrangian minimizer $(A+λ'I)e=−g$, project it onto the ball if $∥e∥>r$, and re evaluate the true objective $U=q(e_w)$. Since $e_w$ is ball feasible, $U$ is a valid upper bound of $q^*$, and $U<−τ$ certifies infeasibility; no optimality of $e_w$ is claimed. All remaining cases return unresolved and are never treated as rejections.

## 7 Numerical verification
### 7.1 Analytical implementation checks
With seed 20261006, 30 random compatible spectral tuples are checked against a 4097 point scalar grid, local refinement and direct evaluation of $s=0,δ,β$. The maximum $|R^2−R_ref^2|$ is $1.11×10^{−16}$. The grid only underestimate reaches $4.12×10^{−5}$ and is reported separately. Random inequality checks show no violation of (7). The three fixed direction cases cover the three branches; hand checks include zero spectral gains and zero input amplitudes.

The direction completion matches the three spectral summaries to machine precision (nine formal cases, and 251 additional repaired constructions covering asymmetric tuples, the upper spectral endpoint, $β=β_0$ and direction endpoints; maximum six scalar discrepancy $2.89×10^{−15}$, objective identity residual $1.25×10^{−14}$). The nominal channel construction matches all six summaries, and the deficit $q(e_ε)−Q_G$ obeys the proven bound (22) and falls by roughly one decade per decade of $ε$. A 36 case degenerate battery covers $α=0$, $β=0$, $r=0$, $I_0=0$, both spectral endpoints, both branch connection points, the trust region hard case equality, repeated minimum eigenvalues and radius equality. For the deterministic flip, the separable bound is −0.0425, the old approximate minimum is −0.024580 and the new approximate minimum is 0.034624; their interval lower bounds at 256 cells are −0.027068 and +0.027262, and the new lower bound at 65536 cells is 0.0345956. The positive certified lower bound, rather than the approximate minimum, establishes acceptance. These are consistency checks of the analytical proofs, not substitutes for them.

### 7.2 Controlled threshold-crossing examples
The scan fixes the flip example's other summaries and varies $β_0$ over 21 equally spaced points from 1.02 to $√2$. Each point has an explicit matrix and nominal channel; the channel norm ranges from 2.22 to 9.05. This is an instance family, not an experiment holding a physical channel fixed. The new envelope exceeds the old one at 20 of 21 points, with maximum gain 0.0749 at $β_0=1.1777$; they agree at $β_0=√2$, as predicted by the spectral endpoint. Six points have $Q_{G,LB}>τ$ while the old lower bound is nonpositive, with $β_0∈[1.0791,1.1777]$. The minimum reference to new interval lower bound gap is 0.00663325.

FIGURE spectral

### 7.3 Fixed-candidate certification study
The candidate model has two base stations with four antennas each, four users, and a 16 element RIS with two bit phases. Each seed generates six phase states (all ones plus five random) and two beam types, MRT and RZF (regularization 0.1), giving 12 candidates; per base station power is 1 with 0.5 per stream. The effective link is $ĥ_{j,u}=s_{j,u}[d_{j,u}+G_j^Hdiag(θ)f_u]$ with same-cell scale 1 and cross-cell scale 0.5; entries of $G_j,d_{j,u}$ have circular complex normal variance 1/4 and entries of $f_u$ variance 1/16. Channels and beams are stacked into dimension eight; noise power is 0.05. These are normalized dimensionless links, not a fitted propagation model.

The original pool (seeds 11–55, $γ=0.5$, $ρ∈\{0.01,0.03,0.06\}$, 720 user problems, 180 units) is an easy baseline: all four methods accept all units, with no false acceptance or fallback. The exploratory grid adds seeds 66–110 and evaluates $γ∈\{0.5,1,2\}$ and $ρ∈\{0.06,0.15,0.30\}$ — 1080 units and 4320 user problems, every cell reported, no selection by the new bound's gain. All references are valid; 472 units are clearly feasible and 608 clearly infeasible. There is no observed false acceptance or lower bound violation. At 256 cells the new bound directly accepts 390 units versus 380 for the old bound (feasible coverage 82.63% versus 80.51%; original seeds contribute six additional accepts, added seeds four).

TABLE sensitivity

The quadratic bound uses

EQ (32): L_quad(r)=c_q−2r∥g∥−γβ_0^2r^2.

Unlike the six-scalar envelopes, it uses full quadratic geometry. In the grid it exceeds the implemented new interval bound on 2903 of 4320 user evaluations; that count does not compare it with the exact continuous envelope. A separate witness has $L_quad=0.198195750531$ and a new envelope objective value of 0.196451952840 at a feasible one dimensional argument, so the continuous minimum is strictly below $L_quad$ there; conversely, the constructed spectral family has $Q_G>L_quad$. The bounds therefore have no universal dominance relation; their maximum is safe and yields the highest direct-certification coverage among tested bounds. The pools overlap on 240 user problems; their combined 5040 evaluation records represent 4800 distinct user problems, and multiple radii or thresholds on one candidate are correlated, not independent samples.

### 7.4 Independent validation
Before computing any result on new seeds, we wrote a fixed configuration file recording the unused seeds 211–215, the grid, the candidate rules, the acceptance rule, the reference validity criteria, the timing protocol and SHA-256 hashes of the code files. The batch adds 540 units and 2160 user problems with zero overlap; the union of all pools contains 6960 distinct user problems. All references are valid and no false acceptance or bound violation occurs. Of 540 units, 159 are clearly feasible; direct accepts are D 111, previous 125, new 126 and $L_quad$ 138. The new bound adds one unit over the previous (coverage 78.62% to 79.24%); combined with the exploratory grid (+10 units, +2.12 points), the certification gain is positive but small and varies between batches. Reference calls fall from 511 (previous pipeline) to 496 (new pipeline) against 1319 for the full reference.

TABLE validation

Table 5 separates interval resolution effects: in the exploratory grid the previous envelope gains one accept at 4096 cells while the new envelope is unchanged; in the validation batch the counts are unchanged at every resolution. Decomposing the gap on the 2160 validation rows, the numerical interval loss (upper estimate minus the 256 cell certified bound) has median $1.22×10^{−3}$ and maximum $2.53×10{−2}$, whereas the theoretical envelope gap (reference minus upper estimate) has median $6.08×10^{−2}$ and maximum 2.68: certification outcomes are limited by the theoretical envelope rather than by the interval discretization.

TABLE resolution

A second validation flavor uses a pool with heterogeneous, distance-dependent scaling (seeds 301–305, same 540-unit structure): two base stations and a 16 element RIS with fixed geometry, path loss $PL(d)=10+30log_{10}(d/10)$ dB on the user links, a 12 dB BS–RIS link, Rician factors $K=3$ on BS–RIS and RIS–UE and Rayleigh BS–UE links, the same candidate rules, and relative balls $r=ρ∥ĥ∥$. The relative radius is interpreted through the worst-case NMSE mapping $ρ=0.06/0.15/0.30↔−24.4/−16.5/−10.5$ dB ($ρ^2$ equals the NMSE when the bound is set from $E∥e∥^2/E∥h∥^2$); the setup is an illustrative heterogeneous scaling, not a calibrated link budget. All references are valid, no bound violation occurs, and the main conclusions reappear: of 540 units (211 feasible), the previous envelope accepts 156, the new envelope 163 (77.25% versus 73.93% of feasible, seven flips), $L_quad$ 186, and the joint maximum 190 (90.05%).

TABLE fa

### 7.5 Complementary certificates and staged screening
The quadratic bound and the six-scalar envelope use different information, and neither dominates. Taking the per user maximum and accepting a candidate only when every user passes the common threshold yields a strictly larger direct acceptance set: 424 of 472 feasible units (89.83%) in the exploratory grid versus 408 (86.44%) and 390 (82.63%); 139 of 159 (87.42%) versus 138 and 126 in the validation batch; and 190 of 211 (90.05%) versus 186 and 163 in the realistic pool. The gain over the quadratic bound comes entirely from units certified only by the six-scalar envelope (16, 1 and 4 units respectively); no unit is accepted by the joint rule while both individual rules fail. The separable bound is dominated by the envelope at every user problem in all pools and is dropped from the flow.

TABLE joint

Staged screening pipelines complete the same final decision task: bounds are evaluated only on users not yet certified, a failed bound is never a rejection, the exact reference ends a candidate at the first clearly infeasible user, and each user's factorization is cached. With one warm up and three interleaved repetitions, the ordering quadratic bound, then envelope, then reference — selected by a fixed rule on the exploratory cohort and frozen, with code hashes and unused seeds recorded before any new computation — attains a median paired speedup of 2.12 over the pure 200 iteration reference (25 of 30 pairs) on the exploratory grid and 2.01 (11 of 15) on five fresh confirmation seeds (221–225), reducing total reference calls by 75% and 72% respectively. The quadratic bound alone is the fastest single bound pipeline; the envelope stage contributes the additional certified units, not speed. All pipelines produce identical final decisions in every run. An adaptive interval refinement (64→256→1024→4096 cells, stopping when a rounding guarded feasible point value shows the threshold cannot be crossed) reproduces fixed 4096 cell acceptance exactly with 3.9× fewer interval cells but is slower in wall clock at this scale and is not used.

A controlled scaling study isolates dimension from load: four users, per base station power 1 (per stream 0.5), $γ=1$, $ρ=0.15$, identical candidate rules; only the antenna count varies (stacked dimension 4 to 32). Interval bound evaluation is dimension independent (0.04–0.05 ms per user problem) because it operates on the summaries, while the exact reference grows from 0.36 to 0.45 ms; accordingly, screening advantage over the pure tolerance stopped reference grows from near parity at dimension 4 to about 2.8× at dimensions 16–32. Acceptance counts are identical across pipelines at every dimension; they are not compared across dimensions.

TABLE fb

### 7.6 Independent exact cross-validation
A preregistered 74 case check set and the full 2160 problem validation batch are solved by the corrected SDP pipeline of Section 6.3. On the batch, the certified lower bound never exceeds the reference (maximum excess $−1.2×10^{−13}$), the witness upper bound never falls below it (maximum deficit $1.3×10^{−12}$), the ordering $L≤U$ holds with no violation, there is no false acceptance and no false rejection, all witnesses lie inside the ball (956 of 2160 required projection), and every user problem is decided (1247 feasible, 913 infeasible, zero unresolved) with witness quality $U−q^*$ of median $6.7×10^{−6}$ and maximum $1.6×10^{−4}$. The maximum deviation of the raw SDP optimum from the reference is $5.6×10^{−8}$. Two independent exact methods thus cross validate each other and every certification statement above.

TABLE sdp

### 7.7 Timing and implementation sensitivity
Table 10 reports the original protocol against the 200 iteration reference, stratified by threshold; Table 11 reports the same final decision task against the validated tolerance stopped reference. With the tolerance stop, the pure reference drops from 0.150 s to 0.041 s per seed and threshold, so the apparent advantage of screening pipelines narrows: the best configuration (new bound with tolerance stopped fallback, 0.032 s) is about 20% below the tolerance stopped pure reference, while the screening pipeline with the original fallback (0.075 s) is slower. The new pipeline is faster than the previous envelope pipeline in 6 of 15 validation pairs and 11–13 of 30 exploratory pairs across runs; no stable runtime advantage over the previous envelope is observed, and part of the apparent saving in early measurements traced to multiplier search over iteration in the reference solver rather than to screening. The results do not establish complexity advantages or superiority over an optimized solver. The generic SDP pipeline is a further 4–13× slower than the tolerance stopped reference; its value is cross validation, and the slowness of the generic conic path is not a statement about full information methods.

TABLE timing
TABLE toltiming

## 8 Interface to configuration lifecycle management
The certificate addresses whether a fixed candidate can be deployed at a prescribed uncertainty radius. Once a candidate is selected, a separate complete configuration reference can determine its robust radius for subsequent online retain or refresh management. Candidate generation, deployment costs and channel drift remain outside the present problem. A positive certificate at one radius is not a maximum radius certificate.

## 9 Discussion
The sharpness statement is intentionally strong and narrow: it holds for the unrestricted information class, whose extremal constructions may require large nominal channel norms that a physical RIS channel may not admit. Restricting the class — for example by retaining the channel norm, for which (28) gives a necessary compatibility condition — can only shrink the admissible set, and the six-scalar envelope need not remain sharp; the complete seven-scalar envelope is open. Empirically, the certificate gain over the previous six-scalar envelope is positive in every batch yet small and batch-dependent, and no stable runtime advantage over the previous envelope is observed; measured speedups of screening pipelines are tied to the tested implementation. The certificate is complementary to, not a replacement for, full-information methods: within the certification hierarchy it is the strongest guarantee available at its information level, and combining levels by maximization is safe.

## 10 Conclusion
Mathematically, joint spectral compatibility of three gains of one interference operator yields an explicit three-branch directional extremum, a safe shared-budget QoS envelope dominating the separately-used-gain envelope, equality with the information-class infimum away from the degenerate boundary, and a separately derived, explicitly attained phase-coupled envelope at that boundary. Practically, the envelope is a fixed-size certificate whose empirical gain over the previous six-scalar bound is strict but modest, and which is complementary to a cheap full-geometry quadratic bound — the two have no universal dominance relation and combine safely by maximization, giving the highest direct-certification coverage among the tested bounds. In scope, the certificate does not replace exact full-information solvers, does not optimize candidate generation, and its sharpness concerns the stated unrestricted information class; natural extensions include additional retained information, physically restricted information classes, validated numerical certificates, and lifecycle management of stored configuration certificates.

## References
[1] G. Zhou, C. Pan, H. Ren, K. Wang, and A. Nallanathan, "A Framework of Robust Transmission Design for IRS-Aided MISO Communications With Imperfect Cascaded Channels," IEEE Transactions on Signal Processing, vol. 68, pp. 5092–5106, 2020. DOI: 10.1109/TSP.2020.3019666.

[2] F. Rendl and H. Wolkowicz, "A Semidefinite Framework for Trust Region Subproblems with Applications to Large Scale Minimization," Mathematical Programming, vol. 77, pp. 273–299, 1997.

[3] C. Fortin and H. Wolkowicz, "The Trust Region Subproblem and Semidefinite Programming," Optimization Methods and Software, vol. 19, no. 1, pp. 41–67, 2004. DOI: 10.1080/10556780410001647186.

[4] A. Zheng, W. Ni, W. Wang, H. Tian, and C. Yuen, "Energy-Efficient Robust Beamforming for Multi-Functional RIS-Aided Wireless Communication Under Imperfect CSI," IEEE Transactions on Communications, vol. 73, no. 8, pp. 6306–6321, 2025. DOI: 10.1109/TCOMM.2024.3516491.

[5] T. Bottazzi and C. Conde, "A Generalized Buzano Inequality," Filomat, vol. 37, 2023 (and M. L. Buzano, Rend. Sem. Mat. Univ. Politech. Torino, vol. 31, 1974).

[6] R. Grone, C. R. Johnson, E. M. Sá, and H. Wolkowicz, "Positive Definite Completions of Partial Hermitian Matrices," Linear Algebra and its Applications, vol. 58, pp. 109–124, 1984.

[7] S. Parrott, "On a Quotient Operator and Related Enumeration Theorems," and the operator-ball completion literature; see also "A strong Parrott theorem" for a modern strengthening of the prescribed-column contraction completion.

[8] C. A. Micchelli and T. J. Rivlin, "A Survey of Optimal Recovery," in Optimal Estimation in Approximation Theory, Plenum Press, 1977, pp. 1–54; and A. A. Melkman and C. A. Micchelli, "Optimal Estimation of Linear Operators in Hilbert Spaces from Inaccurate Data," SIAM Journal on Numerical Analysis, vol. 16, 1979.

[9] H. Sun et al., "Robust Beamforming Design in a NOMA Cognitive Radio Network," IEEE 2018 (representative of parallel/orthogonal error decomposition in bounded-error robust beamforming); exact bibliographic details to be finalized.

## Reproducibility
The supplied code reproduces every number. `run_experiments.py` (original pools), `audit_and_stress.py` (sensitivity study), `run_v41.py` (preregistered validation batch and theory checks), `run_optimization.py` (joint certification, staged screening, adaptive refinement, scaling), `run_sdp_v2.py` (corrected SDP baseline) and `experiments_v6.py` (realistic pool, NMSE mapping and controlled scaling) write into separate immutable result directories and exit nonzero on any failed check. An earlier SDP implementation with an invalid negative-lower-bound rejection rule is retained as history with correction headers and is not cited. Numerical minima are diagnostics; interval lower bounds are used for acceptance.
