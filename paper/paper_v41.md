# Spectral Compatible Shared Error Budget Envelopes for Robust QoS Certification
Fixed beamforming and RIS candidates under a six scalar information model
Version 4.1 • English revision • 6 October 2026

## Abstract
We study the certification of a fixed downlink configuration under a Euclidean channel error ball. The certificate retains six scalar summaries of the desired beam, nominal channel and interference operator. Joint use of the desired direction gain, orthogonal subspace gain and total spectral norm yields a continuous three branch envelope for interference perturbations. Its associated shared error budget QoS bound is safe and dominates the previous envelope based on separate norm inequalities. For a nondegenerate information class, we prove that the bound equals the infimum of the true worst case QoS margin over all instances with the same six summaries. The proof uses a positive definite operator completion and handles zero desired direction interference separately. At the degenerate boundary where the desired direction attains the total interference norm, we give an explicit construction attaining a different envelope that respects signal and interference phase coupling. The boundary envelope can be 0.59 while the general safe bound is −0.21 for the same boundary summaries. Numerical checks reproduce the analytical envelopes and certify six flips in a 21 point constructed spectral scan. In an exploratory fixed parameter candidate study, the implemented bound certifies 390 of 472 feasible units, versus 380 for the previous bound, with no observed false acceptance. The resulting reduction in reference calls does not establish a stable runtime advantage over the previous envelope. The optimality statements concern the unrestricted information class, and the numerical results are limited to the tested candidates and reference implementation.

Keywords: robust QoS certification; shared error budget; spectral compatibility; information class optimality; RIS; fixed configuration.

## 1 Introduction
For a candidate downlink configuration with fixed beamforming vectors and discrete reconfigurable intelligent surface (RIS) phases, deployment requires verification of all user QoS constraints under channel uncertainty. Robust transmission design with imperfect RIS channels is an established problem [1]. The fixed configuration verification considered here is a ball constrained quadratic problem. Trust region methods and the single ball S procedure provide exact optimization formulations [2,3]. These tools are reference methods in this paper rather than contributions.

Our question concerns the information available to a certificate. Suppose that the complete matrix geometry is compressed into six scalar summaries. What uniform guarantee can a certificate obtain from these summaries? A separable bound uses the full error radius both to reduce the desired signal and to increase interference. Shared budget bounds improve this construction by decomposing the same error into its desired direction and orthogonal components. However, a bound that separately uses three interference gains still does not fully exploit their origin in one operator.

We solve the corresponding direction extremum under the joint spectral constraints. The resulting envelope gives a stronger safe QoS bound. We then characterize its optimality over a fixed six scalar information class, with separate treatment of the degenerate boundary. The contributions are the explicit extremum, the completion proof of information class sharpness and the boundary construction. Orthogonal decompositions, positive semidefinite (PSD) constraints, matrix factorization and trust region optimization are established tools used in these arguments.

The numerical study asks whether the formulas are implemented correctly, whether the new envelope crosses certification thresholds that the previous envelope misses, and how a screening and reference fallback pipeline behaves on fixed candidates. It does not address candidate generation optimality, channel evolution or online configuration refresh.

## 2 Fixed configuration model and information class
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

For fixed dimensions $n,m$, let $S_6$ denote the nonempty set of triples $(b,ĥ,C)$ having the prescribed summaries. We assume $p>0$, $a,I_0≥0$. The class does not constrain the nominal channel norm, interference rank, base station support pattern or the set of physically reachable RIS channels. Thus, information class sharpness does not imply optimality over a specific wireless deployment.

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

EQ (12): inf_{(b,ĥ,C)∈S_6} q^*(r)=Q_G(r).

The dimensions are sufficient for the construction; no minimal dimension claim is made.

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

Since $Q_G≤q_ε^*≤q_ε(e_ε)$, the compatible sequence proves that the information class infimum is at most $Q_G$. Proposition 1 proves the reverse inequality. The infimum need not be attained by a finite norm nominal channel.

**Zero alpha and zero radius.** If $α=0$, compatibility gives $β_0=β$. Use $C^H=[0,β]$, $ĥ=(a,I_0/β)^T$ and $b=√p e_1$. For $e=(−x,y)^T$, the signal equals $p(a−x)^2$ and interference equals $(I_0+βy)^2$. Since $R=βy$, choosing a minimizing $x$ attains (9) exactly. This handles $y=0$ without normalizing a zero output. If $r=0$, every compatible instance has the same margin $pa^2−γI_0^2−γσ^2$, equal to $Q_G(0)$. Finally, append zero input coordinates and zero output rows to reach any fixed $n≥5,m≥3$; these operations preserve all six summaries and the evaluated error norm. The proof covers the upper spectral endpoint and $β=β_0$ without excluding them. □

### 4.3 Degenerate boundary and explicit attainment
Suppose $β_0=α>0$. Since $v$ attains the operator norm, $CC^Hv=α^2v$. Thus $d^HC^HP=0$: the parallel and orthogonal interference outputs are orthogonal. Put $s=v^Hĥ$, $d̂=d/α$ and

EQ (23): J_⊥=√{I_0^2−α^2a^2},   c_∥=p−γα^2.

Boundary compatibility requires $I_0≥αa$ and $0≤β≤α$. If $β=0$, it additionally requires $I_0=αa$, so $J_⊥=0$. The margin decomposes exactly as

EQ (24): q(e)=c_∥|s+t|^2−γ∥(I−d̂d̂^H)C^H(ĥ+w)∥^2−γσ^2.

**Theorem 3 — Boundary sharpness.** On a nonempty boundary information class with $β_0=α>0$ and $β>0$, for fixed dimensions $n,m≥2$, the minimum information class margin is

EQ (25): Q_∂(r)=min_{0≤x≤r}{c_∥F_{c_∥}(a,x)−γ[J_⊥+β√{r^2−x^2}]^2−γσ^2},
EQ (26): F_{c_∥}(a,x)=cases{(a−x)^2&c_∥≥0;(a+x)^2&c_∥<0}.

In particular, $inf_{S_6}q^*(r)=Q_∂(r)$ and this value is attained. For $β=0$, the same formula with $J_⊥=0$ reduces to the parallel problem and is attainable in every compatible dimension.

**Proof of safety.** In (24), if $c_∥≥0$, the reverse triangle inequality gives $c_∥|s+t|^2≥c_∥(a−x)^2$. If $c_∥<0$, use $|s+t|≤a+x$ and reverse the inequality upon multiplication. The nominal orthogonal interference has norm $J_⊥$; its perturbation has norm at most $βy$, with $y≤√{r^2−x^2}$. Applying the triangle inequality and minimizing yields $Q_∂≤q^*$.

**Explicit construction.** For $β>0$, take

EQ (27): C^H=mat{α&0;0&β},   ĥ=mat{a;J_⊥/β},   b=mat{√p;0}.

The six summaries are exactly $p,a,I_0,α,β,α$, since $β≤α$ and $α^2a^2+J_⊥^2=I_0^2$. For a minimizing $x$ in (25), put $y=√{r^2−x^2}$ and choose $e=(−x,y)^T$ for $c_∥≥0$, or $e=(x,y)^T$ for $c_∥<0$. The parallel term in (24) then equals $c_∥F_{c_∥}(a,x)$, and the orthogonal interference norm equals $J_⊥+βy$. The error is feasible and its objective equals $Q_∂$. Safety gives equality with $q^*$ for this instance. Zero padding establishes attainment in larger dimensions. If $β=0$, take $C^H=[α,0,…]$, $ĥ=ae_1$ and $b=√p e_1$; the same phase choice attains the parallel formula. At $r=0$, no limiting argument is needed. □

The general safe bound $Q_G$ remains valid on the boundary, while $Q_∂$ incorporates additional phase compatibility forced by that boundary. Comparing them for the same boundary summaries is legitimate. Applying $Q_∂$ to an instance with $β_0>α$ is not justified.

## 5 Boundary behavior and additional channel norm information
Take $p=2$, $a=I_0=α=γ=1$, $β=0.7$, $σ^2=0.05$ and $r=0.2$. The boundary construction (27) has $J_⊥=0$, $c_∥=1$ and $ĥ=(1,0)^T$. Its one dimensional objective is $(1−x)^2−0.49(0.04−x^2)−0.05$, which decreases on $[0,0.2]$. Thus $Q_∂=q^*=0.59$, attained at $x=0.2$. The general safe bound for these boundary summaries is $Q_G=−0.21$.

For comparison, consider a nondegenerate family with $C^H=mat{1&δ&0;0&0&0.7}$, $ĥ=(1,−2/δ,0)^T$ and $b=(√2,0,0)^T$. It preserves the other five summaries while $β_0=√{1+δ^2}$ approaches $α$ from above. For $δ=0.1,0.03,0.01$, the reference and $Q_G$ agree to numerical accuracy at −0.211028, −0.210093 and −0.210010. The nominal channel norms are 20.02, 66.67 and 200.00. The objective at $e=(−0.2,0,0)^T$ is always −0.21.

This also establishes the discontinuity of the information class optimum for this parameter sequence without relying on numerical agreement: the safe lower bound tends to −0.21, while the displayed feasible error gives an upper bound of −0.21 on the family margin. At the boundary, Theorem 3 gives 0.59. The nominal channels do not approach a finite fixed channel. Hence this is a change in the optimum over information classes, not a discontinuity of one physical system's QoS.

Table 1 reports four boundary checks. They vary $p$ and $I_0$, while $a=α=γ=1$, $β=0.7$, $σ^2=0.05$ and $r=0.2$ remain fixed. In particular, the cases with $I_0=1.3$ are not cases with $a=1.3$.

TABLE boundary

If the channel norm $H=∥ĥ∥$ is additionally retained, set $z=C^Hĥ$. A necessary phase compatibility condition is

EQ (28): |d^Hz−α^2s|≤α√{β_0^2−α^2}√{H^2−a^2}.

Indeed, $(CC^H−α^2I)v$ is orthogonal to $v$ and its squared norm is at most $α^2(β_0^2−α^2)$. Cauchy–Schwarz with the orthogonal nominal channel component proves (28). This necessary condition excludes some constructions requiring an unbounded channel norm. An optimal seven scalar envelope is not established here.

## 6 Numerical certification and reference interface
### 6.1 Interval lower bounds
A sampled minimum is an upper estimate of a continuous minimum and must not be used for acceptance. We set $x=r sinθ$, $y=r cosθ$ and partition $[0,π/2]$ into $N_θ$ cells. For a cell midpoint $θ_m$ and width $Δθ$, use

EQ (29): R_U=min{β_0r, R(r sinθ_m,r cosθ_m)+r(α+β)Δθ/2}.

The Lipschitz bound $|ΔR|≤α|Δx|+β|Δy|$ follows by representing the direction envelope as the supremum of the compatible linear output norms; the completion in Theorem 2 approaches the scalar extremum for each direction. Degenerate spectra have explicit compatible diagonal or collinear constructions. Each fixed norm is Lipschitz with these constants, so their supremum is as well.

For the signal, use $p$ times the squared distance from $a$ to $[r sinθ_lo,r sinθ_hi]$. Subtract $γ(I_0+R_U)^2+γσ^2$ and take the minimum cell lower bound. This produces $Q_{G,LB}≤Q_G$ in real arithmetic. The old envelope uses the same partition and cap. Main comparisons use $N_θ=256$; refinements to 1024 and 4096 are reported separately.

### 6.2 Engineering threshold and reference solver
For reference evaluation define $A=bb^H−γCC^H$, $g=Aĥ$ and $c_q=ĥ^HAĥ−γσ^2$. The margin is $c_q+e^HAe+2Re(e^Hg)$. The common engineering threshold is

EQ (30): τ=10^{−8}max{1,|c_q|,2r∥g∥,r^2∥A∥_2}.

The theoretical envelope depends only on six summaries; this experimental threshold and the reference validity checks additionally use full quadratic data. A method directly accepts a candidate only if every user's implemented lower bound exceeds its threshold. Otherwise it returns unknown. Reference statuses are clearly feasible for $q^*>τ$, clearly infeasible for $q^*<−τ$, ambiguous near zero, and invalid if numerical validity checks fail. Invalid references remain explicitly counted and cannot establish a rejection or a safety conclusion.

The reference solves the Hermitian trust region problem by eigenvalue decomposition and a multiplier search, handling PSD interior solutions, singular cases and the indefinite hard case. The global optimality conditions are stationarity, ball feasibility, $A+μI≽0$, $μ≥0$ and complementary slackness [2,3]. Numerical checks inspect the full stationarity residual, shifted minimum eigenvalue, ball violation and complementarity. An additional 106 case audit reports a maximum primal–dual gap of $4.97×10^{−14}$. Ordinary double precision and engineering thresholds do not provide a machine checked formal certificate.

Fallback skips users already certified by the lower bound and stops at the first clearly infeasible user. Each user's eigenvalue factorization is reused across radii. The exact baseline has the same reuse and early stopping. Screening and reference calculations do not read the truth pass cache during timing. These implementation choices affect measured runtime and are not theoretical contributions.

## 7 Numerical verification
### 7.1 Directional and completion checks
With seed 20261006, 30 random compatible spectral tuples are checked against a 4097 point scalar grid, local refinement and direct evaluation of $s=0,δ,β$. The maximum $|R^2−R_ref^2|$ is $1.11×10^{−16}$. The grid only underestimate reaches $4.12×10^{−5}$ and is reported separately. Random inequality checks show no violation of (7). The three fixed direction cases cover the three branches; hand checks include zero spectral gains and zero input amplitudes.

The A2 direction completion has nine cases: three direction points and three contraction levels. It matches the three spectral summaries to machine precision. The separate A3 nominal channel construction matches all six summaries. In the flip example with $p=4$, $a=1.3$, $I_0=0.5$, $α=β=1$, $β_0=1.1$, $γ=1$, $σ^2=1.5$ and $r=0.5$, the approximate new minimum and reference are 0.0346243. The objective at the constructed error exceeds the envelope by approximately $5.1×10^{−4}$, $5.2×10^{−5}$ and $5.2×10^{−6}$ for the three contraction levels.

For the deterministic flip, the separable bound is −0.0425, the old approximate minimum is −0.024580 and the new approximate minimum is 0.034624. Their interval lower bounds at 256 cells are −0.027068 and +0.027262. The new lower bound at 65536 cells is 0.0345956. The positive certified lower bound, rather than the numerical approximate minimum, establishes acceptance.

The repaired V4.1 completion is also checked on 251 constructions, including asymmetric spectral tuples, the upper spectral endpoint, $β=β_0$, and direction endpoints. The maximum six scalar discrepancy is $2.89×10^{−15}$ and the objective identity residual is $1.25×10^{−14}$. Sixty additional explicit boundary instances cover both signs of $c_∥$ and match the constructed objective and reference within $10^{−8}$. These are consistency checks of the analytical proofs, not substitutes for them.

### 7.2 Constructed spectral scan
The scan fixes the flip example's other summaries and varies $β_0$ over 21 equally spaced points from 1.02 to $√2$. Each point has an explicit matrix and nominal channel; the channel norm ranges from 2.22 to 9.05. This is an instance family, not an experiment holding a physical channel fixed.

The new approximate envelope exceeds the old one at 20 of 21 points. The maximum gain is 0.0749 at $β_0=1.1777$. They agree at $β_0=√2$, as predicted by the spectral endpoint. Six points have $Q_{G,LB}>τ$ while the old lower bound is nonpositive, with $β_0∈[1.0791,1.1777]$. The minimum reference to new interval lower bound gap is 0.00663325. Approximate minima differ from the reference by tiny negative roundoff at some points, which is not an acceptance criterion.

FIGURE spectral

### 7.3 Fixed RIS candidates and easy regime
The candidate model has two base stations with four antennas each, four users, and a 16 element RIS with two bit phases. Each seed generates six phase states and two beam types, MRT and RZF, giving 12 candidates. The first phase state is all ones and five states are sampled from $\{1,i,−1,−i\}$. Each base station serves two users with power 0.5 per stream. RZF uses regularization 0.1 before per stream normalization.

The effective link is $ĥ_{j,u}=s_{j,u}[d_{j,u}+G_j^Hdiag(θ)f_u]$, with same cell scale 1 and cross cell scale 0.5. Entries of $G_j,d_{j,u}$ have circular complex normal variance 1/4; entries of $f_u$ have variance 1/16. Real and imaginary parts each have half the stated variance. Channels and beams are stacked into dimension eight. Noise power is 0.05. These are normalized dimensionless links, not a fitted propagation or hardware model.

The original seeds are 11, 22, 33, 44 and 55, with $γ=0.5$ and $r_u=ρ∥ĥ_u∥$, $ρ∈\{0.01,0.03,0.06\}$. The pool contains 720 user ball evaluations and 180 candidate–radius units. All references are valid and all four methods accept all units directly, with no false acceptance or fallback. This is an easy baseline and cannot demonstrate differences between certificates.

### 7.4 Exploratory sensitivity grid
The candidate generation is unchanged. We add seeds 66, 77, 88, 99 and 110 and evaluate $γ∈\{0.5,1,2\}$ and $ρ∈\{0.06,0.15,0.30\}$. The fixed grid is an exploratory response to saturation of the original pool; it is not a preregistered validation. Every grid cell is reported, without selecting candidates by the new bound's gain. There are 120 candidates, 1080 candidate–radius–threshold units and 4320 user ball evaluations. Larger radii are stress settings, not calibrated statistical CSI coverage.

TABLE sensitivity

All references are valid; 472 units are clearly feasible and 608 clearly infeasible, with no ambiguous units. There is no observed false acceptance or lower bound violation. At 256 cells, the new bound directly accepts 390 units, versus 380 for the old bound. Feasible coverage increases from 80.51% to 82.63%, or 2.12 percentage points; the total direct acceptance increase is 0.93 percentage points. Original seeds contribute six additional accepts and added seeds contribute four.

The comparison with the quadratic bound uses

EQ (31): L_quad(r)=c_q−2r∥g∥−γβ_0^2r^2.

Unlike the six scalar envelopes, it uses full quadratic geometry. In the sensitivity study it exceeds the implemented new interval bound on 2903 of 4320 user evaluations. That count does not compare it with the exact continuous envelope. A separate original pool witness has $L_quad=0.198195750531$ and a new envelope objective value of 0.196451952840 at a feasible one dimensional argument; thus the continuous minimum is strictly below $L_quad$ there. Conversely, the constructed spectral family has $Q_G>L_quad$. The bounds therefore have no universal dominance relation. Their maximum can be used as a safe implementation combination, without changing the six scalar optimality statement.

The original and sensitivity pools overlap on 240 user problems at $γ=0.5$, $ρ=0.06$ for the original seeds. Their combined 5040 evaluation records represent 4800 distinct user problems. Multiple radii and thresholds on the same candidate are correlated and do not constitute independent samples.

### 7.5 Resolution and end to end timing
Table 3 separates interval resolution effects from the main 256 cell comparison. Increasing resolution to 4096 lets the old envelope certify one additional unit; the new envelope still certifies 390. The claim of ten additional accepts is therefore tied to the implemented 256 cell setting, not a claim that every such unit is a continuous envelope sign flip.

TABLE resolution

Timing measures the same final candidate decision task for all pipelines. It includes summary or quadratic data construction, all prescribed radii, caching, reference validity checks and early stopping, and excludes candidate generation. Each seed and threshold has one warmup followed by three interleaved repetitions. Table 4 shows the median across ten per seed medians. Each seed contributes 36 candidate units per threshold; parentheses are total reference calls across ten seeds.

TABLE timing

The new pipeline reduces total reference calls from 850 for the old bound to 821, or 3.41%. It is faster than the old pipeline in only 11 of the 30 seed–threshold pairs. Hence, the data show a modest certification improvement without a stable additional runtime advantage. Screening with fallback saves time relative to this serial reference implementation, with nonzero fallback in the sensitivity study. The fixed 200 iteration multiplier search and other implementation details limit broader efficiency interpretations. The results do not establish complexity advantages or superiority over an optimized low rank solver. The 720 original user evaluations are five seeds combined; the corresponding per seed workload is 144.

### 7.6 Scope of numerical evidence
The experiments support formula consistency, constructed attainment or approximation, finite instance safety checks, and the reported acceptance and call counts. Universal safety and information class sharpness rely on Sections 3 and 4, not on finite sample counts. Three timing repetitions do not establish statistical significance. No observed false acceptance is an empirical statement and is not a directed rounding guarantee.

## 8 Interface to configuration lifecycle management
The certificate addresses whether a fixed candidate can be deployed at a prescribed uncertainty radius. Once a candidate is selected, a separate complete configuration reference can determine its robust radius for subsequent online retain or refresh management. Candidate generation, deployment costs and future channel drift remain outside the present optimization problem. A positive certificate at one radius must not be interpreted as a maximum radius certificate.

## 9 Discussion and limitations
The information class optimum is deliberately unrestricted in channel norm, rank and wireless support. Some extremal constructions need large nominal channel norms and may not be reachable in a physical RIS configuration. Additional information may yield a stronger certificate; the channel norm constraint in (28) is one example, but its complete optimal envelope is not derived.

The mathematical conclusions use real or complex exact arithmetic and finite dimensional compactness for one dimensional minima. The implementation uses ordinary floating point arithmetic with engineering safeguards. A formal numerical certificate would require validated interval arithmetic or directed rounding.

Runtime and acceptance gains depend on the candidate pool and solver implementation. A matrix based bound can outperform a six scalar envelope because it retains more information. The present results do not justify universal speed claims. Further empirical assessment should fix the protocol before using additional independent seed batches and should report paired variability rather than treating correlated user evaluations as independent observations.

A complete formula level comparison with operator interpolation, contraction completion and information based optimal recovery literature remains necessary before making a novelty claim. The paper does not claim first discovery of the underlying tools. The explicit proofs establish the stated mathematical results; they do not establish bibliographic novelty.

## 10 Conclusion
Joint spectral compatibility yields a continuous directional envelope and a safe shared budget QoS certificate. The certificate dominates the previous envelope under the same six scalar information budget. A positive definite completion proves sharpness over the nondegenerate information class, with explicit treatment of zero alpha and radius. At the degenerate boundary, an explicit diagonal construction attains the phase coupled boundary envelope. The difference between boundary and nondegenerate limiting optima is an information class effect rather than a physical QoS discontinuity.

The tested sensitivity grid shows a 2.12 percentage point feasible coverage increase at 256 cells and a 3.41% reduction in reference calls, without a stable runtime advantage over the previous envelope. These results support the limited certification interface, while the main conclusions concern the spectral extremum and information class optimality. Further work concerns additional physical information, independent empirical validation and formula level literature positioning.

## References
[1] G. Zhou, C. Pan, H. Ren, K. Wang, and A. Nallanathan, “A Framework of Robust Transmission Design for IRS-Aided MISO Communications With Imperfect Cascaded Channels,” IEEE Transactions on Signal Processing, vol. 68, pp. 5092–5106, 2020. DOI: 10.1109/TSP.2020.3019666. Author manuscript: https://arxiv.org/abs/2001.07054.

[2] F. Rendl and H. Wolkowicz, “A Semidefinite Framework for Trust Region Subproblems with Applications to Large Scale Minimization,” Mathematical Programming, vol. 77, pp. 273–299, 1997.

[3] C. Fortin and H. Wolkowicz, “The Trust Region Subproblem and Semidefinite Programming,” Optimization Methods and Software, vol. 19, no. 1, pp. 41–67, 2004. DOI: 10.1080/10556780410001647186. Author manuscript: https://www.math.uwaterloo.ca/~hwolkowi/henry/reports/2004fortinwolkTRS.pdf.

## Reproducibility
The supplied code runs the original experiments with run_experiments.py, the sensitivity study with audit_and_stress.py, and per unit decision checks with verify_review.py. The repaired construction checks are in verify_v41_proofs.py, and the resolution and dominance checks are in review_v4_checks.py. Original outputs remain in results/; reproduction, sensitivity and supplementary checks are in results_review/. Numerical minima are diagnostics; interval lower bounds are used for acceptance.
