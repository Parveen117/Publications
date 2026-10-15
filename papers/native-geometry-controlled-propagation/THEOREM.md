# GCP-1 — Quantum geometry controls propagation: feedback and the single-metric obstruction

15 October 2026. Repository continuation of CCL-1 and MGH-1.

**Result.** A positive operator-valued lapse determines the principal part of a joint matter–geometry evolution. Its metric branches have their own characteristic cones, shared by all matter species using the same lapse operator. An explicit bounded family includes a flat branch and a genuinely curved Lorentzian branch. The joint evolution is unitary, local and has a finite propagation envelope. Its kinetic coupling produces matter–geometry feedback and entanglement without adding a separate mass-weighted holonomy phase. Compressing geometry to its mean generally destroys the principal operator's second moment. A variance criterion identifies when a single-metric reduction is exact.

The finite geometry register, branch lapse functions, ordering coordinate and geometry-mixing generator are declared model data. The construction is a rigorously specified quantum matter–geometry propagation model. Selecting a gravitational action, implementing its constraints and recovering its interacting continuum physics remain separate obligations. No claim that quantum gravity is completely solved is made.

## 1. Native input and what is newly supplied

Use the existing native cut and turn

\[
C=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
C^2=I,\ R^2=-I,\ CR=-RC,
\]

with `C†=C`, `R†=-R` and central `iota²=-1`. Let a finite independent geometry register have orthogonal cuts `Q_j`, summing to its identity. Supply real functions

\[
0<v_-\le v_j(x)\le v_+<\infty,
\quad v_j\in C_b^2(\mathbb R),\qquad
V(x)=\sum_j v_j(x)Q_j.                                      \tag{1}
\]

The analytic carrier is the positive completion
`L²(R,dx; C² tensor C^N)` already admitted by CCL-1, enlarged by the finite geometry register. This completion is an explicit analytic adapter, not a new primitive axiom of the cut algebra. Let `B†=-B` be a fixed finite matrix acting on geometry and `g` a real dimensionless rest-turn rate.

The source is the geometry observable `V`, with its full quantum state retained. The target is the matter propagation operator. The map is

\[
\boxed{
\mathcal G=(C\otimes V)\partial_x
+\tfrac12 C\otimes V'
+gR\otimes V+I\otimes B,\qquad
\partial_s\Psi=\mathcal G\Psi.}                            \tag{2}
\]

`s` orders operations, as in CCL-1; no laboratory clock or value of `c` is inferred. The factor `V` in the rest term adopts the branch proper-count increment `v_j ds`. It is a chosen minimal coframe coupling, not a consequence of the principal symbol alone. General bounded anti-dagger local terms may be added without changing the principal symbol or the support argument below.

The geometry register now enters the *derivative coefficient*. This is the gap left by CCL-1 Section 7: a geometry-dependent phase alone was a zeroth-order term and could not change the propagation metric. MGH-1's quantum composition and ensemble-consistency contracts are reused. Its Euclidean hinge geometries are not silently identified with the Lorentzian branches supplied here.

## 2. Why the half-derivative term is forced

**GCP-T1 (symmetric transport).** For a smooth Hermitian `A(x)` on compactly supported smooth fields,

\[
(A\partial_x)^\dagger=-A\partial_x-A'.
\]

Consequently `A∂x + A'/2` is formally anti-dagger. Among operators `A∂x + D(x)` with this fixed leading coefficient, anti-daggerness requires

\[
D+D^\dagger=A',
\]

so `D=A'/2+K` with `K†=-K`. In (2), `A=C tensor V` and the remaining displayed coefficients are anti-dagger.

**Proof.** Integration by parts moves the derivative onto both `A` and the test field; the endpoint term vanishes. Comparing the zeroth-order part of `G†=-G` gives the condition. This proves the uniqueness of the Hermitian correction, while retaining the anti-dagger freedom. QED.

Leaving out `V'/2` creates a spurious change of the positive norm. Thus the metric-to-propagator map has a necessary correction, not merely a substitution of a position-dependent speed into the flat equation.

**GCP-T2 (unitary group and finite propagation envelope).** Under (1), (2) has a skew-adjoint realization generating a strongly continuous unitary group. If the initial field has support in `[a,b]`, then for `s>=0`

\[
\operatorname{supp}\Psi(s)\subseteq[a-v_+s,b+v_+s].          \tag{3}
\]

**Proof of existence.** For each branch set

\[
y_j(x)=\int_0^x\frac{du}{v_j(u)},\qquad
(T_jf)(y)=\sqrt{v_j(x_j(y))}\,f(x_j(y)).                    \tag{4}
\]

The positive upper and lower bounds make `y_j` a diffeomorphism of the real line. Substitution in the norm shows `T_j` is unitary from `L²(dx)` to `L²(dy)`. The product and chain rules give exactly

\[
T_j(v_j\partial_x+\tfrac12v'_j)T_j^{-1}=\partial_y.        \tag{5}
\]

The direct sum of these transforms conjugates the leading transport to the self-contained translation group `C∂y`, on its pulled-back `H¹` domain. The remaining operator `W=gR tensor V+I tensor B` is bounded and anti-dagger. In the interaction picture of this translation group, its Dyson series converges in operator norm on finite count intervals, with the `n`th time integral bounded by `(|s| ||W||)^n/n!`. Each bounded-generator interaction evolution is unitary by differentiating its squared norm, and the inverse evolution exists. Multiplication by the free unitary gives the asserted group; on the generator domain its generator is the sum (2).

**Proof of locality.** The free characteristic in branch `j` obeys `|dx/ds|=v_j(x)<=v_+`. The weight in (4) changes amplitudes but not support. Every insertion of the bounded multiplication operator `W` is local in the original `x` address. Every term of the Dyson series therefore has total displacement at most `v_+|s|`; the norm limit preserves that support property. QED.

For regular fields the exact local balance is

\[
\partial_s(\Psi^\dagger\Psi)=\partial_x(\Psi^\dagger A\Psi),
\quad |\Psi^\dagger A\Psi|\le v_+\Psi^\dagger\Psi.          \tag{6}
\]

The flux sign reflects the convention `S_a f(x)=f(x+a)`; physical support speeds are its magnitudes. This theorem is not a construction of a positive-energy quantum field vacuum. The self-adjoint count generator `iota G` need not be bounded below.

## 3. The metric comes from the same principal coefficient

**GCP-T3 (branch cones and universal coupling within a branch).** The principal symbol of (2), using the real covector convention, is

\[
\sigma(\tau,\xi)=\tau I-\xi C\otimes V(x).
\]

In branch `j` its determinant is

\[
\det\sigma_j=\tau^2-v_j(x)^2\xi^2.                         \tag{7}
\]

Its null set is the null set of the inverse of the branch metric

\[
\boxed{e_j^0=v_j(x)ds,\quad e_j^1=dx,\qquad
q_j=-(v_j(x))^2ds^2+dx^2.}                                \tag{8}
\]

Any finite collection of rest rates `g_a` has the same cone in a given branch when all sectors use the same `V`. Their mass terms and local mixing are zeroth order. This is a common metric *within each geometry branch*, not a claim that a superposition of distinct metrics has one cone.

**Proof.** The native cut has eigenvalues `+1,-1`; restricting the symbol to `Q_j` gives the two factors `tau-v_j xi` and `tau+v_j xi`. The inverse metric of (8) gives `-tau²/v_j²+xi²`, proportional by a nonzero factor to (7). Zeroth-order terms do not enter this calculation. QED.

The principal symbol determines a conformal class. Choosing `e_j^1=dx` supplies the spatial comparison unit and selects the representative (8); its full scale is not uniquely reconstructed from the null set. Frame changes must transform the coframe, symbol and source together.

### 3.1 A curved, bounded, exact branch family

Take two branches

\[
v_0(x)=1,\qquad
v_1(x)=1+\frac{a x^2}{1+x^2},\qquad a=3/5.                 \tag{9}
\]

They satisfy `1<=v_j<=8/5`. This is a selected rational test profile, not a fitted observation or a gravitational field equation. In the torsion-free, metric-compatible adapter of (8), the Cartan equations give

\[
\omega^0{}_1=\omega^1{}_0=v'_j ds,\qquad
\Omega^0{}_1=-\frac{v''_j}{v_j}e_j^0\wedge e_j^1,
\qquad \mathcal R_j=-2\frac{v''_j}{v_j}.                   \tag{10}
\]

The convention is `R^rho_{sigma mu nu}=partial_mu Gamma^rho_{nu sigma}-partial_nu Gamma^rho_{mu sigma}+...`. Equation (10) follows directly by differentiating the displayed connection; the torsion equations vanish. This is the explicitly declared metric-curvature adapter, not an assertion that every native commutator equals Riemann curvature.

Here

\[
v'_1=\frac{2ax}{(1+x^2)^2},\qquad
v''_1=\frac{2a(1-3x^2)}{(1+x^2)^3},
\]

so

\[
\mathcal R_0=0,\quad \mathcal R_1(0)=-12/5,\quad
v_1(1)=13/10,\quad\mathcal R_1(1)=6/13.                    \tag{11}
\]

At `x=0` the two cones and their first derivatives coincide while their curvatures differ. At `x=1` both the speed reading and curvature distinguish them. No coordinate rescaling can turn the nonzero curvature branch into the flat branch. This gives genuine curved branch data for the propagation operator, beyond a static two-label phase.

### 3.2 The full branch equation matches its coframe connection

The connection match is stronger than matching a null cone. Reuse the native Clifford realization

\[
\Gamma^0=R,\quad\Gamma^1=CR,\qquad
\{\Gamma^a,\Gamma^b\}=2\operatorname{diag}(-1,1)^{ab}I.
\]

In the declared torsion-free spin-connection adapter, lower the first index of (10). Then `omega_01=-v' ds` and

\[
\mathcal A_s=\tfrac14\omega_{ab,s}\Gamma^a\Gamma^b
=-\tfrac12v'C,\qquad \mathcal A_x=0.
\]

The sum includes both antisymmetric index orders. The branch covariant equation is

\[
\left[\frac{R}{v}\left(\partial_s-\frac{v'}2C\right)
+CR\partial_x+gI\right]\psi=0.                            \tag{11a}
\]

Multiplication on the left by `vR`, using `RCR=C`, gives precisely

\[
\partial_s\psi=(vC\partial_x+\tfrac12v'C+gvR)\psi.
\]

This is the branch restriction of (2) before the configuration mixing. It matches the mass and connection terms as well as the principal symbol. No gravitational field equation is inserted. The coframe and the torsion-free spin adapter are the stated hypotheses, reused from the earlier MG/SM type of construction.

Furthermore the actual derivative commutator is

\[
[\nabla_s,\nabla_x]=\tfrac12v''C
=-\tfrac14v\mathcal R\,C.                                 \tag{11b}
\]

For branch 1 at the origin it equals `(3/5)C`, while branch 0 gives zero. Thus the same specified geometry determines the propagation law, its connection and its operator curvature. Off-diagonal configuration mixing `B` is additional dynamics; this identity does not prove that an arbitrary `B` is compatible with full spacetime covariance.

## 4. Quantum superposition cannot usually be replaced by one average metric

**GCP-T4 (single-metric obstruction).** For two distinct lapse eigenvalues at a point, the full principal determinant is

\[
\det\sigma=(\tau^2-v_0^2\xi^2)(\tau^2-v_1^2\xi^2).         \tag{12}
\]

If both geometry cuts are retained, this cannot be the determinant of two identical copies of one scalar-speed first-order cut operator unless `v_0=v_1` (positivity excludes the opposite-sign alternative).

**Proof.** A repeated scalar-speed determinant would be `(tau²-u²xi²)²`. Equating the two polynomial coefficients gives `2u²=v_0²+v_1²` and `u⁴=v_0²v_1²`. Subtracting yields `(v_0²-v_1²)²=0`, hence equality. QED.

For a normalized geometry vector `chi` and its retained cut `P=I_dir tensor |chi><chi|`, put `A=C tensor V`. Then

\[
\boxed{PA^2P-(PAP)^2
=\operatorname{Var}_\chi(V)P
=PA(I-P)AP\succeq0.}                                     \tag{13}
\]

**Proof.** Insert `I=P+(I-P)` between the two `A` factors. The second term is the dagger square of `(I-P)AP`. Using `C²=I` identifies the scalar with `<V²>-<V>²`. QED.

For the balanced geometry cut and `V=diag(1,2)`, the lost coefficient is exactly `1/4`; the mean speed is `3/2`, while the two actual speeds are `1,2`. At `x=1` in the curved family (9), the analogous loss is `9/400`. It is a lost principal-operator moment, not a numerical error. Averaging `q`, averaging its inverse and averaging `V` are different prescriptions; none removes the exact obstruction (12).

This is a concrete instance of the programme's retained-memory principle. For a **homogeneous** fixed Fourier mode, the full finite generator `G_k` has the exact compression

\[
P(z-G_k)^{-1}P=
\left[z-PG_kP-PG_kQ(z-QG_kQ)^{-1}QG_kP\right]^{-1}_{P},     \tag{14}
\]

where `Q=I-P`, the relevant blocks are restricted to their ranges, and the inverses must exist. This is the existing Schur/return identity applied to the geometry record, not a new general Schur theorem. Dropping the return term is precisely the extra approximation that a mean-only propagation law makes. Equation (14) is not asserted as a fixed-mode formula for the spatially varying case, where Fourier modes couple.

**GCP-T5 (exact single-metric sector).** Let `S` be a geometry subspace such that, for every `x`, `V(x)|S=v(x)I_S`, and suppose `B S` is contained in `S`. Then the joint evolution restricted to `L² tensor S` factorizes into the scalar-lapse matter evolution and `e^{sB|S}`. It generates no new matter–geometry entanglement from product input. For a retained rank-one geometry cut, absence of the principal leakage in (13) is equivalent to `chi` being an eigenvector of `V(x)` at that point.

**Proof.** On the invariant subspace the two generator terms act on independent factors and commute. Their unitary groups factorize. The last assertion is `||(V-<V>)chi||²=0`. QED.

A pointwise zero variance is insufficient for a global reduction if `B` subsequently moves the state into different lapse eigenspaces. A scalar-lapse sector may also be curved: vanishing quantum geometry dispersion and vanishing spacetime curvature are different conditions.

## 5. The same kinetic operator produces feedback

For any fixed finite geometry observable `F`, on the common generator core the exact Heisenberg identity is

\[
\boxed{
\frac{d}{ds}\langle I\otimes F\rangle
=\left\langle
C\otimes[F,V]\partial_x+\tfrac12C\otimes[F,V']
+gR\otimes[F,V]+I\otimes[F,B]
\right\rangle.}                                          \tag{15}
\]

It follows by differentiating `U†(I tensor F)U`, giving `[I tensor F,G]`. Thus quantum geometry observables that do not commute with the lapse receive a source from the matter translation/turn operators. Geometry simultaneously controls those translation/turn coefficients in (2). This is reciprocal coupling within one fixed linear joint law; no state-dependent nonlinear evolution of an unknown density matrix is inserted.

For `F=Q_j`, the first three commutators vanish, but `[Q_j,B]` need not. Populations change through geometry coherences, and the kinetic coupling changes those coherences. The feedback is not visible by examining only the first population derivative.

### 5.1 Exact momentum-dependent geometry response

Use a periodic homogeneous control with `g=0`, `V=diag(1,2)`, `B=omega R_g`, and one exact Fourier mode `k`. This control has flat branch metrics and is distinct from the curved profile (9). In a direction eigenstate `C=r`, `r=+1` or `-1`, the geometry generator is

\[
G_{k,r}=\iota kr\bar v I+\iota kr\delta v Z+\omega R_g,
\quad \bar v=3/2,\quad\delta v=-1/2.
\]

The traceless part squares to `-h_k² I`, with `h_k²=k²(delta v)²+omega²`. Starting in `Q_0`, direct exponentiation gives

\[
p_1(s;k)=\frac{\omega^2}{h_k^2}\sin^2(h_ks),
\quad p_1''(0;k)=2\omega^2,
\quad p_1^{(4)}(0;k)=-8\omega^2h_k^2.                       \tag{16}
\]

At `omega=2/3`, the fourth derivatives are `-128/81` for `k=0` and `-200/81` for `k=1`. The geometry observable therefore responds differently to different matter momentum preparations under the same fixed joint law. The two preparations have the same initial geometry and the same second population derivative. Keeping only that derivative would miss the feedback.

Set `H=iota G_k`, `H_kin=-k C tensor V`, `H_geo=iota omega I tensor R_g`. Then

\[
\dot H_{\rm kin}=2k\omega\delta v C\otimes X_g,
\qquad \dot H_{\rm geo}=-\dot H_{\rm kin}.                 \tag{17}
\]

For `k=1`, direction cut `C=+1` and geometry `X_g=+1`, the initial exchange rates are exactly `-2/3` and `+2/3`. Their sum is zero. These are dimensionless count-generator readings; physical gravitational stress-energy conservation is not inferred from this finite identity.

## 6. Exact entanglement generated by the metric-dependent principal term

For `k=1`, `omega=2/3`, `h=5/6`, choose `s=(6/5) arccos(3/5)`. Remove only the commuting local direction phase `exp(iota k bar(v) s C)`, which cannot change entanglement. The remaining joint unitary is exactly

\[
U=\frac35 I+\frac{16}{25}I\otimes R_g
-\frac{12\iota}{25}C\otimes Z_g.                           \tag{18}
\]

Its unitarity follows from anticommutation of the two anti-dagger generators and `9/25+256/625+144/625=1`. For an initially pure balanced direction cut `(I+X_dir)/2` and geometry `Q_0`, the amplitude array, with the common `1/sqrt(2)` factor outside, is

\[
\Psi=\frac1{\sqrt2}
\begin{pmatrix}(15-12\iota)/25&16/25\\(15+12\iota)/25&16/25\end{pmatrix}. \tag{19}
\]

Every entry of its density matrix is a native rational cut. Exact partial traces give

\[
\boxed{
p_1=\frac{256}{625},\qquad
S_L=1-\tr\rho_g^2=1-\tr\rho_m^2=\frac{73728}{390625}>0.}    \tag{20}
\]

**Proof.** In the general form of (19), write the branch vectors as `(c-iota d,b)` and `(c+iota d,b)`, with `c²+b²+d²=1`. Their overlap has squared magnitude `1-4b²d²`; the reduced purity of the equal pure bipartite preparation is `(1+|overlap|²)/2`. Hence `S_L=2b²d²`. Substitute `b=16/25`, `d=12/25`. QED.

The metric branches are independent quantum registers. Entanglement in (20) is between the direction subsystem at fixed mode and that register. A mixed direction preparation `I/2` gives the same geometry marginal but a separable joint mixture of the two branch products. The positive geometry entropy alone would therefore not certify entanglement without the stated pure preparation. The native certificate checks both cases.

With `delta v=0`, the principal coupling disappears and the joint evolution factorizes. With the displayed initial geometry and `omega=0`, no geometry superposition is created and this preparation remains a product. These are explicit controls, not a statement that every possible initial superposition stays separable when `omega=0`.

## 7. Frame covariance of the complete operator

Let `T(s,x)` be a smooth unitary change of the whole internal frame, write `G=A∂x+A'/2+K`, and set `Psi_new=T† Psi`. Then

\[
A_{\rm new}=T^\dagger A T,
\]
\[
K_{\rm new}=T^\dagger K T
+\tfrac12(T^\dagger A T_x-T_x^\dagger A T)-T^\dagger T_s.    \tag{21}
\]

The transformed equation is `∂s Psi_new=(A_new ∂x+(A_new)_x/2+K_new)Psi_new`. Direct product differentiation proves (21), and each extra zeroth-order term is anti-dagger. Characteristic eigenvalues are preserved under the pointwise conjugation. In particular, a position-dependent frame change cannot be made by conjugating the lapse and discarding its induced connection terms.

This covariance concerns the admitted internal frame; spacetime diffeomorphism constraints and refoliation invariance are not established by it.

## 8. Proof boundary and continuation target

| Obligation | Result here |
|---|---|
| Geometry affects matter's derivative coefficient | Explicit operator map (2) |
| Genuine curved Lorentzian branch | Bounded rational profile, curvature (11) |
| Positive norm, unique joint evolution, finite propagation | Written analytic proof T1–T2 |
| Same geometry for different matter species | Common branch symbol T3; shared coupling remains a premise |
| Replace quantum geometry by its mean | Refuted in general; exact variance and determinant obstruction |
| Exact single-metric reduction | Invariant scalar-lapse subspace criterion T5 |
| Matter-to-geometry feedback and geometry-to-matter transport | One joint generator; exact response/exchange controls |
| Entanglement caused by the principal metric coupling | Exact finite-mode demo (20) |
| Derive the mixing law B and branch family from one gravitational action | Open |
| Match the prior Euclidean hinge to this Lorentzian coframe | Open; no implicit identification |
| Physical constraint algebra, positive quantum field interpretation, 3+1 interacting limit | Open |

The next step is to replace the supplied geometry mixing by a law derived from a selected native gravitational source/action, with an explicitly closed constraint system. It must preserve the same matter–metric map and the retained geometry record established here. Arbitrarily adding another phase would not close that obligation.

No Einstein field equation, numerical `c`, measured mass or fitted physical datum was used in this stage. The continuum completion, coframe representative, torsion-free curvature adapter and finite quantum interpretation are stated where used. Native word certificates cover algebraic identities, with exact rational fixtures and negative controls. The infinite-dimensional existence and support arguments are written proofs, not proof-assistant certificates.
