# AGD-1 — Geometry dynamics from the native cut action

15 October 2026. Continues CID-D2 and GCP-1 without changing their frozen packets.

**Result.** The real-circle sector of the existing native cut action has a fixed inertia and potential. Under an explicit periodic canonical quantization contract it yields an infinite geometry generator, replacing GCP-1's freely supplied finite mixing matrix. The same cut coordinate determines the lapse observable, so GCP-1's two branches are endpoints of one dynamical geometry family. The full joint operator is essentially self-adjoint, preserves norm and has a finite spatial propagation envelope. An explicit geometry-momentum source law and a factorial finite-section error estimate follow. A purported parity constraint fails once matter is coupled; only the stated global clock constraint is implemented, not the gravitational hypersurface-deformation algebra.

## 1. Reuse the action, not a new fitted mixing matrix

CID-D2 selected native anti-dagger controls `(R,J,L)`, information weights `(1,1,2)`, curvature weights `(1,2,2)`, `alpha=2`, and `beta=1`. Its already proved rank-one potential is

\[
\Phi_I=6-2z^2,\qquad \Phi_C=4-2z^2,\qquad
\Phi=10-4z^2.                                             \tag{1}
\]

These control and action choices remain assumptions. They are not adjusted to fit any example below.

Work in CID-D2's parallel reference gauge `A_s=0`, with covariantly constant controls represented by constant matrices. The operation parameter is denoted `s` to match GCP; it is not yet a calibrated physical clock. Thus `D_s P=dot(P)` in the restricted action. A nontrivial retained boundary holonomy or a different parameter normalization would be additional data, not something discarded by this reduction.

On the real cut circle put

\[
P(\theta)=\tfrac12(I+\cos\theta\,Z+\sin\theta\,X),
\quad Z=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
X=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad R=XZ.
\]

`P²=P=P†`, `Tr P=1`, and `Tr((partial_theta P)²)=1/2`. The circle is the fixed locus of the reflection of the other Bloch coordinate. The round cut metric and the potential (1) respect that reflection. Initial position and velocity tangent to this circle therefore remain on it in the classical CID-D2 equation. It is a consistent classical sector, not a parametrization of the full cut sphere.

**AGD-T1 (action reduction).** The restricted action, momentum and Hamiltonian are

\[
S_g=\int\left[\tfrac14\dot\theta^2-(8-2\cos2\theta)\right]ds,
\quad p=\tfrac12\dot\theta,
\quad h_g=p^2+8-2\cos2\theta.                              \tag{2}
\]

Consequently `dot(theta)=2p` and `dot(p)=-4 sin(2 theta)`.

**Proof.** Substitute the cut derivative into the existing kinetic term `Tr(dot(P)²)/2`, giving `dot(theta)²/4`. Substitute `z=cos(theta)` into (1). Legendre transformation with `p=dot(theta)/2` gives (2). The force follows by differentiation. QED.

The inertia `1/2` and coefficients `8,-2` are fixed by this inherited action. Introducing an unrelated two-state rotation frequency at this point would discard that information.

## 2. Quantization contract and the infinite circle carrier

Declare scalar periodic geometry wavefunctions in `L²(S¹,dtheta/(2pi))`, and a positive dimensionless action scale `eta`. With the Fourier basis `|n>=exp(iota n theta)`, `n in Z`, define

\[
U|n\rangle=|n+1\rangle,\quad L|n\rangle=n|n\rangle,
\quad \widehat p=\eta L=-\iota\eta\partial_\theta.
\]

Thus `U†U=UU†=I` and `[L,U]=U` on finite Fourier words. Quantize (2) using this momentum representation and multiplication by its potential:

\[
\boxed{H_g=\eta^2L^2+8I-(U^2+U^{\dagger2}),\qquad
B_g=-\frac{\iota}{\eta}H_g.}                              \tag{3}
\]

The periodic scalar sector, representation and `eta` are declared quantization inputs. The cut action alone does not select this quantum prescription or a physical Planck constant. A twisted circle boundary condition changes the spectrum to `eta(n+zeta)`; the retained boundary phase must be specified before claiming unique quantum dynamics. This packet uses `zeta=0`.

**AGD-T2 (geometry generator and finite-carrier obstruction).** `H_g` is self-adjoint on periodic `H²`, bounded below by `6I`, and has a compact resolvent. Its Fourier diagonal is `eta²n²+8`, and its only off-diagonal matrix elements are `-1` between `n` and `n±2`. No nonzero finite-dimensional carrier can satisfy both exact unitary `U` and `[L,U]=U`.

**Proof.** `eta²L²` is diagonal self-adjoint, and `8I-(U²+U†²)` is a bounded Hermitian multiplication operator in `[6,10]`. Its bounded perturbation has the same domain and is self-adjoint; the bounded perturbation construction can equivalently be obtained by the norm-convergent interaction series for its unitary evolution. The diagonal resolvent is compact because its eigenvalues tend to zero. The resolvent identity preserves compactness for the perturbation. The Fourier entries follow directly from the shift definition. In finite dimension, `[L,U]=U` and unitarity imply `U†LU=L+I`; taking the finite trace gives `0=dimension`, a contradiction. QED.

Finite matrices will therefore be controlled approximations to (3), not exact substitutes for its circle algebra.

## 3. One cut observable supplies the full metric family

Use the fixed reference cut `Q_1=(I-Z)/2`. The native overlap is

\[
r(\theta)=\operatorname{Tr}(P(\theta)Q_1)
=\tfrac12(1-\cos\theta).
\]

This scalar is invariant under simultaneous conjugation of the cut and reference. Fix the same profile and amplitude used by GCP-1,

\[
f(x)=\frac{x^2}{1+x^2},\qquad a=3/5,
\quad v(\theta,x)=1+a f(x)r(\theta).                       \tag{4}
\]

As a geometry observable,

\[
\boxed{V(x)=\left(1+\frac{af(x)}2\right)I
-\frac{af(x)}4(U+U^\dagger),\qquad 1\le V(x)\le8/5.}      \tag{5}
\]

At `theta=0,pi`, its values are exactly the two GCP-1 lapse branches. All intermediate circle cuts now carry their own metric; the geometry is no longer restricted to two labels. The profile `f`, amplitude `a` and coframe interpretation remain supplied comparison data. They are not selected by the action (2).

GCP-1's calculation holds at every fixed `theta`:

\[
e^0=v\,ds,\quad e^1=dx,\quad
q_\theta=-v^2ds^2+dx^2,\quad
\mathcal R(\theta,x)=-2\frac{v_{xx}}{v},\qquad
[\nabla_s,\nabla_x]=\tfrac12 v_{xx} C.                     \tag{6}
\]

These are conditional geometry-configuration readouts; `s`-derivatives of a single substituted classical path `theta(s)` are not being added to (6). The quantum state evolves over configuration labels, rather than being replaced by one classical time-dependent metric. In particular

\[
\mathcal R(\theta,0)=-4a r(\theta),\qquad
\mathcal R(\theta,1)=\frac{a r(\theta)}{1+a r(\theta)/2}.
\]

The second readout is nonlinear: substituting the mean of `r` before evaluation is generally a different observable.

## 4. The joint action-based quantum law

Take the positive completion
`H = L²(R_x times S¹_theta; C²_direction)`. Keep GCP-1's native direction cut `C=Z`, turn `R`, rest rate `g` and symmetric derivative convention. Define

\[
H_m=\iota\eta C\left(v\partial_x+\tfrac12v_x\right)
+\iota\eta g Rv,
\]
\[
\boxed{
\iota\eta\partial_s\Psi=H\Psi,\qquad H=H_g+H_m.}          \tag{7}
\]

Equivalently, `partial_s Psi` has the GCP-1 transport and mass terms plus the geometry generator `B_g` of (3). `B_g` is now fixed once the inherited action and quantization contract are fixed; it is no longer a free finite mixing matrix. Equation (7) is linear on the full joint state and hence respects QGC-1's ensemble-consistency requirement.

The quantum action is

\[
S_Q=\int ds\left[
\frac{\iota\eta}{2}(\langle\Psi,\dot\Psi\rangle-
\langle\dot\Psi,\Psi\rangle)-\langle\Psi,H\Psi\rangle
\right].                                                  \tag{8}
\]

Endpoint-fixed variation in `Psi†` gives (7). This variational identity is not used to pretend that an arbitrary `H` has been derived: its geometry part was already fixed by (2) and the stated quantization rule, while its matter part is the retained GCP coframe prescription.

### 4.1 A full-domain existence proof

**AGD-T3 (essential self-adjointness).** For `eta>0`, fixed real `g`, and the smooth bounded coefficient (4), `H` on the core of smooth compactly supported `x` functions and smooth periodic `theta` functions is essentially self-adjoint. Its closure generates a unique norm-preserving unitary evolution on the full joint carrier. This is an infinite geometry result, not a finite Fourier check.

**Proof.** Let `N=I-partial_x²-partial_theta²`, a positive self-adjoint comparison operator on `H²`. The core is dense in `H²`. Integration by parts makes `H` symmetric. All required coefficient derivatives are bounded. Expanding products gives constants `c_1,c_2<infinity` such that

\[
\|Hu\|\le c_1\|Nu\|,\qquad
|\langle u,[H,N]w\rangle|
\le c_2\|N^{1/2}u\|\|N^{1/2}w\|.                         \tag{9}
\]

To see the order in the second estimate, `-eta² partial_theta²` commutes with `N`. For a coefficient `b`, the commutator of `b partial_x` with either comparison Laplacian consists of terms `2 b_j partial_j partial_x+b_jj partial_x`, up to an overall sign. Integration by parts moves one of the two derivatives onto `u`, leaving bounded coefficients and at most first derivatives on each field. Commutators of multiplication terms with `N` have order at most one. Cauchy–Schwarz gives (9).

Here is the regularization argument completing the proof, rather than assuming that formal symmetry suffices. Put `R_e=(I+eN)^(-1)`. Its range is in `H²`, hence in the closure domain of `H`, and the distributional identity

\[
[H,R_e]=-e R_e[H,N]R_e                                  \tag{10}
\]

extends by (9). If a deficiency vector satisfies `H*psi=+iota psi` or `-iota psi`, the imaginary part of
`<R_e psi,H R_e psi>` is zero. Using (10), its absolute value implies

\[
\|R_e\psi\|^2
\le c_2 e\|N^{1/2}R_e^2\psi\|\|N^{1/2}R_e\psi\|.         \tag{11}
\]

The right side tends to zero: by the spectral resolution of `N`,
`e ||N^(1/2)R_e psi||²` tends to zero by dominated convergence, since
`e lambda/(1+e lambda)²<=1/4`, and the analogous norm with `R_e²` is smaller. The left side tends to `||psi||²`. Both deficiency spaces are therefore zero. This proves essential self-adjointness; its spectral exponential gives the unique unitary group. QED.

The comparison operator in this proof is an analytic domain tool, not an added physical Hamiltonian. The total `H` need not be bounded below because of the matter translation term. A positive-energy quantum field construction is still a separate obligation.

The analytic method is the standard commutator criterion for essential self-adjointness, with its hypotheses checked in (9) and its regularization argument supplied in (10)–(11). For a reference statement, see Lauritz van Luijk, *Energy-limited quantum dynamics*, [arXiv:2405.10259, Corollary B.2](https://arxiv.org/html/2405.10259v2#A2). This is mathematical infrastructure; neither that theorem nor the native finite certificate selects a gravitational action.

### 4.2 Spatial causality survives unbounded geometry dynamics

**AGD-T4 (finite spatial envelope).** If the initial joint field is supported in `x in [b_-,b_+]`, then for `s>=0`

\[
\operatorname{supp}_x\Psi(s)
\subseteq[b_--(8/5)s,b_++(8/5)s].                         \tag{12}
\]

**Proof.** Geometry derivatives and all multiplication terms commute with bounded smooth scalar functions `F(x)`. On the core,

\[
[H,F]=\iota\eta C v F_x.
\]

This bounded commutator extends to the closure domain. Therefore
`d_s<F>=<F_s-Cv F_x>`. Choose a smooth increasing exterior cutoff
`F(s,x)=chi(x-b_+-(8/5)s)`, zero for arguments at most zero. Its derivative is
`-<(8I/5+Cv)chi'> <=0`, because the eigenvalues of `Cv` lie between `-8/5` and `8/5`. The initial exterior probability is zero and remains zero. Left exterior cutoffs give the other side. Domain approximation extends the conclusion to every initial `L²` field with the stated support. QED.

Equivalently, after integrating over the periodic geometry variable,
`rho(s,x)=integral ||Psi||² dtheta/(2pi)` satisfies
`partial_s rho=partial_x integral Psi† C v Psi dtheta/(2pi)`.
The geometry kinetic term contributes a periodic geometry-current divergence whose integral vanishes. No bounded-matrix assumption about `B_g` is used in this proof.

## 5. Geometry momentum receives a derived matter source

On a common smooth core, with expectations taken in the full joint state,

\[
\boxed{
\frac{d}{ds}\langle p\rangle
=-\left\langle4\sin2\theta+\partial_\theta H_m\right\rangle,
}\tag{13}
\]
\[
\partial_\theta H_m=
\iota\eta C\left(v_\theta\partial_x+\tfrac12v_{x\theta}\right)
+\iota\eta g Rv_\theta.
\]

**Proof.** `p=-iota eta partial_theta` gives `[H,p]=iota eta partial_theta(Phi+H_m)`. Multiply by `iota/eta` in the Heisenberg equation. The geometry kinetic term commutes with `p`. QED.

This is a source law from the same joint generator, not an independently postulated mass–phase coupling or a nonlinear function of the unknown source state. Smooth domain-controlled states justify the expectation derivative; the commutator is not asserted bounded on all `L²` states.

### Exact homogeneous source control

For a separate translation-invariant control set `f=1`, `g=0`, fix Fourier mode `k` and direction `C=r`. Then

\[
H_{m,k,r}=-\eta kr\left[1+\tfrac a2-\tfrac a2\cos\theta\right],
\quad
\dot{\langle p\rangle}
=-4\langle\sin2\theta\rangle
+\frac{\eta kra}{2}\langle\sin\theta\rangle.               \tag{14}
\]

At `eta=1`, `k=2`, `r=+1`, `a=3/5`, choose the normalized geometry vector
`(|0>+iota|1>)/sqrt(2)`. Its readings are

\[
\langle\sin\theta\rangle=-1/2,\quad
\langle\sin2\theta\rangle=0,\quad
\boxed{\dot{\langle p\rangle}=-3/10.}                      \tag{15}
\]

For the same preparation the initial geometry-energy rate is `-3/10`, the matter-energy rate is `+3/10`, and their sum is zero. These are dimensionless generator readings in the supplied homogeneous control, not a gravitational experiment.

## 6. Quantum geometry cannot become an exact sharp metric in a normalizable circle state

At any `x` where `af(x)!=0`, the multiplication observable `v(theta,x)` has no nonzero normalizable eigenvector: each fixed-value level set of `cos(theta)` has zero circle measure. Consequently every nonzero normalized geometry wavefunction has strictly positive `Var(V(x))`. This strengthens GCP-T5 for this continuous geometry carrier: an exactly sharp nontrivial branch would be distributional, while ordinary states can only approach it.

For finite momentum variance, positivity of the native inner product gives

\[
\operatorname{Var}(p)\operatorname{Var}(V)
\ge\frac{\eta^2}{4}|\langle v_\theta\rangle|^2.            \tag{16}
\]

**Proof.** Apply Cauchy–Schwarz to `(p-<p>)psi` and `(V-<V>)psi`, and bound their imaginary inner product by half the expectation of `[p,V]=-iota eta v_theta`. The bound may be zero at a symmetric preparation; that does not invalidate the strict-variance level-set argument. QED.

For the state in (15), `Var(p)=1/4`, `<V>=13/10`, `Var(V)=9/200`, and `<v_theta>=-3/20`. The product is `9/800`, above the lower bound `9/1600`. Through GCP's retained-cut identity, `9/200` is exactly the principal second-moment coefficient lost by replacing this geometry state with its mean lapse.

## 7. A rigorous finite-geometry computation budget

The following estimate is for a **homogeneous** `f=1` fixed spatial Fourier mode, not the variable-profile equation. Let `Pi_L` retain geometry indices `|n|<=L`, and define the finite self-adjoint Hamiltonian `H_{k,L}=Pi_L H_k Pi_L`. Split

\[
H_k=D_k+W_k,\quad
D_k=\eta^2L_g^2+8I+(1+a/2)M_k,
\quad M_k=-\eta k C+\iota\eta gR,
\]
\[
W_k=-(U^2+U^{\dagger2})-\tfrac a4 M_k(U+U^\dagger),
\qquad \|W_k\|\le w=2+\tfrac a2\eta\sqrt{k^2+g^2}.        \tag{17}
\]

`L_g` denotes geometry angular momentum here, to distinguish it from the integer truncation radius. Assume `a>=0`.

**AGD-T5 (finite-section error).** For initial geometry support `|n|<=L_0<=L`, put
`r=floor((L-L_0)/2)+1` and `t=w|s|/eta`. With the finite evolution embedded in the full carrier,

\[
\boxed{
\|(e^{-\iota sH_k/\eta}-e^{-\iota sH_{k,L}/\eta})\psi\|
\le 2\sum_{j=r}^\infty\frac{t^j}{j!}\|\psi\|
\le\frac{2e^t t^r}{r!}\|\psi\|.}                         \tag{18}
\]

**Proof.** In the interaction picture of the geometry-diagonal `D_k`, each insertion of `W_k` shifts the geometry index by at most two. All Dyson terms with fewer than `r` insertions stay inside the finite section and agree exactly with the full terms. Each remaining term is bounded by `t^j/j!` in either evolution. The triangle inequality gives the first bound. Using `(r+l)!>=r! l!` gives the second. QED.

For `eta=1`, `g=0`, `k=1`, `a=3/5`, `s=1/10`, `L_0=1`, `L=11`, one has `w=23/10`, `t=23/100`, `r=6`. The elementary series bound `e^t<3` for `0<=t<=1` yields the entirely rational error ceiling

\[
\boxed{\frac{148035889}{120000000000000}<\frac1{800000}.}   \tag{19}
\]

This is a worst-case norm error bound for the stated model and inputs, not an experimental accuracy claim. Finite matrices do not certify the analytic limit by themselves; the path-length argument (18) supplies that link.

## 8. Constraints: what closes and what fails

### 8.1 Kinematic constraints and a global relational clock

The circle parametrization enforces the cut relations pointwise. It does not impose an extra projection onto a special wavefunction. The unitary law preserves the total norm and the admitted domain.

One can parametrize the chosen action by a monotone label `tau`, introduce clock pair `(s,p_s)` and multiplier `N`, and impose the single constraint

\[
\mathcal C=p_s+h=0.
\]

Its quantum history equation is

\[
(-\iota\eta\partial_s+H)\Psi=0.                            \tag{20}
\]

The solution is `Psi(s)=exp(-iota sH/eta)Psi(0)`. Initial data carry a positive conserved inner product, independent of the evaluation slice. On the history domain, a bounded matter–geometry reading transported as
`O_rel(s)=exp(-iota sH/eta) O exp(iota sH/eta)` commutes with the constraint when the required derivatives exist. This follows by direct differentiation.

These histories are generally not normalizable in the additional infinite `s` direction; the physical inner product here is on initial data, not a claimed nonzero `L²(ds dx dtheta)` kernel. There is one global parametrization constraint. Its self-commutator being zero is **not** a proof of local gravitational constraint closure or refoliation invariance.

### 8.2 A concrete constraint rejection

Let `T_pi|n>=(-1)^n|n>`. The derived geometry action has `[T_pi,H_g]=0`, because it only shifts by two. It might therefore be tempting to retain only even geometry modes as a purported physical constraint. Matter invalidates that reduction:

\[
T_\pi U T_\pi=-U,\qquad [T_\pi,V]\ne0
\]

whenever the nonconstant lapse coupling is present. For `eta=1`, `g=0`, `k=2`, `r=+1`, `a=3/5`, the initially even state `|0>` has

\[
\boxed{\|\Pi_{\rm odd}H_k|0\rangle\|^2=9/50.}            \tag{21}
\]

This follows from the two odd amplitudes `ka/4=3/10` at `n=+1,-1`. The symmetry of the uncoupled geometry action is not a constraint preserved by the coupled theory. The map `theta -> theta+pi` also exchanges a cut with its complement, so it cannot be declared a gauge redundancy without transforming its observation/source contract.

### 8.3 Why naive local closure is not obtained

For a diagnostic candidate on spatial half-densities, write `D_f=f partial_x+f'/2` and `A_N=C D_N+gR N`. Direct operator multiplication gives

\[
[D_f,D_h]=D_{fh'-hf'},\quad
[A_N,A_M]=D_{NM'-MN'},
\]
\[
\boxed{[D_f,A_N]=A_{fN'-Nf'}+gRNf'.}                      \tag{22}
\]

The last term is nonzero in general. It is absent in the massless fixed-frame special case but not in the massive one. The coefficient and field transformation laws needed for a gravitational hypersurface-deformation algebra have not been supplied by the global circle model. Equation (22) rejects this *particular naive closure*; it is not a no-go theorem for all possible enlarged constraints. Freezing the coframe while demanding full spacetime gauge closure is precisely an additional assumption that must be repaired.

For completeness, `[D_f,N]=fN'` and the product rule give the first identity. In the second, `C²=I` supplies the derivative bracket; `CR=-RC` makes the two mass cross terms cancel, and the mass–mass bracket vanishes. In the third, the mass contribution is `gR fN'`, whereas the mass part of `A_(fN'-Nf')` is `gR(fN'-Nf')`; their difference is exactly the displayed residual. These identities hold for arbitrary smooth scalar smearings on the common compact-support core, beyond the polynomial fixtures in the certificate.

## 9. Status and the next exact target

| Obligation | AGD-1 result |
|---|---|
| Replace the arbitrary finite geometry mixing | Done under the inherited action and explicit periodic quantization contract |
| Preserve GCP matter–metric map | Same cut overlap supplies the lapse; endpoints reproduce GCP branches |
| Infinite geometry carrier and joint existence | Self-adjointness and locality proved on the full carrier |
| Matter source for evolving geometry | Exact momentum/energy source identities and rational control |
| Finite computation versus infinite geometry | Explicit band-path error estimate, with rational budget |
| Exact sharp classical geometry | Excluded for normalizable states at nonconstant lapse; limiting preparations remain possible |
| Global clock constraint | Explicit history solution and conserved physical inner product |
| Local gravitational constraints | Not closed; parity and naive local-closure obstructions displayed |
| Derivation of the lapse profile, physical action selection, 3+1 gravity | Open |

**Next target:** promote the geometric data needed by local frame changes to dynamical fields and derive their source/action constraints, including the terms that cancel or reinterpret the residual in (22). The current single global cut coordinate does not supply those local degrees of freedom. A second finite phase demo or an auxiliary commuting constraint would not solve that problem.

The inherited action, analytic completion, periodic quantization and coframe dictionary are stated explicitly. Written proofs establish the analytic claims. The native certificate replays exact word identities, rational cut metrics, sparse mode operators, source moments, finite-boundary controls and the rational error ceiling. No empirical fit or completion claim for physical quantum gravity is attached to these results.
