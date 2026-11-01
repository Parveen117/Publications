# YM93 — Local Bochner control and a fourfold Wilson continuation window

Monty Dabas. 6 October 2026. Python 3.12 evidence packet.

YM93 enlarges the actual Wilson gap and corrected-memory window from
0 <= theta <= 1/400 to **0 <= theta <= 1/100**. The new input is a local
Bochner estimate for the logarithm of the actual positive ground state:

\[
\boxed{\sup|\nabla_e\log h_\theta|\le4\theta
       \quad\text{for every finite }\theta\ge0.}                 \tag{1}
\]

The gradient bound is valid beyond the gap window. It is not itself a
positive-curvature or mass-gap theorem at arbitrary coupling. Combining
it with a new one-link second-derivative estimate and a continuation
barrier gives, throughout the enlarged window,

\[
\boxed{\operatorname{gap}(A_\theta-E_0)
       \ge d_{93}(\theta)=\tfrac12-32\theta-520\theta^2
       \ge\frac{16}{125}.}                                  \tag{2}
\]

For straight b-link blocking, b>=2, the actual relative dynamic memory is

\[
\boxed{\Sigma_b(0)\le\eta_{93}(\theta)C_b,\qquad
 \eta_{93}(\theta)=
 \left[\frac{4\theta/9+(1592/3)\theta^2}{d_{93}(\theta)}\right]^2
 \le\frac{418609}{2073600}<\frac14.}                          \tag{3}
\]

These are uniform in finite volume and block factor. They keep the
original Hamiltonian and the entire nonlinear response correction.
The factor four compares endpoints in theta, not gauge coupling or a
physical continuum cutoff. No new weak-coupling continuum construction
is claimed. The methods are standard Bochner, heat-kernel, Duhamel and
continuation methods; the new contribution to this programme is their
explicit combination for the actual vacuum with these certified bounds.

## 1. Model, normalization and dependency boundary

Use the periodic spatial cubic lattice of side N>=4, with native links e
and Wilson plaquettes p. For blocking take N=Lb, L>=4, b>=2. Set

\[
 A_\theta=-\Delta-V,\quad V=\theta\sum_pW_p,\quad
 W_p=\tfrac12\operatorname{Tr}U_p,\quad
 A_\theta h=E_0h,\quad h>0,\quad\chi=\log h,\quad d\mu=h^2d\Phi.
                                                               \tag{4}
\]

The product measure Phi is normalized Haar. Each SU(2) factor is the
round sphere of radius two: Ric=(1/2)g, diameter 2 pi, fundamental
Casimir 3/4. Each link meets four plaquettes. Its potential contribution
V_e and every unit invariant derivative D on that link satisfy

\[
 |V_e|\le M:=4\theta,\qquad |DV|\le M/2,\qquad |D^2V|\le M/4.
                                                               \tag{5}
\]

The last identity uses D^2 W_p=-W_p/4 when e belongs to p. Unit invariant
fields have nabla_D D=0. All Hessians below are covariant Hessians.

On every finite volume, ellipticity, positivity improvement and compactness
give a unique smooth positive normalized h. Its simplicity and smooth
parameter dependence imply continuity in the smooth norms used below.
These finite-volume facts do not assume a volume-uniform interacting gap.

YM88's mixed tensor equation and YM92's plaquette identities are geometric
identities valid at every finite coupling. YM91's conditional-flux principle
is likewise an operator statement. We reuse those identities, while
re-establishing every quantitative response and gap bound needed in the
larger interval. Neither the old estimate beta<=5 theta nor an old numeric
helper is extrapolated beyond theta=1/400. The old spatially weighted walk
majorant is not asserted on the enlarged interval.

## 2. T1 — A factorwise Bochner bound at every finite coupling

The logarithmic ground equation is

\[
 \Delta\chi+|\nabla\chi|^2+V=-E_0.                            \tag{6}
\]

Let U_e=d_e chi, a one-form in the e factor over the full product, and
let L_chi=Delta+2 grad chi dot grad. Product splitting is parallel.
Commuting the differential with the Laplacian in (6) gives

\[
 (-\Delta^{\rm rough}-2\nabla_{\nabla\chi}+\tfrac12)U_e=d_eV.
                                                               \tag{7}
\]

Equivalently, the partial Bochner identity reads

\[
 \tfrac12 L_\chi |U_e|^2
 =|\nabla U_e|^2+\tfrac12|U_e|^2-\langle U_e,d_eV\rangle.
                                                               \tag{8}
\]

The derivative norm in (8) includes all factors; its mixed contributions
are nonnegative. At a maximum of |U_e|^2 on the full compact configuration
space the left side is nonpositive. If this maximum is nonzero, (5) gives
(1/2)|U_e|^2 <= (M/2)|U_e|, hence |U_e|<=M. The zero case is immediate.
This proves (1), with no dimension or volume factor and no small-coupling
assumption. It controls the entire three-component gradient, not just its
coordinate derivatives.

More generally, on a compact Riemannian product with factor Ricci lower
bound rho_e>0, the same calculation gives
sup|d_e log h| <= sup|d_e V|/rho_e. This is an application of the classical
Bochner identity, not a priority claim for a new general geometric method.
The background Ricci tensor in (7) must not be confused with the weighted
Ricci tensor Ric-2 Hess chi used to prove a gap later.

Integrating (1) along one-factor geodesics, and writing
g_e(z)=integral h(x,z) dPhi(x), gives the uniform conditional comparison

\[
 r^{-1}g_e\le h\le rg_e,\qquad r=e^{2\pi M}.
                                                               \tag{9}
\]

This follows from the oscillation of log h at fixed exterior z, without
assuming a particular vacuum density.

## 3. T2 — A sharper diagonal Hessian without assuming a gap

### A complete time-four heat bound

For the native isotropic SU(2) heat kernel, the representation n has
dimension n+1 and electric eigenvalue n(n+2)/4. Character domination gives

\[
 |k_4(x,y)-1|\le\sum_{n\ge1}(n+1)^2e^{-n(n+2)}.
                                                               \tag{10}
\]

The finite positive sums through degree eight for exp(3) and degree seven
for exp(4) are 89641/4480>20 and 16319/315>50, respectively. Separate n=1
and use n(n+2)>=4n for n>=2. With q=1/50, the full infinite sum is at most

\[
 \frac15+\left[\frac{1+q}{(1-q)^3}-1-4q\right]
 =\frac{599222}{2941225}<\frac{21}{100}=:k_*.
                                                               \tag{11}
\]

Thus ||K_4-Pi_e||_(infinity to infinity)<=k_*. No representation tail is
dropped. The older time-nine bound is not simply rescaled to time four.

### The second-derivative Duhamel estimate

Split A=L_e+A_rest-V_e, where A_rest is independent of q_e. Integrating
Ah=E_0h over q_e gives (A_rest-U_e^0)g_e=E_0g_e with
U_e^0=integral V_e h/g_e and |U_e^0|<=M. Positive semigroup comparison yields

\[
 e^{sE_0}e^{-sA_{\rm rest}}g_e\le e^{sM}g_e,\qquad s\ge0.
                                                               \tag{12}
\]

For a unit invariant D on e, set Y=D^2h. The electric Casimir commutes
with D, and the full differentiated equation is

\[
 (L_e+A_{\rm rest}-E_0)Y
 =V_eY+(D^2V)h+2(DV)Dh.                                    \tag{13}
\]

By (1), |Dh|<=Mh. The last two terms in (13) have absolute value at most
M(1/4+M)h. Moreover Pi_e Y=0 by Haar integration. Set
N_2=sup |Y|/g_e, finite by compactness. In the time-four Duhamel formula,
replace K_4 by K_4-Pi_e only on the initial term. Use (9), (11), (12),
and keep the integral of exp(sM) rather than replacing it by 4 exp(4M):

\[
 N_2\le\big[(1+k_*)e^{4M}-1\big]N_2
       +r(\tfrac14+M)(e^{4M}-1).
                                                               \tag{14}
\]

The formula includes M=0 by its continuous limit. Whenever its denominator
is positive, division by h>=g_e/r and then differentiation of log h give

\[
 \|\operatorname{Hess}_{ee}\chi\|\le\delta_*(\theta):=
 \frac{e^{4\pi M}(e^{4M}-1)(1/4+M)}{2-(121/100)e^{4M}}+M^2.
                                                               \tag{15}
\]

Only quadratic forms Hess chi(D,D)=D^2 chi were used. Symmetry of the
covariant Hessian turns their bound into the operator-norm bound.

For a fully rational enclosure, put

\[
 E(x)=1+x+\frac{x^2}{2(1-x/3)},\qquad 0\le x<3.
                                                               \tag{16}
\]

Successive exponential-series ratios from degree two onward are at most
x/3, so exp(x)<=E(x), including the full tail. Use pi<22/7 and M=4 theta
to replace exp(4M) by X=E(16 theta) and exp(4 pi M) by
R_2=E((352/7)theta). The resulting upper bound delta_hat obeys

\[
 2-\frac{121}{100}X\ge\frac{102957}{177500}>0,
 \qquad \widehat\delta(1/100)
       =\frac{4081841509}{28120130625}<\frac3{20}.
                                                               \tag{17}
\]

The ratio delta_hat(theta)/theta is increasing: R_2, (1/4+4 theta),
(X-1)/theta=16+128 theta/(1-16 theta/3), and the reciprocal of the
positive denominator all increase; the remaining contribution is
16 theta. Consequently, including theta=0,

\[
 \boxed{\max_e\sup\|\operatorname{Hess}_{ee}\chi\|\le15\theta
        \quad(0\le\theta\le1/100).}                          \tag{18}
\]

## 4. T3 — Continuation controls the full mixed remainder

Write S_ef=Hess_ef chi, chi_1=(theta/3) sum W_p, and
E_ef=Hess_ef(chi-chi_1), e!=f. Define

\[
 \beta=\max_e\sum_{f\ne e}\sup\|S_{ef}\|,\qquad
 \varepsilon=\max_e\sum_{f\ne e}\sup\|E_{ef}\|.
                                                               \tag{19}
\]

YM88 and YM92 supply the all-coupling identities

\[
 \begin{split}
 \mathcal T_\chi S_{ef}&=\operatorname{Hess}_{ef}V+2\sum_gS_{eg}S_{gf},\\
 \mathcal T_\chi E_{ef}&=2\sum_gS_{eg}S_{gf}
                    +2\nabla_{\nabla\chi}\operatorname{Hess}_{ef}\chi_1,\\
 \mathcal T_\chi&=-\Delta^{\rm rough}-2\nabla_{\nabla\chi}+1.
 \end{split}                                                   \tag{20}
\]

Orthogonal parallel transport makes the drift connection diffusion
contractive in supremum fiber operator norm; Ricci damping one makes
the inverse norm at most one. The plaquette incidence row sum is twelve,
Hess_ef W_p has norm at most 1/4, and its derivative in each of the four
plaquette factors has norm at most 1/8. The repeated-factor case uses
Hess_ee W_p=-(W_p/4)g_e; ordinary scalar component derivatives would not
give the same formula. Therefore the general estimates, valid before
choosing any coupling window, are

\[
 \beta\le\theta+\varepsilon,\qquad
 \varepsilon\le4\delta\beta+2\beta^2+4\theta G.
                                                               \tag{21}
\]

Insert the newly proved G<=4 theta and delta<=15 theta. Then

\[
 0\le P_\theta(\varepsilon):=
 78\theta^2-(1-64\theta)\varepsilon+2\varepsilon^2.
                                                               \tag{22}
\]

This quadratic inequality admits a large branch, so it alone is not an
upper bound. Choose the constant barrier e_*=13/500. Direct exact
arithmetic and monotonicity in theta give, throughout the closed interval,

\[
 P_\theta(e_*)\le P_{1/100}(e_*)=-\frac{13}{62500}<0.
                                                               \tag{23}
\]

For each fixed finite volume epsilon(theta) is continuous: h depends
smoothly on theta, its minimum stays positive on compact parameter
intervals, and finite maxima and sums of uniform tensor norms are
continuous. At zero h is constant and epsilon(0)=0. It cannot cross e_*
because equality would contradict (22)--(23). Hence epsilon<e_* on the
whole interval. Substitution back into (22) yields

\[
 \varepsilon\le\frac{78\theta^2}{1-64\theta-2e_*}
 \le\frac{19500}{77}\theta^2<260\theta^2\quad(\theta>0).
                                                               \tag{24}
\]

At zero every remainder is zero. This proof uses finite-volume continuity
to select the branch and volume-independent constants to bound it. No
uniform spectral gap, truncated perturbation expansion or old-window
mixed-Hessian estimate enters this continuation argument.

## 5. T4 — Full gap and actual corrected memory on the enlarged interval

From (18), (21), (24), the symmetric block row bound on the full covariant
Hessian is

\[
 \|\operatorname{Hess}\chi\|\le h_*=16\theta+260\theta^2.
                                                               \tag{25}
\]

Thus Ric-2 Hess chi >= d_93 g. Weighted Bochner and the ground-state
transform prove (2) for the full finite configuration space. The same
argument works for every coordinate conditional, uniformly in every
exterior value. Gauge invariance of the unique positive vacuum means
the full gap lower bound restricts to the physical Gauss-invariant sector.
A physical Hamiltonian clock factor kappa multiplies d_93.

Here is how the blocking conclusions are justified without importing an
old-window floor. YM90--YM91's gauge-coordinate identities hold for any
smooth positive gauge-invariant h. At fixed coarse products Q, the actual
hidden conditional is nu_Q(z) proportional to h(section(Q),z)^2 Haar_Z.
Its centered generator A_Q has floor at least d=d_93 by (25). The full
hidden compression has floor d too. The coarse compression
C_b=I_b^* L I_b has floor d on centered coarse functions, since I_b is an
isometry and the newly proved full gap applies to I_b f. Thus no earlier
coarse numerical gap estimate is needed even at zero frequency.

For completeness, the actual-flux calculation uses exactly the following
quantities:

\[
 q=\theta+260\theta^2,\quad G=4\theta,\quad
 j=8\theta/9,\quad m=520\theta^2+2h_*j+(8/3)\theta G.
                                                               \tag{26}
\]

Every active fine plaquette has k=2 or 3 hidden links. Its leading source
r^1_(p,u)=(2 theta/3) sum_(e in p retained) grad_e W_p dot (u_R)_e
has hidden Haar eigenvalue 3k/4. The explicit trial flux is
J_u^0=grad_Z sum_p r^1_(p,u)/(3k/4). Prefix-frame transport has
|u_R|^2=b|u|^2 and is independent of z. The incidence and native derivative
calculations of YM92 give |J_u^0|<=sqrt(b) j |u|.

For the actual centered source r_u, retain the full weighted-divergence
residual, including the density drift:

\[
 e_u=r_u+\operatorname{div}_{\nu_Q}J_u^0
     =r_u-r_u^1+2\nabla_Z\chi\cdot J_u^0,
 \quad\nu_Qe_u=0,\quad |\nabla_Ze_u|\le\sqrt b\,m|u|.
                                                               \tag{27}
\]

The three contributions to m are the full mixed-response remainder,
Hess_ZZ chi acting on J^0, and the covariant derivative of J^0 contracted
with grad_Z chi. Conditional Poincare and its full positive inverse give

\[
 K_u=\nabla_ZA_Q^{-1}e_u,\qquad
 \|K_u\|_{L^2(\nu_Q)}\le\sqrt b\,(m/d)|u|.
                                                               \tag{28}
\]

The exact completed flux J^0+K satisfies -div_nu(J^0+K)=r_u. Its Gram
is at most b(j+m/d)^2 I. There is no cutoff in A_Q inverse and every
omitted source mode is included in (27)--(28). Since d+2h_*=1/2,

\[
 j+m/d=\frac{4\theta/9+(1592/3)\theta^2}{d}.
                                                               \tag{29}
\]

YM91's variational conditional-to-dynamic comparison now proves (3);
the factor b cancels against C_b[f]=b integral |grad f|^2 w_b dPhi.
In particular, the conclusion concerns the actual dynamic memory rather
than only a conditional surrogate tensor.

For s>=0, the additional actual source-variance estimate is
sigma^2=4q^2/d. Consequently

\[
 \Sigma_b(s)\le\alpha(s)C_b,\qquad
 \alpha(s)=\min\{\eta_{93},\sigma^2/(d+s)\},
 \quad
 0\le R_b(s)-(C_b+s)^{-1}
 \le\frac{\alpha(s)}{1-\alpha(s)}(C_b+s)^{-1}.             \tag{30}
\]

The s=0 statement is on centered coarse functions, with the floor proved
above. A stronger frequency-decay factor is not inferred just from a
zero-frequency bound.

## 6. Exact endpoints, evidence and remaining scope

| Coupling | Proved bare gap lower | Proved relative memory upper |
| --- | --- | --- |
| 1/400, YM92 | 251/800 = 0.31375 | 2105401/3189425625 < 1/1500 |
| 1/400, YM93 | 1667/4000 = 0.41675 | 635209/5627250225 < 1/8800 |
| 1/100, YM93 | 16/125 = 0.128 | (647/1440)^2 = 418609/2073600 < 1/4 |

At the same old endpoint the memory upper certificate improves by more
than 5.8. These are comparisons of proved bounds, not measurements of
the vacuum gap or memory. At the new endpoint the relative resolvent
error coefficient is 418609/1654991 < 1/3. The simpler constants were
chosen to give explicit closed-interval bounds, not an optimal endpoint.

The exact packet checks one-form and scalar partial-Bochner identities on
one, two and three native factors; detects omitted Ricci, drift and
connection terms; checks arbitrary-direction Casimir commutators and
second-derivative product rules; retains the full heat and exponential
tails; checks barrier signs, rational endpoints, and full residual Gram
controls. Engineered polynomial comparison potentials are explicitly
distinguished from the actual Wilson vacuum. The universal estimates
are analytic proofs, supported rather than replaced by finite controls.
They are not formal machine proofs or independent peer review.

Only (1) is asserted for arbitrary finite theta. In the inherited physical
normalization theta=c_B/(c_E g^4) tends to infinity as g tends to zero.
The diagonal, remainder and gap estimates have not been controlled on
that trajectory. A continuum Hilbert space, nonzero continuum states,
the continuum axioms and a uniform physical mass gap still need proofs.
The 4D continuum problem remains open in this programme.

Classical context: D. Bakry and M. Ledoux, *Levy--Gromov's isoperimetric
inequality for an infinite dimensional diffusion generator*, Invent.
Math. 123 (1996), 259--281, Section 2, for Bochner/diffusion calculus:
<https://www.math.univ-toulouse.fr/~ledoux/LevyGromov.pdf>.
The separate continuum target is the Jaffe--Witten official problem:
<https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf>.
Both consulted 6 October 2026. Neither source is cited as proving the
specific Wilson constants derived here.

Parent: YM92 commit `15f9062ec73d563b5fcb33d50135ec5b2f14994c`, tree
`a7a37917b2e0a201a79c97e07d604392137d81d4`. All 434 upstream inputs remain
byte-identical; the fresh result records 439 input hashes. Only the new
stage's --write command writes its result and pin. Default/--check are
read-only.

```sh
python -B papers/yang-mills-certified-benchmark/certificates/ym93_bochner_continuation.py --check
python -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym93_bochner_continuation.py -v
python -B papers/yang-mills-certified-benchmark/certificates/ym92_wilson_flux_corrector.py --check
python -B papers/yang-mills-certified-benchmark/certificates/ym88_spatial_score_contraction.py --check
```
