# YM88: spatial score contraction for the actual small-coupling vacuum

YM88 closes YM87's conditional energy-transfer gate for the original finite
Wilson Hamiltonian in the explicit regime **0 <= theta <= 1/400**. The main
new program result is a pointwise, spatially summable bound on the actual
vacuum's mixed log-Hessian. It is uniform in the entire exterior configuration
and the total number of links. It gives both a constructive walk majorant for
conditional score coefficients and an untruncated, volume-independent gap:

\[
 \boxed{\operatorname{gap}(A_\theta-E_\theta)
       \ge \frac12-80\theta\ge\frac3{10}.}                 \tag{1}
\]

This extends the endpoint of the earlier program's strict small-coupling
window, theta < 1/8000, by a factor of twenty. It is a sufficient regime,
not a sharp threshold. Small theta is the strong gauge-coupling regime in
the inherited Hamiltonian convention. The large-theta weak gauge-coupling
estimates and the four-dimensional continuum construction remain open.
No novelty claim about the general methods of strong-coupling analysis,
Bochner estimates, matrix comparison, or walk expansions is made.

Base: YM87 commit `b7e7cf5e9d5e92b2bd6a6071a0e23f427475a800`,
tree `213e7b9e30887235ef56c25fa83ade6f6b5b27eb`.
All earlier evidence files remain byte-identical. The analytic proofs below
are accompanied by exact finite controls; they are not formally machine
verified proofs, and no numerical approximation of the Wilson ground state
is substituted for that ground state.

## 1. Native model and statement

Let the spatial lattice be the periodic cubic lattice of side L >= 4.
Let E denote its positively oriented links, P its elementary plaquettes,
and Phi the product normalized Haar measure on SU(2)^E. As in YM62--YM87,

\[
 A_\theta=L_E-V_\theta,\quad
 L_E=-\sum_{e,a}D_{e,a}^2=-\Delta,\quad
 V_\theta=\theta\sum_{p\in P}W_p,\qquad
 W_p=\tfrac12\operatorname{Tr}U_p.                         \tag{2}
\]

The ordered Wilson word contains each of its four links once. A unit
native derivative inserts a quaternion generator with norm 1/2.
Every link meets four plaquettes. The native SU(2) metric is the round
three-sphere of radius two: its fundamental Casimir is 3/4 and its Ricci
tensor is (1/2)g. Use this metric throughout, including in operator norms.

On every finite volume, positivity improvement gives the unique positive
normalized ground state h_theta, eigenvalue E_theta, and actual probability
measure mu_theta = h_theta^2 Phi. Write chi = log h_theta. These are full
compact-manifold objects. The ground state is gauge invariant, so the same
gap lower bound holds after restriction to the Gauss-invariant sector.
For the physical Hamiltonian H = kappa A_theta plus an additive constant,
the lower bound in (1) is multiplied by kappa > 0.

Join two distinct links when they belong to a common plaquette; d(e,f)
is distance in this graph. Define

\[
 S_{ef}=\operatorname{Hess}_{ef}\chi,\qquad
 H_{ef}=\sup_{q\in SU(2)^E}\|S_{ef}(q)\|_{\rm op}.
                                                               \tag{3}
\]

In the specified window the following bounds hold:

\[
 H_{ee}\le35\theta,\qquad
 \max_e\sum_{f\ne e}2^{d(e,f)}H_{ef}\le10\theta,\qquad
 \max_e\sum_{f\ne e}H_{ef}\le5\theta.                     \tag{4}
\]

The last factor of two matters: every off-diagonal distance is at least
one. All suprema in (3) include the entire exterior, not a marginal or an
average over exterior configurations.

## 2. Extend the one-link derivative contraction

This section re-establishes the local part of YM67 at a larger coupling.
It does not invoke YM67's software outside its declared window and does
not assume the earlier small-coupling spectral gap.

Set M = 4 theta. Split A = L_e + A_rest - V_e, where V_e collects exactly
the plaquettes meeting e, |V_e| <= M, and A_rest is independent of q_e.
Set g_e(z) = integral h(x,z) dPhi(x). The full native time-nine heat bound
from YM53 is 19/20 <= k_9(x,y) <= 21/20. Positive Feynman--Kac comparison
then gives, pointwise,

\[
 r_0^{-1}g_e\le h\le r_0g_e,
 \quad r_0=\frac{21}{19}e^{18M},\qquad
 e^{sE}e^{-sA_{\rm rest}}g_e\le e^{sM}g_e.                \tag{5}
\]

For completeness, the last inequality follows by integrating Ah = Eh
over q_e: (A_rest - U_e)g_e = Eg_e with
U_e = integral V_e h / g_e and |U_e| <= M. These comparisons require
only finite volume and bounded V_e, not small coupling.

Every native invariant derivative commutes with L_E. For Y = D_e h,
its reference average in q_e vanishes and
(L_e + A_rest - E)Y = V_e Y + (D_e V)h. In the Duhamel formula over
time nine, replace K_9 by K_9 - Pi_e on the initial term. Its supremum
norm is at most 1/20. Since |D_e V| <= M/2, (5) gives

\[
 N\le c_0N+\frac92e^{9M}M r_0,
 \quad N=\sup\frac{|Y|}{g_e},\quad
 c_0=e^{9M}\left(\frac1{20}+9M\right).                    \tag{6}
\]

Smoothness and strict positivity make N finite before this estimate is
used. If c_0 < 1, solving (6) and dividing by h >= g_e/r_0 proves

\[
 |\nabla_e\chi|\le B,
 \qquad B=\frac{(9/2)M(21/19)^2e^{45M}}{1-c_0}.           \tag{7}
\]

The bound holds for every unit direction, so B bounds the entire
three-component gradient, without a factor sqrt(3).
For any ordered pair of unit native derivatives i,j,

\[
 (A-E)D_iD_jh=(D_iD_jV)h+(D_jV)D_ih+(D_iV)D_jh.
                                                               \tag{8}
\]

Its reference average in the link of i again vanishes. Using
|D_iD_j V| <= M/4 and repeating (6) bounds |D_iD_j h|/h by
B_2 = B/2 + 2B^2. On a single SU(2) factor, an invariant unit vector
field u has nabla_u u = 0. Consequently the quadratic form of the
covariant Hessian satisfies

\[
 |\operatorname{Hess}_{ee}\chi(u,u)|
   =|D_u^2\chi|\le B_2+B^2=B/2+3B^2.                    \tag{9}
\]

The covariant Hessian is symmetric, so this is its operator-norm bound.
One must not identify its off-diagonal frame components with unordered
second Lie derivatives; (9) only uses the correct quadratic form.

Here M <= 1/100. The elementary positive-series bounds

\[
 e^x\le\frac1{1-x}\ (0\le x<1),\qquad
 e^x\le1+x+\frac{x^2}{2(1-x/3)}\ (0\le x<3)             \tag{10}
\]

give e^(9/100) <= 100/91 < 11/10 and
e^(45/100) <= 1067/680 < 8/5. In the second estimate, the ratios of
successive Taylor terms starting with x^2/2 are at most x/3.
Thus c_0 <= 77/500 < 1 and, at the endpoint,
B <= 1764/16967 < 21/200. The quotient B/M increases with M:
its numerator increases and its positive denominator decreases. Hence

\[
 B\le42\theta,\qquad
 B/2+3B^2\le21\theta+5292\theta^2
          \le\frac{3423}{100}\theta<35\theta.             \tag{11}
\]

At theta = 0 the bounds vanish exactly. This proves the diagonal part
of (4) throughout the closed enlarged window.

## 3. The mixed covariant Hessian equation

The ground equation, divided by h, is

\[
 \Delta\chi+|\nabla\chi|^2+V=-E.                         \tag{12}
\]

For e != f, S_ef is a section of T_e^* tensor T_f^*. Denote its connection
Laplacian by Delta^rough. Product curvature and Ric = (1/2)g on each
factor give

\[
 \operatorname{Hess}_{ef}\Delta\chi
       =(\Delta^{\rm rough}-1)S_{ef}.                   \tag{13}
\]

There is one Ricci contribution from each of the two different factors.
In Hess_ef |grad chi|^2 the additional Riemann-curvature terms vanish:
their indices would need to lie in a single factor while e != f.
Differentiating (12) therefore gives the exact tensor identity

\[
 \boxed{(-\Delta^{\rm rough}-2\nabla_{\nabla\chi}+1)S_{ef}
       =\operatorname{Hess}_{ef}V
        +2\sum_g S_{eg}S_{gf}.}                         \tag{14}
\]

Compositions in (14) contract the tangent index of g. The tensors in
the diagonal terms are covariant Hessians as well. Equation (14) is not
an equation for nine scalar components with ordinary scalar Laplacians:
the connection terms are essential. Exact polynomial controls on two
and three native SU(2) factors check all 36 mixed components, including
terms mediated by a third factor. Dropping the Ricci damping or replacing
the rough Laplacian by componentwise scalar Laplacians changes these
controls. The universal justification is the covariant calculation above;
these finite controls check signs and normalization.

The connection diffusion with generator
Delta^rough + 2 nabla_(grad chi) is a contraction in the supremum of
the fiber operator norm. Indeed its stochastic parallel transport is
orthogonal in both factors and preserves the factor splitting; expectation
cannot increase that norm. The same statement follows from the tensor
maximum principle. Compactness makes this diffusion well defined for the
smooth drift at each finite volume. Integrating the semigroup against
e^(-t) gives a resolvent of norm at most one, with no volume constant.

Let N_ef count plaquettes containing e and f, with N_ee = 0, and set
C_ef = theta N_ef/4. A double quaternion insertion bounds
sup ||Hess_ef V||_op <= C_ef. Every row of N sums to twelve: four
incident plaquettes, with three other links each, including multiplicities.
Taking the above norm in (14), and using H_ee <= delta = 35 theta, gives

\[
 (1-4\delta)H_{ef}
 \le C_{ef}+2\sum_{g\ne e,f}H_{eg}H_{gf}\quad(e\ne f).   \tag{15}
\]

## 4. Continuation selects the small branch

For matrices indexed by links define the weighted row norm

\[
 \|X\|_2=\max_e\sum_f2^{d(e,f)}|X_{ef}|.                \tag{16}
\]

The subscript denotes the spatial weight, not a Hilbert-space operator
norm. The triangle inequality for d proves submultiplicativity. Write
J_ef = H_ef off the diagonal and J_ee = 0, and u = ||J||_2. Since C
has support at graph distance one, ||C||_2 = 6 theta. Equation (15)
implies

\[
 u\le6\theta+140\theta u+2u^2.                           \tag{17}
\]

**The quadratic inequality alone does not bound u by its smaller root.**
Its large branch is admissible too. For each fixed finite volume, the
simple positive ground eigenfunction depends continuously on theta in
every required smooth norm. This follows from bounded smooth parameter
perturbation of the simple elliptic eigenvalue and elliptic regularity;
only a local, finite-volume isolated eigenvalue is needed. Thus u(theta)
is continuous, and u(0) = 0. The finite diameter may depend on volume;
it is used only to establish this finite-volume continuity.

Throughout theta in [0,1/400], (17) implies

\[
 0\le\frac3{200}-\frac{13}{20}u+2u^2
     =2\left(u-\frac1{40}\right)\left(u-\frac3{10}\right).
                                                               \tag{18}
\]

The interval (1/40,3/10) is forbidden. Continuity from zero therefore
forces u <= 1/40 throughout the interval. Substituting this into (17)
gives 140 theta + 2u <= 2/5, so (3/5)u <= 6 theta, or u <= 10 theta.
No passage between the two roots is possible. All constants are independent
of volume, completing (4).

## 5. A computable walk majorant with a full tail

Let c = (1 - 140 theta)^(-1), b = 10 theta. In the weighted-norm ball
of radius b, the map F(X) = cC + 2cX^2 is a contraction:

\[
 c(6\theta+2b^2)\le b,\qquad
 \|F(X)-F(Y)\|_2\le4cb\|X-Y\|_2,
 \quad4cb\le\frac2{13}.                                 \tag{19}
\]

It maps the nonnegative ball to itself. Equation (15) says J <= F(J)
entrywise, including the trivial zero diagonal of J. By (4), J is in
the ball. Iteration from J increases entrywise and converges to the
unique fixed point W. Therefore J <= W. This step uses contraction and
the previously established small-branch bound; it does not infer an
upper bound merely from a nonlinear supersolution inequality.

Iteration from zero gives increasing lower approximations W_n with
||W-W_n||_2 <= b(4cb)^n. A spatially organized form is the Catalan series

\[
 W=\sum_{n\ge1}\mathrm{Cat}_{n-1}\,2^{n-1}c^{2n-1}C^n,
 \qquad \mathrm{Cat}_m=\frac1{m+1}{2m\choose m}.           \tag{20}
\]

All entries are nonnegative, and W is symmetric. Its diagonal is generally
positive because walks can return; it must not be silently set to zero
inside the quadratic equation. With lambda = ||C||_2 = 6 theta and
q = 8c^2 lambda = 48c^2 theta <= 48/169, the first m terms obey

\[
 \|W-W^{[m]}\|_2\le\frac{c\lambda q^m}{1-q},\qquad m\ge0.
                                                               \tag{21}
\]

This follows from Cat_(n-1) <= 4^(n-1). At the endpoint the prefactor
is 39/1210. For theta = 0 the entire series and tail vanish. Each C^n
is supported at d <= n, so (20)--(21) are a local walk computation with
a bound on every omitted walk. A truncation alone is a lower approximation
to W, not an upper certificate for the actual response. A norm tail must
also not be assigned independently to every entry and then summed without
accounting for that sum.

## 6. Actual conditional score coefficients

Fix all links other than e arbitrarily. The conditional law is

\[
 \mu_e(dx\mid z)=\frac{h(x,z)^2}{\int h(y,z)^2d\Phi(y)}d\Phi(x).
                                                               \tag{22}
\]

Its weighted Ricci tensor is Ric_e - 2 Hess_ee chi. Equations (11) and
the one-factor integrated Bochner identity give the uniform conditional
Poincare floor

\[
 \rho=\frac12-70\theta\ge\frac{13}{40}.                  \tag{23}
\]

To spell out this step, for the symmetric conditional generator
L_e^mu = Delta_e + 2 grad_e chi dot grad_e, integration of Bochner gives
||L_e^mu f||_2^2 >= rho integral |grad_e f|^2. Its spectrum is discrete
on a compact connected factor and only constants have zero energy.
Applying the inequality to each nonconstant eigenfunction gives (23).
It does not assume a global interacting spectral gap.

For a unit direction v on f != e, differentiating the normalization in
(22) yields the centered score

\[
 s_v=2\big(D_{f,v}\chi-\mu_e(D_{f,v}\chi)\big),\quad
 \nabla_e s_v=2S_{ef}v.                                 \tag{24}
\]

The second term in s_v is constant in the integration variable q_e.
Conditional Poincare first bounds mu_e(s_v^2) by 4H_ef^2/rho; its
inverse-generator bound then gives the energy-dual score estimate

\[
 \langle s_v,(-L_e^\mu)^{-1}s_v\rangle_{\mu_e}
       \le\frac{4H_{ef}^2}{\rho^2}
       \le\left(\frac{2W_{ef}}\rho\right)^2.              \tag{25}
\]

This holds simultaneously for every unit v and every full exterior z.
Consequently the actual YM87 energy-transfer comparison can use

\[
 K_{fe}=\frac{2W_{ef}}\rho\ (e\ne f),\quad K_{ee}=0,
 \qquad
 \max_e\sum_{f\ne e}2^{d(e,f)}K_{fe}
       \le\frac{20\theta}\rho\le\frac2{13}.              \tag{26}
\]

In particular, for integer R >= 1,

\[
 \sup_e\sum_{d(e,f)\ge R}K_{fe}
       \le\frac{20\theta}\rho\,2^{-R},\qquad
 r:=\max_e\sum_fK_{fe}\le\frac{10\theta}\rho\le\frac1{13}.
                                                               \tag{27}
\]

Symmetry gives the same row and column bounds. The diagonal of W is
discarded only when defining the update matrix K. Equations (25)--(27)
discharge actual small-theta response coefficients under full conditioning,
not merely coefficients of a comparison density or an averaged exterior.

## 7. Close the sweep gate and check the gap independently

Define the finite positive definite comparison metric

\[
 G=\rho^{-1}(I-K)^{-1},\qquad
 \|G\|_{\ell^2\to\ell^2}\le\frac1{\rho(1-r)}
                   \le\frac1{1/2-80\theta}.             \tag{28}
\]

For the YM87 update on link i, M_i = I - (I-K)e_i e_i^T, direct
multiplication using K_ii = 0 gives

\[
 G-M_i^TG M_i=\rho^{-1}e_i e_i^T.                        \tag{29}
\]

Telescoping over any permutation of links yields exactly the YM87 sweep
Lyapunov identity G - M_sweep^T G M_sweep = D_sweep. Thus (28) closes
that criterion and gives (1), with a uniform constant and every link
in the sweep. No original-support preservation is assumed for conditional
projections.

There is also an independent direct proof. The symmetric block Hessian
has diagonal norm at most 35 theta and off-diagonal row and column sums
at most 5 theta, by (4). The block Schur bound gives
Hess chi <= 40 theta I. Therefore on the entire product manifold

\[
 \operatorname{Ric}-2\operatorname{Hess}\chi
                  \ge(1/2-80\theta)I.                   \tag{30}
\]

The integrated Bochner identity for L_mu = Delta + 2 grad chi dot grad
now gives its gap at least 1/2 - 80 theta, by the same compact spectral
argument as in Section 6. Ground-state transformation identifies -L_mu
with A-E. This proves (1) independently of the sweep criterion. The
bound is conservative: at theta = 0 the exact unrestricted native gap
is 3/4, whereas (30) gives 1/2.

## 8. Reproduction and limits

Run with Python 3.12, without optimization:

```sh
python -B papers/yang-mills-certified-benchmark/certificates/ym88_spatial_score_contraction.py --check
python -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym88_spatial_score_contraction.py -v
```

The dedicated workflow also replays the pinned YM87 conditional-energy,
YM83 actual-vacuum, and YM67 local-derivative evidence. The manifest pins
389 upstream files; the result records those, four local inputs, and the
manifest, for 394 input hashes. Default and --check never rewrite evidence.
The --write mode explicitly regenerates only the YM88 result and digest.

Controls cover native curvature and Lie brackets, the mixed tensor identity,
independent exponential enclosures, the two quadratic branches, exact
periodic plaquette incidence, walk support and returns, Catalan convolution
and remainder bounds, and every order of a finite symmetric sweep metric.
The polynomial parent densities and finite matrices are labeled comparison
controls. They do not approximate or discretize the actual Wilson vacuum;
the passage to it is the analytic argument in Sections 2--7.

The proven statement concerns the full finite spatial lattice Hamiltonian,
uniformly in its volume, at the declared theta. Uniform finite-volume
estimates do not by themselves construct a continuum vacuum or show that
nonzero states and their normalization survive a joint scaling limit.
Nor may the small-theta bound be continued to large theta through the
forbidden-root argument: its constants cease to establish the required
barrier outside the specified interval. A fixed-lattice bare lower bound
and its clock multiplier are not a four-dimensional continuum mass gap.

The next genuine obstacle is a full-exterior summability estimate that
remains effective along the weak gauge-coupling scaling trajectory, together
with a nontrivial continuum state construction. YM88 resolves this gate in
an explicit enlarged strong-coupling regime only.
