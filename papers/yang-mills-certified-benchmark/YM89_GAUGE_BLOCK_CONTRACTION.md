# YM89: curvature through marginalization and contraction of gauge blocks

YM89 develops the static coarse-graining step from YM88. Integrating out
arbitrarily many hidden links preserves a computable curvature comparison:
the updated matrix is a Schur complement. Gauge invariance then turns this
coordinate-marginal estimate into a bound for the actual ordered product
links, without treating their hidden coordinates as independent.

The resulting actual coarse-vacuum score decreases exponentially with
block length. For the inherited native Wilson Hamiltonian, at
0 <= theta <= 1/400 and every integer block length b >= 2,

\[
 \boxed{\sup|\nabla_E\log w_b|
       \le\frac{13}{110}\,2^{1-b}.}                     \tag{1}
\]

Here w_b is the actual coarse vacuum density with respect to product Haar,
and E is any coarse link. The bound is uniform in the volume and every
value of all other coarse links. In particular a conditional law on m
coarse links, with its whole exterior fixed, approaches product Haar with

\[
 \boxed{\operatorname{TV}(\mu_{b,S\mid S^c=z},\Phi_S)
       \le\min\{1,\tfrac{13\pi m}{220}\,2^{1-b}\}.}      \tag{2}
\]

This is an infrared blocking result in the specified strong gauge-coupling
regime. It neither extends that coupling window to weak gauge coupling
nor constructs a four-dimensional continuum theory. The exact reduced
dynamics retains a Schur-memory term, which is not bounded by a small
relative constant here. Section 8 identifies that remaining requirement
and gives a reversible finite comparison showing why it matters.

Base: YM88 commit `fcff60dc91f9dacd045f181c37e9e9716a0ff407`,
tree `3c8feafa24b66b230c2fc82dbc3bc389e1c423e6`. All earlier files remain
byte-identical. These are written analytic proofs with finite exact
controls, not formal machine verification or an independent review.

## 1. Model and the comparison matrix

Use the periodic three-dimensional cubic spatial lattice, native SU(2)
links with fundamental Casimir 3/4 and Ricci tensor (1/2)g, and

\[
 A_\theta=L_E-\theta\sum_p W_p,\quad
 A_\theta h=E_\theta h,\quad \chi=\log h,
 \quad d\mu=h^2d\Phi.                                   \tag{3}
\]

The finite-volume positive vacuum h is smooth and gauge invariant. For
the block results the fine side is Lb, with coarse side L >= 4. The
paths are the disjoint straight b-link paths of YM63. Interior vertices
of a retained path meet exactly two retained fine links. This geometric
condition is used in the gauge-fiber argument, not merely disjointness
of the edge sets.

Let W be YM88's nonnegative symmetric walk majorant for mixed log-Hessians.
Set J_ef = 2W_ef for e != f and J_ee = 0. YM88 proves

\[
 \rho=\tfrac12-70\theta,\quad
 r=10\theta,\quad r_w=20\theta,\quad
 \gamma=\rho-r=\tfrac12-80\theta\ge\tfrac3{10},           \tag{4}
\]

and the following uniform block bounds for the weighted Ricci tensor:

\[
 \operatorname{Ric}_e-2\operatorname{Hess}_{ee}\chi
       \ge\rho I,\qquad
 2\|\operatorname{Hess}_{ef}\chi\|\le J_{ef},\qquad
 \|J\|_\infty\le r,\quad \|J\|_w\le r_w.                \tag{5}
\]

The weighted row norm is ||B||_w = max_e sum_f 2^{d(e,f)} |B_ef|,
where d is the fine plaquette-link graph distance. It is submultiplicative.
Define the scalar comparison matrix C = rho I - J. Its off-diagonal
entries are nonpositive and C1 >= gamma 1. These component inequalities
imply the tensor lower bound Ric - 2 Hess chi >= gamma I.

The covariance and multiscale methods used below have established
antecedents. See [Menz, *A Brascamp--Lieb type covariance estimate*,
arXiv:1402.5160](https://arxiv.org/abs/1402.5160) for inverse-matrix
covariance comparisons, and [Grunewald, Otto, Villani and Westdickenberg,
*A two-scale approach to logarithmic Sobolev inequalities and the
hydrodynamic limit*](https://www.numdam.org/item/AIHPB_2009__45_2_302_0/)
for coarse-graining criteria. The proofs here specify the compact native
geometry, the actual Wilson-vacuum input and the gauge-product map;
those references are not used to assume a Wilson continuum estimate.

## 2. A conditional covariance comparison

We first state a version valid beyond the numerical window whenever its
geometric hypotheses can be established. On a finite product of compact
connected native factors with smooth positive density exp(2chi), suppose
a symmetric scalar matrix C has strictly positive row margin, nonpositive
off-diagonal entries, and satisfies

\[
 \operatorname{Ric}_i-2\operatorname{Hess}_{ii}\chi
       \ge C_{ii}I,\qquad
 2\|\operatorname{Hess}_{ij}\chi\|\le -C_{ij}\quad(i\ne j).
                                                               \tag{6}
\]

Then C is a positive definite M-matrix; its inverse is entrywise
nonnegative. Fix any complementary coordinates arbitrarily and let Z
be the coordinates still being integrated. For smooth scalar f,g on Z,
write a_i = sup |grad_i f| and b_i = sup |grad_i g|. The conditional law
nu obeys

\[
 \boxed{|\operatorname{Cov}_\nu(f,g)|
                   \le a^T C_{ZZ}^{-1}b.}               \tag{7}
\]

To prove this, let P_t be the conditional diffusion with generator
Delta_Z + 2 grad_Z chi dot grad_Z. Its differentiated equation contains
the weighted Ricci tensor. The componentwise Kato inequality and the
cooperative comparison system from (6) give

\[
 \big(\sup|\nabla_iP_tg|\big)_{i\in Z}
                       \le e^{-tC_{ZZ}}b.               \tag{8}
\]

One may regularize each component norm before applying the maximum
principle and then remove the regularization. The scalar diffusion
terms have nonpositive contribution at a component maximum; diagonal
curvature gives C_ii, and the off-diagonal terms are bounded by -C_ij.
All coefficients are smooth on a finite compact product. The matrix
exponential is entrywise nonnegative because -C_ZZ is cooperative.
Integration of the reversible identity

\[
 \operatorname{Cov}_\nu(f,g)
       =\int_0^\infty\nu\langle\nabla f,\nabla P_tg\rangle dt
                                                               \tag{9}
\]

now proves (7). Compact ergodicity justifies the endpoint at infinity;
the positive matrix bound makes the integral absolutely convergent.
No average over the coordinates held fixed enters (7).

## 3. Marginal curvature is controlled by a Schur complement

Let R be retained coordinates and Z hidden coordinates. Exact pins on
any other coordinates are allowed. Define

\[
 \chi_R(x)=\tfrac12\log\int\exp(2\chi(x,z))d\Phi_Z(z).
                                                               \tag{10}
\]

A constant normalization does not affect derivatives. Differentiating
the normalized conditional integral, with covariant derivatives on R,
gives the exact identity

\[
 \operatorname{Hess}_{ef}\chi_R
 =\nu_x(\operatorname{Hess}_{ef}\chi)
       +2\operatorname{Cov}_{\nu_x}
                    (\nabla_e\chi,\nabla_f\chi).         \tag{11}
\]

Apply (7) to scalar components of these two gradients. Their Z-gradient
norms are bounded by -C_ez/2. Hence the covariance contribution in twice
(11) is bounded by the entries of C_RZ C_ZZ^(-1) C_ZR. This proves (6)
for the marginal, with comparison matrix

\[
 \boxed{C^{\rm eff}
       =C_{RR}-C_{RZ}C_{ZZ}^{-1}C_{ZR}.}                 \tag{12}
\]

For the diagonal this means a lower bound on Ric_e - 2 Hess_ee chi_R;
for distinct links it means 2||Hess_ef chi_R|| <= -C_eff,ef.
The matrix is an upper/lower comparison as specified, not an equality
with the actual marginal Hessian. Equality holds for the appropriate
Gaussian comparison, but is not assumed for the Wilson vacuum.

Both uniform positivity and a row margin survive. Write C_RZ = -J_RZ.
From C1 >= gamma 1 and positivity of C_ZZ^(-1),

\[
 C_{ZZ}^{-1}J_{ZR}1_R
       \le1_Z-\gamma C_{ZZ}^{-1}1_Z,
 \quad
 C^{\rm eff}1_R
       \ge\gamma(1_R+J_{RZ}C_{ZZ}^{-1}1_Z)\ge\gamma1_R.
                                                               \tag{13}
\]

The Schur variational formula also gives C_eff >= gamma I. Taking a
principal submatrix, corresponding to further exact coordinate pins,
only improves the row margin. Repeated elimination is associative:
eliminating Z1 and then Z2 produces the same comparison matrix as
eliminating their union. These facts apply to arbitrary finite sets;
there is no small joint-density denominator or growing pin-count factor.

For the actual YM88 input the inverse in (12) has the convergent full
walk expansion

\[
 C_{ZZ}^{-1}=\rho^{-1}\sum_{n\ge0}(J_{ZZ}/\rho)^n.        \tag{14}
\]

After degree m, the omitted weighted norm is at most
rho^(-1)(r_w/rho)^(m+1)/(1-r_w/rho). At the endpoint r_w/rho <= 2/13.
Together with YM88's certified tail for J, this specifies full-tail
computation of the static comparison; finite walk sums alone are not
upper certificates.

Extending the submatrices by zero and using row norms bounds the
off-diagonal marginal log-Hessian uniformly by

\[
 \begin{split}
 \max_e\sum_{f\ne e}\sup\|\operatorname{Hess}_{ef}\chi_R\|
       &\le T:=\frac{\rho r}{2(\rho-r)}\le\frac{13}{960},\\
 \max_e\sum_{f\ne e}2^{d(e,f)}\sup\|\operatorname{Hess}_{ef}\chi_R\|
       &\le T_w:=\frac{\rho r_w}{2(\rho-r_w)}\le\frac{13}{440}.
 \end{split}                                                \tag{15}
\]

Indeed the relevant nonnegative matrix is J_RR plus the fill term
J_RZ C_ZZ^(-1) J_ZR. Its row norms are at most r+r^2/(rho-r)
and r_w+r_w^2/(rho-r_w), respectively. Its diagonal may be positive;
discarding it for the off-diagonal sums only reduces these bounds.

## 4. Gauge symmetry identifies the actual product-link density

Now retain all fine links on the disjoint straight paths, integrating
out all other fine links without fine-coordinate pins. Let w_R be this
retained density. Gauge invariance of the original vacuum and invariance
of Haar integration make w_R gauge invariant at every retained vertex.

At an internal vertex of a retained path only its two neighboring retained
links are transformed. Gauge transformations there leave the ordered
product Q = q_1 ... q_b unchanged. Their action is transitive on each
fiber of this product map. Explicitly, with the endpoint gauges fixed
to one, choose successive internal gauges g_(i+1) = g_i q_i. The
transformed path is (1,...,1,Q). Thus w_R depends on that path only
through Q. The argument works simultaneously on every retained path.

The product coordinate map preserves product Haar after keeping b-1
factors as hidden variables. Since their density is constant along each
gauge fiber, integration gives exactly

\[
 \boxed{w_R=w_b\circ\pi_b,\qquad
        \chi_R=\psi_b\circ\pi_b,\quad\psi_b=\tfrac12\log w_b.}
                                                               \tag{16}
\]

This is the actual product-link marginal. It is not a fitted Wilson
action on the coarse lattice. Arbitrary fine-coordinate pins can break
the internal gauge invariance and are not inserted into this argument.
After (16) has been established, all ensuing pointwise coarse derivative
bounds remain uniform when the other coarse products are pinned.

## 5. A separated derivative sees the coarse score

For distinct coarse paths E,F, choose any fine links e in E and f in F.
Transport fine derivative directions through the fixed prefixes of
their ordered products. These adjoint maps are orthogonal. The two
different products then give

\[
 \|\operatorname{Hess}_{ef}\chi_R(q)\|
       =\|\operatorname{Hess}_{EF}\psi_b(\pi_bq)\|.       \tag{17}
\]

For two distinct fine links on the same path, choose their transported
directions to be native coarse frame directions u_a,u_b. If e precedes
f, the resulting cross-derivative matrix is

\[
 B_{ab}=D_bD_a\psi_b
       =(\operatorname{Hess}_{EE}\psi_b)_{ab}
           -\tfrac12\epsilon_{abc}D_c\psi_b.             \tag{18}
\]

Reversing path order transposes this expression. One way to verify the
ordering is to vary the two fine factors: the coarse product becomes
exp(-t u_a/2) exp(-s u_b/2) Q. Native fields satisfy
[D_a,D_b] = epsilon_abc D_c. The symmetric part of B is the covariant
coarse Hessian. Its antisymmetric part is one half the cross-product
matrix of the coarse gradient, whose operator norm is |grad psi_b|/2.
Consequently

\[
 \|\operatorname{Hess}_{EE}\psi_b\|
      \le\|\operatorname{Hess}_{ef}\chi_R\|,\qquad
 |\nabla_E\psi_b|\le2\|\operatorname{Hess}_{ef}\chi_R\|.
                                                               \tag{19}
\]

In particular, the score is recovered from the antisymmetric part of a
mixed fine response, not set to zero by treating that matrix as a scalar
or a symmetric coarse Hessian. Exact polynomial controls test (17)--(18)
for noncommuting rational quaternion configurations with b = 2 and b = 3.

Choose the first and last fine links on E. Their distance is at least
b-1: a plaquette-link graph step changes the midpoint's torus Manhattan
distance by at most one. The fine side Lb >= 4b prevents a shorter
wrap around this segment. From (15) and (19),

\[
 \sup\|\operatorname{Hess}_{EE}\psi_b\|\le T_w2^{1-b},
 \quad \sup|\nabla_E\psi_b|\le2T_w2^{1-b},
 \quad \sup|\nabla_E\log w_b|\le4T_w2^{1-b}.             \tag{20}
\]

The endpoint value T_w <= 13/440 proves (1).

## 6. Row bounds and the coarse diagonal generator

Write Delta_E = sup ||Hess_EE psi_b|| and
H^b_EF = sup ||Hess_EF psi_b|| for E != F. For a fixed fine e in E,
sum (17) and (19) over every retained f != e. Equation (15) yields

\[
 (b-1)\Delta_E+b\sum_{F\ne E}H^b_{EF}\le T.              \tag{21}
\]

There is also an exponential row estimate. For any retained straight
path F, the number of its fine links within distance strictly less
than an integer R of e is at most 2R-1. To see this, project midpoints
onto the path's axis. Distance < R implies coordinate distance <= R-1;
an interval of that length contains at most 2R-1 of its unit-spaced
midpoints. The same bound holds on the fine torus in the range used here.
Thus, for 1 <= R <= ceil(b/2),

\[
 \sum_{f\in F}2^{d(e,f)}\ge(b-2R+1)2^R,
 \qquad M_b=\max_{1\le R\le\lceil b/2\rceil}(b-2R+1)2^R.
                                                               \tag{22}
\]

Using (17) in the weighted row bound of (15),
sum_(F != E) H^b_EF <= T_w/M_b. Combining this with (20)--(21) gives

\[
 \max_E\left(\Delta_E+\sum_{F\ne E}H^b_{EF}\right)
 \le\varepsilon_b:=\min\left\{\frac{T}{b-1},
              T_w\left(2^{1-b}+M_b^{-1}\right)\right\}.
                                                               \tag{23}
\]

Here M_2 = 2; for even b = 2k >= 4, M_b = 3 times 2^(k-1);
for odd b = 2k+1 >= 3, M_b = 2^(k+1). Thus the second bound is
O(2^(-b/2)), without a volume factor. The first bound is sharper for
some small b. Both are derived from the same actual marginal.

Block Schur bounds and the native Bochner identity now give

\[
 \operatorname{gap}\big(-w_b^{-1}\operatorname{div}(w_b\nabla)\big)
       \ge\lambda_b:=\tfrac12-2\varepsilon_b
       \ge\frac{227}{480}.                              \tag{24}
\]

The same floor holds for every conditional law on any subset of coarse
links with the other coarse products fixed. At theta = 0 the estimates
give lambda_b = 1/2; the exact unrestricted free native gap is 3/4.
The floor in (24) is conservative, not an assertion that the limiting
gap equals 1/2.

YM64's exact pullback energy identity gives the compressed diagonal
form C_b = b times the generator in (24). Therefore its centered gap
is at least b lambda_b. This is a statement about the diagonal form;
it is not the gap of the exact reduced interacting dynamics.

## 7. Uniform conditional approach to Haar

The native factor has diameter 2pi. The normalized density of any m
target coarse links, given all other coarse products, has log oscillation
at most

\[
 \Omega_{m,b}=2\pi m\,4T_w2^{1-b}.                       \tag{25}
\]

Its derivative in a target link is the same as that of log w_b because
the conditional normalizer is constant in the target variables. A
normalized positive density of log oscillation at most Omega has
total variation from its reference measure at most tanh(Omega/4).
For example, bounding the convex function (w-1)_+ by its chord between
the density's minimum and maximum reduces the bound to maximizing
(B-1)(1-A)/(B-A) under B/A <= exp(Omega); the maximum is attained at
A = exp(-Omega/2), B = exp(Omega/2). Since tanh(x) <= x for x >= 0,

\[
 \operatorname{TV}\le\min\{1,\Omega_{m,b}/4\}
       =\min\{1,2\pi mT_w2^{1-b}\}.                    \tag{26}
\]

This proves (2). At b = 8 its endpoint bound for m = 1 is at most
13/8960, using pi < 22/7. It applies to the worst coarse exterior;
unconditioned one-link marginals are already Haar by gauge invariance.
For fixed m the bound tends to zero along any sequence of admissible
finite volumes as b grows. Averaging over some exterior products also
preserves it. No assertion about existence of an infinite-volume or
continuum vacuum is needed to state this finite-volume uniform limit.

## 8. The dynamic memory requirement survives blocking

The static comparison matrix (12) and the operator Schur complement of
the Hamiltonian serve different purposes. For YM64's vacuum-adapted
projection the exact reduced resolvent has the form

\[
 I_b^*(\mathcal L+s)^{-1}I_b
       =(C_b-\Sigma_b(s)+s)^{-1}.                       \tag{27}
\]

The actual marginal determines C_b. Equations (1)--(26) do not establish
a small relative bound on Sigma_b. YM88's gamma upgrades the inherited
hidden-sector floor D_b >= gamma in the enlarged coupling window, but
the resulting general bound Sigma_b(0) <= C_b alone is insufficient
for transferring the improved diagonal gap to the full dynamics.

A reversible finite example makes this failure explicit. On three states
with uniform stationary law take conductances a > 0 on 0--1 and
epsilon > 0 on 1--2. Project by actual conditional expectation onto the
partition {{0,2},{1}}. On the centered retained and hidden directions,

\[
 C=\tfrac32(a+\epsilon),\quad D=\tfrac12(a+\epsilon),\quad
 |B|^2=\tfrac34(a-\epsilon)^2,\quad
 C-|B|^2/D=\frac{6a\epsilon}{a+\epsilon}.                \tag{28}
\]

The stationary coarse law remains (2/3,1/3). The nonzero full eigenvalues
have sum 2(a+epsilon) and product 3a epsilon. With a = b and
epsilon = b^(-2), both C and D grow, while the full gap is bounded above
by 3a epsilon/(a+epsilon) and tends to zero. Its relative memory is

\[
 \eta=\frac{|B|^2}{CD}
           =\left(\frac{a-\epsilon}{a+\epsilon}\right)^2\longrightarrow1.
                                                               \tag{29}
\]

This comparison is an actual reversible Markov chain, not a Wilson
vacuum or a continuum counterexample. It demonstrates that the coarse
diagonal and hidden-sector floors alone do not control their coupling.

A sufficient quantitative replacement can be stated precisely. On the
centered spaces suppose C_b >= c I, D_b >= d I and
Sigma_b(0) <= eta C_b with 0 <= eta < 1. Energy Cauchy--Schwarz bounds
the off-diagonal form by sqrt(eta) times the two diagonal energy norms.
The elementary two-by-two comparison then gives

\[
 \operatorname{gap}(\mathcal L)
 \ge\frac{c+d-\sqrt{(c-d)^2+4\eta cd}}2
 \ge\frac{cd(1-\eta)}{c+d}.                             \tag{30}
\]

For clarity, if the two diagonal energies are x,y, then the total form
is at least x+y-2sqrt(eta xy), while the vector norm is at most x/c+y/d.
Minimizing this ratio proves the first inequality; determinant divided
by trace proves the second. Domains are the vacuum-adapted form domains
of YM64. The retained centered part and every hidden vector are orthogonal
to the vacuum. No bounded L2 off-diagonal operator is assumed.

To use this route toward weak gauge coupling one still needs verified
actual bounds for the block coefficients and the relative memory along
that scaling trajectory, and a nontrivial continuum-state construction.
The static convergence toward Haar in (2) is infrared behavior at the
declared small theta. It is not a UV continuum limit or a shortcut around
the memory bound.

## 9. Reproduction

Run with Python 3.12, without optimization:

```sh
python -B papers/yang-mills-certified-benchmark/certificates/ym89_gauge_block_contraction.py --check
python -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym89_gauge_block_contraction.py -v
```

The manifest pins 396 upstream files. The result records those, four
new local inputs and the manifest, for 401 input hashes. Default and
--check are read-only; --write explicitly regenerates only the YM89
result and digest. The dedicated workflow also replays YM88 and the
vacuum-adapted blocking evidence of YM64.

Finite controls check static Schur associativity and fill signs, Gaussian
marginal calibration, all hidden-walk tails, quaternion prefix rotations,
the symmetric and antisymmetric same-path derivative parts, exact gauge
fiber flattening, native retained-network incidence, path-ball counts,
conditional total variation, and the reversible memory obstruction.
They test constants and mechanisms supporting the written proof; none
is passed off as a simulation or a certified discretization of the
interacting Wilson vacuum.
