# YM87 — Conditional energy transport and a block-sweep gap certificate

Monty Dabas. 6 October 2026. Python 3.12 evidence packet.

## What is established

YM86 controls omitted response on a declared energy class. To use such a
class after conditional averaging, one must also control the energy of
the conditional projections. This stage supplies that missing type of
estimate and a conditional route from it to a full quantum gap.

1. In the actual finite Wilson vacuum, each native one-link logarithmic
   gradient has an explicit upper bound at every finite coupling,
   independent of total volume. Unlike the derivative contraction in
   YM67, this bound has no small-coupling restriction. Its constants
   deteriorate at large coupling.
2. The energy transferred by a block conditional projection is controlled
   by the conditional score in the dual of the block energy norm. The
   coefficient is an explicit conditional Poisson-resolvent matrix.
   A residual formula provides both lower and upper matrix enclosures;
   a finite variational approximation alone provides only the lower one.
3. Directed energy-transfer bounds give a block-sweep Lyapunov criterion
   for a Poincare inequality and hence the full ground-state gap. This
   criterion uses full-energy estimates; it does not insert an
   energy-constrained response norm into YM83's unrestricted row theorem.
4. A smooth gauge-invariant SU(2) comparison family has fixed density
   ratio and fixed conditional L2 response but unbounded projection
   energy amplification. An explicit Gaussian calibration makes the
   matrix certificate sharp, including when the unweighted row test
   is inconclusive.

Conditional covariance identities, Poisson variational bounds,
functional inequalities and the Gaussian/Gauss--Seidel correspondence
are established methods. The development is their explicit, normalized
combination for this program, including the all-coupling actual score
bound, residual requirement, native gauge example and auditable gap
certificate. No general priority claim is made for those methods.

**The actual weak-coupling transfer coefficients have not been shown to
satisfy the matrix gap criterion. A uniform continuum gap and nonzero
continuum-state survival remain open.**

## 1. Model and notation

On the periodic cubic lattice of side L>=4 keep the full SU(2) link
spaces, native one-link Casimir 3/4, and homogeneous theta>=0:

\[
 A=L_E-V,\qquad V=\theta\sum_pW_p,\qquad
 Ah=E_0h,\qquad d\mu=h^2d\Phi,\qquad
 \mathcal E(f)=\mu\sum_e|\nabla_e f|^2.
 \tag{1}
\]

h is the normalized positive smooth vacuum. Multiplication by h
identifies this Dirichlet form with A-E_0. The physical Hamiltonian
is kappa(A-E_0), so a bare gap gamma corresponds to kappa gamma.
The density is not replaced by a plaquette Gibbs ansatz.

Partition the links into nonempty disjoint blocks B_1,...,B_N. Write

\[
 P_i f=\mu(f\mid X_{\bar B_i}),\quad Q_i=I-P_i,\quad
 \mathcal E_i(f)=\mu|\nabla_{B_i}f|^2,\quad
 v_i(f)=\sqrt{\mathcal E_i(f)}.
 \tag{2}
\]

The actual conditional law of B_i at exterior z is nu_i,z. Its positive
reversible generator is
\(L_{i,z}=-\Delta_{B_i}-2\nabla_{B_i}\log h\cdot\nabla_{B_i}\).
Choose a conditional gap floor rho_i>0 uniform in z. Such floors exist
explicitly at every finite theta by YM83; their large-coupling size is
not asserted to be useful for a continuum limit.

All estimates first concern real smooth functions. Closure gives the
form-domain statements. Conditional projections preserve physical gauge
invariance, so estimates valid on the full space restrict to that sector.

## 2. Theorem 87.1 — An actual all-coupling one-link score bound

Fix one link e. More generally, write
\(A=L_e+A_{\rm rest}-V_e\), with A_rest independent of that link,
\(\|V_e\|_\infty\le M\), and a uniform ratio

\[
 h(x,z)/h(y,z)\le R.
 \tag{3}
\]

For M>0 and any t>0,

\[
 |\nabla_e\log h|
 \le R^2 e^{tM}\left((2t)^{-1/2}+M\sqrt{2t}\right).
 \tag{4}
\]

Taking t=1/(4M) gives the simple rational-constant bound

\[
 \boxed{|\nabla_e\log h|^2\le8MR^4.}
 \tag{5}
\]

If M=0 the derivative is zero. The bounds are independent of the number
of exterior links. They do not require a derivative bound on V_e.

**Proof.** Set g(z)=integral h(x,z)dPhi_e(x). Equation (3) gives
R^(-1)g<=h<=Rg. Integrating the vacuum equation in x gives
\((A_{\rm rest}-U)g=E_0g\), with
\(U=\int V_eh\,d\Phi_e/g\) and |U|<=M. Positivity and bounded-potential
comparison imply

\[
 e^{sE_0}e^{-sA_{\rm rest}}g\le e^{sM}g.
 \tag{6}
\]

The free one-link heat K_s=e^(-s L_e) obeys
\(\|\nabla K_su\|_\infty\le\|u\|_\infty/\sqrt{2s}\).
Indeed the native SU(2) metric is the round three-sphere of radius two,
with nonnegative Ricci curvature. The Bochner gradient inequality and
the heat variance identity give

\[
 2s\Gamma(K_su)\le K_s(u^2)-(K_su)^2\le\|u\|_\infty^2.
 \tag{7}
\]

This is a bound for the whole three-component gradient, not one bound
per coordinate followed by an extra dimension factor.

With C=L_e+A_rest, Duhamel's identity is

\[
 h=e^{tE_0}e^{-tC}h+
       \int_0^t e^{sE_0}e^{-sC}(V_eh)\,ds.
 \tag{8}
\]

The rest semigroup is positive and commutes with K_s. Equations (3),
(6), (7) bound the gradient of its first term by
R e^(tM)g/sqrt(2t), and its source integrand by
MR e^(sM)g/sqrt(2s). Integrate, use e^(sM)<=e^(tM), and divide by
h>=g/R to obtain (4). The integrable s^(-1/2) estimate also justifies
differentiation of the Duhamel integral. At t=1/(4M), the square of
the displayed coefficient is (9/2)sqrt(e) M R^4<8MR^4, since
e<3<256/81. If M=0 the positive ground state is constant on that
decoupled link. QED.

For the original homogeneous Wilson model M=4theta. Define

\[
 \mathfrak D_1=\min\left\{
   (21/19)^2 e^{144\theta},\quad
   64e^{8\pi\sqrt{6\theta}}\right\}\qquad(\theta>0),
 \quad\mathfrak D_1(0)=1.
 \tag{9}
\]

YM83 allows R^2<=mathfrak D_1. Thus a valid squared log-gradient bound is

\[
 \boxed{b^2=32\theta\mathfrak D_1^2.}
 \tag{10}
\]

This is a deliberately conservative all-coupling upper bound. In the
small-coupling region YM67's O(theta) bound on the gradient itself is
stronger; one may retain the better estimate. Neither (10) nor YM67
establishes decay with distance under arbitrary full exterior conditioning.

## 3. Theorem 87.2 — Energy transferred by a conditional projection

For j different from i, use the native orthonormal frame on block j and
define the centered conditional score vector

\[
 s_{j|i}(x,z)=2\left[\nabla_{B_j}\log h(x,z)
                 -P_i\nabla_{B_j}\log h(z)\right].
 \tag{11}
\]

Its conditional mean is zero. On that mean-zero space put

\[
 (\Xi_{j\leftarrow i}(z))_{ab}
  =\langle(s_{j|i})_a,L_{i,z}^{-1}(s_{j|i})_b\rangle_{\nu_{i,z}},
 \qquad
 k_{ji}^2\ge\sup_z\|\Xi_{j\leftarrow i}(z)\|_{\rm op}.
 \tag{12}
\]

The inverse uses the actual conditional generator on the mean-zero
space, not the free electric operator or a truncated inverse. These
coefficients are directed: k_ji need not equal k_ij.

Then

\[
 \boxed{v_i(P_i f)=0,\qquad
        v_j(P_i f)\le v_j(f)+k_{ji}v_i(f)\quad(j\ne i).}
 \tag{13}
\]

**Proof.** Differentiate the actual conditional density, including its
normalizing denominator. The exact identity is

\[
 \nabla_{B_j}P_i f
   =P_i\nabla_{B_j}f+\langle f-P_i f,s_{j|i}\rangle_{\nu_{i,z}}.
 \tag{14}
\]

For a unit vector u in the j-gradient fiber, conditional integration
by parts and energy Cauchy--Schwarz give

\[
 |\langle f-P_i f,u\cdot s_{j|i}\rangle|^2
 \le\nu_{i,z}|\nabla_{B_i}f|^2\,
       \langle u\cdot s_{j|i},L_{i,z}^{-1}(u\cdot s_{j|i})\rangle.
 \tag{15}
\]

Take the supremum over u, integrate z, and apply conditional Jensen to
the first term in (14). Minkowski gives (13). Independence from block i
gives its zero gradient. The finite bounds below extend P_i continuously
on the form domain at each finite regulator. QED.

### Actual bounds, and a global energy consequence

The conditional spectral theorem yields

\[
 \Xi_{j\leftarrow i}
 \le\rho_i^{-1}\nu_{i,z}(s_{j|i}s_{j|i}^{T}).
 \tag{16}
\]

If block j has n_j links, (10) bounds the uncentered vector
2 grad_(B_j) log h in squared norm by 4n_j b^2. Centering only reduces
its covariance. Consequently the actual original model has the valid
all-coupling bound

\[
 \boxed{\sup_z\|\Xi_{j\leftarrow i}(z)\|_{\rm op}\le\frac{4n_jb^2}{\rho_i}
              =\frac{128n_j\theta\mathfrak D_1^2}{\rho_i}.}
 \tag{17}
\]

The right side supplies a valid choice of k_ji^2 in (12).
For fixed block sizes these constants are uniform in total
volume. Assigning this same bound to every distant block gives no
summable influence estimate as the number of blocks grows.

Writing c_i^2=sum_(j!=i) k_ji^2, equation (13) also gives

\[
 \mathcal E(P_i f)
 \le\big(\sqrt{\mathcal E(f)-\mathcal E_i(f)}
                +c_i\sqrt{\mathcal E_i(f)}\big)^2
 \le(1+c_i^2)\mathcal E(f).
 \tag{18}
\]

This closes the fixed-regulator energy-domain issue in YM86. It does
not preserve the original spatial support: P_i f can depend on all
exterior links through h.

### A residual certificate for the inverse in (12)

Fix the exterior. For each score component s_a choose a centered trial
function u_a in the domain of L=L_i,z and set r_a=s_a-Lu_a. Define

\[
 G_{ab}=\langle\nabla u_a,\nabla u_b\rangle,\quad
 C_{ab}=\langle u_a,r_b\rangle,\quad R_{ab}=\langle r_a,r_b\rangle,
 \quad F=G+C+C^T.
 \tag{19}
\]

All inner products are actual conditional ones. Then

\[
 \boxed{\Xi=F+[\langle r_a,L^{-1}r_b\rangle]_{ab},
             \qquad F\le\Xi\le F+\rho_i^{-1}R.}
 \tag{20}
\]

To prove it, substitute s_a=Lu_a+r_a and use self-adjointness.
The last matrix is positive and at most R/rho_i by the spectral theorem.
If an exact Galerkin condition makes C=0, the lower matrix is G. Without
that condition the cross terms must be retained. A finite trial space
alone does not certify an upper bound; R must include the entire
residual. No completeness of the trial space is needed when that full
residual is bounded. Uniformity over z is still necessary for (12).

## 4. Theorem 87.3 — A finite matrix certificate for the full gap

Let k_ji>=0 satisfy (12), with k_ii=0. Define the nonnegative matrices
M_i by replacing row i of I by zero and putting (M_i)_ji=k_ji for j!=i.
Thus M_i deletes component i and distributes it to the other components
according to (13). Fix a sweep order 1,...,N, and define

\[
 R_0=I,\quad R_i=M_iR_{i-1},\quad M=R_N,\quad
 \ell_i=e_i^TR_{i-1},\quad
 D=\sum_i\rho_i^{-1}\ell_i^T\ell_i.
 \tag{21}
\]

D has nonnegative entries and is positive definite. The rows ell_i,
in the sweep order, form a triangular matrix with diagonal one.
Suppose a symmetric matrix G_* and beta>0 satisfy

\[
 \boxed{G_*>0,\qquad
       G_*-M^TG_*M\ge D,\qquad G_*\le\beta I.}
 \tag{22}
\]

Then, on the full space and hence on the physical sector,

\[
 \boxed{\operatorname{Var}_\mu f\le\beta\mathcal E(f),
        \qquad A-E_0\ge\beta^{-1}(I-|h\rangle\langle h|).}
 \tag{23}
\]

The physical gap is at least kappa/beta. A volume-uniform beta is needed
for a volume-uniform conclusion. Equation (17) supplies finite inputs
but does not establish (22) for the weak-coupling Wilson vacuum.

**Proof.** Let T=P_N...P_1, f_0=f and f_i=P_i f_(i-1). Equation (13)
gives v(f_i)<=R_i v(f) componentwise. Orthogonal projection and the
conditional Poincare inequality yield the exact loss and its bound:

\[
 \|f\|_2^2-\|Tf\|_2^2
  =\sum_i\|Q_if_{i-1}\|_2^2
  \le\sum_i\rho_i^{-1}\mathcal E_i(f_{i-1})
  \le v(f)^TDv(f).
 \tag{24}
\]

The first two inequalities in (22), with D positive definite, imply
strict contraction of M in the G_* norm, and thus M^r tends to zero.
Iterating (13) shows E(T^r f) tends to zero. At this fixed compact
regulator mu has a finite comparison Poincare constant relative to
product Haar, so T^r f converges in L2 to mu f. That preliminary
constant is used only for convergence, not in the final bound.

Apply (24) to successive T^r f. Since D and M are entrywise nonnegative,
the componentwise energy bounds can be inserted into its quadratic form.
For centered f, telescoping gives

\[
 \begin{aligned}
 \|f\|_2^2
 &\le v(f)^T\sum_{r=0}^{\infty}(M^r)^TD M^r v(f)\\
 &\le v(f)^TG_*v(f)\le\beta\mathcal E(f).
 \end{aligned}
 \tag{25}
\]

The middle bound follows by telescoping (22); convergence follows from
the already proved contraction. G_* itself need not have nonnegative
entries. Form closure and the ground-state transform prove (23). QED.

If the spectral radius of M is below one, the minimal solution of the
matrix equality is the convergent series in (25). Conversely (22)
certifies that stability without a numerical eigenvalue assertion.
These are finite matrix tests only after the actual conditional
coefficient bounds in (12) have been justified.

### Output cutoff consequence

Even without (22), a finite sequence of updates has
E(Tf)<=||M||_2^2 E(f). If its output is supported on S_out, YM86 gives

\[
 \|(I-P^{\rm out}_{\Lambda})Tf\|_\mu^2
 \le D_{S_{\rm out}}\,\|M\|_2^2\mathcal E(f)/\Lambda.
 \tag{26}
\]

The output cutoff must be complete. In the absence of a proved smaller
support, S_out is the entire lattice. Thus (26) supplies a valid omitted
output estimate but no free locality or continuum cost bound.

## 5. Two comparison controls

### 5.1 Fixed L2 response can conceal unbounded energy transfer

Take two edge-disjoint four-link plaquettes A,B, let a=W_A, b=W_B,
and let T_n be the first-kind Chebyshev polynomial. For 0<epsilon<1
and integer n>=3 define the smooth gauge-invariant density

\[
 d\mu_n=[1+\epsilon aT_n(b)]d\Phi.
 \tag{27}
\]

It is normalized and its maximum/minimum ratio is exactly
(1+epsilon)/(1-epsilon), independently of n. The two block marginals
are Haar. This is a comparison density, not the original Wilson vacuum.
It can be the ground-state density of its separately defined smooth
parent Hamiltonian, whose potential depends on n.

The native character identity T_n=(chi_n-chi_(n-2))/2 and
L_(four-link) chi_m=m(m+2)chi_m give

\[
 \Phi(a^2)=1/4,\quad \Phi(T_n^2)=1/2,\quad
 \mathcal E_\Phi(a)=3/4,\quad
 \mathcal E_\Phi(T_n)=n^2/2.
 \tag{28}
\]

The two-block maximal correlation is the constant
epsilon/(2sqrt(2)): the centered conditional response kernel
epsilon aT_n(b) has rank one. Its extremizing functions are physical.
But for f=a, integrating block A gives

\[
 P_Aa=\frac\epsilon4T_n(b),\qquad
 \mathcal E_{\mu_n}(a)=\frac34,\qquad
 \boxed{\frac{\mathcal E_{\mu_n}(P_Aa)}{\mathcal E_{\mu_n}(a)}
            =\frac{\epsilon^2n^2}{24}.}
 \tag{29}
\]

At b=cos(pi/(2n)), T_n(b)=0 and the conditional A law is Haar.
There s_(B|A)=epsilon a grad_B T_n, L_A a=3a, and
|grad_B T_n|^2=n^2. Hence

\[
 \|\Xi_{B\leftarrow A}\|_{\rm op}
       =\epsilon^2n^2/12
       \quad\hbox{at that exterior.}
 \tag{30}
\]

Equation (15) is saturated there by f=a. The energy-transfer coefficient
detects precisely the information invisible to the fixed L2 correlation
and fixed density ratio. The example itself still has a uniform full
Poincare floor (3/4)(1-epsilon)/(1+epsilon) by density comparison.
Thus failure of a proposed energy-transfer criterion is not a proof
that the physical gap closes.

### 5.2 Gaussian star: an exact matrix calibration

This separate noncompact Gaussian check is not a Yang--Mills vacuum.
For a centered Gaussian with positive precision A, scalar conditional
gaps are rho_i=A_ii and the exact coefficients are

\[
 k_{ji}=|A_{ji}|/A_{ii}.
 \tag{31}
\]

Indeed the conditional score is -A_ji times the centered i-variable,
whose variance is 1/A_ii and conditional generator eigenvalue A_ii.
The same sweep argument applies since a finite Gaussian Poincare
constant is available.

Take a center and four leaves, unit diagonal, center-leaf precision
-3/10 and zero other off-diagonal entries. The conditional L2
correlations are 3/10 on those four edges; their unweighted maximum
row sum is 6/5, so YM83's stated unweighted test is inconclusive.
For the order center then leaves, M has only its first row nonzero:

\[
 M_{0\bullet}=(9/25,3/10,3/10,3/10,3/10).
 \tag{32}
\]

Direct calculation gives

\[
 G_*=A^{-1},\qquad G_*-M^TG_*M=D,
 \qquad\lambda_{\max}(G_*)=5/2.
 \tag{33}
\]

Thus (23) gives the exact Gaussian gap 2/5. The linear observable with
gradient (2,1,1,1,1) attains it. A weighted correlation-row argument
also works for this Gaussian, as does YM82's Hessian matrix criterion.
This example verifies normalization and sharpness of the new certificate;
it does not claim to outperform every previous optimized criterion.

## 6. What remains for the continuum problem

The actual all-coupling score bound and conditional energy-domain
control are established here. The new gap theorem is a sufficient
criterion with an explicit unproved input: useful bounds on the
actual matrices Xi under full exterior conditioning. Equation (17)
is generally too large and lacks spatial decay. A finite Poisson
approximation must retain its full residual and exterior-uniform error.

To obtain a continuum conclusion one still needs a sequence of actual
matrix certificates with the required physical scaling, nonzero
renormalized states and their field identification, and the continuum
construction and covariance properties. No spectral or variance
assumption is discharged by the Gaussian or Chebyshev controls.

The next target is therefore an exterior-uniform, spatially summable
bound on the energy-dual score matrices, or another argument that
controls their accumulated effect without taking their worst separate
entries. Equations (12), (20), (22) make that target concrete.

## 7. Reproduction and references

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym87_conditional_energy_transport.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym87_conditional_energy_transport.py -v
```

Exact controls check native polynomial moments, conditional projection
energy, residual matrix enclosures, directed sweep composition,
Lyapunov inequalities and Gaussian sharpness. Symbolic exponential
budgets avoid floating-point underflow/overflow. The source manifest
preserves the previous evidence bytes. Finite computations do not
formally verify the universal proofs, and no independent review is claimed.

- [YM86](YM86_FINITE_ENERGY_RESPONSE.md): complete-cutoff input tails.
- [YM83](YM83_ACTUAL_VACUUM_CONDITIONAL_RESPONSE.md): actual conditional
  gaps and the unrestricted L2 response row criterion.
- [YM67](YM67_POINTWISE_COARSE_VACUUM.md): one-link rest-semigroup
  comparison and sharper small-coupling derivatives.
- [YM82](YM82_NONLINEAR_LATTICE_SPECTRAL_GAPS.md): established
  conditional-gap principles and the mixed-Hessian criterion.
- F. Otto and M. G. Reznikoff, *A new criterion for the logarithmic
  Sobolev inequality and two applications*, J. Funct. Anal. 243 (2007),
  121--157, [doi:10.1016/j.jfa.2006.10.002](https://doi.org/10.1016/j.jfa.2006.10.002).
- G. Menz, *The approach of Otto--Reznikoff revisited*,
  [arXiv:1309.0862](https://arxiv.org/abs/1309.0862), for conditional
  covariance and functional-inequality methods.
- C. Fox and A. Parker, *Accelerated Gibbs sampling of normal
  distributions using matrix splittings and polynomials*,
  [arXiv:1505.03512](https://arxiv.org/abs/1505.03512), for the established
  Gaussian Gibbs/Gauss--Seidel correspondence.
