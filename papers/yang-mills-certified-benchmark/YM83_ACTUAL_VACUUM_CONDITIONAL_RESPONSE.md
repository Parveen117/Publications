# YM83 — Actual vacuum conditional gaps and global response

Monty Dabas. 5 October 2026. Python 3.12 evidence packet.

## Result and scope

NSB1 rules out obtaining a uniform gap from arbitrarily many local jets.
YM82 instead asks for conditional relaxation estimates for the **actual**
vacuum. This stage supplies one of those estimates for the unchanged
Hamiltonian, at every finite coupling and independently of total volume.
It also gives a global assembly test using conditional correlations of
all block functions, rather than pointwise mixed vacuum Hessians.

The three results are:

1. Every finite link block has an explicit, exterior-uniform conditional
   Poincare lower bound in the actual vacuum. Optimizing a heat comparison
   gives an exponential in the square root of interaction amplitude,
   instead of the linear-amplitude exponent from a fixed heat time.
2. These local bounds imply a full quantum gap if a matrix of conditional
   maximal correlations has row sums below one. The matrix is defined in
   the actual vacuum; its entries are not asserted small in weak coupling.
3. A finite matrix of marked nonlinear responses gives a valid upper
   bound on such a correlation only when accompanied by a certified
   remainder for the omitted functions. There is an exact Hilbert-space
   remainder formula and a finite example showing why the remainder is
   necessary.

The new unconditional estimate is local, not a uniform four-dimensional
mass gap. The global correlation bound, its physical scaling, and a
nontrivial covariant continuum construction remain open. Heat comparison,
bounded-density Poincare comparison and the theory of two projections
are established methods. No priority claim is made for those methods or
for a general conditional functional-inequality principle.

## 1. Unchanged model

On a finite periodic cubic lattice of side at least four, use

\[
 A=L_E-V,\qquad V=\sum_p\Phi_p(W_p),\qquad
 L_E=-\sum_e\Delta_e,\qquad W_p=\tfrac12\operatorname{Tr}U_{\partial p}.
 \tag{1}
\]

The functions \(\Phi_p\) are smooth and real. The original homogeneous
model is \(\Phi_p(w)=\theta w\), \(\theta\ge0\). Native SU(2) has
fundamental Casimir \(3/4\), dimension three, diameter \(2\pi\), and
Ricci tensor \(g/2\). Haar measure is normalized to one.

Let \(h>0\) be the normalized ground state, \(Ah=E_0h\), and put

\[
 H=A-E_0,\qquad d\mu=h^2d\Phi,\qquad
 \mathcal L=h^{-1}Hh,\qquad
 \mathcal E(f)=\mu\sum_e|\nabla_e f|^2.
 \tag{2}
\]

At fixed volume the vacuum is smooth, strictly positive and unique by
ellipticity and positivity improvement. Equation (2) is its exact unitary
ground-state transform. No formula such as \(h=e^{V/2}/Z^{1/2}\) is used.

For a nonempty set B of n links, write the remaining configuration as z.
Let \(V_B\) contain every face potential whose support meets B, counted
once, and define the computable bound

\[
 \Omega_B=\sum_{p:\,p\cap B\ne\varnothing}
                  \operatorname{osc}_{[-1,1]}\Phi_p.
 \tag{3}
\]

Adding constants to the face potentials does not change this quantity.
Let \(m_t=\inf k_t\), \(M_t=\sup k_t\) for the one-link free heat
kernel. They are strictly positive and finite for t>0. Define

\[
 R_B(t)=e^{t\Omega_B}(M_t/m_t)^n.
 \tag{4}
\]

## 2. Theorem 83.1 — An actual conditional gap at every coupling

For every exterior configuration z, every x,y in the block, and t>0,

\[
 \frac{h(x,z)}{h(y,z)}\le R_B(t),\qquad
 \operatorname{Var}_{\mu(\cdot\mid z)} f
 \le\frac{4R_B(t)^2}{3}\,
       \mu\left(|\nabla_B f|^2\mid z\right).
 \tag{5}
\]

Consequently a valid conditional gap floor is
\(\rho_B\ge(3/4)R_B(t)^{-2}\), uniformly in z and total volume.
Two explicit versions are

\[
 \rho_B\ge\frac34(19/21)^{2n}e^{-18\Omega_B},
 \tag{6}
\]

and, for \(\Omega_B>0\),

\[
 \boxed{\rho_B\ge\frac34\,64^{-n}
        \exp[-4\pi\sqrt{3n\Omega_B}].}
 \tag{7}
\]

One may use the larger floor. If \(\Omega_B=0\), the block is free and
its exact conditional gap is 3/4. None of these floors is claimed sharp.

**Proof of (5).** Center each face potential in \(V_B\) by the midpoint
of its range. After an irrelevant scalar shift of A,

\[
 A=L_{E,B}+A_{\bar B}-\widetilde V_B,
 \qquad\|\widetilde V_B\|_\infty\le\Omega_B/2,
 \tag{8}
\]

where \(A_{\bar B}\) is independent of x. Positivity and Feynman--Kac
give a pointwise sandwich, on positive functions,

\[
 e^{-t\Omega_B/2}e^{-t(L_{E,B}+A_{\bar B})}h
 \le e^{-tA}h=e^{-tE_0}h
 \le e^{t\Omega_B/2}e^{-t(L_{E,B}+A_{\bar B})}h.
 \tag{9}
\]

The decoupled semigroup is a product in B and its exterior. If
\(g(z)=\int h(x,z)d\Phi_B(x)\), its value in (9) is between
\(m_t^n e^{-tA_{\bar B}}g(z)\) and
\(M_t^n e^{-tA_{\bar B}}g(z)\). The same strictly positive exterior
factor occurs at x and y and cancels. This proves the ratio bound.
There is no extensive bound on \(\|V\|_\infty\) in this argument.

The conditional density is proportional to \(h(x,z)^2\), so its
maximum/minimum ratio is at most \(R_B(t)^2\). Product Haar has gap
3/4 for any nonempty number of links. Center a test function by its Haar
mean, use the upper conditional density in its variance and the lower
density in its energy. The ratio of these two bounds proves (5).
Closure extends the inequality from smooth functions to the conditional
Sobolev form domain. Positivity on the compact finite product makes the
statement valid at every z, not just almost everywhere. QED.

**Explicit heat times.** The established native bound
\(19/20\le k_9\le21/20\) from YM53/YM67 proves (6).
For (7), use the classical Li--Yau heat Harnack inequality, with the
generator \(\Delta\), on a compact d-dimensional manifold of nonnegative
Ricci curvature:

\[
 u(s,x)\le u(t,y)(t/s)^{d/2}
            \exp[d(x,y)^2/(4(t-s))],\quad 0<s<t.
 \tag{10}
\]

Apply it to \(u(t,x)=k_t(x,z)\). Integrating its early spatial
variable at s=t/2, and separately its late spatial variable at 2t,
uses \(\int k_s(x,z)d\Phi(x)=1\) and gives, for one native link,

\[
 m_t\ge2^{-3/2}e^{-2\pi^2/t},\qquad
 M_t\le2^{3/2}e^{\pi^2/t},\qquad
 M_t/m_t\le8e^{3\pi^2/t}.
 \tag{11}
\]

Thus (5) gives
\(\rho_B\ge(3/4)64^{-n}e^{-2t\Omega_B-6\pi^2n/t}\).
The exponent is minimized at \(t=\pi\sqrt{3n/\Omega_B}\), with value
\(4\pi\sqrt{3n\Omega_B}\). This proves (7). Only the **free** heat
uses nonnegative curvature; no curvature sign for the actual vacuum is
assumed.

### Explicit original-model cases

Let \(N_B\) be the number of faces meeting B. Then
\(\Omega_B=2\theta N_B\). A single link has n=1, \(N_B=4\), hence

\[
 \boxed{\rho_{\rm link}\ge\max\left\{
 \frac34(19/21)^2e^{-144\theta},\quad
 \frac3{256}e^{-8\pi\sqrt{6\theta}}\right\}.}
 \tag{12}
\]

For example, at \(\theta=1/16000\), \(e^{-9/1000}\ge991/1000\)
gives the actual conditional bound
\(\rho_{\rm link}\ge357751/588000>3/5\), at every exterior and
every finite lattice volume. This number is a local conditional gap,
not a global gap or a particle mass.

The four links of one elementary plaquette have n=4, \(N_B=13\):
the selected face is counted once and twelve other faces meet its edges.
The optimized exponent is \(4\pi\sqrt{312\theta}\), with prefactor
\((3/4)64^{-4}\). The evidence checks these incidences on sides 4, 6
and 8. They are properties of the full overlapping lattice, not of an
isolated plaquette replacement.

YM67 already contains the one-link fixed-time ratio mechanism in its
small-coupling setting. The present continuation extracts conditional
Poincare estimates, extends the comparison to arbitrary finite blocks
and amplitudes, and adds the optimized square-root-amplitude estimate.
It does not claim that the underlying Feynman--Kac comparison is new.

## 3. Theorem 83.2 — Global assembly by conditional response

Partition all links into nonempty disjoint blocks \(B_1,\ldots,B_N\).
In \(L^2(\mu)\) define

\[
 P_i f=\mu(f\mid X_{\bar B_i}),\qquad Q_i=I-P_i,
 \qquad K=\sum_iQ_i.
 \tag{13}
\]

K is a bounded positive conditional-resampling generator, with each
block updated at rate one. It is an auxiliary comparison operator,
not the quantum Hamiltonian or a change of physical clock.

For each i<j and fixed exterior \(z=X_{\overline{B_i\cup B_j}}\),
take the actual two-block conditional law \(\nu_{ij,z}\). Its maximal
correlation is

\[
 c_{ij}(z)=\sup\left\{|\nu_{ij,z}(uv)|:
 u=u(X_{B_i}),\ v=v(X_{B_j}),\quad
 \nu u=\nu v=0,\ \nu u^2=\nu v^2=1\right\}.
 \tag{14}
\]

Choose symmetric bounds \(c_{ij}\ge\operatorname*{ess\,sup}_z
c_{ij}(z)\), \(c_{ii}=0\), and set

\[
 r=\max_i\sum_{j\ne i}c_{ij},\qquad
 \rho_* =\min_i\rho_{B_i}>0.
 \tag{15}
\]

If r<1, then

\[
 \boxed{\operatorname{gap}K\ge1-r,\qquad
        H\ge\rho_*(1-r)(I-|h\rangle\langle h|).}
 \tag{16}
\]

The assertion is on the full link space and therefore restricts to the
physical gauge-invariant space. For a two-block partition, the gap of K
is exactly \(1-c_{12}\), where there is no exterior and c is the true
maximal correlation, not merely an upper bound.

**Proof.** On a two-block conditional fiber, P_i and P_j project onto
functions of the other block. Their common range is the constants:
the joint density is everywhere positive relative to product Haar.
On the complement of constants their product norm is c from (14).
The elementary two-projection spectral identity gives spectrum
\(1\pm s\), with \(0\le s\le c\), for the sum of the projections,
plus possible 0 and 1 parts. Equivalently,

\[
 (Q_i+Q_j)^2\ge(1-c)(Q_i+Q_j),
 \quad Q_iQ_j+Q_jQ_i\ge-c(Q_i+Q_j).
 \tag{17}
\]

One way to obtain the identity is to write the generic projection
blocks as \(\operatorname{diag}(I,0)\) and
\(\begin{pmatrix}C^2&C\sqrt{I-C^2}\\
C\sqrt{I-C^2}&I-C^2\end{pmatrix}\), and apply scalar spectral
calculus to the positive contraction C. The norm of C is c. The
orthogonal and common subspaces satisfy (17) directly. This also proves
the stated exact two-block gap. This is standard two-projection theory.

Disintegrate (17) over z. Since \(Q_i^2=Q_i\), summing over pairs gives

\[
 K^2\ge\sum_i\left(1-\sum_{j\ne i}c_{ij}\right)Q_i
        \ge(1-r)K.
 \tag{18}
\]

The spectral theorem excludes spectrum of K in (0,1-r). Its kernel is
the constants: being in every range P_i makes a function independent
of each block; equivalence to product Haar at this fixed volume then
forces it to be constant. Thus
\(\sum_i\|Q_i f\|_2^2\ge(1-r)\operatorname{Var}_\mu f\).
By the conditional inequalities in (5),

\[
 \mathcal E(f)\ge\sum_i\rho_{B_i}\|Q_i f\|_2^2
 \ge\rho_*(1-r)\operatorname{Var}_\mu f.
 \tag{19}
\]

The closed form and (2) prove (16). Gauge transformations act separately
on every link and preserve \(\mu\). Conditional integration over any
fixed subset of links is equivariant, so P_i and K preserve the Gauss
sector. Restricting the full-space inequality is consequently valid.
QED.

### What the theorem does and does not remove

It replaces the pointwise mixed-Hessian assumption of YM82 by bounded
conditional expectation operators. It does not deduce their smallness
from finite jets, unconditioned correlations, or the averaged Ward row
identity. The actual vacuum may have nonlocal conditional dependence:
\(c_{ij}=0\) must not be assigned simply because two blocks share no
bare plaquette. If r>=1, this sufficient test is inconclusive; it does
not establish a closing quantum gap.

## 4. Theorem 83.3 — Marked responses with an explicit omitted tail

Fix a pair and exterior. Let \(\nu\) be its joint law, with marginals
\(\nu_1,\nu_2\), and define

\[
 k(x,y)=\frac{d\nu}{d(\nu_1\otimes\nu_2)}(x,y)-1,
 \qquad\Xi=\int|k|^2d\nu_1d\nu_2.
 \tag{20}
\]

At any fixed finite regulator and exterior the smooth positive compact
density makes \(\Xi\) finite. Its magnitude is not bounded here
uniformly in the regulator. Choose finite real orthonormal centered
families \(a_1,\ldots,a_m\) in \(L^2(\nu_1)\) and
\(b_1,\ldots,b_n\) in \(L^2(\nu_2)\), and measure the marked matrix

\[
 M_{\alpha\beta}=\nu(a_\alpha b_\beta).
 \tag{21}
\]

Then the exact residual and the correlation bounds are

\[
 \|T-E_mTF_n\|_{\rm HS}^2=\Xi-\|M\|_F^2,
 \qquad
 \boxed{\|M\|_{\rm op}\le c\le
 \min\{1,\|M\|_{\rm op}+\sqrt{\Xi-\|M\|_F^2}\}.}
 \tag{22}
\]

Here \(Tv(x)=\int k(x,y)v(y)d\nu_2(y)\) is the centered conditional
expectation operator, c=\(\|T\|\), and E_m,F_n are the indicated
orthogonal projections. To use an inexact matrix \(\widehat M\), an
operator-norm error \(\delta\) and a separately certified omitted
Hilbert--Schmidt norm \(\varepsilon\) suffice:

\[
 c\le\min\{1,\|\widehat M\|_{\rm op}+\delta+\varepsilon\}.
 \tag{23}
\]

A certified upper value for \(\Xi\) and certified lower value for
\(\|M\|_F^2\) give an upper value for the squared residual. Incorrectly
subtracting the norm of an inexact matrix can invalidate that bound.
Every such bound must also be uniform over the required exterior z
before insertion into (15).

**Proof.** The kernel integrates to zero in either marginal variable.
The functions \(a_\alpha(x)b_\beta(y)\) are orthonormal in the product
space, and their kernel coefficients equal (21). Orthogonal projection
and Pythagoras give the residual identity. The operator norm of an
integral kernel is bounded by its Hilbert--Schmidt norm; the triangle
inequality proves the upper bound. Compression proves the lower bound.
Another triangle inequality proves (23). QED.

This is the precise nonlinear-response repair used here. A finite
observed matrix alone supplies a lower bound on maximal correlation,
whereas the gap test needs an upper bound. The omitted-tail estimate
is therefore part of the theorem's input, not an optional numerical
accuracy setting.

### Exact finite control: invisible high-order dependence

For N>=3 binary variables define the positive probability law

\[
 \mu_t(x)=2^{-N}(1+t x_1\cdots x_N),\qquad |t|<1.
 \tag{24}
\]

Every proper marginal is the same independent uniform law, for every t.
Nevertheless, conditioned on the other N-2 variables, any pair has
maximal correlation \(|t|\). Its single-site conditional-resampling
generator has exact gap

\[
 \boxed{\operatorname{gap}K_t=
 \frac{N-\sqrt{(N-2)^2+4t^2(N-1)}}2\longrightarrow0
 \quad\text{as }|t|\uparrow1.}
 \tag{25}
\]

To prove this, write \(\chi_S=\prod_{i\in S}x_i\). For nonempty
proper S of size k,
\(K_t\chi_S=k\chi_S-kt\chi_{S^c}\). The two-dimensional span of
\(\chi_S,\chi_{S^c}\) has eigenvalues
\([N\pm\sqrt{(N-2k)^2+4t^2k(N-k)}]/2\).
The smaller is minimized at k=1 or N-1. Constants and the full parity
give eigenvalues 0 and N. This proves (25). The row test supplies the
valid, possibly weaker floor \(1-(N-1)|t|\) when positive.

For two *blocks* of binary variables, the density
\(1+t\chi_A(x)\chi_B(y)\) has a rank-one centered correlation
operator of norm \(|t|\). All markers orthogonal to either parity see
zero, while the squared omitted norm in (22) is exactly \(t^2\).
Adding the two parity markers recovers it exactly. These finite controls
are comparison probability models, not Yang--Mills vacuum computations.

## 5. Physical scaling and the next actual estimate

The physical Hamiltonian is \(\mathsf H_a=\kappa_a H_a\), with
\(\kappa_a\sqrt{\theta_a}=c/a\) as in YM76--YM82. The sufficient
physical conclusion from (16) would be

\[
 \operatorname{gap}\mathsf H_a\ge
       \kappa_a\rho_{*,a}(1-r_a).
 \tag{26}
\]

Theorem 83.1 proves \(\rho_{*,a}>0\) for bounded-size blocks at every
finite coupling, uniformly in volume and exterior configuration. It
does not prove that this floor has the required continuum scaling.
For example, on the **illustrative**, not identified renormalized,
trajectory \(\sqrt\theta=b\log(1/a)\), the single-link floor (7)
times the physical clock equals

\[
 \frac{3c}{256b\log(1/a)}
       a^{8\pi\sqrt6\,b-1}.
 \tag{27}
\]

It tends to zero whenever \(8\pi\sqrt6\,b>1\), even if a positive
row margin were available. A vanishing **lower bound** is not evidence
that the actual gap vanishes. Larger physical blocks or better local
estimates may be needed; (7) itself worsens with block size.

The concrete next target is to bound the conditional response kernels
of the homogeneous vacuum, including their omitted tails, and find
blocks for which the physical product in (26) stays positive. The
existing local Ward identities do not certify these tails. Theorem
83.1 removes the absence of any explicit all-coupling local conditional
floor; it does not remove the global or scale-uniform obligations.

Nonzero local states must also survive the continuum limit (YM77--YM81),
and their field identification and Euclidean covariance must be proved.
A positive conditional or fixed-lattice gap alone establishes none of
those construction requirements.

## 6. Evidence and reproduction

The written arguments above prove the universal statements subject to
their explicit hypotheses. The Python packet checks constants using
exact rational arithmetic, native face incidences, projection algebra,
response-tail controls and full finite probability matrices. It does
not solve for the actual weak-coupling vacuum or formally verify the
analytic proofs. Tests are a separate algebraic implementation where
possible; no independent human or agent review is claimed.

The source manifest pins the unchanged NSB1 lineage and the new local
inputs. The result's canonical hash also binds the manifest. Default
and --check runs are read-only; --write only regenerates the new result
and its hash. Previous stage evidence is left unchanged.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym83_actual_vacuum_response.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym83_actual_vacuum_response.py -v
```

## References and lineage

- [NSB1](../nonlinear-spectral-blindness/NSB1_NONLINEAR_SPECTRAL_BLINDNESS.md),
  especially the finite-jet gauge counterfamily and Yang--Mills dock.
- [YM67](YM67_POINTWISE_COARSE_VACUUM.md), Theorem 1, for the existing
  one-link fixed-time vacuum-ratio mechanism.
- [YM82](YM82_NONLINEAR_LATTICE_SPECTRAL_GAPS.md), Section 4, for the
  conditional mixed-Hessian criterion and its actual-vacuum obligation.
- P. Li and S.-T. Yau, *On the parabolic kernel of the Schroedinger
  operator*, Acta Math. 156 (1986), 153--201, Section 2 (q=0, nonnegative
  Ricci case). [Primary publication](https://doi.org/10.1007/BF02399203).
- Q. Qin, *Analysis of two-component Gibbs samplers using the theory of
  two projections*, Sections 2.2--2.3 and 4,
  [arXiv:2201.12500v4](https://arxiv.org/html/2201.12500v4).
- F. Otto and M. G. Reznikoff, *A new criterion for the logarithmic
  Sobolev inequality and two applications*, JFA 243 (2007), 121--157;
  Remark 4 describes the earlier spectral-gap principle.
  [Primary preprint](https://webdoc.sub.gwdg.de/ebook/serien/e/sfb611/274.pdf).
- A. Jaffe and E. Witten, *Quantum Yang--Mills Theory*, for the full
  existence and mass-gap target.
  [Official problem statement](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf).
