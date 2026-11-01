# YM91 — Conditional memory tensors and exact flux certificates

YM90 bounds the actual dynamic memory by a worst-case cross-Hessian norm
and a conditional gap. This continuation replaces that sufficient test
by a conditional Poisson response tensor, then gives a flux certificate
for its upper bound. The replacement applies to the actual finite Wilson
vacuum at every finite coupling; a useful uniform numerical bound at weak
coupling still has to be established for that vacuum.

A separate, explicit native SU(2) parent family shows that the replacement
can succeed when the YM90 Hessian test cannot. Its pointwise weighted
Ricci curvature tends to minus infinity, while its certified relative
memory tends to zero. These are engineered comparison Hamiltonians,
not the original cubic Wilson Hamiltonian. The original Wilson numerical
window remains 0 <= theta <= 1/400, inherited from YM90. No continuum
mass gap or continuum-state construction is claimed.

## 1. The actual conditional tensor bounds dynamic memory

Keep the actual model, native metric, block paths, and notation of
[YM90](YM90_GAUGE_MEMORY_BOUND.md). In particular

\[
 d\mu=h^2d\Phi,\quad\chi=\log h,\quad
 \mathcal L=-\Delta-2\nabla\chi\cdot\nabla,
 \qquad\mathcal L=\begin{pmatrix}C_b&B^*\\B&D\end{pmatrix}.
 \tag{1}
\]

The paths have length b >= 2 on a fine periodic cubic lattice of side
Lb, L >= 4. Let q denote all coarse path products, R their retained
fine links, and Z the other fine links. Internal gauge fixing gives

\[
 d\mu(\cdot\mid q)=d\Phi_G(g)\,d\nu_q(z),\qquad
 d\nu_q(z)\propto h(\operatorname{section}(q),z)^2d\Phi_Z(z).
 \tag{2}
\]

No arbitrary fine pins are imposed in constructing (2). For each q let

\[
 A_q=-\nu_q^{-1}\operatorname{div}_Z(\nu_q\nabla_Z)
 \tag{3}
\]

on its conditional mean-zero space. At each fixed finite volume and
finite coupling, smooth positivity and compactness give a positive
conditional floor uniform in q. This existence assertion does not make
that floor uniform in volume or useful in weak coupling.

Let S_E be YM90's transported horizontal log-wave gradient. Concatenate
the three native components for all coarse links, and define

\[
 r(q,z)=2(S-\nu_q S),\qquad Bf=-r\cdot\nabla f(q).
 \tag{4}
\]

The source is internal-gauge invariant. On that invariant hidden sector
write D_0 for D. Its form dominates the direct integral of (3):
\(D_0\ge\mathfrak A:=\int^\oplus A_q\,w_b(q)d\Phi(q)\).
Indeed, retaining just the unretained native link energies gives the
right side; the other fine gradient energies are nonnegative. The
coordinate change is an isometry in z with retained coordinates fixed,
although it is not a product decomposition of the entire fine metric.

**Conditional-memory theorem.** Define the positive tensor

\[
 T_s(q)_{ab}=\langle r_a,(A_q+s)^{-1}r_b\rangle_{\nu_q},
 \qquad s\ge0.                                         \tag{5}
\]

For every smooth coarse f,

\[
 \boxed{0\le\Sigma_b(s)[f]
 \le\int \nabla f(q)^T T_s(q)\nabla f(q)\,w_b(q)d\Phi(q).}
 \tag{6}
\]

**Proof.** Apply the variational formula for an inverse positive form to
D_0+s and to the dominated direct-integral form A+s. Its latter domain
is allowed to be larger: it requires z energy but no q regularity.
The direct-integral inverse acts separately at each q. Inserting (4)
gives (6). Form closure extends the assertion to the relevant energy
domain. This does not replace D_0 by A as an operator identity. QED.

In particular, the sufficient bound

\[
 T_0(q)\le\eta b I\quad\text{for almost every }q,
 \qquad\eta<1,                                        \tag{7}
\]

gives \(\Sigma_b(0)\le\eta C_b\), since
\(C_b[f]=b\int|\nabla f|^2w_b d\Phi\). One may instead prove the
integrated quadratic-form inequality itself for every admissible f;
an unweighted exterior average of T_0 is not that inequality. Its
quadratic form must be controlled on all gradients being certified.
All mixed source components are kept in T_s; no diagonal or commuting
response assumption is made.

This tensor refines YM90's bound. If A_q >= d and the cross-Hessian
matrix has norm a, the conditional variance proof of YM90 gives

\[
 T_s\le(d+s)^{-1}\nu_q(rr^T)
       \le\frac{4ba^2}{d(d+s)}I.                         \tag{8}
\]

Thus (6) recovers its actual 1/144 endpoint. The new option is to bound
T_s without first replacing all its source channels by that worst-case
Hessian estimate.

## 2. Exact fluxes certify an upper bound without a gap constant

All operations in this section are conditional at fixed q. A vector
field J on the unretained link manifold has weighted divergence
\(\operatorname{div}_{\nu}J=\nu^{-1}\operatorname{div}(\nu J)\).
For each component r_a, find a flux J_a satisfying

\[
 -\operatorname{div}_{\nu}J_a=r_a.                        \tag{9}
\]

Then the Gram matrix

\[
 U_{ab}=\int J_a\cdot J_b\,d\nu
 \quad\text{satisfies}\quad \boxed{T_0\le U}.            \tag{10}
\]

To prove it, for arbitrary coefficients u put r_u=sum u_a r_a and
J_u=sum u_a J_a. Integration by parts and Cauchy--Schwarz give
\(\langle r_u,v\rangle^2\le\|J_u\|^2\mathcal E_Z(v)\).
Take the supremum over v. Equality is achieved by the gradient flux
J_a=grad A^(-1)r_a. Thus (10) is the minimum-flux principle, applied
simultaneously to all source directions. A numerical conditional gap
is not needed to verify an exact feasible flux and its norm.

There is a monotone improvement mechanism. Given divergence-free fields
K_1,...,K_m, set G_ij=nu(K_i dot K_j) and L_ia=nu(K_i dot J_a).
For independent K_i, subtract their orthogonal projection:

\[
 J_a^{\rm new}=J_a-\sum_i K_i(G^{-1}L)_{ia},\qquad
 U^{\rm new}=U-L^TG^{-1}L\le U.                          \tag{11}
\]

The divergence equation is unchanged. Enlarging nested spaces of such
fields monotonically lowers the certified Gram. The full orthogonal
projection onto the divergence-free subspace gives the optimal gradient
flux. Harmonic fields, when present, belong to that subspace and are not
silently omitted from the limiting statement.

An approximate flux with a nonzero divergence residual cannot simply
use (10). An alternative all-tail certificate uses centered trial
functions u_a in the domain of A+s. Set

\[
 e_a=r_a-(A+s)u_a,\qquad
 F_{ab}=\langle r_a,u_b\rangle+\langle u_a,r_b\rangle
             -\langle u_a,(A+s)u_b\rangle.
\]

The completion identity already underlying YM87 becomes

\[
 T_s=F+[\langle e_a,(A+s)^{-1}e_b\rangle]_{ab},\qquad
 F\le T_s\le F+(d+s)^{-1}[\langle e_a,e_b\rangle]_{ab}
 \tag{12}
\]

when A >= d. Every omitted mode must be included in the residual norm.
A finite variational solve supplies only the lower conditional tensor;
it supplies no upper bound until that residual is controlled. Moreover,
because T_s itself is an upper comparison for dynamic memory, its lower
bound F is not a lower bound on Sigma_b(s).

## 3. What a certified flux does, and does not, imply

If the flux Gram obeys U <= eta b I with eta < 1, inverse monotonicity
gives the exact projected-resolvent comparison

\[
 0\le I_b^*(\mathcal L+s)^{-1}I_b-(C_b+s)^{-1}
 \le\frac{\eta}{1-\eta}(C_b+s)^{-1}\quad(s>0).            \tag{13}
\]

At s=0 use the centered coarse space, with invertibility supplied by
positive coarse and hidden floors. The source certificate alone does
not supply a full hidden floor. If those independent bounds are c and d,
the physical two-block argument gives the sufficient floor

\[
 \Delta_{\rm physical}\ge\frac{cd(1-\eta)}{c+d}.          \tag{14}
\]

In contrast, an unobserved hidden mode can be arbitrarily slow with
fixed small memory: the matrix with diagonal (1,1,epsilon) and only
H_12=H_21=1/2 has Sigma(0)=1/4 and full gap epsilon for 0<epsilon<1/2.
Its projected zero-frequency resolvent is always 4/3.

The estimate Sigma(0) <= eta C also does not imply
Sigma(s) <= eta d/(d+s) C; YM90's scalar counterexample is replayed.
Such frequency decay requires additional source-energy information.
If D_0 >= d is separately known, one still obtains

\[
 M_n\le\eta d^{-n}C_b,\qquad
 \int_T^\infty K_b(t)dt\le\eta e^{-dT}C_b.               \tag{15}
\]

For t>0 the corresponding pointwise kernel bound is
\(K_b(t)\le\eta\sup_{\lambda\ge d}\lambda e^{-t\lambda}C_b\).
The supremum equals 1/(e t) when dt<=1, and d exp(-dt) when dt>=1.
There is no uniform time-zero source bound from a flux certificate alone.

## 4. A native family with decreasing memory and diverging curvature

Take two disjoint b-link paths Q and R from vertex o to vertex v and a
single further link z from o to v. For b=2 this is exactly the five-link
graph used as a comparison in YM90. All internal path vertices have
degree two. Let U=Qz^(-1), and for an odd integer n>=1 define

\[
 F_n(U)=\frac{\chi_n(U)}{n+1}
       =\frac{U_n(\operatorname{Re}U)}{n+1},\quad
 \rho_n(U)=1+\epsilon F_n(U),\quad h_n=\sqrt{\rho_n(U)},
 \qquad0\le\epsilon<1.                                  \tag{16}
\]

Here chi_n is the spin-n/2 character, U_n the Chebyshev polynomial of the
second kind, and its native Casimir is lambda_n=n(n+2)/4. Unitarity and
Haar character orthogonality give

\[
 |F_n|\le1,\quad\Phi F_n=0,\quad
 \Phi F_n^2=(n+1)^{-2},\quad F_n(-U)=-F_n(U).             \tag{17}
\]

h_n is already normalized in full product Haar. It is gauge invariant
and is the exact positive zero-energy ground state of -Delta-V_n, where

\[
 V_n=(b+1)\left[
   \frac{\epsilon\lambda_nF_n}{2\rho_n}
   +\frac{\epsilon^2\Gamma_U(F_n)}{4\rho_n^2}\right].     \tag{18}
\]

Indeed V_n=-Delta h_n/h_n; ground transformation makes the resulting
operator nonnegative. Equation (18) defines an engineered parent, with
n-dependent nonlinear interactions. It is not a pure Wilson action.

Integrating z leaves exactly product Haar on the coarse pair (Q,R).
At fixed Q, the unretained law is rho_n(Qz^(-1))dPhi(z). Its native gap
has the elementary density-comparison floor

\[
 d_n\ge d_*:=\tfrac34\frac{1-\epsilon}{1+\epsilon}.       \tag{19}
\]

On physical coarse functions, the only relative holonomy is QR^(-1)
up to conjugation. Characters give the exact centered coarse floor
c=3b/2. These facts hold for every n without any positive-curvature
assumption.

### An explicit flux removes the representation-frequency penalty

For a coarse native direction u on Q, choose the right-invariant native
field V_u on z whose variation of U equals the left-native variation
of Q. It has norm |u| and zero Haar divergence. Thus the source from
(4) is r_u=b V_u log rho_n. For any constant c_0, the field

\[
 J_u=-b(1-c_0/\rho_n)V_u
 \quad\text{satisfies}\quad
 -\operatorname{div}_{\nu}J_u=bV_u\log\rho_n=r_u.         \tag{20}
\]

This follows directly by differentiating rho_n J_u; it is valid also
when derivatives of rho_n are large. The source for the other coarse
path R is zero. Write H_n=Phi(1/rho_n). Its flux cost is

\[
 \nu|J_u|^2=b^2(1-2c_0+c_0^2H_n)|u|^2.
\]

Minimizing the scalar c_0 gives c_0=1/H_n and the valid relative bound

\[
 \eta_n\le b(1-H_n^{-1}).                               \tag{21}
\]

This minimizes this explicit family of feasible fluxes, not all fluxes;
the exact conditional tensor, and especially the actual dynamic
memory, may be smaller. Odd parity makes the entire reciprocal series
positive:

\[
 H_n=\sum_{k=0}^{\infty}\epsilon^{2k}\Phi F_n^{2k},\qquad
 H_n-1\le\frac{\epsilon^2}{(1-\epsilon^2)(n+1)^2}=:q_n.
\]

Since x/(1+x) increases for x>=0, (21) yields

\[
 \boxed{\eta_n\le\frac{b\epsilon^2}
 {(1-\epsilon^2)(n+1)^2+\epsilon^2}.}                    \tag{22}
\]

For b=2 and epsilon=1/4 this is **2/[15(n+1)^2+1]**, tending to zero
as n grows. It holds on the entire coarse energy domain. Restricting
to the physical coarse subspace can only improve a valid upper bound.
With d_*=9/20 and c=3, (14) supplies a positive physical comparison gap
for every n; a direct density comparison can give better numerical gap
constants on this fixed graph. The result of interest is the memory
estimate and its separation from a pointwise Hessian criterion.

### The YM90 Hessian sufficient test genuinely fails here

At U=-1 the gradient vanishes and

\[
 \operatorname{Hess}_z\log h_n
       =\frac{\epsilon\lambda_n}{6(1-\epsilon)}I=:k_n I,
 \qquad\operatorname{Ric}_z-2\operatorname{Hess}_z\log h_n
       =(\tfrac12-2k_n)I.                               \tag{23}
\]

Each of the b cross blocks from z to the first retained path has norm
k_n at that configuration. Thus the worst-case cross-Hessian matrix
has a squared >= b k_n squared. The conditional native gap is at most
3/4: each imaginary coordinate of U has mean zero, variance 1/4, and
energy 3/16 under rho_n, by odd parity and centrality. Therefore no
valid floor d can exceed 3/4, and the best possible ratio in the YM90
sufficient test has lower bound

\[
 \frac{4a^2}{d^2}\ge\frac{4bk_n^2}{(3/4)^2}.
 \tag{24}
\]

At epsilon=1/4 this equals b[n(n+2)]^2/729. With b=2 it already
exceeds one for n=5, and grows as n to the fourth power.

| Odd order n | Certified relative memory upper | Weighted Ricci at U=-1 |
| --- | --- | --- |
| 1 | 2/61 | 5/12 |
| 3 | 2/241 | 1/12 |
| 5 | 2/541 | -17/36 |
| 9 | 2/1501 | -9/4 |
| 31 | 2/15361 | -335/12 |

This demonstrates failure of that particular worst-case sufficient
test, not a physical instability. The exact feasible flux and its full
norm succeed despite it. The comparison does not establish that the
original Wilson vacuum has this special translated conditional density.

## 5. Full tails, controls, and the remaining Wilson calculation

All Haar moments in the reciprocal series are exact rational character
integrals. If H_(n,N) stops at k=N, then

\[
 0\le H_n-H_{n,N}\le
 \frac{\epsilon^{2N+2}\Phi F_n^{2N+2}}{1-\epsilon^2}.      \tag{25}
\]

This follows from |F_n|<=1, and encloses every omitted term. With an
upper endpoint H_U, choosing the rational constant c_0=1/H_U in (20)
already certifies the flux cost b squared times (1-1/H_U). A finite
integral computation is thus accompanied by an explicit full-tail
error; no finite-spin approximation to the density or dynamics is made.

Controls verify the native Casimir identity, conditional Rayleigh upper
bound, exact cross-Hessian normalization at U=-1, noncommuting word
derivatives, and gauge invariance of the physical source. Finite operator
controls check the conditional-to-dynamic inequality, exact divergence
constraints, monotone flux corrections, residual Grams and all matrix
cross terms. They also reject three erroneous conclusions:

- a lower conditional tensor is a lower dynamic-memory bound;
- a small exterior average bounds every coarse energy direction;
- small observed memory excludes an uncoupled slow hidden mode.

For the second point, on a unit circle take a tensor equal to M on an
arc of length delta and zero elsewhere. A periodic H1 function with
slopes +1 and -1 on the two half-arcs, and zero slope elsewhere, has
energy-relative response M while the exterior average is delta M.
Taking delta=M^(-2) makes the average small and the response large.
Smooth approximations preserve the strict separation.

Run Python 3.12 without optimization:

```sh
python -B papers/yang-mills-certified-benchmark/certificates/ym91_conditional_memory_flux.py --check
python -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym91_conditional_memory_flux.py -v
```

The workflow also replays unchanged YM90 and YM87 certificates and
tests. The manifest pins 410 upstream inputs; the result includes those,
four local inputs, and the manifest, for 415 input hashes. Default and
--check are read-only. Explicit --write regenerates only this stage's
result and digest. Exact finite controls support the written analytic
proofs; they are not formal verification or independent peer review.

Conditional effective dynamics and variational Poisson estimates are
established methods; see F. Legoll and T. Lelievre,
[Effective dynamics using conditional expectations](https://arxiv.org/abs/0906.4865),
and the residual identity in [YM87](YM87_CONDITIONAL_ENERGY_TRANSPORT.md).
This development applies the conditional inverse order to the gauge
memory source, gives its flux certificate, and demonstrates the native
curvature/memory separation above. No priority claim is made for the
general conditional-expectation or minimum-flux principles.

The next actual Wilson task is to construct certified fluxes for its
full r(q,z), or solve and enclose its conditional Poisson tensor, with a
useful form bound along the desired running coupling. Numerical point
samples, unweighted exterior averages, and a comparison parent's flux
do not discharge that task. A separate physical hidden floor, surviving
continuum observables, and the continuum axioms are still required.
