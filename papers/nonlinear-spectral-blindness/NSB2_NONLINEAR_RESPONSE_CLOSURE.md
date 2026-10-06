# NSB2 — Nonlinear response closure and Hessian-power rigidity

Development: 6 October 2026. Continuation of NSB1, with the EMK carrier
and the YM82 scope retained. This packet supplies proofs and exact finite
controls. It makes no external priority claim and no continuum mass-gap claim.

## 1. What the master equation represents

In a declared Euclidean affine calibration, write

\[
H(x)=D^2U(x)=\frac{a+c}{2}I+\frac{a-c}{2}K+bRK
          =\begin{pmatrix}a&b\\b&c\end{pmatrix}.
\tag{1}
\]

Here R=[[0,-1],[1,0]], K=diag(1,-1), and RK=[[0,1],[1,0]].
The existing EMK1 file uses K0=[[0,1],[1,0]] instead: the exact adapter
is K=-R K0 and RK=K0. It changes the named basis, not the native algebra.
K^2=I, R^2=-I, and KR=-RK hold in both conventions.

Equation (1) spans every symmetric 2-by-2 response. Its coefficients can
already depend nonlinearly on the state. A linear expansion in an algebra
basis does not imply a linear constitutive law. Likewise, a polynomial of
a fixed linear operator is still a linear operator on its argument.

### T1. Powers do not supply an independent response channel

Cayley–Hamilton gives

\[
H^2=(a+c)H-(ac-b^2)I,\qquad
H+\lambda H^2=[1+\lambda\operatorname{tr}H]H
                -\lambda\det(H)I.
\tag{2}
\]

Every scalar polynomial p(H) lies in span{I,H} and commutes with H.
It preserves each H eigenspace; distinct eigenvalues may merge under p.
Thus two fields with identical eigenvalue fields still have identical
eigenvalue fields after every same declared polynomial. Powers alone
cannot repair the spectral blindness in NSB1. They can change eigenvalues
and the geometry of a matrix field; this is not a claim of no effect.

For two symmetric responses Hj=alpha_j I+beta_j K+b_j RK, instead,

\[
[H_1,H_2]=2(b_1\beta_2-\beta_1b_2)R.
\tag{3}
\]

Different channels or response derivatives can carry information absent
from powers of one response. Equations (2) and (3) are elementary matrix
identities, not new general algebra theorems.

## 2. T2 — A potential-preserving square extension is rigid

**Theorem.** Let U be C4 on an open set Omega in R^n and H=D^2U>0.
Fix a nonzero **constant** lambda. The following conditions are equivalent
as identities throughout Omega:

1. H+lambda H^2 is locally a Hessian in the same affine coordinates.
2. H^2 is locally a Hessian in those coordinates.
3. [H,partial_k H]=0 for every k.

If they hold, the Hessian response metric H is flat. Wherever
H+lambda H^2 is positive, that Hessian metric is also flat. On a region
with simple spectrum, the Euclidean eigenlines are constant locally and
U is locally a sum of one-variable potentials in a fixed orthogonal frame,
plus an affine function. Flatness of H alone does not imply these conditions.

**Proof.** Local Hessian integrability of a symmetric field M is
partial_k M_ij=partial_i M_kj. This follows in both directions by two local
applications of the Poincare lemma. Subtract H to obtain the first
equivalence. At any point choose a constant orthogonal coordinate change
diagonalizing H there, with entries h_i. Set C_ijk=U_ijk, fully symmetric.
Then

\[
\partial_k(H^2)_{ij}=(h_i+h_j)C_{ijk},\qquad
\partial_k(H^2)_{ij}-\partial_i(H^2)_{kj}
                    =(h_i-h_k)C_{ijk}.
\tag{4}
\]

By symmetry of C, vanishing of all these differences is equivalent to
(h_i-h_j)C_ijk=0, which is exactly [H,partial_k H]=0. This argument also
works at repeated eigenvalues.

For the flatness conclusion, the condition must hold on the open set.
Differentiate [H,H_i]=0 in direction j and subtract the equation with
i,j exchanged. Equality of mixed second derivatives of H yields
[H_i,H_j]=0. Since H also commutes with every H_i,

\[
[H^{-1}H_i,H^{-1}H_j]=H^{-2}[H_i,H_j]=0.
\tag{5}
\]

NSB1's Hessian Levi–Civita formula F_ij=-[H^-1 H_i,H^-1 H_j]/4
proves flatness. For M=H+lambda H^2, M_i=(I+2lambda H)H_i, so the same
commutations and formula apply. At simple spectrum, differentiating an
eigenvector equation shows that all off-diagonal changes of the
orthonormal eigenframe vanish. The Hessian is diagonal in a locally fixed
orthogonal frame, giving the stated local separation. QED.

No global separation across eigenvalue collisions is asserted. A condition
at one point does not imply the differentiated identity (5). For example,

\[
U=(x^2+y^2+xy^2)/2
\tag{6}
\]

has H(0)=I, so (4) vanishes at zero, but kappa(0)=1/4. Its square fails
integrability at nearby points. Conversely,

\[
U=x^2/2+(x+y)^2/2+(x+y)^3/6
\tag{7}
\]

is strictly convex on 1+x+y>0 and has H_x=H_y, hence zero curvature.
At zero its square has the two Codazzi defects (1,1), so the converse
to the flatness implication is false.

### Strictly convex or concave spectral functions

The same equivalence holds with H^2 replaced by f(H), for a C3 scalar
function f whose derivative is strictly increasing or strictly decreasing
on an interval containing the spectra. If f(H)>0, its metric is flat too.
Indeed the matrix derivative in an eigenbasis is

\[
\partial_k f(H)_{ij}=f^{[1]}(h_i,h_j)C_{ijk},
\quad f^{[1]}(s,t)=\begin{cases}(f(s)-f(t))/(s-t),&s\ne t,\\f'(t),&s=t.
\end{cases}
\tag{8}
\]

For fixed t this divided difference is strictly monotone in s, as follows
from f^[1](s,t)=integral_0^1 f'(t+r(s-t)) dr. Integrability forces each
nonzero C_ijk to have h_i=h_j=h_k, precisely condition 3. Conversely
that condition makes (8) fully symmetric. Differentiating the condition
as above proves flatness, also for the positive transformed metric.
The matrix derivative identity is classical; see source [2].

The nonlinear spectral-function hypothesis matters: f(t)=alpha t+beta
preserves Hessian structure for every H and need not give a flat metric.
Negative lambda is permitted in T2, but positivity of the transformed
metric must be checked separately.

### Variable coefficients and calibration

In dimension two a general field G=u(x)H+v(x)I is locally a Hessian exactly
when

\[
\nabla v=-\operatorname{adj}(H)\nabla u.
\tag{9}
\]

This follows by differentiating G_11,G_12,G_22 and cancelling the Hessian
identities of H. It is the appropriate integrability test for variable
coefficients. T2 does not cover arbitrary state-dependent lambda.

Squaring a covariant Hessian also needs the declared calibration. Under
x=A y, H_y=A^T H_x A, and generally H_y^2 is not A^T H_x^2 A.
Given a fixed positive background metric G0, the covariant version is
H G0^-1 H. T2 applies in a G0-orthonormal affine chart when G0 is constant.
It is not a background-independent tensor operation.

## 3. T3 — A nonlinear extension that preserves the potential

Change the energy itself: U_tilde=F(U), with a declared normalization of U.
Write g=dU, H=D^2U and C=D^3U. The chain rule gives

\[
\boxed{\widetilde H=F'(U)H+F''(U)\,g\otimes g,}
\tag{10}
\]
\[
\widetilde C_{ijk}=F'C_{ijk}
 +F''(H_{ij}g_k+H_{ik}g_j+H_{jk}g_i)+F'''g_i g_j g_k.
\tag{11}
\]

Thus F'>0 and F''>=0 suffice for strict convexity if H>0. The added
rank-one term uses calibrated force data. It need not commute with H
and cannot generally be computed from the eigenvalues of H.
Changing an affine term in U changes this construction, so energy and
force calibration must be part of the model. The chain rule itself is
standard; it is the admissible replacement for an unchecked Hessian square.

For several potentials, the more general exact expression is

\[
D^2F(U_1,\ldots,U_m)=\sum_a F_a D^2U_a
               +\sum_{a,b}F_{ab}\,dU_a\otimes dU_b.
\tag{12}
\]

All cross terms must be retained. Positivity follows, for example, if
the first sum is positive and the matrix (F_ab) is positive semidefinite.

As a control, take U0=(x^2+y^2)/2 and F(s)=s+lambda s^2, lambda>0.
Then

\[
\widetilde H=(1+\lambda r^2)I+2\lambda xx^T\ge I,\quad
\widetilde\kappa=-\frac{2\lambda^2r^2}
 {(1+\lambda r^2)^2(1+3\lambda r^2)^2}.
\tag{13}
\]

At lambda=1 and (x,y)=(1,0), H_tilde=diag(4,2) and kappa=-1/32.
Applying p to H0=I instead gives a constant scalar matrix. This radial
example is the m=2 case already present in NSB1; its role here is a
control of the new composition adapter, not a newly discovered example.

## 4. T4 — Closing the local spectral measurement deficit

Let H=D^2U>0 have **simple eigenvalues at a point**. Fix calibrated unit
eigenvectors e_i there and smooth local eigenvalue labels h_i. Define

\[
G_{ik}=\partial_{e_k}h_i,\qquad
C_{ijk}=D^3U(e_i,e_j,e_k).
\tag{14}
\]

**Reconstruction theorem.** The eigenvalues, calibrated eigenframe and
eigenvalue gradients determine every component with a repeated index:

\[
C_{iik}=G_{ik}.
\tag{15}
\]

Exactly binomial(n,3) independent components of the symmetric cubic tensor
remain: C_ijk with i<j<k. Supplying these marked mixed responses completes
D^3U and therefore the entire response curvature at that point. The count
is sharp for reconstruction of the cubic tensor from these data; it is
not a claim of a minimal measurement count for curvature alone.

**Proof.** Differentiate H e_i=h_i e_i and contract with e_i. Normalization
cancels the eigenvector derivative, giving (15). By Hessian symmetry this
also supplies every permutation of these components. A symmetric cubic
tensor has binomial(n+2,3) components; n^2 have a repeated index and the
remaining binomial(n,3) have three distinct indices. They are independently
realizable by adding sum_(i<j<k) t_ijk x_i x_j x_k to a positive diagonal
quadratic potential. Positivity persists in a small neighborhood. Their
Hessian value and all G_ik at zero are unchanged. Finally H and C determine
F_ij by NSB1's formula. QED.

This is a local reconstruction theorem, not a claim that any prescribed
spectral fields integrate to a potential on a neighborhood.

### Two dimensions: one vector probe closes the deficit

In dimension two there are no three-distinct-index components. Once the
eigenvalue fields and their first derivatives are given, a single response
Hv for any known nonzero vector v supplies the missing H calibration.
For e=v/|v| and w=R e, set

\[
a=e^T H e,\quad b=w^T H e,\quad c=h_1+h_2-a.
\tag{16}
\]

The two components of Hv determine a,b, hence H and its eigenframe up to
irrelevant sign choices. The probe is a **vector** measurement, not one
scalar quadratic form. No probe is insufficient in general by NSB1's
identical-entire-spectrum flat/curved pair.

In that eigenframe, write p=G_11, q=G_12, r=G_21, s=G_22. Then

\[
\boxed{\kappa=-\frac{h_1(qs-r^2)+h_2(pr-q^2)}{4h_1^2h_2^2}.}
\tag{17}
\]

Substitution into the 2D Hessian curvature determinant proves the formula.
For NSB1's pair at c=1,(x,y)=(2,0), h=(1,1/2) and the physical eigenvalue
gradients are (0,0),(1/4,0) for both metrics. The probe v=(1,2) returns
(1/2,2) in the flat case and (1,1) in the curved case. After rotation to
the measured eigenframes, (17) gives 0 and 1/16 respectively.
For comparison, the scalar probe v^T H v with v=(1,1) returns 3/2 in
both cases and leaves that curvature difference undetected.

This repairs the pointwise observability limitation under the stated
Hessian and simple-spectrum assumptions. It does not remove the
finite-jet/neighborhood limitation in NSB1. Eigenvalue collisions require
direct cubic data or additional hypotheses; the algorithm rejects a
repeated-eigenvalue input and does not assert uniform conditioning there.

### Dimension three and four: the remaining mixed data are real

For

\[
U_\tau=(x^2+2y^2+3z^2)/2+\tau xyz
\tag{18}
\]

the origin has H=diag(1,2,3), G=0 for every tau, but
K_12=tau^2/24. For tau=1, H>=I/2 on the cube |x|,|y|,|z|<=1/4,
by the maximum row sum bound on the off-diagonal perturbation.
Thus known H and eigenvalue first derivatives alone are insufficient
already in dimension three. In dimension four the cubic deficit is four
components: C_123,C_124,C_134,C_234. This concerns four **constitutive
state coordinates**, not a reconstruction of a four-dimensional quantum
Yang–Mills vacuum.

## 5. Toolkit and Yang–Mills boundary

The new adapter supplies EMK polynomial reduction, a Hessian-integrability
guard, potential-composition jets, calibrated spectral reconstruction in
rational orthonormal frames, and marked mixed-derivative reconstruction.
It reuses the existing EMK exact matrix operations and NSB1/NSR1 geometry.
Polynomial identities, independent Christoffel calculations and invalid
inputs are checked. Finite checks supplement the proofs; they do not
constitute formal verification or establish a field identity from samples.

YM82's Phi(W)=theta W+lambda W^2 is a nonlinear **plaquette potential**
inside a quantum Hamiltonian. It is not the operation H -> H+lambda H^2
on the Hessian response metric, nor a squaring of the quantum Hamiltonian.
For the homogeneous cubic lattice, its stated sufficient region is
osc(Phi)/2<1/8000. Since osc(Phi)/2>=|theta|, adding this quadratic term
does not extend that criterion to theta->infinity on the weak-coupling
trajectory. Other YM82 block results have their own hypotheses. None
is an unconditional gap result for arbitrary nonlinear coefficients.

NSB2 closes a specified **local response-data** deficit and rules out an
incorrect potential extension. It supplies no new Wilson weak-coupling
gap bound, continuum measure, or continuum state-survival estimate.
The global vacuum and uniform estimates required in YM91 remain open.

## 6. Provenance and replay

Parent commit: ff2b35e8c8749da00d212bdb820217e3faf8c5c2 (YM91).
The source manifest pins the inherited chain and the EMK adapter sources.
Existing source and evidence files are left byte-for-byte unchanged.

Primary background, read 6 October 2026:

1. B. Totaro, *The curvature of a Hessian metric*, International Journal
   of Mathematics 15 (2004), 369–391, Section 2:
   <https://www.math.ucla.edu/~totaro/papers/public_html/hessian.pdf>.
   The third-derivative curvature formula is established Hessian geometry.
2. V. Noferini, *A formula for the Frechet derivative of a generalized
   matrix function*, Manchester MIMS eprint 2016.24 (2016), Section 2.4,
   Theorem 2.11 (the classical Daleckii–Krein theorem):
   <https://eprints.maths.manchester.ac.uk/2544/>.
3. Repository NSB1 Sections 4–7; EMK1 and EMK-C1; YM82 Section 2 and its
   exact quadratic region. The manifest records the exact versions.

The packet combines standard algebra/geometry with the explicitly proved
rigidity and reconstruction statements for this measurement contract.
A literature-wide novelty assessment has not been established.

```bash
python -B papers/nonlinear-spectral-blindness/certificates/nsb2_nonlinear_response.py --check
python -B -m unittest discover -s papers/nonlinear-spectral-blindness/tests -v
python -B papers/nonlinear-spectral-blindness/certificates/nsb1_nonlinear_blindness.py --check
python -B papers/yang-mills-certified-benchmark/certificates/ym82_nonlinear_lattice_gap.py --check
```
