# N4H-1 — Four-dimensional metric gluing, lifted curvature and information response

14 October 2026. Next development after NJC-1.

The previous certificate left an arbitrary cut-even connection. This packet
adds a precise geometric contract that fixes discrete transport, constructs
actual four-dimensional simplices, solves one interior-edge gravitational
equation, and determines when a curvature phase needs information beyond a
compact holonomy. No external time parameter is introduced.

## 1. What is inherited and what is added

NJC-1 T3 and CID-D1 already establish

\[
A_i=\tfrac12 J\partial_iJ+B_i,\qquad [B_i,J]=0.
\]

NJC-1 T8 demonstrates that the complete cut derivative tower does not
identify an independently allowed normal connection. That result remains
valid. Here we restrict to a new adapter: a piecewise Euclidean metric on an
oriented simplicial four-manifold, nondegenerate four-simplices, complete
matching of shared tetrahedra, and the metric transport obtained by
unfolding neighboring cells. These additional geometric data select the
discrete Levi-Civita transport. The continuous connection formula is not
silently imposed on every hinge cut simultaneously.

The native rational cut field and Pāṇinian word compiler are unchanged.
Coordinates, positive Gram forms, oriented cell incidence and analytical
angle completion are explicit interpretation contracts. No primitive
Hilbert-space axiom is used. This is a Euclidean geometric sector; a
Lorentzian spacetime interpretation is a separate extension.

All lengths are divided by a declared reference length L. Write
q_e = length(e)^2/L^2, a_h = area(h)/L^2 and v_sigma = volume(sigma)/L^4.
These and all angles are dimensionless. Rescaling the geometry and L
together changes none of them. Dimensionless variables do not determine
the absolute scale L or a gravitational coupling.

## H1. Complete tetrahedral matching fixes the transport

Let E and F be the 4-by-3 matrices of corresponding edge vectors of a
shared tetrahedron in the two simplex frames. Assume

\[
E^\dagger E=F^\dagger F=G>0.
\]

Let n_s and n_t be the unit normals to the shared tetrahedron pointing into
the source and target simplices. The unique isometry mapping every shared
edge correctly and placing the two interiors on opposite sides is

\[
\boxed{T_{ts}=FG^{-1}E^\dagger-n_tn_s^\dagger.}\tag{1}
\]

**Proof.** It maps E to F and n_s to -n_t. The source decomposes
orthogonally into the image of E and the line spanned by n_s. The formula
is an isometry on each summand, and their images are orthogonal. These four
independent directions span the entire carrier, proving uniqueness.
With consistently oriented simplex frames, the unfolding map lies in SO(4).
The inverse is T_st. Independent frame changes O_s,O_t give
T_ts -> O_t T_ts O_s^dagger. Loop products therefore transform by conjugacy.

An arbitrary extra rotation of a triangle's normal plane would preserve
that triangle but generally fail to match the third direction in the
shared tetrahedron. This is how the new metric data remove the earlier
freedom. The cut alone still does not remove it. Using a plus sign before
n_t n_s^dagger gives another orthogonal map that folds the two interiors
onto the same side; the certificate explicitly rejects it.

## H2. An exact four-dimensional hinge star

Take a triangular hinge ABC with A=0, B=e_1, C=e_2. For a cyclic sequence
D_0,...,D_(n-1), use n four-simplices

\[
\sigma_i=[A,B,C,D_i,D_{i+1}],\quad
D_i=e_3,\quad D_{i+1}=c_i e_3+s_i e_4
\quad\text{in frame }i,\tag{2}
\]

where c_i^2+s_i^2=1 and s_i>0. Indices are cyclic. Coordinates of a named
vertex belong to its particular simplex frame; no common flat embedding
of a curved star is assumed. Every repeated edge length agrees. Shared
tetrahedra have the same Gram matrix, and each simplex has positive
four-volume s_i/24. Thus these are genuine four-simplices, rather than
abstract rotation matrices without compatible cell geometry.

Let P=diag(1,1,0,0), Q=I-P and let R be the oriented quarter-turn on the
normal plane, with R e_3=e_4 and R e_4=-e_3. Applying (1) gives

\[
T_{i+1,i}=P+c_iQ-s_iR.
\]

For the interior dihedral angles theta_i in (0,pi), c_i=cos(theta_i) and
s_i=sin(theta_i). Their loop transport is

\[
\boxed{H_h=P+\cos\delta_h\,Q+\sin\delta_h\,R,
\qquad \delta_h=2\pi-\sum_i\theta_i.}\tag{3}
\]

This is curvature concentrated on a triangular hinge. The formula is
unchanged by a common scaling of all coordinates. Native exact fixtures
use c=(1-t^2)/(1+t^2), s=2t/(1+t^2) for positive rational t.

The invariant defect

\[
\kappa_h=\tfrac14\|H_h-I\|_F^2=1-\cos\delta_h\tag{4}
\]

loses the sign of the oriented angle. The oriented bivector R together
with H also determines sin(delta). Even these two components determine
delta only modulo 2pi. Orientation reversal changes R and the signed
angle together; no external frame is selected by (3).

## H3. Identical compact holonomy can conceal different curvature lifts

Use the same six-cell complex (2), with the same hinge area a_h=1/2.
Compare these positive rational Cayley parameters:

| Geometry | Parameters around the hinge | Sum of interior angles | Deficit | H_h |
|---|---|---|---|---|
| Flat lift | (1/2,1/2,3/4) repeated twice | 2pi | 0 | I |
| Excess-angle lift | (2,2,4/3) repeated twice | 4pi | -2pi | I |

**Exact branch proof.** Put cos(alpha)=3/5, sin(alpha)=4/5. Then
pi/4<alpha<pi/2. The third angle beta has cosine 7/25 and sine 24/25,
so beta=pi-2alpha and 0<beta<pi/2. The first triple therefore sums to pi.
Replacing t by 1/t replaces each angle by its supplement. The second
triple sums to 3pi-pi=2pi. Doubling the triples gives the table.
The native engine checks the rational products and the strict range
inequalities; the interval argument specifies the real angle lifts.

Both stars have identity holonomy and zero defect (4), but the second has
a nonzero conical curvature integral. Its larger total angle is not a
degenerate simplex: every individual simplex remains nondegenerate.
The cell metrics distinguish the two configurations, although the single
coarse loop observation does not. Their boundary metric data differ; this
pair is not presented as two histories with the same physical boundary.

Consequently an information record adequate for integrated curvature must
retain the oriented dihedral-angle sum, or equivalent lift data. A phase
or norm alone is insufficient. This is an application of the previous
observation-factorization principle, not a new claim that all holonomies
lose all curvature.

## H4. A clock-free metric and information action

For a finite triangulated region with fixed boundary edge lengths, select
the dimensionless curvature functional

\[
\mathcal R(q)=\sum_{h\text{ interior}}a_h(q)\delta_h(q)
+\sum_{h\text{ boundary}}a_h(q)\psi_h(q),\qquad
\psi_h=\pi-\sum_{\sigma\supset h}\theta_{\sigma h}.\tag{5}
\]

These are the standard Regge bulk and boundary terms, reused here as a
gravitational adapter. The sign convention in this packet is +R. A
cosmological volume term is set to zero in this first model. The metric
and the cut algebra alone do not uniquely derive this choice of action.

A concrete gauge-invariant information mismatch can use dimensionless
carrier data z_sigma on neighboring cells:

\[
\mathcal I(q,z)=\frac12\sum_{\{\sigma,\tau\}}w_{\sigma\tau}(q)
\|z_\tau-T_{\tau\sigma}(q)z_\sigma\|^2,
\qquad w_{\sigma\tau}>0.\tag{6}
\]

Each unoriented dual edge is counted once. Pairings and state normalization
are declared. This defines an operational mismatch functional; it is not
an unproved identification with thermodynamic or Shannon entropy. Both
normal and retained channels are available. On (2), data z_i=e_3 and
w=1 give I=sum_i(1-c_i); data z_i=e_1 give zero for every choice of angles.
Thus a retained-only readout can miss the normal transport completely.
The certificate verifies (6) under different frame rotations in every cell.

Select

\[
\boxed{\mathcal S(q,z)=\gamma\mathcal R(q)+\alpha\mathcal I(q,z),
\qquad \gamma>0.}\tag{7}
\]

For interior squared-edge variables the stationary equations are

\[
\boxed{\gamma\sum_{h\supset e}\delta_h
\frac{\partial a_h}{\partial q_e}
+\alpha\frac{\partial\mathcal I}{\partial q_e}=0,
\qquad \nabla_{z_\sigma}\mathcal I=0.}\tag{8}
\]

The derivative of I includes the geometry dependence of transport and
weights; holding T fixed while varying q would be a different model.
Vary only unfixed interior data. The carrier equation in the real pairing
is the covariant graph equation

\[
\sum_{\tau\sim\sigma}w_{\sigma\tau}
(z_\sigma-T_{\sigma\tau}z_\tau)=0.
\]

No external time variable occurs: (8) relates geometry and transported
data on the whole region. This is a conditional stationary dynamics law,
not yet a Hamiltonian evolution or a quantum constraint algebra.

### Why the angle-derivative terms cancel

For a nondegenerate Euclidean four-simplex let V be its volume, V_i its
tetrahedral facet volumes, n_i its outward unit normals, and A_ij the
area of the triangle shared by facets i,j. Closure gives
sum_i V_i n_i=0. Base-height geometry gives

\[
4VA_{ij}=3V_iV_j\sin\theta_{ij}.
\]

Since n_i dot n_j=-cos(theta_ij),

\[
\begin{aligned}
\sum_{i<j}A_{ij}\,d\theta_{ij}
&=\frac{3}{4V}\sum_{i<j}V_iV_j\,d(n_i\cdot n_j)\\
&=\frac{3}{4V}\sum_i V_i\,dn_i\cdot
\sum_{j\ne i}V_jn_j=0.
\end{aligned}
\]

The last equality uses dn_i dot n_i=0. This is the flat Schläfli identity.
Apply it to every simplex in the variation of (5). Boundary triangle
areas are fixed, but their angles must be included in the cancellation.
The surviving term is sum_h delta_h da_h, giving (8). This argument does
not require all hinge areas to be independently variable.

The certificate also checks this cancellation separately for all ten
hinges of a four-simplex at four exact edge geometries. If B_ij denotes
the Gram matrix of the five affine barycentric gradients, it verifies

\[
A_{ij}^2=144V^2(B_{ii}B_{jj}-B_{ij}^2),\qquad
\frac{A_{ij}\theta'_{ij}}{12V}
=B'_{ij}-\frac{B_{ij}}2
\left(\frac{B'_{ii}}{B_{ii}}+\frac{B'_{jj}}{B_{jj}}\right).
\]

The sum on the right is exactly zero by the row sums of B and B'. This
avoids approximating square roots or inverse trigonometric functions.

## H5. A solved four-dimensional interior-edge model

The triangulation consists of four four-simplices:

\[
ABCDE,\quad ABCDF,\quad ABCEF,\quad ABDEF.
\]

It is the join of edge AB with the tetrahedral two-sphere CDEF. There is
exactly one interior edge AB, four interior triangles ABC,ABD,ABE,ABF,
six internal tetrahedra and eight boundary tetrahedra. Each interior
triangle belongs to three four-simplices. This model has a genuine bulk
metric variable; the single-hinge stars above served a different purpose.

Fix all eight squared distances from A or B to C,D,E,F at 4, and all six
squared distances among C,D,E,F at 8. Set x=q_AB. In every simplex the
edge Gram matrix based at A is

\[
G(x)=\begin{pmatrix}
x&x/2&x/2&x/2\\
x/2&4&0&0\\
x/2&0&4&0\\
x/2&0&0&4
\end{pmatrix}.
\]

Its transverse Schur complement has eigenvalues 4,4,4-3x/4. Therefore
the admissible open domain is precisely

\[
0<x<16/3.
\]

All four internal triangle areas and their interior dihedral angles agree:

\[
a(x)=\tfrac14\sqrt{x(16-x)},\qquad
\cos\theta(x)=-\frac{x}{2(8-x)},\qquad
\delta(x)=2\pi-3\theta(x),\qquad
v(x)^2=\frac{x(16-3x)}{144}.\tag{9}
\]

**Derivation.** Project the three outer edge vectors perpendicular to AB.
Their Gram matrix has diagonal 4-x/4 and off-diagonal -x/4. Project two
of these vectors perpendicular to the third as well. Their angle has the
cosine in (9). Equivalently use the two corresponding inverse-Gram
facet normals. The two-edge Gram determinant of ABX gives a(x)^2; the
four-edge determinant divided by (4!)^2 gives v(x)^2.

Including the boundary action in (5), Schläfli cancellation gives

\[
\boxed{\mathcal R'(x)=4a'(x)\delta(x)
=\frac{8-x}{\sqrt{x(16-x)}}\delta(x).}\tag{10}
\]

Although the boundary lengths are fixed, omitting their angle terms from
the action would leave extra derivatives and spoil (10).

In vacuum, R'=0. The prefactor in (10) is positive throughout the
admissible domain. Hence delta=0, theta=2pi/3 and cos(theta)=-1/2. Equation
(9) then yields the **unique vacuum solution x=4**. It is a flat
configuration, not a claim that every gravitational vacuum is flat.

An explicit common embedding at this solution is

\[
\begin{array}{lll}
A=(0,0,0,-1),&B=(0,0,0,1),&C=(1,1,1,0),\\
D=(1,-1,-1,0),&E=(-1,1,-1,0),&F=(-1,-1,1,0).
\end{array}
\]

The four simplices triangulate this bipyramid and each has volume 1/3.
The native engine checks all four complete Gram matrices and volumes.

For comparison, x=48/11 is admissible but curved. Its interior cosine is
-3/5 and positive sine is 4/5; a loop in an adapted hinge frame has normal
cosine 117/125 and sine -44/125. It fails the vacuum equation. It could
only be stationary with an appropriate source or a different action.

### Exact local information response

At x=4,

\[
\mathcal R'(4)=0,\qquad \boxed{\mathcal R''(4)=-\tfrac12.}\tag{11}
\]

Indeed cos(theta)'=-4/(8-x)^2. At the flat solution a=sqrt(3),
sin(theta)=sqrt(3)/2 and a sin(theta)=3/2. Differentiate (10), using
delta=0, to obtain (11). The native check computes the cosine derivative
independently by differentiating the inverse Gram matrix.

Let epsilon=alpha/gamma. For a prescribed twice differentiable information
functional I(x), the nonzero derivative (11) and the implicit function
theorem give a unique nearby stationary branch for sufficiently small
epsilon:

\[
\boxed{x=4+2\epsilon\mathcal I'(4)+O(\epsilon^2).}\tag{12}
\]

Equation (12) also applies to a twice differentiable on-shell I_eff(x)
if the carrier equations have been solved on a chosen smooth branch.
Existence and smoothness of that branch must be established for the
chosen data. Unspecified matter dynamics is not assumed solved.

This is a concrete metric response to a dimensionless information force.
The positive coefficient 2 is a property of this boundary geometry and
the sign convention (7), not a universal gravitational constant. The
negative Hessian means the vacuum is a local maximum of +R along this
edge variable; stationarity alone is not an energy-minimum or stability
claim. Adding I or adopting a different overall action sign changes the
stability question.

## H6. A quantum phase detects the discarded winding

After declaring the usual completed central phase exp(iota s), a possible
single-hinge action phase is

\[
\chi_h(\delta_h)=\exp(\iota q_h\delta_h),\qquad q_h=\gamma a_h.
\tag{13}
\]

This is an amplitude prescription from the selected action. An orthogonal
geometric holonomy is not by itself a quantum evolution law.

**Descent theorem.** For a fixed real weight q, (13) depends only on the
SO(2) holonomy, for every real angle lift, if and only if q is an integer.

**Proof.** Equivalent angle lifts differ by 2pi k, k an integer. The phase
changes by exp(2pi iota q k). This is one for all k exactly when q is an
integer, by the kernel 2pi Z of the completed circle exponential. Integer
weights give the Laurent characters u^q of u=exp(iota delta); their inverse
and composition identities have native word certificates.

The two geometries in H3 have the same a_h=1/2 and H_h=I. At gamma=1,
their single-hinge phases are +1 and -1. This is an exact example of
information discarded by a compact holonomy but retained by a selected
action phase. It is not a computed full boundary transition amplitude:
other hinges, the boundary action and the measure have not been integrated.

There are therefore two distinct modeling choices:

* Retain the real curvature lift. No integer condition on gamma a_h follows.
* Demand a single-valued phase on the unlifted SO(2) record alone. Then
  gamma a_h must be integral for each independently identified angle lift.

A specified double cover instead admits half-integer weights. In a full
constrained lattice, independent changes of hinge lifts need not be
allowed; the theorem states its observation domain explicitly. It does
not establish a universal area spectrum or quantization from the cut
alone. This conditional phase criterion is the precise point where the
present circle structure meets a possible geometric quantization rule.

Replacing (5) by a defect sum sum_h a_h(1-cos(delta_h)) would lose the lift
and begin quadratically in small curvature. It is a different action; a
unitary loop and a positive norm do not prove equivalence to gravity.

## 2. Certification boundary and lineage

The executable packet supplies exact rational geometry, complete matching
checks, native rewrite witnesses, independent frame transformations,
all-ten-hinge Schläfli checks and negative controls. The general geometric
proofs, interval choices, analytic implicit-function argument and circle
descent theorem are written proofs with explicit hypotheses. This is not
a Lean/Coq formalization of real analysis or a certificate of quantum
gravity as a physical theory.

| Ingredient | Status in this packet |
|---|---|
| Cut-compatible connection and free B | Reused NJC-1 / CID-D1 |
| Exact circular transport arithmetic | Reused native engine and NJC-1 bridge |
| Shape-matched metric gluing and Regge action | Standard geometric adapter, explicitly added |
| Complete rational 4D hinge fixtures and lift counterexample | New executable development in this research sequence |
| One genuine interior-edge vacuum and source response | Derived and checked for the specified four-cell boundary data |
| Integer phase criterion | Standard circle-character fact applied conditionally to the retained curvature record |
| Quantum measure, physical observables, Lorentzian limit | Not constructed here |
| Full gravitational constraint algebra and continuum recovery | Open; not supplied by commuting auxiliary constraints |

Primary references for the reused Regge variation and its limitations:

1. B. Bahr and B. Dittrich, *Regge calculus from a new angle*,
   [arXiv:0907.4325](https://arxiv.org/abs/0907.4325), equations (2)-(3)
   and Appendix A. The flat Schläfli specialization is used here.
2. B. Bahr and B. Dittrich, *(Broken) Gauge Symmetries and Constraints in
   Regge Calculus*, [arXiv:0905.1670](https://arxiv.org/abs/0905.1670),
   section 3 and Appendix B. A finite curved Regge discretization does not
   automatically retain exact continuum diffeomorphism constraints.

Local predecessors are hash-pinned in SOURCE_PINS.json. Earlier results
are reused, and no priority or literature-wide novelty claim is made.
