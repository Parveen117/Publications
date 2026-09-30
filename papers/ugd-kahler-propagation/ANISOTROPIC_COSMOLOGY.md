# Native spin, rotating shear and an exact anisotropic cosmology

Author: Monty Dabas. Research development R8, 30 September 2026.

This note extends the coupled solution of [CS](COUPLED_COSMOLOGICAL_SOLUTION.md)
to a general spatial metric of Bianchi type I. The native spin source produces
anisotropic stress. Its commutator with the shear rotates the shear axes while
preserving the squared volume-rescaled shear. In the homogeneous rest sector,
the volume, matter field, independent connection and full spatial coframe can
all be written explicitly. Their expanding positive-Lambda branch isotropizes
in expansion rate. Linearization also supplies the five anisotropic metric
modes missing from the homogeneous FLRW analysis in
[PS](PERTURBATION_STABILITY.md).

The action, classical commuting field, zero-Holst torsion law and neutral
coupling g=0 are exactly those of [SM](NATIVE_MASSIVE_MATTER.md). Throughout,
kappa>0, Z>0, m>=0, and the principal solved family has Lambda>0. This is an
extension within that declared theory. Classical Einstein–Dirac Bianchi I
solutions and their nondiagonal stress have established precedents, credited
in BI-8; the general mechanism is not claimed as a first discovery.

## BI-1. General homogeneous coframe and the missing stress

Use cosmic time and a real spatial coframe matrix E(t) with det(E)>0:

\[
e^0=dt,\qquad e^i=E^i{}_j(t)dx^j,\qquad
g=-dt^2+dx^TE^TE\,dx,\qquad v=\det E.
\tag{BI.1}
\]

Choose the spatial orthonormal frame parallel transported along the normal
time lines with respect to the Levi-Civita connection. Its expansion matrix
is symmetric:

\[
K=\dot E E^{-1}=K^T=hI+\sigma,\qquad
h=\tfrac13\operatorname{tr}K=\frac{\dot v}{3v},\qquad
\sigma^T=\sigma,\quad\operatorname{tr}\sigma=0.
\]

This is a rotational frame choice, not a diagonal metric restriction. A
time-dependent orthogonal rotation removes the antisymmetric part of
dot(E) E^(-1). The letter K in this note denotes expansion; the native cut
matrix and the contorsion tensor retain their distinct meanings in SM.
Write contorsion as C_abc here.

The Levi-Civita connection is
omega_LC^{0i}=K_ij e^j, omega_LC^{ij}=0. Therefore

\[
D_i^{\rm LC}\Psi=-\tfrac12\sum_jK_{ji}\mathcal A_j\Psi,\qquad
\gamma^\mu D_\mu^{\rm LC}\Psi
=\gamma^0(\dot\Psi+\tfrac32h\Psi)
\tag{BI.2}
\]

for a spatially homogeneous matter field. Symmetry of K cancels every
trace-free contribution to the slashed equation. In particular PS.1–PS.2
remain valid for chi=sqrt(v) Psi with all eight components retained.

However, the spatial stress has an additional term. Define the antisymmetric
three-by-three matrix

\[
Q_{ij}=\bar\Psi\gamma^{0ij}\Psi\quad(i,j=1,2,3).
\]

Direct substitution into the Hilbert stress CS.8 gives, in this orthonormal
frame,

\[
\boxed{T_{00}=\rho,\qquad T_{0i}=0,\qquad
T_{ij}=p\delta_{ij}+\pi_{ij},\qquad
\pi=\frac Z8[K,Q]=\frac Z8[\sigma,Q].}
\tag{BI.3}
\]

Here rho and p are precisely PS.4, including its contact interaction. For
the sign and coefficient, the spatial derivative part is

\[
-\frac Z4(\bar\Psi\gamma_iD_j\Psi+
                  \bar\Psi\gamma_jD_i\Psi)
=-\frac Z8\sum_k(K_{kj}Q_{ik}+K_{ki}Q_{jk})
=\frac Z8[K,Q]_{ij}.
\]

Repeated spatial indices in the triple give a skew quadratic form and
vanish. The on-shell scalar Lagrangian supplies p delta_ij; the mixed
components vanish by the same skew bilinears as in PS-1. Since K is symmetric
and Q antisymmetric, pi is symmetric and trace-free, and

\[
\operatorname{tr}(\sigma\pi)=0.
\tag{BI.4}
\]

Thus this stress changes the shear orientation but does no shear work in
the homogeneous energy balance. The matter continuity equation remains
dot(rho)+3h(rho+p)=0. Setting pi=0 without checking its commutator would
discard genuine coframe equations.

## BI-2. Conserved spin and a shear equation with constant eigenvalues

Define the volume-rescaled spin matrix and its fixed coefficient

\[
Q_0=vQ,\qquad (Q_0)_{ij}=\chi^TH\gamma^{0ij}\chi,\qquad
\Omega=\frac{\kappa Z}{8}Q_0.
\tag{BI.5}
\]

Each H gamma^{0ij} commutes with both gamma^0 and gamma^0 Jcal. These are
the two skew generators of PS.1. Differentiation consequently proves
dot(Q_0)=0 for every homogeneous solution, including fields outside the
rest sector. Omega is a constant antisymmetric matrix in the chosen
parallel frame. SM's Fierz identity and H gamma^{123}=-Ical Jcal also give

\[
\sum_{i<j}(Q_0)_{ij}^2=\Sigma^2+\Pi^2+B^2,
\tag{BI.6}
\]

with the PS bilinears. In the rest sector this equals N^2, where
N=chi^T chi=2n_0/Z>0.

The curvature of (BI.2) gives
R_00=-tr(dot K+K^2), R_0i=0 and R_ij=dot K_ij+3h K_ij. The Einstein
equations are therefore

\[
3h^2=\Lambda+\kappa\rho+\sigma^2,\qquad
2\dot h+3h^2+\sigma^2=\Lambda-\kappa p,\qquad
\sigma^2:=\tfrac12\operatorname{tr}(\sigma^2_{\rm matrix}),
\tag{BI.7}
\]
\[
\boxed{\dot\sigma+3h\sigma=\kappa\pi
                =\frac1v[\sigma,\Omega].}
\tag{BI.8}
\]

The scalar notation sigma^2 in BI.7 means half the matrix-square trace;
matrix products elsewhere are explicit. The momentum constraints vanish
because the geometry is homogeneous and T_0i=0. The lapse constraint in
BI.7 is retained. With the spatial equations, its defect obeys dot(C)=-3h C;
on the constraint, dot(h)=-kappa(rho+p)/2-sigma^2.

Set S=v sigma and introduce the clock u(t)=integral from t* to t of
d(tau)/v(tau), on any interval with v>0. Equation (BI.8) becomes

\[
\frac{dS}{du}=[S,\Omega],\qquad
\boxed{S(t)=O(u)S_*O(u)^T,\qquad O(u)=e^{-\Omega u}.}
\tag{BI.9}
\]

This is an orthogonal conjugation, so all eigenvalues of S are constant.
In particular, for any symmetric trace-free initial S_*,

\[
\boxed{D^2:=\tfrac12\operatorname{tr}(S_*^2),\qquad
\sigma^2=\frac{D^2}{v^2},\qquad
\|\sigma\|_F=\frac{\|S_*\|_F}{v}.}
\tag{BI.10}
\]

The Frobenius pairing makes the linear map S -> [S,Omega] skew on the
five-dimensional space of symmetric trace-free matrices. These statements
hold before restricting the matter to the rest sector. The remaining
volume and bilinear equations need not have the same elementary solution
for arbitrary nonrest data.

## BI-3. Closed nonlinear volume and native matter phase

Now choose H Xi=Xi, Xi^T Xi=N=2n_0/Z and the invariant homogeneous rest
sector. The matter and isotropic pressure retain their CS form:

\[
\Psi=v^{-1/2}e^{-\mathcal I\theta}\Xi,\qquad
n=\frac{n_0}{v},\qquad
\rho=mn+\frac{3\kappa n^2}{16},\qquad
p=\frac{3\kappa n^2}{16}.
\]

Equation (BI.7) gives

\[
\boxed{\dot v^2=3\Lambda v^2+3\kappa mn_0v
                 +\frac{9\kappa^2n_0^2}{16}+3D^2.}
\tag{BI.11}
\]

For Lambda>0 define

\[
\omega=\sqrt{3\Lambda},\quad
q=\sqrt{\frac{9\kappa^2n_0^2}{16}+3D^2},\quad
d=\frac{\kappa mn_0}{2\Lambda},\quad b=\frac q\omega.
\]

On the expanding branch with past volume zero at t=0, the exact solution is

\[
\boxed{v(t)=d(\cosh\omega t-1)+b\sinh\omega t,\qquad t>0.}
\tag{BI.12}
\]

An arbitrary time origin follows by replacing t with t-t_0. Direct
differentiation gives dot(v)^2 as in BI.11 and
ddot(v)=3 Lambda v+3 kappa m n_0/2. Since d>=0, b>0, both v and dot(v)
are positive for t>0. At D=0 this recovers CS exactly.

With z=exp(omega t) and z*=exp(omega t*), integration gives

\[
\boxed{u(t)=\frac1q\log\!left[
\frac{(z-1)[(d+b)z_*-(d-b)]}
     {(z_*-1)[(d+b)z-(d-b)]}\right],}
\tag{BI.13}
\]
\[
\boxed{\theta(t)-\theta(t_*)=m(t-t_*)+
                         \frac{3\kappa n_0}{8}u(t).}
\tag{BI.14}
\]

Thus the same explicitly integrated inverse-volume clock controls the
contact part of the native phase and the spin-induced shear rotation.
The fields remain real: the phase exponential uses SM's real Ical with
Ical^2=-I.

## BI-4. The full coframe and independent torsion connection

A solved volume and shear must still be lifted to a spatial metric. Let
E_* be any real initial coframe with det(E_*)=v(t*)>0, and set

\[
\boxed{E(t)=\left(\frac{v(t)}{v(t_*)}\right)^{1/3}
 e^{-\Omega u(t)}e^{(S_*+\Omega)u(t)}E_*.}
\tag{BI.15}
\]

No simultaneous diagonalization of S_* and Omega is assumed. To verify
the ordering, put F=exp(-Omega u) exp((S_*+Omega)u). Then

\[
\frac{dF}{du}=O S_* O^T F=S(u)F,
\]

because Omega commutes with O. It follows that
dot(E) E^(-1)=hI+S(u)/v, precisely the symmetric expansion matrix of
BI.1. The exponentials are invertible for every finite u. Their determinants
are one since tr(S_*)=tr(Omega)=0, so det(E)=v. Consequently
g_spatial=E^T E is positive definite for every t>0.

The complete independent connection is

\[
\omega^{ab}=\omega_{\rm LC}^{ab}+C^{ab}{}_c e^c,\qquad
C_{abc}=\frac{\kappa Z}{8}\bar\Psi\gamma_{abc}\Psi.
\tag{BI.16}
\]

Equations BI.12–BI.16, with BI.9, solve the original matter, connection and
coframe equations: BI.2 supplies the matter equation, BI.16 is SM's unique
algebraic connection solution, and BI.7–BI.8 supply every metric component.
SM's elimination equivalence then lifts them to the independent-connection
system. The neutral F=0 Maxwell equations hold as in CS.

This supplies a real analytic nonlinear family with arbitrary finite
symmetric trace-free S_*. No diagonal metric assumption or spin averaging
has been introduced. It is a solved spatially homogeneous rest family,
not a general spatially varying Einstein–matter solution.

## BI-5. Diagonal restrictions and scalar information loss

In a fixed diagonal Bianchi I coframe, the off-diagonal geometric Einstein
components vanish. The omitted equations would require

\[
(h_i-h_j)Q_{ij}=0\quad(i\ne j).
\tag{BI.17}
\]

Equivalently the shear must commute with the constant spin generator.
In the nonzero rest sector BI.6 forbids Q_0=0. Three distinct principal
expansion rates therefore cannot remain diagonal in that parallel frame.
For nonzero Omega, a symmetric trace-free matrix commuting with it has
equal eigenvalues in the plane orthogonal to its rotation axis: the
permitted nonzero diagonal branch is locally rotationally symmetric and
aligned with the spin axis. The isotropic member is also permitted.
Equation BI.15 includes the general rotating, nondiagonal solutions that
this restriction would otherwise remove.

The volume and several scalar readouts also fail to determine the relative
spin/shear alignment. For an explicit pair, take the CS rest vector
Xi=(1,0,0,0,0,0,-1,0), kappa=1, Z=4/3, n_0=4/3. Then

\[
\Omega=\begin{pmatrix}0&0&-1/3\\0&0&0\\1/3&0&0\end{pmatrix}.
\]

The two initial shears

\[
S_A=\operatorname{diag}(1,-2,1),\qquad
S_B=\operatorname{diag}(1,1,-2)
\]

both have D^2=3. At the same m, Lambda and time origin they have identical
v, rho, p, sigma^2, metric scalar curvature, full connection scalar curvature
and torsion contraction. Nevertheless

\[
[S_A,\Omega]=0,\qquad
\operatorname{tr}([S_B,\Omega]^2)=2,
\]
\[
\boxed{\pi_{ij}\pi^{ij}=0\quad\hbox{versus}\quad
                  \frac{2}{\kappa^2v^4}.}
\tag{BI.18}
\]

The squared stress norm is invariant under a simultaneous frame rotation,
so this distinction is not removed by relabelling axes. The compared scalar
curvatures are R_LC=4 Lambda+kappa mn-3 kappa^2 n^2/8 and
R_omega=4 Lambda+kappa mn; T_abc T^{abc}=-3 kappa^2 n^2/2 as in CS.
This is a counterexample to completeness of these selected scalar
readouts, not a claim that all scalar invariants coincide.

## BI-6. Future isotropization and all five linear shear modes

For every fixed t*>0 in BI.12,

\[
v(t)\sim\tfrac12(d+b)e^{\omega t},\qquad
h(t)\longrightarrow\frac\omega3,\qquad
\|\sigma(t)\|_F=O(e^{-\omega t}),\qquad
\|\pi(t)\|_F=O(e^{-2\omega t}).
\tag{BI.19}
\]

The exact shear norm identity proves these rates for arbitrary finite
S_*, not only small anisotropies in the solved rest family. Also

\[
u_\infty=\frac1q\log\!left[
\frac{(d+b)z_*-(d-b)}{(d+b)(z_*-1)}\right]<\infty.
\tag{BI.20}
\]

The normalized spatial coframe E/v^(1/3) tends to an invertible constant
matrix, and g_spatial/v^(2/3) tends to a positive-definite constant matrix.
The expansion becomes isotropic while a constant shape remains. Locally a
constant coordinate change can absorb that shape; on a fixed torus it may
remain a modulus. None of these estimates crosses the singular endpoint
t=0 or proves nonlinear stability for arbitrary nonrest matter data.

Now linearize the full homogeneous system about the isotropic CS rest
solution. Since the background shear is zero, the five shear components
obey

\[
\delta\dot\sigma+3h\delta\sigma
 =\frac1v[\delta\sigma,\Omega_b],\qquad
\boxed{\delta\sigma(t)=\frac1v
 e^{-\Omega_bu}\delta S_*e^{\Omega_bu}.}
\tag{BI.21}
\]

A perturbation of the spin generator multiplies zero background shear
and hence does not enter this equation. Conversely shear has no linear
effect on the scalar constraint, because sigma^2 is quadratic, or on the
homogeneous matter equation, because that equation contains only tr(K).
The eight matter components and scalar metric modes are exactly those
already bounded in PS-2–PS-3.

Writing E=a(I+epsilon B)E_0 locally about the isotropic background, the
trace-free metric perturbation obeys dot(B)_TF=delta sigma in the parallel
frame. Hence

\[
\|B_{\rm TF}(t)-B_{\rm TF}(t_*)\|_F
 \leq\|\delta S_*\|_F\int_{t_*}^t\frac{d\tau}{v(\tau)}
 \leq\|\delta S_*\|_F u_\infty.
\tag{BI.22}
\]

This completes bounded **homogeneous Bianchi I linear perturbations** of
the CS rest solution, including arbitrary symmetric shear and constant
shape modes, in the stated lapse/shift/frame gauge. It extends the earlier
FLRW result. Spatially varying gravitational modes and their constraints
are still outside this proof. The nonlinear result in BI.19 is specifically
isotropization within the exactly solved rest family.

## BI-7. A dimensionless phase relation and a concrete anisotropic benchmark

BI.6 in the rest sector fixes the angular rate of O about the spin axis:

\[
\|\Omega\|_{\rm op}=\frac{\kappa ZN}{8}
                     =\frac{\kappa n_0}{4}.
\]

Let phi=(kappa n_0/4)u denote its nonnegative, unwrapped rotation angle
from t*. Subtract the free mass phase m(t-t*) from BI.14. Then

\[
\boxed{\Delta\theta_{\rm contact}=\frac32\,\phi.}
\tag{BI.23}
\]

The ratio follows from the fixed cubic contact and spin-stress coefficients
of the same action. It is not an independently fitted number. It relates
a classical native phase to the shear rotation parameter; it is not a
particle-spin assignment or a prediction of a measured fundamental
constant. The comparison uses the Levi-Civita parallel frame fixed in BI-1;
it is not asserted to be an additional frame-independent observable.
If S_* commutes with Omega, O rotates within the symmetry of S_*
and does not move its shear axes. For distinct shear eigenvalues, as in
the following example, BI.23 tracks their actual rotation in the parallel
frame. Its direction is set by the signed spin generator.

Choose kappa=1, Z=4/3, m=1/2, Lambda=1/3, n_0=4/3, the Xi and Omega of
BI-5, and S_*=diag(1,0,-1). Then D^2=1, q=2, d=1, b=2, omega=1, and

\[
\boxed{v(t)=\cosh t-1+2\sinh t
 =\frac{(e^t-1)(3e^t+1)}{2e^t}.}
\tag{BI.24}
\]

Starting at t*=log(2), take any E_* with determinant 7/4. With z=exp(t),

\[
u(t)=\tfrac12\log\frac{7(z-1)}{3z+1},\qquad
u_\infty=\tfrac12\log(7/3),\qquad
\phi_\infty=\tfrac16\log(7/3),\qquad
\Delta\theta_{{\rm contact},\infty}=\tfrac14\log(7/3).
\tag{BI.25}
\]

The remaining axis rotation from this starting time is about 0.141216 radians;
this decimal is only a presentation of the exact logarithm. The shear and
anisotropic-stress squares are exactly 1/v^2 and 8/(9v^4), respectively.

| z=exp(t) | Volume v | Mean expansion h | Shear square sigma^2 |
|---|---|---|---|
| 2 | 7/4 | 13/21 | 16/49 |
| 3 | 10/3 | 7/15 | 9/100 |
| 5 | 32/5 | 19/48 | 25/1024 |
| 9 | 112/9 | 61/168 | 81/12544 |

These supplied dimensionless parameters illustrate the exact solution.
They are not an observational cosmological fit. Setting S_*=0 at the same
couplings instead gives the previous v=exp(t)-1 benchmark.

## BI-8. Exact controls, attribution and the next open problem

The module [anisotropic_cosmology.py](anisotropic_cosmology.py) consumes
the unchanged SM, CS and PS modules. It verifies:

- six conserved-spin commutators, all six symmetric expansion directions,
  and 54 complete quadratic-form identities for the anisotropic stress;
- three general spatial coframe jets, including a nonrest matter jet,
  against 48 original coframe, 72 connection, 24 matter and 48 effective
  stress components;
- a negative control: omitting spin torque leaves matter and Cartan
  equations solved but produces nine nonzero coframe residuals;
- 16 exact volume/constraint/clock checks, the CS zero-shear limit and
  the 3/2 phase coefficient relation;
- two noncommuting coframe constructions through sixth-order exact series
  identities, all 25 pairings of the five shear directions, and a failed
  coframe construction with the rotation omitted;
- the diagonal alignment condition, the distinct stress norms for equal
  selected scalar data, and the rational benchmark table.

The universal proofs are the matrix and differential arguments in this
note. Finite jet/series controls do not establish those universal claims
by sampling or replace independent review. Source hashes and parent Git
identities are retained in [SOURCE_PINS.json](SOURCE_PINS.json). All twelve
R7 numerical result groups remain unchanged in the R8 certificate.

Primary precedents are M. Henneaux, *Bianchi type-I cosmologies and spinor
fields*, Physical Review D 21 (1980), 857–863,
[DOI](https://doi.org/10.1103/PhysRevD.21.857); R. T. Jantzen, *Spinor sources
in cosmology*, Journal of Mathematical Physics 23 (1982), 1137–1146,
[author-hosted paper](https://homepage.villanova.edu/robert.jantzen/research/articles/JMathPhys23_1137-1982-spinorsources.pdf);
and B. Saha, *Nonlinear Spinor Fields in Bianchi type-I spacetime: Problems
and Possibilities*, [arXiv:1409.4993](https://arxiv.org/abs/1409.4993).
Jantzen explicitly discusses rotational frame freedom and credits the
general classical Bianchi I solution to Henneaux. Saha emphasizes the
nondiagonal spin stress and the restrictions imposed by a diagonal metric.
The present contribution to this package is the explicit reduction,
coefficient relations, full coframe solution and controls for the particular
source-bound SM contact action, together with its CS/PS limits. General
spinor cosmology, shear precession and future isotropization are not claimed
as new subjects.

The remaining perturbation task is spatially varying Einstein–matter:
coframe and connection perturbations, gauge fixing, lapse/shift constraints,
and a coupled evolution estimate. PS's spatial matter result still holds
only on prescribed geometry. The current development supplies no general
nonlinear Einstein–matter stability theorem, quantum statistics, primitive
selection of the action, or measured values of G, c, Lambda, electromagnetic
alpha, hbar or particle masses.
