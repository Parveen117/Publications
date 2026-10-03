# Coherence first and thermodynamic distinction formation

Monty Dabas · CF-1–CF-8 · 30 September 2026

## Proposed principle and mathematical scope

The motivating proposal is: coherence precedes distinction; a cut gives form
to a distinction within a whole, and its record must retain target-relevant
memory. The centre signifies dissolution of distinguishability, rather than
an externally located creator or an already measured empty state. This is
a proposed foundational interpretation, not an experimentally established law.

We construct one explicit native algebra and thermodynamic realization. It
has vanishing distinction amplitude and finite formation energy at its centre,
divergent response-volume density near that centre, radius-dependent seam
orientation, and calculable thermodynamic curvature. None of its constitutive
parameters is claimed to be uniquely selected by coherence alone.

The source order is the canonical RKF cut field and EMK algebra, then an
explicit cut/record construction, then a calibrated thermodynamic chart.
The shared engine is imported from RKF commit
`3cc5a33b05c16d59c90994ddda69dedc0d392424`; no alternate engine is introduced.
The original thermo-compass model and CID-1 covariance routines are consumed
unchanged. Exact file dependencies are in [SOURCE_PINS.json](SOURCE_PINS.json).

The literal pre-representational whole is not identified with a zero vector,
a probability distribution, or any particular state of this model. The
coordinates below describe distinction records after a chart is admitted.
No private philosophical manuscript is reproduced in this publication.

## CF-1 Native distinction amplitude and dissolution of the cut

Use the source relations

\[
R^2=-1,\quad K^2=1,\quad KR=-RK,\quad
R^\dagger=-R,\quad K^\dagger=K.
\]

Set L=RK. Native multiplication gives L^dagger=L, L^2=1 and KL=-LK.
For radial coefficients a,b define D=aK+bL. Then

\[
D^\dagger=D,\qquad D^2=(a^2+b^2)1.
\tag{1}
\]

For r=sqrt(a^2+b^2)>0, N=D/r is an involution and

\[
P_D=(1+N)/2,\qquad Q_D=(1-N)/2
\tag{2}
\]

are complementary native apertures. Equation (1) follows by cancellation of
the two mixed terms. Equation (2) then follows by direct multiplication.
Radial square roots are taken in the completed cut field; rational witnesses
use Pythagorean coefficient pairs.

There is no direction-independent continuous extension of P_D to a=b=0.
Along (a,b)=(r,0) it equals (1+K)/2; along (-r,0) it equals (1-K)/2. These
are different because the certified normal basis has K nonzero. D itself
does tend to zero. Substituting N=0 at the centre would give P=1/2, which
is not idempotent. Thus disappearance of distinction amplitude does not
select a surviving two-way cut. The centre is a boundary of this cut chart.

## CF-2 Constitutive records and the memory cost of forgetting a branch

On an admitted native module, introduce a tagged record operation

\[
\mathcal C_D(\psi)=(P_D\psi,Q_D\psi).
\tag{3}
\]

The tags distinguish the two record channels. Summing the channels recovers
psi exactly, so this operation is injective. A physical implementation that
creates such records is an additional interaction hypothesis; the formula
does not select a single outcome or assume Born probabilities.

For a later arrow A, the recognized future is

\[
P_DA\psi=P_DAP_D\psi+P_DAQ_D\psi.
\tag{4}
\]

Dropping the second channel is harmless for all states exactly when
P_D A Q_D vanishes on the admitted module. In the faithful regular module,
this is the native algebra condition P_D A Q_D=0. In particular, at
P=(1+K)/2, PRP=0 while PRQRP=-P. A first invisible return is therefore
not evidence of an absent source. This consumes the earlier NC memory law.
For P_K=(1+K)/2 and P_L=(1+L)/2, native multiplication also gives
[P_K,P_L]=-R/2. The order of tagged branch operations can therefore change
their records even when summing all branches reconstructs the original input.

## CF-3 One explicit thermodynamic energy family

Choose positive reference units and translated, dimensionless thermo
coordinates s=(S-S0)/S*, v=(V-V0)/V*. All formulas here use normalized energy
and its conjugates. Declare the chart adapter a=s, b=v in CF-1, so the native
contrast D0=sK+vL has D0^2=r^2*1. This identification is an explicit adapter
postulate, not a claim that primitive EMK coefficients already are entropy
and volume. Let r=sqrt(s^2+v^2), and declare

\[
\mathcal U(s,v)=U_0+T_0s-P_0v+W(s,v),\qquad
W=\kappa r^\nu+\frac{\epsilon}{2}s^2,
\quad \kappa>0,\ 1<\nu<2,\ \epsilon\ge0.
\tag{5}
\]

This is the constitutive postulate, not a consequence of (1) alone. W is the
proposed distinction-formation energy. It is nonnegative, vanishes only at
r=0, and is continuously differentiable at the centre with zero gradient.
It is smooth on the punctured chart, where the source thermo theorems apply.

Write h=kappa*nu*r^(nu-2), e=(s,v)/r and delta=nu-2. Then

\[
T=T_0+hs+\epsilon s,\quad P=P_0-hv,
\qquad H=D^2\mathcal U=h(I+\delta ee^T)+\epsilon
\begin{pmatrix}1&0\\0&0\end{pmatrix}.
\tag{6}
\]

The radial part has eigenvalues h(nu-1),h. Adding the positive semidefinite
epsilon term keeps H positive definite for every r>0. Positive S,V,T,P hold
in a sufficiently small punctured neighbourhood after choosing positive
S0,V0,T0,P0. The singular centre itself is not a regular thermo chart.

For H=[[A,B],[B,C]], the original compass formulas give

\[
C_V=T/A,\quad C_P=TC/\det H,\quad
K_S=VC,\quad K_T=V\det H/A,
\quad C_P/C_V=K_S/K_T.
\tag{7}
\]

Thus all four Legendre potentials, their Maxwell relations and the six
source lambda/z responses follow in their regular charts. The exponent nu,
amplitude kappa and anisotropy epsilon are not those six differently typed
lambda coordinates. No new universal scalar lambda is silently defined.

## CF-4 Divergent local response density with finite energy and finite measure

In the declared affine chart, define the response-volume density
w_H=sqrt(det H). With e=(cos phi,sin phi),

\[
\det H=h^2(\nu-1)+\epsilon h
 [1+(\nu-2)\sin^2\phi].
\tag{8}
\]

Consequently, as r decreases to zero,

\[
w_H\sim\kappa\nu\sqrt{\nu-1}\,r^{\nu-2}\longrightarrow\infty,
\qquad W\longrightarrow0.
\tag{9}
\]

The density is a coordinate density, not a coordinate-invariant scalar.
Its measure w_H ds dv is invariant under legitimate coordinate changes.
For epsilon=0 the exact disk measure and radial metric distance are

\[
\mu_H(B_R)=2\pi\kappa\sqrt{\nu-1}\,R^\nu,
\quad
\ell(0,R)=\frac{2\sqrt{\kappa\nu(\nu-1)}}{\nu}R^{\nu/2}.
\tag{10}
\]

Both are finite and tend to zero. For fixed epsilon>=0, (8) gives the same
leading local integrability. A divergent local coefficient therefore does
not establish an infinite total information content. Nor does this measure
already have the units of physical spatial volume or energy density.

At the exact rational benchmark nu=3/2, kappa=2/3, epsilon=0, r=t^4,

\[
W=\tfrac23 t^6,\quad h=t^{-2},\quad
\det H=\tfrac12t^{-4}.
\tag{11}
\]

This is a calculable realization of small formation energy with large local
response sensitivity. The designation 'information density' additionally
requires the information interface in CF-7 and its stated scope.

## CF-5 Radial seam rotation from the energy Hessian

Define the native response-anisotropy element

\[
B_H=(H_{ss}-H_{vv})K+2H_{sv}L.
\tag{12}
\]

CF-1 gives B_H^2=d_H^2*1, with
d_H^2=(H_ss-H_vv)^2+4H_sv^2. Where d_H>0, its normalized involution defines
a native aperture. It is the cut associated with the principal response
directions in this chosen normalized chart. At d_H=0 the orientation is
undetermined even though H may remain positive definite.

The energy-selected aperture of B_H is generally different from the aperture
of D0. To attach the same vanishing contrast amplitude to that response cut,
set D_H=r B_H/d_H; then D_H^2=r^2*1. Near the centre, B_H/d_H tends to
-(cos(2phi)K+sin(2phi)L), so D_H tends to zero while its aperture again has
different directional limits. CF-2 can be applied to this response cut.

If theta denotes a continuously chosen principal-axis angle on such a patch,

\[
2\theta=\arg\{\epsilon+(\nu-2)h e^{2i\phi}\},
\quad
\left.\frac{\partial\theta}{\partial r}\right|_\phi
=\frac{\epsilon(\nu-2)^2h\sin2\phi}
 {2r[\epsilon^2+2\epsilon(\nu-2)h\cos2\phi+(\nu-2)^2h^2]}.
\tag{13}
\]

The complex expression is a coordinate shorthand for the ordered pair of
coefficients in (12), not a replacement for native UGD arithmetic. Equation
(13) follows from differentiating that pair's angle and h'=(nu-2)h/r.
For epsilon>0 and sin(2phi) nonzero, the seam orientation changes with radius.
On the 45-degree ray s=v the derivative is positive. With epsilon=0 there
is no radial rotation. Thus the new law supplies an explicit load-bearing
coupling for the proposed rotation, rather than assuming every shrinking
circle must rotate its seams. The 45-degree reference uses fixed normalized
S,V axes; it is not unit-independent without their metric/unit transport.

A mere radial rotation of a frame does not create curvature. CF-6 computes
curvature of the energy-derived metric independently of this angle.

## CF-6 Curvature computed from the original thermo response connection

Consume TC-4/TC-5: in the affine Hessian chart,
A_i=(1/2)H^-1 partial_i H and
F_sv=-(1/4)[H^-1 partial_s H,H^-1 partial_v H]. Its Gaussian curvature is

\[
\boxed{\mathcal K_H=
\frac{\epsilon h^2(\nu-2)^2
 [\cos^2\phi-(\nu-1)\sin^2\phi]}
 {4r^2(\det H)^2}.}
\tag{14}
\]

For a derivation, rotate the affine axes at the evaluation point into radial
and tangential directions, without differentiating this constant rotation.
The Hessian entries become
a=h(nu-1)+epsilon cos^2 phi,
b=-epsilon sin phi cos phi,
c=h+epsilon sin^2 phi.
The independent cubic entries are
U_111=(nu-1)B, U_122=B, U_112=U_222=0, where B=h(nu-2)/r.
Substitution in
K=-det([[a,b,c],[U_111,U_112,U_122],[U_112,U_122,U_222]])/(4 det(H)^2)
gives (14); the determinant identity also follows by expanding the displayed
TC commutator. The verifier separately uses full Christoffel differentiation
and the analytic fourth derivatives, not only this simplified expression.

At epsilon=0 the punctured metric is flat despite (9). In that case (10)
puts it in the form d ell^2+gamma^2 ell^2 d phi^2, with
gamma=nu/(2 sqrt(nu-1)). Its apex is a singular metric-completion point;
flatness on r>0 makes no assertion of a smooth or flat extension at r=0.
For epsilon>0, curvature can be positive, zero or negative depending on
direction. This is thermodynamic response curvature, not yet spacetime
curvature, and it is distinct from CID-1's metric-free loop obstruction.

## CF-7 The original CID information ledger applies pointwise

For the benchmark in (11), take e=(c,d), f=(-d,c), c^2+d^2=1. Declare five
outcomes with probabilities (1/4,1/4,1/8,1/8,1/4) and two-component statistic
vectors (e/t,-e/t,2f/t,-2f/t,0). Their mean is zero and their covariance is

\[
\operatorname{Cov}(T)=\frac{ee^T}{2t^2}+\frac{ff^T}{t^2}=H.
\tag{15}
\]

The original CID-1 covariance and coarse_grain functions verify this equality.
Merging the two signs within each axis makes every conditional mean zero:
recognized covariance vanishes and discarded covariance equals H. Retaining
singleton outcomes recognizes all of H. Intermediate partitions give exact
positive recognized/discarded sums. Thus a zero recognized metric need not
mean a zero full response metric.

This is a family of pointwise information realizations with statistics
depending on the chosen point, not a proved global fixed-statistic Fisher
manifold. Its outcome probabilities are independent of t, so its Shannon
entropy is constant even while H diverges. Coarse-graining here is a passive
information operation; it must not be conflated with the constitutive record
interaction hypothesized in CF-2. CID-1's obstruction theorem remains in its
original response algebra and is not relabelled as (14) or gravitational force.

## CF-8 Energy, inward restoring response and the remaining physical bridge

Along an arbitrary differentiable chart path,

\[
\frac{d\mathcal U}{dr}=T\frac{ds}{dr}-P\frac{dv}{dr}.
\tag{16}
\]

At fixed phi, subtracting the affine reference reservoir terms in (5) gives

\[
\frac{dW}{dr}=\kappa\nu r^{\nu-1}+\epsilon r\cos^2\phi>0.
\tag{17}
\]

The negative energy derivative is an inward generalized restoring response
on this state-coordinate path. Its magnitude tends to zero at the centre,
even while the response density diverges. A relaxation law, mobility, clock,
and physical spatial interpretation are not implied by that sign. This is
not a derivation of gravitational acceleration or an inverse-square law.

For the benchmark, W decreases from 2/3 to 1/96 as r decreases from 1 to
1/16, while det(H) increases from 1/2 to 8. These numbers follow from the
declared model and are not measured constants or fitted observations.

Adding a constant to U changes its reference value but not H, seams, response
curvature or any differential in (16). Thus these data alone cannot determine
an absolute total energy or a mass through an uncalibrated E=mc^2 declaration.
The next physical obligation is to derive an interaction/clock/source law and
identify spatial and inertial measurements using the complete native records.

This development makes the proposed centre, seams and energy calculable. It
does not claim all its postulates follow from the phrase 'coherence first',
prove a universal diagram, or certify metaphysical infinity. General written
proofs, finite exact controls, source pins and independent review are distinct.
