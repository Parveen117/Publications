# A common noisy parent for LMD geometry and its actual density-pulse dynamics

Author: Monty. 30 October 2026.
Frozen parent: `349175122103717730fce423ce7728695f566d9f`.

The previous mass-chart measurements retained unresolved configuration modes
sharply. Their complete records could not be jointly refined when those modes'
conjugates were added. Here a different, explicitly declared measurement puts
all modes into one noisy parent from the outset. Its finite marginals agree
exactly. A Wick-calibrated outcome measure reproduces the existing LMD density
operator as a first moment. Smoothing that measure and the conjugate outcomes
gives a fixed-h continuum coarse geometry with compatible cutoff readings.

The countable Gaussian isometry is the established CMG/coherent-measurement
method, now on the actual LMD covariance. It is not a new general positivity
or coherent-state theorem. The density calibration, geometry map, pulse audit
and limits below are derived for LMD rather than transferred from MPQ.

Use 0<h<1/2 throughout the density/geometry results. The extra measurement
noise doubles the logarithmic covariance, and this restriction keeps its
measure in the already proved L2 range. No impossibility is claimed outside
that range. The parent itself exists for all h>0. All results are written
arguments conditional on the frozen inputs; no new numerical computation,
certificate replay, NRQC run, hardware execution or CI dispatch is reported.

## 1. One parent on the existing LMD Hilbert space

Use the orthonormal real basis e0=1, e_n^c=sqrt(2)cos(2 pi n z),
e_n^s=sqrt(2)sin(2 pi n z), with

    lambda_0=1, lambda_(n,c/s)=1/(2n),
    H_h=L2(mu_h),  mu_h=N(0,h C).

Define an ordinary countable product probability kappa_h on labels (u,v):

    u_j ~ N(0,2h lambda_j),  v_j ~ N(0,h/(2lambda_j)),
    all these variables independent,
    zeta_j=u_j/(2 sqrt(h lambda_j))+i sqrt(lambda_j/h) v_j.          (1)

For a normalized finite-multi-index Hermite basis vector H_a of H_h, put

    W H_a=product_j zeta_j^(a_j)/sqrt(a_j!).                         (2)

The circular Gaussian moments E[zeta^a conjugate(zeta)^b]=delta_ab a!
make W an isometry into L2(kappa_h). Let B_h=Ran W and P_B=WW^*. Then

    E(A)=W^* 1_A(u,v) W                                            (3)

is a normalized countably additive POVM on H_h. Every normal-state law is
absolutely continuous with respect to kappa_h. There is no assertion of a
normalizable coherent vector at each typical infinite-dimensional label.

For any finite set of modes, integrating the unused complex Gaussians in
matrix elements of (3) gives their Hermite Kronecker deltas. The marginal
is precisely the corresponding finite coherent POVM tensor the identity
on the unresolved state modes. For one mode its normalized vectors are

    chi_(a,b)(s)=exp[a s/(2h lambda)-a^2/(4h lambda)] exp[-i b s/h],
    dgamma=da db/(2 pi h).

Indeed <chi,H_n>/<chi,1>=zeta^n/sqrt(n!), and |<chi,1>|^2 dgamma is
the one-mode version of (1). Thus these are actual marginals of a single
measurement, not independently normalized cutoff prescriptions.

The u marginal is the sharp configuration distribution convolved with an
independent N(0,h C) field. In particular, for a pure state it has the coupling

    s ~ |psi(s)|^2 mu_h(ds),   xi ~ mu_h independently,   u=s+xi.    (4)

This is a statement about the configuration marginal only, not a joint
measurement of sharp configuration and sharp conjugate fields. The v_j
marginal is the spectral D_j measurement with added independent Gaussian
variance h/(4lambda_j). Here D_j=i h partial_(s_j)-i s_j/(2lambda_j).
These marginal identities follow first on finite modes and then by the
uniqueness of countable cylinder laws.

The old radial and mass-chart POVMs remain unchanged. This parent does not
claim to reproduce their incompatible complete records as marginals.

## 2. Wick calibration preserves the source density first moment

Let K_N(z,z)=1+sum_(n=1)^N 1/n and u_N=P_N u. Under the reference law,
Var u_N(z)=2h K_N(z,z). Define the outcome measures

    Mhat_N(u,dz)=exp[u_N(z)-h K_N(z,z)] dz.                         (5)

Apply the frozen LMD density proof with parameter 2h. For 0<h<1/2,
Mhat_N converges weakly almost surely and in L2 for bounded smearings to
a finite positive nonatomic full-support measure Mhat. Its reference moments
are

    E_kappa Mhat(w)=integral w,
    E_kappa[Mhat(w)Mhat(v)]
                       =integral integral w(z)v(y) exp[2h K(z-y)] dz dy. (6)

The same almost-sure existence holds in every normal-state outcome law by
absolute continuity. Reference moment bounds do not automatically transfer
to every such state.

The subtraction in (5) is essential. In the marginal coupling (4), condition
on s and integrate xi. At finite N,

    E_xi exp[(s+xi)_N-h K_N]=exp[s_N-h K_N/2].                    (7)

The right side is exactly the source Wick density, not exp(s_N). For bounded
w, L2 convergence of Mhat_N(s+xi,w) on mu_h tensor mu_h and contraction of
conditional expectation give

    E_xi Mhat(s+xi,w)=M_h(s,w)                                   (8)

almost surely in s. Positivity and Tonelli extend (8), for w>=0, to arbitrary
normal-state configuration laws, allowing infinite first moments. Consequently
the positive first-moment form of this outcome is exactly the source form:

    integral Mhat(w) |W psi|^2 dkappa_h
                          =integral M_h(s,w)|psi(s)|^2 dmu_h(s).  (9)

Its associated self-adjoint mean operator is J_h(w), on the existing positive
LMD sector. Equation (9) is equality of positive forms with their domains.
It is not equality of the noisy density POVM and the sharp spectral measure
of J_h(w), and it does not assert simultaneous sharp density/conjugate data.
For signed w use its positive and negative parts where the corresponding
moments are finite.

If only the physical subtraction h K_N/2 were applied to u_N, conditional
expectation would give exp(s_N), missing the source Wick normalization.
Calibrating the measurement noise is thus required, not optional notation.

## 3. Compatible continuum coarse geometry at fixed h

Choose a material resolution ell>0 and the positive circle heat kernel k_ell
with Fourier multiplier t_n=exp(-ell^2 n^2), t_0=1. Write R_ell for its
convolution operator. Define

    j_ell(z)=integral k_ell(z-y) Mhat(u,dy),
    d_ell(z)=(R_ell v)(z),      y_ell=d_ell/j_ell,
    theta_ell=Theta(y_ell),     a(theta)=1+(1/2)cos theta,
    Y(theta)=integral_0^theta dt/a(t),   Theta=Y^(-1),
    m_ell=j_ell/a(theta_ell),   P_ell=integral m_ell,
    e_ell(z)=exp(theta_ell(z)) m_ell(z)/P_ell,
    L_ell=integral e_ell.                                      (10)

The overall density scale and the global P term are retained. This is the
original smooth inverse chart applied after the declared coarse measurement.
It smooths the Wick measure itself; it is not exp(R_ell u).
The heat operator is a spatial detector response with nonlocal support, not
time evolution or a construction of a local physical apparatus.

For each outcome in a kappa-full set, j_ell and d_ell are C-infinity,
j_ell>=min(k_ell) Mhat(S1)>0, and all quantities in (10) are finite smooth
profiles with L_ell>0. The v field belongs to H^(-1-epsilon) almost surely
for every epsilon>0, because Var(v_n)=h n on each real nonconstant mode.
Heat smoothing maps it into every Sobolev space. The finite positive measure
Mhat can likewise be convolved with all derivatives of k_ell. These facts
transfer to every normal-state law.

Use the finite-label reconstructions

    j_(ell,N)=k_ell*Mhat_N,    d_(ell,N)=R_ell P_N v,               (11)

and the same inverse (10). Weak convergence of Mhat_N plus bounded total
mass makes j_(ell,N)->j_ell in each C^r norm: apply weak convergence uniformly
to the compact families of continuous tests partial_z^r k_ell(z-.) and use
their uniform continuity. Gaussian series estimates give d_(ell,N)->d_ell
in every H^r almost surely. Positivity of min j_ell and continuity of the
inverse imply convergence of the reconstructed profiles in every finite
C^r norm, on a common probability-one set.

All cutoff reconstructions, and reconstructions at any finite collection of
ell values, are functions of the same (u,v) outcome. Thus they have a common
joint POVM. For every bounded continuous function f of the smooth output,

    A_(ell,N)(f)=W^* f(G_(ell,N)) W -> W^* f(G_ell) W
                                             strongly.           (12)

Here G denotes the profiles in (10). Bounded dominated convergence proves
(12); it yields convergence of every normal-state bounded reading. Exact
compatibility refers to the common parent and its finite label marginals.
It does not mean that simply deleting Fourier coefficients of a nonlinear
coframe recovers its coarser reconstruction.

This closes a fixed-h continuum *coarse* measurement construction in the
stated range and calibration. The spatial resolution ell remains positive;
removing the Fourier cutoff is not removing ell. No pointwise unsmoothed
quantum metric or spacetime curvature is inferred.

The coarse density calibration remains exact as a first moment:

    mean-operator of integral f(z)j_ell(z) dz = J_h(R_ell f)       (13)

for nonnegative f, by (9), symmetry of k_ell and Tonelli.

## 4. Actual source pulses on the common parent

Keep the existing physical pulse U(t)=exp[-i t J_h(w)/h], for real smooth w.
It preserves the configuration probability and hence the entire u marginal
of (3), while it can change v statistics and u-v correlations. There is no
assumption that it is a deterministic classical shift of the outcome labels.

There is an explicit finite-cutoff representation of its generator. Define

    J_(h,N)(w)=integral w(z) exp[s_N(z)-h K_N(z,z)/2] dz,
    g_(N,z),j=sqrt(h lambda_j)e_j(z) for j through N, zero otherwise.

In B_h the creation and annihilation operations are multiplication by zeta_j
and differentiation partial_(zeta_j). On holomorphic cylinder polynomials,
normal ordering gives

    (mathcal J_N(w)F)(zeta)
       =integral w(z) exp[g_(N,z).zeta]
                                    F(zeta+g_(N,z)) dz,
    mathcal J_N(w)=W J_(h,N)(w) W^* on B_h.                       (14)

Indeed S_j becomes sqrt(h lambda_j)(zeta_j+partial_j), and
exp[g.(zeta+partial)-||g||^2/2]=exp(g.zeta)exp(g.partial).
Self-adjointness and the operator domain are inherited from the physical
multiplication operator J_(h,N), not guessed from the displayed expression.
No infinite pointwise shift vector g_z is asserted: its squared norm would
diverge with K_N(z,z).

The integrated source measure supplies the limit. Almost-sure convergence
J_(h,N)(w)->J_h(w) implies strong convergence of the corresponding pulse
unitaries, uniformly for t in a compact interval, by dominated convergence
of their multiplication phases. Thus the exact parent amplitude is

    Psi_t=W U(t)psi
          =lim_(N->infinity) exp[-i t mathcal J_N(w)/h] W psi.    (15)

Its outcome law is |Psi_t|^2 dkappa_h. Equation (14) is an analytic generator
formula and (15) a strong-limit definition, not an executed solver or NRQC
implementation. There is no formal infinite-dimensional Lebesgue integral.

Cutoff dynamics and coarse readouts converge together. With U_N generated
by J_(h,N), for each bounded continuous f,

    U_N(t)^* A_(ell,N)(f) U_N(t)
                           -> U(t)^* A_ell(f) U(t) strongly,      (16)

uniformly on compact t intervals on each fixed input vector. Combine the
uniform strong unitary convergence with (12); uniformly bounded strongly
convergent operators converge uniformly on compact sets of vectors. At each
fixed t all resolutions are still marginals of U(t)^* E U(t). A common
measurement of quantum trajectories at different times is not asserted.

For w>=0, (9) also identifies mathcal J(w) with the positive form compression
of multiplication by Mhat(w) to B_h. Compression is important: evolving by
the scalar multiplier Mhat(w) on the full outcome L2 space would only change
an outcome phase and would leave all probabilities unchanged. It is not (15).

## 5. An exact variance test rejects a tempting classical outcome update

Take the physical vacuum psi=1 and a uniform density pulse with w=1. Define

    I(gamma)=integral integral exp[gamma K(z-y)] dz dy,
    J0=M_h(S1),  D0=D(1).

The measured v0 has the D0 spectral law convolved with N(0,h/4). Direct
differentiation in the physical zero mode gives

    D0 U(t)1=U(t)[t J0-i s0/2].

Both terms are in L2; their squared-norm cross term has zero real part.
The derivative is a Stone-domain statement: constant Cameron--Martin shifts
multiply J0 by an exponential, and their phase difference quotients are
locally bounded by a constant times J0. Its finite L2 vacuum moment justifies
the derivative; no unbounded operator product is being assumed.
Using E J0=1 and E J0^2=I(h), the actual outcome moments are

    E_t v0=t,
    Var_t(v0)=h/2+t^2[I(h)-1].                                (17)

No product-domain assumption for D0^2 is needed; the second moment uses
the norm of D0 U(t)1. In contrast, the proposed update of the vacuum labels

    u unchanged,   v0 -> v0+t Mhat(u,S1)                        (18)

would give the same mean t but variance h/2+t^2[I(2h)-1], since u and v
are independent in the reference law. Both integrals are finite for h<1/2,
and I(2h)>I(h): K(r)>=1-log 2>0 away from the diagonal. Thus (18) is wrong
for every nonzero t. Matching density first moments is insufficient to
identify the quantum pulse law. This rejects this particular additive-label
rule, not every possible stochastic representation of the dynamics.

## 6. Length moments and a declared optional detection threshold

Each L_ell is finite almost surely, but its ungated vacuum mean is infinite
for every h>0 in the stated range and every ell>0. To see this, decompose

    Mhat=e^(u0-h) Mhat_perp,  d_ell=v0+d_ell,perp.

For almost every perpendicular outcome, k_ell*Mhat_perp is bounded above
by a finite positive constant and d_ell,perp is bounded. The independent
Gaussian v0 has positive probability of exceeding that latter bound by 1.
On this event min y_ell>=c exp(-u0), for some c>0 independent of u0.
The frozen inverse satisfies Theta(y)=sqrt(3)y/2+b(y), with bounded b.
Since L is a positive weighted average of exp(theta), it is bounded below
by exp[c' exp(-u0)-||b||_infinity]. Integrating u0~N(0,2h) diverges.
Tonelli proves the assertion without interchanging conditionally convergent
integrals. The positive zero-mode Wick factor is included in c'.

If a finite-mean length operator is required, reuse CMG's explicit additional
contract: choose delta>0, accept only if min j_ell>=delta, and record all
other runs as a rejection symbol. The retained number is

    L_acc=1_(min j_ell>=delta) L_ell,
    L_acc<=exp[||b||_infinity+(sqrt(3)/2)||d_ell||_infinity/delta]. (19)

The flag distinguishes rejection from a zero length. The Gaussian d_ell has
trace-class covariance in H^3 and exponential moments of its norm. Finite
Hermite images W H_a are polynomials in finitely many Gaussian labels.
Holder's inequality proves all finite moments of L_acc on their dense span.
Consequently the positive form integral L_acc |W psi|^2 dkappa_h is closed
and densely defined; it specifies a positive self-adjoint mean operator.
Its spectral PVM is not the original POVM. This is the same established
closed-form method as CMG, with the LMD/Wick-density hypotheses checked here.

The threshold boundary has reference probability zero: min j_ell equals
e^(u0-h) times an independent positive minimum. Therefore the cutoff
rejection-inclusive effects also have the strong bounded-test limit (12).
The threshold and ell are detector calibrations, not derived physical laws.
Invariance of this length-form domain under U(t) is not established.

## 7. Classical profiles, mean length and the pulse comparison

In the LMD coherent state centered at smooth s_*,d_*, the parent labels are
independent Gaussian fields with means s_*,d_* and covariances 2h C and
(h/2)C^(-1). Thus Mhat has law e^(s_*) times the centered 2h Wick measure.
The density and conjugate means after smoothing are

    E j_ell=R_ell(exp s_*),   E d_ell=R_ell d_*.                  (20)

Notice R_ell(exp s_*), not exp(R_ell s_*). Write j_*=exp s_*>0. For
0<h<=h0<1/2, every real Fourier coefficient of Mhat-j_*dz has variance
bounded by C h independently of its index: use the bound on the Fourier
basis and the integrable kernel inequality
|exp(2hK)-1|<=2h K exp(2h0 K). Heat summation then gives

    E||j_ell-R_ell j_*||_H3^2 <= C h ell^(-7),
    E||d_ell-R_ell d_*||_H3^2 <= C' h ell^(-8),     0<ell<=1.    (21)

For the second estimate the exact reference trace is
(h/2)[1+4 sum_(n>=1) n(1+n^2)^3 exp(-2ell^2 n^2)]. Comparing the
Gaussian-decaying sums with their integrals proves both powers. Constants
are independent of h,ell in these ranges; no optimality is asserted.

Choose ell(h)->0 and h ell(h)^(-8)->0. The means in (20) approach the
smooth source centers, and (21) gives H^3 convergence in probability of
(j_ell,d_ell) to (j_*,d_*). Sobolev embedding, positivity of min j_* and
continuity of the full inverse (10) recover the original coframe and length.
All bounded continuous geometry readings converge. For a fixed
0<delta<min j_*, acceptance tends to one. Moreover every finite moment of
L_acc converges to the corresponding classical length moment: the centered
Gaussian d_ell has uniformly bounded H^3 covariance trace, so (19) gives
uniform integrability of every power. These moment claims concern this
unpulsed coherent family.

There is also a controlled actual-pulse comparison for bounded geometry.
Let B_t be a prescribed smooth integrated pulse on a compact t interval,
and psi_h(t)=exp[-i J_h(B_t)/h]psi_h. The u law is unchanged. On these
coherent-state orbits the affine derivative satisfies

    U_h(t)^* D(e_j) U_h(t)=D(e_j)+J_h(B_t e_j)

in its action on the initial coherent vector. The centered D(e_j) action
is a purely imaginary Gaussian multiplier of that vector; the centered J
action is real. Their squared-norm cross term therefore vanishes. The v_j
marginal and LMD second moments give exactly

    E_t v_j=(d_*+j_* B_t)_j,
    Var_t(v_j)=h/(2lambda_j)+Var_(psi_h)[J_h(B_t e_j)].           (22)

The domain assertion follows by differentiating the Cameron--Martin group:
M_h(s+r e_j,dz)=exp(r e_j(z))M_h(s,dz), and e_j is bounded. Its differentiated
pulse term is controlled in L2 by M_h(1) on the smooth coherent center.
This proves the action used in (22) without imposing a common domain for
arbitrary products of GMC operators.

The last variance is at most C_T h uniformly in j,t, by the same bounded-window
kernel estimate. Smoothing and summing add only O_T(h ell^(-7)) to the second
bound in (21). Hence the same scale condition recovers j_*, d_*+j_* B_t,
and y_*+B_t in probability, uniformly in single-time error probabilities.
This is the exact source pulse orbit, not the incorrect update (18).

After a nonlinear pulse the v law is generally not Gaussian. We therefore
do not extend the preceding accepted-length moment convergence or the dense
form-domain invariance claim to these pulse-evolved states. Bounded geometry
recovery and finite mean-length recovery have distinct proved hypotheses.
Nor is uniform single-time concentration a joint quantum trajectory law.

## 8. Material-flow audit and the remaining gate

The fixed heat kernel and Gaussian measurement are calibrated in the chosen
material circle coordinate. No full material-diffeomorphism covariance of
this same calibrated parent is claimed. There is a concrete normalization
test that any proposed label action must pass. The physical LMD flow uses
a_h=1+h/2 in s -> s composed g+a_h log g'. But the outcome measure has
parameter 2h, so applying that same affine formula to u gives

    Mhat(u composed g+a_h log g',dz)
                     =(g')^(-h/2) g^* Mhat(u,dz),               (23)

by the frozen LMD density law with parameter 2h. The coefficient that would
make this *label measure* transform as a pure density is 1+h instead.
Equation (23) audits a proposed classical label transformation; it is not
an identification of the actual quantum action of T_g on this POVM. A
covariant calibration family must handle the measurement noise as well as
the physical mean field. Heat smoothing also fails to commute with a general
non-rigid change of material coordinate. Neither issue can be dropped.

Constructed within the stated assumptions: one countable noisy LMD parent
with exact finite marginals; density first-moment agreement with the original
J_h; compatible fixed-h continuum coarse geometry; source pulse dynamics
with a convergent cutoff realization; a variance test excluding the naive
classical outcome shift; and controlled classical geometry recovery, including
accepted mean length for the stated unpulsed coherent family.

Remaining: all-material-flow covariance of the measurement/calibration family;
fixed-h removal of ell; accepted-length domains after nonlinear pulses; a
local physical apparatus; autonomous gravitational dynamics, spacetime
curvature and physical 3+1 gravity. The common-parent obstruction is addressed
by a different noisy measurement, not by claiming the incompatible older
records suddenly agree. No quantum advantage or NRQC implementation follows.
