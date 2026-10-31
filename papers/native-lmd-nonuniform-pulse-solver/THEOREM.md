# A phase-retaining solver for one nonuniform LMD source

Author: Monty. 31 October 2026.
Frozen research parent: `ca080861d646774ad94a8e48be668194a8086d59`.

This packet implements a first nonuniform-source marginal, conditional on
the existing LMD/common-parent construction. A relative spatial angle is
retained in the conditional source and integrated in its characteristic
function. Higher spatial modes have an explicit error budget. The uniform
full-mass theorem is not imported as a weighted mass law.

## 1. Fixed source and actual measured target

For 0<h<1/2 and |e|<1, choose a positive prescribed profile

    w_(e,beta)(z)=1+e cos(2 pi z-beta),
    J_e=integral w_(e,beta)(z) M_h(dz),
    U_(t,e)=exp(-i t J_e/h).

Use the same vacuum and the same total conjugate-label marginal as in the
frozen pulse-characteristic packet. Constant translation multiplies J_e by
exp(a). Its bounded matrix-element proof therefore gives

    Phi_(t,e)(k)=d E exp(i a J_e),
    d=exp(-h k^2/4),  a=(2t/h)sinh(hk/2),
    p_cos,+=(1+Re Phi)/2,  p_sin,+=(1+Im Phi)/2.        (1)

The pulse is spatially nonuniform but the detector here is still the total
label. A common rotation of a single profile has the same vacuum marginal;
this does not identify local or joint detector laws. The source profile and
the existing detector remain specified inputs, not autonomous field equations.

The preceding next-step note derived

    E J_e=1,
    I_e(h)=E J_e^2=I(h)[1+e^2 h/(2(2-h))],
    I(h)=exp(h) Gamma(1-h)/Gamma(1-h/2)^2.              (2)

For completeness, the two-point kernel is
L_h(r)=exp(h)(2|sin(pi r)|)^(-h). Integrating the product of the two profiles
over their common origin gives 1+(e^2/2)cos(2 pi r). Integration by parts
in the sine-power integral gives its first Fourier coefficient
I(h)h/(2-h). Bounded weights and h<1 justify these integrals. In particular,
at e=1/2,h=3/8 the second moment is (107/104)I(3/8), not I(3/8).

Two derivatives of (1), justified by this L2 moment, give

    E v_0=t,
    Var(v_0)=h/2+t^2[I_e(h)-1].                        (3)

Thus changing this profile at fixed integrated strength leaves the mean
unchanged and changes the variance by t^2 e^2 I(h)h/[2(2-h)]. First-order
mean response alone misses the spatial modulation in this marginal.

## 2. Conditional source with its relative phase

Retain the first pair of raw Fourier coefficients X,Y, each of variance h,
and q of variance h, independent. Write X=R cos(phi), Y=R sin(phi),
S=R^2/(2h). Then S~Exp(1), phi is uniform, and q,S,phi are independent.
Let B_n denote the standard modified Bessel I_n, to avoid confusion with I(h).
The first-mode Wick martingale is

    J_(1,e)=exp(q-h)[B_0(R)+e cos(phi-beta) B_1(R)]
           =Q(q)[W_0(S)+e cos(phi-beta)W_1(S)],
    Q=exp(q-h/2),
    W_n(s)=exp(-h/2) B_n(sqrt(2hs)).                   (4)

Positive circle integration yields (4), and the weighted Wick-martingale
identity gives E[J_e|q,S,phi]=J_(1,e). Since 0<=B_1(r)<B_0(r) for r>=0
(strict except the zero value of B_1), this conditional source is positive.
Deleting the cosine term would change the source compression: it preserves
the first moment but loses a positive second-moment contribution.

The first-mode two-point kernel is exp(h[1+cos(theta-theta')]). Circle
integration of the profile product then gives

    E J_(1,e)^2=exp(h)[B_0(h)+(e^2/2)B_1(h)].          (5)

Equivalently E W_0^2=B_0(h), E W_1^2=B_1(h); the angular average kills the
cross term and averages cos^2 to 1/2. The omitted-mode variance is exactly
I_e minus (5), by conditional-expectation orthogonality.

## 3. Positive radial cell moments, including infinity

Retain the predecessor Gaussian probabilities p_i and tilted first masses
a_i for q bins, including both unbounded tails. For radial cells (u_j,v_j],
including a final unbounded tail, define

    p_j=exp(-u_j)-exp(-v_j),
    b_(n,j)=integral_(u_j)^(v_j) exp(-s) W_n(s) ds,
    S_n=sum_j b_(n,j)^2/p_j,   S_q=sum_i a_i^2/p_i.

The predecessor supplies b_(0,j) by a positive integer-Gamma mixture.
For b_(1,j), set lambda=h/2. The B_1 power series gives a new cell formula

    b_(1,j)=sum_(n>=0) c_n [Q_(n+3/2)(u_j)-Q_(n+3/2)(v_j)],
    c_n=exp(-lambda)lambda^(n+1/2)Gamma(n+3/2)/(n!(n+1)!),
    c_(n+1)/c_n=lambda(n+3/2)/[(n+1)(n+2)],            (6)

where Q_s(x) is the normalized upper Gamma integral, Q_s(0)=1 and
Q_s(infinity)=0. All coefficients and bin probabilities are positive.
The coefficient ratio decreases with n. After n=0,...,N the omitted mass
is at most c_N r_N/(1-r_(N+1)), where r_n is the ratio in (6).

The native implementation uses N=32, computes Q_(1/2) using the predecessor's
certified lower Gamma integral, and increases shape via

    Q_(s+1)(x)=Q_s(x)+exp(-x)x^s/Gamma(s+1).            (7)

The increment is updated by multiplication by x/(s+1). It is not computed
by subtracting two large unnormalized Gamma values. Known bounds [0,1] and
0<=b_1<=b_0 may be intersected with the directed enclosures. No probability
is dropped and no retained mass is renormalized.

## 4. Continuous angle, finite source cells, complete error

Condition on the q and S bins but keep the relative angle theta=phi-beta.
Define P_ij=p_i p_j and m_(n,ij)=(a_i/p_i)(b_(n,j)/p_j). The conditional
source is m_0+e cos(theta)m_1. Its exact captured second moment is

    S_G=S_q[S_0+(e^2/2)S_1],
    delta_G^2=I_e-S_G
      =[I_e-exp(h)(B_0(h)+(e^2/2)B_1(h))]
       +[exp(h)-S_q][B_0(h)+(e^2/2)B_1(h)]
       +S_q[(B_0(h)-S_0)+(e^2/2)(B_1(h)-S_1)].        (8)

These are omitted-mode, Gaussian-bin and radial-bin contributions, each
nonnegative. Tail-cell quantization is included. The earlier conditional
Taylor proof applies to this sigma algebra, even though theta is continuous:

    |E exp(i a J_e)-E exp(i a E[J_e|G])|
          <= min(2,a^2 delta_G^2/2).                   (9)

Let C(x) denote ordinary Bessel J_0. Its circle integral performs the retained
angular expectation exactly:

    chi_G(a)=sum_(i,j) P_ij exp(i a m_(0,ij)) C(a e m_(1,ij)). (10)

Equation (10) includes modulation-dependent angular cancellation. Replacing
C by 1 would average the source before applying its nonlinear phase and
would produce the uniform compression, not this weighted expectation.
This correct integration explains how an angle can be removed from the
numerical grid without removing its contribution to the result.

Each base phase in (10) uses the nine pinned NRQC outputs and their inherited
certified interpolator. Since |C|<=1, its phase remainder is not amplified.
Native intervals enclose C via

    C(x)=sum_(n>=0) (-x^2/4)^n/(n!)^2.

For an argument interval with m=max|x| and last term t_N, an absolute tail
bound is |t_N| (m^2/4)/(N+1)^2 / [1-(m^2/4)/(N+2)^2]. The predeclared
thresholds m<=1,4,16,32 use N=16,24,40,64 respectively. For m>32 the certified
fallback is [-1,1], with the affected source probability recorded; no cell
is skipped. Arithmetic dependency width and all remainders propagate to
the assembled result. Then widen by (9), multiply by d and divide by two.

For implementation checks independent midpoint angular quadrature uses
the derivative bound |(d^2/dtheta^2)cos(x cos(theta))|<=|x|+x^2. Its normalized
error is at most (2 pi/M)^2(|x|+x^2)/24. These checks do not replace the
series or continuum error proofs.

## 5. Frozen scope and interpretation

The request contains four positive e=1/2 cases across three source strengths,
including h=3/8,t=1/2; a stronger e=3/4 modulation; amplitude-sign, time-reversal
and zero-pulse controls; and a uniform overlap reference. Two nested grids
use q bins 64 then 128 on [-2,2] plus tails, and radial bins 16 then 32 on
[0,16] plus the unbounded tail. Probability halfwidth remains 1/200.

No accuracy outcome is presumed. If the ideal first-mode contribution in
(8), after conversion to probability error, already exceeds 1/200, finer
q/S bins cannot remove that obstacle in this certificate. It is not a lower
bound on the actual error. All cases and both refinements are retained.

Native JS performs scientific arithmetic with the inherited rational field.
Python 3.12 only orchestrates and hashes. The protocol is committed before
new probability evaluations. Nine old NRQC phase outputs are reused; no new
calculator, hardware or CI job is dispatched. No full weighted mass law,
joint coframe distribution, autonomous gravity, physical 3+1 closure or
computational advantage is asserted.
