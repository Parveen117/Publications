# The next gate: a nonuniform source marginal

Author: Monty. 31 October 2026. Post-protocol analytic boundary; no additional
numerical campaign is represented by this note.

The full uniform total-mass law removes the need for Fourier truncation in
the present scalar expectation. It does not supply the distribution of a
weighted mass or several masses together. A concrete next target is

    w_e(z)=1+e cos(2 pi z),   0<|e|<1,
    J_e=integral w_e(z) M_h(dz),   U_(t,e)=exp(-i t J_e/h).

This is a nonnegative, nonconstant source in the existing LMD measure. Keep
the same vacuum and the same total conjugate-label detector. Constant
translation still gives J_e(q+a,eta)=exp(a)J_e(q,eta). Therefore the bounded
translation proof in the frozen pulse-characteristic packet applies to this
specific marginal without modification of its Gaussian overlap:

    Phi_(t,e)(k)=exp(-h k^2/4)
       E exp(i [(2t/h)sinh(hk/2)] J_e).

This statement concerns a total-label marginal after a nonuniform pulse;
it does not determine local geometry or the joint detector law.

## Why the uniform mass law cannot simply be reused

There is an exact second-moment witness. Put I(h)=E J_0^2 and
L_h(r)=exp(h)(2|sin(pi r)|)^(-h). Translation invariance and integration over
one circle variable give

    E J_e=1,
    E J_e^2=I(h)+(e^2/2) integral_0^1 L_h(r)cos(2 pi r)dr
           =I(h) [1+ e^2 h/(2(2-h))],       0<h<1.

For the last equality, use cos(2x)=1-2sin^2(x). The ratio of the sine-power
integrals with exponents 2-h and -h is (1-h)/(2-h), by integration by parts
on integral sin^(2-h)(x) (the boundary terms vanish for h<1). Thus the first
Fourier coefficient is I(h)h/(2-h)>0. The linear cross terms vanish because
the cosine has zero mean. All exchanges are justified by LMD's integrable
positive two-point kernel and bounded weights.

For e=1/2,h=3/8 the second moment is exactly (107/104)I(3/8).
The mean is still one, so neither unchanged uniform sampling nor scalar
rescaling can supply this law. This is an analytic mismatch, independent of
any finite solver performance. It also guards against treating the new
auxiliary Exp(1) variable as a replacement for the entire field.

## Proposed development and stopping rule

1. Freeze this nonuniform source, the detector and the required probability
   tolerance before running a new campaign.
2. Construct positive conditional masses with the phase relative to the
   external cosine profile retained, or find and audit an applicable weighted
   mass-law theorem. A common angular orientation can no longer be integrated
   out solely because the unweighted total mass was rotation invariant.
3. Prove the complete variance/error budget. Reproduce the two moments above,
   recover e=0, and check the distributional e -> -e symmetry from a half-circle
   shift. Reject the uniform-law substitution.
4. Evaluate bounded probabilities with native directed arithmetic and the
   pinned phase engine. Retain a wider enclosure if the tolerance is missed.

This is a specific next expectation problem. Joint geometry and currents,
calibration covariance, operator domains, fixed-h removal of geometric
smoothing, autonomous gravity and physical 3+1 dynamics remain open.
