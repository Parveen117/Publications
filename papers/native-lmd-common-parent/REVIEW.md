# Written audit and closure ledger for the common parent

Author: Monty. 30 October 2026.

This is a written analytic audit, not independent peer review. File hashes
and repository checks establish provenance only. There are no new scientific
test counts, numerical fixtures or calculator performance claims.

| Target | Status in this packet | Conditions and limits |
|---|---|---|
| One parent with exact finite adapted-coordinate marginals | Constructed | Matched Gaussian noise on every mode; LMD Hilbert space |
| Retaining the old mass-chart records exactly | Not claimed | Those incompatible sharp-residual records are preserved separately |
| Source density first moment | Exact positive-form equality | Correct outcome Wick subtraction h K_N; 0<h<1/2 |
| Compatible fixed-h continuum coarse geometry | Constructed | Fourier cutoff removed; detector resolution ell remains positive |
| Actual density-pulse evolution | Strong cutoff realization | Physical J_h(w) retained; analytic generator, no executed solver |
| Classical additive update using measured density | Rejected by variance | The specific rule (18); no universal stochastic no-go claimed |
| Ungated vacuum mean length | Infinite | Even at fixed positive ell |
| Accepted-length mean operator | Closed positive form, dense domain | Added threshold delta and retained rejection flag |
| Classical geometry and prescribed pulse orbit | Recovered in probability | ell->0 and h ell^(-8)->0; single-time statements |
| Accepted-length moment recovery | Proved for unpulsed coherent family | Not extended to nonlinear pulse-evolved states |
| Full material-flow covariance and zero-resolution quantum geometry | Open | Noise/calibration transport and fixed-h ell->0 need further work |
| Physical apparatus, autonomous gravity, spacetime curvature, 3+1 | Open | No new gravitational Hamiltonian is introduced |

## Normalization and domain checks

1. **Basis.** This packet uses orthonormal sqrt(2) sine/cosine modes, with
   lambda_n=1/(2n). It does not reuse the unnormalized coefficient variance
   h/n from the finite mass-chart packet without conversion.
2. **Noise.** Source configuration variance is h C, added measurement
   variance is h C, and vacuum u variance is 2h C. The physical v spectral
   variance and added variance are each h/(4lambda_j) in the vacuum.
3. **Range.** The parent Gaussian isometry exists for any h>0. The Wick
   density and geometry claims here use 2h<1, so the existing L2 proof
   applies. This is a dimensionless analytic range, not an experimental
   upper bound on Planck's constant or a no-go outside the range.
4. **Conditional calibration.** At finite N, averaging the xi part of
   exp[(s+xi)_N-h K_N] leaves exp[s_N-h K_N/2]. The limit uses an L2
   contraction. Positivity then extends the first-moment identity to all
   normal states as an equality of extended nonnegative forms.
5. **Mean versus measurement.** Equality of first-moment forms gives the
   original mean operator J_h(w), not its sharp outcome law. The noisy
   outcome covariance is exp(2h K), while the source covariance is exp(h K).
6. **Density reconstruction.** The coarse density is k_ell*Mhat, not
   exp(R_ell u). Its classical mean is R_ell(exp s_*). Positivity of the
   circle heat kernel gives a strictly positive sample minimum at ell>0.
7. **Compatibility.** All labels and all reconstruction functions live
   on one kappa_h space. Gaussian marginal consistency is exact. Nonlinear
   reconstructed coframes are jointly readable functions of that parent,
   but are not claimed to be related by deleting coframe Fourier coefficients.
8. **Generator.** With [partial,zeta]=1, the physical Wick exponential is
   exp(g.zeta) exp(g.partial). Its action shifts F to F(zeta+g), with the
   factor exp(g.zeta). Self-adjointness comes from W J_N W^* on Ran W.
   No pointwise infinite vector g_z or uncontrolled exponential limit is used.
9. **Pulse variance.** D0 U(t)1=U(t)(t J0-i s0/2), so the real/imaginary
   norm cross term vanishes. Adding the detector's h/4 variance gives
   h/2+t^2[I(h)-1]. Substituting the noisy density as a classical shift
   instead uses I(2h), which is strictly larger.
10. **No scalar phase substitute.** Multiplication by exp[-it Mhat(w)/h]
    on the whole outcome L2 space leaves every outcome probability fixed.
    The physical generator is its positive form compression for w>=0,
    and its unitary acts in the proper analytic range.

## Limit and moment checks

Weak convergence of positive measures plus their bounded sample total masses
gives C^r convergence after convolution with a fixed smooth kernel. Gaussian
heat-filtered v series converge in every Sobolev norm. This controls the
inverse chart before compressing bounded readouts by W. The statement holds
at fixed ell; it makes no assertion of a raw quantum metric at ell=0.

The simultaneous cutoff dynamics/readout limit uses strong unitary convergence
and uniform boundedness of the readout effects. Uniformity on compact pulse
intervals is on each fixed input vector. It is not an operator-norm estimate
for the unbounded source generator and not a joint trajectory measurement.

For classical recovery, the density Fourier variance bound is uniform in
frequency and O(h). Summing the H^3 heat weights costs ell^(-7). The v
variance grows linearly in frequency and costs ell^(-8). These estimates
give a sufficient, nonoptimal simultaneous scale condition. The accepted
length moment proof additionally uses a Gaussian exponential bound for v;
that hypothesis is available for the unpulsed coherent family, not generally
after the nonlinear source pulse. Only second moments are used for the
pulsed bounded-geometry recovery. Those two conclusions are kept separate.

The rejection boundary is null because the independent u0 multiplies the
coarse density by e^(u0-h). Dense Hermite-domain length moments follow from
the threshold bound and Gaussian v norm moments. The ungated mean diverges
by the same zero-mode double-exponential mechanism. Domain invariance under
the pulse is a further question, not a consequence of finite sample lengths.

## Next concrete task

Construct and audit a material-flow-covariant *calibration family* for this
same noisy parent: specify how the detector Gaussian covariance, background
shift, Wick subtraction and smoothing kernel transform. Test both the
configuration first moment and the density first moment against the actual
T_g, and intertwine the same source J_h(w) pulse. The mismatch in (23) shows
why applying the physical affine coefficient to noisy u is insufficient.

Keep separate the fixed-h removal of ell and accepted-length domains after
nonlinear pulses. Neither is needed to state the existing coarse result,
but both must be addressed before stronger quantum geometry claims. A local
apparatus and an autonomous gravitational generator remain independent gates.
