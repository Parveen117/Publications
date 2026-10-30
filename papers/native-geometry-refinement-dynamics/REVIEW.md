# Written review and closure decisions

Author: Monty. 30 October 2026.

This is an analytic audit of the new derivation, not independent peer review.
No new numerical results, certificate counts or calculator benchmarks are
reported. Source hashes and repository consistency checks are provenance/I/O,
not evidence for a continuum theorem.

| Target | Decision | Boundary |
|---|---|---|
| Finite Fourier mass chart and normalized POVM | Constructed conditionally on LMD | Finite family; full residual Hilbert space retained |
| Exact refinement of complete predecessor records | Obstructed | Adding a conjugate while retaining its old sharp residual configuration |
| Compatibility of coframe-only marginals | Unresolved | The no-go uses labels that the geometry map can discard |
| Bounded length response under resolved density pulses | Derived | Exact operator-norm derivative and remainder for F(L)=L/(1+L) |
| Raw mean length made finite by such a pulse | Ruled out when initially infinite | Two-sided multiplicative bounds preserve integrability status |
| Arbitrary smooth source geometry and pulse orbit | Recovered in an ordered limit | First h->0 at fixed family, then Fourier resolution->infinity |
| Fixed-h continuum geometry and all-flow covariance | Unresolved | No interchange of limits or common noisy parent established here |
| Autonomous gravitational dynamics and 3+1 theory | Unresolved | The generator remains the existing density-probe interaction |

## Analytic checks

1. **Normalization conventions.** Unnormalized sine/cosine coefficients have
   variance h/n, while the orthonormal sqrt(2) modes have variance h/(2n).
   The full removed variance is h sum alpha_a f_a^2. The mass-coordinate
   density uses det H, not the radial x0 Jacobian.
2. **Moment-cone boundary.** Positive H alone gives a local chart. The convex
   minimization proves onto the moment-cone interior; the supporting-face
   argument supplies compact cone support after H localization. These are
   separate necessary steps.
3. **Characteristic sign.** With the plus Fourier kernel, exp(i t p_b/h)
   yields g(x-t e_b). Therefore q_b(x)>q_b(x-t e_b) along an interior segment.
   The positive Gaussian overlap cannot remove the resulting commutator.
4. **Scope of incompatibility.** The coarse record contains q_b sharply;
   the refined p_b effects fail to commute with it. The proof cannot be
   applied without checking that this sharp record is still retained.
5. **Feedback sign.** For pi proportional to nu/a, P'/P=-integral a' w dpi.
   Hence subtracting log P in the coframe contributes the plus sign in (11).
6. **Bounded evolution.** The bounds |r'|<=5M/2, |r''|<=5M^2/2 imply
   |(log L)''|<=35M^2/4. Applying the logistic function gives
   |F''(L_t)| along the path <=15M^2/4, and the Taylor remainder <=15t^2M^2/8.
   Here the notation refers to the second time derivative of F(L_t), not
   the second derivative of F as a function of L.
7. **Arbitrary source center.** At fixed A the classical conjugate is
   (H*)^(-1)D*, which is the weighted projection of d_*/j_*. It equals the
   full profile only when that profile lies in the resolved span.
8. **Arbitrary pulse.** Its finite-h phase gradient is
   H^(-1)[D*+M(B_t f)]. A constant Fourier shift is asserted only for B_t
   in the resolved span. The general pulse is treated through its exact
   phase and LMD concentration, not a false finite-mode covariance identity.
9. **Order of limits.** Constants in the inverse-H, Gaussian and filter
   estimates depend on the selected finite family. No uniform-in-N estimate
   is supplied. The limit is explicitly classical first, refinement second.
10. **Time statement.** Uniform single-time concentration is the supremum of
    error probabilities over t, not a probability on a joint quantum path.

## Source sanity checks, by written calculation

For A={1}, the construction has an empty configuration filter and reduces
to the frozen radial measurement. For A={1,cos(2 pi z)}, it reduces to the
two-window measurement, with the same Wick variance and Fourier sign.

For a flat classical density at y=0 and w=c=cos(2 pi z), (11) gives
r'_0=3c/2. Its derivative gives r''_0=3c^2/4-3/8. Therefore

    L'_0=integral 3c/2=0,
    L''_0=integral [(3c^2/4-3/8)+(9c^2/4)]=9/8.

This recovers the predecessor's quadratic length coefficient 9/16 using
the new evolution equation with its feedback term. No executable replay
was needed or performed for this identity.

For w=0, the generator, B_F and Taylor remainder vanish. Constant positive w
has b(z)=a(theta(z))w-a'(theta(z))w+integral a'(theta)w dpi; it simplifies
to a(theta)w when theta is constant, agreeing with Theta'(y)=a(theta).

## Next task contract

A candidate common-parent measurement must specify its outcome space,
normalization and actual finite marginals on the LMD Hilbert space. It must
avoid the sharp residual/conjugate conflict identified here. Its geometry
map must retain P and distinguish bounded output laws from unbounded length
moments. Its density-pulse response must be derived from the existing J(w),
and material-flow covariance must be checked for that same prescription.

Use the existing CMG noisy-parent method as a credited design reference if
helpful, with an explicit LMD covariance and Wick-density audit. A new list of
independently normalized cutoff POVMs, or the ordered classical limit alone,
does not discharge this common-parent requirement. The autonomous gravity
generator remains a separate source problem.
