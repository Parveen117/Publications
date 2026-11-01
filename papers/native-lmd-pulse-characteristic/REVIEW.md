# Proof and implementation review

Author: Monty. 30 October 2026.

This is an internal source/domain review, not independent peer review.

1. **Translation sign.** D=i h partial_q-iq/2 becomes +i h partial_q.
   exp(i k D) translates the flattened ket by -hk. The bra's pulse phase
   minus the ket's phase is +2t sinh(hk/2) J/h, hence E v_0=+t.
2. **Two Gaussian factors.** The amplitude overlap contributes exp(-hk^2/8).
   The common parent's added variance h/4 supplies another exp(-hk^2/8).
   Dropping either one gives the wrong zero-pulse variance h/2.
3. **Bounded proof before moments.** The characteristic formula uses unitary
   translations and bounded phases. Its proof never interchanges an infinite
   moment expansion of the lognormal/GMC mass with an integral. The optional
   formal third-order fixture uses arbitrary finite moment symbols and is not
   a claim that an arbitrary GMC moment exists.
4. **Source and noisy mass.** The actual law uses chi_h at a nonlinear
   argument. The rejected additive label rule uses chi_(2h) at tk. Both
   discrepancies matter: even substituting chi_h in the additive rule leaves
   the bounded linear-response mismatch.
5. **Cutoff theorem.** J_N is conditional expectation of the L2 limit J.
   Therefore ||J-J_N||_2^2=E J^2-E J_N^2. Positivity of that difference is
   not inferred from a false pointwise monotonicity of K_N.
6. **Density regularity.** Gaussian domination of the characteristic function
   justifies Fourier inversion and every spatial derivative. The additional
   sinh factor in the cutoff error is integrable against that Gaussian.
   Smooth bounded probability density does not imply finite coframe length
   moments or settle any earlier geometric domain issue.
7. **Source interval arithmetic.** Positive zeta integral tails enclose each
   retained series coefficient; the n>80 tail is bounded by
   2g^81/[81(1-g)]. Exponential Taylor tails use a geometric ratio bound.
   All interval arithmetic rounds outward on an 80-bit grid. Matrix control
   coefficients are 28-bit dyadics, with explicitly checked source errors.
8. **Cutoff arithmetic.** All Laurent coefficients of the exponent are
   nonnegative; rounding every convolution contribution preserves enclosure.
   The 48-term tail is bounded by the scalar exponential at h H_N. The
   independent N=1 Bessel series checks a separate scalar implementation.
9. **Actual calculator interface.** Inputs contain only declared gates and
   signed Pauli words. Weighted controls and time shears are assembled outside
   NRQC by the native exact engine. Four source-free state vectors provide
   the independent matrix baseline; only one state carries source moments.
10. **Finite representation boundary.** Trace([D_f,J_f])=0 rules out the
    positive-affine identity with trace(J_f)=4. No finite source time evolution
    is claimed from a supplied covariance shear. J_f is not asserted positive.
11. **Negative controls.** Wrong phase sign, missing detector noise, the naive
    density shift, finite affine identification and naive bounded response
    must all fail. Zero pulse is a positive equality control, not counted as
    evidence of a separation.
12. **Provenance and limits.** The calculator commit is pinned and its tracked
    tree must be clean. Result/source hashes are recorded. Runtime values are
    process observations only; there is no speed or memory comparison claim,
    hardware job, circuit-scale inference or CI dispatch.

Initial integration exposed two output-parser assumptions: NRQC's success
label is `completed`, and mean coefficients are nested under `exact`.
Those schema adapters were corrected before the final clean run. No physical
formula or calculator result was changed to make a failed comparison pass.

Still open: evaluation of the oscillatory source characteristic expectation,
joint/nonuniform pulse measurements, material calibration covariance, removal
of smoothing at fixed h, post-pulse accepted-length domains, autonomous
gravity and 3+1 physical recovery.
