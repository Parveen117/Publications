# Source, error-budget and implementation review

Author: Monty. 30 October 2026. Internal review, not independent peer review.

1. **Actual target.** The evaluated object is the frozen uniform-vacuum
   characteristic expectation, and the reported probabilities are explicitly
   the cosine/sine binary postprocessings. They are not a threshold event,
   the full (u,v) law, a coframe distribution or experimental frequencies.
2. **Conditional mean, not an arbitrary quadrature node.** Each node is
   E[J|bin], obtained by Gaussian tilting. Thus the first Taylor term cancels
   exactly, allowing the quadratic continuum-error bound. Using the central
   q-value instead would invalidate that cancellation and is not implemented.
3. **Tails included.** Both unbounded q tails have their own probability and
   conditional mass. No discarded tail mass or postselection normalization
   is hidden in the computation.
4. **Full source budget.** I(gamma)-sum A_b^2/P_b includes all nonzero modes
   and finite-bin variation. It is not just numerical integration error.
   The limit I(gamma)-exp(gamma) is an irreducible budget for this projection;
   fine central bins cannot establish arbitrarily small errors.
5. **Different noise strengths.** The actual characteristic uses gamma=h and
   the nonlinear phase argument; the naive comparator uses gamma=2h and tk.
   Each source has its own second-moment and projection error. All selected
   gamma values remain below the L2 threshold 1.
6. **Complex error applies to each component.** The Taylor remainder bounds
   the modulus. It therefore bounds the real and imaginary errors separately.
   The Gaussian damping and factor 1/2 are retained in probability intervals.
7. **NRQC execution boundary.** Nine actual calculator calls produce exact
   phase X/Y expectations. Native code performs source conditioning and cubic
   interpolation. Negative interpolation coefficients are postprocessing,
   not probabilities or source amplitudes. No infinite-field circuit or
   quantum computational advantage is claimed.
8. **Independent comparison.** NRQC phase nodes and their assembled discrete
   expectations are compared with separately implemented rational Taylor
   sine/cosine intervals. Time reversal and exact zero-pulse readings are
   additional analytic checks. Interval overlap is a consistency check;
   certification itself rests on the separately stated remainder bounds.
9. **Phase remainder.** The four interpolation nodes are spaced pi/32;
   period wrapping uses exact sign identities. The phase coordinate interval
   is checked before applying the polynomial remainder bound. The finite
   decimal display of a calculator answer is used only with its outward
   enclosure, never as an exact real number.
10. **Gaussian and elementary functions.** Signed directed rounding is used
    throughout. CDF series can have cancellation; 192-bit interval arithmetic
    encloses it. The omitted term is in the decreasing alternating tail in
    the asserted range. Pi, square roots, exp and trig have explicit tails
    or integer bounds. The code rejects requests outside its declared domain.
11. **Probability positivity.** The known original-source probability lies in
    [0,1], so intersecting its rigorous interval with [0,1] is sound. The
    interpolated full k-function is not asserted positive-definite or claimed
    to describe a separately prepared state. Piecewise constant source
    replacement would not retain the translation covariance needed for (1).
12. **Scope of success.** PASS means the checks passed and the enclosures were
    computed. A separate per-case status says whether the requested accuracy
    was achieved. Disjoint actual/naive intervals certify that comparison;
    overlapping intervals leave it unresolved in this execution.

The frozen characteristic/common-parent/LMD packets and NRQC implementation
remain unchanged. The new work is an executable expectation algorithm and
its error transfer, not an independent validation of those continuum axioms.
The immediate refinement is to include nonzero modes in the conditional
source compression for cases whose source-error budget remains too large.
