# LDC-1 retained actual-source execution

The inputs are the complete outward source jets already certified by
LAA-1 at original s=29,31,33. Their certificate and source pins are
retained unchanged. This stage reuses those jets rather than adding
another pair of derivatives. The integrated replay independently
recomputes the predecessor source evaluations.

| Source label s | Anchor a | Available carrier transitions | Exact source record lengths |
|---|---|---|---|
| 29 | 3249/4 | N=1 through 9, each to N+1 | 0 through 19 |
| 31 | 3721/4 | N=1 through 9, each to N+1 | 0 through 19 |
| 33 | 4225/4 | N=1 through 9, each to N+1 | 0 through 19 |

For all 27 transitions, native outward arithmetic checks:

- Positive monic boundary norm h_N, with the LDL pivot overlapping
  its independently evaluated Schur-norm expression.
- Orthogonality of the monic new polynomial to every old monomial.
- Coupling beta_N=h_N/h_(N-1) strictly between zero and one quarter.
- Positive response increment at the same held-out original s=35,
  with agreement of direct subtraction and the determinant formula.
- Agreement of the residual and Schur expressions used in the proof.

The response increment uses normalized determinant ratios before
multiplication. This prevents avoidable loss of tiny products in the
unchanged 110-place outward interval engine. No precision adjustment,
floating-point fallback, midpoint selection or source correction is used.

At depth ten, the held-out responses and positive source errors overlap
the retained LAA results. The law describes increasing carrier depth
at each fixed source anchor; it does not equate different anchors.

Every record row is positive and its probability sum contains one.
All available falling-factorial readings overlap the same source c_r.
For every adjacent record length 1->2 through 18->19, the explicit
deletion coupling has squared discrepancy overlapping
(c_1-c_2)/(m(m+1)), giving 54 finite coupling checks. Length-zero
normalization and deletion from length one are also checked.

The exact toy fixtures separately verify the universal operator
identities, determinant increments for negative, zero and positive
lambda, the first new even moment, and a rational two-mode cut
compression with defect (9/100)I. The formal noncommuting projection
identities replay in the frozen rewrite engine. These fixtures test
the algebra; the written proofs establish its quantified statements.

The positive actual h_N values invoke LDC-T4 to certify nonzero cut
compression defects without calculating matrix square roots. Exact
stabilization applies to the specified branch probabilities, not to
arbitrary coherent amplitudes or a fixed finite matrix at every depth.

The finite record statements are unconditional for the certified
source bank. The infinite L2 statement requires the actual all-order
positive source condition. Physical interpretation remains open and
is not selected by the finite examples.
