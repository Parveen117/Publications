# Lineage — EMK recognition geometry

## Origin

The rotational seam metric lives in the vault's
`emk_ugd_recognition_geometry` tree as LaTeX: definitions, theorems and
proofs, never executed. As with EMK-1, EMK-2 and UGD-1, no executable
version existed anywhere in the corpus. EMK-G1 is the first
machine-checkable realization.

## What EMK-G1 adds

The source's mathematics is reproduced, not extended. What is added:

- an arithmetic reformulation that makes the whole geometry exact:
  carrying W = A^2 turns the connection, the Riemann component and the
  Gaussian curvature into rational functions, so K(0) = -kappa is
  certified with no square root and no floating point;
- a two-route check of the curvature (R_{uvuv}/det g against the closed
  form) and of the Christoffels (W-form against Levi-Civita computed
  independently from the metric components);
- a control for the geodesic seam criterion: an odd perturbation of the
  warp breaks reflection symmetry and destroys the geodesic property,
  showing the symmetry is load-bearing rather than decorative;
- a certified SEPARATION of metric curvature from recognition curvature
  (T6), turning the source's stated principle into a machine-checkable
  guard.

## Relationship to the other capsules

- **EMK-1 / EMK-2** (`../emk-ugd-algebra/`): the algebra identifies the
  seam and grades memory as the odd sector; this capsule gives the seam
  its metric. T6 certifies that the two curvatures — Gaussian here, odd
  grade there — are independent, so the layers may be composed but never
  conflated.
- **Cut-First Equivalence** (`../cut-first-equivalence/`): the
  recognition two-form of CFE is one of the two objects separated in T6.
  Nothing in this capsule licenses identifying it with the metric
  curvature.

## Claim boundary

The metric family is declared, as in the source. K(0) = -kappa is a
theorem for that family and not a universal identification of flow
parameters with curvature scalars. Global, completeness and analytic
statements are not claimed. Untouched: RH / K0 / L0 / YM continuum
gates, and quantum gravity.
