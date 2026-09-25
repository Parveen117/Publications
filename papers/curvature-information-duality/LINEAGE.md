# Lineage — CID-1

- **Consumes (pinned):** `papers/cut-first-equivalence/certificates/`
  `cfeu_uniqueness.py` (memoryless space ker d₁, residue functional,
  response complex, potential reconstruction) and through it
  `cfe1_cut_first_equivalence.py` (witness family, BETA, LOOP_RECT,
  circulation). Nothing is re-declared; a change in either pinned
  module breaks this capsule's tests by construction.
- **Companion context:** `papers/information-invariance/` (D-series:
  δI = 0, exactness of ω, selection principle) — the conceptual
  neighbour, not a code dependency.
- **Vault cross-reference — DIFF DONE (Aug 21).** The vault appendix
  `APPENDIX_GENERALIZED_EULER_INFORMATION_CURVATURE_DUALITY.tex` (259
  lines) has now been read and compared against CID-1. **Result: the
  appendix is largely declared structure, and CID-1's theorems are not
  in it.** The appendix declares a generalized Euler phase-seam space,
  a recognition connection, an information-curvature tensor
  I^GE_ab = <grad rho, grad rho> + <grad theta, grad theta> + M^Sigma,
  a curvature-residue norm, and a duality *map* D_Sigma sending
  curvature components to closure-residue components; its theorem is
  the qualitative statement that visible return is not recognition
  return when any active component is unclosed. CID-1 proves something
  the appendix does not state: that the obstruction is **invariant**
  under change of information metric while the Hodge split is not, and
  that the g-orthogonal complement has rank exactly 1 for every
  positive-definite metric. So CID-1 stands as written; the appendix's
  tensor I^GE remains **declared, not certified**, and is recorded as
  such in the EMK-TOP-1 claim boundary. The appendix's sheet-rotation
  class and visible-return statement *are* certified — in EMK-TOP-1
  T5 and T1.
- **Prior executable version:** NONE — first bridge capsule between
  the curvature and information layers.

## QTH-1 scope correction, 25 September 2026

The original T3 commuting and noncommuting witnesses remain valid. The
overgeneralized converse has been removed: a nonzero SLD commutator can
have zero state average, explicitly at rho=I/2 with Pauli x/y derivatives.
QTH1_RESULT.json and its expected digest are regenerated with this negative
control. A difference defined as QFI-CFI always gives an algebraic equality;
its positivity and equality conditions are the substantive measurement claims.
This correction does not derive a physical quantum theory from native inputs.
