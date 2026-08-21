# RST-1 — Recognition-Seam Topology: the exact spine, certified

Source: **"Recognition-Seam Topology and Daggered Index Geometry"**
(M. Dabas, July 2026, arXiv source v1), sections 2–7. The paper is the
mature front door of the RH programme (MP PR #94 lineage) and carries its
own claim ledger; this capsule certifies the ledger's *proved* finite,
circle-operator and model statements in exact rational arithmetic — no
floats, no transcendental evaluation, no truncation posing as a limit —
and records the analytic layer as declared-not-certified.

## Blocks

| Block | Content | Key exact fact |
|---|---|---|
| T1 | Cut-stable topology, Eye quotient, seam, curvature defect (§2) | All 64 subsets checked; τ = {∅, {0..3}, {4,5}, X} is a genuine non-discrete topology; seam = zero-cost locus with an off-seam point of *infinite* cost; the defect fires exactly where an unlawful declared composite breaks recognized path equality |
| T2 | Rotation instead of gluing (§3.2, §4) | Orientation target factors through **no** map on the glued quotient; one bit of sheet label repairs it; recognition equality ≠ state identity; reflected side target is target-null exactly on the fixed locus |
| T3 | Rotation-induced involution (§3.3) | J_U self-adjoint involution with orthogonal projections for the rational 3-4-5 rotation; **bridge:** J_{U=I} is the cut swap — fourth certified presentation of the derived involution (EMK-1 T1, EMK-G1 T1, EMK-T1 T1) |
| T4 | Flow generator (§3.4, §5) | Cayley loop as an exact polynomial matrix in h: **L_h = I + h²[D₁,D₂] exactly** (the nilpotent loop terminates at degree 2 — stronger than the paper's O(h³)); commuting control gives L_h ≡ I; eigenmode law certified term-by-term with no exponential evaluated; covariance defect ⟺ conjugation invariance |
| T5 | Derived dagger algebra (§6) | R₀\* = −R₀ but R₀♯ = +R₀ over Gaussian rationals: dagger self-adjointness is **not** Hilbert self-adjointness; leakage three-way equivalence on all 256 grid matrices |
| T6 | Circle carrier (§7) | **sf = Wind = ind = q** for q ∈ {−3,−1,0,1,2,5}: rational crossing times with a proof the window is exhaustive for the *full* Fourier spectrum; winding by exact signed crossings on the Pythagorean parametrization; Toeplitz kernel/cokernel read off the shift |
| T7 | Staircase negative control (§5.8) | Length exactly 2 for every n while sup-distance² ≤ 1/(2n²): limit of lengths ≠ length of limit; not evidence of curvature |

## Findings

- **F1 (derived involution, fourth presentation).** The rotation-induced
  J_U at U = I *is* the cut swap. Geometry involution (EMK-G1), algebra
  seam reflection (EMK-1), tensor cut-swap (EMK-T1) and now the
  sheet-exchange of rotation transport are one derived object. The
  involution is never primitive — the paper's derived-involution
  principle is now a certificate, not a remark.
- **F2 (curvature = non-closure, without limits).** The commutator-loop
  residue is a polynomial-coefficient identity, and for nilpotent
  generators the loop *equals* I + h²[D₁,D₂] with no higher terms.
- **F3 (index without truncation).** The circle identity
  sf = Wind = ind = q is certified with every count exact and exhaustive.
- **F4 (memory is one bit, minimally).** Gluing destroys orientation
  targets; the minimal repair is a single binary sheet label.

## Claim boundary

Analytic theorems of the source (damped trace-class family Thm 8.8,
undamped strip classification Thm 9.1, Schatten thresholds, archimedean
defect) are **declared, not certified** — they are analytic, outside
exact rational verification, and are not consumed. No arithmetic meaning
is assigned to q. The represented seam is **not** identified with the
critical line. Model zeros are **not** zeros of ξi. **Riemann
hypothesis: ABSTAIN** — exactly as the source paper states.

## Next

RST-2 will certify the finite arithmetic/transfer/fusion spine
(§8–14): Möbius pressure extraction, radical recognition current,
prime-shell transfer and asymmetric leakage, zero-identification
boundary, and the Surya/Weil sign fusion on rational-eigenvalue models.

## Reproduce

```
python certificates/rst1_recognition_seam_topology.py
python -m pytest tests -v
```

The certificate regenerates byte-identically; its sha256 is pinned in
`certificates/EXPECTED_RST1.sha256`.
