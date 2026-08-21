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

---

# RST-2 — Arithmetic and fusion spine (§8–14), certified

Same folder, second capsule, 35 further tests (70 total in this folder).

| Block | Content | Key exact fact |
|---|---|---|
| T1 | Sampled reflection, spectral blindness (§8.2–8.3) | ρ(s)=1−s swaps the rational sample pairs exactly; seam projection provably independent of the anti-seam component |
| T2 | Weights, leakage norm, crossings (§8.4, §8.7) | Equal pair weights ⟺ swap isometric (both directions); finite leakage norm² = max aₖ² attained; avoided family has det = −((η−c)²+µ²) ≤ −µ² < 0 — zero never in spectrum |
| T3 | Arithmetic recognition (§10) | **The formal-log move:** log n carried as the vector of prime exponents, so (L1)∗µ = Λ is an exact integer-vector identity for all n ≤ 60 — no logarithm evaluated. Λ-current 0 on every anti-seam class; contamination → exactly ϵe₆ |
| T4 | Prime-shell transfer (§11) | U A = A P, JA = A, P₋A = 0, JU\*J = U exactly; asymmetric leakage = ((w_R−w_L)/2)(e_R−e_L) per shell; **converse-fails witness:** a common seam-sector error is invisible to the anti-seam current |
| T5 | Zero-identification boundary (§13) | Paired determinant: zeros ±λⱼ **with multiplicity by exact division** (double level → double zero); log-derivative identity in cross-multiplied polynomial form; nonvanishing factor preserves the divisor; working Z1/Z2/Z3 obstruction validator |
| T6 | Surya/Weil fusion (§14) | Feshbach congruence M = Sᵀdiag(A,F)S as an exact identity with negative witnesses transported through it; **Surya without square roots** (squared amplitude x²/(1+x²) on x<0): A_SS = 0 ⟺ W₀ ⪰ 0; seam charge additive, = −spectral flow, every crossing downward; µ = 1 strictly outside the charge; principal holonomy blindness (m mod 1 = 0 ∀m) |

**Findings.** F5: the projection-blindness thread now spans four capsules
(UGD-1 T4, EMK-T1 T5, RST-2 T1, principal holonomy) — the visible
coordinate never carries the full state. F6: the anti-seam current
detects only the *antisymmetric* discrepancy — the seam sector needs its
own audit (converse-fails witness). F7: the Surya sign problem is
decidable entirely in squared-amplitude arithmetic — no functional
calculus square root is ever needed for the vanishing question.

**Claim boundary.** Gaussian damping form and 1/√2 normalization are
declared; W₃, Weil↔RH, the completed parity/Feshbach reduction and the
even compact operator are imported by the source and NOT certified here;
fusion gates F1–F5 remain OPEN; Zmodel→Zξ forbidden. **RH: ABSTAIN.**

```
python certificates/rst2_arithmetic_fusion_spine.py
python -m pytest tests -v
```
