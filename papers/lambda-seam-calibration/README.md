# Lambda-Seam Calibration (LAM-1)

Calibration and audit of the **legacy Λ-analytic layer** — the
Λ-Eisenstein / Λ-mock functional-equation machinery that grew out of the
T-V-S-P compass diagram of *Geometric Completion of Thermodynamic
Response* (arXiv:2603.20773) — against exact, executable arithmetic.

The legacy program carried the compass response block

```
L = [[λp, Lpv], [Lvp, λv]],   φ = arg(Lpv / Lvp)
```

(the Onsager non-reciprocity bias phase) into automorphic machinery,
producing "Λ-functional equations" of the form
`Ξ_φ(s) = e^{iφ} Ξ_φ(1−s)`. That layer is exactly where the legacy
framework went silent on RH and Yang–Mills. This capsule certifies the
finite content of that layer and names the mechanism of the silence.

## Certified blocks (LAM-1)

| Block | Statement | Verdict |
|---|---|---|
| T1 | Ramanujan two-route: exact cyclotomic reduction of `c_c(n)` in `Z[x]/Φ_c(x)` equals `Σ_{d|gcd(c,n)} d·μ(c/d)` on a 1476-pair grid — the layer the Λ-twist multiplies but never enters | PASS |
| T2 | Memory-lattice multiplicativity: `χ_φ(n)=e^{iφΩ(n)}` is exact integer bookkeeping — Ω additive on all `mn ≤ 400`; formal Euler product reproduces every `n ≤ 300` with marker exactly `Ω(n)`, multiplicity exactly 1 | PASS |
| T3 | Theta seam certified: `θ(1/t)=√t·θ(t)` at `t=2` with two-sided directed rational enclosures, widths `< 1e-20` (actual agreement 30+ digits), **no floats in any verdict path** | PASS |
| T4 | Declared-phase audit: under the exact S-reindexing bijection `(c,d)→(d,−c)` on the primitive region `|c|,|d| ≤ 60` (8812 pairs), the coset-weight layer has **zero net change** (16512 = 16512) | PASS |

Controls: coprimality tamper separates (T1); non-additive fake winding
detected (T2); `√3` tamper separates the enclosures with certified gap
`> 0.319` (T3); asymmetric region separates (T4).

## Finding LAM1-F1

The `e^{iφ}` functional-equation factor of the legacy Λ-construction is
a **declared** right S-multiplier, not a derived quantity. The legacy
canvas itself records that no homomorphism `SL(2,Z) → S¹` can realize
arbitrary φ (finite abelianization), and T4 certifies that the
coset-weight bookkeeping contributes zero net phase across the S-move.
The legacy analytic layer therefore decorated the classical substrate
with an *input* phase — which is the precise mechanism of its silence on
RH/YM: **the phase was never load-bearing**. Any native functional
equation (obligation N2) must *derive* its seam-crossing phase; the
theta seam of T3 is the classical shape of that single crossing.

## Claim boundary

- Legacy Λ-FE: declared multiplier over the classical substrate.
- N1 (native continuation), N2 (native functional equation),
  N3 (identification): **OPEN**.
- K₀ / L₀ / RH: **OPEN**. YM continuum gates: **OPEN**.
- No claims about zeros of any L-function.
- Classical anchors (Poisson summation, Jacobi theta inversion,
  SL(2,Z) coset structure, Ramanujan-sum closed form) are pinned named
  dependencies, cited and never rederived natively.

## Reproduce

```
python papers/lambda-seam-calibration/certificates/lam1_seam_interface.py
python -m pytest papers/lambda-seam-calibration/tests -v
```

The generator is deterministic; CI regenerates the certificate and
fails on any pin drift (`EXPECTED_LAM1.sha256`).

## Arithmetic discipline

Integers, exact rationals, exact integer polynomial arithmetic, and
directed rational intervals with relative outward rounding to 300-bit
significands. No floating point in any verdict path.
