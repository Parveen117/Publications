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

## Certified blocks (LAM-2: derived seam phase)

The counterpoint to LAM1-F1 — for a *genuine* twist the seam phase is
**derived**, as finite arithmetic on the cut's residue classes:

| Block | Statement | Verdict |
|---|---|---|
| T1 | Gauss magnitude exact: `|τ(χ)|² = q` in `ℤ[ζ_L]` for all 34 primitive characters, moduli {3,4,5,7,9,11,13} (cyclic convolution + reduction mod `Φ_L`; no floats, no complex numbers) | PASS |
| T1b | Quarter-turn instance: `τ(χ₄)² = −4` exactly — the derived seam phase at q=4 satisfies the cut relation `ι² = −1` | PASS |
| T2 | The derivation step exact: separated-Gauss identity `Σ_a χ̄(a)ζ_q^{an} = χ(n)τ(χ̄)` at **every** residue n for every primitive χ (335 identities), including `gcd(n,q)>1` where both sides vanish — the vanishing *is* primitivity | PASS |
| T3 | Derived phase certified analytically: twisted theta seam `θ_χ(1/2) = √8·θ_χ(2)` at q=4 with derived multiplier `τ/(i√q) = 1`, enclosures agreeing to 30+ digits | PASS |
| T4 | No conductor for the legacy twist: for every modulus `q ≤ 60`, explicit congruent pair with different `Ω` — the winding twist is a Dirichlet character of **no** modulus | PASS |

Controls: imprimitive character mod 9 gives `|τ|² = 0` and fails the
derivation step at residue n=3 (both exact); multiplier tamper
(`√2` for `√8`) separates the certified enclosures with gap `> 0.58`.

## Findings LAM2-F1/F2/F3

**F1 (template):** for primitive-character twists the FE seam phase is
finite cut arithmetic — exact magnitude, exact derivation step, and a
certified analytic seam confirming the derived multiplier. A derived
phase is *load-bearing*: tampering it separates certified enclosures.

**F2 (native hook, remark-level):** `τ(χ₄)² = −4` — the derived phase at
q=4 is exactly the quarter turn: `(τ/√q)² = −1`, the cut's defining
relation, obtained from seam arithmetic rather than declared.

**F3 (impossibility half):** the legacy winding twist
`χ_φ(n) = e^{iφΩ(n)}` has **no conductor** (witness for every q ≤ 60),
so the Gauss-seam derivation route does not exist for it: the declared
phase of LAM1-F1 was *forced* by the twist's own arithmetic. Together,
LAM1-F1 + LAM2-F1 state the N2 design constraint: a native functional
equation must ride a twist whose phase is derivable as finite seam
arithmetic.

## Certified blocks (LAM-3: the crossing, derived)

The derivation of the `s ↔ 1−s` crossing itself, executed and certified
at machine level — theta seam → split-flip identity → manifestly
symmetric completed kernel = Euler-product side — at points where the
native `ζ_Σ` exists. The FE phase for ζ is the derived `+1` of the
untwisted seam: an *output* of the crossing, never an input.

| Block | Statement | Verdict |
|---|---|---|
| T1 | Memory-grammar values: `γ`, `ln 2`, `ln π`, `ζ(3)`, `ζ(5)` as exact-rational enclosures (recognized body + Bernoulli correction window + declared bracketed tail), widths `< 1e-24`, two truncation depths overlap | PASS |
| T2 | Split-flip identity at s=3: `∫₀¹ t^{1/2}ω dt = 1/6 + ∫₁^∞ t^{-2}ω dt`, two independent routes (incomplete-gamma series vs `E₁` route with exact rational cancellation) agreeing to 28+ digits — the crossing step, certified | PASS |
| T3 | Completed identity at s=3: `ζ(3)/(2π) = 1/6 + ∫₁^∞ (t^{1/2}+t^{-2})ω dt` | PASS |
| T4 | Second point s=5: `3ζ(5)/(4π²) = 1/20 + ∫₁^∞ (t^{3/2}+t^{-3})ω dt` — kills point-luck | PASS |

Controls: polar-drop separates (gap > 1/10); seam-weight tamper
`ω(1/t) → t·ω(t)+(t−1)/2` separates (gap > 2/5); reflected-kernel
exponent tamper separates (gap > 1/200); body-only ζ separates
(memory grammar load-bearing); Bernoulli-window tamper separates.

**Note on the E₁ route:** at `c = 16π` the convergent-series route
suffers ~22 digits of cancellation — fatal in floats, *exact* in
rationals. The capsule's positivity and width assertions at that point
are the machine witness that the arithmetic discipline is load-bearing.

## Findings LAM3-F1/F2

**F1:** the crossing is derived end-to-end and machine-certified on the
`Re s > 1` side. Contrast chain complete across the paper:
LAM-1 — declared phase, never load-bearing; LAM-2 — derived phase,
load-bearing finite seam arithmetic; LAM-3 — the derivation itself
executed, with every tamper separating certified enclosures.

**F2:** the memory grammar is load-bearing: body-only ζ(3) separates by
more than `1e-5` at widths below `1e-24` — the public echo of "without
memory terms there is no critical strip."

## Honest boundary (LAM-3)

Continuation into the critical strip via the symmetric kernel is
**definitional** here, not a native theorem. The exact remaining native
obligation for N2 is now pinned: *reproduce the split-flip step (T2)
inside the cut grammar, where the theta seam is a native theorem rather
than a pinned classical anchor.* T3 is T2 composed with exact termwise
Mellin evaluation (recorded per TAUT-1); the tails on both sides share
the Euler–Maclaurin utility.

## Reproduce

```
python papers/lambda-seam-calibration/certificates/lam1_seam_interface.py
python papers/lambda-seam-calibration/certificates/lam2_derived_phase.py
python papers/lambda-seam-calibration/certificates/lam3_crossing_derivation.py
python -m pytest papers/lambda-seam-calibration/tests -v
```

The generator is deterministic; CI regenerates the certificate and
fails on any pin drift (`EXPECTED_LAM1.sha256`, `EXPECTED_LAM2.sha256`, `EXPECTED_LAM3.sha256`).

## Arithmetic discipline

Integers, exact rationals, exact integer polynomial arithmetic, and
directed rational intervals with relative outward rounding to 300-bit
significands. No floating point in any verdict path.
