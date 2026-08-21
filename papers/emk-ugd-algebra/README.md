# EMK / UGD algebra — the determinant seam ladder (EMK-1)

The framework's **own algebraic core**, brought from the vault into the
public certified corpus. Until now the Cut-First Equivalence capsules
used an explicit rational equation of state as their witness — a
scaffold. This is the theory's native geometry instead.

## The primitive algebra

```
K² = I      seam reflection
R² = −I     circular phase transport (the quarter turn)
RK = −KR    mixed rotation-cut parity,  (RK)² = I
[R, K] = 2RK                             the UGD bracket
```

A primitive block is `M = aI + bK + cR + dRK`.

## The determinant identity

> **det M = (a² − b²) + (c² − d²) = Δ∥ + Δ⊥**

The determinant splits **additively** into a seam-compatible channel and
a rotational / anti-seam channel, and the two close or fail
**independently**.

## Certified blocks

| Block | Statement | Verdict |
|---|---|---|
| T1 | Primitive relations K²=I, R²=−I, RK=−KR, (RK)²=I, and the UGD bracket [R,K]=2RK | PASS |
| T2 | The additive determinant identity on the full rational grid (2401 quadruples), two independent routes; **channel independence** exhibited — each channel can vanish alone, and both can vanish while the operator is nonzero | PASS |
| T3 | Seam projectors P±=(I±K)/2 idempotent/orthogonal/complete; `[A,K]=0` **iff** the mixing blocks vanish — the commutator *is* the seam-mixing detector | PASS |
| T4 | Exact Schur identity + the quadratic determinant correction; **control:** both diagonal blocks have det 1 yet the quadratic residue is 1/6 — first-order closure is provably incomplete | PASS |
| T5 | Cut-commutator curvature `C_iC_x − C_xC_i = −Δx·Δi·[D_x,D_i]` exactly; path equality ⟺ vanishing commutator channel | PASS |
| T6 | Winding index an exact integer, homotopy invariant (1, 1, 0 across scaled and off-origin loops) — **non-healable curvature**, via rational circle parametrization and exact crossing counts, no trigonometry | PASS |
| T7 | **Bridge to CFE:** a block is memoryless (commutes with K) **exactly** when rotation-free (c=d=0), and then det M = Δ∥ alone | PASS |

## EMK-2 — CFE rebuilt on this algebra, and the UGD multiplicative law

EMK-1 T7 gave the licence; EMK-2 performs the replacement and retires the
scaffold equation of state.

**The algebra is ℤ/2-graded:** even (seam) sector `{I, K}`, odd
(rotational) sector `{R, RK}`, with `even·even = even`, `odd·odd = even`,
`even·odd = odd`. The determinant channels read the two grades — so
**memory is the odd grade**, not a dial bolted onto a model.

**The UGD multiplicative law lives in the seam channel.** In the Cayley
coordinate `y = b/a`, even-sector composition is exactly

```
y₁ ⊕ y₂ = (y₁ + y₂) / (1 + y₁y₂)
```

and under the native chart `x = (1+y)/(1−y)` this is precisely
`Log_Σ(x₁x₂) = Log_Σ(x₁) + Log_Σ(x₂)` — the framework's own
multiplicative-to-additive law (Recognition-Kernel-Framework, F00G
Theorem 7.1). Certified with **no series, no transcendental evaluation
and no floats**. The seam determinant channel is correspondingly
multiplicative, `Δ∥(MN) = Δ∥(M)·Δ∥(N)`, while its coordinate is additive:
the determinant split and the multiplicative law are one grading seen
twice.

| Block | Statement | Verdict |
|---|---|---|
| T1 | The ℤ/2 grading, with explicit composition laws for each sector | PASS |
| T2 | The multiplicative law in the Cayley chart, exact; seam channel multiplicative; naive-addition control separates | PASS |
| T3 | CFE rebuilt natively: memoryless ⟺ even grade; zero curvature, zero residue | PASS |
| T4 | Residue on the native carrier: ∮ω=∬Ω exact, zero iff even grade, strictly monotone in the odd amplitude | PASS |
| T5 | **The dial was the odd grade:** under ρ = χ−1 the scaffold and native densities are identical — the residues reproduce CFE-1's certified −32, −16, 0, +16, +32 exactly | PASS |

## UGD-1 — UGD number

The most original layer: a **seam-aware analogue of positional notation**.
A classical digit records a scalar at a scale. A UGD digit records
**phase, scale index, and local seam charge**:

```
d_n = (φ_a, n, s)  ∈  Φ_K × ℤ × {−1, 0, 1}
```

Phase is carried as an **exponent in ℤ/K** — never evaluated as a root of
unity — and the scale base is an exact rational, so no transcendental
arithmetic enters any verdict.

| Block | Statement | Verdict |
|---|---|---|
| T1 | Phase exponents form ℤ/K with the generator of full order; **at K=4 the generator has order 4 and its square is the half turn** — the quarter turn, at the numeral level | PASS |
| T2 | **Carry conservation:** phase overflow carries to the next scale, seam overflow returns to {−1,0,1} and is pushed to the ledger; total seam charge is conserved exactly (1296 digit pairs). Control: a clamping policy that drops overflow loses charge and separates | PASS |
| T3 | Multiplication is **scale convolution** — scale indices *add*, phases multiply mod K, so π_scale is multiplicative: the multiplicative-to-additive law at the numeral level | PASS |
| T4 | The classical projection is **blind to seam**: null-seam recovers ordinary positional notation exactly, while three numerals share the projection 7 with total seam charges 0, 0, 2 | PASS |
| T5 | **Cut-zero:** the neutral digit's neutrality is *derived*, its projection and ledger are null — yet it occupies a scale and carries a phase exponent, so it is distinct from the *absence* of a numeral. Zero is the first cut | PASS |
| T6 | λ-logic: negation `N(φ,σ,k)=(φ,−σ,−k)` is an involution; the seam predicate marks exactly `k≠0`; the null-seam trivial-phase limit recovers classical behaviour (36-point domain, exhaustive) | PASS |

**UGD1-F1** — the quarter turn now stands certified in **four** presentations:
arithmetic (LAM-2 ι²=−1), quantum geometry (CFE-Q J²=−I), native algebra
(EMK-1 R²=−I), and number (UGD-1 at K=4).

**UGD1-F2** — the carry discipline is a *conservation law*; the ledger is
load-bearing, not bookkeeping.

**UGD1-F3** — the classical projection is provably lossy: UGD numbers carry
strictly more than their classical shadows.

## Findings

**F1** — the additive split is *structure, not notation*: the channels are
genuinely independent.

**F2** — first-order closure is provably incomplete; diagonal-block tests
cannot certify recognition closure.

**F3** — the CFE memory dial and the EMK rotational determinant channel are
the same object in two presentations, so CFE's witness EOS can be
replaced by this native algebra.

**EMK2-F1** — memory *is* the odd grade; the memory dial was never an extra
parameter.

**EMK2-F2** — multiplication of native scalars is addition of seam
coordinates; the determinant split and the UGD multiplicative law are one
grading.

**EMK2-F3** — the scaffold retires without cost: every certified curvature
and residue value is unchanged under ρ = χ−1.

## Provenance

Source: vault LaTeX — the EMK algebraic core ladder I–VII and the UGD
algebraic-operators appendix. **No executable version existed anywhere**
in the corpus (checked Recognition-Kernel-Framework and RH-Framework:
neither contains a UGD algebra implementation). This capsule is the first
machine-checkable realization.

## Claim boundary

Certified: the finite-dimensional theorem spine. **Not claimed:** the
infinite-dimensional trace-class, essential-spectrum, Hurwitz-zeta,
ξ-function or Hilbert–Pólya extensions — the vault source explicitly
defers those, and so does this capsule. CFE's (U) uniqueness stays open.
No RH / K₀ / L₀ / YM continuum gate is touched; quantum gravity is not
touched.

## Reproduce

```
python papers/emk-ugd-algebra/certificates/emk1_determinant_seam_ladder.py
python papers/emk-ugd-algebra/certificates/emk2_native_carrier.py
python papers/emk-ugd-algebra/certificates/ugd1_numerals.py
python -m pytest papers/emk-ugd-algebra/tests -v
```

CI regenerates the certificate and fails on pin drift
(`EXPECTED_EMK1.sha256`, `EXPECTED_EMK2.sha256`, `EXPECTED_UGD1.sha256`).
