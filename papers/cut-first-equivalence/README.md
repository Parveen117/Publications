# Cut-First Equivalence — certified core (CFE-1)

The candidate **flagship theorem** of the recognition / cut-flow
framework (arXiv:2603.20773), with its memoryless-limit direction and
residue identity certified in exact rational arithmetic.

**Read [`THEOREM.md`](THEOREM.md) first** — it states the theorem, marks
certified vs. open, and maps every obligation to its scaffolding.

## The statement

> Classical thermodynamics is the **memoryless sector** of recognition
> thermodynamics, and the loop residue **∮_∂D ω = ∬_D Ω** is the
> obstruction. Equilibrium is flat (Ω=0); irreversibility is curvature
> (Ω≠0).

A bridge identity of the Noether type: recognition/memory structure ↔
thermodynamic response geometry, with a computable residue.

## Certified blocks (CFE-1)

| Block | Statement | Verdict |
|---|---|---|
| T1 | Memoryless limit is classical: at χ=1 curvature vanishes on the whole grid, every Maxwell loop closes, Onsager asymmetry is exactly 0 | PASS |
| T2 | Residue is exact: ∮_∂D ω = ∬_D Ω in ℚ (boundary walk vs area sum, independent), both = 32 at χ=7/5 | PASS |
| T3 | Obstruction is faithful: Ω=0 iff χ=1; residue strictly monotone in μ=χ−1 (−32,−16,0,+16,+32), zero exactly at μ=0 | PASS |
| T4 | Equilibrium invariant 𝓘=Γ_c·Γ_m=1 exactly at χ=1, departs when χ≠1 | PASS |

All exact rationals, no floating point in any verdict path.

## Surjectivity — DISCHARGED (CFE-2)

**(S)** the cut grammar maps *onto* the classical response algebra: all
four Maxwell relations as cut-closedness, the response-coefficient
identities (γ ratio law + Mayer), and a rank-equality surjectivity
certificate (5 = 5, cokernel 0) built by cut-native finite differences
independent of the classical partials — renaming test passed, non-closed
control separates.

## Quantum witness (CFE-Q)

The same bridge holds quantum-mechanically, certified on a dissipative
qubit: memory dial χ ↔ non-Markovianity, ∮ω=∬Ω ↔ Berry holonomy, χ=1 ↔
CP-divisible (Markovian) limit, and the Bloch quarter-turn J²=−I echoes
the seam phase ι²=−1. A witness, not a general theorem; quantum gravity
is a horizon, not a claim.

## Open — the single remaining obligation

- **(U) Uniqueness** — ∬Ω is the *unique* obstruction to χ→1. Finite,
  gate-free, attackable by the same exact-rank machinery (H²(response)
  rank 1). Touches no RH / K0 / L0 / YM continuum gate.

## Reproduce

```
python papers/cut-first-equivalence/certificates/cfe1_cut_first_equivalence.py
python papers/cut-first-equivalence/certificates/cfe2_surjectivity.py
python papers/cut-first-equivalence/certificates/cfeq_quantum_witness.py
python -m pytest papers/cut-first-equivalence/tests -v
```

CI regenerates the certificate and fails on any pin drift
(`EXPECTED_CFE1.sha256`, `EXPECTED_CFE2.sha256`, `EXPECTED_CFEQ.sha256`).

## Relationship to the rest of the corpus

CFE-1 is the **parent** capsule. Its downstream corollaries live in
[`../lambda-seam-calibration/`](../lambda-seam-calibration/): the
seam-integer quantization (LAM-1/LAM-2) becomes "the obstruction takes
quantized seam values," and the "no memory, no critical strip" finding
(LAM-3) becomes "∬Ω is non-trivial iff the memory channel is retained."

## Capsule CFE-U: uniqueness — the flagship closed

`certificates/cfeu_uniqueness.py` discharges Part (U) with the same
exact-rank machinery that closed (S):

| Block | Certifies | Result |
|---|---|---|
| **T1** | The response complex Λ⁰→Λ¹→Λ² has exact ranks 5, 1 and d₁∘d₀=0; im d₀ **is** CFE-2's rank-5 classical space | PASS |
| **T2** | H¹ = 0: closed ⟺ gradient, constructively — the explicit potential integrates every closed form; zero curvature leaves no residual obstruction | PASS |
| **T3** | **The uniqueness theorem:** the annihilator of the memoryless space is EXACTLY 1-dimensional; the loop residue lies in it, is nonzero, hence spans — every obstruction is c·∮ω; and R = Area·d₁ on the whole basis | PASS |
| **T4** | The dial's direction 5p dv excites the generator; the residue ladder −32…+32 reproduced from the **pinned** CFE-1 module (two routes agree in ℚ); equal residue ⟹ equal obstruction class with explicit potential witness | PASS |
| **T5** | The 16-cell residue matrix has rank exactly 1 with kernel = the memoryless space: one obstruction repeated, not many | PASS |
| **T6** | Control: quadratic-coefficient algebra has rank 3 — uniqueness is a theorem about the **declared** algebra, and the boundary itself is certified | PASS |

**With CFE-1 + CFE-2 + CFE-Q + CFE-U, every declared obligation of the
Cut-First Equivalence theorem is machine-certified** on the witness EOS
and its declared response algebra.
