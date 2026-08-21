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
