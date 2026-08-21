# The Cut-First Equivalence Theorem

*Candidate flagship statement of the recognition / cut-flow framework
(arXiv:2603.20773). This document states the theorem, marks precisely
which parts are certified and which are open, and maps every proof
obligation to where its scaffolding lives.*

---

## Statement (plain language)

> Recognition/cut response reproduces classical thermodynamic response
> **exactly**. The classical laws — Maxwell relations, Onsager
> reciprocity, the equilibrium invariant 𝓘 = 1 — are the **memoryless
> limit** (χ → 1, equivalently dω = 0). The deviation from that limit is
> **exactly** the loop residue
>
> ∮_∂D ω = ∬_D Ω,
>
> the integrated cut-loop curvature. Equilibrium is the flat face
> (Ω = 0); irreversibility is curvature (Ω ≠ 0); and the residue is the
> complete obstruction to flatness.

In one line, of the Noether type — *a bridge identity between two
independently computed quantities*:

> **Classical thermodynamics is the memoryless sector of recognition
> thermodynamics, and ∬Ω is the obstruction.**

## Statement (precise)

Let (𝓜, J) be a regular response manifold carrying a cut involution
J (J² = 1) and a memory dial χ, with entropy-weighted response
one-form ω = λ_p dp + λ_v dv and curvature two-form Ω = dω. Then:

1. **(Memoryless limit.)** In the limit χ → 1 (equivalently dω = 0), the
   cut-flow response coincides with classical thermodynamic response:
   every Maxwell relation holds, the Onsager block is symmetric
   (L_pv = L_vp), and the invariant 𝓘 = Γ_c · Γ_m = 1.

2. **(Residue.)** For any 2-cell D, the circulation of ω around ∂D
   equals the curvature flux through D: ∮_∂D ω = ∬_D Ω, and this residue
   vanishes iff the response is memoryless on D.

3. **(Completeness — the two open halves.)**
   - **(S) Surjectivity.** Every classical response identity is realized
     by some cut protocol: the cut grammar maps *onto* the classical
     response algebra.
   - **(U) Uniqueness.** ∬Ω is the *unique* obstruction to the χ → 1
     limit: no second, independent obstruction exists.

Parts (1) and (2) are the **certified core** (see below). Parts (S) and
(U) are **open** and constitute the theorem's remaining proof
obligation.

---

## What is certified: capsule CFE-1

The capsule `certificates/cfe1_cut_first_equivalence.py` proves, in
**exact rational arithmetic** on an explicit witness equation of state,
four blocks (pin `EXPECTED_CFE1.sha256`, all tampers separate):

| Block | Certifies | Result |
|---|---|---|
| **T1** | Memoryless limit is classical: at χ=1 the curvature vanishes on the whole grid, every Maxwell loop sums to exactly 0, the Onsager cross-asymmetry is exactly 0 | PASS |
| **T2** | The residue is exact: ∮_∂D ω = ∬_D Ω as an identity in ℚ (boundary walk and area sum computed independently), witnessed at χ=7/5 with both sides = 32 | PASS |
| **T3** | The obstruction is faithful: Ω = 0 on the grid **iff** χ = 1; the residue is strictly monotone in the memory dial μ = χ−1 (−32, −16, 0, +16, +32 across χ = 3/5 … 7/5), zero exactly at μ=0, sign tracking sign(μ) | PASS |
| **T4** | The equilibrium invariant: 𝓘 = Γ_c·Γ_m = 1 exactly at χ=1, departing exactly when χ≠1 | PASS |

This maps Part (1) → T1+T4 and Part (2) → T2+T3. The residue's exact
linearity in the memory dial (T3) is the machine witness that ∬Ω is a
*genuine* order parameter for irreversibility, not an artefact of the
model.

## Why this is the right flagship

Discourse-changing theorems are **bridge identities**: two quantities
computed by different means that must be equal. Noether bridges symmetry
↔ conserved current; Einstein bridges geometry ↔ stress-energy; Onsager
bridges microscopic reversibility ↔ macroscopic transport symmetry;
Landauer bridges information erasure ↔ heat. Cut-First Equivalence
bridges **recognition/memory structure ↔ thermodynamic response
geometry**, with a computable residue ∬Ω. It upgrades the second law's
*inequality* slot to an *identity with a computable residue* — on the
witness EOS, certified.

It is also the **parent** statement. Once established:
- the **recognition-Noether** direction (bias phase φ as the generator
  of a conserved recognition current) becomes "the conserved current of
  the memoryless symmetry";
- the **memory-completion** finding (no memory, no critical strip;
  certified at s=3,5 in `../lambda-seam-calibration/`) becomes "the
  obstruction ∬Ω is non-trivial iff the memory channel is retained";
- the **seam-quantization** results (|τ(χ)|²=q, τ(χ₄)²=−4 ⇒ ι²=−1, in
  `../lambda-seam-calibration/`) become "the obstruction takes quantized
  seam values."

One theorem organizes four, instead of four competing claims.

## The remaining obligation (open, precisely stated)

To promote the certified core to the full theorem, prove — **within the
regular / certified regime, touching no RH / K0 / L0 / YM continuum
gate**:

- **(S) Surjectivity / completeness.** That the cut grammar reproduces
  *every* classical response identity, not only the Maxwell / Onsager /
  invariant content verified on the witness EOS. This is the true
  novelty half: it says recognition thermodynamics *loses nothing*
  classical.
- **(U) Uniqueness.** That ∬Ω is the *only* obstruction to χ → 1 — that
  no independent second obstruction can arise. This is the harder half.

Both are finite, well-posed, gate-free obligations. Neither requires the
analytic-continuation machinery behind the framework's open number-theory
gates.

## Proof-scaffolding map

| Obligation | Scaffolding already in the corpus |
|---|---|
| Part (1), memoryless limit | this capsule T1/T4; paper Corollary .1.1 (χ=1 ⇒ Onsager), contact-form Maxwell derivation, 𝓘=1 |
| Part (2), residue = curvature | this capsule T2/T3; paper Stokes identity ∮ω=∬Ω (Eqs. for Ω=dω) |
| Curvature = micro-correlator | paper Weighted Reciprocity (Eq. 22) + Microscopic Irreversibility theorem |
| (S) surjectivity | **to build** — the cut grammar's image in the classical response algebra |
| (U) uniqueness | **to build** — obstruction-theoretic argument that H²(response) is rank 1 |
| quantized residue values (downstream) | `../lambda-seam-calibration/` LAM-1/LAM-2 (seam integers, Gauss-sum cuts) |

## Honest boundary

The witness equation of state is an explicit rational model chosen so the
arithmetic is exact. **The theorem is a claim about the structure**
(one-form, curvature two-form, memory dial), **not about this one EOS** —
the EOS is the witness, not the theorem. The certified core establishes
the memoryless-limit direction and the residue identity; it does **not**
establish (S) or (U), and it makes no claim about the framework's open
number-theory or Yang–Mills continuum gates. What the framework's genuine
originality rests on is the **bridge to finite seam arithmetic** developed
in the neighbouring Λ-seam capsules; this document's flagship should be
read together with those.
