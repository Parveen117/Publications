# Lineage — Cut-First Equivalence

## Origin

CFE-1 formalizes the candidate flagship theorem identified for
arXiv:2603.20773 ("Geometric Completion of Thermodynamic Response").
The paper already proves the memoryless-limit direction as corollaries
(χ=1 ⇒ Onsager; the contact-form origin of Maxwell relations; 𝓘=1) and
states the Stokes residue ∮ω=∬Ω. CFE-1 lifts these from prose and
scattered lemmas into a single certified capsule with an explicit
witness equation of state, so the flagship rests on machine-checkable
arithmetic.

## What is new here vs. the paper

Nothing mathematically new is *claimed*. What is added is:
- an explicit rational witness EOS on which the residue is computed two
  independent ways (boundary walk, area sum) and shown exactly equal;
- the faithfulness result T3 (residue strictly monotone in the memory
  dial, zero iff equilibrium) which makes ∬Ω a certified order
  parameter, not merely a defined quantity;
- an explicit, pinned statement of the two open obligations (S)
  surjectivity and (U) uniqueness that separate the certified core from
  the full theorem.

## Relationship to prior capsules

- The classical substrate (Stokes, contact-form Maxwell, Onsager in the
  reversible limit) is a **pinned named dependency**, cited and not
  rederived — the same discipline used in the Λ-seam and Yang–Mills
  capsules.
- The downstream quantization of the residue lives in
  `../lambda-seam-calibration/`: LAM-1/LAM-2's seam integers and
  Gauss-sum cuts are the arithmetic side of the bridge, and LAM-3's
  "no memory, no critical strip" is the memory-channel necessity. CFE-1
  is the thermodynamic parent under which those become corollaries.

## Claim boundary

Certified: the memoryless-limit direction and the residue identity on an
explicit EOS, in exact arithmetic. Open: (S), (U), and generality beyond
the witness EOS. Untouched: RH / K0 / L0 / YM continuum gates. The EOS is
the witness, not the theorem.

## CFE-2 addendum

CFE-2 discharges obligation (S), surjectivity, on the witness EOS: the
cut grammar generates the entire classical response algebra. The cut
basis is built by cut-native centered finite differences of the
potential — a procedure independent of the analytic partials used for
the classical frame — and the two span the same exact 5-dimensional
space (cokernel 0). The renaming test (TAUT-1) was applied and passed:
random potentials of the same shape also reach full rank (structural,
not point-luck), and a non-closed shear control separates. Only (U)
uniqueness remains open for the full theorem; no RH/K0/L0/YM continuum
gate is touched.

## CFE-Q addendum

CFE-Q is step (2) of the quantum sequence: a certified quantum-
thermodynamics witness on a dissipative qubit. It shows the Cut-First
Equivalence structure realizes quantum-mechanically — Markovian
(CP-divisible) limit as the flat χ=1 face, closed-cycle geometric
holonomy as a faithful non-Markovianity measure, and the Bloch
quarter-turn J²=−I as the Berry-phase echo of the Λ-seam ι²=−1. It is a
witness, not a general open-system theorem; generality and the full
theorem are open. Quantum gravity is deliberately NOT addressed — it
remains a horizon, to be attempted only by connecting a certified piece
to an established bridge (e.g. Jacobson 1995), never asserted as
foundation. No RH/K0/L0/YM continuum gate is touched.
