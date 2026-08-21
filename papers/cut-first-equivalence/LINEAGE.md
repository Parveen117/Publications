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
