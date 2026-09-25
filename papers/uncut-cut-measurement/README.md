# Physics from Cuts: Uncut Ground, Measurement and Lawful Continuation

**Monty Dabas — working research programme, version 0.1, 25 September 2026.**

The starting proposal is reflexive: definitions create distinctions, so the
mathematical framework used to describe the uncut is itself a cut. A measurement
must therefore retain its declared observational context and unresolved
information. Spacetime is a possible derived representation of the programme;
it is not an input to the finite results here.

Read [the manuscript](MANUSCRIPT.md), then the [claim ledger](CLAIMS.md).
The [research programme](RESEARCH_PROGRAMME.md) specifies the next derivations.

## What this first version establishes

Within an explicitly declared finite presentation, the manuscript proves:

1. A target can be recovered from a cut exactly when it is constant on every
   readout fibre. A deterministic correction cannot create missing distinctions.
2. The minimum extra memory alphabet is the maximum number of distinct target
   values hidden behind one readout. This is a target-relative cost.
3. Refinement reduces target ambiguity; sequential repair obeys a product bound,
   which need not be an equality.
4. A clock-free transition descends to observed dynamics exactly when equivalent
   observations have equivalent next observations.
5. Repeated distinction by future readouts terminates on a finite presentation
   and constructs the coarsest dynamically closed refinement of a cut.
6. Common scalar invariants may be trivial even when several cuts jointly
   distinguish all presented states. Pairwise consistency need not give a global
   state compatible with all observations.

The inherited RSC cut-corner composition identity is restated and checked
separately as a downstream linear representation. These are elementary finite
derivations and a connection to the existing native framework, not a claim of
priority for the underlying mathematical facts.

## Reproduce

Python 3.11 or 3.12; standard library only. From the repository root:

```bash
python -m unittest discover -s papers/uncut-cut-measurement/tests -v
python papers/uncut-cut-measurement/certificate.py --check
```

The second command is read-only. It regenerates exact exhaustive checks in
memory and compares both the committed certificate bytes and SHA-256 pin.
It also binds the local manuscript, implementation, tests, source manifest,
claim ledger and research programme. An intentional revision uses
`certificate.py --write`, followed by review of the changed source and pin.

`CERTIFICATE.json` records **PASS_FINITE_CHECKS** with the enumerated domains.
The written proofs cover the hypotheses stated in the manuscript; Python
enumeration is not formal verification of those proofs. Source pins identify
the inspected repository versions, not an independent certification of them.

## Scope of the physical programme

Uncut ground is a proposed primitive, not identified with the finite carrier.
The carrier and its functions are declared mathematical representations.
No spacetime metric, external clock, probability law, quantum state space,
Hamiltonian or classical field equation is assumed by these finite propositions.
No physical quantum/classical correspondence, spacetime reconstruction, minimum
action law, universal curvature-information identity or quantum-gravity closure
has yet been derived by this package.

This addition is independent of pending [correction PR #4](https://github.com/Parveen117/Publications/pull/4).
It uses pinned native sources and does not treat archived publication claims as
new premises. See [SOURCE_PINS.json](SOURCE_PINS.json).
