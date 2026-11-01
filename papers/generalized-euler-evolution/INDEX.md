# Generalized Euler evolution: certified stages

| Stage | Result | Evidence |
| --- | --- | --- |
| [GE2: clock-free recognition generator](GE2_CLOCK_FREE_RECOGNITION_GENERATOR.md) | Rational finite arrows recover the phase derivation; symmetric retained records recover native heat through recognition-memory density and an intrinsic history ledger. Unequal-step errors and native-curvature clock normalization are explicit. | Six written results, 18 exact/adversarial tests; [audit](GE2_SOURCE_AUDIT.md), [result](certificates/GE2_RESULT.json). |
| [GE1: domains and stable evolution](README.md) | Native generator domains, stable phase/resolvent limits, counterexamples to unrestricted Euler/Taylor use and a nonuniform-gap control. | Six written results, 17 exact/adversarial tests. Frozen predecessor files remain unchanged. |

Runtime: **Python 3.12 only**. Run both packets read-only:

```bash
python3.12 -B papers/generalized-euler-evolution/certificates/ge1_evolution.py --check
python3.12 -B papers/generalized-euler-evolution/certificates/ge2_clock_free.py --check
python3.12 -B -m unittest discover -s papers/generalized-euler-evolution/tests -v
```

Certification means scoped written proofs and exact finite controls, not
formal proof-assistant or independent expert certification. The selected
native free heat is covered; physical clock/protocol selection, general
feedback and interacting/4D/Clay obligations remain open.
