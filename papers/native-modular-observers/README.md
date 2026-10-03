# Native modular observers — NM1

Author: Monty Dabas. Development edition: 3 October 2026.

This folder starts a mathematical route from the existing associative EMK
operator algebra to modular symmetry and observer-dependent loss of information.
It works in Publications only; no physics repository is updated by this stage.

**Result:** a chosen rational lattice adapter realizes SL(2,Z) in the native
algebra. Its projective observer hides one central sign. An exact sign ledger
restores operator composition, and a proof shows that no choice of representatives
can remove the correction everywhere. This sign loss persists under all common
left/right continuations in that observer; coarser trace loss need not persist.

| Written result | Scope |
| --- | --- |
| NM1-T1 | Native rotation/reflection, derived shear, explicit lattice, Euclidean modular generation |
| NM1-T2 | Minimal two-valued lift, sign cocycle, associativity and nonsplitting |
| NM1-T3 | Stable projective blindness versus recoverable trace blindness; endpoint versus history |
| NM1-T4 | Automorphy factors, integer-weight parity and central-sign obstruction |
| NM1-T5 | Functional transformation defects and the ambiguity of a completion |

Read [THEOREMS.md](THEOREMS.md) for proofs and
[SOURCE_AUDIT.md](SOURCE_AUDIT.md) for lineage and source contracts.
[The certificate](certificates/NM1_RESULT.json) provides exact finite controls;
[the pins](certificates/NM1_SOURCE_PINS.json) bind inputs and unchanged sources.
Written general proofs and executable controls are separate evidence. Neither
formal proof-assistant verification nor independent expert certification is claimed.

## Reproduce

From the repository root, with Python **3.12 only** and no third-party packages:

```sh
python3.12 -B papers/native-modular-observers/certificates/nm1_modular_observers.py --check
python3.12 -B -m unittest discover -s papers/native-modular-observers/tests -v
```

Default execution and `--check` are read-only. `--write` regenerates only this
stage's result and digest. The narrow CI workflow runs both checks and verifies
that evidence and consumed source files remain unchanged.

## What this does and does not unify

Rotation, reflection and shear coexist in the same admitted algebra. A circle's
periodic endpoint does not retain its entire path. The native matrix sign,
projective action, integer winding and a coarser scalar readout are distinct
targets. None is silently renamed connection curvature.

This is a concrete bridge to classical modular algebra, not a claim to have
invented modular groups, cocycles, covering spaces or modular-form theory.
The lattice and observer are declared choices. No canonical physical selection,
new mock theta function, shadow, nonholomorphic completion, RH theorem or
Yang–Mills mass gap is derived here.

The next research gate is specific: select a native history generating function
and an admitted analytic function space, derive its transformation defect, then
prove existence and normalization of an appropriate correction. A sign ledger
alone does not supply a half-integral-weight multiplier or a mock-modular shadow.
