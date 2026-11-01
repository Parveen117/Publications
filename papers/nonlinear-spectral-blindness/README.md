# Nonlinear spectral blindness and response flatness

NSB1 develops the spectral-observation programme in nonlinear dynamics
and constitutive physics, then returns to the native gauge model.

The strongest physics example consists of two strictly convex energies
with **the same Hessian eigenvalues at every state in a common region**.
One response metric is flat; the other has positive curvature, equal
to **1/16** at the stated reference state. A separate electromagnetic
example has the same signed wave speeds at every bias but different
calibrated curvature. Impedance supplies missing information there.

The nonlinear theorem proves that every finite jet-pullback
characteristic polynomial, including mixed unmarked spectral words,
can omit all nonlinear coefficients. Coordinate response markers
recover those coefficients exactly on the declared finite jet model.
An explicit quartic mechanical system shows why this distinction
matters for stability and energy recovery.

The T–V–S–P connection is precise: a strictly convex potential gives
two dual flat affine connections, while their metric-compatible mean
can be curved. Arbitrary physical symmetry does not imply that every
curvature vanishes. Mixed-control process loops and metric transport
loops are also distinct.

## Read the development

- [NSB1 manuscript: proofs and exact examples](NSB1_NONLINEAR_SPECTRAL_BLINDNESS.md).
- [Source audit and correction to SPECTRAL_2 Section 16](SOURCE_AUDIT.md).
- [Yang–Mills dock and the remaining vacuum estimate](YANG_MILLS_DOCK.md).
- [Exact result](certificates/NSB1_RESULT.json),
  [source manifest](certificates/NSB1_SOURCE_MANIFEST.json), and
  [canonical hash](certificates/EXPECTED_NSB1.sha256).

The source PDFs remain unchanged. Their general negative-compression
bound is not inherited: the audit records a counterexample and the
corrected `max(0, mu_n) + tail` bound.

## Reproduce

Python **3.12 only**, standard library, with assertions enabled:

```bash
python3.12 -B papers/nonlinear-spectral-blindness/certificates/nsb1_nonlinear_blindness.py --check
python3.12 -B -m unittest discover -s papers/nonlinear-spectral-blindness/tests -v
```

There are **48 new tests**: 26 nonlinear jet and mechanical tests,
10 independently differentiated geometry tests, and 12 native-gauge,
source-correction and evidence tests. The dedicated workflow also
replays NSR1, NSR2 and YM82. Normal reproduction is read-only;
`--write` regenerates only the NSB1 result and its hash after an
intentional reviewed change.

The written proofs carry the universal statements. Finite exact
controls and independent tests check formulas, examples, evidence
integrity and claim boundaries; they are not formal proofs of every
analytic assertion.

## Scope of the gauge consequence

The native single-plaquette construction matches arbitrarily many
local potential derivatives of the free operator while its physical
gauge-invariant gap can approach zero. It therefore rules out a
positive uniform gap conclusion from those finite local jets alone.
This is a separately specified parent family. It does not have an
identical full Hamiltonian spectrum, and does not prove a gap for
the original homogeneous four-dimensional continuum Yang–Mills model.

Carleman triangularity and dual-flat Hessian geometry are established
inputs, acknowledged in the manuscript. No exhaustive external
priority claim is made.
