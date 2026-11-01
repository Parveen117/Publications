# Native symmetry response: electromagnetic properties and higher susceptibilities

Monty Dabas. 5 October 2026. **NSR1, Python 3.12 only.**

The thermodynamic response construction transfers to a reversible electromagnetic
constitutive law. This folder develops six useful consequences, with written
proofs, exact examples, and an explicit route to reconstructing the response
geometry from wave measurements. It also treats general linear symmetries and
equilibrated internal modes.

The main electromagnetic result is the oriented response-curvature two-form

\[
\boxed{f=-\frac{dr\wedge d\log Z}{2(1-r^2)^{3/2}},\qquad
r=\frac{b}{\sqrt{ac}},\quad Z=\sqrt{a/c},\quad
H=\begin{pmatrix}a&b\\b&c\end{pmatrix}>0.}
\]

Here H is the energy Hessian in declared electric-displacement and magnetic-flux
coordinates. Z is normalized differential wave impedance in the closed
one-polarization sector specified in the paper. The formula describes the
**selected native response connection**, not electromagnetic gauge curvature or
an automatically accumulated optical phase.

| Result | Useful consequence | Status |
| --- | --- | --- |
| NSR1-T1 | Releasing a conjugate constraint enhances both responses by the same factor 1/(1-r²) | Standard Schur identity, explicitly transferred from the thermo construction |
| NSR1-T2 | Coupling and impedance gradients determine oriented response curvature; a common stiffness multiplier cancels | Native electromagnetic synthesis with a written proof |
| NSR1-T3 | Two signed wave speeds and differential impedance reconstruct the response matrix; bias derivatives reconstruct its curvature | Exact result for the declared closed wave sector |
| NSR1-T4 | Symmetry selects the allowed higher responses; parity can hide curvature at zero bias while rotation symmetry permits it away from zero | General covariance identities and exact examples |
| NSR1-T5 | Born–Infeld has nonzero response curvature on an explicit stable patch, although a quadratic law can match all point linear responses | Worked use of an established nonlinear electromagnetic model |
| NSR1-T6 | Relaxing an internal mode can create response curvature absent in the clamped family | Exact equilibrium Schur reduction and a local stable example |

These are additions to this research programme. No literature-priority claim
or new experimentally established electromagnetic law is made. Symmetry
constrains response functions; their coefficients still require a constitutive
law or measurement.

- [Derivations and proofs](THEOREMS.md)
- [Measurement protocol and useful next targets](MEASUREMENT.md)
- [Source map and novelty boundary](SOURCE_MAP.md)
- [Exact executable evidence](certificates/NSR1_RESULT.json)
- [Pinned dependencies](certificates/NSR1_SOURCE_PINS.json)

## Reproduce

From the Publications root:

```bash
python3.12 -B papers/native-symmetry-response/certificates/nsr1_response.py --check
python3.12 -B -m unittest discover -s papers/native-symmetry-response/tests -v
```

Only the Python standard library is required. The verifier reuses the unchanged
NT independent Christoffel implementation, compares it with the native
connection's represented commutator, and checks rational Taylor coefficients
through fourth order. It does not install a new native algebra engine.
`--check` is read-only; `--write` regenerates this stage's result and digest
after source verification. General statements have written proofs; the finite
checks test their implementations and failure controls.

The Yang–Mills continuum problem is not advanced by an asserted mass-gap claim
here. These response diagnostics may be applied to a *separately established*
gauge-response family; no such identification is silently supplied.
