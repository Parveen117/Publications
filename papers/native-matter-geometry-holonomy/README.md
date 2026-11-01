# MGH-1 — Matter at an observation cut: geometric holonomy and quantum correlations

14 October 2026. Continuation of N4H-1, CRC-1 and QGC-1.

This packet realizes the proposed matter–cut coupling using the actual
holonomy of two four-dimensional simplicial metrics. The configurations
share a six-cell triangulation and hinge area. One is flat at the selected
hinge; the other has an exact nonzero deficit. A declared quantum register
represents these complete metric configurations, including their boundary
data. The original native arithmetic and proof engines are unchanged.

For a local unit-mass occupancy cut, pure matter input with occupancy
probability `p`, flat initial geometry and the specified coherent
interaction, the exact information–curvature relation is

\[
\boxed{S_L=p(1-p)\kappa_h,\qquad
\kappa_h=1-\cos\delta_h,\qquad
S_L=1-\operatorname{Tr}\rho_m^2.}
\]

The geometric phase fixes the half-angle in the preceding chat model.
For the exact fixture, `exp(i delta_h)=(-7+24 i)/25` and
`exp(i delta_h/2)=(3+4 i)/5`. With `p=1/2`:

| Quantity | Exact result |
|---|---:|
| Hinge defect `kappa_h` of the curved configuration | `32/25` |
| Joint state purity | `1` |
| Either reduced-state purity | `17/25` |
| Entanglement linear entropy | `8/25` |
| Final curved-configuration probability | `8/25` |
| Mean compact hinge defect | `256/625` |
| Five-setting entanglement witness | `-2/5` |
| White-noise visibility sufficient for this witness | `v > 5/13` |

The geometry is computed from positive four-simplex Gram matrices and
whole-tetrahedron unfolding. Quantum configuration coherence, its mixing
operation, and the mass-weighted phase coupling are explicit model choices.
The result is an exact finite matter–geometry model and a reproducible
quantum demo. Its physical gravitational coupling, continuum limit and
relation to the earlier Regge stationary equations remain to be derived.
The calculations are predictions, not hardware or gravitational observations.

- [THEOREM.md](THEOREM.md): constructions, proofs, assumptions and lineage.
- [DEMO.md](DEMO.md): register meaning, operation, exact outputs, witness
  measurements and rejection controls.
- [CERTIFICATE.json](CERTIFICATE.json): exact geometric and quantum witnesses.
- [SOURCE_PINS.json](SOURCE_PINS.json): engine and predecessor provenance.

From the repository root:

```bash
python3.12 papers/native-matter-geometry-holonomy/run.py --check
```

Python 3.12 performs I/O only. The pinned native JavaScript engines perform
all scientific arithmetic and polynomial proof replay. The read-only check
recomputes the certificate bytes and checksum. Analytic exponential/angle
arguments and general positivity statements have written proofs; this is
not a proof-assistant formalization of all real analysis.

The construction reuses established controlled-unitary and entanglement
witness methods. The specific contribution is the exact geometric
realization, the resulting conditional information–hinge relation, and
the tests showing what fails when source coherence, higher moments or
the curvature lift is omitted.
