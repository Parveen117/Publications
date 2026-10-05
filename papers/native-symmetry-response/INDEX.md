# Native symmetry-response developments

| Stage | Question | Main result |
| --- | --- | --- |
| [NSR2: conjugate-diagram properties](NSR2_CONJUGATE_DIAGRAM_PROPERTIES.md) | Do lambda_P, lambda_V and the other T–V–S–P properties have direct counterparts in other diagrams? | Yes: a six-property dictionary transfers with the two conjugate pairs and four Legendre potentials, retaining constraints, signs and units. |
| [NSR1: electromagnetic response curvature](README.md) | What information lives in the higher response after a constitutive law is specified? | Coupling and impedance derivatives reconstruct a selected native response curvature. |

NSR2 is the direct property transfer. Its [diagram comparison](figures/conjugate_diagrams.svg)
shows the thermodynamic, electromagnetic and charge–flux versions side by
side. NSR1 develops an additional geometric diagnostic on such response
families. Both stages retain the existing native thermo connection and
make their physical adapters explicit.

## Reproduce the direct transfer

```bash
python3.12 -B papers/native-symmetry-response/certificates/nsr2_diagram_properties.py --check
python3.12 -B -m unittest discover -s papers/native-symmetry-response/tests -v
```

The calculations require only the standard library. The optional figure
builder uses Matplotlib under Python 3.12. Written proofs, exact finite
checks, and empirical validation are distinct evidence levels.
