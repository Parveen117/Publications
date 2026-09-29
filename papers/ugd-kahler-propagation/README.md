# UGD coefficient geometry and a variational propagation sector

Author: Monty Dabas. Research development: R2, 30 September 2026.

This package preserves the earlier λ-Hessian/Kähler–Einstein development and
extends it with an explicit action for a selected native-coefficient propagation
law. It belongs to Publications PR #5 on
`research/uncut-cut-measurement-2026-09-25`; inclusion on this branch is not a
claim of a merged paper, peer review or an archived release.

The new result is precise: two admitted EMK coefficient blocks realize a real
Clifford propagation system. Its native volume element supplies the unique
constant skew variational multiplier, up to scale, for the stated four-component
real field. A weighted conservation law fixes the symmetric lower-order drift.
The resulting quadratic action has exactly that first-order equation as its
Euler–Lagrange equation. Choosing the propagation law, continuum, density and
constitutive coefficients remains input. No Einstein–Hilbert action or measured
physical constants have been derived from the native axioms alone.

**R2: thermodynamic gauge completion.** [CP-1–CP-8](THERMO_GAUGE_COMPLETION.md)
derive the Lorentz cone of positive native response operators, turn the field's
phase symmetry into its own conserved source current, and give an explicit
thermodynamic potential whose phase connection has invertible Darboux
coordinates. Varying the thermo state maps recovers every sourced Maxwell
equation on rank-four patches of the declared coupled action. Four curved
thermo channels are necessary and sufficient for this regular representation
to include zero curvature. A two-channel counterexample demonstrates why
field representation alone does not establish variational completeness.

| Result | Location | Scope |
|---|---|---|
| Native complex sector, transport compatibility, corrected H/S derivatives, Hessian lift and determinant–Einstein equivalence | [Recovered bridge](history/RKF_UGD_Kahler_Einstein_Bridge.md) | Original note preserved byte for byte; its historical manuscript references are identified as such. |
| Even-index signature obstruction; hyperbolicity criterion; smooth Monge–Ampère benchmark with convex Dirichlet uniqueness and an inverse variational functional | [Geometry and signature](GEOMETRY_AND_SIGNATURE.md), GS-1–GS-4 | Written deductions with explicit hypotheses; not a general existence theorem. |
| An explicit Lorentzian Einstein product and a failed gradient-reflection control | [Geometry and signature](GEOMETRY_AND_SIGNATURE.md), GS-5 | A standard product geometry with a freely chosen length scale. |
| Clifford symbol, native variational multiplier, action, density drift, and the precise lower-order compatibility obstruction | [Propagation and action](NATIVE_PROPAGATION_ACTION.md), NP-1–NP-6 | A selected coefficient-sector field law, not selection of nature's dynamics. |
| Positive response cone, native source current, explicit thermo connection chart and full gauge variations | [Thermo gauge completion](THERMO_GAUGE_COMPLETION.md), CP-1–CP-8 | Coupled equations follow from the declared gauge kinetic law and fixed backgrounds; no metric dynamics or physical constants are selected. |
| Exact reproducible controls and source hashes | [Certificate](CERTIFICATE.json), [verifier](verify.py), [source pins](SOURCE_PINS.json) | Finite rational checks support the written proofs; they are not a proof assistant or empirical validation. |

The current thermo programme starts with a stable fundamental relation U(S,V)
and its declared Hessian response metric; see [TC-1–TC-6](../thermo-compass-foundations/THEOREM.md)
and [PF-1–PF-7](../thermo-phase-field/THEOREM.md). The older H/S Hessian discussed
in the recovered note is a distinct construction. This package does not equate
their potentials, state spaces or metrics with physical spacetime.

The shared native operator engine remains in
[Recognition-Kernel-Framework/operator_foundation](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/main/operator_foundation).
`model.py` is a publication-specific exact fixture implementation, not another
native arithmetic engine. The finite coefficient representation used here does
not faithfully represent a proper one-sided seam TS=I, ST≠I.

## Reproduce

From the repository root, with Python 3.11 or 3.12 and only the standard library:

```bash
python papers/ugd-kahler-propagation/verify.py --check
```

`--check` is read-only and is the default. Explicit `--write` regenerates the
certificate and digest after a reviewed source change. The verifier also reruns
the recovered bridge script. CI checks this package alongside the existing
thermo compass and phase-field packages.

## Next unresolved derivation

The potential-to-gauge-connection map and its admissible gauge variations now
have a concrete realization in CP-5–CP-8. Derive or select the continuum
propagation coefficients and their own dynamics from an explicit native
process law. In particular, a λ-to-coframe map, a rule
for admissible metric variations, and the coefficient of a curvature action
are still missing. A prescribed Lorentzian principal symbol and a matter-field
action do not supply those ingredients. G, c, Λ and α remain unpredicted.
