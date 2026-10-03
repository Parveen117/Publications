# Thermo phase field: electromagnetic bridge R1

Monty Dabas · 29 September 2026

This develops the [thermo-compass foundations](../thermo-compass-foundations/README.md) into an explicit Abelian phase-field model. The [PF-1–PF-7 proof](THEOREM.md) constructs the phase connection, identifies the required cut contract, completes the local curvature representation with two thermo channels, and derives a discrete field equation with energy and propagation controls.

| Result | Status |
|---|---|
| Explicit scalar connection A=(db-(b/a)da)/(2 sqrt(det H)) | Derived from the declared thermo response metric; curvature agrees with TC-5. |
| Continuous phase and readout compatibility | The phase partner must survive the readout; a fixed reflection stabilizer is discrete, and a parallel reflection forces zero curvature on that line. |
| One versus two thermo channels | A single equilibrium pullback has rank at most two. Two independent phase channels represent every closed nondegenerate two-form locally in four dimensions, by Darboux's theorem. |
| Maxwell-form field equations | Derived from a declared quadratic link-phase action on an oriented cell complex; includes source continuity and both divergence constraints. |
| Positive field energy | Exact conservation without sources and exact source-work balance; rational midpoint transport preserves the energy and constraints. |
| Propagation | Two transverse modes per nonzero Fourier character; model speed follows from selected response weights and clock/ruler calibration. |
| Rest energy and inertial mass | E0=m_in c^2 follows for a separate gapped mode under the stated relativistic dispersion and energy/momentum calibration. Gap and calibration are not derived. |

The information connection is explicit: a calibrated response catalogue supplies a Gram weight, and its recognized/discarded parts add exactly. Choosing the catalogue's electric/magnetic interpretation, the cell complex, source units and process law remains part of the model contract. Local curvature completeness does not itself select the physical electromagnetic field.

The two-channel construction gives one tensor-product U(1) connection. It does not claim two photon species. A finite lattice has dispersive corrections; exact continuum Lorentz symmetry, numerical c, hbar, alpha and particle masses are not claimed. The standard geometric and discrete-electrodynamics identities are credited in the proof.

## Reproduce

Python 3.11 or 3.12, standard library only, from repository root:

```bash
python papers/thermo-phase-field/verify.py --check
python papers/thermo-compass-foundations/verify.py --check
```

Default verification is read-only. For intentional evidence regeneration after source review:

```bash
python papers/thermo-phase-field/verify.py --write
python papers/thermo-phase-field/verify.py --check
```

The [certificate](CERTIFICATE.json) records independent scalar/full-Christoffel curvature comparisons, compact native-phase holonomies, boundary-of-boundary identities, action variations, source work, charge continuity, exact finite transport, Fourier polarization ranks, negative controls and input refusals. The general arguments are the written proofs; the sample counts do not prove physical identification.

[Source pins](SOURCE_PINS.json) identify the existing results consumed. [Earlier physics inventory](../thermo-compass-foundations/RESULTS_INDEX.md) locates the rest of the programme.

## Next development target

Select the response catalogue and reversible process law from the native transport/recognition rules, then compare the resulting propagation and source coupling across independent sectors. The completed R1 model provides exact quantities to constrain that selection: phase curvature, two transverse modes, conserved charge, positive energy and the dispersion relation. Matter-gap selection follows as a separate problem.
