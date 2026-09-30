# UGD coefficient geometry and a variational propagation sector

Author: Monty Dabas. Research development: R9, vacuum speed and physical calibration.

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

**R3: metric dynamics.** [MG-1–MG-8](NATIVE_METRIC_DYNAMICS.md) classify a
declared native polynomial trace-action class into Palatini–Holst, cosmological
and bulk-inactive terms. For nondegenerate coframes and nonzero Palatini
coupling, connection/coframe variations give vacuum Einstein equations for
every real Holst coefficient. Sixteen stable thermo channels explicitly supply
all coframe variations and all ten metric variations. Adding the declared
Maxwell law with independent gauge variations yields Einstein–Maxwell equations.
The continuum, action class, constitutive adapter and coupling values remain
inputs. This uses established first-order gravity, with primary references and
explicit coefficient, rank and failure controls.

**R4: native loop curvature.** [LR-1–LR-8](NATIVE_LOOP_CURVATURE.md) derive a
curvature-square continuum action from a declared oriented quadratic loop
readout. Its compact matrix stationary law is D_A K(F)=0; for the unprojected
readout it reduces to D_A{J,F}=0. The even/odd equations recover MG's torsion
and Einstein equations. The coefficient relations are Lambda=-12 sigma u^2
and kappa Lambda=6/beta. The original real-four module permits the negative-
Lambda sector; one extra native R factor with R^2=-I supplies a real positive-
Lambda sector. Explicit de Sitter/anti-de Sitter patches have flat combined
Cartan transport and nonzero metric curvature. Exact readout-symmetry and
identical-return scale controls keep primitive selection and physical
calibration distinct from this conditional construction. The established
MacDowell–Mansouri/Cartan mechanism is credited.

**R5: massive classical matter.** [SM-1–SM-8](NATIVE_MASSIVE_MATTER.md)
repair NP's massive action obstruction by doubling the fixed real module.
The additional native return gives an invertible mass-compatible variational
multiplier. A declared covariant minimal action supplies the field equation,
phase current, coframe stress and spin source, including the off-shell
torsion-trace term. In LR's zero-Holst bulk sector the full connection equation
has a unique algebraic torsion solution. Eliminating it gives an explicit
quartic contact action and matching cubic matter equation. This is a native
coefficient realization of the established Einstein–Cartan mechanism, with
all 24 connection directions checked. Two universal commuting-field Fierz
identities and a failed cosmological-mass shortcut make the scope precise:
particle masses, quantum statistics and a primitive action selection remain
open.

**R6: an explicit coupled solution.** [CS-1–CS-8](COUPLED_COSMOLOGICAL_SOLUTION.md)
solve the neutral homogeneous rest sector of R5 with an FLRW metric and
nonzero polarized torsion. The derived source has rho=m n+3 kappa n^2/16,
p=3 kappa n^2/16 and n=n_0/a^3. A closed hyperbolic volume law and logarithmic
native spinor phase solve the full first-order equations in LR's finite-scale
positive-Lambda sector. An exact benchmark has a^3=exp(t)-1. The lapse
constraint, all coframe/connection/matter components and effective stress are
checked independently. In the massless member the full connection scalar is
constant even while a torsion invariant diverges at the past endpoint. The
model therefore supplies a solved family, not an automatic bounce, quantum
cosmology or measured parameter prediction. General Cauchy existence remains
open; the R7 extension below addresses two specified linear stability problems.

**R7: perturbation bounds.** [PS-1–PS-8](PERTURBATION_STABILITY.md) extend the
homogeneous equations to all eight matter components and prove bounded linear
perturbations of the coupled expanding positive-Lambda FLRW sector. A separate
linear matter problem on the prescribed CS geometry has a wavelength-independent
Sobolev energy bound. For the benchmark starting at t*=log(2), its rescaled
future norm gain is at most sqrt(2). Spatial gradients leave the rest sector,
so all eight components are retained. The massless zero-Lambda limit supplies
an explicit logarithmic relative-growth control. R8 below adds homogeneous
anisotropic metric modes; full inhomogeneous Einstein–matter and general
nonlinear stability remain open.

**R8: exact anisotropic cosmology.** [BI-1–BI-8](ANISOTROPIC_COSMOLOGY.md)
retain a general Bianchi I spatial metric and derive its spin-shear stress
pi=(Z/8)[sigma,Q]. The volume-rescaled shear evolves by orthogonal
conjugation, so its eigenvalues are conserved. An explicit nonlinear rest
family solves the volume, native phase, full coframe and torsion equations;
physical shear decays as 1/v on its expanding positive-Lambda branch.
Linearization supplies all five homogeneous shear modes alongside the R7
matter/scalar modes. A diagonal-metric restriction requires spin alignment;
equal volume and selected curvature scalars can hide different anisotropic
stress. The fixed parallel-frame contact phase/rotation ratio is 3/2 within
this action. These results extend the source-bound model and credit classical
Einstein–Dirac precedents; they do not select a physical action or constants.

**R9: vacuum speed and calibration.** [SC-1–SC-8](SPEED_AND_CALIBRATION.md)
derive the physical characteristic kernels of the existing common-metric
vacuum action: two photon and two classical metric polarizations share its
null cone. Thus c_GW/c_gamma=1 in the geometric-optics limit, consistent with
the specified published GW170817 bound. This is an inherited Einstein–Maxwell
prediction, not a new discriminator against GR or a primitive-only selection.
Native matter shares the wavefront while a massive mode has a smaller group
speed. An explicit admissible speed family and a different-metric control
expose the unselected constitutive/clock/length inputs. The SI value of c is
defined exactly, so choosing units to reproduce it is not an independent
prediction. [Reference data](SPEED_REFERENCE_DATA.json) retain the external
definition and the observation's emission assumptions separately from the
model. No quantum graviton or numerical dimensional constant is derived.

| Result | Location | Scope |
|---|---|---|
| Native complex sector, transport compatibility, corrected H/S derivatives, Hessian lift and determinant–Einstein equivalence | [Recovered bridge](history/RKF_UGD_Kahler_Einstein_Bridge.md) | Original note preserved byte for byte; its historical manuscript references are identified as such. |
| Even-index signature obstruction; hyperbolicity criterion; smooth Monge–Ampère benchmark with convex Dirichlet uniqueness and an inverse variational functional | [Geometry and signature](GEOMETRY_AND_SIGNATURE.md), GS-1–GS-4 | Written deductions with explicit hypotheses; not a general existence theorem. |
| An explicit Lorentzian Einstein product and a failed gradient-reflection control | [Geometry and signature](GEOMETRY_AND_SIGNATURE.md), GS-5 | A standard product geometry with a freely chosen length scale. |
| Clifford symbol, native variational multiplier, action, density drift, and the precise lower-order compatibility obstruction | [Propagation and action](NATIVE_PROPAGATION_ACTION.md), NP-1–NP-6 | A selected coefficient-sector field law, not selection of nature's dynamics. |
| Positive response cone, native source current, explicit thermo connection chart and full gauge variations | [Thermo gauge completion](THERMO_GAUGE_COMPLETION.md), CP-1–CP-8 | Coupled equations follow from the declared gauge kinetic law and fixed backgrounds; no metric dynamics or physical constants are selected. |
| Native curvature-action classification, full thermo metric variations, conditional Einstein and Einstein–Maxwell equations | [Metric dynamics](NATIVE_METRIC_DYNAMICS.md), MG-1–MG-8 | Explicit polynomial action class, independent Lorentz connection and a local constitutive adapter; not primitive-only selection or a prediction of couplings. |
| Quadratic native loop limit, compact matrix field equation, readout classification and cosmological sign/scale relations | [Loop curvature](NATIVE_LOOP_CURVATURE.md), LR-1–LR-8 | Declared local loop process, module and readout; explicit curved vacua and a scale ambiguity, not measured constants or a unique primitive process. |
| Doubled massive action, covariant matter sources and full algebraic torsion elimination | [Massive matter](NATIVE_MASSIVE_MATTER.md), SM-1–SM-8 | Declared commuting classical field and minimal coupling in the zero-Holst sector; Einstein–Cartan reduction with the contact term retained, not quantum fermions or selected physical masses. |
| Analytic coupled spinor/coframe/torsion family with full field-equation controls | [Coupled cosmological solution](COUPLED_COSMOLOGICAL_SOLUTION.md), CS-1–CS-8 | Neutral homogeneous rest sector with an isotropic metric and polarized torsion; exact positive-Lambda solution and scalar-readout counterexample, not general PDE existence or empirical cosmology. |
| Coupled homogeneous linear stability, gradient leakage, and an exact spatial matter gain bound | [Perturbation stability](PERTURBATION_STABILITY.md), PS-1–PS-8 | All eight homogeneous matter components retained; spatial estimate holds on prescribed geometry. BI below adds homogeneous shear; full inhomogeneous gravity and general nonlinear stability remain open. |
| Spin-shear stress, exact nondiagonal coframe, nonlinear rest-family isotropization and five homogeneous shear modes | [Anisotropic cosmology](ANISOTROPIC_COSMOLOGY.md), BI-1–BI-8 | General Bianchi I spatial metric; nonlinear solution restricted to homogeneous rest matter, with homogeneous linear completion around CS. Spatially varying gravity and general nonlinear stability remain open. |
| Vacuum characteristic quotient, common speed ratio and calibration nonselection | [Speed and calibration](SPEED_AND_CALIBRATION.md), SC-1–SC-8 | Existing common-metric Einstein vacuum; inherited ratio 1 is consistent with a specified published bound. No primitive-only SI speed, new distinction from GR or quantum graviton is derived. |
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

CP-5–CP-8 supply the gauge variation interface and MG-2–MG-5 the metric
interface. LR supplies an explicit loop-process candidate, its continuum
action and the resulting coefficient relations; SM supplies a conditional
massive matter action and its current, stress and torsion derivation. CS now
constructs an analytic self-consistent family of these coupled fields. PS
adds its coupled homogeneous linear stability and a prescribed-geometry
spatial matter bound. BI extends the solved family to rotating Bianchi I
shear and completes the homogeneous anisotropic linear modes. What remains
is selection of the module, readout, adapter, matter couplings,
inverse-length scale and action normalization by an independently specified
native process principle. Complete return records alone retain the LR-8 scale
ambiguity, and the naive positive-Lambda Cartan coupling does not produce SM's
mass term. SC now proves the conditional vacuum speed ratio and identifies
the independent constitutive/calibration input needed for an operational
native speed relation. G, c in SI, Λ, α, hbar and particle masses remain
unpredicted from primitive laws; c's SI numeral is itself a definition. Quantum
statistics, the generic spatially varying Einstein–matter Cauchy problem,
full spatially inhomogeneous stability and nonlinear stability for general
matter data remain open. Sixteen thermo
channels are sufficient here, not proved minimal.
