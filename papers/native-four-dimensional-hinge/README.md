# N4H-1 — Four-dimensional hinge dynamics

14 October 2026. Continuation of the native cut/operator gravity program.

The next step after NJC-1 is now concrete: complete simplex metrics fix
the previously free normal transport, a four-cell model has a genuine
interior-edge equation, and the curvature phase exposes a missing angle
lift in compact holonomy observations.

Read [THEOREM.md](THEOREM.md) for assumptions, derivations and scope.

## Results

* Unique whole-tetrahedron gluing:
  `T = F (E^dagger E)^(-1) E^dagger - n_target n_source^dagger`.
* Two nondegenerate six-cell stars have identical compact holonomy but
  lifted deficits `0` and `-2*pi`. Complete metric data distinguish them.
* Four four-simplices with one interior squared edge `x` give the unique
  vacuum solution `x=4` within `0<x<16/3`, under the selected Regge action
  and specified fixed boundary data.
* The exact curvature-functional Hessian is `R''(4)=-1/2`.
  A small information source gives
  `x-4 = 2 (alpha/gamma) I'(4) + O((alpha/gamma)^2)` under the stated
  smoothness conditions. This is a stationary response, not a stability claim.
* A phase `exp(iota*gamma*a*delta)` descends to an unlifted SO(2) record
  for all angle lifts iff `gamma*a` is integral. Retaining the lift avoids
  that restriction. This is a conditional quantization criterion.

All variables and couplings in the equations are dimensionless. No
external time, physical measurement, path-integral measure or fundamental
Hilbert-space axiom is introduced. Euclidean metric geometry and the
action are explicit additional contracts. This is an intermediate step
toward quantum gravity, not its completion.

## Reproduce

From the repository root, with Node.js and Python 3.12:

```sh
python3.12 papers/native-four-dimensional-hinge/run.py --check
```

The launcher uses Python only for file handling, hash validation and
invocation. It extracts the two unchanged native JS engines from the
existing `papers/Native_Phase_Word_Reconstruction_Source.zip` into a
temporary directory. All exact mathematical operations run in those
engines. No numerical or symbolic mathematics package is required.

`--check` verifies source pins, reconstructs every certificate and result,
and compares the full payload byte for byte with [CERTIFICATE.json](CERTIFICATE.json)
and [EXPECTED.sha256](EXPECTED.sha256). To regenerate deliberately after
updating source pins, use `--write`.

The certificate reports algebraic witnesses and exact finite controls.
Its symbolic proof/replay checks are included in the exact-check count;
these counts must not be added as independent theorems. General real-angle
and analytic statements retain their written hypotheses.

## Next research gates

1. Couple several interior metric variables and solve the carrier equations
   with controlled boundary data; test nontrivial curvature and constraints.
2. Construct the Lorentzian and quantum amplitude/measure contracts, choosing
   explicitly whether curvature lifts or a specified cover are retained.
3. Establish refinement/continuum behavior, gravitational constraint closure,
   and measurable predictions. The present finite certificate supplies none
   of these by extrapolation.

Only this new packet was added for this development. The native engines,
uploaded manuscripts and previous certificate files are unchanged.
