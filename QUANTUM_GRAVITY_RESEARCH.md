# Quantum gravity: canonical research record

Updated 15 October 2026. Current stage: **GCP-1**.

## Working agreement

The research belongs in this repository: derivations, code, source pins, certificates, rejected stronger claims and the next exact obligation. Chat reports summarize repository work; standalone documents do not substitute for it. Continue autonomously from the strongest verified endpoint, reuse previous results, and avoid parallel duplicate research. The goal is a rigorous quantum-gravity construction, not a growing count of loosely connected examples.

Preserve frozen predecessor packets. Distinguish an algebraic proof, a written analytic theorem, a selected model and a physical identification. No certificate count by itself closes a continuum or physical claim. Use Python 3.12 for orchestration and the existing native arithmetic/rewrite engines for the scientific fixtures.

## Current result

[GCP-1](papers/native-geometry-controlled-propagation/README.md) places an independent quantum geometry observable in matter's derivative coefficient. A bounded, genuinely curved Lorentzian branch family supplies the same propagation cone, connection and operator curvature. A joint linear law gives reciprocal matter–geometry coupling. Exact results identify the geometry variance lost by replacing the joint model with a single mean metric, and the invariant subspace condition under which that replacement becomes exact.

This closes the **construction of a causal quantum geometry-dependent propagation model**. The geometry profiles and mixing law remain selected. It does not yet close gravitational dynamical selection, refoliation constraints or the physical interacting continuum.

## Lineage to reuse

| Stage | Established content | Boundary to preserve |
|---|---|---|
| [CID-D2](papers/native-cut-transport-feedback/README.md) | Native cut and transport feedback from a selected conservative action | Finite carrier; not quantum gravitational constraints |
| [NJC-1](papers/native-jacobian-cut-certification/README.md) | Native Jacobian/cut/connection identities and retained free connection | A cut does not fix all ambient connection data |
| [EGR-1](papers/native-exponential-gravity-response/README.md) / [EGR-2](papers/native-nonlinear-gravity-continuation/README.md) | Regge stationary response and local nonlinear branch | Selected action, boundary and local hypotheses |
| [CRC-1](papers/native-cut-record-completeness/README.md) | Operation-dependent record completeness and source-law controls | A short invariant list need not close an arbitrary operation |
| [QGC-1](papers/native-quantum-gravity-consistency/README.md) | Ensemble obstruction and consistent quantum adapters | Old nonlinear relaxation is not an autonomous quantum channel |
| [N4H-1](papers/native-four-dimensional-hinge/README.md) | Actual four-dimensional Euclidean hinge metrics and lifted curvature | Lorentzian quantum dynamics not obtained by relabelling |
| [MGH-1](papers/native-matter-geometry-holonomy/README.md) | Independent metric register, geometric phase and exact quantum witness | Phase interaction and metric mixing selected |
| [CCL-1](papers/native-common-cone-limit/README.md) | Shared-refinement cone, local mixing and a conditional curved symbol | Shared scaling and physical dictionary explicit |
| [GCP-1](papers/native-geometry-controlled-propagation/README.md) | Geometry controls derivative, connection and curvature; reciprocal coupling and exact reduction obstruction | Geometry action and constraints remain open |

The earlier classical MG/SM action results are pinned and included as source snapshots in CCL-1. Reuse their existing classical coframe/source variation instead of proposing it as new. They do not supply the missing positive interacting quantum theory.

## Next proof target

**Derive a compatible geometry generator from the source/action and close its constraints while retaining the GCP matter–metric interface.**

The next packet should specify the following before claiming a dynamical completion:

1. A geometry carrier and a source functional/action that determine the geometry transition law, including its free connection and boundary data.
2. One joint state-independent quantum operation or generator, with a declared physical constraint subspace and proof that the dynamics preserves it. Auxiliary commuting constraints must not be called the gravitational constraint algebra.
3. A controlled regime in which the geometry readout recovers the existing classical field equations, with a quantified error. The two-metric demo alone is not this limit.
4. A clear account of locality and common geometry for all tested matter sectors; retain the return term when eliminating quantum geometry.

A complete physical claim additionally requires the intended 3+1 interacting limit, physically meaningful observables, stability/positivity on the physical sector and calibrated predictions. Work on one concrete missing theorem at a time; do not replace these gates by a completion percentage.

## Replay status

At this update all eight newly integrated predecessor packets (NJC-1, angular calibration, EGR-1, EGR-2, CRC-1, QGC-1, N4H-1, MGH-1) replayed unchanged. CCL-1's 69 exact fixture checks and 3 controls pass. GCP-1 has 115 exact checks, 6 native word proofs/replays and 9 controls. CI replays the source-bound chain using the repository runner below.

```sh
python3 scripts/replay_quantum_geometry.py
```
