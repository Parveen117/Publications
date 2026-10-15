# Quantum gravity: canonical research record

Updated 15 October 2026. Current stage: **ZCR-1**.

## Working agreement

The research belongs in this repository: derivations, code, source pins, certificates, rejected stronger claims and the next exact obligation. Chat reports summarize repository work; standalone documents do not substitute for it. Continue autonomously from the strongest verified endpoint, reuse previous results, and avoid parallel duplicate research. The goal is a rigorous quantum-gravity construction, not a growing count of loosely connected examples.

Preserve frozen predecessor packets. Distinguish an algebraic proof, a written analytic theorem, a selected model and a physical identification. No certificate count by itself closes a continuum or physical claim. Use Python 3.12 for orchestration and the existing native arithmetic/rewrite engines for the scientific fixtures.

## Current result

[ZCR-1](papers/native-zero-reference-recovery/README.md) makes the proposed fixed-reference recovery condition quantitative. It reuses GCP's exact invariant scalar-sector criterion, supplies finite-carrier leakage/readout bounds, and proves finite-time flat-metric recovery for the actual spatially varying AGD operator. An explicit normalizable circle preparation has exact moments; scaling its width and dimensionless action parameter together makes the matter reading approach the flat-metric law while geometry continues evolving.

This closes a **controlled flat-metric recovery estimate within the selected AGD model**. It does not establish classical matter trajectories, a unique quantum measurement outcome, physical action-scale calibration, or local gravitational constraint closure. At fixed action scale, indefinite narrowing raises the preparation's kinetic energy. A change of the reference called zero does not remove the retained-record defect.

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
| [GCP-1](papers/native-geometry-controlled-propagation/README.md) | Geometry controls derivative, connection and curvature; reciprocal coupling and exact reduction obstruction | Geometry mixing and profiles selected |
| [AGD-1](papers/native-action-geometry-dynamics/README.md) | Inherited cut action fixes an infinite geometry generator under an explicit quantization contract; joint existence, source law and finite-section bound | Global circle geometry; local gravitational constraints, physical action selection and lapse profile remain open |
| [ZCR-1](papers/native-zero-reference-recovery/README.md) | Exact packet moments, finite observation-sector error and variable-profile flat-metric recovery | Joint preparation/action-scale limit; matter remains quantum; general classical gravitational field limit open |

The earlier classical MG/SM action results are pinned and included as source snapshots in CCL-1. Reuse their existing classical coframe/source variation instead of proposing it as new. They do not supply the missing positive interacting quantum theory.

## Next proof target

**Derive a local dynamical coframe/source action and determine whether its constraints account for the residual `g R N f'` in AGD-1's naive local algebra.**

The inherited geometry generator is now available. The next packet must address the actual missing local structure:

1. Specify local geometry fields and their canonical data, rather than treating the single global circle coordinate as all gravitational degrees of freedom.
2. Derive their transformation laws and constraints from a declared action. Show that the coupled dynamics preserves those constraints; neither the solved global clock constraint nor the uncoupled parity symmetry supplies this result.
3. Retain the GCP matter–metric interface, the geometry source response and common spatial cone, or state precisely why their forms must change.
4. Preserve or extend ZCR's controlled flat-metric recovery to a general curved evolving classical field. ZCR controls matter propagation within the selected model; it does not yet recover classical gravitational dynamics. AGD's homogeneous Fourier-truncation estimate separately controls computation.

A complete physical claim additionally requires the intended 3+1 interacting limit, physically meaningful observables, stability/positivity on the physical sector and calibrated predictions. Work on one concrete missing theorem at a time; do not replace these gates by a completion percentage.

## Replay status

The eight integrated predecessor packets (NJC-1, angular calibration, EGR-1, EGR-2, CRC-1, QGC-1, N4H-1, MGH-1) remain frozen. CCL-1 has 69 exact fixture checks and 3 controls. GCP-1 has 115 exact checks, 6 native word proofs/replays and 9 controls. AGD-1 adds 92 exact checks, 8 native word proofs/replays and 8 controls. ZCR-1 adds 169 exact checks, 11 native word proofs/replays and 13 controls. The repository runner replays all 12 stages read-only; CI separately replays the source CID-D2 action packet. Analytic proofs and their assumptions are documented beside the finite certificates.

```sh
python3 scripts/replay_quantum_geometry.py
```
