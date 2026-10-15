# Quantum gravity: canonical research record

Updated 15 October 2026. Current stage: **LCF-1**.

## Working agreement

The research belongs in this repository: derivations, code, source pins, certificates, rejected stronger claims and the next exact obligation. Chat reports summarize repository work; standalone documents do not substitute for it. Continue autonomously from the strongest verified endpoint, reuse previous results, and avoid parallel duplicate research. The goal is a rigorous quantum-gravity construction, not a growing count of loosely connected examples.

Preserve frozen predecessor packets. Distinguish an algebraic proof, a written analytic theorem, a selected model and a physical identification. No certificate count by itself closes a continuum or physical claim. Use Python 3.12 for orchestration and the existing native arithmetic/rewrite engines for the scientific fixtures.

## Current result

[LCF-1](papers/native-local-coframe-covariance/README.md) addresses AGD's actual mixed local obstruction. RCG/CID's metric coefficient yields a coframe density on regular spatial cut patches; for the circle, `e_cut=|partial_x theta|/2`. Its derived transformation `delta_f e=(f e)'` cancels the mixed residual in `A_N(e)=C D_N+gRNe`. EMK-C1's full covariance contract is used directly. The UGD-A0 carrier restriction and UGD-G1 phase/scale/sheet distinctions are retained rather than replacing native associative operator arithmetic by raw carry addition.

The enlarged spatial/spatial and spatial/normal brackets close exactly. The normal/normal bracket still differs by the coframe-field direction `V_(NM'-MN')`. This identifies the next missing term, rather than declaring a complete constraint algebra. A spatially covariant matter action supplies `J=i eta g psi*R psi`; on an oriented circle-cut patch its fixed-N variation gives matter force `partial_x(NJ)/2` on the cut momentum.

The local cut-field interpretation is additional structure: a spatially constant cut has zero sensitivity metric and cannot supply this nondegenerate coframe. AGD's global quantum geometry is not silently promoted to a local quantum field. Its established dynamics, IRT reflection invariant and ZCR recovery estimates remain frozen and scoped to their original models.

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
| [IRT-1](papers/native-invariant-ratio-transport/README.md) | Coupled reflection-weight ratio, classification of the common geometry invariant algebra, and calibrated rate/diagonal theorem | Specified invariant algebra is not a complete state record or local gravitational constraint algebra |
| [CLD-1](papers/native-cut-local-dynamics/README.md) | Changing-cut memory, finite covariant closure/conditioning budget, and AGD moving-reflection criterion | Cut-relative locality is not spatial locality; infinite memory norm estimates and local gravitational constraints are not proved |
| [RCG-1](papers/native-relational-cut-space/README.md) | Based cut geometry, rank-one overlap/cycle reconstruction, common-transport residual and its exact geometry-change law | Standard reconstruction method credited; operational geometry does not select physical spacetime or an interacting gravitational action |
| [LCF-1](papers/native-local-coframe-covariance/README.md) | Cut-gradient coframe density, cancellation of mixed local residual, covariant matter action and local cut-source variation | Normal geometry transport, full constraints and quantum local-field completion remain open; regular spatial patch required |

The earlier classical MG/SM action results are pinned and included as source snapshots in CCL-1. Reuse their existing classical coframe/source variation instead of proposing it as new. They do not supply the missing positive interacting quantum theory.

## Next proof target

**Derive normal evolution of the local geometry from an admitted action so that it supplies LCF's missing `V_(NM'-MN')`, and prove preservation of the full constraints.**

LCF resolves the specified mixed covariance defect; it does not yet provide the normal geometry generator. The next work must:

1. Specify the local cut-field canonical/quantum data and a geometry action, retaining the derived fixed-N matter source and distinguishing independent coframes from composite cut coframes.
2. Compute all brackets with the resulting geometry-field normal action. Remove or consistently reinterpret the explicit normal obstruction without dropping an observable field direction.
3. Retain the GCP coframe/matter interface and RCG comparison/phase records. Treat degenerate cut-gradient patches and global UGD sheet records explicitly where they enter.
4. Establish domain, positivity and controlled classical-field recovery for the enlarged model; the prior global AGD/ZCR results do not automatically prove these local-field claims.

A complete physical claim still requires the intended 3+1 interacting limit, physically meaningful observables and calibrated predictions. Work on the specific missing normal-evolution theorem rather than reporting a completion percentage.

## Replay status

The eight integrated predecessor packets (NJC-1, angular calibration, EGR-1, EGR-2, CRC-1, QGC-1, N4H-1, MGH-1) remain frozen. CCL-1 has 69 exact fixture checks and 3 controls. GCP-1 has 115 exact checks, 6 native word proofs/replays and 9 controls. AGD-1 adds 92 exact checks, 8 native word proofs/replays and 8 controls. ZCR-1 adds 169 exact checks, 11 native word proofs/replays and 13 controls. IRT-1 adds 103 exact checks, 9 native word proofs/replays and 16 controls. CLD-1 adds 167 exact checks, 15 native word proofs/replays and 14 controls. RCG-1 adds 108 exact checks, 9 native word proofs/replays and 10 controls. LCF-1 adds 95 exact checks, 8 native word proofs/replays and 17 controls. The repository runner replays all 16 stages read-only; CI separately replays the source CID-D2 action packet. Analytic proofs and their assumptions are documented beside the finite certificates.

```sh
python3 scripts/replay_quantum_geometry.py
```
