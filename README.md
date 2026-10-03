# Publications

Latest Yang–Mills continuation: [YM-54](papers/yang-mills-certified-benchmark/YM54_RESPONSE_CURVATURE_PROTOCOL_BRIDGE.md)
constructs a counted compact heat protocol from two native shape-response
directions. Its second tensor eigenvalue obeys beta>=4 f^2/tau, where f
is the existing response-curvature marker and tau is the shape-Gram trace.
Uniform |f|>=f0>0 and tau<=T feed YM53's interacting window
abs(theta)<f0^2/(600 T). The original curved thermo fixture gives an
exact example. Seven written results, twenty new and 171 related tests
use **Python 3.12 only**; 110 upstream sources retain their bytes.
Equal signed counts, frozen response jets, directions and heat clock are
declared choices. Physical selection, anisotropic joint/volume limits and
4D/Clay remain open.

Previous interacting continuation: [YM-53](papers/yang-mills-certified-benchmark/YM53_ANISOTROPIC_INTERACTING_GAP.md)
extends the interacting chain to native anisotropic and rank-two heat.
If every site's second tensor eigenvalue is at least beta>0,
abs(theta)<beta/2400 suffices for a gap uniform in finite width and
fine heat step. The proof permits unequal, differently oriented tensors,
handles inadmissible short bridges, and carries the gap into each
fixed-width time limit. Interacting energy uses the derived ground-source
weight and obeys the same lower bound. Seven written results, eighteen
new exact-control tests and 151 related tests use **Python 3.12 only**.
The anisotropic joint/volume limit, physical selection and 4D/Clay remain open.

Previous energy continuation: [YM-52](papers/yang-mills-certified-benchmark/YM52_ENERGY_NOISE_AND_OBSERVER_GAP.md)
connects the declared native heat protocol to derivative energy, discarded
record variance and bounded-density entropy dissipation. For tensor
eigenvalues a<=b<=c, its exact centered rate is min(tr(C)/4,a+b);
two noncommuting active directions suffice on this compact carrier.
A sign-insensitive observer has rate a+b. At fixed trace, optimizing
that even-observer rate selects isotropy with an explicit stability bound,
whereas optimizing the full-observer rate leaves an anisotropic plateau.
The known spectral formula is credited to Lauret; the native record
bridge and observer selection conditions are explicit.
Six written results and sixteen new tests use **Python 3.12 only**.
Physical energy/time calibration, selection of the optimization objective
remain open; YM53 develops the declared anisotropic chain extension.

Previous tool audit: [YM-51](papers/yang-mills-certified-benchmark/YM51_TOOL_DEPENDENCY_AND_HEAT_SELECTION.md)
audits 29 tool dependencies and identifies the heat protocol's symmetric
second moment as a remaining selection choice. The same native brackets,
positive reference and all linear-coefficient decay can coexist with
arbitrarily slow quadratic decay. Curvature alone does not choose the heat
law or its gap. Fixed-law internal invariance selects isotropy up to a
scale; covariance alone does not. Six quadratic channels recover the
general symmetric protocol tensor.
Seventeen new and 117 related tests pass. YM-50's native reference bridge
and the earlier fixed-protocol chain results remain unchanged.
Physical isotropy/clock/action selection, general UGD integration, the actual
row defect, NCG quantum measure and 4D/Clay remain open.
Current YM verification uses **Python 3.12 only**.

[YM-43](papers/yang-mills-certified-benchmark/YM43_NATIVE_TIME_TRANSFER_BOUND.md)
supplies the preceding width-uniform time bound at three coarse heat-kernel cells.
[YM-42](papers/yang-mills-certified-benchmark/YM42_SILENCE_TIME_LADDER.md)
supplies the preceding structural silence ball and spatial memory-tail result.

New mathematical development: [Native compact gauge completion](papers/native-compact-gauge/README.md)
derives the action-compatible U(2)/SU(2) sector from independent EMK factors,
with native curvature, source and Gauss identities and scoped exact controls.

For common mathematical development, start with [Morphic/EMK/UGD mathematical foundation](MATHEMATICAL_FOUNDATION.md) and [EMK-C1: native connection and curvature calculus](papers/emk-ugd-algebra/CONNECTION_CALCULUS.md). Mathematical chapters develop here; current physical applications continue in [extra-ideas](https://github.com/Parveen117/extra-ideas). The shared engine is the single canonical RKF source linked below.

[Native thermodynamic curvature, NT-1 to NT-8](papers/native-thermodynamic-curvature/README.md) now connects the thermodynamic response paper and theorem register to that calculus: an exact EMK-to-Hessian intertwiner, ordered Smriti balance, observer-return curvature and the scoped Onsager-flatness criterion. The native edition includes written proofs, a compiled paper and independent exact controls.

Citable research outputs of **Monty Dabas / Celextrix Pvt Ltd**. One
research foundation — the recognition-kernel theorem ladder — feeds
several application programmes; this repository holds the papers and
certified evidence common to all of them. **If you arrived here from a
specific application, jump straight to your programme below.**

## Programme index — find your lane

| If you are evaluating... | Go to | Headline |
| --- | --- | --- |
| **RNKE-Q / quantum verification** (SINE, IIT Bombay) | [`papers/quantum-certified-verdicts/`](papers/quantum-certified-verdicts/) + [Quantum-Classical-public](https://github.com/Parveen117/Quantum-Classical-public) | Real 156-qubit IBM Heron QPU: auditor certified an error-suppressed GHZ state (z = -2.8), refused the uncorrected one (z = +40.8), cut corrupted circuits every time; cross-platform concordance with fault localization at z = +12.2; and a seven-qubit multi-ledger operator-memory benchmark with **CROSS_BACKEND_REPLICATION certified on three real IBM QPUs** (ibm_fez, ibm_kingston, ibm_marrakesh — 3/3 runs, all declared gates) |
| **Battery health certification** (IITM Pravartak) | [Energy repository](https://github.com/Parveen117/Energy) | Preregistered aging-trend certification on six CALCE cells, combined p = 3.1e-11, with an Ed25519 independent-evaluator harness |
| **ATHENA navigation / SETU communications** (FITT, TIDES) | Patent estate + product repos | These tracks are patent- and product-led (PCT/IB2025/060887, PCT/IB2026/051695, PCT/IB2026/058465); this repository is their shared research foundation |
| **The mathematics itself** | [`book/recognition-kernel-collected-volume/`](book/recognition-kernel-collected-volume/) + [Recognition-Kernel-Framework](https://github.com/Parveen117/Recognition-Kernel-Framework) | 600+ page collected theorem volume; 60+ certified theorems with stated claim boundaries |
| **Thermodynamic response papers** | [Native thermodynamic curvature](papers/native-thermodynamic-curvature/README.md) + [Thermodynamics-Reproducibility](https://github.com/Parveen117/Thermodynamics-Reproducibility) | arXiv:2603.20773 core converted into native response/curvature language, with an explicit source map and finite certificate |

**One discipline everywhere:** preregistration before data, SHA-256
pins on every certificate, failed predictions published with
mechanisms named, and CI that re-runs any result on demand.

## Publication layout

```text
papers/
  representation-complete-native-seam-integer/
    source/        LaTeX source tree
    pdf/           final compiled PDF
    metadata/      arXiv and release metadata
```

## Current paper

**A Representation-Complete Reduction of the Riemann Hypothesis to a Native Seam Integer: Bilateral Recognition, Spectral Blindness, and a Finite Birman--Schwinger Obstruction**

Author: Monty Dabas  
ORCID: 0009-0005-6948-209X

Claim boundary:

```text
SOURCE FLOOR S_(.01,-) >= .008 I               PROVED
SOURCE INVERSE NORM <= 125                      PROVED
RANK-AT-MOST-FIVE HERMITIAN REDUCTION           PROVED
FINAL OUTWARD ACCEPTANCE THEOREM                PROVED
ACTUAL OUTWARD FIVE-BY-FIVE VALUE               OPEN
CONTROLLED eta_j -> 0                           OPEN
RIEMANN HYPOTHESIS                              ABSTAIN
```

The source manuscript is developed in `Parveen117/MP` on branch
`agent/paper-source-floor-five-by-five-theorem` and draft PR #244.

## Zenodo release rule

Create a GitHub release only after:

1. the final source tree is copied here;
2. the clean compiled PDF is included;
3. metadata and claim boundaries match the PDF;
4. the repository is public;
5. the repository is enabled in Zenodo's GitHub integration.

Each GitHub release is intended to be archived by Zenodo and assigned a version DOI.
