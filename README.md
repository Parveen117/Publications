# Publications

Citable research outputs of **Monty Dabas / Celextrix Pvt Ltd**. One
research foundation — the recognition-kernel theorem ladder — feeds
several application programmes; this repository holds the papers and
certified evidence common to all of them. **If you arrived here from a
specific application, jump straight to your programme below.**

## Corrections and version scope

[Corrections dated 2026-09-25](papers/corrections-2026-09-25/README.md) supersede
identified statements in T50/T54, SPECTRAL 2 Corollary 16.2, information
invariance, Morphic manuscripts and the thermodynamic Pluecker test.
The [corrected thermo source edition](papers/thermodynamic-response-corrections/)
includes the complete source and exact provenance. Archived PDFs and the
collected volume retain their original edition identity; read them with the
correction record. A finite certificate PASS is distinct from a general proof,
formal verification, empirical support and independent review.

## Programme index — find your lane

| If you are evaluating... | Go to | Headline |
| --- | --- | --- |
| **RNKE-Q / quantum verification** (SINE, IIT Bombay) | [`papers/quantum-certified-verdicts/`](papers/quantum-certified-verdicts/) + [Quantum-Classical-public](https://github.com/Parveen117/Quantum-Classical-public) | Real 156-qubit IBM Heron QPU: auditor certified an error-suppressed GHZ state (z = -2.8), refused the uncorrected one (z = +40.8), cut corrupted circuits every time; cross-platform concordance with fault localization at z = +12.2; and a seven-qubit multi-ledger operator-memory benchmark with **CROSS_BACKEND_REPLICATION certified on three real IBM QPUs** (ibm_fez, ibm_kingston, ibm_marrakesh — 3/3 runs, all declared gates) |
| **Battery health certification** (IITM Pravartak) | [Energy repository](https://github.com/Parveen117/Energy) | Preregistered aging-trend certification on six CALCE cells, combined p = 3.1e-11, with an Ed25519 independent-evaluator harness |
| **ATHENA navigation / SETU communications** (FITT, TIDES) | Patent estate + product repos | These tracks are patent- and product-led (PCT/IB2025/060887, PCT/IB2026/051695, PCT/IB2026/058465); this repository is their shared research foundation |
| **The mathematics itself** | [`book/recognition-kernel-collected-volume/`](book/recognition-kernel-collected-volume/) + [Recognition-Kernel-Framework](https://github.com/Parveen117/Recognition-Kernel-Framework) | 600+ page collected theorem volume; 60+ certified theorems with stated claim boundaries |
| **Thermodynamic response papers** | [`papers/`](papers/) + [Thermodynamics-Reproducibility](https://github.com/Parveen117/Thermodynamics-Reproducibility) | arXiv:2603.20773 and the T01 audit discipline every programme follows |

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

Claim boundary reported by the latest source (not a fresh independent
revalidation of its full numerical certificate):

```text
SOURCE FLOOR S_(.01,-) >= .008 I               PROVED
SOURCE INVERSE NORM <= 125                      PROVED
RANK-AT-MOST-FIVE HERMITIAN REDUCTION           PROVED
FINAL OUTWARD ACCEPTANCE THEOREM                PROVED
ACTUAL FIXED-SHIFT FIVE-BY-FIVE VALUE           REPORTED CERTIFIED BY LATEST MP COMPANION
UNSHIFTED ENDPOINT REDUCTION                    PROVED AS REDUCTION IN LATEST SOURCE
ACTUAL ENDPOINT SIGN K0 >= 0                    OPEN
RIEMANN HYPOTHESIS                              ABSTAIN
```

The source manuscript is developed in `Parveen117/MP` on branch
`agent/paper-source-floor-five-by-five-theorem` and draft PR #244. The source
version used for this status reconciliation is commit
`686b5506c7c2fedd53641dbf550218b60c996a13`, especially section 07l and
A_claim_ledger. It reports the fixed-shift bound 0.8388901637561339 < 1;
that bound alone does not establish the unshifted K0 sign or RH.

## Zenodo release rule

Create a GitHub release only after:

1. the final source tree is copied here;
2. the clean compiled PDF is included;
3. metadata and claim boundaries match the PDF;
4. the repository is public;
5. the repository is enabled in Zenodo's GitHub integration.

Each GitHub release is intended to be archived by Zenodo and assigned a version DOI.
