# Publications

Public release repository for citable research outputs by Monty Dabas.

## Latest certified result — real quantum hardware (Aug 29, 2026)

**Cross-platform typed-verdict certification, executed from this
repository's CI at zero cost** —
[`papers/quantum-certified-verdicts/AZQ1/`](papers/quantum-certified-verdicts/AZQ1/):
the same preregistered audit ran on a Quantinuum H2 trapped-ion
emulator, Rigetti QVM, and a **real 156-qubit IBM Heron QPU
(ibm_kingston)**. On the real machine the auditor **refused** to
certify an uncorrected GHZ state (z = +40.8), **certified** the
error-suppressed one (z = −2.8) under the identical unchanged
contract, and cut the deliberately corrupted circuit both times
(z = +73.4 / +66.6). Five pinned certificates, four preregistrations,
single-shot no-retry discipline throughout. The instrument measures
execution quality — it does not rubber-stamp it.

Sibling certified programmes: battery-health falsification ladder with
an independent-evaluator harness in
[Energy](https://github.com/Parveen117/Energy) (flagship: preregistered
aging-trend certification on six CALCE cells, combined p = 3.1e-11).

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
