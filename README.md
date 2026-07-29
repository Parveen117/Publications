# Publications

Public release repository for citable research outputs by Monty Dabas.

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
