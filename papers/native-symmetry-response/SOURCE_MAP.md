# NSR1 source and contribution map

This is an additive folder in Publications, based on the YM69 publication
tree `7f8162638e9ac75e9add4157c640135ac9a4f79a`, published in commit
`7bc141c21aa409b7a4718593603a2348730c00b7`. Earlier chapters are preserved.

| Source | Consumed result | Role in NSR1 |
| --- | --- | --- |
| `native-thermodynamic-curvature/THEOREM.md`, NT-1 | positive native response element, typed derivatives | electromagnetic constitutive adapter with its own units and signs |
| Same, NT-3/4 | A=H⁻¹dH/2, F=-[X_i,X_j]/4, oriented phase form | T2 derives the coupling/impedance representation |
| Same, NT-8 | covariant response tower and its slot actions | T4 retains the native tower while deriving symmetry selection rules |
| `native-thermodynamic-curvature/controls.py` | independent full Christoffel calculation | reused unchanged to check rational Taylor fixtures |
| `emk-ugd-algebra/CONNECTION_CALCULUS.md` | native connection, covariance and derivative contract | foundational dependency; not a newly implemented algebra |
| `native-thermodynamic-curvature/SOURCE_MAP.md` | original thermodynamic source order and limitations | distinguishes unit-typed thermo response from electromagnetic quantities |
| Mezincescu–Russo–Townsend, arXiv:2311.04278v3 | Hamiltonian Maxwell equations and Born–Infeld energy | established input for T3/T5, with new readout calculations derived here |
| Gibbons–Rasheed, hep-th/9506035 | nonlinear electromagnetic duality context | prevents treating duality in nonlinear electrodynamics as a new discovery |

The historical six-response summary was also read in Publications at
`f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56`, path
`papers/thermodynamic-response-corrections/source/sections/01a_response_summary.tex`.
It fixes the earlier temperature/pressure/volume typing. It is not copied
or silently merged into this branch; the current native edition is the
direct mathematical dependency.

## Scope of the contribution

- Standard inputs: potential reciprocity, inverse/Schur response identities,
  Hessian geometry, differentiation of symmetry invariance, Maxwell's
  equations, and the Born–Infeld constitutive model.
- Developed here: a native electromagnetic response adapter; the explicit
  curvature formula in normalized coupling and impedance; reconstruction
  from signed wave speeds; bias-selection and internal-mode comparisons;
  exact controls that distinguish point response from its higher derivatives.
- Proposed physical use: reconstruct the diagnostic from constitutive or
  wave data and compare clamped and relaxed preparations.
- Not established: literature priority, a newly observed effect, automatic
  equality to optical/Berry phase, a photon mass, a universal coupling
  constant, or a continuum Yang–Mills result.

## Verification contract

The source manifest pins the exact local dependencies and this chapter's
proof, protocol, verifier, tests, and CI workflow. The result records rational
Taylor calculations, independent Christoffel agreement, wave eigenvectors,
symmetry controls, and Schur reduction. All comparisons are finite exact
checks; the six general results are written proofs, not a proof-assistant
formalization. The existing NT rational controls are also replayed unchanged.

External references are attributed mathematical/physical sources, not
re-certified by hashing this repository. The public literature check is
targeted and does not establish global novelty.
