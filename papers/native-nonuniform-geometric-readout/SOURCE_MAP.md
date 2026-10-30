# Sources, choices and verification boundary

Author: Monty. 30 October 2026.
Frozen parent: `3b487951ccf5c4939a8c0249b435c904f47fedfe`.

## Frozen source inputs

- [LMD](../native-logarithmic-material-density/THEOREM.md), sections 2–6:
  Gaussian mode variances, positive full-support Wick measure, affine
  conjugates, coherent states and density covariance. All continuum results
  here are conditional on those frozen inputs, not a new independent audit.
- [MPQ](../native-material-profile-quantization/THEOREM.md), section 2:
  the original y=d/j chart, a(theta), Y, the inverse geometry including P,
  and the exact density-pulse source solution. Its different quantum Gaussian
  representation is not imported.
- [RMG](../native-resolved-material-geometry/THEOREM.md), section 1:
  original spatial coframe and the need to distinguish pointwise fields from
  measure-valued readings.
- [Radial predecessor](../native-radial-geometric-readout/THEOREM.md):
  total-mass factorization, plus-sign Fourier convention and the limitation
  to a constant classical conjugate ratio. Its POVM is preserved unchanged.
- [Earlier prescribed-angle reading](../native-nrqc-gravity-bridge/BOUNDED_COFRAME.md):
  normalized coframe formula and outcome-wise bounds. Here the finite profile
  and density-shape coefficient are outcomes of a different joint POVM.
- [Native rational engine](../native-common-cone-limit/native_engine/native_operator.cjs):
  unchanged exact F arithmetic and finite matrices for the Taylor diagnostics.

Exact blob and SHA-256 pins are in SOURCE_PINS.json. The executable and runner
are hashed as well. No external Fable result is a new dependency.

## Established methods, credited separately

1. Y. Shikano and A. Hosoya, *Optimal Covariant Measurement of Momentum on a
   Half Line in Quantum Mechanics*, J. Math. Phys. 49, 052104 (2008),
   [arXiv:0710.1724](https://arxiv.org/abs/0710.1724),
   [full text](https://arxiv.org/html/0710.1724v3).
   Context for generalized momentum measurements and larger-space dilation.
   The two-coordinate cone compression is specified and normalized directly
   here; no optimality or apparatus theorem is imported.
2. C. Carmeli, T. Heinosaari and A. Toigo, *Sequential measurements of conjugate
   observables*, J. Phys. A: Math. Theor. 44, 285304 (2011),
   [arXiv:1105.4976](https://arxiv.org/abs/1105.4976),
   [DOI](https://doi.org/10.1088/1751-8113/44/28/285304).
   Established framework for sequential generalized measurements of conjugate
   variables. Our q1 filter is nonlinear in the mass coordinates; we do not
   assert it is a Weyl-Heisenberg covariant phase-space observable. The exact
   covariance claimed here is only the momentum translation under the
   specified density pulses, proved in THEOREM.md.

Other tools used are elementary exponential tilting on a compact support,
change of variables in L2, Plancherel, Gaussian filtering, PVM compression by
an isometry, and localized semiclassical estimates. These are standard methods,
not new general results. The source-specific construction is the mass chart,
its selected conjugate/shape measurement, the nonlinear coframe pushforward
and its limited nonuniform source dynamics and classical recovery.

## Choices and nonclaims

- The physical Hilbert space stays LMD; no finite-dimensional replacement.
- Retaining all residual eta modes does not make the conjugate resolution
  complete: only span{1,cos(2 pi z)} is resolved.
- The Gaussian shape filter, sigma and cone Fourier POVM are added ideal
  measurement choices. They are not uniquely derived from constraints.
- The sharp density/conjugate joint-measurement obstruction is not bypassed;
  Q has added noise and the p marginal is a generalized measurement.
- The full finite-h configuration law is unchanged by a density pulse.
  Phase sensitivity comes from the p outcomes, not hidden density evolution.
- A fixed cosine calibration is not covariant under all material flows.
- Coframes are measure valued. No raw finite-mean length, pointwise metric,
  spacetime curvature, local apparatus, speed advantage or 3+1 closure follows.

## Validation performed

Run `python3.12 papers/native-nonuniform-geometric-readout/run_validation.py`.
Python performs only process orchestration, hashes and I/O. The existing native
JS exact arithmetic computes the scientific fixtures. The run verifies source
pins and reproduces VALIDATION.json; `--write` explicitly replaces that file.

The finite diagnostic solves Theta'=a(Theta) as a rational Taylor jet, checks
the independent half-angle denominator identity, integrates exact cosine
moments, and compares the circumference's second derivative with a separate
differential calculation. It also records the error from dropping P feedback,
the failure of a radial-only profile and the omitted second harmonic.

Validation: 16 exact finite checks and 5 negative controls pass. The initial
draft check for Theta's fifth derivative had an incorrect expected constant;
the ODE coefficient equation 5 a5=63/256 fixes it to 189/32. This correction
is included before the recorded successful run; no predecessor result changed.

Written review checked: the first cosine variance h; the full removed Wick
variance h(1+c^2); det H rather than x0 as Jacobian; cone surjectivity; the
Fourier pulse sign; the Gaussian filter normalization and noise derivative;
the distinction between a projection and a full field inverse; compact-cone
localization before inverse-H estimates; and the nonlinear P reconstruction.

Finite tests do not prove the continuum chart, POVM normalization, classical
limit, all-beta strict circumference inequality or Taylor remainder bounds.
Those statements have written arguments in THEOREM.md. No independent peer
review, fresh NRQC execution, hardware run or CI dispatch accompanies this packet.
