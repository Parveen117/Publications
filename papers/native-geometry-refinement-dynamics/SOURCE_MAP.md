# Source map and scope of the new arguments

Author: Monty. 30 October 2026.
Frozen parent: `eb7916e6e4f4bf30db4f98f3792b46d60b4ba636`.

## Frozen repository inputs

- [LMD](../native-logarithmic-material-density/THEOREM.md), sections 2–6:
  the actual Hilbert space, logarithmic Fourier covariance, Wick measure,
  smooth density pulses and coherent-state moment estimates. This packet
  relies on those written results; it does not independently re-prove them.
- [Two-window geometric readout](../native-nonuniform-geometric-readout/THEOREM.md):
  the finite mass-coordinate chart, Gaussian shape filter, Fourier compression,
  phase convention, localizing argument and the stated refinement obligation.
- [Radial predecessor](../native-radial-geometric-readout/THEOREM.md):
  the one-mass special case, positive-sector sharp-ratio boundary and the
  distinction between bounded readings and divergent raw length moments.
- [MPQ](../native-material-profile-quantization/THEOREM.md), section 2:
  the source inverse y=d/j, theta=Theta(y), m=j/a(theta), P=integral m,
  and exact density-pulse evolution. Only these source identities are reused;
  the older quantum Gaussian measure class is not substituted for LMD.
- [RMG](../native-resolved-material-geometry/THEOREM.md):
  the interpretation of the spatial coframe, its insufficiency to specify
  spacetime curvature, and measurement-prescription-specific obstructions.
- [CMG](../native-coarse-material-geometry/THEOREM.md), sections 1–3:
  an existing countable noisy parent measurement with consistent finite
  marginals, spatial smoothing and a detection threshold in the older MPQ
  representation. This is an alternative design precedent, not a newly
  claimed construction or an automatic transfer of its geometry to LMD.
- [MDD](../native-material-detector-dilation/THEOREM.md), sections 1–2:
  the distinction between a declared generalized measurement, an actual source
  apparatus and the original density-only interaction class.

Git blob IDs and SHA-256 digests of these unchanged inputs are recorded in
SOURCE_PINS.json. No external residual-gauge/Fable theorem is imported.

## Established methods

T. Heinosaari, D. Reitzner and P. Stano, *Notes on Joint Measurability of Quantum
Observables*, Foundations of Physics 38, 1133–1147 (2008),
[arXiv:0811.0783v2](https://arxiv.org/abs/0811.0783v2),
[full text](https://arxiv.org/html/0811.0783v2),
[DOI](https://doi.org/10.1007/s10701-008-9256-7).
Proposition 8 supplies the standard sharp-marginal commutation criterion.
THEOREM.md repeats its short effect-support argument and applies it to the
explicit characteristic operator of the new conjugate mode. The general
criterion is not an original theorem of this project.

C. Carmeli, T. Heinosaari and A. Toigo, *Sequential measurements of conjugate
observables*, J. Phys. A: Math. Theor. 44, 285304 (2011),
[arXiv:1105.4976](https://arxiv.org/abs/1105.4976),
[DOI](https://doi.org/10.1088/1751-8113/44/28/285304).
Background for sequential unsharp measurements. The nonlinear mass-coordinate
filter used here is normalized directly; a full Weyl-covariant phase-space
classification is not assumed to apply to it.

Finite exponential-family convexity, L2 coordinate transformations, Plancherel,
Gaussian filter overlaps, orthogonal projection, finite Fourier norm estimates
and dominated convergence are standard mathematical methods. Their proofs or
required inequalities are given directly for this model. No new general
positivity, Fourier, uncertainty or convexity theorem is claimed.

## What this packet changes

It resolves the specific complete-record refinement question negatively for
the predecessor's sharp-residual prescription. It extends the finite-window
construction and classical comparison to arbitrary finite Fourier families,
derives a bounded observable's exact norm-differentiable pulse response, and
connects their ordered classical/refinement limit to the full smooth source
profile and its existing probe-driven evolution.

It does not silently replace a rejected strategy with a new physical law.
The Gaussian filter and cone Fourier POVM remain declared choices. The
incompatibility theorem does not settle the geometry-only marginals. The
ordered limit does not provide a fixed-h quantum continuum or autonomous
gravitational dynamics. REVIEW.md separates these outcomes explicitly.
