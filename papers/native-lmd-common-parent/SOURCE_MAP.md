# Sources, inherited methods and added choices

Author: Monty. 30 October 2026.
Frozen parent: `349175122103717730fce423ce7728695f566d9f`.

## Frozen inputs

- [LMD](../native-logarithmic-material-density/THEOREM.md): the actual
  logarithmic covariance, Gaussian physical space, Wick measure in the L2
  range, coherent states, affine conjugates, density pulses and corrected
  material-flow coefficient. The source measure uses parameter h; the new
  noisy outcome measure uses 2h. Both are kept distinct.
- [CMG](../native-coarse-material-geometry/THEOREM.md): the countable Hermite
  isometry with exact finite marginals; absolute continuity of all normal
  outcome laws; smoothing, rejection-inclusive measurements and closed mean
  forms. These are inherited methods, not newly discovered constructions.
  CMG's old MPQ covariance and its exp(smoothed log density) geometry are
  not silently transferred: here the density is a smoothed Wick measure.
- [RMG](../native-resolved-material-geometry/THEOREM.md): finite coherent
  marginal formulas, their added variances, the inverse-chart interpretation
  of spatial geometry and the bounded-periodic decomposition of Theta.
- [MPQ](../native-material-profile-quantization/THEOREM.md), section 2:
  exact source inverse, global P feedback and prescribed density-pulse orbit.
  Its Gaussian measure is not identified with LMD.
- [Refinement and dynamics](../native-geometry-refinement-dynamics/THEOREM.md):
  the precise sharp-residual incompatibility and the requirement for a common
  noisy parent. Its mass-chart measurements remain frozen and are not claimed
  as marginals of the new parent.
- [Nonuniform geometry](../native-nonuniform-geometric-readout/THEOREM.md):
  the two-window classical comparison and the need to preserve the original
  nonlinear denominator. Its simple conjugate-shift covariance is not
  transferred to the new noisy adapted-coordinate measurement.
- [MDD](../native-material-detector-dilation/THEOREM.md): distinction between
  an ideal generalized measurement and a source apparatus, including the
  phase blindness of a density-only interaction protocol.

SOURCE_PINS.json records the unchanged input blobs and SHA-256 digests.
There is no new dependency on an external residual-gauge/Fable result.

## Established external methods

1. B. Chakraborty, J.-P. Gazeau and A. Youssef, *Coherent state quantization
   of angle, time, and more irregular functions and distributions*,
   [arXiv:0805.1847v2](https://arxiv.org/abs/0805.1847v2).
   Established coherent-state/anti-Wick methods, already credited in CMG.
   The countable isometry and finite marginals needed here are proved directly
   in the source normalization; no novel general quantization theorem is claimed.
2. N. Berestycki, *An elementary approach to Gaussian multiplicative chaos*,
   [arXiv:1506.09113v2](https://arxiv.org/abs/1506.09113v2),
   [full text](https://arxiv.org/html/1506.09113v2).
   Standard renormalized measure construction and regularization methods,
   already used by LMD. We use the frozen L2 proof at parameter 2h<1 and
   do not import a critical/supercritical extension or identify this source
   model with Liouville gravity.

Other standard tools are conditional expectation as an L2 contraction,
Tonelli, Gaussian moments, Sobolev embedding, circle heat-kernel smoothing,
normal ordering in the Hermite/Bargmann representation, the spectral theorem,
and representation of a closed positive form by a self-adjoint operator.
The relevant constants, domains and convergence arguments are stated directly.

## New contract and exact scope

The parent measurement adds matched Gaussian noise to all adapted-coordinate
modes. The outcome Wick measure subtracts h K_N, accounting for both source
and measurement fluctuations. The geometric map smooths this measure and v
at a declared ell>0 before applying the original inverse. Optional delta>0
records rejection when the coarse density is too small; discarded outcomes
are not silently removed or renormalized.

The source-specific results are density first-moment calibration, fixed-h
compatible coarse geometry, the exact compressed pulse generator and its
cutoff limit, the variance mismatch of a naive classical label update, and
the LMD smoothing budget h ell^(-8) for classical recovery. These do not
derive a unique detector, a local apparatus, a zero-resolution quantum metric,
autonomous gravity or 3+1 dynamics. Material-flow covariance needs an audit
of the detector noise and calibration, not only of the physical density.

This packet contains written analytic arguments only. No numerical execution,
certificate replay, NRQC benchmark, hardware experiment or CI dispatch was
performed. REVIEW.md records the proof checks and remaining obligations.
