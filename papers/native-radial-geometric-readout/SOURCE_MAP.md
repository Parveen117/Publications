# Sources, assumptions and review boundary

Author: Monty. 30 October 2026.

Frozen repository parent: `cf72a467522b3516d3c3c5531c7666ff2da96a2c`.

## Internal inputs

- [LMD](../native-logarithmic-material-density/THEOREM.md), sections 2–6:
  Gaussian zero-mode factorization; the positive Wick measure; Cameron--Martin
  conjugates; corrected material-flow unitaries; coherent-state density moments.
- [MPQ](../native-material-profile-quantization/THEOREM.md), section 2:
  y=d/j, the exact inverse chart, and the density-pulse law. Classical source
  identities are reused without asserting equivalence of quantum measure classes.
- [RMG](../native-resolved-material-geometry/THEOREM.md), section 1:
  spatial coframe exp(r)dz and the linear-plus-periodic behavior of Theta(y).
- [HQS-T4](../native-homogeneous-quantum-sector/THEOREM.md), section 5, and
  MPQ section 3 already record the standard half-line and positive-spectrum
  canonical-observable obstructions. We reuse that distinction; the new work
  is the LMD radial factorization and the specified source-geometric readout.
- [Previous coframe reading](../native-nrqc-gravity-bridge/BOUNDED_COFRAME.md):
  the restriction that the angle was supplied externally and density-only
  readings could not distinguish conjugate phases.

## Established external methods

1. Y. Shikano and A. Hosoya, *Optimal Covariant Measurement of Momentum on a
   Half Line in Quantum Mechanics*, J. Math. Phys. 49, 052104 (2008),
   [arXiv:0710.1724](https://arxiv.org/abs/0710.1724),
   [DOI](https://doi.org/10.1063/1.2917064).
   Half-line momentum requires a generalized measurement; a larger-space
   construction gives a covariant POVM. Our zero-extension/Fourier isometry is
   specified directly. No optimality or apparatus claim from that paper is
   imported into this source model.
2. M. H. Al-Hashimi and U.-J. Wiese, *Canonical quantization on the half-line
   and in an interval based upon a new concept for the momentum in a space
   with boundaries*, Phys. Rev. Research 3, 033079 (2021),
   [arXiv:2103.01715](https://arxiv.org/abs/2103.01715),
   [DOI](https://doi.org/10.1103/PhysRevResearch.3.033079).
   Different momentum concepts and enlarged constructions illustrate why the
   no-extension result must be restricted to the specified minimal derivative.
   Such alternatives are not silently substituted here.
3. L. Gouba, *Affine Quantization on the Half Line* (2021),
   [arXiv:2005.08696](https://arxiv.org/abs/2005.08696).
   Context for the standard positive-coordinate/dilation pair. No gravitational
   closure or independent source theorem is taken from this reference.

The proofs additionally use change of variables in L2, the one-dimensional
Sobolev trace, Plancherel, the deficiency-index criterion, compression of a PVM
by an isometry, and elementary Gaussian and Cauchy integrals. They are standard
tools, not new discoveries. THEOREM.md provides the source-specific formulas.

## Written audit performed

- Keep D0=D(1) distinct from the material reference charge R[1].
- Keep the complete Wick zero-mode factor exp(-h/2).
- Check the radial Jacobian p_h(q)/x and the resulting dilation term +1/2.
- Check deficiency signs for +i h partial_x: (0,1) on a scalar fiber.
- Check the Fourier sign against the uniform pulse: y shifts by +beta.
- Require a zero Sobolev boundary trace only for finite second y moments;
  do not apply it to all normalized measurement states.
- Use a compact mass cutoff in the classical-limit proof, avoiding an
  unsupported inverse-GMC-moment assumption.
- Keep the Cauchy location distinct from an undefined expectation.
- Check the raw exponential-length moment divergence instead of inferring
  moment convergence from weak convergence.
- Restrict the source inverse identification to constant classical y.

These are written analytic checks. No fresh numerical tests, certificate
replay, independent peer review, CI dispatch or hardware validation was
performed for this packet. The earlier NRQC pilot and its evidence are frozen.
