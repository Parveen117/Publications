# Sources and method boundary

Author: Monty. 31 October 2026. Frozen parent: `ca080861d646774ad94a8e48be668194a8086d59`.

- [LMD-T3](../native-logarithmic-material-density/THEOREM.md) supplies the
  positive Wick martingale and two-point kernel for bounded weighted masses.
- [Common parent](../native-lmd-common-parent/THEOREM.md) fixes the detector.
  [Pulse characteristic](../native-lmd-pulse-characteristic/THEOREM.md) supplies
  its bounded translation proof; constant translation scales the weighted
  source too, as proved in this packet.
- [Nonuniform next-step derivation](../native-lmd-exact-mass-marginal/NEXT_STEP.md)
  supplies the weighted second moment and the uniform-law mismatch witness.
  Rémy's uniform mass theorem is not assumed to describe this weighted law.
- [First-mode solver](../native-lmd-first-mode-pulse-solver/THEOREM.md) supplies
  the q/radial representation and the integer-Gamma cells. Its raw q and
  radial engine is reused; its nine earlier phase outputs remain pinned.
- [Exact-mass Gamma implementation](../native-lmd-exact-mass-marginal/THEOREM.md)
  supplies native half-shape Gamma integrals. Using this integration code
  does not apply its uniform probability-law theorem to the weighted model.
- Standard Bessel identities are credited to NIST DLMF:
  [10.25.2](https://dlmf.nist.gov/10.25.E2), modified power series;
  [10.32.3](https://dlmf.nist.gov/10.32.E3), integer-order circle integral;
  [10.2.2](https://dlmf.nist.gov/10.2.E2), ordinary power series;
  [10.9.1](https://dlmf.nist.gov/10.9.E1), ordinary order-zero circle integral.
  The weighted conditional moments, positive half-Gamma cells and complete
  error decomposition are derived explicitly here. No new Bessel or general
  conditional-expectation theorem is claimed.

The written proofs carry continuum claims. Executed checks test finite
implementation, bounds and normalization. Native rational arithmetic is the
scientific engine; Python is orchestration only. This is a prescribed
nonuniform source and one global marginal, not a field-state preparation,
joint geometry solver, hardware run or autonomous gravity construction.
