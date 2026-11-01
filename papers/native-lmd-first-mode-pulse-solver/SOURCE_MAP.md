# Sources and method boundary

Author: Monty. 31 October 2026.

- Frozen Publications parent: `13389614767542792097c189dcf6e11747806b42`.
- [LMD source](../native-logarithmic-material-density/THEOREM.md): Fourier
  covariance convention, positive Wick martingale and L2 source mass.
- [Common parent](../native-lmd-common-parent/THEOREM.md): specified noisy
  detector and its calibration; unchanged in this refinement.
- [Characteristic law](../native-lmd-pulse-characteristic/THEOREM.md): actual
  uniform-vacuum marginal, phase direction and Gaussian damping.
- [Conditional solver](../native-lmd-certified-pulse-solver/THEOREM.md):
  conditional-mean quadratic remainder, Gaussian bins, phase interpolation
  and native directed rational elementary functions. `phase.cjs` adapts its
  interpolation into a module without changing the error theorem.
- [Frozen moments](../native-lmd-pulse-characteristic/results/source_intervals.json):
  source I(gamma) bounds and a separately evaluated N=1 Laurent moment.
- [Old output](../native-lmd-certified-pulse-solver/results/solver.results.json):
  four wide baseline intervals and eight earlier successes.
- The nine pinned predecessor NRQC outputs are reused as phase data, with
  their exact-backend statuses and certified decimal enclosures checked.
  There are **zero fresh calculator requests** in this packet.
- Native scientific arithmetic: the frozen `interval.cjs` and its native
  `F` field. Python 3.12 performs orchestration and hashes only.
- Standard modified Bessel identities: NIST DLMF
  [10.25.2](https://dlmf.nist.gov/10.25.E2) (power series) and
  [10.32.1](https://dlmf.nist.gov/10.32.E1) (circle integral). The needed
  specialization, Gaussian second moment and Poisson/Gamma cell identity
  are derived explicitly in THEOREM.md. They are not claimed as new
  special-function or positivity theorems.

SOURCE_PINS.json checks inherited files. The run's execution.json binds
the new implementation and request; manifest.json binds all its outputs.
The original field, cutoff measure, pulse and parent detector have not been
replaced by a two-variable physical system. The compression evaluates an
expectation in the already derived measurement law.
