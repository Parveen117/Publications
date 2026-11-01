# Source and algorithm map

Author: Monty. 30 October 2026.
Parent: `ce60e932bb46142eafeef304c099c15ea89ad088`.

- [LMD](../native-logarithmic-material-density/THEOREM.md): source Gaussian
  space, zero-mode factorization, Wick martingale, mean-one source density
  and finite second moment for 0<gamma<1.
- [Common parent](../native-lmd-common-parent/THEOREM.md): the measurement
  marginal, added conjugate noise and doubled configuration covariance.
- [Characteristic law](../native-lmd-pulse-characteristic/THEOREM.md): the
  exact target, phase sign, Gaussian damping and the previous cutoff bound.
- [Frozen source intervals](../native-lmd-pulse-characteristic/results/source_intervals.json):
  enclosures of I(gamma), derived using the credited log-Gamma/zeta method.
  The native solver reuses these values; the file is pinned and unchanged.
- [Native rational engine](../native-common-cone-limit/native_engine/native_operator.cjs):
  exact F arithmetic. The new adapter adds directed interval operations and
  bounded elementary functions; it does not replace the source arithmetic
  with floating-point library functions.
- Pinned Parveen117/quantum-classical-provisional commit
  `fdd39048a6ddaa2dd78aa3197c274d7766f31a0e`: calculator.py,
  cyclotomic64.py and phase_backend.py. Actual p64/X/Y requests use the
  documented exact phase backend with no fallback or engine modification.

Standard methods used here are L2 conditional-expectation orthogonality,
Gaussian completion of the square, the integral Taylor remainder, cubic
Lagrange interpolation, alternating-series remainders and integer square
roots. Their application and needed constants are proved in THEOREM.md.
No new external continuum theorem, source measure or residual-gauge result
is assumed. Python 3.12 is I/O orchestration only.

SOURCE_PINS.json binds inherited sources; results/execution.json binds the
executed implementation; results/manifest.json binds all retained outputs.
