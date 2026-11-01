# Source and execution map

Author: Monty. 30 October 2026.

Frozen Publications parent: `70aef2535a1c124a410e4a4bf96780236ca25c23`.
Frozen calculator repository: Parveen117/quantum-classical-provisional,
commit `fdd39048a6ddaa2dd78aa3197c274d7766f31a0e`.

- LMD's THEOREM.md supplies the logarithmic Gaussian representation,
  source Wick martingale, L2 domain, zero-mode scaling and density pulse.
- native-lmd-common-parent/THEOREM.md supplies the selected measurement,
  h/4 added conjugate noise, doubled density covariance and the variance
  test. Its detector, state and source assumptions are preserved.
- native-nrqc-gravity-bridge/source_interval.cjs supplies the already-used
  outward rational log-Gamma/zeta/exponential method, adapted here to I(h).
  Its frozen pilot is unchanged.
- native-common-cone-limit/native_engine/native_operator.cjs supplies all
  exact rational/complex finite matrix operations. Python 3.12 only launches
  commands, records files, hashes and runtimes.
- NRQC calculator.py and emk_tensor.py define the actual supported signed
  Pauli geometry interface. observer and clifford_t are explicitly selected;
  there is no backend fallback and no change to the calculator repository.

Standard analytic ingredients are the Gaussian unitary transformation to
Lebesgue measure, the translation group for i h partial_q, the spectral
theorem, L2 martingale orthogonality, and Fourier inversion. Each is applied
explicitly in THEOREM.md. No external residual-gauge result is imported.

Primary special-function references, consulted 30 October 2026:

- [NIST DLMF 5.12.1-2](https://dlmf.nist.gov/5.12): Euler beta and sine-power
  integrals used in I(g).
- [NIST DLMF 5.5.5](https://dlmf.nist.gov/5.5#E5): gamma duplication.
- [NIST DLMF 5.7.3](https://dlmf.nist.gov/5.7#E3): convergent log-Gamma series.
  In the ratio the Euler constant cancels, giving
  log I(g)=g+sum_(n>=2) zeta(n)(1-2^(1-n))g^n/n.

The exact characteristic formula, bounded response difference, cutoff
error transfer and their scope are written derivations in this packet.
The finite checks do not constitute a proof assistant verification of those
continuum statements. Raw executions and machine-readable verification are
retained in results/, with a reproducible command in EXECUTION.md.
