# Quantum Certified Verdicts — RNKE-Q evidence ladder

Typed verdicts (COMMITTED / CUT, CONCORDANT / DISCORDANT) for quantum
executions, preregistered before every submission, pinned, and
re-runnable from this repository's CI at zero cost.

## Certified results

- **Real hardware (AZQ1):** on a real 156-qubit IBM Heron QPU
  (ibm_kingston), the auditor certified an error-suppressed GHZ state
  (z = -2.8), refused an uncorrected one (z = +40.8), and cut a
  deliberately corrupted circuit both times — same unchanged contract.
  Matching runs on a Quantinuum H2 trapped-ion emulator and Rigetti
  QVM. Five pinned certificates.
- **Cross-platform concordance (AZQ2):** Bell observables agree across
  two independent platforms within a declared tolerance (certified),
  and an injected phase-calibration fault is caught at z = +12.2 and
  localized to exactly the phase-sensitive observable.

The instrument measures execution quality across platforms — it does
not rubber-stamp it. Every preregistration is git-timestamped before
its run; every certificate carries a SHA-256 pin; workflows `azq1.yml`
and `azq1_ibm.yml` re-execute any rung on demand.
