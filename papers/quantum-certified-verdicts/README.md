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

**Sensor-layer lineage:** the programme's verdict engine was first
proven at the quantum-sensor layer — an NV/ODMR magnetometer
verification cartridge (digital-twin campaign, seven cases, all
predeclared statuses matched, certificate pinned 1ec177be..., eleven
tests) and a phase-trust rule separating clean phase walk / correctable
decoherence / leakage into ALLOW / CORRECT / REJECT — the same
three-region behaviour later demonstrated on a real QPU in AZQ1
(refuse uncorrected, certify error-suppressed, cut corrupted). The
sensor-layer sources live in the private hardware repository and are
available for technical diligence.

**Sibling evidence:** the [Quantum-Classical-public](https://github.com/Parveen117/Quantum-Classical-public) repository carries the programme's heaviest hardware result — a seven-qubit multi-ledger operator-memory benchmark with certified CROSS_BACKEND_REPLICATION across three real IBM QPUs (ibm_fez, ibm_kingston, ibm_marrakesh; 3/3 runs, all six per-run gates on every backend), with SHA-256 manifests and claim-scope documentation.

The instrument measures execution quality across platforms — it does
not rubber-stamp it. Every preregistration is git-timestamped before
its run; every certificate carries a SHA-256 pin; workflows `azq1.yml`
and `azq1_ibm.yml` re-execute any rung on demand.
