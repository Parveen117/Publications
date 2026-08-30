# AZQ1 — Response to External Review (Aug 2026)

An external review raised five objections to the AZQ1 GHZ parity-audit
campaign. Each is answered below with pointers to files in this folder.
Two acknowledged gaps are closed by the AZQ1-I3 rung (see
`PREREGISTRATION_I3.md`), which upgrades the audit to a
calibration-derived noise budget with full histograms.

## 1. "Single shot on a GHZ state has enormous variance"

Misreading. "Single shot per rung" in the preregistration means **one
job submission, no re-tries**. Every submission measured **512 shots
per basis** (`PREREGISTRATION.md`, "Verdict rule"; every certificate
carries `"shots": 512`). The variance argument does not apply.

## 2. "The noise bound was fitted post-hoc to produce the verdicts"

The bound was fixed **before any job was submitted**
(`PREREGISTRATION.md`, first line: "Committed before any Azure job is
submitted"): p0 = 0.05, binomial SE, threshold z = 5. The same numbers
appear unchanged in every certificate.

The stronger answer is the **teeth control**: under the *same*
contract, the clean GHZ passed (z = −2.75 on `ibm_kingston`) and a
circuit corrupted by construction (X inserted mid-chain,
`run_azq1_i2.py`) failed at z = +66.6. A single fixed threshold that
simultaneously admits the clean state and rejects the corrupted one is
not a threshold that can be tuned to make "any execution pass" or
"any execution fail". The reviewer read the corrupted control as an
"uncorrected GHZ"; it is a deliberate falsification target, published
as such.

## 3. "Show raw counts, not z-scores"

Partially valid. Certificates publish the forbidden-outcome weights
per basis (W_Z, W_X), the readout-calibration control fractions, and
the shot count — from which the forbidden-outcome counts are exact
(W × shots). Full per-bitstring histograms were not pinned. AZQ1-I3
pins the complete counts dictionary for every circuit.

## 4. "Was the bound derived from the device's published error rates?"

No — and this is the one methodological upgrade the review correctly
identifies. p0 = 0.05 was a declared generic floor, not derived from
device calibration. It was fixed a priori (so it is not post-hoc), but
it was not physically motivated per-device. AZQ1-I3 replaces it with a
**declared formula**: p0 computed from the live calibration data of
the selected qubit chain (two-qubit gate errors + readout errors) at
submission time, written to a hash-pinned precommit file **before**
any job is submitted, so the derivation order is verifiable.

## 5. "Cross-platform means 1 real QPU + 2 simulators"

Correct for AZQ1, and the wording is amended to "multi-backend (one
hardware QPU, two emulators)". Hardware cross-backend replication is
established in a separate campaign: the Multi-Ledger Operator-Memory
Benchmark passed its declared release gate on **three physical IBM
backends** (`ibm_fez`, `ibm_kingston`, `ibm_marrakesh`) with
manifests, duplicate-job rejection, and cross-run stability checks —
see the `Quantum-Classical-public` repository.

## Standing invitation

The evaluator harness and certificates are public. Independent
replication — running the declared contracts on your own account and
comparing pins — is invited and is the intended use of this release.
