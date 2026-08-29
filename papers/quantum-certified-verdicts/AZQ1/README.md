# AZQ1 — Cross-platform typed-verdict certification (GHZ parity chart)

First certified rung of the RNKE-Q programme: the same preregistered
audit — a GHZ state's parity chart, tested with a declared noise
contract (p0 = 5%, 5 sigma) and a deliberately corrupted control —
executed across three quantum backend families at zero cost, driven
entirely from this repository's CI.

Every run below was preregistered (statistic, thresholds, circuits
committed to git BEFORE submission), single-shot with no re-tries, and
its certificate is pinned. Read the four certificates in order — they
tell one story.

| Certificate | Backend | Clean GHZ | Corrupted GHZ |
| --- | --- | --- | --- |
| `AZQ1_CERTIFICATE.json` | Quantinuum H2 emulator (trapped-ion noise model) | NOT_FALSIFIED (W=0.003) | FALSIFIED z=+65.9 |
| `AZQ1_CERTIFICATE.json` | Rigetti QVM | NOT_FALSIFIED (W=0.000) | FALSIFIED z=+157.6 |
| `AZQ1_R_CERTIFICATE.json` | Quantinuum H2-1 real QPU | NOT_EXECUTABLE — the workspace plan is emulator-only; recorded verbatim, not hidden | — |
| `AZQ1_I_CERTIFICATE.json` | **ibm_kingston, real 156-qubit Heron QPU** (uncorrected execution) | **FALSIFIED z=+40.8 — the auditor REFUSED to certify** a state the hardware had genuinely degraded (clean-prediction failure, published as it fell) | FALSIFIED z=+73.4 |
| `AZQ1_I2_CERTIFICATE.json` | **ibm_kingston** (declared error suppression: best calibrated chain, DD XY4, gate+measure twirling; readout controls 99.2%/100%) | **NOT_FALSIFIED z=−2.8 — CERTIFIED** | FALSIFIED z=+66.6 |

## What the pair AZQ1-I / AZQ1-I2 demonstrates

Same machine, same circuits, same unchanged contract. Uncorrected
execution: verdict CUT. Error-suppressed execution: verdict COMMITTED.
Corrupted circuit: CUT both times. The auditor measures execution
quality — it neither rubber-stamps a degraded state nor blocks a good
one, and the contract was never loosened to manufacture a pass.

## Reproduce / re-run

GitHub Actions workflows `azq1.yml` (Azure targets) and `azq1_ibm.yml`
(IBM, script selectable) re-execute any rung; outputs land in this
folder as they fall. Preregistrations: `PREREGISTRATION.md`, `_R`,
`_I`, `_I2`. Total cost of everything above: zero.
