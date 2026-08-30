# AZQ1-I3 Preregistration — calibration-derived GHZ parity audit (IBM hardware)

Committed before any I3 job is submitted. This rung upgrades AZQ1-I2 on
the two points identified in external review: the noise budget is now
**derived from the device's live calibration data before submission**,
and **full per-bitstring histograms are pinned** in the certificate.

## Observable and identity (unchanged from AZQ1)

GHZ_3: H on q0, CX chain. Parity chart: (a) Z basis — only 000/111;
(b) X basis — even total parity. Statistic per basis: forbidden-outcome
weight W. Audit statistic: W = (W_Z + W_X)/2.

## Noise budget (declared formula — fixed now)

Let (a,b,c) be the qubit chain selected by the same lowest-error rule
as I2, from the backend target at submission time. Define, from that
same calibration snapshot:

    p0 = e2(a,b) + e2(b,c) + r(a) + r(b) + r(c)

where e2 is the reported two-qubit gate error of the chain pairs and
r the reported readout error of each chain qubit. This first-order
error budget is the declared bound; no other tuning parameter exists.
A floor of p0 >= 0.005 applies (guards against zero-reported errors).

**Sequencing guarantee:** the runner computes p0, records the chain,
the individual calibration values, and the timestamp into
`AZQ1_I3_PRECOMMIT.json`, prints its SHA-256, and only then submits
jobs. The certificate must reference that precommit hash. Any
mismatch voids the run.

## Verdict rule (fixed now)

- shots = 8192 per basis per rung.
- SE = sqrt(p0(1−p0)/(2·shots)); z = (W − p0)/SE.
- z < 5 → NOT_FALSIFIED (COMMITTED chart);
  z ≥ 5 → FALSIFIED_MEASUREMENT_CONTRACT (CUT).

## Rungs and predictions

- I3-A clean GHZ with declared error suppression (DD XY4 + twirling,
  same options as I2) → predicted NOT_FALSIFIED.
- I3-B corrupted circuit (X on q1 mid-chain, breaking the chart by
  construction) → predicted FALSIFIED. Teeth control: if it is not
  rejected, the audit has no power and that is reported.
- One submission per rung per basis, no re-tries; full counts
  dictionaries for every circuit (calibration controls included) are
  embedded in the certificate and pinned by SHA-256.

## Claim boundary

Chart-consistency verdicts on the declared backend only. No quantum
advantage, supremacy, or hardware-quality claim.
