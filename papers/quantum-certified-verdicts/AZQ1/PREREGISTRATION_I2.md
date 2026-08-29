# AZQ1-I2 Preregistration — real IBM hardware with declared error suppression

Committed before submission. Statistic, thresholds (p0 = 0.05, 5 sigma),
shots (512/basis) and circuits IDENTICAL to AZQ1-I — the contract is
not loosened. What changes is experimental technique, declared here:
(1) qubit selection: the connected 3-qubit chain minimizing the sum of
two-qubit gate error and readout error from the device's live
calibration data, pinned in the certificate; (2) dynamical decoupling
XY4; (3) gate and measurement twirling; (4) transpiler optimization
level 3; (5) a descriptive readout-calibration control (prepare 000
and 111, report flip rates — no threshold attached).

Predictions: clean -> NOT_FALSIFIED; corrupted -> FALSIFIED at
z >= 5. Single shot per rung, no re-tries, results as they fall into
AZQ1_I2_CERTIFICATE.json. Paired with AZQ1-I this tests whether the
auditor tracks the QUALITY DIFFERENCE between uncorrected and
error-suppressed execution on the same machine.
