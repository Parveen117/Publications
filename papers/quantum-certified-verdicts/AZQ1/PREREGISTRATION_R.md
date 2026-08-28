# AZQ1-R Preregistration — real Quantinuum H2-1 hardware rung

Committed before any real-hardware job is submitted. Identical
statistic, thresholds and circuits to AZQ1 (see PREREGISTRATION.md),
now on the real trapped-ion machine `quantinuum.qpu.h2-1`, after the
exact circuits were validated on its emulator (AZQ1 certificate pin
9cffe27a45dcdc12a1e28b41dc62eb1b81000d2f690dfebfcd654d4fb469cd40).

- Target: quantinuum.qpu.h2-1. shots = 512 per basis. p0 = 0.05
  (same as the noise-modeled emulator).
- Rungs: clean (prediction NOT_FALSIFIED) and corrupted (prediction
  FALSIFIED at z >= 5). Four jobs total, single shot each, no
  re-tries.
- Cost boundary: if the workspace plan's real-hardware quota is
  insufficient, the job error is recorded as NOT_EXECUTABLE with the
  quota message verbatim — an honest outcome, not a failure to hide.
- Results published as they fall into AZQ1_R_CERTIFICATE.json.
