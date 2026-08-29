# AZQ2 Preregistration — cross-platform concordance certificate

Committed before submission. AZQ1 certified states per platform; AZQ2
certifies AGREEMENT BETWEEN PLATFORMS — the core RNKE-Q product claim:
a computation's observable is platform-independent within a declared
contract, or the pair is cut.

## Observables and statistic

Bell state |Phi+> on 2 qubits. Observables E_ZZ and E_XX (parity
expectations from 512-shot counts, Z and X bases). Platforms:
quantinuum.sim.h2-1e (trapped-ion noise model) and rigetti.sim.qvm.

Per observable O: d = E_O(platform A) - E_O(platform B);
SE_d = sqrt(SE_A^2 + SE_B^2), SE = sqrt((1 - E^2)/shots).
Declared concordance tolerance tau = 0.10 (covers the emulator's
honest noise gap against an ideal simulator).
Verdict per observable: z = (|d| - tau)/SE_d;
z < 3 -> CONCORDANT; z >= 3 -> DISCORDANT.

## Rungs and predictions (fixed now)

- AZQ2-A clean Bell pair on both platforms: prediction CONCORDANT on
  both E_ZZ and E_XX.
- AZQ2-B miscalibration teeth: an Rz(0.60 rad) is inserted before
  measurement on qubit 0, on the Quantinuum side only (a realistic
  phase-calibration fault). Prediction: E_XX DISCORDANT (z >= 3);
  E_ZZ remains CONCORDANT (Rz is Z-diagonal). This asymmetry is the
  falsification control: the auditor must localize the discordance to
  the phase-sensitive observable only.
- 512 shots per basis per platform; single shot per rung; no
  re-tries; results as they fall into AZQ2_CERTIFICATE.json.

## Claim boundary

Concordance verdicts on the named backends under the declared
tolerance only. No hardware-quality ranking, no advantage claim.
