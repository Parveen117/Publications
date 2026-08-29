# AZQ2 — Cross-platform concordance certificate (results as they fell)

Preregistered before submission; single shot, no re-tries. Bell-state
observables compared across quantinuum.sim.h2-1e and rigetti.sim.qvm
under a declared tolerance tau = 0.10 at 3 sigma, 512 shots.

- **AZQ2-A clean: CONCORDANT on both E_ZZ and E_XX** (|d| = 0.008) —
  the platforms agree within contract; prediction passed.
- **AZQ2-B localization: E_ZZ CONCORDANT under the injected phase
  fault** — the Rz miscalibration is Z-diagonal and the auditor
  correctly left the phase-insensitive observable uncut; prediction
  passed.
- **AZQ2-B teeth: prediction FAILED.** The injected Rz(0.60) fault
  produced |d| = 0.164 on E_XX — clearly beyond tau — but the declared
  statistic landed at z = +2.6, under the 3-sigma cut. The
  preregistered fault angle sat exactly at the resolution boundary of
  512 shots (expected z ~ 3.0): a design error by the implementer,
  published as it fell. The verdict machine was honest about its own
  power; the fault magnitude/shot budget pairing was the mistake. A
  future AZQ2-B2 preregistration may test a fault sized to be
  resolvable (e.g. larger angle or larger shot budget) — committed
  before running, never retrofitted to this data.

Certificate: `AZQ2_CERTIFICATE.json` (pinned
40e16187e0c73b121e9027147171051c08c1f474f35209f5769d785a58245fdb).
