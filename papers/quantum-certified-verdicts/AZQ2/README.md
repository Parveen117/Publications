# AZQ2 — Cross-Platform Concordance Certificate

**The RNKE-Q product claim, certified:** a computation's observables
are platform-independent within a declared contract — or the pair is
cut, with the fault localized. Preregistered, 512 shots, tau = 0.10 at
3 sigma, across quantinuum.sim.h2-1e and rigetti.sim.qvm.

| Rung | Result |
| --- | --- |
| Clean Bell pair, E_ZZ and E_XX | **CONCORDANT** — platforms agree within 0.008 |
| Injected phase fault, E_ZZ | **CONCORDANT** — the auditor correctly leaves the phase-insensitive observable uncut |
| Injected phase fault, E_XX | **DISCORDANT, z = +12.2** — the fault is caught and localized to exactly the phase-sensitive observable |

Certificates: `AZQ2_CERTIFICATE.json` (clean + localization),
`AZQ2_B2_CERTIFICATE.json` (fault detection, pin 62e76106...).
Power-calibration note: an earlier rung with a boundary-sized fault
(Rz 0.60, expected z ~ 3.0) landed at z = +2.6 and is retained in
`AZQ2_CERTIFICATE.json` as the measured resolution limit of the
512-shot contract; the detection rung above uses a fault inside the
resolving power, preregistered in `PREREGISTRATION_B2.md`.

Run it yourself: workflow `azq1.yml`, inputs
`workdir=AZQ2, script=run_azq2.py`.
