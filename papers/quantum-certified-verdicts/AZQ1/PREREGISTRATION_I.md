# AZQ1-I Preregistration — real IBM superconducting hardware rung

Committed before any IBM job is submitted. Identical circuits,
statistic and thresholds to AZQ1; target = least-busy operational IBM
QPU on the open plan (name recorded in the certificate). shots = 512
per basis; p0 = 0.05 (real superconducting error floor, same as the
noisy-emulator allowance). Rungs: clean -> NOT_FALSIFIED predicted;
corrupted -> FALSIFIED at z >= 5 predicted. Four jobs, single shot
each, no re-tries; quota or access errors recorded verbatim as
NOT_EXECUTABLE. Results as they fall into AZQ1_I_CERTIFICATE.json.
Together with AZQ1 (Quantinuum H2 emulator + Rigetti QVM) this
completes three backend families for the cross-platform claim.
