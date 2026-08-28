# AZQ1 Preregistration — GHZ parity-chart audit on Azure Quantum backends

Committed before any Azure job is submitted. First rung of the RNKE-Q
cross-platform certification programme: the same typed-verdict audit,
identical statistic, run on two independent backend families
(Quantinuum H-series emulator with hardware noise model; Rigetti QVM).
IBM-hardware rung follows separately.

## Observable and identity

Prepare GHZ_n (n = 3): H on q0, CX chain. A GHZ state satisfies the
parity chart: (a) Z-basis: only all-0 / all-1 outcomes; (b) X-basis
(H on all before measure): even total parity only. Declared statistic
per basis: excess weight W = (observed forbidden-outcome probability).
Acceptance identity: W_Z + W_X = 0 up to shot noise.

## Verdict rule (fixed now)

- shots = 512 per basis per backend.
- Noise allowance: binomial SE at p0 = the backend's declared/typical
  error floor; threshold: z = (W - p0)/SE with p0 = 0.05 for the
  noisy Quantinuum emulator and p0 = 0.01 for the ideal Rigetti QVM.
- Verdicts: z < 5 -> NOT_FALSIFIED (COMMITTED chart);
  z >= 5 -> FALSIFIED_MEASUREMENT_CONTRACT (CUT).

## Rungs and predictions

- AZQ1-A clean GHZ, both backends -> NOT_FALSIFIED.
- AZQ1-B corrupted circuit (an X inserted on q1 mid-circuit, breaking
  the parity chart by construction) -> FALSIFIED on both backends.
  This is the teeth control: if the corrupted circuit is not rejected,
  the audit has no power and that is reported.
- Single shot per rung per backend; results published as they fall;
  certificate JSON pinned into this folder with SHA-256.

## Claim boundary

Chart-consistency verdicts on declared backends only. No quantum
advantage, supremacy, or hardware-quality claim.
