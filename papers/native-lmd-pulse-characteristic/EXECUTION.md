# Fresh calculator execution and source-cutoff evidence

Author: Monty. 30 October 2026.

Final clean run: **PASS — 379 checks and 25 negative controls**.
Eight fresh ideal-circuit NRQC requests completed, using two qubits and
the pinned existing calculator. There were no hardware jobs or shots.

The four states are vacuum, Bell, controlled-S on |++> and its controlled-S
inverse; the engines are observer and clifford_t. Both engines agree with
independently specified exact 4x4 state/operator matrices. Only the vacuum
is a source-moment encoding. The phase circuits exercise non-Clifford gates
but do not implement the source pulse.

## Pulse-moment comparison

Three source h values and five prescribed times (-1,-1/2,0,1/2,1) were
checked. Weighted means/covariances are assembled with the existing native
JS rational/complex engine from the Pauli output; no weighted-generator
support is invented in the NRQC API. The table shows conditional source
predictions at t=1, rounded to six decimals. The encoded source Gram entries
are each certified within 10^-6; exact error bounds are in the JSON.

| h | Actual measured variance | Naive density-shift variance | Difference |
| --- | ---: | ---: | ---: |
| 1/8 | 0.203696 | 0.387700 | 0.184004 |
| 1/4 | 0.450200 | 1.071053 | 0.620852 |
| 3/8 | 0.764428 | 2.917333 | 2.152905 |

Actual: h/2+t^2[I(h)-1]. Naive: h/2+t^2[I(2h)-1]. Every nonzero declared
time separates their certified intervals; t=0 is the expected equality.
Wrong cross-phase sign and omitted detector noise are detected separately.

The finite pair fails [D_f,J_f]=ih J_f by trace, as it must. Therefore this
is an executed source-moment contract, not an executed source time evolution.
The analytic characteristic law is proved separately in THEOREM.md.

## Native source refinement and bounded response

Positive Laurent-series interval calculations evaluated the actual source
cutoff second moments at N=0,1,2,4,8,16. The table reports the squared L2
mass error delta_N^2=I(h)-I_N(h), rounded to six decimals.

| h | N=0 | N=1 | N=2 | N=4 | N=8 | N=16 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1/8 | 0.008048 | 0.003617 | 0.002367 | 0.001431 | 0.000824 | 0.000462 |
| 1/4 | 0.041175 | 0.021033 | 0.014671 | 0.009550 | 0.005956 | 0.003629 |
| 3/8 | 0.121937 | 0.070334 | 0.052189 | 0.036526 | 0.024654 | 0.016315 |

These are squared errors, not errors in a measured probability or proof that
the last cutoff meets a small requested tolerance. The proved characteristic
error bound uses delta_N, the square root of this quantity.

The exact new bounded witness uses sin(2v_0). At zero pulse, its actual
response slope exceeds the naive slope by approximately 0.004600, 0.016276
and 0.032444 for the three h values, respectively. Native series arithmetic
certifies the sign and intervals. No NRQC computation of the full infinite
dimensional characteristic expectation is claimed.

## Reproduce

Use the unmodified quantum-classical-provisional checkout at
`fdd39048a6ddaa2dd78aa3197c274d7766f31a0e`, Python 3.12 and Node:

```sh
python3.12 papers/native-lmd-pulse-characteristic/run_validation.py \
  --nrqc /path/to/pinned/quantum-classical-provisional \
  --out /tmp/lmd-pulse-validation
```

The orchestrator checks the calculator commit and tracked-tree cleanliness,
builds certified source intervals, runs all eight calculator requests and
executes the independent native verifier. To replay just the saved finite
comparison without new calculator executions:

```sh
node papers/native-lmd-pulse-characteristic/verify.cjs \
  papers/native-lmd-pulse-characteristic/results
```

- [Verification and all checked cases](results/verification.json)
- [Source intervals, cutoff error budgets and coefficients](results/source_intervals.json)
- [Execution metadata and source hashes](results/execution.json)
- [Result manifest](results/manifest.json)
- [Proof, source/measurement assumptions and next implementation gate](THEOREM.md)

The final environment was Python 3.12.14 and Node v24.22.1. Source preparation
took approximately 2.284 seconds; individual NRQC command wall times were
0.457-0.631 seconds. These are local process observations, not a paired
benchmark or scaling result. No speed, memory or quantum advantage is claimed.
No hardware/API execution or CI dispatch was performed.
