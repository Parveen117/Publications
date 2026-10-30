# Executed certified marginal solver

Author: Monty. 30 October 2026.

The first certified solver run completed: **286 checks and 4 negative
controls passed**, with **9 fresh NRQC phase requests**. It computes intervals
for actual continuum-model bounded measurement probabilities, conditional on
the frozen LMD/common-parent theorem. Its source approximation error includes
the omitted modes; it is not just a finite moment match.

## Accuracy and separation

The declared request uses h=1/8,1/4,3/8; t=-1/4,0,1/4,1/2; k=2; and
32,64,128 central Gaussian bins with two retained tail bins. The requested
maximum probability halfwidth is 1/200, or 0.5 percentage point.

At the finest partition, **8 of 12 cases meet that tolerance for both
cosine and sine readings**: five nonzero-pulse cases and all three exact
zero-pulse baselines. Four cases return valid but wider enclosures and are
explicitly marked `certified_enclosure_tolerance_not_met`.

The following are probabilities of the declared '+' postprocessing outcome,
not probabilities that the theory is true or that an experiment succeeds.
Every printed endpoint is rounded outwards.

| h | t | Actual cosine '+' interval | Actual sine '+' interval | Requested tolerance |
| --- | --- | --- | --- | --- |
| 1/8 | -1/4 | [0.880368744, 0.881274763] | [0.291428288, 0.292334307] | Met |
| 1/8 | 0 | [0.941248451, 0.941248452] | [0.500000000, 0.500000000] | Met |
| 1/8 | 1/4 | [0.880368744, 0.881274763] | [0.707665693, 0.708571712] | Met |
| 1/8 | 1/2 | [0.723750532, 0.727364357] | [0.844641372, 0.848255197] | Met |
| 1/4 | -1/4 | [0.827928675, 0.832037081] | [0.317734286, 0.321842692] | Met |
| 1/4 | 0 | [0.889400391, 0.889400392] | [0.500000000, 0.500000000] | Met |
| 1/4 | 1/4 | [0.827928675, 0.832037081] | [0.678157308, 0.682265714] | Met |
| 1/4 | 1/2 | [0.684321785, 0.700746361] | [0.774610662, 0.791035239] | Not met |
| 3/8 | -1/4 | [0.780280412, 0.791352807] | [0.338979191, 0.350051587] | Not met |
| 3/8 | 0 | [0.843644639, 0.843644640] | [0.500000000, 0.500000000] | Met |
| 3/8 | 1/4 | [0.780280412, 0.791352807] | [0.649948413, 0.661020809] | Not met |
| 3/8 | 1/2 | [0.645785189, 0.690066787] | [0.709046229, 0.753327826] | Not met |

For h=1/8,t=1/4,k=2, the actual cosine '+' interval is
**[0.880368744, 0.881274763]**, while the naive density-shift comparator is
**[0.872931701, 0.877491417]**. The complete certified intervals are disjoint;
this is a finite-time probability separation, extending the earlier variance
and linear-response diagnostics. Sine also separates at that point and at
the time-reversed point. At h=1/8,t=1/2 only the sine comparison separates.
Other nonzero cases remain unresolved against the naive comparator with this
compression budget. Zero pulse correctly agrees.

## What was actually executed

NRQC's pinned cyclotomic64 engine computed 34 exact Pauli expectations at
17 phase nodes in nine one/two-qubit requests. Native rational interval
arithmetic constructed Gaussian conditional masses, evaluated the certified
cubic phase interpolation and compared it with independent Taylor trig
evaluation. The pipeline included both Gaussian tails, checked nested
partition progress, exact baseline readings and time reversal, and detected
the wrong phase direction and omission of the combined Gaussian damping.

The scientifically new algorithmic step is the conditional-mean cancellation:
the complete source error is bounded by (a^2/2) times its residual variance.
That enables certified continuum expectation intervals without representing
every field mode. The NRQC phase probes are not a prepared field state or a
full pulse-evolution circuit. Their numerical contribution is explicit and
replaceable by the separately checked native phase calculation; no quantum
speed or memory advantage is claimed.

The native solver took approximately **272.916 seconds** on this host. The
nine NRQC command wall times ranged from approximately **0.656 to 11.069
seconds**. This is a serial local execution record, not a comparative or
scaling benchmark. Python 3.12.14 orchestrated the run; Node v24.22.1 executed
the native scientific arithmetic. No hardware, shots or CI dispatch occurred.

## Reproduce and inspect

Use the unmodified quantum-classical-provisional repository at commit
`fdd39048a6ddaa2dd78aa3197c274d7766f31a0e`:

```sh
python3.12 papers/native-lmd-certified-pulse-solver/run_solver.py \
  --nrqc /path/to/pinned/quantum-classical-provisional \
  --out /tmp/lmd-certified-pulse
```

The optional `--request` accepts a JSON request with the documented schema.
Supported h values are 1/8,1/4,3/8; |t|<=1, |k|<=4, radius=2, and central
bin counts from 32,64,128. The default request is retained in request.json.
Invalid/unsupported choices fail instead of silently changing the problem.

To re-evaluate native processing using the saved calculator nodes:

```sh
node papers/native-lmd-certified-pulse-solver/solve.cjs \
  papers/native-lmd-certified-pulse-solver/results
```

- [Checks and compact results](results/verification.json)
- [Every interval, error budget and partition refinement](results/solver.results.json)
- [Execution metadata and source hashes](results/execution.json)
- [Artifact manifest](results/manifest.json)
- [Conditional continuum proof and algorithm](THEOREM.md)
- [Review and limits](REVIEW.md)

For the four wide cases, more central q bins cannot remove the unresolved
nonzero-mode variance. The next accuracy step is a source compression that
includes nonzero Fourier modes. General states, nonuniform/joint geometry
measurements and autonomous physical gravity are still separate tasks.
