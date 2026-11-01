# Executed first-spatial-mode refinement

Author: Monty. 31 October 2026.

Three of the four previously unresolved uniform-vacuum pulse cases now meet
the unchanged probability halfwidth target **1/200**. Combined with the
eight retained predecessor successes, **11 of the original 12 cases** have
certificates at that accuracy. Only four cases were newly evaluated here.
The result is conditional on the same LMD/common-parent source theorem;
the source, pulse and detector definitions are unchanged.

The new run passes **234 checks and 8 negative controls**. Its probability
error includes all omitted field modes, q/radial quantization, all tails,
phase interpolation and arithmetic enclosures. The last case retains a
valid but wider enclosure and is explicitly marked tolerance-not-met.

## Complete probability intervals

The '+' events are the existing bounded cosine/sine postprocessings of v_0.
They are not theory-confidence scores, threshold probabilities P(v_0>0),
hardware frequencies or a proof of full quantum gravity. All cases use k=2.

| h | t | New cosine '+' interval | New sine '+' interval | Halfwidth, percentage points (approximately) | Target met? |
| --- | --- | --- | --- | --- | --- |
| 1/4 | 1/2 | [0.688453864, 0.697017226] | [0.775845580, 0.784408942] | 0.42817 | Yes |
| 3/8 | -1/4 | [0.781790228, 0.788300514] | [0.342806983, 0.349317269] | 0.32552 | Yes |
| 3/8 | 1/4 | [0.781790228, 0.788300514] | [0.650682731, 0.657193017] | 0.32552 | Yes |
| 3/8 | 1/2 | [0.656330954, 0.682364115] | [0.713878700, 0.739911861] | 1.30166 | No |

Interval endpoints are rounded outward. Exact rational halfwidths are retained
in [solver.results.json](results/solver.results.json). The previous halfwidths
were respectively 0.82123, 0.55362, 0.55362 and 2.21408 percentage points.
The new error widths fall by about 47.86%, 41.20%, 41.20% and 41.21%.
These are reductions of certified uncertainty widths, not measurements of
prediction accuracy against experimental truth.

## What changed and what was reused

The first spatial sine/cosine pair is represented by its radius. Its absolute
orientation integrates out exactly for uniform total mass; q plus one radial
variable therefore suffice for this conditional expectation. A positive
Poisson/Gamma series gives exact conditional radial-cell masses. The two
declared product partitions contain 1,122 and 4,290 cells, including both q
tails and the radial tail. Both refinements and all cell masses are retained.

The implementation, theorem and fixed request were committed at
[`917d68f`](https://github.com/Parveen117/Publications/commit/917d68f28da80bb7d4a8c8690ce43b686c0b4d7c)
before the probability evaluation. The old packets remain frozen at
`13389614767542792097c189dcf6e11747806b42`.

There were **zero fresh NRQC requests**. The nine previous calculator outputs
at NRQC commit `fdd39048a6ddaa2dd78aa3197c274d7766f31a0e` are reused, including
their exact-backend statuses and outward enclosures. Seventeen phase nodes
are checked again against native Taylor trig; seven selected cell phases
per case/refinement are checked independently. The whole new assembled sum
was not independently replayed by a second solver. Certification rests on
the stated conditional remainder and interval arithmetic.

Native Node v24.22.1 scientific processing took **82.73 seconds**. Python
3.12.14 performed orchestration and hashing only. Runtime is not a controlled
comparison with the predecessor: the source dimension, set of cases and
independent-trig checking workload differ. No speed or memory advantage is
claimed. No field state, full pulse circuit, hardware job or CI dispatch was
executed. The output hash manifest and final execution-source hashes passed
the final audit.

## The remaining boundary

For h=3/8,t=1/2,k=2, even ideal first-mode conditioning leaves a term above
0.01266209 in this quadratic error certificate, versus target 0.005.
Finer bins alone cannot remove that term. It is a certificate limitation,
not a lower bound on the true unknown oscillatory error.

The post-run [next-step derivation](NEXT_STEP.md) uses previously certified
mode moments: N=8 is the first tested Fourier cutoff in {0,1,2,4,8,16} with
room for source-bin error. It also proves that conditioning on **all mode
amplitudes but discarding their angular relations** leaves a budget above
0.00635767 for this same certificate. A two-harmonic mass witness explicitly
depends on the invariant relative phase alpha_2-2 alpha_1. The next source
compression must retain such relations or supply a sharper error theorem.
An executable 8-mode expectation solver is still open.

The prior actual-versus-naive probability separation remains a predecessor
result. The naive comparator was not recalculated here. Nonuniform/joint
measurements, general states, detector-flow covariance, fixed-h smoothing
removal, pulse-evolved length domains, autonomous gravity and physical 3+1
recovery retain their original open status.

## Reproduction

```sh
python3.12 papers/native-lmd-first-mode-pulse-solver/run.py \
  --out /tmp/lmd-first-mode-reproduction \
  --protocol-commit 917d68f28da80bb7d4a8c8690ce43b686c0b4d7c
node papers/native-lmd-first-mode-pulse-solver/verify.cjs
node papers/native-lmd-first-mode-pulse-solver/frontier.cjs
```

The output directory must not exist. Source-pin mismatches and incomplete
execution fail closed. The third command prints the separate post-run
mode-budget diagnostic, retained in [MODE_BUDGET.json](MODE_BUDGET.json).

- [Conditional first-mode theorem](THEOREM.md)
- [All cells, refinements and probability enclosures](results/solver.results.json)
- [Checks, controls and achieved accuracies](results/verification.json)
- [Execution provenance](results/execution.json) and [hash manifest](results/manifest.json)
- [Source map](SOURCE_MAP.md) and [internal review](REVIEW.md)
