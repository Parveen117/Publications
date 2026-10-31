# Certified uniform-pulse marginal from the full source mass law

Author: Monty. 31 October 2026.

The remaining h=3/8,t=1/2,k=2 probability case now meets the frozen 1/200
halfwidth. Its final cosine '+' enclosure is [0.670953808,0.671586342] and
sine '+' enclosure is [0.721602998,0.722235532]. The halfwidth is at most
0.031626657 percentage points, compared with about 1.30166 in the first-mode
predecessor. All decimal endpoints are rounded outward from rational bounds.

The new route uses Rémy's full circular GMC mass theorem as an explicitly
credited analytic dependency, with the source covariance and normalization
audited in [SOURCE_AUDIT.md](SOURCE_AUDIT.md). It retains the existing zero
mode, source model, pulse, detector and tolerance. The theorem supplies a
scalar law including every spatial mode's contribution to total mass.
It does not implement an eight-mode circuit or a full field-state solver.

- [Derivation and complete error budget](THEOREM.md)
- [Executed results, controls and provenance](EXECUTION.md)
- [Frozen request](request.json), [inherited hashes](SOURCE_PINS.json),
  [raw results](results/solver.results.json), [verification](results/verification.json)
- [Nonuniform-source next gate and exact mismatch witness](NEXT_STEP.md)

The frozen campaign passed 239 checks and 10 negative controls in 87.976 s
of native solver wall time. Four cases and two refinements were evaluated.
Eleven original cases were already certified in two frozen packets; this
packet adds the last one, giving 12/12 at the original accuracy. That count
combines retained evidence, rather than twelve fresh executions.

## Reproduce

Use Python 3.12 and Node, from the repository root. The output directory must
not already exist. No third-party scientific arithmetic package is used.

```sh
python3.12 papers/native-lmd-exact-mass-marginal/run.py \
  --out /tmp/lmd-exact-mass-replay \
  --protocol-commit 5fde12e85080a1f086295c8234d2922d1d27271f
```

The protocol commit is the recorded source provenance; the runner also
checks the inherited hashes and records its own execution-source hashes.
Compare these and the output manifest with the retained run. Scientific
arithmetic uses the existing native rational engine with directed 192-bit
intervals. Python handles orchestration, file integrity and timing only.

Nine pinned NRQC phase outputs are reused; there are zero fresh NRQC or
hardware calls. The certificate is conditional on the stated mathematical
dependencies. Nonuniform/joint readouts, autonomous gravity and physical
3+1 dynamics are not established by these selected marginal probabilities.
