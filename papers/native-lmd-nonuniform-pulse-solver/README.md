# First certified nonuniform-source pulse marginal

Author: Monty. 31 October 2026.

The native solver now handles a positive spatially varying profile
w(z)=1+e cos(2 pi z), retaining its relative angle with the first source
mode. A modified-Bessel source term and an ordinary-Bessel angular average
preserve the modulation's contribution. Higher spatial modes remain in the
complete error budget.

The frozen run passed **266 checks and 11 negative controls** in **351.998 s**.
Eight of nine cases meet the probability halfwidth 1/200 at both refinements.
This includes four of five primary nonuniform cases and four controls or
references. The h=3/8,t=1/2,e=1/2 case remains wider and is retained explicitly.

At h=1/8,t=1/2,k=2, the strong e=3/4 profile is distinguished from e=0 by
certified sine '+' probabilities:

| Source | Probability enclosure |
|---|---|
| Uniform e=0 | [0.844775485,0.846462258] |
| Nonuniform e=3/4 | [0.840676900,0.843121654] |

These intervals are disjoint after the complete source and evaluation error
budgets. The total integrated source strength and mean label response agree;
its spatial distribution changes the bounded reading. This is a result in
the specified mathematical model, not a hardware observation.

- [Derivation and error budget](THEOREM.md)
- [Execution, both refinements and unresolved case](EXECUTION.md)
- [Source map and credit](SOURCE_MAP.md), [inherited hashes](SOURCE_PINS.json)
- [Frozen request](request.json), [raw results](results/solver.results.json),
  [verification](results/verification.json)
- [Post-run exact joint-law derivation and next implementation](EXACT_JOINT_LAW.md)

The last link adds a separate analytic route using the credited
Chhaibi--Najnudel measure-valued theorem. It preserves the dependence between
mass and first harmonic and includes all spatial modes for this source
family. That representation's finite-cell solver has not been executed here;
it does not retrospectively close the numerical failure in this campaign.

## Reproduce

Use Python 3.12 and Node from the repository root, with a new output directory:

```sh
python3.12 papers/native-lmd-nonuniform-pulse-solver/run.py \
  --out /tmp/lmd-nonuniform-replay \
  --protocol-commit bb04408a16c52b165872eaea95c915fad7795c76
```

The protocol was published before probability evaluation. The runner checks
inherited hashes and records execution-source hashes and output manifests.
Scientific arithmetic uses the existing native rational engine; Python only
orchestrates. Nine pinned NRQC phase outputs are reused, with zero new NRQC,
hardware or CI requests. The previous uniform 12/12 certificate is unchanged.
Local or joint geometry, autonomous gravity and physical 3+1 dynamics remain
separate development obligations.
