# Frozen nonuniform-source execution

Author: Monty. 31 October 2026.

**Result:** eight of nine cases meet the original 1/200 probability-halfwidth
target at both tested refinements. Four of the five primary nonuniform
cases pass; the four controls/references also pass. The difficult
h=3/8,t=1/2,k=2,e=1/2 case returns a valid wider enclosure.

## Scope and provenance

| Item | Recorded value |
|---|---|
| Frozen source parent | `ca080861d646774ad94a8e48be668194a8086d59` |
| Protocol published before probability evaluation | `bb04408a16c52b165872eaea95c915fad7795c76` |
| Original NRQC phase-data commit | `fdd39048a6ddaa2dd78aa3197c274d7766f31a0e` |
| Implementation checks / negative controls | 266 / 11, all pass |
| Cases / refinements per case | 9 / 2 |
| Product cells per case | 1,122 coarse; 4,290 fine |
| Complete native solver wall time | 351.998 s |
| Environment | Python 3.12.14; Node v24.22.1 |
| Fresh / reused NRQC requests | 0 / 9 |
| Hardware / CI requests | 0 / 0 |
| Inherited files pinned | 37 |

This is a first-mode conditional expectation with a retained angular integral
and a complete continuum error bound under the frozen LMD/common-parent
contract. The wider result is a declared tolerance limitation, not an
arithmetic-check failure. Runtime applies to this campaign only; no comparison
or speed advantage is inferred.

## Fine-grid probabilities

All rows use k=2. Decimal endpoints are outward-rounded from retained rational
enclosures. Halfwidths below are upward-rounded bounds in percentage points;
the required limit is 0.5 pp. Both refinements and all intermediate source
cells remain in [solver.results.json](results/solver.results.json).

| Case | h | t | e | Cosine '+' | Sine '+' | Halfwidth, pp |
|---|---|---|---|---|---|---|
| Small h | 1/8 | 1/2 | 1/2 | [0.723965769,0.725989422] | [0.842946228,0.844969882] | <=0.101182605 |
| Middle h | 1/4 | 1/4 | 1/2 | [0.827648313,0.830159456] | [0.677588833,0.680099976] | <=0.125557093 |
| Large h | 3/8 | 1/4 | 1/2 | [0.780841882,0.788362647] | [0.649157348,0.656678113] | <=0.376038208 |
| Difficult, tolerance unmet | 3/8 | 1/2 | 1/2 | [0.655315575,0.685390653] | [0.709256644,0.739331722] | <=1.503753890 |
| Strong modulation | 1/8 | 1/2 | 3/4 | [0.723286528,0.725731281] | [0.840676900,0.843121654] | <=0.122237628 |
| Amplitude sign | 1/8 | 1/2 | -1/2 | [0.723965769,0.725989422] | [0.842946228,0.844969882] | <=0.101182605 |
| Reverse pulse | 1/8 | -1/2 | 1/2 | [0.723965769,0.725989422] | [0.155030118,0.157053772] | <=0.101182605 |
| Zero pulse | 1/8 | 0 | 1/2 | [0.941248451,0.941248452] | Exactly 1/2 | <0.000000001 |
| Uniform reference | 1/8 | 1/2 | 0 | [0.724505485,0.726192258] | [0.844775485,0.846462258] | <=0.084338593 |

The sine '+' difference between the strong-modulation and uniform rows is
[-0.005785357,-0.001653832]. Thus the prescribed spatial modulation is
detected by a bounded reading after including every error term. Their cosine
intervals overlap. At e=1/2 both sine and cosine comparisons overlap; no
separation is claimed for that weaker modulation. The uniform row overlaps
its frozen predecessor certificate. These are deliberately retained positive
and unresolved comparisons, not a selection of successful observables.

## Error diagnosis and controls

The difficult case's halfwidth improves from at most 1.572129375 pp on the
coarse grid to 1.503753890 pp on the fine grid. Its *ideal first-mode*
probability-budget term is enclosed by [0.014603524,0.014603528], already
above the target 0.005. Refining only the q/radial cells cannot remove this
obstacle in the present variance certificate. It is not a lower bound on
the unknown actual approximation error.

Every source probability and tail is retained. The ordinary-Bessel
large-argument fallback was available by protocol but used by **zero cells**
in every case and refinement. The full moment, captured moment, decomposition,
phase remainder, Bessel series remainder and source weights are recorded.

Checks cover half-Gamma normalization and recurrence; independent Bessel
angular integration; positive weighted cells and tails; nested captured
moments; calculator nodes and selected base phases against Taylor trig;
source normalization; amplitude sign; pulse reversal; zero pulse; and
uniform predecessor overlap. The eleven negative controls detect a dropped
tail, fixed aligned angle, deleted relative-angle contribution, confusion
between ordinary and modified Bessel functions, a missing angular factor,
an incorrect Gamma shape, uniform-law substitution at the moment level,
unsupported partition/endpoint, pulse reversal and the explicit weighted
second-moment mismatch.

The theoretical same-mean/different-variance identity follows from the
written model derivation. The tests above validate its implementation
ingredients and the selected bounded probabilities; they do not prove the
continuum construction or a physical apparatus.

## Post-run analytic next step

The separate [exact joint-law derivation](EXACT_JOINT_LAW.md) audits and applies
Chhaibi--Najnudel's measure-valued theorem to this source. It derives a
joint auxiliary representation of total mass J and its first Fourier mass L,
preserving their dependence, plus a concrete positive-cell error formula.
This is a new credited analytic dependency introduced after the frozen run;
the recorded code and probabilities do not use it. Its implementation and
runtime qualification are pending. It therefore supplies a route past the
first-mode budget without relabeling the difficult run as successful.

[Metadata](results/execution.json), [verification](results/verification.json)
and the [manifest](results/manifest.json) retain the protocol and evidence.
All executed source hashes, inherited hashes and output hashes were checked
after completion. The previous uniform certificates remain unchanged.
Joint detector/geometry laws, calibration, operator domains, autonomous
gravity, physical 3+1 dynamics and hardware implementation remain open.
