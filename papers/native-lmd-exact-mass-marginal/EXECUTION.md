# Frozen exact-mass-law marginal execution

Author: Monty. 31 October 2026.

**Result:** all four frozen cases meet the 1/200 probability-halfwidth target
at both tested refinements. The last unresolved original accuracy case is
closed under the existing LMD/common-parent contract plus the explicitly
credited circular GMC law in [SOURCE_AUDIT.md](SOURCE_AUDIT.md).

## Provenance and numerical scope

| Item | Recorded value |
|---|---|
| Frozen source parent | `ace56a0a3faa21708e4a6779c8e4a8f02dab70e0` |
| New protocol commit, published before target evaluation | `5fde12e85080a1f086295c8234d2922d1d27271f` |
| Original phase-result calculator commit | `fdd39048a6ddaa2dd78aa3197c274d7766f31a0e` |
| Native checks / negative controls | 239 / 10, all pass |
| Runtime | 87.976 s for this complete native campaign |
| Environment | Python 3.12.14; Node v24.22.1 |
| New / reused NRQC requests | 0 / 9 |
| Hardware calls / CI dispatch | 0 / 0 |
| Coarse / fine product cells per case | 1,428 / 5,412 |
| Requested probability halfwidth | 1/200 = 0.5 percentage points |

The run uses an exact scalar distribution identity, then certified finite
conditional expectations. It does not encode the field state in the calculator.
The wall time is for this workload only and establishes no speed advantage
over another simulator or physical device.

## Retained probability enclosures

All endpoints below are outward decimal displays of the full rational
enclosures in [solver.results.json](results/solver.results.json). The halfwidth
column is an upward-rounded bound, in percentage points, for both readings.

| Case (h,t,k) | Partition | Cosine '+' | Sine '+' | Halfwidth, pp |
|---|---|---|---|---|
| Remaining (3/8,1/2,2) | Coarse | [0.670565593,0.671897269] | [0.721417451,0.722749127] | <=0.066583740 |
| Remaining (3/8,1/2,2) | Fine | [0.670953808,0.671586342] | [0.721602998,0.722235532] | <=0.031626657 |
| Reverse (3/8,-1/2,2) | Fine | [0.670953808,0.671586342] | [0.277764468,0.278397002] | <=0.031626657 |
| Predecessor overlap (1/8,1/4,2) | Fine | [0.880472012,0.880517910] | [0.707837104,0.707883003] | <=0.002294868 |
| Zero pulse (3/8,0,2) | Fine | [0.843644639,0.843644640] | Exactly 1/2 | <0.000000001 |

The remaining case's first-mode predecessor gave cosine
[0.656330954,0.682364115] and sine [0.713878700,0.739911861], with about
1.30166 pp halfwidth. The new halfwidth is approximately 41.157 times smaller.
Both new refinements overlap that retained certificate. The independent
small-h result overlaps its older zero-mode certificate. Reverse time has
the same cosine and opposite sine characteristic component; the zero-pulse
reading agrees with its analytic Gaussian value.

Eight original certifications and three first-mode certifications remain
unchanged. Adding the one remaining original case gives **12/12** original
cases at the declared accuracy. Four cases were executed here, including
controls and overlap cases; the combined count is not a fresh twelve-case run.

## Complete error and controls

For h=3/8 the residual source variance shrinks from
[0.003691111,0.003691112] to [0.001749364,0.001749365]. At the fine target,
the probability-budget contributions are bounded by 0.000270743 from q
binning and 0.000044195 from E binning, including both variables' retained
tails. The inherited phase-component interpolation remainder is below
0.000003870690 before damping and probability conversion. Native arithmetic,
Gamma integration and phase-node enclosures propagate through the final
rational endpoints. No spatial cutoff error remains for this scalar law,
conditional on the external theorem and its audited applicability.

Checks cover independent Gamma and incomplete-Gamma recurrences, analytic
shape-one values, mass/probability normalization, an independent frozen
second moment, nested conditional variances, retained calculator phase
nodes, selected cell phases against Taylor trig, predecessor overlap, time
reversal and the zero pulse. They are numerical implementation checks;
the written theorems carry the continuum claims.

The ten negative controls reject dropping the singular E endpoint or its
unbounded tail; missing Gamma normalization or zero-mode Wick subtraction;
the wrong coupling exponent; invalid Gamma shape, fractional power or
unfrozen partition; reversal of pulse direction; and using only two matched
moments to infer the full probability law. In the last control, a two-point
law with the same first two mass moments gives cosine '+' probability
[0.616279225,0.616280913], disjoint from the certified target. Therefore the
external full-law input cannot be replaced by a moment fit.

## Integrity and open boundary

[Execution metadata](results/execution.json) records protocol identity,
environment, wall time and the executed source hashes. [SOURCE_PINS.json](SOURCE_PINS.json)
pins 34 inherited inputs. [Manifest](results/manifest.json) covers the retained
JSON results; all source and output hashes were checked after completion.
No predecessor packet was changed, and no campaign failure was discarded.

This closes the selected uniform-vacuum marginal accuracy gate. The
[next-step note](NEXT_STEP.md) gives a concrete nonuniform pulse and proves
why its mass law differs. It has not been numerically executed here. Joint
geometry, current laws, calibration covariance, domains, fixed-h smoothing
removal, autonomous gravity and physical 3+1 dynamics remain open.
