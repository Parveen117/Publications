# LAA-1 actual-source evaluation

The completed response, original coordinate s, and relation
a=(s-1/2)^2 are unchanged. The pinned LAT jet at s=29 is reused.
The unchanged ILC evaluator is executed afresh at s=31 and s=33
through derivative order 19, using the same complete-source protocol:
prime cutoff 4096, Gamma shift 512, 32 Euler--Maclaurin terms, and
the inherited 110-place outward rational arithmetic. Every derivative
tail and source component is retained in the certificate.

Both 10x10 source forms are positive at all three anchors. All 210
binomial readings through total degree 19 are positive at each one.
Every moment enclosure has width below 1e-62. The new s=31 value
and derivative independently overlap the previously unfitted PIW
source calculation. The composite source label s=33 is valid: these
confluent derivative forms do not require a prime anchor label.

## Whole-interval source certificates

The general-anchor analytic majorants Q_w are exactly 221/228,
225/244 and 229/260 at s=29,31,33 respectively. The normalized
majorant uses the lower endpoint of the corresponding positive q(a).

| Original source anchor | a | Order 3, d=3/20 | Order 4, d=1/12 |
|---|---|---|---|
| 29 | 3249/4 | [55233/80,74727/80] | [11913/16,14079/16] |
| 31 | 3721/4 | [63257/80,85583/80] | [40931/48,48373/48] |
| 33 | 4225/4 | [14365/16,19435/16] | [46475/48,54925/48] |

Each row is certified by positive interval LDL factors for
A_m-tau_A I and B_m-tau_B I, plus the normalization floor.
Both orders have strictly overlapping neighboring intervals.

| Fixed order | Exact connected a interval | Included real s interval | Inherited cut channels |
|---|---|---|---|
| 3 | [55233/80,19435/16] | [27,35] | 6 |
| 4 | [11913/16,54925/48] | [28,34] | 8 |

The overlap certificates refer to the same actual source, not to
equality of independently fitted finite rational functions. The
native checks verify endpoints and strict overlaps. LAA-T3 supplies
the uniform analytic argument between them. Ten-dimensional
positivity throughout these intervals is not claimed.

## Transport error and legitimate restart

The ledger 29->31->33 gives exactly the same polynomial enclosure
and budget as direct 29->33, retaining the s=29 source and degree 19.
All twenty normalized forecast intervals contain the fresh s=33
source moments. The new exact tail formula gives bounds even for
the highest retained derivative; a bounded interval need not be
precise enough to certify a form.

For the unnormalized fifth coefficient at s=33:

| Source used for forecast | Uniform remainder upper bound |
|---|---|
| Original s=29 jet, either path | <0.356 |
| Fresh s=31 jet | <8.251e-7 |

All twenty budgets decrease after that legitimate fresh-source
restart. The order-19 bounds still exceed one in both cases. This
is the explicit precision cost: retaining degree 19 algebraically
does not preserve useful accuracy at every layer. The shifted
polynomial's value interval is disjoint from the actual s=33 value,
so resetting its source error to zero would be demonstrably wrong.

An inverse path returns to the original coordinate and zero analytic
remainder, because the ledger retains the original source rather
than successively rounding and replacing it. A path of four locally
small factors 5/4 is rejected when its cumulative ratio exits the
original analytic disk. A missing twenty-first layer is rejected.

## Independent held-out source at s=35

The unchanged PIW full_source.cjs evaluates q and q' at s=35 without
using either in the three finite fits. Both readings lie in the
s=33 source forecast. All three ten-direction responses have
positive errors and satisfy the inherited unconditional ILC bound.

| Anchor of finite rational fit | Approximate positive error at s=35 |
|---|---|
| 29 | 7.7505380453e-26 |
| 31 | 2.2148786438e-29 |
| 33 | 1.9706735470e-35 |

The certificate stores full outward intervals; the displayed values
are only rounded descriptions. The final error is certified below
1.971e-35. The three response intervals are disjoint. Hence fresh
local fits are not three anchor descriptions of one fixed rational
model. LAT transport preserves the response of a fixed model;
reevaluating and fitting the actual source changes that model.

The all-order actual-source inequality remains open. No extrapolation
from the finite atlas proves RH or identifies a physical wave, mass,
time evolution, Maxwell field or quantum-gravity theory.
