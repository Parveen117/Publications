# RH positivity at L=3/4, 4/5, 17/20, 9/10, 19/20 and 1: reproducible fixed-source certificates

Monty · 23 October 2026, updated 24 October 2026

This packet supports [BRT-T64--T74](../../RH_BOUNDARY_REGULAR_TRIAL.md#brt-t64-a-width-dependent-complement-comparison-with-boundary-energy-allocation).
For the real even closed source family declared there, with every
active prime power carrying its actual weight Lambda(n)/sqrt(n)
(2, 3, 4 at L=3/4 and 4/5; 2, 3, 4, 5 at L=17/20, 9/10 and 19/20;
2, 3, 4, 5, 7 at L=1), the results are

\[
0<L\le L_0 \ \Longrightarrow\ q_L[g]\ge\mu_{L_0}\|g\|^2
\quad(g\in D_L),\qquad
L_0\in\{3/4,\ 4/5,\ 17/20,\ 9/10,\ 19/20,\ 1\},
\]

with the explicit margins recorded in each certificate
(`full_family_lower_bound_computed`; `full_family_lower_bound` is the
declared value the verifier asserts against it). It is a bounded-support, computer-assisted
working proof. All-width positivity and RH remain open. The packet is
not an independent external review or a formal proof-assistant
development.

## Layout

One verifier, several packet directories:

| Directory | Run | Role |
|---|---|---|
| `t53_route/` | `--packet t53_route` | L=3/4, N=8 retained modes, complement hierarchy from BRT-T53/T54 (`declared` rule) |
| `boundary_route/` | `--packet boundary_route` | L=3/4, N=14 retained modes, complement hierarchy from BRT-T64 (`boundary_split` rule); independent of the T41--T53 prime-graph comparison |
| `near_null/` | `--packet near_null --rayleigh` | L=3/4, two-sided enclosure of the Rayleigh quotient of three fixed profiles (BRT-T65); not a positivity statement |
| `l08_route/` | `--packet l08_route` | L=4/5, N=16 retained modes, thirty-mode trial, kernel degree 300, `boundary_split` rule (BRT-T69) |
| `near_null_08/` | `--packet near_null_08 --rayleigh` | L=4/5, two-sided enclosure of three fixed profiles (BRT-T69); not a positivity statement |
| `l085_route/` | `--packet l085_route` | L=17/20, four active prime powers, N=24 retained modes, thirty-six-mode trial, kernel degree 400, `boundary_split` rule (BRT-T70) |
| `near_null_085/` | `--packet near_null_085 --rayleigh` | L=17/20, two-sided enclosure of three fixed profiles (BRT-T70); not a positivity statement |
| `l09_route/` | `--packet l09_route` (Arb engine) | L=9/10, N=30 retained modes, fifty-two-mode trial, kernel degree 700, `boundary_split` rule (BRT-T71) |
| `near_null_09/` | `--packet near_null_09 --rayleigh` | L=9/10, two-sided enclosure of three fixed profiles (BRT-T71); not a positivity statement |
| `l095_route/` | `--packet l095_route` (Arb engine) | L=19/20, N=40 retained modes, sixty-four-mode trial, kernel degree 300 with the radius-5/2 remainder bound of BRT-T72, `boundary_split` rule (BRT-T73) |
| `near_null_095/` | `--packet near_null_095 --rayleigh` (Arb engine) | L=19/20, two-sided enclosure of three fixed profiles (BRT-T73); not a positivity statement |
| `l1_route/` | `--packet l1_route` (Arb engine) | L=1, five active prime powers 2, 3, 4, 5, 7, N=56 retained modes, eighty-mode trial (denominator 10^18), kernel degree 400 with the BRT-T72 bound, fourteen cosh terms, `boundary_split` rule (BRT-T74) |
| `near_null_1/` | `--packet near_null_1 --rayleigh` (Arb engine) | L=1, two-sided enclosure of three fixed profiles (BRT-T74); not a positivity statement |
| `l0625_route/` | `--packet l0625_route` | L=5/8, N=5, `boundary_split` rule; two-prime control for the same rule |
| `l055_control/` | `--packet l055_control` | L=11/20, N=5, the `rh_l055` trial and declared T57 constant re-run through this verifier; reproduces the exact rational witness of that packet |

Each directory holds `config.json` (width, active prime powers, retained
modes, kernel and cosh degrees, complement rule, margins, and from
L=19/20 on the `kernel_tail` pair of BRT-T72), `trial.json`
(fixed rational columns) and the committed `result-<precision>.json`
(or `rayleigh-<precision>.json`) from `verify.py` and
`result-arb-<precision>.json` (or `rayleigh-arb-<precision>.json`)
from `verify_arb.py`.

## Reproduce

The evidence was generated with Python 3.13 using only its standard
library (the L=11/20 packet used 3.12; the arithmetic is specified by
the `decimal` module, not by the interpreter version). Run from the
repository root:

~~~sh
python3 proof_lab/rh_l075/verify.py 220 --packet proof_lab/rh_l075/t53_route
python3 proof_lab/rh_l075/verify.py 220 --packet proof_lab/rh_l075/boundary_route
python3 proof_lab/rh_l075/verify.py 160 --packet proof_lab/rh_l075/near_null --rayleigh
python3 proof_lab/rh_l075/verify.py 160 --packet proof_lab/rh_l075/l0625_route
python3 proof_lab/rh_l075/verify.py 220 --packet proof_lab/rh_l075/l055_control
python3 proof_lab/rh_l075/verify.py 220 --packet proof_lab/rh_l075/l08_route
python3 proof_lab/rh_l075/verify.py 160 --packet proof_lab/rh_l075/near_null_08 --rayleigh
python3 proof_lab/rh_l075/verify.py 230 --packet proof_lab/rh_l075/l085_route
python3 proof_lab/rh_l075/verify.py 160 --packet proof_lab/rh_l075/near_null_085 --rayleigh
python3 proof_lab/rh_l075/verify.py 160 --packet proof_lab/rh_l075/near_null_09 --rayleigh
~~~

The L=9/10, 19/20 and 1 certificates and the 19/20 and 1 near-null
enclosures were produced with the Arb engine (`pip install python-flint`):

~~~sh
python3 proof_lab/rh_l075/verify_arb.py 260 --packet proof_lab/rh_l075/l09_route
python3 proof_lab/rh_l075/verify_arb.py 320 --packet proof_lab/rh_l075/l095_route
python3 proof_lab/rh_l075/verify_arb.py 200 --packet proof_lab/rh_l075/near_null_095 --rayleigh
python3 proof_lab/rh_l075/verify_arb.py 420 --packet proof_lab/rh_l075/l1_route
python3 proof_lab/rh_l075/verify_arb.py 300 --packet proof_lab/rh_l075/near_null_1 --rayleigh
~~~

`verify.py` accepts the `l09_route` packet (it would need of the order
of days). It does not read the `kernel_tail` key: on `l095_route` and
`near_null_095` it would apply T60's radius-2 remainder bound, whose
value at degree 300 and ratio 19/20 is of the order 10^(-5), far above
the 10^(-30) shift, and it asserts L<1, which rejects `l1_route` and
`near_null_1`; those four packets are certified by the Arb engine
only. Every other packet runs on both engines.

Run without -O or PYTHONOPTIMIZE: assertions are part of the verifier.
Each command reads the packet's `config.json` and `trial.json` and writes
the result file next to them. Copy the packet to a temporary directory
first if you want to preserve the committed evidence byte for byte.
The two L=3/4 certificate runs take of the order of half an hour each
on one core, the L=4/5 run of the order of two hours and the L=17/20
run of the order of five hours; the others take minutes. On the Arb
engine the L=9/10 run takes about sixteen minutes, the L=19/20 run
about ten minutes at 320 digits and the L=1 run about twenty-three
minutes at 420 digits; the rest take seconds.

A successful certificate run exits with status zero, prints
PASS_FIXED_WINDOW_LOWER_MATRIX, reports both the interval shifted sign
test and the exact rational shifted sign test as true, and reports all
rejection controls as true. A failed assertion is failure or
inconclusive verification, not evidence that the actual mathematical
source is negative.

## Fixed choices

The code SHA-256 is recorded in every result file as `code_sha256`;
the trial and config hashes as `trial_sha256` and `config_sha256`.

Trials: `t53_route` and `boundary_route` both use 24 unnormalized even
Legendre modes P_0, P_2, ..., P_46 of x/L, denominator 10^14; the first
N rows are the identity. `l08_route` uses 30 modes P_0, ..., P_58 with
N=16, `l085_route` 36 modes P_0, ..., P_70 with N=24, `l09_route` 52
modes P_0, ..., P_102 with N=30, `l095_route` 64 modes
P_0, ..., P_126 with N=40, and `l1_route` 80 modes P_0, ..., P_158
with N=56 and denominator 10^18 (at the 10^(-30) near-null scale of
that width a 10^(-14) rounding of the tail would cost more than the
margin). The remaining
rows are a fixed rational rounding of the truncated minimizing tail of
the 24-, 30-, 36-, 52-, 64- or 80-mode operator matrix found by
preliminary numerical candidate exploration. That
exploration is not a proof dependency: the verifier accepts the fixed
trial and charges every residual mode. A different trial can fail the
sufficient certificate without disproving positivity.

Coordinates use P_(2j)(x/L) without normalization; basis squared norms
are 2L/(4j+1). The certified lower-matrix margin (`shift`) is in those
coordinates; its sign transfers by congruence to the normalized matrix.
The all-function bound `full_family_lower_bound` is stated for
||g||_2^2 and is the quantity to quote.

## What is enclosed

The operator is BRT-T60's (236) with the subtraction
sum_n a_n (S_(log n)+S_(-log n))/2 over every active prime power. Only
the regular kernel j and cosh(x/2) are replaced, by a degree-200
(L=3/4), degree-300 (L=4/5 and 19/20), degree-400 (L=17/20 and 1) or
degree-700 (L=9/10) Taylor polynomial and twelve (fourteen at L=1)
even Taylor terms; BRT-T66, T69, T70, T71, T73 and T74 bound the
combined operator error epsilon. The kernel remainder uses T60's
radius-2 Cauchy bound (geometric ratio L) up to 9/10 and, from 19/20
on, T72's radius-5/2 bound |j|<=16/5 (ratio 4L/5), selected by the
packet key `kernel_tail` = {radius, sup} with the two admitted pairs
(2, 5) and (5/2, 16/5).
Logarithmic endpoint integrals, the square of the boundary multiplier,
and every mixed prime term are integrated analytically on the pieces
of y=(x+L)/(2L) cut at s_n'=log n/(2L) and 1-s_n' for every active n:
seven pieces at L=3/4 and 4/5, nine at L=17/20, 9/10 and 19/20,
eleven at L=1. The prime-2 translation caps overlap at every width;
the piecewise integration does not use cap disjointness. The verifier
classifies every prime power below 64 against 2L with directed
arithmetic: at L=4/5 log 5 > 8/5, so the prime power 5 is not active;
at L=17/20, 9/10 and 19/20 log 5 < 2L < log 7, so exactly 2, 3, 4, 5
are active; at L=1 log 7 < 2 < log 8, so exactly 2, 3, 4, 5, 7 are
active.

The residual Gram is the complete one,
(bar T Z)^*(bar T Z) - (P_N bar T Z)^*(P_N bar T Z). The complement
inverse is charged level by level: explicit levels N, ..., N+h-1 with the
proved weights 1/d_j, and every higher level with the single weight
1/d_(N+h), using the monotone hierarchy d_(j+1) > d_j. For `t53_route`,
d_j = 2429/128000 + (H_(2j)-H_16)/2 (BRT-T53 (213)--(216)); for
`boundary_route` and `l0625_route`, d_j is computed by the verifier from
BRT-T64 with directed arithmetic and recorded in `complement_report`.

The directed entry enclosure gives a rational symmetric matrix Q; exact
Fraction LDL^T proves Q-(rounding allowance + shift)I positive. A
separate interval LDL^T test also passes. The evidence includes Q's
integer numerators and the exact rational pivots; no floating eigenvalue
is accepted as a sign proof.

## Arithmetic and constants

As in the L=11/20 packet: directed Decimal endpoints with ROUND_FLOOR /
ROUND_CEILING; correctly rounded ln, exp and sqrt widened by one
adjacent representable number; Machin's identity for pi with
alternating-series bounds; Euler's constant from the digamma expansion
at 1024 with 32 Bernoulli corrections and the next correction as
remainder ([NIST DLMF 5.11(ii)](https://dlmf.nist.gov/5.11));
arctan of e^(-L) by its alternating series with a first-omitted-term
bound; all polynomial coefficients and the final sign check in exact
rational arithmetic.

The rejection controls cover an indefinite positive-diagonal matrix,
omitted residual debt, a zero pivot, a negative direction, inflated
uncertainty, an interval inverse containing zero, binary float input,
and a check that the residual debt actually charged is nonzero. A
further line records whether the uniform allowance 1/d_N alone would
also have passed; it is recorded, not asserted, in either direction.

## Two engines: directed Decimal and Arb ball arithmetic

`verify.py` is the standard-library reference engine (directed Decimal
endpoints, hand-proved enclosures of pi, Euler's constant and the
arctangent series). `verify_arb.py` is a second implementation of the
same certificate on FLINT/Arb ball arithmetic through `python-flint`
(0.9.0): every quantity is a midpoint-radius ball with outward rounding
done by the library, polynomial products run in Arb's compiled
polynomial arithmetic, and pi, Euler's constant, logarithms, exp, sqrt
and arctan are Arb's certified enclosures. The exact lower and upper
bounds handed to the exact rational LDL^T are exact dyadic numbers
below and above the whole ball. The analytic pipeline, the complement
rule, the level-by-level residual charge, the rational rounding
allowance, the shift, the all-function transfer and all controls are
the same code, line for line where the arithmetic allows.

Both engines were run on every packet directory and both results are
committed (`result-<digits>.json` from verify.py, `result-arb-<digits>.json`
from verify_arb.py; `rayleigh-*.json` likewise). For all six
certificates the exact rational witness matrices and the exact shifted
pivots are byte-for-byte identical between the two engines; the ball
enclosures are strictly inside the Decimal enclosures; the four
Rayleigh-quotient enclosures overlap.

| Packet | L | N, M, K | verify.py | verify_arb.py | witness identical |
|---|---|---|---|---|---|
| l055_control | 11/20 | 5, 12, 48 | 17 s | 0.4 s | yes |
| l0625_route | 5/8 | 5, 16, 120 | ~4 min | 1.3 s | yes |
| t53_route | 3/4 | 8, 24, 200 | ~12 min | 9 s | yes |
| boundary_route | 3/4 | 14, 24, 200 | ~22 min | 16 s | yes |
| l08_route | 4/5 | 16, 30, 300 | ~82 min | 50 s | yes |
| l085_route | 17/20 | 24, 36, 400 | ~7.5 h | 169 s | yes |
| l09_route | 9/10 | 30, 52, 700 | not run (days) | ~16 min | — |
| l095_route | 19/20 | 40, 64, 300 (T72 bound) | not run (T72 bound is Arb-only) | ~10 min | — |
| l1_route | 1 | 56, 80, 400 (T72 bound) | not run (asserts L<1) | ~23 min | — |

Reproduce with the Arb engine from the repository root, for example

~~~sh
pip install python-flint
python3 proof_lab/rh_l075/verify_arb.py 230 --packet proof_lab/rh_l075/l085_route
~~~

The L=9/10 certificate (`l09_route`, BRT-T71) was produced with the Arb
engine; the Decimal engine accepts the same packet and would need of
the order of a few days on one core for it. The L=19/20 and L=1
certificates (`l095_route`, BRT-T73; `l1_route`, BRT-T74) and their
near-null enclosures use the `kernel_tail` option (radius 5/2, sup
16/5; BRT-T72), which the frozen Decimal engine does not implement,
so they are Arb-only; the option changes one scalar (epsilon_j) and
nothing else in the pipeline.

## Independent cross-check of the analytic pipeline

`crosscheck_quadrature.py` (SciPy, binary floats; not part of the
certificate) evaluates the source formula ASM (2) directly by adaptive
quadrature on the unnormalized Legendre profiles P_0, ..., P_10 and
compares with the verifier's own operator entries
(`crosscheck_reference.json`, first 6-by-6 blocks of the exploration
runs). It shares no code path with the verifier: no Taylor kernel, no
Legendre energy identity, no boundary decomposition, no piecewise
formula. Agreement at the five widths 11/20, 5/8, 3/4, 4/5, 17/20 is
at machine precision:

| L | active prime powers | max entry deviation |
|---|---|---|
| 11/20 | 2, 3 | 4.8e-16 |
| 5/8 | 2, 3 | 2.9e-16 |
| 3/4 | 2, 3, 4 | 7.9e-16 |
| 4/5 | 2, 3, 4 | 8.4e-16 |
| 17/20 | 2, 3, 4, 5 | 1.1e-15 |

This checks the analytic source identities and the prime, pole, scalar
and boundary terms end to end on the low modes; the directed-interval
enclosures and the exact rational sign check remain the certificate.

## Scope

Fresh deterministic certificates for the fixed actual source at L=3/4,
4/5, 17/20, 9/10, 19/20 and 1. No earlier certificate is replayed. The
result is positivity on the bounded supports (0,3/4], (0,4/5],
(0,17/20], (0,9/10], (0,19/20] and (0,1]; widths beyond 1 (where the
prime power 8 activates at L=1.040) and all-width positivity/RH remain
separate open tasks. No percentage completion of an
RH proof is inferred.
