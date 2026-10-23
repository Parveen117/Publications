# RH positivity at L=3/4 and L=4/5: reproducible fixed-source certificates

Monty · 23 October 2026

This packet supports [BRT-T64--T69](../../RH_BOUNDARY_REGULAR_TRIAL.md#brt-t64-a-width-dependent-complement-comparison-with-boundary-energy-allocation).
For the real even closed source family declared there, with the three
active prime powers 2, 3, 4 carrying their actual weights
log(2)/sqrt(2), log(3)/sqrt(3), log(2)/2, the results are

\[
0<L\le3/4 \quad\Longrightarrow\quad
q_L[g]\ge\mu_{3/4}\,\|g\|^2,\qquad
0<L\le4/5 \quad\Longrightarrow\quad
q_L[g]\ge\mu_{4/5}\,\|g\|^2\qquad(g\in D_L),
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
| `l0625_route/` | `--packet l0625_route` | L=5/8, N=5, `boundary_split` rule; two-prime control for the same rule |
| `l055_control/` | `--packet l055_control` | L=11/20, N=5, the `rh_l055` trial and declared T57 constant re-run through this verifier; reproduces the exact rational witness of that packet |

Each directory holds `config.json` (width, active prime powers, retained
modes, kernel and cosh degrees, complement rule, margins), `trial.json`
(fixed rational columns) and the committed `result-<precision>.json`
(or `rayleigh-<precision>.json`).

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
~~~

Run without -O or PYTHONOPTIMIZE: assertions are part of the verifier.
Each command reads the packet's `config.json` and `trial.json` and writes
the result file next to them. Copy the packet to a temporary directory
first if you want to preserve the committed evidence byte for byte.
The two L=3/4 certificate runs take of the order of half an hour each
on one core, the L=4/5 run of the order of two hours; the others take
minutes.

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
N=16. The remaining rows are a fixed rational rounding of the truncated
minimizing tail of the 24- or 30-mode operator matrix found by
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

The operator is BRT-T60's (236) with the three-prime subtraction
P_234 = sum_n a_n (S_(log n)+S_(-log n))/2. Only the regular kernel j and
cosh(x/2) are replaced, by a degree-200 (L=3/4) or degree-300 (L=4/5)
Taylor polynomial and twelve even Taylor terms; BRT-T66 and T69 bound
the combined operator error epsilon.
Logarithmic endpoint integrals, the square of the boundary multiplier,
and every mixed prime term are integrated analytically on the seven
pieces 0, 1-w', 1-v', u', 1-u', v', w', 1 of y=(x+L)/(2L), where
u'=log2/(2L), v'=log3/(2L), w'=log4/(2L). The prime-2 translation caps
overlap at both widths; the piecewise integration does not use cap
disjointness. At L=4/5 the verifier checks log 5 > 2L = 8/5 with
directed arithmetic, so the prime power 5 is not active.

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

## Scope

Fresh deterministic certificates for the fixed actual source at L=3/4
and L=4/5. No earlier certificate is replayed. The result is positivity
on the bounded supports (0,3/4] and (0,4/5]; the sliver up to
(log 5)/2, the activation of the prime power 5 at log 5 < 2L, and
all-width positivity/RH remain separate open tasks. No percentage completion of an
RH proof is inferred.
