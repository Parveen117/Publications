# RH positivity at L=11/20: reproducible fixed-source certificate

Monty · 23 October 2026

This packet supports [BRT-T60--T63](../../RH_BOUNDARY_REGULAR_TRIAL.md#brt-t60-an-exact-source-column-with-a-bounded-polynomial-remainder).
For the real even closed source family declared there, the result is

\[
0<L\le11/20 \quad\Longrightarrow\quad
q_L[g]\ge10^{-11}\|g\|^2\qquad(g\in D_L).
\]

It is a bounded-support, computer-assisted working proof. All-width
positivity and RH remain open. It depends on the inherited source,
operator-domain and positive-complement theorems as well as the
analytic enclosure argument in BRT-T60--T63. This packet is not an
independent external review or a formal proof-assistant development.

## Reproduce

The evidence was generated with Python 3.12.14 using only its standard
library. Run from any working directory:

~~~sh
python3.12 proof_lab/rh_l055/verify.py 220
~~~

Run without -O or PYTHONOPTIMIZE: assertions are part of the verifier.
The command reads the adjacent trial.json and writes the adjacent
result-220.json. Copy the packet to a temporary directory first if you
want to preserve the committed evidence byte for byte.

A successful run exits with status zero, prints PASS_FIXED_WINDOW_LOWER_MATRIX,
reports both the interval shifted sign test and exact rational shifted
sign test as true, and reports all seven rejection controls as true.
A failed assertion is failure or inconclusive verification, not evidence
that the actual mathematical source is negative.

The default precision is 180; the committed current evidence uses the
explicit 220 argument above. Changing precision or code is a new
verification run. No external package, network access, or downloaded
zero table is required.

## Files and fixed choices

| File | Role |
|---|---|
| verify.py | Directed interval integration, exact rational matrix sign check, and full-family bound checks |
| trial.json | Fixed rational 12-by-5 even Legendre trial; denominator 10^14 |
| result-220.json | Complete interval matrix, scalar allowances, rational witness, exact pivots, controls, and hashes |
| README.md | Reproduction instructions and dependency map |

The code SHA-256 is
9505fdab80e29313cb9f4d2466bdeafb51c84e5ce9b6d5f998f178eac0633b66.

The trial SHA-256 is
d36cceac9070318a92868b6d30e45eeb4b5222013e0db140aa26be7fbe7bb2e0.

The 12 trial modes have degrees 0,2,...,22. The first five retained
coefficients form the identity, so P_5 Z=(P_0,P_2,...,P_8).
The remaining seven rows are a fixed approximation to the minimizing
tail. They were selected by preliminary numerical candidate exploration,
then rounded to rationals. That exploration is not a proof dependency:
the verifier accepts the fixed trial and checks all residual debt.
A different trial can fail the sufficient certificate without
disproving positivity.

Coordinates use P_(2j)(x/L) without normalization. The certified
returned-matrix margin 10^(-10) is in these coefficient coordinates,
whose basis squared norms are 2L/(4j+1). Its sign transfers by
congruence to the normalized matrix; the numerical margin must not be
silently relabeled as a normalized-basis eigenvalue.

## What is enclosed

At L=11/20 only prime powers 2 and 3 are active, with their actual
weights log(2)/sqrt(2) and log(3)/sqrt(3). The operator retains its
scalar constant, singular Legendre energy, logarithmic boundary
multiplier, pole rank-one term and both translated prime terms.

Only the regular kernel j is replaced by a degree-48 polynomial and
cosh(x/2) by a degree-16 polynomial. BRT-T60 bounds their combined
operator error by epsilon < 1.159 * 10^(-12). Logarithmic endpoint
integrals are evaluated analytically, including the square of the
boundary multiplier and all mixed prime terms. Integration uses
polynomial and log primitives, not spatial quadrature.

T57 supplies C_5 >= (1871/6300)I on the entire infinite complement.
The residual Gram is computed as

\[
\overline R^*\overline R
=(\overline T Z)^*(\overline T Z)
-(P_5\overline T Z)^*(P_5\overline T Z).
\]

Thus every unretained residual mode is included. It is not the Gram of
a finite list of sampled higher coefficients. With kappa=6300/1871,

\[
\widetilde S_5
\succeq \overline A-\kappa\overline R^*\overline R-\Delta I,
\]

where Delta < 2.485 * 10^(-12) accounts for the analytic operator
approximation, including its effect on the residual Gram.

The directed entry enclosure gives a rational symmetric matrix Q
with ||M-Q|| < 10^(-17). Exact Fraction LDL^T proves
Q-(10^(-17)+10^(-10))I positive. A separate interval LDL^T test also
passes. The evidence includes Q's integer numerators and the exact
rational pivots; no floating eigenvalue is accepted as a sign proof.

The final norm estimates verify ||Z-C_5^(-1)R|| < 8/5. The exact
completed-square identity gives q[g] >= 10^(-11)||g||^2 at the
endpoint. Source inclusion and form closure extend this to every
smaller L, as written in T63.

## Arithmetic and constants

Basic interval arithmetic uses decimal endpoints rounded downward
and upward. Logarithms and square roots are correctly rounded before
expansion by one adjacent representable endpoint. Python's decimal
behavior is specified in the
[Python 3.12 documentation](https://docs.python.org/3.12/library/decimal.html).

Pi uses Machin's arctangent identity with alternating-series bounds.
Euler's constant uses the harmonic/digamma expansion at 1024 with
32 exact Bernoulli corrections and the next correction as a rigorous
absolute remainder. The positive-real first-neglected-term bound is
given in [NIST DLMF 5.11(ii)](https://dlmf.nist.gov/5.11).
All polynomial coefficients and the final second sign check use exact
rational arithmetic.

The seven controls reject an indefinite positive-diagonal matrix,
an example with omitted residual debt, a zero pivot, a negative
direction, an inflated uncertainty allowance, an interval inverse
containing zero, and binary float input. These are targeted checks of
the acceptance logic; they are not a substitute for auditing the
source identities or the interval implementation.

## Scope

This is a fresh deterministic certificate for the fixed actual source.
It neither replays an older certificate nor performs a physical or
Monte Carlo simulation. Earlier written-only sections are historical.
The actual returned sign at L=3/4 and an all-width argument remain
separate open tasks. No percentage completion of an RH proof is inferred.
