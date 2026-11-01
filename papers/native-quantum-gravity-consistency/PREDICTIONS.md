# QGC-1 — Predictions, calibration and rejection protocols

14 October 2026. All numerical values below are exact model predictions
or chosen protocol settings. No new experimental dataset has been acquired
or analysed in this packet. Certificate fixtures are not observations.

## 1. Quantum operation: predict the whole outcome record

Prepare the normalized cut `P=(I+(3/5)Z+(4/5)X)/2`. Calibrate the fixed
filter independently to `a=3/5`, `b=4/5`, with the outcome instrument
defined in THEOREM.md, equations (6)–(7). Retain both outcomes and the
number of input preparations. The predictions are:

| Reading | Exact prediction |
|---|---:|
| Success fraction | `109/125` |
| Failure fraction | `16/125` |
| Successful branch, conditional `Z` mean | `91/109` |
| Successful branch, conditional `X` mean | `60/109` |
| Failure branch, conditional `Z` mean | `-1` |
| Unconditional `Z` mean | `3/5` |
| Unconditional `X` mean | `12/25` |
| All-success fraction after `m` repetitions | `4/5+(1/5)(3/5)^(2m)` |

Estimate the branch means on separately prepared repetitions when the
readouts are incompatible. Do not treat one event as full state tomography.
Changing source amount from `N=1` to a known `N` multiplies all unnormalized
readings by `N`; it leaves the displayed probabilities unchanged.

The source-moment relation to the **linear** geometry readout is
`y1_out=y1_in`, `y2_out=(3/5)y2_in`, provided the same `alpha/gamma`
and readout adapter are used. A normalized successful branch alone
does not predict the total source. Testing this instrument on a quantum
processor would test its quantum implementation, not demonstrate that
physical gravity uses this channel.

**Rejection:** after accounting for independently assessed preparation,
measurement and filter errors, a discrepancy rejects the specified
instrument/adapter. It does not reject the full native operator algebra.
Calibration should not use the same output statistics later presented
as validation. Store all outcomes, attempted preparations and losses;
an unrecorded failure category can imitate concentration.

## 2. Joint response readout: a bound and a saturating design

At fixed nonzero `chi=alpha N/gamma`, the declared first-response
observables are `A1=chi Z/10` and `A2=chi X/47`. For a commuting joint
pointer that is unbiased for **every** source state, the added noise
operators obey

\[
 \|\mathsf N_1\|\|\mathsf N_2\|\ge\chi^4/220900.
\]

The four-effect construction in THEOREM.md predicts the probabilities

\[
 p_{st}=\tfrac14(1+s\tfrac35 r_z+t\tfrac45 r_x),\quad s,t=\pm1,
\]

with outcome values `b1=s chi/6`, `b2=5t chi/188`. On the source in
section 1 its ordered probabilities are `(1/2,9/50,8/25,0)`.
Its added noises are constant operators
`N1=(4 chi^2/225)I` and `N2=(9 chi^2/35344)I`, saturating the bound.

This is an operator-norm bound under universal unbiasedness. A variance
product measured on a single source is not its experimental estimator.
One must reconstruct the effects, first moments and second moments on
a source set spanning the Hermitian algebra, with calibration uncertainty.
The `I,Z,X,Y` record from CRC-1 specifies such a set of directions.
Subtract intrinsic response variance from pointer variance in accordance
with the definitions of `N_i`; do not count the same noise twice.

**Rejection:** failure of the exact saturating statistics rejects that
chosen readout design. An apparent violation of the universal bound first
requires auditing unbiasedness, positivity, commuting pointer assumptions
and uncertainty propagation. There is no claim of a new universal
gravitational noise floor independent of the specified response adapter.

## 3. Nonlinear geometry: two calibrations, one held-out slope

Use the exact CRC-1 triangulation, fixed source direction `(3/5,4/5)`,
and source law `F_kappa`, with finite `kappa != 1`. Control a small
source parameter with `epsilon(0)=0`. Its unknown initial scale is
`c=epsilon'(0) != 0`. Record cell volumes relative to their zero-source
values; the common volume unit cancels in logarithms.

Estimate the derivatives at zero of the log volumes for the outer, cap
and inner classes. These are `v_o,v_cap,v_in`. The protocol is:

1. Fix the cell registration, source preparation, control range and
   derivative-estimation procedure before examining the held-out class.
2. Use only the outer and cap slopes to estimate

   \[
   c=\frac{10}{3}\left(v_o+\frac{456869}{17643476}v_{\rm cap}\right),
   \qquad \kappa=1-\frac{4410869c}{1060320v_{\rm cap}}.
   \]

3. Predict the inner slope without another fitted parameter:

   \[
   v_{\rm in}^{\rm pred}=\frac{640}{141}v_o
                  -\frac{1860931}{621932529}v_{\rm cap}.
   \]

4. Also check equality of the four slopes within each symmetric class.
   These checks require the declared source and geometric symmetry.

Examples at `c=1` are:

| Coupling | Outer slope | Cap slope | Inner slope |
|---|---:|---:|---:|
| `kappa=0` | `163103/848256` | `4410869/1060320` | `912209/1060320` |
| `kappa=2` | `1729253/4241280` | `-4410869/1060320` | `1975471/1060320` |

These reproduce the predecessor results. They are not two competing
experimental datasets. The coefficient `kappa` remains a parameter until
it is independently fixed by an action principle or identified from data.

Define the held-out residual

\[
 D=v_{\rm in}-\frac{640}{141}v_o
                   +\frac{1860931}{621932529}v_{\rm cap}.
\]

For a slope covariance matrix `Sigma` in order `(outer,cap,inner)`, its
statistical variance is `ell^T Sigma ell`, where
`ell=(-640/141,1860931/621932529,1)`. Account separately for bias from
finite control range, volume reconstruction, source drift and departure
from the assumed action. Specify the rejection threshold and treatment
of multiple tests before unblinding the inner data. Agreement of two
calibrating slopes alone is not evidence for the family.

**Domain limits:** this is an exact derivative prediction at zero source.
No certified finite neighbourhood or derivative-remainder bound is
provided. Measurements at one nonzero control value cannot be directly
substituted for these slopes. Zero fitted scale or zero cap slope blocks
finite-parameter identification as detailed in THEOREM.md. Near zero
signal, ratio estimates can be poorly conditioned even when the symbolic
inverse exists. Approaching `kappa=1` is not licensed by this local branch.

An independent geometric simulation can test the declared discretization.
A physical test additionally needs an independently specified relation
between these labelled cell volumes and experimentally observed geometry.
Passing an analog simulation is not validation of quantum gravity.

## 4. What is required for a physical gravitational prediction

The current exact packets do not fix the conversion of the native source
amount or cut components to a measured stress-energy source, the unit of
squared edge length, the physical action normalization, or the relation
between the operation parameter and relational/laboratory clocks. They
also do not derive a joint quantum state of matter and spacetime geometry.

Consequently this packet does not supply a numerical gravitational phase,
force correction, spectrum, collapse time or entanglement rate. Inserting
the measured `G`, `hbar` or a chosen time scale into an unconstrained
dimensionless coefficient and fitting the output would not derive one.

For orientation, Bose et al.,
[A Spin Entanglement Witness for Quantum Gravity](https://arxiv.org/abs/1707.06050),
and Marletto and Vedral,
[Gravitationally-induced entanglement between two massive particles is sufficient evidence of quantum effects in gravity](https://arxiv.org/abs/1707.06036),
give concrete gravity-mediated entanglement proposals under explicit
assumptions. They are references for the type of physical observable
needed, not evidence that this native model predicts their outcomes or
that those proposed experiments have been performed.

The decisive next physical calculation would be a joint matter–geometry
transition amplitude or relational correlator, derived from one fixed
quantum dynamics, with its physical calibration and a prediction that
cannot be removed by refitting the source law. None of the finite
instrument or volume identities substitutes for that calculation.
