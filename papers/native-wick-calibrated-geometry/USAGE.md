# WCG-1 reading contract and native replay

## What the detector reports

| Object | Definition and role | Domain or calibration |
|---|---|---|
| N_u | Measured Wick mass from noisy configuration u | Gaussian parameter 2h; elementary L2 range 0<h<1/2 |
| N_u(w) | Positive smeared density output, w>=0 | First-moment form is exactly that of system M(w) |
| v(k) | Conjugate output | First moment D(k); second-moment domain exactly Dom D(k) |
| q_ell | Positive heat convolution of N_u | Smooth, strictly positive at every ell>0 |
| d_ell | Heat-smoothed conjugate output | Smooth at every ell>0 |
| r_ell,m_ell | Original nonlinear material inverse of epsilon q_ell,d_ell | Entire profile, including global P=integral m, is reconstructed |
| g/(1+g) | Bounded reading of g=exp(2r) | Bounded self-adjoint operator on the whole Hilbert space |
| 1_A L | Accepted circumference with rejection recorded separately | Positive self-adjoint mean defined by its closed form |

This is one joint noisy measurement, not a joint sharp density/current
measurement. The conditional identity E_probe N_(s+probe)=M_s is what
makes the density reading source-calibrated. Its full probability law
still has extra detector fluctuations. The system remains in the LMD
sector; the Gaussian probe is an explicitly declared apparatus resource.

Smoothing the positive measured mass is essential. Neither exp(H_ell u)
nor a Fourier truncation of q_ell is the stated prescription. Positive
cutoffs use H_ell*[exp(u_N-h c_N)dz], including the zero-mode counterterm.

## Exact native adapter

From the repository root:

```javascript
const WCG = require('./papers/native-wick-calibrated-geometry/wick_geometry_adapter.cjs');

const calibration = WCG.calibration('1/4', 4, '3/5', 8);
const mode = WCG.modeCovariance('1/4', 3, '1/7');
const compressed = WCG.zeroModeCompression('1/4', '1/2');
const profile = WCG.zeroAngleProfile([1, 2], [0, 0], {
  sign: 1, threshold: '1/2'
});
```

`calibration(h,cutoff,mean,order)` derives finite polynomial coefficients
by actual Gaussian moments. Its conditional coefficients agree with
exp(t s-h c_N t²/2), and its vacuum coefficients with one. The finite
coefficient record is not a sampled Gaussian field, numerical exponential,
or a proof of convergence of a continuum measure.

`modeCovariance(h,n,phaseMassVariance)` returns the exact system/probe
configuration and current noise factors using lambda_0=1 and
lambda_n=1/(2n). The last argument represents an independently established
Var M(kb). The example's 1/7 is an illustrative supplied variance, not
a computed continuum value.

`zeroModeCompression(h,rootH)` computes the actual first-two-Hermite
matrix elements of the one-mode coherent measurement, provided rootH
is rational and its square equals h. It gives the calibrated density
matrix and current first/second moments. A compressed current matrix
does not represent the full current commutation algebra; compressed
second moments differ from the square of compressed first moments.

`zeroAngleProfile` uses the exact inverse value Y^(-1)(0)=0 and a=3/2
for accepted profiles with zero current. It retains the global P and
returns rational metric and bounded-metric labels. In the example these
are [4/9,16/9] and [4/13,16/25], with circumference one under normalized
two-cell quadrature. Rejected records retain their flag. This finite
comparison does not evaluate Y inverse or cosine at arbitrary arguments.

The other helpers expose exact Gaussian/Wick coefficients, calibrated
finite-mode noise, normalized profile weights, threshold factors,
isometry compression and rational error budgets. No floating-point
scientific arithmetic or new operator engine is used.

## Resolution and error budget

For the actual PCR phase-pulsed coherent states, WCG-T4 supplies constants
C_rho,C_d from integrable covariance bounds and proves

    E||q_ell-H_ell rho_*||H3² <=2085 C_rho h ell^-7,
    E||d_ell-H_ell d_*||H3² <=(8283/2+2085 C_d)h ell^-8.

`rateBudget` accepts h, ell, `densityTraceConstant`,
`currentTraceConstant` and an optional tolerance. Supplying valid trace
constants gives the two squared-error bounds and their joint Markov
probability bound. Its scalar arithmetic cannot infer those constants
from an unspecified continuum state. The bound concerns the smoothed
density/current errors about their smoothed means; total classical
error also includes deterministic heat bias and the nonlinear inverse's
local continuity. A Markov upper bound can exceed one.

`sourceRateBudget` accepts the theorem's `densityKernelConstant` C_rho
and nonnegative `currentPhaseConstant` C_d, and derives both trace
constants above. `heatTraceSplit(ell,L)` exposes the rational low/high
bounds for an integer split satisfying 1/ell<=L<=2/ell. These helpers
keep the analytical input constants distinct from their exact arithmetic.

`heatMajorant` records the rational bound exp(-x)<=m!/x^m for
x=2ell²n². It does not numerically evaluate the heat multiplier. The
explicit split into low modes and this summable high-mode tail gives
the conservative constants 2085 and 8283 in the theorem.

The sufficient classical schedule ell=h^(1/9) makes h ell^-8=ell
tend to zero. It is not an optimal resolution or a dimensional hardware
specification. At fixed h, reducing ell to zero is a different limit:
the raw point r output is proved non-tight.

## Threshold, instruments and replay

Acceptance is min q_ell>=delta. Keep the rejection probability and its
distinct outcome flag. The accepted-length mean operator is an
unnormalized first moment; a conditional accepted mean additionally
divides by the state-dependent acceptance probability. Without the gate,
the vacuum mean circumference is infinite even at fixed ell.

Every finite Hermite vector belongs to the accepted-length operator's
domain and has all finite accepted POVM moments. Such higher moments
are not automatically powers of the mean operator. These domain estimates
do not imply finite length moments for every state or convergence of
mean length in PCR's non-Gaussian phase-pulsed family.

For a finite detector partition the chosen square-root instrument has
K_b=sqrt(E(B_b)); its cutoff operations converge in trace norm on each
state. It need not be the apparatus-induced instrument of MDD, since
effects alone do not determine the postmeasurement state.

```bash
python3.12 papers/native-wick-calibrated-geometry/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

Python 3.12 handles I/O and hashes. Frozen native F/Cut and Paninian code
handles 354 exact checks, 15 symbolic certificates and 50 negative controls.
The integrated research
chain contains 48 stages. The written proofs establish the continuum
measure, domains, cutoff convergence, classical limit and obstructions.
