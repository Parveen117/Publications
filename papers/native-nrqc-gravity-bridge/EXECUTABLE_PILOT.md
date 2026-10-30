# Executable source-covariance pilot and comparison

Author: Monty. 30 October 2026. Parent analytic bridge at Publications
`4de11ba53035125c3910d9505d711126a0705df3`; unchanged NRQC at
`fdd39048a6ddaa2dd78aa3197c274d7766f31a0e`.

## Closed form for the source coefficient

The earlier analytic bridge left B(h) as a double integral. Translation
invariance and the cosine convolution reduce it to

    B(h)=(exp(h)/2) integral_0^1 cos(2 pi r)
                                   (2 sin(pi r))^(-h) dr.

Write I0 for the same integral without the cosine and I1 for the integral
with it. Integrating the derivative of sin(theta)(2 sin(theta/2))^(-h)
over a period gives (1-h/2)I1=(h/2)I0. The boundary term vanishes like
theta^(1-h), since 0<h<1. Euler's beta integral and the gamma duplication
identity give I0=Gamma(1-h)/Gamma(1-h/2)^2. Therefore

    B(h)= exp(h) h / [2(2-h)]
                    * Gamma(1-h)/Gamma(1-h/2)^2.         (1)

The beta/gamma identities are standard, not new results:
[NIST DLMF 5.12](https://dlmf.nist.gov/5.12) and
[5.5](https://dlmf.nist.gov/5.5). The source-specific result is their
application with the normalized cosine and the full Wick zero-mode factor.

## Rigorous arithmetic and finite encoding

All new source arithmetic uses the frozen native JS rational engine. From
the standard log-gamma series ([DLMF 5.7](https://dlmf.nist.gov/5.7)),

    log I0 = sum_(n>=2) zeta(n) h^n (1-2^(1-n))/n.

At h=1/4,1/2,3/4, the implementation takes 80 outer terms. For each zeta it
sums k^(-n), k=1..1024, with the integral-test lower and upper tails
1/[(n-1)1025^(n-1)] and 1/[(n-1)1024^(n-1)]. Every summand is rounded
outward on a 2^(-80) grid using exact integers. The remaining outer tail is
bounded by 2 h^81/[81(1-h)]. Positive exponential Taylor sums through 120
are bounded by the next term divided by 1-x/122, with x<4 checked.
These are rational enclosures, not ordinary floating-point gamma outputs.

For each source point, choose 22-bit dyadic alpha, beta and gamma in the
three-level construction from ANALYTIC_BRIDGE.md. The actual encoded matrix
has diagonal alpha^2 and beta^2+gamma^2, and cross term i alpha beta.
Compare each of those exact rational entries against the source enclosure.
The predeclared absolute tolerance is 10^(-5) per entry. No error on full
source dynamics or on an uncomputed observable is inferred.

The two weighted controls expand into six signed Pauli words. The existing
NRQC geometry query computes their exact covariance matrix G. Native JS then
forms W G W^dagger with the declared dyadic coefficient matrix W. No NRQC
code is changed, and no unsupported weighted-geometry API is invented.

## Independent comparison and controls

The baseline uses explicit 4x4 matrices from their three-level definitions
and computes Tr(P H_a H_b)-Tr(P H_a)Tr(P H_b). It does not use the NRQC
Pauli product implementation. Four finite states are tested: vacuum,
flipped first qubit, Bell state and a phase-sensitive complex product state.
Both existing NRQC backends, stabilizer and observer, must agree exactly
with the independent matrices for all three coefficient choices.

Only the finite vacuum comparison is mapped to the LMD vacuum moments.
The other three states validate the finite algebra adapter; they are not
asserted to represent three new continuum LMD states.

Negative controls remove the second tangent, use the leading one-qubit
tensor, corrupt an ordered XY phase, reverse the cross-response sign, and
omit a projector factor. The executable checks reject these claims. This
is not a replay of the frozen continuum certificates.

The two-qubit model is supplied with coefficients calculated from (1).
Consequently this pilot verifies an encoding; it does not show that NRQC
computes the continuum coefficient more efficiently. Its output agrees
exactly with the finite baseline and approximately with the source moments,
with the approximation error separately bounded.

## Run and evidence

Use Python 3.12 for orchestration and the user-authorized existing NRQC engine;
new source arithmetic, weighted assembly and baseline calculations are native JS.
From the Publications root, with a clean checkout of NRQC at the pinned commit:

```bash
python3 papers/native-nrqc-gravity-bridge/run_pilot.py \
  --nrqc /path/to/pinned/quantum-classical-provisional \
  --out /tmp/nrqc-gravity-pilot
```

The committed `pilot_results/` contains requests, actual outputs, rational
source intervals, execution metadata, verification and hashes. No hardware,
paid compute or GitHub workflow dispatch is involved. Memory consumption was
not compared, so there is no memory-advantage claim.

`verification.json` is the authoritative measured summary. Timing uses five
alternating cold-process repetitions. The native command computes all three
weighted matrices; the NRQC command computes one common six-Pauli matrix
which is reusable for all three points. Weighted assembly is additional for
NRQC. Both share the same separately reported source preparation. CLI startup,
language runtime and JSON overhead are included, so this is an end-to-end
local check, not an isolated algorithmic benchmark or a scaling conclusion.

The separate BOUNDED_COFRAME.md constructs a restricted source-geometric
reading. That observable has NOT been simulated by this covariance pilot.
The unrestricted geometric inverse, full field dynamics and 3+1 gravity
remain open. No computational advantage is established.
