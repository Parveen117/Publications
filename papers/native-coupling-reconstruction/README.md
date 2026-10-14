# Native coupling reconstruction — CID-D3

The information and curvature coefficients in CID-D2 can now be separated
by a specified observation protocol. Its original symmetric helix cannot
do this: the entire dynamics depends only on `alpha+2 beta`. The new
protocol changes one calibrated control and starts from a suitable cut
at rest. Two scalar readings of its covariant acceleration recover
`beta/alpha`, even under an unknown regular nonlinear reparameterization.

[THEOREM.md](THEOREM.md) contains the proof, force signatures, observation
formula, exact protocol and independent-preparation prediction test.

## Results

| Result | Meaning |
| --- | --- |
| `T_C=2 T_I` in the original protocol | Distinct coefficient pairs yield identical cut, generator and full transport histories; five-cycle memory does not identify the action weights |
| Initial-rest acceleration | `[D_tau²P,P]=f'²(alpha T_I+beta T_C)`; the unknown parameter speed is one positive scalar |
| Two native force readings | A nondegenerate force Gram matrix recovers the coefficient ratio, with no physical clock calibration |
| Exact protocol criterion | For controls `(s R,J,L)`, `det G=256(s²-1)²(s²-2)²x²y²z²`; both blind control settings and blind cuts are explicit |
| Independent prediction | A fixed recovered ratio predicts force direction on another preparation; the residual is exactly `c²(alpha rho-beta)² det G` |

The rational example uses `s=2`, cut coordinates `(2/3,2/3,1/3)` and
`det G=16384/81`. Synthetic source coefficients `(2,1)` return ratio
`1/2`; the previously indistinguishable pair `(3,1/2)` returns `1/6`.
These are exact reconstruction fixtures, not experimental measurements
or predictions of universal constants.

## Reproduce

```bash
python3.12 papers/native-coupling-reconstruction/verify.py --check
```

Python **3.12 only** orchestrates; the pinned unchanged Node **24.19.0**
RKF engine executes all arithmetic, matrices and native word proofs.
The command replays CID-D2 and its full CID-D1/QG-R2/QG-R1/phase/six-packet
chain, checks direct source pins and reproduces the new report byte for
byte without modifying committed evidence. `--write` explicitly generates
only this continuation's report and checksum.

New evidence: **25 algebraic certificates and serialized replays, 91 exact
checks and 16 negative controls**. This includes 16 parameter-change cases,
three frame cases, and predictions on three independent preparations
without refitting the coupling ratio.

| Coverage | Evidence |
| --- | --- |
| Force extraction and cut constraints | Native noncommutative word proofs |
| Gram recovery, scale cancellation, predictive residual | Native central-polynomial proofs |
| Symmetric ambiguity and calibrated control criterion | Native sphere-cut matrix calculus, with confluence and dagger checks |
| Regular nonlinear parameter and variable-frame changes | Written covariant chain rule and exact rational cases |
| Equal complete histories for equal symmetric kappa | Equality of their full vector fields, proved by native force identity; written flow implication |
| Physical measurement and action-family selection | Open; ideal prepared cut and covariant derivatives are the observation inputs |

## Admission rules and controls

The controls and weights must be calibrated together. A mere change of
units or frame does not create another experiment. Initial rest must be
covariant rest. The Gram determinant must be positive; zero signals,
negative recovered coefficients, singular parameters and zero prediction
denominators are rejected. The prediction needs positive orientation,
not only a squared alignment check. Pure-curvature data is represented
without dividing by a zero information coefficient.

Away from rest an exact counterexample changes the apparent ratio from
`1/2` to `5/8` solely through nonlinear reparameterization. That is why
the rest condition is part of the theorem. No claim about the stability
of noisy second derivatives is made.

The predecessor's printed `Omega_ab^{,2}` exponent is corrected to
`Omega_ab^2`; its source pin and report are regenerated. Its 39 native
proofs, 150 checks, 11 controls and mathematical implementation are
unchanged.
