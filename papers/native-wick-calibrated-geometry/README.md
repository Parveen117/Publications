# WCG-1 — Wick-calibrated joint readings give resolved material geometry

PCR's phase/current information now has an actual nonlinear material
measurement in the same LMD Hilbert sector. A countable coherent detector
produces a positive measured Wick mass N_u and conjugate labels v. Its
calibration is exact:

    E_probe[N_(s+probe)(w) | s]=M_s(w),
    current first moment=D(k),
    current second-moment domain=Dom D(k).

The measured mass has doubled Gaussian covariance. Its half-variance
subtraction therefore differs from the system's. In the elementary
range 0<h<1/2, conditional L2 convergence proves the calibration for
the full positive density form, not only a finite-mode approximation.

At every fixed positive spatial resolution ell, smoothing the measured
mass and current gives smooth original material profiles. Their normalized
POVM supplies bounded metric operators on every state and, with a recorded
density threshold, a positive self-adjoint circumference mean with a dense
domain. The global normalization in the original inverse is retained.

PCR's actual density-phase-pulsed coherent states recover smooth classical
geometry in H3/C2 probability under the sufficient joint limit

    ell -> 0, h ell^(-8) -> 0.

The density and current squared-error rates are O(h ell^-7) and
O(h ell^-8). The conjugate law is generally non-Gaussian, so this does
not imply convergence of unbounded length expectations. At fixed h,
raw point r is proved non-tight as ell vanishes, and ungated length
has infinite vacuum mean even at fixed ell.

- [Five theorems and proofs](THEOREM.md)
- [Reading contract, API and replay](USAGE.md)
- [Sources and evidence boundary](SOURCE_MAP.md)

```bash
python3.12 papers/native-wick-calibrated-geometry/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

Validation: **354 exact checks, 15 symbolic certificates and 50 negative
controls**. Certificate SHA256:
`5c42b770cd54bd083ace48863399d48c8a5e319a3976a92093c82cc2bd93f6b7`.

The detector's exact first moment does not establish full material-flow
covariance of its outcome law. The apparatus, spatial resolution and
acceptance threshold remain specified measurement data. This constructs
resolved material-slice observables, not a pointwise quantum spacetime.
