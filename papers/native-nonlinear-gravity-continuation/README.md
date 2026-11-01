# EGR-2 — Nonlinear cut-source geometry and local focusing

14 October 2026. Continuation of EGR-1 using its same nine four-simplices,
ten interior squared edges, boundary data and unchanged native engines.

Read [THEOREM.md](THEOREM.md) for the model contract and proofs.

The new result is a local analytic solution of **all ten nonlinear metric
equations**, with an exact second-order response calculation. Freezing
the eight flat-background displacement directions does not produce that
solution: nonlinear compatibility selects their first-order correction.

For the specified source continuation, complementary rank-one cuts give:

| Fixed readout, with positive source strength | Metric action near the branch | Initial volumes |
|---|---|---|
| `(3/5,4/5)` | Saddle: 2 positive and 8 negative Hessian directions | All nine increase |
| `(-3/5,-4/5)` | Strict local minimum: 10 positive Hessian directions | All nine decrease |

Both cuts have the same trace, rank and spectrum. This makes source
direction and the coupling law necessary parts of the concentration
statement. The stability here is variational stability of a finite
Euclidean metric action, not dynamical or quantum stability.

The stationary metric branch supplies a nonlinear cut force
`H_eff=-alpha*N*(y1*Z+y2*X)`. Exponential cut transport survives, while
the simple exponential law for the source-channel ratio receives a
nonzero nonlinear correction. The earlier ratio law remains the leading
metric approximation.

Exact polynomial cell-volume conditions identify the boundary of valid
metrics. Initial contraction is established on the stable branch; arrival
at zero volume, singularity transport and a finite collapse threshold
are not established. The local theorem supplies no certified numerical
radius in source strength.

## Reproduce

From the repository root:

```sh
python3.12 papers/native-nonlinear-gravity-continuation/run.py --check
```

Python performs I/O only. Node runs the exact rational geometry and the
Pāṇinian word certificates, using the pinned original engines extracted
from `Native_Phase_Word_Reconstruction_Source.zip`.

[CERTIFICATE.json](CERTIFICATE.json) contains the full third-derivative
tensor, both branches' finite jets, the nonsingular continuation Jacobian,
volume identities, cut feedback coefficients and negative controls.
The check recomputes the payload and compares it byte for byte. Use
`--write` only after deliberately updating the source pins.

Analytic existence, inertia for sufficiently small source, and local flow
follow the written proofs and their hypotheses. The native certificate
checks their finite algebraic inputs. Counts include symbolic proof and
replay checks; they are not counts of independent theorems. No generic
proof-assistant certification or complete quantum-gravity result is claimed.

All prior packets remain unchanged. The source coupling is an explicit
affine continuation in squared lengths, and is only invariant under the
old vertex displacements to first order at the flat background. A different
nonlinear coupling can change the stability and focusing result.
