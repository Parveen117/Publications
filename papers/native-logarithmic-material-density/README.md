# LMD-1 — Logarithmic material density and local flow

An explicit alternative reduced representation after MRM-1:
logarithmic Gaussian covariance supports every smooth material-circle
flow, and a Wick exponential supplies a positive smeared density in the
L2 range 0<hbar<1 (dimensionless units).

The exact density transformation requires the quantum shift
a_h=1+hbar/2 in s -> s composed with g + a_h log g'.
Keeping the classical coefficient 1 gives a different density weight.
The old smooth MPQ representation is preserved, not silently replaced.

Read [THEOREM.md](THEOREM.md), [SOURCE_MAP.md](SOURCE_MAP.md), and
[FABLE_REVIEW.md](FABLE_REVIEW.md). The latter audits Fable's bridge at
97bb7fe and gives the precise follow-up corrections.

    python3.12 papers/native-logarithmic-material-density/run.py --check

Scientific finite arithmetic uses the unchanged native F/Cut and
Paninian engines. Python 3.12 performs I/O only. Analytic proofs,
external theorem inputs and finite checks have separate scopes.

This constructs local material-flow generators and a renormalized
density, not the original pointwise metric or a 3+1 quantum-gravity
theory. An exact three-dimensional shear test proves that this scalar
Gaussian prescription does not automatically extend to full spatial
diffeomorphism covariance.
