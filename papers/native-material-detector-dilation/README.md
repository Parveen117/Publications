# MDD-1 — Material detector dilation and calibration

Frozen parent: CMG-1 at `3bec523b8b89a7c98b204da012920d0f8b2dea47`.

This packet derives an explicit additional apparatus for the existing
RMG/CMG joint measurement. A matched Gaussian probe, canonical rotation
and two commuting output quadratures realize exactly its POVM. Finite
versions have a classical first-class source lift; the reduced quantum
rotations, clock norm and finite-bin instruments have a countable limit.

It also proves why the old J-only apparatus cannot produce this joint
readout, computes the added noise, and calibrates resolution-independent
smooth adapted readings and the original J-density first moment. It does
not remove the CMG point-geometry or ungated-length obstructions.

Read [THEOREM.md](THEOREM.md) and [SOURCE_MAP.md](SOURCE_MAP.md).

Replay: `python3.12 papers/native-material-detector-dilation/run.py --check`
from the repository root. Python performs I/O only; the frozen native JS
rational/Cut and Paninian engines perform the scientific witnesses.
General analytic theorems are written proofs, not certified by finite counts.

New apparatus, probe preparation, ideal quadrature readout and pulse area
are declared inputs. This is not a local material hardware model, a new
unreduced quantum ordering, a detector-selected covariance, or physical
3+1 quantum gravity. Fable's ergodicity branch is not a dependency.
