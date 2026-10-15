# SRT-1 — Axis information through a singular response chart

A diagonal invariant does not determine the full axis response. This packet
retains the missing blocks and proves a constructive singular-chart theorem:
for `A(t)=diag(t^w_i)`, the transformed response `C=A^-T H A^-1` extends
analytically across zero exactly when `H_ij` vanishes to order at least `w_i+w_j`.
Its boundary value is the coefficient of that order. Higher boundary derivatives
require correspondingly higher jets.

For `H=[[a,bt],[bt,ct^2]]`, the collapsed point retains only `a`, while first and
second jets recover `b,c` and the full symmetric response. Native exact controls
show why symmetry alone, diagonal readings, or unsigned ratios are insufficient.

The signed thermo/EM potential mapping and the separate differential-chain
condition for transporting electromagnetic field equations are explicit.
This is a response-chart theorem with declared scaling weights, not a proof of
physical singularity resolution or completed quantum gravity.

Replay: **143 exact checks, 5 native symbolic certificates/replays and 21 negative controls**.

- [Theorem, reconstruction and scope](THEOREM.md)
- [Native certificate](CERTIFICATE.json)
- [Source pins](SOURCE_PINS.json)
- Replay: `python3 papers/native-singular-response-transport/run.py --check`
