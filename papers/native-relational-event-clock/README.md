# REC-1 — Dynamics from relational event readings

LNC's homogeneous matter-free sector admits phi as a local event-clock when
L>0 and pi is nonzero. Its reduced Hamiltonian is `H_red=L V(phi)/pi`.
The equations derive the conserved comparison `L/pi`, together with
`pi^2+2 W(phi)` where `W'=V`. Complete observables at a clock reading form a
canonical pair within this sector.

Normalized length and momentum readings lie on the same diagonal. This does not
select the potential or rate: two explicit potentials produce identical endpoint
readings and the same diagonal, but different interior readings. Dense
differentiable readings recover `V=-(pi^2)'/2`; two readings recover only its integral.

The packet reuses IRT's 45-degree calibration result and LNC's selected action.
It includes clock-failure, global-winding and noncommutative-closure controls.
This is a classical homogeneous reduction, not a completed local quantum theory.

Replay: **115 exact checks, 4 native symbolic certificates/replays and 13 negative controls**.

- [Derivation, domain, examples and scope](THEOREM.md)
- [Native exact certificate](CERTIFICATE.json)
- [Source pins](SOURCE_PINS.json)
- Replay: `python3 papers/native-relational-event-clock/run.py --check`

The continuum proof is written explicitly. Native checks and algebraic witness
replays certify their stated finite and symbolic identities, not all analytic claims.
