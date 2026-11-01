# DCD-1 — Density and coherent cuts generate TVSP transport

QRT's quantum response transport has the exact source generator

    B=[Q',Q]=[Pi+Q/2, dot(theta).X].

Pi projects onto the filtered base state, Q onto its two density-response
directions, and X are the existing LMD density observables. The same
unitary transports both the base state and the TVSP response frame.
The cut weights 1 and 1/2 are forced within this two-cut ansatz.

Smooth preparations supported at bounded total density retain
1/8<kappa<3/8 and give preservation of the existing material-current and
Cameron--Martin momentum domains. Alternating original density phases
with declared coherent cut phases implements the transport in operator
norm, with explicit finite-cycle and moving-path error bounds.

The old configuration-normalizing density/current flows cannot reproduce the
small-loop base-state and response action. This is proved using the
nonzero cubic response moment, not inferred from a finite fixture.

The material-metric comparison also has a precise result. A conjugate
phase changes momentum while leaving all density-response geometry
unchanged; the original material metric depends on that missing data.
Thus its reconstruction requires a conjugate-sensitive extension.

- [Six theorems and proofs](THEOREM.md)
- [Control meaning, API and replay](USAGE.md)
- [Sources and evidence boundary](SOURCE_MAP.md)

```bash
python3.12 papers/native-density-cut-dynamics/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

Validation: **143 exact checks, 13 symbolic certificates and 23 negative
controls**. Certificate SHA256:
`834d9785f01c86f178ab2b041e50ecb90c5957c75bb2b76439846324288a4a6a`.

The coherent phase controls are an explicit operational assumption.
This packet does not derive a laboratory actuator, a physical time law,
or an identification of response curvature with spacetime curvature.
