# PCR-1 — Phase and current complete the density-response coordinates

The existing density pulse supplies the phase information missing from
DCD's density-only response geometry. In the fixed prepared family

    psi_(theta,beta)=exp(-i beta.X/h) psi_theta,

the original Cameron--Martin currents P_i=D(k_i) have exact means

    m_i=<X_i>, p_i=<P_i>=beta_i m_i, beta_i=p_i/m_i.

The four ensemble means identify this entire declared family: m fixes
theta and p/m fixes beta. Density statistics alone remain phase-blind.
No extra Hilbert-space factor or finite matrix current is introduced.

The response tensor is Kähler. Its real metric has rank four while its
complex quantum geometric tensor has rank two. With A=i<psi,dpsi>,

    F_Berry=(1/h)sum_i d(log m_i) wedge dp_i=-2 omega.

The exact source lift retains TVSP response transport and has controlled
original current domains. Mixed amplitude/phase loops have a nonzero
Berry phase; the uncorrected source lift's dynamical phase cancels it.
A separate phase-pulsed coherent family gives smeared semiclassical
material data and an explicit homogeneous classical-label reconstruction.

- [Five theorems and proofs](THEOREM.md)
- [Meaning, API and replay](USAGE.md)
- [Sources and evidence boundary](SOURCE_MAP.md)

```bash
python3.12 papers/native-phase-current-response/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

Validation: **247 exact checks, 14 symbolic certificates and 35 negative
controls**. Certificate SHA256:
`7e57cb1ff8971775d559de29e099887262cbffb8b4adafe02fe1803bb855a737`.

The four coordinates describe a quantum control family. They do not
define a Lorentzian spacetime, a sharp pointwise quantum material metric,
or a completion of the full 3+1 quantum-gravity source.
