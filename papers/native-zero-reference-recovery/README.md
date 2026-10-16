# ZCR-1: Fixed reference and controlled classical-metric recovery

The proposed fixed-zero idea now has an operational version: the geometry reading must remain definite to the required precision under evolution. Naming a reference zero does not by itself remove geometry dispersion or observation-sector leakage.

The exact scalar-sector criterion is reused from GCP-1. The new main result is an explicit error bound comparing the **full spatially varying AGD model** with the flat-metric matter equation. The geometry wavefunction evolves in the comparison; it is neither frozen nor repeatedly measured/reset.

For the normalizable circle preparation of width parameter `m`, with `d=||partial_x chi||`,

\[
\|\Psi(s)-\Psi_{\rm ref}(s)\|
\le a(d+1/2+|g|)\left[s\sqrt{q_m}+\frac{s^2}{2}\sqrt{\epsilon_m}\right],
\]

where both preparation quantities have exact formulas:

\[
q_m=\frac1{4m+2},\qquad
\epsilon_m=\eta^2\frac{m^2}{4m-1}
+\frac{2(4m+1)}{(m+1)(2m+1)}.
\]

Choosing `m=n²` and `eta=1/m` gives convergence for every fixed finite operation interval and bounded matter derivative norm. This is an explicit family of preparation/action scales, not a change of units within a fixed physical model.

| Exact result or certified budget | Value |
|---|---|
| Retained second-moment defect | `(QHP)†(QHP)` |
| Finite-carrier retained-amplitude error | At most `||QHP||² s²/(2 eta²)` |
| Packet `m=2`: lapse mean at `x=1` | `103/100` |
| Same packet: lapse variance at `x=1` | `27/20000` |
| Recovery example: joint state-norm error | At most `279/200000` |
| Same example: reading error for observable norm at most one | At most `279/100000` |

The recovery example uses `a=3/5`, `n=1000`, `m=1000000`, `eta=1/1000000`, `g=0`, `0<=s<=1`, and `d<=1`. The million-mode preparation parameter is evaluated analytically, not by claiming a million-mode simulation. The bounds are mathematical budgets, not experimental accuracy.

**Executed native certificate:** 169 exact checks, 11 symbolic certificates, 11 serialized replays and 13 negative controls. The general packet moments follow from an integration recurrence; finite Fourier fixtures independently check them. Written domain and evolution estimates establish the infinite-dimensional recovery theorem. These analytic proofs are not proof-assistant certified by the finite checks.

The limit recovers a **classical metric for matter propagation**. The matter equation remains quantum. Classical trajectories, a unique measurement-outcome rule, physical constants, local gravitational constraints and a complete gravitational classical limit remain separate obligations. Narrowing at fixed `eta` raises kinetic energy, so the stated controlled limit cannot be replaced by narrowing alone.

- [THEOREM.md](THEOREM.md): full derivations, domains, recovery bound and scope.
- [CERTIFICATE.json](CERTIFICATE.json): exact arithmetic, source-bound word proofs and controls.
- [SOURCE_PINS.json](SOURCE_PINS.json): predecessor and unchanged-engine hashes.
- [Canonical research record](../../QUANTUM_GRAVITY_RESEARCH.md): programme status and next target.

Run with Python 3.12 and Node.js from the repository root:

```sh
python3 papers/native-zero-reference-recovery/run.py --check
python3 scripts/replay_quantum_geometry.py
```

Python handles I/O and hash checks only. The unchanged native engines perform all scientific certificate arithmetic and symbolic replay.
