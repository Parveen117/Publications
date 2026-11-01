# AGD-1: Geometry dynamics from the native cut action

The existing CID-D2 action now determines the geometry generator, under an explicit periodic canonical quantization contract. It replaces GCP-1's arbitrary finite mixing matrix with an infinite circle operator:

\[
H_g=\eta^2 L^2+8I-(U^2+U^{\dagger2}),\qquad
V(x)=\left(1+\frac{af(x)}2\right)I-\frac{af(x)}4(U+U^\dagger).
\]

The same cut coordinate supplies the lapse observable. With the inherited matter prescription, the joint law is

\[
\iota\eta\partial_s\Psi=
\left[H_g+\iota\eta C(v\partial_x+v_x/2)+\iota\eta gRv\right]\Psi.
\]

The written proofs establish the action reduction, essential self-adjointness of this full joint operator, finite spatial propagation, the matter source for geometry momentum, and a controlled finite Fourier approximation. No finite-dimensional matrix can preserve both the exact unitary circle shift and its canonical momentum commutator.

Concrete exact results in the specified homogeneous controls:

| Reading | Result |
|---|---|
| Initial geometry-momentum rate | `-3/10` |
| Geometry and matter energy rates | `-3/10`, `+3/10` |
| Lapse variance lost by a mean replacement | `9/200` |
| Squared leakage out of an initially even geometry sector | `9/50` |
| Finite-section state-norm error, with the inputs in AGD-T5 | At most `148035889/120000000000000 < 1/800000` |

**Executed certificate:** 92 exact checks, 8 native symbolic certificates, 8 serialized replays and 8 negative controls. The analytic domain and approximation proofs are written proofs; the native certificate verifies the specified algebra and rational controls. These counts do not certify physical quantum gravity.

The action, periodic quantum representation, scale `eta`, lapse profile and coframe dictionary remain declared choices. The global clock constraint is solved. Local gravitational constraints remain open: the proposed even-mode restriction fails under matter coupling, and a naive local commutator leaves the explicit residual `g R N f'`. The next target is a local dynamical coframe/source action that accounts for this residual and preserves its constraints.

- [THEOREM.md](THEOREM.md): derivations, domain proofs, error estimate and constraint obstructions.
- [CERTIFICATE.json](CERTIFICATE.json): native witnesses, exact fixtures and rejected claims.
- [SOURCE_PINS.json](SOURCE_PINS.json): unchanged predecessor and engine provenance.
- [Canonical research record](../../QUANTUM_GRAVITY_RESEARCH.md): lineage and next obligation.

From the repository root, using Python 3.12 and Node.js:

```sh
python3 papers/native-action-geometry-dynamics/run.py --check
python3 scripts/replay_quantum_geometry.py
```

Python handles I/O and source verification only. Scientific arithmetic and word replay use the unchanged native engines. The check command recomputes the certificate without writing and compares its exact bytes and checksum.
