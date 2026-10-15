# GCP-1: Quantum geometry controls propagation

The geometry observable now enters the derivative coefficient of the matter equation. This closes the specific operator-construction gap left by CCL-1's local phase attachment.

\[
\boxed{\partial_s\Psi=\left[(C\otimes V)\partial_x+\frac12 C\otimes V'+gR\otimes V+I\otimes B\right]\Psi.}
\]

The full theorem proves:

- Unique unitary evolution and finite propagation for bounded positive lapse branches and bounded anti-dagger local coupling.
- A shared characteristic cone for matter species within each geometry branch.
- An exact curved-coframe realization, including the connection term and its nonzero derivative-commutator curvature.
- Reciprocal kinetic coupling: geometry changes matter propagation, and matter momentum changes geometry-transition dynamics.
- An obstruction to replacing quantum geometry by its mean: the missing principal coefficient is precisely `Var(V)`.
- An exact scalar-metric sector criterion and a finite quantum entanglement example.

For the supplied curved branch `v_1(x)=1+(3/5)x²/(1+x²)`, the curvature is `R_1(0)=-12/5` and the spin-transport commutator at the origin is `(3/5)C`. For the separate homogeneous quantum demo, the metric-transition probability is `256/625` and the entanglement linear entropy is `73728/390625`.

**Executed verification:** 115 exact checks, 6 native symbolic certificates, 6 serialized replays and 9 negative controls. The analytic existence/locality proofs are written proofs; the certificate covers the declared native algebra and finite exact fixtures.

The model supplies its branch profiles and geometry mixing `B`. Deriving them from a gravitational action with a closed constraint algebra remains the next target. This is not a completed physical quantum-gravity theory or an experimental result.

- [THEOREM.md](THEOREM.md): assumptions, constructions, proofs and boundaries.
- [CERTIFICATE.json](CERTIFICATE.json): exact results and serialized proof witnesses.
- [SOURCE_PINS.json](SOURCE_PINS.json): source and unchanged engine hashes.
- [Research status](../../QUANTUM_GRAVITY_RESEARCH.md): canonical continuation target.

From the repository root, with Python 3.12 and Node.js:

```sh
python3 papers/native-geometry-controlled-propagation/run.py --check
```

Python performs I/O and provenance checks only. All scientific fixture arithmetic and word-proof replay use the unchanged native JavaScript engines. The command recomputes the certificate and verifies its exact bytes and checksum.
