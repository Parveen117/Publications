# IRT-1: Dimensionless invariant ratios under transport

The coupled AGD model has an exact reflection invariant. With `S phi(theta)=phi(-theta)` and `Pi_±=(I±S)/2`, the weights of these two geometry sectors are conserved even with the full spatially varying matter coupling:

\[
\boxed{\mathcal I_{+:-}=\frac{\operatorname{Tr}(\Pi_+\rho)}{\operatorname{Tr}(\Pi_-\rho)}=\text{constant}.}
\]

The denominator must be nonzero. This reflection is different from the half-turn translation whose conservation failed in AGD-1.

The written infinite-dimensional proof also classifies the bounded geometry observables commuting separately with the geometry generator and lapse:

\[
\mathcal A_{\rm inv}=\operatorname{span}\{I,S\}.
\]

Thus the two sector weights give the complete record for this specified invariant algebra. They do not determine the complete quantum state or all future nonconserved readings.

For finite linear transport, a nontrivial ratio is conserved for every initial source exactly when its two observables have a common dual generator rate. For unitary transport this rate must be zero, so both observables separately commute with the generator. A separate exact attenuation channel illustrates nonzero common rates with a retained erasure record.

The diagonal statement is calibrated precisely:

\[
\frac{n'}n=\frac{d'}d
\quad\Longleftrightarrow\quad
\log\frac{n(s)}{n(0)}=\log\frac{d(s)}{d(0)}.
\]

This is a 45-degree line on equally scaled Euclidean log axes. Independent axis rescaling changes its plotted slope. Constant-rate exponential scaling additionally requires a continuous multiplicative semigroup; the diagonal alone does not impose that rate law.

| Exact result | Value |
|---|---|
| AGD example: reflection-sector weights | `1/3`, `2/3` |
| Preserved ratio | `1/2` |
| Dimension of the common bounded geometry observable algebra | `2` |
| Two states with the same invariant sector record: lapse variances | `9/200`, `27/400` |
| Attenuation example: retained population / survival | `(36/125)/(9/25)=4/5` |
| Separately retained erasure weight | `16/25` |

**Executed verification:** 103 exact checks, 9 native word certificates, 9 serialized replays and 16 negative controls. Finite rank checks use compressed `L²` and cosine. An explicit control rejects extrapolation from an accidental three-mode Hamiltonian symmetry to the infinite operator. The infinite classification and general ratio theorem have written proofs; they are not proved by finite sampling.

The result identifies conserved information within the selected model. It does not derive a physical gravitational rate, a universal physical 45-degree angle, or local gravitational constraint closure. The action, quantum adapter and coframe dictionary retain their declared assumptions.

- [THEOREM.md](THEOREM.md): derivations, classification and calibrated-angle scope.
- [CERTIFICATE.json](CERTIFICATE.json): exact witnesses and rejected stronger claims.
- [SOURCE_PINS.json](SOURCE_PINS.json): unchanged engines and source provenance.
- [Canonical research record](../../QUANTUM_GRAVITY_RESEARCH.md): next gravity obligation.

From the repository root, using Python 3.12 and Node.js:

```sh
python3 papers/native-invariant-ratio-transport/run.py --check
python3 scripts/replay_quantum_geometry.py
```

Python performs I/O and hash verification. Scientific arithmetic and symbolic replay use the unchanged native engines.
