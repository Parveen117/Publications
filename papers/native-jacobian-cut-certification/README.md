# NJC-1 — Native certification of the Jacobian uploads

The eight uploaded files have been read and audited. Their useful algebraic
core is recorded in [THEOREM.md](THEOREM.md); every substantive source claim
is classified in [CLAIM_AUDIT.md](CLAIM_AUDIT.md). Original manuscripts are
unchanged. The two circular–hyperbolic uploads are duplicate pointers to one
MP commit; the tower pointer resolves to RKF T42.

This packet reuses the existing cut and connection calculus. It certifies
corrected formulas, records explicit counterexamples, and identifies the
normal-plane data still needed by the gravity construction.

| Retained result | Immediate use |
|---|---|
| Reflection residual `Jz-z=-2P_-z` and exact sine form | Dimensionless mismatch observable |
| Typed projector/tensor decomposition | Consistent observer, retained and transverse components |
| CID-D1 connection `A_i=J partial_i J/2+B_i`, `[B_i,J]=0` | Separates cut-determined transport from free sector transport |
| Minimal-connection curvature `[partial_i J,partial_j J]/4` | Exact curvature from two changing cut directions |
| Quadratic-generator circular/parabolic/hyperbolic bridge | Correct normal rotation and boost branches |
| Response tower and corrected blind-target criterion | Tests whether additional derivatives actually reveal the target |
| Bounded relaxation and corrected memory propagator | Supplies lawful dynamics for the declared damping model |
| Constant-cut/nonzero-normal-curvature counterexample | Shows why a cut tower alone cannot recover independently allowed normal transport |

The source Einstein equation, universal entropy/irreversibility claims,
quantized dissipation and Millennium conclusions do not pass as derived
theorems. The audit states the concrete algebraic, type or missing-hypothesis
reason in each case.

## Reproduce

From the Publications repository root:

```bash
python3.12 papers/native-jacobian-cut-certification/run.py --check
```

The launcher uses only standard-library I/O, hashes and process execution.
It extracts the two pinned, unchanged native JS modules from the source ZIP
already present in this repository. All new certificate arithmetic and
word proofs run in those modules. Node must be available.

`--check` regenerates the mathematical payload in memory, replays its
serialized witnesses, and compares it byte-for-byte with `CERTIFICATE.json`
and `EXPECTED.sha256`. `--write` is an explicit certificate-generation mode.
Packet documents and implementation are hash-pinned by `SOURCE_PINS.json`.

The original RKF tower was separately replayed with Python 3.12: ten tests
passed and the historical payload hash matched
`3d8db4f085a3b624192c45b9f6e0356b600834193967bc89478381b53ae44a3a`.
That Fraction-based historical check is not counted as a new native proof.
The old MP SymPy suite was inspected; the new native port supplies the
certificate used here.

## Meaning of PASS

`PASS_CERTIFIED_CORE_WITH_RECORDED_REJECTIONS` means that the declared native
algebra identities, replay witnesses and exact finite controls passed,
including rejection of specified false inferences. It is not blanket
approval of the source manuscripts. The general analytic statements retain
their written completion/domain hypotheses, and physical applications need
their stated adapters. This is not a general-purpose proof-assistant or
experimental quantum-gravity certificate.

The exact counts, proof objects, presentations, negative controls and
four-dimensional curvature example are in [CERTIFICATE.json](CERTIFICATE.json).
