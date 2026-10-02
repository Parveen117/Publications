# Native compact gauge completion

The [eight written results](THEOREM.md) connect the existing native massive
module and its U(1) phase to an explicit U(2) internal symmetry with an SU(2)
sector. They supply the positive invariant pairing, native non-Abelian
curvature, classical source/Noether/Gauss identities, and a complete
thermodynamic coefficient-variation interface.

The result is conditional on a declared replicated module, normalized massive
multiplier and local action class. It does not select nature's gauge group,
fermion statistics, coupling constants, continuum measure or mass gap.
The [source map](SOURCE_MAP.md) distinguishes existing results from this bridge.

Run from a Publications checkout with Python 3.12 and Node 24.19.0:

```bash
bash papers/native-compact-gauge/verify.sh \
  --rkf-root ../Recognition-Kernel-Framework \
  --physics-root ../Publications-physics --check
```

Use RKF commit `3cc5a33b05c16d59c90994ddda69dedc0d392424` and the separate
Publications physics checkout at `f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56`.
Every consumed file and local input is bound in [SOURCE_PINS.json](SOURCE_PINS.json).
The verifier checks pins before executing the independent controls. Its
read-only result must equal [CERTIFICATE.json](CERTIFICATE.json) and
[EXPECTED.sha256](EXPECTED.sha256). `--write` is the explicit maintainer path
for a newly reviewed packet, never a substitute for a failing check.

The canonical RKF engine is consumed unchanged. `presentation.json` only
declares two independent EMK factors. `controls.py` is an independent rational
matrix/coordinate verifier, not a second native operator engine. The certificate
separates exhaustive finite linear problems, sampled variable-field controls,
native derivation replays and the eight written proofs. It does not claim a new
proof-assistant certificate or physical validation.
