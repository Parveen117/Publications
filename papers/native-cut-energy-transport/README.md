# Native cut transport and internal response clock

Monty Dabas · 30 September 2026

[CT-1–CT-7](THEOREM.md) continue the coherence-first thermodynamic model with
one explicit process hypothesis: multiply both reference-subtracted conjugate
responses by the same positive native scalar a≤1.

From that hypothesis and the previous energy we derive:

- A unique global state semigroup F^-1(aF(x)), where F=grad W.
- A native additive clock -LogΣ(a), still uncalibrated to physical seconds.
- Strict energy descent and a positive inward radial rate.
- A native amplitude-plus-rotation transport and exact moving-cut memory law.
- Radius-dependent angular drift under anisotropy, with the principal-response
  seam's path derivative kept separate from the contrast angle.
- A local variational characterization of the derived differential law.

At ν=3/2, κ=2/3 and ε=0, halving conjugate responses takes r from 1 to 1/4,
W from 2/3 to 1/12, and det H from 1/2 to 2. The interval is native LogΣ(2),
whose rational Cayley series is enclosed with a retained positive tail.

This is a conditional model and mathematical existence/uniqueness proof, not
a unique physical law selected by coherence. A pure dagger-unitary conjugation
cannot reduce D²; the amplitude-changing part is explicitly nonunitary, but
invertible for every finite positive radius. It is not implicitly a closed
quantum-system evolution or a measured energy reservoir. Gravity, physical
seconds, spatial distance and a universal constant remain unestablished.

## Reproduce

With Python 3.11/3.12 (without -O), Node, and the canonical RKF checkout at
`3cc5a33b05c16d59c90994ddda69dedc0d392424`:

```bash
node papers/native-cut-energy-transport/verify.cjs --rkf-root /path/to/Recognition-Kernel-Framework --check
```

Default/check mode is read-only. Explicit --write generates this packet's
source-bound certificate. The verifier imports the canonical native engine
and original thermo model; no alternate arithmetic engine is introduced.
Source pins include the CF energy theorem, native propagation/memory theorem,
native logarithm and the Morphic clock definition. No upstream files change.

The finite controls check algebraic transports, composition, winding parity,
aligned and mismatched cut memory, Hessian-derived trajectories, energy descent,
angular drift, principal seam derivatives and rational clock-tail ledgers.
They supplement the written general proofs; no proof-assistant or independent
review status is claimed. The pre-existing upstream master-pin limitation in
[native propagation](../native-ugd-propagation/README.md) remains unchanged.
Configured GitHub CI does not itself establish a hosted PASS.
