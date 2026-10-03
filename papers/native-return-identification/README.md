# Native interaction-strength identification

Monty Dabas · NI-1–NI-4 · 30 September 2026

[The proof](THEOREM.md) turns the unresolved parameter b of the exact-cut R2
family into an identifiable response quantity. Two normalized symmetric
probe readings u>1 and 0<v<1 with u+v<2 determine both b and the unknown
probe amplitude uniquely. An exact quadratic and rational bisection produce
certified parameter intervals. The unchanged native R2 solver supplies the
finite-aperture uncertainty; it is carried into the inferred b interval.

Example: readings 6/5 and 2/5 give b=5/7 and probe amplitude 3/4. These are
chosen exact controls, not empirical data. A separate differential shape
observable cancels multiplicative output gain and affine input gain.

The symmetry audit finds that identical adjacent cells cannot realize finite
exact-cut closure in this family, while equal opening/closing scalar amplitudes
leave b free. Physical reciprocity is not inferred from scalar balance.

**Identification is not a prediction from primitives. Alpha is NOT DERIVED.**
The required symmetric physical probe, electron/photon source identification
and action normalization are not established by this mathematical experiment.

## Reproduce

```bash
node papers/native-return-identification/verify.cjs --rkf-root /path/to/Recognition-Kernel-Framework --check
```

Pin: 3cc5a33b05c16d59c90994ddda69dedc0d392424. Default/check mode is read-only;
--write explicitly generates evidence for this packet. The canonical native
arithmetic and R2 solver are imported unchanged. No extra package is needed.
Finite PASS does not establish the general proof, physical identification or
independent review. The previously documented upstream master-pin limitation
remains unchanged; no full-engine recertification is claimed.
