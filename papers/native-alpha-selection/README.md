# Native fine-structure selection audit

Monty Dabas · AS-1–AS-5 · 30 September 2026

**Alpha is not numerically derived.** [The written result](THEOREM.md) tests
whether actual native aperture closure and the newly developed internal clock
select the electromagnetic coupling. It imports the original R2 depth-return
solver and canonical native operator/weighted-completion engines.

The exact-cut family (a,b)=(b+1,b) has the same aperture for every b>0, but
its response slope is 1/(1+2b). For periods (2,1) and (3,2), the slopes are
1/3 and 1/5. Exact jets and retained-tail aperture enclosures distinguish them.
Thus exact cut closure alone has not fixed interaction sensitivity.

The source-complete reduction keeps both Z=A-BC^-1B† and the dressed source
qeff=qv-BC^-1qh. Full inversion verifies the pole residue and separate contact
term of an admitted quadratic model. The native noncommutative WC5 API also
checks an example where erasing the hidden source changes the visible state.
These are not automatically an electron/photon identification.

The missing law is a selection of interaction strength relative to charged
source and action normalization. The theorem states a concrete acceptance
contract for an alpha prediction, without using its observed value or tuning
parameters to it. A negative result for the audited premises is not a claim
that the full framework can never select a coupling.

## Reproduce

```bash
node papers/native-alpha-selection/verify.cjs --rkf-root /path/to/Recognition-Kernel-Framework --check
```

Use the pinned RKF commit 3cc5a33b05c16d59c90994ddda69dedc0d392424. Check mode
is read-only; --write explicitly generates this packet's source-bound evidence.
No packages or alternate arithmetic engine are installed. Written proofs,
finite exact controls, physical matching and independent review are separate.
The documented upstream master-pin limitation is unchanged; this scoped
certificate is not a recertification of the whole upstream framework.
