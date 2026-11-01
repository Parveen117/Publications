# ILC-1 — Single-point infinite layers and a positive cut

17 October 2026. Continuation of PIW-1 in Publications.

At one regular source point a, define dimensionless derivative layers
c_n=(-a)^n q^(n)(a)/(n!q(a)). Their two Hankel forms are the
confluent limits of the retained L/J observations. The all-order
condition has two equivalent readings for the actual source:

- both Hankel forms stay positive at every depth;
- all alternating binomial readings D_(r,k) are nonnegative.

The latter readings give normalized binary-record weights. Under
the universal condition the acquired polynomial pairing determines
a positive contraction Y. Its finite compressions yield cuts that
converge strongly to the two-channel projection

\[
\Pi=\begin{pmatrix}Y&\sqrt{Y(I-Y)}\\
\sqrt{Y(I-Y)}&I-Y\end{pmatrix}.
\]

Two reflections give an exact norm-preserving rotation and a wave
recurrence. A co-transported observer retains its reading. Infinite
derivative depth, repeated wave transport and averaging are separately
defined operations; their limits must not be conflated.

The actual execution evaluates 16 complete derivative layers at
original s=29. Both 8x8 forms and all 136 available binomial readings
are positive. At independently evaluated, unfitted s=31 the local
response error is positive and below 6e-28. Every omitted analytic
tail is retained. This finite result does not establish the universal
source condition, which is equivalent to RH.

- [Five theorems and proofs](THEOREM.md)
- [Complete higher-source evaluation](SOURCE_EVALUATION.md)
- [Source and claim ledger](SOURCE_MAP.md)
- [Usage and replay contract](USAGE.md)
- [Source pins](SOURCE_PINS.json)
- [Exact certificate](CERTIFICATE.json)

```bash
python3.12 papers/native-infinite-layer-cut/run.py --check
```

Classical compact moment, projection-dilation and spectral methods
are explicitly credited and rederived for this source. The wave is
not identified with an electromagnetic field, and no physical mass
or quantum-gravity closure is asserted.
