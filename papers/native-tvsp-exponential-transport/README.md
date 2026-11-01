# TET-1 — TVSP layers generate exponential observer transport

The existing TVSP source now supplies an exact bridge from its retained
response geometry to ordered transport, a moving observer and a positive
observation cut. This continues TLM-1 rather than introducing a separate
layer model or operator engine.

- Source H layers through depth N determine U layers through depth N.
- U^T H U=H_0 and O=O_0 U^(-1) preserve the co-moving pairing.
- The positive cut P=U P_0 U^(-1) satisfies the existing EGR contract after
  conversion to an orthonormal frame.
- Straight paths of the TLM source have U=(I+t H_0^(-1)D)^(-1/2).
- On gamma(t)=(t,t²), the first error from discarding exponential ordering
  is t³ F_sv(0)/6. Small-loop return also recovers the same curvature.
- A fixed base observer can see that return; a co-moving observer cancels it.
- Actual stable TVSP sources show why finite layers alone do not determine
  a nonzero endpoint. A residual bound gives controlled approximation when
  interval bounds are supplied.

Read [the six theorems and proofs](THEOREM.md), [usage and observer
interpretation](USAGE.md), and [source/claim ledger](SOURCE_MAP.md).

Replay from the repository root:

```bash
python3.12 papers/native-tvsp-exponential-transport/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

The certificate passes 302 exact native checks, 7 symbolic identity replays
and 43 negative controls. Universal analytic claims rest on the written proofs;
the checks are not a general proof-assistant certificate. The physical path,
clock and connection to a continuum quantum geometry state space remain
separate obligations. No physical space quantization is inferred.
