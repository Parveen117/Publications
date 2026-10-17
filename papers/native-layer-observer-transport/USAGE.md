# LOT-1 replay

From the repository root:

```bash
python3.12 papers/native-layer-observer-transport/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

The first command checks every source pin, recreates all exact native
and outward interval fixtures, replays the symbolic certificates, and
requires byte-for-byte equality with CERTIFICATE.json and EXPECTED.sha256.
Python 3.12 performs I/O and process orchestration only.

`transport_adapter.cjs` exports:

| Function | Contract |
| --- | --- |
| `weightedCut(raw,N)` | Exact rational G_N, rational cut, channel reflection and wave; both pairings must certify strict positivity |
| `exactTransport(raw,N,M)` | Strict M>N, weighted inclusion and acquired adjoint, compressed cut and wave, monic boundary and rank-one correction |
| `sourceTransport(raw,N)` | Outward h_N, W_N, canonical position, cut defect, wave loss and transport discrepancy for one certified actual transition |
| `inclusion(N,M)` | Coefficient inclusion for 1<=N<=M |
| `adjoint(T,G_old,G_new)` | Acquired adjoint G_old^-1 T^dagger G_new |

The returned `boundary` is always the first rank-one boundary
projection. It equals the full-cut difference only when M=N+1.
For M>N+1 the full difference has rank M-N. Its first boundary's
compression to the original N carrier still equals K_N. The theorem
and fixtures distinguish this first boundary from the entire deeper
cut increment.

Actual-source input remains the pinned LDC moments. Exact toy moments
are illustrative algebra fixtures, not replacements for those jets.
To intentionally regenerate the local certificate after a reviewed
change, update source pins and run `run.py --write` before replay.
