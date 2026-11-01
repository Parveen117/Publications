# LWP-1 replay

From the repository root:

```bash
python3.12 papers/native-layer-wave-propagation/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

The packet runner checks source hashes, executes native scientific
arithmetic and symbolic replays, and requires exact certificate bytes
and SHA256 equality. Python3.12 performs I/O and orchestration only.

`propagation_adapter.cjs` exports these main operations:

| Function | Contract |
| --- | --- |
| `action(pair,kind,steps)` | Exact real polynomial action for wave, inverse wave or generator; degree rises at most one per step |
| `pairColumn(pair,N)` | Embed two coefficient polynomials into a finite carrier, rejecting an insufficient dimension |
| `correlation(u,v)` | Polynomial y*u_1*v_1+(1-y)*u_2*v_2 for real coefficient pairs |
| `functional(raw,poly)` | Outward evaluation using retained source moments; unavailable coefficients rejected |
| `finiteReturns(raw,M,K,kind)` | Strict positive source forms and native LDL solves; returns for wave or generator, including power0 |
| `returnPolynomial(k)` | Exact polynomial y*T_k(2y-1) |
| `budgets(N,M,t)` | Exact factorial-tail bounds with d=M-N; each tail ratio must be below one |
| `exponentialReading(raw,M,t)` | Common first-channel exponential return enclosure with source uncertainty and the real cosine remainder retained |

Polynomial fixture arithmetic is real rational. The written theorems
extend correlations to complex states by conjugating the left input.
The finite matrix engine retains its acquired adjoint convention.

An enclosure for exp(tL) is not an assertion that exp(L)=W. Nor does
finite source certification prove the existence of every deeper source
carrier. The corresponding scope flags are retained in the certificate.
