# QRT-1 — Actual quantum states realize the TVSP transport

This continues TET on LMD's existing physical Hilbert space. Two actual
commuting density observables generate normalized filtered states.
Their response vectors E have covariance Gram matrix H and satisfy

    H=E*E,
    A_i=H^(-1)E*partial_i E=(1/2)H^(-1)partial_i H,
    K(t)E(0)=E(t)U_TET(t).

K is a genuine full-space unitary obtained from the response projection
Q=EH^(-1)E*. Its bounded instantaneous generator and the preserved density
operator domains are explicit. No extra two-state carrier is attached.

An actual prepared state of the unbinned density channels has response
curvature strictly between 1/8 and 3/8. This follows from uniform moment
bounds for three small positive-probability spectral cells; the native
certificate does not pretend to evaluate the continuum density integral.
A fixed response detector after a loop has probability sin²(Phi), with
Phi the integrated oriented response curvature.
A second fixed detector at 45 degrees has probability
1/2+(1/2)sin(2Phi), resolving the orientation sign for small loops.

- [Five theorems and proofs](THEOREM.md)
- [Interpretation and replay](USAGE.md)
- [Sources and evidence boundary](SOURCE_MAP.md)

```bash
python3.12 papers/native-quantum-response-transport/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

The packet replays **123 exact checks, 7 symbolic certificates and
28 negative controls**. Its SHA256 certificate is
`0da65242c759497be548a8678f92ae1cb1d1b4be4e511613a45505dbbb136ed1`.

The connection is now realized on actual quantum response states. Its
identification with the original material geometry and a source current
that implements the chosen transport remain further obligations.
