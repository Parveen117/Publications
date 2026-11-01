# POC-1 — Equal-reading centres and prime-sampled source recovery

16 October 2026. Continues POS-1 in the same Publications research chain.

Three noncollinear points have a unique equal-squared-distance centre
in the retained positive TET metric. Prime-content inputs need not give
a prime-content centre: the explicit inputs (2,0), (0,3), (5,0) give
(145/54,44/27), outside the adapted integer lattice. The correct local
rule is determinant compatibility: every prime not dividing det A
has one compatible centre at every prime-power depth. Bad-prime first
layers can fail to lift.

The actual RH source supports a separate positive development. The
SE48 spectral pencil and SE49 innovation/convergence theorem extend to
observations at odd-prime original s values. Prime-reciprocal divergence
gives sufficient analytic uniqueness density, with compact-uniform
Schwarz-product decay bounded by a power of 1/log(p_N+1). Consequently
RH is equivalent to positivity of both actual pencils at every prime
bank order, or both actual Schur innovations at every step. Those
universal signs remain open.

An exact outward replay selects s=3,5,7,11,13,17,19,23 from SE49's
complete 24-direction source archive. Both 8x8 pencils are positive.
At the independent held-out s=27 reading, the response error is less
than 10^-15. This is a new selection and replay of existing source
data, not a larger source space, a global eta improvement, or RH proof.

- [Theorems and analytic proofs](THEOREM.md)
- [Source citations and claim ledger](SOURCE_MAP.md)
- [Replay and observation contract](USAGE.md)
- [Source pins](SOURCE_PINS.json)
- [Exact certificate](CERTIFICATE.json)
- [Native adapter](observation_adapter.cjs)

```bash
python3.12 papers/native-prime-observation-closure/run.py --check
```

The implementation uses the unchanged native rational/Cut matrix and
Paninian rewrite engines. Python 3.12 handles file integrity and process
orchestration only. The finite arithmetic witnesses supplement the
written universal and analytic arguments; they do not prove those
arguments by enumeration. No map from a centre to a zeta zero, a
source-selected physical mass, or quantized spacetime is supplied.
