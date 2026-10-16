# POC-1 observation contract and replay

Run from the Publications repository root:

```bash
python3.12 papers/native-prime-observation-closure/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

The packet command checks SHA-256 source pins and the imported Git blob
IDs, then reproduces its certificate without rewriting it. `--write`
deliberately regenerates this packet's certificate after its own source
manifest has been reviewed. Neither mode edits predecessor packets.

| Quantity | Contract |
|---|---|
| Prime-content point | Integer adapted coordinates with prime gcd; a chosen point model |
| Equal-reading centre | Three equal squared distances in one retained positive source metric |
| Prime local observer | The ambient POS module modulo p; the centre does not choose p |
| Depth r | A solution modulo p^r compatible with reduction; bad primes need lifting checks |
| Prime RH observation | Actual value and derivative at t=(p-1/2)^2, original s=p |
| Actual source | Completed xi logarithmic derivative, including all original Gamma, prime and pole contributions |
| L/J pivots | Coordinate residual signs; not spectral gaps or mass scales |
| Observation density | Divergence sufficient for analytic uniqueness in a stated half-plane; not geometric point density |
| Finite q_N pole | Pole of a fitted response; not a certified actual xi zero |

The adapter `centre(points,H)` works with native exact real rational
coordinates and positive symmetric 2x2 H. `uniqueResidue(A,b,n)` requires
unit determinant modulo n. `residues` enumerates only the declared small
moduli 2 through 256; the universal prime-power theorem is written,
not established by enumeration.

`Interval` uses native F arithmetic and a fixed 110-decimal-place grid.
Every constructor and operation rounds lower endpoints down and upper
endpoints up using integer division, including negative fractions.
`intervalLDL` accepts only symmetric interval matrices and certifies
positivity only if every pivot lower endpoint is strictly positive.
A failed interval sign can mean insufficient precision; it is not
automatically a negative actual eigenvalue.

`response(L,J,u,z)` encloses u^T(J+zL)^-1u at positive real z by an
outward LDL solve. Its dimensions must agree and its pencil must be
certified positive. The source archive is pinned independently; checks
of interval-overlap identities establish consistency, not a fresh
derivation of the special-function boxes. No binary floating arithmetic
is used for scientific values. JavaScript numbers are limited to exact
indices, dimensions and bounded exponent bookkeeping.

Two rational toy fixtures test interpolation and finite-sign masking.
Their residues and poles are declared in the certificate. They must
not be substituted for actual xi data. Repeated reading labels do not
add independent analytic conditions unless further source derivatives
are supplied and their own sign contract is proved.

The next RH obligation is both actual innovation signs at every new
prime, including all mixed entries and tails. A finite extension is
useful evidence, but the all-order statement is the remaining theorem.
The separate mass/quantum-gravity route still needs its own physical
coupling and observation contract.
