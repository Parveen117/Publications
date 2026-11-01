# PIW-1 replay and source contract

```bash
python3.12 papers/native-prime-innovation-window/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

`--check` verifies source hashes and reproduces the exact certificate
without rewriting it. `--write` deliberately regenerates this packet's
certificate after its own manifest has been reviewed. Both retain all
predecessors. The integrated runner appends PIW to POC.

| Object | Meaning |
|---|---|
| delta | Actual new source value minus the old finite response |
| kappa | Positive coefficient calculated from the retained L/J forms and the new probe |
| varphi=kappa delta | Value position in the strict acceptance interval |
| vartheta | Allocation of the positive innovation budget to x alpha; the remaining fraction is beta |
| g_N | Square of the determinant transfer factor, carrying the known interpolation zeros |
| R_N | Actual response residual after those already-known zeros are divided out |
| Prime 29 | New actual retained direction beyond the old source archive |
| Prime 31 | Independent held-out check; not part of the nine-prime fit |

`window_adapter.cjs` expects a compatible divided-difference source bank
whose intervals enclose the retained exact L,J,u and nodes. Positivity
of arbitrary unrelated matrices is not a replacement for that source
contract. The adapter certifies the old forms, solves their positive
pencils outwardly, and reports both original and residual Schur signs.
It refuses a repeated next node, a nonpositive next node, or uncertified
old forms. It does not repair failed signs by adding a positive shift.

`full_source.cjs` freshly evaluates odd integer s in its declared finite
range. The replay uses 4096 as the prime cutoff, 128 logarithm terms,
96 arctangent terms, digamma shift 256 and 24 Bernoulli terms. Exact
tail formulas retain everything beyond those truncations. Every prime
power for a retained prime is summed in its exact rational weight.
The same native 110-place interval grid is used throughout.

The certificate records s=23 as an independent old-source check,
s=29 as the new retained jet, and s=31 as the unfitted check. It also
contains the resulting nine-prime bank for later source-bound reuse.
The four-atom rational fixtures are separately labelled toys and must
not be substituted for actual xi. No zero catalogue or RH premise is
used to generate the complete source boxes.

An interval that contains zero is inconclusive, not proof of a negative
actual sign. For the accepted new prime both the full LDL and the two
Schur calculations have strictly positive lower endpoints. The next
mathematical obligation is control of the actual value and derivative
windows at all prime steps; the finite execution leaves that open.
