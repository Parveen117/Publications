# POS-1 observer contract and replay

Run from the repository root:

```bash
python3.12 papers/native-prime-observer-sectors/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

The first command verifies all source pins and replays the exact native
certificate without rewriting it. `--write` is the deliberate certificate
regeneration mode, used only after updating this packet's own reviewed
source pins. It never edits predecessors. The integrated runner retains
all earlier stages and adds POS-1.

## Coordinates and interpretation

| Object | Meaning |
|---|---|
| Native Lambda=S Z^2 | The lattice already chosen in NM1; S=I+R |
| Adapted (a,b) | Integral coordinates in the S basis; reduce these modulo n |
| V_n | Finite additive observer with n^2 labels |
| Prime n | No proper nonzero additive subobserver invariant under every modular generator |
| p^k resolution | Nested congruence depth; each refinement adds p^2 possibilities |
| Coordinate content | gcd(a,b), invariant under modular transport; prime-depth zero readings expose its valuations |
| Finite complex register | A declared permutation realization; not an LMD physical-sector identification |
| Co-moving TET lattice | U_gamma Z^2 in the tangent frame; S U_gamma Z^2 in original native coordinates |
| Positive scale a | Chosen energy in the countermodel; not derived mass |

Transport stability here is persistence of the observer's equivalence
relation under common continuation. The reading may be permuted. This
contract does not assert state-vector stationarity, dynamical robustness,
physical particle stability, or constant energy.

## Implementation surface

`observer_adapter.cjs` uses exact native F/Cut matrices and exact integer
residues. Matrix inputs with fractional/complex entries are rejected where
integral adapted action is required. Moduli, prime bases, digit depths and
finite enumeration sizes have explicit guards; resource bounds limit only
the executable fixtures, not the written theorems.

The adapter supplies the frozen-sign source matrices, intrinsic modular
action, exact lattice-invariance tests, bounded HNF enumeration, modular
closure, coprime CRT maps, prime digits with signed carries, permutation
registers, coherent coarse records, uniform layer embeddings, selected
Hamiltonians, and the corrected TET frame bridge. Inspect its exported
functions and certificate for the exact finite ranges exercised.

No rational matrix is reduced modulo a prime by silently inverting a
denominator. In particular S is singular in ambient characteristic two,
although the intrinsic Lambda/2Lambda quotient is valid and has four
elements. The change of basis is performed before residue reduction.

The supplied CRT inverse requires coprime moduli. Digit carrying uses
Euclidean residues, including when the rotation produces negative
coordinates. Independent digit transport is deliberately rejected.

## What can be reused next

The elementary-observer theorem and its co-moving hierarchy are available
as exact algebraic inputs. A proposed physical mass law must specify how
these records couple to the existing material sector and how their energy
is selected. A fixed-lattice interpretation in the curved source must
address the actual holonomy obstruction. The packet does not prescribe
either missing rule by relabeling an invariant or a spectral parameter.
