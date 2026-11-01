# QGC-1 — Quantum consistency and falsifiable cut–geometry response

14 October 2026. A continuation of EGR-1, EGR-2 and CRC-1 in the native
cut/operator programme. The original arithmetic and proof engines are unchanged.

The central result is an obstruction: the selected EGR gradient relaxation
cannot be the affine quantum evolution of an unknown source while retaining
all its pure-cut trajectories. Two ensembles representing exactly the same
density operator give different averaged initial velocities. Preserving a
cut's spectrum does not resolve that inconsistency.

This packet also supplies a consistent, explicitly selected operation:
a fixed exponential filter with a unitary dilation and a complete outcome
ledger. It derives the quantum readout cost of simultaneously observing the
two linear response modes, and makes CRC-1's free nonlinear source coefficient
identifiable from two volume slopes with a third, held-out prediction.

| Result | What has been established |
|---|---|
| Ensemble obstruction | Same density `diag(4/5,1/5)`; average old velocities `0` and `(888/29375) Z` |
| Quantum instrument | Exact joint unitary; trace preservation, affine ensemble evolution and unchanged remote marginal when all outcomes are included |
| Outcome ledger | For the displayed input, success probability `109/125`; normalized successful cut `(91/109,60/109)`; unconditional source moments remain different |
| Joint readout | Added-noise bound `||N1|| ||N2|| >= chi^4/220900`, with an exact saturating finite instrument |
| Coupling reconstruction | Two initial volume slopes identify source scale and `kappa`; a third slope obeys an independent linear relation |

These results establish finite operational consistency and conditional model
predictions. They do **not** establish a quantum spacetime dynamics, a continuum
limit, a mapping to measured gravitational units, or experimental confirmation.
The certificate's scope is part of its payload. The audit does not mark the
physical quantum-gravity programme complete.

- [THEOREM.md](THEOREM.md): premises, constructions and proofs.
- [PREDICTIONS.md](PREDICTIONS.md): equations, calibration, held-out observations
  and rejection conditions; no observations are fabricated.
- [CHAPTER.md](CHAPTER.md): consolidated scientific argument and the exact
  remaining barriers to a completed quantum-gravity claim.
- [CERTIFICATE.json](CERTIFICATE.json): native witnesses and serialized proofs.
- [SOURCE_PINS.json](SOURCE_PINS.json): engine, source and predecessor hashes.

From the repository root:

```bash
python3.12 papers/native-quantum-gravity-consistency/run.py --check
```

`run.py` performs I/O and invokes the unchanged native JavaScript engines
extracted from the pinned source archive. `--check` recomputes the entire
certificate and compares both its bytes and checksum. No floating-point
arithmetic or replacement symbolic algebra is used. Written proofs cover
arbitrary finite ancillary dimension, operator-norm bounds and analytic
limits; the executable certificate checks their polynomial identities,
exact fixtures and predecessor-dependent coefficients. This is not a claim
of formal verification in a general-purpose proof assistant.

The standard ideas of quantum instruments, dilation, mixture equivalence and
joint-measurement noise are not claimed as new discoveries. The additions
are their explicit application to the existing native gravity response,
the exact counterexample to an overextended relaxation interpretation,
and the inverse/held-out calculation for the declared source-law family.
