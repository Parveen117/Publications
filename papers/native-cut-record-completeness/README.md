# CRC-1 — Which cut information must a carrier retain?

14 October 2026. Continuation of the invariant/rotation discussion and
EGR-2's nonlinear geometry. This packet reuses the existing native
row-closure implementation, observer-kernel criteria and history results.

Read [THEOREM.md](THEOREM.md) for the precise targets, proofs and lineage.

Three concrete results:

1. For the declared real-plane cut readings, the record
   `(N,z,x)=(Tr rho, Tr Z rho, Tr X rho)` closes under every admitted
   rotation. Amount plus the first cut reading is insufficient; one
   exact continuation exposes the missing direction. Two calibrated
   cut orientations supply the three readings when both outcomes of
   the first cut are retained.
2. If an unknown additional native phase component is admitted and a
   transport mixes it into the observed channel, the required linear
   record becomes `(N,z,x,y)`. Native closure ranks change from `2 -> 3`
   to `2 -> 4`. If that component is known to vanish initially, no new
   unknown initial parameter is created by transport.
3. In the same nine-cell gravity model, two explicit source laws have
   the **same source cut, amount, source gradient and first-order curvature
   response**, but different second-order response and opposite metric
   stability. The missing law datum is a quadratic source-coupling term.

Cut count, ordering and lifted angle are separate records. They are needed
for targets that depend on them. An identical endpoint state gives identical
readings after identical memoryless continuations, so retaining every past
operation is not a universal requirement for predicting an endpoint reading.

## Reproduce

From the repository root:

```sh
python3.12 papers/native-cut-record-completeness/run.py --check
```

Python performs I/O and starts Node. The original native operator engine
and Paaninian compiler are extracted from the pinned repository archive.
No new arithmetic implementation or third-party mathematical package is
used. The native `rowClosure` and `descendedAction` routines are reused.

[CERTIFICATE.json](CERTIFICATE.json) records the closure bases and operation
witnesses, exact source reconstructions, frame comparisons, history controls,
both gravitational source-law jets, symbolic proofs and serialized replays.
Read-only `--check` reconstructs it byte for byte; `--write` generates it
after deliberate source-pin updates.

Minimality is among linear records for the declared source/operation class.
Readings are exact native traces, not one-shot experimental state recovery.
The general all-word criterion, analytic branch existence and local inertia
use the written proofs. Exact checks do not amount to a generic formal
proof-assistant certificate or a complete quantum-gravity theory.

No older packet or engine is changed. The comparison does not choose a
physical coupling or replace EGR-2's original source action. It identifies
which additional law information is needed before that choice can be tested.
