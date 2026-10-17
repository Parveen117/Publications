# LAA-1 native replay

From the repository root:

```sh
python3.12 papers/native-layer-anchor-atlas/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

run.py is Python 3.12 I/O/orchestration only. It verifies every source
pin before invoking the native JavaScript verifier. All scientific
arithmetic uses the unchanged F/Cut engine and inherited outward
Interval implementation. There is no floating-point scientific path.

The verifier recomputes the fresh source jets at s=31 and s=33 and
the independent held-out value/slope at s=35. It checks exact tail
identities, polynomial matrix composition, replayed rewrite proofs,
source forms and tails, all six uniform chart certificates, overlap
geometry, retained-origin forecasts, and rejected stronger claims.
The old s=29 jet is loaded from the pinned LAT certificate; the full
chain independently replays that predecessor evaluation.

For intentional certificate generation after reviewing source edits:

```sh
python3.12 papers/native-layer-anchor-atlas/run.py --write
```

SOURCE_PINS.json covers this packet's sources and its frozen lineage.
CERTIFICATE.json is deterministic output; EXPECTED.sha256 records
its digest. --check requires both to match exactly. A finite positive
certificate does not prove the open all-depth or physical claims.
