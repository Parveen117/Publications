# LDC-1 native replay

From the repository root:

```sh
python3.12 papers/native-layer-depth-compatibility/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

Python 3.12 verifies file hashes and invokes the native verifier.
Scientific arithmetic uses the frozen F/Cut and inherited outward
Interval engines. The verifier reuses pinned complete source jets;
the integrated chain also reruns their original source evaluations.

The replay checks exact operator fixtures, positive response
innovations, fixed-length observer records, deletion couplings,
formal projection identities and rejected stronger claims.
CERTIFICATE.json records deterministic output and EXPECTED.sha256
records its digest. Source files and retained dependencies are pinned.

After intentional reviewed changes, regenerate with:

```sh
python3.12 papers/native-layer-depth-compatibility/run.py --write
```

The finite certificate does not establish the all-order premise of
the conditional infinite record limit. Its physical interpretation
remains open.
