# Execution contract

Scientific execution for this continuation runs on GitHub Actions only,
as requested on 17 October 2026. Local work is limited to derivation,
editing, source inspection, Git operations and artifact/hash I/O.
Do not duplicate passing GitHub replay locally or launch an independent
audit agent without a new request.

The integrated runner has a narrowly scoped bootstrap option:

```bash
python3 scripts/replay_quantum_geometry.py --bootstrap-stage native-layer-analytic-return
```

On GitHub only, if both new certificate files are absent, it checks all
58 predecessors and executes this stage with `--write` once. It then
exports the certificate, expected hash and CI provenance. If both files
exist, the same invocation uses the usual `--check` mode. One missing
file is an error; an existing certificate is never silently replaced.
Without this option the integrated runner remains read-only.

The source PR triggers one GitHub run. Its immutable artifact is imported
with an output-only `[skip ci]` commit, avoiding another identical replay.
The provenance identifies the PR head, checked-out commit/tree, run ID,
59-stage reports and output hashes. The final import commit is not itself
reported as a freshly tested head. Future substantive changes still
require their normal GitHub checks; this procedure does not waive merge
requirements or mark a skipped check as passed.

Native arithmetic stays in JavaScript; Python 3.12 performs file hashing,
process orchestration and provenance I/O only. `SOURCE_PINS.json` covers
the fixed packet and inherited dependencies. Read the theorem for every
infinite or positive-realization assumption.
