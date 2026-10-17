# GitHub execution and predecessor reuse

No local scientific replay or independent audit agent is used.
The source branch's GitHub job invokes ci.py. It verifies the successful
predecessor run, its exact source commit/tree, the output-only import lineage,
the 59-stage report, and the absence of changes to every predecessor source.

Only then does it execute this new stage. Thus the result is **59 reused
stages plus one new stage**, not a claim that all 60 stages were rerun.
The ordinary integrated runner still supports a full 60-stage read-only
replay. Changed predecessor sources are rejected by this incremental route.

PIN_INPUTS.json declares the files to pin. Initial pins are generated in
GitHub from the immutable checked-out source; run.py checks them before the
native JS verifier emits its certificate. Existing output sets are checked,
never silently replaced. Partial output sets are rejected.

The GitHub artifact contains SOURCE_PINS.json, CERTIFICATE.json,
EXPECTED.sha256 and CI_PROVENANCE.json. Importing them is file/hash I/O only
and uses an output-only commit without another scientific run. Provenance
identifies the tested source commit; it does not relabel an import commit
as a freshly tested head.

Python 3.12 performs orchestration, hashing and GitHub API/file I/O only.
All scientific fixture arithmetic is performed by the existing native
JavaScript engines. Infinite claims are written proofs in THEOREM.md.
