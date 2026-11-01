# TLM-1 — TVSP lambda-layered matrices

Every entry carries a source-generated local response record and an observer
reading. Its layers come from the existing TVSP lambda/Jacobian tower.
The normalized scalar jets support associative matrix composition and exact
ordinary-matrix recovery; declared context reconstructs the covariant tower.

The concrete source example keeps the observer reading **5** while the
selected curvature is **5 mu²/324**. At mu=5 it is **125/324**.
Even all derivatives of that one reading cannot recover the curvature;
the first layer of the surrounding TVSP response supplies the missing data.

A single six-by-six matrix displays the point response and its first layer.
Its entries distinguish the deformation although its eigenvalues do not.
This is an algebraic record, not a positive self-adjoint quantum observable.
Discrete layer order alone does not establish physical space quantization.

- [Definitions and proofs](THEOREM.md)
- [Source and claim ledger](SOURCE_MAP.md)
- [Reading and using an entry](USAGE.md)
- [Native response adapter](tvsp_adapter.cjs)

Replay with Python 3.12 and Node:

    python3.12 papers/native-tvsp-lambda-matrices/run.py --check

Python performs I/O and hash checks only. Exact scientific operations use
the unchanged frozen native engines. No engine is forked or modified.

This packet implements the user's corrected TVSP proposal. The paused
generic C direct_sum M3 recursion was not published and is not used here.
