# LAR-1 — Entire source return

The LAA source Cauchy bound makes the LWP exponential seed series an
entire function without assuming all-order positivity. This packet
derives a rational error bound for that actual infinite source function,
then uses every retained coefficient c_1,...,c_19 in a degree-36 polynomial
and bounds all remaining contributions. No fresh source derivative is used.

- [Theorems and proofs](THEOREM.md): entire return, coarse and refined
  remainder bounds, conditional positive realization and a counterexample.
- [Source evaluation](SOURCE_EVALUATION.md): three retained actual anchors.
- [Source map](SOURCE_MAP.md): inherited hypotheses and classical reference.
- [Execution](USAGE.md): GitHub-only generation and subsequent replay.
- `analytic_adapter.cjs`, `verify.cjs`, `SOURCE_PINS.json` and `run.py`:
  native rational/interval execution and byte-pinned dependencies.
- `CERTIFICATE.json`, `EXPECTED.sha256`, `CI_PROVENANCE.json`: outputs
  generated on GitHub and imported without local scientific execution.

The analytic function is unconditional under the retained source bounds.
Its unitary observer interpretation retains the inherited all-order
positivity condition. RH and physical interpretation remain open.
