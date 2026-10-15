# CCL-1: Common cone, rest invariant and the geometry interface

15 October 2026. Conditional continuation of PR1 / PR2 / IN1, reusing MG / SM and the local MGH-1 example.

The written theorem establishes a band-limited convergence estimate for the shared cut-shift scaling, a common support cone under bounded local sector mixing, and a conditional curved-coframe principal-symbol identity. Physical universality, source-to-coframe reconstruction and interacting quantum gravity remain open.

The repository source is [THEOREM.tex](THEOREM.tex). Build with `pdflatex THEOREM.tex` if a PDF is needed. The check report, unchanged native engines and pinned source snapshots live beside the theorem.

The next stage is [GCP-1](../native-geometry-controlled-propagation/README.md): geometry enters the derivative coefficient and determines the branch connection and curvature, with an exact quantum-to-single-metric reduction criterion.

Replay supporting exact checks with Python 3.12 and Node.js:

```sh
python3 run.py --check
```

`CHECKS.json` separates 69 exact native-arithmetic fixtures and 3 negative controls from floating-point dispersion illustrations. It is not a proof-assistant certificate for the analytic theorems. The unchanged native engines are bundled and hash-pinned. Source snapshots and local predecessor hashes document lineage; previous packets are not altered. Literature-wide novelty is not claimed.

The PDF is a research draft. Its main new interface is shared refinement plus local coupling; the existing boost algebra, classical gravitational action and finite geometry entanglement demo are reused rather than restated as new discoveries.
