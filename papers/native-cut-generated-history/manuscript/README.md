# Finalized manuscript v1.0

**Cut-Generated Information Transport and Complete Observer Verification**

*A finite mathematical framework for quantum-gravity research*

16 October 2026. This manuscript closes the finite quadratic-observer
representation and verification study. Its scientific scope is the integrated
CHD-1/2/3 theorem, not a completed physical theory of quantum gravity.
It is a finalized research manuscript, not a claim of journal acceptance,
peer review, or independent experimental validation. No submission was made.

- [Read the eight-page paper](paper.pdf)
- [Editable LaTeX source](paper.tex)
- [Evidence and artifact hashes](MANUSCRIPT_MANIFEST.json)
- [Complete source theorem](../CIRCULAR_OBSERVER_COMPLETION.md)

The paper includes the observation contract, coherent-history proof, minimum
weighted amplitude record, complete circular observer reconstruction,
minimum fixed probe count, iff transport certificate, decoder, common-frame
covariance, deterministic error bound, exact thermo example and bibliography.
It explicitly distinguishes standard mathematical ingredients, native
integration, written general proofs, finite executable certificates and
physical interpretation. The historical QGC chapter and certified predecessor
sources remain unchanged.

Evidence revision: `ac291f7aefea5b6f3d0eaf7f7679e83625341285`.
Combined CHD replay: 594 exact checks, 12 symbolic replays and 89 negative
controls. All 23 stages and both research CI runs passed at that revision.

## Reproduce

From the repository root, with Python 3.12:

```sh
python3 papers/native-cut-generated-history/manuscript/build.py --check
python3 papers/native-cut-generated-history/run.py --check
python3 scripts/replay_quantum_geometry.py
```

The first command checks manuscript/proof-source hashes only; it does not
prove a theorem. The latter commands execute the existing native certificates.

To rebuild the PDF, install pdfLaTeX with the packages used by `paper.tex`:

```sh
python3 papers/native-cut-generated-history/manuscript/build.py --build
```

The build uses a temporary directory for auxiliary files and copies only the
resulting PDF here. A changed typesetting environment can change PDF bytes;
review any changed artifact before updating its manifest hash. PDF production
does not modify scientific certificates. Python is used only for I/O and
orchestration; scientific calculations remain in the native engines.

## Closure statement

Complete within the declared finite quadratic/linear observation contract:
representation, retained information, observable reconstruction and invariant
transport verification. Independent physical gravitational dynamics, local
quantum constraints and physical FTL remain open research claims. They are
not counted as completed by this manuscript.
