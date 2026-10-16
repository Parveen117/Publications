# LSC-1 source map

Frozen predecessor commit: 064f593ea7aeca54ee0cb2c3822ef3f990de4362.

- LQA-1: source, reality/adjoint contract, shared zero modes, reversed
  oscillator domain and unmodified local constraints.
- LPS-1/2: auxiliary c=7 factors, equal-weight descendant pairing,
  regular physical norm, observable matching and local classical reduction.
- HQS-1: the specific physical Q expression, full-line matching and the
  established self-adjoint observable representation.
- oscillator_fixture.cjs: frozen operator routines copied verbatim from
  the contiguous helper section of LPS verify.cjs, with exports added.
  This avoids modifying or re-running a predecessor's certificate as an
  implicit module side effect. Its provenance is byte-pinned and checked.
- Native CCL arithmetic and Paninian engines: unchanged.

The positive-level identity, cyclicity induction, dual-word reduction,
axis-distribution decomposition and observable-compatibility calculation
are proved explicitly here. Hadamard division, the positive-distribution
measure representation and full-line boundary forms use standard analysis,
with the needed arguments provided. No Kac determinant or external
irreducibility theorem is assumed. No novelty priority is asserted for
free-field irreducibility or the established pairing mechanism.

Standard analysis reference checked for the positive-distribution step:
Patrick Gerard, *Distributions and Partial Differential Equations*,
lecture notes, Fall 2019, Proposition 1.9.11 and Remark 1.9.12:
https://www.imo.universite-paris-saclay.fr/~patrick.gerard/Distributions2019_Chap1,2,3.pdf
This supports the measure representation step only; it is not evidence for
the source-specific oscillator or local-constraint classification.
