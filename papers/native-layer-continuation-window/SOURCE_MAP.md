# LCW-1 source map

Parent: Publications commit 1e62ad12a5a2e4c61643ae1f636ca694a3abd0a0,
ILC-1, draft PR #81. Every predecessor file remains unchanged.

| Source | Use | Boundary |
|---|---|---|
| ILC-1 theorem and certificate | Actual single-point normalization, old 16 layers, finite response, cut dilation, all-order RH equivalence | Universal source signs were not proved |
| ILC source_layers.cjs | Fresh complete derivatives 0..17 at s=29 | Retains the original Euler/Gamma/pole/tail protocol |
| ILC layer_adapter.cjs | Both Hankel forms, binomial readings, normalized finite response | Finite positive measure is a moment model |
| PIW window_adapter.cjs | Outward LDL solves and inverse quadratic forms | No source positivity correction |
| POC observation_adapter.cjs | 110-place outward rational interval engine | Inconclusive positivity is rejected |
| Frozen native_operator.cjs | Exact F/Cut arithmetic and matrices | Scientific fixtures use no floating point |
| Frozen paninian_operator.cjs | Audited scalar rewrite identities and replay | Symbolic checks supplement the written general proof |
| Dette and Nagel (2012), Section 2, equations (2.1)–(2.2) | Historical attribution of canonical moment positions | Not a source-specific xi inequality |

Primary reference: H. Dette and J. Nagel, *Distributions on unbounded
moment spaces and random moment sequences*, Annals of Probability 40,
2690–2704, DOI 10.1214/11-AOP693,
<https://arxiv.org/pdf/1007.3369>.

LCW-T1–T3 supply the finite Gram/Schur derivation in full. ILC's
analytic source evaluation and error theorem are reused with degree 17;
SOURCE_PINS.json pins their text, executable dependencies and original
certificates. No source integral with an unproved positive xi-zero
measure is used to justify the new actual signs.
