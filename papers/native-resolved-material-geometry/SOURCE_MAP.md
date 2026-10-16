# RMG-1 source and scope map

Frozen parent: d9c70208520f53f31d93bd94c3190f48b34acb4c.

| Input | Used here | Boundary |
|---|---|---|
| LNC-1 theorem, canonical q=log e | Material-slice coframe interpretation | No spacetime curvature from one slice |
| MRF-T3 | Exact physical r,m and gauge representative | Regular, definite-sign patch |
| MPQ-T2/T3/T6 | Inverse adapted chart, Gaussian physical space, coherent states | No pre-existing nonlinear point operators inferred |
| RGR-T5 | Selected physical norm and quotient intertwiner | No classification of all invariant distributions |
| Existing CCL native_operator.cjs and paninian_operator.cjs | Exact rational/Cut arithmetic and replayable rewrites | Finite witnesses, not machine-checked analytic theorems |

Standard method: B. Chakraborty, J.-P. Gazeau and A. Youssef,
*Coherent state quantization of angle, time, and more irregular functions
and distributions*, arXiv:0805.1847v2, section 2, equations (1)-(5):
https://arxiv.org/html/0805.1847v2 . Accessed 16 October 2026.
It describes the coherent resolution of identity, operator-valued
integration and Gaussian lower symbols. RMG derives its own covariance,
normalization and source-dependent estimates explicitly. The closed
positive form representation theorem, Gaussian moments, Sobolev algebra
property and continuous mapping theorem are standard functional analysis.
No novelty is claimed for those methods.

New written results: this source's finite-resolution geometry POVM and
domains, noise versus the frozen sharp probe, conjugate-phase distinction,
fixed-hbar non-tightness, infinite vacuum mean length and a sufficient
simultaneous classical/resolution limit. Exact native fixtures check
finite Gaussian moments, noise scaling, inverse-chart identities and
POVM dilation algebra. They do not prove measure theory, domain density,
infinite limits, curvature recovery or physical calibration.

Coordination: the external branch claude/residual-gauge-invariant-dual,
report commit 60fd188b0b8888dfe6b948d845cf800bc8a24c07, owns the separate
full residual-Diff uniqueness/completeness question. Its abelianized
uniqueness is not used as a full-group theorem here.
