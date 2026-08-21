# CID-1 — Curvature–information duality: the obstruction is metric-free

A **bridge capsule**. CFE-U certified that the obstruction to the
memoryless limit is unique (a 1-dimensional obstruction space spanned by
the loop residue). This capsule asks whether the *information* side sees
the same obstruction — or a second, independent one.

**The answer is a separation with an invariance:** the Hodge/Helmholtz
split of a response *depends* on the information metric; the obstruction
*does not*. The metric is a gauge on the decomposition and an invariant
on the obstruction.

| Block | Certifies | Result |
|---|---|---|
| **T1** | The information metric at a declared rational point is g = Cov_p(T), exact — **no exponential or logarithm evaluated anywhere**. PD ⟺ statistics affinely independent; degenerate with explicit kernel otherwise; constants always in the kernel (matching CFE-U's ker d₀) | PASS |
| **T2** | **Recognition ledger = law of total covariance:** g = g_recognized + g_discarded exactly (between-block + within-block), both PSD. Monotonicity follows as an *identity*, not an inequality: recognition can lose information, never manufacture it. Sufficient partition discards exactly nothing; total collapse recognizes exactly nothing | PASS |
| **T3** | **The obstruction is metric-free:** curvature is d, not ∇ — no Christoffel of any metric enters. Euclidean, anisotropic, sheared and the information metric Cov_p(T) itself all give the same residue, reproducing CFE-1's ladder −32…+32 | PASS |
| **T4** | **The split is not:** exact rational least squares gives ω = dΦ_g + h_g, genuine for each metric (closed part memoryless with an explicit potential, remainder g-orthogonal), and the splits **genuinely differ** — yet ∮h_g = 32 for every metric | PASS |
| **T5** | For every metric the g-orthogonal complement of the memoryless space has rank **exactly 1** — CFE-U's uniqueness seen through the information metric. Generators differ as vectors, all land in the same obstruction class. Control: a degenerate PSD form gives no well-posed split | PASS |

## Findings

- **CID-F1.** Curvature and information are **dual descriptions of one
  object**, not two obstructions. Every positive-definite information
  metric reproduces CFE-U's rank-1 obstruction space exactly.
- **CID-F2.** The law of total covariance *is* a recognition ledger:
  what the Eye map keeps and what it discards sum exactly to the whole,
  with both parts PSD.
- **CID-F3.** Monotonicity of information under coarse-graining is here
  an **identity**, which is strictly stronger than the usual inequality.

## Claim boundary

The exponential family's form is **declared**; only rational-point
arithmetic is certified. Continuum information geometry, α-connections,
dually flat structures and quantum (Fubini–Study / Bures / SLD) metrics
are **not claimed**. The vault appendix named for Euler
information–curvature duality **has not been read and is not the
source** — this is a bridge between two certified capsules of this
corpus. Uniqueness inherits CFE-U's scope (the declared response
algebra). RH / K0 / L0 / YM / QG untouched.

```
python certificates/cid1_curvature_information_duality.py
python -m pytest tests -v
```
