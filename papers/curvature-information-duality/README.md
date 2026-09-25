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
information–curvature duality **is not the source** (the later comparison is recorded in LINEAGE.md) — this is a bridge between two certified capsules of this
corpus. Uniqueness inherits CFE-U's scope (the declared response
algebra). RH / K0 / L0 / YM / QG untouched.

```
python certificates/cid1_curvature_information_duality.py
python -m pytest tests -v
```

---

# QTH-1 — Quantum recognition information: the ledger becomes an inequality

The question CID-1 left open, answered. CID-1 certified the **classical**
recognition ledger g = g_recognized + g_discarded as an **exact identity**
(law of total covariance). Does that shape survive quantum-ly?

**The measurement gap is a separate object from classical conditional variance.** 30 QTH tests; 50 tests with CID-1.

> QFI = CFI(measurement) + discard,  discard ≥ 0,
> vanishing **exactly** for an optimal measurement.

Here discard is defined as QFI minus CFI, so the displayed equality is
an identity. The substantive result is its nonnegative gap, with exact
saturation in the stated witnesses. Classical total covariance remains
valid inside outcome data. The declared commuting family and eigenbasis
measurement have zero measurement gap.

| Block | Certifies | Result |
|---|---|---|
| **T1** | Rational Bloch states; the SLD from an exact 4×4 rational solve, its defining equation verified as a matrix identity; QFI by **two independent routes** (solve vs closed Bloch form) agreeing in ℚ | PASS |
| **T2** | CFI ≤ QFI on the sampled rational measurement directions, **discard exactly 0** at the aligned direction and **exactly 1 (everything)** at the orthogonal one. The quantum ledger is an inequality with exact saturation | PASS |
| **T3** | **One tensor, two certified halves.** Q_ij = Tr(ρL_iL_j) splits exactly into a real symmetric SLD information matrix (diagonal = QFI) and an imaginary antisymmetric part (state-averaged SLD commutator divided by 2i, exactly −1/2 here, zero diagonal). The antisymmetric half is the **state-averaged commutator divided by 2i**. Noncommutativity is necessary but not sufficient for a nonzero average; rho=I/2 supplies an exact counterexample | PASS |
| **T4** | The classical bookkeeping inside the outcome data still holds exactly, while the quantum discard is nonzero by an exactly computed amount **no partition of outcomes can see**. The remainder is strictly quantum | PASS |
| **T5** | **The declared commuting family saturates the measurement bound**: SLDs commute, antisymmetric SLD pairing exactly 0, eigenbasis measurement optimal, discard exactly 0, and QFI equals CID-1's classical covariance route | PASS |
| **T6** | A rational unitary rotation preserves QFI **exactly**; the declared tangential depolarizing witness scales it by exactly λ². These witnesses preserve or lose information | PASS |

## Findings

- **QTH-F1.** Classical total covariance and the quantum measurement gap
  are different ledgers. QFI >= CFI is checked on the stated witnesses;
  defining discard as their difference also gives an exact equality.
- **QTH-F2.** One SLD Gram tensor has symmetric information and
  antisymmetric state-average parts. A connection-curvature identification
  needs an explicit adapter and normalization.
- **QTH-F3.** The selected commuting family with eigenbasis measurement has zero measurement gap.

**Companion, not duplicate:** CFE-Q certified the Bloch/holonomy side of
quantum CFE; this is the information side.

Witness capsule (one qubit, explicit rational families) like CFE-Q — not
a general open-system or quantum-estimation theorem. Kubo–Mori/Bogoliubov
metrics need logarithms and are **never evaluated**. **Not claimed:**
entropy production, work, heat, or that the antisymmetric part is a
physical Berry phase.

## Correction audit, 25 September 2026

The SLD tensor is a declared quantum comparison model, not a derivation of
quantum states or probabilities from EMK primitives. Its antisymmetric part
is a state average: zero does not imply commuting operators. At rho=I/2,
Lx=sigma_x and Ly=sigma_y have [Lx,Ly]=2i sigma_z but zero average. The new
regression and regenerated QTH1 certificate bind this control. The nonzero
witness at the original R0 remains valid. A connection, normalization, and
explicit adapter are needed to identify this pairing with a chosen curvature.

The native construction and cut-ledger comparison are in
[Physics from Cuts v0.8](../uncut-cut-measurement/NATIVE_RESPONSE_TENSOR_AND_CUT_LEDGER.md).
