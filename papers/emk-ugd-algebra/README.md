# EMK / UGD algebra — the determinant seam ladder (EMK-1)

The framework's **own algebraic core**, brought from the vault into the
public certified corpus. Until now the Cut-First Equivalence capsules
used an explicit rational equation of state as their witness — a
scaffold. This is the theory's native geometry instead.

## The primitive algebra

```
K² = I      seam reflection
R² = −I     circular phase transport (the quarter turn)
RK = −KR    mixed rotation-cut parity,  (RK)² = I
[R, K] = 2RK                             the UGD bracket
```

A primitive block is `M = aI + bK + cR + dRK`.

## The determinant identity

> **det M = (a² − b²) + (c² − d²) = Δ∥ + Δ⊥**

The determinant splits **additively** into a seam-compatible channel and
a rotational / anti-seam channel, and the two close or fail
**independently**.

## Certified blocks

| Block | Statement | Verdict |
|---|---|---|
| T1 | Primitive relations K²=I, R²=−I, RK=−KR, (RK)²=I, and the UGD bracket [R,K]=2RK | PASS |
| T2 | The additive determinant identity on the full rational grid (2401 quadruples), two independent routes; **channel independence** exhibited — each channel can vanish alone, and both can vanish while the operator is nonzero | PASS |
| T3 | Seam projectors P±=(I±K)/2 idempotent/orthogonal/complete; `[A,K]=0` **iff** the mixing blocks vanish — the commutator *is* the seam-mixing detector | PASS |
| T4 | Exact Schur identity + the quadratic determinant correction; **control:** both diagonal blocks have det 1 yet the quadratic residue is 1/6 — first-order closure is provably incomplete | PASS |
| T5 | Cut-commutator curvature `C_iC_x − C_xC_i = −Δx·Δi·[D_x,D_i]` exactly; path equality ⟺ vanishing commutator channel | PASS |
| T6 | Winding index an exact integer, homotopy invariant (1, 1, 0 across scaled and off-origin loops) — **non-healable curvature**, via rational circle parametrization and exact crossing counts, no trigonometry | PASS |
| T7 | **Bridge to CFE:** a block is memoryless (commutes with K) **exactly** when rotation-free (c=d=0), and then det M = Δ∥ alone | PASS |

## Findings

**F1** — the additive split is *structure, not notation*: the channels are
genuinely independent.

**F2** — first-order closure is provably incomplete; diagonal-block tests
cannot certify recognition closure.

**F3** — the CFE memory dial and the EMK rotational determinant channel are
the same object in two presentations, so CFE's witness EOS can be
replaced by this native algebra.

## Provenance

Source: vault LaTeX — the EMK algebraic core ladder I–VII and the UGD
algebraic-operators appendix. **No executable version existed anywhere**
in the corpus (checked Recognition-Kernel-Framework and RH-Framework:
neither contains a UGD algebra implementation). This capsule is the first
machine-checkable realization.

## Claim boundary

Certified: the finite-dimensional theorem spine. **Not claimed:** the
infinite-dimensional trace-class, essential-spectrum, Hurwitz-zeta,
ξ-function or Hilbert–Pólya extensions — the vault source explicitly
defers those, and so does this capsule. CFE's (U) uniqueness stays open.
No RH / K₀ / L₀ / YM continuum gate is touched; quantum gravity is not
touched.

## Reproduce

```
python papers/emk-ugd-algebra/certificates/emk1_determinant_seam_ladder.py
python -m pytest papers/emk-ugd-algebra/tests -v
```

CI regenerates the certificate and fails on pin drift
(`EXPECTED_EMK1.sha256`).
