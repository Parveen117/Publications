# EMK recognition geometry — the rotational seam metric (EMK-G1)

The geometry layer. EMK-1/EMK-2 identify the seam algebraically; this
capsule certifies the metric that gives it lengths, geodesics and
curvature.

## The construction

The local involution `J(z) = i·conj(z)` acts on real coordinates as
`(x,y) ↦ (y,x)`. In seam-adapted coordinates

```
u = (x+y)/√2   tangent to the seam
v = (y−x)/√2   normal to it        J : (u,v) ↦ (u,−v)
```

so seam reflection is ordinary reflection of the normal coordinate and
the fixed seam is the **45° line** `v = 0`. The reflection-symmetric seam
metric is

```
g_A = A(v)² du² + dv²      det g_A = A(v)²      K_A = −A″/A
```

## The move that keeps it exact

Carrying **W = A²** removes every square root:

```
Γᵘ_uv = W′/(2W)      Γᵛ_uu = −W′/2      K = (W′)²/(4W²) − W″/(2W)
```

For the quadratic warp `W = 1 + κv²` this gives, exactly,

```
K_κ(v) = −κ / (1 + κv²)²        K_κ(0) = −κ
```

Every quantity is a rational function of rational data — the whole seam
geometry is certified with **no square root and no float**. The `√2` in
the coordinate normalization is a positive constant that changes no line,
no fixed set and no curvature sign; the certified arithmetic uses
unnormalized seam coordinates and records the normalization explicitly.

## Certified blocks

| Block | Statement | Verdict |
|---|---|---|
| T1 | The involution acts as `(x,y)↦(y,x)`; `u` invariant, `v` negated, exactly on an integer grid; fixed set is precisely the line `y=x`; `J²=id` | PASS |
| T2 | `det g_A = W`; positive definite exactly where `W>0`, degenerate exactly at the zeros of `A`; for `κ<0` the locus is the exact condition `1+κv²=0` | PASS |
| T3 | Christoffels in W-form agree with Levi-Civita computed independently; curvature via `R_uvuv/det g` agrees with the closed form — two routes, 25 grid points; **K(0) = −κ exactly** | PASS |
| T4 | **The 45° seam is an intrinsic geodesic** of every reflection-symmetric warp (`W` even ⇒ `W′(0)=0`). Control: an odd perturbation breaks the symmetry, `W′(0)≠0`, and the seam is no longer geodesic | PASS |
| T5 | `κ>0` → globally definite, `K<0`; `κ=0` → Euclidean, `K≡0`; `κ<0` → `K>0` strictly inside `v² < −1/κ`, degenerate exactly at the boundary | PASS |
| T6 | **The two curvatures are independent** — see below | PASS |

## The guard (T6)

The source states as a *principle* that metric (Gaussian) curvature and
the declared connection curvature are distinct. This capsule certifies
that as a **separation** rather than asserting it:

| configuration | Gaussian curvature | recognition curvature |
|---|---|---|
| Euclidean seam metric (κ=0) | 0 | 2 |
| curved metric, rotation-free | −3 | 0 |

Neither determines the other. Any future identification must be **proved**,
never inherited from the fact that both are called "curvature". This is
the most tempting over-claim in the geometry layer, and it is now blocked
by a certificate.

## Findings

**F1** — the 45° seam is an intrinsic geodesic, not a visual bisector, with
the reflection symmetry certified load-bearing.

**F2** — the entire seam geometry is exactly rational under `W = A²`, so
`K(0) = −κ` needs no square root and no float.

**F3** — metric curvature and recognition curvature are certified
independent.

## Claim boundary

The metric family is **declared** (a reflection-symmetric warp), as in the
source; `K(0) = −κ` is a theorem for this family, **not** a universal law
equating flow parameters with curvature scalars. Global and completeness
statements are not claimed. No identification with the recognition or
thermodynamic two-form is claimed — T6 certifies the opposite. CFE's
(U) uniqueness stays open. No RH / K₀ / L₀ / YM continuum gate is
touched; quantum gravity is not touched.

## Reproduce

```
python papers/emk-recognition-geometry/certificates/emkg1_rotational_seam_metric.py
python -m pytest papers/emk-recognition-geometry/tests -v
```

CI regenerates the certificate and fails on pin drift
(`EXPECTED_EMKG1.sha256`).

---

# EMK-G2 — Global quotient, completeness, and seam holonomy

Source: vault `emk_ugd_recognition_geometry/sections/global_quotient_holonomy.tex`
(281 lines) — the global continuation of EMK-G1: same seam, now on the
periodic cylinder. 29 tests; 49 in this folder.

**The section's real content is a second guard**, and this capsule
certifies it in both directions: *trivial Levi-Civita holonomy does not
imply recognition closure.*

| Block | Certifies | Result |
|---|---|---|
| **T1** | Since A > 0, the geodesic criterion A′(0) = 0 ⟺ W′(0) = 0, so EMK-G1's local criterion transfers unchanged; seam length² = L²·W(0), exactly L² for the family. Control: an odd perturbation kills geodesy while the curve still **closes** — closure and geodesy are different properties | PASS |
| **T2** | Orbits of ȷ(u,v) = (u,−v) have size exactly 2 off the seam, 1 on it ⇒ the quotient is the half-cylinder with a **mirror** boundary; counting witness (double = 2× quotient off-seam); the side target sgn(v) factors through **no** map on the quotient | PASS |
| **T3** | The jet criterion in **W-form** (equivalent since A > 0): even warps glue smoothly to themselves; an odd term fails at n = 1 exactly yet glues to its own reflection; a mismatched pair fails at the identified jet. **The seam shift α is metric-invisible** while a declared ledger sees it | PASS |
| **T4** | κ > 0 ⇒ W ≥ 1 everywhere; κ = 0 flat; κ < 0 ⇒ positivity exactly on v² < −1/κ with degeneracy at **normal distance² exactly −1/κ** (finite ⇒ incomplete). Curvature blow-up as an **exact growth witness**, not a limit claim; the polar-cap condition fails | PASS |
| **T5** | Holonomy by **two independent exact routes** — boundary walk of ω¹₂ = A′du, and polynomial integration of K·A where K·A = −A″ is certified by **exact polynomial division** (remainder zero), not evaluation. Seam loop: Θ_LC = 0 **exactly**, not merely mod 2π | PASS |
| **T6** | **The guard, both directions:** (a) Θ_LC = 0 on the seam while recognition monodromy = −I; (b) Θ_LC ≠ 0 on a rectangle whose monodromy is exactly the identity. Neither determines the other; the ledger closes the residue; global closure fails when any active sector stays open | PASS |

**Guard thread.** EMK-G1 T6 separated *metric curvature* from
*recognition curvature*. EMK-G2 T6 separates *Levi-Civita holonomy*
from *recognition monodromy* — the same guard at the global level.

Declared: both warp families, the auxiliary EMK bundle, ρ_Σ,
tolerances. Identification of the two holonomies is **not claimed** —
T6 certifies the opposite.
