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

---

# EMK-G3 — The helical sheet-memory lift: flat is not memoryless

Source: vault `emk_ugd_recognition_geometry/sections/helical_sheet_memory.tex`
(245 lines) — the door EMK-G2's closing remark left open. 27 tests;
76 in this folder.

**The sharpest principle in the geometry layer**, certified as an exact
separation: *a flat lifted connection may carry nontrivial global
monodromy.* Local curvature detects infinitesimal nonclosure; monodromy
detects global return data.

| Block | Certifies | Result |
|---|---|---|
| **T1** | The helix is the **orbit structure of a deck transformation**, not a drawn curve: ρ is a bijection with exact inverse, monodromy after n circuits is exactly ρⁿ (n positive, zero, negative), the action is **free** and additive, and every orbit meets 0 ≤ u < L exactly once | PASS |
| **T2** | The compensated variables φ − p_φu, σ − p_σu are **exactly invariant** under the deck map; the uncompensated ones fail by exactly α and β; a wrong pitch fails by an exactly computed amount | PASS |
| **T3** | **The principle.** A = p du has curvature **exactly zero** as a bivariate polynomial identity and **exactly zero** holonomy on every contractible rectangle — yet one circuit of the base circle gives **exactly α ≠ 0**, and n circuits exactly nα. A curved companion (curvature exactly −1, holonomy = exact flux) completes the separation in both directions | PASS |
| **T4** | Visible return is exactly α ≡ 0 (mod K); full lifted return also needs β = 0 and q = 0 — exhaustive on a grid, with lifted ⟹ visible and never the converse. After n circuits the displacement is exactly (nα, nβ, nq): **a sheet increment is never undone by more circuits** — sheet memory is not a phase that eventually wraps | PASS |
| **T5** | The four return classes with a witness each, mutually exclusive. **Memory-bearing return is distinct from exact return**: remove the ledger and the same state becomes an open obstruction; de-admit the sheet sector and it reads as ordinary lawful transport | PASS |
| **T6** | Consuming **pinned EMK-G2**: base return holds, Θ_LC = 0 exactly, RTC curvature zero, lift exactly flat — **and the transition is still open**. Strong form: for **all 31 proper subsets** of the five sectors, a residue vector closes exactly that subset while the whole stays open. A **non-faithful projection** dropping the sheet sector maps an open lift to a closed image — the converse of the projection theorem fails exactly | PASS |

**Guard thread, three levels.** EMK-G1 T6: metric curvature vs
recognition curvature. EMK-G2 T6: Levi-Civita holonomy vs recognition
monodromy. EMK-G3 T3: local flatness vs global monodromy. Same guard,
each time one level higher.

Declared: the fibre model, ρ, the pitches, A^EMK, ledgers, tolerances,
the thermodynamic presentation map. Not certified: the non-abelian
lift, the cut-localized sheet term. **Not claimed:** that α, β, q are
phase lag, scale drift or hysteretic branch count.

---

# UGD-G1 — Phase–scale–seam geometry: the number layer meets the surface

Source: vault `emk_ugd_recognition_geometry/sections/ugd_phase_scale_seam_geometry.tex`
(213 lines). **This closes a loop.** EMK-G3's fibre ℝ_φ × ℝ_σ × ℤ_k is
not an arbitrary choice — it *is* the UGD state, whose numeral
presentation UGD-1 certified independently. 28 tests; 104 in this folder.

| Block | Certifies | Result |
|---|---|---|
| **T1** | Three sectors, three topologies, **three different closure rules** (phase mod K, scale by equality in ℚ, sheet by equality in ℤ). No collapsing audit works: a mod-K audit **passes** an open sheet, an exact-equality audit **rejects** a genuine phase return. The residue vector is irreducibly a vector | PASS |
| **T2** | **The cocycle theorem, sharpened.** Under any chart re-choice n→n+m_a−m_b the triple sum is **exactly invariant** — certified exhaustively *and* by the identity making the shifts cancel in pairs — so a nonzero integer sheet defect is a genuine class that must enter the ledger and **cannot be erased by chart choice**. Scale: exact composed-identity condition with the closing chart constructed | PASS |
| **T3** | Abelian curvature componentwise; flat-with-period **consumed** from pinned EMK-G3, not repeated. New: the sheet defect isn't a curvature component at all — its invariance involves no continuous data, so no relaxation (flat or curved) changes it. F_UGD = 0 ⇏ Δk = 0 for a **structural** reason | PASS |
| **T4** | The two named non-implications with witnesses, strengthened to **pairwise independence** across all 12 ordered pairs of the four sectors | PASS |
| **T5** | Faithful ⟹ preserves **and** reflects closure; lossy preserves but provably fails to reflect. Third mode separated: a projection that **identifies** two sectors discards nothing yet still fails — phase +3 with sheet −3 merges to zero and reads closed. **Injectivity, not non-discarding, is load-bearing** | PASS |
| **T6** | **The bridge.** Consuming pinned UGD-1: numerals with identical classical projection and phase content but different seam charge give a transition whose phase and scale residues are **exactly zero** and whose sheet residue is **exactly the seam-charge difference**. UGD-1 T4 ("classical projection blind to seam") and this section's "lossy projection reports closure" are **one statement at two layers** | PASS |

The number layer and the geometry layer are the same object seen twice.

Declared: state space, transition, ledgers, tolerances, chart cover,
P_UGD, and the thermodynamic dictionary (which the section itself calls
a constitutive presentation, not a universal identity — its
σ = log λ_s option is **never evaluated**). Phase closure certified in
exponent-mod-K form throughout.
