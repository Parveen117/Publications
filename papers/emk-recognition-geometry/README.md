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
