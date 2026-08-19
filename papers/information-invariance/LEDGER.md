# Ledger — Information Invariance paper

Status vocabulary: **PROVED** (finite theorem + executable certificate) ·
**OPEN** (obligation named, not discharged) · **CONJECTURE** (stated, no proof
attempted here) · **RETRACTED** (present in v1, removed in v2 with reason).

## PROVED
| id | statement | certificate |
|---|---|---|
| A | δ𝓘 = 0 ⟺ dω = 0 ⟺ ω exact, on a finite simply-connected chart; on the λ-chart this is the Onsager-flat sector Γ_cΓ_m = 1 | `certificates/l1_certificates.py` → `L1_certificate.json` |
| B | A single cycle's loop area does not witness curvature: explicit flat-with-offset and curved-without-offset constructions with equal first-cycle area, separated only at cycle 2 (reopen ratio 0 vs 1) | same |

Pin: `certificates/EXPECTED_L1.sha256`, regenerated in CI on every push.

## OPEN (named obligations, v2 §V)
| id | obligation |
|---|---|
| D1 | Derive the null-coordinate energy dynamics *from* Definition 1, rather than positing it; show why α_G ℰ² is the unique quadratic term under δ𝓘 = 0. |
| D2 | Exhibit the rescaling and coarse-graining that take the linearised dynamics to the Schrödinger equation, including the origin of the complex phase from a real energy density, and the identification of the correlation scale with ħ/m. |
| D3 | Derive L_P (which contains ħ) inside a dynamics declared deterministic and non-quantised, instead of inserting it as a regulator; then derive the emergent Einstein equations rather than citing the analogy. |
| D4 | Relate Γ_grav ~ Gm²/ħV to the Diósi–Penrose rate: either reduce to it, or state the regime where the two differ observably. |

## RETRACTED from v1 (with reason)
- "Validated Onsager violations" / "absence of defect-driven entropy generation" → the measured first-cycle loops **do not close** at 300 K; by Theorem B the observation is compatible with irreversible drift. Replaced by the honest statement plus the predeclared cycle-2 rule (T02).
- "Reproducibility across material platforms" → one platform (CVD graphene on SiO₂), two spots, one cycle. Rewritten.
- Numerical decoherence example (10⁻³ s⁻¹ for a 10⁻¹⁷ kg, 10⁻²¹ m³ resonator) → arithmetic not reproduced; removed pending D4.
- "Triple invariance" section incl. the Riemann-Hypothesis remark → removed; no obligation for it is discharged anywhere, and RH is OPEN in the framework's own ledger.

## Patent boundary
The theory content here is covered by, and later than, PCT/IB2025/060887
(filed 26 Oct 2025, priority 6 Aug 2025). Publication of the theory does not
affect the device/method claims. No unfiled device detail appears in this paper.
