# Ledger — Information Invariance paper

Status vocabulary: **PROVED** (finite theorem + executable certificate) ·
**OPEN** (obligation named, not discharged) · **CONJECTURE** (stated, no proof
attempted here) · **RETRACTED** (present in v1, removed in v2 with reason).

## PROVED
| id | statement | certificate |
|---|---|---|
| A | δ𝓘 = 0 ⟺ dω = 0 ⟺ ω exact, on a finite simply-connected chart; on the λ-chart this is the Onsager-flat sector Γ_cΓ_m = 1 | `certificates/l1_certificates.py` → `L1_certificate.json` |
| B | A single cycle's loop area does not witness curvature: explicit flat-with-offset and curved-without-offset constructions with equal first-cycle area, separated only at cycle 2 (reopen ratio 0 vs 1) | same |
| C | A real linear flow x' = Ax is exactly a Schrödinger flow iff it carries a complex structure J with A antisymmetric and AJ = JA; then H = −ħJA is Hermitian and iħψ' = Hψ | `certificates/d2_phase_origin.py` → `D2_certificate.json` |
| E | On a balanced cut, a real generator is a Schrödinger generator iff its diagonal (cut-even) blocks are equal and antisymmetric and its transport (cut-odd) block is symmetric with C = −B; phase lives in the odd transporting part, damping in the symmetric part | `certificates/d2a_cut_supplies_structure.py` → `D2a_certificate.json` |
| D | One-channel obstruction: an odd-dimensional real flow — in particular a single scalar channel — admits no complex structure (det J² = (−1)^N), so no rescaling of the scalar linearisation can yield Schrödinger evolution | same |

Pins: `EXPECTED_L1.sha256`, `EXPECTED_D2.sha256`, `EXPECTED_D2A.sha256` (in `certificates/`) — regenerated in CI on every push.

## OPEN (named obligations, v2 §V)
| id | obligation |
|---|---|
| D1 | Derive the null-coordinate energy dynamics *from* Definition 1, rather than positing it; show why α_G ℰ² is the unique quadratic term under δ𝓘 = 0. |
| D2 | **Reduced by Theorems C, D** (not discharged). The scalar linearisation d(δℰ)/du = 2α_G Ē δℰ + … is one real channel and therefore cannot be Schrödinger evolution. What remains: |
| ~~D2a~~ | **DISCHARGED by Theorem E** — the conjugate channel is the image of the fluctuation on the *opposite sheet of the oriented cut*, and the complex structure is the framework's derived quarter turn ι (F00E: ι²=−1, ι†=−ι; grading from theorum/41, pin 34afc445…). Nothing new is assumed: both are already PROVED/certified in Recognition-Kernel-Framework. |
| ~~D2b~~ | **EXECUTED — negative result** (`certificates/d2b_linearisation_verdict.py`). On a transverse lattice with the parity cut: the quadratic term contributes pointwise multiplication (symmetric, cut-even → pure growth/decay); a diffusive term is symmetric (→ diffusion); only a first-difference transport term is antisymmetric and cut-odd (→ phase). For the full generator the symmetric part is exactly M + L, which does not vanish at a nonzero background, and rescaling cannot remove it. **Therefore the linearised dynamics as written is not a Schrödinger generator; it is a damped or amplified rotation.** Unitarity requires M + L = 0 at the background — a concrete modelling requirement on the dynamics, now stated rather than assumed. |
| ~~D3~~ | **EXECUTED — negative result** (`certificates/d3_regulator_scale.py`). Theorem F: ħ is dimensionally independent of {c, G}; no length can be built from c and G alone, so L_P = √(ħG/c³) requires ħ as an *input* — inserting it presupposes the constant the framework claims to derive. Theorem G: the length the dynamics does generate is L_E = c²/√(G ℰ) ≡ Λ^(−1/2), of order the Hubble scale (~10⁶¹ L_P). **Hence the model's own regularisation is cosmological, not Planckian, and the ultraviolet-finiteness claim does not follow.** The emergent-Einstein half of the conjecture remains untouched and OPEN. |
| ~~D4~~ | **EXECUTED — corrected** (`certificates/d4_decoherence_audit.py`). Theorem H: Gm²/(ħV) has dimensions M⁰L⁻²T⁻¹ — a rate *per unit area*, not a rate; the unique length-based repair Gm²/(ħR) is exactly the Diósi–Penrose rate E_G/ħ. Numerically the draft's 10⁻³ s⁻¹ is recovered (6.3×10⁻⁴ s⁻¹) **only** after the repair, i.e. it is a DP number. No regime of departure from DP has been exhibited, so this must be stated as a consistency requirement, not an independent prediction. Three concrete ways it could become one are listed in the certificate. |

## DISCHARGED WITH NEGATIVE RESULT
- **D3 (Planck-scale regularisation).** Proved that the regulator was inserted, not derived, and identified the scale the theory actually produces (cosmological). The UV-finiteness claim of v1 is withdrawn.
- **D2 (emergent Schrödinger step).** Theorems C, D, E fixed the criterion; the computation (D2b) shows the proposed energy dynamics fails it. The emergent-quantum claim of v1 is therefore withdrawn as a result and retained only as a *conditional* statement: it holds if and only if the dynamics is amended so the cut-even symmetric part vanishes at the background. Theorems A and B are independent and unaffected.

- **D4 (decoherence prediction).** Formula corrected (dimensional error), identified with Diósi–Penrose, and demoted from prediction to consistency check.

- **D1 (origin of the dynamics).** Negative half: invariance alone derives nothing. Positive half: the quadratic law follows from a *stated* selection principle (one new length, used once), which also repairs the coupling's dimensions.

## RETRACTED from v1 (with reason)
- "Validated Onsager violations" / "absence of defect-driven entropy generation" → the measured first-cycle loops **do not close** at 300 K; by Theorem B the observation is compatible with irreversible drift. Replaced by the honest statement plus the predeclared cycle-2 rule (T02).
- "Reproducibility across material platforms" → one platform (CVD graphene on SiO₂), two spots, one cycle. Rewritten.
- Numerical decoherence example → the printed formula is dimensionally invalid (Theorem H); the quoted value is recovered only under the Diósi–Penrose repair, and is presented as such.
- "the fluctuation dynamics yield the Schrödinger equation" → shown false for the stated dynamics by D2b; replaced with the conditional statement and its exact condition.
- "Triple invariance" section incl. the Riemann-Hypothesis remark → removed; no obligation for it is discharged anywhere, and RH is OPEN in the framework's own ledger.

## Patent boundary
The theory content here is covered by, and later than, PCT/IB2025/060887
(filed 26 Oct 2025, priority 6 Aug 2025). Publication of the theory does not
affect the device/method claims. No unfiled device detail appears in this paper.
