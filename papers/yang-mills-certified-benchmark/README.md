# Yang–Mills Certified Benchmark (YM-1)

First rung of resuming the MP Yang–Mills adapter with the certified RNKE
machinery: the reduced finite-cutoff gap of the SU(2) one-holonomy Wilson
benchmark, previously a float, is now an exact two-sided rational enclosure.

## Certified statement

For the normalized SU(2) Wilson convolution operator on class functions
(source: `MP/adapters/yang_mills/cutoff_gap_benchmark.tex`), with exact
character spectrum `lambda_j(beta) = I_{2j+1}(beta)/I_1(beta)`:

```text
Delta_red(a=1, beta=2) = -log( I_2(2) / I_1(2) )
  in [0.83672330623158891305006780105454768321... ]  (width < 1e-30, exact ℚ)
```

Verdict engine: directed rational interval arithmetic only. No floats in any
verdict. Controls: (C1) `lambda_1` strictly inside (0,1) hence gap strictly
positive; (C2) tampering `I_2 -> I_3` produces a separated bracket; (C3)
width budget.

## Claim boundary (fail-closed)

- CERTIFIED: the reduced finite-cutoff benchmark number only.
- NOT CERTIFIED / OPEN: full-lattice gap, theta-graph interacting transfer,
  ultraviolet/infinite-volume uniformity, OS continuum reconstruction, the
  Clay existence and mass-gap predicate. Adapter verdict remains `hold`.

## Ledger

| Item | Status |
|---|---|
| YM-1 reduced gap enclosure (a=1, beta=2) | PASS (pinned) |
| YM-2 interacting theta-graph gap, certified small coupling (kappa < Delta_red/6) | PASS (pinned) |
| YM-3 first-order crossing direction: rank-one, A<->B transported, slope lambda_half/4 | PASS (pinned) |
| YM-4 symmetry-protection theorem [S,T_kappa]=0 + finite-kappa certified lower bounds | PASS (pinned) |
| YM-5 certified TWO-SIDED gap at finite kappa, beyond the sandwich threshold (Kill-Lemma realization) | PASS (pinned) |
| YM-6 seam-integer dock: EXACT eigenvalue counts on the native 5x5, certified to kappa=1/2 (D03 congruence) | PASS (pinned) |
| YM-7 V7 carrier: kappa=7/10 unlocked (exact count), certified 7-eigenvalue crossing curves | PASS (pinned) |
| YM-8 CAPSTONE: theta-graph gap is a THEOREM at every coupling (certified kernel floors + pinned Jentzsch anchor); remaining Millennium content = uniformity in the cutoff | PASS (pinned) |
| YM-9 FIRST UNIFORMITY: exact heat-kernel refinement family; free gap exactly cutoff-independent; interacting uniform bound 3/8 along a declared trajectory | PASS (pinned) |
| YM-10 blindness ledger (SPECTRAL-1/2 framework): exact multiplicity law, exact probe counts, composition ledger, YM-9 T4 correction | PASS (pinned) |
| YM-11 three gate verdicts on the toy carrier: gauge CLOSED; volume/IR SPLIT (free closed, sandwich route proven insufficient, critical volume exact); universality SPLIT (counting closed, metric open) | PASS (pinned) |
| YM-12 the T03 move: positivity SQUARE-SOURCED (classical anchor demoted to simplicity-only); volume+cutoff-uniform interacting gap 3/8 on the factorized chain; governance capsule | PASS (pinned) |
| YM-13 verification-overlap beta audit: beta 4/5 -> 3/5 via independent witnesses (CF Bessel, compound-interest exp, det/trace inertia); three lineage overlaps retained with refusal notes; RH remains lowest at 1/5 | PASS (pinned) |
| YM-14 interaction-overlap dock, first carrier: bridged pair — overlap seam rank-one and swap-transported (same theta integral), certified k=1 to kappa=1/2, overlap price kappa*lam/4 on one line; kappa=1 carrier limit recorded | PASS (pinned) |
| YM-15 1D block-transfer dock opened: exact closed form of the chain compression for EVERY m (vacuum f0^{m-1}, A-line f0^{m-1}·KMS_m(r), r=I2/I1; SU(2) convolution lemma, machine-checked as a formal polynomial identity m=2..8, no truncation remainder); m-UNIFORM ratio bracket λ(1−r)/(1+r)<ρ_m<λ(1+r)/(1−r), certified positive normalised gap for all m at κ≤1; κ=2 fails closed. A-line carrier ONLY — complement not bounded | PASS (pinned) |
| YM-16 chain dock: exact B-factorisation ⇒ gap(S_m) ≤ Δ_red for ALL m,κ (m-uniform upper bound, exact eigenvector); vacuum volume bracket f0^{m−1} ≤ λ₁ ≤ e^{κ(m−1)}; certified two-sided gap (ratio ≤ 9/10) of S_m for m ≤ m*(κ) — m*=13 (κ=1/8), 7 (1/4), 4 (1/2) — with m* computed exactly and matched by the dock; m*+1 refused. Sup route provably non-uniform: death volume has coordinates | PASS (pinned) |
| YM-17 interleaving seam: m_κ = B_odd·B_even; each dressed half-chain is a tensor product of the two-site chain at 2κ ⇒ product spectrum ⇒ gap ratio EXACTLY m-uniform (pair ratio certified ≤ 9/10 at κ≤1/2; bracket [1,3] honest at κ=1); vacuum bracket narrowed to rate κ² (per-bridge upper end 0.125→0.030 at κ=1/8); whole volume problem reduced to ONE inequality: λ₁(T) vs σ₁(X)σ₁(Y) (vacuum tracking of interleaved halves) | PASS (pinned) |
| YM-18 CALIBRATION + DESIGN: enlarged exact carrier W'_m (vacuum + sites + bridge-pairs, dim 2m) gives certified Rayleigh lower bounds on λ₁(T_A,m), m=2..6; per-bridge vacuum rate φ_m^lo is m-independent to ~1e-5 (0.002177 @κ=1/8, 0.008698 @1/4, 0.03461 @1/2) — vacuum tracking behaves as a clean per-bridge rate; bracket width = looseness of the upper end. Cluster dock DESIGNED (monomers = odd pairs at 2κ, activities = even bridges, KP on a line) — no claim | PASS (calibration, pinned) |
| YM-19 Dobrushin dock — **DEMOTED to ANCHORED** (governance rule: no classical borrowing; T3 rests on a cited Dobrushin/Föllmer anchor, same fate as YM-8 Jentzsch). T1 (coefficient bound) and T2 (α<1 arithmetic) are native and stand; the m-uniform GAP claim is WITHDRAWN until derived natively. Was: on the bounded-overlap chain — Dobrushin coefficient C_ij ≤ 1−1/δ² (elementary, proved), δ_s = e^{2κ} exact, δ_t ≤ (1+S_a)/(1−S_a) certified; α = 2(1−1/δ_t²)+2(1−1/δ_s²) < 1 ⇒ (pinned anchor: Dobrushin/Föllmer/Künsch covariance decay; extraction lemma proved) λ₂/λ₁ ≤ α for EVERY m. Certified gap/step ≥ 0.263 (a=6,κ=1/16), ≥ 0.145 (a=8,κ=1/8, ON YM-9's trajectory θ=1/64), ≥ 0.235 (a=12,κ=1/8); refused at (6,1/8),(4,1/16),(10,1/4). Coupling ceiling of the route exact: κ < log2/4. STANDING CORRECTION to YM-17 T3 (vacuum tracking insufficient; correlation-decay route replaces it). NOT uniform in the cutoff (δ_t→∞ as a→0) | ANCHORED (pinned, not a native theorem) |
| YM-20 NATIVE ORIGIN AUDIT: machine-checked origin ledger for YM-1..19 in the framework's vocabulary (carrier = COORDINATE_SHADOW throughout; counting = NATIVE_DERIVED throughout; anchors: YM-8 and YM-19 demoted imports, YM-17 SV-Weyl as technique). CAYLEY FORM certified: YM-15 ceiling (1+r)/(1−r) and YM-19 δ_t are Exp_Σ(A_Σ(y)) (F00-G Thm 5.1 pinned) with seam coordinates r, S_a; native odd-series log agrees with YM-1 route to 1e-20; uniform gap in native form = −Log_Σ(λ) − A_Σ(r). Bessel coefficients = native factorial series (F00-E Lemma 2.1 shape); Haar identification = the shadow. NAMED YM-21: bridge-local recognition-energy contraction under the seam-involution flow (F00-E Thm 6.2) | PASS (pinned) |
| **YM-F1 THE CHAIN AS A RECOGNITION FABRIC (native carrier)**: T1 SU(2) = unit sphere of the EMK block over C_Σ — span{I, R, ιK, ιRK} is Hamilton's quaternion algebra, det = Δ∥+Δ⊥ becomes the quaternion norm (twisted seam channels), χ½ = 2×identity-sector coefficient; T2 chain = ladder fabric, bridge faces = plaquette holonomies H_l = A_i⁻¹A_{i+1}; T3 Recognition-Stokes EXACT on the chain (interior rungs cancel, E_int = 0, orientation tamper bites); T4 action = sum of face residues ρ_l = 1−Tr H_l/2 ∈ [0,2], weight = Exp_Σ(κ(m−1))·Exp_Σ(−κΣρ). Source: vault fabric appendix + EMK-1 + F00-E. Remaining shadow: Haar ↔ Φ_Σ (declared). Native form of YM-21 posed on the fabric | PASS (pinned) |
| YM-21 THE TILING LAW (native, on the fabric): STANDING CORRECTION to YM-F1's posed question (boundary residue cannot control face residue — witness g, g⁻¹, certified); the right object is the tiling. T2 Stokes sub-telescoping ⇒ YM-15 entry = r^{#faces crossed}: SPATIAL decay along the chain exactly geometric, rate −Log r, every m. T3 leading-order TIME two-point = (λ·KMS_m(r))^t, bracketed λ^t((1∓r)/(1±r))^t for every m ⇒ YM-15/20's uniform rate −Log λ − A_Σ(r) DERIVED as a tiling count (t=1..6, m=2..8 inside bracket). T4 remainder named natively: higher-content/branching tilings, expansion parameter = face-coefficient ladder r_j (certified decreasing: r_½ ≫ r_1 ≫ r_{3/2}); uniform convergence = native cluster statement | PASS (pinned) |
| **YM-22 THE TILING RATE SURVIVES THE CUTOFF**: on the heat-kernel family along κ=θa, the time-face log is EXACTLY 3a/4 (no transcendental), so the leading tiling rate per unit time γ(a) = 3/4 − A_Σ(r(θa))/a. Certified over SIX DECADES a=1…1e-6 at θ=1/16: **23/32 ≤ γ(a) ≤ 3/4 for every a** — the limit 3/4 − θ/2 is a uniform LOWER bound in the cutoff (approached from above); KMS bracket holds with λ_a for m=2..8 at every a ⇒ leading-order rate uniform in volume AND cutoff simultaneously. Beats YM-9's sandwich 3/4 − 6θ by factor 12 in the θ-coefficient; trajectory ceiling θ < 3/2 (θ=2 refused). Higher face coefficients per unit time r_1/a, r_{3/2}/a → 0: in the cutoff limit only j=½ branching survives. SCOPE: leading tiling order | PASS (pinned) |
| **YM-23 THE WEAK-COUPLING TURN**: T1 tiling route has an EXACT weak-coupling boundary — ladder collapses (r_½,r_1,r_{3/2} = 0.970, 0.922, 0.859 at κ=50: no small parameter), leading rate sign change at r_c(a) = tanh(3a/8) in native Cayley form, κ_c(a) bracketed (1.574 @a=1, 0.759 @1/2, 0.376 @1/4, 0.188 @1/8), and on a DECLARED AF-shaped trajectory κ = ¼log(1/a) the route is left at a = 1/8 — the continuum lies outside the tiling route. T2 native linearization: ρ = 2s/(1+s) in the stereographic chart; odd(H₁H₂) = p₁+p₂+½[p₁,p₂] EXACTLY — the commutator is the only non-additive term = EMK-1 T5 / RST-1 T4 curvature on the fabric. T3 EXACT NEGATIVE: the abelianized fabric's softest spatial stiffness ≤ 12κ/(m(m+1)) (Rayleigh identity certified m=2..40) — gapless in volume. T4 native statement: a weak-coupling gap lives entirely in the commutator sector (non-healable curvature) — the mass gap is non-perturbative, said in the framework's words | PASS (pinned) |
| YM-24 ABELIAN SUBFABRIC + NON-ABELIAN RESIDUE (exact algebra, no gap claim): T1 commuting faces = F00-E's circular Euler orbit, compose by the native addition law, Stokes additive in the Euler coordinate (the rotor chain), m=2..8 on rational circle points; T2 non-abelian residue = Gram defect |p×q|² = |p|²|q|² − ⟨p,q⟩² (zero iff one subfabric); T3 exact second-order Stokes odd₂ = Σp_i + Σ_{i<j} p_i×p_j, even₂ = 1 − Σ⟨p_i,p_j⟩, remainder degree ≥ 3 by exact interpolation; T4 YM-23's soft mode lies in an EXACT abelian subfabric, E_na = |Σ_{i<j}p_i×p_j|² vanishes there, positive generically, rotation-invariant; order tamper = 2 p_i×p_j (EMK-T2: order is content). Consequence: a weak-coupling gap must come from transitions BETWEEN abelian subfabrics (direction changes) | PASS (pinned) |
| YM-25 THE DIRECTION-CHANGE SECTOR: two exact negatives — the time kernel is direction-blind (bi-invariant: one scalar per spin-j block, d_j² coefficients, refinement-invariant) and the strong-coupling carriers are direction-blind (class functions, rotation-invariant: structural SPECTRAL-1 blindness) — and one exact positive: direction change between adjacent faces IS the 6j recoupling (MP gold/01 seed), squared overlaps between (12)3 and 1(23) schemes exactly [[1/4,3/4],[3/4,1/4]], unitary, no √ evaluated (RST-2 discipline), singlet tamper breaks unitarity. REDUCTION: a weak-coupling gap can live only in the intertwiner transfer | PASS (pinned) |
| **YM-26 DOCK TO theorum/28 (Recognition-Complete Finite-to-Infinite Cut Theorem)** — owner's correction: the framework's own convergence machinery. Dictionary: P_n = A-line carrier (YM-15), Q_n = complement, floor λ₁ ≥ f₀^{m−1}, threshold μ = νf₀^{m−1}, e_n = ‖QSQ‖/μ, N₊(B_n−I) = Haynsworth inertia; outward certificate e_n < 1 ∧ inertia count = 1 ⇒ λ₂ < νλ₁. T1: LAWFUL (e_n present, not a shadow certificate) for every m ≤ m*(κ) = 13, 7, 4 at κ = 1/8, 1/4, 1/2. T2: e_n(m) certified increasing, crosses 1 exactly at m*+1 — theorum/28's hypothesis 3 (memory-channel Cauchy bound uniform in m) is NOT YET BUILT; full hypothesis ledger 1–7 written in the theorem's vocabulary. T3: via the tiling law, item 3 reduces to a Cauchy bound uniform in m for the multi-insertion content-½ memory channel (the j=½ branching sector) | PASS (pinned) |
| **YM-27 DOCK TO RH-FRAMEWORK T01** (owner: use the RH journey's proved results, cite them): T1 the sup route IS T01-E4C (bounded native simple multiplier, U_Σ(w) = e^κ exactly; YM-16's e_n reproduced to the last digit) ⇒ theorum/28 hypothesis 3 = a memory-channel bound beating T01-E4C by a factor uniform in m — measured slack per face 56.5 / 27.8 / 13.5 at κ = 1/8, 1/4, 1/2 (YM-18 calibration vs E4C price). T2 the abelian subfabric is a T01-E5A rational UGD packet: wrap law = rotor-chain composition, winding = seam memory μ; (g,g⁻¹) witness vs full-turn fabric have the same visible boundary but μ = 0 vs 1 — the native form of EMK-1 winding. T3 the Haar shadow of YM-F1 is exactly RH's OPEN T01-E5C/E6 — YM inherits RH's item, invents none; aliasing structure of class functions on E1 grids certified (exact for n > 2j). T4 T01-B/C recognition energy = chain vacuum form, floor f₀^{m−1} | PASS (pinned) |
| **YM-28 T01-E4D ON THE CHAIN** (RH-Framework PR #19 applied): T1 the mean law is EXACT with zero deviation on single insertions — I_rec(W·χ½(A_i)) = f₀(2κ)^{m−1} for every m (formal identity, χ₁ branch dies at the open end); T2 two insertions: deviation = r₁(2κ)^{|i−j|}, a geometric correlator; T3 |S| ≤ 4, m ≤ 7: deviation ≤ (1+2r₁)^{|S|−1} — PER INSERTION, not per face (m does not enter); T4 E4D per-face price log(√f₀(2κ)/f₀(κ)) = 0.001949 / 0.007742 / 0.03016 accounts for 90 / 89 / 87 % of YM-18's measured true rate — the slack of YM-27 is now DERIVED (E4C/E4D = 63 / 31 / 16); T5 conditional projection m*_E4D ≥ 400 / 200 / 50 vs sup-route 13 / 7 / 4 IF E4D-C held on superpositions | PASS (pinned) |
| YM-29 CLOSED-FORM VACUUM FLOOR FOR EVERY m: native face-independence lemma (bridge faces independent contents, W a product over faces ⇒ face-product integrals factorise; fusion-rule certified: ∫wχ½ = 2f_½, ∫wχ½² = f₀+3f₁); T0^{1/2}χ½(H_ℓ) = λχ½(H_ℓ); ansatz ψ = 1 + cΣχ½(H_ℓ) gives λ₁ ≥ f₀^{m−1}·Q_m(c*) in closed form for ALL m, Q_m > 1 always, ~4λ²r²(m−1) for m ≫ 1/(4λ²r²) (Q₂₀₀₀ = 2.58 / 6.99 / 24.1). Honest: linear floor moves the dock's m* by 0 — the floor at the TRUE rate (0.000228/face margin, YM-28) needs the product ansatz, whose T-expectation is the one-time-step ladder = 6j recoupling contraction | PASS (pinned) |
| **YM-30 THE 2-ROW RECOUPLING ENGINE — BUILT AND VALIDATED**: exact surd arithmetic (Q adjoined √), Clebsch–Gordan by Racah's closed form (no floats), 3j, three-matrix-coefficient vertex integral, ε-removal of conjugates, ladder (one-time-step fabric) evaluated as a transfer over cut indices with label sums folded in. V1 every evaluation rational (surds cancel — the built-in consistency check); V2 one plaquette = 1/d_a²; V3 single-row limit recovers the YM-15 chain; V4 two plaquettes: Σ_c d_c N = 1/4 and the split d₀N₀ : d₁N₁ = 1/4 : 3/4 = YM-25's recoupling squares, now from the engine; V5 rail-swap and reversal symmetry. FIRST USE (calibration, κ=1/8, content ≤ 1, point coefficients): the W^{1/2}-dressed vacuum for T'' = W^{1/2}T0W^{1/2} beats the old floor f₀^{m−1} by a GEOMETRIC factor 1.000732 per face (m=2..6) — floor rate 0.002685/face > YM-18's bound 0.002177 > E4D excitation rate 0.001949: margin 0.000736/face; single-face excitation (c ≠ 0) does not help | PASS (pinned) |
| **YM-31 OUTWARD m-UNIFORM DRESSED VACUUM FLOOR — theorum/28 ITEM 6 AT THE TRUE RATE**: trial φ = W^{1/2}·(W_pt/W) with W_pt a product of rational truncated faces (no error by definition); numerator = ladder ⟨W_pt, T0 W_pt⟩ (YM-30 engine), rung truncation is a rigorous LOWER bound (T0 = Π_i Σ_c λ_c Π_c: a sum of nonnegative projection terms, monotone in λ_c); denominator ∫W_pt²/W = [Σ d_jd_kd_l f_j^pt f_k^pt (−1)^{2l} f_l N_jkl]^{m−1} exact with outward l-tail; NATIVE LOG-CONVEXITY Z_m² ≤ Z_{m−1}Z_{m+1} (Z_m = ⟨g, S^{m−2}g⟩, S = C^{1/2}M_K C^{1/2} symmetric PSD, cut-square inequality) ⇒ ratio nondecreasing ⇒ **λ₁(S_m) ≥ FLOOR_{m₀}·(ρ_lo/I_hi)^{m−m₀} for every m ≥ 6**, rung tail at m₀ bounded (5.9e-9). Per-face floor rate **0.002685 / 0.010722 / 0.042607** at κ = 1/8, 1/4, 1/2 vs old floor 0.001952 / 0.007802 / 0.031089 and E4D excitation rate 0.001949 / 0.007742 / 0.030161 — margin **0.000736 / 0.002980 / 0.012446 per face**, uniformly in m | PASS (pinned) |
| YM-32 DRESSED SINGLE-INSERTION ENERGY, m-UNIFORM (engine extended to rung insertions via CG reduction of 4-valent vertices; validated: j=0 reproduces YM-30, one-rail insertion vanishes, reversal identity, mixed-rung tamper smaller): E_m(p) = ⟨W_ptχ½(A_p), T0 W_ptχ½(A_p)⟩/⟨W_pt,T0W_pt⟩ exact for m=2..7, all p — bulk value m-independent to < 1e-6, ends differ; E_bulk = 0.43368 / 0.43532 / 0.44171 at κ = 1/8, 1/4, 1/2, i.e. λ·(1 + 1.27e-3 / 5.05e-3 / 1.98e-2), inside YM-15's bracket (I first wrote 'below λ'; the engine corrected me). Dressed denominator factorises (rung insertion free) ⇒ the J=½ sector's dressed excitation/vacuum ratio is exactly E_bulk, uniformly in m — a certified LOWER bound on that sector's top, not an upper bound on λ₂ | PASS (pinned) |
| YM-33 t-ROW FABRIC ENGINE + EXACT TIME TWO-POINT FUNCTION (owner's framing: time carries no infinity, only the flow λ^τ = Exp_Σ(−τ·gap); expand the space coupling only): general cut transfer over t rails with rungs contracted up the column (interior vertices 5-leg: two rungs, two faces, insertion — two CG fusions); validated: t=2 reproduces YM-32 exactly, t=1 the chain, 3-row fabrics rational with exact time reflection. Dressed sequence a_τ = ⟨φ, T''_pt^τ φ⟩ and vacuum z_τ exact for τ=1,2, m=5..8: ratio sequence ρ₁ = E_bulk (YM-32), ρ₂ ≥ ρ₁, and the SECOND-ORDER connected correction δ₂ = ρ₂ − ρ₁ is a BULK CONSTANT — 0.00076998 / 0.00305278 / 0.01178845 at κ = 1/8, 1/4, 1/2 — with boundary influence decaying geometrically: moving an end from distance 2→3 shifts δ₂ by 1.3e-9 (either end, identically), 3→4 by 1.6e-12, ratio certified in [r²/2, 2r²]. The framework's own strong-coupling series for the gap, order by order as exact rationals, m-uniform with exponentially localised edges. Lower-bound side (Rayleigh data) — item 3's upper bound still open | PASS (pinned) |
| YM-34 E4D-C RESTATED (native Laplace principle): the iterated multiplier WITHOUT smoothing, q_n = ∫wⁿ/∫wⁿ⁻¹, is strictly increasing to sup w = e^{2κ} (certified n ≤ 60) — so E4D-C as written in the RH ledger (deviation control of the multiplier iterated on one state) is unprovable, and the exponential of YM-16/26 is this growth. WITH the time kernel between multiplications the ratio is nondecreasing (log-convexity) and converges to λ₁(K^{1/2}WK^{1/2}) = f₀(1+δ), δ = 0.003582 / 0.014032 / 0.051911 (κ=1/8,1/4,1/2) vs sup price 0.274 / 0.598 / 1.405. E4D-C is restated for the alternating product (K^{1/2}M_wK^{1/2})ⁿ: λ₂/λ₁ ≤ 1−δ_gap uniform in m — the vacuum side of this is YM-31; only the excitation side is open | PASS (pinned) |
| Remaining: theorum/28 item 3 — the memory channel's OPERATOR norm (not basis-state norms) uniform in m, via the same engine (E4D-C quadratic form on content-½); with YM-31's floor the dock's e_n becomes (‖QTQ‖/FLOOR) and the sup numerator is the last exponential; (interval coefficients, content tail, Perron bound on the column transfer → m-uniform floor at rate ≥ 0.00268) and the E4D-C quadratic form on content-½ via the same engine; (exact rational ladder evaluations) — closes item 6 at the true rate AND item 3 (E4D-C quadratic form on content-½) in one build; E4D-C on SUPERPOSITIONS of content-½ basis states (a quadratic-form statement on the multi-insertion channel) = theorum/28 item 3 (content-½ multi-insertion memory channel, Cauchy-uniform in m); (weak) YM-27 the intertwiner (spin-network) transfer on the chain fabric — recoupling matrices as weak-coupling tiling weights, floor uniform in m? (declared: chain ~ O(4)-type rotor model in 1+1D, continuum gap at the frontier for every method); (strong) FULL-ORDER tiling convergence uniform in m and a (j=½ branching); (weak) the COMMUTATOR-RESIDUE TRANSFER on the fabric — does the non-abelian sector alone carry a gap; AF trajectory derivation; (convergence of the face-coefficient tiling expansion) — the cluster statement in fabric language; then (replace the borrowed decay anchor by program machinery — square-sourcing/seam count/flow generator; generalized Euler flow and EMK/UGD sections to be read first); then uniformity in the CUTOFF a→0 of the volume-uniform gap (Dobrushin degrades as the heat kernel sharpens) — this is now the sharp wall; AF trajectory, tightness, OS, non-triviality, metric universality, 2D lattice, Clay | OPEN |
| (superseded as the main route) cluster dock YM-18 proper: certified Kotecký–Preiss radius for the interleaving polymer gas with two-site pairs as monomers; per-bridge/local control in the content basis — certified strong-coupling cluster radius), 2D lattice, AF trajectory, tightness, OS reconstruction, non-triviality, metric universality, Clay predicate | OPEN |
| Continuum existence / mass gap | OPEN |

## Reproduce

New in YM-2: for beta=2 the interacting theta-graph transfer
`T_kappa = M^{1/2}(K x K)M^{1/2}` has a certified positive spectral gap for
every coupling `0 <= kappa < kappa_0 = Delta_red/6 = 0.13945388...`, via the
min-max sandwich `m_- T_0 <= T_kappa <= m_+ T_0` with
`m_+/m_- = e^{6 kappa}` from `|Tr U| <= 2`. At `kappa = 1/8` the reduced gap
is certified `>= 0.08672330623...`. Fail-closed above threshold.

```bash
python papers/yang-mills-certified-benchmark/certificates/ym1_certified_gap.py
python papers/yang-mills-certified-benchmark/certificates/ym2_theta_interacting_gap.py
python -m pytest papers/yang-mills-certified-benchmark/tests -v
```

New in YM-3 (exact, zero coupling): the free top-excited eigenspace is exactly
two-dimensional; the exact theta integral `Int chi12(A)chi12(B)chi12(AB^-1) = 1/2`
splits it at first order with derivatives `+-lambda_half/4`; the vacuum derivative
is exactly 0. Hence `r'(0) = lambda_half/4` (certified `0.10828185668...`), the
first-order CROSSING DIRECTION is the single symmetric line
`chi12(A)+chi12(B)` — rank-one and transported by the graph symmetry — and
YM-2's sandwich slope is slack by the exact factor 24.

Single-command reproduction (regenerates and pins all three):

```bash
python papers/yang-mills-certified-benchmark/certificates/ym12_square_sourced.py
```

YM-12 transplants the RH journey's decisive step (T03: re-source
positivity as a square with declared defects) to Yang-Mills:

- **T1** manifest square: `T_kappa(a) = S*S` with
  `S = [K_{a/2} x K_{a/2}] m_{kappa/2}` — both factors are the program's
  own halved-parameter objects (YM-9 semigroup + YM-5 doubling).
  Positivity is now **square-sourced**; YM-8's Jentzsch anchor is demoted
  to top-eigenvalue *simplicity* only. The compressed defect
  `A - N*G^{-1}N` is exactly a Gram of `(I-P)S phi` — a source-restriction
  object, the T03 shape.
- **T2** on the theta chain `L_m` (vertex-glued, no shared faces) the
  transfer tensor-factorizes exactly, so `lambda_2/lambda_1` is
  `m`-independent, and with YM-9: `Delta(L_m, a, a/16) >= 3/8` for
  **every** `m` and **every** `a` — an interacting gap uniform in volume
  and cutoff simultaneously. This closes YM-11 gate 2's interacting half
  on the factorized family and localizes the remaining volume obstruction
  to **interaction overlap** (shared faces): the chain (no sharing) is
  uniform, `B_n` (complete sharing) kills the sup route, and a physical
  lattice's bounded sharing sits strictly between the two certified
  extremes.
- **T3** governance capsule (RH 07-20 analog): claim/consumption chain for
  YM-1..12, standing corrections (YM-9 T4 amendment, YM-8 anchor
  demotion), and a do-not-reopen list — machine-checked for completeness.

YM-11 takes three dependency-map gates and gives each an exact verdict on
the toy carrier — the generalized theta ("n-banana") graphs `B_n`, two
vertices joined by `n` edges, so `b_1 = n-1` holonomies and
`F(n) = n(n-1)/2` plaquettes; `B_3` is the theta graph of YM-1..10.
**These are toy-carrier verdicts, not Clay-sense closures**, and two of
the three are partly negative.

- **Gauge — CLOSED.** Tree gauge-fixing on a connected graph is exact (no
  Gribov obstruction), leaving `b_1 = E-V+1` holonomies and a residual
  diagonal `G` acting by conjugation, so gauge-invariant states are
  exactly `L^2(G^{b_1})^Ad`. For `B_3` that **is** the carrier every
  capsule has used — now certified rather than assumed, with Wilson loops
  spanning it and YM-10's multiplicity law counting them.
- **Volume / IR — SPLIT.** The *free* reduced gap is `C_(1/2) = 3/4` for
  every `n` and every `a`: volume- **and** cutoff-uniform, closed. The
  *interacting* sandwich gives `Delta >= C_(1/2) - 2 F(n) theta`, which
  degrades quadratically in `n`: at `theta = 1/16` the bounds are `3/8`
  (n=3), exactly `0` (n=4), `-1/2` (n=5). The route is **proven
  insufficient** for volume-uniformity, with the critical volume computed
  exactly. What a replacement must do is named: control the interaction
  per-plaquette (dock route, YM-6) rather than by a global sup (envelope
  route, YM-5) — the same lesson one level up.
- **Universality — SPLIT.** For the whole regulator class
  `K = sum_j d_j c_j chi_j` with `c_j > 0`, the **counting layer**
  (content grading, multiplicity law, carrier dimensions, blindness
  ledger, probe counts) depends only on SU(2) representation theory and
  is therefore regulator-independent — closed, and it means YM-10's
  ledger is kinematics, not a choice of action. The **metric layer** (the
  gap value) genuinely differs between regulators (witnessed) and remains
  the real universality gate.

YM-10 applies the source-restriction framework of the companion
manuscripts *When Spectra Forget Order* and *Stable Recovery Beyond
Spectral Blindness* to our own compressions. Every dock computes on a
compressed carrier `V`; in that language `V` is a source restriction, so
the honest question is not whether the complement bound is tight but
whether the carrier is **blind** to the target.

- **T1** exact multiplicity law `m(j1,j2) = min(2j1,2j2)+1` on the
  Ad-invariant sector — it reproduces `dim V5 = 5` and `dim V7 = 7`
  independently, and re-derives YM-6's `(1/2,1/2)` two-dimensionality from
  representation theory rather than from the Gram matrix.
- **T2 correction to YM-9 T4** (self-audit): YM-9's uniform count counted
  *contents*; `k_Sigma` counts *eigenvalues with multiplicity*. At `s = 2`
  that is 4 versus 5. Both are `a`-independent, so YM-9's uniformity
  conclusion is unaffected — only the arithmetic. YM-9 is amended to
  report both.
- **T3/T4** blindness ledger and exact probe count
  `rank(Pi_s | ker E) = b(V,s)`. Result worth noting: `b(V5, s) = 0` for
  `s <= 2` — the YM-6 carrier is **exactly blindness-free in its own
  threshold window**, so that choice was correct, not lucky. At `s = 3`,
  V5 needs 6 probes, V7 needs 4, V9 needs 0 more than V7 minus 4.
- **T5** composition ledger: blindness is nonincreasing along
  `V5 ⊂ V7 ⊂ V9` with nonnegative stage increments.
- **T6** the whole ledger depends only on contents and Casimir levels, so
  it is exactly `a`-independent — one ledger covers the entire YM-9
  refinement family.

Honest remainder: the ledger is exact for the **free** target; the
interacting faces mix the content grading, so `b(V,s)` is a lower bound
there and YM-6/7's declared truncation remainder is the extra obligation.

YM-9 — the program's first statement that is uniform in the cutoff rather
than at a fixed cutoff. Switching the link action from Wilson to the
HEAT KERNEL `K_a(g) = sum_j d_j e^{-a C_j} chi_j(g)` makes refinement
EXACT (`K_a * K_b = K_{a+b}`; n links of spacing a/n compose to one link
of spacing a with no discretization error), and the cutoff dependence then
cancels identically:

- **T2** free reduced gap `= C_(1/2) = 3/4` EXACTLY for every `a > 0`;
- **T3** along the declared trajectory `kappa(a) = theta*a`,
  `Delta(a, kappa(a)) >= 3/4 - 6*theta` for EVERY `a > 0` — at
  `theta = 1/16` a cutoff-independent gap of exactly **3/8**; fail-closed
  at `theta >= 1/8`;
- **T4** the seam count `k(a, e^{-as})` is exactly `a`-independent
  (combinatorial), so the dock's threshold grammar transports across the
  whole family at once;
- **T5** Wilson contrast, computed: at fixed `beta` the reduced gap grows
  like `1/a` and diverges — no cutoff-independent gap without a running
  coupling.

Honest remainder (in-certificate): the trajectory is DECLARED not derived
(the physical one is fixed by asymptotic freedom, `beta ~ log(1/a)`, not
modelled); the graph is FIXED (no infinite-volume growth); the heat-kernel
action is a CHOICE (universality gate open); `Delta` is a reduced graph
gap, not a reconstructed physical mass gap. Net movement: one gate from
"untouched" to "touched on a toy carrier".

Capstone (YM-8): the interacting theta-graph transfer has a strictly
positive spectral gap at EVERY coupling — certified pointwise kernel
floors (`k >= e^{-3 kappa} [c0^{-1} e^{-beta}]^2 > 0`, exact enclosures)
plus the pinned classical Jentzsch/Krein-Rutman anchor (simple Perron
eigenvalue; CIRC-1 discipline: cited, never rederived). Quantified window
`kappa in [0, 7/10]` carries exact counts and enclosures (YM-2..7). The
honest remainder is named exactly: nothing here is uniform in the lattice;
a gap at every fixed cutoff is compatible with the gap closing in the
continuum limit — the fixed-regulator lesson. The Millennium content is
precisely the open uniformity/OS/vacuum gates of the MP dependency map.

New in YM-7: carrier enlarged to V7 (adds `chi_1(A), chi_1(B)`; Gram stays
block-diagonal, all seven exact T0 eigenvectors), dropping the complement
top from `lambda_1` to `lambda_1*lambda_half` (certified ordering) — which
UNLOCKS the previously refused `kappa = 7/10` cell: exact seam count
`k_Sigma = 1` at `(7/10, 13/10)`. Also ships certified brackets for all
seven compressed eigenvalue curves (inertia bisection; every step an exact
LDL count), published with honest kappa-drift — adopting the RH-line's
post-pause lesson that limited-denominator "exact constants" are never
promoted from bisection.

New in YM-6 (framework-native): the gap becomes an exact threshold COUNT
`k_Sigma(kappa,mu) = N(spec(T_kappa) > mu)` via a Haynsworth congruence —
the D03 "infinite ko finite" move — on the native 5x5 carrier
`{1, chi12(A), chi12(B), chi12(A)chi12(B), chi12(AB^-1)}` (all exact T0
eigenvectors; exact Gram has a single 1/2 overlap). Certified exact counts:
`k=1` at (kappa,mu) = (1/8,3/5), (1/4,3/5), (1/2,1) — so
`lambda_2 <= mu < lambda_1` EXACTLY, past YM-5's reach — and constant
count across the window = no seam spectral flow (thermodynamics/02
instance). The (7/10,5/4) cell is honestly REFUSED (bracket [1,5]).

New in YM-5: certified TWO-SIDED spectral gap of the interacting theta
transfer at finite coupling — including BEYOND YM-2's sandwich threshold
`kappa_0 = 0.1394...`. Complement control uses the exact doubling identity
`m_kappa^2 = m_{2 kappa}`: `|PMQ|^2 <= lammax(P m_{2k} P - (P m_k P)^2)`,
with `|QSQ| <= lambda_half^2 e^{3k}` and a Weyl bound. Certified rows
(ratio upper bound, gap lower bound): kappa=1/8: 0.5003, gap>=0.6924;
kappa=1/4: 0.5610, gap>=0.5781; kappa=3/10: 0.5867, gap>=0.5332.
Fail-closed at kappa=1. Lineage: realizes the old lambda-manuscripts'
Kill-Lemma skeleton unconditionally (see LINEAGE.md); no claim imported.

New in YM-4: (T1) the graph swap S:(A,B)->(B,A) commutes with T_kappa for
EVERY coupling — `Tr(BA^-1) = Tr(AB^-1)` on SU(2) — so YM-3's rank-one
crossing line `chi12(A)+chi12(B)` is symmetry-protected at all kappa: a
theorem, not a first-order accident. (T2) first certified finite-kappa
spectral data from the exact character expansion
`exp[(k/2)TrU] = sum d_j f_j chi_j`, `f_j = 2 I_{2j+1}(k)/k`, with exact
Clebsch-Gordan ring arithmetic, the exact pairing tensor
`Int chi_p(A)chi_q(B)chi_r(AB^-1) = delta_{p=q=r}/d_p`, and a certified
truncation remainder. Sector-resolved lower bounds on the grid show the
swap-even branch rising and the swap-odd branch falling with coupling.

CI runs exactly that one command plus pin diff and pytest.

## YM-35 — T54′ energy version + commutator remainder (E4D-C with the RKF operator ladder)

Framework inputs: theorum/54 (constant-μ chain is exactly its uncovered case — `sum_mu_diverges`), theorum/50 A4 (polarization visibility = YM-28's superposition gap), theorum/53 (energy, not mass), theorum/51 + 61 (commutator calculus). On the fabric, with the YM-30 engine as the only evaluator:

- **T1 locality** — `[K^{1/2}, B_l]` vanishes on states with no content at sites l, l+1; non-adjacent face pairs do not mix (site content exactly ½). Commutator support is `{l−1, l+1}` — bridge-local, every m.
- **T2 exact adjacent commutator energy** — `E_Σ([K^{1/2},B_l] B_{l+1}) = λ²·(3/4)·(1−√λ₁)²`, κ- and m-independent (site weights 1/4 : 3/4 = YM-25 recoupling squares, from the engine).
- **T3 T54′ ratio theorem (commuting model, exact m ≤ 8)** — vacuum and excitation carry the same per-face factor; ratio exactly m-independent; `Σμ<∞` not needed for a ratio.
- **T4 second-order remainder** — relative size ρ = μ²·(3/4)(1−√λ₁)² = 2.9e-4 / 1.2e-3 / 4.6e-3 at κ = 1/8, 1/4, 1/2 (a = 1); kill criterion ρ < 1/10 passed.

**Verdict:** E4D-C's excitation bound opens at the first non-commutative order — the first obstruction is an exact local number, not a sup. **Not done:** the full face expansion (all words, contents j ≥ 1) — the worldline cluster problem (YM-33). E4D-C OPEN. Pin `EXPECTED_YM35.sha256`.
