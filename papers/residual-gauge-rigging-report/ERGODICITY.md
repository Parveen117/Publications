# Residual-Diff ergodicity on (X_c, nu_2) and the classification of the physical rigging

16 October 2026. Separate branch `claude/residual-gauge-invariant-dual`,
continuing from 943182d. RGR-1 at d9c70208520f53f31d93bd94c3190f48b34acb4c
is a **stated dependency** and is neither reconstructed nor modified. The
shared endpoint 059d4b3 was read for coordination only: its new packet
(RMG-1, nonlinear material geometry) does not overlap this note. Nonlinear
physical geometry, fixed-hbar composites and length domains are the other
track's, and are not touched here.

Labels: **[P]** proved here; **[D]** RGR-1 statement used as a dependency
(ledger in Section 7); **[E]** external theorem; **[C]** open. Finite exact
witnesses (Section 8) are regression checks only and certify no continuum
statement.

## Result

**Theorem E [P, using D and E].** The residual group G_c acts ergodically
on (X_c, nu_2). Every G_c-invariant F in L^∞(nu_2) is a.e. constant.

Together with Theorem U of `UNIQUENESS.md` (unchanged, ledger below):

**Corollary (classification) [P].**

1. Every linear functional on RGR's Psi that is continuous for (19) and
   G_c-invariant is a multiple of `ell`.
2. Every separately continuous, Hermitian, positive, G_c-invariant form on
   Psi equals `c · conj(ell) ⊗ ell` with c >= 0. The physical gauge
   quotient is C (when c > 0). RGR-T5's channel is the only one.
3. With the MPQ spectator factor, under the explicit contract (BC) of
   Section 5, every admissible positive form is
   `conj(ell(psi)) ell(phi) <v, (c_+ P_+ + c_- P_-) u>`. It is unitarily
   equivalent to RGR's (29) through `W = sqrt(c_+) P_+ + sqrt(c_-) P_-`,
   and it equals (29) exactly under RGR's isometric normalization.

What is classified, and what is not, is listed in Section 6.

The route is a Gaussian divergence closure built from genuine generator
identities (Section 3). It uses **de la Vallée Poussin** weights on the
paired cos/sin Fourier modes of the field-dependent gauge parameter. With
sharp Dirichlet cutoffs, the D^{1,2} estimate fails at the Cameron–Martin
borderline. Whether the sharp-sum version of Condition (L) holds in
L^1(nu_2) is left undecided, and nothing here depends on it (Section 4).

---

## 1. Setting (one chirality)

- X = S^1 x H^3_0 ∋ (alpha, sigma). kappa = e^sigma / ∫ e^sigma is the
  density, `I = I_sigma(t) = ∫_0^t kappa`, and `chi = alpha + I`.
- `T_g(alpha, sigma) = (g(alpha), sigma + P_0 log g'(chi))`.
- nu_2 = d alpha ⊗ mu_2, where mu_2 is the centred Gaussian with
  covariance 2 lambda_n, lambda_n = (1+n^2)^-4. Its Cameron–Martin space
  is `H_2 = H^4_0`, with `||u||^2_{H_2} = sum_n |u_n|^2 / (2 lambda_n) = (1/2) ||u||^2_{H^4}`.
  W_2(h) is the Gaussian linear functional with variance ||h||^2_{H_2}.
  Theorem O's notation is unchanged.
- For a real smooth f, RGR gives the field
  `X_f = (f(alpha), v_f)` with `v_f = P_0 f'(chi)`, the score A_f, and the
  generator `B_f = -X_f + A_f/2` (23).
- A is RGR's algebra of smooth functions whose Fréchet derivatives on H^3
  charts are bounded by `C_r (1+||sigma||)^{N_r} exp(C_r ||sigma||)`.
- Fourier conventions: `<n> = (1+n^2)^{1/2}`,
  `||w||^2_{H^s} = sum <n>^{2s} |w_n|^2`, and
  `||w||_{A^r} = sum |n|^r |w_n|` (with A^0 using |n|^0 = 1).
- On the product X_c = X x X, everything below is applied to each chiral
  factor in turn.

## 2. The invariance identity used

**Lemma 1 [P from D].** Let F in L^∞(nu_2) satisfy F o T_g = F a.e. for
every g in G_c. Then `∫ F B_f a dnu_2 = 0` for every real smooth f and
every a in A.

*Proof.*

- Change of variables with the nu_2 cocycle (11), (21) gives
  `∫ F U_g a dnu_2 = ∫ (F o T_g) a dnu_2 = ∫ F a dnu_2`.
- By (26), `U_t a - a = ∫_0^t U_s B_f a ds` as a Bochner integral in
  L^1(nu_2). Pairing with F, which is bounded, and applying the first
  step to `B_f a` (which lies in A),

      0 = ∫ F (U_t a - a) dnu_2 = t ∫ F B_f a dnu_2. ∎

For a finite family of real smooth e_k and coefficients b_k in A, the
Leibniz rule gives the pointwise identity

    sum_k B_{e_k}(b_k psi) = -Y psi - psi div Y,    Y = sum_k b_k X_{e_k},   (2.1)

where `div Y := sum_k ( X_{e_k} b_k - b_k A_{e_k}/2 )`.

Since b_k psi lies in A, Lemma 1 gives `∫ F (-Y psi - psi div Y) dnu_2 = 0`.

## 3. Theorem E: divergence closure with de la Vallée Poussin weights

**The field-dependent gauge parameter.** Fix a real smooth h in H^3_0 and
put `m = ∫ h kappa`.

- Set `G(z) = DI_sigma[h](I^-1(z))`, where
  `DI_sigma[h](t) = ∫_0^t h kappa - m I(t)`, and set `f_x(y) = G(y - alpha)`.
- Then `Phi := G' = (h - m) o I^-1`, and `f_x(alpha) = G(0) = 0`.
- The field X_{f_x} at x is (0, h): the angular part vanishes, and
  `P_0 (f_x' o chi) = P_0 (h - m) = h`.
- The Fourier coefficients of f_x are `c_k(x) = e^{-2 pi i k alpha} G_k(sigma)`, with

      G_k(sigma) = ∫ DI_sigma[h](t) e_{-k}(I(t)) kappa(t) dt.

  This is an explicit smooth function of sigma, with no inversion of
  diffeomorphisms, and each derivative is bounded by
  `poly(|k|, ||sigma||_{H^3}) exp(C ||sigma||_∞)`. So c_k is in A, and so
  are its real cos/sin versions.

**Paired VP truncation.**

- V_K is the Fourier multiplier with symbol rho_K(q): it equals 1 for
  |q| <= K, equals (2K - |q|)/K for K < |q| < 2K, and equals 0 beyond.
- rho_K is even, so the cos_k and sin_k modes always carry the same
  weight. It is Lipschitz with constant 1/K, satisfies |rho_K| <= 1, and
  fixes constants.
- The truncated field is `Y_K = sum_{|k|<2K} rho_K(k) c_k X_{e_k}`, a
  finite real combination.
- By translation covariance, its angular part is
  `(V_K f_x)(alpha) = (V_K G)(0)`, which is independent of alpha. Its
  shape part is `u_K = P_0 [ (V_K Phi) o I ]`, also independent of alpha.

**Lemma 2 (divergence identification) [P from D, E].** For each K,
u_K is in D^{1,2}(mu_2; H_2), and

    div Y_K = -delta_2(u_K)      nu_2-a.e.,

where delta_2 is the Skorokhod divergence of mu_2.

*Proof.*

- u_K is a finite sum `sum rho_K(k) G_k(sigma) P_0[e_k o I]` of products
  of A-functions with smooth maps H^3 → H^4. Its derivative is
  Hilbert–Schmidt, by the bound in Lemma 3 below with K fixed, and has
  integrable bounds. So u_K is in D^{1,2}(H_2) [E: Fréchet-C^1 maps with
  HS H-derivative and integrable bounds, e.g. Nualart §1.2 / Bogachev
  §5.7].
- Test against a bounded smooth cylinder a. The left side, using (27)
  with c = 2 applied to c_k a in A, gives

      ∫ div(Y_K) a = sum rho_K(k) [ ∫ X_{e_k}(c_k a) - ∫ c_k X_{e_k} a - (1/2)∫ c_k a A_{e_k} ] = -∫ Y_K a.

- The right side: the angular part of Y_K does not depend on alpha, so
  the alpha integration by parts produces no term. Gaussian duality gives
  `∫ <D a, u_K>_{H_2} = ∫ a delta_2(u_K)`. So the right side also equals
  `-∫ Y_K a`.
- Bounded smooth cylinders determine L^1(nu_2) classes. ∎

**Lemma 3 (commutator) [P].** Let D and w be functions on S^1 and K >= 2.
Then

    ||[M_D, V_K] w||_{H^4} <= 48 ||D||_{A^1} ||1_{K/2<|q|<5K/2} w||_{H^3}
                             + 7^4 ||P_{>K/2} D||_{H^4} ||w||_{A^0}.

*Proof.*

- The p-th coefficient is `sum_l d_l w_q (rho(q) - rho(p))`, with
  q = p - l, and `|rho(q) - rho(p)| <= min(1, |l|/K)`.
- The coefficient vanishes when |q|, |p| <= K, and also when both exceed
  2K.
- **Case |l| <= K/2.** A nonzero term needs K/2 < |q| < 5K/2, so
  `<p> <= 2<q>` and `<q> <= 3K`. Hence
  `<p>^4 |rho(q) - rho(p)| <= 16 <q>^4 |l|/K <= 48 <q>^3 |l|`. Young's
  inequality, l^1 in l against l^2 in q, gives the first term.
- **Case |l| > K/2.** One of |p|, |q| is below 2K, so `<p> <= 7 <l>`.
  Young's inequality, l^2 in l against l^1 in q, gives the second term. ∎

The sharp symbol has jumps instead: `|rho(q) - rho(p)| = 1` on every
crossing pair. So its first case carries `<q>^4` rather than `<q>^3`
(see Section 4).

**Lemma 4 (gradient formula) [P].** Write `D~_e := (D_e I) o I^-1`, where
`D_e I(t) = ∫_0^t e kappa - I(t) ∫ e kappa`. For e in H_2,

    D_e u_K = P_0 [ ( [M_{D~_e}, V_K] Phi' ) o I ].

*Proof.*

- Differentiate `Phi_k = ∫ (h - m)(t) e_{-k}(I(t)) kappa(t) dt` and
  change variables y = I(t), using `(D_e I)' dt = D~_e'(y) dy`. This gives

      D_e Phi_k = -(Phi' D~_e)_k - (D_e m) delta_{k0}.

- Differentiating `P_0[e_k o I]` gives `P_0[e_k'(I) D_e I]`.
- Summing with weights rho_K(k), and using that V_K fixes constants and
  commutes with d/dy, the result is
  `P_0[((V_K Phi') D~_e - V_K(Phi' D~_e)) o I]`. ∎

**Proposition 5 (D^{1,2} convergence) [P, using E].** As K → ∞,
u_K → h in D^{1,2}(mu_2; H_2).

*Proof.* Fix sigma in H^3.

- Then I and I^-1 are H^4 diffeomorphisms, so Phi is in H^4, w := Phi'
  is in H^3, and right composition with I is bounded on H^4 with norm
  `C(sigma)` [E: composition in D^s(S^1), s > 3/2; e.g. Inci–Kappeler–Topalov,
  Mem. AMS 2013].
- (i) `||u_K - h||_{H^4} <= C(sigma) ||P_{>K} Phi||_{H^4} → 0`.
- (ii) Let e_n be an orthonormal basis of H_2, so e_n is of size n^-4 in
  C^0. Then:
  - `||D~_{e_n}||_{A^1} <= C ||D~_{e_n}||_{H^2} <= C(sigma) n^-3`, and
  - `||D~_{e_n}||_{H^4} <= C(sigma) ||e_n kappa||_{H^3} <= C(sigma) n^-1`.
- By Lemmas 3 and 4,

      sum_n ||D_{e_n} u_K||^2_{H_2} <= C(sigma)^2 sum_n [ 48^2 ||D~_{e_n}||^2_{A^1} ||band_K w||^2_{H^3}
                                                     + 7^8 ||P_{>K/2} D~_{e_n}||^2_{H^4} ||w||^2_{A^0} ].

  Each bracket tends to 0, because `band_K w → 0` in H^3 and each tail
  tends to 0. The sum in n is dominated by `sum_n n^-2`.
- So `||D u_K||_{HS} → 0` for every sigma in H^3. All constants have the
  form `(1+||sigma||_{H^3})^N exp(C ||sigma||_{H^3})`.
- Fernique [E] makes these constants integrable, and dominated
  convergence gives `E||u_K - h||^2 + E||D u_K||^2_{HS} → 0`. The limit h
  is deterministic, so Dh = 0. ∎

**Proof of Theorem E.**

- By (2.1), Lemma 1 and Lemma 2, for every psi in A,

      0 = ∫ F [ -Y_K psi + psi delta_2(u_K) ] dnu_2.

- Now let K → ∞. First, `Y_K psi = (V_K G)(0) ∂_alpha psi + <D psi, u_K>_{H_2}`,
  with `|(V_K G)(0) - G(0)| <= ||P_{>K} G||_{A^0} → 0`. Second, by the
  Skorokhod isometry bound [E: Nualart Prop. 1.3.1],

      E|delta_2(u_K - h)|^2 <= E||u_K - h||^2_{H_2} + E||D u_K||^2_{HS} → 0.

  Since psi and its derivatives have all moments, the bracket converges
  in L^1(nu_2) to `-∂_h psi + psi W_2(h) = delta_2(psi h)`. **This is the
  actual L^1(nu_2) convergence of the complete divergence expression**
  (2.1), not of the vector fields alone.
- Since F is bounded, `∫ F delta_2(psi h) dnu_2 = 0` for all psi in A and
  all smooth h.
- Take `psi = a(alpha) b(sigma)` with a smooth and b a multi-Hermite
  polynomial in W_2(e_i), for smooth orthonormal e_i.
  - The identity `delta_2(He_{k}(W_2(e_j)) prod_{i≠j} He_{k_i}(W_2(e_i)) e_j)`
    equals the same product with `He_{k+1}` in slot j.
  - So these values span the orthogonal complement of the constants in
    L^2(mu_2) [E: Wiener chaos].
  - Hence, for a.e. alpha, F(alpha, ·) is mu_2-a.e. a constant c(alpha).
- Rotations lie in G_c and act by alpha → alpha + beta with sigma fixed,
  so c is constant.
- On X_c, apply this to each chiral factor: the + and - factors' shapes
  and angles independently. ∎

**Why this is not the excluded route.**

- No finite Cameron–Martin shift is ever claimed to be a gauge
  transformation; Theorem O still holds.
- Orbit density is not used.
- The abelian 0–1 law is not transferred through orbits. The only
  Gaussian input is the L^2 orthogonality statement of Wiener chaos,
  applied to a weak identity that is derived from finitely many genuine
  generator identities followed by a proved L^1 limit.
- Ergodicity enters only through Lemma 1's exact invariance.

## 4. Sharp Dirichlet sums: the precise defect

Condition (L) in `UNIQUENESS.md` uses sharp partial sums |k| <= K. The
proof above goes through for sharp sums **if and only if** u_K^sharp → h
in a topology that makes delta_2 converge. In the sharp case:

- The gradient is `D_e u_K^sharp = P_0[([M_{D~_e}, P_K] Phi') o I]`.
- Every crossing pair (|q| <= K < |p|) contributes `d_l w_q` with weight
  `<p>^4 ~ K^4` and **no 1/K gain**.
- So `||D u_K^sharp||_{HS}` is controlled by `K^4 |w_q|` for |q| ≈ K.
- Since w = Phi' lies in H^{7/2-} almost surely and not better, this is
  O(1) at best and is not small.
- The exact witness confirms this on the borderline profile w_q = |q|^-4:
  - the sharp band H^4-norm squared is 1.55, 1.14, 1.03, 1.009, 1.002 for
    K = 8, ..., 128, which stays at or above 1;
  - the VP norm decays: 4.6e-2, 7.5e-3, 1.4e-3, 3.0e-4, 6.8e-5.

Status:

- **[P]** The D^{1,2} route fails for sharp sums in the sense that Lemma
  3 has no sharp analogue with a vanishing right side.
- **[C]** A lower bound `E||D u_K^sharp||^2_HS >= c > 0` would need
  `E|w_K|^2 >= c K^-8`. That is a nondegeneracy of the nonlinear
  coefficient field, and it is not proved.
- **[C]** Even then, L^1 convergence of the sharp divergence could hold
  through a different cancellation. A uniform bound
  `sup_K E||D u_K^sharp||^2_HS < ∞` would give weak L^2 convergence,
  which would already suffice for ergodicity. That bound is also not
  proved.
- So sharp (L) is **neither proved nor refuted**. The obstruction is that
  the jump of the Dirichlet symbol meets coefficients exactly at the H_2
  borderline. VP weights remove it with half a derivative to spare.
  Failure of sharp (L) would not have implied nonergodicity, and no
  invariant set has been or needs to be produced.

## 5. Theorem S repaired: domain, topology, upgrade, and the sector reading

**Domain and topology.**

- A_0 is the span of polynomials in finitely many S(f), with f real
  smooth. It is contained in MPQ's domain A.
- `D_red = A_0 ⊕ A_0` inside `H_red = H_+ ⊕ H_-` (MPQ (8)), with the
  **norm topology of H_red**.
- The form is defined on `Psi ⊗_alg D_red`. It is required to be
  positive and Hermitian, separately continuous in the Psi slot for (19),
  and G_c-invariant in that slot.

**Compatibility contract (BC).**

- (BC1) The reduced form `B(v, u) := eta(1 ⊗ v)[1 ⊗ u]` is bounded:
  `|B(v, u)| <= C ||v|| ||u||`. So `B = <·, T ·>` with T bounded and
  T >= 0.
- (BC2) Weak symmetry: `B(O v, u) = B(v, O u)` on D_red for
  `O = S(f), D(h)` (smooth f, h), and for J(f_0) for one smooth f_0 > 0.

**Lemma S1 (core) [P, using E: Nelson's analytic vector theorem].**

- A_0 lies in a finite Wiener chaos, and it is invariant under S(f) and
  D(h).
- On chaos of order at most N, `||S(f)^n p||` and `||D(h)^n p||` are at
  most `C^n ((N+n)!/N!)^{1/2} ||p||`. So A_0 consists of analytic
  vectors, and both restrictions are essentially self-adjoint.
- The closure of `D(h)|A_0` is a self-adjoint restriction of MPQ's
  self-adjoint D(h), so the two are equal. Hence A_0 is a core for each
  S(f) and D(h).
- For J(f_0) I do **not** prove a core statement. (BC2) for J is used
  only to remove off-diagonal blocks, and (BC3) below may replace it.

**Lemma S2 (upgrade) [P].** Let O be self-adjoint, let C be a core for O,
and let T be bounded with `<O v, T u> = <v, T O u>` for v, u in C. Then T
commutes with every bounded Borel function of O, and in particular with
`exp(-itO)`.

*Proof.*

- For u in C, the map v → <Ov, Tu> is bounded on the core. So Tu lies in
  Dom (O|C)^* = Dom O, with `OTu = TOu`.
- For general u in Dom O, take core approximants u_n → u with Ou_n → Ou.
  Closedness gives `TO ⊂ OT`.
- Hence `T (O - z)^-1 = (O - z)^-1 T` for z not real, and the spectral
  theorem extends this to all bounded Borel functions. ∎

Weak symmetry **without** (BC1), or without the core property, does not
give this. That is why boundedness is an explicit contract, and not a
derived property. An unbounded closed form must instead be required to
be invariant under the unitary groups `exp(-itO)` directly (contract
BC1').

**Theorem S' [P, using Theorem E, (BC), E].** Under (BC1), (BC2), every
admissible form is

    eta(psi ⊗ v)[phi ⊗ u] = conj(ell(psi)) ell(phi) <v, (c_+ P_+ + c_- P_-) u>,    c_± >= 0.

*Proof.*

- **Gauge slot.** For fixed v, u, the map (psi, phi) → eta is separately
  continuous, invariant and sesquilinear. Two applications of Theorem
  U(a) together with Theorem E give `conj(ell(psi)) ell(phi) B(v, u)`.
- **Within a sector.** By Lemmas S1 and S2, T commutes with the spectral
  projections of every S(f), which generate the maximal abelian algebra
  L^∞(mu). So T acts within each sector as multiplication by some tau.
  Commuting with `U_h(t)` gives `tau(s + t h) = tau(s)`, a.e., for a
  dense set of Cameron–Martin directions. Hence tau is constant
  [E: ergodicity of a Gaussian measure under a dense subgroup of its
  Cameron–Martin space].
- **Across sectors.** The off-diagonal block X: H_- → H_+ commutes in the
  same way, so X = x·1. (BC2) for J(f_0) reads
  `x J_- = J_+ x`, where `J_± = ±∫ f_0 e^s`, so x = 0. ∎

Without J or the sign in the declared algebra, the two sectors carry
equivalent representations. The commutant is then M_2(C), and the
off-diagonal weights are not fixed. Witness W4 shows the commutant
dimension dropping from 2 to 1 when a sign-mixing observable is added.

**(BC3) Alternative.** Declare the sign epsilon itself, or E_J((0, ∞))
for one J(f_0) with f_0 > 0, as a bounded observable. In H_+ the operator
J(f_0) is positive, and in H_- it is negative (MPQ-T3). So
`P_+ = E_{J(f_0)}((0, ∞))` is a spectral projection of a declared
observable. Bounded commutation with P_+ kills the off-diagonal block
without any core statement for J.

**Corrected sector interpretation.** This replaces the reading in
`UNIQUENESS.md` §3, which said the ratio is "not observable through any
declared MPQ observable". That reading was wrong.

- **Populations are measured.** P_+ is a spectral projection of the
  declared J(f_0). For a fixed reduced vector `v = v_+ + v_-`, the
  predicted population of the + sector is

      c_+ ||v_+||^2 / (c_+ ||v_+||^2 + c_- ||v_-||^2),

  so it depends on c_+/c_-. Measuring it needs no sign-mixing observable.
  Witness: 1/2 at equal weights and 1/5 at (1, 4) for the same vector.
- **Relative phase is not measured.** Every declared observable commutes
  with P_±, so `<v, O v>_B` is unchanged under `v_- → e^{i theta} v_-`.
  This is the superselection statement, and it is the only thing that
  sign preservation alone implies.
- **The weights are a unitary rescaling, not new physics.**
  `W = sqrt(c_+) P_+ + sqrt(c_-) P_-` satisfies `W^* W = T`, and W
  commutes with every declared observable. So the theories with weights
  (c_+, c_-) and (1, 1) are unitarily equivalent as representations of
  the declared algebra. The ratio changes only **which reduced vector
  represents a given physical state**, that is, the normalization of the
  map R in (29).
- RGR's isometric normalization ("the generalized state F → <v, RF> has
  norm ||v||_{H_red}") is exactly the extra condition `c_+ = c_- = 1`
  (witness: it passes at equal weights and fails at (1, 4)).
- So the full physical inner product is **unique** under
  (BC1) + (BC2 or BC3) + isometric normalization, and **unique up to
  unitary equivalence** under (BC1) + (BC2 or BC3) alone.

## 6. What is classified

Classified (given Theorem E and the ledger):

- the continuous G_c-invariant dual of RGR's Psi for the topology (19):
  it is C·ell;
- all separately continuous, positive, invariant Hermitian forms on Psi:
  they are c·conj(ell)⊗ell, with a one-dimensional quotient;
- all G_c-invariant probability measures absolutely continuous with
  respect to nu_2: nu_2 is the only one;
- with the spectator, all forms obeying (BC): Theorem S'.

Not classified:

- functionals that are not continuous for (19), or that live on other
  test spaces;
- invariant measures singular to nu_2. Other Gaussian covariances are
  already excluded by RGR (27), but non-Gaussian singular measures are
  not addressed;
- spectator forms that violate (BC1) (unbounded, without unitary
  invariance), or declared algebras without J or sign;
- anything about nonlinear geometry observables.

## 7. Dependency ledger

**RGR-1 statements used [D]**, at d9c7020, file
`papers/native-residual-gauge-rigging/THEOREM.md`:

| Item | Used for | Audit here |
|---|---|---|
| (10)–(11) continuous cocycles, `sqrt(R^(1)) = R^(2)` | Lemma 1, Theorem U | agrees with REPORT.md (C5) |
| (19) Psi and its topology; (21) L^1 isometry | Theorem U, Lemma 1 | read; T_g bounded on balls checked |
| (22) separating bump | Theorem U step 2 | cylinder bump in UC_b checked |
| (23) X_f, A_f, B_f; A and `B_f A ⊂ A` | Lemmas 1, 2 | the c_k-in-A step is proved here |
| (26) Bochner identity in L^1(nu_2) for a in A | Lemma 1 | not re-derived |
| (27) at c = 2: `∫ X_f a = (1/2) ∫ A_f a` | Lemma 2 | not re-derived |
| (17) u_g, v_f in Gaussian D^{1,2} | Lemma 2 | replaced by the explicit u_K bound |
| CGD-T2 step 3 (no invariant kinematic vector) | Theorem U step 2 | as in UNIQUENESS.md |
| MPQ-T3 (8)–(12), and A domain | Section 5 | A_0 ⊂ A only |

**External theorems [E]:**

- Skorokhod divergence: D^{1,2}(H) ⊂ Dom delta, with
  `E delta(u)^2 <= E||u||^2 + E||Du||^2_HS` (Nualart, *The Malliavin
  Calculus and Related Topics*, Prop. 1.3.1).
- Product rule and Gaussian duality (ibid. §1.3).
- Wiener chaos density (ibid. Thm 1.1.1).
- D^{1,2} membership of Fréchet-C^1 maps with HS H-derivative (ibid.
  §1.2; Bogachev, *Gaussian Measures*, §5.7).
- Fernique's theorem.
- Composition and inversion bounds in D^s(S^1), s > 3/2
  (Inci–Kappeler–Topalov, Mem. AMS 226 (2013)).
- Nelson's analytic vector theorem.
- Ergodicity of a Gaussian measure under a dense subgroup of its
  Cameron–Martin space (Bogachev, op. cit., §2.5).
- The spectral theorem.

**Proved here [P]:** Lemmas 1–4, Proposition 5, Theorem E, Lemmas S1–S2,
Theorem S', and the sector reading. Theorem U, Theorem O and their proofs
are unchanged in `UNIQUENESS.md`.

## 8. Finite exact witnesses (`ergodicity/`)

These are native JS exact rational arithmetic, run with
`python3.12 ergodicity/run.py --check`. They check identities and specific
finite instances only.

- **W1.** Commutator Fourier formula
  `[M_D, S]_p = sum d_l w_q (rho(q) - rho(p))` for sharp and VP symbols,
  for two trigonometric D, two w, and K = 4, 6. Also: commutation with
  d/dy, the cancellation identity of Lemma 4, and the fixing of
  constants.
  - *Negative controls:* the reversed sign, and the identity without the
    commuting derivative.
- **W2.** Band H^4 norms on the borderline profile `|q|^-4`, for
  K = 8, ..., 128. VP decays by at least a factor 3/5 per doubling, and
  `K · norm <= 2`. Lemma 3's low-frequency bound holds with the stated
  constant 48. Sharp stays at or above 1.
  - *Negative control:* sharp fails the `K · norm <= 2` smallness bound
    that VP satisfies, at every K.
- **W3.** Finite Gaussian closure. The constraints
  `E[F(-∂_j psi + x_j psi / v_j)] = 0` on a cubic space give nullity 1
  (constants only) when imposed in all directions.
  - *Negative controls:* one direction leaves nullity 4 (functions of
    x_2), and the undoubled score leaves nullity 0 (even constants fail).
- **W4.** Sector model on C^2 ⊕ C^2. The declared sign-preserving
  algebra plus P_+ has a commutant of dimension 2; adding a sign-mixing
  observable makes it 1. Populations 1/2 against 1/5. W^*W = T, and W
  intertwines the declared algebra.
  - *Negative controls:* W does not intertwine the mixing observable,
    and unequal weights fail isometric normalization.

The results block of `CERTIFICATE.json` records
`continuum_proofs_machine_certified: false`.
