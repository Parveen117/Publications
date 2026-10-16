# Residual-Diff ergodicity on (X_c, nu_2) and the classification of the physical rigging

16 October 2026. Separate branch `claude/residual-gauge-invariant-dual`,
continuing from 943182d. **Revised the same day** after a five-point audit;
the changes and the audit verdict are listed in `AUDIT.md`. RGR-1 at
d9c70208520f53f31d93bd94c3190f48b34acb4c is a **stated dependency** and is
neither reconstructed nor modified. The shared endpoint 059d4b3 was read for
coordination only. Nonlinear physical geometry, fixed-hbar composites,
length domains, detector calibration and apparatus compatibility are the
other track's, and are not touched here.

Labels: **[P]** proved here; **[D]** RGR-1 statement used as a dependency
(ledger in Section 7); **[E]** external theorem; **[C]** open. Finite exact
witnesses (Section 8) are regression checks only and certify no continuum
statement.

## Result

**Theorem E [P, using D and E].** The action of G_c on (X_c, nu_2) is
nonsingular [D], and it is **ergodic as a nonsingular action**: every Borel
set B with `nu_2(T_g^{-1} B Δ B) = 0` for all g in G_c is nu_2-null or
conull; equivalently, every F in L^∞(nu_2) with `F o T_g = F` nu_2-a.e. for
each g is a.e. constant.

nu_2 itself is quasi-invariant and **not invariant** (Proposition M); no
invariance of the measure is asserted anywhere. Ergodicity of the
nonsingular action is exactly the hypothesis Theorem U(b),(d) of
`UNIQUENESS.md` consumes, so with Theorem U (unchanged, ledger below):

**Corollary (classification) [P].**

1. Every linear functional on RGR's Psi that is continuous for (19) and
   G_c-invariant is a multiple of `ell`.
2. Every separately continuous, Hermitian, positive form on Psi that is
   G_c-invariant **in each argument** equals `c · conj(ell) ⊗ ell` with
   c >= 0. The physical gauge quotient is C (when c > 0). RGR-T5's channel
   is the only one.
3. With the MPQ spectator factor, under the explicit contract (BC) of
   Section 5 and slotwise gauge invariance, every admissible positive form
   is `conj(ell(psi)) ell(phi) <v, (c_+ P_+ + c_- P_-) u>` with
   c_± >= 0. When **both weights are strictly positive** it is unitarily
   equivalent to RGR's (29) through `W = sqrt(c_+) P_+ + sqrt(c_-) P_-`,
   and equals (29) exactly under RGR's isometric normalization. When
   exactly one weight vanishes, the physical quotient retains only the
   supported sector H_+ or H_-; when both vanish, the form is the zero
   form with trivial quotient (Theorem S', trichotomy).

What is classified, and what is not, is listed in Section 6, together with
what is actually proved about ordinary invariant measures (Proposition M).

The route is a Gaussian divergence closure built from genuine generator
identities (Section 3). It uses **de la Vallée Poussin** weights on the
paired cos/sin Fourier modes of the field-dependent gauge parameter. With
sharp Dirichlet cutoffs, the D^{1,2} estimate fails at the Cameron–Martin
borderline. Whether the sharp-sum version of Condition (L) holds in
L^1(nu_2) is left undecided, and nothing here depends on it (Section 4).

---

## 1. Setting and conventions (one chirality)

- X = S^1 x H^3_0 ∋ (alpha, sigma). kappa = e^sigma / ∫ e^sigma is the
  density, `I = I_sigma(t) = ∫_0^t kappa`, and `chi = alpha + I`.
- `T_g(alpha, sigma) = (g(alpha), sigma + P_0 log g'(chi))`.
- nu_2 = d alpha ⊗ mu_2, where mu_2 is the centred Gaussian with
  covariance 2 lambda_n, lambda_n = (1+n^2)^-4. Its Cameron–Martin space
  is `H_2 = H^4_0`, with `||u||^2_{H_2} = sum_n |u_n|^2 / (2 lambda_n) = (1/2) ||u||^2_{H^4}`.
  W_2(h) is the Gaussian linear functional with variance ||h||^2_{H_2}.
  Theorem O's notation is unchanged.
- **Fourier conventions.** `e_k(y) = e^{2 pi i k y}`, `u_k = ∫ u e_{-k}`,
  so `(u')_k = 2 pi i k u_k` and `e_k' = 2 pi i k e_k`. `<n> = (1+n^2)^{1/2}`,
  `||u||^2_{H^s} = sum <n>^{2s} |u_n|^2`, `||u||_{A^r} = sum |n|^r |u_n|`
  (A^0 uses |n|^0 = 1). Real gauge parameters use the paired real modes
  cos_k, sin_k; every multiplier below has an even symbol, so it is
  well defined on the real pairs, and all identities are written in the
  complex basis with this understood. The witness file computes "up to the
  common factor 2 pi i", which is stated there.
- For a real smooth f, RGR gives the field
  `X_f = (f(alpha), v_f)` with `v_f = P_0[f'(chi)]`, the score A_f, and
  the generator `B_f = -X_f + A_f/2` (23).
- A is RGR's algebra of smooth functions whose Fréchet derivatives on H^3
  charts are bounded by `C_r (1+||sigma||)^{N_r} exp(C_r ||sigma||)`,
  uniformly in the angles.
- On the product X_c = X x X, a smearing f acts on one chirality (say +)
  and everything below is applied to each chiral factor in turn; test
  functions always live on the full X_c.

## 2. The invariance identity used

**Lemma 1 [P from D].** Let F in L^∞(nu_2) satisfy F o T_g = F a.e. for
every g in G_c. Then `∫ F B_f a dnu_2 = 0` for every real smooth f and
every a in A.

*Proof.*

- By (11) and (21), R^{(2)}_g is exactly the Radon–Nikodym density of the
  pushforward of nu_2, so the change of variables gives
  `∫ F U_g a dnu_2 = ∫ (F o T_g) a dnu_2 = ∫ F a dnu_2`. Only
  quasi-invariance is used.
- By (26), `U_t a - a = ∫_0^t U_s B_f a ds` as a Bochner integral in
  L^1(nu_2). Pairing with F, which is bounded, and applying the first
  step to `B_f a` (which lies in A),

      0 = ∫ F (U_t a - a) dnu_2 = t ∫ F B_f a dnu_2. ∎

For a finite family of real smooth e-modes and real coefficients b_k in A,
the Leibniz rule gives the pointwise identity

    sum_k B_{e_k}(b_k psi) = -Y psi - psi div Y,    Y = sum_k b_k X_{e_k},   (2.1)

where `div Y := sum_k ( X_{e_k} b_k - b_k A_{e_k}/2 )`. (For the complex
modes below this is applied to the real and imaginary pairs.)

Since b_k psi lies in A, Lemma 1 gives `∫ F (-Y psi - psi div Y) dnu_2 = 0`.

## 3. Theorem E: divergence closure with de la Vallée Poussin weights

**The field-dependent gauge parameter.** Fix a real smooth h in H^3_0 and
put `m = ∫ h kappa`.

- Set `G(z) = DI_sigma[h](I^-1(z))`, where
  `DI_sigma[h](t) = ∫_0^t h kappa - m I(t)`, and set `f_x(y) = G(y - alpha)`.
- Then `Phi := G' = (h - m) o I^-1`, and `f_x(alpha) = G(0) = 0`.
- The field X_{f_x} at x is (0, h): the angular part vanishes, and
  `P_0[f_x'(chi)] = P_0[(h - m)] = h`.
- The Fourier coefficients of f_x are
  `c_k(x) = e^{-2 pi i k alpha} G_k(sigma)`, with

      G_k(sigma) = ∫ DI_sigma[h](t) e_{-k}(I(t)) kappa(t) dt,

  and the coefficients of Phi are `Phi_k = 2 pi i k G_k` (integrate by
  parts, or use (G')_k = 2 pi i k G_k). These are explicit smooth
  functions of sigma, with no inversion of diffeomorphisms.

**Lemma 0 (coefficient membership) [P].** For each k, G_k, Phi_k, c_k and
their real and imaginary parts lie in A, with r-th Fréchet derivative
bounded by `C_r (1+|k|)^{r+1} (1+||sigma||_{H^3})^{N_r} exp(C_r ||sigma||_{H^3})`,
uniformly in alpha.

*Proof.* Each sigma-derivative of kappa = e^sigma/∫e^sigma is a finite
combination of products of e^sigma, inverse powers of ∫e^sigma and the
increments, bounded by `C exp(C||sigma||_∞)`; the same holds for I and
DI[h]. A derivative falling on `e_{-k}(I)` brings down a factor
`2 pi k` times derivatives of I. Since `||sigma||_∞ <= C ||sigma||_{H^3}`,
every r-th derivative of the integrand is bounded by
`C_r (1+|k|)^r exp(C_r ||sigma||_{H^3})`, and integration over t keeps the
bound. The factor e^{-2 pi i k alpha} is bounded with alpha-derivatives of
size (2 pi k)^j, a constant for each fixed k. Multiplying by 2 pi i k for
Phi_k gives the stated power. ∎

The k-growth of these coefficientwise bounds is harmless: no K-uniform
estimate below is run through coefficient sums; all K-uniform estimates
use the composed representation of u_K.

**Paired VP truncation.**

- V_K is the Fourier multiplier with symbol rho_K(q): it equals 1 for
  |q| <= K, equals (2K - |q|)/K for K < |q| < 2K, and equals 0 beyond.
- rho_K is even, so the cos_k and sin_k modes always carry the same
  weight. It is Lipschitz with constant 1/K, satisfies |rho_K| <= 1, and
  fixes constants.
- The truncated field is `Y_K = sum_{|k|<2K} rho_K(k) c_k X_{e_k}`, a
  finite combination which is real because `c_{-k} = conj(c_k)` and
  rho_K is even; in the real basis it is
  `sum_{k>=0} rho_K(k) (a_k X_{cos_k} + b_k X_{sin_k})` with a_k, b_k in A.
- By translation covariance, its angular part is
  `(V_K f_x)(alpha) = (V_K G)(0)`, which is independent of alpha. Its
  shape part is

      u_K = P_0 [ (V_K Phi) o I ] = sum_{|k|<2K} rho_K(k) Phi_k(sigma) P_0[e_k o I],

  also independent of alpha. The second expression follows from the
  first because V_K Phi is the trigonometric polynomial
  `sum rho_K(k) Phi_k e_k`; equivalently it is
  `sum_k rho_K(k) G_k(sigma) P_0[e_k' o I]`, the derivative factor coming
  from `v_{e_k} = P_0[e_k'(chi)]` and `e_k'(chi) = 2 pi i k e_k(alpha) e_k(I)`,
  whose angular phase cancels against the one in c_k. As K → ∞,
  `P_0[Phi o I] = P_0[h - m] = h`.

**Lemma A (quantitative composition and inversion) [P].** There are
absolute N_0, C_0 such that, with
`C_0(sigma) := C_0 (1+||sigma||_{H^3})^{N_0} exp(C_0 ||sigma||_{H^3})`,
for every sigma in H^3_0:

(a) I is an H^4 diffeomorphism of S^1 (degree one on the circle of
circumference one), `||I||_{H^4} + ||I^{-1}||_{H^4} <= C_0(sigma)`, and
`inf I' = inf kappa >= C_0(sigma)^{-1}`.

(b) For integer 0 <= s <= 4 and u in H^s,
`||u o I||_{H^s} + ||u o I^{-1}||_{H^s} <= C_0(sigma) ||u||_{H^s}`.

*Proof.* (a) H^3 is an algebra, so `||e^sigma||_{H^3} <= e^{C||sigma||_{H^3}}`
by the power series, `∫ e^sigma >= e^{-||sigma||_∞}`, and
`inf kappa >= e^{-2||sigma||_∞}`; I' = kappa gives the I bound. For
I^{-1}: `(I^{-1})' = 1/(kappa o I^{-1})`, and higher derivatives follow by
the chain rule from derivatives of kappa composed with I^{-1}; each
composition is estimated in L^2 by the change of variables
`∫ |v(I^{-1}(z))|^2 dz = ∫ |v|^2 kappa <= ||kappa||_∞ ||v||^2`, and
reciprocals cost powers of (inf kappa)^{-1}.

(b) Faà di Bruno for s <= 4: each term of (u o I)^{(s)} is
`u^{(j)} o I` times a product of derivatives of I of total order s. In
every term, all factors except at most one (the highest u-derivative or
the highest I-derivative) are in L^∞ with norms controlled by
`||u||_{H^s}` and `||I||_{H^4}` (one-dimensional Sobolev embedding), and
the remaining factor is in L^2; compositions are estimated as in (a).
Same for I^{-1}. ∎

(Continuity of composition and inversion in H^s/D^s for s > 3/2 is also
in Inci–Kappeler–Topalov [E]; it is cited as context only, since the
quantitative bounds used below are (a)–(b), proved here, so no external
hypothesis is left unverified on this point.)

**Lemma B (Gaussian Sobolev membership) [P, using E].** Each of the maps
u_K and the constant map h belongs to D^{1,2}(mu_2; H_2), with Malliavin
derivative the restriction to H_2 of the Fréchet derivative; likewise
every a in A belongs to D^{1,2}(mu_2) (for each fixed alpha, with bounds
uniform in alpha).

*Proof.* u_K is a finite sum of products of A-functions (Lemma 0) with
the maps `sigma → P_0[e_k o I_sigma]`, which are smooth from H^3 to H^4
with derivatives bounded as in Lemma A. So u_K: H^3 → H_2 is Fréchet-C^1
with `||u_K(sigma)||_{H_2} + ||Du_K(sigma)||_{HS(H_2,H_2)} <= G_K(sigma)`
of poly-exp type (HS-summability from Lemma 3 and the basis decay in
Proposition 5). Take a smooth orthonormal basis (eps_m) of H_2 and the
scalar coordinates `u_m = <u_K, eps_m>_{H_2}`. The finite-dimensional
functions `u_m o Pi_n` (Pi_n the projection onto the first n modes) are
C^1 with poly-exp growth of value and gradient, hence lie in the
finite-dimensional Gaussian Sobolev space W^{1,2}(gamma_n) [E: Bogachev,
§5.2], i.e. in D^{1,2}; `u_m o Pi_n → u_m` in L^2(mu_2) by continuity and
domination (`||Pi_n sigma||_{H^3} <= ||sigma||_{H^3}` keeps one dominating
function), and `sup_n E||D(u_m o Pi_n)||^2 < ∞` likewise. Nualart's
closability lemma [E: Lemma 1.2.3] gives u_m in D^{1,2} with the Fréchet
derivative as Malliavin derivative. The vector statement follows since
`D^{1,2}(H_2) = { sum u_m eps_m : sum_m (E u_m^2 + E ||D u_m||^2) < ∞ }`,
and that sum is E||u_K||^2 + E||Du_K||^2_{HS} < ∞ by Fernique [E]. The
same argument applies to a in A. ∎

**Lemma 2 (divergence identification) [P from D, E].** For each K,

    div Y_K = -delta_2(u_K)      nu_2-a.e.,

where delta_2 is the Skorokhod divergence of mu_2 (acting on the shape
variable of the smeared chirality; both sides are independent of alpha).

*Proof.*

- u_K in D^{1,2}(H_2) is Lemma B, so delta_2(u_K) is defined and lies in
  L^2 [E: Nualart Prop. 1.3.1].
- Test against a bounded smooth cylinder a. On the left, (27) with c = 2
  applied to the real and imaginary parts of c_k a (in A, Lemma 0), plus
  the Leibniz rule `X_{e_k}(c_k a) = (X_{e_k} c_k) a + c_k X_{e_k} a`,
  gives

      ∫ div(Y_K) a dnu_2 = sum_k rho_K(k) [ ∫ X_{e_k}(c_k a) - ∫ c_k X_{e_k} a - (1/2) ∫ c_k a A_{e_k} ] = - ∫ Y_K a dnu_2.

- On the right: `Y_K a = (V_K G)(0) ∂_alpha a + <Da, u_K>_{H_2}`. The
  angular term integrates to zero in alpha for each sigma, because
  (V_K G)(0) is alpha-independent and a is periodic. Gaussian duality
  [E: Nualart §1.3] gives `∫ <Da, u_K> dmu_2 = ∫ a delta_2(u_K) dmu_2`
  for each alpha. So the right side also equals `-∫ Y_K a dnu_2`.
- Bounded smooth cylinders determine L^1(nu_2) classes (full support,
  RGR-T3 [D]), so the two L^1 functions agree a.e. ∎

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

In particular, K-uniformly,

    ||[M_D, V_K] w||_{H^4} <= (48 ||D||_{A^1} + 7^4 ||D||_{H^4}) (||w||_{H^3} + ||w||_{A^0}).   (3.1)

The sharp symbol has jumps instead: `|rho(q) - rho(p)| = 1` on every
crossing pair. So its first case carries `<q>^4` rather than `<q>^3`
(see Section 4).

**Lemma 4 (gradient formula) [P].** Write `D~_e := (D_e I) o I^-1`, where
`D_e I(t) = ∫_0^t e kappa - I(t) ∫ e kappa`. For e in H_2,

    D_e u_K = P_0 [ ( [M_{D~_e}, V_K] Phi' ) o I ].

*Proof.*

- Differentiating Phi = (h - m) o I^{-1} in the direction e gives
  `D_e Phi = -Phi' D~_e - D_e m` (chain rule for I^{-1}, with
  `D_e(I^{-1}) = -D~_e/(kappa o I^{-1})` and
  `(h-m)' o I^{-1} = Phi' · (kappa o I^{-1})`); the constant -D_e m is
  annihilated later by P_0.
- Differentiating u_K = P_0[(V_K Phi) o I]:

      D_e u_K = P_0[ ((V_K Phi)' o I) · D_e I + (V_K D_e Phi) o I ]
              = P_0[ ( (V_K Phi') D~_e - V_K(Phi' D~_e) ) o I ]
              = P_0[ ( [M_{D~_e}, V_K] Phi' ) o I ],

  using `D_e I = D~_e o I`, that V_K commutes with d/dy, and that V_K
  fixes the constant D_e m. ∎

**Uniform domination.** Fix the smooth h. With C_0(sigma) from Lemma A
there are poly-exp C(sigma) (depending on h) such that, for the smooth
orthonormal basis eps_n of H_2 (so `||eps_n||_{H^s} = sqrt(2) <n>^{s-4}`):

- `||Phi||_{H^4} + ||w||_{H^3} + ||w||_{A^0} <= C(sigma)`, where
  w := Phi' (Lemma A(b) applied to h - m and I^{-1}; `||w||_{A^0} <= C ||w||_{H^1}`);
- `||D~_{eps_n}||_{A^1} <= C ||D~_{eps_n}||_{H^2} <= C(sigma) <n>^{-3}`
  and `||D~_{eps_n}||_{H^4} <= C(sigma) ||eps_n kappa||_{H^3} <= C(sigma) <n>^{-1}`
  (Lemma A(b) at s = 2, 4; H^3 and H^1 are algebras);
- hence, by Lemma 4 and (3.1),

      sup_K [ ||u_K||_{H_2} + ||Du_K||_{HS} ] <= G(sigma) := C_h · C_0(sigma)^{p},

  a single K-independent dominating function with `E[G^2] < ∞` by
  Fernique [E]. Every K-limit below is dominated by G.

**Proposition 5 (D^{1,2} convergence) [P, using E].** As K → ∞,
u_K → h in D^{1,2}(mu_2; H_2).

*Proof.* Fix sigma in H^3_0.

- (i) `||u_K - h||_{H^4} = ||P_0[((V_K - I)Phi) o I]||_{H^4}
  <= C_0(sigma) ||(I - V_K) Phi||_{H^4} → 0`, since Phi in H^4
  (Lemma A) and |1 - rho_K| <= 1 is supported in |q| > K.
- (ii) By Lemmas 3 and 4 and Lemma A(b) at s = 4,

      sum_n ||D_{eps_n} u_K||^2_{H_2} <= C(sigma)^2 sum_n [ 48^2 ||D~_{eps_n}||^2_{A^1} ||1_{K/2<|q|<5K/2} w||^2_{H^3}
                                                     + 7^8 ||P_{>K/2} D~_{eps_n}||^2_{H^4} ||w||^2_{A^0} ].

  The first bracket is `<= C(sigma)^2 <n>^{-6} · o_K(1)` (tail of the
  convergent H^3 series of w); the second is `<= C(sigma)^2 <n>^{-2}`
  termwise and each tail `||P_{>K/2} D~_{eps_n}||_{H^4} → 0`. The sum in
  n is dominated by `sum <n>^{-2}`, so dominated convergence in n gives
  `||Du_K||_{HS} → 0` pointwise in sigma.
- Both quantities are dominated K-uniformly by G(sigma) with E[G^2] < ∞,
  so dominated convergence in sigma gives
  `E||u_K - h||^2_{H_2} + E||Du_K||^2_{HS} → 0`. The limit h is
  deterministic, so Dh = 0, and u_K - h ∈ D^{1,2}(H_2) by Lemma B. ∎

**Proof of Theorem E.**

Fix the chirality the smearing acts on (say +); write delta_2 for the
Skorokhod divergence in sigma_+.

- By (2.1), Lemma 1 and Lemma 2, for every psi in A and every K,

      0 = ∫ F [ -Y_K psi + psi delta_2(u_K) ] dnu_2.

- **L^1(nu_2) limit of the complete expression.** Three pieces:
  1. Angular: `(V_K G)(0) ∂_alpha psi → G(0) ∂_alpha psi = 0`, since
     `|(V_K G)(0) - G(0)| <= ||P_{>K} G||_{A^0} → 0` pointwise with
     poly-exp domination, G(0) = 0, and ∂_alpha psi has all moments.
  2. Gradient: `E|<D psi, u_K - h>| <= (E ||D psi||^2)^{1/2} (E ||u_K - h||^2)^{1/2} → 0`
     by Proposition 5 (Dpsi in L^2 by Lemma B for A).
  3. Divergence: by the Skorokhod isometry bound [E: Nualart Prop. 1.3.1]
     applied to u_K - h ∈ D^{1,2}(H_2),

         E|delta_2(u_K) - W_2(h)|^2 <= E||u_K - h||^2_{H_2} + E||D u_K||^2_{HS} → 0,

     and `E|psi (delta_2(u_K) - W_2(h))| <= ||psi||_{L^2} ||...||_{L^2} → 0`.

  Hence `-Y_K psi + psi delta_2(u_K) → -∂_h psi + psi W_2(h) = delta_2(psi h)`
  in L^1(nu_2), the last identity by the divergence product rule
  [E: Nualart Prop. 1.3.5; hypotheses: psi ∈ D^{1,2} with all moments
  (Lemma B), h deterministic in H_2]. **This is the actual L^1(nu_2)
  convergence of the complete divergence expression** (2.1), not of the
  vector fields alone.
- Since F is bounded, `∫ F delta_2(psi h) dnu_2 = 0` for all psi in A and
  all smooth h in H^3_0.
- **Chaos endgame.** Fix a smooth orthonormal basis (eps_j) of H_2 and
  take the countable test family
  `psi = a(alpha_+) b(sigma_+) psi'(x_-)`, where a runs over a countable
  dense set of smooth functions of alpha_+, psi' over a countable
  L^2-dense subset of bounded smooth cylinders of x_-, b over finite
  Hermite products `prod_i He_{k_i}(W_2(eps_i))`, and h over the eps_j.
  Each such psi is in A (Gaussian linear functionals and their
  polynomials have poly-exp bounds). The recurrence
  `delta_2(He_k(W(eps_j)) prod_{i≠j} He_{k_i}(W(eps_i)) · eps_j)
   = He_{k+1}(W(eps_j)) prod_{i≠j} He_{k_i}(W(eps_i))`
  [E: Hermite recurrence plus the product rule] shows that the values
  `delta_2(b eps_j)` span a dense subspace of the orthogonal complement
  of the constants in L^2(mu_2) [E: Wiener chaos completeness,
  Nualart Thm 1.1.1].
- Therefore, for nu-a.e. (alpha_+, x_-), the function
  `F(alpha_+, ·, x_-)` is mu_2-orthogonal to every element of that dense
  subspace, hence mu_2-a.e. equal to a constant: F does not depend on
  sigma_+ almost everywhere.
- Rotations (r_beta, id) lie in G_c, act by alpha_+ → alpha_+ + beta with
  sigma unchanged and cocycle 1; invariance over rational beta plus
  continuity in measure (the Fubini argument of Theorem U(b)) removes the
  alpha_+ dependence.
- Repeating both steps on the - factor leaves F a.e. constant. ∎

**Why this is not the excluded route.**

- No finite Cameron–Martin shift is ever claimed to be a gauge
  transformation; Theorem O still holds.
- Orbit density is not used.
- The abelian 0–1 law is not transferred through orbits. The only
  Gaussian input is the L^2 completeness of Wiener chaos, applied to a
  weak identity derived from finitely many genuine generator identities
  followed by a proved L^1 limit.
- Invariance of F enters only through Lemma 1's exact identity; no
  invariance of the measure nu_2 is used anywhere (Lemma 1 uses
  quasi-invariance).

## 4. Sharp Dirichlet sums: the precise defect

Condition (L) in `UNIQUENESS.md` uses sharp partial sums |k| <= K. The
proof above goes through for sharp sums **if and only if** u_K^sharp → h
in a topology that makes delta_2 converge. In the sharp case:

- The gradient is `D_e u_K^sharp = P_0[([M_{D~_e}, P_K] Phi') o I]`.
- Every crossing pair (|q| <= K < |p|) contributes `d_l w_q` with weight
  `<p>^4 ~ K^4` and **no 1/K gain**.
- So `||D u_K^sharp||_{HS}` is controlled by `K^4 |w_q|` for |q| ≈ K.
- Since w = Phi' lies in H^{7/2-} almost surely (mu_2 samples are
  H^{7/2-}, so I^{-1} is H^{9/2-}) and not better, this is O(1) at best
  and is not small.
- The exact witness confirms this on the borderline profile w_q = |q|^-4:
  - the sharp band H^4-norm squared is 1.55, 1.14, 1.03, 1.009, 1.002 for
    K = 8, ..., 128, which stays at or above 1;
  - the VP norm decays: 4.6e-2, 7.5e-3, 1.4e-3, 3.0e-4, 6.8e-5.

Status (unchanged by the audit):

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

## 5. Theorem S repaired: domain, contract, upgrade, and the sector reading

**Domain and topology.**

- A_0 is the span of polynomials in finitely many S(f), with f real
  smooth. It is contained in MPQ's domain A, and it is invariant under
  every S(f) and every D(h) with smooth h (for smooth h,
  beta_h = S(C^{-1}h) with C^{-1}h smooth). It is **not** invariant under
  J(f): J(f_0) v is in H_red (Fernique moments) but not in A_0.
- `D_red = A_0 ⊕ A_0` inside `H_red = H_+ ⊕ H_-` (MPQ (8)), with the
  **norm topology of H_red**.
- The form eta is defined on `Psi ⊗_alg D_red`, positive and Hermitian.

**Gauge invariance hypothesis (slotwise).** eta is required to be
G_c-invariant **separately in each gauge slot**:
`eta(U_g F)[G] = eta(F)[G]` for all g (by Hermiticity this implies the
same in the second slot). Separate continuity in the Psi slot for (19) is
required as before.

**Remark N [P, using D].** Slotwise invariance cannot be weakened to
invariance under the simultaneous action
`eta(U_g F)[U_g G] = eta(F)[G]`. The kinematic pairing
`eta_1(psi)[phi] = <psi, phi>_{L^2(nu_1)}` on Psi is positive, Hermitian,
separately continuous for (19), and simultaneously invariant (each U_g is
nu_1-unitary [D]); but it is not slotwise invariant and is not a multiple
of conj(ell) ⊗ ell. Indeed if it were slotwise invariant, Theorem U(d)
would give `eta_1 = c conj(ell) ⊗ ell`; testing on RGR's (22) sequence
(||psi_N||_{L^2(nu_1)} → 0, ell(psi_N) → 1) forces c = 0, contradicting
eta_1 ≠ 0. So simultaneous invariance alone does not imply factorization
through the invariant dual, and the class of simultaneously invariant
forms is strictly larger; it is not classified here.

**Compatibility contract (BC).**

- (BC1) **Boundedness.** The reduced form
  `B(v, u) := eta(1 ⊗ v)[1 ⊗ u]` satisfies `|B(v, u)| <= C ||v|| ||u||`
  on D_red. Equivalently B extends continuously to H_red, so `B = <·, T ·>`
  with T bounded, T >= 0.
- (BC2a) **S/D weak symmetry on D_red.**
  `B(O v, u) = B(v, O u)` for v, u in D_red and O = S(f), D(h), f, h real
  smooth. (Well defined: A_0 is O-invariant.)
- (BC2b) **J compatibility through the extension.** For one smooth
  f_0 > 0 and all v, u in D_red,
  `<J(f_0) v, T u> = <v, T J(f_0) u>`.
  This is stated with the bounded T from (BC1) and the Hilbert inner
  product, because J(f_0) maps D_red outside itself, so "B(Jv, u)" is
  undefined; the continuous extension is the meaning of the contract.
- (BC3) **Alternative to (BC2b) only.** T commutes with
  `P_+ = E_{J(f_0)}((0, ∞))` (a spectral projection of the declared
  J(f_0), which is positive on H_+ and negative on H_- by MPQ-T3), or
  with the sign operator epsilon. (BC2a) is retained under either
  variant.

An **unbounded** closed form in place of (BC1) is not admissible here;
the substitute contract would be direct invariance of the form domain and
values under the unitaries exp(-itO), which is a strictly stronger
assumption than weak symmetry and is not analyzed further.

**Lemma S1 (cores) [P, using E: Nelson].** A_0 consists of analytic
vectors for every S(f) and every D(h) (h smooth): on chaos of order <= N
both operators raise the order by at most one with
`||O^n p|| <= C^n ((N+n)!/N!)^{1/2} ||p||`. Hence each restriction is
essentially self-adjoint on A_0, and since its closure is a self-adjoint
restriction of MPQ's self-adjoint operator, A_0 is a core for each S(f)
and D(h). No core statement is made for J; none is needed under (BC2b)
or (BC3).

**Lemma S2 (upgrade) [P].** Let O be self-adjoint, C a core for O, T
bounded with `<O v, T u> = <v, T O u>` for v, u in C. Then T commutes with
every bounded Borel function of O.

*Proof.* For u in C, v → <Ov, Tu> is bounded on C, so Tu ∈ Dom(O|_C)^* =
Dom O and OTu = TOu. For u in Dom O take graph-norm approximants from the
core; closedness gives TO ⊂ OT, hence T(O-z)^{-1} = (O-z)^{-1}T for
nonreal z, and the spectral theorem finishes. ∎

Weak symmetry **without** (BC1), or without the core property, does not
give this; that is why boundedness and the cores are explicit contract
items, not derived properties.

**Lemma S3 (dense-shift upgrade) [P, using E].** Let tau in L^∞(mu) with
`tau(· + t h) = tau` mu-a.e. for all t and all smooth h. Then tau is
constant. *Proof.* For bounded tau the map h → tau(· + h) ∈ L^1(mu) is
continuous on H_C (approximate tau in L^1 by bounded continuous s; the
shifted measures mu_h have uniformly L^1-convergent densities for h in
compact sets of H_C). Smooth h are dense in H_C = H^4, so tau is
invariant under **all** Cameron–Martin shifts, and the Gaussian 0–1 law
[E] makes it constant. ∎

**Theorem S' [P, using Theorem E, (BC), E].** Under slotwise invariance,
(BC1), (BC2a) and ((BC2b) or (BC3)), every admissible form is

    eta(psi ⊗ v)[phi ⊗ u] = conj(ell(psi)) ell(phi) <v, (c_+ P_+ + c_- P_-) u>,    c_± >= 0,

the zero form (c_± = 0) included.

*Proof.*

- **Gauge slot.** For fixed v, u, the form (psi, phi) → eta(psi⊗v)[phi⊗u]
  is separately continuous, sesquilinear and slotwise invariant. Two
  applications of Theorem U(a) with Theorem E give
  `conj(ell(psi)) ell(phi) B(v, u)`.
- **Within and across sectors, diagonal part.** The S(f) act identically
  on the two summands; their bounded Borel functions generate the
  diagonal algebra `{M_phi ⊕ M_phi : phi in L^∞(mu)}` (smooth f separate
  points of H^3 and generate its Borel sigma-algebra). By Lemmas S1 and
  S2 (and (BC2a)), T commutes with it, so T is a 2x2 matrix of
  multiplication operators `(M_{tau_{ij}})`, since L^∞(mu) is maximal
  abelian in B(L^2(mu)) [E].
- Commutation with U_h(t) = exp(-itD(h)) (Lemmas S1, S2 again) conjugates
  each M_{tau_{ij}} to M_{tau_{ij}(· + t h)} (the Cameron–Martin weight
  cancels in the conjugation), so each tau_{ij} is invariant under all
  smooth shifts, hence constant by Lemma S3. Thus
  `T = [[a, x], [conj(x), b]]` with scalars a, b >= 0.
- **Across sectors.** Under (BC2b): testing v = (0, v_-), u = (u_+, 0)
  with J = M_j ⊕ M_{-j}, j = ∫ f_0 e^s > 0, gives
  `-conj(x) <v_-, M_j u_+> = conj(x) <v_-, M_j u_+>` for all
  v_-, u_+ in A_0; choosing v_- = u_+ = 1 gives `<1, j·1> = E[j] > 0`,
  so x = 0. Under (BC3): T commutes with P_+, so x = 0 directly. Set
  c_+ = a, c_- = b. ∎

**Trichotomy (zero weights included) [P].** Let B = <·, (c_+P_+ + c_-P_-)·>
as above, and let H_B be the completion of D_red/ker B.

- **c_+, c_- > 0.** `W = sqrt(c_+)P_+ + sqrt(c_-)P_-` is bounded,
  **invertible**, commutes with every declared observable, and
  `W^*W = T`, so `B(v, u) = <Wv, Wu>`. The map v → Wv extends to a
  unitary from H_B onto H_red intertwining the declared algebra: the
  theory with weights (c_+, c_-) is unitarily equivalent to RGR's (29),
  and **equal** to it iff c_+ = c_- = 1 (RGR's isometric normalization).
  The equivalence claim **requires both weights strictly positive**; W is
  not invertible otherwise.
- **c_- = 0 < c_+** (symmetrically c_+ = 0 < c_-). ker B ⊇ D_red ∩ H_-,
  and H_B ≅ H_+ carrying only the + sector representation: the quotient
  theory retains the supported sector alone (MPQ's "one sign treated
  separately"). It is **not** unitarily equivalent to the two-sector
  (29): every H_- state has zero physical norm.
- **c_+ = c_- = 0.** B = 0 and eta = 0 (or eta is the zero multiple of
  the gauge factor): the physical quotient is {0}. The zero form
  satisfies all hypotheses and is part of the classification.

**Corrected sector interpretation.** This replaces the reading in
`UNIQUENESS.md` §3, which said the ratio is "not observable through any
declared MPQ observable". That reading was wrong.

- **Populations are measured.** P_+ is a spectral projection of the
  declared J(f_0). For a fixed reduced vector `v = v_+ + v_-` (and
  c_+ + c_- > 0, v not B-null), the predicted population of the + sector
  is

      c_+ ||v_+||^2 / (c_+ ||v_+||^2 + c_- ||v_-||^2),

  so it depends on c_+/c_-. Measuring it needs no sign-mixing observable.
  Witness: 1/2 at equal weights and 1/5 at (1, 4) for the same vector.
- **Relative phase is not measured.** Every declared observable commutes
  with P_±, so `<v, O v>_B` is unchanged under `v_- → e^{i theta} v_-`.
  This is the superselection statement, and it is the only thing that
  sign preservation alone implies.
- **Strictly positive weights are a unitary rescaling, not new physics.**
  For c_± > 0, W above commutes with every declared observable, so the
  theories (c_+, c_-) and (1, 1) are unitarily equivalent as
  representations of the declared algebra: the ratio changes only which
  reduced vector represents a given physical state, i.e. the
  normalization of the map R in (29). A vanishing weight is **not** a
  rescaling: it changes the theory to the one-sector quotient.
- RGR's isometric normalization ("the generalized state F → <v, RF> has
  norm ||v||_{H_red}") is exactly the extra condition c_+ = c_- = 1
  (witness: it passes at equal weights and fails at (1, 4), and it
  excludes the degenerate cases, which are not isometric).
- So the full physical inner product is **unique** under
  (BC1) + (BC2a) + (BC2b or BC3) + isometric normalization, and unique
  up to the trichotomy above without the normalization.

## 6. Invariant measures, and what is classified

**Proposition M (ordinary invariant measures) [P, using D].**

(a) nu_2 is quasi-invariant under G_c with continuous positive cocycle
R^{(2)} [D: (11), (18), (21)], and it is **not invariant**: if it were,
every R^{(2)}_g would equal 1 a.e., hence Delta_g = 0 a.e. for every g;
differentiating `t → ∫ Delta_{g_f(t)} a dnu_2` at t = 0 (dominated by
RGR-T4's derivative estimates [D]) gives `∫ A_f a dnu_2 = 0` for all
a in A, so A_f = 0 a.e., contradicting `∫ A_f^2 dnu_2 > 0` for
nonconstant f [D: (27) block].

(b) There is **no nonzero finite G_c-invariant measure absolutely
continuous with respect to nu_2**. In particular nu_2 is not equivalent
to any invariant probability measure.

(c) Had one existed, it would have been unique: for two invariant
probabilities m_i = F_i nu_2, the ratio F_1/F_2 is G_c-invariant on
{F_2 > 0} (the cocycles cancel), and Theorem E applies.

*Proof of (b).* Let m = F nu_2 be invariant, finite, nonzero. The set
{F > 0} is invariant mod nu_2-null because the Radon–Nikodym cocycle is
strictly positive, so by Theorem E it is conull, and m can be normalized
to an invariant probability. Theorem U's Step 1 provides
g_n = (gamma_n^{-1}, id) with: (i) `inf_{x in B_R} ||shape(T_{g_n}^{-1}x)||_{H^3} =: c_n → ∞`
(the escape estimate), and (ii) `T_{g_n}^{-1}(B_R) ⊂ {||sigma_+||_{H^3} <= S_n}`
for some finite S_n [D2: T_g^{±1} map balls to bounded sets]. So
`T_{g_n}^{-1}(B_R)` lies in the shell `{c_n <= ||sigma_+|| <= S_n}`.
Choose n_1 < n_2 < ... with c_{n_{k+1}} > S_{n_k}: the shells are
pairwise disjoint, while invariance gives
`m(T_{g_{n_k}}^{-1} B_R) = m(B_R)` for every k. Finiteness forces
m(B_R) = 0 for every R, and the B_R exhaust X_c (angles are compact),
so m = 0. ∎

The invariant object in this theory is therefore the **functional** ell
on the cocycle-weighted action, RGR (20)–(21), not a measure: (21) is an
isometry statement for U_g, and Lemma 1 uses only quasi-invariance. The
sigma-finite infinite case is **not addressed**: the shell argument needs
finiteness.

**Classified** (given Theorem E and the ledger):

- the continuous G_c-invariant dual of RGR's Psi for the topology (19):
  it is C·ell;
- all separately continuous, positive, Hermitian forms on Psi invariant
  in each argument: they are c·conj(ell)⊗ell, with a one-dimensional
  quotient for c > 0;
- ordinary invariant measures as in Proposition M: none absolutely
  continuous and finite, nu_2 itself quasi-invariant but not invariant;
- with the spectator, all forms obeying slotwise invariance and (BC):
  Theorem S' with the trichotomy.

**Not classified:**

- functionals that are not continuous for (19), or on other test spaces;
- forms invariant only under the simultaneous action (Remark N shows the
  class is strictly larger: the kinematic pairing is in it);
- sigma-finite infinite invariant measures absolutely continuous with
  respect to nu_2;
- invariant measures singular to nu_2. (Within the Gaussian scale
  family, RGR (27) selects c = 2 as the only constraint-solving
  functional [D]; that is a statement about the constraints, not about
  invariant measures.);
- spectator forms violating (BC1) (unbounded, without unitary
  invariance), or declared algebras without J and sign;
- anything about nonlinear geometry observables.

## 7. Dependency ledger

**RGR-1 statements used [D]**, at d9c7020, file
`papers/native-residual-gauge-rigging/THEOREM.md`:

| Item | Used for | Audit here |
|---|---|---|
| (10)–(11) continuous cocycles, `sqrt(R^(1)) = R^(2)` | Lemma 1, Theorem U | agrees with REPORT.md (C5) |
| (18)+(21): R^(2) is the nu_2 pushforward density / L^1 isometry | Lemma 1, Prop. M(a) | read; only quasi-invariance used |
| (19) Psi and its topology | Theorem U, Remark N | read |
| (22) separating bump | Theorem U step 2, Remark N | cylinder bump in UC_b checked |
| (23) X_f, A_f, B_f; A and `B_f A ⊂ A`; `d/dt Delta|_0 = 2A_f` | Lemmas 1, 2, Prop. M(a) | c_k, Phi_k in A proved here (Lemma 0) |
| (26) Bochner identity in L^1(nu_2) for a in A | Lemma 1 | not re-derived |
| (27) at c = 2: `∫ X_f a = (1/2) ∫ A_f a`; `∫ A_f^2 > 0` | Lemma 2, Prop. M(a) | not re-derived |
| RGR-T4 derivative estimates along flows | Prop. M(a) differentiation | not re-derived |
| D2: T_g^{±1} map H^3 balls to bounded sets | Prop. M(b) | as in UNIQUENESS.md |
| CGD-T2 step 3 (no invariant kinematic vector) | Theorem U step 2 | as in UNIQUENESS.md |
| MPQ-T3 (8)–(12): H_red, S, J, D(h), sign contract; A domain | Section 5 | A_0 ⊂ A; J positivity per sector read |

**External theorems [E], hypotheses checked where used:**

- Nualart, *The Malliavin Calculus and Related Topics* (2nd ed.):
  Thm 1.1.1 (chaos completeness); Lemma 1.2.3 (closability limit:
  F_n → F in L^2, sup E||DF_n||^2 < ∞ ⟹ F ∈ D^{1,2}); Prop 1.3.1
  (D^{1,2}(H) ⊂ Dom delta with `E delta(u)^2 <= E||u||^2 + E||Du||^2_HS`);
  Prop 1.3.5 (`delta(psi h) = psi W(h) - <Dpsi, h>` for psi ∈ D^{1,2},
  h ∈ H, with the L^2 integrability supplied by A's moments).
- Bogachev, *Gaussian Measures*: finite-dimensional Sobolev membership of
  C^1 functions with square-integrable value and gradient (§5.2); the
  Cameron–Martin 0–1 law (invariance under all H-shifts ⟹ constant);
  Fernique's theorem. The dense-shift upgrade needed on top of the 0–1
  law is Lemma S3 [P].
- Hermite recurrence He_{k+1}(x) = x He_k(x) - k He_{k-1}(x) (standard).
- Nelson's analytic vector theorem (Lemma S1).
- L^∞(mu) is maximal abelian in B(L^2(mu)) (standard).
- Inci–Kappeler–Topalov, Mem. AMS 226 (2013): cited as context for H^s
  composition/inversion; the quantitative bounds actually used are
  Lemma A [P].

**Proved here [P]:** Lemmas 0, A, B, 1–4, Proposition 5, Theorem E,
Proposition M, Remark N, Lemmas S1–S3, Theorem S' with the trichotomy,
and the sector reading. Theorem U, Theorem O and their proofs are
unchanged in `UNIQUENESS.md`.

## 8. Finite exact witnesses (`ergodicity/`)

Native JS exact rational arithmetic, run with
`python3.12 ergodicity/run.py --check`. They check identities and specific
finite instances only; they certify no continuum statement, and in
particular not ergodicity.

- **W1.** Commutator Fourier formula
  `[M_D, S]_p = sum d_l w_q (rho(q) - rho(p))` for sharp and VP symbols,
  for two trigonometric D, two w, and K = 4, 6. Also: commutation with
  d/dy, the cancellation identity of Lemma 4, and the fixing of
  constants.
  - *Negative controls:* the reversed sign, and the identity without the
    commuting derivative.
- **W1b (audit, point 3).** Coefficient assembly of u_K at the exactly
  representable angle alpha = 1/2 and I = id: summing
  `rho_K(k) c_k · (e_k' factor k)` reproduces V_K(Phi) with
  Phi_k = k G_k (up to the common 2 pi i), i.e. the derivative factor in
  the corrected Lemma 2 expansion.
  - *Negative control:* the assembly without the derivative factor (the
    pre-audit displayed formula) does not.
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
- **W3b (audit, point 1).** Quasi-invariance versus invariance in the
  finite model: the weighted divergence functional has zero mean against
  every monomial (the invariant **functional**), while the plain flow
  does not preserve the Gaussian (`E[∂_1 psi] ≠ 0` for psi = x_1, x_1^3:
  Liouville fails, the **measure** is not invariant).
- **W4.** Sector model on C^2 ⊕ C^2. The declared sign-preserving
  algebra plus P_+ has a commutant of dimension 2; adding a sign-mixing
  observable makes it 1. Populations 1/2 against 1/5. W^*W = T, and W
  intertwines the declared algebra.
  - *Negative controls:* W does not intertwine the mixing observable,
    and unequal weights fail isometric normalization.
- **W5 (audit, point 2).** Zero-weight sectors: rank(T) is 4, 2, 0 for
  weights (1,4), (1,0), (0,0) — full theory, one-sector quotient, zero
  form.
  - *Negative control:* the zero-weight W = diag(1,1,0,0) is not
    invertible (rank 2), so no unitary equivalence claim survives there.
- **W6 (audit, point 4).** Slotwise versus simultaneous invariance in a
  Z_2 swap model: the identity Gram is invariant under simultaneous
  conjugation but not slotwise, and it has rank 2 while the factored form
  conj(ell)⊗ell has rank 1.
  - *Negative control:* "simultaneous invariance forces rank one" is
    rejected.

The results block of `CERTIFICATE.json` records
`continuum_proofs_machine_certified: false`.
