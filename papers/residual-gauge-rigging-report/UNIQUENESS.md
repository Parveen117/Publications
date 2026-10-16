# Residual rigging uniqueness for G_c: reduction to ergodicity and an orbit obstruction

16 October 2026. Separate branch `claude/residual-gauge-invariant-dual`.
This note builds on the shared research branch at d9c70208520f53f31d93bd94c3190f48b34acb4c,
whose packet `papers/native-residual-gauge-rigging` (RGR-1) is treated as a
**stated dependency**. RGR-1 is not reconstructed or modified here, and
neither is any other frozen packet. Nonlinear geometry observables and
curvature are out of scope; that is the other researcher's track.

> **Update, 16 October 2026 (same day, later).** `ERGODICITY.md` now proves
> ergodicity of G_c on (X_c, nu_2) (Theorem E) by a divergence closure with
> de la Vallee Poussin weights, so U(b),(d) apply unconditionally. Condition
> (L) below, stated with sharp partial sums, is superseded as a route and
> remains undecided as stated (`ERGODICITY.md` Section 4). Theorem S is
> superseded by Theorem S' there, with an explicit domain/boundedness
> contract (BC); and Section 3's claim that the sector ratio is "not
> observable through any declared MPQ observable" is **withdrawn**: sign
> projections do measure sector populations. See `ERGODICITY.md` Section 5.
> `ERGODICITY.md` was audited and revised later the same day (quasi-invariance
> versus invariant measures, zero-weight sectors, coefficient conventions,
> the (BC) contract and slotwise invariance); the change list is `AUDIT.md`.

Labels: **[P]** proved here; **[D]** taken from RGR-1, with the step audited
below; **[E]** external theorem; **[C]** open or conjectural.

## Result in one paragraph

Take any linear functional on RGR's test space Psi that is continuous for
RGR's own topology (19) and invariant under the actual residual group G_c.
Every such functional is integration against the doubled-covariance measure
nu_2 with a bounded G_c-invariant density F (Theorem U). So uniqueness of
`ell`, and of the positive rigging up to scale, holds **if and only if G_c
acts ergodically on (X_c, nu_2)**. That is the exact additional condition
asked for. If it fails, each nontrivial invariant set gives an explicit
additional channel, which is a counterexample to uniqueness.

Whether ergodicity holds is left open. The obvious route through the
Cameron–Martin (abelian) 0–1 law is *proved* to be unavailable. Almost
surely, no fixed nonconstant smooth Cameron–Martin direction is tangent to
a G_c-orbit, and this stays true after completing G_c to C^6. To move a
configuration along such a direction one needs a diffeomorphism that is
one derivative below the C^6 regularity at which quasi-invariance is
verified (Theorem O). A sufficient analytic condition that would close
the gap is stated precisely as (L).

With the MPQ spectator factor, even assuming ergodicity, the full physical
inner product is unique only **up to independent weights on the two sign
superselection sectors** H_+ and H_- (Theorem S). RGR-T5's choice of equal
weights is a normalization, not a consequence of the constraints.

---

## 1. Dependency ledger: RGR-1 facts used, and their audit

The following are used, with the same notation as RGR-1: X_c = X x X,
X = S^1 x H^3_0, nu_c = d alpha ⊗ mu_c, kinematic nu_1, T_g, U_g.

- **(D1) Pointwise cocycle identity, RGR-T2 (10)–(11).**
  `R_g^(1)(x)^(1/2) = R_g^(2)(x)` at every point, with continuous strictly
  positive versions. *Audit:* this is consistent with, and independently
  agrees with, my earlier determinant identity (`REPORT.md`, Theorem C.2,
  (C5)). RGR-T1 (4) and my (C5) state the same formula,
  `det_2(I + D u_g) = exp(-∫ log g')`, and were derived by different
  arguments. The continuous extension of the energy difference (9) and the
  limit (14)–(16) were read but not re-derived here. Only the
  *consequences* (D1)–(D3) are used.
- **(D2) Test space and continuity, RGR (19) and RGR-T3.** Psi is
  `UC_b ∩ L^2(nu_1) ∩ L^1(nu_2)` with seminorms p_R, ||.||_L2(nu1) and
  ||.||_L1(nu2). Each U_g acts continuously, with
  `p_R(U_g psi) = sup_{B_R} |R_g^(1)^(1/2) psi o T_g^-1|`. U_g is unitary
  in L^2(nu_1) and, by (21), an isometry in L^1(nu_2). *Audit:* the step
  "T_g^±1 take H^3 balls to bounded sets" holds because
  `||u_g||_{H^3} <= C_g(1 + ||sigma||_{H^3})` on balls.
- **(D3) Separating bump, RGR (22).** The empirical variance
  `V_N = N^-1 sum_{j<=N} xi_j^2` tends to c almost surely under mu_c.
  *Audit:* for a smooth bump beta, the cylinder function `beta(V_N)` is
  Lipschitz on each H^3 ball for fixed N, so it lies in UC_b. This is all
  that is used.
- **(D4) Full support** of nu_1 and nu_2 in H^3, so a continuous function
  is determined by its class (RGR-T3).

No gap was found in the steps above. RGR's continuum analysis outside
(D1)–(D4), namely T4's common core and T5's MPQ map, is not needed for
Theorem U and was not re-audited.

---

## 2. Classification of all continuous invariant functionals

**Theorem U [P, using D1–D4].** Let phi be a linear functional on Psi that
is continuous for the topology (19) and satisfies `phi(U_g psi) = phi(psi)`
for all g in G_c and psi in Psi.

(a) There is a unique F in L^∞(X_c, nu_2) with `F o T_g = F` nu_2-a.e. for
every g, such that

    phi(psi) = ∫ psi F dnu_2      for all psi in Psi.

Conversely, every such F defines a continuous invariant functional. Hence

    (Psi')^{G_c} ≅ L^∞(X_c, nu_2)^{G_c}.

(b) `(Psi')^{G_c} = C · ell` if and only if the action of G_c on
(X_c, nu_2) is ergodic. Ergodicity of the product action is equivalent to
ergodicity of Diff^∞_+(S^1) on each chiral factor (X, nu_2).

(c) If B is an invariant Borel set with 0 < nu_2(B) < 1, then
`ell_B(psi) = ∫_B psi dnu_2` is a continuous invariant functional not
proportional to `ell`. The form `|ell_B|^2 + |ell_{B^c}|^2` is a positive,
separately continuous, invariant Hermitian pairing whose physical quotient
is two-dimensional. That is an explicit additional channel.

(d) Assume ergodicity. Let eta be any separately continuous, Hermitian,
positive form on Psi that is G_c-invariant in each argument. Then
`eta(psi)[phi] = c · conj(ell(psi)) ell(phi)` for some c >= 0. The
physical quotient has dimension at most one, and exactly one when eta is
nonzero.

*Proof.*

**Step 1 (escape to infinity).** Let Psi_b be the set of elements of Psi
whose support is bounded in both shape H^3 norms.

- Psi_b is dense in Psi. Let theta be smooth, equal to 1 on [0, 1] and 0
  on [2, ∞). Set `psi_M = psi · theta(||sigma_+||/M) theta(||sigma_-||/M)`.
  The H^3 norm is Lipschitz, so psi_M is in Psi. Moreover psi_M = psi on
  B_R once M > R, and psi_M → psi in L^2(nu_1) and L^1(nu_2) by dominated
  convergence.
- Fix psi in Psi_b supported in the set where both norms are <= M, and fix
  R. Let `phi_n(y) = n^(-5/2) sin(2 pi n y)` and let gamma_n be the smooth
  diffeomorphism with `log gamma_n' = phi_n - log ∫ e^{phi_n}`.
- On B_R, sigma lies in a bounded subset of C^2, so `chi' = k` is bounded
  above and below, and chi'', chi''' are bounded. Hence

      ||(phi_n o chi)'''||_L2 >= c_R ||phi_n'''||_L2 - C_R (||phi_n''||_∞ + ||phi_n'||_∞)
                              >= c_R' n^(1/2) - C_R' n^(-1/2).

  (The first term follows by changing variables in `phi_n'''(chi) chi'^3`.)
- Therefore `inf_{x in B_R} ||sigma_+ + u_{gamma_n}(x)||_{H^3} → ∞`. Take
  `g = (gamma_n^-1, id)` with n large. Then `T_g^-1 (B_R)` is disjoint from
  supp psi, and so `p_R(U_g psi) = 0`.
- Continuity of phi gives
  `|phi(chi)| <= C (p_R(chi) + ||chi||_{L2(nu1)} + ||chi||_{L1(nu2)})`
  for some R and C. Apply it to `chi = U_g psi` and use invariance,
  unitarity and (21):

      |phi(psi)| = |phi(U_g psi)| <= C (||psi||_{L2(nu1)} + ||psi||_{L1(nu2)}).

- So phi is continuous on Psi_b for the norm
  `N = ||.||_{L2(nu1)} + ||.||_{L1(nu2)}`. Its N-continuous extension and
  phi are both continuous on Psi and agree on the dense subspace Psi_b, so
  they agree on Psi. Thus **every** continuous invariant functional is
  N-continuous.

**Step 2 (splitting along the two singular measures).** The map
`psi → ([psi]_{nu1}, [psi]_{nu2})` is an isometry from (Psi, N) into
L^2(nu_1) ⊕_1 L^1(nu_2).

- Its image is dense. Let a be a bounded continuous cylinder function and
  `theta_N = beta(V_N(sigma_+))` as in (D3). Then
  `a theta_N → (0, a)` and `a(1 - theta_N) → (a, 0)`, and such a are dense
  in both summands.
- By (D1), U_g acts diagonally: as the Koopman unitary
  `R^(1)^(1/2) · o T_g^-1` on the first summand, and as the L^1 isometry
  `R^(2) · o T_g^-1` on the second.
- Hence `phi = <v, ·>_{L2(nu1)} + ∫ · F dnu_2`, with v in L^2(nu_1) and F in
  L^∞(nu_2). By density, each summand is separately invariant.
- An invariant v is a G_c-invariant vector of the kinematic Hilbert space,
  so v = 0 (CGD-T2, step 3).
- Invariance of the second summand reads
  `∫ R^(2)_g psi(T_g^-1 x) F(x) dnu_2 = ∫ psi (F o T_g) dnu_2 = ∫ psi F dnu_2`
  for all psi. Hence `F o T_g = F` almost everywhere.
- Conversely, such an F is invariant by the same identity, and N-continuous
  because F is bounded. This proves (a).

**(b)** follows from (a). The functional `ell` corresponds to F = 1, and
ergodicity means exactly that the invariant F are constants. For the
product: if each factor is ergodic and F(x, y) is invariant, then for
almost every y, F(·, y) is invariant under the first factor and hence
constant (apply Fubini over a countable dense subgroup, using continuity in
measure). So F depends on y alone and is therefore constant. The converse
is immediate.

**(c)** The functional is N-continuous and invariant by (a). If it were
`c · ell`, then `∫ psi (1_B - c) dnu_2 = 0` for all psi in Psi. Psi is
dense in L^1(nu_2) by Step 2, so 1_B = c almost everywhere, which
contradicts 0 < nu_2(B) < 1. Positivity of the form is clear, and its two
coordinate functionals are linearly independent.

**(d)** For each psi, the functional `eta(psi)[·]` is continuous and
invariant, so it equals `a(psi) ell`. Hermiticity gives
`a(psi) ell(phi) = conj(a(phi) ell(psi))`. Set phi = 1, using ell(1) = 1:
then a(psi) = conj(a(1)) conj(ell(psi)), and a(1) = eta(1)[1] >= 0. ∎

*Remark.* Step 1 is what makes the classification complete. RGR's sup
seminorms on balls cannot carry an invariant functional, because the
residual group pushes every bounded set of shapes out of every ball.
This holds for arbitrary continuous functionals, not only for those given
by measures.

---

## 3. Spectator factor: what the constraints do and do not fix

**Theorem S [P, using E].** Assume ergodicity. Let eta be a separately
continuous positive Hermitian form on `Psi ⊗_alg D_red` that is invariant
under G_c in the gauge slot, and symmetric for every MPQ observable O in
the sense `eta(F)[(1 ⊗ O)G] = eta((1 ⊗ O)F)[G]` on the declared domains.
Then

    eta(psi ⊗ v)[phi ⊗ u] = conj(ell(psi)) ell(phi) · ( c_+ <P_+ v, P_+ u> + c_- <P_- v, P_- u> ),

with c_± >= 0 independent. Here P_± are the sign projections of
`H_red = H_+ ⊕ H_-`.

*Proof.*

- By Theorem U(d), applied slotwise with v and u fixed,
  `eta = conj(ell) ⊗ ell ⊗ B(v, u)`, where B is a positive form on D_red.
- Symmetry for O means B commutes with O. On each H_ε the MPQ algebra
  contains the multiplication operators S(f), whose spectral functions
  generate L^∞(mu), a maximal abelian algebra. It also contains the
  Cameron–Martin unitaries `exp(-it D(h))`.
- An operator commuting with both is multiplication by a function
  invariant under all Cameron–Martin shifts, which is constant by the
  Gaussian 0–1 law [E]. So each H_ε is irreducible.
- The sign operator commutes with the whole declared algebra (MPQ-T3), so
  B can only be a sum of scalars on the two sectors. ∎

So the physical inner product is fixed up to the ratio c_+/c_-. That ratio
is not observable through any declared MPQ observable, because all of them
preserve the sign. It does affect the weights of states that superpose the
two sectors. Fixing it needs one of:

- an additional declared sign-mixing observable;
- an explicit normalization convention (RGR's (29) uses c_+ = c_-);
- restricting to one sign sector, as MPQ allows.

This is the second "additional condition" for uniqueness.

---

## 4. Why abelian uniqueness does not transfer: an orbit obstruction

Per chirality, write `chi_sigma = alpha + I_sigma` and

    W_sigma = P_0 { phi o I_sigma : phi in C^∞(S^1) } ⊂ H.

**Theorem O [P, using IGQ-T2's path statistic].**

(a) The G_c-orbit of (alpha, sigma) is exactly `S^1 x (sigma + W_sigma)`,
and `W_{sigma'} = W_sigma` along the orbit. Each W_sigma is dense in H.

(b) Fix a nonconstant smooth h in H. For nu_2-almost every sigma, h is not
in W_sigma. More strongly, h is not in `P_0 { phi o I_sigma : phi in C^5 }`.
So the Cameron–Martin shift `sigma → sigma + t h`, with t ≠ 0, leaves the
orbit under the C^6 completion of G_c, almost surely.

(c) The unique diffeomorphism g with `g(alpha) = alpha` and
`T_g(alpha, sigma) = (alpha, sigma + h)` has
`log g' = h o I_sigma^-1 + const`. It lies in H^5 but almost surely not in
C^6. That is one derivative below the C^6 regularity at which CGD-T3 and
IGQ-T3 verify quasi-invariance.

*Proof.*

- **(a)** `T_g(alpha, sigma) = (g(alpha), sigma + P_0 log g'(chi))`, where
  `log g'` ranges over all smooth functions, up to a constant, as g ranges
  over Diff^∞_+. Rotations move alpha freely. Since
  `I_{sigma'} = g o I_sigma`, up to rotation, with g smooth,
  `C^∞ o I_{sigma'} = C^∞ o I_sigma`. Density holds because composition
  with a fixed H^4 diffeomorphism is an isomorphism of H^4, and C^∞ is
  dense in it.
- **(b)** Suppose `h = phi o I_sigma + c` with phi in C^5. Choose an
  interval J on which h' ≠ 0. Then phi' ≠ 0 on I_sigma(J), so on J
  `I_sigma = phi^-1 o (h - c)` is C^5. Hence `sigma = log I_sigma' + const`
  is C^4 on J, and sigma''' has zero quadratic variation on J. But the
  argument of IGQ-T2 ((11)–(12), with the variance scaled by c = 2 and
  dyadic meshes of J) gives, almost surely and simultaneously for all
  rational intervals J, a quadratic variation of sigma''' on J equal to a
  positive constant times |J|. Contradiction.
- **(c)** `P_0 log g'(chi) = h` forces `log g' = h o chi^-1 + const`. Since
  chi is in H^4, this is H^4, so g is H^5. If g were C^6, then
  log g' would be C^5, which (b) excludes almost surely.

Consequence: the Gaussian 0–1 law for Cameron–Martin shifts (the abelian
uniqueness of `REPORT.md`, Theorem C.1(c)) cannot be transferred to G_c by
inclusion of orbits, even after completing G_c to the largest group for
which quasi-invariance is currently verified. Any proof of ergodicity must
be genuinely nonlinear.

---

## 5. A sufficient condition that would close the gap

Let `f_{x,h} = (DI_sigma[h]) o chi^-1`. This is the field-dependent
gauge parameter realizing the shift h at x, and it lies in H^5. Write
`f_{x,h} = sum_k c_k(x) e_k` in the Fourier basis in y.

RGR's algebra A and generators `B_f = -X_f + A_f/2` from (23) are used.
For an invariant F in L^∞(nu_2), differentiating
`∫ F (V_t a - a) dnu_2 = 0` shows `∫ F B_f a dnu_2 = 0` for all a in A
[P, by the argument of RGR (26) in L^1(nu_2)].

**Condition (L) [C].** For every smooth h and every psi in A, the
coefficients c_k psi lie in A, and

    sum_{|k|<=K} B_{e_k}(c_k psi)  →  -∂_h psi + (1/2) s_h psi    in L^1(nu_2),

where `∫ ∂_h a dnu_2 = (1/2) ∫ s_h a dnu_2` defines the Gaussian score s_h.

**Proposition L [P].** (L) implies ergodicity, and hence the uniqueness in
Theorem U(b),(d).

*Proof.* (L) gives `∫ F (-∂_h psi + (1/2) s_h psi) dnu_2 = 0` for all psi
and all smooth h. So F is weakly invariant under every Cameron–Martin
direction. That makes F constant, by the Malliavin-derivative
characterization (DF = 0 implies F is constant) [E], and rotation
invariance removes the angular dependence. ∎

**Where (L) is delicate.** Expanding,

    sum_k B_{e_k}(c_k psi) = -Y_K psi + psi [ (1/2) sum_k c_k A_{e_k} - sum_k X_{e_k} c_k ],

where `Y_K = sum_k c_k X_{e_k} → (0, h)` in H, because f_{x,h} is in H^5.

- The first term converges.
- The bracket is a regularized Leibniz identity. Term by term,
  X_{e_k} c_k has an O(k) part, which cancels because ∫ f' = 0, and an
  O(1) part, of size about `∫ k_sigma h`. So the second sum is
  conditionally convergent at best, like the Fourier traces of RGR-T1,
  and it must cancel against the renormalized pairings A_{e_k}.

Proving (L) means proving that this cancellation holds in L^1(nu_2) under
the paired-sine/cosine summation. That is the precise remaining
convergence statement. Abelian uniqueness is never assumed in it.

---

## 6. Status

| Statement | Status |
|---|---|
| Every continuous G_c-invariant functional on Psi is ∫ · F dnu_2 with F invariant | [P] Theorem U(a) |
| Uniqueness of ell, and of the positive rigging up to scale, iff G_c is ergodic on (X_c, nu_2) | [P] U(b),(d) |
| Explicit extra channel from any nontrivial invariant set | [P] U(c) |
| With spectator: unique only up to sign-sector weights c_+, c_- | [P]+[E] Theorem S |
| Abelian 0–1 law does not transfer by orbit inclusion, even with a C^6 completion | [P] Theorem O |
| Ergodicity of Diff^∞_+(S^1) on (S^1 x H^3_0, nu_2) | **open** [C] |
| Regularized Leibniz condition (L), which implies ergodicity | [C], sufficient by Prop. L |

**Recommended next step:** prove (L). It is a single conditional-trace
computation of the same kind as RGR-T1 (5), for the operator
`v → D_sigma f_{x,h}[v]` on Diff directions. A finite exact witness of the
O(k) cancellation and of the diagonal value at sigma = 0 would be a useful
regression check before the continuum estimate. Alternatively, one could
look for a nontrivial invariant set, for example one built from
fine-scale statistics of sigma''' along the orbit-invariant smooth
structure W_sigma. That would yield a genuine extra channel by U(c).
