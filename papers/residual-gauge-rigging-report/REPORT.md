# Residual gauge rigging after central integer averaging

Independent report on the CGD-1 endpoint. 16 October 2026.

Repository `Parveen117/Publications`, branch
`codex/quantum-geometry-propagation-20261015`. Verified HEAD:
`fea9ef35e1357e8fa1afee24104472ccc2846bce` ("Derive central gauge rigging
and common-core source constraint algebra"), matching the stated start
commit. Nothing was committed, pushed or commented on the branch or PR #71.
All frozen packets are untouched.

Labels used throughout:

- **[P]** proved in this report (proof written here);
- **[E]** externally invoked theorem (cited; not reproved here);
- **[C]** conjecture or precisely stated open problem;
- **[W]** finite exact witness (regression evidence, not a proof).

---

## 0. Sources actually inspected and assumptions

Read in full at the verified HEAD:

- `papers/native-central-gauge-domain/THEOREM.md`, `SOURCE_MAP.md`,
  `README.md`, `SOURCE_PINS.json`, `run.py`, `verify.cjs` (CGD-1).
- `papers/native-inverse-gauge-quantization/THEOREM.md` (IGQ-1).
- `papers/native-material-profile-quantization/THEOREM.md`, sections 1 to 4
  (MPQ-1).
- The CGD and IGQ diffs of `RESEARCH_CLOSURES.md` and
  `QUANTUM_GRAVITY_RESEARCH.md` (commits `fea9ef3` and `fa98e69`), and the
  log of the 25 most recent commits on the branch.

Not read: the full text of Shimomura (2001), which CGD itself only cites
by its abstract, and the MRF packet beyond what IGQ/MPQ restate. No claim
below depends on reading them, except where it is marked [E].

Replays: `python3.12 papers/native-central-gauge-domain/run.py --check`
passes with certificate hash `9704439d…5094`, which equals
`EXPECTED.sha256`. The IGQ and MPQ `verify.cjs --check` runs also pass.
(Note: the container's default `python3` is 3.13; `run.py` rejects it by
design, so 3.12 was invoked explicitly.)

Assumptions inherited without change (from CGD/IGQ/MPQ):

- (A1) Per chirality, `H_one = L^2(S^1_alpha x H^3_0, d alpha d mu_0)`, with
  `mu_0` the mean-zero Gaussian of Fourier variances
  `lambda_n=(1+n^2)^(-4)` (cos and sin), Cameron–Martin space `H^4_0`.
  `H_c = H_one tensor H_one` (CGD (5),(12)).
- (A2) Configuration `chi = alpha + I_sigma`, `I_sigma(x)=int_0^x k`,
  `k = e^sigma / int e^sigma` (IGQ (14)).
- (A3) Residual action of `g` in `Diff^inf_+(S^1)` on one factor:
  `Phi_g(alpha,sigma) = (g(alpha), sigma + u_g)`, with
  `u_g = P_0 log g'(alpha + I_sigma)` (IGQ (16)). Unitaries
  `U^c_g psi = R_g^(1/2) psi o Phi_g^(-1)` (IGQ (19), descended by CGD-T1).
- (A4) IGQ-T3: Ramer hypotheses hold, `D_sigma u_g` is Hilbert–Schmidt on
  `H^4_0` (IGQ (17),(18)).

---

## A. Audit of the zero-invariant-dual theorem (CGD-T2)

**Verdict: CGD-T2 is correct as stated, and its scope statement is
accurate.** Checked steps:

1. The tail inverse (10)–(11) is right. `||h_n||^2 = O(n^(-2r-3)) p_(r+2)(v)^2`
   multiplied by `(1+n^2)^r` is `O(n^-3)`, which is summable. Together with
   the right inverse `E`, this gives the topological isomorphism
   `s(Z;H_c)/(I-S)s(Z;H_c) = H_c`.
2. Factorisation is right. A Z_c-invariant continuous functional kills
   `(I-S)s`, so it factors through `P_Z`. It is then Hilbert-continuous on
   `H_c` and given by a Riesz vector `v`.
3. There are no G_c-invariant vectors in `H_c`. The argument through the
   finite invariant measure `|v|^2 m` and its `alpha = chi(0)` marginal is
   sound, because `(g o chi)(0) = g(alpha)` is equivariant. It works with
   spectators.
4. The extension to separately continuous forms is right.

What the proof really uses is that, on `Phi_Z`, every invariant functional
is automatically bounded in the `H_c` norm. Stated that way, it is a
special case of a trivial but useful general fact.

**Lemma A.1 [P] (Hilbert-bounded no-go).** Let `Psi` be any
G_c-invariant subspace of `H_c`, with any topology. If an invariant linear
functional on `Psi` is bounded in the `H_c` norm, it is zero.

*Proof.* It extends by continuity to `closure(Psi)`. Its Riesz vector is
then G_c-invariant, and by step 3 that vector is 0. ∎

So a nonzero residual physical state must be a functional that is
*unbounded* on `H_c`. A stronger obstruction holds for the abelian
reduction introduced in B: such a functional cannot even be integration
against a measure absolutely continuous with respect to the kinematic
measure.

**Lemma A.2 [P] (no invariant absolutely continuous measure under
Cameron–Martin shifts).** Let `e_j` be the Fourier orthonormal basis of the
Cameron–Martin space, and `y_j` the coordinates (i.i.d. N(0,1) under
`mu_0`). Suppose `rho >= 0` is Borel and `lambda = rho mu_0` satisfies
`lambda(B + t e_j) = lambda(B)` for every Borel set B, every real t and
every j. Then `rho = 0` almost everywhere. The same holds on
`S^1 x H^3_0`, with rotations included, and for complex measures of
σ-finite variation.

*Proof.* Translating gives
`rho(sigma + t e_j) exp(-t y_j - t^2/2) = rho(sigma)` a.e., for each t. So
`rho exp(-y_j^2/2)` is invariant under translations in the `e_j`
direction. By Fubini, with `mu_0 = N(0,1) x mu'`, it is a.e. a function
of the complementary coordinates. Iterating over `j <= n` gives

    rho = exp( (1/2) sum_(j<=n) y_j^2 ) kappa_n(Q_n sigma),

where `kappa_n` is measurable with respect to the coordinates `y_j`,
`j > n`. The set `{rho > 0} = {kappa_n > 0}` is then a tail event, so its
measure is 0 or 1 by Kolmogorov's 0–1 law. Suppose it is 1, and put
`f = log rho`, which is finite a.s. Then `f = (1/2) chi^2_n + G_n` with
`G_n` independent of `chi^2_n`. Lévy's concentration function satisfies
`Q(X+Y; 1) <= Q(X; 1)` for independent X and Y. Hence
`Q(f; 1) <= Q(chi^2_n/2; 1) = O(n^(-1/2)) -> 0`. But for a finite real
random variable, `Q(f; 1) > 0`, since countably many unit intervals cover
R. This is a contradiction. For the complex case, the total variation of
an invariant complex measure is itself invariant. ∎

Scope remarks (audit). (i) The zero dual in CGD-T2 is a statement about
the Hilbert residual factor. Raising the integer decay weights cannot
change it, as CGD itself says. (ii) CGD-T3's density step rests on
Shimomura (2001), and its proof was not accessed (CGD's own access limit).
Nothing below uses CGD-T3, except the formal identification in B.3.
(iii) By Lemma A.2, any nonzero invariant distribution must be *singular*
to the kinematic measure. That is exactly what the construction below
produces.

---

## B. Proposed space and operators: an abelianised residual reduction

### B.1 Classical equivalence

The chart `chi <-> (alpha, sigma)`, with `alpha = chi(0)` and
`sigma = log chi' - mean`, is a *global* diffeomorphism of the smooth
inverse-embedding space onto `S^1 x C^inf_0`. Its cotangent lift gives
momenta `(pi_alpha, pi_sigma)`. The residual moment map is

    Q(f) = <Pi, X_f> = f(alpha) pi_alpha + pi_sigma( P_0 f' o chi ).      (B1)

**Proposition B.1 [P].** On the smooth regular source, the constraint sets
`{T_eps[f] = 0 for all f}` and `{pi_alpha = 0, pi_sigma = 0}` coincide.
They generate the same ideal, and their gauge orbits on the constraint
surface coincide: each is the whole smooth configuration space, so the
reduced phase space of the embedding sector is a point.

*Proof.* `f -> X_f(chi)` is a linear bijection from `C^inf(S^1)` onto
`R x C^inf_0`. Given `(a, h)`, put
`delta chi = a + DI_sigma[h]`, which is periodic because
`DI_sigma[h](1) = 0`. Then `f = delta chi o chi^-1` gives
`delta alpha = a` and `delta sigma = P_0[(delta chi)'/chi'] = h`. Hence
`pi` is a combination of the `Q(f)` with smooth field-dependent `f`, and
conversely by (B1). Both groups act transitively on smooth
configurations. ∎

In the measure completion the two orbit relations differ. Diff orbits lie
inside `S^1 x (sigma + H^4_0)`, but they are not equal to it. This is
where quantum inequivalence could enter, and it is addressed in C.3.

### B.2 Abelian quantum constraints and the stronger test topology

Per chirality, let `A = S^1 x C^inf_0` act on `H_one` by rotations
`alpha -> alpha + t` and Cameron–Martin Weyl shifts, as in MPQ (10):

    V_h psi(alpha,sigma) = exp(-beta_h(sigma)/2 - ||h||_C^2/4) psi(alpha, sigma + h).

The rotations are genuine elements of `G_c`. The shifts are the quantised
abelian constraints `pi_sigma`.

**Proposed residual test space.** Pass to white-noise coordinates
`x = C^(-1/2) sigma`. Under `mu_0`, x is mean-zero white noise on the
circle, and `C^(-1/2) = (1 - Delta/4pi^2)^2` is an isomorphism of `C^inf_0`.
Let `(S)` be the Hida test space over the nuclear triple
`C^inf_0 ⊂ L^2_0 ⊂ D'_0`, built with `A_0 = 1 - Delta/4pi^2`, and set

    Psi_A = [ C^inf(S^1_alpha) ⊗̂ (S) ]^{⊗̂ 2}   (two chiralities),
    Phi_Z^Hida = J^-1 s(Z; Psi_A)               (central factor, as in CGD (6)).

`(S)` embeds continuously and densely in `L^2(mu_0)`, and its topology is
strictly stronger than the Hilbert topology. `Phi_Z^Hida` is invariant
under Z_c, under rotations and under `V_h` for smooth h. For the shifts
this uses [E]: on `(S)`, translations by elements of `C^inf_0` and
multiplication by `exp<., k>` are continuous (Kuo, *White Noise
Distribution Theory*, 1996, ch. 4–5; Obata, LNM 1577, 1994). In white-noise
coordinates, `V_h = e^{-|k|^2/4} M_{exp(-<.,k>/2)} T_k` with
`k = C^(-1/2) h`.

### B.3 Relation to the actual residual generators (formal, identified exactly)

Write `u = P_0 f' o chi` and `D_j` for the Cameron–Martin generators.
Differentiating (A3) gives `Q(f) = -i(X_f + (1/2) div_m X_f)`. Since
`-iD_j = d_j - beta_j/2`, this equals

    Q(f) = f(alpha) pi_alpha  +  (1/2) sum_j { u_j(sigma), D_j }   (+ the α analogue),   (B2)

where `{ , }` is the symmetrised product. The real Diff constraints are
therefore the *Weyl-symmetrised field-dependent combinations* of the
abelian constraints. The only ordering correction, `(i/2)` times the
regularised divergence of `X_f` with respect to the formal flat measure
`D alpha D sigma`, is computed in Theorem C.2. The sum over j in (B2) is
formal; (B2) is used only for motivation.

---

## C. A nonzero invariant distribution

### C.1 The abelian construction (decisive construction)

Let `nu = d alpha ⊗ gamma_(2C)`, where `gamma_(2C)` is the Gaussian
measure with **doubled covariance** `2C`. Define `Omega` on `(S)` by its
S-transform:

    S Omega(eta) = <<Omega, :exp<., eta>:>> = exp( (1/2) <eta, eta> ),
    equivalently  Omega(psi) = ∫ psi dgamma_(2C)  (continuous versions).    (C1)

Formally, `Omega` is the flat half-density `(D sigma / dmu_0)^(1/2)`,
renormalised by `2^(-1/2)` per mode. The geometric mean of "Lebesgue" and
`gamma_C` is `gamma_(2C)`.

**Theorem C.1 [P, using E].**
(a) `Omega` is a nonzero element of `(S)*`, with `Omega(1) = 1`. It is a
positive distribution and is singular to `mu_0`; in particular it is not
in `L^2` (Feldman–Hájek).
(b) `Omega(V_h psi) = Omega(psi)` for every `psi` in `(S)` and every
smooth `h`.
(c) **Uniqueness:** every `Phi` in `(S)*` with `Phi o V_h = Phi` for all
smooth h equals `Phi(1) Omega`.
(d) On `Psi_A`, the distributions invariant under rotations and shifts in
both chiralities are exactly `C · (1_alpha ⊗ Omega)^{⊗2}`. Composing with
`P_Z` gives the unique, up to scale, nonzero invariant functional on
`Phi_Z^Hida`:

    r(f) = (Omega_+ ⊗ Omega_-)(P_Z J f).

*Proof.* Write `psi_eta = :exp<x,eta>: = exp(<x,eta> - <eta,eta>/2)`.
Direct substitution gives

    V_k psi_eta = exp( <k,eta>/2 - <k,k>/8 ) psi_(eta - k/2).                   (C2)

So invariance of Phi is equivalent to
`F(eta) = e^{<k,eta>/2 - <k,k>/8} F(eta - k/2)` for `F = S Phi`.

- (b) `e^{<eta,eta>/2}` satisfies this identically. The exponents cancel,
  as witnessed exactly in [W2]. Exponential vectors span a dense subspace
  of `(S)`, and `V_k` and `Omega` are continuous [E], so invariance holds
  on all of `(S)`.
- (c) Take `k = 2 eta` with eta real. Then
  `F(eta) = e^{<eta,eta>/2} F(0)`. F is entire, so the identity extends to
  complex eta. The S-transform is injective [E], so `Phi = F(0) Omega`.
- (a) Existence follows from the Potthoff–Streit characterisation [E],
  since `|e^{<eta,eta>/2}| <= e^{|eta|_0^2/2}`. The Gaussian
  `N(0, 2 I)` integrates `psi_eta` to `e^{<eta,eta>/2}`, which gives (C1)
  and positivity. Singularity holds because `C^(-1/2)(2C)C^(-1/2) - I = I`
  is not Hilbert–Schmidt.
- (d) Expand in α-Fourier modes. A rotation-invariant component must be
  the zero mode. Apply (c) to each chirality, using nuclearity of the
  completed tensor products. Z_c is handled by CGD-T2's coinvariants. ∎

This answers B and C for the abelianised reduction. The test topology is
explicit and strictly stronger than Hilbert, and it supports a nonzero
continuous invariant distribution. That distribution is unique and is
not reached by Lemma A.1 or A.2: it is unbounded on `H_c` and singular to
the kinematic measure.

### C.2 Invariance under the actual residual group: what is proved

**Theorem C.2 [P] (regularised Liouville identity).** Fix g in
`Diff^inf_+(S^1)` and a point `(alpha, sigma)` in `S^1 x H^3_0`.

(i) The Fourier-cutoff trace of `D_sigma u_g`,
`tau_g = lim_N sum_(n<=N) [<c_n, D u c_n> + <s_n, D u s_n>]`, exists for
**every** `(alpha, sigma)`, and

    tau_g = -log g'(alpha) - kappa_g ,   kappa_g := -∫_(S^1) log g' >= 0.        (C3)

The infinitesimal version, for a vector field f, is
`tau(D_sigma X_f^sigma) = -f'(alpha)`. Hence the regularised flat
divergence vanishes identically:

    div_reg X_f = f'(alpha) + tau(D_sigma X_f^sigma) = 0.                      (C4)

(ii) `J_g := g'(alpha) det_2(I + D_sigma u_g) e^{tau_g}` is a strict
cocycle: `J_(gh)(x) = J_g(Phi_h x) J_h(x)` at every point.

(iii) `J_g ≡ 1`. Equivalently, for every g and every configuration,

    det_2( I + D_sigma u_g(alpha, sigma) ) = exp( -∫_(S^1) log g' ).            (C5)

So the formal flat measure `D alpha D sigma` is exactly left-invariant
under the residual group, with Carleman-regularised determinant.

*Proof.*

(i) Write `D_sigma u_g[h] = P_0[a · DI[h]]` with `a = (g''/g') o chi`. The
kernel of DI is `K(x,t) = k(t)(1_{t<x} - I(x))`. The real Fourier
projector kernel is `D_N(x-t) - 1`, with `D_N` the Dirichlet kernel. Then

    tau_N = ∫∫ a(x) K(x,t) [D_N(x-t) - 1] dt dx.

For the `I(x)` part, `S_N k -> k` in `L^2`, which gives
`∫ a I k - ∫ a I`. For the `1_{t<x}` part,
`∫_0^x k(t) D_N(x-t) dt = S_N(k 1_[0,x])(x)`. At the jump point x this
tends to `k(x)/2` (Dirichlet), with uniform bounds for bounded-variation
data, so dominated convergence gives `(1/2) ∫ a k - ∫ a I`. Hence
`tau = ∫ a k (1/2 - I) dx`. Since `a k = d/dx log g'(chi(x))`,
integration by parts gives the boundary term `-log g'(alpha)` (using
`chi(1) = chi(0) + 1`) plus `∫ log g'(chi) chi' dx = ∫_(S^1) log g'`. The
regularity used is `a` and `k` in C^1, which follows from `chi` in H^4
and `sigma` in H^3. The vector-field case is the derivative at t = 0.

(ii) `alpha' = g(alpha)` is independent of sigma, so the Jacobian is block
triangular. Its σ-block for `gh` is `(I+A)(I+B)`, with
`A = D u_g(Phi_h x)` and `B = D u_h(x)`. `AB` is trace class, so its
cutoff trace equals its trace, and the cutoff trace is linear where the
limits exist. Combining this with
`det_2((I+A)(I+B)) = det_2(I+A) det_2(I+B) e^{-tr AB}` [E: Simon, *Trace
Ideals*, Thm 9.2] gives the cocycle.

(iii) Fix a flow `g_t`. From (ii),
`d/dt J_(g_t)(x) = J_(g_t)(x) · d/ds|_0 J_(g_s)(Phi_(g_t) x)`. At s = 0,
`det_2(I + K_s) = 1 + O(||K_s||_HS^2)` and `kappa_(g_s) = O(s^2)`, and the
`g'` factor cancels the `-log g'` in (C3). So the derivative is zero and
`J_(g_t) ≡ 1`. Every element of `Diff^inf_+(S^1)` is a finite product of
time-one flows, because the subgroup they generate is normal and nontrivial
and the group is simple [E: Herman 1971; Mather; Thurston 1974]. The
cocycle then gives `J ≡ 1`, and (C5) follows from (C3). ∎

[W1] gives an exact check of the second-order term of (C5) at `sigma = 0`.
There `log det_2 = -(1/2) tr K_1^2`, and (C5) requires
`tr K_1^2 = -∫ psi^2`. With the native rational engine and a rational
trig polynomial psi, the truncated traces converge to `-106549/88200`,
which is exactly `-∫ psi^2`, with error times L approaching a constant
(about 6.6). This is evidence only.

**Proposition C.3 [P, using E] (matching Radon–Nikodym versions).**
`nu` is quasi-invariant under every `Phi_g`; IGQ-T3's Ramer check applies
verbatim, since the Cameron–Martin space is the same with norm scaled by
`2^(-1/2)`. Take the Ramer–Kusuoka densities
`Lambda = g'(alpha) |det_2| exp(-delta(u) - (1/2)|u|^2)` for `m` and for
`nu` [E: Ramer 1974; Kusuoka 1982; Üstünel–Zakai 2000, ch. 3]. Their
Skorokhod integrals are the limits `delta(P_N u) = L_N^(.) - tau_N`, with
`L_N^m = sum_(j<=N) u_j sigma_j / lambda_j` and `L_N^nu = L_N^m / 2`.

Because `tau_N -> tau_g` *everywhere* (C3), a single subsequence makes
`L_N^m` converge a.s. under both measures. So there is one Borel set of
full `m`- and `nu`-measure on which both densities are given by the same
Borel formula, and on it

    Lambda^nu_g = (Lambda^m_g)^(1/2) · J_g^(1/2) = (Lambda^m_g)^(1/2).          (C6)

Equation (C6) is exactly the statement that integration against `nu` is
formally invariant under the half-density unitaries `U^c_g`.

### C.3 The remaining gap

`nu` is singular to `m`, so `U^c_g psi` is an `m`-class while `Omega` needs
its values on an `m`-null set. Theorem C.2 and Proposition C.3 prove
invariance for the *canonical Borel versions*. They do not prove that
those versions are the ones the class `U^c_g psi` determines, and they do
not prove that `(S)` is `U^c`-invariant. `U^c_g` composes with
non-analytic smooth g, while `(S)` consists of functions with entire
extensions; the non-invariance of `(S)` was not proved either. A
formulation that avoids versions uses the absolutely continuous
approximants

    nu_N = gamma_(2C)|_(modes <= N) ⊗ gamma_C|_(modes > N) ⊗ d alpha,
    Omega_N(phi) = ∫ phi (d nu_N / dm) dm .

These are genuine functionals on `m`-classes, with `dnu_N/dm` in `L^q` for
`q < 2`. Also `Omega_N -> Omega` in `(S)*`, because
`S Omega_N(eta) = exp((1/2) sum_(j<=N) eta_j^2)` [E: the
Potthoff–Streit convergence theorem].

**Conjecture R [C] (residual transfer).** For every g in `G_c` and every
`psi` in `Psi_A`:

- (R1) `U^c_g psi` lies in `L^p(m)` for some `p > 2`;
- (R2) `lim_N Omega_N(U^c_g psi) = Omega(psi)`.

**Consequence of R [P, conditional on C].** Define `Omega~ = lim Omega_N`
on its natural domain. It is well defined on classes and contains
`Psi_G = span U^c(G_c) Psi_A`, which is G_c-invariant and dense. On it,
`Omega~(sum_i U_(g_i) psi_i) = sum_i Omega(psi_i)`. `Omega~` is
G_c-invariant, nonzero (`Omega~(1) = 1`), and continuous for the final
topology induced by the maps `psi -> U_g psi` from `Psi_A`. That is a
nonzero continuous generalized physical state for the full residual
group.

---

## D. Positive physical pairing, null space and channel multiplicity

**Theorem D.1 [P] (abelianised reduction).** Let `eta` be a rigging map
on `Phi_Z^Hida`: antilinear, real, positive, and commuting with the
constraint group (Giulini–Marolf axioms). Then

    eta(f)[h] = c · conj(r(f)) r(h),   c > 0,

so the pairing is unique up to scale. Its null space is
`ker r = { f : (Omega_+ ⊗ Omega_-)(P_Z J f) = 0 }`, a closed hyperplane.

Explicitly, with chaos kernels `f_n`, the condition per chirality is
`sum_n (2n-1)!! Tr^n f_(2n) = 0`, where `Tr` contracts pairs of variables
along the diagonal. This follows from expanding `e^{<eta,eta>/2}`.

The physical space is `C`, with **exactly one surviving channel**. This
was derived from the uniqueness in Theorem C.1(c),(d), not assumed from
"pure gauge".

*Proof.* `eta(f)` is an invariant distribution, so by C.1(d)
`eta(f) = a(f) r`. Commutation with gauge makes a an invariant antilinear
functional. `V_h` preserves real functions, so `conj(r)` is the only one,
up to scale. Reality makes c real, positivity makes `c >= 0`, and
nontriviality makes `c > 0`. ∎

Positivity is not inherited from kinematic positivity here. The form is
the rank-one square `|r(f)|^2`, and its sign is fixed by the rigging
axioms. (It happens to be a positive distribution, `nu`, but nothing uses
that.)

**Corollary D.2 [P].** The hypotheses of IGQ-T6/CGD-T5 hold for the
abelianised reduction, with amplitude map
`r_g = r = ∫ (P_Z J f) d(nu_+ ⊗ nu_-)`. The product rigging
`conj(r f) r h <v,u>_MPQ` completes unitarily to `H_red`. The MPQ
observables, pulse and pointer/clock instrument act unchanged.

**For the actual residual group G_c [C].** If Conjecture R holds, the same
rank-one form built on `Omega~` is a nonzero positive pairing. Uniqueness
over `Psi_G`, and hence the channel count, is **not** established. A
G_c-invariant distribution that is not A-invariant on `Psi_A` is not
excluded. So for G_c itself we can only say that the multiplicity is at
least 1, conditionally on R.

---

## E. Exact obstruction and precise unresolved statements

1. **Proved obstructions.** Invariant functionals bounded in the Hilbert
   norm vanish on any invariant test space (A.1), which strictly contains
   CGD-T2. Under the abelian constraints, no invariant measure absolutely
   continuous with respect to the kinematic measure exists (A.2). Any
   physical state must therefore be singular to the kinematic measure.
   The construction meets this: `Omega` is the doubled-covariance
   Gaussian `nu`.
2. **Proved positive results.** For the classically equivalent abelian
   constraints, there is a unique nonzero invariant Hida distribution, a
   unique positive rigging, one channel, and recovery of MPQ (C.1, D.1,
   D.2). The flat measure is exactly left-invariant at the level of
   regularised determinants (C.2, unconditional and pointwise), and the
   Radon–Nikodym versions match (C.3).
3. **The single analytic gap** is Conjecture R. It is version consistency
   on the `m`-null support of `nu`, together with an integrability bound.
   Formally the anomaly in (B2) is `(i/2) div_reg X_f`, and it is
   identically zero by (C4).
4. **Open beyond R:** (a) the channel count for G_c alone; (b) whether a
   G_c-invariant σ-finite measure absolutely continuous with respect to
   `m` can exist (A.2 for G_c rather than A; conjectured not to); (c) all
   obligations CGD lists as separate (nonlinear geometry, other orderings,
   measure selection, 3+1).

| Item | Status |
|---|---|
| CGD-T2 correct, scope accurate | [P] audit |
| Hilbert-bounded no-go on any invariant test space | [P] A.1 |
| No invariant measure abs. continuous to m (abelian constraints) | [P] A.2 |
| Classical equivalence of abelianised constraints | [P] B.1 |
| Unique invariant Hida distribution `Omega = nu` | [P]+[E] C.1 |
| Regularised Liouville identity (C5) for Diff | [P]+[E: simplicity, det_2] C.2 |
| Matching Radon–Nikodym versions (C6) | [P]+[E: Ramer–Kusuoka] C.3 |
| Unique positive rigging, one channel, MPQ recovery (abelianised) | [P] D.1–D.2 |
| Nonzero continuous G_c-invariant state | [C] conditional on R |
| Channel count for G_c alone | open |

---

## Recommendation for the next research step

Prove Conjecture R. It is now a concrete Gaussian-analysis problem with
no group-theoretic content left, because (C5) is already proved.

1. **Integrability (R1).** Bound the Ramer density `R_g` in `L^(p/2)(m)`
   for some `p > 2`. Use IGQ's estimate `||D u_g||_HS <= C exp(c ||sigma||_(H^3))`
   with Fernique, and a Kusuoka-type `L^p` bound for `det_2 e^{-delta u}`.
2. **Convergence (R2).** Write `Omega_N(U^c_g psi)` as
   `∫ psi d nu_N^(g)`, where `nu_N^(g)` is the geometric mean of `m` and
   the `Phi_g`-transport of the partially flat measure
   `Leb_(<=N) ⊗ gamma_C|_(>N)`. Its density with respect to `nu_N` tends
   to `J_g = 1`: the tail Cameron–Martin norm and tail trace of `u_g`
   vanish and the head trace tends to `tau_g`. What remains is showing
   uniform integrability against `Psi_A`.
3. **Channel count for G_c.** Try to show that every G_c-invariant
   distribution on `Psi_G` is A-invariant on `Psi_A`. The plan is to
   integrate (B2) against field-dependent smearings `f = delta chi o chi^-1`
   that reproduce constant `h`. If that works, D.1 transfers verbatim and
   the MPQ equivalence becomes unconditional.
4. **Witnesses to add (separate branch only, after your go-ahead).** Exact
   native checks of the third-order term of (C5) at `sigma = 0`, and of the
   cutoff trace (C3) on a rational-coefficient family. These would remain
   [W] evidence.

---

## Appendix: local witnesses (not committed)

The script `rl_witness.cjs`, kept next to this report, uses only the
repository's `native_operator.cjs` (exact Q+iQ arithmetic). Python was not
used for any calculation. Output:

- [W1] Truncated `tr K_1^2` at L = 8, 16, 32 and 64, for
  `psi_hat = {1: 1/2 + i/3, 2: -2/5 + i/7, 3: 1/4}`. The limit is
  `-106549/88200`, which equals `-∫ psi^2`. The error times L is 7.14,
  6.83, 6.70 and 6.64.
- [W2] The exponent identity behind (C2) vanishes exactly in 50 random
  rational five-dimensional trials.

Run: `node rl_witness.cjs /path/to/Publications`
