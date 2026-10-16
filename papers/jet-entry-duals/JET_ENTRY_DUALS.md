# JED-1 — What can read a layer: duals, invariants and the finite-transport boundary of the jet-entry algebra

16 October 2026. Separate branch `claude/jet-entry-development`, based on
the chain tip `4b4a1c4` (POS-1), which contains TLM-1 (`babb4fb`) and
TET-1 (`9a06216`) unchanged. Nothing on the shared branch, PR #71, any
Codex branch or any frozen packet is modified. This develops the
operator tool of TLM-1 (jet-entry matrices) and its transport use in
TET-1; the quantum realization QRT-1 and the later packets are used only
as interfaces (Section 7).

Labels: **[P]** proved here; **[D]** statement of a repository packet used
as a dependency (ledger in Section 9); **[E]** external theorem; **[C]**
conjecture or open; **[W]** finite exact witness (Section 8). Finite
witnesses are regression checks only and certify no statement about
arbitrary depth, smooth sources or the continuum.

## 0. Result

TLM-1 stores an observer's reading together with its source-generated
layers in the algebra `A_N = M_m(J_N)`, `J_N = C[eps_s, eps_v]/(eps)^(N+1)`,
with the degree-zero map `Pi` as the observer. TLM-T2 shows the spectrum
of the one-matrix encoding cannot see the layers and that the jet algebra
has no faithful positive bounded *-representation. This packet answers the
general question those two facts raise: **which linear readouts of a
layered entry exist at all, and what do the layers determine.**

**Theorem (reading a layer) [P].** Let `d_N = (N+1)(N+2)/2 = dim J_N`.

1. *Representations.* Every irreducible representation of `A_N` factors
   through `Pi`; it is the ordinary representation on `C^m`. The layers
   act in any representation only through a nilpotent ideal
   (Theorem 1).
2. *States.* Every positive linear functional on `A_N` is
   `M -> tr(rho M_0)` with `rho >= 0`. No state, and no vector state of
   any *-representation, bounded or not, reads a layer (Theorem 2).
3. *Frames.* The linear functionals invariant under TLM's matched frame
   changes `M -> G^(-1) * M * G` (constant unitaries, unitary jets, or all
   invertible jets: the same answer) are exactly `lambda o Tr_J`, where
   `Tr_J(M) = sum_alpha tr(M_alpha) eps^alpha` is the `J_N`-valued trace.
   They form a `d_N`-dimensional space; the only positive one is the
   layer-blind trace state. For the TLM source, `Tr_J(H_mu) = 7 + mu eps_v`
   reads the deformation with its sign, which neither the curvature
   `5 mu^2/324` nor the spectrum of `Reg_1` does (Theorem 3).
4. *Charts.* Chart changes act on `J_N` by its automorphism group, the
   N-jet group of diffeomorphisms fixing the observer point. Under any
   chart group containing the linear volume-preserving changes `SL(2,R)`,
   the invariant linear functionals on `J_N` are `C * pi_N` alone; on
   `A_N`, `C * tr o Pi`. Under rotations alone, exactly the Laplacian
   readings `(Delta^k f)(p)`, `2k <= N`, survive. The full dual `J_N'` is
   the N-th order rigging of the observer by its derivative distributions
   `d^alpha delta_p`, and the observer is its only positive, only
   chart-invariant and only multiplicative element (Theorem 4).
5. *Transport.* The complete formal tower (all layers, `N = infinity`) of
   an admissible smooth TVSP source does not determine TET's finite
   transport, nor the curvature along the path: an explicit flat
   deformation of the TET source has the same tower as the undeformed
   source at the observer point, a different transport at every `t != 0`,
   and a different curvature there. For real-analytic sources the tower
   determines the transport along the whole path (Theorem 5).

So the layers are invisible to every state, every irreducible
representation and every chart-invariant linear functional, and are
visible exactly to the derivative distributions of the observer; the
frame group leaves the layer traces; and the tower is complete for
transport in the analytic class and incomplete in the smooth class at
every depth including infinite depth. Nonlinear invariants (the curvature
and its covariant derivatives, TLM-T4) are not classified here; that is
Codex's curvature-reconstruction track (Section 7).

## 1. Setting, from TLM-1 and TET-1 [D]

`J_N = C[eps_s, eps_v]/(eps_s, eps_v)^(N+1)` with the product (4) of
TLM-T1, basis `eps^alpha = eps_s^a eps_v^b`, `|alpha| = a + b <= N`, and
the coefficientwise involution `x -> x*` (both coordinates self-adjoint).
`A_N = M_m(J_N)` carries the `star` product and the dagger `M^dagger`
of TLM-T2, `(M^dagger)_ij = (M_ji)*`. The observer is
`Pi(M) = M_0`, a unital *-homomorphism `A_N -> M_m(C)` (TLM (7)).

`I_N = ker Pi = {M : M_0 = 0}` is a two-sided ideal with
`I_N^(N+1) = 0` and `I_N^dagger = I_N`, since every product raises the
total degree and the dagger preserves it. For `m = 1`, `I_N = J_N^+` is the
positive-degree ideal. Products act on column vectors from the left.

TLM's frame contract (10) is `M -> G_out^dagger * M * G_in` for unitary
jets `G`; matched frames (`G_out = G_in = G`) give the conjugation
`M -> G^(-1) * M * G`. TLM's chart changes act by the jet chain rule and
are separate from frame changes; Theorem 4 makes that action explicit for
scalar entries.

TET-T1 is the transport `U' = -A U`, `U(0) = I`, `A = (1/2) H^(-1) H'`
along a path, with the recurrence (2) producing the transport jet through
degree N from the source jets through degree N. TET-T6 shows that for
every finite N two admissible sources with equal jets through N can have
different transport at every positive endpoint.

The adapter's coefficient field is the native `Q + iota Q`. All proofs
below use only field operations and conjugation, so every statement holds
over `C` and verbatim over `Q + iota Q`.

## 2. Theorem 1 — characters and irreducible representations [P]

**Theorem 1.**

(a) The only unital algebra homomorphism `J_N -> C` is the observer
`pi_N`.

(b) Let `rho: A_N -> End(V)` be a unital representation on a nonzero
complex vector space with no invariant subspace other than `0` and `V`.
Then `rho(I_N) = 0`, so `rho` factors through `Pi`, and `rho` is
equivalent to the standard representation of `M_m(C)` on `C^m`. For
`m = 1`, `J_N` has exactly one irreducible representation, the
one-dimensional observer.

(c) For every finite-dimensional unital representation, the chain
`V ⊇ rho(I_N)V ⊇ rho(I_N)^2 V ⊇ ... ⊇ rho(I_N)^(N+1) V = 0` consists of
subrepresentations, and each quotient is a direct sum of copies of the
standard `C^m`. The layers act only through the nilpotent ideal
`rho(I_N)`.

**Proof.** (a) A homomorphism `chi` satisfies
`chi(eps_s)^(N+1) = chi(eps_s^(N+1)) = 0` in `C`, so `chi(eps_s) = 0`;
likewise `chi(eps_v) = 0`; hence `chi(x) = x_0`.

(b) Put `W = rho(I_N) V`, the span of all `rho(y) xi` with `y` in `I_N`.
It is invariant: `rho(x) rho(y) xi = rho(x * y) xi` and `x * y` lies in
the ideal. If `W = V` then
`V = rho(I_N) V = rho(I_N) rho(I_N) V = ... = rho(I_N^(N+1)) V = 0`,
a contradiction. So `W = 0`, that is `rho(I_N) = 0`, and `rho` is a
representation of `A_N / I_N`. The map `Pi` induces
`A_N / I_N = M_m(J_N) / M_m(J_N^+) ≅ M_m(J_N / J_N^+) = M_m(C)`.

Let `V` be an irreducible unital module over `M_m(C)` with matrix units
`e_ij`. Choose `xi != 0`; since `sum_i e_ii = 1` some `eta = e_ii xi` is
nonzero. The map `C^m -> V`, `(c_k) -> sum_k c_k e_ki eta`, intertwines
the actions because `e_jl e_ki = delta_lk e_ji`. It is injective: if
`sum_k c_k e_ki eta = 0`, applying `e_jj` gives `c_j e_ji eta = 0`, and
`e_ji eta != 0` because `e_ij e_ji eta = e_ii eta = eta != 0`. Its image
is a nonzero invariant subspace, hence all of `V`. So `V ≅ C^m`.

(c) Each `rho(I_N)^k V` is invariant as in (b), and the quotient
`rho(I_N)^k V / rho(I_N)^(k+1) V` is annihilated by `rho(I_N)`, hence a
finite-dimensional module over `M_m(C)`. Every such module is a direct
sum of copies of `C^m` [E: Wedderburn; elementary form: `V' ≅ C^m ⊗ e_11 V'`
by the map in (b)]. ∎

TLM's one-matrix encoding `Reg_N` (TLM (8)) is such a representation:
its chain is the total-degree filtration and its diagonal blocks `M_0`
are the standard quotients, which is TLM (9). Theorem 1 says every
representation has this shape; the "visibility" of layers in `Reg_N`
is the visibility of a nilpotent ideal, never of an irreducible.

## 3. Theorem 2 — the positive dual collapses to the observer [P]

A linear functional `phi` on `A_N` is **positive** if
`phi(M^dagger * M)` is real and nonnegative for every `M`.

**Lemma 2.1 (Cauchy–Schwarz) [E, proof included].** For positive `phi`,
`phi(B^dagger * A) = conj(phi(A^dagger * B))` and
`|phi(B^dagger * A)|^2 <= phi(A^dagger * A) phi(B^dagger * B)`. In
particular `|phi(A)|^2 <= phi(1) phi(A^dagger * A)`.

*Proof.* Expanding `phi((A + B)^dagger (A + B)) >= 0` and
`phi((A + iB)^dagger (A + iB)) >= 0` shows that
`x = phi(A^dagger B)`, `y = phi(B^dagger A)` satisfy `x + y` real and
`x - y` imaginary, so `y = conj(x)`. Thus `<A, B> = phi(B^dagger A)` is a
positive semidefinite Hermitian form, and the inequality is the usual
discriminant argument on `phi((A + t B)^dagger (A + t B)) >= 0`. ∎

**Theorem 2.** Every positive `phi` on `A_N` vanishes on `I_N`. Hence
`phi(M) = tr(rho_phi M_0)` for a unique positive semidefinite `rho_phi`,
and conversely every such functional is positive. The states of `A_N`
are exactly the density matrices on `C^m` read through `Pi`.

**Proof.** Let `a` be in `I_N` and put `a_0 = a`,
`a_(k+1) = a_k^dagger * a_k`. By induction `a_k` lies in `I_N^(2^k)`,
because `I_N` is dagger-closed. Choose `k` with `2^k >= N + 1`; then
`a_k = 0` and `phi(a_k) = 0`. Lemma 2.1 gives
`|phi(a_j)|^2 <= phi(1) phi(a_(j+1))` for every `j`, so downward
induction yields `phi(a_j) = 0` for all `j`, in particular `phi(a) = 0`.
Therefore `phi = phi_0 o Pi` with `phi_0(X) = phi(X)` on constant
matrices, which is positive on `M_m(C)`, hence of the form `tr(rho X)`
with `rho >= 0`; `rho` is unique because `tr(rho X)` determines `rho`.
Conversely `Pi` is a unital *-homomorphism, so
`tr(rho Pi(M^dagger M)) = tr(rho M_0^* M_0) >= 0`. ∎

**Corollary 2.2 [P].** Let `rho` be a *-representation of `A_N` by
operators on a pre-Hilbert space `D` (each `rho(a)` defined on all of
`D`, `rho(a) D ⊆ D`, `<rho(a) xi, eta> = <xi, rho(a^dagger) eta>`),
bounded or not. Then `rho(I_N) = 0`.

*Proof.* For `xi` in `D`, `phi_xi(a) = <xi, rho(a) xi>` is positive since
`phi_xi(a^dagger a) = ||rho(a) xi||^2`. For `a` in `I_N`, Theorem 2 gives
`||rho(a) xi||^2 = phi_xi(a^dagger a) = 0`. ∎

This strengthens TLM-T2's bounded statement (a self-adjoint nilpotent
bounded operator is zero) to unbounded *-representations and to
functionals that come from no representation at all. The GNS
representation of every state is a multiple of the standard one, in
agreement with Theorem 1.

What Theorem 2 does **not** say: that no quantum experiment can read a
layer. It says that a reading of a layer is never a single expectation
value of an element of the layered algebra. QRT-1 reads the first layer
of `H` as derivatives of a family of states (third central moments,
QRT (3)); those are derivative functionals of Theorem 4(d) applied to the
family, not states of `A_N`. See Section 7.

## 4. Theorem 3 — the frame-invariant linear dual is the layer trace [P]

Define the `J_N`-valued trace `Tr_J: A_N -> J_N`,
`Tr_J(M) = sum_i M_ii = sum_alpha tr(M_alpha) eps^alpha`.

**Theorem 3.**

(a) `Tr_J(X * Y) = Tr_J(Y * X)` for every `X` in `M_(m x n)(J_N)` and
`Y` in `M_(n x m)(J_N)`. Hence `Tr_J(G^(-1) * M * G) = Tr_J(M)` for every
invertible `G` in `A_N`.

(b) For a linear functional `ell` on `A_N` the following are
equivalent: (i) `ell` is invariant under conjugation by the signed
permutation matrices (constant, a finite group); (ii) under all constant
unitaries; (iii) under all unitary jets `G^dagger * G = 1`; (iv) under
all invertible jets; (v) `ell = lambda o Tr_J` for a unique linear
`lambda: J_N -> C`, i.e. `ell(M) = sum_alpha c_alpha tr(M_alpha)`.

(c) The invariant linear dual therefore has dimension `d_N`, spanned by
the layer traces `M -> tr(M_alpha)`. Its positive elements are the
multiples of `tr o Pi`: the unique invariant state is the normalized
trace state, which is layer-blind.

(d) Under the two-frame action `M -> G_out^dagger * M * G_in` with
independent frames, the only invariant linear functional is zero.

**Proof.** (a) Coefficientwise,
`(Tr_J(X*Y))_alpha = sum_i sum_k sum_(beta+gamma=alpha) (X_ik)_beta (Y_ki)_gamma`,
and the coefficient field is commutative, so the two orderings agree.
Then `Tr_J(G^(-1) M G) = Tr_J(M G G^(-1)) = Tr_J(M)`.

(b) (v) implies (iv) by (a), and (iv) ⇒ (iii) ⇒ (ii) ⇒ (i) by restriction.
Assume (i). For each `alpha` put `ell_alpha(X) = ell(X eps^alpha)` on
`M_m(C)`. With `D = diag(e_1, ..., e_m)`, `e_i = -1` and all other
entries `+1`, `D^(-1) e_ij D = -e_ij` for `j != i`, so
`ell_alpha(e_ij) = -ell_alpha(e_ij) = 0`. With the transposition
`P = P_(ij)`, `P^(-1) e_ii P = e_jj`, so `ell_alpha(e_ii) = ell_alpha(e_jj)`.
Hence `ell_alpha = c_alpha tr` and
`ell(M) = sum_alpha c_alpha tr(M_alpha) = lambda(Tr_J M)` with
`lambda(eps^alpha) = c_alpha`.

(c) `Tr_J` is onto `J_N` (`Tr_J(x e_11) = x`), so `lambda` is unique and
the dimension is `d_N`. If `ell` is also positive, Theorem 2 gives
`ell = tr(rho ·) o Pi`; invariance under the signed permutations forces
`rho` diagonal with equal entries, so `ell = c tr o Pi`, `c >= 0`.

(d) Take `G_out = -1`, `G_in = 1`: `ell(M) = ell(-M)`. ∎

**Remarks.**

1. The spectrum of `Reg_N(M)` is invariant under conjugation by all of
   `GL(m d_N, C)`, a far larger group than `Reg_N` of the frame group;
   that is why it sees only `M_0` (TLM (9)). `Tr_J` is the invariant that
   survives for the frame group and not for `GL(m d_N)`. The ordinary
   trace of the encoding is `tr Reg_N(M) = d_N tr(M_0)` [W]. The Jordan
   type of `Reg_N(M)` is a `GL(m d_N)`-invariant that detects whether a
   layer is nonzero (nilpotent rank) but not its value [W].
2. *TLM's example.* For the source (15) of TLM-T5 at depth one,
   `Tr_J(H_mu) = 7 + mu eps_v`: a frame-invariant **linear** reading of
   the first layer that recovers `mu` with its sign [W]. The curvature
   `kappa(p) = 5 mu^2/324` forgets the sign (TLM-T5), the characteristic
   polynomial of `Reg_1(H_mu)` is `(t^2 - 7t + 9)^3` for every `mu`
   (TLM (20)). The three readouts are consistent: they are invariants of
   three different groups.
3. *Metric records.* Under the congruence `H -> G^T * H * G` with
   `G = c I`, `ell(c^2 H) = ell(H)` forces `ell = 0`. So there is no
   nonzero linear invariant of a metric record under constant rescaling
   of the frame [W]; every frame-covariant readout of a metric record,
   `kappa` included (TLM (14)), is nonlinear. Under orthogonal frames,
   congruence is conjugation and Theorem 3 applies.
4. The theorem is about TLM's frame contract (10) with the chart fixed. It
   is not chart invariance; Theorem 4 treats that.

## 5. Theorem 4 — chart jets, the observer and the Laplacian tower [P]

For `p, q` in `J_N^+` let `theta_(p,q)` be the substitution
`eps_s -> p`, `eps_v -> q`, extended multiplicatively:
`theta(x) = sum_alpha x_alpha p^a q^b`. Entrywise, it acts on scalar
records in `A_N`. This is the action of an N-jet of a chart change fixing
the observer point on scalar entries; for tensor records such as `H` the
Jacobian jet acts in addition (TLM-T2, "the actual jet chain rule"),
which is not needed below except in Corollary 4.3.

**Theorem 4.**

(a) `theta_(p,q)` is a well-defined unital homomorphism of `J_N` fixing
`pi_N`, a *-homomorphism when `p, q` have real coefficients, and an
automorphism exactly when the degree-one parts of `p` and `q` are
linearly independent. These are all the unital automorphisms of `J_N`;
they form the N-jet group of local diffeomorphisms of `(R^2, 0)` (for
real coefficients) and act on `J_N` by the N-jet of the chart change.

(b) Let `Gamma` be any group of such substitutions containing the three
linear ones `(eps_s, eps_v) -> (2 eps_s, eps_v/2)`,
`(eps_s + eps_v, eps_v)`, `(eps_s, eps_v + eps_s)`; for example
`Gamma ⊇ SL(2,R)` (linear volume-preserving chart changes) or the full
jet group. Then the `Gamma`-invariant linear functionals on `J_N` are
`C * pi_N`, and the functionals on `A_N` invariant under `Gamma` acting
entrywise together with constant unitary frame changes are `C * tr o Pi`.

(c) Under the rotations `(eps_s, eps_v) -> (c eps_s - s eps_v, s eps_s + c eps_v)`,
`c^2 + s^2 = 1`, alone, the invariant linear functionals on `J_N` are
spanned by `pi_N` and the Laplacian readings

    ell_k(x) = sum_(j=0)^k C(k,j) (2j)! (2k-2j)! x_(2j, 2k-2j) = (Delta^k f)(p),
    1 <= k <= N/2,

where `x = j^N f`; the invariant dual has dimension `floor(N/2) + 1`.

(d) *(Rigging identification.)* Let `Germ_p` be the algebra of `C^inf`
germs at `p`. The entry map `j^N_p: Germ_p -> J_N` of TLM (3) is a
surjective unital *-homomorphism with kernel `m_p^(N+1)`, `m_p` the
germs vanishing at `p` [E: Hadamard's lemma]. Consequently

    J_N' = { ell in Germ_p' : ell(m_p^(N+1)) = 0 } = span { d^alpha delta_p : |alpha| <= N },

with the coefficient functional `x -> x_alpha` equal to
`d^alpha delta_p / alpha!`. In words: the dual of the layered entry
algebra is the N-th order rigging of the observer `delta_p` by its
derivative distributions. In this rigging, the observer is the only
positive element up to scale (Theorem 2), the only chart-invariant
element up to scale (part (b)), and the only multiplicative element
(Theorem 1(a)).

**Proof.** (a) The polynomial algebra maps to `J_N` by substitution; the
map descends to `J_N` iff `p^a q^b = 0` whenever `a + b >= N + 1`, which
holds iff `p, q` have zero constant term (if `p_0 != 0` then
`p^(N+1) != 0`). Then `theta(x)_0 = x_0`, and real coefficients commute
with the coefficientwise involution. `theta` preserves the degree
filtration and acts on the degree-`k` graded piece by the `k`-th
symmetric power of the linear part `L` of `(p, q)`; it is bijective iff
each symmetric power is, iff `L` is invertible (a block-triangular
argument). Conversely any unital automorphism sends `eps_s, eps_v` to
nilpotents, hence to `J_N^+`, so it is some `theta_(p,q)`, and
invertibility forces `L` invertible. Composition of substitutions is the
truncated composition of the maps `(p, q)`, which is the N-jet group.

(b) Let `ell` be invariant and write `ell_alpha = ell(eps^alpha)`.
Under `(2 eps_s, eps_v/2)`, `eps^alpha -> 2^(a-b) eps^alpha`, so
`(2^(a-b) - 1) ell_alpha = 0` and `ell_alpha = 0` unless `a = b`. For
`a = b >= 1` apply the shear `(eps_s + eps_v, eps_v)` to
`y = eps_s^(a+1) eps_v^(a-1)`:

    theta(y) = sum_(j=0)^(a+1) C(a+1, j) eps_s^(a+1-j) eps_v^(a-1+j).

Only the `j = 1` term has equal exponents, so
`ell(theta(y)) - ell(y) = (a+1) ell_(a,a)`, and invariance gives
`ell_(a,a) = 0`. Hence `ell = ell_0 pi_N`. For `A_N`, invariance under
the constant signed permutations gives `ell = sum_alpha c_alpha tr(M_alpha)`
by Theorem 3(b); the chart action commutes with the matrix part, so the
scalar functional `x -> sum_alpha c_alpha x_alpha` is `Gamma`-invariant,
hence `c_alpha = 0` for `alpha != 0`.

(c) On the degree-`n` graded piece use `z = eps_s + i eps_v`,
`zbar = eps_s - i eps_v`. A rotation by angle `t` sends
`z^j zbar^(n-j)` to `e^(i(n-2j)t) z^j zbar^(n-j)` (up to the sign
convention). An invariant functional must vanish on every eigenvector
whose eigenvalue is not identically `1`, i.e. on all `j != n/2`. So the
invariant functionals of degree `n` form a space of dimension at most one,
and zero for odd `n`. For even `n = 2k`, `ell_k` is invariant because
`Delta` commutes with rotations and `ell_k(j^N f) = (Delta^k f)(p)` by
`d_s^(2j) d_v^(2k-2j) f(p) = (2j)! (2k-2j)! x_(2j, 2k-2j)`; it is nonzero
on `(eps_s^2 + eps_v^2)^k`. Hence exactly one invariant per even degree
`<= N`, which is `floor(N/2) + 1` in all.

(d) Hadamard's lemma: a smooth `f` equals its degree-N Taylor polynomial
at `p` plus `sum_(|beta| = N+1) (x-p)^beta g_beta` with smooth `g_beta`;
so `ker j^N ⊆ m_p^(N+1)`, and the reverse inclusion holds because monomials
of degree `N+1` have zero N-jet. Surjectivity is clear on polynomials.
Thus `ell` in `Germ_p'` kills `m_p^(N+1)` iff `ell = lambda o j^N` for a
unique `lambda` in `J_N'`, and `lambda = sum_alpha lambda(eps^alpha)
(coefficient_alpha)`, where `coefficient_alpha(j^N f) = d^alpha f(p)/alpha!`. ∎

**Remarks.**

1. Which chart group is "physical" is a contract choice (TLM uses a
   normalized affine oriented chart; calibrated `S, V`). Theorem 4 gives
   the invariant linear dual for each natural choice: trivial group, all
   `d_N` layer coefficients; rotations, the Laplacian tower; any group with
   a non-compact linear part as in (b), the observer alone. The collapse
   needs a non-compact element (a scaling or a shear); a compact group
   leaves averaged readings. This is the finite-dimensional shadow of the
   mechanism in `ERGODICITY.md`: compact groups average, `Diff` does not.
2. No topology on `Germ_p` is used in (d); "distribution" means a linear
   functional of the displayed form. With Schwartz's theorem that
   distributions supported at a point are finite sums of derivatives of
   `delta_p` [E], the identification also holds in the topological sense.
   That theorem is not needed for anything here.
3. Part (d) is the precise sense in which a TLM entry `[C | j^N f]` is
   "the observer together with its rigging": the pairing of `f` with
   `span{d^alpha delta_p}`.

**Corollary 4.3 (metric records) [P].** Under linear chart changes
acting on a metric record by `H -> R^T (theta_R H) R`, the scaling
`R = c I` multiplies the layer `alpha` by `c^(2 + |alpha|)`, so every
chart-invariant linear functional of a metric record is zero. All
chart-invariant readouts of `H`, in particular `kappa` and its covariant
derivatives, are nonlinear in the layers. Their classification is the
invariant theory of metric jets and is not done here; see Section 7.

## 6. Theorem 5 — the complete tower versus finite transport

TET-T6 [D] exhibits, for each finite `N`, admissible sources with equal
jets through `N` and different transport at every positive endpoint.
Here the depth is infinite.

**Theorem 5.**

(a) [P] There is a one-parameter family of admissible TVSP sources,
`|eta| <= 1`, on TLM's chart `|s|, |v| < 1/20`, whose complete formal
towers at `p = 0` (all jets of `u`, `H`, the six responses `Lambda`, and
every covariant layer) coincide for all `eta`, while along the path
`gamma(t) = (t, 0)` TET's transport and the curvature satisfy, for every
`0 < t < 1/20` and `eta != 0`,

    U_eta(t) = diag(1, (1 + eta theta(t))^(-1/2)) != I = U_0(t),
    kappa_eta(t, 0) = eta^2 theta'(t)^2 / (4 (1 + eta theta(t))^2) != 0 = kappa_0(t, 0),

with `theta(s) = exp(-1/s^2)` for `s != 0`, `theta(0) = 0`.

(b) [P from E] If the supplied energy `u` is real-analytic on a
connected patch, then for every analytic path in the patch the complete
tower at `p` determines `H` along the path and hence `U(t)` for every `t`;
the coefficients are those of TET's recurrence (2), and `U` is analytic
on the whole parameter interval.

(c) Consequently the bridge from infinitesimal data (layers at the
observer point) to finite transport is complete in the real-analytic
class and incomplete in the smooth class at every finite depth (TET-T6)
and at infinite depth (a). Analyticity is a sufficient condition; the
exact class boundary is not determined here [C: the natural candidate is
a quasi-analytic class in the sense of Denjoy–Carleman [E], not proved].

**Proof.** (a) Put

    u_eta(s, v) = 10 s - 10 v + (s^2 + v^2)/2 + eta theta(s) v^2 / 2.

Then `T = u_s = 10 + s + eta theta'(s) v^2/2`, `P = -u_v = 10 - v - eta theta(s) v`,

    H_eta = [[1 + eta theta''(s) v^2/2,  eta theta'(s) v],
             [eta theta'(s) v,            1 + eta theta(s)]].

On `|s| <= 1/20` one has `0 <= theta <= e^(-400)` and
`|theta'|, |theta''| <= 10^9 e^(-400)` (both `s^(-k) e^(-1/s^2)` factors increase on the chart), so for `|eta| <= 1` the matrix
`H_eta` is within `10^(-160)` of the identity: `a > 0`, `Delta > 0`,
`T, P > 0`, and all six responses (1) of TLM are smooth on the chart,
their denominators `a, c, Delta` being nonzero. The sources are
admissible in the sense of TLM/TET.

Every derivative of `theta` vanishes at `0` [E: standard]. By the Leibniz
rule every partial derivative of `u_eta - u_0 = eta theta(s) v^2 / 2` at
`p = 0` vanishes, so `u_eta` and `u_0` have the same formal Taylor series
at `p`. The formal series of a smooth function of several smooth functions
depends only on their formal series (differentiate the composition; every
term is a finite sum of products of derivatives at `p`), so `H`, `Delta`,
`H^(-1)`, the six responses, the connection `(1/2) H^(-1) dH` and every
covariant layer (TLM (12)) have `eta`-independent formal towers at `p`.

Along `gamma(t) = (t, 0)`: `H_eta(t, 0) = diag(1, 1 + eta theta(t))`,
`A(t) = (1/2) H^(-1) H' = diag(0, eta theta'(t) / (2 (1 + eta theta(t))))`,
and `U_eta(t) = diag(1, (1 + eta theta(t))^(-1/2))` satisfies
`U' = -A U`, `U(0) = I` by direct differentiation; it is the transport by
uniqueness (TET-T1). For `t != 0`, `theta(t) > 0`, so `U_eta(t) != I`.

For the curvature use TLM (14) at `(t, 0)`:
`X_s = H^(-1) d_s H = diag(0, eta theta'/(1 + eta theta))`,
`X_v = H^(-1) d_v H = [[0, eta theta'], [eta theta'/(1 + eta theta), 0]]`
(the `v`-derivative of the `(1,1)` entry vanishes at `v = 0`). Then
`[X_s, X_v] = (eta theta')^2 [[0, -1/a], [1/a^2, 0]]` with
`a = 1 + eta theta`, `F_sv = -[X_s, X_v]/4`,
`(H F_sv)_12 = (eta theta')^2 / (4 a)` and `Delta = a`, giving the
displayed `kappa`. Since `theta'(t) = 2 t^(-3) e^(-1/t^2) != 0` for
`t != 0`, `kappa_eta(t, 0) != 0` for `eta != 0`, while `kappa_0 ≡ 0`
(constant `H_0`).

(b) The formal tower at `p` of an analytic `u` is its Taylor series,
which determines `u` on the connected patch [E: identity theorem]. Hence
`H`, `A` along the path and the unique solution `U` of `U' = -A U` are
determined. `A(t)` is real-analytic in `t`, so `U` is real-analytic on
the whole interval [E: Cauchy's theorem for linear ODEs with analytic
coefficients], and its Taylor coefficients at `0` are TET (2) by TET-T1.
∎

**Remarks.**

1. In (a) the small-loop return at the observer point is the same for all
   `eta` (the towers agree, so `F_sv(p)` and all its derivatives agree):
   the difference is invisible to every loop based at `p` of every
   finite order, visible to open-path transport and to loops based at
   `t != 0`. This refines TET-T5/T6: loop recovery at `p` reads
   `F_sv(p)`, which the tower contains; the tower does not contain
   `kappa(t, 0)` for `t != 0`.
2. Nothing in (a) is rational; no finite certificate evaluates it, and
   none could. The witnesses for Theorem 5 replay the TET-T6 family at
   rational endpoints and the analytic positive control of TET-T4.
3. The analogy with `BRIDGE.md` (Theorem K): there, infinitesimal
   constraint solutions give finite gauge invariance exactly on the
   closure of the dense core; here, infinitesimal layers give finite
   transport exactly on a class where the germ determines the function.
   In both, the class is the whole content of the bridge.

## 7. What this means for the chain, and the interfaces

| Item | Consequence here | Track |
|---|---|---|
| TLM-T2 (no positive bounded *-rep of `J_N`) | Strengthened: no positive functional and no *-representation, bounded or not, reads a layer (Theorem 2, Corollary 2.2); no irreducible representation does (Theorem 1) | Claude |
| TLM (9) (spectrum blind) | Explained: spectrum is a `GL(m d_N)`-invariant; the frame-group invariants are the layer traces `Tr_J`, which do read layers (Theorem 3) | Claude |
| TLM-T5 (`mu` vs `-mu`) | `Tr_J(H_mu) = 7 + mu eps_v` is a linear frame-invariant that keeps the sign the curvature forgets [W] | Claude |
| TLM chart changes | Made explicit as the jet group; chart-invariant linear dual is the observer alone, or the Laplacian tower under rotations (Theorem 4) | Claude |
| TET-T6 (finite depth) | Extended to infinite depth and to the curvature along the path; analytic class is complete (Theorem 5) | Claude (bridge) |
| QRT-T2: `d_i H_jk = -<Y_i Y_j Y_k>` | The quantum sector reads layers as derivatives of a **family** of states `theta -> omega_theta`, i.e. by `d^alpha delta_theta` applied to the family (Theorem 4(d) type functionals), never by one state (Theorem 2). This is consistent and is the mechanism by which QRT's `K` carries layer information: through the connection, not through a representation of `A_N` | Interface; QRT is Codex's |
| Curvature `kappa` and its covariant derivatives | Nonlinear invariants; their classification (invariant theory of metric jets) is curvature reconstruction | Codex; not developed |
| Physical path/time selection, path-independence rule (TET, POS-T4) | Theorem 5 says which class makes the tower complete; which paths are physical is not addressed | Codex; interface |
| POS-1 | Untouched | — |

Possible next step in Claude's track, if wanted: apply Theorem 4(d) to
QRT's state family. The functionals on the family algebra that are
second-order derivative distributions and are invariant under
reparametrization of `theta` should be classifiable in the same way;
the curvature-reading functional of QRT-T4 would then be characterized
by invariance rather than by construction. That is a statement about
uniqueness of a readout, which is this track; it would consume QRT only
as a dependency.

## 8. Finite exact witnesses (protocol JED1_JET_ENTRY_DUALS_V1)

`witnesses/verify.cjs` uses the unchanged native engine over
`Q + iota Q` through the unchanged TLM-1 jet adapter and TET-1 path
adapter: **112 exact checks and 46 negative controls**, at depths
`N = 1, 2, 3` unless stated. Run `python3.12 witnesses/run.py --check`
(Python 3.12 does I/O and pin checks only).

| Block | Checks | Negative controls |
|---|---|---|
| T1 | `eps^(N+1) = 0`, `eps != 0`; on the regular module the layer ideal acts nilpotently with index exactly `N+1`; its common kernel is the top degree (dimension `N+1`); every top-degree vector is a common eigenvector with the observer eigenvalue | a nonzero value of `eps` is not a character; a layer coefficient is not multiplicative; the layers do not act trivially; a degree-zero vector is not a common eigenvector |
| T2 | three point states `tr(rho M_0)` are nonnegative on squares and Hermitian; **every** sampled layer-reading functional (each single layer with four weights, and three mixed ones) is violated by an explicit `1 + t eps^alpha + r eps^(2 alpha)`; the observer is not | a reading functional nonnegative on all squares of the TLM fixture records still fails positivity: sampling certifies nothing |
| T3 | `Tr_J` cyclic (including `2 x 3` by `3 x 2`); sampled unitary jets `exp(eps X) exp(eps' Y) U_0` are unitary with `G^dagger = G^(-1)`; invertible jets inverted; `Tr_J` invariant under all of them while the record moves; signed permutations force the trace on `M_2(C)` (difference span rank 3); the difference span for `M_2(J_N)` has rank exactly `3 d_N` with and without jet frames, and the `d_N` layer traces annihilate it; `Tr_J(H_mu) = 7 + mu eps_v`, trace state blind, `mu` separated from `-mu`; Jordan type detects a nonzero layer | an off-diagonal layer coefficient is not invariant; `Tr_J` is not invariant under independent in/out frames; jet frames do not cut the invariant dual below `d_N`; `tr Reg_N` reads no layer; Jordan type does not read the layer value; no linear congruence invariant of the metric record |
| T4 | six substitutions (shears, scaling, two rational rotations, a curved one) are unital *-homomorphisms fixing the observer and of full rank; shear and curved-chart inverses compose to the identity through depth `N`; the SL(2) differences span codimension exactly 1 (the observer); the rotation differences span codimension exactly `floor(N/2)+1` and the Laplacian readings annihilate them and are independent; layer coefficients equal `d^alpha / alpha!` at `p` | a singular linear part is not an automorphism; a substitution moving the point breaks the quotient relation; the Laplacian reading is not shear-invariant (`N >= 2`); rotations do not collapse the dual to the observer (`N >= 2`) |
| T5 | TET-T6 family `diag(1 + eta t^(N+1), 1)`: transport jets equal through `N`, first difference `-eta/2`, exact rational endpoint `diag(2/3, 1)` at `t = 1/2`; TET-T4 affine path: retained tower equals the binomial closed form | the truncated tower evaluated at the endpoint is not the transport |

The rank checks are two-sided: an explicit invariant family annihilates
the difference span (lower bound), and the exact rank of the span gives
the upper bound. For rotations a single rational rotation of infinite
order suffices, because the span of `theta^k x - x` over all `k` is the
image of `theta - 1`. These are checks of the written proofs at three
depths; they certify nothing at other depths or in the continuum.

## 9. Dependency ledger

| Used | Source | Label |
|---|---|---|
| `J_N`, products (4), involution, `A_N`, dagger, `Pi`, (7), `Reg_N` (8)–(9), frame contract (10), curvature readout (14), source (15)–(20) | TLM-1 at `babb4fb` (`papers/native-tvsp-lambda-matrices/`), unchanged | [D] |
| Transport (1), recurrence (2), affine path (9), finite-depth obstruction TET-T6 (15)–(16), bound (17) | TET-1 at `9a06216` (`papers/native-tvsp-exponential-transport/`), unchanged | [D] |
| Response Hessian and third-moment identity (3), connection (4), intertwiner (8) | QRT-1 at `59c68ca`, read only as an interface in Section 7 | [D], not used in any proof |
| Wedderburn structure of `M_m(C)`-modules (Theorem 1(c) only) | standard | [E] |
| Cauchy–Schwarz for positive functionals | standard; proof included (Lemma 2.1) | [E] |
| Hadamard's lemma (Taylor with smooth remainder) | standard | [E] |
| `exp(-1/s^2)` is smooth and flat at `0` | standard | [E] |
| Identity theorem; Cauchy's theorem for linear analytic ODEs | standard | [E] |
| Denjoy–Carleman quasi-analytic classes | named only in the open item of Theorem 5(c) | [E], not used |
| Schwartz's point-support theorem | remark only | [E], not used |
| Mechanism comparison with `ERGODICITY.md`, `BRIDGE.md` | my earlier reports (a9f44d5, 97bb7fe) | remarks only, not used |

Nothing from the shared branch after `35d63e4`, from POS-1, WCG-1,
PCR-1 or DCD-1 is used in a proof.

## 10. Open items

1. [C] The exact function class in which the complete tower determines
   transport (Theorem 5(c)).
2. [C] Classification of the polynomial frame- and chart-invariants of
   the layered metric record (where `kappa` sits); classical invariant
   theory would answer it, but it is Codex's curvature-reconstruction
   track and is left as an interface.
3. The QRT state-family uniqueness statement sketched at the end of
   Section 7.
