# The bridge from infinitesimal Dirac constraints to finite gauge invariance on RGR's test space

16 October 2026. Separate branch `claude/residual-gauge-invariant-dual`,
continuing from a9f44d5. RGR-1 at d9c7020 is a stated dependency, read and
audited, never modified. Nonlinear material geometry, reference-momentum
observables and detector questions are the other track's. Sharp Condition
(L), sigma-finite invariant measures and unbounded-form generalizations
are not touched; nothing here needs them.

Labels: **[P]** proved here; **[D]** RGR-1 dependency (ledger in Section
8); **[E]** external theorem; **[C]** open. Finite witnesses (Section 9)
certify no continuum statement.

## 0. The object and the result

Throughout: Psi, its Frechet topology (19) (seminorms p_R, ||.||_{L2(nu_1)},
||.||_{L1(nu_2)}), the algebra A, the generators `Q_eps(f) = i B_f` for
real smooth f acting on chirality eps, and `D = span{U_g a : a in A, g in G_c}`
are **exactly RGR's** (RGR (19), (23), (24)) [D]. D is invariant under
every U_g (by its definition and the group law) and under every Q_eps(f)
[D: RGR-T4]. Define the **infinitesimal solution space**

    K_inf = { F in Psi' : F(Q_eps(f) psi) = 0 for all real smooth f,
              both chiralities eps, and all psi in D },

where Psi' is the space of linear functionals continuous for (19). The
condition is well posed because Q_eps(f) psi lies in D ⊂ Psi [D]. Let
`Dbar` denote the closure of D in (19).

**Theorem K (structure of the infinitesimal solution space) [P over D, E].**

    K_inf = C·ell ⊕ { F in Psi' : F = 0 on Dbar }.

In particular:

(a) Every F in K_inf agrees with `F(1)·ell` on Dbar. On the constraint
    core and its closure there is exactly one infinitesimal solution
    channel: **finite-group uniqueness (Theorems U + E) does imply
    completeness of the infinitesimal physical solution space there,
    unconditionally.**
(b) `K_inf = C·ell` **if and only if** D is dense in Psi for the full
    topology (19), and `dim(K_inf / C·ell) = dim (Psi/Dbar)'`.
(c) Any hypothetical extra element of K_inf vanishes identically on
    Dbar, which contains A, every U_g-translate, every bounded uniformly
    continuous cylinder function in Psi, and all of RGR's constructed
    states (Section 5); it is invisible to every pairing RGR writes down.

The density of D in (19) is **not proved in RGR-1** (its proof gives
Hilbert-space density only) and is **not decided here**: Section 5 shows
cylinder approximation provably cannot give it, that no obstruction
exists at first or second order (Lasry–Lions), and that the question is
equivalent to a sharp uniform-approximation problem on Hilbert balls,
left [C]. Section 6 states the revised test-space contract (explicitly
distinct from the frozen one) under which `K_inf = C·ell` holds
unconditionally. Section 7 transfers the classification to positive
slotwise constraint-annihilating rigging forms under the corrected
spectator contract, with zero weights and normalization explicit.

Plain answer to the assignment's final question: **yes on the constraint
core Dbar; on all of the frozen Psi it is exactly equivalent to the
unproved (19)-density of D**, by Theorem K(b).

## 1. Hypotheses, stated explicitly

- Topology: F ranges over (19)-continuous functionals; nothing weaker
  (Hilbert continuity) and nothing stronger is assumed. Every F in Psi'
  is bounded by finitely many seminorms:
  `|F(psi)| <= C (p_{R_0}(psi) + ||psi||_{L2(nu_1)} + ||psi||_{L1(nu_2)})`
  for some R_0, C (Frechet structure; p_R increase in R).
- Domain: the constraints are tested on D and only on D, matching RGR-T4,
  where the generators are defined and (25) holds [D].
- Both chiralities are smeared; rotations are the constant smearings.

`ell ∈ K_inf` is RGR (25) [D]. So `C·ell ⊆ K_inf`, and the question is
whether anything else solves all constraints.

## 2. The bridge, part 1: infinitesimal to one-parameter flows

**Lemma B1 (orbit differentiability in Psi on D) [P from D].** For every
real smooth f, every g in G_c and every a in A, the curve
`t → U_{g_f(t)} U_g a` is differentiable in the topology (19) at every t,
with derivative `U_{g_f(t)} U_g B_{f~} a`, where `f~ = (g^{-1})_* f` is
the pullback vector field, and

    U_{g_f(t)} U_g B_{f~} a = Q-image of an element of D:
    U_{g_f(t)} U_g B_{f~} a = -i Q(f) [ U_{g_f(t)} U_g a ].

*Proof.* RGR-T4 proves that `t → U_{g_h(t)} a` is Psi-differentiable for
a in A and real smooth h, with derivative `U_{g_h(t)} B_h a`, and that
B_h preserves A [D: (26) and its uniform differentiation on bounded H^3
sets]. Conjugation of flows: `g^{-1} g_f(t) g = g_{f~}(t)` with
`f~ = (g^{-1})_* f` smooth [D: "conjugating a flow by g is another smooth
flow"], so `U_{g_f(t)} U_g = U_g U_{g_{f~}(t)}` by the representation
property. Since U_g is (19)-continuous [D: RGR-T3], the chain rule for a
continuous linear map gives differentiability with derivative
`U_g U_{g_{f~}(t)} B_{f~} a`. Finally, on the Stone domain the group
commutes with its generator and intertwines under U_g
(`U_g B_{f~} = B_{g_* f~} U_g = B_f U_g` [D: T4]), which rewrites the
derivative as `B_f` applied to `U_{g_f(t)} U_g a ∈ D`. ∎

**Lemma B2 (integration) [P].** Let F ∈ K_inf. Then for every real
smooth f, every t, and every chi in D:

    F(U_{g_f(t)} chi) = F(chi).

*Proof.* By linearity take chi = U_g a. Set `phi(t) = F(U_{g_f(t)} chi)`.
By Lemma B1 and continuity of F, phi is differentiable at every real t
with `phi'(t) = -i F( Q(f) [U_{g_f(t)} chi] ) = 0`, because
`U_{g_f(t)} chi ∈ D` and F ∈ K_inf. A real- (or complex-) valued
function of one real variable, differentiable everywhere with vanishing
derivative, is constant (mean value theorem on real and imaginary
parts). ∎

No continuity of t → phi'(t) is needed, and no exponentiation of an
unbounded operator is invoked: only RGR's own Psi-differentiability and
one-variable calculus.

## 3. The bridge, part 2: flows to the whole residual group

Not every orientation-preserving smooth circle diffeomorphism is a
single autonomous flow map; that is not assumed. What is used:

**Lemma B3 (flows generate) [P using E].** Let G_0 be the subgroup of
`Diff^∞_+(S^1)` generated by all time-t maps g_f(t) of smooth vector
fields. Then `G_0 = Diff^∞_+(S^1)`.

*Proof.* The generating set is closed under inversion
(`g_f(t)^{-1} = g_f(-t)`) and under conjugation: for any h,
`h g_f(t) h^{-1} = g_{h_* f}(t)`, the flow of the pushforward field,
which is again smooth. A subgroup generated by a conjugation-invariant
set is normal. G_0 contains all rigid rotations (flows of constant
fields), so G_0 ≠ {id}. `Diff^∞_+(S^1)` is a simple group
[E: Herman, building on Mather; see Banyaga, *The Structure of Classical
Diffeomorphism Groups*]. A nontrivial normal subgroup of a simple group
is the whole group. ∎

This is where the actual group enters; only the abstract simplicity of
the C^∞ group is used, and representation continuity does the rest:

**Theorem B (the bridge) [P over D, E].** For F in Psi', the following
are equivalent:

1. F ∈ K_inf, i.e. F annihilates `Q_eps(f) D` for all smooth f and both
   chiralities;
2. `F(U_g psi) = F(psi)` for every g in G_c and every psi in Dbar.

*Proof.* (2 ⟹ 1): for a ∈ D, `(U_{g_f(t)} a - a)/t → B_f a` in (19)
[D: Lemma B1 at t = 0], so
`F(B_f a) = lim t^{-1}[F(U_{g_f(t)} a) - F(a)] = 0`, and Q = iB.

(1 ⟹ 2): By Lemma B2, F is invariant under every flow map on D. Let
`H = {g ∈ Diff^∞_+(S^1) : F(U_{(g, id)} psi) = F(psi) for all psi ∈ Dbar}`
(chirality +; the − case is symmetric). H is a subgroup: closure under
composition uses that Dbar is U-invariant (D is, and each U_g is
(19)-continuous [D]), so `F(U_{g1 g2} psi) = F(U_{g1}(U_{g2} psi)) = F(U_{g2} psi) = F(psi)`.
H contains every flow map: invariance holds on D by Lemma B2, and
extends to Dbar because F and U_g are continuous. By Lemma B3, H is all
of `Diff^∞_+(S^1)`. For a general `g = (g_+, g_-) ∈ G_c`, write
`U_g = U_{(g_+, id)} U_{(id, g_-)}` and apply both chiralities. ∎

## 4. Classification on the closure of the constraint core

**Theorem U' [P over D, E].** Let F ∈ Psi' satisfy `F(U_g psi) = F(psi)`
for all g ∈ G_c and psi ∈ Dbar. Then `F = F(1)·ell` on Dbar.

This is Theorem U of `UNIQUENESS.md` with the invariance hypothesis
weakened from Psi to Dbar; the proof is the same argument run inside
Dbar, and every auxiliary object it uses lies in D. The three
ingredients are checked:

- **Stability under bumps [P].** Let `theta_M(x) = theta(||sigma_+||^2_{H^3}/M^2) theta(||sigma_-||^2_{H^3}/M^2)`
  with theta smooth, 1 on [0,1], 0 beyond 4. Then theta_M ∈ A, and for
  every g, `theta_M ∘ T_g ∈ A`: the map `x → u_g(x) = P_0 log g'(chi)`
  is smooth with poly-exp Fréchet bounds (the Lemma A machinery of
  `ERGODICITY.md` applied to chi = alpha + I_sigma and the fixed smooth
  log g'), and `||sigma + u_g||^2` composed with theta keeps A bounds.
  Hence, by the exact identity

      theta_M · U_g a = U_g ( (theta_M ∘ T_g) · a ),

  D is stable under multiplication by theta_M, and so is Dbar
  (multiplication by a fixed bounded continuous function is
  (19)-continuous). Therefore `Dbar_b := {psi ∈ Dbar, bounded support}`
  is (19)-dense in Dbar (theta_M psi → psi as in Theorem U, Step 1).
- **Escape (Step 1 of U) inside Dbar [P from D].** For psi ∈ Dbar_b and
  the finitely many seminorms bounding F, the escape elements
  g_n = (gamma_n^{-1}, id) of Theorem U give `p_{R_0}(U_{g_n} psi) = 0`
  for large n, while `U_{g_n} psi ∈ Dbar` and the L2 and L1 norms are
  preserved [D: unitarity and (21)]. Invariance on Dbar then gives
  `|F(psi)| <= C·N(psi)` with `N = ||.||_{L2(nu_1)} + ||.||_{L1(nu_2)}`,
  first on Dbar_b, then on Dbar by density and `theta <= 1`.
- **Splitting (Step 2 of U) inside Dbar [P from D, E].** The separating
  sequences of Step 2 are `a·theta_N` and `a(1 - theta_N)` with a a
  bounded smooth cylinder and theta_N the RGR (22) bump — all elements
  of A ⊂ D. So the image of D under `psi → ([psi]_{nu_1}, [psi]_{nu_2})`
  is already N-dense in `L^2(nu_1) ⊕_1 L^1(nu_2)`, F's N-continuous
  restriction to Dbar extends uniquely to the sum, the extension is
  invariant under the diagonal action by density and the two isometries
  [D], its Hilbert part vanishes [D: CGD-T2 step 3], and its L^∞(nu_2)
  density is G_c-invariant, hence **constant by Theorem E**. Evaluating
  at 1 ∈ A fixes the constant: `F = F(1)·ell` on Dbar. ∎

**Proof of Theorem K.** If F ∈ K_inf, Theorem B gives invariance on
Dbar, Theorem U' gives `F - F(1) ell = 0` on Dbar, so
`F ∈ C·ell + Ann(Dbar)` where `Ann(Dbar) = {G ∈ Psi' : G|_{Dbar} = 0}`.
Conversely `ell ∈ K_inf` [D: (25)], and any G ∈ Ann(Dbar) is in K_inf
because `Q_eps(f) psi ∈ D ⊂ Dbar` for psi ∈ D. The sum is direct since
`ell(1) = 1` with 1 ∈ D. Part (b): in a locally convex space, a subspace
is dense iff its annihilator in the dual is trivial (Hahn–Banach), and
`K_inf / C·ell ≅ Ann(Dbar) = (Psi/Dbar)'`. ∎

## 5. Density audit (assignment step 1)

**(a) What the frozen file proves.** RGR-T4 asserts "D is a dense
gauge-invariant subspace of Psi and H_c". Its proof paragraph establishes
the Stone domains, essential self-adjointness and (25); the only density
proved anywhere in RGR-1 is that bounded smooth cylinder functions are
dense in **H_c** (RGR-T3's proof). No argument for density in the
topology (19) appears. The frozen files are preserved; this is an audit
finding, not an edit.

**(b) Hilbert density is provably insufficient, and cylinders provably
cannot close the gap [P].** Let beta be smooth, 0 <= beta <= 1,
beta(0) = 0, beta(1) = 1. Then `psi* = beta(||sigma_+||^2_{H^3})` lies
in A ⊂ Psi. For **every** cylinder function c (a function of finitely
many modes and the angles): pick tau in the orthogonal complement of
those modes with `||tau||_{H^3} = 1`; the points x_0 = (alpha, 0, x_-)
and x_1 = (alpha, tau, x_-) satisfy c(x_0) = c(x_1) while psi*(x_0) = 0,
psi*(x_1) = 1. Hence `p_1(psi* - c) >= 1/2`. So the cylinder algebra,
which **is** dense in H_c [D], is at (19)-distance >= 1/2 from an
element of A. (The same psi* is approximated in L^2(nu_1) by cylinders
to any accuracy, which exhibits the gap between the two densities.)
Witness W9 replays the two-point mechanism exactly.

**(c) What Dbar provably contains [P].** (i) A, and every U_g A. (ii)
Every bounded uniformly continuous cylinder element of Psi: mollify the
finite-dimensional profile; the mollified cylinders lie in A and
converge uniformly on all of X_c by the profile's modulus, hence in
every p_R, and in L^2, L^1 by dominated convergence. (iii) Products of
elements of Dbar with the bumps theta_M and theta_M ∘ T_g (Section 4).
(iv) All of RGR's constructed objects: the (22) sequences, the
separating bumps, 1, and the common core D itself. So an element of
`Ann(Dbar)` pairs to zero with every state RGR-1 constructs or uses.

**(d) No first- or second-order obstruction [P using E].** Every
psi ∈ Psi with bounded support is the uniform-on-X_c limit of C^{1,1}
functions: inf-sup (Lasry–Lions) regularization on the Hilbert space
H^3 x H^3 (angles handled by finite-dimensional mollification) converges
uniformly with C^{1,1} regularizers [E: Lasry–Lions 1986, bounded
uniformly continuous functions on a Hilbert space]. Consequently no
counterexample to density can be detected by first- or second-order
(Taylor/Hessian-type) functionals; an obstruction, if any, must occur at
third order or beyond.

**(e) The exact remaining question [C].** By (b)–(d) and the bump
reduction, density of D in (19) holds **iff** every bounded uniformly
continuous function on a ball of H^3_0 x H^3_0 can be uniformly
approximated by smooth functions all of whose Fréchet derivatives are
bounded on that ball (A-regularity; each order may have its own
constant). Sufficient for `K_inf = C·ell`; its failure for some psi*
would, conversely, produce by Hahn–Banach a nonzero functional in
`Ann(Abar)` — which refutes `K_inf = C·ell` provided the functional also
annihilates the possibly larger Dbar (the U_g-translates add cocycle
factors whose A-regularity is plausible from RGR's derivative bounds but
is not needed or claimed here). Classical methods stop short of
settling it: smooth partitions of unity give C^∞ approximants with no
derivative bounds; Lasry–Lions gives bounded second differences but only
C^{1,1}; finite-dimensional mollification does not act on genuinely
infinite-dimensional arguments. I have not located a theorem or a
counterexample covering the bounded-all-orders case, and I do not decide
it here. This is the precise obstruction between "uniqueness on the
constraint core" and "uniqueness on all of the frozen Psi".

**(f) Consequences either way [P].** If density holds, `K_inf = C·ell`
on the frozen contract. If it fails, `K_inf = C·ell ⊕ (Psi/Dbar)'` with
the second summand nonzero — but every extra solution vanishes on Dbar,
so it assigns the value 0 to every state in (c), carries no physical
pairing against RGR's rigging, and cannot be normalized against ell
(no element of Ann(Dbar) takes the value 1 at 1). Completeness **on the
core** is unconditional either way (Theorem K(a)).

## 6. A revised test-space contract (explicitly not the frozen one)

Define `Psi* := Dbar` with the topology induced by (19). This is a
**revision**, distinguished from RGR's frozen contract (19)-on-Psi:

- Psi* is a complete Frechet space (closed subspace), gauge invariant,
  contains A and D, and D is dense in Psi* **by construction**.
- ell restricts continuously, ell(1) = 1; RGR-T3's continuity and
  invariance statements, T4's core statements (D ⊂ Psi*), and T5's
  pairing statements all restrict verbatim, because every object in
  their statements and proofs lies in D or Dbar (Section 5(c)).
- On Psi*: `K_inf(Psi*) = C·ell` unconditionally (Theorem K applied
  with Psi replaced by Psi*, where the annihilator term is trivial).

Nothing frozen is altered: RGR-1's files and claims stand as stated, with
the (19)-density of D in the **frozen** Psi recorded as unproved there
and undecided here.

## 7. Transfer to constraint-annihilating rigging forms (assignment step 4)

**Theorem R [P over D, E].** Let eta be a positive, Hermitian form on
`Psi ⊗_alg D_red` (or on `Psi* ⊗_alg D_red`), separately (19)-continuous
in the gauge slot, satisfying the corrected spectator contract of
`ERGODICITY.md` Section 5 — (BC1) bounded reduced form, (BC2a) S/D weak
symmetry on D_red, and (BC2b) J-compatibility through the bounded
extension, or (BC3) in place of (BC2b) — and **constraint-annihilating
slotwise**:

    eta( Q_eps(f) psi ⊗ v )[ phi ⊗ u ] = 0
    for all real smooth f, both chiralities, psi, phi ∈ D, v, u ∈ D_red

(by Hermiticity the second slot is annihilated as well). Then on
`Dbar ⊗ D_red`:

    eta(psi ⊗ v)[phi ⊗ u] = conj(ell(psi)) ell(phi) · <v, (c_+ P_+ + c_- P_-) u>,   c_± >= 0,

with exactly the trichotomy of the audited Theorem S': both weights
positive — unitarily equivalent to RGR's (29) via
`W = sqrt(c_+)P_+ + sqrt(c_-)P_-`, equal to (29) iff c_+ = c_- = 1
(isometric normalization); one weight zero — the quotient retains only
the supported sector; both zero — the zero form. Positivity of c_± comes
from eta's positivity at psi = phi = 1 ∈ D.

*Proof.* Fix v, u and phi ∈ D. The first slot is conjugate linear, so
`psi → conj( eta(psi ⊗ v)[phi ⊗ u] )` is a linear, (19)-continuous
functional on Psi, and it annihilates `Q_eps(f) D` by the slotwise
annihilation hypothesis. Theorems B and U' give, for psi ∈ Dbar,
`eta(psi ⊗ v)[phi ⊗ u] = conj(ell(psi)) · eta(1 ⊗ v)[phi ⊗ u]`.
The second slot is linear, continuous, and constraint-annihilating by
Hermiticity, so the same argument gives, for psi, phi ∈ Dbar,
`eta = conj(ell(psi)) ell(phi) B(v, u)` with
`B(v, u) = eta(1 ⊗ v)[1 ⊗ u]`. B inherits positivity and (BC); the
sector analysis of Theorem S' (Lemmas S1–S3 and the across-sector step)
applies verbatim and yields `B = <·, (c_+P_+ + c_-P_-)·>`. ∎

**Scope.** On arguments outside `Dbar ⊗ D_red` the constraints test
nothing, so no classification is asserted there (on the frozen Psi this
is again exactly the density question). The physical (GNS) quotient of
the constraint core under any admissible eta is therefore: {0} (zero
form), one sector's H_±, or H_red up to the unitary W — never a second
gauge channel. Combined with Theorem K(a): **the infinitesimal physical
solution space over the core is complete, and its positive pairings are
classified, with the normalization fixed exactly by RGR's isometric
convention.**

## 8. Dependency ledger

**RGR-1 [D]** (file papers/native-residual-gauge-rigging/THEOREM.md at
d9c7020): (19) Psi and topology; (21) L^1 isometry and unitarity; (22)
separating bumps; (23) B_f, Q(f), A, `B_f A ⊂ A`; (24) D and its
invariance under U_g and all generators; (25) `ell(Q_eps(f) psi) = 0` on
D; (26) with its uniform differentiation on bounded sets (the
Psi-differentiability of orbits on A); flow conjugation and the
intertwining `U_g B_{f~} = B_f U_g` (T4); T_g^{±1} map balls to bounded
sets (D2); CGD-T2 step 3. From my packets: Theorem U's Steps 1–2
(UNIQUENESS.md), Theorem E and Lemma A machinery and the (BC) contract
with Lemmas S1–S3 (ERGODICITY.md).

**External [E]:** simplicity of `Diff^∞_+(S^1)` (Herman, with Mather's
commutator theorems; exposition in Banyaga, *The Structure of Classical
Diffeomorphism Groups*); Lasry–Lions regularization on Hilbert spaces
(C^{1,1} uniform approximation of bounded uniformly continuous
functions); Hahn–Banach; and the [E] items already in ERGODICITY.md's
ledger (used through Theorems E, U, S').

**Open [C]:** uniform approximation, on balls of H^3_0 x H^3_0, of
bounded uniformly continuous functions by smooth functions with
ball-bounded Fréchet derivatives of every order — equivalently (Section
5(e)) the (19)-density of D in the frozen Psi, equivalently
`K_inf = C·ell` on the frozen contract.

## 9. Finite exact witnesses (`bridge/`)

Native JS exact arithmetic; run `python3.12 bridge/run.py --check`.
Instances and algebraic identities only; no continuum claim, and in
particular neither density nor its failure is certified.

- **W7 (flows generate via simplicity).** In the simple group A_5
  (5x5 exact permutation matrices), the conjugation-closed set of
  3-cycles generates the whole group (closure size 60). Negative
  control: in the non-simple A_4, the conjugation-closed class of double
  transpositions generates the proper normal subgroup V_4 (size 4, not
  12) — simplicity, not normality alone, is the needed input.
- **W8 (Theorem K's structure at finite level).** In a 4-dimensional
  model with core W, constraints Q_i mapping W into ker(ell|_W):
  the solution space of `F(Q_i w) = 0` has dimension exactly
  1 + codim(W), and every solution restricted to W is proportional to
  ell|_W. Negative control: a functional with F(w_2) ≠ 0 is not a
  solution.
- **W9 (cylinder insufficiency mechanism).** The two-point argument of
  Section 5(b): any function of the first coordinate is at sup-distance
  >= 1/2 from the two prescribed values at (0,0), (0,r) (exact:
  max(|t|, |1-t|) >= 1/2 since |t| + |1-t| >= 1), while a function of
  both coordinates achieves error 0.
- **W10 (infinitesimal to flow, exactly).** For the nilpotent 4x4 model
  e^{sN} (a polynomial flow, exact at rational s): F annihilating
  {N v, N^2 v, N^3 v} satisfies F(e^{sN} v) = F(v) identically in s;
  negative control: dropping the N^3 condition leaves a nonconstant
  cubic in s.
- **W11 (slotwise constraint-annihilating forms factor).** On the model
  (W ⊕ complement) ⊗ C^2: sesquilinear forms annihilating the
  constraint images in both slots and commuting with the declared sector
  observables form exactly the 2-parameter family
  `ell^† ell ⊗ diag(c_+, c_-)`; adding the sign-mixing observable cuts
  the dimension to 1. The ell-block forms have rank 2, 1, 0 for weights
  both positive, one zero, both zero (the trichotomy at finite level).

`CERTIFICATE.json` records `continuum_proofs_machine_certified: false`
and that the density question is not decided by any finite check.
