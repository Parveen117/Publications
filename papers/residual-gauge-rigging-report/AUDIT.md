# Audit note: ERGODICITY.md revision of 16 October 2026

Response to the five-point audit of the first version (commit 50ccc41).
Changes by point, then remaining gaps, then the verdict. Labels as in
`ERGODICITY.md`.

## 1. Invariant-measure inference (stated error; removed and replaced)

The old Section 6 listed "all G_c-invariant probability measures
absolutely continuous with respect to nu_2: nu_2 is the only one". That
presumed nu_2 is invariant. It is not: RGR gives quasi-invariance with
the cocycle R^(2) and a weighted invariant **functional** (20)–(21), not
an invariant measure. The bullet is removed. What is actually proved now
stands separately as Proposition M [P]: (a) nu_2 is quasi-invariant and
**not** invariant (via `d/dt Delta|_0 = 2 A_f` and `∫ A_f^2 > 0`, RGR
(23)/(27) [D]); (b) **no** nonzero finite invariant measure absolutely
continuous with respect to nu_2 exists at all (Theorem E makes a
candidate equivalent to nu_2; the escape-to-infinity shells of Theorem
U's Step 1 then force total mass zero); (c) uniqueness, had one existed,
by the cocycle-ratio argument. The sigma-finite infinite case is not
addressed. Theorem E is restated as ergodicity of the **nonsingular**
action (invariant sets are null or conull; equivalently invariant L^∞
functions are constant), which is exactly what Theorem U(b),(d) consume,
so the classification corollary is unaffected. Lemma 1 is annotated: it
uses quasi-invariance only.

## 2. Zero-weight sectors (gap; fixed)

The unitary-equivalence claim through `W = sqrt(c_+)P_+ + sqrt(c_-)P_-`
now carries the hypothesis c_+, c_- > 0, where W is invertible. Theorem
S' gains the explicit trichotomy: both weights positive, unitary
equivalence to RGR's (29) and equality iff c_+ = c_- = 1; exactly one
weight zero, the B-null space swallows the unsupported sector and the
physical quotient is the one-sector theory (not unitarily equivalent to
(29)); both zero, the zero form with trivial quotient, included in the
classification. Witness W5.

## 3. Lemma 2 coefficient expression (derivative factor; fixed)

The displayed finite expansion of u_K read `sum rho_K(k) G_k P_0[e_k o I]`;
correct is `sum rho_K(k) Phi_k P_0[e_k o I]` with `Phi_k = 2 pi i k G_k`
(equivalently `G_k` paired with `e_k'`), the factor coming from
`v_{e_k} = P_0[e_k'(chi)]` with the angular phase of e_k'(chi) cancelling
against the phase of c_k. The composed representation
`u_K = P_0[(V_K Phi) o I]`, from which Lemma 4, Lemma 3's application and
Proposition 5 were derived, was and is correct, so no estimate changes.
Fourier conventions (e_k = e^{2 pi i k y}, derivative factors, real
cos/sin pairing, the witness file's "up to 2 pi i") are now fixed in
Section 1, and both expansions are displayed and reconciled in Section 3.
Witness W1b checks the corrected assembly exactly and rejects the
pre-audit formula.

## 4. Spectator compatibility contract (formulation tightened)

- (BC2) is split. (BC2a): weak symmetry for S(f) and D(h) on D_red,
  which these operators preserve. (BC2b): the J(f_0) condition is stated
  through (BC1)'s continuous extension, `<J v, T u> = <v, T J u>` on
  A_0 ⊂ Dom J, because J(f_0) maps D_red outside itself and "B(Jv, u)"
  is undefined; the bounded operator T is the meaning of the contract.
- (BC3) may replace **(BC2b) only**; (BC2a) is retained in both
  variants, as required.
- Gauge invariance is required **slotwise** (first slot plus Hermiticity
  suffices). Remark N [P] shows this cannot be weakened: the kinematic
  L^2(nu_1) pairing is positive, separately continuous and invariant
  under the **simultaneous** action (U_g is nu_1-unitary), yet not
  slotwise invariant and not a multiple of conj(ell)⊗ell (by RGR (22)).
  So simultaneous sesquilinear invariance alone does not imply
  factorization through the invariant dual, and that larger class is
  explicitly left unclassified. Witness W6.
- Newly explicit step: the within-sector constancy used the Gaussian 0–1
  law for **all** Cameron–Martin shifts while the contract only yields
  smooth-direction shifts. Lemma S3 [P] supplies the dense-shift upgrade
  (L^1 continuity of translations on bounded functions). The unbounded
  alternative to (BC1) is named (direct invariance under exp(-itO)) and
  not analyzed.

## 5. Continuum-proof audit (survives; implicit steps made explicit)

Chain audited: Lemma 1 → (2.1) → Lemma 2 → Lemmas 3–4 → Proposition 5 →
Theorem E. Findings:

- **Coefficient membership in A**: now Lemma 0 [P], with the
  (1+|k|)^{r+1} poly-exp bound, uniform in angles. The k-growth is
  harmless because no K-uniform estimate runs through coefficient sums;
  all run through the composed form of u_K. (This is now stated.)
- **H^4 composition/inversion**: now Lemma A [P], quantitative bounds
  for I and I^{-1} at integer s <= 4 with poly-exp constants, proved by
  Faà di Bruno and change of variables. Inci–Kappeler–Topalov is demoted
  to a context citation; no external hypothesis left unverified there.
- **Gaussian D^{1,2} membership**: now Lemma B [P using E], via scalar
  coordinates, finite-dimensional Gaussian Sobolev membership (Bogachev
  §5.2), Nualart's closability Lemma 1.2.3, and the coordinate-sum
  characterization of D^{1,2}(H_2). Also covers A ⊂ D^{1,2}, used for
  the test functions.
- **Hilbert–Schmidt summability**: the basis decay
  `||D~_{eps_n}||_{A^1} <= C(sigma)<n>^{-3}`,
  `||D~_{eps_n}||_{H^4} <= C(sigma)<n>^{-1}` is derived from Lemma A;
  the n-sum is dominated by sum <n>^{-2}.
- **K-independent domination**: a single G(sigma) with E[G^2] < ∞
  (Fernique) now dominates `||u_K||_{H_2} + ||Du_K||_{HS}` uniformly in
  K, via the K-uniform commutator bound (3.1). All limits are dominated
  by it.
- **Complete divergence expression**: the L^1(nu_2) limit is assembled
  from three pieces (angular, with (V_K G)(0) → G(0) = 0; gradient, by
  Cauchy–Schwarz against Proposition 5; divergence, by the Skorokhod
  isometry bound on u_K − h), then identified with delta_2(psi h) by the
  product rule, whose hypotheses (psi ∈ D^{1,2} with moments, h
  deterministic) are checked. Pointwise vector-field convergence is
  nowhere used as a substitute.
- **Endgame**: the chaos/rotation/two-factor step is written out with a
  countable test family, the Hermite recurrence, Fubini and rational
  rotations; the per-chirality bookkeeping on X_c is stated.
- **External references with hypotheses**: Nualart Thm 1.1.1, Lemma
  1.2.3, Prop 1.3.1, Prop 1.3.5; Bogachev §5.2, CM 0–1 law, Fernique;
  Nelson; maximal abelian L^∞; Hermite recurrence. Each is cited at its
  point of use with its hypotheses discharged there; the ledger in
  Section 7 separates them from the RGR dependencies [D], which are
  (10)–(11), (18)/(21), (19), (22), (23), (26), (27), T4's derivative
  estimates, D2, CGD-T2 step 3, MPQ-T3.

## Remaining gaps and assumptions

- Sharp Condition (L): undecided, unchanged (Section 4); the defect
  remains the absent 1/K gain at the H_2 borderline. Not independently
  resolved, so left open as instructed.
- Sigma-finite infinite invariant measures ≪ nu_2: open.
- Invariant measures singular to nu_2: not addressed.
- Forms with only simultaneous gauge invariance: strictly larger class,
  unclassified (Remark N).
- J(f_0) core on A_0: not proved; not needed under (BC2b)/(BC3).
- RGR-1 is assumed as stated (ledger items above); it is not re-proved.
- Finite witnesses certify algebraic identities and instances only;
  no continuum theorem, and not ergodicity.

## Verdict

The main ergodicity argument **survives the audit**. The two substantive
repairs — the Section 6 invariant-measure inference (removed, replaced by
Proposition M) and the displayed coefficient expansion in Lemma 2 (the
derivative factor) — do not touch the proof chain, which ran through the
composed form of u_K and quasi-invariance throughout. The remaining
changes add previously implicit lemmas (0, A, B, S3), the K-independent
domination, the slotwise-invariance hypothesis with its counterexample,
and the zero-weight trichotomy.
