# Lineage: the lambda-framework Yang-Mills manuscripts

Two pre-framework manuscripts ("The Yang-Mills Mass Gap as a
lambda-Topological Positivity Condition", v1/v2) proposed a conditional
route: IF a vacuum state omega_0 exists on the observable algebra, THEN
sector decorrelation (ACA) + local orthogonality (LO) + a Schur contraction
estimate + kernel positivity force a mass gap ("Kill Lemma").

Ruling (mine, don't merge — same doctrine as the morphic-calculus chapters):

- NOT imported: the vacuum conjecture, the ACA/LO "first-principles" proofs
  (prose, tied to the retired lambda-RH citations), and any continuum claim.
  Those documents remain uncertified lineage.
- MINED: the Kill-Lemma SKELETON — (block Schur estimate) + (positivity)
  => gap. YM-5 realizes exactly this skeleton on the theta-graph carrier,
  unconditionally and machine-certified:
    their Schur estimate  ->  U1/U2 block + off-block interval bounds
                              (off-block via the exact doubling m_k^2=m_2k)
    their kernel positivity -> certified character-coefficient positivity
                              (f_j = 2 I_{2j+1}/kappa > 0, positive series)
    their Kill Lemma       -> U3 Weyl bound + interlacing => two-sided gap.
- The old manuscripts' honest reduction ("vacuum exists => gap") maps onto
  the still-open gates of MP adapters/yang_mills/DEPENDENCY_MAP.md; nothing
  here advances those gates.

## YM-15 (Aug 22 2026): the 1D block-transfer dock, opened

- CONSUMES: YM-14 T4 (the named object), YM-14 carrier conventions
  (beta=2 Wilson T0, A-line basis, B-parity inertness), YM-3/4/14's exact
  SU(2) pairing integral generalised to the convolution lemma
  Int chi_a(XA^{-1})chi_b(AY) = delta_ab chi_a(XY)/d_a, YM-6 exact LDL
  inertia, YM-7 inertia-bisection pattern.
- SESSION LINEAGE: an earlier working session reached the chain topology,
  the convolution formula (after two self-corrections) and exact spot
  values at m = 4, 6, 8. This capsule replaces the spot checks by the
  closed form for every m, with the formula machine-checked as a formal
  polynomial identity (engine vs closed form, m = 2..8, all entries).
- IMPORTS ZERO numerical claims; the m=4,6,8 spot values are superseded.
- HONEST: A-line carrier only; the complement of W_{m+1} grows with m and
  is untouched. YM-16 (m-uniform complement bound) is the named next dock.

## YM-16 (Aug 22 2026): the chain dock — finite volume certified, death volume computed

- CONSUMES: YM-15 closed form (both M_kappa and M_{2kappa} on the chain
  carrier), YM-6/14 Haynsworth dock, YM-5 doubling identity (survives the
  bridge product verbatim), YM-1 Delta_red (reappears as the exact
  m-uniform UPPER bound on the chain gap via B-factorisation).
- DIRECTION: the Millennium predicate is volume-uniformity. This capsule
  certifies the chain gap for m <= m*(kappa) and PROVES the sup route
  cannot go further (per-bridge price e^kappa/f0 > 1 is exponential in m).
  The wall is now a formula, not a mood: the obstruction is the use of
  |m_kappa|_inf on the complement; its replacement (local/per-bridge
  control in the content basis = certified strong-coupling cluster
  radius, Osterwalder-Seiler shape) is YM-17.
- IMPORTS ZERO external claims. Jentzsch anchor not used.

## YM-17 (Aug 22 2026): the interleaving seam

- CONSUMES: YM-16 dock at m = 2 (the pair), YM-5 doubling (pairwise),
  YM-15 closed form, singular-value Weyl inequality (classical, cited).
- RESULT: half-chains are exactly volume-uniform (product spectrum); the
  A-chain's entire volume problem is the vacuum-tracking inequality
  lambda_1(T) vs sigma_1(X) sigma_1(Y). Vacuum bracket per bridge now
  O(kappa^2) at both ends. No m-uniform chain gap claimed.
- YM-18 named: certified Kotecky-Preiss radius for the interleaving
  polymer gas (two-site pairs as monomers). This is the first capsule
  whose next object is a classical convergence criterion rather than a
  new carrier — the program has reached the cluster-expansion wall with
  its own coordinates.

## YM-18 calibration (Aug 22 2026): vacuum tracking measured, cluster dock designed

- CONSUMES: YM-15 engine (generalised to multiple insertions per site),
  YM-16 f0, YM-17 upper end. PC-1 precedent for "calibration, not claim".
- MEASURED: the enlarged-carrier vacuum rate is m-independent to ~1e-5
  for m = 2..6 at every grid kappa. Evidence, not proof, that the
  interleaved vacuum tracks a per-bridge rate.
- DESIGNED, NOT EXECUTED: the Kotecky-Preiss dock (T3). Every line is an
  obligation. YM-18 proper = discharging them.

## YM-19 (Aug 22 2026): the Dobrushin dock — volume-uniform gap certified

- CONSUMES: YM-9 heat-kernel family (time step a), YM-8 Perron floors
  (simple positive vacuum), YM-16 B-factorisation (consistency control),
  YM-5 exp two-route, YM-13 compound exp route (e^{-x} = (e^{-x/n})^n).
- ANCHOR (CIRC-1, cited not rederived): Dobrushin uniqueness and the
  Foellmer/Kuensch covariance decay (Georgii ch. 8). Same governance as
  YM-8's Jentzsch anchor, which YM-12 later demoted; the program should
  look for a square-sourced replacement of this anchor too.
- STANDING CORRECTION: YM-17 T3's "vacuum-tracking inequality" is not a
  sufficient reduction; recorded in the certificate. YM-18 retained as
  calibration. The KP design (YM-18 T3) is superseded as the main route
  but remains valid as an alternative.
- CLOSES: volume-uniformity of the interacting gap on the bounded-degree
  chain at strong coupling, fixed time step. Does NOT close cutoff
  uniformity — the Dobrushin condition degrades as a -> 0. That is now
  the sharp wall, one level deeper than YM-9 left it.

## GOVERNANCE (Aug 22 2026, owner's rule): NO CLASSICAL BORROWING

- Owner's standing rule: the program proves with its own machinery
  (EMK/UGD algebra, generalized Euler flow, seam calculus); classical
  results are not imported as load-bearing anchors.
- Applied: YM-19 T3 DEMOTED to "anchored"; the m-uniform gap claim is
  withdrawn until derived natively. YM-19 T1/T2 (exact, self-derived)
  stand. YM-8's Jentzsch anchor was already demoted by YM-12 the same
  way. Audit of YM-1..19: these were the only two borrowed anchors.
- Standing obligation: a native replacement for correlation decay.

## YM-20 (Aug 22 2026): native origin audit + Cayley form

- READ for this capsule: Recognition-Kernel-Framework F00-E (native
  Euler, circular-hyperbolic via sheet involution), F00-G (Cayley
  coordinate, odd-series log, Exp(A(y))=(1+y)/(1-y)), RH-Framework T01
  (cut-tail mass, recognition energy, native action bound) and its
  NATIVE_DERIVATION_CONTAMINATION_AUDIT (vocabulary adopted verbatim).
- FINDING: the YM carrier is a coordinate shadow by the framework's own
  standard (same class as RH T03-A); the counting layer is native; the
  program's uniform bounds are Cayley exponentials of seam coordinates.
- Framework-internal pins (F00-G 5.1/7.1, F00-E 2.1/6.2, T01-C) are
  citations INSIDE the program, permitted as in EMK-2.
- YM-21 named: native contraction replacing the decay anchor.

## YM-F1 (Aug 22 2026): the chain as a recognition fabric — native carrier

- SOURCE (vault, read this session): APPENDIX_DISCRETE_HOLONOMY_FABRICS_AND_
  RECOGNITION_STOKES — the program's lattice objects (plaquettes, Wilson
  action, gauge-fixed tree) are its fabric, face holonomy, face residue
  and rails. Also read: generalized_euler_phase_ratio_space and
  EULER_INFORMATION_CURVATURE_DUALITY (no decay mechanism — recorded),
  mp_gold/01 (YM-1..8's ancestor: MP PR#30), gold/05 (curvature-to-seam-
  index spectral-flow law), gold/06 (same-projection two-block bound).
- CONSUMES: EMK-1 (Publications) determinant seam ladder; F00-E iota_Sigma.
- NEW NATIVE FACT: SU(2) is the unit sphere of the EMK primitive block
  with iota-twisted seam channels — the gauge group was inside the
  framework's own algebra. EMK-1's channel witness (2,2,3,1) has split
  determinant 8 untwisted and norm 18 twisted (test).
- Carrier layer of YM-1..18 is now NATIVE except the declared Haar <->
  Phi_Sigma identification. No number changed. No gap claimed.

## YM-21 (Aug 22 2026): the tiling law

- STANDING CORRECTION: YM-F1's "native form of YM-21" was ill-posed
  (boundary closure never bounds interior residue — the fabric
  appendix's own principle); killed by witness (g, g^{-1}) and replaced.
- DERIVED, not docked: YM-15's compressed A-block is the leading-order
  tiling transfer of the 2D fabric; YM-15/20's m-uniform bound is a
  tiling count. Spatial decay exact via Stokes sub-telescoping.
- OPEN: full-order time decay uniform in m = convergence of the
  face-coefficient tiling expansion (native cluster statement).
  Expansion parameter certified small and decreasing (r_j ladder).

## YM-22 (Aug 22 2026): the tiling rate survives the cutoff

- CONSUMES: YM-9 heat-kernel family + declared trajectory; YM-21 T3
  tiling transfer (scalar in lambda); YM-20 native odd-series log.
- RESULT: at leading tiling order the chain's rate per unit time is
  uniform in m AND a along kappa = theta a, infimum exactly 3/4 - theta/2.
  First native touch of the cutoff wall. Full order OPEN: in the limit
  only j = 1/2 branching tilings survive (r_j/a -> 0 for j >= 1).
- Build note: I first wrote the convergence as approach-from-below; the
  arithmetic showed approach-from-above (the limit is a uniform lower
  bound) — controls corrected to what the arithmetic says (YM-14 C4
  lesson again).
- CI LESSON (b2291d1 failed on hosted Tests, pins passed): YM-F1/YM-21
  seeded the shared `random` module at import; in the full pytest run
  the test order changed the stream and YM-F1's pinned demo drifted.
  Fix: reseed inside run(). Pins unchanged. Rule: a pinned certificate
  must be a pure function of its inputs, never of module import order.

## YM-23 (Aug 22 2026): the weak-coupling turn

- CONSUMES: YM-F1 stereographic chart + quaternion product, YM-21/22
  tiling rate, EMK-1 T5 cut-commutator curvature, RST-1 T4 flow loop,
  F00-G Cayley/odd-series log (native tanh and log 2 = A(1/3)).
- RESULT: exact boundary of the strong-coupling route; native
  linearization with the commutator as the only curvature; exact
  negative — abelianized fabric has no volume-uniform gap. The
  weak-coupling gap problem is now the commutator-residue transfer.
- DECLARED, not certified: AF trajectory shape; harmonic excitation
  relation. No weak-coupling gap claimed.

## YM-24 (Aug 22 2026): the abelian subfabric and the non-abelian residue

- CONSUMES: YM-F1 quaternion chart, YM-23 soft mode, F00-E Thm 3.1/5.2
  (native addition law / circular orbit), EMK-1 rotational determinant
  channel, EMK-T2 (order is content).
- RESULT: the commutator sector's exact structure — abelian subfabrics
  are copies of the native circle with additive Stokes; the non-abelian
  residue is a Gram defect; the gapless soft mode is exactly abelian.
  Weak-coupling gap, if any, lives in direction changes between
  subfabrics. YM-25 named. No gap claimed.

## YM-25 (Aug 22 2026): the direction-change sector

- CONSUMES: YM-10 multiplicity law, YM-9 semigroup, YM-F1 chart, YM-24
  rotation, RST-2 squared-amplitude discipline, MP gold/01 (PR #30
  spin-network seed — reopened here at weak coupling).
- RESULT: time kernel and class-function carriers are both
  direction-blind (exact); direction change = 6j recoupling with exact
  rational squares 1/4 : 3/4. The weak-coupling gap is isolated to the
  intertwiner transfer. YM-26 named. Frontier status declared.

## YM-26 (Aug 22 2026): dock to theorum/28

- OWNER'S CORRECTION recorded: the framework already has the finite-to-
  infinite machinery (recognition-Cauchy + Smriti tails + outward
  certificate u_n + e_n < 1). YM-16 was that theorem without its name.
- RESULT: YM-16's dock re-certified as a lawful theorum/28 outward
  certificate for m <= m*; the remaining obligation is theorum/28
  hypothesis 3 (memory-channel Cauchy bound uniform in m), reduced via
  the tiling law to the content-1/2 multi-insertion channel.
- RULE going forward: before declaring any object "infinite-dimensional
  with no finite truncation", check theorum/28's hypothesis list first.
