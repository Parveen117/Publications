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
