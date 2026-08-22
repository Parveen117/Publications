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
| Remaining: m-uniform complement bound for the chain (YM-16, named), 2D lattice, AF trajectory, tightness, OS reconstruction, non-triviality, metric universality, Clay predicate | OPEN |
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
