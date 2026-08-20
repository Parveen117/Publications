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
| Remaining: physical (asymptotically-free) trajectory, growing lattice, universality, OS reconstruction, Clay predicate | OPEN (the real wall) |
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
python papers/yang-mills-certified-benchmark/certificates/ym9_uniform_heat_kernel.py
```

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
