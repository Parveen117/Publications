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
| YM-6 larger compression V / theta-network basis; push certified region toward strong coupling | OPEN |
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
python papers/yang-mills-certified-benchmark/certificates/ym5_two_sided_gap.py
```

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
