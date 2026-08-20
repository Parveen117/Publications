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
| YM-3 exact character transfer at general coupling + gap-closing direction | OPEN (next) |
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

CI regenerates the certificate and fails on any pin drift.
