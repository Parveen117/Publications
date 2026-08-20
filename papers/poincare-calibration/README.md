# Poincaré Calibration — PC-1

**This is a calibration, not a proof.** The Poincaré conjecture is a proved
theorem (Hamilton–Perelman, 2003). Nothing here claims any part of it and
nothing here is new mathematics. The capsule instantiates the dictionary of
`MP adapters/poincare/adapter.tex` v0.2 ("Poincaré Calibration Adapter"),
whose stated purpose is to test whether the framework can carry native
hypotheses, singularity control, surgery, extinction and topological
reconstruction **without hiding any step inside a decorative closure
predicate**. PC-1 answers that for one piece — exactly.

## Carrier

Left-invariant metrics on **SU(2)** (deliberately the same group as the
Yang–Mills line), diagonal in a Milnor frame `g = diag(A,B,C)`. Ricci flow
is an exact rational ODE: `A' = 2[(B−C)² − A²]/(BC)` and cyclically.

## Certified

| | Statement |
|---|---|
| **T1** | `λ(g) = min spec(−4Δ + R) = R(g)` on any homogeneous metric; `R = [2(AB+BC+CA) − (A²+B²+C²)]/(ABC)` verified |
| **T2** | **`dR/dt − 2‖Ric‖² ≡ 0`** in `ℚ(A,B,C)` — Perelman's λ-monotonicity as a *polynomial identity*, with an explicit 3-term sum-of-squares witness. No analysis, no tolerance, no approximation |
| **T3** | Round solution exact: `a(t) = a₀ − 2t`, `λ = 3/a`, extinction `t* = a₀/2` exactly |
| **T4** | Berger rounding exact: `du/dt = −4(x−z)/(xz)` for `u = x/z`, so `sign(du/dt) = −sign(u−1)` — monotone convergence to round |
| **T5** | Seam readout: λ is a bottom-of-spectrum threshold object, strictly increasing ⇒ threshold crossings are irreversible (one-directional seam flow; same primitive as the YM seam-integer dock, opposite regime) |

Controls: R-formula tamper and flow tamper both break the identity;
round-case specialization consistent; `‖Ric‖² > 0` strictly; collapse
boundary sampled and recorded as an obstruction-tower entry, not claimed.

## Honest remainder

Untouched, restated from the adapter's obstruction tower: surgery and
neck/cap data; non-collapsing, κ-solutions, canonical neighbourhoods;
singularity classification off the homogeneous carrier; finite-time
extinction in general; topological reconstruction; **and the Poincaré
conjecture itself — NOT CLAIMED**. Hamilton–Perelman stays a pinned named
dependency, cited never rederived (CIRC-1).

The homogeneous carrier is precisely where surgery and collapse are
absent. That is why this is a calibration.

## Reproduce

```bash
python papers/poincare-calibration/certificates/pc1_lambda_monotonicity.py
python -m pytest papers/poincare-calibration/tests -v
```
