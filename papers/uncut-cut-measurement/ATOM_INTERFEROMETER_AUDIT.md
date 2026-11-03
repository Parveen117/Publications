# v1.7: native source response meets an atom-interferometer benchmark

**Outcome:** an executable conditional source-to-phase interface, real published summary observations, and a numerical apparatus reference. The reference explicitly uses Newtonian gravity. **No native gravity prediction has yet been compared with data.** This is a retrospective audit; the observed values were read before the calculation and are not a held-out test.

## 1. Evidence intake

[Rosi et al. (2014)](https://arxiv.org/abs/1412.7954v1) measured a source-configuration double-difference phase of `0.547870(63) rad` using two rubidium clouds and 24 tungsten cylinders. Their nominal geometry and pulse timings are transcribed in [the input manifest](experiments/atom_interferometer_2014_inputs.json); the cylinder mass is approximately 516 kg. They attribute about 97% of the signal to cylinders, with the remainder from other moving apparatus. The [authors' detailed analysis](https://doi.org/10.1098/rsta.2014.0030) identifies the quoted phase error as a rescaled aggregate statistical error and reports a full simulation of `0.548105(1) rad` using historical CODATA G. That simulation's quoted error covers Monte Carlo sampling, not all systematic uncertainty.

[Rosi et al. (2015)](https://arxiv.org/abs/1501.01500v1) use three clouds and report `0.5533(6) rad` for a different differential observable, with a simulated `0.5528 rad`. Their gravity-field curvature is a spatial derivative of acceleration, not directly a Riemann tensor. The different geometry and observable prevent treating these two phase numbers as a train/test pair with one unchanged transfer factor.

The papers, publisher landing page and targeted data searches were inspected. No machine-readable shot table, complete source-mass metrology/covariance bundle or authors' full simulation implementation was located in this audit. This is a search outcome, not a claim that those data do not exist. Neither a figure digitization nor the published G estimate is treated as an independent raw observation.

## 2. AI-1: the native source-to-phase interface

U36 supplies a stationary field `phi_s=H^(-1)b_s` after **declaring** a positive grounded graph, weights and source map. Keep that graph fixed as the physical source configuration `s` changes; moving matter changes `b_s`. If a candidate changes the graph, use its separately specified `H_s` and do not apply the fixed-H simplification.

Supply a physical trajectory observation row `D_j(t)` and an independently calibrated factor `kappa`, so that the model's vertical acceleration is

`a_(j,s)(t)=kappa D_j(t) H^(-1) b_s`.

`D` includes the graph-to-laboratory position map and the signed spatial difference/gradient reading; `kappa` includes energy, source/probe and inertial calibration. This is the missing physical adapter, not a consequence of renaming a graph coordinate. Physical times, lengths, inertia and laser readout are explicitly added downstream; they are **not imported into the uncut primitive**.

For an ideal point-probe acceleration readout with instantaneous pulses at `0,T,2T`, define the conventional instrument sensitivity

`w_T(t)=t` for `0<=t<=T`, and `w_T(t)=2T-t` for `T<=t<=2T`.

The linearized sensor model gives `phi_(j,s)=k_eff integral w_T(t) a_(j,s)(t) dt + nu_(j,s)`, where `nu` denotes independently budgeted nuisance phase. This is a **supplied measurement law**. It is not derived here from the native Born-form construction. Finite pulses, recoil, cloud distributions and trajectory feedback need corrections or uncertainty bounds before precision comparison.

**Kinematic identity supporting this adapter.** If `z''=a`, integrate twice to obtain `z(t)=z(0)+t z'(0)+integral_0^t (t-u)a(u)du`. The combination `z(0)-2z(T)+z(2T)` cancels position and velocity and equals `integral_0^(2T) w_T(u)a(u)du`. A linear laser phase readout multiplies this by `k_eff`. This identity alone does not establish the full quantum interferometer dynamics.

With the signed cloud difference **lower minus upper**, form

`r_s=integral w_T(t)[D_(L,s)(t)-D_(U,s)(t)]dt`.

The measurable source contrast is therefore

`Phi_CF=k_eff kappa [r_C H^(-1)b_C-r_F H^(-1)b_F]+Delta nu`.

If independently justified reference paths give `r_C=r_F=r`, this reduces to

`Phi_CF=k_eff kappa r H^(-1)(b_C-b_F)+Delta nu`.

This is the explicit experiment-facing prediction interface. Source-independent common contributions cancel only under the stated shared-response assumptions. A change in trajectory or calibration must propagate through `r_s`, not disappear inside an adjusted gain. A fixed-H linear model also predicts proportional source response and sign reversal under C/F exchange, conditional on these controls.

The exact certificate checks the sensitivity moments `T², T³, 7T⁴/6`, common-response cancellation, source linearity, and a changed-trajectory control that defeats an unjustified common-row simplification. Its small rational matrices are algebraic controls, not experimental predictions.

## 3. Numerical apparatus reference without fitting the observed phase

To test the apparatus integration machinery, [the reference script](experiments/atom_interferometer_reference.py) computes the external Newtonian field of finite uniform cylinders. It integrates the signed vertical field over each cylinder with Gauss-Legendre volume quadrature, samples two fixed on-axis ballistic trajectories, and applies the triangular time sensitivity. Six cylinders are placed at each of two radial distances on each of two platforms. The native graph is not used in this reference.

Inputs are nominal and explicit: equal cylinder masses, uniform density, point clouds, instantaneous pulses, fixed background `9.81 m/s²`, nominal laser wavelength `780 nm`, and historical independent reference `G=6.67384e-11 m³ kg⁻¹ s⁻²`. No parameter or global amplitude is fitted to the observed phase. Cloud averaging, recoil, support masses, other moving masses and the full systematic ledger are omitted.

The signed C-minus-F, lower-minus-upper result is approximately **0.5304493433 rad**. It is about **3.18% below** the measured total. Increasing quadrature orders from 20 to 28 changes it by less than `1e-8 rad`; this establishes numerical stability for the simplified model, not a bound on physical model error. The omitted non-cylinder signal is substantial, and further differences remain from the approximations. We do **not** divide by the reported cylinder fraction to manufacture a fitted total or assign a significance using only the experimental statistical error.

This reference recovers the signal's scale from apparatus geometry but is not a precision reproduction, native-theory success, independent discovery, or new gravity law. Its main purpose is to make the input-to-observable calculation inspectable before inserting a native field.

## 4. Exact stopping point and next validation gate

The native field cannot yet replace the external reference because a physical graph/weight/boundary selection, mass-source deposition, and force-to-acceleration calibration have not been supplied independently of this signal. The original source-coupled law is still a candidate. A single aggregate scalar cannot both fit its free amplitude and validate that amplitude.

The next comparison must freeze those inputs, include the apparatus corrections with covariance, and reserve additional source configurations or pulse times before inspecting their outcomes. Evaluate both a conventional baseline and the native prediction using the same geometry, instrument response and nuisance budget. A disagreement then tests the declared native model and adapter jointly. A native match supports that bounded physical model; it does not uniquely select all of nature.

Current status is **NOT_EVALUATED_NATIVE_MODEL**. Public summary values have been ingested, but no reserved measurements, raw-shot reanalysis or native empirical PASS exists. No outside researcher was contacted and no hardware run was submitted.

## Reproduction and evidence separation

From `papers/uncut-cut-measurement`:

```bash
python certificates/atom_interferometer_bridge.py
python experiments/atom_interferometer_reference.py --check
python certificate.py --check
```

The reference script requires NumPy and uses floating-point quadrature. Its committed [result](experiments/ATOM_INTERFEROMETER_REFERENCE.json) is explicitly separate from the exact rational certificate. The master certificate checks the rational adapter controls and binds all new file bytes; it does not silently promote floating-point reference calculations into exact physical certification. Run the reference check separately to reproduce the numerical calculation.
