# NSB1 — Return to the Yang–Mills programme

The physics-first development establishes a concrete reason to retain
nonlinear marked response: spectra of a compressed nonlinear representation
can be exact while omitting the target. It does not say that an exact
Hamiltonian spectrum fails to determine its own gap.

## What can now be used

1. **Choose the target before choosing the spectrum.** A linearized mode,
   a finite jet spectrum, a selected correlation function, and the full
   physical Hamiltonian gap are different targets.
2. **Retain marked nonlinear coefficients.** Coordinate response markers
   reconstruct a declared finite jet. In the two-channel response model,
   spectra plus calibrated impedance and bias derivatives reconstruct
   the response geometry.
3. **Keep a quantitative margin.** Local algebraic recovery is insufficient
   if its conditioning degenerates with volume, cutoff or proximity to
   a blind state.
4. **Test for invisible nonlinear barriers.** The exact compact-gauge
   family below proves that arbitrarily deep local jets alone cannot
   establish a coupling-uniform gap.

## The new actual compact-gauge theorem

On one SU(2) plaquette, let

\[
 s=1-W^2,\quad f_N=1-s^{N+1},\qquad
 K_{t,N}=L_E-\frac t2L_Ef_N+\frac{t^2}{4}\Gamma(f_N).
\]

The vacuum is exactly proportional to exp(t f_N/2), with energy zero.
For every prescribed finite local degree r, choose 2N>r.
In exponential link coordinates at the trivial configuration,
the multiplication potential has zero derivatives through degree r
for every t>=0. Thus these local operator jets agree with the free
electric model.

But the physical gauge-invariant gap obeys

\[
 \Delta_{t,N}\le
 162\exp\{-t[(3/4)^{N+1}-(1/4)^{N+1}]\}\longrightarrow0.
\]

The proof and assumptions are in Section 7 of
[NSB1](NSB1_NONLINEAR_SPECTRAL_BLINDNESS.md).
This extends YM82's quartic parent obstruction while matching
arbitrarily many local jets. It is a separate polynomial model;
overlapping active plaquettes do not have this asserted vacuum.

## The admissible next target for the original model

For the actual ground state h of the original gauge Hamiltonian,
YM82's conditional method requires estimates on

\[
 S=2\log h,\qquad
 \rho_B=\hbox{a uniform conditional Poincare lower bound},\qquad
 c_{BC}\ge\|\operatorname{Hess}_{BC}S\|_{\rm op}.
\]

A positive comparison matrix
diag(rho_B)-(c_BC) would prove a full gap under the stated hypotheses.
NSB1 explains why pointwise spectral data or a finite local response
jet do not automatically supply those estimates. Candidate markers
must be shown to control the actual vacuum quantities, not a chosen
comparison potential.

YM81's actual source-resolved moments and continuum-survival targets
remain relevant. A source must retain nonzero norm in a fixed physical
energy window; an observed microscopic band that escapes as the cutoff
is removed does not do that.

No new estimate for rho_B or c_BC along the homogeneous weak-coupling
continuum trajectory is proved in NSB1. Nontrivial continuum fields,
physical scale control, covariance, and a positive continuum mass gap
remain obligations. The next work should address one of those actual
quantities using the repaired observations, rather than infer a gap
from response curvature alone.
