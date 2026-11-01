# Measuring the electromagnetic response diagnostic

The immediate experimental target is a **map of nonlinear constitutive
response over two bias controls**. The map can be useful without claiming
that an optical beam automatically follows the selected native connection.

## Direct response route

1. Declare a reversible preparation and a local energy-density convention.
   Record units, temperature, mechanical clamping, frequency window, and
   internal-state preparation. Use independent D_y and B_z controls or a
   calibrated conversion from available controls.
2. At each bias, estimate a=∂e/∂x, b=∂e/∂y=∂h/∂x, and c=∂h/∂y.
   Check a>0 and ac-b²>0. Check the mixed-response agreement rather than
   forcing it in a fit. Across the map check a_y=b_x and b_y=c_x, the
   additional local integrability conditions for an energy Hessian.
3. Compute r=b/√(ac) and Z=√(a/c), then estimate their derivatives in
   both bias directions. Evaluate the two-form in NSR1-T2 with orientation
   +x,+y,-x,-y fixed in advance.
4. Compare against a constant-H quadratic control fitted at one bias.
   This control matches that point's linear response but predicts f=0.
   A significant nonzero reconstructed f distinguishes higher response
   structure within the stated model.

These are measurements proposed for future work. No acquired data are
included in this folder. Near |r|=1 the reconstruction is poorly conditioned;
a large inferred curvature there must be accompanied by uncertainty on r
and its derivatives, not interpreted as an unqualified divergence signal.

## Wave route in a verified closed polarization sector

Measure both signed characteristic velocities and differential impedance
at each bias. Normalize by v_*=U_*/(D_*B_*) and impedance scale B_*/D_*.
Equation (11) reconstructs a,b,c. Then follow the same two-bias derivative
procedure. The two directions of propagation are not two polarization
measurements. Coupling to other polarization components invalidates the
two-by-two reconstruction and requires the full response matrix.

Born–Infeld provides a theoretical validation fixture, not a supplied
material fit. At x=y=1/2 its exact values are:

| Quantity | Value |
| --- | --- |
| a, c | 4/5, 4/5 |
| b | 1/5 |
| relaxed/clamped gain | 16/15 |
| r, Z | 1/4, 1 |
| normalized signed speeds | 1, -3/5 |
| r_x, r_y | 1/2, 1/2 |
| (log Z)_x, (log Z)_y | -4/5, 4/5 |
| f_xy | -128/(75√15) |
| Gaussian response curvature κ | 128/225 |

Z equals one at this point, while its derivatives are nonzero. Both bias
directions are therefore essential. A zero-bias test alone also misses this
model's curvature because symmetry forces it to vanish at the origin.

## Distinguishing physical mechanisms

| Comparison | What it can reveal | Additional requirement |
| --- | --- | --- |
| One-bias linear fit versus two-bias derivative map | higher constitutive response absent from the linear fit | calibrated independent bias derivatives |
| Clamped versus quasistatically relaxed internal mode | static Schur contribution to response curvature | resolved equilibrium branch and relaxation time |
| Low versus higher probe frequency | possible retarded/internal-mode response | a causal dispersive constitutive model; NSR1 alone is insufficient |
| Reconstructed ∫f versus measured return of a probe | whether physical probe transport realizes the selected connection | derive or calibrate the probe's actual transport equation |

For elasticity, spin response, or another symmetry sector, first identify
the physical potential, its conjugate variables and their units, and a
stable response patch. NSR1-T4 then determines the symmetry restrictions on
its response tower. A named symmetry by itself does not select the potential
or prove that this electromagnetic adapter applies.

For Yang–Mills, the relevant next interface would be an actual family of
gauge-invariant source susceptibilities with its scale dependence retained.
The classical electromagnetic calculations here neither extend the YM62
coupling window nor control the YM64 energy-dependent Schur operator.
