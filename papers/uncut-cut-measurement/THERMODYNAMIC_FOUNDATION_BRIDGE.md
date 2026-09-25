# From native source cost to the thermodynamic compass

Status: **conditional derivation with a written proof and an independent exact finite witness**. This is an equilibrium bridge on a cut. It does not assign entropy or volume to the uncut ground, or derive a physical clock, the gravitational action, or quantum measurement.

## The typed foundation

The starting point is U36's grounded connected weighted graph: after fixing boundary values, its scalar aligned carrier has a positive definite stiffness matrix `K` and source `b`. The **candidate cost law**, chosen in U36 rather than implied by the word “uncut”, is

`E(x;b) = x^T K x / 2 - b^T x`.

Choose two independently readable carrier coordinates `y=(s,v)` and write `x=(y,z)`, with `z` hidden by the cut. The physical assertions that calibrated `s` is entropy `S`, calibrated `v` is volume `V`, and measured energy is `a E+c` for `a>0` are **adapter assumptions**. We initially set `a=1` in native units. A physical equilibrium chart also requires `S,V>0`, `T>0`, and the appropriate preparation. Neither `s` nor `v` is a property of the literal uncut primitive before this identification.

### Theorem TF-1 (exact equilibrium projection)

Let `K` be symmetric positive definite, with block partition `K_yy,K_yz,K_zz`, and partition `b=(b_y,b_z)`. Define the constrained equilibrium cost `U(y)=min_z E((y,z);b)`. Then

`z*(y)=K_zz^{-1}(b_z-K_zy y)`,

`U(y)=y^T D y/2 - d^T y - b_z^T K_zz^{-1} b_z/2`,

`D=K_yy-K_yz K_zz^{-1} K_zy`, `d=b_y-K_yz K_zz^{-1} b_z`, and `D` is positive definite. The full stationary minimizer and the minimizer of `U` have the same visible coordinates. The result is invariant under an invertible change of hidden basis. It is **constrained equilibrium elimination**, not a closed dynamics for `y`.

**Proof.** Complete the square in `z`:

`E(y,z)=U(y)+(z-z*(y))^T K_zz (z-z*(y))/2`.

Positivity of `K_zz` gives unique `z*(y)`. For `y≠0`, apply positivity of `K` to `(y,-K_zz^{-1}K_zy y)`: its quadratic form is `y^T D y>0`, proving `D>0`. Minimizing first over `z` and then over `y` equals unconstrained minimization. A hidden-basis change changes coordinates of the minimizing `z` but leaves its minimum value unchanged. □

### Corollary TF-2 (thermodynamic chart, after calibration)

On any open region where calibrated `S,V` are valid independent positive coordinates and `T:=U_S>0`, define `P:=-U_V`. Then

`dU=T dS-P dV`, `dT∧dS=dP∧dV`, `T_V=-P_S`,

and the three other Maxwell relations follow on the locally regular Legendre charts. The response Hessian is exactly `D` in these units. The thermodynamic compass is therefore a chart representation of the **same projected cost**. In particular, loop integrals of `dU` vanish on this chart; a nonzero response or memory loop requires another, explicitly specified response one-form or an extended state/ledger. No irreversible entropy production follows from TF-1.

**Proof.** Differentiate `U`; `D=D^T` gives `U_SV=U_VS`. Positive definiteness implies `U_SS>0`, `U_VV>0` and `det D>0`, so the local Legendre changes are regular. Pulling back the thermodynamic contact form `dU-T dS+P dV` gives zero. The remaining Maxwell relations follow from the exact transformed potentials. □

### Boundary of the derivation

This derivation takes a stationary cut of a **supplied** finite graph and supplied quadratic source cost. Physical entropy/volume/energy units, temperature and pressure measurements, an equation of state beyond the positive quadratic sector, and universality across preparations are not consequences of positive definiteness. The conservative affine Cayley continuation in U36 preserves excess cost; it does **not** force a disturbed field to relax to the minimum. Stationary preparation or a reservoir is an independent ingredient. A physical metric and Einstein-type equations require the further propagation symbol, action and identification hypotheses in the thermodynamic paper. A Gaussian quantum theory there additionally uses its stated `hbar`, action and positivity hypotheses.

## Exact rational witness and falsification controls

Three free graph vertices form a triangle with unit links; ground vertex 1 with stiffness 1 and vertex 2 with stiffness 2. In the order `(S,V,z)`, take source `b=(0,10,0)` and

`K=[[3,-1,-1],[-1,4,-1],[-1,-1,2]]`.

It is a U36 grounded graph with positive weights, not an arbitrary matrix. Eliminating `z` gives

`D=[[5/2,-3/2],[-3/2,7/2]]`, `d=(0,10)`, `det D=13/2>0`,

`U(S,V)=5S²/4-3SV/2+7V²/4-10V`.

At `(S,V)=(1,1)`, `z*=1`, `T=1`, `P=8`; an open neighborhood has positive `S,V,T,P`. Exactly,

`C_V=2/5`, `C_P=7/13`, `K_S=7/2`, `K_T=13/5`,

so `C_P/C_V=K_S/K_T=35/26`. Here `K_S=-V(∂P/∂V)_S` and `K_T=-V(∂P/∂V)_T`. The response equalities are consequences of this specified `U`; they are consistency checks, not independent observations of nature.

The executable [exact witness and negative controls](certificates/thermodynamic_foundation_bridge.py) check the graph decomposition, rational Schur complement, response identities, hidden-coordinate invariance, and a perturbed nonsymmetric response map whose Maxwell identity fails. The general proof above does not depend on enumeration.

## Source and dependency boundary

- U36 and W51: [native graph source law and stationary-only hidden elimination](SOURCE_RESPONSE_AND_GRAVITY_TESTS.md); U28 separately governs dynamical closure.
- [Thermodynamic paper, equilibrium foundation](https://github.com/Parveen117/Publications/blob/a3ff8e9f8e6afdea44c55557346e8b5fea0c5d22/papers/thermodynamic-response-corrections/source/sections/01_foundations.tex): starts with an admissible `U(S,V)` and derives contact/Maxwell/response results. TF-1 supplies **one explicit conditional constructor** for that starting object.
- [Thermodynamic action lift](https://github.com/Parveen117/Publications/blob/a3ff8e9f8e6afdea44c55557346e8b5fea0c5d22/papers/thermodynamic-response-corrections/source/sections/06b_classical_field_equations.tex): its action is separately declared, so TF-1 does not silently imply its metric field equations.
- [Thermodynamic spacetime identification](https://github.com/Parveen117/Publications/blob/a3ff8e9f8e6afdea44c55557346e8b5fea0c5d22/papers/thermodynamic-response-corrections/source/sections/06a00_response_propagation_metric.tex): Lorentzian principal symbol/scale are further data, not properties of the positive equilibrium Hessian.

This bridge proves a mathematical *route* from U36's native cost to the thermo paper's equilibrium input. It is the foundation on which a physical measurement adapter can now be tested without importing spacetime into the native start.
