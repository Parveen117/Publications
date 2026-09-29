# Native curvature action and complete thermodynamic metric variations

Author: Monty Dabas. Development R3. Results MG-1–MG-8.

**Result.** In a declared four-dimensional polynomial trace-action class, the
admitted native coefficient representation reduces the bulk gravitational
action to Palatini–Holst plus a cosmological term. For an invertible coframe,
nonzero Palatini coefficient, real Holst coefficient and no connection-coupled
matter, independent coframe/connection variations give vacuum Einstein
equations. An explicit local construction using sixteen stable thermo channels
supplies every coframe variation, hence every metric variation. Adding the
declared Maxwell kinetic law gives the source-free Einstein–Maxwell system,
with electromagnetic stress sourcing the metric.

This is a conditional action-class reduction and a complete variation interface.
The continuum, its dimension, the representation, action class, constitutive
extraction and coupling constants are inputs. It does not establish that native
primitive axioms uniquely select general relativity, or compute G, c, Λ or α.
Palatini–Holst equivalence and the role of topological terms are established
results [1,2]. The contribution here is their explicit native-trace realization,
connection to the CP thermo chart, and reproducible variation/counterexample
controls. No priority claim is made for those established gravity results.

All statements concern local bulk equations and compactly supported smooth
variations. Boundary charges, global bundles and boundary conditions are outside
the theorem. Differentiating curvature requires the usual smoothness; no
low-regularity or Cauchy existence theorem is inferred from these calculations.

## MG-1. Native coframe, connection and trace identities

Use the unchanged real coefficient representation of NP-1. Its four matrices
Γ_a and native volume element J obey

\[
\{\Gamma_a,\Gamma_b\}=2\eta_{ab}I,\quad
\eta=\operatorname{diag}(-1,1,1,1),\quad
J=\Gamma_0\Gamma_1\Gamma_2\Gamma_3,\quad J^2=-I.
\]

J anticommutes with vectors and commutes with bivectors. The trace functional
is τ(A)=tr_R(A)/4. Fix ε_0123=+1. The Clifford relations give

\[
\tau(\Gamma_a\Gamma_b)=\eta_{ab},\qquad
\tau(\Gamma_a\Gamma_b\Gamma_c\Gamma_d)
=\eta_{ab}\eta_{cd}-\eta_{ac}\eta_{bd}+\eta_{ad}\eta_{bc},
\]
\[
\tau(J\Gamma_a\Gamma_b\Gamma_c\Gamma_d)=-\epsilon_{abcd}.
\tag{1}
\]

For repeated indices, anticommutation reduces the ordinary trace to the first
formula and the J trace vanishes. For four distinct indices, their product is
ε_abcd J, whose J trace is −ε_abcd. This proves (1) for all indices.
The six Γ_aΓ_b, a<b, span the Lorentz bivectors, since

\[
[\Gamma_a\Gamma_b,\Gamma_c]
=2(\eta_{bc}\Gamma_a-\eta_{ac}\Gamma_b).
\]

On an admitted oriented four-manifold, let e^a be an invertible coframe and
ω^{ab}=−ω^{ba} an **independent Lorentz connection**. Define

\[
E=e^a\Gamma_a,\quad \Omega=\frac14\omega^{ab}\Gamma_a\Gamma_b,
\quad F_\Omega=d\Omega+\Omega\wedge\Omega
=\frac14 R^{ab}\Gamma_a\Gamma_b,
\]
\[
T=DE=dE+\Omega\wedge E+E\wedge\Omega=T^a\Gamma_a,
\quad R^{ab}=d\omega^{ab}+\omega^a{}_c\wedge\omega^{cb}.
\tag{2}
\]

Repeated internal indices run over all four values, including both orders of
antisymmetric pairs. The spacetime metric is q=η_ab e^a⊗e^b, with
det q=−(det e)^2 and volume vol_e=e^0∧e^1∧e^2∧e^3. Choose the orientation
det e>0, so vol_e=√(−det q) d^4x. This produces Lorentz signature from the
admitted Clifford/coframe realization; it does not obtain Lorentz signature by
turning a pseudo-Kähler tangent metric into a Lorentzian one. GS-1's even-index
obstruction remains intact. The finite coefficient module is still not a
faithful representation of a proper one-sided seam.

Under local proper spin transformations, E and F_Ω transform by conjugation,
and J is invariant. The traces below are thus local Lorentz invariants.

## MG-2. A finite action-class reduction

Declare the action class to consist of real **constant** linear combinations
of τ(W) and τ(JW), where W is a wedge-word of total differential-form degree
four in E (degree one), T (degree two) and F_Ω (degree two). Fixed insertions
are only I and J. There are no inverse coframes, Hodge stars, additional fields,
field-dependent coefficients, extra internal tensors or further derivatives in
this gravitational class. This is an assumption about admissible actions,
not a consequence of the coefficient algebra alone.

There are exactly eleven words:

\[
EEEE;\quad EET,ETE,TEE;\quad EEF,EFE,FEE;\quad TT,TF,FT,FF.
\]

Here F abbreviates F_Ω. Odd Clifford traces remove EET, ETE, TEE, TF and FT,
with either insertion. Formula (1) and graded trace cyclicity leave

| Trace | Differential form |
|---|---|
| τ(JE^4) | −ε_abcd e^a∧e^b∧e^c∧e^d = −24 vol_e |
| τ(E^4) | 0 |
| τ(JE∧E∧F_Ω) | −(1/4) ε_abcd e^a∧e^b∧R^{cd} |
| τ(E∧E∧F_Ω) | −(1/2) e^a∧e^b∧R_ab |
| τ(T∧T) | T_a∧T^a |
| τ(JT∧T) | 0 |
| τ(F_Ω∧F_Ω), τ(JF_Ω∧F_Ω) | The two curvature characteristic forms |

In particular τ(EFE)=−τ(EEF), but **τ(JEFE)=+τ(JEEF)**:
moving a vector past J introduces an extra minus sign. Neither ordering gives
a new invariant. Multiple allowed J insertions, if included, reduce to the
same list using JE=−EJ, JT=−TJ, JF=FJ and J²=−I.

The covariant product rule and D²e^a=R^a{}_b∧e^b prove the Nieh–Yan identity

\[
d(e_a\wedge T^a)=T_a\wedge T^a-e^a\wedge e^b\wedge R_{ab}.
\tag{3}
\]

Thus the torsion-square term changes the coefficient of the Holst form modulo
an exact form. It is not equal to that form pointwise. Also DF_Ω=0 and DJ=0
give, for Q=I or J,

\[
\delta\tau(QF_\Omega\wedge F_\Omega)
=2d\tau(Q\,\delta\Omega\wedge F_\Omega).
\]

These constant-coefficient characteristic terms do not change the local bulk
equations. Their integrals need not vanish, and may affect global or boundary
physics. Modulo these terms, the whole declared class therefore reduces to
Palatini, Holst and volume. If the Palatini coefficient is nonzero, write it as

\[
\boxed{
S_g=-\frac1\kappa\int\tau[(J+\gamma I)E\wedge E\wedge F_\Omega]
+\frac\Lambda{24\kappa}\int\tau(JE^4),
\qquad \kappa\ne0,\quad \gamma\in\mathbb R.
}
\tag{4}
\]

Our γ is the coefficient of the Holst term relative to the displayed Palatini
normalization; it is not the usual reciprocal-parameter convention. Equation
(4) is, in real tensor notation,

\[
S_g=\frac1{4\kappa}\int\epsilon_{abcd}e^a\wedge e^b\wedge R^{cd}
+\frac\gamma{2\kappa}\int e^a\wedge e^b\wedge R_{ab}
-\frac\Lambda{24\kappa}\int\epsilon_{abcd}e^a\wedge e^b\wedge e^c\wedge e^d.
\tag{5}
\]

The algebra fixes the displayed conversion factors. It does not fix κ, γ or Λ,
their physical units, or whether nature selects this action class.

## MG-3. Connection variation removes torsion

Since δF_Ω=DδΩ, integration by parts in (4) gives the connection equation
D[(J+γI)E∧E]=0. The trace pairing on the six bivectors is nondegenerate;
multiplication by J preserves that space. For every real γ,

\[
(J+\gamma I)^{-1}=\frac{\gamma I-J}{1+\gamma^2}.
\]

Hence D(E∧E)=0, or, in components,

\[
e^a\wedge T^b-e^b\wedge T^a=0.
\tag{6}
\]

This implies T^a=0 when e is invertible in dimension four. For completeness,
let ι_b contract with the dual coframe vector, and put V=Σ_b ι_bT^b. Contract
(6) with ι_b and sum over b. Using Σ_b e^b∧ι_bT^a=2T^a gives
(3−4)T^a−e^a∧V=0, so T^a=−e^a∧V. Contract again to obtain V=−3V.
Thus V=0 and T^a=0. Because ω is Lorentz-valued, it is now the unique
Levi-Civita spin connection of e. No torsion-free condition was imposed before
varying the action.

The hypotheses matter: a degenerate coframe has no dual frame; a complex
γ=±i makes the multiplier singular; connection-coupled spin matter changes
the right-hand side. None of these cases is included in this theorem.

## MG-4. Complete coframe variation gives Einstein equations

The Holst contribution to the coframe equation is proportional to
e^b∧R_ab. It vanishes when T=0 by the first Bianchi identity. Varying the
other terms in (5) gives

\[
\epsilon_{abcd}e^b\wedge
\left(R^{cd}-\frac\Lambda3e^c\wedge e^d\right)=0.
\tag{7}
\]

For an invertible coframe and its torsion-free connection this is equivalent to

\[
\boxed{G_{\mu\nu}(q)+\Lambda q_{\mu\nu}=0.}
\tag{8}
\]

One normalization check is ε_abcd e^a∧e^b∧R^{cd}=2R(q) vol_e; (5) consequently
reduces to (1/2κ)∫(R−2Λ)vol_e after solving the connection equation. More
explicitly, the coframe Euler density is

\[
\mathcal E_a{}^\mu=-\frac{\sqrt{-q}}\kappa
\bigl(G^{\mu\nu}+\Lambda q^{\mu\nu}\bigr)\eta_{ab}e^b{}_\nu.
\tag{9}
\]

There are sixteen coframe components and ten metric components. The variation

\[
\delta q=e^T\eta\,\delta e+\delta e^T\eta e
\]

is onto all symmetric tensors. For any symmetric h, an explicit right inverse
is δe=(1/2)e q^{-1}h. Its kernel consists of δe=L e with L^Tη+ηL=0,
the six Lorentz gauge directions. Thus (9) really tests the full Einstein
tensor, rather than a trace or another projection. The determinant also fixes
the volume density used here; it is no longer an independent background input
in this metric/gauge action.

## MG-5. A stable thermo realization with all metric variations

Use sixteen independent copies, indexed by (a,μ), of the CP-5 thermo potential

\[
U(s,v)=10s-v+s^2+sv+\frac52v^2+\frac12s^2v,
\quad H=\begin{pmatrix}2+v&1+s\\1+s&5\end{pmatrix},
\quad \Delta=5(2+v)-(1+s)^2.
\]

Here s,v are local dimensionless coordinates. Physical scales and positive
reference S,V are supplied if a material interpretation is wanted. Around
(s,v)=(0,0), H is positive definite and the reference temperature and pressure
are positive. In this stable patch the phase connection has the CP-5 chart

\[
\alpha=P\,dQ+d\chi,\quad Q=s,\quad
P=\frac16-\frac1{2\sqrt\Delta},\quad
P_v=\frac5{4\Delta^{3/2}}>0,
\tag{10}
\]

with χ=arcsin[(1+s)/√(5(2+v))]−s/6. Its inverse is

\[
s=Q,\qquad v=\frac{(1+Q)^2+[4(P-1/6)^2]^{-1}}5-2,
\qquad P<1/6.
\tag{11}
\]

**Declare the constitutive extraction**

\[
e^a=\sum_{\mu=0}^3(\delta^a_\mu+P_{a\mu})\,dQ_{a\mu}.
\tag{12}
\]

This chooses the Darboux part α−dχ, its normalization and an offset. Removing
dχ from a coframe is not a U(1) gauge invariance of a metric; the choice in (12)
is an explicit additional constitutive rule. The claim is existence of a
complete local thermo adapter, not a unique gauge-independent λ-to-spacetime
law forced by the Hessian.

To realize any smooth oriented coframe sufficiently close to the identity in
a coordinate patch, choose

\[
Q_{a\mu}(x)=x^\mu,\qquad P_{a\mu}(x)=e^a{}_\mu(x)-\delta^a_\mu,
\tag{13}
\]

and recover v_{aμ} by (11). The patch can be centered so x=0, e=I; continuity
keeps all sixteen channels inside the positive-temperature, positive-pressure,
strictly stable neighborhood. An arbitrary smooth nondegenerate coframe can
be put in this form near a point by a coordinate/calibration choice. This is
a local representation statement, with neither a global extension nor a
preferred physical calibration asserted.

Hold all Q fields fixed in (13) and vary the sixteen v fields independently.
Then, at every point of this realization,

\[
\delta e^a{}_\mu=P_{v,a\mu}\,\delta v_{a\mu}.
\tag{14}
\]

All diagonal coefficients are strictly positive. Small compactly supported
variations preserve the strict inequalities and coframe invertibility. If a
coframe action has first variation ∫𝓔_a{}^μ δe^a{}_μ, its v equations on this
chart are exactly P_v,aμ 𝓔_a{}^μ=0, equivalent to all sixteen coframe equations.
Conversely, when all coframe equations vanish, the remaining thermo variations
also vanish after integration by parts. More generally the same implication
holds whenever the four dQ_{aμ} are independent for each fixed a, an open rank
condition. The independent ω equations are retained unchanged.

The thermo-to-coframe derivative has rank sixteen and its composition with
coframe-to-metric has rank ten. Consequently the pullback of (4) through (12)
has the full Einstein equations on these regular patches. Sixteen channels
are a **sufficient construction**, with no minimality claim. The independent
spin connection and action class still require physical selection.

## MG-6. Three controls against incomplete gravity claims

**Coordinate-only coframe.** If e^a=dX^a for four independent scalar functions,
q=X*η is locally flat. Such an ansatz cannot represent a curved metric.
Variations of X are coordinate variations and test differential identities,
not all ten independent metric equations. For example X^1=(x^1)^2 on x^1>0,
with the other coordinates unchanged, gives a nonconstant but flat metric.

**Conformal-only variations.** For q=e^{2f}\bar q, δq=2q δf. The vacuum metric
action therefore yields only R=4Λ. The explicit metric

\[
q=t^2(-dt^2+dx^2+dy^2+dz^2),\qquad t>0,
\]
\[
\operatorname{Ric}_{\mu\nu}=t^{-2}\operatorname{diag}(3,1,1,1),\qquad R=0
\]

satisfies this trace equation at Λ=0 but fails the full vacuum Einstein
equations. This is a direct real-coordinate curvature calculation, not an
inference from a numerical fit.

**Pure Holst exception.** If the Palatini coefficient is zero, it cannot be
normalized as (4). A nonzero pure Holst action without a cosmological term
still imposes T=0, but its coframe equation then reduces to the first Bianchi
identity; it does not impose Einstein equations. Thus nonzero Palatini coupling
is substantive. Degenerate coframes provide a separate failure: the tested
rank-three frame has a torsion-equation map of rank eighteen instead of
twenty-four, so the torsion conclusion itself fails.

## MG-7. Dynamical metric with the Maxwell sector

Add the previously declared Maxwell law, now using q(e):

\[
S_{\rm EM}=-\frac1{4e_{\rm em}^{\,2}}\int\sqrt{-q}\,
F_{\mu\nu}F^{\mu\nu}\,d^4x,\quad F=da,\quad e_{\rm em}^{\,2}>0.
\tag{15}
\]

The gauge potential and coframe have independent variations. They can instead
be parametrized by **independent** thermo-channel collections: (12) for e, and
the regular four-channel CP-8 construction for a. The latter recovers all
Maxwell equations even at F=0; the gauge Noether divergence identity holds
off shell in this source-free sector. Sharing state variables without proving
the joint variation rank would not justify this conclusion.

Since (15) does not depend on ω, MG-3 still gives T=0. With the definition
δS_EM=(1/2)∫√(−q) T^{μν}δq_μν, the full equations are

\[
\boxed{G_{\mu\nu}+\Lambda q_{\mu\nu}=\kappa T^{\rm EM}_{\mu\nu}},
\quad
T^{\rm EM}_{\mu\nu}=\frac1{e_{\rm em}^{\,2}}
\left(F_{\mu\rho}F_\nu{}^\rho-\frac14q_{\mu\nu}F_{\rho\sigma}F^{\rho\sigma}\right),
\tag{16}
\]
\[
\partial_\mu(\sqrt{-q}\,F^{\mu\nu})=0,\qquad dF=0.
\]

These are source-free electromagnetic equations coupled to gravity, not a
claim to have coupled the NP commuting matter field to an arbitrary curved
spin connection. That matter action's local Lorentz covariance, stress and
spin-current/torsion variation need a separate treatment. The gauge kinetic
law, coupling and common metric choice in (15) remain constitutive input.

## MG-8. Verification and the remaining physical selection

`metric_dynamics.py` implements rational matrix-valued exterior forms, independent
index contractions, first jets and action variations. The source-bound R3
certificate records:

- 512 fourth-trace entries, 24 spin/vector commutators and all 22 degree-four
  trace choices on three coframe/curvature/torsion fixtures;
- a directly differentiated Nieh–Yan identity on three nontrivial first-jet
  fixtures, including nonzero exact-form density controls;
- complete metric-variation rank ten, Lorentz kernel dimension six and torsion
  rank twenty-four on three invertible frames; thirty explicit metric right
  inverses and a degenerate-frame rejection;
- 432 exact native curvature-action coframe variations against the independent
  Einstein-tensor contraction, with zero, positive and negative real Holst
  coefficients and both Einstein and non-Einstein curvature fixtures;
- 48 independently differentiated Maxwell coframe variations against its stress
  tensor, using exact first-order jets through determinant and matrix inverse;
- sixteen explicit stable thermo inverses and variation ranks sixteen/ten;
- the conformal trace-only counterexample at three rational times and a
  nonlinear coordinate-only flat control.

The gravity controls include a nonzero algebraic Weyl tensor with Ricci=0 and
R_abcd R^{abcd}=48. It is a pointwise curvature witness, not a new global vacuum
solution. Centered coframe differences in the polynomial gravitational density
are exact, because each particular coframe entry occurs at most linearly.
Maxwell's rational dependence on the metric is differentiated with jets rather
than a finite-difference approximation. Finite controls support the written
proofs above; they do not establish a universal theorem by sampling.

This step replaces an unspecified λ-to-coframe variation gap with a concrete
conditional adapter and derives the Einstein sector of an explicit native trace
action class. The next physical question is why a native process selects that
class, this independent connection, this adapter and its coupling scales.
Admitting inverse-coframe/Hodge curvature terms, variable coefficients or
additional fields broadens the class and can change the equations. The earlier
`06a1_dynamical_nonselection.tex` result is therefore preserved, not evaded.
G, c, Λ and α are not numerical predictions of this construction.

## Primary references

1. S. Holst, “Barbero's Hamiltonian derived from a generalized Hilbert-Palatini
   action,” *Physical Review D* **53**, 5966–5969 (1996).
   [arXiv:gr-qc/9511026](https://arxiv.org/abs/gr-qc/9511026),
   [doi:10.1103/PhysRevD.53.5966](https://doi.org/10.1103/PhysRevD.53.5966).
2. D. J. Rezende and A. Perez, “Four-dimensional Lorentzian Holst action with
   topological terms,” *Physical Review D* **79**, 064026 (2009).
   [arXiv:0902.3416](https://arxiv.org/abs/0902.3416),
   [doi:10.1103/PhysRevD.79.064026](https://doi.org/10.1103/PhysRevD.79.064026).

The source-pinned local dependencies are GS/NP, CP-5–CP-8, the corrected
spacetime/action datum and classical field equations, and the dynamical
nonselection theorem. This note adds a concrete action-class calculation to
those conditional results; it does not silently strengthen their physical
premises.
