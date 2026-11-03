# Thermodynamic response, native phase charge and complete gauge variations

Author: Monty Dabas. Development edition: 30 September 2026.

This extension joins [NP-1–NP-6](NATIVE_PROPAGATION_ACTION.md) to the
[thermodynamic phase connection PF-1–PF-7](../thermo-phase-field/THEOREM.md).
The main advance is an exact criterion for recovering the **full sourced
Maxwell equation by varying thermodynamic state maps**, rather than treating
the gauge potential as an independently variable primitive. An explicit
curved thermodynamic potential realizes the required local coordinates.

The coefficient cone and source current are derived inside the admitted native
representation. The four-dimensional continuum, propagation coefficients,
background densities, field statistics and quadratic gauge action remain
specified model data. In particular the thermodynamic maps below enter the
connection only: the background metric and densities are held fixed under
their variations. Their own dynamics, and a λ-to-spacetime-coframe law, remain
separate targets.

## CP-1. A Lorentz cone in the native positive response sector

Retain the real four-dimensional representation, A_i and J of NP-1–NP-3.
Write τ(X)=tr_R(X)/4. The subspace

\[
\mathcal H_J=\{X:X^T=X,\ [X,J]=0\}
=\operatorname{span}_{\mathbb R}\{I,A_1,A_2,A_3\}
\]

has real dimension four. One proof identifies ψ=(p,q,r,s) with
ζ=(p+iq,r−is): J acts as multiplication by i and A_i as the three Pauli
matrices. J-linear self-adjoint maps are precisely two-by-two complex
Hermitian matrices. Equivalently the same dimension and basis follow by
solving the real symmetry and commutator equations.

For X=tI+x_iA_i define the native quadratic response

\[
\boxed{\Delta_J(X)=2\tau(X)^2-\tau(X^2)=t^2-|x|^2.}
\]

The eigenvalues of X are t±|x|, each repeated twice in the real
representation. Thus det_R X=Δ_J(X)², and

\[
X>0\iff t>|x|,\qquad
X\ge0\iff t\ge|x|.
\]

Consequently Δ_J has signature (+,−,−,−); positivity selects its future
component. The order X≼Y defined by Y−X≥0 is the corresponding cone order.
This is a Lorentzian form on the **space of admissible response operators**,
not a Lorentzian Kähler tangent metric. It therefore respects GS-1. Nor does
this algebraic statement identify operator coordinates with spacetime events.

The larger real EMK coefficient space without this J-compatible self-adjoint
gate has a different determinant: for aI+bR+cK+dL in one block it is
a²+b²−c²−d², of signature (2,2). The positive response gate and selected
two-block representation are essential hypotheses of the Lorentz-cone result.
This is the established Hermitian-matrix realization of the Lorentz cone,
expressed and checked in the native coefficient basis; its general priority
is not claimed.

## CP-2. A stable thermo Hessian is a sum of two native null responses

For a real field ψ define

\[
n=\psi^T\psi,\qquad b_i=\psi^TA_i\psi,\qquad
P_\psi=\psi\psi^T+(J\psi)(J\psi)^T.
\]

Direct expansion, or the rank-one complex matrix ζζ†, proves

\[
P_\psi=\tfrac12(nI+b_iA_i),\qquad |b|^2=n^2,
\qquad \Delta_J(P_\psi)=0.
\]

Every nonzero P_ψ is positive semidefinite and has real rank two (complex
rank one). Sums are positive and lie in the future cone; a positive-definite
X requires at least two summands and two suffice by complex Cholesky
factorization. This is a minimal decomposition of response data, not a
claim that a chosen decomposition solves the subsequent field equations.

In particular let H=[[a,b],[b,c]]=D²U>0 be the normalized thermo Hessian,
Δ_H=ac−b². An admitted relative phase φ gives

\[
\widehat H_\phi=\tfrac{a+c}{2}I
+b\cos\phi\,A_1+b\sin\phi\,A_2+\tfrac{a-c}{2}A_3,
\qquad \Delta_J(\widehat H_\phi)=\Delta_H>0.
\]

An explicit decomposition uses the complex columns
(√a,b e^(iφ)/√a) and (0,√Δ_H/√a), with the preceding real identification.
For the exact rational example H=[[4,2],[2,10]], φ=0, the real fields
ψ₁=(2,0,1,0) and ψ₂=(0,0,3,0) give Ĥ=P_ψ₁+P_ψ₂ and Δ_J(Ĥ)=36.

There are two different phases here. The central field phase exp(θJ)
leaves P_ψ unchanged. The relative phase in Ĥ_φ is generated instead by
exp(−φJA₃/2), which rotates A₁ into cosφ A₁+sinφ A₂. A real thermo
Hessian alone supplies the three-dimensional slice x₂=0; central phase
invariance cannot supply the missing coordinate. Admitting and specifying
relative phase is an additional response/readout contract.

## CP-3. The propagation norm becomes the source current of phase symmetry

Keep the background M^i, ρ and Q from NP-2 and NP-4, and a nonzero real
coupling g. For a real one-form a=a_μ dx^μ set

\[
\nabla_\mu=\partial_\mu+g a_\mu J,\qquad
B_a=\nabla_t-M^i\nabla_i+Q.
\]

Its zero-order gauge term is C=gJ(a_tI−M^ia_i). It is skew and commutes
with J, so the exact variational compatibility condition of NP-5 holds.
Under

\[
a'_\mu=a_\mu+\partial_\mu\chi,\qquad
\psi'=e^{-g\chi J}\psi,
\]

one has B_(a')ψ'=e^(−gχJ)B_aψ. The field action

\[
S_\psi=\frac Z2\int\rho\,\psi^TJ B_a\psi\,d^4x
\]

is therefore locally gauge invariant. Expansion of the coupling gives
S_ψ=S_ψ,0−∫j^μ a_μ d⁴x, with the **current density**

\[
\boxed{j^0=\frac{Zg\rho}{2}\psi^T\psi,\qquad
j^i=-\frac{Zg\rho}{2}\psi^TM^i\psi.}
\]

The field equation B_aψ=0 implies ∂_μj^μ=0 by NP-4's local balance.
Thus that conserved norm has an exact role as phase charge in this gauged
model. It does not become the action Hamiltonian; NP-5's Hamiltonian-sign
counterexample is unchanged. For several fields with the same coupling and
normalization, sum their currents.

There is also a causal check. For the metric q from NP-2,

\[
q(v,v)=-(v^0)^2+(v^{\rm sp}+\beta v^0)^T
(V^TV)^{-1}(v^{\rm sp}+\beta v^0).
\]

Using |b|²=n² shows q(j,j)=0 for one field. If the summed response is Ĥ,

\[
\boxed{q(j,j)=-(Zg\rho)^2\Delta_J(\widehat H).}
\]

The density weight on both sides matches; j/ρ is an ordinary vector in the
specified chart. If Zg>0 the time component is positive. At Zgρ=1, β=0,
V=I, the preceding H example gives j=(7,−2,0,3), q(j,j)=−36.
This is an exact invariant of the supplied response, not a particle mass
prediction or an identification of thermodynamic entropy with charge.

## CP-4. A closed matter–gauge variational system

Choose e²>0 and, as a common-cone constitutive hypothesis, use the same fixed
Lorentzian metric q in the gauge action. Put w=√(−det q), F=da and

\[
S[a,\psi]=S_\psi-\frac1{4e^2}\int w F_{\mu\nu}F^{\mu\nu}\,d^4x.
\]

Indices on F are raised with q; ρ and w need not be identified. Independent
compact-support variations give

\[
\boxed{B_a\psi=0,\qquad
\frac1{e^2}\partial_\nu(w F^{\nu\mu})=j^\mu,\qquad dF=0.}
\]

Antisymmetry makes the divergence of the gauge equation consistent with
CP-3's current conservation. The source is now the same native field whose
phase is transported, rather than a prescribed external current.

The quadratic gauge kinetic law, common metric and constants Z,g,e remain
model choices. This step establishes their coupled Euler–Lagrange system;
it does not select the gauge action from native axioms or prove well-posedness
of the full coupled system, a bounded matter Hamiltonian, or quantum QED.

## CP-5. An explicit λ/Hessian-to-connection chart

Use the stable curved thermo potential, in normalized affine coordinates,

\[
U(s,v)=10s-v+s^2+sv+\tfrac52v^2+\tfrac12s^2v.
\]

Then a=U_ss=2+v, b=U_sv=1+s, c=U_vv=5 and
Δ=5(2+v)−(1+s)². The affine terms set the reference values T=10 and P=1;
they leave the Hessian and connection unchanged. Use translated entropy and
volume coordinates S=S₀+s, V=V₀+v with positive S₀,V₀ and a sufficiently
small patch. On a neighborhood of (0,0), a>0 and Δ>0. The existing
PF-1 phase connection specializes to

\[
\alpha=\frac{ds-(b/a)dv}{2\sqrt\Delta},\qquad
d\alpha=-\frac5{4\Delta^{3/2}}\,ds\wedge dv.
\]

This is a direct function of the chosen potential's second and third
derivatives. Its curvature never vanishes on this chart; at the origin it is
−5/108, recovering the earlier thermo fixture.

More strongly, define

\[
Q=s,\quad P=\frac16-\frac1{2\sqrt\Delta},\quad
\chi=\arcsin\!\left(\frac{1+s}{\sqrt{5(2+v)}}\right)-\frac s6.
\]

Direct differentiation gives the exact Darboux form

\[
\boxed{\alpha=P\,dQ+d\chi,\qquad d\alpha=dP\wedge dQ.}
\]

Here P_v=5/(4Δ^(3/2))>0. Its explicit local inverse is

\[
s=Q,\qquad
v=\frac15\left[(1+Q)^2+\frac1{4(P-1/6)^2}\right]-2,
\qquad P<1/6.
\]

This supplies a concrete thermo channel in place of an abstract Darboux
existence assertion. U is an explicitly chosen stable prototype; no material
equation of state is predicted. It is distinct from the historical enthalpy/S
potential discussed in the recovered bridge.

## CP-6. Exactly when thermodynamic variations recover every gauge equation

Take N independent curved thermo channels, each with fixed profile U_A and
nonzero phase curvature. Their 2N-dimensional target T has one-form
α=Σα_A and nondegenerate σ=dα. Admit smooth state maps
Φ:M⁴→T and define a=Φ*α. This identifies the pulled connection with that
in CP-3. No independent variation of a is now assumed.

For a variation δΦ, the chain rule, equivalently Cartan's formula, gives

\[
\delta a_\mu=
\underbrace{\sigma_{AB}(\Phi)\partial_\mu\Phi^B}_{C_{A\mu}}
\delta\Phi^A+\partial_\mu(\alpha_A\delta\Phi^A).
\]

Let E^μ=e^(−2)∂_ν(wF^(νμ))−j^μ denote the full gauge residual. The state-map
Euler–Lagrange equations are exactly

\[
C_{A\mu}E^\mu-\alpha_A\partial_\mu E^\mu=0.
\]

On solutions of the matter equation, ∂_μj^μ=0 and hence ∂_μE^μ=0.
The map equations therefore reduce to CE=0. Since σ is invertible,

\[
\boxed{\operatorname{rank}C=\operatorname{rank}d\Phi;
\quad \operatorname{rank}d\Phi=4\ \Longrightarrow\ E=0.}
\]

Conversely E=0 and the matter equation imply the state-map equations at
every rank. Thus on an open rank-four patch the pulled-back action and the
independent-potential action have the same field equations after a is
identified with Φ*α. At smaller rank only that many independent combinations
of E are tested; the omitted directions must not be declared solved.

This is a **local equation-equivalence theorem**, not a bijection of gauge
classes, a global parametrization, or a proof of hyperbolicity of the redundant
scalar equations. It requires unrestricted small state-map variations inside
the stable target patch. If q, M, ρ or the target potential profiles also vary
with Φ, their additional variations must be included and the displayed
equivalence does not automatically hold.

## CP-7. Four thermo channels are necessary and sufficient at F=0

At a rank-four point, let W=im dΦ⊂T_ΦT. Then F is the restriction of σ to
W. Since dim W^σ=2N−4,

\[
\operatorname{rank}F
=4-\dim(W\cap W^\sigma)\ge 8-2N.
\]

In particular F=0 makes W an isotropic four-plane, so N≥4. This yields:

| Curved thermo channels N | Rank-four variations possible? | Smallest rank of F at such a point |
|---|---|---|
| 1 | No | Not applicable |
| 2 | Yes | 4 |
| 3 | Yes | 2 |
| 4 | Yes | 0 |

The bound is attained in canonical coordinates (Q_A,P_A), σ=ΣdP_A∧dQ_A.
For N=4, any local one-form a=a_μ(x)dx^μ has the representation

\[
Q_\mu=x^\mu,\quad P_\mu=a_\mu(x),\qquad
\sum_\mu P_\mu dQ_\mu=a.
\]

The derivative dΦ contains the identity in its Q rows, so it has rank four
even when F vanishes. Variations of P at fixed Q give arbitrary δa directly.
CP-5's inverse converts each pair to actual (s_A,v_A) values in its thermo
patch. The resulting connection differs by Σdχ_A, an exact gauge term
absorbed by the field transformation in CP-3. An ordinary local gauge can
set a_μ to zero at a chosen point; shrinking the coordinate patch keeps P
and Q inside the prototype chart. Thus the construction covers arbitrary
local smooth potentials, including the vacuum, with regular variations.

The number four here is a minimal channel count for this **symplectic
pullback representation on a supplied four-dimensional base**. It is not a
derivation of spacetime dimension or the number of fundamental interactions.

## CP-8. An explicit spurious solution excluded by the rank test

Take flat backgrounds M^i=A_i, ρ=1, Z=g=e=1, and the constant real field
ψ=(1,0,0,1). At a=0 it solves B_aψ=0 and has j=(1,0,1,0).

With just two canonical channels choose

\[
Q_1=x^1,\quad Q_2=x^3,\quad P_1=P_2=0.
\]

Then a=F=0, dΦ has rank two, and CE=0 for E=−j≠0. All matter and
state-map equations hold, while the full Maxwell equation fails. This is
an actual local counterexample in the coupled model, not merely a count of
unknowns. In the raw thermo gauge apply the corresponding phase change to ψ.

The four-channel realization Q_μ=x^μ, P_μ=0 detects the failure: its P
equations are precisely E^μ=0 and reject this nonzero source at zero field.
The verifier checks this example separately from its general chain-rule tests.

## What has advanced and what remains

The native response algebra now supplies an explicit Lorentz cone and a
causal source-current identity. The previously conserved propagation norm is
identified as the phase charge of a concrete gauge coupling. Most directly,
the full sourced gauge equations can be obtained from variations of explicit
thermodynamic state maps, with a sharp regularity/channel criterion and a
failure witness when it is omitted.

The metric/coframe and density are still backgrounds. Their dynamics,
Einstein–Hilbert selection, physical normalization, global topology and the
stability/quantization of the coupled matter model remain open. No claim
derives G, c, Λ, α, ℏ or a particle mass. The prior massive-multiplier
obstruction and the distinction between positive norm and Hamiltonian remain.

## Established mathematics and source lineage

The Hermitian Lorentz cone, rank-one spinor bilinears, minimal phase coupling,
Noether current, and Darboux/Cartan identities are established mechanisms.
The contribution here is their explicit assembly with the pinned native
representation and thermodynamic phase connection, the variation-rank
criterion, and the constructive thermo chart and counterexample. The labels
CP-1–CP-8 are not independent claims of historical priority.

* S. Başkal, E. Georgieva and Y. S. Kim, [Wigner's new physics frontier: Physics of two-by-two matrices, including the Lorentz group and optical instruments](https://arxiv.org/abs/math-ph/0310068): established two-by-two matrix Lorentz representation.
* S. E. Gralla and T. Jacobson, [Spacetime approach to force-free magnetospheres](https://arxiv.org/abs/1401.6159v4), §3.3: varying restricted Euler potentials in a Maxwell action can produce projected/force-free equations. CP-6 checks the precise analogous risk instead of assuming equation equivalence.
* A. Stern, Y. Tong, M. Desbrun and J. E. Marsden, [Variational Integrators for Maxwell's Equations with Sources](https://arxiv.org/abs/0803.2070): existing variational source and conservation methods; no new general Maxwell variational principle is asserted here.
* [PF-1–PF-3](../thermo-phase-field/THEOREM.md) supply the thermo connection and the earlier two-channel representation on nondegenerate-curvature patches. CP-6–CP-8 strengthen the variation analysis and extend regular coverage to F=0.
