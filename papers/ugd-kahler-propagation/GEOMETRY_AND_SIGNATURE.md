# λ geometry, signature and an explicit Einstein benchmark

Author: Monty Dabas. Development edition: 29 September 2026.

This note consolidates the geometric development following the
[recovered bridge](history/RKF_UGD_Kahler_Einstein_Bridge.md). All coordinates
and potentials in the calculations are nondimensionalized. Physical units and
their calibration would be additional data.

## GS-1. A Kähler metric cannot have Lorentzian signature

Let a real nondegenerate symmetric form g admit J with J²=−I and
g(Ju,Jv)=g(u,v). For any non-null v, g(v,Jv)=0 and
g(Jv,Jv)=g(v,v). Thus span{v,Jv} is a definite two-plane; its orthogonal
complement is nondegenerate and J-invariant. Induction pairs both the positive
and the negative directions. The real signature is (2p,2q).

Consequently a pseudo-Kähler four-manifold can have signatures (4,0), (2,2) or
(0,4), but not (1,3) or (3,1). The original suggestion that ordinary
pseudo-Kähler compatibility permits (1,3) is false. Choosing a time orientation,
requiring entropy increase, or imposing C_P>C_V cannot defeat this pointwise
linear-algebra obstruction. Thermodynamic stability constrains the relevant
state-space response, not automatically the spacetime metric.

Likewise, if a native involution L anticommutes with an invertible R, then R
is an isomorphism between the +1 and −1 eigenspaces of L. On a real
four-dimensional carrier the projector (I+L)/2 therefore has rank two, not one.
Flipping this whole cut in an admitted positive metric gives a two-plus-two
split when L is self-adjoint. It does not select one time direction.

The tensor g+iΩ is Hermitian when Ω is real antisymmetric. This alone does not
give J²=−I, integrability or dΩ=0. Also, the 2n-by-2n matrix G+iΩ on the
complexified real tangent of a compatible Kähler metric has a kernel. The
determinant below is that of the n-by-n complex metric h, not this doubled
matrix.

## GS-2. What a quadratic hyperbolicity criterion actually selects

Let P(ξ)=Q(ξ,ξ) be a nondegenerate real quadratic form and let P(n)>0.
Assume that every polynomial s↦P(k+sn) has only real roots. Decompose k as
an n-parallel component plus w with Q(w,n)=0. For such w,

\[
P(w+sn)=P(w)+s^2P(n).
\]

Real roots require P(w)≤0. The restriction to n-perpendicular is
nondegenerate, so it is negative definite. Hence P has signature (1,d−1).
Conversely that signature with P(n)>0 gives real roots by this decomposition.
An overall sign convention reverses the signature labels.

This proves a conditional signature criterion; it does not derive the
hyperbolic process law from entropy. Merely having nonzero null vectors is
insufficient: diag(1,1,−1,−1) has them, but for n=e₀ and w=e₁ the polynomial
is s²+1. The explicit field symbol in NP-2 supplies a Lorentzian example under
its continuum and constitutive hypotheses.

## GS-3. The Hessian lift and its Einstein equation

Fix a flat torsion-free affine connection and affine coordinates x on a
connected chart. Let h_ab=∂a∂b λ be smooth and nondegenerate. On an independent
companion copy y, set z=x+iy and

\[
\mathcal G=2h_{ab}(dx^a dx^b+dy^a dy^b),\qquad
J\partial_{x^a}=\partial_{y^a},\qquad
\Omega=2h_{ab}\,dx^a\wedge dy^b.
\]

Third-derivative symmetry gives dΩ=0; this J is integrable. Thus the lift is
pseudo-Kähler, and Kähler when h>0. This is the established Hessian/affine
r-map construction, not a new general theorem about thermodynamic spacetime.
Our normalization has Kähler potential 4λ and complex components h_ab.

With constant χ, the complex Ricci tensor is

\[
\operatorname{Ric}_{a\bar b}
=-\tfrac14\partial_a\partial_b\log|\det D^2\lambda|.
\]

Therefore, locally on the connected affine chart,

\[
\boxed{\operatorname{Ric}(\mathcal G)=\chi\mathcal G
\iff \log|\det D^2\lambda|+4\chi\lambda=a\cdot x+b.}
\]

In real dimension 2n the associated vacuum Einstein equation has
Λ_geom=(n−1)χ. This is the Einstein condition on this restricted family of
metrics. It is not the general sourced Einstein equation, and a variation
restricted to λ need not be equivalent to arbitrary metric variation.

## GS-4. A solved boundary problem and a constructed scalar functional

Choose Ω₀=B_(1/4)((1,1))⊂R² and boundary data from

\[
\lambda_*=-\log x_1-\log x_2.
\]

Then D²λ_*=diag(x₁⁻²,x₂⁻²)>0 and

\[
\det D^2\lambda_*=\frac1{x_1^2x_2^2}=e^{2\lambda_*}.
\]

This supplies a real-analytic solution on a neighborhood of the closed domain
with χ=−1/2 and a=b=0 in GS-3. It also proves uniqueness among
C²(Ω₀)∩C⁰(cl Ω₀) solutions with positive-definite Hessian and the same boundary
data. If u−v has a positive interior maximum, then D²u≤D²v there. Positivity
implies det D²u≤det D²v, contradicting e^(2u)>e^(2v). Interchanging u and v
rules out a negative minimum. Thus u=v.

This is existence, uniqueness and smoothness for this explicitly supplied
positive-branch boundary problem. It is not an existence theorem for arbitrary
thermodynamic data or the indefinite Monge–Ampère equation. Yau's compact
Kähler theorem does not automatically apply to that different problem.

For compactly supported variations in two affine dimensions, define

\[
\mathcal I[\lambda]=\int_{\Omega_0}
\left(\frac13\lambda\det D^2\lambda-\frac12e^{2\lambda}\right)d^2x.
\]

Let C^ij be the cofactor matrix of D²λ. Direct differentiation gives
∂i C^ij=0 and C^ij λ_ij=2 det D²λ. Two integrations by parts therefore yield

\[
\delta\mathcal I
=\int_{\Omega_0}(\det D^2\lambda-e^{2\lambda})\,\delta\lambda\,d^2x.
\]

This is an inverse variational construction for the benchmark PDE. No native
principle has selected this functional or identified it with a four-dimensional
gravitational action. Only compact-support variations are asserted; no
unwritten higher-derivative boundary conditions are imposed.

## GS-5. An explicit Lorentzian Einstein sector, with a failed alternative

After making an additional choice of one companion direction as t and
reversing its sign, define, for ℓ>0 and x₁,x₂>0,

\[
q=2\ell^2\left[
\frac{-dt^2+dx_1^2}{x_1^2}
+\frac{dx_2^2+dy^2}{x_2^2}\right].
\]

Each two-dimensional factor has Gaussian curvature −1/(2ℓ²). The product
connection preserves both factors and has no mixed curvature, so

\[
\operatorname{Ric}(q)=-\frac{q}{2\ell^2},\qquad
R(q)=-\frac2{\ell^2},\qquad
\Lambda=-\frac1{2\ell^2}.
\]

Thus G(q)+Λq=0. At ℓ=x₁=x₂=1, in order (t,x₁,x₂,y),
q=diag(−2,2,2,2), Ric=diag(1,−1,−1,−1), R=−2 and Λ=−1/2.
This is a standard AdS₂×H² local product. The chosen Lorentzian metric is no
longer Kähler with its original J. ℓ is free, so this is not a prediction of
the observed cosmological constant.

A signature flip alone need not preserve Einstein geometry. For the same
logarithmic λ and the positive lift G, take u=dλ/||dλ||_G and
q_grad=G−2u⊗u. In order (x₁,y₁,x₂,y₂), its only nonzero components are

\[
(q_{\rm grad})_{13}=(q_{\rm grad})_{31}=-\frac2{x_1x_2},\quad
(q_{\rm grad})_{22}=\frac2{x_1^2},\quad
(q_{\rm grad})_{44}=\frac2{x_2^2},
\]

where displayed component labels are one-based. At x₁=x₂=1 direct
Levi–Civita calculation gives Ric=diag(−1,1,−1,1) and scalar curvature 1.
The nonzero diagonal Ricci entries where q_grad has zero diagonal exclude
Ric=χq_grad for every χ. The exact verifier independently checks both metrics.

## Sources and relationship to established work

* P. Osipov, [Selfsimilar Hessian and conformally Kähler manifolds](https://arxiv.org/abs/2012.03791v2), §2: the established Hessian-to-Kähler tangent-bundle construction.
* S.-T. Yau, *On the Ricci curvature of a compact Kähler manifold and the complex Monge–Ampère equation, I*, CPAM 31 (1978), 339–411, DOI 10.1002/cpa.3160310304: a distinct compact Kähler existence problem; no invocation of its theorem is needed for GS-4.
* T. Jacobson, [Thermodynamics of Spacetime: The Einstein Equation of State](https://arxiv.org/abs/gr-qc/9504004v2): a distinct route requiring entropy/area proportionality, local causal horizons, Unruh temperature and a Clausius relation.
* [Corrected thermo propagation metric](../thermodynamic-response-corrections/source/sections/06a00_response_propagation_metric.tex), [classical action](../thermodynamic-response-corrections/source/sections/06b_classical_field_equations.tex), and [dynamical nonselection](../thermodynamic-response-corrections/source/sections/06a1_dynamical_nonselection.tex): existing action-conditional reductions remain conditional.

Priority for these general mathematical mechanisms is not claimed. The
framework-specific contribution here is the explicit coefficient compatibility,
the hypotheses and counterexamples connecting the separate sectors, and their
reproducible verification.
