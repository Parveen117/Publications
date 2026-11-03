# Thermo compass to response geometry and transport

Monty Dabas research programme · TC-1–TC-6 · 29 September 2026

**Result.** In a stable two-dimensional thermodynamic chart, the energy response Hessian determines the caloric and mechanical responses. Its positive form and an orientation give a local representation of the native quarter-turn scalar. Declaring that Hessian as a response metric gives a unique torsion-free metric connection. In the natural affine chart its curvature depends on third derivatives of the fundamental relation through an exact commutator; fourth derivatives cancel. This constructs a response-plane transport geometry before attempting constants.

This is a conditional mathematical synthesis of established thermodynamic and Hessian-geometric identities, not a claim of new general geometric identities or a completed microscopic physics theory. The fundamental relation, stable regime, orientation, unit interface and use of its Hessian as a metric are the inputs. The resulting geometry belongs to thermodynamic state space. It is not automatically EMK recognition curvature, spacetime curvature, a physical clock, or a quantum Hamiltonian.

## 1. Source order and hypotheses

The source is the corrected [thermodynamic manuscript](../thermodynamic-response-corrections/README.md), especially its foundations, unit-typed response section, intrinsic bracket, and Cut–Flow–Jet capstone. The exact source edition was preserved from Publications commit `c1a35e0f509cd137d20f9504543f97b5ec7479fc`; [SOURCE_PINS.json](SOURCE_PINS.json) records its hashes and the RKF foundation dependencies.

Use a simple compressible equilibrium chart with fixed amount of substance, positive temperature and volume, and a smooth specific or molar fundamental relation U(S,V). Assume U is C^4 and its Hessian is positive definite on the chart. A fully extensive U(S,V,N) before fixing N is not the positive-definite two-variable object used here.

T,V,S,P are overlapping thermodynamic chart variables, not four independent spacetime coordinates. The corrected source has six constrained response components in a unit-typed direct sum. Its covariant tower may be used after its connection and transport hypotheses are met; an ordinary four-variable derivative matrix is not assumed skew.

TC-1 below keeps physical units. For TC-2 onward, choose constant positive reference scales S_*, V_*, U_* and set s=S/S_*, v=V/V_*, u=U/U_*. Coordinate translations are allowed. Let

\[
H=D^2_{s,v}u=\begin{pmatrix}a&b\\b&c\end{pmatrix},
\quad \Delta=ac-b^2>0,
\quad R=\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]

Choosing the ordered (s,v) orientation fixes R; reversing orientation reverses the quarter-turn below. The letter H in these formulas is the response Hessian, not thermodynamic enthalpy.

## TC-1. Reconstruct the four basic responses from one Hessian

In dimensional (S,V) variables temporarily write a=U_SS, b=U_SV, c=U_VV and Delta=ac-b^2. The compass foundation gives

\[
T=U_S,\qquad P=-U_V,\qquad dU=T\,dS-P\,dV,
\]

\[
dT=a\,dS+b\,dV,\qquad dP=-b\,dS-c\,dV.
\]

At fixed V, (dT/dS)_V=a. At fixed P, dV=-(b/c)dS, so (dT/dS)_P=Delta/c. At fixed T, dS=-(b/a)dV, so (dP/dV)_T=-Delta/a. Thus

\[
\boxed{C_V=T/a,\quad C_P=Tc/\Delta,\quad K_S=Vc,\quad K_T=V\Delta/a.}
\]

It follows that

\[
\frac{C_P}{C_V}=\frac{K_S}{K_T}=\frac{ac}{\Delta},
\quad \Gamma_c\Gamma_m=1,
\quad 1-\frac{C_V}{C_P}=\frac{b^2}{ac}\in[0,1).
\]

The last ratio measures the normalized cross-response in the distinguished thermodynamic S,V directions. It is invariant under separate constant rescalings of these axes, not under arbitrary mixing of what the constraints mean. It is neither an electromagnetic coupling nor a return-sector sign.

The six source responses follow without new constants:

\[
\lambda_P=-S\Delta/c,\quad\lambda_V=-Sa,\quad
\lambda_S=-Vc,\quad z_T=-V\Delta/a,\quad
z_S=P/c,\quad\lambda_T=Pa/\Delta.
\]

These are exact constrained-derivative identities inside the supplied stable equation of state.

## TC-2. Represent the native quarter-turn on the response plane

Return to the normalized Hessian H and define

\[
\boxed{J_H=\frac{RH}{\sqrt\Delta}}.
\]

Direct multiplication gives RH RH=-Delta I and H R H=Delta R. Hence

\[
J_H^2=-I,\qquad J_H^T H=-HJ_H,\qquad J_H^THJ_H=H.
\]

For native scalar z=x+iota_Sigma y, define rho_H(z)=xI+yJ_H. The quarter-turn identity proves multiplicativity, and the H-adjoint relation proves

\[
\rho_H(z)^{{\dagger}_H}=\rho_H(z^\dagger),
\quad \rho_H(z)^{{\dagger}_H}\rho_H(z)=(x^2+y^2)I.
\]

This is faithful: taking the trace of xI+yJ_H=0 gives x=0, and J_H is invertible, so y=0. The native scalar is inherited from F00/F00-E. This theorem constructs its representation on a specific thermodynamic tangent plane; it does not derive the primitive scalar a second time or identify this plane with an electron state.

For a constant orientation-preserving change x'=Ax, transport the metric as H'=A^{-T}HA^{-1}. The identity R A^{-T}=(A R)/det A gives J_H'=A J_H A^{-1}. Constant positive rescaling H->lambda H also leaves J_H unchanged. Under orientation reversal one must transport the orientation too; silently holding R fixed changes the sign.

The raw coordinate Hessian of a scalar is not a tensor under arbitrary nonlinear coordinate changes. After defining g=H_ij dx^i dx^j in the affine thermodynamic chart, transport g as a metric tensor. Recomputing raw second derivatives in a nonlinear or Legendre chart is a different operation.

## TC-3. First response jet of the quarter-turn

For a differentiable path of positive Hessians, Jacobi's determinant derivative gives

\[
\boxed{\partial_iJ_H=
\frac{R\,\partial_iH}{\sqrt\Delta}
-\frac12\operatorname{tr}(H^{-1}\partial_iH)J_H.}
\]

This also follows by differentiating Delta=ac-b^2. Differentiating J_H^2=-I yields

\[
J_H(\partial_iJ_H)+(\partial_iJ_H)J_H=0.
\]

Thus the derivative of the represented phase structure is constrained by the third-order energy response. It is not an independent parameter to be fitted. Higher derivatives follow by differentiating this identity on the same regular chart. At Delta=0 this normalization fails; such a point is outside the stable chart, not an invitation to divide by zero.

## TC-4. The Hessian response metric supplies a transport connection

Declare the mathematical response metric

\[
g=H_{ij}\,dx^i dx^j.
\]

This is the positive second-variation geometry of the supplied fundamental relation. Impose metric compatibility and zero torsion. These select its Levi-Civita connection; they are the specified selection conditions here, not a claim that the compass uniquely selects every physical transport law.

The general connection formula is

\[
\Gamma^k{}_{ij}=\frac12H^{k\ell}
(\partial_iH_{\ell j}+\partial_jH_{\ell i}-\partial_\ell H_{ij}).
\]

Since H is a Hessian, its third derivative tensor is completely symmetric. Therefore, with (A_i)^k{}_j=Gamma^k_ij,

\[
\boxed{A_i=\frac12H^{-1}\partial_iH.}
\]

The coordinate expression holds in the affine Hessian chart. Under a general coordinate change, the connection transforms with its usual derivative term.

Metric compatibility follows immediately:

\[
\partial_iH=A_i^TH+HA_i.
\]

Moreover

\[
\boxed{\partial_iJ_H+[A_i,J_H]=0.}
\]

To check the last identity, put X=H^{-1}partial_iH. The two-dimensional identity X+J_H X J_H^{-1}=tr(X)I for H-self-adjoint X, combined with TC-3, gives cancellation. Equivalently, the metric connection preserves the oriented metric area form, and therefore its compatible quarter-turn. This is parallel transport of the representation already constructed in TC-2.

## TC-5. Curvature from the third-order response tower

Use the convention F_ij=[nabla_i,nabla_j] on tangent vectors, so

\[
F_{ij}=\partial_iA_j-\partial_jA_i+[A_i,A_j].
\]

Let X_i=H^{-1}partial_iH. Differentiating the inverse gives

\[
\partial_iX_j=-X_iX_j+H^{-1}\partial_i\partial_jH.
\]

The symmetric second derivatives of H cancel in the antisymmetrization. Hence

\[
\boxed{F_{ij}=-\frac14[X_i,X_j]
=-\frac14[H^{-1}\partial_iH,H^{-1}\partial_jH].}
\]

This is an exact formula for the curvature of the declared response metric. It depends on U'' and U''' at the point; U'''' is needed for the differentiated expression but cancels. The derivative matrices X_i need not be skew. Curvature here arises from their order of composition, not from incorrectly declaring every Jacobian antisymmetric.

In two dimensions define the Gaussian curvature by

\[
K=\frac{(HF_{sv})_{12}}{\det H}.
\]

With the above orientation, F_sv=-K sqrt(Delta) J_H. This is thermodynamic metric curvature. Its identification with an EMK metric or recognition connection would require the corresponding metric map or intertwiner; none is asserted here.

### Curved exact example

In shifted, dimensionless affine coordinates, take

\[
u(s,v)=10s-10v+s^2+sv+\frac52v^2+\frac12s^2v.
\]

Let physical entropy and volume origins be positive, for example S/S_*=1+s and V/V_*=1+v. At the origin temperature is positive and H is positive definite, so an admissible neighbourhood exists. At that point

\[
H=\begin{pmatrix}2&1\\1&5\end{pmatrix},\quad
H_s=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
H_v=\begin{pmatrix}1&0\\0&0\end{pmatrix}.
\]

Then

\[
J_H=\frac13\begin{pmatrix}-1&-5\\2&1\end{pmatrix},\quad
F_{sv}=\begin{pmatrix}5/324&25/324\\-5/162&-5/324\end{pmatrix},
\quad \boxed{K=5/324}.
\]

The verifier obtains this curvature both from the commutator formula and from differentiating the full Christoffel formula. A quadratic u with the same point Hessian has zero third derivatives and zero curvature. This comparison explains why retaining the higher response tower matters; a point susceptibility matrix alone is insufficient to calculate neighbourhood transport geometry.

## TC-6. Lawful path transport and its exact control

Given a C^1 path x(t) in the stable chart and initial tangent response w_0, solve

\[
\dot w+A_i(x(t))\dot x^i w=0.
\]

Continuous coefficients give a unique solution along every compact path segment within the chart. TC-4 implies

\[
\frac{d}{dt}(w^THw)=0,
\qquad U_\gamma J_H(x_0)=J_H(x_1)U_\gamma.
\]

For the norm identity, substitute the transport equation and partial_iH=A_i^TH+HA_i. For the intertwiner, TC-4 says that J_Hw solves the same transport equation as w, with initial value J_H(x_0)w_0, and uniqueness completes the proof.

An exact nonconstant example is the Hessian of

\[
u=10s-10v+\frac{s^2}{2}+\frac{s^3}{3}+\frac{s^4}{12}+\frac{v^2}{2},
\quad H=\operatorname{diag}((1+s)^2,1),\quad s>-1,
\]

on a local positive-temperature chart. Along v=0, from s=0,

\[
U(s)=\operatorname{diag}((1+s)^{-1},1)
\]

satisfies the transport equation exactly and preserves both the metric norm and J_H intertwining. Parameter s, or a chosen path parameter t, is not thereby a physical clock.

## 2. How this enters the Jacobian tower

The source response section Lambda is a collection of constrained scalar responses with fixed unit scales. TC-4 provides a tangent/cotangent connection. On tensors, extend it by the product rule; on the nondimensional scalar-response components use a declared flat internal connection as a first local realization. Then

\[
L_1=\widehat\Lambda,\qquad L_{n+1}=\nabla L_n
\]

has the source paper's tensor typing. Other response-fibre connections are possible and must be declared. This tangent construction does not silently replace the thermodynamic quotient connection or the EMK recognition connection.

The finite jet raiser and process transport of the source capstone can consume this connection once their holonomic and cut-compatibility hypotheses hold. In particular, do not assume [T_X,R_nabla]=0 merely because both were built from derivatives. Where it does not vanish, retain the commutator defect instead of applying the commuting exponential factorization.

The current result supplies equilibrium response, its local phase representation, a selected metric-compatible path transport, and an explicit higher-jet curvature. A physical process law selecting paths, kinetic coefficients, a clock, electromagnetic matter, or spacetime geometry remains a subsequent theorem/adapter problem. No value of alpha, hbar, c or G is claimed.

## 3. Evidence and lineage

The proofs above establish the stated identities for the specified class. [verify.py](verify.py) checks exact rational implementations, the direct Christoffel route, coordinate covariance, nonconstant transport, inadmissible inputs and negative controls. Its generated [CERTIFICATE.json](CERTIFICATE.json) binds the local source bytes; default verification never rewrites approved evidence.

The recovered helical result has its own earlier proof and controls and is invoked separately by the new verification entry point. The earlier Physics-from-Cuts master certificate remains unchanged. The downstream alpha note is retained as a later-stage conditional reduction, not the starting axiom of this programme.

These are written proofs and reproducible finite controls. Proof-assistant formalization, independent peer review and new empirical validation are not claimed. The theorem selection is new to this programme's current development step; priority for the underlying standard identities is not claimed.
