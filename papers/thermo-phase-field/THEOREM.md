# Thermodynamic phase transport to an Abelian field model

Monty Dabas · PF-1–PF-7 · 29 September 2026

**Result.** The oriented response metric of TC-1–TC-6 has an explicit scalar phase connection. A single two-variable equilibrium pullback supplies only a decomposable field strength; two independent thermo phase channels remove that local rank restriction. A declared quadratic phase dynamics on a cell complex then gives a complete discrete Maxwell system, charge continuity, positive conserved field energy and two transverse propagation modes. A separate conditional dispersion theorem identifies rest energy with inertial mass times the squared propagation speed.

The contribution here is the explicit interface between the existing thermo construction, native phase arithmetic, cut compatibility and a reproducible field model. The underlying connection, Darboux and discrete-electrodynamics identities are established mathematics. This is not a priority claim for those identities. Neither the physical vacuum, a universal clock/ruler, charge normalization nor numerical fundamental constants are selected by this construction.

## Inputs, lineage and notation

Use the stable normalized thermodynamic chart of [TC-1–TC-6](../thermo-compass-foundations/THEOREM.md):

\[
H=\begin{pmatrix}a&b\\b&c\end{pmatrix}>0,
\quad \Delta=ac-b^2,
\quad R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},
\quad J_H=RH/\sqrt\Delta.
\]

Here H is the Hessian of the chosen fundamental relation, not enthalpy. Assume the fundamental relation is smooth for the local Darboux statement below; the explicit connection/curvature calculations need only the TC differentiability hypotheses. The connection on this tangent plane is Gamma_i = H^{-1} partial_i H / 2. TC-4 proves metric compatibility and parallelness of J_H. Its native scalar representation is inherited, not reconstructed from a pre-assumed electromagnetic field.

The prior [response tensor and cut ledger](../uncut-cut-measurement/NATIVE_RESPONSE_TENSOR_AND_CUT_LEDGER.md) supplies the phase-descent condition for a readout. The corrected thermo manuscript already has a [quotient-bundle interface](../thermodynamic-response-corrections/source/sections/06a0_thermodynamic_recognition_square.tex), a [conditional propagation-symbol metric](../thermodynamic-response-corrections/source/sections/06a00_response_propagation_metric.tex) and [action-lifted gauge equations](../thermodynamic-response-corrections/source/sections/06b_classical_field_equations.tex). PF-1–PF-6 give a concrete phase sector and finite dynamical realization of that programme. [SOURCE_PINS.json](SOURCE_PINS.json) binds the versions consumed.

## PF-1. Scalar phase connection of the thermo response plane

The orientation-preserving H-isometries of a two-dimensional fibre form SO(H), isomorphic to U(1). Every such isometry is

\[
G_\theta=\cos\theta\,I+\sin\theta\,J_H.
\]

**Proof.** Use an oriented orthonormal coframe B with B^T B=H. It conjugates J_H to R and every orientation-preserving H-isometry to an ordinary plane rotation. TC-2's representation x+iota y maps exactly to xI+yJ_H. Thus the group statement consumes the prior native phase representation. It does not identify this group with electric charge.

An explicit choice is

\[
B=\begin{pmatrix}\sqrt a&b/\sqrt a\\0&\sqrt\Delta/\sqrt a\end{pmatrix}.
\]

For orthonormal response components z=Bw, the connection becomes

\[
\Omega=B\Gamma B^{-1}-dB\,B^{-1}=R\,\mathcal A,
\qquad
\boxed{\mathcal A=\frac{db-(b/a)\,da}{2\sqrt\Delta}}.
\]

**Proof.** Metric compatibility makes Omega skew. Its (2,1) component is sqrt(Delta) Gamma_21/a, since dB B^{-1} is upper triangular. Substitute Gamma=H^{-1}dH/2. The resulting component is (a db-b da)/(2a sqrt(Delta)), giving the formula. Here da and db differentiate the entries of H, so the phase connection uses third derivatives of the thermodynamic potential.

For the component convention z'=exp(-chi R)z,

\[
\mathcal A'=\mathcal A+d\chi,
\qquad f=d\mathcal A,
\qquad f'=f.
\]

The scalar curvature is related to TC-5 by

\[
\boxed{f=-K\sqrt\Delta\,ds\wedge dv}.
\]

Indeed, Omega wedge Omega=0, so the transformed matrix curvature is R dA; TC-5 gives B F_Gamma B^{-1}=-K sqrt(Delta) R ds wedge dv. For a positively oriented contractible loop, parallel transport in this convention is exp(-R integral A)=exp(-R integral f). On noncontractible loops the flat holonomy data must also be retained.

At the curved TC fixture H=[[2,1],[1,5]], H_s=[[0,1],[1,0]], H_v=[[1,0],[0,0]],

\[
\mathcal A_s=1/6,\quad \mathcal A_v=-1/12,
\quad \boxed{f_{sv}=-5/108}.
\]

The derivative of A, not merely its value at the reference point, gives f. The verifier differentiates the scalar formula with nonzero fourth-order jets and checks agreement with the independently differentiated full Christoffel curvature.

## PF-2. Which cuts retain continuous phase?

For an orthogonal readout projector P, all phases descend through the readout precisely when [P,J_H]=0. This is the prior U31 criterion, now applied to the thermo plane. In real dimension two an invariant projector has rank zero or two: a nonzero real vector and its J_H image are linearly independent. A single real response channel therefore needs its phase partner for a faithful complex-line readout.

A different object is a reflection cut C satisfying C^2=I, C^T H C=H and CJ_H=-J_H C. Then

\[
C G_\theta C=G_{-\theta}.
\]

If the cut matrix itself must remain fixed, its stabilizer inside this phase group is only {I,-I}. A covariant change of readout, C'=G C G^{-1}, retains the full phase freedom of the datum. These are different operational contracts; the second must be stated rather than silently treating a fixed real cut as U(1)-invariant.

There is also a transport consequence. If C is a globally defined parallel reflection and the same connection has F=fJ_H, then

\[
0=\nabla^2 C=[F,C]=f[J_H,C]=2fJ_HC
\quad\Longrightarrow\quad f=0.
\]

Nonzero phase curvature therefore excludes imposing a parallel reflection on that same line bundle. It is compatible with covariantly described, nonparallel measurement cuts. This constrains the intended phase sector without claiming that every native cut has to carry electromagnetism.

## PF-3. Two independent thermo channels give local field-strength completeness

Let Theta:M -> Sigma pull back one two-dimensional thermo chart to an event manifold. No dimension or Lorentzian metric on M is inferred from Sigma. Its induced phase curvature is

\[
F_1=\Theta^*f
=q(\Theta)\,d(s\circ\Theta)\wedge d(v\circ\Theta).
\]

Hence rank(F_1)<=2 and F_1 wedge F_1=0. In a four-dimensional electromagnetic identification this restricts the electric and magnetic fields to E dot B=0. This is a restriction on a **single equilibrium pullback**, not on every U(1) connection or the full response-jet bundle.

Use two independent maps Theta_1, Theta_2 and the tensor product of their complex phase lines. Tensor-product parallel transport multiplies phases, so

\[
\mathcal A=\Theta_1^*\mathcal A_1+\Theta_2^*\mathcal A_2,
\qquad F=F_1+F_2.
\]

Both curvatures are closed. Their sum can have rank four: in local coordinates,

\[
F_1=dx^0\wedge dx^1,\qquad F_2=dx^2\wedge dx^3,
\qquad F\wedge F=2\,dx^0\wedge dx^1\wedge dx^2\wedge dx^3.
\]

The two scalar connections give one tensor-product U(1) connection; this construction does not assert two photon species.

**Local completeness theorem.** Every smooth closed nondegenerate two-form on a four-dimensional patch can, after restricting the patch, be represented as the sum of pullbacks of two copies of any thermo phase-curvature patch with f nonzero. Thus two independent channels are necessary for rank four and sufficient locally in that nondegenerate class.

**Proof.** On the thermo patch write f=q(s,v)ds wedge dv with q nonzero. Set p=s and r=integral from v0 to v of q(s,t)dt. Then dp wedge dr=f and (p,r) are local coordinates. Darboux's theorem for the target closed nondegenerate form gives F=dP1 wedge dQ1+dP2 wedge dQ2. Define Theta_i by composing (Pi,Qi) with the inverse thermo coordinate chart, after coordinate translations and restriction. The pullback identities give the result. The curved TC example supplies a nonzero patch because f_sv=-5/108 at its reference state. [Darboux reference](https://ocw.mit.edu/courses/18-969-topics-in-geometry-dirac-geometry-fall-2006/0977bd35075df6dbd9ab50071f0c86c3_lecture1.pdf).

This is an existence/representation theorem. It does not select the maps from native evolution. Rank-changing loci, global bundles, charge lattices and flux quantization need additional analysis. On a contractible patch two potentials with the same curvature differ by an exact form, so the representation also matches local connections up to gauge.

## PF-4. A declared quadratic phase law yields discrete Maxwell dynamics

For dynamics choose a finite oriented three-dimensional cell complex and a process parameter tau. The complex and its dimension are model data. Let d0, d1, d2 be vertex-to-edge, edge-to-face and face-to-cell coboundaries, with d1 d0=d2 d1=0. Compact unit native phases on oriented links multiply along paths; a reversed link has its dagger phase. Vertex phase changes cancel around each face, and the oriented product of face holonomies around a cell equals one.

In a local unwrapped phase chart let a be edge connection coordinates and phi a vertex potential. Define

\[
e=-\dot a-d_0\phi,\qquad b=d_1 a.
\]

The gauge change a'=a+d0 chi, phi'=phi-dot chi leaves both invariant. Compact face phases are exp(-iota b_f). The unwrapped quadratic theory is the local phase response model; periodic cosine loop costs agree with it only to quadratic order near a flat reference. Nontrivial global holonomies are not removed by declaring a chart.

Choose constant positive response matrices epsilon and nu. One explicit constitutive construction is epsilon=Y_e^T Y_e and nu=Y_b^T Y_b for full-column-rank calibrated response catalogues. Their channels and units must be declared. For any orthogonal recognition projector P,

\[
Y^T Y=Y^T P Y+Y^T(I-P)Y,
\]

with both terms positive semidefinite. This is how the earlier real information ledger enters the quadratic field response. A wholly discarded channel can destroy strict positivity; the full catalogue, not its incomplete cut, defines the model here.

Select the reversible local quadratic process action

\[
L=\tfrac12 e^T\epsilon e-\tfrac12 b^T\nu b+j^Ta-\rho^T\phi.
\]

This action, time parameter, source calibration and response weights are explicit additional dynamical data. The existence theorem PF-3 does not force their selection. Independent link variations here enlarge the dynamical model beyond one fixed equilibrium pullback. No pre-existing thermodynamic connection is silently varied as an unconstrained physical field.

Put d=epsilon e, h=nu b. The equations are

\[
\boxed{\dot b=-d_1e,\quad d_2b=0,\quad
\dot d=d_1^Th-j,\quad -d_0^Td=\rho.}
\]

**Proof.** The first two equations follow from the definitions and boundary-of-boundary identities. Variation in a gives d/dtau(-epsilon e)-(-d1^T nu b+j)=0. Variation in phi gives -d0^T epsilon e-rho=0. These give the second pair. We define discrete divergence as -d0^T, fixing the Gauss-law sign.

Differentiating Gauss and using d1 d0=0 gives

\[
\dot\rho=d_0^Tj,
\quad\text{equivalently}\quad \dot\rho+\operatorname{Div}j=0.
\]

The same condition makes the source action gauge invariant up to an endpoint term: its change integrates to integral chi^T(d0^T j-dot rho). Thus charge conservation and the propagated Gauss constraint are compatible consequences of the same finite model.

These are standard discrete Maxwell equations, derived here for the declared thermo-response phase law. Geometric variational electrodynamics predates this package; see Stern, Tong, Desbrun and Marsden [0707.4470v3](https://arxiv.org/abs/0707.4470v3) and [0803.2070v1](https://arxiv.org/abs/0803.2070v1). The result does not prove that this cell complex or these weights are the physical vacuum.

## PF-5. Positive field energy and exact finite transport

For constant epsilon,nu on a closed or periodic complex define

\[
\mathcal E=\tfrac12d^T\epsilon^{-1}d+\tfrac12b^T\nu b.
\]

Then

\[
\boxed{\dot{\mathcal E}=-e^Tj}.
\]

**Proof.** Differentiate, substitute PF-4, and cancel e^T d1^T h against h^T d1 e. Energy is nonnegative and is zero only when d=b=0. For source-free motion it is constant. Boundary ports or evolving constitutive weights require their corresponding energy terms.

With y=(d,b), source-free evolution is dot y=K y with

\[
K=\begin{pmatrix}0&d_1^T\nu\\-d_1\epsilon^{-1}&0\end{pmatrix},
\qquad W=\operatorname{diag}(\epsilon^{-1},\nu),
\quad K^TW+WK=0.
\]

The implicit-midpoint map U_h=(I-hK/2)^{-1}(I+hK/2) preserves W exactly for every real h. To see invertibility, a hypothetical (I-hK/2)v=0 gives v^T Wv=0 after multiplication by v^T W, so v=0. Expanding (I+hK/2)^T W(I+hK/2) and its minus counterpart proves their equality, hence U_h^T W U_h=W. Left null rows giving Gauss and magnetic-divergence constraints are also preserved. This is a time-discretization control, not the exact continuous exponential.

The certificate checks nonconstant positive response weights on a periodic 2-by-2-by-2 complex: 8 vertices, 24 links, 24 faces, 8 cells. It verifies the energy and both constraints through exact rational midpoint steps. It also checks source work and charge continuity independently at sampled states.

## PF-6. Two transverse modes and the propagation scale

For a homogeneous cubic lattice set epsilon=epsilon0 I, nu=nu0 I with positive scalars. In temporal gauge phi=0, source-free propagation obeys

\[
\ddot a=-\frac{\nu_0}{\epsilon_0}d_1^Td_1a.
\]

For a Fourier character with unit phases z_i=exp(i q_i), set delta_i=z_i-1. The curl symbol is C(delta)v=delta cross v, and direct multiplication gives

\[
C^\dagger C=\Lambda I-\delta\delta^\dagger,
\qquad \Lambda=\sum_i|\delta_i|^2=4\sum_i\sin^2(q_i/2).
\]

For delta nonzero the longitudinal direction delta has eigenvalue zero; its Hermitian orthogonal complement has dimension two and eigenvalue Lambda. Gauss selects that transverse space. Therefore the model has **two propagating transverse modes** per nonzero wavevector, with

\[
\boxed{\omega_\tau^2=4\frac{\nu_0}{\epsilon_0}\sum_i\sin^2(q_i/2)}.
\]

The constant mode delta=0 is a separate zero-frequency sector. On the finite side-two periodic complex, the exact ranks are rank d0=7, rank d1=14, rank d2=7, agreeing with two curl modes for each of seven nonconstant Fourier characters. Harmonic/global modes are retained.

With physical lattice spacing ell and time calibration t=t_* tau, q_i=ell k_i, the long-wavelength expansion is

\[
\omega^2=c_{\rm eff}^2|k|^2+O\!\left(\frac{\nu_0\ell^4}{\epsilon_0t_*^2}\sum_i k_i^4\right),
\qquad
\boxed{c_{\rm eff}^2=\frac{\nu_0}{\epsilon_0}\frac{\ell^2}{t_*^2}}.
\]

This computes a propagation scale from the selected response weights and calibrations. A finite lattice has dispersive corrections; exact Lorentz symmetry or universal physical c is not inferred. Scaling nu while holding the phase geometry fixed changes this speed, demonstrating why the dynamical input matters.

## PF-7. Conditional rest-energy/inertial-mass bridge

This corollary concerns a separate gapped matter mode, not a mass term inserted into the gauge potential. Suppose its continuum dispersion and energy/momentum calibration are

\[
\omega(k)^2=\Omega_0^2+c^2|k|^2,
\qquad E=\eta\omega,\quad p=\eta k,
\qquad \Omega_0,c,\eta>0.
\]

A scalar quadratic process law with a positive local restoring term can realize such a dispersion in its continuum limit. Selecting that mode, its gap and the universal status of eta is further work. Eta is an action scale here and has not been identified with or used to derive Planck's constant.

Define inertial mass operationally through the low-momentum energy curvature,

\[
\left.\frac{\partial^2E}{\partial p_i\partial p_j}\right|_{p=0}
=\frac{\delta_{ij}}{m_{\rm in}}.
\]

Substitution gives E=sqrt(eta^2 Omega0^2+c^2|p|^2), so

\[
\boxed{m_{\rm in}=\frac{\eta\Omega_0}{c^2},\qquad
E_0=m_{\rm in}c^2,\qquad E^2-c^2|p|^2=m_{\rm in}^2c^4.}
\]

Equivalently E=E0+|p|^2/(2m_in)+O(|p|^4). Thus the mass is read from the quadratic response to momentum, and the rest-energy relation follows under the dispersion/calibration hypotheses. This is not a numerical mass prediction or a derivation of those hypotheses. For the gapless gauge mode E=c|p| in the continuum regime; nonzero travelling energy does not require nonzero rest mass.

## What the next theorem must select

The open physical step is to select the two-channel maps or broader jet sector, the event complex and process law from native response rules, then identify charge, clocks, rulers and matter gaps. The present result provides explicit phase geometry, a local field-strength representation, a field equation, conservation laws and propagation controls to test such a selector against. It does not substitute the old name “information” for measured energy: the response ledger fixes a positive quadratic structure, while energy units and dynamics enter through declared interfaces.
