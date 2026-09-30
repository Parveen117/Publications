# Native massive matter, covariant sources and algebraic torsion

Author: Monty Dabas. Research development R5, 30 September 2026.

This is a conditional classical extension of [NP](NATIVE_PROPAGATION_ACTION.md),
[CP](THERMO_GAUGE_COMPLETION.md), [MG](NATIVE_METRIC_DYNAMICS.md) and
[LR](NATIVE_LOOP_CURVATURE.md). It repairs NP's massive variational obstruction
by doubling the admitted real coefficient module, constructs a local Lorentz-
and phase-invariant matter action, and varies its field, gauge, coframe and
connection arguments. In LR's unprojected bulk sector, those variations yield
Einstein–Cartan–Maxwell equations with an explicit massive classical source.
Eliminating the independent connection gives Einstein equations with a
specific contact interaction in the effective matter action.

The continuum, spin structure, nondegenerate coframe, classical commuting
field, minimal coupling prescription, module and constants are hypotheses.
This does not derive a quantum fermion, measured particle mass, or the choice
of this action from a primitive native process. The mechanism of a spin source
producing algebraic torsion and a contact interaction is established
Einstein–Cartan theory, not a new discovery of gravity; see the primary
references below. The contribution here is its explicit realization and
normalization in the already source-bound native coefficient sector.

## Conventions and hypotheses

Work locally on an oriented, time-oriented Lorentzian four-manifold admitting
a spin structure, or on a contractible patch with a chosen spin frame. All
variations have compact support. Set the coordinate propagation speed to one;
this is a unit convention, not a derivation of physical c or hbar. Let

\[
\eta_{ab}=\operatorname{diag}(-1,1,1,1),\qquad
e^a=e^a{}_{\mu}dx^{\mu},\qquad w=\det(e^a{}_{\mu})>0,
\qquad g_{\mu\nu}=\eta_{ab}e^a{}_{\mu}e^b{}_{\nu}.
\]

The inverse is denoted e_a{}^mu, and vol_e=w d^4x. Lorentz indices are raised
and lowered with eta. The independent metric-compatible connection has
omega^{ab}=-omega^{ba}, with

\[
T^a=de^a+\omega^a{}_b\wedge e^b,\qquad
R^{ab}=d\omega^{ab}+\omega^a{}_c\wedge\omega^{cb},\qquad
\epsilon_{0123}=+1.
\tag{SM.1}
\]

The electromagnetic potential is a, with F=da. In particular, F is not
defined by antisymmetrizing a torsionful affine derivative. This keeps the
phase gauge symmetry intact. Parameters Z>0, kappa!=0, m and g are supplied;
g is the matter charge coupling, not the metric. The Maxwell kinetic
normalization e_em^2>0 is another supplied parameter. We restrict MG's real
Holst coefficient to **gamma_H=0**, the bulk value of LR's unprojected readout.
Nonzero Holst coefficients with spin sources require a different elimination.

NP's four numerical matrices are now denoted Gamma^a (upper Lorentz index):

\[
\Gamma^0=R\otimes I,\quad \Gamma^1=K\otimes K,\quad
\Gamma^2=K\otimes L,\quad \Gamma^3=L\otimes I,\quad
A_i=\Gamma^0\Gamma^i,\quad J=\Gamma^0\Gamma^1\Gamma^2\Gamma^3.
\]

Here R^2=-I, K^2=L^2=I, L=KR and KR=-RK, exactly as before. This is a
coefficient-module calculation, not a faithful finite-dimensional realization
of a proper one-sided seam TS=I, ST!=I.

## SM-1. Minimal doubling removes the constant-multiplier obstruction

On two copies of the four-dimensional real module define

\[
\gamma^a=I_2\otimes\Gamma^a,\quad
\mathcal I=R\otimes I_4,\quad \mathcal J=I_2\otimes J,\quad
\mathcal A_i=I_2\otimes A_i,\quad \Psi=\begin{pmatrix}u\\v\end{pmatrix}.
\]

The external native return Ical satisfies Ical^2=-I, Ical^T=-I and commutes
with every gamma^a. It is distinct from the Clifford volume Jcal: Jcal
anticommutes with the gamma^a. For the NP massive equation

\[
B_m\Psi=\left(\partial_t-\mathcal A_i\partial_i+m\gamma^0\right)\Psi=0,
\qquad S_{\rm flat}=\frac Z2\int\Psi^T\mathcal I B_m\Psi\,d^4x,
\tag{SM.2}
\]

Ical B_m is formally self-adjoint. Its variation is Z Ical B_m Psi; the mass
term survives. Indeed B_m is formally skew-adjoint and commutes with Ical.
The mass matrix Ical gamma^0 is symmetric and nonzero. In contrast, the
single-copy J Gamma^0 in NP is skew and contributes zero to a commuting
quadratic action. That earlier obstruction remains valid on its original
module.

**Minimality in the stated class.** Gamma^0 and the A_i generate all M_4(R):
Gamma^i=-Gamma^0 A_i, and the sixteen Clifford words span the full algebra.
On r identical copies their commutant is M_r(R) tensor I_4. A constant skew
multiplier commuting with the kinetic and massive matrices therefore has the
form C tensor I_4 with C^T=-C. Such a multiplier is invertible precisely when
r is even. The smallest choice is r=2, where C is a nonzero multiple of R.
The exact unrestricted 64-variable matrix calculation gives commutant rank
60 and commutant-plus-skew rank 63, independently confirming this conclusion.
This is not a minimality theorem for arbitrary modules, derivative-dependent
multipliers, auxiliary fields or Grassmann-valued variables.

The evolution generator G=Acal_i partial_i-m gamma^0 obeys
G^2=Delta-m^2. Thus every component satisfies
(partial_t^2-Delta+m^2)Psi=0. The positive norm remains conserved under
suitable boundary conditions, but the classical Hamiltonian is indefinite.
The exact rest-state control has equal-norm states with opposite energies.
No positive-energy quantum Hilbert space or spin-statistics result follows.

## SM-2. A covariant action and the off-shell torsion term

Set

\[
H=-\mathcal I\gamma^0,\quad \bar\Psi=\Psi^T H,\quad
\gamma^\mu=e_a{}^\mu\gamma^a,\quad
\Omega_\mu=\frac14\omega_{ab\mu}\gamma^a\gamma^b,
\quad D_\mu\Psi=(\partial_\mu+\Omega_\mu+g a_\mu\mathcal I)\Psi.
\tag{SM.3}
\]

Here H^T=H, H^2=I, (H gamma^a)^T=-H gamma^a, and
Omega^T H+H Omega=Ical^T H+H Ical=0. The lowered omega_ab in (SM.3)
matters. Relative to MG's convention omega^{ab} Gamma_a Gamma_b/4, this
matter representation is the constant similarity by I_2 tensor Gamma^0;
its boosts have the corresponding reversed matrix sign. It uses the same
real coframe and connection as MG.

Declare the minimal curved completion

\[
S_\Psi=\frac Z2\int w\left(\bar\Psi\gamma^\mu D_\mu\Psi
                         -m\bar\Psi\Psi\right)d^4x.
\tag{SM.4}
\]

For commuting real fields, this equals the symmetrized kinetic action
Z/4 times [bar(Psi) gamma^mu D_mu Psi-(D_mu bar(Psi)) gamma^mu Psi],
with the same mass term. The equality follows from the skew matrix H gamma^mu;
it does not require imposing the connection equation. At e^a_mu=delta^a_mu,
omega=a=0, H(gamma^a partial_a-m)=Ical B_m, so (SM.4) is exactly (SM.2).

Local spin transformations S obey S^T H S=H and
S gamma^a S^{-1}=(Lambda^{-1})^a_b gamma^b. Together with
e'=Lambda e and Omega'=S Omega S^{-1}-(dS)S^{-1}, these make (SM.4)
Lorentz invariant. Phase transformations

\[
a'=a+d\chi,\qquad \Psi'=\exp(-g\chi\mathcal I)\Psi
\tag{SM.5}
\]

also leave it invariant. The two actions commute. Exact controls include both
a rational boost and a rational rotation with nonzero connection
inhomogeneities, plus local phase transformations.

**Field variation.** Write delta S_Psi=int delta Psi^T E_Psi d^4x. Direct
integration by parts gives

\[
\frac{E_\Psi}{Zw}
=H\left[\gamma^\mu\partial_\mu
 +\frac{1}{2w}\partial_\mu(w\gamma^\mu)
 +\frac12\{\gamma^\mu,\Omega_\mu+g a_\mu\mathcal I\}-m\right]\Psi.
\tag{SM.6}
\]

Let T^rho_mu nu be the spacetime torsion corresponding to (SM.1), and define
t_mu=T^nu_mu nu. The tetrad compatibility identity implies

\[
\frac1w\partial_\mu(w\gamma^\mu)+[\Omega_\mu,\gamma^\mu]
=t_\mu\gamma^\mu.
\]

For example, compatibility gives partial_mu log w=Gamma^nu_mu nu and
the contracted derivative of gamma^mu contains -Gamma^mu_mu nu gamma^nu;
their difference is t_nu gamma^nu. Therefore the actual independent-connection
Euler equation is

\[
\boxed{\left(\gamma^\mu D_\mu-m+\tfrac12t_\mu\gamma^\mu\right)\Psi=0.}
\tag{SM.7}
\]

Dropping the torsion-trace term before varying changes the off-shell theory.
It vanishes for the on-shell axial torsion derived below, not for an arbitrary
independent connection. The verifier computes partial L/partial Psi and the
divergence of partial L/partial(partial Psi) separately on rational jets;
omitting t_mu/2 fails both nontrivial fixtures.

## SM-3. Conserved gauge source and a massive causal-current identity

Varying the gauge potential defines

\[
\frac{\delta S_\Psi}{\delta a_\mu}=-w j^\mu,\qquad
j^\mu=-\frac{Zg}{2}\bar\Psi\gamma^\mu\mathcal I\Psi.
\tag{SM.8}
\]

The local phase variation (SM.5), integrated by parts, gives the off-shell
identity

\[
\partial_\mu(wj^\mu)=g E_\Psi^T\mathcal I\Psi.
\tag{SM.9}
\]

Thus the source is conserved on the matter equation, including with an
independent torsionful connection. In the orthonormal frame write

\[
S=\bar\Psi\Psi,\quad P=\bar\Psi\mathcal J\Psi,\quad
V^a=-\bar\Psi\gamma^a\mathcal I\Psi
=\left(\Psi^T\Psi,-\Psi^T\mathcal A_i\Psi\right).
\]

Then

\[
\boxed{\eta_{ab}V^aV^b=-S^2-P^2,\qquad
       g_{\mu\nu}j^\mu j^\nu=-\frac{Z^2g^2}{4}(S^2+P^2).}
\tag{SM.10}
\]

**Native algebra proof.** For Psi=(u,v), S=2u^T Gamma^0 v and
P=2u^T Gamma^0 J v. CP's pure-response identity is
uu^T+(Ju)(Ju)^T=(|u|^2 I+(u^T A_i u)A_i)/2. Since Gamma^0 anticommutes
with every A_i, applying this identity to Gamma^0 u gives

\[
|u|^2I-(u^TA_i u)A_i
=2\left[(\Gamma^0u)(\Gamma^0u)^T
        +(\Gamma^0Ju)(\Gamma^0Ju)^T\right].
\]

Sandwich by v and use CP's pure null identities for u and v. The result is
(SM.10). Hence V is future causal; it is timelike precisely when S^2+P^2>0.
For Z>0 the number current j/g is future causal when g!=0; the sign of the
charge determines the orientation of the charge current. Timelikeness of this
bilinear is not a determination of the supplied parameter m.

Add the declared Maxwell action

\[
S_{\rm EM}=-\frac1{4e_{\rm em}^2}\int w F_{\mu\nu}F^{\mu\nu}d^4x.
\]

Independent a variations give

\[
\frac1{e_{\rm em}^2}\partial_\nu(wF^{\nu\mu})=wj^\mu,\qquad dF=0.
\tag{SM.11}
\]

This is a sourced Maxwell reduction of the declared action. It does not
derive the electromagnetic kinetic law or its measured coupling.

## SM-4. The same action supplies coframe stress and spin

Let l_Psi be the bracketed scalar density in (SM.4), including Z/2, so that
S_Psi=int w l_Psi d^4x. At fixed independent omega, a and Psi,

\[
\tau_a{}^\mu:=\frac1w\frac{\delta S_\Psi}{\delta e^a{}_{\mu}}
=e_a{}^\mu l_\Psi
 -\frac Z2 e_a{}^\nu\bar\Psi\gamma^\mu D_\nu\Psi.
\tag{SM.12}
\]

This follows from delta w=w e_a^mu delta e^a_mu and
delta e_b^nu=-e_b^mu e_a^nu delta e^a_mu. It is the coframe source of the
first-order theory, generally nonsymmetric; it must not prematurely be called
the symmetric metric stress tensor with omega already set to Levi–Civita.

Define gamma^{abc}=gamma^{[a}gamma^b gamma^{c]} and
A^{abc}=bar(Psi) gamma^{abc} Psi, with unit-weight antisymmetrization. The
spin source is totally antisymmetric because H gamma^a is skew:
the vector terms in gamma^c gamma^a gamma^b have zero commuting bilinear.
With all a,b summed, define the three-forms by

\[
\delta S_\Psi=\int\delta e^a\wedge\mathcal T_a
 +\frac12\delta\omega^{ab}\wedge\mathcal S_{ab}
 -\delta a\wedge\boldsymbol j+\delta\Psi^T E_\Psi\,d^4x,
\]
\[
\mathcal T_a=w\tau_a{}^\mu\iota_{\partial_\mu}d^4x,\qquad
\boldsymbol j=w j^\mu\iota_{\partial_\mu}d^4x,\qquad
\boxed{\mathcal S_{ab}=\frac Z4\eta_{aa}\eta_{bb}
                 A^{cab}\iota_{e_c}\mathrm{vol}_e.}
\tag{SM.13}
\]

There is no sum in eta_aa eta_bb here; c is summed. Equivalently, for each
independent a<b, delta S_Psi/delta omega^{ab}_mu is
Zw eta_aa eta_bb e_c^mu A^{cab}/4. Direct variations of all 16 coframe and
24 independent connection coefficients verify both source formulas.

## SM-5. Conditional Einstein–Cartan–Maxwell equations

Use the gamma_H=0 gravitational bulk action already obtained in LR:

\[
S_g=\frac1{4\kappa}\int\epsilon_{abcd}e^a\wedge e^b\wedge R^{cd}
 -\frac{\Lambda}{24\kappa}\int\epsilon_{abcd}e^a\wedge e^b\wedge e^c\wedge e^d.
\tag{SM.14}
\]

On LR's chosen unprojected loop sector, kappa=-1/(2 beta sigma u^2) and
Lambda=-12 sigma u^2. The scalar loop scale u here is unrelated to the
four-component field block called u in SM-1–SM-3. Characteristic and
Nieh–Yan boundary terms do not affect compactly supported bulk variations.

For S_total=S_g+S_EM+S_Psi, the independent coframe and connection equations are

\[
\boxed{\epsilon_{abcd}e^b\wedge
  \left(R^{cd}-\frac\Lambda3e^c\wedge e^d\right)
  =-2\kappa\left(\mathcal T_a^\Psi+\mathcal T_a^{\rm EM}\right),}
\tag{SM.15}
\]
\[
\boxed{\epsilon_{abcd}T^a\wedge e^b=-\kappa\mathcal S_{cd}.}
\tag{SM.16}
\]

For signs and factors, the connection part of delta S_g is
(1/(2 kappa)) int delta omega^{cd} wedge epsilon_abcd T^a wedge e^b.
It follows by integrating epsilon ee D(delta omega) by parts, then moving
the one-form delta omega to the left. Adding (SM.13) gives (SM.16).
The coframe part is (1/(2 kappa)) int delta e^a wedge the left-hand side of
(SM.15). Maxwell has zero independent spin source because F=da.

Equations (SM.7), (SM.11), (SM.15) and (SM.16) form the coupled classical
system. They follow from full independent variations. CP's four regular
thermo gauge channels and MG's sixteen coframe channels can implement those
variations using **independent collections of state maps**. Their proven
surjective variation interfaces survive the addition of this matter action.
The scalar gauge-potential chart is unchanged; its action on the doubled
matter field uses Ical instead of CP's single-copy phase generator J.
Identifying the collections or holding omega fixed by an unvaried
Levi–Civita substitution would change the problem. No new channel-minimality
claim or primitive selection of that constitutive adapter is made.

## SM-6. Unique algebraic torsion and the effective contact term

Write omega=omega_LC+K, with K_abc=K_ab,mu e_c^mu. For invertible e the
linear map T -> epsilon_abcd T^a wedge e^b is invertible: its zero-source
equation is exactly MG's nondegenerate Palatini connection equation, whose
kernel is T=0. This proves uniqueness for any supplied spin source, point by
point; it is not a global existence theorem for the coupled PDE system.

Substitution of the totally antisymmetric source (SM.13) into (SM.16) gives

\[
\boxed{K_{abc}=\frac{\kappa Z}{8}A_{abc},\qquad
       T_{abc}=-\frac{\kappa Z}{4}A_{abc},\qquad t_\mu=0.}
\tag{SM.17}
\]

The sign T_abc=-2K_abc follows directly from T^a=K^a_b wedge e^b for a
totally antisymmetric K. Epsilon contraction of this candidate verifies
(SM.16); invertibility proves it is the full solution, not merely a solution
within an assumed axial subspace. The verifier also forms the quadratic
Hessian in all 24 independent connection coefficients, obtains rank 24,
and checks all 24 stationarity components at (SM.17).

**Eliminating the connection.** The curvature expansion is
R(omega)=R(omega_LC)+D_LC K+K wedge K. The term linear in D_LC K is a
boundary term because D_LC e=0. On the axial solution the remaining
K-dependent scalar density is

\[
\frac{\mathcal L_K}{w}
=-\frac1{2\kappa}K_{abc}K^{abc}+\frac Z8 K_{abc}A^{abc}.
\]

One must first solve the full connection equation; varying only this axial
restriction would not establish the discarded equations. Substituting
(SM.17), which already solves all of them, yields

\[
\boxed{S_{\rm eff}=S_g[e,\omega_{\rm LC}]+S_{\rm EM}
 +S_\Psi[e,\omega_{\rm LC},a,\Psi]
 +\frac{\kappa Z^2}{128}\int w A_{abc}A^{abc}\,d^4x.}
\tag{SM.18}
\]

Since H gamma^{abc} is symmetric, delta A^{abc}=2 delta Psi^T H
gamma^{abc} Psi. Thus the effective matter equation is

\[
\boxed{(\gamma^\mu D_\mu^{\rm LC}-m)\Psi
 +\frac{\kappa Z}{32}A_{abc}\gamma^{abc}\Psi=0.}
\tag{SM.19}
\]

It agrees with inserting (SM.17) directly into (SM.7). This equality is
checked by independent quartic variation and direct connection multiplication.
Setting torsion to zero while omitting the contact term would lose this
interaction. In the effective metric formulation,

\[
G_{\mu\nu}(g)+\Lambda g_{\mu\nu}
 =\kappa\left(T_{\mu\nu}^{\rm EM}+T_{\mu\nu}^{\Psi,\rm eff}\right),
\tag{SM.20}
\]

where T_eff is the metric stress of the eliminated action, including the
contact term and the variation of omega_LC(e). Equivalently it is the
symmetric coframe source after elimination, on the matter equations. Vacuum
Psi=0 returns MG/LR's Einstein–Maxwell system; with F=0 it returns their
vacuum Einstein system. These are conditional action reductions, not proof
that the native axioms uniquely select physical general relativity.

## SM-7. Exact bilinear numbers and a failed cosmological-mass shortcut

The same native matrices give the second commuting Fierz identity

\[
A_{abc}A^{abc}=6\eta_{ab}V^aV^b=-6(S^2+P^2).
\tag{SM.21}
\]

For an explicit algebra check, A^{123}=2u^T Jv and
A^{0ij}=-2u^T Gamma^i Gamma^j v. The coefficient identity

\[
\sum_{i<j}(\Gamma^i\Gamma^j u)(\Gamma^i\Gamma^j u)^T
 -(Ju)(Ju)^T
 =(\Gamma^0u)(\Gamma^0u)^T+(\Gamma^0Ju)(\Gamma^0Ju)^T,
\quad i,j\in\{1,2,3\},
\]

obtained by substituting the native R,K,L matrices, proves (SM.21) after
sandwiching by v and summing the six permutations of each triple. The code
also exhausts all 330 degree-four monomials in the eight real variables for
each of (SM.10) and (SM.21). These are universal finite polynomial identity
checks, distinct from the sampled variable-coefficient field controls.

For example, take

\[
\Psi=(2,0,1,0,\;0,0,3,0)^T.
\]

Then V=(14,-4,0,6), S=-12, P=0, eta(V,V)=-144 and A_abc A^{abc}=-864.
The only nonzero ordered triple with a<b<c is A^{013}=12. Consequently

\[
j^a=Zg(7,-2,0,3),\qquad K_{013}=-\frac32\kappa Z,
\qquad \frac{\mathcal L_{\rm contact}}w=-\frac{27}{4}\kappa Z^2.
\]

These are exact rational benchmark values for a chosen field, not measured
constants or a self-consistent spacetime solution. The complete coupled
Einstein–matter PDE initial/boundary-value problem is not solved here.

There is also a useful negative control. LR's positive-Lambda representation
has the Cartan translation generators R tensor Gamma^a. Under the fixed
similarity needed in SM-2, its naive matter connection contains
-ell^{-1} Ical gamma_a e^a. Slashing it produces -4 ell^{-1} Ical, but
H Ical=gamma^0 is skew, so

\[
\Psi^T H\mathcal I\Psi=0.
\]

The resulting algebraic term disappears from the commuting quadratic action.
Therefore simply calling 4 ell^{-1} the mass does **not** derive the massive
action above in that positive-Lambda sector. The extra return factor permits
both the LR sign extension and SM's massive multiplier, but their existence
does not establish m^2 proportional to Lambda. A further coupling/selection
principle would be needed.

## SM-8. Verification scope and next unresolved selection

`spin_matter.py` imports the unchanged publication fixtures and uses only
Python 3.11/3.12 standard-library rational arithmetic. `verify.py --check`
reruns it along with every earlier package control. New source hashes are
bound into `CERTIFICATE.json`; LR and its implementation are additionally
pinned to the pre-R5 parent commit. Earlier GS/NP/CP/MG/LR proof bytes and
result groups are retained unchanged.

The added checks cover the exhaustive doubled multiplier problem, two full
quartic coefficient identities, independent field/gauge/coframe/spin action
variations, the off-shell Noether identity, local Lorentz and phase covariance,
the full connection Hessian and stationarity equation, torsion signs, and
equality of the two effective matter equations. These support the written
deductions; they are not a proof-assistant certificate, independent review or
experimental validation.

The precise advance is a consistent sourced classical matter completion of
the selected native coefficient action, including the interaction required
when its connection is eliminated. What remains is a native process principle
that selects the module, matter action, adapter, readout and scales; physical
calibration of G, c, Lambda, alpha, hbar and masses; and a quantum/statistical
construction if physical fermions are intended. Additional nonminimal or
pseudoscalar mass couplings have not been excluded by a classification theorem.
Nothing here changes GS's obstruction to a Lorentzian real metric compatible
with an ordinary pseudo-Kähler complex structure on the same tangent space.

## Primary references and attribution

1. A. Perez and C. Rovelli, *Physical effects of the Immirzi parameter in loop
   quantum gravity*, Physical Review D **73**, 044013 (2006),
   [doi:10.1103/PhysRevD.73.044013](https://doi.org/10.1103/PhysRevD.73.044013),
   [arXiv:gr-qc/0505081](https://arxiv.org/abs/gr-qc/0505081). This establishes
   why vacuum arguments for dropping a Holst contribution do not extend
   automatically to a theory with spin sources.
2. L. Freidel, D. Minic and T. Takeuchi, *Quantum gravity, torsion, parity
   violation, and all that*, Physical Review D **72**, 104002 (2005),
   [doi:10.1103/PhysRevD.72.104002](https://doi.org/10.1103/PhysRevD.72.104002),
   [arXiv:hep-th/0507253](https://arxiv.org/abs/hep-th/0507253). This treats
   torsion-induced effective interactions and the dependence on coupling
   choices. Our explicit real commuting-field normalization is derived above;
   it is not imported as a quantum-fermion coefficient from that paper.
