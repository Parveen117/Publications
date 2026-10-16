# MRF-1 — Dynamical material windows and reduced physical profiles

16 October 2026. Frozen parent: RAL-1 at
af76ce357626e1cf6cf221f00e6599841de3b55e. This packet adds a specified
classical canonical reference field to the regular periodic LNC source.
It resolves RAL's fixed-window obstruction within that extension and
derives its physical phase space. It does not assume that the vacuum
quantum representation remains anomaly free after adding a physical field.

## 1. Canonical material reference and full local closure

Keep h=pk+q'phi'-phi'', d=pq'+k phi'-p', with canonical pairs (q,p),(phi,k).
Write v=phi', u=v/p, A_epsilon=p+epsilon v, c=q-log|A_+ A_-|/2.
Use the smooth periodic regular patch A_+ A_->0 with both currents
nowhere zero, x in R/Z. Thus p and P=integral p have one nonzero sign,
and |u|<1.

Add a degree-one orientation-preserving material map Z:S^1_x -> S^1_z,
with lift Z(x+1)=Z(x)+1, Z'>0, and periodic conjugate density rho,
{Z(x),rho(y)}=delta(x-y). Variations of Z are periodic. Set

    J=rho Z',
    H_tilde[N]=integral N(h+u J), D_tilde[f]=integral f(d+J),
    T_tilde_epsilon=(h+uJ+epsilon(d+J))/2
                   =T_epsilon+epsilon A_epsilon J/(2p).     (1)

The reference is dynamical: rho is not constrained to zero. Its canonical
contract and the coefficient u are selected inputs. This is not the
Brown--Kuchar dust action or an established positive-energy matter model.
The densities in (1) use fields and finitely many derivatives at the same
x, on the regular patch; this is the locality asserted for this source.

**MRF-T1 (all-smearing classical closure).** For arbitrary smooth periodic
smearings, the full off-shell identities are

    {D_tilde[f],D_tilde[l]}=D_tilde[fl'-lf'],
    {D_tilde[f],H_tilde[N]}=H_tilde[fN'-Nf'],
    {H_tilde[N],H_tilde[M]}=D_tilde[NM'-MN'].                 (2)

The chiral families therefore commute with one another and each has the
inherited signed Witt bracket.

Proof. D_Z[f]=integral f J satisfies {D_Z[f],D_Z[l]}=D_Z[fl'-lf'],
delta_f Z=fZ', delta_f rho=(f rho)'. The old D makes u a scalar:
delta_f u=f u'. The old/reference mixed spatial-normal terms cancel
the unwanted N f u' J, leaving the weight-two density uJ. This proves
DD and DH using the old LNC identities; the two spatial generators commute.

For HH the old normal transformations give

    delta_N p=(Nv)', delta_N v=(Np)',
    delta_N u=N'(1-u^2)+N[(1-u^2)p'/p-u u'].                (3)

Since u depends only on the mutually commuting p and phi', {u(x),u(y)}=0.
With omega=NM'-MN', the antisymmetrized old/reference bracket from (3)
is integral omega(1-u^2)J: every NM term cancels. The reference-pair
bracket is D_Z[u^2 omega]. Their sum is D_Z[omega], which with the old
D[omega] proves HH. No reference source term is dropped.

The coefficient u makes Z track the same gauge transport as p. In the
more general local ansatz h+b(u)J, d+J, normal closure would require
b^2+(1-u^2)b_u=1. It allows b=u and b=+1 or -1. Requiring p and Z to
have the same normal transport velocity selects b=u within this ansatz.
Closure alone does not select the measurement contract; b=1 generally
fails the window test.

## 2. Exact cancellation of the window residual

Let w(z) be any real smooth periodic material profile, possibly supported
in a proper material interval. For real smooth a(c,P), P!=0, put

    G_(a,w)=integral p(x) a(c(x),P) w(Z(x)) dx,
    Zcal=(A_-/A_+)T_tilde_+ +(A_+/A_-)T_tilde_-.             (4)

w is prescribed as a function of the dynamical material coordinate;
Z and rho are varied canonically. This differs from fixing w(x).

**MRF-T2 (material-window first-class coupling).** Exactly,

    {G_(a,w),T_tilde_epsilon[f]}
                     =(1/2)integral f w(Z) a_c Zcal.         (5)

No w_z residual remains. With the global clock/pointer pairs (T,E),(X,Pi),
C_T=E+g(T)Pi G_(a,w) closes with all the new local constraints. The action is

    S=integral ds [integral(p dot q+k dot phi+rho dot Z)dx
        +E dot T+Pi dot X
        -H_tilde[N]-D_tilde[beta]-nu C_T].                 (6)

All its secondary constraints are first class on the patch. The reference
is part of the canonical source, not an external label transformation.

Proof. Put V_epsilon=epsilon f A_epsilon/(2p). The added term in (1)
leaves the old transformations of p,phi,A_+,A_- unchanged. It changes q
by -f vJ/(2p^2) relative to the old chiral flow. Direct calculation gives

    delta_epsilon p=(V_epsilon p)',
    delta_epsilon Z=V_epsilon Z', delta_epsilon P=0,
    c'=T_tilde_+/A_+ - T_tilde_-/A_- - J/p,
    delta_epsilon c=V_epsilon c'+f Zcal/(2p).               (7)

Consequently

    delta_epsilon[p a(c,P)w(Z)]
       =partial_x[V_epsilon p a(c,P)w(Z)]
                              +(f/2)w(Z)a_c Zcal.           (8)

The periodic derivative integrates to zero. Its w_z part is canceled by
the canonical variation of Z. Equations (2), (5) and {C_T,C_T}=0 prove
closure and conditional constraint preservation along smooth solutions
within the patch.

The new source feedback, in addition to the old LNC field equations, is

    dot q|_ref=-N vJ/p^2, dot k|_ref=partial_x(NJ/p),
    dot Z|_ref=(Nu+beta)Z',
    dot rho|_ref=partial_x[(Nu+beta)rho].                    (9)

The apparatus term additionally gives

    dot Z|_app=0, dot rho|_app=-nu g Pi p a w_z(Z).          (10)

Thus imposing rho=0 as an extra constraint would generally destroy this
construction: (10) does not preserve it for a nonconstant window. The
reference momentum carries required backreaction.

## 3. Complete regular material-coordinate reduction

The reference introduces one physical canonical pair per material point
after the two local gauge families are removed. It is not a pure gauge
addition leaving the vacuum physical space unchanged.

For z in S^1, let x_z=Z^{-1}(z), and define on the local constraint surface

    m(z)=[p/Z'](x_z),
    c_mat(z)=c(x_z), r(z)=c_mat(z)+log|m(z)|.                (11)

The logarithms use inherited dimensionless units. m has one nonzero sign,
and integral m dz=P.

**MRF-T3 (physical profiles and canonical quotient).** The profiles r,m
are invariant under both local gauge families. Every regular constraint
orbit in the degree-one Z patch has a representative

    Z=x mod 1, phi=tau, k=0,
    p=m(x), q=r(x), rho=m'-m r'.                            (12)

The profiles distinguish these orbits. At fixed tau, the quotient one-form
and smeared canonical brackets are

    Theta_red=integral m delta r dz,
    {integral f r, integral l m}_red=integral f l dz.        (13)

The r-r and m-m brackets vanish. This is a classical physical phase space;
it supplies neither an infinite-dimensional Lebesgue measure nor a quantum
physical inner product.

Proof. On shell (7) transports c as a scalar and p,Z' as densities with
the same velocity V_epsilon. Hence p/Z' is a scalar with that velocity.
The inverse-coordinate variation delta x_z=-V_epsilon(x_z) cancels
both transports in (11), proving actual gauge invariance.

For existence of (12), the currents retain exactly the LPS independent
chiral density transformations. Straighten both to A_+=A_-=P using its
degree-one circle diffeomorphisms. Along those smooth finite gauge paths,
the current patch stays regular; Z and rho obey scalar and density
transport with smooth velocity. Z therefore remains a degree-one monotone
map. The remaining field equations along these specified gauge paths have
smooth coefficients and sources, so the straightening still applies.
Then phi'=0 and H_tilde=Pk enforces k=0. Residual constant normal smearing
sets the additive phi to tau. Finally a common spatial diffeomorphism
puts Z=x, preserving phi'=k=0. D_tilde=0 now gives rho=p'-pq'.

Conversely arbitrary smooth r and definite-sign m inserted in (12) solve
both constraints, and (11) returns exactly r,m. Representatives with the
same profiles coincide; invariance forbids gauge equivalence between
different profiles. This proves classification under the connected smooth
gauge flows used here, including circle rotations. Material labels are
fixed calibration labels; there is no extra quotient by target relabelings.

The restricted one-form is (13), because delta Z=delta phi=0 and k=0.
This gives the bracket. The gauge is also locally transverse: there
delta Z=beta and delta phi=Np, so preserving both gauges forces N=beta=0.
Without the apparatus pulse the scalar-clock reduced Hamiltonian is zero
for this reference model. The profiles remain physical labels, not gauge.

## 4. Reduced apparatus, feedback and precise locality boundary

The reduced coupling is

    G_red[r,m]=integral m(z) w(z) a(r(z)-log|m(z)|,P) dz,
    P=integral m dz.                                       (14)

Write J_P=integral m w a_P dz, with a_P holding the first argument fixed.
Set the new clock T=s, h_T=g(T)Pi, and E+g Pi G_red=0.

**MRF-T4 (physical feedback and material localization).** The exact equations are

    dot r=h_T[w(a-a_c)+J_P], dot m=-h_T m w a_c,
    dot X=g G_red, dot Pi=0.                               (15)

They preserve the reconstructed reference constraint:

    d/dT(m'-m r')=-h_T m a w'.                             (16)

This equals (10) in the material gauge. The nonconstant window is
consistent because reference momentum responds.

Proof. At fixed z, delta c_mat=delta r-delta m/m and delta P=integral
delta m. The functional derivatives of (14) are m w a_c and
w(a-a_c)+J_P, proving (15). Differentiate rho=m'-m r'. Because
a'=a_c(r'-m'/m) and J_P'=0, the a_c' terms cancel and (16) remains.
No field equation is discarded.

For the RPM measuring function a(c,P)=1+(1/2)cos(c+log|P|), write

    theta=r+log|P/m|, a=1+(1/2)cos theta, b=-(1/2)sin theta,
    B=integral m w b, J_P=B/P.                             (17)

This realizes the original measuring function with a material window.
It still depends on global P. Outside supp w, dot m=0 but dot r=h_T B/P
can be nonzero. A material window does not make this entire coupling local.

Alternatively select a=a(c) independent of P. Then J_P=0 and both rates
in (15) vanish outside supp w. The density p a(c)w(Z) has only finite-jet
field dependence on the patch; a(c)=1+(1/2)cos c is one example.
This is a different observable from RPM's function. The pointer and
pulse clock remain global pairs. No propagation speed, positive-energy
apparatus or spacetime causal completion follows from the support statement.

## 5. Controlled classical comparison with RPM

**MRF-T5 (readout recovery and homogeneity test).** For (17), start at
r(z)=Q, m(z)=P!=0 and normalize integral w dz=1. The exact classical
pointer output through a finite pulse is

    X_out-X_in=kappa P[1+(1/2)cos Q], kappa=integral g dT,   (18)

with initial Q,P on the right. This matches the classical RPM readout
even for nonconstant normalized material windows. It does not identify
the full field dynamics or the quantum theories.

Proof. Initially theta=Q, and (14) equals the probe on the right of (18).
G_red has no explicit T dependence and is the Hamiltonian's sole system
factor, so dot G_red={G_red,g Pi G_red}=0. Integrating dot X proves (18)
including all field backreaction.

For w=1 the homogeneous profiles form an invariant classical submanifold:
J_P=b and (15) give Qdot=g Pi a, Pdot=-g Pi P b. With Z=x,rho=0 this
is exactly the homogeneous RAL/RPM extension. Applying its inherited
homogeneous quantization gives the old RPM model, but does not construct
an invariant quantum subspace or quantization of the new full field theory.

For nonconstant w and g Pi!=0, homogeneity fails already in the initial
vector field: if b!=0, dot m=-g Pi P b w is nonconstant; if b=0,
dot r=g Pi a w is nonconstant because a>0. The reference therefore cannot
be discarded merely because an initial pointer reading agrees.

For an explicit short-pulse witness choose Q=0,P=2 and
w(z)=1+(1/2)cos(2 pi z), giving a=3/2,b=0 initially. With pulse area
s=integral g and constant Pi,

    partial_s r|_0=(3/2)Pi w, partial_s m|_0=0,
    partial_s rho|_0=-3 Pi w'.                             (19)

The spatial response starts at first order, although the full classical
pointer shift is exactly 3 kappa. The native Fourier fixture displays
the common 2 pi spatial-derivative normalization.

For completeness these reduced comparisons exist through finite pulses
for all smooth initial definite-sign m and smooth r. Since |b|<=1/2,

    |log(|m(T,z)|/|m(T_0,z)|)|
       <=(|Pi|/2)||w||_infinity integral |g|dT.              (20)

Thus m and P stay away from zero on finite intervals. Also
|B/P|<=(1/2)||w||_infinity and
|dot r|<=(5/2)|g Pi| ||w||_infinity. The vector field (15),(17) is locally
Lipschitz on each periodic C^k profile space, k>=0, away from m=0,P=0.
These bounds give C^0 continuation. Differentiation in z gives equations
linear in the highest spatial derivatives with bounded coefficients and
lower-order sources. Induction and Gronwall bound every finite C^k norm
on finite intervals. Consistent C^k solutions therefore yield smooth
solutions from smooth initial profiles. This is a reduced nonlocal ODE
argument, not a general existence theorem for arbitrary unreduced gauges.

## 6. New endpoint and quantum obligation

Closed here: the specified local classical reference-field action and
arbitrary-smearing first-class algebra; canonical window-residual
cancellation; complete material-coordinate reduction and physical profile
brackets; rod/geometry feedback; exact classical pointer comparison and
the distinction between readout recovery and homogeneous dynamics.

The source gained physical field degrees of freedom. LQA's old anomaly
cancellation and LPS/CCS/COS measures cannot simply be reused as its quantum
theory. Next: construct a positive quantum representation and domains for
the reduced profiles, establish a controlled continuum/regularization
contract, and compare to quantization of the extended constraints.
The old homogeneous model is a comparison target, not a proof of that
continuum construction. Critical patches, fully local apparatus dynamics,
physical 3+1 recovery and calibration remain open. Predecessors are frozen.
