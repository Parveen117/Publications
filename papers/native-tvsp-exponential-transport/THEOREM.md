# TET-1 — TVSP layers generate exponential observer transport

16 October 2026. Continuation of TLM-1, NT-1/3/4/8 and EGR-1.

## Contract and result

Retain TLM's normalized affine chart (s,v), supplied smooth constitutive
potential u, positive Hessian H=D²u, six actual response components Lambda,
and the selected connection A_i=(1/2)H^(-1) partial_i H. The connection is
selected in that source trivialization and subsequently transformed as a
connection. It is not selected afresh after an arbitrary moving frame change.
All paths below remain in the regular positive chart. The movement parameter
t is independent of response lambda, deformation mu and derivative depth N;
it has not been identified with a physical clock.

This packet connects retained source layers to transport layers, to a moving
observer and positive cut, and to a curvature-sensitive return measurement.
The exponential, parallel transport, metric compatibility and jet calculus
are established mathematical constructions reused here. The contribution to
this programme is their source-bound assembly and explicit TVSP examples,
including the cubic ordering defect F_sv(0)/6 on a curved source path.
No new general theory of exponentials or physical space quantization is claimed.

The continuum statements have direct written proofs below. Exact native
fixtures and serialized word identities are reproducible checks, not a
proof-assistant verification of all smooth sources.

## TET-T1. Source-jet to transport-jet reconstruction

Let gamma be a supplied smooth path with gamma(0)=p. Write H(t)=H(gamma(t))
and define

    A(t)=dot(gamma)^i A_i(gamma(t))=(1/2)H(t)^(-1)H'(t).
    U'=-A U,                         U(0)=I.                       (1)

Then U is invertible and is the ordered exponential of -A. Later path
segments multiply on the left. If A(t) commutes with A(r) at every pair
of parameter values, U=exp(-integral A). In particular constant A gives
U=exp(-tA). No exponential coordinate trajectory is assumed by this formula.

For N>=1, the path jet and H source jet through degree N determine the
transport jet through degree N. With normalized coefficients
H(t)=sum H_n t^n and H(t)^(-1)=sum B_n t^n, the exact recurrences are

    B_0=H_0^(-1),
    B_n=-H_0^(-1) sum_(j=1)^n H_j B_(n-j),
    A_n=(1/2) sum_(j=0)^n B_j (n-j+1) H_(n-j+1),
    U_0=I,
    (n+1)U_(n+1)=-sum_(j=0)^n A_j U_(n-j).                       (2)

Only A_0,...,A_(N-1) are needed. Products retain the printed matrix order.
Degree zero alone determines U(0), but not its first derivative. Unknown
higher coefficients are never set to zero and called exact derivatives.

**Proof.** On a compact parameter interval where A is continuous, successive
substitution in U=I-integral A U gives the ordered integral series. If
||A||<=M in an induced norm, the nth term has norm at most (Mt)^n/n!,
so the series converges uniformly and solves (1). The same construction
for V'=V A, V(0)=I gives (VU)'=0, proving invertibility and uniqueness.
Pairwise commutativity permits the integration simplexes to combine into
the ordinary exponential series. Equating coefficients in HB=I,
A=(1/2)BH' and (1) gives (2). Taylor composition along a supplied path
uses only its and the source's retained coefficients through the requested
degree. This proves determination, not convergence of an arbitrary smooth
Taylor series at a nonzero endpoint. TET-T6 states the endpoint boundary.

## TET-T2. Moving observers, metric transport and positive cuts

For (1), the following identities hold throughout the path:

    U(t)^T H(t) U(t)=H_0,
    w(t)=U(t)w_0,      O(t)=O_0 U(t)^(-1),
    O(t)w(t)=O_0 w_0.                                             (3)

Here w is a tangent carrier vector and O its dual observer. A fixed
observer O_0 need not give a constant reading on U(t)w_0.

If P_0²=P_0 and P_0^T H_0=H_0 P_0, then

    P(t)=U(t)P_0U(t)^(-1)                                        (4)

is an H(t)-orthogonal cut. Its rank, trace and spectrum are unchanged,
and P'=-[A,P]. For a constant source weight q>=0, rho=qP retains that
weight and its discrete cut spectrum. In rank one Tr(rho)=q, as in EGR
section 1; in general Tr(rho)=q rank(P_0). Neither q nor a physical energy
coupling is derived by this transport statement.

**Proof.** H'=A^T H+HA follows from symmetry of H and H'. Differentiating
U^T H U therefore gives zero. Observer cancellation proves the other
identity in (3). Similarity proves idempotence, rank and spectrum in (4).
Substituting H(t)=U^(-T)H_0U^(-1) proves H(t)-self-adjointness.
Differentiating (4) proves its commutator evolution.

To connect exactly to EGR's fixed positive pairing, choose a smooth positive
factor S(t) with S(t)^T S(t)=H(t), and put

    W(t)=S(t)U(t)S(0)^(-1),
    V(t)=W'(t)W(t)^(-1)=S'S^(-1)-S A S^(-1).                     (5)

Equation (3) gives W^T W=I, hence V^T=-V. The cut in this orthonormal
frame is W [S(0)P_0S(0)^(-1)] W^T. In rank one this is precisely EGR's
normalized positive-cut transport contract. The original U generally is not Euclidean
orthogonal. This construction uses the positive pointwise tangent pairing,
not a positive representation of TLM's nilpotent jet algebra.

There is an additional simplification in this oriented rank-two sector.
Take the positive-diagonal triangular factor
S=[[sqrt(a),b/sqrt(a)],[0,sqrt(Delta/a)]]. NT-4 directly computes

    S A S^(-1)-dS S^(-1)=R alpha,
    R=[[0,-1],[1,0]],
    alpha=[db-(b/a) da]/(2 sqrt(Delta)),
    W_gamma=exp[-R integral_gamma alpha].                        (5a)

Equivalently the full source-coordinate transport is

    U_gamma(t)=S(gamma(t))^(-1)
               exp[-R integral_gamma alpha] S(gamma(0)).         (5b)

Indeed the transformed connection is skew, hence a multiple of the single
constant generator R; its lower-left entry gives alpha. Consequently the
orthonormal-frame exponential is ordinary in rank two. Coordinate-frame
generators may still fail to commute: omitting ordering for the original
A is not the frame change (5). This is EGR's transport mechanism on the
declared tangent carrier, not an identification with EGR's Regge coupling.

For an arbitrary smooth invertible change of tangent frame g(t),
A_tilde=g^(-1)Ag+g^(-1)g' and U_tilde=g(t)^(-1)U(t)g(0).
Direct differentiation proves the law. Replacing A_tilde by the bare
strain formula in the new frame would drop the frame-motion term.

## TET-T3. Transport of the actual response layers

TLM stores the six-response section in a rank-six fibre E with its declared
flat connection. Its kth covariant derivative has type E tensor (T*)^k.
The two-by-two U in (1) acts on the tangent slots; it is not a six-by-six
action on Lambda. Parallel transport of an initial covariant tensor T_0 is

    T_parallel(t)(z_1,...,z_k)
       =T_0(U(t)^(-1)z_1,...,U(t)^(-1)z_k).                       (6)

For an actual source response tensor T, evaluated along gamma,

    d/dt [T(gamma(t))(U v_1,...,U v_k)]
      =(nabla_dot(gamma) T)(gamma(t))(U v_1,...,U v_k).             (7)

Thus source variation survives removal of frame motion. An actual scalar
response is not forced to be constant; at k=0, (7) is its ordinary path
derivative. This also specifies how the retained response tower is used.

**Proof.** Differentiate each slot in (6), using (U^(-1))'=U^(-1)A.
Differentiate the left side of (7), using (Uv)'=-A Uv. The negative
connection terms are exactly those of the covariant derivative on the
lower slots. The E fibre has no extra connection term in the chosen frame.

More explicitly, use TLM's indexing L_1=Lambda and L_(k+1)=nabla L_k.
Define the pulled-back layer Ltilde_k=L_k(gamma) composed with U in its
k-1 lower slots, and xi=U^(-1) dot(gamma). Then (7) gives

    (d Ltilde_k/dt)(v_1,...,v_(k-1))
       =Ltilde_(k+1)(xi,v_1,...,v_(k-1)).                         (7a)

Here the first slot is the newest derivative; the remaining slots are
the displayed tensor's arguments. This connects actual tower layers to
observer motion. A finite tower does not close its own dynamics unless
its next layer is supplied, determined by a further law, or bounded.

For a smooth endomorphism field M, j^N(U^(-1) M U) is the ordered product
(j^N U)^(-1) (j^N M) (j^N U) in the truncated jet algebra. TLM's faithful
regular multiplication encoding preserves this identity. A metric uses
congruence U^T H U instead. Applying similarity to a scalar cannot uncover
its source context, and the encoded formal adjoint is not ordinary matrix
transpose. These are consequences of TLM-T1/T2, not a new positive *-model.

## TET-T4. Exact exponential paths of the existing TVSP source

Use the unchanged TLM source, with |s|,|v|<1/20 and |mu|<=5:

    u_mu=10s-10v+s²+sv+(5/2)v²+(mu/2)s²v,
    H_mu=[[2+mu v,1+mu s],[1+mu s,5]],
    H_0=[[2,1],[1,5]],
    H_s=[[0,mu],[mu,0]],     H_v=[[mu,0],[0,0]],
    f=-lambdaS=5(1+v),      Delta=5(2+mu v)-(1+mu s)².             (8)

For a straight path gamma(t)=(a t,b t), set D=aH_s+bH_v and B=H_0^(-1)D.
While this path stays in the patch,

    A(t)=(1/2)(I+tB)^(-1)B,
    U(t)=(I+tB)^(-1/2)=exp[-(1/2)log(I+tB)].                       (9)

All generators on this particular path commute. The inverse square root
and logarithm are the positive H_0-self-adjoint branches. Near t=0 the
normalized transport coefficients are binom(-1/2,n) B^n.

**Proof.** H=H_0(I+tB), and B^T H_0=H_0 B. Conjugation by a positive
factor of H_0 makes B real symmetric. Positivity of H on the path gives
positive eigenvalues of I+tB, which defines the branches in (9).
Differentiating each scalar spectral component proves (1) and (9), or
equivalently the commuting power-series identity proves it locally and
the positive branch extends it along the path. Formula (3) follows too.

On the v-axis, introduce Q=[[1,0],[-1/5,0]], Q²=Q. This is an algebraic
idempotent, not an ordinary orthogonal observation cut. Then

    Delta(t)=9+5mu t,
    U(t)=exp[-(1/2)log(Delta(t)/9) Q]
        =I+(sqrt(9/Delta(t))-1)Q.                                (10)

For mu=5 and t=369/10000, Delta/9=(21/20)² and U=I-Q/21, an exact
rational endpoint inside the certified patch. With the H_0-orthogonal cut
P_0=[[1,1/2],[0,0]], (4) gives a concrete positive transported cut.
The actual mechanical response on this path is 5(1+t), not a preserved 5.

On the s-axis instead, f=5 at every point, while the ambient curvature is

    kappa(t)=5mu² / [4(10-(1+mu t)²)²].                           (11)

This is the two-dimensional source curvature sampled along the path, not
intrinsic curvature of a one-dimensional path. Equations (9) and (11)
exhibit moving source frames and changing geometry with the same displayed
scalar reading. Constancy of that one reading does not determine transport.

## TET-T5. Ordered exponential correction and loop recovery

Write C=H_0^(-1)H_s and D=H_0^(-1)H_v. The retained first source layers
give the already established selected curvature

    F_0=F_sv(0)=-[C,D]/4
       =mu²/324 [[5,25],[-10,-5]].                               (12)

For the curved path gamma(t)=(t,t²) in the same source, let
A(t)=A_0+tA_1+t²A_2+.... Then

    A_0=C/2,     A_1=D-C²/2,     [A_0,A_1]=-2F_0,
    U(t)-exp(-integral_0^t A(r)dr)=t³ F_0/6+O(t^4).              (13)

Thus even this single holonomic source supplies a nonzero error from
discarding ordering. For mu=1 the cubic coefficient is
[[5,25],[-10,-5]]/1944. The ordinary exponential in (13) uses the full
integrated generator, not just its constant term.

**Proof.** Along the path H=H_0(I+tC+t²D), so expansion of (1) gives
the stated A_0,A_1. Recurrence (2) gives

    U_1=-A_0,
    U_2=(A_0²-A_1)/2,
    U_3=-A_0³/6+A_0A_1/6+A_1A_0/3-A_2/3.

Expanding exp(-tA_0-t²A_1/2-t³A_2/3+...) instead gives mixed coefficients
(A_0A_1+A_1A_0)/4. Their difference is -[A_0,A_1]/12=F_0/6.
This is an exact coefficient identity; the fixture checks the source-derived
coefficients rather than fitting sampled numerical trajectories.

For a small rectangle based at p, following +s,+v,-s,-v with side lengths
h,k and later segments acting on the left, the full return is

    U_box=I-hk F_sv(p)+O(|h|²|k|+|h||k|²).                       (14)

**Proof.** Set X=A_s(p), Y=A_v(p), and retain each edge's mixed-degree
terms. The four edges in chronological order have terms

    E1=I-hX,
    E2=I-kY-hk partial_s Y,
    E3=I+hX+hk partial_v X,
    E4=I+kY.

Multiplying E4 E3 E2 E1, the mixed coefficient is
-partial_s Y+partial_v X-[X,Y]=-F_sv. Pure edge-square terms cancel
in the full expansion and cannot affect the mixed coefficient. Smoothness
and U_box=I when either side vanishes give the stated remainder locally.
This is NT-2's sign convention, now recovered from the actual source jets.
Reversal reverses the leading sign. Keeping only constant A_s,A_v would
discard the derivative terms and give the wrong source curvature.

A fixed base observer O_0 on a returned w_0 reads the leading difference
-hk O_0 F_sv w_0. Sufficient independent base readings recover F_sv, and
the known oriented H and Delta give kappa=(HF_sv)_12/Delta. A single
observer may annihilate the signal. In particular, co-transporting BOTH
observer and vector makes their pairing constant by (3) and hides this
return. Curvature recovery compares the returned frame/vector with the
retained base frame; it does not contradict moving-observer invariance.

There is also a finite-loop formula in this rank-two sector. For a
rectangle bounding a surface Sigma inside the positive source chart,
with the same orientation and base p, equations (5a) and NT-4 give

    d alpha=-kappa sqrt(Delta) ds wedge dv,
    U_box=S(p)^(-1) exp[R integral_Sigma kappa sqrt(Delta) ds dv] S(p).
                                                                    (14a)

Integrating partial_s alpha_v-partial_v alpha_s over the rectangle and
using the one-dimensional fundamental theorem on its two coordinates
gives integral_boundary alpha=integral_Sigma d alpha, proving (14a).
For the source (8), direct differentiation gives
d alpha=-5mu²/(4 Delta^(3/2)) ds wedge dv, agreeing with its kappa.
The formula recovers the finite return from the full supplied source;
an endpoint rotation alone determines its phase only modulo 2pi. Its
leading term agrees with (14), since S F_sv S^(-1)=-R kappa sqrt(Delta).

## TET-T6. Finite-depth limits and controlled finite transport

For every finite N>=0, the entire retained H and six-response jets through
degree N need not determine transport at any specified positive endpoint.
An actual TVSP counterexample is the pair eta=0,1:

    u_eta=10s-10v+(s²+v²)/2
          +eta s^(N+3)/[(N+2)(N+3)],
    H_eta=diag(1+eta s^(N+1),1),      S=1+s, V=1+v.                (15)

These sources are positive near zero. Their H and actual six-response
jets through degree N coincide. Along gamma(t)=(t,0),

    U_eta(t)=diag((1+eta t^(N+1))^(-1/2),1).                      (16)

The two U jets through N coincide and their first differing coefficient
is -eta/2 at degree N+1 in the first channel. Their finite positive
endpoint values differ, however small that endpoint is.

**Proof.** The added potential begins at order N+3, so its Hessian begins
at N+1; its temperature change starts at N+2 and its pressure is unchanged.
All denominators in NT's six response formulas are nonzero at zero.
Products and inversion therefore preserve equality of the degree-N jets.
Direct integration of the diagonal generator gives (16).
This supplies a holonomic, stable-source obstruction, not arbitrary formal
coefficients. It also makes the derivative-depth boundary sharp.

A finite approximation can nevertheless have a proved error bound when
interval information is supplied. Let U_N=sum_(n=0)^N U_n t^n, define its
residual R_N=U_N'+A U_N, and assume, on 0<=t<=T in an induced norm,

    ||A(t)||<=M,             ||R_N(t)||<=C t^N.

Then

    ||U(t)-U_N(t)|| <= C exp(Mt) t^(N+1)/(N+1).                  (17)

**Proof.** E=U-U_N satisfies E'=-AE-R_N, E(0)=0. The evolution between
r and t has norm at most exp(M(t-r)), from the same ordered series as
TET-T1. Variation of constants bounds E by
integral_0^t exp(M(t-r)) C r^N dr, which implies (17).
Finite retained jets alone do not supply the interval-wide constants M,C.
For the full explicit source of TET-T4 the closed formula supplies actual
finite transport; it is not inferred from an unspecified truncated record.

## Research boundary and next obligation

This closes the declared source-layer/exponential-observer transport bridge.
It does not select the physical path gamma, establish a clock, quantize
spatial distance, or construct a 3+1 gravitational quantum state space.
The EGR cut transport is realized on a positive tangent carrier; no
intertwiner to the distinct LMD/MPQ continuum quantum representation has
been constructed. Its constraints, domains and geometry obligations remain
at their existing endpoints.

The next physical bridge must identify a source-selected trajectory or
observable in that existing quantum sector, supply its domain and show that
this transport or its controlled limit acts there. Adding formal layers
or asserting that an exponential is automatically unitary cannot replace
that identification.
