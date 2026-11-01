# MPQ-1 — Positive continuum quantization of the reduced material profiles

16 October 2026. Frozen parent: MRF-1 at
839f4cf99a3f32c34a88e0c3138b4d7fff301815. The material action, regular
patch, reference field and RPM measuring function are unchanged. This
packet supplies a selected quantization **after** MRF's proved reduction.
It closes the reduced state/observable/pulse problem under the contract
below, not quantization of the original unreduced local constraints.

The new ingredients are an explicit source canonical transformation, a
global pulse-adapted physical chart, a specified Gaussian representation
of its affine field algebra, and a convergent Fourier cutoff. The Gaussian
covariance and observable contract are choices, not consequences of the
classical constraints. In particular this representation does not inherit
GOA's orientation-mixing matching. The comparison obstruction is proved.

## 1. Source audit and explicit canonical Abelianization

All fields are periodic on the unit circle except the lift of the
degree-one monotone material map Z, for which Z(x+1)=Z(x)+1. Keep MRF's
p of one nonzero sign, v=phi', u=v/p with |u|<1, and Z'>0. Put

    c=q-(1/2)log(p^2-v^2),
    K_phi=k-partial_x atanh(v/p),
    C_Z=p c'+rho Z'.                                      (1)

**MPQ-T1 (canonical source bridge).** The change (q,p,phi,k) to
(c,p,phi,K_phi) is canonical. MRF's actual constraints become exactly

    H_tilde=p K_phi+u C_Z, D_tilde=v K_phi+C_Z.             (2)

Next set c_m(z)=c(Z^{-1}(z)), m(z)=(p/Z')(Z^{-1}(z)) and
R_Z(x)=rho(x)+p(x)c'(x)/Z'(x). Then (c_m,m), (phi,K_phi),
(Z,R_Z) are canonical pairs. On this open chart the original constraints
are equivalent to K_phi=R_Z=0, whose Poisson brackets vanish strongly.
This is a classical statement about the exact source, not a quantum
anomaly calculation or a global translation group on the restricted chart.

Proof. Use the type-two functional

    F_2=integral [q p+phi K_phi+F(p,phi')] dx,
    F(p,v)=-(p/2)log(p^2-v^2)+p-v atanh(v/p).

Its derivatives are F_p=-(1/2)log(p^2-v^2) and F_v=-atanh(v/p).
Thus the new coordinate is c, the old k is K_phi+(atanh u)', and
the transformation is canonical, with the usual exact one-form term.
Substitution gives (2); its determinant is p(1-u^2), nowhere zero.
At fixed old x, delta c(x)=delta c_m(Z(x))+c_m'(Z(x))delta Z(x).
Changing variables x to z in integral p delta c consequently gives

    integral(p delta c+rho delta Z)dx
       =integral m delta c_m dz+integral R_Z delta Z dx.   (3)

The phi pair is unaffected by this second transformation. The new gauge
momenta therefore commute, and C_Z=Z'R_Z. MRF's physical r is
r=c_m+log|m|, so integral m delta r=integral m delta c_m+delta P,
P=integral m. This recovers exactly its reduced symplectic form.

The chart restrictions still couple the source fields: arbitrary finite
translations of phi or Z can leave |phi'|<|p| or Z'>0. It would be wrong
to replace this restricted gauge chart by an unrestricted flat gauge
Hilbert space and claim an unreduced quantum equivalence from (2) alone.

## 2. Global physical chart and exact pulse solution

Fix the existing measuring function, without changing its global P term:

    theta(z)=r(z)+log|P/m(z)|,
    a(theta)=1+(1/2)cos(theta),
    Y(theta)=integral_0^theta dt/a(t),
    y(z)=Y(theta(z)), j(z)=m(z)a(theta(z)).                (4)

Here theta and y are real lifts, not angles identified modulo 2 pi.
Since 1/2<=a<=3/2, Y is a smooth increasing diffeomorphism R -> R.

**MPQ-T2 (global symplectic linearization).** On either sign component,
(4) is a global bijection of the smooth physical profile phase space,
and its canonical one-form is exactly

    integral m delta r=integral m delta theta
                      =integral j delta y.               (5)

Every material-window coupling in the frozen family is

    G_w=integral w(z)j(z)dz.                              (6)

Thus all G_w Poisson commute. For pulse area alpha(T)=integral g dT,
fixed pointer Pi, and an initial slice with alpha=0, the entire reduced
classical solution is

    j_T=j_0, y_T=y_0+Pi alpha(T)w,
    theta_T=Y^{-1}(y_T), m_T=j_0/a(theta_T),
    P_T=integral m_T, r_T=theta_T+log|m_T/P_T|,
    X_T=X_0+alpha(T)G_w.                                 (7)

Proof. In delta theta the term -delta m/m integrates against m to
-delta P, while delta log|P| integrates to +delta P. They cancel,
proving the first equality in (5); the second uses dY=dtheta/a.
Conversely choose arbitrary smooth y and definite-sign j. Then
theta=Y^{-1}(y), m=j/a(theta), P=integral m and
r=theta+log|m/P| give the unique inverse. No zero of m or P occurs.
Equations (6),(7) now follow from Hamilton's equations for
C_T=E+g(T)Pi G_w. In particular the old apparently nonlocal feedback has
not been deleted: reconstructing r through P_T restores it exactly.
Nonconstant w still excites the original profiles, as MRF proves.

For quantization put epsilon=sign(j)=sign(m), and use the further chart

    s(z)=log|j(z)|, d(z)=j(z)y(z).

It obeys {d(h),s(f)}=integral h f, {d(h),j(f)}=j(hf), with
d(h)=integral h d and j(f)=integral f j. The d-d and j-j brackets vanish.
Indeed j delta y=delta d-d delta s, so d is the coordinate conjugate
to the logarithmic momentum s. The pair s,d separates every classical
physical profile: j=epsilon exp(s), y=d/j and the inverse above applies.
It is this complete separating chart, not a claim to quantize every
nonlinear classical function, that defines the observable contract below.

## 3. A genuine continuum state space and physical observables

Set hbar=1 until section 6. In the real orthonormal Fourier basis
e_0=1, e_n^c=sqrt(2)cos(2 pi n z), e_n^s=sqrt(2)sin(2 pi n z),
select covariance C with lambda_0=1 and lambda_n=(1+n^2)^(-4) for
both n modes. Let mu be the probability law of independent real Gaussian
coefficients s_j of variances lambda_j. Our reduced Hilbert spaces are

    H_epsilon=L^2(mu), H_red=H_+ direct_sum H_-.            (8)

Both signs use the same chosen covariance and conventional unit weights.
One sign can also be treated separately. This is a correlated continuum
measure, not a formal product of Lebesgue measures at material points.

**MPQ-T3 (positive continuum affine representation).** The series for s
converges mu-almost surely in H^3(S^1), hence in C^2. In particular
j=epsilon exp(s) is a continuous, nowhere-zero physical density with
one sign on every sample. For real smooth f, the operators

    S(f)psi=(integral f s)psi,
    J(f)psi=epsilon(integral f exp(s))psi                 (9)

are self-adjoint on their maximal multiplication domains, respectively
{psi:S(f)psi in L^2} and {psi:J(f)psi in L^2}. For real h in the
Cameron--Martin space H_C=H^4, define

    beta_h(s)=sum_j h_j s_j/lambda_j,
    ||h||_C^2=sum_j h_j^2/lambda_j,
    U_h(t)psi(s)=exp[-t beta_h(s)/2-t^2||h||_C^2/4]
                  psi(s+t h).                            (10)

These are strongly continuous unitary groups. Their self-adjoint
generators, U_h(t)=exp(-it D(h)), have the explicit domain

    Dom D(h)={psi:lim_(t->0) i[U_h(t)psi-psi]/t
                         exists in L^2(mu)}.             (11)

For smooth h their action on the common invariant domain below is
D(h)=i[partial_h-beta_h/2]. On that domain

    [D(h),S(f)]=i integral h f,
    [D(h),J(f)]=i J(hf), [D(h),D(l)]=[J(f),J(g)]=0.       (12)

The sign operator acts as epsilon on H_epsilon and commutes with the
declared observable algebra. The covariance, this sign contract and the
chart/polarization are selected quantization data.

Proof. The expected squared H^3 norm is
1+2 sum_(n>=1)(1+n^2)^(-1)<infinity. Thus the Fourier series converges
in H^3 almost surely, and the one-dimensional Sobolev embedding supplies
C^2 representatives. The pushforward of the countable Gaussian product
is a countably additive probability measure on H^3. In particular its
L^2 norm is positive and complete; psi=1 is a normalized nonzero state.
Real measurable multipliers have the maximal self-adjoint domains (9).

Cameron--Martin translation gives (10). One may also prove it from the
finite Gaussian densities: the square multiplier cancels precisely the
density ratio after s -> s+t h. The Gaussian linear functional beta_h
has variance ||h||_C^2. Finite projections converge in that variance,
which yields (10) by L^2 limits and strong continuity first on cylindrical
tests, then on H. The group law follows from
beta_h(s+t h)=beta_h(s)+t||h||_C^2. Stone's theorem proves (11).

A common dense invariant domain A for smooth smearings is the linear
span of products of polynomials in finitely many S(f)'s and finitely
many J(f)'s, with all f smooth. It contains the finite Fourier Hermite
polynomials and is dense. Every such product and all its derivatives
are square integrable: Gaussian linear functionals have all moments,
and Fernique integrability on H^3 bounds exp(b||s||_C2) for every finite
b. For smooth h, beta_h=S(C^{-1}h) also has a smooth smearing.
Directional differentiation sends J(f) to J(hf), and sends a Gaussian
linear functional to a constant. The formulas therefore preserve A.
Differentiating (10) in L^2 on A gives (12). No point product of two
operator-valued distributions is used. The unitary definition, not an
unproved essential self-adjointness assertion for A, specifies D(h).

For later use an explicitly invariant pulse domain in pointer Fourier
variables is the span of exp[-i pi J(v)] a(s)b(pi), where v is real
smooth, a is in A, and b is smooth of compact support. It is dense,
and the pulse simply changes v to v+kappa w. Acting with D(h), X or
Pi adds only J(vh), J(v) or bounded pointer-momentum factors and
derivatives of a,b. The same moment bounds justify these operations.

An additive self-adjoint canonical partner of a positive J(f), f>=0
and not identically zero,
cannot satisfy the full Weyl translation relation: it would conjugate
J(f) to J(f)+tI for all real t, whereas J(f) has spectrum [0,infinity)
in H_+ (and the negative half-line in H_-). Unitary conjugation preserves
spectrum. Our affine relation respects this boundary. This no-go concerns
the full Weyl contract, not a formal commutator on a restricted domain.

The representation supplies S(f), J(f), D(h) and their spectral
measurements. It does **not** assert self-adjoint point operators r(z),
m(z), rho(z), or a quantization of arbitrary nonlinear inverses of (4).
Those would require additional composite-operator definitions. Complete
classical coordinates and a specified quantum observable algebra do not
make such definitions automatic.

## 4. Controlled finite-mode approximation

Let P_N be the Fourier projection through n=N, s_N=P_Ns, and
mu_N its (2N+1)-dimensional marginal. The embeddings
I_N:L^2(mu_N)->H, I_N psi(s)=psi(s_N), are isometric. Their ranges
increase and have dense union. They use the actual marginal measures,
not a cell-product ansatz unrelated to (8). Set

    G_w,N=epsilon integral w exp(s_N), G_w=J(w).          (13)

**MPQ-T4 (continuum and pulse convergence).** G_w,N -> G_w almost
surely and in every finite L^p(mu). As multiplication operators on H,
G_w,N converges in strong resolvent sense to G_w. With the canonical
pointer Pi=-i partial_X and its usual H^1 domain, the joint-spectral
products A_N=G_w,N tensor Pi and A=G_w tensor Pi are self-adjoint and

    exp(-it A_N) -> exp(-it A) strongly,

uniformly for t in compact intervals. The finite-mode dynamics restrict
to I_N L^2(mu_N) tensor L^2(pointer). Thus this is a compatible
approximation of the full continuum dynamics, not just coefficientwise
agreement of a few matrix elements. D(P_N h) -> D(h) in strong resolvent
sense for h in H_C, through their unitary groups.

Proof. Fourier convergence in H^3 implies uniform convergence of s_N,
and |G_w,N|<=||w||_L1 exp(c||s||_H3), with a constant independent of N.
Fernique integrability and dominated convergence give every finite L^p
statement. Bounded resolvents converge pointwise almost surely; dominated
convergence proves their strong limits. In pointer Fourier variables
(s,pi), the propagators are multiplication by exp[-it G_w,N(s) pi].
Their difference is bounded by min(2,|t pi| |G_w,N-G_w|), giving strong
convergence uniformly on compact t intervals for every vector. The same
Cameron--Martin argument used in (10), with P_Nh -> h in H_C, proves
strong convergence of its unitaries and hence the generator resolvents.

The cutoff is in log density. Exponentiating a finite Fourier profile
already produces higher harmonics; this construction is not a claim
that exp(s_N) itself has Fourier support <=N. Infinite covariance modes
remain present in the limit. The Gaussian weight and its derivative term
in D cannot be dropped while preserving these adjoints and norms.

## 5. Clock constraint, physical norm and material pointer instrument

Reuse RPM's ideal clock/pointer construction with the newly constructed
system H_red and G_w, rather than assuming the old system Hilbert space.
Let alpha(T)=integral_(-infinity)^T g(t)dt for real smooth compactly
supported g. On H_0=H_red tensor L^2(R_X) let A=G_w tensor Pi, with
joint-spectral domain integral |G_w pi|^2 |psi_hat(s,pi)|^2 dmu dpi<infinity.
Put W(T)=exp[-i alpha(T)A] and let mathcal W multiply histories by W(T).

**MPQ-T5 (derived positive physical completion).** On
L^2(R_T;H_0), define

    C_T=mathcal W(-i partial_T)mathcal W^*,
    Dom C_T=mathcal W H^1(R_T;H_0).                       (14)

This is self-adjoint. On its spectral test core it is the expected
-i partial_T+g(T)A. For D_A the union of bounded spectral subspaces of A,
use Phi=mathcal W[S(R_T) tensor_alg D_A] and

    Rf=integral W(T)^* f(T)dT.

The group-averaged form and its completion are

    integral dt <f,exp(it C_T)h>=<Rf,Rh>_(H_0),
    H_phys=H_red tensor L^2(R_X).                         (15)

Every v in H_0 supplies a generalized solution Psi_v(T)=W(T)v with
physical inner product <Psi_v,Psi_u>_phys=<v,u>. There are no nonzero
solutions in the kinematic L^2 history space. This derives the norm for
the clock constraint conditional on the chosen reduced representation;
it does not derive mu from the source constraints.

Proof. Unitary conjugation proves self-adjointness and the domain in
(14). Removing W turns the test space into Schwartz vector-valued
histories with bounded spectral vectors. Integrating the translation
kernel of -i partial_T gives the two integrals in <Rf,Rh>; equivalently
Fourier transformation gives 2 pi times the zero-frequency pairing.
R has dense image D_A (choose a Schwartz clock function of integral one).
Its null quotient therefore completes to H_0. A solution constant in the
rotated frame has infinite clock L^2 norm unless it is zero. This is the
same general averaging argument already proved in RPM, now with (8),(13)
as its actual system inputs. No infinite-dimensional gauge Haar measure
has been introduced.

For any self-adjoint system observable O its relational realization on a
clock slice is W(tau)^* O W(tau), on W(tau)^* Dom O. The corresponding
history operator is mathcal W[W(tau)^* O W(tau)]mathcal W^* and strongly
commutes with C_T. This defines domains even where a formal sum below
would not by itself specify a self-adjoint operator.

At pulse end kappa=integral g, U=exp(-i kappa A). On the invariant test
domain of section 3 the exact Heisenberg identities are

    U^* X U=X+kappa G_w,
    U^* D(h) U=D(h)+kappa Pi J(hw),
    U^* J(f) U=J(f).                                     (16)

Their self-adjoint realizations are the unitary conjugates, with their
transported domains. Differentiating the multiplication phase in (s,pi)
proves (16), including the sign of the second term. The reference profile
has quantum backreaction; the interaction is not just a classical readout.

For a centered Gaussian pointer phi_sigma with position variance sigma^2,
the Kraus operators and their normalization are

    K_x=phi_sigma(x-kappa G_w), integral K_x^*K_x dx=I.    (17)

This positive normalized instrument follows from the spectral theorem
and the unitary pointer coupling. For a normalized system state psi the
outcome density is

    p_psi(x)=sum_epsilon integral |psi_epsilon(s)|^2
                  |phi_sigma(x-kappa G_w(s))|^2 dmu(s).  (18)

Where the moments exist, its mean is kappa <G_w> and its variance is
sigma^2+kappa^2 Var(G_w). These statements and (15) reuse RPM's theorem,
not its special phase-packet numerical prediction. The finite-mode clock
operators from (14) converge in strong resolvent sense: their resolvents
are W_N(-i partial_T-z)^(-1)W_N^*, and both W_N and W_N^* converge strongly.
The Gaussian-pointer isometries converge strongly too. Hence their output
states converge in trace norm for trace-class inputs, and all bounded
pointer outcome tests converge. This supplies an operational cutoff limit.

## 6. Information, semiclassical control and the homogeneous comparison

**MPQ-T6 (observable information and controlled classical recovery).**
The probe family G_w reads j, not all the conjugate profile information.
For example, in either sign sector the normalized states

    psi_+(s)=exp(i b s_0), psi_-(s)=exp(-i b s_0), b!=0,

have identical G_w statistics for every material window. But

    <D(1)>_+=-b, <D(1)>_-=b, Var_+(D(1))=Var_-(D(1))=1/4. (19)

Thus phase information exists and a specified self-adjoint conjugate
observable distinguishes it. Equation (19) is a statement about its
spectral measurement; this packet does not add a second source apparatus
coupled to D(1). The frozen G_w pointer alone does not measure that phase.

For a controlled classical limit reinstate hbar>0 and choose covariance
hbar C. Write beta_h and ||h||_C using C, not hbar C. For smooth centers
s_*,d_* the normalized states

    psi_(s_*,d_*)(s)=exp[beta_(s_*)(s)/(2 hbar)
                         -||s_*||_C^2/(4 hbar)]
                    exp[-i integral d_* s/hbar]           (20)

have log-density distribution N(s_*,hbar C). With
D_hbar(h)=i hbar partial_h-i beta_h/2 and
K(z,z')=sum lambda_j e_j(z)e_j(z'),

    <D_hbar(h)>=integral h d_*,
    Var(D_hbar(h))=(hbar/4)||h||_C^2,
    <J(f)>=epsilon integral f exp[s_*+(hbar/2)K(z,z)]dz.  (21)

The J covariance is the double integral against f(z)l(z') of

    exp[s_*(z)+s_*(z')+(hbar/2)(K(z,z)+K(z',z'))]
                         (exp[hbar K(z,z')]-1).          (22)

It tends to zero, and the means tend to the chosen classical profiles.
For all finite sets of these smearings the construction therefore has
the correct peaked classical data, not merely a positive abstract norm.
The inverse of (4) assigns every smooth regular MRF classical profile
such centers. It does not define its nonlinear quantum composites.

Proof. In (19), derivative of the phase contributes respectively -b,+b,
and the Gaussian derivative term has mean zero and variance 1/4 because
lambda_0=1. In (20), the squared real factor is the Gaussian translation
density. Acting with D_hbar leaves the mean integral h d_* and the
imaginary centered Gaussian term -i beta_h(s-s_*)/2, of variance
hbar ||h||_C^2/4. Gaussian exponential moments prove (21),(22).
Continuity and boundedness of K on the circle justify the integrations
and hbar ->0 limits. Under the pulse, (16) with restored hbar is exact;
in particular its means follow j'=0 and d'=Pi w j. Taking an independent
pointer packet with momentum mean pi_* and variance tending to zero
gives the classical profile backreaction and pointer mean (7).
One can let its position and momentum variances both be O(hbar).

**MPQ-T7 (precise failure of the inherited quantum identification).**
Classically the w=1 homogeneous surface is invariant and (7) is the
frozen homogeneous RPM flow. Quantum mechanically the present sign-
preserving theory is not unitarily equivalent to GOA/RPM under a map
that identifies both the classical orientation label and the measured
probe. In a positive sector G_1 has spectrum [0,infinity), whereas
the old full-line F is unitarily equivalent to -i partial_y and has
spectrum R. Even with both new signs, [epsilon,G_1]=0, while the old
F does not commute with sign(P).

Proof. G_1>0 pointwise in H_+, and its essential range is (0,infinity):
the independent Gaussian constant mode rescales it by arbitrary positive
factors. There is no zero eigenvector, although zero belongs to its
spectrum. The negative sector is analogous. GOA's full-line unitary
straightening supplies the old spectrum. Its explicit unit positive-
and negative-P packets have <psi_-,F psi_+>=-1/32, so the sign(P)
commutator is nonzero. A unitary intertwining both operators and their
domains cannot change either fact. This is a conditional non-equivalence,
not a contradiction in GOA or a no-go for every other MRF quantization.

The N=0 cutoff is an exact one-degree sign-preserving affine model:
after removing the Gaussian density it is L^2(ds)=L^2(d|j|/|j|),
J=epsilon exp(s), D=i partial_s. It is a controlled homogeneous
quantization choice, but its one-coordinate image is **not** an
invariant homogeneous sector of the full continuum G_1 dynamics:
exp[-i kappa pi integral exp(s)] depends on nonzero modes. The set of
exactly homogeneous s profiles has mu measure zero. No delta-supported
homogeneous vector has been declared normalizable in (8).

## 7. Closure boundary and the next source obligation

Closed under the explicit representation contract: the canonical source
bridge; exact reduced pulse dynamics; a positive continuum Hilbert space
for the complete separating adapted profile chart; self-adjoint smeared
observables and their common invariant calculation domain; a compatible
finite-mode limit; generalized clock-constraint solutions and their
positive physical inner product; normalized pointer measurement and
backreaction; coherent semiclassical recovery; and the precise failure
of the inherited orientation/probe quantum identification.

There is no source uniqueness theorem for C. For example replacing every
lambda_j by t lambda_j with t>0 gives other constructions. In infinitely
many independent modes, the average of s_j^2/lambda_j converges almost
surely to 1 under mu and to t under the rescaled measure. The two measures
are mutually singular for t!=1. No translation-invariant field Lebesgue
measure or arbitrary unitary equivalence has been presumed.

Still open in the **same frozen model**: quantize the original local
H_tilde,D_tilde on a justified common domain, including inverse p, the
reference modes and the coupling; then compare its physical quotient
with (8),(15), or prove an obstruction for an explicitly stated ordering.
The classical Abelianization (1)-(3) is the concrete starting point, but
nonlinear canonical transformations do not automatically lift to quantum
unitaries, and the restricted source gauge chart is not a vector space.
LQA's vacuum anomaly cancellation is not reused for this extension.
Self-adjoint/composite quantum geometry in the original r,m variables,
patches with zeros, source selection of the measure, a fully local
apparatus, physical 3+1 gravity and calibration are also not closed here.

Standard Gaussian, Cameron--Martin, spectral and affine methods are
credited in SOURCE_MAP.md. The source-specific contribution is their
explicit integration with the unchanged MRF action and probe, the global
chart, cutoff control and comparison obstruction. All infinite-mode,
domain, measure and spectral claims above have written proofs; the finite
native certificates are separate regression witnesses, not a general
proof-assistant certification or experimental evidence.
