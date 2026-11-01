# RGR-1 — A positive residual gauge rigging by covariance doubling

16 October 2026. Frozen parent: CGD-1 at
fea9ef35e1357e8fa1afee24104472ccc2846bce. The MRF source, IGQ canonical
moment maps, inverse-embedding kinematic representation and selected
strongly invariant probe are unchanged. This packet constructs a
nonzero full gauge-invariant test functional and a selected positive
physical quotient. It does not classify every invariant distribution or
prove uniqueness of physical quantization.

The essential calculation is a cancellation between the angular
Jacobian and a conditional Fredholm determinant in the shape variables.
It gives an explicit continuous version of the existing Radon--Nikodym
cocycle. Its square root is the ordinary density cocycle for a second,
mutually singular Gaussian measure. Integrating continuous tests against
that second measure supplies the physical functional. No Haar measure on
an infinite-dimensional group, divergent Gaussian energy, or arbitrary
extension of an almost-everywhere density is used.

## 1. Fixed source, covariance family and notation

After CGD's proved central integer average, one chiral configuration is

    X=S1_alpha x E, E=H^3_0(S1),
    chi(x)=alpha+I_sigma(x),
    k_sigma=exp(sigma)/integral exp(sigma), I_sigma(x)=integral_0^x k_sigma.

The circle has length one. Let mu_c be the centered real Fourier Gaussian
on E with variance c lambda_n in each sine/cosine coordinate,
lambda_n=(1+n^2)^(-4), n>=1, c>0. Put nu_c=d alpha dmu_c, with normalized
angular Haar measure. The *kinematic* measure is still nu_1. The
Cameron--Martin space is H=H^4_0 with

    ||h||_H^2=sum_(n,real modes) (1+n^2)^4 h_n^2
             =sum_(j=0)^4 b_j integral (partial_x^j h)^2,
    b_j=binomial(4,j)(2 pi)^(-2j).                         (1)

For a smooth orientation-preserving circle diffeomorphism g, its action is

    T_g(alpha,sigma)=(g(alpha),sigma+u_g),
    u_g=P_0 log g'(chi), P_0 v=v-integral v.               (2)

Circle lifts can be used in formulas; every expression is lift-independent.
The two chiral factors use X x X and nu_c tensor nu_c. Symbols below
initially refer to one factor. IGQ's Ramer hypotheses hold for every
c>0: a scalar covariance change preserves H, Hilbert--Schmidt derivatives
and invertibility. Quasi-invariance is therefore already available.

## 2. Conditional shape trace and the actual Jacobian cancellation

Let P_N be the complete sine/cosine Fourier cutoff 1<=n<=N. The trace
used below is specifically

    tr_F A=lim_N Tr_H(P_N A P_N).                         (3)

This is not an assertion that A is trace class or that every ordering of
its diagonal is allowed. Define det_F(I+A)=det_2(I+A) exp(tr_F A)
only for the operators proved to have this trace. The ordinary
Carleman--Fredholm determinant det_2 is the one in Ramer's formula.

**RGR-T1 (shape determinant and angular cancellation).** For
A_g=D_sigma u_g on H,

    tr_F A_g=integral_0^1 log g'(y)dy-log g'(alpha),
    det_F(I+A_g)=1/g'(alpha),
    det_2(I+A_g)=exp[-integral_0^1 log g'(y)dy].           (4)

Thus g'(alpha) det_F(I+A_g)=1. The determinant is positive.

Proof. For h in H, set

    B_sigma h(x)=integral_0^x k(t)h(t)dt
                         -I(x) integral_0^1 k(t)h(t)dt.

Then A_g=P_0 a B_sigma, a=(g''/g')(chi). IGQ proves this operator
Hilbert--Schmidt on H. Its Fourier diagonal is unchanged by the H
normalization of each Fourier vector. The kernel of B_sigma is
k(t)[1_(t<x)-I(x)], and B_sigma 1=0. Symmetric Fourier summation of its
jump takes the average of the two one-sided values. Consequently

    tr_F(P_0 a B_sigma)=integral a(x)k(x)[1/2-I(x)]dx.     (5)

Here is a direct check that does not assume a trace-class kernel. Let R_0
be the periodic mean-zero primitive, let M_k mean multiplication by k,
and put j(x)=x-I(x), a periodic H^4 function. On the full Fourier space,

    B_sigma=R_0 M_k-1 tensor (ev_0 R_0 M_k)+j tensor <k,.>.

The conditional trace of M_a R_0 M_k is zero. In complex Fourier
coordinates its partial trace is

    sum_l a_(-l) k_l sum_(|n|<=N,n+l!=0) 1/[2 pi i(n+l)].

For each l the inner symmetric harmonic sum tends to zero, and is
bounded uniformly in N by C[1+log(1+|l|)]. The H^3 coefficient bounds
make the outer sum absolutely dominated, also for parameter derivatives.
The remaining two terms are finite rank and have total trace
-R_0(ka)(0)+integral ka(x-I). The Fourier identity
-R_0(q)(0)=integral q(1/2-x) gives (5). Finally the constant mode has
zero diagonal because B_sigma 1=0; applying P_0 therefore leaves the
same trace on H_0. This argument also proves the locally uniform
convergence and parameter differentiation used below. Substituting
y=I(x) into (5), and integrating by parts, gives the first identity (4).

For clarity the finite determinant is not inferred just by exponentiating
that trace. Join g to the identity by a smooth path of increasing lifts
g_t, for example the linear interpolation of the lifts. Set
f_t=dot g_t composed with g_t^-1. The shape velocity is
v_f=P_0 f'(chi), and its derivative K_f=P_0 f''(chi)B_sigma has

    tr_F K_f=-f'(alpha).                                (6)

Indeed (5) becomes integral_0^1 f''(alpha+y)(1/2-y)dy.
Writing B_t=I+A_(g_t), the variational equation at fixed initial alpha
is dot B_t=K_(f_t)(T_(g_t)x) B_t. All increments and their derivatives
are Hilbert--Schmidt. The differentiation formula for det_2, together
with (3), gives

    d/dt log det_F B_t
      =Tr[(B_t^-1-I) dot B_t]+tr_F dot B_t
      =tr_F(dot B_t B_t^-1)=-f_t'(g_t(alpha)).           (7)

Products of two Hilbert--Schmidt operators in this calculation are
trace class, so their trace is absolutely convergent and cyclic.
The conditional traces and their derivatives exist by the kernel
calculation above; hence the cancellation in (7) is legitimate.
Integrating (7) proves det_F B_t=1/g_t'(alpha). Combining this with
the first identity (4) proves the last. Invertibility and the path
from the identity fix the positive sign.

## 3. A finite energy difference on rough configurations

Typical sigma has infinite H norm. Expressions such as
exp(||sigma||_H^2/4) are therefore not definitions of states.
We instead extend the *difference* of two smooth energies.

For h(y)=log g'(y), u=P_0 h(chi), w=sigma''', set a=u'=h'(chi)k.
Differentiating k'=k sigma' gives

    u''''=a w+b,
    b=h''''(chi)k^4+6h'''(chi)k^3 sigma'
       +h''(chi)k^2[4 sigma''+7(sigma')^2]
       +h'(chi)k[3 sigma' sigma''+(sigma')^3].          (8)

Here a has a bounded derivative and b is H^1 when sigma is H^3.
Define

    C_g(sigma)=sum_(j=0)^3 b_j integral sigma^(j) u^(j)
             +b_4[-(1/2)integral a' w^2-integral w b'],
    Delta_g(alpha,sigma)=2 C_g(sigma)+||u_g||_H^2.       (9)

The letters b_j in (1) are fixed coefficients; the function b in (8)
is a different, explicitly displayed remainder.

**RGR-T2 (continuous energy cocycle and Gaussian density).** Delta_g is
continuous on X, smooth for smooth g, and bounded and locally Lipschitz
on each H^3 ball. It agrees on smooth profiles with

    ||sigma+u_g||_H^2-||sigma||_H^2.

It obeys Delta_(gh)(x)=Delta_g(T_h x)+Delta_h(x). The Radon--Nikodym
derivative for nu_c has the following strictly positive continuous version:

    R_g^(c)(T_g x)=exp[Delta_g(x)/(2c)],
    R_g^(c)(x)=exp[-Delta_(g^-1)(x)/(2c)].               (10)

In particular

    sqrt(R_g^(1)(x))=R_g^(2)(x)                         (11)

at every point of X, not merely on an unspecified full-measure set.

Proof of extension and the crucial Fourier limit. Integration by parts
on smooth profiles gives

    integral sigma'''' u''''
       =-(1/2)integral a' w^2-integral w b'.             (12)

Every term on the right is continuous in H^3. Product estimates in one
dimension prove the asserted bounds and smoothness. Smooth profiles are
dense in H^3, and (2) is continuous there. The smooth energy-difference
identity therefore proves the cocycle identity everywhere.

To use Ramer's formula we must also identify the stochastic contraction;
continuity of (9) alone would not do this. Fix 0<epsilon<1/2. Under every
mu_c, almost surely

    M_epsilon=sup_(n,real modes) n^(4-epsilon)|sigma_n|<infinity,          (13)

and M_epsilon has moments of every order. This follows directly from
Gaussian tails and summing exp(-C n^(2epsilon)). Thus the complex Fourier
coefficients of w are bounded by C M_epsilon |n|^(-1+epsilon).
With P_N as above, the only delicate term in the truncated contraction is

    integral (P_N w)' a w
      =-(1/2)integral a'(P_N w)^2
                       +integral (P_N w)' a (I-P_N)w. (14)

The last term is bounded by

    C M_epsilon^2 N^(-1+2epsilon) sum_l |l| |a_l|.      (15)

For each frequency difference l, at most 2|l| pairs cross the cutoff;
for |n|<=N<|m|, |n||w_n||w_m|<=C M_epsilon^2
N^(-1+2epsilon). The weighted Fourier sum of a is finite since a is
H^3. This proves (15). The b term integrates by parts and converges
using P_N w -> w in L2. Lower derivative terms converge directly.
Hence

    lim_N sum_(n<=N,real modes) sigma_n (u_g)_n/lambda_n=C_g.             (16)

The convergence also holds in L2(mu_c): all coefficients and their
derivatives used here are bounded by a polynomial in ||sigma||_H3
times exp(C||sigma||_H3); (13), Gaussian exponential integrability and
the same estimates provide the necessary domination.

For fixed alpha, u_g is in the Gaussian first Sobolev space with values
in H: its H norm and the Hilbert--Schmidt norm of its H derivative
have the same integrable bounds. The divergence identity consequently is

    delta_c u_g=C_g/c-tr_F A_g.                         (17)

Indeed divergences of P_N u_g converge in L2, and their two finite sums
are exactly (16)/c and (3). This checks, rather than assumes, the
conditional contraction in the change-of-variables theorem.

Ramer's formula for the shape pullback density is

    det_2(I+A_g) exp[-delta_c u_g-||u_g||_H^2/(2c)].    (18)

The angular map contributes g'(alpha). Substituting (17), (9) and
RGR-T1 into (18) leaves precisely exp[-Delta_g/(2c)] for the full
pullback density. Its reciprocal at the image point is (10).
This is the already existing nonsingular representation, now with an
explicit density; no limit of finite probability densities is substituted
for the infinite-dimensional change-of-variables theorem. Both nu_1
and nu_2 have full H^3 support, so continuous versions are unique.
Equation (11) is thus fixed by the source action and the derived cocycle,
not chosen freely on a nu_1-null set.

## 4. A genuine stronger test space and a nonzero physical functional

Use two chiral copies, X_c=X x X, nu_c^(2)=nu_c tensor nu_c. In this
section R_g denotes the product cocycle. Let UC_b(X_c) denote continuous
complex functions uniformly continuous and bounded on each set
with both shape H^3 norms at most R. It does not mean globally bounded.
Define

    Psi=UC_b(X_c) intersect L2(nu_1^(2)) intersect L1(nu_2^(2)),
    p_R(psi)=sup_(||sigma_+||,||sigma_-||<=R)|psi|,
    topology: p_R (R=1,2,...), ||.||_L2(nu1), ||.||_L1(nu2).             (19)

**RGR-T3 (continuous invariant test functional).** Psi is a complete
Frechet space continuously and injectively embedded, with dense image,
in the unchanged H_c=L2(nu_1^(2)). Each residual gauge transformation
acts continuously on Psi by the same kinematic unitary

    U_g psi=R_g^(1)^(1/2) psi composed with T_(g^-1).

The continuous linear functional

    ell(psi)=integral psi dnu_2^(2)                     (20)

is nonzero, ell(1)=1, and ell(U_g psi)=ell(psi) for every residual g.

Proof. Uniform convergence on bounded sets preserves continuity and
uniform continuity on such sets. A Cauchy sequence for (19) has this
pointwise limit as well as L2(nu_1) and L1(nu_2) limits; subsequences
identify the limits almost everywhere. This proves completeness.
Full support makes a continuous function determined by its nu_1 class,
so the embedding is injective. Bounded smooth cylinder functions are
in Psi and dense in H_c.

The maps T_g and T_g^-1 take H^3 balls to bounded sets and are locally
Lipschitz with uniform bounds on those balls. The same holds for the
continuous density factor (10). Thus
p_R(U_g psi)<=C_(g,R) p_(S(g,R))(psi). Unitarity gives its L2 bound.
Equation (11) and ordinary change of variables give the exact L1 isometry

    integral |U_g psi| dnu_2^(2)=integral |psi| dnu_2^(2).                (21)

They also prove invariance of (20), which is continuous by the L1
seminorm. Strong continuity on Psi follows on bounded sets from uniform
continuity, and on L2 and L1 from the strongly continuous nonsingular
representations. For L1, first use bounded continuous tests and
R_g^(2)->1 pointwise with integral one (hence L1 convergence), then
density and (21). All these claims are for actual pointwise continuous
representatives; an arbitrary L2 class cannot simply be integrated
against the singular measure nu_2.

The distinction from CGD-T2 is explicit. If xi_j are the real shape
coefficients divided by sqrt(lambda_j), then N^-1 sum_(j<=N) xi_j^2
converges almost surely to c under mu_c. Thus mu_1 and mu_2 are mutually
singular. A smooth bounded bump equal to one near 2 and zero near 1,
applied to that empirical variance, gives cylinder tests psi_N with

    ||psi_N||_L2(nu1)->0, ell(psi_N)->1.                 (22)

The functional is not Hilbert continuous and is not a kinematic vector.
There is no contradiction with the zero-dual theorem on CGD's weaker
topology. The kinematic covariance has not been doubled: only the
explicit physical test functional uses the second measure.

## 5. Differentiated constraints on an explicit common test core

For one real smooth smearing f, the configuration vector field is

    X_f alpha=f(alpha), X_f sigma=v_f=P_0 f'(chi).

Define A_f by the expression C_g in (9), with u_g replaced by v_f
and h replaced by f'. Equivalently it is the continuous extension of
<sigma,v_f>_H. Then

    d/dt Delta_(g_f(t))|_0=2 A_f,
    B_f=-X_f+A_f/2,
    Q(f)=i B_f, U_(g_f(t))=exp[-it Q(f)].               (23)

Let A be the algebra of smooth functions on X_c whose every Frechet
derivative on H^3 charts is bounded by
C_r(1+||sigma||)^N_r exp(C_r||sigma||), uniformly in the angles.
It contains constants and bounded smooth cylinders with bounded
derivatives. The fields X_f and scores A_f have these bounds;
product/chain rules show B_f A is contained in A. Every element of A
is in Psi and in all finite Lp(nu_c^(2)), by Gaussian integrability.
Set

    D=span{U_g a : a in A, g in G_c}.                  (24)

**RGR-T4 (full constraint solutions and common domains).** D is a
dense gauge-invariant subspace of Psi and H_c, invariant under every
smoothly smeared source generator and every finite mixed product.
The real smeared generators are essentially self-adjoint on D, with
the same Stone closures as IGQ/CGD. On D,

    ell(Q_epsilon(f) psi)=0,
    ell(T_epsilon[f] psi)=ell(H[f] psi)=ell(D_source[f] psi)=0.           (25)

The inherited chiral and H/D algebra closes there with the same zero
central term. D_source denotes the momentum constraint, not test space D.

Proof. Pointwise differentiation of (10) gives (23). For a in A,
the fundamental theorem of calculus along the flow yields

    U_t a-a=integral_0^t U_s B_f a ds.                  (26)

The right side is a Bochner integral both in L2(nu_1) and L1(nu_2),
using strong continuity and the respective isometries. Fubini and the
pointwise formula identify these integrals. On each bounded H^3 set
the same identity differentiates uniformly. Hence a is in the Stone
domain, with (23), and the derivative exists in Psi. Iteration gives
all mixed products because B_f preserves A. Conjugating a flow by g
is another smooth flow, so (24) is also invariant under all generators.
The derivative estimates are continuous in the smooth smearings;
the same argument gives the smooth orbit derivatives needed for the
finite-group-law commutators. Alternatively differentiating the
continuous cocycle (10) and its pullback gives them directly on D.
CGD's source signs T_epsilon=-epsilon Q_epsilon are unchanged.

For essential self-adjointness, restrict any real Stone generator
(including a real H or D_source smearing) to D. It is dense, symmetric
and invariant under its entire one-parameter unitary group. A deficiency
vector would have bounded unitary matrix coefficients satisfying an
exponentially growing scalar ODE in one time direction. The coefficients
must vanish, and density kills the deficiency vector. This is the same
core argument used in CGD, now on the explicit rigging-compatible space.
Differentiating ell(U_t psi)=ell(psi) in the Psi topology proves (25).

Covariance doubling is selected within the Gaussian scale family, not
an arbitrary lucky normalization. Integration by parts from (17),(6)
gives, for a in A and any c>0,

    integral X_f a dnu_c=(1/c) integral A_f a dnu_c,
    integral Q(f)a dnu_c=i(1/2-1/c) integral A_f a dnu_c.                 (27)

For nonconstant f, A_f is not the zero function: its first shape
derivative at sigma=0 pairs h with the nonzero profile f'(alpha+x).
Taking a=A_f, full support makes integral A_f^2 positive. Thus c=2 is
the only member of this Gaussian covariance family whose integration
functional solves every unchanged constraint. This is not uniqueness
among all distributions or all admissible pairings.

## 6. Positive pairing and actual reduced-observable recovery

**RGR-T5 (a selected positive full gauge sector).** On Psi define

    eta(psi)[phi]=overline(ell(psi)) ell(phi).            (28)

It is separately continuous, Hermitian, positive, and invariant in
each argument under the whole residual group. Its null space is
ker ell and its quotient is exactly C, with physical norm |ell(psi)|.
The same statements hold on D, since 1 belongs to D and ell(1)=1.

This channel count is derived for the explicit pairing (28). We do
not infer that the entire invariant dual is one-dimensional or that
all other positive forms are excluded.

There is also a full finite-source-gauge functional before central
reduction. In CGD notation take

    Phi_full={f in Phi_Z : P_Z Jf in Psi}

with the graph topology of Phi_Z and Psi, and set
ell_full(f)=ell(P_Z Jf). This is a complete dense test space: P_Z is
Hilbert continuous and Psi is complete; finite sequences with entries
in Psi are dense in the kinematic sequence Hilbert space. Its full
source-group action is continuous by CGD's bounded integer-shift
estimate and its exact intertwiner. Consequently ell_full is nonzero
and invariant under every finite element of G, including the central
subgroup. This does not assert a new common unreduced infinitesimal
domain across the integer-section seams; (25) is on the proved central
quotient, as in CGD.

Tensor with the unchanged MPQ physical representation. On
Psi tensor_alg H_red define

    R(sum_j psi_j tensor v_j)=sum_j ell(psi_j)v_j,
    eta_red(F)[G]=<R F,R G>_Hred.                       (29)

R is onto because R(1 tensor v)=v. Its positive null quotient therefore
completes to H_red by an explicitly defined unitary, not an abstract
Hilbert-space dimension argument. Generalized physical states are
F -> <v,R F> for v in H_red. Their norm is exactly ||v||_Hred.

For every selected MPQ self-adjoint observable O and its spectral
measure E_O, R(I tensor E_O(B))=E_O(B)R. Thus the descended operator
has precisely MPQ's full maximal spectral domain Dom O; its action,
adjoints and all bounded spectral functions are intertwined. On
D tensor_alg D_red the source generators annihilate the rigging and
commute with the physical calculation algebra. In particular the
declared source probe G_w^* becomes J(w), and the clock/pointer unitaries,
domains and normalized instrument of MPQ pass through the same map.
This discharges IGQ-T6's existence condition for this *selected* rigging
and observable contract. It is not a derivation for every operator in
the full kinematic commutant or an alternative original-field ordering.

## 7. Scope and next obligation

Established here: one nonzero generalized solution functional for the
full finite gauge group; an explicit stronger complete test topology;
a compatible common constraint core; a positive selected pairing; and
an actual norm-and-spectral-domain intertwiner with MPQ. The physical
pairing is derived from (20),(28), never from kinematic positivity.

Still open: classification of the full invariant dual and admissible
positive forms, uniqueness under an appropriately complete observable
contract, nonlinear quantum geometry beyond MPQ's adapted spectral
observables, original-field ordering comparisons, critical-current
patches, fully local apparatus, matter/potentials and physical 3+1
identification. No equivalence to the older orientation-mixing HQS/GOA
quantization is asserted; MPQ's existing obstruction remains intact.

The Gaussian change-of-variables theorem is an external analytic input.
The determinant cancellation, rough-profile energy extension and rigging
construction are written proofs tied to the fixed source. Native exact
fixtures replay their algebraic identities and reject specified wrong
formulas; they do not machine-prove the continuum analysis. Independent
review of the conditional trace and contraction estimates is especially
valuable. All predecessor packets remain frozen.
