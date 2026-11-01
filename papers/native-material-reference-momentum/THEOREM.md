# MRM-1 — Physical material-reference momentum

16 October 2026. Frozen parent: MDD-1 at
112d3d012f7f4189abf4e9e92903c543c3af0817. Keep the actual MRF source,
MPQ chart, selected RGR physical pairing, covariance, and MDD apparatus.
Write h=hbar>0; the unit material circle has Fourier frequency 2 pi n.
The two orientation sectors are treated identically and are not mixed.

This packet addresses the original reference momentum rho, a nonlinear
source-geometric quantity beyond the configuration observable J. Its
total has a self-adjoint physical realization with a controlled cutoff.
The noisy MDD product estimator of that same total has no finite vacuum
continuum law. Nonconstant local smearings face a further canonical-flow
obstruction in this chosen measure class. These are different statements.

## 1. Source current, physical transformation and constraint lift

MRF's representative of the regular physical quotient is

    Z=x, phi=tau, k=0, p=m, q=r, rho=m'-m r'.

MPQ/RMG give j=epsilon exp(s), y=d/j, theta=Y^{-1}(y),
a(theta)=1+cos(theta)/2, m=j/a(theta),
P=integral m, and r=theta+log|m/P|. The physical brackets are
{d(z),s(z')}=delta(z-z') and the canonical potential is
-integral d delta s up to an exact variation. All classical fields and
test functions in this section are smooth and periodic; m has one sign.

**MRM-T1 (source reference current).** The exact current and its smear are

    rho=-j y'=-d'+d s',
    R[w]=integral w rho=integral d(w'+w s'), R=R[1].       (1)

Its Hamiltonian action and algebra are

    {s,R[w]}=-(w'+w s'), {d,R[w]}=-(w d)',
    {R[w],R[v]}=R[v w'-w v'].                            (2)

If g_t is the circle flow of the vector field -w, its finite physical
action is

    s_t=s composed with g_t+log g_t',
    d_t=(d composed with g_t)g_t'.                       (3)

It follows that j,m are densities, y,theta are scalars, and

    r_t=r composed with g_t+log g_t',
    rho_t=(rho composed with g_t)(g_t')^2, P_t=P.         (4)

The reconstructed circumference integral exp(r) is unchanged. These
are transformations of physical target-material profiles, not an
additional quotient or the inverse-embedding source gauge action.
The target rod labels have not been declared redundant.

Proof. Differentiating r gives rho=-m theta'=-j y', since
theta'=a(theta)y'. This proves (1). Functional differentiation gives
delta R/d d=w'+w s' and delta R/d s=-(w d)'; use the stated bracket
and periodic integration by parts to obtain (2). In particular its
smearing sign is the negative of w v'-v w'. Solving the transport
equations gives (3). Substitute in the inverse chart; P is preserved
by change of variables. The g'' terms cancel in m_t'-m_t r_t',
proving (4).

There is a local classical weak-Dirac source lift for every such smooth
physical functional. In MPQ's exact canonical source coordinates
(c_m,m),(phi,K_phi),(Z,R_Z), set r=c_m+log|m| and use (1).
This lift commutes with K_phi and R_Z. Since

    H_tilde=p K_phi+u Z'R_Z, D_tilde=v K_phi+Z'R_Z,

its brackets with the actual source constraints lie in their smooth
constraint ideal. Adjoin a canonical ideal clock and pointer and
C_T=E+g(T)Pi R[w]; the extended constraints are first class locally
while the source remains in the regular chart. This is a classical
lift. The physical quantum realization below specifies an ordering
after reduction, not equality with every unreduced nonlinear ordering.

## 2. A total-charge operator with an exact physical domain

Use H_epsilon=L^2(mu_h), where mu_h=N(0,h C) and
lambda_0=1, lambda_n=(1+n^2)^(-4) on each real sine/cosine pair.
Let S_j be multiplication by s_j and

    D_j=i h partial_(s_j)-i s_j/(2 lambda_j),
    [D_j,S_k]=i h delta_jk.                              (5)

This is MPQ's selected positive physical representation; no new inner
product is inferred from kinematic positivity.

For material translations define

    T(t)psi(s)=psi(s( . +t)), T(t+1)=T(t).                (6)

The paired covariances make T(t) unitary. It is strongly continuous:
first use finite-coordinate continuous tests and then density and
unitarity. Write T(t)=exp(-it Rhat/h).

**MRM-T2 (total physical charge and cutoff).** Rhat is self-adjoint. On
finite Fourier polynomials its action is

    Rhat_N=sum_(n=1)^N 2 pi n
            (S_n^s D_n^c-S_n^c D_n^s)
           =i h sum_(n=1)^N 2 pi n
            (s_n^s partial_(s_n^c)-s_n^c partial_(s_n^s)). (7)

The Gaussian score terms cancel because the two variances are equal.
This is a proved quadratic composite, not substitution of a random
smearing s' into MPQ's operator D(h).

The orthogonal spectral projections and full domain are

    E_k=integral_0^1 exp(2 pi i k t) T(t) dt, k in Z,
    Dom Rhat={psi:sum_k (2 pi h k)^2 ||E_k psi||^2<infinity},
    Rhat psi=sum_k 2 pi h k E_k psi.                    (8)

Its spectrum is exactly 2 pi h Z, in both directions. The vacuum 1 is
in its kernel. No physical energy interpretation or lower bound follows.

Proof. Differentiating (6) on a polynomial gives (7) and Stone's theorem
gives self-adjointness. Fourier analysis of the period-one unitary group
gives (8); completeness of the E_k follows by the Fejer approximate
identity for continuous Hilbert-valued orbits. Every k occurs:
(s_1^c+i s_1^s)^k for k>=0, and the corresponding conjugate for k<0,
are nonzero L^2 eigenvectors. All other spectrum is excluded by (8).

For a graph core, project onto polynomials of bounded Gaussian chaos
degree involving only the first M full sine/cosine pairs and the zero
mode. These finite-dimensional spaces reduce T(t); their increasing
orthogonal projections tend strongly to identity. They commute with
Rhat and therefore also approximate every vector in its graph norm.
Thus the finite Fourier polynomials form a core, not merely a dense
subset of the Hilbert space.

Let T_N rotate only the first N pairs and define Rhat_N by its Stone
generator on all of H. On each finite-coordinate vector T_N(t) and
T(t) agree once N is large enough, for every t. Density gives strong
convergence uniformly in t, hence strong resolvent convergence
Rhat_N -> Rhat. The formula (7) is its action on the polynomial core.

## 3. Domains, source feedback, spectral readout and classical recovery

Let A be MPQ's common algebra of polynomial products of smooth S(f)
and J(f)=epsilon integral f exp(s). It is contained in Dom Rhat and
invariant under Rhat. Translation differentiates the smooth smearing,
while the inherited Gaussian exponential moment bounds justify every
derivative in L^2. In particular,

    [Rhat,S(f)]=-i h S(f'),
    [Rhat,D(g)]=-i h D(g'),
    [Rhat,J(f)]=-i h J(f')                              (9)

on A. The full self-adjoint operator transport, including domains, is

    T(t)^* O(f) T(t)=O(f( . +t)), O=S,D,J.              (10)

For D, conjugate the Cameron--Martin unitary group of MPQ, using
translation invariance of C, and then its Stone generator. This proves
the full-domain statement; a formal commutator alone would not.

**MRM-T3 (feedback and measurement).** The unchanged density pulse
U=exp[-i kappa J(w) tensor Pi/h] has the exact balance

    U^*(Rhat tensor I)U=Rhat tensor I-kappa J(w') tensor Pi. (11)

The right-hand side has the self-adjoint realization with domain
U^*Dom(Rhat tensor I). Equality on a common calculation domain is
justified on the span of exp[-i Pi J(v)/h] a(s)b(Pi), with smooth real
v, a in A, and b smooth of compact support. This domain is dense and
stable under the relevant differentiations and pulses; the same
Gaussian bounds apply. The single commutator in (9) gives (11), with
all higher commutators zero. It agrees with MRF's source balance
d rho/dT=-g Pi j w', including the sign and reference feedback.

The joint-spectral product Rhat tensor Pi is self-adjoint. For the
additional ideal spectral pointer choose real smooth compactly
supported g(T), real pulse area kappa, and take

    W(T)=exp[-i alpha(T) Rhat tensor Pi/h], alpha'=g,
    C_T=W(-i h partial_T)W^*,
    Dom C_T=W H^1(R_T;H_red tensor L^2(R_X)).             (12)

On its differentiable core C_T=E+g(T)Rhat tensor Pi.
The standard free-clock average transported by W gives the physical
norm ||psi_0||^2 on H_red tensor L^2(R_X); histories are W(T)psi_0.
Thus this norm is derived for the declared apparatus, as in MPQ/RPM.

A centered Gaussian pointer with position variance sigma^2 gives
Kraus density M_x=phi_sigma(x-kappa Rhat). Its effects integrate to I,
and for normalized psi in Dom Rhat,

    E X=kappa <Rhat>, Var X=sigma^2+kappa^2 Var Rhat.     (13)

The corresponding cutoff pulses converge strongly: use the T_N
convergence in each pointer-momentum fiber, then truncate that
momentum and use dominated convergence. This is an added ideal
spectral measurement contract. MDD's joint-quadrature measurement is
a different apparatus and need not give this law.

For smooth classical centers s_*,d_*, use the normalized MPQ coherent
state

    psi_*(s)=exp[beta_(s_*)(s)/(2h)-||s_*||_C^2/(4h)]
              exp[-i integral d_* s/h].                (14)

It belongs to Dom Rhat, the cutoffs converge on it in graph action,
and its charge mean and variance are

    <Rhat>=integral d_* s_*',
    Var Rhat=h [C(d_*',d_*')+(1/4)||s_*'||_C^2].         (15)

Here C(f,f)=<f,Cf> and ||f||_C^2=<f,C^{-1}f>. To prove (15),
differentiate the translation orbit of (14) in L^2:

    Rhat psi_*=[integral d_* s'
                      +(i/2)<C^{-1}s_*,s'>] psi_*.

All coefficients are smooth. Under |psi_*|^2 mu_h, s has law
N(s_*,hC), and <C^{-1}s_*,s_*'>=0 because C commutes with the
skew derivative. Gaussian variances of the real and imaginary linear
parts give (15); their cross term in the squared modulus is zero.
Fourier tails of these smooth coefficients give the stated graph-action
limit. Chebyshev then proves spectral concentration at the classical
charge as h -> 0. An actual pointer also needs sigma -> 0 to remove
its preparation noise.

For example s_*=a e_1^c and d_*=b e_1^s give mean -2 pi a b.
The pair d_*=+b e_1^s and d_*=-b e_1^s has the same complete J law
but opposite charge means, separated by 4 pi |a b| when ab!=0.
This is physical conjugate-phase information within the selected
model, not a calibrated gravitational experiment.

## 4. Why the raw joint-detector estimate fails

MDD/CMG's balanced vacuum labels are independent real Gaussians,

    u_j ~ N(0,2h lambda_j), v_j ~ N(0,h/(2 lambda_j)).    (16)

Their smooth finite reconstruction of the total current is

    Z_N=integral v_N u_N'
       =sum_(n=1)^N 2 pi n (u_n^s v_n^c-u_n^c v_n^s).  (17)

**MRM-T4 (operator versus noisy product).** At each finite N the first
moment operator of (17), on the common polynomial domain, is Rhat_N.
Nevertheless its vacuum laws are not tight as N -> infinity, even
after any deterministic scalar subtraction. The spectral Rhat law
of that same vacuum is delta_0.

Proof. The joint dilation reads u_j=S_j+A_j, v_j=D_j-B_j with centered
independent matched probes. Every mixed term in (17) involves distinct
sine/cosine indices, so probe compression leaves S_n^s D_n^c and
S_n^c D_n^s with no same-index ordering correction. This proves the
finite first-moment identity, not an identity of measurements.

Set X_n=u_n^s v_n^c-u_n^c v_n^s. Conditional on the two u's, X_n is
Gaussian with variance h[(u_n^s)^2+(u_n^c)^2]/(2 lambda_n).
Gaussian integration gives

    E exp(it X_n)=(1+h^2 t^2)^(-1).

Thus the independent X_n have centered Laplace density
exp(-|x|/h)/(2h). The last summand 2 pi N X_N has density bounded by
1/(4 pi h N), and convolution cannot increase that bound. For every
L>0 and every real a_N,

    Pr(|Z_N-a_N|<=L)<=L/(2 pi h N) -> 0.                (18)

In particular Var Z_N=2h^2 sum_(n=1)^N(2 pi n)^2 diverges.
Variance divergence alone would not prove non-tightness; (18) does.
This excludes the raw full-record product prescription, not every
renormalized estimator or every measurement of the valid operator
Rhat. It also does not remove CMG's scale and threshold obstructions.

## 5. A local Weyl current fails on the inherited vacuum domain

For a real smooth window w let P_N include modes -N,...,N. In a real
orthonormal basis set A_ij=<e_i,w e_j'> and b_i=<e_i,w'>. The natural
symmetric cutoff of (1) is

    Rhat_w,N=D(P_N w')
                 +(1/2)sum_ij A_ij (D_i S_j+S_j D_i).  (19)

Each finite operator has a self-adjoint realization from the Gaussian
half-density action of the complete affine flow with velocity b+A s.
Conjugating by the square root of the finite Gaussian density gives
i h[(b+A s).partial+Tr A/2] on Lebesgue L^2. Compact smooth tests are
an invariant core for this affine unitary group; Gaussian polynomials,
in particular the vacuum in Gaussian representation, lie in its
generator domain. This finite statement is not its continuum limit.

**MRM-T5 (vacuum-domain obstruction).** If w is nonconstant, the vectors
Rhat_w,N 1 have unbounded norms, with a second-chaos component growing
at least as a positive constant times sqrt(N). Scalar counterterms
cannot make them converge. Additive linear S/D counterterms cannot
cancel that second-chaos component either.

Proof. Put B=C_N^{-1/2} A C_N^{1/2}, S=(B+B^*)/2,
b_C=C_N^{-1/2}P_N w', and x=C_N^{-1/2}s/sqrt(h). Direct application
of (19) gives

    Rhat_w,N 1=-(i sqrt(h)/2)b_C.x
                   -(i h/2)(x^T S x-Tr S),
    ||Rhat_w,N 1||^2=(h/4)||b_C||^2+(h^2/2)||S||_HS^2. (20)

The two terms are orthogonal Gaussian chaos degrees one and two.
These formulas may be evaluated in the complex Fourier basis without
changing the real operator's Hilbert--Schmidt norm. Write
sigma_n=sqrt(lambda_n)=(1+n^2)^(-2), including n=0. Then

    S_mn=pi i w_(m-n)
           [n sigma_n/sigma_m-m sigma_m/sigma_n].        (21)

Choose k!=0 with w_k!=0 and set m=n+k. As n -> +infinity,

    n sigma_n/sigma_(n+k)
       -(n+k) sigma_(n+k)/sigma_n -> 7k.                (22)

Indeed sigma_n/sigma_(n+k)=1+4k/n+O(n^-2) and its inverse has
the opposite first-order coefficient. There are order N entries of
(21) bounded away from zero along this diagonal, so ||S||_HS^2>=c_w N
for all sufficiently large N. The smoothness of w makes the linear
norm in (20) bounded, but cannot cancel an orthogonal chaos component.

For an exact illustration w=sin(2 pi z), a_n=1+n^2,
b_n=1+(n+1)^2, the positive adjacent entries satisfy, for n>=1,

    S_(n+1,n)/(2 pi)
       =(1/4)[n b_n^2/a_n^2-(n+1)a_n^2/b_n^2]>=1.       (23)

Clear denominators and put x=n-1. The numerator of the bracket minus
4 becomes

    n b_n^4-(n+1)a_n^4-4 a_n^2 b_n^2
       =193+1041x+2084x^2+2292x^3+1570x^4
                    +702x^5+204x^6+36x^7+3x^8>0.       (24)

The transpose entries give ||S||_HS^2>=(2 pi)^2 2(N-1), and hence
||Rhat_w,N 1||^2>=(2 pi)^2 h^2(N-1). This last bound is only an
illustration; (22) proves the assertion for every nonconstant smooth w.

Subtracting a scalar changes only chaos degree zero. A linear S/D
term on the vacuum has degree at most one. Neither removes degree
two. For constant w the matrix S and b_C vanish, consistently with
the total-charge vacuum kernel in section 2.

Failure on this domain alone does NOT exclude all self-adjoint
operators, strong resolvent limits on other domains, or other
renormalizations. The next theorem has a stronger but different
hypothesis: implementing the exact classical finite configuration flow.

## 6. No canonical nonconstant local flow in this measure class

For a smooth orientation-preserving circle diffeomorphism g define

    F_g(s)=s composed with g+log g'.                     (25)

**MRM-T6 (physical configuration implementability obstruction).** For
every non-rigid g, (F_g)_*mu_h and mu_h are mutually singular. Thus no
unitary on this configuration multiplication representation, including
any multiplicity/spectator extension or equivalent strictly positive
density reweighting, implements (25) on all its spectral projections.
In particular no self-adjoint local R[w] for nonconstant w can generate
the exact finite canonical action (3) under this contract.

Proof. Reuse the explicit path-statistic proof of frozen IGQ-T2,
now on the physical log density, not on the forward gauge variable.
The extra zero mode and the spatial mean projection in IGQ are
annihilated by the third derivative. With the factor h restored,

    B_s=(2 pi)^(-4) partial_z^3 s
       has continuous quadratic-variation measure h dz.

The Fourier amplitude ratio to sqrt(h) times a periodic Brownian
bridge is n^4/(1+n^2)^2; its defect is O(n^-2), so the difference
is almost surely C^1. The full details, including the measurable
dyadic path statistic, are in IGQ-T2 and are unchanged. The chain
rule gives

    B_(F_g s)=(g')^3 B_s composed with g + a C^1 term,
    d[B_(F_g s)]=h(g')^7 dz.                            (26)

Since integral g'=1, strict convexity gives integral(g')^7>1 unless
g'=1. The total dyadic variations h and h integral(g')^7 distinguish
disjoint full-measure sets. For g'=1+cos(2 pi z)/2 the latter factor
is exactly 4587/1024, as already proved in IGQ.

Unitary conjugation preserves which configuration spectral
projections are zero. Singularity contradicts that requirement, also
with arbitrary multiplicity; an equivalent positive density preserves
the same null sets. For nonconstant w, g_t'=1-t w'+O(t^2) is
nonconstant for sufficiently small nonzero t, so a self-adjoint
generator implementing all these flows is impossible. Rigid
translations are precisely the exception constructed in section 2.

This does not obstruct the earlier residual source gauge action: RGR
acts on inverse-embedding gauge variables in a different measure
class, and the physical s factor is a spectator there. No new
R[w]=0 constraint, central anomaly or gauge-rigging change is inferred.
Other physical representations, and local composites not implementing
(25), remain open.

## 7. Scope, native replay and the next obligation

MRM adds a source-derived, phase-sensitive, self-adjoint total charge,
its domain, spectrum, clock/pointer realization, exact density-pulse
balance and classical recovery. It distinguishes measurement noise
from operator existence, and locates a nonconstant local-flow
obstruction using a frozen general proof rather than duplicate work.

The native packet checks Weyl identities, transport cancellation,
exact finite Gaussian moments, covariance rotations, coherent-state
coefficients, the all-n polynomial identity (24), local finite-matrix
norms and the path-statistic coefficient. Stone/domain, measure
singularity, tightness and continuum limit arguments are written
proofs; finite checks do not replace them.

The next physical-geometry decision is now constrained: retaining
this Gaussian configuration measure cannot supply self-adjoint local
reference currents with the full canonical finite-flow contract.
A new representation/observable contract must be justified, or a
specific renormalized geometric observable must be defined without
claiming that flow. Neither an arbitrary scalar subtraction nor a
formal commutator settles it. Full local geometry, physical dynamics,
matter/potentials and 3+1 gravity remain open.

Separately, Fable's corrected ergodicity packet at a9f44d5 has been
read and replayed. Its next task is the source-test-topology bridge
between infinitesimal constraint solutions and finite-group
invariant duals, as specified in FABLE_NEXT_PROMPT.md. That unmerged
work is not a premise of any theorem here.
