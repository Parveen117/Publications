# Geometry refinement: a compatibility obstruction and controlled source dynamics

Author: Monty. 30 October 2026.
Frozen parent: `eb7916e6e4f4bf30db4f98f3792b46d60b4ba636`.

The two-window predecessor left two precise questions: whether increasing its
resolved windows defines consistent quantum measurements, and whether the
source geometry and its pulse evolution survive refinement. They have different
answers. The complete measurement records do not form exact refinements at
fixed h: a previously sharp residual configuration cannot remain sharp while
its noncommuting conjugate is added. Nevertheless, every finite Fourier family
has a normalized measurement, an exact bounded-length evolution for pulses
in that family, and a controlled classical comparison. Taking the classical
limit first and then increasing the family recovers arbitrary smooth source
profiles and their prescribed density-pulse evolution.

All claims here are written analytic arguments conditional on the frozen LMD
inputs and the declared measurement rules. No scientific numerical calculation,
certificate replay, new NRQC run, hardware execution or CI dispatch accompanies
this packet. The dynamics below is the existing source probe interaction;
no autonomous gravitational Hamiltonian is derived.

## 1. Finite smooth-mode construction and the precise refinement question

Keep the positive LMD space H_h=L2(mu_h), 0<h<1. Let A be any finite collection
of distinct real Fourier modes including f0=1, with the nonconstant modes
cos(2 pi n z) or sin(2 pi n z), without sqrt(2) normalization. Write

    s=sum_(a in A) q_a f_a+eta_A,
    Var(q_a)=h alpha_a,  alpha_0=1, alpha_(n,c/s)=1/n.   (1)

The q_a are independent of one another and eta_A. Put m=|A|. The complete
removed Wick variance and the residual base measure are

    V_A(z)=h sum_(a in A) alpha_a f_a(z)^2,
    B_eta(dz)=exp[-V_A(z)/2] M_perp,A(eta,dz),
    M(q,eta,dz)=exp[sum q_a f_a(z)] B_eta(dz).           (2)

B_eta is finite, positive, nonatomic and full support almost surely, as in
the predecessor's finite-mode factorization. All unresolved continuum modes
remain in the state space. Define

    F_eta(q)=M(1),     x_a=partial_(q_a)F_eta=M(f_a),
    H_ab=partial_(q_b)x_a=M(f_a f_b).                   (3)

H is strictly positive because the f_a are linearly independent and M has
full support. Let f' be the vector of nonconstant modes and

    K_A=conv{f'(z):z in S1},
    C_A={x0>0, x'/x0 in int K_A}.                      (4)

For m=1 use C_A=R_+. Linear independence with f0 makes K_A full dimensional.
For every good eta, q -> x is a smooth bijection onto C_A with Jacobian
det H. This includes the two-window cone |x1|<x0.

Proof of the onto claim. At fixed eta write
Z(t)=integral exp[t.f'(z)] B_eta(dz). For b in int K_A, the function
log Z(t)-t.b is strictly convex. It is coercive: uniformly in unit u,
max_z u.(f'(z)-b) is bounded below by a positive number. Uniform continuity
of f' and full support of B_eta give a positive lower mass for fixed-size
neighborhoods of these maximizers, uniformly in u (use a finite cover of the
circle). Thus log Z(r u)-r u.b >= c r-C for r large. Its unique minimum
solves grad log Z(t)=b; its Hessian is the positive covariance matrix of f'.
Take t=q' and q0=log[x0/Z(t)]. This proves the claimed inverse and smoothness.
The inverse is measurable in eta, for example by the unique strictly convex
minimization and rational approximations to its minimizer.

If p_h^A is the Gaussian density in (1), then

    (U_A psi)(x,eta)=sqrt[p_h^A(q(x,eta))/det H(q(x,eta),eta)]
                              psi(q(x,eta),eta)         (5)

is an onto unitary to L2(C_A,dx) tensor L2(mu_perp,A). This is a coordinate
change within LMD, not a change to the older MPQ Gaussian representation.

Fix sigma>0 and let k_sigma(t)=(2 pi sigma^2)^(-1/4)exp[-t^2/(4sigma^2)].
For Q in R^(m-1), filter all resolved nonconstant coefficients:

    K_Q(x,eta)=product_(a!=0) k_sigma(Q_a-q_a(x,eta)),
    W_(A,sigma) psi(Q,p,eta)
      =(2 pi h)^(-m/2) integral_(C_A) exp[i p.x/h]
                                  K_Q(x,eta) U_A psi(x,eta) dx.       (6)

Zero extension precedes the L2 Fourier transform; no convergence of an
absolute integral is assumed. Integrating K_Q^2 over Q gives one, so
Plancherel proves that W_(A,sigma) is an isometry. Its compressed multiplication
PVM is the normalized joint POVM E_(A,sigma) on the full H_h. For m=1 the
empty filter is one and this recovers the radial construction.

The complete record includes eta_A sharply. The Q marginal is q' with
independent N(0,sigma^2) errors. It is important not to confuse the exact
residual record with the noisy resolved coefficients.

For each outcome use the original source inverse:

    nu_(Q,eta)(dz)=exp[sum_(a!=0) Q_a f_a(z)] B_eta(dz)
                            /B_eta(exp[sum_(a!=0) Q_a f_a]),
    y_p(z)=sum_(a in A) p_a f_a(z),  theta_p=Theta(y_p),
    a(theta)=1+(1/2)cos(theta),  Theta=Y^(-1),
    Y(theta)=integral_0^theta dt/a(t),
    e[y_p,nu](dz)=[exp(theta_p)/a(theta_p)] nu(dz)
                                              /nu(1/a(theta_p)),
    L[y_p,nu]=integral e[y_p,nu].                       (7)

Each outcome has finite positive L, though its probability law can have
infinite raw length moments. The coframe is measure valued. It is not a
pointwise metric or a spacetime curvature operator.

Exact refinement would require, for A subset B, that the complete E_A record
be obtained as a classical, state-independent postprocessing of E_B. A weaker
necessary condition is a common joint POVM having E_A and E_B as marginals.
The next section disproves even this weaker condition for this family.

## 2. The complete records are not compatible across refinement

Take a nonconstant Fourier mode b in B\A. Its Gaussian coefficient q_b is
part of eta_A. Consequently the sharp spectral measure of q_b is a marginal
of E_(A,sigma_A), independently of sigma_A.

In the B chart define the bounded characteristic operator of p_b by

    T_t=integral exp[i t p_b/h] E_(B,sigma_B)(dQ dp deta).

Let e_b denote the unit vector in the mass coordinate x_b. Fourier inversion
and the overlap of the real Gaussian filters give, in U_B coordinates,

    (T_t g)(x,eta)=1_(C_B)(x-t e_b) R_t(x,eta) g(x-t e_b,eta),
    R_t=exp[-sum_(a in B,a!=0)
           (q_a(x,eta)-q_a(x-t e_b,eta))^2/(8 sigma_B^2)].             (8)

This formula is for x in C_B; R_t is needed only where both points lie in C_B.
It is strictly positive on that overlap for every finite sigma_B>0.

The cone is open and convex. For t>0 and a segment contained in it,

    q_b(x,eta)-q_b(x-t e_b,eta)
       =integral_0^t (H^(-1))_bb(x-u e_b,eta) du >0.    (9)

Choose any sufficiently small t>0 and a nonempty interior overlap. The
bounded strictly increasing function arctan(q_b) therefore does not commute
with T_t: their commutator has a nonzero multiplier on this overlap. This
also proves failure of commutation with at least one spectral projection of
q_b and at least one effect of the p_b marginal. The argument holds fiberwise
for almost every eta and hence on H_h.

Use the standard sharp-marginal joint-measurement criterion: if a PVM P and
a POVM F have a joint observable G, then their effects commute. Indeed
0<=G(X,Y)<=P(X) implies G(X,Y)=P(X)G(X,Y)=G(X,Y)P(X); splitting F(Y) into
X and its complement gives P(X)F(Y)=F(Y)P(X). Attribution is in SOURCE_MAP.md.
Equation (9) contradicts this necessary condition.

Thus E_(A,sigma_A) and E_(B,sigma_B) have no common joint POVM, and E_A cannot
be recovered by classical postprocessing of E_B, for a proper finite Fourier
refinement adding a nonconstant mode and finite positive resolutions.

This is a source-specific application of an established measurement theorem,
not a new general uncertainty principle. Its cause is the decision to retain
every unresolved configuration mode sharply. It does not prove that no
compatible quantum geometry exists. In particular it does not by itself
settle compatibility of the coframe-only pushforwards, which discard record
information. It also does not prevent convergence of other bounded marginal
readings, a differently designed noisy parent measurement, or a classical
limit. No inconsistency of the LMD Hilbert space or its unitaries follows.

## 3. Exact finite-family pulse law and a bounded geometric evolution equation

For a real w=sum_(a in A) beta_a f_a, the unchanged density pulse is
U_t=exp[-i t J(w)/h]. In (5) it multiplies by exp[-i t beta.x/h]. Hence

    W_(A,sigma) U_t psi(Q,p,eta)
                    =W_(A,sigma) psi(Q,p-t beta,eta).   (10)

The outcome density shape is fixed and y evolves exactly to y+t w. All
claims below hold at every fixed 0<h<1 and finite A,sigma; they use no
classical approximation. They concern the ideal readout under this probe.

Write theta_t=Theta(y+t w), and define the normalized source-mass weight

    pi_t(dz)=[nu(dz)/a(theta_t(z))]/nu(1/a(theta_t)).

Differentiating the full inverse (7), including its denominator, gives

    d/dt log(de_t/dnu)(z)
       =[a(theta_t(z))-a'(theta_t(z))] w(z)
                                +integral a'(theta_t) w dpi_t.       (11)

The second, spatially constant term is the original global P feedback.
Omitting it changes the evolution law. The formula is meaningful for the
positive Radon--Nikodym density relative to nu, including singular nu.

Set M=||w||_infinity. Since 1/2<=a<=3/2 and |a'|,|a''|<=1/2,

    |d log(de_t/dnu)/dt|<=5M/2,
    |d log L_t/dt|<=5M/2.                              (12)

In particular e^(-5M|t|/2)L_0<=L_t<=e^(5M|t|/2)L_0. Thus a finite pulse
in this span neither repairs nor creates a finite raw mean length from an
infinite one: integrability of L is equivalent before and after this pulse.
This statement compares the outcome laws through (10), not sharp length PVMs.

For a bounded reading choose F(L)=L/(1+L), and define

    A_F=W^* F(L) W,
    B_F=W^* { [L/(1+L)^2] integral b(z) de(z)/L } W,
    b(z)=[a(theta(z))-a'(theta(z))]w(z)+integral a'(theta)w dpi.       (13)

Here W=W_(A,sigma), and all quantities inside braces are outcome functions
at t=0. Both operators are defined on all H_h, and ||B_F||<=5M/8. More
strongly the Heisenberg evolution has an operator-norm derivative:

    d/dt [U_t^* A_F U_t] at t=0 = B_F,
    ||U_t^* A_F U_t-A_F-t B_F|| <= (15/8)t^2 M^2.       (14)

For completeness, let r_t=log(de_t/dnu), R_t=log L_t. Differentiating the
last term in (11) gives

    d/dt integral a'(theta_t) w dpi_t
       =integral a(theta_t)a''(theta_t)w^2 dpi_t
                                      -Var_(pi_t)(a'(theta_t)w).

Therefore |r_t''|<=5M^2/2 and
|R_t''|=|integral r_t'' de_t/L_t+Var_(e_t/L_t)(r_t')|<=35M^2/4.
The first two derivatives of R -> exp(R)/(1+exp(R)) have absolute value
at most 1/4. Combining with |R_t'|<=5M/2 gives |d^2F(L_t)/dt^2|<=15M^2/4.
Taylor's integral remainder is uniform over every outcome; compressing by
the isometry proves (14). On Dom J(w), its derivative agrees with the
commutator form (i/h)[J(w),A_F], which consequently has the bounded extension
B_F. No derivative of an infinite raw expectation is used.

This provides a genuine bounded observable evolution under the specified
source interaction. The parameter t is a probe pulse parameter. It has not
been identified with a lapse-generated gravitational time or a gravitational
Hamiltonian, and it supplies no spacetime curvature equation.

## 4. Fixed-family classical recovery for an arbitrary smooth source center

Let s_*,d_* be any smooth classical profiles with j_*=exp(s_*)>0 and
y_*=d_*/j_*. Define the source Gram matrix, moments and weighted projection

    H*_ab=integral f_a f_b j_*,  D*_a=integral f_a d_*,
    b_A=(H*)^(-1)D*,    P_A^(j_*) y_*=sum_(a in A) (b_A)_a f_a.     (15)

This is the L2(j_* dz) orthogonal projection onto V_A=span{f_a}. There is
no assumption that y_* already belongs to V_A. The finite-family classical
target is its weighted projection, not its full pointwise value.

Take the frozen LMD coherent state psi_h centered at s_*,d_*. For fixed A,
let sigma_h->0 and h/sigma_h->0. The law of (p,nu_(Q,eta)) under E_(A,sigma_h)
then concentrates at (b_A,nu_*), where nu_*=j_* dz/integral j_*. Thus the
coframe (7) converges against continuous test windows to
e[P_A^(j_*) y_*,nu_*], and bounded continuous functions of those readings
converge in expectation.

Here are the domain and error details extending the two-window proof. Under
the coherent configuration law q has covariance h diag(alpha_a) and mean q_*;
each H_ab=M(f_a f_b) has mean H*_ab and variance O_A(h) by LMD. Localize with
a fixed smooth cutoff chi(q,H), equal to 1 near (q_*,H*) and supported where
q is bounded and delta I<=H<=L I. Its norm loss tends to zero.

The support is contained in a compact subset of C_A. To check the boundary
point carefully, if x'/x0 approached a supporting face, there would be unit u
and k=max_z u.f'(z) with integral[k-u.f'] M tending to zero. The nonnegative
function g=k-u.f' is uniformly bounded at fixed A, so integral g^2 M also
tends to zero. But this is the H quadratic form of (k,-u), whose norm is
at least one, contradicting H>=delta I. x0 is bounded above and below as well.
This proves a zero boundary trace after localization; it is not an assumed
inverse-GMC moment estimate.

The inverse coordinate derivative is H^(-1), bounded on this support.
Also |partial_(q_a) log det H|<=m ||f_a||_infinity because
-||f_a||_infinity H<=partial_(q_a)H<=||f_a||_infinity H.
At fixed eta the coherent amplitude in (5) is proportional to
sqrt[p_h^A(q-q_*)/det H], with phase exp[-i D*.q/h]. Consequently

    ||(i h grad_x-b_A)U_A(chi psi_h)||^2<=C_A h,         (16)

using H^(-1)D*-b_A=H^(-1)(H*-H)b_A and the finite Gaussian coefficient
variances. C_A is not asserted uniform as A increases. Filter differentiation
adds exactly

    h^2/(4sigma_h^2) integral |U_A chi psi_h|^2
                                  sum_(a!=0)|grad_x q_a|^2 dx dmu_perp.

The cross term vanishes since integral K_Q grad_x K_Q dQ=0. Plancherel gives
the localized outcome bound C_A h+C'_A h^2/sigma_h^2 for |p-b_A|^2.
Chebyshev and isometry remove the cutoff in event probabilities. This proves
concentration without asserting a finite second moment for the full law.

For the shape marginal, Q_a-q_a are independent N(0,sigma_h^2) errors in a
coupling to the original configuration law. Their finite sum times f_a tends
uniformly to zero. The reconstructed normalized measure is therefore a
reweighting of M/J(1) by a factor tending uniformly to one. LMD concentration
and positivity of the limiting total mass give nu_(Q,eta)->nu_* weakly in
probability. The geometry map (7) is continuous in p and weak nu at each
finite p, with a denominator bounded below by 2/3. This completes the claim.

## 5. Arbitrary smooth density pulses: exact state evolution and its limit

Let B_t(z) be a prescribed real smooth pulse profile with B_0=0, continuous
in t in the smooth-function topology on a compact time interval. Set

    U_h(t)=exp[-i J(B_t)/h],   psi_h(t)=U_h(t)psi_h.      (17)

The J windows commute, so this is the exact integrated source pulse. It is
strongly continuous by dominated convergence. If B_t is differentiable in t,
it solves i h partial_t psi_h(t)=J(dot B_t)psi_h(t) on these coherent-state
orbits; J(1) has finite second moment for 0<h<1 and bounds the smooth windows.
No new evolution generator is added.

For B_t outside V_A there is no claim of a constant p translation at finite h.
Instead in (6) the exact phase multiplier is
exp[-i J(B_t)(q(x,eta),eta)/h], and its q derivative is M(B_t f_a).
The x phase gradient of the evolved coherent state is therefore

    H^(-1)[D*+ (M(B_t f_a))_a].                        (18)

The configuration law, cutoff and amplitude estimates are unchanged.
LMD's second-moment bound gives uniformly on the compact t interval

    E |M(B_t f_a)-integral B_t f_a j_*|^2=O_(A,T)(h).

The argument of section 4 now yields concentration at

    y_(A,t)=P_A^(j_*)[y_*+B_t],    nu_*,                (19)

uniformly in the sense of the supremum over t of each single-time error
probability. This is not the probability of a supremum of a jointly measured
quantum trajectory. No joint time-law or sequential instrument has been
assumed for arbitrary B_t.

## 6. Smooth-profile refinement after the classical limit

Now take A_N to contain the constant and both real modes through frequency N.
Let S_N be ordinary Fourier partial summation. For any smooth real y and
positive smooth j with 0<j_min<=j<=j_max, weighted projection and the finite
trigonometric norm bound give

    ||P_N^(j)y-y||_infinity
       <=||S_N y-y||_infinity
          +sqrt(2N+1) sqrt(j_max/j_min) ||y-S_N y||_L2.  (20)

Indeed P_N^(j)(y-S_Ny) belongs to V_N and its ordinary L2 norm is at most
sqrt(j_max/j_min)||y-S_Ny||_L2. For a degree-N trigonometric polynomial,
the supremum norm is bounded by sqrt(2N+1) times its L2 norm. Smoothness
makes the right side tend to zero; for smooth families on compact t intervals
it does so uniformly in t. This is a deterministic source estimate.

The nonlinear geometric inverse has a useful explicit stability bound. For
two continuous profiles y1,y2 on the same probability measure nu, let
delta=||y1-y2||_infinity. Since Theta'=a,

    |d/dy log[exp(Theta(y))/a(Theta(y))]|=|a-a'|<=2,
    |d/dy log[1/a(Theta(y))]|=|a'|<=1/2.

The denominator in (7) must be varied too. These estimates give pointwise

    exp(-5delta/2) <= de[y1,nu]/de[y2,nu] <= exp(5delta/2),
    ||e[y1,nu]-e[y2,nu]||_TV
                           <=(exp(5delta/2)-1)L[y2,nu].             (21)

Here TV is the total variation norm of the signed coframe measure, with no
factor 1/2. The same multiplicative bounds hold for the two circumferences.
Equations (20),(21) prove uniform recovery of the full smooth source coframe
along y_t=y_*+B_t after increasing N. The full P feedback is retained.

Combining with sections 4–5 establishes the iterated limit

    lim_(N->infinity) lim_(h->0)
       Law_(psi_h(t), E_(A_N,sigma_h)) [bounded geometry readings]
          = point law of e[y_*+B_t,nu_*],               (22)

uniformly in the stated single-time sense on compact t intervals. At each
fixed N, sigma_h->0 and h/sigma_h->0; for example sigma_h=sqrt(h).
There is no assertion that h and N limits commute, no uniform-in-N quantum
noise estimate, and no fixed-h continuum geometry inferred from (22).
Classical recovery of finer spatial structure and compatibility of quantum
measurements are distinct requirements.

## 7. Closure decisions and the next concrete design requirement

Established within the selected contract, by written derivation:

- Every finite Fourier family has the normalized mass-chart measurement.
- Its complete records cannot be nested by a common joint POVM when a new
  conjugate mode is added while the old sharp residual record is preserved.
- For pulses in the resolved span, F(L)=L/(1+L) has the exact bounded
  Heisenberg derivative (13) and the operator-norm remainder (14).
- Arbitrary smooth source centers and prescribed density-pulse profiles are
  recovered in the ordered classical/refinement limit (22).

Not established: a compatible fixed-h continuum geometry POVM, compatibility
or incompatibility of the coframe-only marginals, a material-diffeomorphism
covariant refinement rule, a finite-mean length domain, a local apparatus,
autonomous gravitational dynamics, spacetime curvature or 3+1 gravity.

Next measurement design requirement: do not record future conjugate modes as
sharp residual configurations in a proposed common parent. Either discard
those labels and analyze the geometry-only marginals, or replace them by a
single noisy parent measurement whose finite marginals are fixed in advance.
Then check both the source inverse and the actual density-pulse evolution on
that parent. The older CMG countable noisy parent is an established example
in the different MPQ representation; it does not automatically transfer its
geometric or Wick-density claims to LMD. Adding independent filters at each
cutoff without a common-parent proof is insufficient.

This identifies an exact failure in one refinement strategy and a controlled
source-dynamics bridge that survives it. It does not close the physical
quantum-gravity problem or improve the earlier NRQC performance comparison.
