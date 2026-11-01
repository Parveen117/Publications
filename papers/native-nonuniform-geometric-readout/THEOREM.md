# Nonuniform source geometry from a two-window LMD measurement

Author: Monty. 30 October 2026.
Frozen parent: `3b487951ccf5c4939a8c0249b435c904f47fedfe`.

The radial predecessor reconstructs the source geometry only when y=d/j is
constant. This packet first derives its exact response to an arbitrary density
pulse and exhibits the geometry it misses. It then constructs a different,
explicitly chosen joint POVM on the same positive LMD Hilbert space, resolving
the two conjugate coefficients of y(z)=p0+p1 cos(2 pi z). The density shape is
read with a declared finite resolution. The original inverse chart, including
its global P feedback, turns the outcomes into a coframe measure.

This is a written analytic construction conditional on frozen LMD. The new
POVM and its resolution are measurement choices, not consequences of the
gravitational constraints. Native rational calculations check only a finite
Taylor witness and algebraic diagnostics. No new NRQC run, continuum numerical
simulation, hardware result or physical gravity identification is reported.

## 1. What the previous radial measurement actually sees

Use the predecessor's coordinates (x,eta_r), where x=J(1)>0 and
nu=M/J(1) depends only on eta_r. Its Fourier-compression isometry is W_r.
For every real smooth window w and real pulse area beta,

    J(w)=x nu(w),     U_(beta,w)=exp[-i beta J(w)/h],
    W_r U_(beta,w) psi(y,eta_r)
       =W_r psi(y-beta nu(w),eta_r).                       (1)

This follows fiber by fiber from multiplication by exp[-i beta x nu(w)/h].
It is a shape-dependent shift, valid for arbitrary normal states and without
moment assumptions. It is not a spatially resolved angle measurement. In the
source classical limit it sees only y_bar -> y_bar+beta nu_*(w).

For example j_*=1, y_*=0, w(z)=cos(2 pi z) gives nu_*(w)=0. The classical
radial coframe stays dz, while the actual source y_beta(z)=beta cos(2 pi z)
is nonconstant. Section 7 proves that its circumference exceeds 1 for every
nonzero beta. At finite h, (1) can still change the old quantum law because
nu(w) fluctuates; no exact finite-h invariance is asserted for that example.

## 2. An onto two-moment coordinate chart within LMD

Put c(z)=cos(2 pi z), f0=1, f1=c. Decompose

    s=q0+q1 c+eta,       (q0,q1) ~ N(0,h I_2),           (2)

independently of eta, which omits precisely the constant and first cosine
modes. The variance of q1 is h: LMD's orthonormal sqrt(2)cos coefficient
has variance h/2. Let M_perp be the Wick measure of the remaining modes and

    B_eta(dz)=exp[-h(1+c(z)^2)/2] M_perp(eta,dz),
    M(q,eta,dz)=exp[q0+q1 c(z)] B_eta(dz).               (3)

Removing finitely many smooth Gaussian modes preserves positivity, finiteness,
absence of atoms and full support almost surely: (3) is an invertible bounded
positive reweighting for each finite q. In particular B_eta has these
properties. We retain all remaining continuum modes; this is not a truncation
of the LMD state space to two quantum degrees of freedom.

For a fixed good eta define

    F_eta(q)=integral exp(q0+q1 c) B_eta,
    x_a=partial_(q_a) F_eta=M(f_a),
    H_ab=partial_(q_b) x_a=M(f_a f_b),  a,b=0,1.        (4)

H is strictly positive: u^T H u=integral (u0+u1 c)^2 M>0 for u!=0,
by full support. The map q -> x is a smooth diffeomorphism onto the cone

    C={(x0,x1): x0>0, |x1|<x0}.                         (5)

Proof of surjectivity and uniqueness. Set Z(t)=B_eta(exp(t c)). Then
x1/x0=Z'(q1)/Z(q1). Its derivative is the strictly positive variance of c
under the tilted probability measure. As t tends to +/-infinity its mean
tends to +/-1, since B_eta gives positive mass to every neighborhood of the
maxima/minima of c. Thus each ratio in (-1,1) determines exactly one q1;
q0=log(x0/Z(q1)) then supplies each x0>0. The positive Jacobian det H gives
the smooth inverse. These formulas also give joint measurability in eta.

Let p_h(q) be the two-dimensional N(0,h I_2) density. The change of coordinates

    (U psi)(x,eta)=sqrt[p_h(q(x,eta))/det H(q(x,eta),eta)]
                         psi(q(x,eta),eta)              (6)

is an onto unitary H_h -> L2(C,d^2x) tensor L2(mu_perp,h).
Indeed d^2x=det H d^2q and the positive reciprocal multiplier is the inverse.
It is a unitary coordinate change within the selected LMD representation,
not an equivalence to MPQ's different smooth Gaussian measure class.

J(f_a) becomes multiplication by x_a. Classically the phase momenta conjugate
to x are p=H^(-1)D, with D_a=integral f_a d. If d=j(b0+b1 c), then D=Hb
and p=b exactly. For an arbitrary profile, H^(-1)D gives its weighted
least-squares projection onto span{1,c}, not its full pointwise value.
The operators D(1),D(c) are affine conjugates, not reference charges R[v].

## 3. A normalized joint readout with an explicit resolution

A Fourier measurement of both x coordinates alone discards q1, which changes
the normalized density shape. Reusing the old sharp shape reading would
therefore be unjustified. Instead choose sigma>0 and the real Gaussian filter

    k_sigma(t)=(2 pi sigma^2)^(-1/4) exp[-t^2/(4 sigma^2)],
    K_Q(x,eta)=k_sigma(Q-q1(x,eta)),     Q in R.          (7)

For each x,eta, integral K_Q^2 dQ=1. The filter first gives an unsharp reading
of q1. Follow it with zero extension V_C from C to R^2 and the plus-sign
Fourier transform in x. Define the isometry

    (W_sigma psi)(Q,p,eta)
       =(2 pi h)^(-1) integral_C exp[i p.x/h]
                       K_Q(x,eta)(U psi)(x,eta) d^2x.   (8)

The integral is understood in the L2 Fourier sense when it is not absolutely
convergent. Plancherel and (7) give ||W_sigma psi||=||psi||. Hence

    E_sigma(A)=W_sigma^* 1_A(Q,p,eta) W_sigma             (9)

is a normalized positive operator-valued measure on all of H_h. It extends to
mixed normal states by the trace. Every bounded real outcome function gives
a bounded self-adjoint effect integral. No self-adjoint derivative on C or
joint sharp mass/momentum PVM is assumed.

The (Q,eta) marginal is exactly

    integral_C k_sigma(Q-q1)^2 |U psi(x,eta)|^2 d^2x.

Thus Q=q1+N(0,sigma^2) in this configuration marginal. The eta coordinate
retains its configuration spectral reading; it commutes with the operations
in x. This is an ideal continuum measurement prescription, not a construction
of a spatially local apparatus. It changes the previous radial POVM: equality
of its marginals with the radial packet is neither assumed nor needed.

## 4. Endogenous profile and the original nonlinear coframe

For an outcome (Q,p,eta), define

    nu_(Q,eta)(dz)=exp(Q c(z)) B_eta(dz)/B_eta(exp(Q c)),
    y_p(z)=p0+p1 c(z),       theta_p(z)=Theta(y_p(z)),
    a(theta)=1+(1/2)cos(theta),
    Y(theta)=integral_0^theta dt/a(t),     Theta=Y^(-1),
    e_(Q,p,eta)(dz)
      =[exp(theta_p(z))/a(theta_p(z))] nu_(Q,eta)(dz)
                              /nu_(Q,eta)(1/a(theta_p)),
    ell_(Q,p,eta)=integral e_(Q,p,eta).                  (10)

Both coefficients of y are measurement outcomes, not supplied angle values.
Q is also an outcome; sigma and the chosen mode c are calibration choices.
The missing overall density scale exp(q0) cancels from the source coframe.
The denominator in (10) is P/J(1) for the corresponding normalized source
density. It retains the entire P feedback relevant to geometry; the absolute
value of P itself has not been jointly measured.

Equation (10) pushes (9) to a normalized coframe/circumference POVM. Each
outcome is a positive finite nonatomic full-support measure and

    exp(min_z theta_p) <= ell_(Q,p,eta) <= exp(max_z theta_p).

These bounds are outcome dependent, not uniform in p. Raw length moments
need not exist. Bounded functions and threshold events remain well-defined
for every normal state. The earlier infinite-mean warning is not removed.
This is a measure-valued spatial coframe reading, not a pointwise metric or
a spacetime curvature operator.

## 5. Exact nonuniform density-pulse evolution

For beta=(beta0,beta1), the existing source probe with window
beta0+beta1 c has unitary

    U_beta=exp[-i(beta0 J(1)+beta1 J(c))/h].

It commutes with each K_Q and multiplies U psi by exp[-i beta.x/h]. Therefore

    W_sigma U_beta psi(Q,p,eta)=W_sigma psi(Q,p-beta,eta). (11)

The outcome law is pushed forward by (Q,p,eta)->(Q,p+beta,eta). In particular,
the density-shape outcomes are unchanged and

    y_p(z) -> y_p(z)+beta0+beta1 c(z).                  (12)

Applying (10) after this translation reconstructs theta, the varying 1/a
factor and the changed denominator together. This is an exact measurement-law
identity for all normal states, including beta1!=0, not an equality between
nonexistent sharp local metric operators. A quantum pointer momentum is
included by spectral decomposition of beta_a=kappa_a Pi; no normalizable
exact momentum eigenvector is needed.

The full configuration probability law is preserved by U_beta, although the
p law shifts. A normalized probability law on R^2 cannot be invariant under
a nonzero translation. Thus the new readout detects preparation phases that
all configuration-only readings miss. It is not an NRQC implementation.

For windows outside span{1,c}, their J(w) is generally nonlinear in x at
fixed eta, so (11) does not give a constant two-vector shift. Arbitrary-window
local geometry is not closed. Nor does the fixed cosine calibration remain
the same under every material diffeomorphism: c composed with g generally
leaves this span. An all-flow covariant family requires further work.

## 6. Source classical limit with both readout errors controlled

Take the frozen LMD coherent state centered at smooth s_*,d_*, with

    j_*=exp(s_*),    d_*=j_* (b0+b1 c).

Write q_* for the two removed coefficients of s_* and
H*_ab=integral f_a f_b j_*, D*_a=integral f_a d_*. Then H*>0 and D*=H*b.
Under the coherent configuration law, q-q_* is N(0,h I_2) and

    E ||H-H*||^2=O(h),      E |q-q_*|^2=2h.            (13)

The first estimate follows from LMD's Wick two-point formula for the finite
family 1,c,c^2 and the integrable log-kernel estimate used in the radial
packet. In (6), the modulus of the coherent wavefunction is proportional,
at fixed eta, to sqrt[p_h(q-q_*)/det H]. Its q-dependent phase is
exp[-i D*.q/h]. The eta-dependent phase is independent of x.

To avoid unjustified inverse-moment bounds, multiply the state by a smooth
configuration cutoff chi(q,H), equal to 1 near (q_*,H*) and supported where
q is bounded and delta I <= H <= L I for fixed positive delta,L. Its norm
loss tends to zero by (13). Its x derivatives are bounded. This support
lies in a compact subset of C: x0 is bounded above and below, while
det H <= x0^2-x1^2 bounds the distance from the cone edges below.

On this support, partial_x q=H^(-1) is bounded. Moreover
partial_(q0) log det H=2 and |partial_(q1) log det H|<=2; the latter follows
from -H <= partial_(q1)H <= H, since |c|<=1. Differentiating the real
amplitude in (6) therefore gives an O(h) bound on its squared h-gradient
norm. The phase gradient is H^(-1)D*, and

    H^(-1)D*-b = H^(-1)(H*-H)b.

Its localized squared mean is O(h) by (13). Combining these estimates gives

    ||(i h grad_x-b) U(chi psi_h)||^2=O(h).             (14)

The cutoff ensures zero boundary trace. For the filter there are the exact
pointwise identities

    integral K_Q grad_x K_Q dQ=0,
    integral |grad_x K_Q|^2 dQ
                  =|grad_x q1|^2/(4 sigma^2).          (15)

Apply Plancherel to (8) and use (15) to cancel the cross term. The unnormalized
cutoff outcome law consequently obeys

    integral |p-b|^2 |W_sigma(chi psi_h)|^2 dQ d^2p dmu_perp
                    <= C1 h+C2 h^2/sigma^2.            (16)

Constants depend on the fixed source center and cutoff, not h,sigma. If

    sigma_h -> 0,       h/sigma_h -> 0,

for example sigma_h=sqrt(h), (16), Chebyshev and the vanishing norm loss
imply p->b in probability in the full outcome law. This does not assert a
finite second moment for the unlocalized law. Isometry makes the difference
of the two event probabilities tend to zero uniformly in the event.

The (Q,eta) marginal has a particularly direct coupling: draw a configuration
from |psi_h|^2 mu_h and add independent Gaussian noise to q1. For that draw,

    nu_(Q,eta)(dz)
       =exp[(Q-q1)c(z)] nu_(q1,eta)(dz)
                           /nu_(q1,eta)(exp[(Q-q1)c]),
    nu_(q1,eta)=M/J(1).                                  (17)

Since |c|<=1 and Q-q1->0 in probability, this reweighting tends to identity
in total variation. LMD density concentration therefore gives weak
concentration of nu_(Q,eta) at nu_*=j_* dz/integral j_*. Union bounds combine
this with p->b; independence between p and shape is not needed. Compactness
of the circle and continuity of (10) now give the original source coframe
in probability against each finite family of continuous test windows.
Bounded continuous functions of these readings converge in expectation;
raw unbounded moment convergence is not inferred.

Thus the classical comparison sector includes genuinely nonconstant
y_*(z)=b0+b1 c(z), arbitrary smooth positive j_*, and its entire pulse orbit
under windows in span{1,c}. This is strictly larger than the predecessor's
constant-y sector, but does not include every smooth y.

## 7. A source-specific witness that the radial average misses

Set j_*=1, y_*=0 and apply beta c. Write theta_beta(z)=Theta(beta c(z)).
The actual source circumference is

    L(beta)=[integral exp(theta_beta)/a(theta_beta) dz]
                                      /[integral 1/a(theta_beta) dz]. (18)

a is even, Y and Theta are odd, and c(z+1/2)=-c(z). Pairing the two half
circles replaces the numerator by integral cosh(theta_beta)/a(theta_beta).
Consequently L is even and L(beta)>1 for every beta!=0: cosh(theta)>1
except at the zero set of c. In contrast, the classical radial ratio is
integral j y/integral j=0 for all beta and predicts circumference 1.

For an explicit local diagnostic, Theta'(0)=3/2, Theta''(0)=0 and
Theta'''(0)=-9/8. With c=c(z), the source expansions are

    P(beta)=integral 1/a(theta_beta)=2/3+beta^2/8+O(beta^4),
    e_beta(dz)=[1+(3/2)beta c
                    +beta^2((3/2)c^2-3/16)+O(beta^3)] dz,
    L(beta)=1+(9/16)beta^2+(27/1024)beta^4+O(beta^6),
    integral c e_beta=(3/4)beta+O(beta^3).               (19)

Freezing P at P(0) would instead give the incorrect quadratic coefficient
3/4 in L. The native exact jet calculation checks through degree six, with
coefficients recorded in VALIDATION.json. As a separate analytic cross-check,
the tangent half-angle solution of Theta'=a(Theta) gives globally

    1/a(Theta(t))=4/3-(2/3)cos(sqrt(3)t/2).

The rational jet of this identity independently checks the denominator.
Analyticity near beta=0 justifies the displayed Taylor remainders; the finite
calculation does not supply a numerical bound on a remainder. These are classical source
identities and bounded-reading classical-limit targets. They are not finite-h
quantum expected lengths, sampled detector data or calculator gravity outputs.

## 8. Closure boundary and the next obligation

Constructed under the frozen LMD inputs and stated ideal measurement:
the exact arbitrary-window response of the old radial POVM; an onto two-mass
chart; a normalized joint conjugate/shape POVM; a positive coframe pushforward;
exact nonuniform pulse translation for span{1,c}; and the corresponding
source classical limit with the two measurement errors controlled.

New choices are the cosine calibration, Gaussian q1 filter and its sigma,
and Fourier compression on the mass cone. No uniqueness or local apparatus
derivation is supplied. A density shape and conjugate profile have a joint
unsharp measurement; no joint sharp canonical fields have been constructed.

Next: determine whether increasing smooth-window families can be made
consistent with one another and material-flow covariance, with a controlled
limit of the reconstructed bounded geometry laws. Merely adding modes does
not prove such a limit. Finite-mean length domains, a local apparatus, an
independently ordered unreduced source, spacetime curvature and 3+1 gravity
remain open. The NRQC continuum/state-preparation gap remains open as well.
