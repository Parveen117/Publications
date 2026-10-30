# State-dependent radial geometry from the LMD affine pair

Author: Monty. 30 October 2026.
Frozen parent: `cf72a467522b3516d3c3c5531c7666ff2da96a2c`.

The previous bounded coframe required a supplied angle profile. This packet
replaces the **global constant angle** by an outcome of an explicit normalized
POVM on the existing LMD physical space. It derives its uniform-probe pulse
law and a restricted source classical limit. It does not reconstruct a local
angle field. The measurement prescription is a stated addition, not a unique
consequence of the gravitational constraints.

All continuum statements below are written arguments conditional on the frozen
LMD results. Half-line momentum and its Fourier-compression POVM are standard
constructions; SOURCE_MAP.md credits them. No new numerical, NRQC or hardware
run is reported, and no finite check is offered as a continuum proof.

## 1. The existing source pair and a genuine radial factorization

Use the positive orientation sector H_h=L2(mu_h), 0<h<1, in LMD's dimensionless
canonical units. Define

    J0=M(S1)>0,             D0=D(1),
    [D0,J0]=i h J0.                                      (1)

D0 is the total affine conjugate integral d. It is **not** the material
reference charge R[1]. Their source meanings and generators remain separate.

Write s=q+eta, where q is the constant Fourier coefficient, q~N(0,h), and eta
has zero spatial mean and is independent of q. Let M_perp be the Wick measure
using only the nonconstant modes. The full variance subtraction gives

    M(s,dz)=exp(q-h/2) M_perp(eta,dz),
    Z(eta)=exp(-h/2) M_perp(eta,S1),
    x=J0=exp(q) Z(eta),
    nu_eta(dz)=M_perp(eta,dz)/M_perp(eta,S1).             (2)

LMD-T3 implies 0<Z<infinity almost surely. Its positivity, absence of atoms
and full support transfer to nu_eta. Let p_h be the N(0,h) density and put
K_h=L2(mu_perp,h). The map

    (U psi)(x,eta)=sqrt[p_h(log(x/Z(eta)))/x]
                     psi(log(x/Z(eta))+eta), x>0,        (3)

is an onto unitary H_h -> L2(R_+,dx) tensor K_h. Indeed dq=dx/x at fixed eta;
substitution proves norm equality, and the reciprocal multiplier gives an
inverse on the entire L2 space. This is a spectral coordinate transformation
within LMD, not an equivalence with the old smooth MPQ measure class.

Under U the two operators are

    U J0 U^*=x,
    U D0 U^*=i h(x partial_x+1/2).                       (4)

For the second formula, LMD gives D0=i h partial_q-i q/2. Differentiating
sqrt[x/p_h(q)] in the inverse of (3) cancels the Gaussian q term and leaves
the 1/2 dilation term. Equivalently, conjugate the actual Cameron--Martin
shift group; it becomes the unitary dilation group on dx. This specifies
the self-adjoint generator and its Stone domain, beyond a formal expression.

The normalized configuration shape M(w)/J0 is multiplication by nu_eta(w)
on K_h alone. Thus it is independent of the radial coordinate x.

## 2. The sharp-ratio obstruction, with its exact scope

Classically, if y(z)=y_* is constant then D0=y_* J0. A natural symmetric
candidate for the global ratio is

    A=J0^(-1/2) D0 J0^(-1/2).

On U^*[C_c^infinity(R_+) tensor K_h], it becomes i h partial_x. Its closure
has domain H1_0(R_+;K_h), with the zero boundary trace. Its adjoint is the
same derivative on H1(R_+;K_h). On one scalar radial fiber the equations

    (A^*-i)f=0 -> f(x)=c exp(x/h),
    (A^*+i)f=0 -> f(x)=c exp(-x/h)

give deficiency indices (0,1). With multiplicity K_h they are (0,dim K_h).
There is no self-adjoint extension of this minimal derivative on the same
positive-orientation space. Consequently a sharp self-adjoint ratio cannot
be obtained merely by selecting another boundary condition for this operator.

This is the standard half-line obstruction, already recognized in HQS-T4
and in MPQ's positive-spectrum Weyl boundary. Its use here concerns the
specific LMD total-mass ratio, not a new general no-go theorem. It does not rule out a different
observable, a different ordering, an enlarged space, or boundary couplings
between orientation sectors. None of those modifications is made here.

## 3. A normalized endogenous-angle measurement

Let V0 extend a half-line wavefunction by zero into L2(R), and let

    (F_h^+ f)(y)=(2 pi h)^(-1/2) integral_R exp(i y x/h) f(x) dx.

This sign convention diagonalizes i h partial_x as multiplication by y.
The isometry

    W=(F_h^+ V0 tensor I_K) U                            (5)

defines a joint POVM of the radial readout y and the configuration shape:

    E(B)=W^* 1_B(y,eta) W,
    Prob_psi(dy,deta)=|(W psi)(y,eta)|^2 dy mu_perp,h(deta). (6)

Positivity and E(all)=I follow directly from the isometry. Mixed normal states
are treated by the trace, or a positive trace-class decomposition. For every
bounded real outcome function f, W^* f W is bounded self-adjoint on all H_h.

The y marginal is generally not a PVM: the range of zero extension is not
reducing for the full-line momentum spectral projections, since its generated
translations do not preserve support on the half-line in both directions.
On the finite-second-moment domain, W^* y W equals A with domain
U^*H1_0(R_+;K_h), by Plancherel and the zero-extension Sobolev trace condition.
No larger maximal first-moment domain is asserted.

Normalized shape and this y measurement have a joint law. Sharp total mass
and y do not thereby acquire a joint PVM. The Fourier compression is an
explicit ideal measurement choice, not a construction of a local laboratory
apparatus or a claim of uniqueness among covariant POVMs.

## 4. Geometry obtained from the outcome, not a supplied global angle

Keep the frozen source function

    a(theta)=1+(1/2)cos theta,
    Y(theta)=integral_0^theta dt/a(t), Theta=Y^(-1).

For each joint outcome (y,eta), define

    theta_out=Theta(y),
    ell_out=exp(Theta(y)),
    e_out(dz)=ell_out nu_eta(dz).                        (7)

Push (6) forward by this measurable map. It is a normalized coframe and
circumference POVM: every outcome has a finite positive circumference,
but unbounded outcome moments need not exist. All bounded functions and
threshold events of (7) are well-defined effects on the full physical space.
For example, ell/(1+ell) and the event ell>1 are bounded readings.

Equation (7) is the original source inverse in the **constant-y classical
comparison sector**: when d(z)=y_* j(z), theta=Theta(y_*) is constant,
m=j/a(theta), P=(integral j)/a(theta), and exp(r)dz=exp(theta)j dz/integral j.
The global P feedback has not been discarded. No external angle value is
supplied in the quantum readout (7).

For nonconstant classical y(z), the ratio integral d/integral j is a
mass-weighted average. Applying Theta to this average does not reconstruct
the original local theta(z) or coframe. The present result is therefore a
global radial geometry reading, not the unrestricted inverse chart. A spatial
circle metric and its circumference also do not specify spacetime curvature.

## 5. Existing material-flow covariance and uniform-probe dynamics

For LMD's corrected material flow, total mass is invariant. Its underlying
configuration map commutes with s->s+t. The corresponding positive half-density
unitaries commute as well, by their exact Radon--Nikodym chain rule. Therefore
T_g commutes strongly with both J0 and the dilation group generated by D0.

In (3), an operator commuting with all spectral functions of x is decomposable
in x. Commutation with every dilation makes its operator-valued multiplier
constant in x (apply translation invariance in log x to matrix elements).
Consequently U T_g U^*=I_radial tensor V_g for a unitary on K_h. Thus the y
marginal and circumference are invariant, while the shape obeys

    V_g^* nu_eta(w) V_g=nu_eta(w composed g).            (8)

This proves covariance of the joint coframe readout under the already defined
material flows. It imports no new residual-gauge completeness result.

For a real uniform pulse parameter beta, the existing density interaction is

    U_beta=exp(-i beta J0/h).

In (3) it multiplies by exp(-i beta x/h); hence exactly

    W U_beta psi(y,eta)=W psi(y-beta,eta),
    U_beta^* E_y(B) U_beta=E_y(B-beta).                 (9)

The readout shifts y->y+beta, so ell becomes exp(Theta(y+beta)). This is the
frozen uniform-window source pulse law y_T=y_0+Pi alpha(T). A quantum pointer
is included by spectral decomposition of Pi, giving the conditional shift
beta=kappa pi with its pointer distribution. No normalizable exact pointer
momentum eigenstate is required for that direct-integral statement.

These are dynamics under an already specified probe coupling, not a newly
derived gravitational Hamiltonian. A nonconstant window produces a local
y profile and leaves the constant-y classical sector; its full geometry
cannot be inferred from (9).

## 6. Restricted classical limit, without inverse-mass moment assumptions

Take LMD coherent states centered at smooth s_* and
d_*(z)=exp(s_*(z))y_*. Put J_*=integral exp(s_*)>0 and D_*=y_* J_*.
LMD's Gaussian calculation gives

    ||(D0-D_*)psi_h||^2=h/4,
    ||(J0-J_*)psi_h||^2=O(h).                           (10)

For the second estimate, the Wick covariance kernel obeys
|exp(h K)-1|<=h |K| exp(h0 max(K,0)), whose diagonal singularity is integrable
for any fixed h0<1. The same reasoning gives concentration of each M(w).

Choose a smooth compact cutoff chi in (0,infinity), equal to one near J_*.
Then phi_h=chi(J0)psi_h belongs to the domain of A, and
||psi_h-phi_h||->0. On this localized domain,

    (A-y_*)phi_h
      =chi(J0)J0^(-1)(D0-D_*)psi_h
       +y_* chi(J0)J0^(-1)(J_*-J0)psi_h
       +i h chi'(J0)psi_h-(i h/2)chi(J0)J0^(-1)psi_h.   (11)

All inverse factors here are bounded by the fixed cutoff. Equation (10)
therefore makes the norm in (11) O(sqrt(h)). Plancherel identifies its square
with the y second centered moment of W phi_h. Chebyshev, followed by the
isometry estimate for psi_h-phi_h, gives concentration of the full y POVM
at y_*. No negative-moment theorem for total GMC mass has been assumed.

The commuting shape readings concentrate at
nu_*(w)=integral w exp(s_*)/J_*. The joint law consequently concentrates at
(y_*,nu_*), for every finite family of such windows. Thus (7) recovers the
original constant-y source coframe in probability, and bounded continuous
functions of its readings converge in expectation. Events require their
threshold to avoid the limiting value. Unbounded length-moment convergence
does not follow and is not claimed.

## 7. An exact phase-sensitive witness and a moment failure

Choose any normalized chi in K_h and the following radial state, a>0:

    (U psi_(a,b))(x,eta)=sqrt(2a) exp(-a x) exp(-i b x/h) chi(eta).

All b give exactly the same full configuration law, hence the same statistics
for every density-only reading, because their wavefunction moduli agree.
The Fourier integral in (5), evaluated by the elementary Laplace integral,
gives the normalized y density

    p_(a,b)(y)= (a h)/[pi ((y-b)^2+(a h)^2)].           (12)

This is a Cauchy law. The parameter b is its location, **not its mean**;
the ordinary y expectation does not exist. The state is normalizable and has
finite moments of total mass, while its radial boundary value is nonzero.
It is not in the finite-second-y-moment domain H1_0.

Since Theta is strictly increasing and Theta(0)=0, ell>1 iff y>0. Therefore

    Prob(ell>1)=1/2+(1/pi) arctan(b/(a h)).             (13)

| b/(a h) | Probability that the dimensionless circumference exceeds 1 |
|---|---|
| -1 | 1/4 |
| 0 | 1/2 |
| +1 | 3/4 |

These are analytic predictions for the selected ideal readout. They are not
NRQC outputs, hardware observations, or experimentally measured gravity.
The measurement kernel is fixed across these states; b changes the state
preparation phase and is not an angle inserted into the readout rule.
They explicitly demonstrate phase sensitivity missing from a prescribed-angle
configuration-only construction. Density laws are compared on separate
ensembles; simultaneous sharp density and angle measurement is not asserted.

More generally, (9) preserves every configuration probability but shifts the
y law. A probability measure cannot be invariant under a nonzero translation
on R. The bounded strictly increasing reading ell/(1+ell) increases strictly
under beta>0 in every normal state, by the pointwise shift formula.

The same witness also limits the result. RMG's source inverse has
Theta(y)=sqrt(3)y/2+b_0(y), with bounded periodic b_0. Thus exp(Theta(y)) grows
exponentially as y->+infinity, whereas (12) decays quadratically. Its raw mean
circumference is infinite. Threshold probabilities and bounded readouts remain
valid; this example does not establish a finite-mean length operator.

## 8. What is established and what remains open

Established, conditional on the frozen LMD inputs and the stated ideal POVM:
an onto radial factorization; the scope of the sharp-ratio obstruction; a
normalized state-dependent global-angle/coframe measurement; exact phase
sensitivity; material-flow covariance; uniform-probe pulse evolution; and
the constant-y source classical limit for bounded readings.

New choice: the Fourier-compression measurement prescription. It replaces a
supplied global angle by an outcome, but is not uniquely selected by gravity.
The original physical space, positive orientation, source inverse and uniform
probe coupling are retained. The construction is not a new general theorem
about half-line momentum, POVMs, positivity or Gaussian multiplicative chaos.

Open: a local angle field from the conjugate distribution; geometry after a
nonuniform density pulse; a finite-mean length domain for an appropriate
measurement/state family; a physical local apparatus; an independently ordered
unreduced source; spacetime curvature and 3+1 gravity. The next concrete source
test is the nonuniform pulse, where the single radial angle loses information.
