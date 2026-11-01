# Bounded coframe readings with a prescribed angle profile

Author: Monty. 30 October 2026. Written derivation conditional on LMD-T3/T4/T6
and the frozen MPQ/RMG inverse chart. This is an additional readout contract,
not a quantum constraint on the conjugate field and not full reconstruction.

## Source identification before quantization

The original source spatial coframe is exp(r(z)) dz. In MPQ's exact inverse,

    y=d/j, theta=Y^{-1}(y), m=j/a(theta),
    P=integral m, r=theta+log|m/P|.

Restrict the *classical comparison sector* to constant y=Y(theta0), with
theta0 a fixed real calibration parameter. In the positive orientation sector,

    exp(r(z)) dz = exp(theta0) j(z) dz / integral j.       (1)

The a(theta0) factors cancel between m and P. The global feedback P has
therefore been retained, not frozen independently. This formula is false for
general nonconstant theta. A condition d=Y(theta0)j is imposed only on
classical centers used for comparison; no simultaneous sharp quantum values
of the conjugate fields, or normalizable joint constraint kernel, are assumed.

## Definition, domains and convergence

LMD gives a finite, nonzero, nonatomic full-support Wick measure M almost
surely on the already selected positive physical space L2(mu_h), 0<h<1.
Define the readout coframe measure and its smeared multiplication operator:

    nu_theta0(s,dz)=exp(theta0) M(s,dz)/M(s,S1),
    E_theta0(w) psi(s)=nu_theta0(s,w) psi(s).              (2)

For real bounded Borel w this operator is bounded self-adjoint, defined on the
entire Hilbert space, with norm at most exp(theta0)||w||_infinity. For an arc I,
0<=E_theta0(1_I)<=exp(theta0)I. Its spectral law is a well-defined probability
law for every normal state. The family has a joint configuration PVM because
it consists of measurable multiplication functions. This is not a joint PVM
for density and current, and it is not the old RMG joint detector.

No pointwise division by a distribution occurs in (2); the denominator is a
strictly positive scalar total mass. For bounded continuous w, replacing M
by the LMD Fourier-cutoff M_N gives almost-sure convergence of the ratio.
The same holds for each fixed arc because M has no atoms at its endpoints.
The ratios are uniformly bounded. Dominated convergence then gives strong
operator convergence E_theta0,N(w)->E_theta0(w) on every Hilbert-space vector.
This proves existence, not a quantitative cutoff error rate for simulation.

For each fixed sample, normalized cumulative mass is a continuous strictly
increasing circle coordinate, because M has no atoms and has full support.
For two distinct material points x,y, let I_xy be the oriented connecting arc.
The distance

    D_theta0(s;x,y)=min(nu_theta0(I_xy),
                       exp(theta0)-nu_theta0(I_xy))       (3)

is the pullback of the usual circle distance by that cumulative coordinate.
It is a metric and a bounded multiplication observable for fixed endpoints.
All its sample metric spaces are isometric to a circle of circumference
exp(theta0). The material-coordinate readings fluctuate, but this construction
does not create intrinsic one-dimensional curvature or spacetime dynamics.
The total circumference is a declared parameter, not a predicted gravitational
length. This limitation is part of the result.

## Compatibility with the existing material flow

For LMD's corrected coefficient a_h=1+h/2, T_g^* J_h(w) T_g=J_h(w composed g)
and total mass is invariant. Joint functional calculus therefore gives

    T_g^* E_theta0(w) T_g=E_theta0(w composed g).         (4)

The identity holds for bounded Borel w by the spectral measure identity; it
does not require an unbounded inverse-mass operator on a guessed domain.
For an indicator it sends I to g^{-1}(I) in the conjugated observable.
Samplewise, nu_theta0(F_g s,B)=nu_theta0(s,g(B)) for each fixed g.

In the old classical flow the adapted variable y transports as a scalar,
so the constant-y comparison sector is preserved. Equation (4) supplies the
quantum covariance of the declared reading. It does not independently
quantize the nonlinear source constraints or prove a representation of all
source canonical transformations.

## Classical recovery and its precise restriction

For a smooth coherent center s_* and d_*=exp(s_*)Y(theta0), LMD-T6 gives
joint convergence in probability of M_h(w) and M_h(S1) to integral w exp(s_*)
and integral exp(s_*), respectively, as h->0. The limiting denominator is
strictly positive. Thus (2) converges in probability, and in every finite
moment by boundedness, to the source coframe reading (1). Arc distances in
(3) have the corresponding limit.

For centers with other d_* the same configuration reading has the same law,
although the original inverse can have different theta and r. Hence (2)
cannot be the full phase-sensitive source geometry. It is a bounded
configuration-only extension with the correct restricted classical limit.
Calling it a solution of the unrestricted inverse problem would be false.

## Extension to a prescribed smooth nonconstant angle

The same source inverse gives more than a constant circumference if the
readout contract supplies a fixed smooth real profile theta(z). Define

    nu_theta(s,dz) = [exp(theta(z))/a(theta(z))] M(s,dz)
                    / M(s,1/a(theta)),                 (5)
    E_theta(w)=nu_theta(w), L_theta=E_theta(1).

The denominator is at least (2/3)M(S1)>0. Equivalently, L_theta is the
average of exp(theta) with respect to the probability measure proportional
to M/a(theta). Consequently

    exp(min theta)<=L_theta<=exp(max theta),
    ||E_theta(w)||<=exp(max theta)||w||_infinity.        (6)

All these multiplication observables are bounded self-adjoint on the full
physical Hilbert space. The same ratio argument proves their strong cutoff
limit, including fixed arcs. L_theta is a state-dependent random reading in
general; its allowed range is set by the prescribed profile. No universal
fixed value or gravitational length scale is predicted.

For smooth coherent centers with d_*=exp(s_*)Y(theta(z)), the classical limit
of (5) is exactly exp(r_*)dz from the original MPQ inverse, including its
global P. Thus this is a source-geometric reading on a larger stated family
of classical comparison sectors. The profile theta is supplied; it is not
obtained by dividing a quantum conjugate distribution by M.

Its flow law is a covariance of the *whole family* of calibrated readings:

    T_g^* E_theta(w) T_g=E_(theta composed g)(w composed g). (7)

Both numerator and denominator follow from LMD's corrected density law.
Holding a nonconstant profile fixed would not give the constant-profile
formula (4). The calibration profile has to transform as stated. This
distinction prevents an unjustified claim of symmetry at a fixed detector
setting. Constant theta0 reduces exactly to (2)–(4).

This extends the bounded geometric readout but leaves its main limitation:
the quantum conjugate-dependent angle and its dynamics remain unconstructed.
The fixed theta profile does not distinguish quantum states that differ
only in their conjugate phase.

## Relation to the executable work

The separate pilot tests the current–density covariance tensor, not the
spectral distribution of (2). Its two-qubit moment encoding must not be used
to claim an implementation of this normalized random-measure readout.
A genuine NRQC computation of (2) still needs a state/measure preparation,
finite-bin approximation and quantitative error control for the normalized
ratio. The existence proof above does not provide their computational cost.

This closes only the construction and classical-limit question for the
declared prescribed-angle coframe readouts. Endogenous conjugate-dependent
geometry, spacetime curvature, 3+1 constraints and independent physical
validation remain open. No numerical certificate is offered as a proof of
the measure, domain or continuum assertions.
