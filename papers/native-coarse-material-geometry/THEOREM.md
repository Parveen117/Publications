# CMG-1 — Continuum coarse geometry with a declared detection threshold

16 October 2026. Frozen parent: RMG-1 at
059d4b35d4fc306c5d9e65413970c2d446613edc. Work on exactly the selected
RGR/MPQ positive physical space. The source, gauge representation,
covariance, sign sectors and original J-pointer coupling are unchanged.
The additional measurement contract has two calibration parameters:
a material smoothing scale ell>0 and a minimum accepted density delta>0.
Neither parameter is derived as a law of gravity.

RMG showed that its raw point geometry has no fixed-hbar continuum limit
and its vacuum mean length diverges even at N=0. Here a countable
coherent-state dilation, spatial smoothing and an explicitly recorded
rejection outcome give a normalized continuum measurement. Its accepted
length has a dense closed positive form and a self-adjoint mean operator.
Both obstructions remain if the relevant calibration is removed.

## 1. The countable physical measurement, before nonlinear reconstruction

Write h=hbar>0 and H_epsilon=L^2(mu_h), mu_h=N(0,h C), with the exact
MPQ Fourier covariance lambda_0=1 and lambda_n=(1+n^2)^(-4), twice for
n>=1. The physical space is H_+ direct_sum H_-. This h-family is MPQ's
existing semiclassical convention, not a change to the residual gauge
measure in RGR. All constructions at fixed h include h=1.

Let kappa_h be the probability law of countably many independent pairs

    u_j ~ N(0,2h lambda_j), v_j ~ N(0,h/(2lambda_j)).       (1)

This is an ordinary countable product probability measure, not an
infinite-dimensional Lebesgue measure. Put

    zeta_j=u_j/(2 sqrt(h lambda_j))+i sqrt(lambda_j/h) v_j. (2)

The real and imaginary parts have variance 1/2. If H_a is the normalized
product of probabilists' Hermite polynomials in s_j/sqrt(h lambda_j),
with finite multi-index a, define

    W H_a=product_j zeta_j^(a_j)/sqrt(a_j!).              (3)

**CMG-T1 (countable POVM with the correct finite marginals).** W extends
to an isometry from H_epsilon into L^2(kappa_h). Thus

    E(B)=W^* 1_B W                                       (4)

is a normalized countably additive POVM on the full sequence of labels.
Its marginal on modes through N is exactly RMG-T1's coherent-state POVM
on H_N tensor I_H>N. Record epsilon as a central discrete outcome.
For every normal physical state, the label law is absolutely continuous
with respect to kappa_h. No individual normalized continuum coherent
vector at a typical outcome is asserted.

Proof. The circular Gaussian moments are
E[zeta^a conjugate(zeta)^b]=delta_ab a!. Product independence makes (3)
an orthonormal family, and the H_a form a complete basis of H_epsilon.
The moment identity follows by expanding
E exp(t zeta+s conjugate(zeta))=exp[((t+s)^2-(t-s)^2)/4]=exp(ts).
This proves the isometry. Multiplication projections on L^2(kappa_h)
then prove positivity, weak countable additivity and E(all)=I.

For one mode, direct Gaussian integration in RMG's vector chi_(u,v)
gives <chi,H_a>/<chi,1>=zeta^a/sqrt(a!), with
|<chi,1>|^2 du dv/(2pi h) equal to (1). Integrating unused modes in
matrix elements of (4) gives their Hermite Kronecker deltas. Hence the
finite marginals agree on a dense basis and therefore everywhere.
For a pure state psi its law has density |W psi|^2. For
rho=sum p_i |psi_i><psi_i| it has density sum p_i |W psi_i|^2, of integral
Tr rho. On the two-sign direct sum use the diagonal blocks of rho in
this formula. The measurement preserves the existing superselection
contract. This proves absolute continuity for all normal states.

This is a continuum outcome measurement built on the derived physical
norm. It does not require a continuum limit of RMG's measure-and-prepare
postmeasurement states. Finite readout instruments are given in T3.

## 2. An explicit filter and the exact coarse source profiles

Let Lambda e_n=n^2 e_n; equivalently Lambda=-partial_z^2/(2pi)^2.
For ell>0 select the resolvent filter

    R_ell=(I+ell^2 Lambda)^(-4),
    t_n=(1+ell^2 n^2)^(-4), t_0=1,
    u_ell=R_ell u, v_ell=R_ell v.                         (5)

It is a declared smooth material response with nonlocal support, not
a new time evolution or a local source apparatus. Use the unmodified
MPQ inverse chart on these filtered labels:

    j_ell=epsilon exp(u_ell), y_ell=v_ell/j_ell,
    theta_ell=Y^{-1}(y_ell), a(theta)=1+(1/2)cos theta,
    m_ell=j_ell/a(theta_ell), P_ell=integral m_ell,
    r_ell=theta_ell+log|m_ell/P_ell|,
    e_ell=exp(r_ell), L_ell=integral e_ell.               (6)

As in RMG, Y(theta)=integral_0^theta dt/a(t), and
Y^{-1}(y)=omega y+b(y), omega=sqrt(3)/2, B=||b||_infinity<infinity.
The global P term is retained. The outcome satisfies exactly

    rho_ell=m_ell'-m_ell r_ell'=-v_ell'+v_ell u_ell'.      (7)

Thus the reconstructed material-slice source constraints hold at the
regularity available below. Filtering adapted coordinates before their
inverse is a measurement prescription, not a canonical transformation,
nor the same operation as filtering the original r,m fields afterward.

**CMG-T2 (fixed-hbar continuum profiles and cutoff control).** Under
kappa_h, u belongs to H^3 and v to H^-5 almost surely. For each ell>0,
u_ell and v_ell belong to H^3 almost surely, hence (6) gives definite-sign
regular H^3 profiles and C^2 coframe data. These assertions hold almost
surely in every normal physical state. Fourier approximants
u_(ell,N)=P_N u_ell, v_(ell,N)=P_N v_ell converge in H^3, and their
reconstructed profiles converge in H^3 and C^2.

Proof. The raw expected norms are

    E||u||_H3^2=2h[1+2 sum_(n>=1)(1+n^2)^(-1)]<infinity,
    E||v||_H-5^2=(h/2)[1+2 sum_(n>=1)(1+n^2)^(-1)]<infinity.

The map R_ell is bounded H^-5 -> H^3 because
sup_n (1+n^2)^4/(1+ell^2 n^2)^4 is finite. Its exact filtered noise
budgets, with the same paired real Fourier normalization as MPQ, are

    E||u_ell||_H3^2
      =2h[1+2 sum_(n>=1)(1+n^2)^(-1)/(1+ell^2 n^2)^8],
    E||v_ell||_H3^2
      =(h/2)[1+2 sum_(n>=1)(1+n^2)^7/(1+ell^2 n^2)^8].   (8)

For N>=1 their tails obey respectively

    E||(I-P_N)u_ell||_H3^2 <=4h/N,
    E||(I-P_N)v_ell||_H3^2 <=128h ell^(-16)/N.            (9)

Indeed the summands are bounded by n^-2 and
128 ell^-16 n^-2, and sum_(n>N)n^-2<=1/N. Gaussian series convergence
gives the almost-sure statement and the Fourier limits. H^3 is an
algebra embedded in C^2; exp(u) stays strictly positive on the compact
circle. Compositions, division by j and a, and the global logarithm in
(6) are continuous H^3 operations at every realized profile. This
proves reconstruction convergence. T1's absolute continuity transfers
the almost-sure result to every normal state; it does not transfer
vacuum moment bounds to all states.

## 3. The detection threshold, normalized effects and finite instruments

Accept a geometry readout precisely on

    A_(ell,delta)={min_z exp(u_ell(z))>=delta}.           (10)

The reported outcome space is the disjoint union of accepted profile
pairs (epsilon,r_ell,m_ell) and an isolated symbol bottom meaning
"no geometry reading: density threshold not met". Keep this symbol and
its probability. The state is not deleted, renormalized or moved out
of the physical Hilbert space. Let F_(ell,delta) be this measurable
outcome map and E_(ell,delta) its pushforward of (4).

**CMG-T3 (normalized continuum detector and strong bounded-test limit).**
E_(ell,delta) is normalized on the original physical space. Every bounded
real output function gives a bounded self-adjoint operator on all H_red.
Replacing the filtered labels by their Fourier approximants gives, for
every bounded continuous output function f,

    W^* f(F_(ell,delta,N)) W -> W^* f(F_(ell,delta)) W
                                         strongly.     (11)

All normal-state outcome laws converge weakly. For any finite Borel
partition into detector bins whose boundaries have kappa_h-null inverse
images, the bin effects also converge strongly. The operators
K_b=sqrt(E_(ell,delta)(B_b)) define a normalized finite-bin instrument
rho -> K_b rho K_b; its cutoff instruments converge in trace norm on
each trace-class input. These are declared readout instruments, not
an implementation by the frozen J apparatus.

Proof. Normalization and bounded domains follow from T1. The threshold
boundary has kappa_h measure zero: u_ell=u_0+u_ell^perp, with u_0 an
independent nondegenerate Gaussian, and
min u_ell=u_0+min u_ell^perp. Conditioning on the perpendicular modes
makes equality to log delta a single zero-probability value of u_0.
T2 therefore implies almost-sure convergence of the outcome maps in
the disjoint-union topology, with H^3 on profile pairs. Bounded dominated
convergence against |W psi|^2 proves strong convergence of the
multiplication operators, hence their compressions in (11). For bins
the same proof uses the null-boundary assumption. Functional calculus
on [0,1] gives strong convergence of their square roots. Approximating
trace-class inputs by finite-rank operators proves trace-norm convergence
of K_b rho K_b. Finally sum_b K_b^*K_b=sum_b E(B_b)=I.

In RGR's selected quotient, the effects lift as I_gauge tensor E(B_b)
on the gauge-test tensor physical space and intertwine the quotient map.
Existing MPQ clock/pulse unitary conjugation gives the corresponding
relational measurements. None of these bounded constructions depends
on the still-open classification of additional residual gauge channels.

## 4. A positive length mean operator with a dense stated domain

The exact circumference can be written as a normalized weighted mean:

    L_ell=integral exp(theta_ell(z)) w_ell(z) dz,
    w_ell=|m_ell|/|P_ell|>0, integral w_ell=1.            (12)

This uses the definite sign and the global P feedback. On the accepted
event it gives the particularly useful bound

    0<L_ell<=exp[B+omega ||v_ell||_infinity/delta].        (13)

Define the nonnegative recorded number
f_L=1_A L_ell, assigning zero only as a bookkeeping value on bottom.
The rejection flag distinguishes this value from a zero circumference.

**CMG-T4 (length form, moments and physical domain).** The form

    q_L(psi,phi)=integral f_L conjugate(W psi) W phi dkappa_h,
    D(q_L)={psi:integral f_L |W psi|^2 dkappa_h<infinity}  (14)

is closed, nonnegative and densely defined. It contains every finite
Fourier Hermite polynomial and the vacuum. In fact these states have
all finite moments of f_L. The uniquely associated positive
self-adjoint mean operator M_(ell,delta) has domain

    {psi in D(q_L): there is eta in H with
       q_L(phi,psi)=<phi,eta> for every phi in D(q_L)},
    M_(ell,delta)psi=eta.                                (15)

On the sufficient second-moment domain {psi:f_L W psi in L^2(kappa_h)}
it acts as W^*(f_L W psi). Its spectral measure need not equal the
original length POVM; first-moment agreement is not a sharp joint
geometry construction. In particular no self-adjoint point r field
or canonical r,m commutator is inferred.

Proof. Formula (12) follows from e^r=e^theta |m/P|. On A, |y_ell|
is at most ||v_ell||_infinity/delta; (2) of RMG then gives (13).
The Gaussian H^3-valued v_ell has trace-class covariance by (8).
Consequently exp(c||v_ell||_H3) has finite expectation for every c.
One direct proof diagonalizes its covariance Q: for b||Q||<1,
E exp(b||v_ell||^2/2)=product(1-b q_j)^(-1/2)<infinity;
use a x<=b x^2/2+a^2/(2b) for any finite a. H^3 -> C^0 supplies
the same assertion for the supremum norm. The finite Hermite images
in (3) are polynomials in finitely many complex Gaussian coordinates.
Holder's inequality with (13) proves all their length moments finite.
They form a dense subspace of H. Multiplication by sqrt(f_L) is closed
on L^2(kappa_h), and its composition with the isometry W is closed.
Its norm square is (14). The closed positive form representation
theorem gives precisely (15), including the stated sufficient domain.

For fixed ell,delta the Fourier-approximant f_(L,N) converge to f_L in
all finite moments on that dense Hermite span: P_N is a contraction
on H^3, so all accepted lengths are dominated by the same exponential
of ||v_ell||_H3 in (13), up to the fixed Sobolev embedding constant.
This is a moment statement on an explicit domain, not an unproved
strong-resolvent statement about every unbounded cutoff operator.

If p_rho=Tr rho E(A)>0, the mean among accepted runs is
integral f_L dPr_rho / p_rho, whenever its numerator is finite.
This state-dependent ratio is not the expectation of a single linear
observable in general. Both the probability p_rho and the rejection
rate must accompany a conditional mean. Finite length moments are
proved on the specified states/domains, not for every normal state.

## 5. What removal of each calibration does

**CMG-T5 (the two independent obstructions survive).**

(a) Fix h>0 and ell>0. In the physical vacuum, the ungated smoothed
length still has infinite expectation. Thus

    E[1_(A_(ell,delta)) L_ell] -> infinity as delta downarrow 0. (16)

(b) Fix h>0 and delta>0. As ell downarrow 0, accepted point readings
r_ell(z) are not tight even conditionally on acceptance. For any finite R,

    Pr(A_(ell,delta) and |r_ell(z)|<=R) -> 0,
    Pr(A_(ell,delta)) -> p_delta with 0<p_delta<1.        (17)

Hence the fixed-hbar construction retains a material response scale
and a reported density threshold. It is not a removal of all regulators
or a proof that nature contains either calibration parameter.

Proof of (a). Write u_ell=u_0+u_perp and v_ell=v_0+v_perp, with the
two constant Gaussian modes independent of the perpendicular modes.
Choose M such that the event ||u_perp||_infinity,||v_perp||_infinity<=M
has positive probability. Intersect with epsilon v_0>=M+1, an
independent event of positive probability. Then at every z,
y_ell(z)>=exp(-u_0-M). By (12),

    L_ell>=exp[-B+omega exp(-M) exp(-u_0)].              (18)

The expectation of the right side over u_0~N(0,2h) is infinite:
as u_0 goes to minus infinity its double-exponential growth dominates
the Gaussian density. Tonelli proves the claim. Each profile has
positive minimum exp(u_ell), so the events A increase to the whole
label space as delta decreases to zero. Monotone convergence proves (16).

For (b), u_ell -> u in H^3 almost surely and ||u_ell||_H3<=||u||_H3.
The same constant-mode argument as T3 removes the threshold boundary
and gives the probability limit. Gaussian full support and an
independent constant mode give strictly positive acceptance and
rejection probabilities for every fixed delta>0. Independently,

    Var(v_ell(z))=(h/2)[1+2 sum_(n>=1)
                    (1+n^2)^4/(1+ell^2 n^2)^8]
                       is of order h ell^(-9).          (19)

The lower bound sums n between 1/(2ell) and 1/ell. The upper bound
splits at 1/ell, using n^8 below and ell^-16 n^-8 above. In particular
the variance diverges. RMG's bound on log|m/P| is still
2||u_ell||_infinity+log 3. On ||u_ell||_infinity<=M, bounded r therefore
forces v_ell(z) into a fixed bounded interval. Its Gaussian probability
tends to zero, independently of u and acceptance. The remaining event
has uniformly small probability for large M. This proves (17).

## 6. Classical profiles and actual length moments are recovered

**CMG-T6 (controlled classical limit including mean length).** Fix a
smooth regular classical profile r_*,m_* and its exact MPQ adapted
centers s_*,d_*. Choose

    0<delta<min_z exp(s_*(z)), ell(h) -> 0,
    h ell(h)^(-15) -> 0 as h -> 0.                      (20)

In the corresponding MPQ coherent states, acceptance tends to one,
the accepted profiles recover r_*,m_* in probability in H^3 and C^2,
and, for every finite p>0,

    E|1_A L_ell-L_*|^p -> 0,
    E[1_A L_ell]/Pr(A) -> L_*.                          (21)

Thus mean circumference as well as bounded readings recover the
classical value under this calibration contract. The limiting statement
is about outcome moments, not a unitary identification of quantum
representations at different h.

Proof. The coherent law shifts the means in (1) to s_*,d_* with the same
independent covariances, as follows either from T1's marginals or RMG's
overlap formula. The filtered means tend to the smooth centers in H^3.
For 0<ell<=1, (8) gives configuration noise O(h), while

    sum_(n>=1) (1+n^2)^7/(1+ell^2 n^2)^8
      <=ell^-14 sum_(n>=1) 1/(1+ell^2 n^2)
      <=(pi/2) ell^-15.                                (22)

The last inequality compares a decreasing function with its integral.
Hence momentum noise is O(h ell^-15). Markov's inequality, (20) and
continuity of (6) prove profile convergence; the strict threshold
margin gives acceptance probability tending to one.

It remains to justify unbounded moments. The centered Gaussian H^3
noise in v_ell has covariance trace tending to zero by (22), and thus
bounded by a finite T on a sufficiently small h interval. For such a
Gaussian Z, choosing b=1/(2T) in the product estimate of T4 gives
E exp(a||Z||)<=exp(a^2 T+1/2). The filtered mean H^3 norm is bounded.
Therefore (13) uniformly bounds every finite positive moment of
1_A L_ell, with fixed delta, including a strictly higher moment than
any chosen p. Convergence in probability to L_*, together with this
uniform integrability, proves (21). Unlike the ungated RMG length,
the necessary moment hypothesis is now established. No claim of an
optimal smoothing exponent or threshold is made.

## 7. Endpoint, scope and division of work

Constructed: a countable continuum measurement on the existing physical
space; fixed-hbar coarse C^2 geometry with cutoff control for all normal
states; normalized rejection-inclusive detector probabilities and finite
instruments; a positive self-adjoint accepted-length mean operator with
dense domain; and controlled classical recovery of actual length moments.
The threshold acts on measurement outcomes and preserves their original
inverse-chart values on acceptance. It is not a new source constraint.

Still open: calibration-independent sharp geometry, fixed-hbar removal
of the response scale/threshold, a spacetime curvature operator, a
derived local apparatus, source selection of covariance, critical patches,
matter/potential extensions and physical 3+1 quantum gravity. Coarse
material coframe data do not determine spacetime curvature. Standard
coherent-state, Gaussian, Sobolev and closed-form methods are credited;
their integration and the two separate obstruction/recovery results are
specific to the frozen source. Finite native witnesses are not analytic
proof-assistant certification or physical experiments.

External coordination: Fable's 943182d report on the separate
claude/residual-gauge-invariant-dual branch reduces the full continuous
invariant-functional question to residual-Diff ergodicity, conditional
on its stated dependencies. Its next assignment is Condition (L) and
repair of the spectator-form domain/weight claims. We do not assume
that classification or import its unreviewed stronger claims. This
packet uses only RGR's selected physical space and leaves that branch
untouched. Every frozen predecessor packet is preserved.
