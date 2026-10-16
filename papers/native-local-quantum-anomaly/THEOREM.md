# LQA-1 — Inhomogeneous constraint quantization: anomaly and conditional repair

## 1. Target and normalization

Start from LNC's actual canonical chart q=log e, p=e pi, with scalar pair
(phi,k), where k=p_phi is a field density here. Select V=0 and no matter:

\[
h=pk+q'\phi'-\phi'',\qquad
d=pq'+k\phi'-p'.                                      \tag{1}
\]

Use REC's periodic independent-coframe boundary choice, x in R/Z, smooth
real periodic q,phi, and periodic momenta. The Poisson delta integrates to
one over this circle. A global fixed-winding angle lift is not imposed.
The action scale is hbar=1 except where restored explicitly. This is a
selected 1+1 model; the extra scalar and canonical structure are inherited
choices, not deductions from a diagonal drawing.

For chirality epsilon=+1 or -1 define

\[
A_\epsilon=p+\epsilon\phi',\quad B_\epsilon=k+\epsilon q',
\qquad T_\epsilon=(h+\epsilon d)/2.
\]

**LQA-T1 (source factorization).**

\[
T_\epsilon=\tfrac12 A_\epsilon B_\epsilon
                   -\tfrac\epsilon2 A_\epsilon',\quad
\{A_\epsilon(x),B_\epsilon(y)\}=2\epsilon\partial_x\delta(x-y).
                                                               \tag{2}
\]

Same-species brackets vanish and opposite chiralities commute. To prove (2),
expand the product: its derivative term is -phi''/2-epsilon p'/2.
The two cross canonical brackets each contribute epsilon partial_x delta;
for opposite chiralities they cancel. Thus the derivative improvement is
along a **null current**: its self-bracket is zero.

The following asymmetric Fourier normalization avoids square-root constants:

\[
a_m^\epsilon=\tfrac12\int_0^1 e^{2\pi i\epsilon m x}A_\epsilon(x)dx,
\quad b_m^\epsilon=\tfrac1{2\pi}\int_0^1e^{2\pi i\epsilon m x}B_\epsilon(x)dx,
\quad L_m^\epsilon=\tfrac1{2\pi}\int_0^1e^{2\pi i\epsilon m x}T_\epsilon(x)dx.
                                                               \tag{3}
\]

Quantize the current brackets by [F,G]=i{F,G}. Integration by parts gives

\[
[a_m^\epsilon,b_n^\eta]=\delta_{\epsilon\eta}m\delta_{m+n,0}I,
\quad[a_m^\epsilon,a_n^\eta]=[b_m^\epsilon,b_n^\eta]=0,
\quad (a_m^\epsilon)^*=a_{-m}^\epsilon,\ (b_m^\epsilon)^*=b_{-m}^\epsilon.
                                                               \tag{4}
\]

Specifically, the unscaled double Fourier integral of the quantum cross
bracket is 4 pi m delta; the product of the two prefactors in (3) is
1/(4 pi). Opposite Fourier orientations give the same sign in both copies.
Periodicity makes a_0^+=a_0^-=a_0 and b_0^+=b_0^-=b_0. The zero modes
are shared, not doubled independent canonical momenta.

## 2. Representation and common domain

For each chirality let all a_n,b_n with n>0 annihilate a cyclic vacuum.
Let a_-n,b_-n be commuting polynomial creation variables. On their polynomial
space, a_n=n partial/partial b_-n and b_n=n partial/partial a_-n.
Both zero modes act by real multiplication on Schwartz functions of
(a_0,b_0). Define D as that Schwartz space tensored algebraically with
finite polynomials in the two chiral families. Equivalently calculations
can be done in a fixed-real-zero-mode fiber. No finite upper mode is removed
from the representation. The auxiliary zero-mode measure is a chosen
kinematic measure, not HQS's physical group-averaged measure.

The vacuum functional with (4) defines a nondegenerate indefinite Hermitian
form. It is not a positive Hilbert norm. A Hilbert majorant can be chosen
using the sum/difference oscillators; D is a dense algebraic domain, but
star below denotes the indefinite-form adjoint relation on D. No assertion
of self-adjoint field operators on a positive Hilbert space is being made.

Normal ordering puts every positive-index annihilator to the right of
negative-index creators; zero modes commute. Set

\[
S_m=\sum_{r\in\mathbb Z}:a_{m-r}b_r:,
\qquad L_m=S_m+i m a_m.                              \tag{5}
\]

Equation (3) gives precisely (5): the quadratic term has coefficient one,
and integrating the derivative gives i m a_m. The same definitions apply
to both chiralities.

Each sum in (5) acts finitely on each vector of D. A normal-ordered term
with an annihilator vanishes unless its index occurs among the finitely
many occupied modes; terms with two creators require m<0 and have only
finitely many partitions of -m. Multiplication by zero modes preserves
Schwartz space. Hence all S_m,L_m and their finite products preserve D;
commutators have an honest common invariant domain. This does not establish
convergence of arbitrary infinite smearing sums, exponentials of q(x),
or a local physical inner product. Finite trigonometric smearings suffice
for the obstruction below.

## 3. Exact all-mode anomaly

**LQA-T2 (selected oscillator algebra).** On D,

\[
[S_m,a_n]=-n a_{m+n},\quad [S_m,b_n]=-n b_{m+n},
\]
\[
\boxed{[L_m,L_n]=(m-n)L_{m+n}
             +\frac{m^3-m}{6}\delta_{m+n,0}I.}         \tag{6}
\]

The two chiral copies commute. In standard Virasoro notation (6) has c=2
per copy, including the inherited derivative term.

**Proof.** Commuting a current through the normal-ordered quadratic sum
leaves a single nonzero contraction, giving the first line. Jacobi then
shows that [S_m,S_n]-(m-n)S_{m+n} commutes with every nonzero current.
At fixed zero modes, the polynomial oscillator module has scalar commutant:
commuting with creators and their cross derivatives fixes an operator by
its action on the vacuum, and that action is constant. Mode grading makes
this scalar vanish unless m+n=0.

For m>0, S_m annihilates the vacuum and

\[
S_{-m}\Omega=a_0b_{-m}\Omega+b_0a_{-m}\Omega
       +\sum_{j=1}^{m-1}a_{-j}b_{-(m-j)}\Omega.
\]

The vacuum coefficient of S_m S_-m is
2m a_0 b_0 + sum_{j=1}^{m-1}j(m-j). Subtracting 2mS_0 removes the
zero-mode term. Elementary finite sums yield

\[
\sum_{j=1}^{m-1}j(m-j)=\frac{m^3-m}{6}.               \tag{7}
\]

This also follows by induction: the sum vanishes at m=1, and its increment
is m(m+1)/2, the increment of the displayed cubic. Antisymmetry handles
negative m. The result is independent of zero-mode values and so holds
on D with the shared multiplication zero modes.

For K_m=i m a_m, [K_m,K_n]=0 by nullness, whereas

\[
[S_m,K_n]+[K_m,S_n]=i(m^2-n^2)a_{m+n}=(m-n)K_{m+n}.
\]

Thus K changes the noncentral generator but cannot change (7). This proves
(6), including the infinite-mode assertion, without taking a finite-matrix
commutator limit. It is the standard free-current Virasoro calculation
applied to the actual LNC improvement; it is not a new central-charge law.

To locate the obstruction in the inherited normal/spatial generators, set

\[
\mathsf H_m=(2\pi)^{-1}\int e^{2\pi imx}h=L_m^++L_{-m}^-,\quad
\mathsf D_m=(2\pi)^{-1}\int e^{2\pi imx}d=L_m^+-L_{-m}^-.
\]

With F(m)=(m^3-m)/6, direct substitution in (6) gives

\[
[\mathsf D_m,\mathsf D_n]=(m-n)\mathsf D_{m+n},\quad
[\mathsf H_m,\mathsf H_n]=(m-n)\mathsf D_{m+n},
\]
\[
\boxed{[\mathsf D_m,\mathsf H_n]=(m-n)\mathsf H_{m+n}
                    +2F(m)\delta_{m+n,0}I.}          \tag{8}
\]

The central terms cancel in two brackets and add in the mixed one. Equal
left/right central charges therefore do not restore the full constraint
algebra. At m=2,n=-2 the mixed scalar defect is 2.

## 4. Strong-kernel and positivity obstructions

**LQA-T3.** No constant gamma makes all unchanged operators
C_m=L_m+gamma delta_m0 admit a nonzero simultaneous zero vector in a module
on which (6) and the relevant products are defined.

Indeed, on such a vector, [C_m,C_-m] implies
[F(m)-2m gamma]v=0. For m=1 this requires gamma=0; for m=2 it requires
gamma=1/4. Equivalently m=2 and m=3 require 1/4 and 2/3. These cannot
both hold. This includes dual/distributional realizations only if the
same operator identities are well defined there. It is not a statement
about a BRST construction, weak/half constraints or another representation.

There is also a separate kinematic positivity obstruction. In any positive
inner product with the adjoints (4), positive vacuum norm and the chosen
annihilators, for n>0,

\[
\|(a_{-n}-b_{-n})\Omega\|^2
=\langle\Omega|(a_n-b_n)(a_{-n}-b_{-n})|\Omega\rangle
=-2n\|\Omega\|^2.                                  \tag{9}
\]

Thus this vacuum contract forces an indefinite form. In diagonal current
coordinates u=(a+b)/sqrt(2), v=(a-b)/sqrt(2), the current metrics are +1
and -1. Both give +1 to the central charge with this same positive-index
vacuum: the double contraction contains the metric sign twice. Changing
the negative-sign vacuum is a distinct quantization; the known alternative
must be checked with this source's improvement before use. Neither (9)
nor (6) is a no-go theorem for all quantizations of LNC.

## 5. Exact classification of a possible linear correction

Keep the quadratic ordering and module fixed, but examine the explicitly
different operators

\[
L_m(r,s,\gamma)=S_m+i m(r a_m+s b_m)+\gamma\delta_{m,0}I,
\qquad r,s,\gamma\in\mathbb R.                       \tag{10}
\]

They still preserve D and obey L_m^*=L_-m there. The inherited model is
(r,s,gamma)=(1,0,0).

**LQA-T4 (necessary and sufficient within (10)).**

\[
[L_m,L_n]=(m-n)L_{m+n}+\mathcal A(m)\delta_{m+n,0}I,
\quad
\mathcal A(m)=\left(\tfrac16+2rs\right)m^3
                  -\left(\tfrac16+2\gamma\right)m.   \tag{11}
\]

Thus every central coefficient vanishes if and only if

\[
\boxed{rs=-\tfrac1{12},\qquad\gamma=-\tfrac1{12}.}    \tag{12}
\]

**Proof.** Q_m=r a_m+s b_m satisfies [Q_m,Q_n]=2rs m delta.
The cross terms with S again reproduce (m-n)i(m+n)Q_{m+n}.
The improvement commutator is -mn[Q_m,Q_n]=2rs m^3 delta. Replacing
S+imQ by L subtracts 2m gamma delta from the central part. This proves
(11). A polynomial A m^3+B m vanishes for all integer m exactly when
A=B=0 (m=1,2 already suffice). This gives (12).

Fixing r=1 to retain the inherited A derivative makes s=gamma=-1/12
the unique choice *inside this family*. The original null family s=0
cannot work for any r or gamma. There is no asserted uniqueness among
general counterterms, representations or constraint reductions.

Restoring [a_m,b_n]=hbar m delta while leaving S quadratic gives

\[
[L_m,L_n]=\hbar(m-n)L_{m+n}
 +\left[\frac{\hbar^2}{6}(m^3-m)+2\hbar rs m^3
                    -2\hbar\gamma m\right]\delta I.
\]

For hbar>0, cancellation requires rs=gamma=-hbar/12. With r=1 this
is an order-hbar correction preserving the original principal classical
expression; it is still an added quantum prescription. If r,s were held
fixed and non-null as hbar tends to zero, even the classical bracket would
change. A scalar intercept by itself never cancels the cubic term.

## 6. Homogeneous comparison and the next decision

**LQA-T5 (zero-mode compatibility test).** In REC's period-one homogeneous
chart write P=L pi and K for the total scalar momentum, using capital
letters to distinguish them from local densities. Then

\[
a_0=P/2,\quad b_0=K/(2\pi),\quad
4\pi a_0b_0=PK=C_{\rm HQS}.                          \tag{13}
\]

With no oscillators S_0=a_0b_0. The original zero-mode equation therefore
agrees with HQS up to the nonzero multiplier 4 pi. In contrast, (12) gives

\[
4\pi L_0=PK-\pi\hbar/3=0.                          \tag{14}
\]

This is a different zero-mode constraint. A field derivative improvement
alone cannot undo its constant shift. On the circle the corresponding
formal corrected densities, with r=1, are

\[
\widetilde T_\epsilon=T_\epsilon^{\rm ordered}
          +\frac{\epsilon\hbar}{24\pi}B_\epsilon'
          -\frac{\pi\hbar}{6},
\]
\[
\widetilde h=h^{\rm ordered}+\frac{\hbar}{12\pi}q''-\frac{\pi\hbar}{3},
\quad
\widetilde d=d^{\rm ordered}+\frac{\hbar}{12\pi}k'.   \tag{15}
\]

These density expressions mean their finite Fourier modes (10); no
operator-valued pointwise product or smeared self-adjointness is asserted.
Direct substitution of (3) verifies all factors in (13)–(15).

The oscillator vacuum subspace is not invariant under negative constraint
modes. Consequently (13)–(14) are a comparison of zero-mode equations,
not a claim that HQS already embeds as a full local physical subspace.
Quantizing the homogeneous truncation and imposing all local constraints
are distinct operations. No positive physical norm or nonzero all-mode
kernel for (12) has been constructed here; (9) remains a kinematic fact.

The counterterm is not adopted. The following explicit comparison instead
retains the original improvement and the homogeneous zero-mode equation.
HQS and its positive physical sector stay frozen.

## 7. A source-preserving alternative with positive kinematics

Use rational diagonal currents, avoiding any numerical square roots:

\[
u_m=a_m+b_m,\quad v_m=a_m-b_m,\qquad
[u_m,u_n]=2m\delta_{m+n,0},\quad
[v_m,v_n]=-2m\delta_{m+n,0},\quad[u_m,v_n]=0.          \tag{16}
\]

Now choose a different vacuum: u_n and v_-n annihilate it for n>0;
u_-n and v_n create. This explicitly replaces the vacuum of section 2.
It is a known representation choice for an indefinite-sign current pair,
not an identity that changes the previous c=2 result retroactively.

Let U_n,V_n be polynomial creation variables. The realization

\[
u_{-n}=U_n,\ u_n=2n\partial_{U_n},\quad
v_n=V_n,\ v_{-n}=2n\partial_{V_n}\quad(n>0)           \tag{17}
\]

has a positive bosonic Fock inner product: the squared norm of a monomial
is the product of (2n)^j j! for every occupied variable and its occupation
j. Tensor this with L^2(R^2,da_0 db_0), setting u_0=a_0+b_0 and
v_0=a_0-b_0. The algebraic domain D_rev of Schwartz zero-mode functions
times finite oscillator polynomials is dense and invariant under the
operators below. The two chiral copies share only the zero modes.

Normal order separately with the new annihilators and define

\[
S_m^{\rm rev}=\tfrac14\sum_j:u_{m-j}u_j:
                   -\tfrac14\sum_j:v_{m-j}v_j:,
\qquad
L_m^{\rm rev}=S_m^{\rm rev}+\tfrac{im}{2}(u_m+v_m).  \tag{18}
\]

Each mode acts finitely on D_rev. Both-creator u terms have finitely many
partitions when m<0; both-creator v terms have finitely many when m>0;
all remaining nonzero terms contract with an occupied mode. This gives a
common invariant domain for finite products and establishes
<f,L_m g>=<L_-m f,g> on that **positive kinematic** domain. Essential
self-adjointness of real smearings and integration to a gauge group are
not asserted.

**LQA-T6 (anomaly-free current constraints, with an explicit boundary).**
For every integer m,n on D_rev,

\[
\boxed{[L_m^{\rm rev},L_n^{\rm rev}]=(m-n)L_{m+n}^{\rm rev}.} \tag{19}
\]

The two chiral copies commute, so all three H/D brackets now close with
zero central terms. No b-current improvement or scalar intercept is added.
The zero-oscillator zero-mode equation remains PK=0, as in (13).

**Proof.** For a current u of metric 2, the coefficient 1/4 gives a
Sugawara generator with c=1, by the same contraction argument as (7), now
with coefficient (m^3-m)/12. Set w_n=v_-n; then [w_m,w_n]=2m delta and
w_n for n>0 are annihilators. The v part of (18) is minus the ordinary
w Sugawara generator at index -m. Its commutator has the same Witt term
and the opposite central coefficient -(m^3-m)/12. Thus S_rev has no
central term. In both sectors [S_m^rev,u_n]=-n u_{m+n} and likewise for
v. The improvement uses the null combination (u+v)/2=a, so its
self-commutator is zero and its cross terms reproduce the improvement
in (m-n)L_{m+n}. This proves (19). Positivity follows directly from
the monomial norms of (17); adjoints follow by creation/derivative
pairing and the reality of the zero modes.

This choice has an unbounded mode grading:

\[
S_0^{\rm rev}=a_0b_0+\sum_{n>0}n(N_{U,n}-N_{V,n}).   \tag{20}
\]

It is not a positive-energy highest-weight representation. No physical
time or physical energy has been identified with S_0, so (20) is recorded
as a representation property rather than silently ruled out or silently
called stable physical dynamics.

There is **no nonzero strong common kernel inside D_rev itself**. To see
this, take a nonzero polynomial with all occupied indices at most B and
choose m>2B+2. Choose j,m-j>B, distinct. In L_m^rev its term
-V_j V_{m-j}/2 multiplies that polynomial. Neither the u part, the
linear improvement, the zero-mode terms nor a contraction can produce
these two new high-index variables. Its nonzero coefficient cannot cancel.
This argument applies coefficientwise to finite sums with Schwartz
zero-mode coefficients. It does not exclude infinite Fock vectors outside
D_rev or generalized states in a chosen dual; their domains require work.

Thus (19) closes the **local constraint algebra** in a positive kinematic
representation, not the physical quantum theory. A physical state space
requires solving the constraints on a justified extension, constructing a
positive rigging/reduced inner product, defining observables and comparing
the result with HQS. Zero-mode agreement alone is not that equivalence.
This source-preserving representation is the next working candidate;
the counterterm branch stays a documented comparison, not a model change.

## 8. Verification contract and closure verdict

The written proofs above establish both all-integer algebras on their
declared domains and classify (10). The native certificate computes finite oscillator
actions without a fixed truncation: each application includes every term
that can act on the input polynomial, recalculating the bound for products.
It checks central contractions, mixed-bracket signs, source factorization,
indefinite Gram witnesses, wrong repairs, and the counterterm mismatch.
Paninian word certificates check finite Weyl and commuting-polynomial
identities, including the cubic-sum induction step. These are not a machine
verification of the full domain argument or of a continuum field theory.

Closed here: the two selected local current-algebra calculations, a positive
kinematic anomaly-free representation, and the linear-counterterm
classification. Failed here: unchanged strong all-mode closure in the first
vacuum, and a nonzero strong kernel on the second algebraic test domain.
Open: an acceptable local physical quantization, potentials, quantum matter, singular clocks, global cut
gluing, and physical 3+1 gravity. The result identifies an actual obstacle
and its algebraic repair cost; it does not turn certificate success into
success of the rejected physical claim.
