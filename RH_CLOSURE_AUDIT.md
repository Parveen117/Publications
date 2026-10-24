# CLA-1 — Closure audit of the shunya-potential route

17 October 2026. Continuation at Publications commit
550ed8e6a61d8680bbb1ddd062befaa458c0bea5, following
[SPT-1](RH_SHUNYA_POTENTIAL.md).

**Decision: the current argument does not prove RH.** The native cut,
the complete polynomial invariant algebra, potential-preserving
exponential transport, and the repeated-cut limit do not by themselves
imply positivity of the arithmetic source. This note proves that failure
of implication with an explicit faithful comparison source. It also
proves that an arbitrarily long finite moment record can agree with a
positive source while a higher-degree negative direction remains.

The comparison sources are not the Riemann source and are not
counterexamples to RH. They identify precisely which source-independent
arguments cannot close the remaining arithmetic step. The actual source
and all previously proved support-window bounds are retained.

For the exponential deformation, a negative direction, if present,
survives at every finite cut count and every finite deformation parameter,
with an explicit amplitude choice. Thus taking a limit that erases its
unscaled reading cannot establish the missing sign.

These are written algebraic and analytic proofs. No source evaluation,
numerical experiment, scientific replay, independent agent or workflow
dispatch is performed. The elementary methods used here are not claimed
as a new general theory.

## Retained source and the precise target

Fix one finite anchor a>=3249/4. Retain CEM's actual source
\[
\begin{gathered}
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\\
\mu=(\gamma-i\delta)^2,\quad \delta=\beta-\tfrac12,\qquad
u_\mu=\frac a{a+\mu},\qquad
w_\mu=\frac{\nu_\mu}{q(a)(a+\mu)},\\
\mathcal L_a[f]=\sum_\mu w_\mu f(u_\mu),\qquad
\chi(y)=y(1-y),\\
B_a(f,h)=\mathcal L_a[fh],\qquad
M_a(f,h)=\mathcal L_a[\chi fh].
\end{gathered}
\tag{1}
\]
The sum is over distinct source values with their multiplicities.
The inherited bounds are
\[
\mathcal L_a[1]=1,\quad \sum_\mu |w_\mu|\le2,\quad
0<|u_\mu|<1,\quad u_\mu\ne1.
\tag{2}
\]
The support is infinite, conjugation symmetric, and has only 0 as an
accumulation point. Products in (1) are bilinear even at nonreal points.

In the retained/complementary coordinates, the existing objects are
\[
\Pi=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
D_y=\begin{pmatrix}1&0\\0&\chi\end{pmatrix},\qquad
A_0=\begin{pmatrix}0&2\chi\\-2&0\end{pmatrix},\qquad
\phi_y(r,h)=\frac12(r^2+\chi h^2).
\tag{3}
\]
They obey
\[
\Pi A_0\Pi=0,\qquad A_0^TD_y+D_yA_0=0,\qquad
A_0^2=-4\chi I.
\tag{4}
\]
SPT proves that phi_y generates every pointwise polynomial invariant
of this specified flow. Conservation in (4) is not an inequality.

## CLA-T1. The exact closure requirement still has an all-profile quantifier

Put c_j(a)=L_a[y^j]. The retained criterion is
\[
\boxed{
\mathrm{RH}
\Longleftrightarrow M_a(p,p)\ge0\quad\hbox{for every }p\in\mathbb R[y]
\Longleftrightarrow
D_N(a):=(c_{i+j+1}-c_{i+j+2})_{0\le i,j<N}\succeq0
\quad\hbox{for every }N\ge1 .
}
\tag{5}
\]
This is CEM-T5, not a new proof of its right-hand side.

**Proof and dependency check.** Under RH, all actual u_mu lie in (0,1)
and all w_mu are positive. Each source summand of M_a(p,p) is then
nonnegative. Conversely, if RH fails, the actual support contains a
nonreal conjugate pair z,conjugate(z). CEM's finite detector annihilates
the finitely many other support points with modulus at least |z|.
Multiplication by successive powers of y suppresses the remaining
support. The two-observer Gram determinant has leading term
\[
-4|A_m|^2(\operatorname{Im}z)^2<0,
\]
and an error of strictly smaller exponential order. Some finite
polynomial bank is therefore indefinite. All modified weights
w_mu chi(u_mu) are nonzero, so this detector applies to M_a.
Expanding any polynomial in monomials gives the matrix equivalence.
QED.

At a fixed cut count n, CEM already identifies the full matrix of
rescaled second derivatives with D_N(a). Increasing n does not
increase N or prove any missing sign. An invariant first derivative,
a constant unscaled limit, or a positive finite bank cannot replace
the two universal quantifiers in (5).

The separate waveform formulation remains
\[
\mathrm{RH}
\Longleftrightarrow
\mathscr Q_\xi[g]\ge0
\quad\hbox{for every }g\in\mathcal W_L
\hbox{ and every finite }L>0,
\tag{6}
\]
where the fixed family is
\[
\mathcal W_L=
\left\{g\in C_c^\infty((-L,L),\mathbb R):
g(-x)=g(x),\ \int g=0,\ \int g(x)\cosh(x/2)\,dx=0\right\}.
\]
FPW/TPW prove this equivalence. The actual established TPW bound is
\[
0<L\le11/20
\quad\Longrightarrow\quad
\mathscr Q_\xi[g]\ge\|g\|_2^2/500
\quad(g\in\mathcal W_L).
\tag{7}
\]
It does not assert (6) at all widths.

## CLA-T2. A faithful source can preserve the entire native flow and still have negative memory

Here is an explicit comparison model. It is deliberately given a
different functional symbol and is not substituted into (1).

Set
\[
z_\pm=\frac{1\pm i}{2},\qquad
x_j=2^{-j-2}\quad(j\ge1),\qquad
C=z_++z_-+\sum_{j\ge1}x_j=\frac54 ,
\]
and define
\[
\boxed{
\mathcal L_*[f]=\frac1C
\left(z_+f(z_+)+z_-f(z_-)+\sum_{j\ge1}x_jf(x_j)\right).
}
\tag{8}
\]
This source has all the following properties:

- It is real on real-coefficient polynomials and has L_*[1]=1.
- Its infinite support has only 0 as an accumulation point.
- Every source value is nonzero, different from 1, and inside |y|<1.
- Every weight is nonzero; weights at conjugate points are conjugate.
- Its absolute weight sum is (4 sqrt(2)+1)/5<2.
- Both B_*(f,h)=L_*[fh] and M_*(f,h)=L_*[chi fh] have zero radical.
- All the pointwise identities (3)--(4), the full potential conservation,
  and the repeated-cut freezing limit hold for this same source.
- Nevertheless, for p(y)=y-1/2,
\[
\boxed{M_*(p,p)<-\frac7{80}<0.}
\tag{9}
\]

**Proof.** The geometric series sum is 1/4, proving normalization
and the weight budget. Conjugation, support and convergence properties
follow directly from (8). In particular the usual compact-source
bounds used by CEM and SPT remain valid.

For faithfulness, let
\[
Q(y)=(y-z_+)(y-z_-)=y^2-y+\tfrac12.
\]
It vanishes at the two nonreal nodes and is strictly positive on the
real axis. For every nonzero real polynomial p,
\[
\begin{aligned}
B_*(p,Qp)&=\frac1C\sum_{j\ge1}x_jQ(x_j)p(x_j)^2>0,\\
M_*(p,Qp)&=\frac1C\sum_{j\ge1}x_j\chi(x_j)Q(x_j)p(x_j)^2>0 .
\end{aligned}
\tag{10}
\]
The strict inequalities follow from the infinite distinct real support.
Thus no nonzero p belongs to either radical. Faithfulness of a
bilinear form does not mean it is positive.

For the negative direction, chi(z_+)=chi(z_-)=1/2 and
p(z_+)^2=p(z_-)^2=-1/4. The pair contributes -1/8 before
division by C. At each real node,
0<chi(x_j)<1/4 and p(x_j)^2<1/4. Consequently
\[
\begin{aligned}
M_*(p,p)
&=\frac1C\left(-\frac18+
\sum_{j\ge1}x_j\chi(x_j)(x_j-\tfrac12)^2\right)\\
&<\frac1C\left(-\frac18+\frac1{16}\sum_{j\ge1}x_j\right)
=\frac45\left(-\frac18+\frac1{64}\right)
=-\frac7{80}.
\end{aligned}
\]
This proves (9).

The invariant algebra and generator identities are polynomial identities
before any source evaluation. They therefore hold unchanged for (8).
The support and absolute weight budget justify evaluating the entire
flow and taking the locally uniform repeated-cut limit. In particular
\[
[\Pi e^{tA_0/n}\Pi]^n\longrightarrow\Pi
\]
also holds here, while
\[
\mathscr S_*[0,p]=\tfrac12M_*(p,p)<0 .
\]
Thus invariant potential, faithful readings and cut freezing coexist
with a negative source direction. QED.

The weights in this particular model even have w/u=1/C at every node.
That extra algebraic resemblance still does not identify its nodes
with actual zeta zeros. The model is not asserted to satisfy the
zeta functional equation, the prime explicit formula, its sharper
zero-location constraints, or the actual-source estimates in (7).

The implication refuted is precisely:
"the listed cut/invariance/support identities alone imply positivity."
No impossibility theorem for RH or for future arithmetic arguments is
claimed.

## CLA-T3. Every prescribed finite moment depth admits a positive and an indefinite comparison source with identical readings

For every integer d>=0 there exist two normalized functionals L_0
and L_epsilon with the following properties:

1. Both have infinite conjugation-symmetric support inside 0<|y|<1,
   only 0 as an accumulation point, and absolute weight sum 1.
2. All their weights are positive real numbers. L_0 has only real
   support; L_epsilon also has a nonreal conjugate pair.
3. Their moments agree exactly through degree d:
\[
\mathcal L_\epsilon[y^k]=\mathcal L_0[y^k]
\quad(0\le k\le d).
\tag{11}
\]
4. B_0 and M_0 are positive definite on real polynomials.
   B_epsilon and M_epsilon have zero radical, but M_epsilon
   has a negative polynomial direction.

**Proof.** Use the positive baseline
\[
t_j=2^{-j},\qquad \alpha_j=2^{-j},\qquad
\mathcal L_0[f]=\sum_{j\ge1}\alpha_j f(t_j).
\]
Its weights sum to one. Both unweighted and chi-weighted polynomial
squares are strictly positive for nonzero polynomials.

Fix z=(1+i)/2. Select the first d+1 real nodes. Their Vandermonde
matrix is invertible, so there are unique real numbers b_j satisfying
\[
\sum_{j=1}^{d+1}b_jt_j^k=-\operatorname{Re}(z^k),
\qquad 0\le k\le d.
\tag{12}
\]
Define
\[
\mathcal L_\epsilon[f]=\mathcal L_0[f]
+\epsilon\left(\frac{f(z)+f(\bar z)}2
+\sum_{j=1}^{d+1}b_jf(t_j)\right).
\tag{13}
\]
Choose epsilon>0 small enough that every modified real weight
alpha_j+epsilon b_j is positive. Only finitely many are modified,
so such a choice exists. Equation (12) gives (11), including total
weight one. All weights in (13), including epsilon/2 at each nonreal
node, are now positive real numbers. Positive coefficients at complex
evaluation points do not make a bilinear polynomial-square form positive.

Faithfulness follows exactly as in (10), using
Q(y)=(y-z)(y-bar(z)) to kill the nonreal pair.

To display the negative direction without assuming a general moment
theorem, choose real a_m,b_m so that
\[
a_m z+b_m=i\sqrt2\,z^{-m},\qquad
q_m(y)=y^m(a_my+b_m).
\tag{14}
\]
Real linear interpolation at a nonreal point is invertible, so these
coefficients exist. Their absolute values are bounded by a constant
times |z|^{-m}, with that constant independent of m. Thus
q_m(z)=i sqrt(2), q_m(bar(z))=-i sqrt(2), and the two nonreal
nodes contribute exactly -epsilon to M_epsilon(q_m,q_m).

For every real node t_j<=1/2,
\[
|q_m(t_j)|\le C_1(1/\sqrt2)^m
\]
with C_1 independent of j,m. The real weights have total 1-epsilon
and chi(t_j)<=1/4. Their nonnegative contribution to M_epsilon
therefore tends to zero as m increases. For all sufficiently large m,
\[
M_\epsilon(q_m,q_m)<0.
\]
This completes the construction. QED.

Taking d=2N makes the entire memory bank D_N identical for the two
sources, because its entries use only moments through degree 2N.
That bank is positive definite for both, although the second source
has a negative direction at some larger degree. Any finite collection
of polynomial source readings can similarly be included by taking d
at least their maximum degree.

The quantifier is "for every finite d there exists a comparison source."
The comparison source can change with d. This does not assert that a
single indefinite source passes every degree: CLA-T1's detector would
exclude that. Nor does the construction assert that an alternative
source with the same finite readings is an actual zeta source.

For the linear native flow with polynomial initial data, finitely many
time derivatives are polynomial source readings, so the same limitation
applies to a prescribed finite list of those derivatives. Full analytic
functions, all derivatives simultaneously, and the separately defined
TVSP covariant Lambda tower are not silently identified with such a
finite list.

## CLA-T4. A negative direction survives every exponential rate and every finite cut count

For any conjugation-symmetric discrete source satisfying (2), with
conjugate weights at conjugate points, set
\[
\mathscr N_\eta[r]=
\mathcal L[\chi r^2 e^{\eta r^2}],\qquad \eta\in\mathbb R.
\tag{15}
\]
This is SPT's self-response for
F_eta(s)=(exp(eta s)-1)/eta, with F_0(s)=s.

Suppose some real polynomial p satisfies M(p,p)=-d<0.
Then there exists epsilon_0>0 such that the explicit profiles
\[
\boxed{
r_\eta=\frac{\epsilon_0}{\sqrt{1+|\eta|}}\,p
}
\tag{16}
\]
satisfy
\[
\boxed{
\mathscr N_\eta[r_\eta]\le
-\frac{d\,\epsilon_0^2}{2(1+|\eta|)}<0
\quad\hbox{for every finite real }\eta .
}
\tag{17}
\]
The same inequality describes the exact rescaled self-response at
every finite cut count n>=1.

**Proof.** Let K=max_{|y|<=1}|p(y)|>0. On the closed unit disk,
|chi|<=2 and the absolute weight sum is at most 2. The elementary
bound |exp(z)-1|<=|z| exp(|z|) gives, for real epsilon!=0,
\[
\left|
\frac{\mathscr N_\eta[\epsilon p]}{\epsilon^2}
-M(p,p)\right|
\le4|\eta|\epsilon^2K^4
\exp(|\eta|\epsilon^2K^2).
\tag{18}
\]
Choose epsilon_0 small enough that
\[
4\epsilon_0^2K^4\exp(\epsilon_0^2K^2)<d/2.
\]
Substitute epsilon=epsilon_0/sqrt(1+|eta|) in (18).
Its right-hand side is strictly less than d/2, uniformly in eta.
This proves (17).

For SPT's repeated nonlinear cuts the exact identity is
\[
-\frac n4\partial_t^2R_{n,F_\eta}(0;r,r)=\mathscr N_\eta[r].
\tag{19}
\]
Apply (19) with r=r_eta. No n-dependent error or limit is involved.
QED.

For the actual source, this is conditional on a negative direction
existing; it does not produce an actual negative witness. For the
comparison source (8), the explicit p in (9) supplies one.

If n increases without rescaling, the second response in (19) tends
to zero even when its rescaled value is negative. If |eta| increases,
the witness amplitude in (16) shrinks. A zero unscaled endpoint can
therefore erase sign information. It cannot prove the all-amplitude
inequality. This is a direct quantifier issue, not an invocation of
incompleteness or a prohibition on every possible limiting proof.

The symbols n and eta retain these meanings. They are not identified
with every native lambda coordinate, and no literal equation between
lambda=0 and lambda=infinity is introduced.

## CLA-T5. Changing a regular frame or imposing Maxwell exactness does not select the missing sign

Two elementary facts delimit the geometric argument.

First, for a real vector space carrying a quadratic form Q and any
bijective linear change of frame T,
\[
Q_T(v):=Q(Tv)\ge0\ \hbox{for every }v
\quad\Longleftrightarrow\quad
Q(w)\ge0\ \hbox{for every }w.
\tag{20}
\]
A negative direction w is carried to T^{-1}w. This proves the
equivalence directly, including each finite invertible congruence
of a source bank. The same statement holds on any full profile space
on which the declared transformation is actually bijective.
A singular cut can discard a negative direction, but positivity
on the resulting restricted space does not imply positivity on the
original space.

Second, Maxwell relations follow from exact mixed derivatives and
do not by themselves imply thermodynamic stability. For example,
\[
U(S,V)=\tfrac12(S^2-V^2),\qquad T=U_S=S,\qquad P=-U_V=V
\tag{21}
\]
has T_V=-P_S=0 and regular Legendre coordinate changes. All four
potential differentials and their Maxwell relations are valid.
Nevertheless its Hessian is diag(1,-1), which is indefinite.

The stable real-fiber adapter in SPT-T7 remains valid under its stated
hypotheses y in (0,1) and eta>=0. That choice makes a real thermodynamic
test Hessian positive. It is not a derivation that the actual u_mu in
(1) are real. Assuming that identification would insert the RH content
needed by (5).

Equations (20)--(21) settle the two stated claims. They do not assert
that contact and Hessian structures have identical transformation
rules, or that the arithmetic source transport is the TVSP tangent
transport. The distinctions already recorded by SPT remain in force.

## Closure decision and the remaining arithmetic obligation

The following parts of this route have written closure within their
specified assumptions:

| Statement | Status |
|---|---|
| Classification of the native polynomial invariant algebra | Proved in SPT |
| Potential-preserving scalar Hamiltonian flow and exponential family | Proved in SPT |
| Exact repeated-cut response and its retained sign condition | Proved in CEM/SPT |
| Implication from these source-independent identities to positivity | Refuted by CLA-T2 |
| Inference of all-degree positivity from a prescribed finite moment prefix alone | Refuted by CLA-T3 in the stated comparison class |
| Removal of a negative direction by the specified exponential rates or cut count | Excluded by CLA-T4, with an amplitude-adjusted witness |
| RH positivity for the actual source at every degree or every width | Unproved |

A proof of RH in this program must establish, for the actual functional
(1), either the all-degree inequality (5) or an equivalent arithmetic
statement such as (6). In TPW's closed-energy notation the latter is
\[
\boxed{\|K_L^\circ\|\le1\qquad\hbox{for every finite }L>0.}
\tag{22}
\]
The currently proved contraction bound is
\[
\|K_L^\circ\|\le2000/2001
\qquad(0<L\le11/20).
\tag{23}
\]
The earlier stronger BPC margin at L<=1/2 and the unrestricted
PTE window at L<=1/32 remain unchanged. No extension of a positivity
window or finite-bank range is claimed here.

"An equivalent inequality remains" should not be mistaken for a routine
last technicality: proving that universal inequality is RH itself in
this formulation. New input must exploit the actual arithmetic source
beyond the identities retained by the comparison models.

The closure of this audit is therefore a definite limitation of the
current argument, not a proof or disproof of RH. It also does not prove
that RH is undecidable, that a cut explains a physical speed, or that
the source model establishes a black-hole information law.

## References and verification scope

- [SPT-1](RH_SHUNYA_POTENTIAL.md): invariant classification, exponential
  potential, nonlinear cuts and the TVSP source distinction.
- [CEM-1](RH_CUT_EQUILIBRIUM_MEMORY.md): actual source, memory matrices,
  finite detector, bilinear conservation and repeated-cut limit.
- [OQD-1](RH_OBSERVER_QUOTIENT_DESCENT.md): faithful observer pairing.
- [FPW-1](RH_FIRST_PRIME_WINDOW.md): fixed neutral-family RH equivalence.
- [TPW-1](RH_TWO_PRIME_WINDOW.md): the retained two-prime window,
  compact norm and remaining all-width condition.
- [Clay Mathematics Institute: Riemann Hypothesis](https://www.claymath.org/millennium/riemann-hypothesis/):
  official problem statement and unsolved status, checked 17 October 2026.

Vandermonde interpolation, preservation of quadratic-form sign under
invertible changes of variables, and elementary exponential estimates
are standard tools. Their explicit application above is fully derived;
no general novelty claim is made.

Verification concerns the written proofs, source definitions, document
structure, links, exact immutable-head content and draft status. It is
not a proof-assistant certificate. No scientific computation, numerical
source evaluation, independent agent, replay or CI run is added.
