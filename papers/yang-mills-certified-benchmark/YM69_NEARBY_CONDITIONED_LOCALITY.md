# YM-69: exact nearby pins and an infinite distant exterior

Monty Dabas. 5 October 2026. **Python 3.12 only.**

YM68 averages over a distant exterior while integrating out all nearby
unobserved links. This continuation allows a finite nearby set to be
pinned exactly. The spatial estimate is uniform over **every value of
those nearby pins**, and its constants do not depend on the cardinality
of the distant exterior, which may be countably infinite.

The averaging measure must be the exterior law conditioned on the nearby
pins. An exact Bayes identity cancels the small joint denominator against
that measure. Only the finite nearby marginal's positive density floor
remains. The resulting estimate controls conditional mutual information,
native score energy and the conditional target diagonal drift.

The nearby set may also grow slowly with separation. A sufficient regime
is n=o(sqrt(R)); a quantified n~c sqrt(R) regime is given below. This does
not cover an arbitrary filled spatial collar, a supremum over the distant
exterior, or the full Schur-memory kernel.

## Setup and proof inputs

Keep the actual vacuum and independent couplings of YM66,

\[
 A=L_E-\sum_p j_pW_p,\qquad |j_p|\le\eta<1/8000.              \tag{1}
\]

Use the same disjoint product-link paths of length b. Let Y be m>=1
target coarse links and N be n>=0 other coarse links to be pinned.
Let S=Y union N, s=m+n. The word "nearby" distinguishes N from the
distant set; no adjacency within S is required. Every complete fine
path in the distant exterior Z is at link-graph distance at least R>=1
from all complete fine paths in S. In particular, distance from Y alone
does not suffice if some pinned link of N approaches Z.

Initially Z and the volume are finite. The same bounds pass through
YM66's unique vacuum limit and then to a countably infinite Z.
Write a for a configuration of N and z for a configuration of Z. Set

\[
 Q(y,a)=w_S(y,a),\quad Q_N(a)=w_N(a),\quad
 P_z(y,a)=k_{S|Z}(y,a|z),\quad P_{N,z}(a)=\int P_z(y,a)d\Phi_Y.
                                                               \tag{2}
\]

Use YM67's constants B, a1=2B, a2=B+6B^2, ell2=B+10B^2 and
K=exp(4pi B). In particular, Q_N(a)>=K^(-n). The target conditionals
given N alone or N and Z have density between K^(-m) and K^m,
per-link score norm at most a1, and ordered second log derivatives
bounded by ell2. These are the native pivot estimates of YM67, not a
new fitted potential or a representation cutoff.

For context, [Rebeschini and van Handel, arXiv:1308.4117](https://arxiv.org/abs/1308.4117)
develop comparison estimates from local conditional specifications;
[Kuelske and Opoku, arXiv:0711.3764](https://arxiv.org/abs/0711.3764)
treat conditional continuity under transformations with Dobrushin
control. The present proof instead combines YM68's smoothing argument
with the explicit native finite-pin density floor. It makes no claim of
uniformity over the size of an arbitrary nearby collar.

## YM69-T1: joint pointwise control before conditioning nearby values

Let rho be YM66's spatial rate and put

\[
 c_{s,b}(R)=1+\frac{2bs}{3}+\frac1{3\rho R},\qquad
 \delta_{s,b}(R)=\min\{1,c_{s,b}(R)e^{-\rho R}\}.
\]

At eta=0 use delta=0. YM65 bounds the covariance of any norm-one
readout of S and any norm-one readout of Z by delta, using the complete
fine support b s as an upper bound on the smaller support size.

YM68-T1's symmetric-ball argument, **before its final integration over
the target point**, proves

\[
 \boxed{\sup_{(y,a)}\mathbb E_Z
            |P_Z(y,a)-Q(y,a)|\le U_{s,b}(R),}              \tag{3}
\]
where a valid bound is
\[
 U_{s,b}(R)=\min\left\{K^s-K^{-s},
  \inf_{0<r\le1}\left[
    \frac{40^s\delta_{s,b}(R)}{r^{3s}}+s a_2K^s r^2
  \right]\right\}.                                        \tag{4}
\]

Indeed the reference volume of a product of s native balls is at least
r^(3s)/40^s, symmetric averaging has error at most s a2 K^s r^2,
and covariance duality bounds the expected absolute ball integral.
This is uniform in the center (y,a).

The order of operations in (3) matters: the expectation is inside the
pointwise supremum. No expectation of a supremum is claimed. Also the
cap 2 for a normalized L1 difference cannot be used in (4), which is
a pointwise density estimate. A later conditional L1 bound does have
that cap. For fixed s,b, (4) is O(delta^(2/(3s+2))).

## YM69-T2: posterior normalization cancels the exterior denominator

Fix **any** a. Define the target marginal under this exact pin and
the target conditional after additionally observing z by

\[
 q_a(y)=\frac{Q(y,a)}{Q_N(a)},\qquad
 p_{a,z}(y)=\frac{P_z(y,a)}{P_{N,z}(a)}.
\]

If mu_Z is the unconditioned exterior law, Bayes' rule gives the
correct averaging law

\[
 d\nu_a(z)=d\mu_{Z|N=a}(z)
          =\frac{P_{N,z}(a)}{Q_N(a)}d\mu_Z(z).              \tag{5}
\]

Put F_a(y,z)=P_z(y,a)-Q(y,a) and H_a(z)=integral_Y F_a dPhi_Y.
Subtracting the target conditionals yields the exact identity

\[
 p_{a,z}-q_a=\frac{F_a-q_a H_a}{P_{N,z}(a)}.
\]

Multiplying by (5) cancels the denominator P_(N,z)(a), so

\[
 \begin{split}
 \mathbb E_{\nu_a}\|p_{a,Z}-q_a\|_1
 &=\frac1{Q_N(a)}\int\!\int |F_a-q_aH_a|d\Phi_Yd\mu_Z\\
 &\le\frac2{Q_N(a)}\int\!\int |F_a|d\Phi_Yd\mu_Z
 \le 2K^n U_{s,b}(R).                                     \tag{6}
 \end{split}
\]

The factor two uses integral q_a=1 and |H_a|<=integral |F_a|.
Integrating (3) only in y suffices; the near point a remains fixed.
Combining (6) with the target density range proves

\[
 \boxed{\sup_a\mathbb E_{Z|N=a}\|p_{a,Z}-q_a\|_1
 \le\tau_{m,n,b}(R):=
 \min\{2,K^m-K^{-m},2K^n U_{m+n,b}(R)\}.}                  \tag{7}
\]

For n=0, H_a=0 exactly, so the factor two is unnecessary and YM68's
bound is recovered. At eta=0 all the differences vanish exactly.
Expected total variation is at most tau/2. For t>0,

\[
 \sup_a\nu_a\{z:\|p_{a,z}-q_a\|_1>t\}
                    \le\min\{1,\tau_{m,n,b}(R)/t\}.        \tag{8}
\]

For each fixed n, the decay rate is 2rho/(3(m+n)+2). Uniformity is
over the nearby values a, the distant set's cardinality and its values
in the averaged sense of (7). It is not uniform over all n or over
every exact distant configuration z.

## YM69-T3: an infinite distant exterior with a common exceptional set

In the unique coarse vacuum, take a countably infinite Z and exhaust
it by finite Z_j. Apply YM68-T2 with the finite target S=Y union N.
Outside a single mu_Z-null set, the joint conditional densities P_z
exist and are C1 in **all** coordinates (y,a); their density and
derivative bounds are uniform. Set P_(N,z)=integral_Y P_z dPhi_Y.
It is positive and C1 in a on that same exterior full-measure set.

The identities integral_Z P_z(y,a)dmu_Z=Q(y,a) and their N marginals
hold initially almost everywhere, then for every (y,a) by continuity
and uniform bounds. Consequently (5) defines a normalized exterior
probability law for every exact a. Moreover

\[
 K^{-2n}\le\frac{d\nu_a}{d\mu_Z}\le K^{2n}.               \tag{9}
\]

The exceptional set is therefore null for **every** nu_a, with one
common set independent of a. The ratio P_z(y,a)/P_(N,z)(a) defines
the target conditional there and is C1 jointly in (y,a). Its sharper
target-only bounds follow from the finite-pin approximants in YM67.
The limits along finite exterior exhaustions give the same versions
almost surely. In particular, for each finite cylinder readout of Z,
(5) agrees at every a with the continuous finite-dimensional Bayes
formula. These identities determine the exterior probability kernel.

There is no uncountable intersection of separately chosen null sets.
The construction first obtains one C1 joint density on the compact
finite S, then divides its positive N marginal. Bounded convergence
passes (3) and (6)--(8) to this construction. Thus (7) holds for an
infinite distant exterior, still uniformly in every nearby value a.
No particular value of P_z is assigned on the common exceptional set.

## YM69-T4: conditional information, native scores and diagonal drift

Both p_(a,z) and q_a satisfy YM68-T3's target regularity hypotheses
with m links, independently of n and |Z|. Let

\[
 C_{F,m}=6m\ell_2+2m a_1^2,\qquad
 J(p\Vert q)=\int p\sum_{e\in Y,c}|D_{e,c}\log(p/q)|^2d\Phi_Y.
\]

The symmetrized entropy identity and native integration by parts from
YM68 therefore yield, for every a,

\[
 \boxed{I(Y:Z\mid N=a)\le
      \min\{2m\log K,2m\log K\,\tau_{m,n,b}(R)\},}         \tag{10}
\]
\[
 \boxed{\mathbb E_{Z|N=a}J(p_{a,Z}\Vert q_a)
      \le\min\{4m a_1^2,C_{F,m}\tau_{m,n,b}(R)\}.}         \tag{11}
\]

Here (10) means mutual information in the full conditional probability
law at the exact value N=a. Averaging it against the N marginal gives
the usual conditional mutual information I(Y:Z|N) with the same bound.
The measure in (11) is the posterior nu_a, not the unconditioned mu_Z.

For a smooth target function F, define the diagonal expressions
C_(b,p)F=b L_YF-b sum Dlog p DF. Pointwise Cauchy--Schwarz gives

\[
 \boxed{\sup_a\mathbb E_{Y,Z|N=a}
  |b^{-1}(C_{b,p_{a,Z}}-C_{b,q_a})F(Y)|^2
 \le\|\Gamma_Y(F)\|_\infty
       \min\{4m a_1^2,C_{F,m}\tau_{m,n,b}(R)\}.}            \tag{12}
\]

This extends the previous drift estimate to an exact nearby background.
Links outside Y,N,Z are still integrated out. Identifying the result
with a full Schur-reduced dynamical generator would require additional
hidden-resolvent and coupling estimates.

## YM69-T5: controlled growth of the nearby set, and the collar limit

Fix eta>0, b and the target size m. Allow n=n(R) and let s=m+n.
Set A=40^s delta_(s,b)(R) and C=s a2 K^s. If A<=C, the radius
r=(A/C)^(1/(3s+2)) in (4) gives

\[
 U_{s,b}(R)\le2 A^{2/(3s+2)}C^{3s/(3s+2)},\qquad
 \tau_{m,n,b}(R)\le
 4K^n A^{2/(3s+2)}C^{3s/(3s+2)}.                           \tag{13}
\]

With delta<=c_(s,b)(R)exp(-rho R), the logarithm of the last expression
is at most

\[
 \log4+n\log K+
 \frac2{3s+2}\bigl[s\log40+\log c_{s,b}(R)-\rho R\bigr]
 +\frac{3s}{3s+2}\bigl[\log(s a_2)+s\log K\bigr].          \tag{14}
\]

If s^2/R->0, then A<=C eventually and the negative term of order R/s
dominates the positive O(s+log s) terms. Hence tau->0. In particular,
every n=o(sqrt(R)) is sufficient. A sharper boundary regime follows
by taking n/sqrt(R)->c>0:

\[
 \limsup_{R\to\infty}\frac{\log\tau_{m,n,b}(R)}{\sqrt R}
 \le 2c\log K-\frac{2\rho}{3c}<0
 \quad\hbox{if}\quad
 \boxed{c^2<\frac{\rho}{3\log K}.}                         \tag{15}
\]

If tau vanishes exactly, interpret its logarithm as minus infinity.
At eta=0, independence is exact and the separate free case needs no
growth restriction. These are sufficient regimes; failure of (15)
does not demonstrate failure of actual locality.

At eta=1/16000, directed intervals give rho/(3 log K)>3.48724.
For the four-link target with b=2, n=sqrt(R) is therefore admissible.
The proof concerns the number of pinned links and their separation
from Z. It does not assert that a filled three-dimensional collar of
radius R, containing order R^3 links, lies in this regime. Passing
to arbitrary growing collars remains an unresolved step.

## Exact normalization checks and numerical bounds

A positive finite comparison law with weights (1,1,1,1,1,3,3,9) on
binary (Y,N,Z), in lexicographic order, has at N=0

\[
 \mathbb E_{Z|N=0}\|p_{0,Z}-q_0\|_1=2/9,\qquad
 Q_N(0)^{-1}\mathbb E_Z\|F_0(\cdot,Z)\|_1=1/5.
\]

Thus simply omitting the normalization term H_a would invalidate a
factor-one bound. Averaging the conditional difference under the wrong
unconditioned exterior law gives 13/60 instead of 2/9.

An independent smooth native example is
w(x,a,z)=1+alpha a z+epsilon x a z, using real quaternion coordinates,
alpha=3/5 and epsilon=1/10. Its floor is 3/10. At a=1 the posterior
exterior density is 1+alpha z, and the conditional target density is
1+[epsilon z/(1+alpha z)]x. Correct posterior averaging gives target
mean zero; unweighted exterior averaging gives -1/216. The native
moment formula for q0 and a rigorous series tail independently verify
the resolvent integral used in that value. These are comparison laws,
not finite approximations to the interacting YM vacuum.

The certificate also evaluates the theorem's actual inherited constants
with directed arithmetic. To avoid magnifying a tiny absolute rounding
floor, for a dyadic radius r=2^(-k) it combines exponents **before**
evaluating the covariance/ball term:

\[
 40^s c_{s,b}(R)2^{3sk}e^{-\rho R}
 =c_{s,b}(R)\exp[-\rho R+s\log40+3sk\log2].                 \tag{16}
\]

Lower bounds on the positive decay exponent and upper bounds on the
prefactor give upper error bounds. Negative-decay-exponent candidates
may be omitted, leaving the density-range cap. The finite dyadic grid
certifies individual displayed cases; (13)--(15), with unrestricted
radii, establish the asymptotic statements.

For eta=1/16000, m=n=4 and b=2, upward-rounded bounds are:

| Separation R from the entire target-plus-near set | Expected TV, uniformly in near values | Conditional mutual information | Conditional relative Fisher energy |
|---:|---:|---:|---:|
| 256 | 0.0147665 | 0.00435449 | 0.000137671 |
| 512 | 0.000230779 | 0.0000680541 | 0.0000165168 |
| 1024 | 2.10974e-7 | 6.22141e-8 | 1.50994e-8 |
| 2048 | 5.37202e-14 | 1.58416e-14 | 3.84475e-15 |

For the same eta,b and fixed m=4, the growing-near regime has:

| Nearby pin count n | Separation R=n^2 | Uniform-near expected L1 bound |
|---:|---:|---:|
| 128 | 16,384 | 0.000179606 |
| 256 | 65,536 | 2.36105e-9 |
| 512 | 262,144 | 2.08713e-19 |

These conservative analytic budgets apply to finite or countably infinite
distant exteriors. They are not measured correlations or physical masses.

Worst-case distant configurations, arbitrary growing nearby collars,
full Schur-memory locality, time-zero Markov identification and a
nontrivial four-dimensional continuum mass gap remain open. The exact
and interval tests support the calculations; they are not formal proof
verification or independent peer review.

The dedicated workflow runs 25 new tests and replays the unchanged YM68
and YM67 certificates with their 26 and 30 tests: 81 targeted tests in
Python 3.12. The source manifest preserves 218 upstream files.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym69_nearby_conditioned_locality.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym69_nearby_conditioned_locality.py -v
~~~
