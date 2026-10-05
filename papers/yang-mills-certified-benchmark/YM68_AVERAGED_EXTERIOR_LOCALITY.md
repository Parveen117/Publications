# YM-68: averaged locality for an arbitrarily large or infinite exterior

Monty Dabas. 5 October 2026. **Python 3.12 only.**

YM67 controls every exact finite set of coarse pins, but its pointwise
spatial estimate deteriorates with the number of pins. This continuation
removes that size dependence in an **averaged** norm. The constants below
depend on the fixed observed region and its separation from the pins,
even when the pinned exterior is countably infinite.

The main step is symmetric smoothing in the target coordinates alone.
No exterior density floor is used. The result controls total variation,
mutual information, the native logarithmic derivative in mean square,
and the corresponding finite-target diagonal drift. Infinite-exterior
conditional densities are C1 in the target coordinates for almost every
exterior configuration.
Uniform control at every exterior configuration, mixing after arbitrary
nearby conditioning, and full Schur-memory locality remain open.

## Setup and inherited inputs

Retain the actual YM62/YM66 vacuum for

\[
 A=L_E-\sum_p j_pW_p,\qquad |j_p|\le\eta<1/8000,
 \qquad d\mu=h^2d\Phi.                                      \tag{1}
\]

The coarse links are YM64's disjoint ordered product paths of fixed
length b. Work first in a finite admissible volume, then in YM66's
unique infinite-volume state. Let Y be m>=1 target coarse links, and Z
any disjoint finite collection of coarse links whose **complete fine
paths** have distance at least R>=1 from Y's complete fine paths.
Unobserved links between Y and Z are integrated out. Write

\[
 q(y)=w_Y(y),\qquad p_z(y)=k_{Y|Z}(y|z).
\]

Use B, a1=2B, a2=B+6B^2, ell2=B+10B^2 and K=exp(4pi B) from
YM67. Its pivot argument, including conditioning on all non-pivot
variables before averaging, gives for every finite z

\[
 K^{-m}\le p_z,q\le K^m,\qquad
 \|\nabla_e\log p_z\|,\|\nabla_e\log q\|\le a_1,\qquad
 |D_iD_j p_z|,|D_iD_jq|\le a_2K^m,\qquad
 |D_iD_j\log p_z|,|D_iD_j\log q|\le\ell_2.                  \tag{2}
\]

Directions are native unit directions and keep their order. These
target bounds do not contain |Z|. At infinite volume with finite Z,
second derivatives are understood through uniform finite-volume bounds
or as weak derivatives; first derivatives are continuous.

Let gamma and rho be YM66's common gap and spatial rate. YM65's
equal-time covariance estimate, with its smaller-support factor bounded
by b m, gives for every real bounded target readout f and exterior
readout g, of norms at most one,

\[
 |\operatorname{Cov}(f(Y),g(Z))|\le\delta_R,
 \qquad
 \delta_R=\min\left\{1,
 \left(1+\frac{2bm}{3}+\frac1{3\rho R}\right)e^{-\rho R}\right\}.
                                                               \tag{3}
\]

At eta=0 use delta_R=0 exactly. The complete exterior observable g
may involve every link of Z; (3) is not merely a pairwise correlation
bound. No sum over the size of Z is required.

The distinction from uniform conditional continuity is substantive.
For context, [Rebeschini and van Handel, arXiv:1308.4117](https://arxiv.org/abs/1308.4117)
develop comparison estimates from local conditional specifications;
[Kuelske and Opoku, arXiv:0711.3764](https://arxiv.org/abs/0711.3764)
study conditional continuity under local transformations with Dobrushin
control. No such uniform conditional comparison hypothesis is asserted
here. The estimates below are direct consequences of (2)--(3).

## YM68-T1: target-only smoothing removes the exterior-size factor

For a target point y, let B_r(y) be a product of m native balls of
radius 0<r<=1. YM67's reference-volume bound is

\[
 v_r=\Phi_Y(B_r(y))\ge r^{3m}/40^m.                         \tag{4}
\]

Each one-link ball is invariant under reversal of a radial native
geodesic. For a function whose second derivative along that geodesic
is bounded by H, symmetric Taylor averaging differs from its central
value by at most H r^2/2. Apply this to p_z-q, whose bound from (2)
is H=2 a2 K^m. Average one target link at a time. Each averaging
operator is a sup-norm contraction, so a telescoping sum gives

\[
 \left|p_z(y)-q(y)-\frac1{v_r}\int_{B_r(y)}(p_z-q)\,d\Phi_Y\right|
                 \le m a_2 K^m r^2.                       \tag{5}
\]

Only same-direction second derivatives enter the Taylor remainder;
no interchange of noncommuting derivatives is made.

Fix y and r and set H_y(z)=integral_(B_r(y))(p_z-q)dPhi_Y.
Choose g(z)=sgn H_y(z), a measurable norm-one exterior readout.
The tower property gives

\[
 \mathbb E_Z|H_y(Z)|
 =\operatorname{Cov}(\mathbf1_{B_r(y)}(Y),g(Z))
 \le\delta_R.                                             \tag{6}
\]

Combining (4)--(6), then integrating y against normalized reference,
proves

\[
 \boxed{\mathbb E_Z\|p_Z-q\|_1\le\varepsilon_{m,b}(R),}
                                                               \tag{7}
\]
where
\[
 \varepsilon_{m,b}(R)=\min\left\{2,K^m-K^{-m},
 \inf_{0<r\le1}\left[
 \frac{40^m\delta_R}{r^{3m}}+m a_2K^m r^2\right]\right\}.    \tag{8}
\]

There is no division by the probability or density of the exterior.
The range cap follows pointwise from (2). At eta=0 the bound is zero.
For eta>0, the uncapped stationary radius is

\[
 r_*^{3m+2}=\frac{3m\,40^m\delta_R}{2m a_2K^m}.
\]

Use r=min(1,r_*) to see that (8) is
O(delta_R^(2/(3m+2))). For fixed eta,m,b, its asymptotic exponential
rate is 2rho/(3m+2). The observed region determines this rate;
the number of exterior pins does not. Total variation is half-L1,
so its expectation is at most epsilon/2.

For any threshold t>0, Markov's inequality also gives

\[
 \mu_Z\{z:\|p_z-q\|_1>t\}\le\min\{1,\varepsilon_{m,b}(R)/t\}.
                                                               \tag{9}
\]

The supremum over target test functions is inside the expectation in
(7). This is stronger than controlling each fixed test function
separately. It is still an expectation over the exterior, not a
supremum over exterior configurations.

## YM68-T2: infinite exterior conditionals are C1 almost surely

In the unique infinite-volume coarse vacuum, let Z be any countable
set disjoint from Y, with no separation needed for existence. Choose
finite Z_n increasing to Z. The positive continuous finite-pin densities
p_n(y)=k_(Y|Z_n)(y|Z_n) obey (2) uniformly in n.

For every fixed y, p_n(y) is a bounded martingale in the increasing
sigma-fields generated by Z_n. Indeed integrating p_(n+1) over the
additional conditional coordinates gives p_n; the identity initially
holds almost everywhere in y and extends to every y by continuity.
On a countable dense subset of the compact target space, martingale
convergence holds simultaneously outside one null exterior set.
The uniform first-derivative bound makes p_n equicontinuous, so this
dense-set convergence extends to uniform convergence in y.

The ordered second-derivative bounds in (2) give equicontinuous first
derivatives in any of finitely many native coordinate charts. Thus
the family is relatively compact in C1. Every C1 cluster limit has
the already determined uniform function limit; consequently the whole
sequence converges in C1. Its limit p_Z obeys the same positive density
floor, gradient bounds, and locally Lipschitz first derivatives.
Bounded convergence against continuous target test functions identifies
p_Z dPhi_Y as the regular conditional law given all of Z.

Two different finite exhaustions of Z give the same conditional
probability measure almost surely. Their continuous densities therefore
coincide everywhere in y on a common full-probability exterior set.
This provides an exhaustion-independent **almost-sure C1 version**.
Values on null exterior configurations are not canonically fixed.

If every link of Z has fine-path distance at least R from Y, (7)--(9)
pass to the limit by bounded convergence. Thus the same numerical bound
applies to finite exteriors of any cardinality and to the whole
countably infinite distant exterior.

## YM68-T3: entropy and native score locality without pin-count constants

For positive normalized target densities p,q satisfying (2), put
f=log(p/q) and define relative Fisher energy using the native derivatives,

\[
 J(p\Vert q)=\int p\sum_{e\in Y,a}|D_{e,a}f|^2\,d\Phi_Y.
\]

The log ratio obeys |f|<=2m log K. Nonnegativity of both relative
entropies and the symmetrized entropy identity give

\[
 D(p\Vert q)\le D(p\Vert q)+D(q\Vert p)
 =\int(p-q)\log(p/q)\,d\Phi_Y
 \le2m\log K\,\|p-q\|_1.                                 \tag{10}
\]

Let Delta=sum D_(e,a)^2, the negative-spectrum native Laplacian.
Integration by parts with r=p/q yields the useful exact identity

\[
 \begin{split}
 J(p\Vert q)
 &=\int q\sum D r\,D f\,d\Phi_Y\\
 &=-\int(p-q)\left(\Delta f+
            \sum D\log q\,D f\right)d\Phi_Y.
 \end{split}                                               \tag{11}
\]

The constant 1 in r-1 is permitted because an integrated divergence
vanishes. There are 3m diagonal second derivatives, each bounded by
2ell2. Per-link gradient norms are at most a1 for log q and 2a1
for f. Consequently

\[
 \boxed{J(p\Vert q)\le
 (6m\ell_2+2m a_1^2)\|p-q\|_1.}                           \tag{12}
\]

The independent pointwise cap is J<=4m a1^2. In C1 limits, (11)
can be obtained by passing through the smooth approximants or using
weak second derivatives; (12) passes directly by C1 convergence.

Use p=p_Z and average. If C_F=6m ell2+2m a1^2, then for finite or
countably infinite distant Z,

\[
 \boxed{I(Y:Z)\le\min\{2m\log K,\,2m\log K\,\varepsilon_{m,b}(R)\},}
                                                               \tag{13}
\]
\[
 \boxed{\mathbb E_Z J(p_Z\Vert q)
        \le\min\{4m a_1^2,\,C_F\varepsilon_{m,b}(R)\}.}     \tag{14}
\]

For an infinite exterior, bounded positive conditional densities make
the entropy formula well-defined, and the martingale C1 limit passes
both quantities to the limit. The square root of (14) is a bound on
the score difference in L2 of the actual joint vacuum law. It is not
a bound on its value at every exact pin. The conditional Fisher energy
exceeds t^2 with probability at most the right side of (14) divided
by t^2, capped at one.

For a smooth target function F, define the conditional and marginal
diagonal differential expressions by YM64/YM67's formula

\[
 C_{b,p}F=bL_YF-b\sum D\log p\,DF.
\]

Pointwise Cauchy--Schwarz and (14) give

\[
 \boxed{\mathbb E_{Y,Z}
  |b^{-1}(C_{b,p_Z}-C_{b,q})F(Y)|^2
 \le\|\Gamma_Y(F)\|_\infty
          \min\{4m a_1^2,C_F\varepsilon_{m,b}(R)\}.}        \tag{15}
\]

The clock-normalized coefficient has no extra factor b. The covariance
input still depends on the full fine target size b m. These are
finite-target conditional diagonal generators; hidden near links have
been averaged. No equality with the full dynamically reduced generator
or with the Schur correction is asserted.

## YM68-T4: receding exterior information vanishes and the spatial tail is trivial

Let F_R be the sigma-field of coarse coordinates outside expanding
finite boxes. For each fixed target Y, its fine-path distance from
these coordinates tends to infinity. Equations (7) and (13) show that
the information about Y in the whole receding exterior tends to zero.
In particular, for every bounded target readout F,

\[
 \|\mathbb E[F(Y)\mid F_R]-\mathbb E F(Y)\|_{L^1(\mu)}
              \le\|F\|_\infty\varepsilon_{m,b}(d_R)\to0. \tag{16}
\]

The reverse martingale theorem identifies the limit as the conditional
expectation on the spatial tail T=intersection_R F_R. Thus every
bounded local readout is independent of T. Local cylinder events
generate the whole coarse configuration sigma-field; a monotone-class
argument makes every tail event independent of itself. Its probability
is therefore zero or one. This is tail triviality of the constructed
coarse vacuum law, not uniqueness among all possible algebraic ground
states or a time-Markov statement.

For exteriors at each integer fine distance R, the bounds are summable
in R. Tonelli and (7) give almost-sure L1 convergence p_(F_R)->q.
Uniform target C1 compactness from (2) upgrades it to C1 convergence
for almost every vacuum configuration. This typical-configuration
convergence has no uniform assertion on the exceptional null set.

## YM68-T5: why average control cannot be upgraded by changing a quantifier

Here is an explicit smooth native comparison family. It is not the
interacting Yang--Mills vacuum. On Q x Q^n, let x and z_i denote the
real quaternion coordinates and put

\[
 a(z)=(1+z)/2,\quad
 g_n(Z)=\prod_{i=1}^n a(z_i)-2^{-n},\quad
 w_n(x,Z)=1+\epsilon x g_n(Z),\qquad0<\epsilon<1.           \tag{17}
\]

Both marginals are normalized product reference because E x=0 and
E a=1/2. The density floor is at least 1-epsilon, uniformly in n;
target first and second derivatives have n-independent bounds.
Yet

\[
 \mathbb E_Z\|k_{x|Z}-1\|_1\le2\epsilon\,2^{-n}\to0,
 \qquad
 k_{x|z_1=\cdots=z_n=1}=1+\epsilon(1-2^{-n})x.             \tag{18}
\]

At the target quaternion i, where x=0 and D_1 x=1/2, the exact-pin
score is epsilon(1-2^(-n))/2. It stays away from zero. By continuity
the same obstruction persists in sufficiently small positive-measure
neighbourhoods of those pins. The native moment E a^2=5/16 also gives

\[
 \chi^2(w_n\Vert1)=\frac{\epsilon^2}{4}
       [(5/16)^n-4^{-n}],\qquad
 \mathbb E_Z J(k_{x|Z}\Vert1)
 \le\frac{3\epsilon^2}{16(1-\epsilon)}[(5/16)^n-4^{-n}].   \tag{19}
\]

Thus even exponentially small average information and score energy,
positive densities, and uniform target regularity do not force a
vanishing worst-case pin effect.

A second distinction is necessary for full drift locality. On Q^3,
w(x,y,z)=1+epsilon xyz makes x and z independent after integrating y.
Pinning y=1 changes their density to 1+epsilon xz, with covariance
epsilon/16. Therefore unconditional target/exterior mixing does not
establish mixing after a nearby variable has been pinned. The next
stronger gate is a spatial estimate uniform under this nearby
conditioning, followed by the hidden resolvent and unbounded Schur
coupling. Neither gate follows from (7)--(15).

## Reproducible evidence and scope

For eta=1/16000, a four-link target and b=2, upward-rounded bounds are:

| Fine-path distance R | Expected total variation | Mutual information | Expected relative Fisher energy |
|---:|---:|---:|---:|
| 128 | 0.00334084 | 0.000985180 | 0.000137671 |
| 256 | 0.000199146 | 0.000058727 | 0.0000142529 |
| 512 | 0.000000134436 | 0.000000039644 | 0.0000000096216 |
| 1024 | 6.15886e-14 | 1.81619e-14 | 4.40790e-15 |

Every row applies to any finite or countably infinite exterior at the
stated separation. These are conservative analytic bounds, not measured
correlations or a physical mass. At R=512, the probability that the
conditional L1 change exceeds 0.001 is at most 0.000268871.

The certificate uses exact rational arithmetic and YM67's directed
120-place decay arithmetic. For displayed bounds it chooses among
explicit dyadic radii; these are valid upper bounds, not numerical
vacuum solutions. The unrestricted infimum in (8) proves the analytic
asymptotic rate; the finite grid is used only for finite displayed cases.
It independently checks native Haar moments, the Laplacian sign and
weighted integration-by-parts identity, affine-density Fisher integrals
with geometric-series tail bounds, and both conditioning obstructions.
There are 26 new tests. The dedicated Python 3.12 workflow also replays
YM67's 30 native derivative tests and YM65's 30 spatial tests: 86 total.
The source manifest preserves 211 upstream files.

The arithmetic checks accompany the written proofs and are not formal
verification or peer review. The actual fixed-spacing lattice vacuum
is retained throughout the theorems; comparison densities are labelled.
Worst-case exterior locality, full Schur-memory locality, time-zero
Markov identification, and a nontrivial four-dimensional continuum
Yang--Mills mass gap are not established here.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym68_averaged_exterior_locality.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym68_averaged_exterior_locality.py -v
~~~
