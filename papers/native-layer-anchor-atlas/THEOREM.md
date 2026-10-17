# LAA-1 — Source-anchor atlas with retained remainder

17 October 2026. Continue LAT-1 on the same completed source
Xcal(w^2)=xi(1/2+w), q=Xcal'/Xcal. No predecessor packet or engine
is changed. The new results are an exact differentiated tail budget,
a source-preserving finite-jet ledger, and overlapping actual-source
positivity charts. A transported Taylor polynomial, an actual source
jet, and a rebuilt finite positive rational model remain distinct.

At a real anchor a=w^2, q(a)>0, write

\[
c_n(a)=\frac{(-a)^n q^{(n)}(a)}{n!q(a)},\qquad
f_a(\lambda)=\frac{q(a(1+\lambda))}{q(a)}.
\]

Thus f_a(lambda)=sum_n (-1)^n c_n(a) lambda^n near zero.
Here lambda is displacement from an anchor, n is derivative depth,
and log(a/a_0) is the LAT anchor parameter. None is physical time.

## LAA-T1. General source anchor and exact differentiated tail

Let w>4 and R=15/16. The retained unconditional SE49 bound gives

\[
|q(a(1+\lambda))|\le Q_w:=\frac{41}{2w}+\frac14
\quad (|\lambda|\le R),\qquad
|c_n(a)|\le M R^{-n},\quad M=Q_w/q(a).
\tag{1}
\]

For integers 0<=j<=K and 0<=z<1 define

\[
\boxed{S_{K,j}(z)=\sum_{\ell=0}^{j}
\binom{K+1}{j-\ell}
\frac{z^{K+1-j+\ell}}{(1-z)^{\ell+1}}
=\sum_{n=K+1}^{\infty}\binom nj z^{n-j}.}
\tag{2}
\]

This is a finite positive rational expression, including at z=0.
It is valid for every retained derivative, without LAT's additional
first-term ratio restriction. Wherever that older restriction holds,
S_(K,j) is no larger than LAT's T_(K,j).

**Proof.** The principal square root obeys
Re(w sqrt(1+lambda))>=w sqrt(1-R)=w/4>1. SE49 states
|w' q(w'^2)|<41/8+|w'|/4 in this half-plane. Divide by
|w'|>=w/4 to get Q_w. The response is analytic on a neighborhood
of the closed disk; Cauchy's coefficient estimate gives (1).
For (2), differentiate z^(K+1)/(1-z) exactly j times and divide
by j!. Leibniz's rule gives the displayed positive sum; termwise
differentiation of the convergent power series gives its infinite
tail interpretation. The ratio of successive terms on the right is
z(n+1)/(n+1-j), decreasing with n. Bounding it by its first value
gives exactly the former LAT bound. QED.

The execution uses the certified lower endpoint of q(a) in M. Thus
all tail bounds are outward, even when the source jet is an interval.
Cauchy's estimate is classical, cited in SOURCE_MAP.md; the retained
SE49 source domain supplies the needed analytic hypothesis.

## LAA-T2. Finite polynomial transport and the source error ledger

Let p_(a,K)(lambda)=sum_(n=0)^K (-1)^n c_n(a) lambda^n.
For rho>0 set delta=rho-1 and define

\[
(P_{\rho,K}c)_j=
\rho^j\sum_{n=j}^K\binom nj c_n(1-\rho)^{n-j},
\qquad 0\le j\le K.
\tag{3}
\]

These are coefficients of p_(a,K)(delta+rho lambda), still
normalized by the ORIGINAL q(a). They satisfy

\[
\boxed{P_{\sigma,K}P_{\rho,K}=P_{\rho\sigma,K},\qquad
P_{\rho,K}^{-1}=P_{1/\rho,K}.}
\tag{4}
\]

The true unnormalized coefficients at b=rho a are
b_j=rho^j (-1)^j f_a^(j)(delta)/j!. Whenever d=|rho-1|<R,

\[
\boxed{|b_j-(P_{\rho,K}c)_j|\le E^{\rm src}_j(d):=
(1+d)^j M R^{-j} S_{K,j}(d/R).}
\tag{5}
\]

If the resulting interval for b_0 is positive, then
q(b)=q(a)b_0 and c_j(b)=b_j/b_0 have certified outward enclosures.
In particular c_0(b)=1 exactly; it is not evaluated by interval
self-division.

**Proof.** Substitution of the affine variable in the polynomial
gives (3). Two substitutions obey
(rho-1)+rho[(sigma-1)+sigma lambda]=rho sigma-1+rho sigma lambda.
This proves (4) for every polynomial of degree at most K. Applying
(1) term by term to the differentiated missing tail and then (2)
gives (5), since rho^j<=(1+d)^j. Positive interval division proves
the normalized assertions. QED.

For a path rho_1,...,rho_h the ledger retains the original source
jet, its original q(a), degree K, and the product rho=rho_1...rho_h.
It evaluates (3) and (5) from that origin. Thus two paths with the
same final ratio have the same polynomial and the same error budget.
The implementation checks every cumulative prefix is inside the
original disk because it treats the ledger as a sequence of valid
source readings. This is stronger than needed for the final reading
alone. Small individual steps cannot be used to evade the cumulative
disk condition. An inverse path can restore the original reading;
the ledger does not manufacture a growing error merely by rounding
successive coordinate changes.

All K+1 retained layers now admit bounds after a valid move, but the
high-derivative bounds can be very wide. No layer K+1 is created.
If an intermediate polynomial is truncated to degree J<K, its
discarded terms must also be carried: at a subsequent ratio sigma,
their exact contribution to coefficient j<=J is
sigma^j sum_(n=J+1)^K binom(n,j) d_n (1-sigma)^(n-j), where d_n
are the intermediate unnormalized polynomial coefficients. Omitting
that contribution destroys the full-degree identity (4).

### Why a shifted polynomial cannot be relabelled as a fresh source

Take f(lambda)=1/(1+lambda/2), K=3 and rho=5/4. Its exact polynomial
is p_3=1-lambda/2+lambda^2/4-lambda^3/8, and

\[
f(\lambda)-p_3(\lambda)
=\frac{(-\lambda/2)^4}{1+\lambda/2}.
\tag{6}
\]

At delta=1/4 this is positive. Both the shifted value and the
normalized first derivative differ from those of the source. Native
rational checks exhibit both differences. Normalizing the polynomial
does have its own projective semigroup wherever the denominator is
nonzero; the error is identifying that polynomial with f. A new
complete source evaluation is a legitimate restart. A general
restart could also carry its inherited error explicitly, but this
packet does not substitute forecast boxes for fresh exact source jets.

## LAA-T3. General-anchor finite positivity corridor

Suppose K=2N-1 and the actual anchor forms A_N,B_N are positive,
with H_N=A_N+B_N, in the ILC notation
A_m=(c_(i+j+1)), B_m=(c_(i+j)-c_(i+j+1)). Their positive finite
model matches c_0,...,c_K and has all later moments in [0,1]. For
0<d<R put

\[
E_j(d)=(1+d)^j\left[M R^{-j}S_{K,j}(d/R)+S_{K,j}(d)\right].
\tag{7}
\]

This bounds the difference between the true and finite-model
unnormalized moved moments on every rho in [1-d,1+d]. For m<=N set

\[
\begin{split}
\epsilon_A&=\max_{1\le j\le2m-1}E_j,&
\epsilon_B&=\max_{0\le j\le2m-2}(E_j+E_{j+1}),\\
\tau_A&=m\epsilon_A\frac{(1+d)^{2m}}{(1-d)^{2m-1}},&
\tau_B&=m\epsilon_B\frac{(1+d)^{2m}}{(1-d)^{2m-2}}.
\end{split}
\tag{8}
\]

If A_m-tau_A I>0, B_m-tau_B I>0, and
1/(1+d)-E_0>0, then both actual m-dimensional forms are strictly
positive for EVERY anchor b in [a(1-d),a(1+d)].

**Proof.** The missing actual/model coefficient differences have
absolute value at most M R^(-n)+1. Equations (2) and (5) give (7).
LAT-T4's polynomial pullback and inverse l1-norm bound then apply
without change; inserting these smaller error bounds gives (8).
An entrywise form error eps costs at most m eps times the squared
Euclidean coefficient norm. The displayed shifted positive forms
dominate that error uniformly. The positive finite-model constant
is at least 1/(1+d), so the final condition ensures positive actual
normalization. Division gives the stated source forms. QED.

No source values are inferred by a sampling interpolation. Each
chart certifies a full interval using its analytic majorant. Failure
of these sufficient inequalities remains inconclusive.

## LAA-T4. Overlapping charts and the retained cut

For one actual source, let finitely many certified order-m charts
have anchor intervals I_i. If their union is connected, the two
actual m-dimensional forms are positive throughout that union.
On overlaps they are the same forms: each is defined by the same
q and derivatives at the same actual anchor. Different outward
enclosures need only contain those same exact values.

**Proof.** Every point of the union belongs to at least one certified
chart, so LAA-T3 gives positivity there. Source differentiation and
normalization identify the exact moments on overlaps. QED.

Consequently the inherited ILC construction is available at every
anchor in the union: Y_m=H_m^(-1)A_m is a strict positive contraction
in its acquired H_m metric, its positive square root
S_m=[Y_m(I-Y_m)]^(1/2) defines the 2m-channel projection, and two
reflections define the corresponding metric-unitary cut wave.
These canonical same-anchor constructions agree on overlaps.
This is a family of finite constructions, not a newly supplied
physical time propagation between their different metrics.

The local rational models fitted independently at different anchors
are generally different rational functions. Overlap of the actual
source charts does NOT identify them with the LAT transport of one
fixed finite model. The executed held-out response comparisons below
certify that distinction for the present source.

## LAA-T5. Executed finite atlas and remaining obligation

Reuse the pinned LAT source jet through order 19 at original s=29,
and evaluate the complete source afresh through order 19 at s=31
and s=33. At all three anchors both 10x10 forms are positive and all
210 available binomial readings are positive. The source label s=33 is
composite; primality is not a hypothesis of this derivative-source
criterion. Each jet retains all Gamma, Euler and derivative tails.

At EACH anchor the certified radii are d=3/20 for m=3 and d=1/12
for m=4. The exact connected unions in the common a-coordinate are

\[
\boxed{m=3:\quad a\in[55233/80,\ 19435/16],}
\qquad
\boxed{m=4:\quad a\in[11913/16,\ 54925/48].}
\tag{9}
\]

Because a=(s-1/2)^2 is increasing here, these contain respectively
every real s in [27,35] and [28,34]. Exact interval endpoints and
overlaps are recorded in CERTIFICATE.json. This supplies six-channel
cuts throughout the first interval and eight-channel cuts throughout
the second, using the inherited algebraic theorem. Ten-dimensional
forms are only certified at the three printed anchors.

The full-degree source ledger for 29->31->33 equals the direct
29->33 ledger. Its enclosures contain the independently evaluated
s=33 moments. A legitimate fresh restart at s=31 gives smaller
remainder budgets for that target. A separate complete source value
and slope at unfitted s=35 check the s=33 forecast. The three local
ten-direction rational fits at this same held-out point have
certifiably different responses; their source errors satisfy the
inherited unconditional ILC analytic estimate.

Neither finitely many charts nor finer anchor subdivision proves
positivity for every derivative order. The next exact obligation
for the infinite-layer claim is control of the ACTUAL source at
unbounded degree, including the LCW next-layer window membership
or equivalent ILC all-order inequalities. Extending the anchor
domain at fixed m is a separate obligation and cannot replace it.
RH, physical mass, Maxwell dynamics, quantum gravity and the
global eta bound are not closed or improved here.
