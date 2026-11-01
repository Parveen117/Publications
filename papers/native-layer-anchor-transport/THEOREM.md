# LAT-1 — Moving anchors, exact layer flow and continuous source positivity

17 October 2026. Continues LCW-1 and ILC-1 on the same completed
xi response. This stage derives anchor motion directly from source
derivatives, then certifies actual finite-order positivity throughout
two continuous anchor intervals. The intervals are proved using a
uniform error bound, not a finite sampling grid or an assumed positive
xi-zero measure. Actual derivatives 18 and 19 are also evaluated.

Retain Xcal(w^2)=xi(1/2+w), q=Xcal'/Xcal. At a regular positive
anchor a with q(a)>0 define c_n(a)=(-a)^n q^(n)(a)/(n!q(a)).
The executed reference anchor is a_0=3249/4, w_0=57/2, original s=29.

## LAT-T1. The actual source has an exact anchor hierarchy

Along a=a_0 exp(tau), wherever q(a) is nonzero and regular,

\[
\boxed{\dot c_n=(n+c_1)c_n-(n+1)c_{n+1},\qquad \dot c_0=0.}
\tag{1}
\]

For D_(r,k)=sum_(j=0)^k (-1)^j binom(k,j)c_(r+j),

\[
\boxed{\dot D_{r,k}=(r+c_1)D_{r,k}-(r+k+1)D_{r+1,k}.}
\tag{2}
\]

**Proof.** Differentiate the definition with a partial_a. The power
a^n contributes n c_n, q^(n) contributes -(n+1)c_(n+1), and
1/q contributes c_1 c_n because a q'/q=-c_1. This proves (1).
Substitute into the binomial sum and collect coefficients using
j binom(k,j)=k binom(k-1,j-1); this gives (2). QED.

These are source identities without an RH assumption. They do not
close at a finite depth: the derivative of the final available layer
requires the next layer. A boundary sign is therefore not determined
by an autonomous finite evolution. Tau is the logarithm of anchor
scale, not a newly identified physical time.

## LAT-T2. Exact positive-model transport and retained invariants

Let nu be any probability measure on [0,1], and
f(lambda)=integral (1+lambda y)^-1 dnu(y). This may be the finite
positive model constructed from a certified source bank. For rho>0,
put delta=rho-1 and

\[
\phi_\rho(y)=\frac{\rho y}{1+\delta y},\qquad
d\nu_\rho(\phi_\rho(y))=
\frac{d\nu(y)}{(1+\delta y)f(\delta)}.
\tag{3}
\]

This is a probability measure on [0,1]. Successive moves rho and
sigma equal the move rho sigma, including the normalizing weights.
The corresponding response satisfies

\[
\frac{f(\delta+\rho\lambda)}{f(\delta)}
=\int\frac{d\nu_\rho(z)}{1+\lambda z}.
\tag{4}
\]

For finite interior atoms y_j, weights v_j, and source scale q(a),
the rational response has poles -mu_j and residues R_j, where

\[
\mu_j=a\frac{1-y_j}{y_j},\qquad R_j=aq(a)\frac{v_j}{y_j}.
\tag{5}
\]

Both quantities are unchanged by (3) with the new anchor rho a and
new source value q(a)f(delta). The represented rational function is
unchanged. These are finite-model poles, not identified xi zeros.

**Proof.** The denominator in (3) is positive. Its integral is
f(delta), so the weights normalize. Also
1-phi_rho(y)=(1-y)/(1+delta y), proving interval preservation.
Substitution proves (4). Composing the fractional transformations
gives phi_(rho sigma), and their denominator factors multiply to
1+(rho sigma-1)y; the normalizations then telescope. Substitution
in (5) proves both invariants. QED.

Along rho=exp(tau), dot y=y(1-y), while a normalized atom weight
satisfies dot v_j=(c_1-y_j)v_j. The moments therefore obey (1).
The binomial readings obey the explicit formula

\[
D_{r,k}(\rho a)=\frac{\rho^r}{f(\delta)}
\int\frac{y^r(1-y)^k}{(1+\delta y)^{r+k+1}}\,d\nu(y).
\tag{6}
\]

For the finite model every such reading stays positive when the
support is interior. Applying (3) to an actual all-order positive
source would require ILC's all-order premise. The finite theorem
below instead bounds the actual source against this model.

## LAT-T3. A complete source bound for anchor-transfer errors

At a_0 put f(lambda)=q(a_0(1+lambda))/q(a_0). In the closed disk
|lambda|<=R=15/16, choose the principal w=w_0 sqrt(1+lambda).
Then Re w>=w_0/4=57/8>1. The retained SE49 bound
|wq(w^2)|<41/8+|w|/4 gives

\[
|q(a_0(1+\lambda))|<221/228,\qquad
|c_n(a_0)|\le M R^{-n},\quad M=\frac{221}{228q(a_0)}.
\tag{7}
\]

In numerical certification replace q(a_0) in M by its positive lower
interval endpoint. This enlarges M safely. There is no midpoint model.

**Proof.** Re(1+lambda)>=1/16, and the principal-root formula gives
Re sqrt(1+lambda)>=1/4. Dividing the inherited majorant by |w|
gives 41/57+1/4=221/228. The actual source is holomorphic here
because the corresponding s=1/2+w lies in the Euler-product
zero-free half-plane. Cauchy's derivative formula on this circle
bounds each Taylor coefficient by M/R^n. This standard analytic
tool is credited to NIST DLMF 1.9.31,
<https://dlmf.nist.gov/1.9#E31>. QED.

Suppose K=2N-1 and both N-dimensional source forms are positive.
ILC gives an exact finite model f_N matching c_0,...,c_K, with
all its later moments between zero and one. For b=rho a_0 define
the unnormalized layers

\[
b_j(\rho)=\rho^j\frac{(-1)^j f^{(j)}(\rho-1)}{j!},
\quad \widetilde b_j(\rho)=
\rho^j\frac{(-1)^j f_N^{(j)}(\rho-1)}{j!}.
\tag{8}
\]

These are normalized by the old q(a_0). The true new c_j(b) is
b_j/b_0 once b_0>0. For |rho-1|<=d<R set

\[
T_{K,j}(z)=\frac{\binom{K+1}{j}z^{K+1-j}}
 {1-z(K+2)/(K+2-j)},
\]
\[
E_j=(1+d)^j\{M R^{-j}T_{K,j}(d/R)+T_{K,j}(d)\}.
\tag{9}
\]

Whenever 0<=j<=K and the displayed denominators are positive,
|b_j-tilde b_j|<=E_j uniformly throughout that interval.

**Proof.** The first K+1 coefficients of f-f_N cancel exactly.
The remaining coefficient of order n is bounded by M R^-n+1.
After j derivatives, the tail starts at k=K+1-j and has multiplier
binom(j+k,j). The ratio of successive nonnegative tail terms is
z(j+k+1)/(k+1), decreasing with k and bounded by its first value
z(K+2)/(K+2-j). The geometric sum gives (9). The bound on rho^j
is (1+d)^j. QED.

For a forecast from the finite actual Taylor polynomial alone, omit
the second term in (9). It bounds the actual remainder rather than
actual-minus-model. This separate forecast is checked against the
independent s=31 value and derivative.

## LAT-T4. Uniform actual positivity on an anchor interval

Choose m<=N and let A_m,B_m be the retained reference source forms.
Set

\[
\epsilon_A=\max_{1\le j\le2m-1}E_j,\qquad
\epsilon_B=\max_{0\le j\le2m-2}(E_j+E_{j+1}),
\]
\[
\tau_A=\frac{m\epsilon_A(1+d)^{2m}}{(1-d)^{2m-1}},\qquad
\tau_B=\frac{m\epsilon_B(1+d)^{2m}}{(1-d)^{2m-2}}.
\tag{10}
\]

If

\[
\boxed{A_m-\tau_A I>0,\qquad B_m-\tau_B I>0,}
\tag{11}
\]

then the actual normalized m-dimensional forms are positive at
every b=rho a_0, 1-d<=rho<=1+d. This is a finite sufficient test
with a universal quantifier over an interval of actual anchors.

**Proof.** For a polynomial p of degree below m put
P(y)=(1+delta y)^(m-1) p(phi_rho(y)). Its coefficient matrix is

\[
(T_\rho)_{ij}=\binom{m-1-j}{i-j}\rho^j\delta^{i-j},\quad i\ge j,
\tag{12}
\]

and zero otherwise. Direct composition gives
T_rho T_sigma=T_(rho sigma), so T_rho^-1=T_(1/rho).
The sum of absolute values in column j of this inverse is
rho^-j(1+|rho^-1-1|)^(m-1-j), at most (1-d)^(-(m-1)).

Writing the model forms in the common denominator (1+delta y)^(2m)
and using 1+delta y<=1+d gives their matrix lower bounds

\[
\widetilde A_m(\rho)\ \ge\frac{\rho}{(1+d)^{2m}}
T_\rho^T A_mT_\rho,\qquad
\widetilde B_m(\rho)\ \ge\frac1{(1+d)^{2m}}
T_\rho^T B_mT_\rho.
\]

Every actual-minus-model A entry has absolute value <=epsilon_A;
every B entry has absolute value <=epsilon_B. After congruence by
T_rho^-1, their quadratic forms have absolute error at most
m epsilon (1-d)^(-2m+2) times the Euclidean norm squared:
use |v^T E v|<=epsilon||v||_1^2 and ||v||_1<=sqrt(m)||v||_2.
The lower bounds and rho>=1-d now give exactly (10)–(11).
Thus both unnormalized actual forms are positive. Their sum is
the actual H form, whose (0,0) entry b_0 is positive. Division by
b_0 establishes the normalized claim. QED.

The implementation additionally checks b_0>=1/(1+d)-E_0>0.
Failure of (11) is inconclusive: it is not a proof of a negative
actual form. A fixed finite order over a continuous interval is
different from positivity at all derivative orders at one point.

## LAT-T5. Executed continuation, source intervals and cut interpretation

The unchanged complete evaluator at s=29 computes c_0,...,c_19.
Actual c_18,c_19 satisfy the LCW windows; the predicted old c_18
endpoints agree. Independent direct LDL certifies both 10x10 forms.
All 210 binomial readings r+k<=19 are positive. Every omitted
Gamma/Euler tail remains in the source enclosures.

Using K=19, the native replay certifies (11) for:

| Actual finite order | rho=b/a_0 interval |
|---|---|
| m=3 | [17/20,23/20] |
| m=4 | [11/12,13/12] |

The first interval includes the original-source points s=27 and s=31,
and every real anchor between them. The proof is uniform; no grid
interpolation is used. At the reference anchor there are ten positive
directions; ten-direction positivity is not claimed throughout these
intervals. The attempted m=3 interval [5/6,7/6] fails this sufficient
certificate and is explicitly recorded as inconclusive.

At each certified anchor, ILC supplies the finite positive contraction,
two-channel cut and two-reflection wave for that order. The larger
finite rational model transports by (3), but it does not equal the
whole actual source. The held-out response has a strictly positive
residual. The source flow (1) still requires higher layers, so neither
this interval theorem nor the semigroup closes the all-depth RH gate.
No Maxwell equation, physical time, mass gap, global eta improvement,
or physical quantum-gravity identification follows here.
