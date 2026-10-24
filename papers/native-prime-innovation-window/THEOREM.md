# PIW-1 — Sharp innovation windows and a new complete prime reading

17 October 2026 (Asia/Kolkata). Continues POC-1 in Publications.

POC-1 reduced the prime-sampled RH route to both actual Schur signs at
every step. This packet factors that acceptance test into a sharp value
window and a dimensionless slope window. It also computes a new complete
source jet at original s=29 and certifies its addition to the eight-prime
bank. The all-order sign statement remains open.

## Retained source and notation

Use the unchanged SE48/SE49 completed response

\[
q(t)=\mathcal X'(t)/\mathcal X(t),\qquad
\mathcal X(w^2)=\xi(1/2+w).
\]

At distinct positive nodes t_i retain actual u_i=q(t_i), q'(t_i), and
the real symmetric positive definite forms

\[
L_{ij}=-q[t_i,t_j],\qquad J_{ij}=(tq(t))[t_i,t_j].
\]

Thus J+LT=u1^T and J+TL=1u^T, with T=diag(t_i). Positivity is an
explicit premise at the old bank. For a new x>0 different from every
t_i, set

\[
P=J+xL,\quad v=P^{-1}u,\quad
r_0=q_N(x)=u^TP^{-1}u,\quad d_0=q_N'(x)=-v^TLv,
\]
\[
a_i=(x-t_i)^{-1},\quad
A=a^TL^{-1}a>0,\quad C=a^TJ^{-1}a>0,\quad b=a^Tv,
\qquad \kappa=x(A+xC)>0.
\tag{1}
\]

A proposed next source jet is (r,d). It is the actual pair
(q(x),q'(x)) when applying the result to RH. Define its residuals
delta=r-r_0 and delta'=d-d_0. The prime route uses x=(p-1/2)^2;
the algebra below applies to any distinct positive new node.

## PIW-T1. Exact value window and two-sign factorization

The actual next cross columns are

\[
\ell_i=(u_i-r)/(x-t_i),\qquad
j_i=(xr-t_iu_i)/(x-t_i).
\]

The SE49 innovations admit the exact residual expressions

\[
\boxed{\alpha=-\delta'+2b\delta-A\delta^2,\qquad
\beta=\delta+x\delta'-2xb\delta-x^2C\delta^2.}
\tag{2}
\]

In particular their combined budget has the factorization

\[
\boxed{x\alpha+\beta=\delta(1-\kappa\delta).}
\tag{3}
\]

There exists a derivative d making both enlarged forms strictly positive
if and only if

\[
\boxed{q_N(x)<r<q_N(x)+\kappa^{-1}.}
\tag{4}
\]

For a value in that open interval, the precise derivative-residual window is

\[
\boxed{
2b\delta-\delta/x+xC\delta^2
<\delta'<2b\delta-A\delta^2.}
\tag{5}
\]

Its width is delta(1-kappa delta)/x. Equivalently, putting

\[
\varphi=\kappa\delta,\qquad
\vartheta=\frac{x\alpha}{\delta(1-\kappa\delta)},
\tag{6}
\]

the actual next jet is accepted exactly when 0<varphi<1 and
0<vartheta<1. The second expression is used only where its denominator
is positive. Both coordinates are dimensionless.

**Proof.** The displacement relation gives
P=(xI-T)L+1u^T. Applying this to v yields
(xI-T)Lv=u-r_0 1. Consequently

\[
\ell=Lv-\delta a,\qquad j=Jv+x\delta a.
\tag{7}
\]

Expand -d-ell^T L^-1 ell and r+xd-j^T J^-1j. Use
r_0=v^TJv+xv^TLv and d_0=-v^TLv to obtain (2). The terms
containing delta' and b cancel in x alpha+beta, giving (3).
The two strict signs are equivalent to (5). Its interval is nonempty
exactly when delta(1-kappa delta)>0. Since kappa>0, this is
0<delta<1/kappa, proving (4). Finally beta equals the budget times
(1-vartheta), while x alpha equals the budget times vartheta. QED.

Both endpoints in (4) have zero budget, so no derivative at an endpoint
can make both signs strict. A value inside the interval alone is also
insufficient: a derivative outside (5) makes at least one innovation
nonpositive. Thus the value and slope are separate actual observations.
Equation (4) is an acceptance envelope, not an unconditional bound on an
arbitrary actual response whose next signs have not been established.

## PIW-T2. Covariant residual normalization and the slope cone

Define the rational transfer factor and normalized residual by

\[
h_N(x)=\frac{\det L\,\prod_i(x-t_i)}{\det(J+xL)},\qquad
g_N(x)=h_N(x)^2,\qquad
R_N(x)=\frac{q(x)-q_N(x)}{g_N(x)}.
\tag{8}
\]

Away from the retained nodes, g_N(x)>0 for x>0 and

\[
\frac{g_N'(x)}{g_N(x)}=2b(x).
\tag{9}
\]

Where the actual response is holomorphic at an old node, R_N has a
removable singularity there. This normalization removes the double
interpolation zeros already forced by the old observations. It does
not remove an unobserved source pole.

The next actual two-sign test is equivalently

\[
\boxed{
0<R_N<\frac1{\kappa g_N},\qquad
-\frac{R_N}{x}+xCg_NR_N^2<R_N'<-Ag_NR_N^2.}
\tag{10}
\]

In particular every accepted actual jet has R_N'<0. This local
consequence does not assert that R_N is decreasing everywhere from
positivity at a finite collection of samples.

Under any invertible real change of old source basis B, transport
L,J to B^TLB,B^TJB and u,a to B^Tu,B^Ta. Then r_0,d_0,A,C,b,
kappa,h_N,g_N,delta and both normalized coordinates (6) are unchanged.
The cross columns must be transported too. Merely changing matrix
coordinates while leaving the probes fixed does not have this property.

**Proof.** Left-multiplying P=(xI-T)L+1u^T by D=(xI-T)^-1 gives
DP=L+au^T. Taking the trace after right multiplication by P^-1
yields

\[
b=\sum_i(x-t_i)^{-1}-\operatorname{tr}(P^{-1}L).
\]

Jacobi's determinant identity proves (9); also h_N tends to 1 as
x tends to positive infinity. The determinant in (8) is nonzero at
every positive x, whereas each retained-node numerator zero is simple.
SE48 interpolation makes q-q_N vanish to order at least two at that
node, proving removability. At any other point where q has a pole and
g_N is finite and nonzero, that pole remains in R_N.

Substitute delta=g_NR_N and
delta'=g_N(R_N'+2bR_N) into (2). Then

\[
\alpha=-g_NR_N'-Ag_N^2R_N^2,\qquad
\beta=g_NR_N+xg_NR_N'-x^2Cg_N^2R_N^2.
\tag{11}
\]

Their strict signs give (10). The upper bound in (10) also follows
from (3). Since A,g_N,R_N^2 are positive, the upper derivative bound
is negative. Finally inversion under congruence gives v -> B^-1 v;
all displayed contractions are invariant. The determinant ratio is
unchanged because det(B)^2 cancels. This proves covariance. QED.

Combining POC-T5 and PIW-T1 gives a scalar form of the existing RH
criterion: after the initial positive prime bank, both coordinates in
(6) must lie in (0,1) at every subsequent actual prime. The identity
and criterion are proved; their universal source inequalities are not.
The residual slope cone makes that remaining obligation explicit.

## PIW-T3. Complete outward source jets without an assumed zero measure

For odd integer s>=3, w=s-1/2 and t=w^2, the full source jet is

\[
Q(s)=\frac1s+\frac1{s-1}-\frac12\log\pi
      +\frac12\psi(s/2)+\frac{\zeta'}{\zeta}(s),
\]
\[
Q'(s)=-\frac1{s^2}-\frac1{(s-1)^2}
       +\frac14\psi'(s/2)+\left(\frac{\zeta'}{\zeta}\right)'(s),
\]
\[
\boxed{q(t)=\frac{Q(s)}{2w},\qquad
q'(t)=\frac{wQ'(s)-Q(s)}{4w^3}.}
\tag{12}
\]

These follow by differentiating the completed xi definition and changing
the source coordinate. They hold in the Euler-product half-plane and
do not assume RH. The exact outward construction in
[SOURCE_EVALUATION.md](SOURCE_EVALUATION.md) includes every component:

- pi from Machin's formula with alternating-series remainder bounds;
- rational logarithms from a convergent atanh series with its geometric tail;
- digamma from a shifted Bernoulli expansion with a signed first-omitted-term
  bound, followed by exact recurrence;
- half-integer trigamma from pi squared and a finite rational sum;
- the complete Euler logarithmic derivative and its derivative, with all
  retained prime powers summed and every omitted prime bounded by explicit
  logarithmic integrals.

For cutoff M>=2, write a=s-1, ell=log M and c=1-(M+1)^-s. The positive
omitted Euler sums are bounded above by

\[
T_0=\frac{M^{-a}}c\left(\frac\ell a+\frac1{a^2}\right),
\qquad
T_1=\frac{M^{-a}}{c^2}
\left(\frac{\ell^2}a+\frac{2\ell}{a^2}+\frac2{a^3}\right).
\tag{13}
\]

The retained negative log derivative receives a tail in [-T_0,0];
its positive derivative receives one in [0,T_1]. All new numerical
values are native exact rational intervals rounded outwards to a
110-decimal-place grid. The written source identities and remainder
proofs justify those enclosures; floating-point agreement is not used
as certification.

## PIW-T4. Ninth prime accepted, with an independent held-out response

Keep POC's old eight prime directions s=3,5,7,11,13,17,19,23 and
their unchanged complete SE49 enclosures. Compute the new full jet at
s=29 using PIW-T3, with M=4096, digamma shift 256 and 24 Bernoulli
terms. Include all eight new mixed entries and the two new diagonals.
The resulting actual 9x9 L and J are strictly positive by outward
native LDL, independently agreeing with the direct and residual
versions of the Schur signs.

Conservative strict displays of the new signs and normalized readings are

\[
\alpha_9>1.36\times10^{-18},\qquad
\beta_9>6.68\times10^{-15},
\]
\[
0.140<\varphi_9<0.141,\qquad
0.142<\vartheta_9<0.143.
\tag{14}
\]

The previously unobserved original s=31 is separately evaluated with
the same complete-source method and is not inserted into this fit. At
that point the nine-prime response satisfies

\[
\boxed{0<q((61/2)^2)-q_9((61/2)^2)<3\times10^{-17}.}
\tag{15}
\]

The archive stores the full rational enclosures, including all source
tails, new mixed entries, pivots, acceptance coordinates and the separate
finite SE49 analytic error bound. A fresh s=23 jet agrees with the old
independently generated source box in both value and derivative.

This advances the prime bank from eight to nine directions and adds an
actual source sample beyond the old s<=26 archive. It does not certify
the union of all 24 old directions with the new direction as a positive
25-dimensional space. It also does not prove the next prime's two signs
merely by having computed its held-out value and derivative.

## Attribution, exact obligation and limits

SE48 supplies the spectral pencil and double interpolation; SE49 supplies
the original Schur acceptance and finite analytic bounds. POC supplies
the odd-prime sampling criterion. Completion of squares, displacement
identities, determinant differentiation and the special-function tools
are credited standard mathematics. This packet contributes their explicit
residual factorization and covariant normalization for the retained
source, plus the new complete native prime-reading execution. No general
priority claim for interpolation theory or Schur complements is made.

The next RH task is to control the actual normalized residual and its
slope in (10), or equivalently both coordinates (6), at all successive
prime steps with full source terms retained. A repeated finite check or
a freely chosen positive model cannot supply that universal theorem.
The global shift frontier remains 1079/2000. No RH proof, certified new
zeta zero, physical mass or quantum-gravity closure is asserted.
