# LWP-1 — Exact wave horizons and certified exponential readings

17 October 2026. Continues LOT-1, LDC-1 and ILC-1 in the same
source track. The physical interpretation remains open. This stage
distinguishes an exact horizon for whole states from a longer exact
horizon for scalar correlations, then derives continuous unitary
transport with explicit finite-depth errors.

Retain the source forms H_N, A_N, B_N, Y_N=H_N^(-1)A_N and
the weighted pairing G_N=diag(A_N,B_N). All finite statements
assume strict A and B at the largest carrier used. Let P_N, Z_N,
W_N=(2P_N-I)Z_N and the weighted isometries V_(M,N) be those of
LOT-1. Adjoint sharp and norms use these acquired pairings.
An old state is a pair of polynomials of degree below N. No global
positive source measure is assumed in the finite theorems.

Carrier dimension, source derivative degree, discrete wave step k,
and exponential parameter tau are different quantities. Tau here is
neither source displacement, LAT's logarithmic anchor coordinate,
nor an identified physical time.

## LWP-T1. The source wave and its bounded skew generator

The weighted wave and its inverse are matrix polynomials of degree
one in Y_N:

\[
W_N=\begin{pmatrix}2Y_N-I&-2(I-Y_N)\\2Y_N&2Y_N-I\end{pmatrix},
\qquad
W_N^{\sharp}=W_N^{-1}
=\begin{pmatrix}2Y_N-I&2(I-Y_N)\\-2Y_N&2Y_N-I\end{pmatrix}.
\tag{1}
\]

Define their self-adjoint and skew-adjoint parts by

\[
\mathcal C_N=\frac{W_N+W_N^{\sharp}}2
=\operatorname{diag}(2Y_N-I,2Y_N-I),
\]
\[
\boxed{\mathcal L_N=[P_N,Z_N]
=\frac{W_N-W_N^{\sharp}}2
=\begin{pmatrix}0&-2(I-Y_N)\\2Y_N&0\end{pmatrix}.}
\tag{2}
\]

These commute, Lcal_N^sharp=-Lcal_N, and

\[
\mathcal L_N^2
=-4\operatorname{diag}(Y_N(I-Y_N),Y_N(I-Y_N)),\qquad
\|\mathcal L_N\|\le1.
\tag{3}
\]

For k>=1, with T_k and U_(k-1) the classical Chebyshev
polynomials of the first and second kinds,

\[
\boxed{W_N^{\pm k}=T_k(\mathcal C_N)
\ \pm\ U_{k-1}(\mathcal C_N)\mathcal L_N.}
\tag{4}
\]

**Proof.** Substitute LOT's rational P_N into the two products
(2P_N-I)Z_N and Z_N(2P_N-I), obtaining (1)-(2). They are
mutual inverses because P_N and Z_N supply two reflections.
Block multiplication proves (3). Alternatively the norm bound
follows from (2) and the unit norms of W_N and its inverse.
Since Ccal_N and Lcal_N are polynomials in W_N and W_N^(-1),
they commute. The identity W_N+W_N^(-1)=2Ccal_N gives
W_N^(k+1)=2Ccal_N W_N^k-W_N^(k-1). Its initial values are I
and Ccal_N+Lcal_N. The Chebyshev recurrence and its initial
values prove the plus case of (4); inversion changes the sign of
Lcal_N and proves the minus case. QED.

The corresponding polynomial symbols are

\[
\mathsf W(y)=\begin{pmatrix}2y-1&-2(1-y)\\2y&2y-1\end{pmatrix},
\qquad
\mathsf L(y)=\begin{pmatrix}0&-2(1-y)\\2y&0\end{pmatrix}.
\tag{5}
\]

Both W_N^k and Lcal_N^j are matrix polynomials in Y_N, of degree
at most |k| and j respectively. The Chebyshev recurrence is
classical; see DLMF 18.9.1 and Table 18.9.1,
<https://dlmf.nist.gov/18.9>, and the representations in
<https://dlmf.nist.gov/18.5>. The source specialization is (1)-(5).

## LWP-T2. Whole-state and scalar-reading horizons differ

Fix N<=M<=M' and put d=M-N. For every old polynomial-pair state z,

\[
\boxed{W_{M'}^kV_{M',N}z
=V_{M',M}W_M^kV_{M,N}z\quad (|k|\le d).}
\tag{6}
\]

The common value is the coefficient vector of Wsf(y)^k z(y).
The analogous equality for Lcal^j holds for 0<=j<=d.

For old test and input states u,z, define
gamma_M(k)=<V_(M,N)u,W_M^k V_(M,N)z>_(G_M). Then

\[
\boxed{\gamma_{M'}(k)=\gamma_M(k)\quad (|k|\le2d).}
\tag{7}
\]

The same scalar equality for Lcal^j holds for 0<=j<=2d. In
particular a scalar reading can be exact across depths after the
whole propagated state has reached a new boundary direction.

**Proof.** Multiplication by y is exact on degree-below-M
polynomials as long as its output still has degree below M.
Applying |k| degree-one wave factors to z produces degree at
most N-1+|k|. The inverse wave in (1) has the same degree bound.
Thus every intermediate and final multiplication is exact for
|k|<=M-N, proving (6). The same reasoning applies to Lsf^j.

For the scalar assertion, let nu_M be the positive finite model
measure supplied by ILC for carrier M. It matches source moments
c_0,...,c_(2M-1); it is not an assumed infinite source measure.
With g(y)=diag(y,1-y), finite spectral evaluation gives

\[
\gamma_M(k)=\int u(y)^*g(y)\mathsf W(y)^kz(y)\,d\nu_M(y).
\tag{8}
\]

The integrand has degree at most
(N-1)+1+|k|+(N-1)=2N+|k|-1. When |k|<=2(M-N), this is at
most 2M-1. Its integral is therefore the coefficient functional
that replaces y^r by the actual finite source moment c_r.
Every deeper certified model gives exactly the same value. This
proves (7), and the identical degree count proves its Lcal version.
QED.

These are sufficient uniform horizons, with a strict first failure
available. For z=(y^(N-1),0) and k=d+1, the polynomial Wsf(y)^k z
has degree M with nonzero leading coefficient
2^(2k-1)(1,1). The M+1 carrier represents it exactly, while the
M carrier has no such coefficient. Its positive next pairing
prevents this nonzero difference from being a null state. Thus
the vector horizon cannot generally be enlarged by one.

A polynomial observable F(y) of matrix degree q gives the same
finite source identification for <u,F W^k z> when
q+|k|<=2d. For the corresponding same-step reading
<W^k u,F W^k z>, the sufficient bound is q+2|k|<=2d.
Only an F that is self-adjoint in the specified pairing is a
self-adjoint observable; the algebraic correlation identity itself
does not require that additional property.

## LWP-T3. The first unresolved return is a known boundary norm

Let e_+=(e_0,0) and define the real return correlation

\[
a_M(k)=\langle e_+,W_M^k e_+\rangle_{G_M}
=\int yT_{|k|}(2y-1)\,d\nu_M(y).
\tag{9}
\]

Its initial value is c_1, the squared norm of e_+; it is not a
probability. The normalized return is a_M(k)/c_1. Whenever
|k|<=2M-2, (9) is exactly determined by the retained source
moments through degree 2M-1. For a strict M to M+1 extension,

\[
\boxed{a_{M+1}(2M-1)-a_M(2M-1)
=2^{4M-3}h_M>0,\qquad
h_M=\frac{\det H_{M+1}}{\det H_M}.}
\tag{10}
\]

**Proof.** The top-left block of (4) is T_k(2Y_M-I), proving
(9). The polynomial yT_k(2y-1) has degree k+1. For k>=1 its
leading coefficient is 2^(2k-1). At k=2M-1 all lower coefficients
multiply moments already shared by the two models. LDC-T2 gives
the first new moment difference c_(2M)-<e_0,Y_M^(2M)e_0>_(H_M)
=h_M. Multiplying it by the leading coefficient proves (10).
QED.

More generally, the test and input pair (y^(N-1),0) has its first
possible scalar failure at k=2(M-N)+1. The same argument gives
a strictly positive difference 2^(2k-1)h_M between carriers M+1
and M. Thus the doubled scalar horizon in (7) is also a genuine
uniform boundary, rather than a claim of indefinite stabilization.

## LWP-T4. Continuous transport with rational depth-error bounds

For a real, dimensionless parameter tau define

\[
\boxed{\mathscr U_M(\tau)=\exp(\tau\mathcal L_M).}
\tag{11}
\]

This is a unitary group in G_M. The generator is selected directly
from the existing cut and channel reflection by (2). It is not
a logarithm of W_M, and Uscript_M(1)=W_M is not asserted.
Different choices of generator can define different continuous
groups; no unique physical interpolation is claimed.

For r>=1 and s>=0 set

\[
\mathcal E_r(s)=\sum_{j=r}^{\infty}\frac{s^j}{j!},\qquad
\boxed{\mathcal E_r(s)\le B_r(s):=
\frac{s^r/r!}{1-s/(r+1)}\quad(0\le s<r+1).}
\tag{12}
\]

For any two certified depths M'>=M>=N, with d=M-N,

\[
\boxed{\left\|\mathscr U_{M'}(\tau)V_{M',N}z
-V_{M',M}\mathscr U_M(\tau)V_{M,N}z\right\|
\le2\|z\|\mathcal E_{d+1}(|\tau|).}
\tag{13}
\]

For the scalar correlation
b_M(tau;u,z)=<V_(M,N)u,Uscript_M(tau)V_(M,N)z>_(G_M),

\[
\boxed{|b_{M'}(\tau;u,z)-b_M(\tau;u,z)|
\le2\|u\|\|z\|\mathcal E_{2d+1}(|\tau|).}
\tag{14}
\]

Replace Ecal by B in each formula whenever its displayed domain
condition holds. After division by the indicated nonzero state
norms these are relative state and normalized correlation bounds.
Each normalized right side may also be capped at two by unitarity.

**Proof.** The exponential series converges absolutely in operator
norm because Lcal_M is bounded. Its adjoint is exp(-tau Lcal_M),
and multiplication of the convergent series gives the inverse and
group properties. Expand the two sides of (13). By LWP-T2 the
coefficients through j=d cancel exactly. Each remaining term has
norm at most 2||z|| |tau|^j/j! because both generators have norm
at most one and the inclusions are isometries. Summation proves
(13). In the scalar expansion, LWP-T2 cancels coefficients
through j=2d; the remaining correlations are bounded by
2||u||||z|| |tau|^j/j!, proving (14). Finally successive terms
in Ecal_r have ratio at most s/(r+1), so summing that geometric
majorant proves (12). QED.

The tail series bounds hold for every finite real tau. The rational
formula B_r has the explicit additional condition s<r+1; it is
not continued through its pole. The source moments and native
rational arithmetic suffice to compute the finite expressions and
these bounds, without a new numerical exponential or root engine.

## LWP-T5. A sharper common exponential return interval

For the same seed e_+ define b_M(tau)=b_M(tau;e_+,e_+). Then

\[
b_M(\tau)=\int y\cos\!\left(2\tau\sqrt{y(1-y)}\right)d\nu_M(y).
\tag{15}
\]

Given source moments through c_(2M-1), form the finite polynomial

\[
\boxed{S_M(\tau)=\sum_{j=0}^{M-1}
\frac{(-4\tau^2)^j}{(2j)!}\,D_{j+1,j},\qquad
D_{r,k}=\sum_{l=0}^k(-1)^l\binom kl c_{r+l}.}
\tag{16}
\]

For EVERY compatible certified carrier M'>=M, including M itself,

\[
\boxed{|b_{M'}(\tau)-S_M(\tau)|
\le c_1\frac{|\tau|^{2M}}{(2M)!}.}
\tag{17}
\]

Thus the normalized return b_(M')/c_1 belongs to the common
interval centered at S_M/c_1 with radius |tau|^(2M)/(2M)!.
This is a real correlation interval, not a probability interval.

**Proof.** Equation (3) implies that the even powers of Lsf(y)
are (-4y(1-y))^j I and its odd powers have zero diagonal.
Expanding (11) proves (15). The Taylor polynomial of cosine
through degree 2M-1 has only the even terms j=0,...,M-1.
Its real remainder has absolute value at most
|x|^(2M)/(2M)!, since every real derivative of cosine is bounded
by one. Here x=2tau sqrt(y(1-y)) and 2sqrt(y(1-y))<=1.

After multiplication by y and integration against nu_(M'),
the remainder is therefore bounded by
c_1 |tau|^(2M)/(2M)!. The jth retained integral is
integral y^(j+1)(1-y)^j dnu_(M')=D_(j+1,j), whose degree
2j+1 is at most 2M-1. Every such moment is a shared actual
source reading. This gives exactly (16)-(17). QED.

In particular twenty source layers give M=10 and the normalized
common radius 1/20! at tau=1. Rational tau, the finite differences
in (16), and the factorial remainder produce outward enclosures
entirely in the retained native arithmetic. The square root in
(15) is a proof representation; it is not numerically evaluated.
An outward interval for c_1 must be kept positive when dividing
S_M, and the resulting interval arithmetic preserves source
uncertainty rather than replacing it by a midpoint.

The universal quantifier in (17) concerns compatible finite positive
carriers. It does not assert that the actual source has all those
extensions. Under the additional all-order premise below, the
same interval also contains the infinite-source return.

## LWP-T6. Conditional strong limit, uniform on bounded parameter sets

Assume the additional all-order actual-source positivity premise
of ILC-T2. Use LOT-T6's intrinsic completion
Gcal_infinity=L2(y dnu) direct_sum L2((1-y)dnu), and its weighted
polynomial inclusions. Define Lcal by the multiplication symbol
Lsf(y) in (5), and Uscript(tau)=exp(tau Lcal).

Then Lcal is bounded skew-adjoint with norm at most one, and
the finite exponential actions converge strongly to Uscript(tau),
uniformly for tau in every compact real interval. For old
polynomial-pair states the bounds (13)-(14) remain valid with
the larger carrier replaced by this completion. Equation (17)
also holds for the corresponding infinite-source seed return.

**Proof.** On the weighted space, direct multiplication gives
Lcal^sharp=-Lcal and

\[
\|\mathcal L(f,g)\|^2
=\int4y(1-y)\bigl(y|f|^2+(1-y)|g|^2\bigr)d\nu
\le\|(f,g)\|^2.
\]

Endpoint equivalence classes are respected: an undefined
channel value at a zero weight is multiplied by that weight
in the other channel. No inverse endpoint weight is used.
The bounded exponential therefore defines a unitary group.

For a fixed polynomial pair, finite powers agree with the
infinite multiplication powers through the same degree horizons.
The proofs of (13)-(14) apply unchanged. For |tau|<=T, their
tails tend to zero as M increases, uniformly in tau. Polynomial
pairs are dense in the weighted completion by LOT-T6; approximate
an arbitrary state by such a pair and use the unit norms of
both exponential actions to extend strong convergence to it.
One may extend each finite generator by zero off its carrier;
its exponential is then identity on that complement, making the
preceding comparison an ordinary strong limit on one space.
The proof of (17) applies to nu in place of nu_(M'). QED.

This conditional limit does not establish the existence of all
positive actual-source extensions. A finite wave horizon, a
finite common return interval and an all-order source theorem
remain different conclusions.

## Source execution and classical attribution

The pinned source jets at s=29,31,33 provide M=10 returns through
|k|=18, and the data for nine strict first-return discrepancies
at each anchor. No fresh source derivative is required. The
twenty layers also supply every coefficient in S_10 for the
finite exponential-return intervals.

The Chebyshev recurrence, exponential series and real Taylor
remainder are classical tools. Besides the references after
LWP-T1, see DLMF 4.2.19 for the exponential series,
<https://dlmf.nist.gov/4.2#E19>; DLMF 4.19.2 for cosine,
<https://dlmf.nist.gov/4.19#E2>; and DLMF 1.4.35-36 for the
real Taylor formula and remainder,
<https://dlmf.nist.gov/1.4#E36>. No priority claim is made for
these tools. The source-specific horizons, first discrepancies
and common finite intervals are derived from the retained
ILC moment matching, LDC boundary norm and LOT weighted carrier.
