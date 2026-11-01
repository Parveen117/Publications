# RBP-1 — Finite-anchor positivity for relational observer banks

17 October 2026. Written continuation at Publications commit
e29ce1a9a3089796da22bf99871abb399c8f429a, following
[RLO-1](RH_RELATIONAL_OBSERVER_CONFLUENCE.md).

This continuation proves an actual finite-anchor range for RLO's
fixed dyadic relational matrices. It does not stop at their positive
infinite-anchor limit.

Put
\[
b=\sqrt a,\qquad \ell=\log\frac b{2\pi},\qquad
\kappa=\frac{e^{-4}}{2\pi},\qquad
N(a)=\left\lfloor\frac{\kappa b}{2}\right\rfloor,\qquad
m_+(a)=\left\lfloor\frac{N(a)}8\right\rfloor .
\]
For N(a)>=128 and every 1<=m<=m_+(a), the actual first m dyadic
observers x_j=2^(-j-1) satisfy
\[
\boxed{
\mathsf K_m(a)\succ0,\qquad
\mathfrak D_m(a)
=\frac{\det\mathsf K_m(a)}{\prod_{i<j}(x_j-x_i)^2}
\ge\left(\frac1{24\ell\,49^m}\right)^m>0.
}
\]
Thus an explicit unbounded initial segment of the actual relational
criterion is proved at each covered finite anchor. Its size grows
as \(e^{-4}\sqrt a/(32\pi)\).

The proof first controls the wave Toeplitz matrix directly. This
improves the retained square-root constant and establishes both
shifted cut forms A_N and B_N, as well as H_N. Their acquired positive
pairing constructs a finite source carrier whose first 2N moments
match the actual source. The resulting positive rational kernel is
compared with the actual one, retaining the entire omitted-moment
error.

**Status: RH remains unproved.** The established sizes are bounded
by one explicit multiple of sqrt(a), not every size at one fixed
finite anchor or every square-root width. No larger support window
L is claimed. The written result is not a numerical certificate
run; no scientific computation, independent agent or replay is used.

## 1. Source and observer definitions

Retain the actual CWR/SQR source at a>=a_0=3249/4:
\[
\begin{gathered}
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\\
u_\mu=\frac a{a+\mu},\qquad
w_\mu=\frac{\nu_\mu}{q(a)(a+\mu)},\qquad
\mathcal L_a[f]=\sum_\mu w_\mu f(u_\mu),\\
q(a)>0,\qquad 0<|u_\mu|<1,\qquad
\sum_\mu |w_\mu|\le2,\qquad c_j(a)=\mathcal L_a[y^j],\quad c_0=1.
\end{gathered}
\tag{1}
\]
The actual source values and weights may be complex. Its products
remain bilinear, with reality supplied by conjugation symmetry.

Define the source moment forms
\[
H_N=(c_{i+j})_{i,j<N},\qquad
A_N=(c_{i+j+1})_{i,j<N},\qquad
B_N=H_N-A_N.
\tag{2}
\]
Use the arcsine reference only for estimates:
\[
d\nu_*(y)=\frac{dy}{\pi\sqrt{y(1-y)}},\qquad
H_N^*=\left(\int y^{i+j}\,d\nu_*\right),\quad
A_N^*=\left(\int y^{i+j+1}\,d\nu_*\right),\quad
B_N^*=H_N^*-A_N^*.
\tag{3}
\]

The native wave and its Toeplitz bank are
\[
r_k(a)=\mathcal L_a[T_k(2y-1)],\quad r_{-k}=r_k,\quad r_0=1,\qquad
\mathsf T_M=(r_{i-j})_{0\le i,j<M}.
\tag{4}
\]
RLO's separate relational kernel is
\[
F_a(x)=\frac{q(a(1+x))}{q(a)},\qquad
\mathsf K_a(x,z)=\mathcal L_a\!\left[
\frac{y}{(1+xy)(1+zy)}\right].
\tag{5}
\]
Its diagonal is -F_a'(x); off the diagonal it is
(F_a(x)-F_a(z))/(z-x). The symbols T_M and K_m therefore
refer to different, explicitly related observer matrices.

| Parameter | Role |
|---|---|
| a and b=sqrt(a) | Actual source anchor and its square root |
| ell | The logarithm log(b/(2 pi)) |
| M | Number of native wave directions |
| N | Number of polynomial/carrier directions |
| m | Number of relational observer positions |
| x_j | Fixed dyadic positions 2^(-j-1) |
| t_j and omega_j | Nodes and weights of a constructed finite carrier |

Carrier nodes are not asserted to be the actual spectral source
values u_mu. Finite matching and the retained error distinguish them.

## RBP-T1. A direct wave estimate improves the square-root positivity constant

For
\[
b\ge64,\qquad 2\le M\le b/8,
\]
define
\[
d_M^{\rm wave}(a)=
\frac{\log M+2+384M/b}{\ell}.
\]
Then
\[
\boxed{
\|\mathsf T_M-I\|\le d_M^{\rm wave}(a).
}
\tag{6}
\]
In particular, if
\[
2\le M\le\kappa b,\qquad \kappa=e^{-4}/(2\pi),
\]
then
\[
\boxed{
\frac1{2\ell}I\preceq\mathsf T_M\preceq2I.
}
\tag{7}
\]

**Proof of the wave bound.** SQR-T2 retains the complete rational,
Gamma and prime-power source in
\[
C_a(z)=1+2\sum_{k\ge1}r_kz^k
=1+\frac1\ell\log\frac{1-z}{1+z}+E_a(z).
\]
On the circle varrho=1-1/M it proves
\[
\sup_{|z|\le\varrho}|E_a(z)|
\le\frac{128M}{b\ell},\qquad
E_a(z)=2\sum_{k\ge1}e_kz^k.
\tag{8}
\]
Consequently
r_k=-1/(k ell)+e_k for odd k, and r_k=e_k for positive
even k. The leading off-diagonal Toeplitz correction has absolute
row sum at most
\[
\frac2\ell\sum_{\substack{1\le k<M\\k\ {\rm odd}}}\frac1k
\le\frac{\log M+2}{\ell}.
\tag{9}
\]

For completeness, the lower triangular matrix with entries
e_(i-j) varrho^(i-j), i>j, is the compression of multiplication
by E_a(varrho z)/2 on the unit circle. Parseval gives norm at
most sup|E_a|/2. Conjugation by
diag(1,varrho,...,varrho^(M-1)) removes the radial factors.
Adding its transpose therefore bounds the Toeplitz error by
\[
\varrho^{-(M-1)}\sup|E_a|
\le e\,\sup|E_a|
<\frac{384M}{b\ell}.
\]
Together with (9) this proves (6). The reference circle controls
analytic coefficients; no positivity of the arithmetic source is
assumed in this argument.

**Proof of the explicit range.** The condition M<=kappa b gives
log M<=ell-4. It also ensures b>64 and M<b/8 when M>=2.
Using e>8/3 and pi>3,
\[
384\kappa<\frac{81}{64},\qquad
1-d_M^{\rm wave}
\ge\frac{2-384\kappa}{\ell}
>\frac{47}{64\ell}>\frac1{2\ell}.
\]
Equation (6) now gives (7). QED.

This estimate uses the Toeplitz bank directly, avoiding the extra
Hankel part present in SQR's unshifted Chebyshev-bank bound.

## RBP-T2. Both cut channels are positive through a growing carrier size

If
\[
1\le N\le\left\lfloor\frac{\kappa b}{2}\right\rfloor,
\]
then
\[
\boxed{
\begin{aligned}
\frac1{2\ell}A_N^*&\preceq A_N\preceq2A_N^*,\\
\frac1{2\ell}B_N^*&\preceq B_N\preceq2B_N^*,\\
\frac1{2\ell}H_N^*&\preceq H_N\preceq2H_N^*.
\end{aligned}
}
\tag{10}
\]
All three forms are strictly positive.

**Proof.** Introduce a formal real Laurent functional
R_a[z^k]=r_|k|, with the involution z^*=z^(-1). Put
\[
y(z)=\frac{2+z+z^{-1}}4.
\]
The Chebyshev identity
T_k(2y(z)-1)=(z^k+z^(-k))/2 gives
R_a[p(y(z))]=L_a[p(y)] for every real polynomial p.

For degree p<N, the two expressions
\[
P_\pm(z)=\frac{1\pm z}{2}\,z^{N-1}p(y(z))
\]
are ordinary polynomials of degree at most 2N-1. Their formal
products satisfy
\[
P_+^*P_+=y(z)p(y(z))^2,\qquad
P_-^*P_-=(1-y(z))p(y(z))^2.
\tag{11}
\]
If v_+ and v_- are their coefficient vectors, (4) gives
R_a[P_\pm^*P_\pm]=v_\pm^T T_(2N) v_\pm.

For the reference unit circle, Parseval and
y=(1+cos(theta))/2 give
\[
\|v_+\|^2=\int_0^1 y p(y)^2\,d\nu_*,\qquad
\|v_-\|^2=\int_0^1(1-y)p(y)^2\,d\nu_*.
\]
Apply (7) with M=2N to prove the first two comparisons in
(10). Their sum gives the third. The reference weighted
integrals are positive for every nonzero polynomial. QED.

The Laurent involution in (11) is an algebraic identity.
Absolute squares on the unit circle are used only to estimate
the reference coefficient norm. No absolute-square substitution
is made at the actual complex arithmetic source points.

SQR established unshifted positivity through coefficient
e^(-9)/(2 pi) times sqrt(a). Here the simultaneous carrier
coefficient is e^(-4)/(4 pi), larger by e^5/2.
This compares the retained analytic bounds, not an optimal range.
Unlike SQR's stated unshifted result, (10) supplies both A_N and
B_N, the finite positivity premises required by LOT's cut model.

## RBP-T3. The acquired positive forms give an exact finite probability carrier

At any size covered by (10), define
\[
Y_N=H_N^{-1/2}A_NH_N^{-1/2},\qquad
v_N=H_N^{1/2}e_0.
\tag{12}
\]
Then
\[
0\prec Y_N\prec I,\qquad \|v_N\|^2=1.
\]
There are N distinct nodes t_j in (0,1) and positive weights
omega_j summing to one such that
\[
\boxed{
c_k(a)=v_N^TY_N^kv_N
=\sum_{j=1}^N\omega_jt_j^k
\qquad(0\le k\le2N-1).
}
\tag{13}
\]

**Proof.** Positivity of A_N and B_N gives the strict operator
inequalities for Y_N, while c_0=1 gives the seed norm.
On degree-below-N polynomials, H_N is the acquired positive
pairing. Write its compressed multiplication as
\(\mathcal Y_N=H_N^{-1}A_N\), which is conjugate to Y_N
by H_N^(1/2). For degree p<=N-2 this compression equals
actual multiplication by y, because the product still belongs
to the carrier.
Thus the powers of the seed through N-1 represent
1,y,...,y^(N-1).

For k<=2N-2, split k=i+j with i,j<=N-1 and use the H_N
pairing to get v_N^T Y_N^k v_N=c_k. For the last moment,
\[
v_N^TY_N^{2N-1}v_N
=\langle y^{N-1},\mathcal Y_N y^{N-1}\rangle_{H_N}
=(A_N)_{N-1,N-1}=c_{2N-1},
\]
where the middle expression uses polynomial carrier coordinates.

The seed is cyclic because its first N powers span that carrier.
A real symmetric cyclic operator has simple eigenvalues, and the
seed has nonzero projection onto every eigenvector. Its spectral
weights are therefore strictly positive. The spectral theorem
now proves (13). QED.

The positive inner product has been obtained from the proved
finite source forms. An infinite positive completion is not
assumed. Equation (13) matches exactly 2N moments; it does not
replace all subsequent actual-source moments.

The native cut has an explicit angular realization. In the
Euclidean carrier coordinates set
\[
\widehat P_N=
\begin{pmatrix}
Y_N&\sqrt{Y_N(I-Y_N)}\\
\sqrt{Y_N(I-Y_N)}&I-Y_N
\end{pmatrix}.
\tag{14}
\]
It is an orthogonal projection. On an eigenchannel t_j, writing
theta_j=arccos(sqrt(t_j)) gives its block
\[
\begin{pmatrix}
\cos^2\theta_j&\cos\theta_j\sin\theta_j\\
\cos\theta_j\sin\theta_j&\sin^2\theta_j
\end{pmatrix},\qquad 0<\theta_j<\pi/2.
\tag{15}
\]
This is LOT's weighted cut after normalizing its acquired metric.
The probabilities follow from the verified finite positivity,
with the coupling retained.

## RBP-T4. A positive rational model retains an explicit actual-source error

Use the finite carrier (13) to define
\[
\begin{aligned}
F_N(x)&=v_N^T(I+xY_N)^{-1}v_N
       =\sum_j\frac{\omega_j}{1+xt_j},\\
\mathsf K_N^{\rm mod}(x,z)
&=v_N^T(I+xY_N)^{-1}Y_N(I+zY_N)^{-1}v_N\\
&=\sum_j\frac{\omega_jt_j}{(1+xt_j)(1+zt_j)} .
\end{aligned}
\tag{16}
\]
This is a positive Gram kernel. On any m<=N distinct real
positions in [0,1), its matrix is strictly positive.
For |x|,|z|<=h<1, the actual source obeys
\[
\boxed{
\begin{aligned}
|F_a(x)-F_N(x)|
&\le\frac{3h^{2N}}{1-h},\\
|\mathsf K_a(x,z)-\mathsf K_N^{\rm mod}(x,z)|
&\le
\epsilon_N(h):=
\frac{3h^{2N-1}\bigl(2N-(2N-1)h\bigr)}{(1-h)^2}.
\end{aligned}
}
\tag{17}
\]
For an m-position matrix this gives the finite enclosure
\[
\boxed{
\mathsf K_m^{\rm mod}-m\epsilon_N(h)I
\preceq\mathsf K_m(a)
\preceq
\mathsf K_m^{\rm mod}+m\epsilon_N(h)I .
}
\tag{18}
\]

**Proof.** Positivity in (16) is the ordinary Euclidean Gram
identity with vectors
Y_N^(1/2)(I+xY_N)^(-1)v_N. Distinct nodes t_j, positive
weights and the Cauchy determinant identity make every
m<=N distinct-position matrix strictly positive.

Write d_k=c_k-sum_j omega_j t_j^k. Equation (13) gives d_k=0
for k<2N. The actual budget two and model budget one give
|d_k|<=3 at every remaining order. Expand the resolvents:
\[
F_a-F_N=\sum_{k\ge2N}(-x)^k d_k,\qquad
\mathsf K_a-\mathsf K_N^{\rm mod}
=\sum_{i,j\ge0}(-x)^i(-z)^j d_{i+j+1}.
\]
The absolute geometric sums yield (17); specifically
sum_(k>=2N-1)(k+1)h^k
=h^(2N-1)(2N-(2N-1)h)/(1-h)^2.
The row-sum norm bound for the real symmetric error matrix
proves (18). QED.

There is an exact determinant comparison for the positive model.
For m<=N distinct positions x_i in [0,h],
\[
\boxed{
(1+h)^{-2m^2}\det A_m(a)
\le
\frac{\det\mathsf K_m^{\rm mod}}
     {\prod_{i<j}(x_j-x_i)^2}
\le\det A_m(a).
}
\tag{19}
\]
Indeed Cauchy--Binet and the Cauchy determinant formula make the
middle expression
\[
\sum_{\substack{J\subseteq\{1,\ldots,N\}\\|J|=m}}
\frac{
 \prod_{j\in J}\omega_jt_j
 \prod_{\substack{j<k\\j,k\in J}}(t_k-t_j)^2}
{\prod_{i=0}^{m-1}\prod_{j\in J}(1+x_it_j)^2}.
\]
Every term is positive. Removing its denominator gives det A_m,
because the required moments are among those matched in (13).
Each denominator lies between one and (1+h)^(2m^2), proving
(19). Its coincident extension is the same normalized relation
used in RLO.

A positive model alone does not establish the sign of the actual
kernel. Formula (18) retains the missing tail. For the canonical
dyadic sequence, the next theorem controls that tail directly
after the Newton observer change, avoiding small-spacing losses.

## RBP-T5. The dyadic Newton matrices have a uniform comparison bound

Take x_i=2^(-i-1), i=0,...,m-1, with m<=N. Retain RLO's
Newton observers
\[
b_i(y)=\frac{y^i}{\prod_{k=0}^i(1+x_ky)}.
\]
Let
\[
G_m=(\mathcal L_a[y b_i b_j])_{i,j<m},\qquad
G_m^{\rm mod}=\left(\sum_{\ell=1}^N
\omega_\ell t_\ell b_i(t_\ell)b_j(t_\ell)\right)_{i,j<m}.
\]
If the carrier size N is covered by (10), then
\[
\boxed{
G_m^{\rm mod}\succeq\frac1{12\ell\,49^m}I,\qquad
\|G_m-G_m^{\rm mod}\|
\le2304m\left(\frac32\right)^{-(2N-2m+1)}.
}
\tag{20}
\]
RLO's exact congruence gives
\[
G_m=D_m\mathsf K_m(a)D_m^T,\qquad
\det G_m=\mathfrak D_m(a).
\tag{21}
\]

**A reference eigenvalue bound.** First,
\[
\boxed{\lambda_{\min}(A_m^*)\ge24\,49^{-m}.}
\tag{22}
\]
To prove it, put y=(1+cos(theta))/2 and define
\[
\chi_j(y)=\sqrt2\,
\frac{\cos((j+1/2)\theta)}{\cos(\theta/2)} .
\]
The quotient is a polynomial of degree j: its first two
unnormalized polynomials are 1 and 4y-3, and its recurrence is
v_(j+1)=(4y-2)v_j-v_(j-1). They are orthonormal for y dnu_*,
because multiplication by y cancels the squared cosine
denominator and half-integer cosines are orthogonal on (0,pi).

The coefficient l1 norm of v_j is at most 7^j. This follows
from the initial norms 1,7 and the recurrence bound
B_(j+1)<=6B_j+B_(j-1)<=7^(j+1).
If Q contains the coefficient columns of chi_0,...,chi_(m-1),
then
\[
Q^TA_m^*Q=I,\qquad
\|Q\|^2\le\|Q\|_F^2
\le2\sum_{j<m}49^j<49^m/24 .
\]
Inverting the congruence proves (22).

**The positive-model bound.** For a vector v in R^m write
\[
\sum_{i<m}v_i b_i(y)=\frac{p_v(y)}{D(y)},\qquad
D(y)=\prod_{k=0}^{m-1}(1+x_ky),
\]
\[
p_v(y)=\sum_{i<m}v_i y^i
\prod_{k=i+1}^{m-1}(1+x_ky).
\]
Let T map v to the coefficient vector of p_v.
It is triangular with diagonal one. For any column, the sum
of off-diagonal coefficients is at most
exp(sum_(k>i)x_k)-1<=sqrt(e)-1.
For a row, the coefficient of order r is at most (1/2)^r/r!,
so the same bound holds. Therefore
\[
\|T-I\|\le\sqrt e-1<3/4,\qquad
\|Tv\|\ge\|v\|/4.
\tag{23}
\]
Here e<3 and sqrt(3)<7/4 suffice for the strict bound.

At the positive carrier nodes, 1<=D(t)<=e<3.
Positivity is used at these finite model nodes, not at the
possibly complex actual source. Moment matching gives
\[
\begin{aligned}
v^TG_m^{\rm mod}v
&\ge\frac19\sum_j\omega_jt_j p_v(t_j)^2\\
&=\frac19(Tv)^TA_m(a)(Tv)\\
&\ge\frac1{18\ell}(Tv)^TA_m^*(Tv)
\ge\frac1{12\ell\,49^m}\|v\|^2 .
\end{aligned}
\]
The equality needs degree at most 2m-1<=2N-1, precisely the
available matching range. This proves the first part of (20).

**The actual-source remainder.** On |y|=R=3/2,
\[
\prod_{k=0}^i|1+x_ky|
\ge\prod_{k=0}^i(1-Rx_k)\ge\frac1{16}.
\tag{24}
\]
Indeed the first factor is 1/4, while the remaining finite
product is at least one minus its sum, and
sum_(k>=1)R x_k=3/4.
Thus y b_i(y)b_j(y) has modulus at most
256 R^(i+j+1) on this circle.

Its Taylor coefficient of order k is consequently bounded by
256 R^(i+j+1-k). The first 2N actual and model moments agree,
and the rest differ by at most three. Summing the tail gives
an entrywise error at most
\[
3\cdot256\sum_{k\ge2N}R^{i+j+1-k}
=2304 R^{i+j+1-2N}.
\]
Since i+j+1<=2m-1, the symmetric row-sum bound proves the
second part of (20). Equation (21) is RLO-T2. QED.

This estimate uses the summability of the fixed dyadic positions.
It controls the Newton matrices directly; it does not multiply a
raw reading error by a potentially large divided-difference
matrix norm.

## RBP-T6. Actual dyadic relational banks are positive through square-root size

For a fixed anchor set
\[
N=N(a)=\left\lfloor\frac{\kappa\sqrt a}{2}\right\rfloor.
\]
If N>=128 and 1<=m<=floor(N/8), then
\[
\boxed{
G_m(a)\succeq\frac1{24\ell\,49^m}I,\qquad
\mathsf K_m(a)\succ0,\qquad
\mathfrak D_m(a)\ge
\left(\frac1{24\ell\,49^m}\right)^m.
}
\tag{25}
\]
A useful more general finite sufficient condition, with any
carrier N covered by (10), is
\[
\boxed{
55296\,\ell\,m\,49^m
\left(\frac32\right)^{-(2N-2m+1)}\le1.
}
\tag{26}
\]
Whenever (26) holds, the conclusions of (25) hold for that m.

**Proof.** Equation (26) makes the error in (20) at most half
its positive-model lower bound, proving the first inequality
in (25). Congruence (21) and the product of eigenvalue lower
bounds give the other two conclusions.

It remains to prove (26) throughout the explicit schedule.
Put R=3/2. The elementary inequality R^10>49 follows, for
example, from R^8>25 and R^2=9/4. For m<=N/8, the left side
of (26) is less than
\[
55296\ell m R^{-(2N-12m+1)}
\le6912\ell N R^{-N/2}.
\tag{27}
\]
The definition of the floor gives
\[
\ell<4+\log(2(N+1))<N/8\qquad(N\ge128).
\tag{28}
\]
For the second inequality, at N=128 use log(258)<9, which
follows from e>2. The difference
N/8-4-log(2(N+1)) increases thereafter.

Hence the last side of (27) is below
864 N^2 R^(-N/2). The function N^2 R^(-N/2) decreases
for N>=128, since log R>1/3. At N=128, R^64>2^32
because R^2>2, while 864(128)^2<2^24. Thus (27) is
strictly below 2^(-8)<1. This proves (26), and therefore
(25), throughout the stated range. QED.

All threshold estimates are written inequalities; no source
values or spectra have been numerically evaluated.

For each covered anchor, (25) establishes the first m_+(a)
members of RLO's actual fixed-sequence criterion:
\[
\mathsf K_m(a)
=\bigl(\mathsf K_a(2^{-i-1},2^{-j-1})\bigr)_{0\le i,j<m}.
\]
These are its initial prefixes, not only late clusters near
zero. Their growing range is explicit:
\[
m_+(a)\sim\frac{e^{-4}}{32\pi}\sqrt a
\quad(a\to\infty).
\tag{29}
\]
For anchors or sizes outside the stated schedule, (6), (18)
and (26) remain sufficient gates when their inequalities can
be established. This note does not evaluate them there.

## 2. What advances, and what remains open

| Object | Retained result | Result here |
|---|---|---|
| Unshifted H_N | SQR: coefficient e^(-9)/(2 pi) times sqrt(a) | Coefficient e^(-4)/(4 pi) times sqrt(a), with the same type of reference comparison |
| Both shifted forms A_N and B_N | Not supplied by SQR's unshifted theorem | Simultaneously positive throughout (10) |
| Native finite cut | Required A_N,B_N positivity as a premise | That premise is proved on the new growing range |
| Positive finite probability carrier | Finite positivity was needed | Constructed from the actual forms, matching c_0 through c_(2N-1) |
| Actual dyadic relational prefixes | RLO: each fixed size eventually positive, without a size threshold | Explicit square-root range and determinant margin (25) |
| Difference from a positive model | Must be retained | Scalar enclosure (17)--(18) and direct Newton bound (20) |
| All-size RH sign at one finite anchor | Open | Still open |

The finite carrier is useful because its probability and angular
structure are now consequences of proved source inequalities.
It is not an assumed positive model substituted for the unknown
whole source. The omitted moments enter explicitly in every
comparison with actual rational readings.

RLO proves that every dyadic relational prefix at one fixed
finite anchor is positive exactly when RH holds. This note
establishes an explicit growing but finite prefix at each
covered anchor. It does not establish all prefixes at that
same anchor.

Likewise, BRR's sufficient source-return criterion needs
positivity at arbitrarily large square-root widths. The
improved coefficient here is still one fixed width.
The newly proved comparisons do not satisfy that universal
width requirement.

PRW's full even support window L<=1/5 and FCP's two-moment
restricted window L<=111/200 remain unchanged. The anchor a,
carrier size N, relational size m and support half-width L
must not be identified.

## Dependencies, attribution and verification scope

- [RLO-1](RH_RELATIONAL_OBSERVER_CONFLUENCE.md): actual
  relational kernel, Newton congruence and fixed-sequence RH criterion.
- [SQR-1](RH_SQUARE_ROOT_BANKS.md): complete arithmetic
  wave remainder on the adapted analytic circle.
- [CWR-1](RH_CUT_WAVE_RETURN.md): wave coefficients and
  actual-source absolute budget.
- [QRI-1](RH_REFLECTION_RATIO_GATE.md): source/reference
  distinction and the first wave correction.
- [LOT-1](papers/native-layer-observer-transport/THEOREM.md):
  weighted finite cut and its two positive-form premises.
- [NIF-1](RH_INFINITY_CUT_FIXED_POINT.md): arcsine reference
  and finite carrier context.
- [BRR-1](RH_BOUNDARY_OBSERVER_RETURN.md): the remaining
  all-width return requirement.
- [PRW-1](RH_POLE_RECOVERY_WINDOW.md) and
  [FCP-1](RH_FINITE_CUT_POSITIVITY.md): unchanged support windows.

Toeplitz compression, polynomial filters, finite spectral
quadrature, Cauchy--Binet and analytic coefficient bounds are
standard tools. Their actual-source estimates, simultaneous
cut positivity and growing relational range are proved above.
No literature-wide novelty, optimal constant, all-order positive
completion or RH proof is claimed.

Only this note and its README entry change. The documentation
commit includes [skip ci] and is stacked on draft PR #118.
Verification is written argument review, document structure
and links, exact immutable-head content and PR metadata.
No scientific computation, numerical source evaluation,
independent agent, replay or workflow dispatch is performed.
The PR remains an unmerged draft.
