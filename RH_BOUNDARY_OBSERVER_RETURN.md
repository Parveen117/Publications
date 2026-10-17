# BRR-1 — Boundary matching and square-root observer return

17 October 2026. Written continuation of the actual source and native
cut-wave programme, based on Publications commit
0304aeecc6e0eedb64d19d804252e05afcba9523.

The local zero chart can be recovered exactly from the infinite-anchor
wave boundary, provided the observation coordinate and normalization
move with the anchor. Its unscaled boundary reading tends to zero;
the normalized ratio retains the entire local source germ.

That matching suggests a sharper finite observer. A polynomial
resolvent filter made from the existing cut-wave polynomials has a
uniform relative error on the actual complex source. Integrating
these filters gives a polynomial approximation to the square-root
weight in the exact transported observer. The source weights cancel
the observer's large size in the error estimate.

As a result, fixed-bank, fixed-accuracy return now needs at most order
sqrt(rho) polynomial directions at anchor ratio rho. This removes the
log(rho) factor from the previous CWR sufficient return bound. If RH
fails, the first negative unshifted bank must occur at square-root
depth, between explicit lower and witness-dependent upper constants.

**Status: RH remains open in this programme.** SQR supplies positivity
up to one fixed multiple of sqrt(a). The new return theorem does not
prove positivity at every required multiple. All results here are
written analytic derivations, not new machine-checked certificates.

## Retained source and coordinates

Let
\[
a_0=3249/4,\qquad b_0=57/2,\qquad
a=\rho a_0,\quad \rho\ge1,\quad b=\sqrt a,\quad
r=1-\rho^{-1},\quad f_\rho=\frac{q(a)}{q(a_0)}.
\]
Retain the source
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\qquad
\mu=(\gamma-i\delta)^2,\quad \delta=\beta-\tfrac12,
\]
\[
y_a(\mu)=\frac a{a+\mu},\qquad
w_\mu(a)=\frac{\nu_\mu}{q(a)(a+\mu)},\qquad
\mathcal L_a[p]=\sum_\mu w_\mu(a)p(y_a(\mu)).
\tag{1}
\]
Here nu_mu is the retained multiplicity. MOR gives, for every
a>=a_0,
\[
q(a)>0,\quad |y_a(\mu)|<1,\quad
\sum_\mu|w_\mu(a)|\le2,\quad
U_\mu:=\Re\mu>0,\quad |\Im\mu|<a_0+U_\mu.
\tag{2}
\]
The last inequality follows directly from the same low-zero and
high-zero estimates in MOR-T1. No new zero-height input is used.

The moments and actual unshifted banks are
\[
c_n(a)=\mathcal L_a[y^n],\qquad
H_N(a)=(c_{i+j}(a))_{0\le i,j<N}.
\]
They are real by conjugation. Products inside L_a remain bilinear;
the complex source weights are not assumed positive.

The retained moving observer is
\[
(\mathcal O_\rho p)(y)
=(1-r y)^{-1/2}p\!\left(\frac{y}{\rho(1-r y)}\right),
\]
\[
\boxed{
f_\rho\mathcal L_a[(\mathcal O_\rho p)(\mathcal O_\rho h)]
=\mathcal L_{a_0}[ph].
}
\tag{3}
\]
The square root is the analytic branch equal to one at y=0.
For rho=1 this is the identity map.

Use the same native wave polynomials
\[
\tau_K(y)=T_K(2y-1),\qquad
R_a=\frac{b+1/2}{b-1/2}.
\]
CWR-T2 gives, on every actual source point,
\[
|\tau_K(y_a(\mu))|\le R_a^K,\qquad
\log R_a\le\frac2b=\frac4{57\sqrt\rho}<\frac1{4\sqrt\rho}.
\tag{4}
\]
These estimates include any hypothetical off-line source values.

| Coordinate | Role |
|---|---|
| lambda_0 | Local displacement at the fixed reference anchor a_0 |
| Lambda=rho-1 | Outward source-anchor displacement |
| z | Complex wave-generating coordinate |
| K | Degree parameter of the polynomial resolvent filter |
| t>=0 | Auxiliary integration variable in that filter |
| N, M | Reference bank size and returning bank size |

In particular t is not physical time, and lambda_0 is not silently
identified with the reciprocal anchor coordinate or a thermodynamic
lambda component.

## BRR-T1. The local zero chart survives in a shrinking boundary window

Write
\[
\mathfrak F(w)=wq(w^2),\qquad
C_a(z)=1+2\sum_{k\ge1}r_k(a)z^k
=\frac{\mathfrak F(b(1-z)/(1+z))}{\mathfrak F(b)},
\]
where r_k=L_a[tau_k], and use its actual-source continuation where
needed. For |lambda_0|<=1/2 define
\[
w(\lambda_0)=b_0\sqrt{1+\lambda_0},\qquad
z_b(\lambda_0)=\frac{b-w(\lambda_0)}{b+w(\lambda_0)},\qquad
F_{a_0}(\lambda_0)=\frac{q(a_0(1+\lambda_0))}{q(a_0)}.
\tag{5}
\]
All roots are the branches continued from lambda_0=0. Then, for
every finite b>=b_0,
\[
\boxed{
F_{a_0}(\lambda_0)
=\frac1{\sqrt{1+\lambda_0}}\,
\frac{C_a(z_b(\lambda_0))}{C_a(z_b(0))}.
}
\tag{6}
\]
The denominator is nonzero. This recovers every derivative of the
reference germ, not only its value F_(a_0)(0)=1.

Let ell=log(b/(2pi)). Uniformly on this lambda_0 disk as b tends
to infinity,
\[
\boxed{
1-z_b(\lambda_0)=\frac{2w(\lambda_0)}{b+w(\lambda_0)}
=O(b^{-1}),\qquad
\ell C_a(z_b(\lambda_0))\longrightarrow4\mathfrak F(w(\lambda_0)).
}
\tag{7}
\]
In particular C_a(z_b(lambda_0))->0, whereas C_a(z)->1 at each
fixed |z|<1. The limit is not uniform up to this moving boundary.

**Proof.** The Cayley identity gives
b(1-z_b)/(1+z_b)=w(lambda_0). Consequently
\[
C_a(z_b(\lambda_0))
=\frac{\mathfrak F(w(\lambda_0))}{\mathfrak F(b)}.
\]
Dividing by its lambda_0=0 value and using
mathfrak F(w)=w q(w^2) proves (6).
Its denominator is mathfrak F(b_0)/mathfrak F(b)>0 by (2).
The whole lambda_0 disk has Re(w)>2, so the source expressions
are analytic there and the identity also identifies all derivatives.

SQR-T1 gives mathfrak F(b)=ell/4+O(b^(-1)), proving (7)
uniformly since the numerator is bounded on the fixed disk.
For fixed |z|<1, u=(1-z)/(1+z) has Re(u)>0. Applying that same
source bound at bu and at b gives
C_a(z)=1+log(u)/ell+o(1), proving the fixed-z limit. QED.

Equation (6) is a transported, normalized matching of charts.
It does not set the coordinate values zero and infinity equal.
Taking only the two unscaled limits in its ratio would replace
a recoverable germ by the indeterminate expression 0/0.

The entire arithmetic source remains in mathfrak F(w(lambda_0)),
including its prime-power terms. No new invariant is assigned an
assumed value at infinity. The following construction makes the
observer return finite and quantitative.

## BRR-T2. Source-weighted relative errors have an anchor-independent budget

Put
\[
x_\rho(y)=1-r y,\qquad x_\mu=x_\rho(y_a(\mu)).
\]
Then
\[
\boxed{
x_\mu=\frac{a_0+\mu}{a+\mu},\qquad
\frac{y_a(\mu)}{\rho x_\mu}=y_{a_0}(\mu),\qquad
\Re x_\mu\ge\rho^{-1},\qquad |\arg x_\mu|<\pi/4.
}
\tag{8}
\]
Moreover the exact absolute-weight identity gives
\[
\boxed{
f_\rho\sum_\mu |w_\mu(a)|\,|x_\mu|^{-1}
=\frac1{q(a_0)}\sum_\mu\frac{\nu_\mu}{|a_0+\mu|}
\le2.
}
\tag{9}
\]

**Proof.** The first two identities are direct substitution.
Since |y_a|<1 and r>=0, Re(1-r y_a)>=1-r=rho^(-1).

Write U=Re(mu)>0 and V=Im(mu). Direct multiplication by
the conjugate denominator gives
\[
\frac{|\Im x_\mu|}{\Re x_\mu}
=\frac{|V|(a-a_0)}
{(a_0+U)(a+U)+V^2}
<1
\]
when V is nonzero, by |V|<a_0+U; for V=0 the ratio is zero.
This proves the argument bound. Finally
|a+mu| |x_mu|=|a_0+mu|, so the left side of (9) equals
the middle expression term by term. MOR's absolute budget at
a_0 proves its bound. QED.

This budget is the useful cancellation: the transported observer
has size proportional to |x_mu|^(-1/2), but its product is
summed against source weights containing exactly the compensating
factor. A uniform error bound that first discards this dependence
can introduce avoidable anchor growth.

## BRR-T3. The native wave gives a relative polynomial resolvent filter

For rho>1, integer K>=1 and real t>=0, define
\[
d_{\rho,t}=\frac{2(1+t^2)}r-1>1,\qquad
\boxed{
\mathcal B_{\rho,K,t}(y)=
\frac{1-\tau_K(y)/T_K(d_{\rho,t})}{1-r y+t^2}.
}
\tag{10}
\]
Despite its displayed quotient, B_(rho,K,t) is a real polynomial
of degree at most K-1. Set
\[
\varepsilon_{\rho,K}
=2\exp\!\left(-\frac{7K}{4\sqrt\rho}\right).
\tag{11}
\]
At every actual source point,
\[
\boxed{
\left|1-(x_\mu+t^2)\mathcal B_{\rho,K,t}(y_a(\mu))\right|
\le\varepsilon_{\rho,K}.
}
\tag{12}
\]
The bound is uniform over t>=0 and the full complex source.

**Proof.** The denominator in (10) vanishes at
y=(1+t^2)/r. At that point 2y-1=d_(rho,t), so the polynomial
numerator also vanishes. Polynomial division proves the degree
and reality assertion.

For d>1 the wave recurrence gives
T_K(d)=cosh(K arcosh(d))>=exp(K arcosh(d))/2. Since
\[
d_{\rho,t}\ge d_{\rho,0}=\frac{\rho+1}{\rho-1},\qquad
\operatorname{arcosh}d_{\rho,0}
=\log\frac{\sqrt\rho+1}{\sqrt\rho-1}
=2\operatorname{artanh}(\rho^{-1/2})
\ge\frac2{\sqrt\rho},
\]
equation (4) yields
\[
\left|\frac{\tau_K(y_a)}{T_K(d_{\rho,t})}\right|
\le2\exp\left[K\log R_a-K\operatorname{arcosh}d_{\rho,t}\right]
\le2\exp\left(-\frac{7K}{4\sqrt\rho}\right).
\]
This ratio is exactly the residual in (12). QED.

For rho=1 define B_(1,K,t)=1/(1+t^2), a constant polynomial.
Its resolvent residual is zero, and all later estimates continue
to hold. Thus no singular limiting convention is needed at the
unchanged anchor.

The Chebyshev representation used here is the standard identity
T_K((z+z^(-1))/2)=(z^K+z^(-K))/2, recorded in
[NIST DLMF 18.5.1](https://dlmf.nist.gov/18.5#E1) and already
derived from the native wave in CWR. The divided polynomial
and its source-specific residual bound are given explicitly above.

## BRR-T4. Integrating the filters preserves finite polynomial degree

Define
\[
\boxed{
\mathcal S_{\rho,K}(y)=\frac2\pi
\int_0^\infty\mathcal B_{\rho,K,t}(y)\,dt.
}
\tag{13}
\]
The integral defines a real polynomial of degree at most K-1.
At every actual source point,
\[
\boxed{
\left|\mathcal S_{\rho,K}(y_a)-x_\mu^{-1/2}\right|
\le2\varepsilon_{\rho,K}|x_\mu|^{-1/2}.
}
\tag{14}
\]
For rho=1 the polynomial is exactly one.

**Proof.** First justify the polynomial integral without a source
positivity assumption. At each fixed y in [0,1],
|tau_K(y)|<=1 and T_K(d_(rho,t))>=1. Hence
\[
|\mathcal B_{\rho,K,t}(y)|
\le\frac2{\rho^{-1}+t^2}.
\]
Choose any K distinct points in [0,1]. Their values determine the
K polynomial coefficients by a fixed invertible interpolation
matrix. Each value is absolutely integrable in t by the displayed
bound, so every coefficient is absolutely integrable. Integrating
the coefficients defines the polynomial in (13); its integral
representation then holds at every complex y as well.

For Re(x)>0 the elementary resolvent integral is
\[
x^{-1/2}=\frac2\pi\int_0^\infty\frac{dt}{x+t^2}.
\tag{15}
\]
For positive real x, substitute t=sqrt(x)u. Holomorphic continuation
on the right half-plane proves the principal-branch formula;
the integral converges locally uniformly there.

Subtract (15) from (13) at a source point and use (12):
\[
\left|\mathcal S_{\rho,K}(y_a)-x_\mu^{-1/2}\right|
\le\frac{2\varepsilon_{\rho,K}}\pi
\int_0^\infty\frac{dt}{|x_\mu+t^2|}
\le\frac{\varepsilon_{\rho,K}}{\sqrt{\Re x_\mu}}.
\]
By (8), |x_mu|/Re(x_mu)<sqrt(2), so the last quantity is
at most 2 epsilon_(rho,K) |x_mu|^(-1/2). This proves (14).
For rho=1, (15) at x=1 gives S=1. QED.

The auxiliary integration variable builds polynomial coefficients.
It introduces no new actual source measure and changes none of the
source moments. No coefficient integral is numerically evaluated
in this continuation.

## BRR-T5. Relative polynomial observers return whole banks at square-root depth

For a real polynomial p of degree d, put
\[
B_4(p)=\sum_{j=0}^d4^j|p_j|,\qquad
\boxed{
\mathcal P_{\rho,K}[p](y)
=\mathcal S_{\rho,K}(y)\,
p\!\left(\frac{y\,\mathcal B_{\rho,K,0}(y)}{\rho}\right).
}
\tag{16}
\]
This map is linear in p and its output has degree at most
(d+1)K-1. If epsilon_(rho,K)<=1, then on the actual source
\[
\boxed{
|\mathcal P_{\rho,K}[p](y_a)-\mathcal O_\rho p(y_a)|
\le3\varepsilon_{\rho,K}B_4(p)|x_\mu|^{-1/2}.
}
\tag{17}
\]
Under the same condition, for any two real polynomials p,h,
\[
\boxed{
\left|\mathcal L_{a_0}[ph]
-f_\rho\mathcal L_a[
\mathcal P_{\rho,K}[p]\mathcal P_{\rho,K}[h]]\right|
\le30\varepsilon_{\rho,K}B_4(p)B_4(h).
}
\tag{18}
\]

For a reference bank of size N and epsilon_(rho,K)<=1, take M=NK. Let V_(rho,K,N)
be the real M-by-N coefficient matrix whose j-th column is
P_(rho,K)[y^j], for 0<=j<N. Define
\[
G_{\rho,K,N}=
f_\rho V_{\rho,K,N}^TH_{NK}(a)V_{\rho,K,N},\qquad
S_N=\sum_{j=0}^{N-1}16^j=\frac{16^N-1}{15}.
\]
Then
\[
\boxed{
\|H_N(a_0)-G_{\rho,K,N}\|_2
\le60S_N\exp\!\left(-\frac{7K}{4\sqrt\rho}\right).
}
\tag{19}
\]
In particular, for any fixed N and 0<eta<=1, choose
\[
\boxed{
K=\left\lceil\frac{4\sqrt\rho}{7}
\log\frac{60S_N}{\eta}\right\rceil,\qquad M=NK.
}
\tag{20}
\]
The error is at most eta. For fixed N and eta, M is
O(sqrt(rho)), with the explicit constants in (20).

**Proof.** Let epsilon=epsilon_(rho,K), x=x_mu and
v=y_(a_0)(mu)=y_a/(rho x). Equation (12) at t=0 says
x B_(rho,K,0)(y_a)=1+e_mu, with |e_mu|<=epsilon.
Thus the argument of p in (16) is v(1+e_mu), with |v|<1.
For epsilon<=1,
\[
|p(v(1+e_\mu))|\le\sum_j2^j|p_j|\le B_4(p),
\]
\[
|p(v(1+e_\mu))-p(v)|
\le\varepsilon\sum_{j\ge1}j2^{j-1}|p_j|
\le\varepsilon B_4(p).
\]
The first difference estimate follows by the finite geometric
identity for (1+e_mu)^j-1. The last inequality uses
j2^(j-1)<=4^j.

Apply (14) to the root factor and the preceding difference
estimate to p. Their sum gives (17).
Also |O_rho p(y_a)|<=B_4(p)|x|^(-1/2).
Writing P_p=O_p+E_p and expanding P_p P_h-O_p O_h therefore gives
\[
|P_pP_h-O_pO_h|
\le(6\varepsilon+9\varepsilon^2)
B_4(p)B_4(h)|x|^{-1}
\le15\varepsilon B_4(p)B_4(h)|x|^{-1}.
\]
Summing against f_rho times the absolute source weights,
using (9), and then using the exact bilinear return (3),
proves (18). This is an absolute error estimate on the original
bilinear products, not a replacement by absolute squares.

The degree of S is at most K-1 and the degree of y B_(rho,K,0)
is at most K. This proves the degree claim and the size M=NK.
For a real vector alpha in R^N, apply (18) to the polynomial
sum_(j<N) alpha_j y^j on both sides:
\[
|\alpha^T(H_N-G_{\rho,K,N})\alpha|
\le30\varepsilon
\left(\sum_{j<N}4^j|\alpha_j|\right)^2
\le30\varepsilon S_N\|\alpha\|_2^2.
\]
This proves (19). The choice (20) makes
epsilon<=eta/(30S_N)<=1 and makes that error at most eta.
At rho=1 the polynomial observer is exactly p, with zero
return error; padding its coefficient vector gives the same
M=NK formulation. QED.

The return is a genuine congruence of an actual moment bank.
If H_(NK)(a) is positive semidefinite, so is G_(rho,K,N).
The construction does not presume that positivity.

Compared with CWR's uniform absolute observer estimate, (18)
retains the x_mu dependence until the source weights are summed.
This removes the anchor prefactor responsible for the old
log(rho) in the fixed-bank, fixed-accuracy sufficient degree.
The dependence on reference size N in (20) is explicit and is
not claimed optimal.

## BRR-T6. A hypothetical negative direction returns at square-root depth

Suppose a real reference polynomial p of degree d satisfies
\[
\mathcal L_{a_0}[p^2]=-\Delta<0,\qquad B=B_4(p)>0.
\]
For every rho>=1 define
\[
\boxed{
K_\rho(p)=
\left\lceil\frac{4\sqrt\rho}{7}
\log\frac{240B^2}{\Delta}\right\rceil,\qquad
M_\rho(p)=(d+1)K_\rho(p).
}
\tag{21}
\]
Then the polynomial observer in (16) satisfies
\[
\boxed{
f_\rho\mathcal L_a[
\mathcal P_{\rho,K_\rho(p)}[p]^2]
\le-\frac{3\Delta}{4}<0.
}
\tag{22}
\]
In particular the actual bank H_(M_rho(p))(rho a_0) is indefinite.

If RH fails, let n_-(a) denote the first size for which
H_N(a) is not positive semidefinite. Then
\[
\boxed{
n_-(\rho a_0)
\le C_p\sqrt\rho+d+1,\qquad
C_p=\frac{4(d+1)}7\log\frac{240B^2}{\Delta}.
}
\tag{23}
\]
These constants depend on a hypothetical witness. None of p,
Delta or C_p is asserted to have been found for the actual source.

**Proof.** The retained absolute source budget gives
Delta<=2(sum |p_j|)^2<=2B^2. Hence the logarithm in (21)
is positive, and
\[
\varepsilon_{\rho,K_\rho(p)}
\le\frac{\Delta}{120B^2}\le\frac1{60}<1.
\]
Applying (18) with h=p bounds the return error by Delta/4,
proving (22). The returned polynomial has degree below
M_rho(p), so it is a negative direction of the displayed bank.
Its constant direction has reading c_0=1; thus it is indefinite.
The ceiling inequality proves (23). At rho=1 the exact return
already has the original negative value. QED.

Combine this result with SQR-T5 and the retained finite detector.
If RH fails, the same actual source has
\[
\boxed{
0<\kappa\sqrt{a_0}
\le\liminf_{\rho\to\infty}\frac{n_-(\rho a_0)}{\sqrt\rho}
\le\limsup_{\rho\to\infty}\frac{n_-(\rho a_0)}{\sqrt\rho}
\le C_p<\infty,\qquad
\kappa=\frac{e^{-9}}{2\pi}.
}
\tag{24}
\]
Thus any hypothetical first failure has square-root growth in
anchor, bounded between positive constants. This improves
CWR's square-root-log upper bound and uses SQR's square-root
lower bound. It supplies neither a negative witness nor an
all-order positivity proof.

## BRR-T7. A sharper sufficient schedule and an equivalent family of widths

Let rho_j tend to infinity and let M_j be positive integers with
\[
\boxed{\frac{M_j}{\sqrt{\rho_j}}\longrightarrow+\infty.}
\tag{25}
\]
If H_(M_j)(rho_j a_0) is positive semidefinite along this sequence,
then RH holds. Conversely RH makes all these actual banks
positive definite.

For example, the smaller explicit schedule
\[
\boxed{
\mathrm{RH}\quad\Longleftrightarrow\quad
H_{\left\lceil\sqrt m\,\log(2+\log(m+1))\right\rceil}(m a_0)
\succeq0
\quad\text{for arbitrarily large integers }m
}
\tag{26}
\]
is valid. More generally, any multiplier tending to infinity
can replace log(2+log(m+1)).

There is also the exact fixed-width family criterion
\[
\boxed{
\mathrm{RH}\quad\Longleftrightarrow\quad
\text{for every }L>0,\quad
H_{\lceil L\sqrt m\rceil}(m a_0)\succeq0
\text{ for arbitrarily large integers }m.
}
\tag{27}
\]
The quantifier "for every L" is essential to this stated criterion.

**Proof.** Fix a finite reference size N. For all sufficiently
large j set K_j=floor(M_j/N)>=1. Then
K_j/sqrt(rho_j)->infinity by (25), and N K_j<=M_j.
Positivity of the supplied bank implies positivity of its
leading H_(N K_j) submatrix. Its returned congruence G_j from
(19) is positive semidefinite. The error bound in (19) tends
to zero, so H_N(a_0) is a limit of positive semidefinite matrices.
This works for every fixed N. The retained finite detector
therefore implies RH. Under RH, the actual positive weights
and infinitely many distinct real source points make every
finite H_N positive definite, proving the converse.

For the schedule in (26), the size divided by sqrt(m) tends
to infinity. An unbounded sequence of indices with the stated
positivity therefore satisfies (25). This proves (26) and its
generalization.

To prove (27), its RH direction is immediate. If RH fails,
the detector supplies one finite negative polynomial p.
Choose a constant L>C_p from (23). For all sufficiently
large m,
\[
\lceil L\sqrt m\rceil\ge C_p\sqrt m+d+1
\ge n_-(m a_0).
\]
The bank at that size contains a leading submatrix that is not
positive semidefinite, so the larger bank is not positive semidefinite.
There cannot be arbitrarily large positive banks for this L.
This proves the contrapositive and completes the equivalence.
QED.

## What matching zero and infinity has achieved

A fixed-depth observer at increasing anchor sees the arcsine
limit and can lose finite-source information. Equation (6)
retains the entire zero-chart germ in a moving wave window of
width O(1/sqrt(a)). Equations (16)-(20) give finite polynomial
observers that return each fixed reference bank at that same
square-root scale, for any fixed requested accuracy.

The sharpened conditional result (24) is a quantitative account
of a possible hidden sign: it cannot occur in SQR's covered
small-width banks, but one hypothetical witness would reappear
at another finite square-root width. The moving observer does
not create or erase that sign.

| Quantity | Previous retained estimate | This continuation |
|---|---|---|
| Positive actual bank sizes | SQR: one explicit multiple of sqrt(a) | Unchanged |
| Fixed-bank, fixed-accuracy return | CWR: O(sqrt(rho) log(rho)) | O(sqrt(rho)), with explicit N and accuracy dependence |
| First failing bank if RH fails | Between sqrt(rho) and sqrt(rho) log(rho) scales | Between positive constant multiples of sqrt(rho) |
| Sufficient growing-bank condition | M/sqrt(rho)-2log(rho) tends to infinity | M/sqrt(rho) tends to infinity |

The remaining actual-source task can now be stated without the
old mandatory logarithmic factor: establish the positivity in
(25), or establish the family of widths in (27), or prove all
fixed-anchor Hankel signs directly.

SQR supplies N<=max(1,floor(kappa sqrt(a))). Its ratio
N/sqrt(rho) tends to the single finite constant kappa sqrt(a_0).
It therefore does not satisfy (25), and it does not establish
(27) for every L. Removing an avoidable logarithm in a sufficient
return estimate is not removal of this positivity obligation.

## Evidence, attribution and execution

| Statement | Evidence or status |
|---|---|
| Exact zero-chart recovery from the moving boundary | BRR-T1; CWR generating identity and SQR source asymptotic |
| Relative-error budget independent of anchor | BRR-T2; retained source geometry and exact weight cancellation |
| Polynomial resolvent and square-root filters | BRR-T3 to T4; wave recurrence, polynomial division and a derived resolvent integral |
| Square-root whole-bank return | BRR-T5; explicit real congruence and source-weighted error estimate |
| Conditional square-root first-failure bounds | BRR-T6 plus SQR-T5; no actual witness asserted |
| Sharper RH-equivalent schedules | BRR-T7; all-order premise remains open |
| New numerical certificate or physical invariant value | Not supplied |

Source lineage:
[finite detector](RH_FINITE_OBSERVER_DETECTOR.md),
[NIF](RH_INFINITY_CUT_FIXED_POINT.md),
[MOR](RH_MOVING_OBSERVER_RETURN.md),
[CWR](RH_CUT_WAVE_RETURN.md),
[RFT](RH_REFLECTION_TRANSFER.md), and
[SQR](RH_SQUARE_ROOT_BANKS.md).

The polynomial residual construction, the resolvent integral,
congruence and norm estimates use classical tools applied to
the retained native wave and source. Their required steps are
given here. No literature-wide novelty or optimality claim is made.

Only this note and its README entry are added. No local scientific
computation, independent agent, numerical source evaluation or
scientific replay is performed. In a future numerical implementation,
the polynomial coefficients and source normalization would require
their own rigorous enclosures; none is evaluated here. Existing data,
engines, certificates and workflows remain unchanged. These are
written proofs, not new machine-checked certificates. The commit
includes [skip ci]. The change remains a draft stacked on PR #103
without a merge.
