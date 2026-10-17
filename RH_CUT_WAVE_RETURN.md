# CWR-1 — Native cut-wave detection and accelerated observer return

17 October 2026. Written continuation of the same cut, reflection,
source and observer rules, based on Publications commit
9404e5e21eb180d78203b029a72655216ea8bbe1.

The retained two-reflection wave supplies an exact scalar detector
for the actual xi source. Combining its two seed channels gives a
wave return whose generating function is an explicit transformed
xi'/xi. Every off-line zero produces an uncancelled pole and hence
exponential growth along a subsequence of these readings.

The same wave polynomials improve the moving-observer approximation.
Returning any fixed reference bank now needs at most order
sqrt(rho) log rho directions at anchor ratio rho, improving MOR's
linear-depth construction. An explicit smaller growing-bank
schedule is again equivalent to RH.

These conclusions use the existing source product, source budget,
cut-wave recurrence and complex coefficient calculus. No new physical
law or external scientific assumption is introduced. The all-order
source bound that would conclude RH is displayed separately.

## Retained source, wave and observer conventions

Keep
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\qquad
a_0=3249/4,\quad a=\rho a_0,\quad \rho\ge1.
\]
For a positive-height zero beta+i gamma write
\[
\delta=\beta-\tfrac12,\qquad v_\mu=\gamma-i\delta,\qquad
\mu=v_\mu^2,\qquad y_a(\mu)=\frac a{a+\mu}.
\]
Multiplicity is denoted by nu_mu. Retain MOR's absolutely convergent
source functional
\[
q(a)=\sum_\mu\frac{\nu_\mu}{a+\mu},\qquad
\mathcal L_a[p]=\frac1{q(a)}
  \sum_\mu\frac{\nu_\mu}{a+\mu}p(y_a(\mu)).
\tag{1}
\]
Its real polynomial moments and banks are
\[
c_n(a)=\mathcal L_a[y^n],\qquad
H_N(a)=(c_{i+j}(a))_{0\le i,j<N}.
\tag{2}
\]
For all a>=a_0, MOR-T1 gives
\[
q(a)>0,\quad |y_a(\mu)|<1,\quad
\frac1{q(a)}\sum_\mu\frac{\nu_\mu}{|a+\mu|}\le2.
\tag{3}
\]
The weights and evaluations in (1) remain the actual, possibly complex
ones. Products in L_a are bilinear products, not absolute squares.

Put r=1-rho^-1 and f_rho=q(a)/q(a_0). The carried observer is
\[
(\mathcal O_\rho p)(y)
=(1-r y)^{-1/2}p\!\left(\frac{y}{\rho(1-r y)}\right),
\]
\[
\boxed{
f_\rho\,\mathcal L_a[(\mathcal O_\rho p)(\mathcal O_\rho h)]
=\mathcal L_{a_0}[ph].
}
\tag{4}
\]
The branch is one at y=0; O_1 is the identity. MOR-T2 proves (4)
without all-order positivity.

We also have 0<f_rho<=1. Indeed the proof of MOR-T1 establishes
a+Re(mu)>|Im(mu)| for every source term at a>=a_0. For
x=a+Re(mu), v=Im(mu),
\[
\frac d{da}\Re\frac1{a+\mu}
=\frac{v^2-x^2}{(x^2+v^2)^2}<0.
\]
Termwise differentiation is justified on compact anchor intervals
by the convergent inverse fourth-power zero tail. Hence q'(a)<0,
which proves the asserted normalization bound.

Wave step k, polynomial degree, carrier size, rho and the local
source displacement are kept distinct throughout.

## CWR-T1. Both cut channels give one exact reflection-even return

Use the existing LOT/LWP symbols
\[
P(y)=\begin{pmatrix}y&1-y\\y&1-y\end{pmatrix},\quad
Z=\operatorname{diag}(1,-1),\quad
g(y)=\operatorname{diag}(y,1-y),
\]
\[
\boxed{
W(y)=(2P(y)-I)Z
=\begin{pmatrix}2y-1&-2(1-y)\\2y&2y-1\end{pmatrix}.
}
\tag{5}
\]
These obey, as polynomial identities,
\[
P^2=P,\quad Z^2=I,\quad
W^TgW=g,\quad \det W=1,\quad
\frac{W+W^{-1}}2=(2y-1)I.
\tag{6}
\]
No source positivity is required for these algebraic statements.

Define the native wave polynomials by
\[
\tau_0(y)=1,\quad \tau_1(y)=2y-1,\quad
\tau_{k+1}(y)=2(2y-1)\tau_k(y)-\tau_{k-1}(y).
\]
Thus tau_k(y)=T_k(2y-1), as already derived in LWP-T1, and
\[
\frac{W^k+W^{-k}}2=\tau_k(y)I.
\tag{7}
\]
For the two constant seeds e_+=(1,0)^T and e_-=(0,1)^T, define
\[
\boxed{
r_k(a)=\sum_{\sigma\in\{+,-\}}
 \mathcal L_a[e_\sigma^T g W^k e_\sigma]
=\mathcal L_a[\tau_{|k|}(y)],\qquad k\in\mathbb Z.
}
\tag{8}
\]
Then r_0=1 and r_-k=r_k. The balanced wave Gram bank is
\[
\boxed{
\mathsf T_N(a)=(r_{j-i}(a))_{0\le i,j<N}.
}
\tag{9}
\]
It is exactly the sum of the two source Gram banks of
W^i e_+ and W^i e_-, respectively.

**Proof.** Direct multiplication proves (6). The recurrence obtained
from W+W^-1=2(2y-1)I gives (7) and shows that each diagonal entry
of W^k is tau_|k|. Since tr g=1, (8) follows. The source
functional is real on real polynomials by conjugation symmetry.
Finally (W^i)^T g W^j=g W^(j-i), proving (9). QED.

At a certified finite carrier of size M with both acquired forms
A_M,B_M positive, its finite model matches source moments through
degree 2M-1. The sum of the two seed correlations therefore equals
(8) for |k|<=2M-1, since its combined integrand is tau_|k| of
degree |k|. This gives one more combined step than the individual
seed horizon 2M-2 in LWP. No infinite positive completion is used
to define the unconditional readings in (8).

Changing the observer frame while carrying g, W and both seeds
preserves every pairing. Reversing the wave sends W to W^-1 and
leaves r_k unchanged. These are the precise reflection invariances.
Equation W^TgW=g is a bilinear conservation law; its interpretation
as conservation of a positive norm requires the source positivity.

## CWR-T2. Exact generating function of the actual balanced wave

For each source value let b=sqrt(a) and define
\[
\omega_a(\mu)=\frac{b+i v_\mu}{b-i v_\mu},\qquad
R_a=\frac{b+1/2}{b-1/2}.
\tag{10}
\]
Then
\[
\frac{\omega_a+\omega_a^{-1}}2=2y_a-1,\qquad
\tau_k(y_a)=\frac{\omega_a^k+\omega_a^{-k}}2,
\]
\[
\boxed{
\max(|\omega_a|,|\omega_a|^{-1})
=\sqrt{\frac{(b+|\delta|)^2+\gamma^2}
             {(b-|\delta|)^2+\gamma^2}}
\le R_a.
}
\tag{11}
\]
In particular |r_k(a)|<=2R_a^k, so its ordinary generating
function converges at least on |z|<R_a^-1.

Put u(z)=(1-z)/(1+z). On that initial disk,
\[
\boxed{
\mathfrak R_a(z):=\sum_{k=0}^\infty r_k(a)z^k
=\frac12+\frac{u(z)}{2q(a)}q(a u(z)^2)
=\frac12+
\frac{1}{4\sqrt a\,q(a)}
\frac{\xi'}{\xi}\!\left(\frac12+\sqrt a\,u(z)\right).
}
\tag{12}
\]
The right side gives its meromorphic continuation away from z=-1.

**Proof.** Substitution of y_a=a/(a+v_mu^2) proves the first
identity in (11); the wave recurrence proves the second.
The squared modulus of omega is
((b+delta)^2+gamma^2)/((b-delta)^2+gamma^2).
Taking the larger reciprocal and using |delta|<1/2 proves the
displayed bound. The total absolute source budget (3) then gives
the coefficient bound.

The generating identity for the recurrence is obtained by summing it:
\[
\sum_{k=0}^\infty T_k(w)z^k
=\frac{1-zw}{1-2zw+z^2}
=\frac12+\frac{1-z^2}{2(1-2zw+z^2)}.
\]
At w=2y_a-1 the denominator equals
\[
1-2z(2y_a-1)+z^2
=\frac{a(1-z)^2+\mu(1+z)^2}{a+\mu}.
\]
Insert this into (1). On a compact subdisk of |z|<R_a^-1,
the preceding absolute bounds justify interchanging the source
sum and the wave series. The result is the first equality in
(12). Differentiating X(w^2)=xi(1/2+w) gives
2w q(w^2)=xi'/xi(1/2+w), proving the second.
It then continues meromorphically by the retained entire source.
QED.

The first readings are, explicitly,
\[
r_0=1,\qquad r_1=2c_1-1,\qquad
r_2=8c_2-8c_1+1.
\tag{13}
\]
They are actual source combinations, without a fitted spectral model.

## CWR-T3. An uncancelled pole is an exponential observer detector

Define the reflection-invariant gain
\[
\boxed{
\mathfrak A(a)=
\sup_\mu\max(|\omega_a(\mu)|,|\omega_a(\mu)|^{-1})
\in[1,R_a].
}
\tag{14}
\]
Then the exact actual-source growth law is
\[
\boxed{
\limsup_{k\to\infty}|r_k(a)|^{1/k}=\mathfrak A(a).
}
\tag{15}
\]
Consequently, at any fixed a>=a_0, the following conditions
are equivalent:
\[
\boxed{
\begin{aligned}
&\mathrm{RH};\\
&\mathfrak A(a)=1;\\
&\limsup_{k\to\infty}|r_k(a)|^{1/k}\le1;\\
&|r_k(a)|\le1\quad\hbox{for every integer }k\ge0;\\
&\mathsf T_N(a)\succeq0\quad\hbox{for every }N.
\end{aligned}
}
\tag{16}
\]
Under RH all of the finite Toeplitz banks are in fact strictly
positive. If RH fails, there is a finite wave step k such that
\[
\boxed{
\det\begin{pmatrix}1&r_k(a)\\r_k(a)&1\end{pmatrix}
=1-r_k(a)^2<0.
}
\tag{17}
\]
This is a two-observer compression of the balanced wave pairing.
No actual negative reading is asserted to have been found.

**Proof.** The source pole q(t) at t=-mu has residue nu_mu.
Its two preimages in (12) are z=omega_a(mu) and omega_a(mu)^-1.
They are distinct and finite: gamma>0 excludes mu=0 and
omega=+1, while b>0 excludes omega=-1. At either preimage z_*,
d(a u^2)/dz=-4a u/(1+z)^2 is nonzero. The residue in (12) is
\[
\boxed{
\operatorname*{Res}_{z=z_*}\mathfrak R_a(z)
=-\frac{\nu_\mu(1+z_*)^2}{8a\,q(a)}\ne0.
}
\tag{18}
\]
Different mu cannot have the same preimage, since a u(z)^2
determines mu. Repeated zeros add their multiplicities. Thus
no such pole disappears by phase or reflection cancellation.

As gamma tends to infinity, omega tends to -1 in modulus one,
uniformly for |delta|<=1/2. There are finitely many zeros at
bounded height. If an off-line zero exists, A(a)>1 and the
largest gain is attained, so the nearest generating pole has
modulus 1/A(a)<1. There are no nearer singularities of (12).
If RH holds, all these poles lie on the unit circle; at least
one exists, and (12) has no singularity strictly inside it.
In both cases the Taylor radius is exactly 1/A(a).
The coefficient formula for the radius proves (15).

Equation (11) makes A(a)=1 equivalent to delta=0 for every
source zero. Under that condition all y_a lie in (0,1) and
all normalized source weights are positive. Since |tau_k(y)|<=1
there, |r_k|<=1. Conversely bounded readings imply the
subexponential condition in (16), so (15) gives RH.

Under RH, split each positive source weight equally between
omega and its conjugate inverse on the unit circle. These
weights sum to one and have moments r_k. Hence the Toeplitz
quadratic form is a positive sum of absolute squares of a
finite polynomial on the circle. Infinitely many distinct
source points make it strictly positive for every nonzero
polynomial. Conversely positivity of every T_N gives
1-r_k^2>=0 by its rows and columns 0,k. This completes (16).
If RH fails, (15) has growth greater than one, so some finite
k has |r_k|>1, proving (17). QED.

The limsup allows oscillating readings. It does not say that
every sufficiently late reading increases monotonically.

The retained verified-height input gives an additional honest
bound. With T=3*10^12, all possible off-line zeros have gamma>T,
so (11) gives
\[
\boxed{
0\le\log\mathfrak A(a)\le
\operatorname{artanh}
\left(\frac{\sqrt a}{a+T^2+1/4}\right).
}
\tag{19}
\]
To see this, write the logarithm of the gain as
artanh(2 sqrt(a)|delta|/(a+gamma^2+delta^2)).
Its argument increases with |delta|<=1/2 and decreases with
gamma, giving (19). This uses the same previously retained
finite-height result; no new zero verification is performed.
The upper bound is not an assertion that the gain equals one.

## CWR-T4. The renormalized infinite-anchor wave retains the horizontal defect

Let the source-specific maximal horizontal displacement be
\[
D_\xi=\sup_{\beta+i\gamma,\ \gamma>0}|\beta-\tfrac12|
\in[0,1/2].
\]
The actual wave detector satisfies the exact limit
\[
\boxed{
\lim_{a\to\infty}\sqrt a\,\log\mathfrak A(a)=2D_\xi.
}
\tag{20}
\]
Equivalently it is defined entirely by the source observer readings:
\[
\boxed{
D_\xi=
\frac12\lim_{a\to\infty}\sqrt a\,
\log\!\left(\limsup_{k\to\infty}|r_k(a)|^{1/k}\right).
}
\tag{21}
\]
This number is unchanged by reflecting delta to -delta or
reversing the wave, and RH is equivalent to D_xi=0.

**Proof.** For every fixed source zero, (11) gives
\[
\sqrt a\,
\log\max(|\omega_a|,|\omega_a|^{-1})
\longrightarrow2|\delta|.
\]
Taking a zero with |delta| arbitrarily close to the supremum
therefore proves that the lower limit in (20) is at least 2D_xi.
For the upper bound, (11) implies
\[
\mathfrak A(a)\le
\frac{\sqrt a+D_\xi}{\sqrt a-D_\xi}.
\]
Multiply the logarithm by sqrt(a) and let a increase. The
limit of the resulting bound is 2D_xi, proving equality.
No attainment of the supremum is needed. Equation (15) yields
(21), and absolute delta and reciprocal omega establish the
reflection statements. QED.

The normalization and the order of limits matter. The unscaled
A(a) tends to one whether or not RH holds. QRI-T4, meanwhile,
gives r_k(a)->0 at every fixed k>=1, because the arcsine average
of tau_k is zero. Thus
\[
\lim_{a\to\infty}\limsup_{k\to\infty}|r_k(a)|^{1/k}=1,\qquad
\limsup_{k\to\infty}
 \left|\lim_{a\to\infty}r_k(a)\right|^{1/k}=0.
\tag{22}
\]
The finite-depth infinity limit discards the very signal that
(21) retains. The theorem identifies the invariant and its
source meaning; it does not evaluate it as zero.

## CWR-T5. Wave-polynomial observers resolve the moving cut at square-root depth

For a real polynomial p(y)=sum_j p_j y^j put
\[
B_4(p)=\sum_j4^j|p_j|.
\]
For every rho>=1 the carried observer in (4) has a real wave
expansion
\[
\mathcal O_\rho p(y)=
\alpha_0+\sum_{k=1}^\infty\alpha_k\tau_k(y)
\]
on a region containing every actual y_a(mu). Define its finite
wave-polynomial observer by
\[
\boxed{
\mathcal P_{\rho,K}[p](y)=
\alpha_0+\sum_{k=1}^K\alpha_k\tau_k(y),\qquad K\ge0.
}
\tag{23}
\]
This is a linear map in p with output degree at most K.
Its coefficients are defined by the analytic wave expansion,
without replacing any source moment.

Set
\[
t_{\rho,K}=
12\sqrt\rho\,\exp\!\left(-\frac{3(K+1)}{4\sqrt\rho}\right).
\]
The two bounds on the entire actual source support are
\[
\boxed{
|\mathcal O_\rho p(y_a)|\le\sqrt\rho\,B_4(p),\qquad
|\mathcal O_\rho p(y_a)-\mathcal P_{\rho,K}[p](y_a)|
\le\sqrt\rho\,B_4(p)t_{\rho,K}.
}
\tag{24}
\]
In particular the exponent is controlled by K/sqrt(rho).

**Proof.** The wave variable v gives the two-sided parametrization
\[
y(v)=\frac12+\frac{v+v^{-1}}4.
\]
Let R=exp(1/sqrt(rho)). Its outer circle |v|=R maps to an
ellipse enclosing [0,1]. On and inside it,
\[
|y|\le\cosh^2(1/(2\sqrt\rho))
\le1+\frac9{16\rho}\le\frac{25}{16}.
\]
Indeed for 0<=s<=1/2, cosh s<3/2 and hence
sinh s<=3s/2. It follows that
\[
\Re(1-r y)\ge1-r|y|\ge\frac7{16\rho}>0,\qquad
\left|\frac{y}{\rho(1-r y)}\right|\le\frac{25}{7}<4.
\]
Thus O_rho p is holomorphic on a neighborhood of that closed
ellipse, with the same branch as (4), and its supremum there
is at most
\[
M_\rho(p)=2\sqrt\rho\,B_4(p).
\tag{25}
\]

Compose O_rho p with y(v). The result is holomorphic on a
neighborhood of the closed annulus R^-1<=|v|<=R and invariant
under v->v^-1. Its Laurent coefficients therefore satisfy
b_-k=b_k. Cauchy's coefficient formula on |v|=R gives
|b_k|<=M_rho(p)R^-k for k>=0. Since
tau_k(y(v))=(v^k+v^-k)/2, this proves the stated wave expansion,
with alpha_0=b_0 and alpha_k=2b_k for k>=1. Conjugation
symmetry makes these coefficients real.

For the actual support, (11) puts its wave variables inside the
smaller annulus with outer radius R_a. The following uniform
separation holds because sqrt(a_0)=57/2:
\[
\log R_a
=2\operatorname{artanh}(1/(2\sqrt a))
\le\frac2{\sqrt a}
=\frac4{57\sqrt\rho}
<\frac1{4\sqrt\rho}.
\]
Consequently R_a/R<=exp(-3/(4 sqrt(rho))). The wave-series
tail on that support is at most
\[
2M_\rho(p)\,
\frac{\exp(-3(K+1)/(4\sqrt\rho))}
     {1-\exp(-3/(4\sqrt\rho))}
\le\frac{32}{3}\rho\,B_4(p)
     \exp(-3(K+1)/(4\sqrt\rho)).
\]
For the last inequality use 1-exp(-x)>=x/2 when 0<=x<=1.
Increasing 32/3 to 12 proves the error in (24).

Finally, on the actual support |y_a|<1 and
y_a/[rho(1-r y_a)]=y_(a_0), whose modulus is also below one.
Thus |1-r y_a|>=1/rho and
|O_rho p(y_a)|<=sqrt(rho) sum_j|p_j|<=sqrt(rho)B_4(p).
This proves the first part of (24). QED.

The coefficients can equivalently be recovered by the ordinary
circle/Fourier integrals of O_rho p((1+cos theta)/2). These are
exact observer definitions, not numerical values computed in this
note. A future numerical use must enclose those coefficients and
the required actual moments. No extra root, exponential or
scientific execution engine is installed or run here.

The proof uses the existing cut-wave recurrence and the same
Cauchy coefficient principle used in QRI and NIF. The ellipse
is the image of that two-sided wave circle; it does not introduce
another physical geometry or an assumed positive zero measure.

## CWR-T6. Whole-bank return and a smaller sufficient RH schedule

For real polynomials p,h, let t=t_(rho,K). Then the finite
observer in (23) satisfies
\[
\boxed{
\left|\mathcal L_{a_0}[ph]
-f_\rho\mathcal L_a[
\mathcal P_{\rho,K}[p]\mathcal P_{\rho,K}[h]]\right|
\le 2f_\rho\rho(2t+t^2)B_4(p)B_4(h).
}
\tag{26}
\]

Put M=K+1, and let T_(rho,K,N) be the real M-by-N coefficient
matrix whose j-th column is P_(rho,K)[y^j] in the degree-below-M
monomial basis, for 0<=j<N. Set
\[
G_{\rho,K,N}=
f_\rho\,T_{\rho,K,N}^T H_M(a)T_{\rho,K,N},\qquad
S_N=\sum_{j=0}^{N-1}16^j=\frac{16^N-1}{15}.
\]
Then
\[
\boxed{
\|H_N(a_0)-G_{\rho,K,N}\|_2
\le 2f_\rho\rho(2t+t^2)S_N.
}
\tag{27}
\]
When t<=1, this is at most
\[
72f_\rho\rho^{3/2}S_N
 \exp\!\left(-\frac{3M}{4\sqrt\rho}\right).
\tag{28}
\]

For any fixed target size N and 0<epsilon<=1, the sufficient
depth
\[
\boxed{
M\ge
\left\lceil
\frac{4\sqrt\rho}{3}
\left(\frac32\log\rho+\log\frac{72S_N}{\epsilon}\right)
\right\rceil
}
\tag{29}
\]
makes the matrix return error at most epsilon. In particular
fixed-size, fixed-accuracy return needs at most order
sqrt(rho) log rho directions.

More generally, let rho_j tend to infinity and let M_j be
positive integers such that
\[
\frac{M_j}{\sqrt{\rho_j}}-2\log\rho_j\longrightarrow+\infty.
\tag{30}
\]
Positivity of the actual H_(M_j)(rho_j a_0) along that sequence
implies RH. Conversely RH makes all such actual banks positive.

A convenient explicit equivalent target is
\[
\boxed{
\mathrm{RH}
\ \Longleftrightarrow\
H_{\lceil4\sqrt m\,\log(2m)\rceil}(m a_0)\succeq0
\quad\hbox{for arbitrarily large integers }m.
}
\tag{31}
\]

**Proof.** Write F_p=O_rho p and P_p=P_(rho,K)[p]. Equation
(24) gives |F_p|<=sqrt(rho)B_4(p),
|P_p|<=sqrt(rho)B_4(p)(1+t), and
|F_p-P_p|<=sqrt(rho)B_4(p)t on the actual support.
Expanding F_p F_h-P_p P_h bounds its modulus by
rho B_4(p)B_4(h)(2t+t^2). The absolute source budget (3),
followed by the exact return (4), proves (26).

For monomials y^i,y^j, B_4 is 4^i,4^j. Hence for a real
coefficient vector x the error form is bounded by
\[
2f_\rho\rho(2t+t^2)
\left(\sum_{j<N}4^j|x_j|\right)^2
\le2f_\rho\rho(2t+t^2)S_N\|x\|_2^2.
\]
This proves (27). For t<=1, 2t+t^2<=3t, giving (28).
The choice (29) makes t<=epsilon/(6 rho S_N)<=1, and (28)
is at most f_rho epsilon<=epsilon.

Now suppose (30) holds and the larger banks are positive.
Fix any finite N and choose K_j=M_j-1. The returned G matrices
are positive semidefinite by congruence. Condition (30) forces
t_(rho_j,K_j)->0, and the bound in (28) is at most
\[
72S_N\exp\left[
-\frac34\left(\frac{M_j}{\sqrt{\rho_j}}-2\log\rho_j\right)
\right]\longrightarrow0.
\]
Thus H_N(a_0) is a limit of positive semidefinite matrices and
is positive semidefinite. This holds for every N. The retained
actual-source detector D3 then implies RH. Under RH, positive
weights and infinitely many distinct real support points make
every finite source bank strictly positive, proving the converse.

For the schedule in (31), with rho>=2, one has
\[
\exp(-3M/(4\sqrt\rho))\le(2\rho)^{-3},
\qquad t\le\tfrac32\rho^{-5/2}<1.
\]
The explicit returned error therefore satisfies
\[
\boxed{
\|H_N(a_0)-G_{\rho,M-1,N}\|_2
\le9S_N\rho^{-3/2}.
}
\tag{32}
\]
It tends to zero for each fixed N, proving (31) along any
unbounded integer subsequence. QED.

All moments through degree 2M-2 in the new actual bank are
required for this construction. An approximation to an observer
does not supply absent source moments. In any future numerical
implementation the coefficient matrix T and f_rho need their
own rigorous enclosures. Neither is evaluated here.

The explicit return (29)-(32) improves the sufficient resolution
from the previous MOR construction. It proves how positive
large-anchor banks would return the reference signs. It does
not assume or establish their positivity at the required sizes.

## CWR-T7. A hypothetical negative direction persists at square-root-log depth

Suppose a real polynomial p at the reference anchor satisfies
\[
\mathcal L_{a_0}[p^2]=-\delta<0,\qquad B=B_4(p)>0.
\]
Define
\[
\boxed{
M_\rho(p)=
\left\lceil
\frac{4\sqrt\rho}{3}
\left(\frac32\log\rho+\log\frac{288B^2}{\delta}\right)
\right\rceil,\qquad K_\rho=M_\rho(p)-1.
}
\tag{33}
\]
Then for every rho>=1,
\[
\boxed{
f_\rho\,
\mathcal L_a[\mathcal P_{\rho,K_\rho}[p]^2]
\le-\frac{3\delta}{4}<0.
}
\tag{34}
\]
Thus the actual H_(M_rho(p))(rho a_0) is indefinite. In
particular, if RH fails and n_-(a) denotes the first size whose
actual unshifted Hankel bank is not positive semidefinite,
\[
\boxed{
n_-(\rho a_0)\le
2\sqrt\rho\log\rho+
\frac43\sqrt\rho\log\frac{288B^2}{\delta}+1.
}
\tag{35}
\]
The constants depend on one hypothetical witness. No witness
or numerical value of its margin is asserted.

**Proof.** The source budget at a_0 gives
delta<=2(sum |p_j|)^2<=2B^2. Therefore the logarithm in (33)
is positive, and
\[
t_{\rho,K_\rho}\le
\frac{\delta}{24\rho B^2}\le\frac1{12\rho}<1.
\]
Apply (26) with h=p and use (28)'s same scalar estimate, with
B^2 in place of S_N. The error is at most
\[
72f_\rho\rho^{3/2}B^2
\left(\rho^{-3/2}\frac{\delta}{288B^2}\right)
\le\delta/4.
\]
The exact old reading is -delta, proving (34). Since f_rho>0,
the degree-below-M_rho polynomial is a negative direction of
H_M_rho. The constant direction has reading c_0=1, so the
bank is indefinite. The ceiling bound gives (35). QED.

Combined with QRI-T7, this sharpens the conditional visibility
interval:
\[
\boxed{
\liminf_{a\to\infty}\frac{n_-(a)}{\log a}
\ge\frac1{4\log10},\qquad
\limsup_{\rho\to\infty}
\frac{n_-(\rho a_0)}{\sqrt\rho\log\rho}\le2
\quad\hbox{if RH fails}.
}
\tag{36}
\]
The lower bound is unchanged; the old linear upper bound is
now replaced by the square-root-log bound. Both use the actual
same source and keep a hypothetical obstruction visible.

## What is now closed, and what would finish RH

The native balanced-observer detector, its uncancelled-pole
proof, the scaled infinite-anchor defect formula, and the
accelerated finite observer return are established by the
written derivations above. They all use the two-reflection
diagram rule and the already retained actual source.

The remaining source theorem can now be stated in any of
three exact forms:

| Route | Remaining actual-source assertion |
|---|---|
| Single balanced wave observer | abs(r_k(a_0))<=1 for every k |
| Infinite-depth scaled invariant | D_xi in (21) equals zero |
| Smaller moving banks | The actual banks in (31) are positive along an unbounded sequence |

Each assertion would finish RH through the proved implications.
None follows merely from the algebraic reflection identities.
For example P(2)^2=P(2), Z^2=I and W(2)^T g(2) W(2)=g(2),
but
\[
W(2)=\begin{pmatrix}3&2\\4&3\end{pmatrix},\qquad
g(2)=\operatorname{diag}(2,-1),
\]
and the reciprocal wave gains are 3+sqrt(8) and 3-sqrt(8).
This is an algebraic example of the same cut identities with
an indefinite pairing, not an actual xi source point.
A positive source norm is the condition that must be established.

The currently proved actual positive banks from QRI grow as
log a. The new sufficient schedule grows as sqrt(a/a_0) log(a/a_0).
The return theorem narrows this gap substantially but does not
close it. Declaring every formal cut wave unitary before proving
its positive source pairing would assume the missing result.

## Evidence, lineage and execution

| Result | Evidence |
|---|---|
| Balanced two-channel wave and carried reflection invariance | CWR-T1; retained LOT/LWP polynomial identities |
| Exact actual wave generating function | CWR-T2; recurrence and the full source product |
| Growth marker, bounded-wave RH equivalence and two-observer detection | CWR-T3; explicit nonzero residues |
| Scaled infinite-anchor defect equals maximal horizontal displacement | CWR-T4; exact wave gains and a supremum limit |
| Wave-polynomial observer approximation | CWR-T5; explicit nested wave circles and Cauchy bounds |
| Square-root-log bank return and equivalent smaller RH schedule | CWR-T6; actual source budget and congruence |
| Persistence of a hypothetical negative source direction | CWR-T7; explicit return error |
| The bounded-wave condition or required growing-bank positivity | Open |

Retained source lineage:
[ILC](papers/native-infinite-layer-cut/THEOREM.md),
[LOT](papers/native-layer-observer-transport/THEOREM.md),
[LWP](papers/native-layer-wave-propagation/THEOREM.md),
[finite detector](RH_FINITE_OBSERVER_DETECTOR.md),
[NIF](RH_INFINITY_CUT_FIXED_POINT.md),
[MOR](RH_MOVING_OBSERVER_RETURN.md), and
[QRI](RH_REFLECTION_RATIO_GATE.md).

The circle parametrization, Chebyshev recurrence, Laurent/Cauchy
coefficient formula, Taylor-radius theorem and congruence calculus
are standard mathematics already used in this source track.
Their required steps and their precise actual-source application
are given here. No literature-wide novelty or optimality claim
is made for the reformulations or constants. The existing zero
product, critical strip and verified-height inputs retain their
attribution in MOR; no new external scientific input is added.

Only this note and its README entry are added. No local scientific
computation, independent agent, numerical source evaluation or
scientific replay is performed. No workflow is changed. These are
written proofs, not new machine-checked certificates. The commit
contains [skip ci], and the change remains a draft stacked on
PR #99 without a merge.
