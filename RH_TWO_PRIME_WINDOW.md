# TPW-1 — Quadratic kernel cancellation and a two-prime window

17 October 2026. Written continuation at Publications commit
8d08991385f0511c6064960284e4d768617159ec, following
[BPC-1](RH_BOUNDARY_PRIME_CAPS.md).

This continuation crosses the prime-3 threshold while retaining
exactly the same two neutral moments as FPW and BPC. The resulting
written estimate is
\[
\boxed{
0<L\le\frac{11}{20}
\quad\Longrightarrow\quad
\mathscr Q_\xi(g,g)\ge\frac1{500}\|g\|_2^2,
\qquad g\in\mathcal W_L .
}
\]
At L=11/20 both primes 2 and 3 can contribute. No further prime
power is active. BPC's stronger margin 1/80 for L<=1/2 remains
available in its earlier range.

Two ingredients make the extension possible. A quadratic lower
bound for the Gamma kernel contributes an exact scalar on even
mean-zero tests: its two correlation moments vanish. A joint
prime estimate then treats the small prime-3 edge caps separately
from the rest of the window. The worst prime-2 midpoint estimate
need not be added to a worst prime-3 estimate at every position.

A general polynomial-kernel identity records the finite moment
residual that remains at higher degree. Only the quadratic case
cancels under the present conditions; no additional moments are
silently imposed.

**Status: RH remains unproved.** Positivity on this same family
for every finite L is RH-equivalent by FPW-T1. This note proves
the stated two-prime range, not that all-width condition. These
are written analytic proofs, not machine-checked source
certificates. No scientific computation, source evaluation,
independent agent or numerical replay is performed.

## Retained family and exact source

Let I=(-L,L), extend functions by zero, and retain
\[
\begin{aligned}
\mathcal W_L=\{g\in C_c^\infty(I,\mathbb R):\;&g(-x)=g(x),\\
&\int g(x)\,dx=0,\quad
\int g(x)\cosh(x/2)\,dx=0\}.
\end{aligned}
\tag{1}
\]
The two conditions are fixed for every L.

Write
\[
N_g=\|g\|_2^2,\qquad
C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,\qquad
a_n=\frac{\Lambda(n)}{\sqrt n},\quad u_n=\log n,
\]
and let N_L denote the finite set of prime powers n<e^(2L).
On W_L the pole square vanishes and the actual source form is
\[
\boxed{
\mathscr Q_\xi(g,g)=\mathscr A[g]
-\sum_{n\in\mathcal N_L}a_n C_g(u_n).
}
\tag{2}
\]
Retain the exact Gamma identity
\[
\begin{gathered}
\mathscr A[g]=\frac{\Omega_0}{2}N_g+
\int_0^\infty k(u)\bigl(N_g-C_g(u)\bigr)\,du,\\
\Omega_0=-\gamma_E-\pi/2-3\log2-\log\pi,\qquad
k(u)=\frac{e^{u/2}}{2\sinh u},\\
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}).
\end{gathered}
\tag{3}
\]
The Fourier normalization and proof are those of PTE-T2--T3.
In particular the negative scalar Omega_0/2 is retained.

Use the boundary weight and normalized quadratic Legendre
coefficient
\[
w_L(x)=-\frac14\log(1-(x/L)^2),\qquad
b_2=\left\langle g,\sqrt{\frac5{2L}}P_2(x/L)\right\rangle .
\]
FPW-T2--T3 prove
\[
\begin{gathered}
J_L[g]:=\frac18\int_I\int_I
\frac{|g(x)-g(y)|^2}{|x-y|}\,dy\,dx
\ge\frac{25}{24}N_g-\frac7{24}|b_2|^2,\\
|b_2|\le\frac{\sqrt5}{96}L^2\cosh(L/2)\sqrt{N_g}.
\end{gathered}
\tag{4}
\]

## TPW-T1. Polynomial kernel corrections have an exact finite moment residual

Let P be an even real polynomial and suppose
\[
k(u)\ge\frac1{2u}+P(u)\qquad(0<u\le2L).
\tag{5}
\]
For real smooth g supported in I define
\[
M_j(g)=\int_{\mathbb R}x^j g(x)\,dx,\qquad
\mathcal R_P[g]=\int_0^{2L}P(u)C_g(u)\,du.
\]
If P(u)=sum_(j=0)^d p_j u^(2j), then
\[
\boxed{
\mathcal R_P[g]
=\frac12\sum_{j=0}^d p_j
\sum_{r=0}^{2j}(-1)^r\binom{2j}{r}
M_{2j-r}(g)M_r(g).
}
\tag{6}
\]
The Gamma form has the lower bound
\[
\boxed{
\begin{aligned}
\mathscr A[g]\ge{}&
\left(\frac{\Omega_0}{2}+F(2L)+\frac12\log2+
\int_0^{2L}P(u)\,du\right)N_g\\
&+J_L[g]+\int_I w_L(x)|g(x)|^2\,dx-\mathcal R_P[g].
\end{aligned}
}
\tag{7}
\]

For P(u)=A-Bu^2 and even mean-zero g, the residual vanishes:
\[
\boxed{
\mathcal R_P[g]=0,\qquad
\int_0^{2L}P(u)\,du=2AL-\frac83BL^3.
}
\tag{8}
\]
Thus for g in W_L the corresponding consequence of (7) is
\[
\begin{aligned}
\mathscr A[g]\ge{}&
\left(\frac{\Omega_0}{2}+F(2L)+\frac12\log2+
2AL-\frac83BL^3+\frac{25}{24}\right)N_g\\
&+\int_I w_L(x)|g(x)|^2\,dx-\frac7{24}|b_2|^2.
\end{aligned}
\tag{9}
\]

**Proof.** The correlation is even and supported in [-2L,2L].
Changing variables in its integral gives
\[
\int_0^{2L}P(u)C_g(u)\,du
=\frac12\int_I\int_I P(x-y)g(x)g(y)\,dy\,dx.
\]
Expanding P(x-y) proves (6). In particular
\[
\int_0^{2L}C_g(u)\,du=\frac12M_0^2,\qquad
\int_0^{2L}u^2C_g(u)\,du=M_0M_2-M_1^2.
\tag{10}
\]
For even mean-zero g, M_0=M_1=0, proving the cancellation (8).

By Cauchy--Schwarz, N_g-C_g(u)>=0. Consequently (5) can
be integrated against this nonnegative function, regardless of
whether P itself has either sign. The tail u>=2L contributes
exactly F(2L)N_g. For the truncated logarithmic kernel the same
inside/outside decomposition as in BPC-T2 gives
\[
\int_0^{2L}\frac{N_g-C_g(u)}{2u}\,du
=J_L[g]+\int_I\left(\frac12\log2+w_L(x)\right)|g(x)|^2\,dx.
\tag{11}
\]
For example its boundary term is one quarter of
log(2L/(L-x))+log(2L/(L+x)), which is the displayed weight.
The singularity at zero is integrable because the difference
N_g-C_g(u) is O(u^2). No divergent terms at zero are separated.

The polynomial contribution is
N_g integral_0^(2L) P(u)du-R_P[g]. These observations prove (7).
Finally combine (4), (7) and (8) to obtain (9). QED.

The M_j in (6) are ordinary moments of the test profile.
For a quartic term p_2 u^4, evenness and mean zero leave the
residual 3p_2 M_2(g)^2. The cosh moment does not imply M_2=0.
Higher polynomial corrections therefore require their residuals
to be retained or bounded; they do not cancel automatically.

## TPW-T2. A quadratic lower bound holds through u=11/10

The elementary kernel inequality needed here is
\[
\boxed{
k(u)\ge\frac1{2u}+\frac6{25}-\frac{u^2}{25},
\qquad 0<u\le\frac{11}{10}.
}
\tag{12}
\]

**Proof.** The same series comparison as in FPW gives
\[
\frac{\sinh u}{u}
\le1+\frac{u^2}{6(1-u^2/20)}.
\tag{13}
\]
The successive ratios after the u^2/6 term are at most u^2/20,
and 1-u^2/20>0 on this range. Also
\[
e^{u/2}\ge E_4(u):=
1+\frac u2+\frac{u^2}{8}+\frac{u^3}{48}+\frac{u^4}{384}.
\]
Put
\[
B(u)=1+\frac{12u}{25}-\frac{2u^3}{25}>0 .
\]
Exact multiplication yields
\[
\begin{aligned}
&(1-u^2/20)E_4(u)
-B(u)(1-u^2/20+u^2/6)=u\,P_*(u),\\
P_*(u)={}&
\frac1{50}-\frac u{24}+\frac{119u^2}{6000}
-\frac{7u^3}{1920}+\frac{199u^4}{24000}
-\frac{u^5}{7680}.
\end{aligned}
\tag{14}
\]
Using u<=11/10 to bound the two negative higher powers gives
\[
\begin{aligned}
P_*(u)
&\ge\frac1{50}-\frac u{24}
+\frac{1519u^2}{96000}+\frac{1043u^4}{128000}\\
&\ge\frac1{50}-\frac u{24}
+\frac{79u^2}{5000}+\frac{u^4}{125}.
\end{aligned}
\]
The square (u^2-5/8)^2>=0 implies
u^4>=5u^2/4-25/64. Hence
\[
P_*(u)\ge
\frac{129u^2}{5000}-\frac u{24}+\frac{27}{1600}>0.
\tag{15}
\]
The last quadratic has positive leading coefficient and
discriminant
\[
\frac1{576}-4\frac{129}{5000}\frac{27}{1600}
=-\frac{97}{18000000}<0.
\]
Equations (13)--(15), with B(u)>0, therefore imply
\[
e^{u/2}\ge
\left(1+\frac{12u}{25}-\frac{2u^3}{25}\right)\frac{\sinh u}{u}.
\]
Divide by 2 sinh(u) to obtain (12). QED.

The linear-in-u estimate obtained by taking a constant remainder
was not needed here. The quadratic remainder has the exact
cancellation in (8), which preserves a larger scalar contribution.

## TPW-T3. Rational constants for the extended Gamma estimate

The following bounds will be used:
\[
\boxed{
\frac{6931}{10000}<\log2<\frac{1733}{2500},\qquad
\frac{\Omega_0}{2}>-\frac{2689}{1000},\qquad
F(11/10)>\frac{59}{50}.
}
\tag{16}
\]
These are elementary series inequalities; no source values or
spectra are evaluated.

**Logarithm of 2.** Integrating the geometric series gives
artanh(z)=sum_(j>=0) z^(2j+1)/(2j+1) for |z|<1.
Its positive terms at z=1/3 give
\[
\log2>
\frac23+\frac2{81}+\frac2{1215}+\frac2{15309}
=\frac{53056}{76545}>\frac{6931}{10000}.
\]
For the tail after the first three terms, replace its denominators
by 7 and sum the geometric series:
\[
\log2<
\frac{842}{1215}+\frac1{6804}
=\frac{23581}{34020}<\frac{1733}{2500}.
\]
The underlying inverse-hyperbolic integral is also recorded in
[NIST DLMF 4.37.3](https://dlmf.nist.gov/4.37.E3).

**Gamma constant.** The midpoint harmonic estimate used in FPW
gives gamma_E<H_4-log(9/2), where H_4=25/12. The same positive
series yields
\[
\log(9/2)=2\log2+\log(9/8)>
\frac{693}{500}+\frac2{17}>\frac{451}{300}.
\]
Thus gamma_E<29/50. Further, pi/2<11/7 and
\[
e^{573/500}>
\sum_{j=0}^6\frac{(573/500)^j}{j!}>
\frac{3143}{1000}>\frac{22}{7}>\pi ,
\]
so log pi<573/500. For the middle rational inequality it suffices
to bound the terms after 1+573/500 below, respectively, by
1313/2000, 1/4, 143/2000, 2/125 and 3/1000.
Their sum with 1+573/500 is 3143/1000.
Consequently
\[
\begin{aligned}
\gamma_E+\pi/2+3\log2+\log\pi
&<
\frac{29}{50}+\frac{11}{7}+\frac{5199}{2500}
+\frac{573}{500}\\
&=\frac{94098}{17500}<\frac{2689}{500}.
\end{aligned}
\]
This proves the second inequality in (16).

**Boundary tail.** The alternating exponential series gives
\[
e^{-11/20}\ge
\sum_{j=0}^5\frac{(-11/20)^j}{j!}
=\frac{221535049}{384000000}>
r:=\frac{5769}{10000}.
\tag{17}
\]
Put t=1155/32693. The logarithmic identity for artanh gives
\[
\operatorname{artanh}r=\log2-\operatorname{artanh}t.
\]
Indeed t=(3-5r)/(5-3r). Since
\[
0<t<\frac4{113},\qquad
\operatorname{artanh}t
\le\frac{t}{1-t^2}<
\frac{452}{12753}<\frac{71}{2000},
\]
we obtain
\[
\operatorname{artanh}r>
\frac{6931}{10000}-\frac{71}{2000}=\frac{411}{625}.
\tag{18}
\]
Also r<1/sqrt3<26/45. The derivative of arctan is at most one,
and the classical Archimedean lower bound gives pi>157/50.
It follows that
\[
\begin{aligned}
\arctan r
&\ge\frac{\pi}{6}-(1/\sqrt3-r)\\
&>\frac{157}{300}-\frac{26}{45}+\frac{5769}{10000}
=\frac{47021}{90000}.
\end{aligned}
\tag{19}
\]
By monotonicity in (17)--(19),
\[
F(11/10)>
\frac{411}{625}+\frac{47021}{90000}
=\frac{21241}{18000}>\frac{59}{50}.
\]
This completes the proof of (16). QED.

## TPW-T4. The stronger kernel supplies a positive weight at L=11/20

Set L_*=11/20 and define
\[
\boxed{
d_*(x)=\frac18-\frac14\log\bigl(1-(x/L_*)^2\bigr),
\qquad |x|<L_* .
}
\tag{20}
\]
For every g in W_(L_*),
\[
\boxed{
\mathscr A[g]\ge\int_{-L_*}^{L_*}d_*(x)|g(x)|^2\,dx.
}
\tag{21}
\]
This weight is even and log-convex, and has minimum 1/8.

**Proof.** Use (12) in (9), with A=6/25 and B=1/25.
At the endpoint the exact polynomial contribution is
\[
2AL_*-\frac83BL_*^3
=\frac{33}{125}-\frac{1331}{75000}
=\frac{18469}{75000}.
\tag{22}
\]
Since sqrt5<3 and cosh(11/40)<2, (4) gives
\[
|b_2|\le\frac{L_*^2}{16}\sqrt{N_g}
=\frac{121}{6400}\sqrt{N_g}.
\]
The resulting loss has coefficient
\[
\frac7{24}\left(\frac{121}{6400}\right)^2
=\frac{102487}{983040000}<\frac1{8000}.
\tag{23}
\]
Use (16), the weaker log2>693/1000, and (22)--(23) in (9).
The constant remaining before w_(L_*) is greater than
\[
\begin{aligned}
&-\frac{2689}{1000}+\frac{25}{24}+\frac{59}{50}
+\frac{693}{2000}+\frac{18469}{75000}-\frac1{8000}\\
&\hspace{2cm}=\frac{6271}{50000}-\frac1{8000}>\frac18.
\end{aligned}
\tag{24}
\]
This proves (21).

For completeness, the log-convexity is unchanged under scaling
from BPC-T4. For
D(s)=1/8-(1/4)log(1-s^2), |s|<1, one has
\[
D''(s)D(s)-D'(s)^2
\ge\frac{1-s^2+2s^4}{16(1-s^2)^2}>0.
\]
Thus log D is convex. Since d_*(x)=D(x/L_*), the claimed
properties follow. QED.

No positivity of the source was assumed to construct this weight.
The negative Gamma constant and the possible quadratic Legendre
component were both included in (24).

## TPW-T5. Prime 2 and prime 3 obey a joint row bound

At L_*=11/20 the only active prime powers are 2 and 3.
For the weight (20), let R(x) be the weighted row function from
BPC-T3:
\[
R(x)=\frac12\sum_{n=2,3}a_n
\left(
\frac{\mathbf1_I(x-u_n)}{\sqrt{d_*(x)d_*(x-u_n)}}+
\frac{\mathbf1_I(x+u_n)}{\sqrt{d_*(x)d_*(x+u_n)}}
\right).
\tag{25}
\]
Only terms with the shifted point inside I are evaluated.
Then
\[
\boxed{
\operatorname*{ess\,sup}_{x\in I}R(x)\le\frac{123}{125}<1.
}
\tag{26}
\]

**Active translations.** The positive atanh series gives
\[
\log3>1+\frac1{12}+\frac1{80}=\frac{263}{240}.
\tag{27}
\]
The positive exponential series gives
\[
e^{11/10}>
\sum_{j=0}^5\frac{(11/10)^j}{j!}
=\frac{36015101}{12000000}>3.
\]
Thus log3<11/10. Also 2log2>11/10 by (16), so 4 is
outside the prime-power cutoff. This proves that precisely
2 and 3 are active. Both u_2 and u_3 exceed L_*.
Consequently each translation supplies at most one neighbor
to a given point in I.

**Prime-2 bound on the whole window.** Since sqrt2>141/100,
(16) gives
\[
\frac{a_2}{2}
<\frac{1733}{7050}<\frac{123}{500}.
\tag{28}
\]
Even log-convexity, applied to x and -(x+u_2), gives
\[
d_*(x)d_*(x+u_2)\ge d_*(u_2/2)^2
\tag{29}
\]
whenever both points are inside I.
Moreover
\[
\frac{u_2}{2L_*}>\frac{63}{100},\qquad
1-\left(\frac{63}{100}\right)^2
=\frac{6031}{10000}<\frac{29}{48}\le e^{-1/2}.
\]
Hence w_(L_*)(u_2/2)>1/8 and d_*(u_2/2)>1/4.
The prime-2 part of (25), at every point, is therefore at most
\[
\frac{123/500}{1/4}=\frac{123}{125}.
\tag{30}
\]

**The extra energy on the prime-3 caps.** Write
\[
\mathcal C_3=\{x\in I: |x|>u_3-L_*\}.
\tag{31}
\]
Up to endpoints of measure zero, these are exactly the two
parts of the support on which a prime-3 neighbor exists.
For x in C_3, (27) gives
\[
\frac{|x|}{L_*}>
\frac{u_3}{L_*}-1>
\frac{131}{132}.
\]
Thus
\[
1-(x/L_*)^2<\frac{263}{17424}<\frac1{64},
\qquad
w_{L_*}(x)>\frac14\log64=\frac32\log2>1.
\]
Consequently
\[
\boxed{d_*(x)>\frac98\quad(x\in\mathcal C_3).}
\tag{32}
\]
Both points of a prime-3 pair belong to these caps.

At a point of C_3, the prime-2 neighbor, if present, still has
weight at least 1/8. Its contribution to (25) is therefore
bounded by
\[
\frac{123/500}{\sqrt{(9/8)(1/8)}}=\frac{82}{125}.
\tag{33}
\]
Using log3<11/10 and sqrt3>5/3 gives
a_3/2<33/100<1/3. Both weights in a prime-3 pair exceed 9/8,
so its contribution is at most
\[
\frac{1/3}{9/8}=\frac8{27}.
\tag{34}
\]
The joint bound on these caps is
\[
\frac{82}{125}+\frac8{27}
=\frac{3214}{3375}<\frac{123}{125}.
\tag{35}
\]
Outside C_3 there is no prime-3 contribution, and (30) applies.
Equations (30) and (35) prove (26). QED.

The bound is joint in position. It does not add the independent
global maxima of the two prime terms. The cap estimate (32)
simultaneously improves the prime-2 contribution at precisely
the points where prime 3 can contribute.

## TPW-T6. The complete two-prime source has a positive margin

For every 0<L<=11/20 and every g in W_L,
\[
\boxed{
\mathscr Q_\xi(g,g)\ge\frac1{500}\|g\|^2.
}
\tag{36}
\]

**Proof at L_*.** The weighted-square inequality proved in BPC-T3 is
\[
a|g(x)g(y)|
\le\frac{a}{2\sqrt{d_*(x)d_*(y)}}
\bigl(d_*(x)|g(x)|^2+d_*(y)|g(y)|^2\bigr).
\]
Integrating each translation once and changing variables gives
\[
\sum_{n=2,3}a_n|C_g(u_n)|
\le\int_I R(x)d_*(x)|g(x)|^2\,dx.
\]
Combine the actual source identity (2), the Gamma bound (21)
and the joint row bound (26):
\[
\begin{aligned}
\mathscr Q_\xi(g,g)
&\ge\frac2{125}\int_I d_*(x)|g(x)|^2\,dx\\
&\ge\frac2{125}\frac18N_g=\frac1{500}N_g .
\end{aligned}
\]
No sign is assumed for either individual correlation.

**Smaller widths.** Since both neutral moments are fixed,
W_L is contained in W_(L_*) whenever L<=L_*.
The source form and the L2 norm are unchanged by zero extension;
extra translations with disjoint support have zero correlation.
The endpoint estimate therefore proves (36) for every smaller L.
QED.

Prime 3 can contribute nontrivially on this family. Take two
sufficiently narrow even-reflected bumps centered at +/-log3/2,
inside I. Their correlation at log3 is nonzero. Correct the two
moments using two even central smooth profiles with independent
moment vectors, supported strictly inside |x|<log3-L_*.
These central corrections have no overlap at shift log3 with any
profile supported in I, so that correlation remains nonzero.
This is the same exact-moment construction used for prime 2 in FPW.
The proof above therefore handles an active term, rather than
forcing it to vanish by an additional restriction.

The larger window's margin 1/500 does not replace the stronger
BPC margin: for L<=1/2 one still has Q_xi[g]>=||g||^2/80.
Both statements concern the same restricted family.

## TPW-T7. The compact operator remains a strict contraction in the new window

Use the restricted closed energy G_L^circ from FPW/BPC, with
the same smooth neutral form core. Its exact source identity is
\[
\mathscr Q_\xi[g]
=\|(G_L^\circ)^{1/2}g\|^2-c_L\|g\|^2,\qquad
K_L^\circ=c_L(G_L^\circ)^{-1},
\]
\[
c_L=\sum_{n\in\mathcal N_L}a_n-\Omega_0/2>0.
\tag{37}
\]
For L<=11/20, the inequalities above give
\[
c_L<
\frac{123}{250}+\frac23+\frac{2689}{1000}<4.
\]
The form-core estimate (36) extends by closure and implies
G_L^circ>=(c_L+1/500)I. Therefore
\[
\boxed{
\|K_L^\circ\|
\le\frac{c_L}{c_L+1/500}\le\frac{2000}{2001}<1.
}
\tag{38}
\]
For the diagnostic coupling tau(lambda)=lambda/(1+lambda),
including tau(infinity)=1,
\[
\boxed{
I-\tau(\lambda)K_L^\circ\ge\frac1{2001}I
\quad(0\le\lambda\le\infty,\ 0<L\le11/20).
}
\tag{39}
\]
Its negative count is zero throughout this range.
BPC's stronger contraction bound for L<=1/2 is also retained.

The all-width criterion is unchanged:
\[
\boxed{
\mathrm{RH}
\Longleftrightarrow
\mathscr Q_\xi\ge0\text{ on }\mathcal W_L\text{ for every }L
\Longleftrightarrow
\|K_L^\circ\|\le1\text{ for every }L .
}
\tag{40}
\]
This equivalence is the retained FPW-T1/T5 result.
The present proof supplies (38) only in the stated two-prime window.
The restricted inverse is still formed after imposing the moments;
nothing here bounds the unrestricted K_L on every test at this width.

## Results and the remaining problem

| Reading | Retained BPC range | New TPW range |
|---|---|---|
| Fixed family | The two moments in (1) | The same two moments |
| Half-width | L<=1/2 | L<=11/20 |
| Source margin | At least 1/80 | At least 1/500 |
| Active primes at endpoint | 2 | 2 and 3 |
| Restricted compact norm | At most 256/257 | At most 2000/2001 |
| All-width RH positivity | Unproved | Unproved |

The general polynomial identity (6)--(7) gives a way to improve
the Gamma estimate without changing the test family. Its residual
must be controlled explicitly at each degree. The joint comparison
(30)--(35) gives a second tool: when another translation enters,
its location can improve the energy bound available to the existing
translations there.

Neither mechanism proves domination at every width. The current
kernel estimate is proved only through separation 11/10, and
the present joint row bound uses the specific prime-3 caps.
Further widening requires new bounds for the full active sum.
The prime power 4 enters when L>log2 and must eventually be included.

The positive result is for W_L. The earlier unrestricted estimate
L<=1/32 remains the established unrestricted range in this lineage.
No larger unrestricted range, uniform positive finite-bank size,
or new lower constant for n_-(a)/sqrt(a) follows from this note.
CWF's return to the native banks remains a fixed-profile statement.

No actual off-line zero or negative source witness is asserted.
The gap left in (40) is the positivity statement at all widths,
not a missing identification of the test family.

## Dependencies and verification scope

- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): exact Gamma energy,
  the complete prime-power source and closed-form construction.
- [FPW-1](RH_FIRST_PRIME_WINDOW.md): fixed-family RH equivalence,
  Legendre estimate and control of the quadratic component.
- [BPC-1](RH_BOUNDARY_PRIME_CAPS.md): weighted pair comparison,
  log-convexity and the stronger retained one-prime margin.
- [CWF-1](RH_CRITICAL_WAVE_FORMS.md): matching to the actual source
  and the fixed-profile return to native moment banks.
- [NIST DLMF 5.9.16](https://dlmf.nist.gov/5.9.E16):
  retained digamma identity underlying (3).
- [NIST DLMF 4.37.3](https://dlmf.nist.gov/4.37.E3):
  the elementary inverse-hyperbolic integral used for the series bounds.

The new claims are accompanied by analytic proofs, including
the kernel polynomial, moment residual, rational constants and
position-dependent two-prime estimate. Repository verification
concerns document structure, links, exact immutable-head content
and draft status. No scientific test, source computation, numerical
spectrum, independent agent or replay is requested or performed.
