# BRT-1 — A boundary-regular quadratic trial for the first-prime source

Monty · 18 October 2026.

Written continuation at 31df5fe6b71be4b134d5a15677fe910912a4d63b.
This note chooses one explicit nonzero correction to the constant
observer, rather than introducing another unspecified trial. The
correction cancels the trial's endpoint values, removes the divergent
logarithmic boundary term in its source column, and makes its actual
first-prime correlation cubic in the small overlap width.

At L=7/20 the retained prime-2 penalty on this trial is strictly less
than 3/1000000. The complete source value and the residual entering
GRR's certificate are given explicitly. A follow-up written estimate
proves q_L[g_*] > 34733/24000000 > 1/700 and improves the residual
penalty separately on degree-two and higher modes. The sign of the
complete certificate is NOT established. The further degree-four
residual bound in T8 proves that the sufficient threshold (22), using
the established lower trial margin, cannot hold for this fixed trial. T9 also proves that no nonzero finite polynomial trial can cancel the complete residual; T10 checks an endpoint-preserving quartic correction. There is no new unrestricted positive interval,
first-prime closure, RH proof, or claim that this trial is optimal.

## Actual source and domain

Use the source, closed domain D_L and operator T_L of
[OPC-1](RH_FIRST_PRIME_ONE_CHANNEL.md),
[ASM-1](RH_ARITHMETIC_SEAM_REDUCTION.md) and
[GRR-1](RH_GAMMA_REMAINDER_RETURN.md). All functions are real and
even on I=(-L,L), extended by zero. Write
\[
\begin{gathered}
N_g=\|g\|_2^2,\quad C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,\quad
\ell(g)=\int_I g(x)\cosh(x/2)\,dx,\\
k(u)=\frac{e^{-u/2}}{1-e^{-2u}},\quad
F(v)=\int_v^\infty k(u)\,du,\quad
\Omega_0=-\gamma_E-\frac\pi2-3\log2-\log\pi.
\end{gathered}
\]
Fix L=7/20 except where a formula explicitly permits general L.
Only prime power 2 is active. Put u_2=log 2 and a_2=(log 2)/sqrt(2).
Then the unchanged source is
\[
q_L[g]=\ell(g)^2+\frac{\Omega_0}{2}N_g
+\int_0^\infty k(u)(N_g-C_g(u))\,du
-a_2 C_g(u_2).
\tag{1}
\]
The constant and degree-two normalized observers are
\[
\phi_0=(2L)^{-1/2},\qquad
\phi_2=\sqrt{\frac5{2L}}P_2(x/L),\qquad
P_2(s)=\frac{3s^2-1}{2}.
\tag{2}
\]
They are orthonormal. ASM establishes that zero-extended polynomials
belong to the source operator domain. This domain statement, not an
assumption of smoothness across the endpoints, admits the trial below.

## BRT-T1. The degree-two endpoint correction is fixed exactly

Among trials g_t=phi_0-t phi_2, there is exactly one with both
endpoint values zero:
\[
\boxed{
t_*=\frac1{\sqrt5},\quad
h_*=\frac{\phi_2}{\sqrt5},\quad
g_*=\phi_0-h_*=
A_L(1-(x/L)^2),\quad A_L=\frac3{2\sqrt{2L}}.
}
\tag{3}
\]
The trial retains the unit constant reading and has norm
\[
\boxed{\langle g_*,\phi_0\rangle=1,\qquad N_{g_*}=\frac65.}
\tag{4}
\]

**Proof.** Since P_2(1)=P_2(-1)=1,
g_t(+-L)=phi_0(1-t sqrt(5)). Endpoint cancellation is therefore
equivalent to t=1/sqrt(5). Substitution gives (3).
Orthogonality gives both statements in (4):
1+t_*^2=6/5. In particular the residual test must not normalize
g_* to unit L2 norm and then forget the changed constant reading.
QED.

This coefficient is selected by boundary cancellation. It is not
asserted to minimize q_L, the residual norm, or the certificate.

## BRT-T2. Exact overlap and cubic first-prime suppression

For 0<=u<=2L, put z=u/L. The exact autocorrelation is
\[
\boxed{
C_{g_*}(u)=
\frac65-\frac32z^2+\frac34z^3-\frac3{80}z^5.
}
\tag{5}
\]
It is zero for u>=2L and even for negative u. Equivalently,
with d=(2L-u)/L in [0,2],
\[
\boxed{
C_{g_*}(u)=
\frac34d^3\left(1-\frac d2+\frac{d^2}{20}\right).
}
\tag{6}
\]
At the specified endpoint,
\[
\boxed{0<a_2 C_{g_*}(\log2)<\frac3{1000000}.}
\tag{7}
\]

**Proof.** For u>=0 the overlap is (u-L,L). With x=Ls,
\[
C_{g_*}(u)=\frac98
\int_{z-1}^{1}(1-s^2)(1-(s-z)^2)\,ds.
\]
Expanding the polynomial and integrating gives (5); substituting
z=2-d gives (6). The bracket in (6) is positive on [0,2].
It is at most one there, since -d/2+d^2/20<=0.

OPC-T2 retains
\[
\frac{6931}{10000}<\log2<\frac7{10},\qquad a_2<\frac12.
\]
Thus at L=7/20,
\[
0<d_2=\frac{7/10-\log2}{7/20}
<\frac{69}{3500}<\frac1{50}.
\]
Equation (6) now yields
\[
0<a_2 C_{g_*}(\log2)
<\frac12\frac34\left(\frac1{50}\right)^3
=\frac3{1000000}.
\]
This is the actual signed prime term's magnitude for the selected
profile, not a norm bound applying to all source profiles. QED.

For comparison, the constant profile has
C_(phi_0)(u)=1-u/(2L) on this range, which vanishes linearly at
the overlap endpoint. The cubic order in (6) results from the
trial vanishing linearly at both endpoints. No prime contribution
has been deleted or declared positive in (1).

## BRT-T3. The full trial source value has an explicit finite formula

Define the finite kernel moments
\[
M_j(L)=\int_0^{2L}u^j k(u)\,du\qquad(j=2,3,5),
\]
and
\[
p_*=\ell(g_*)=
A_L\left(\frac{16}{L}\cosh(L/2)
-\frac{32}{L^2}\sinh(L/2)\right).
\tag{8}
\]
Then at L=7/20,
\[
\boxed{
\begin{aligned}
q_L[g_*]={}&p_*^2+\frac35\Omega_0+\frac65 F(2L)\\
&+\frac3{2L^2}M_2(L)-\frac3{4L^3}M_3(L)
+\frac3{80L^5}M_5(L)-a_2 C_{g_*}(u_2).
\end{aligned}
}
\tag{9}
\]

**Proof.** Integrating (3) against cosh(x/2) gives (8).
Split the Gamma integral in (1) at 2L. On the finite interval,
(4)--(5) give
\[
N_{g_*}-C_{g_*}(u)
=\frac32(u/L)^2-\frac34(u/L)^3+\frac3{80}(u/L)^5.
\]
On the tail C=0, giving (6/5)F(2L). Substitute into (1).
Since k(u)=1/(2u)+O(1) at zero, every moment in (9) is finite.
The negative cubic coefficient is retained exactly. QED.

A source value for this one trial is an upper bound for the
returned scalar, not by itself a proof of full positivity.

## BRT-T4. The source column and complete residual are bounded

Define j(u)=k(u)-1/(2u), continuously extended by j(0)=1/4.
For x in (-L,L), the Gamma part of the source column is exactly
\[
\boxed{
\begin{aligned}
(G_\Gamma g_*)(x)={}&
\frac{g_*(x)}2\{F(L-x)+F(L+x)\}\\
&+\frac{A_L}{4}\left(1-3(x/L)^2\right)\\
&+\frac{A_L}{2L^2}\int_{-L}^{L}
 j(|x-y|)(y^2-x^2)\,dy .
\end{aligned}
}
\tag{10}
\]
The complete column is
\[
\boxed{
W(x):=(T_Lg_*)(x)
=p_*\cosh(x/2)+\frac{\Omega_0}{2}g_*(x)
+(G_\Gamma g_*)(x)
-\frac{a_2}{2}\{g_*(x-u_2)+g_*(x+u_2)\},
}
\tag{11}
\]
where the translated profiles use zero extension.

Both W and the mean-zero residual
\[
\boxed{
r_*(x)=W(x)-\frac1{2L}\int_{-L}^{L}W(y)\,dy
}
\tag{12}
\]
extend continuously to the closed interval and are bounded.

**Proof of the formula.** The symmetric Gamma operator associated
with (1) acts as
\[
G_\Gamma g(x)=\frac{g(x)}2\{F(L-x)+F(L+x)\}
+\frac12\int_I k(|x-y|)(g(x)-g(y))\,dy.
\]
For the selected polynomial,
g_*(x)-g_*(y)=A_L(y^2-x^2)/L^2. Split k into 1/(2u)+j.
The singular part simplifies exactly because
\[
\int_{-L}^{L}\frac{y^2-x^2}{|y-x|}\,dy=L^2-3x^2.
\]
This proves (10). Polarization of the pole and the real
translation correlation yields the remaining terms in (11).

**Boundary proof.** Since F(t)=-log(t)/2+O(1) as t decreases to
zero, while g_* vanishes linearly at either endpoint, the first
line of (10) tends to zero at that endpoint. The polynomial
term is bounded. The function j is continuous on [0,2L], so
its integral term is continuous on the closed interval.
The zero extension of g_* is continuous, hence the prime terms
are also continuous. This proves the claim for W and then (12).
ASM's polynomial-domain result ensures these columns are columns
of the retained source operator. QED.

The boundary cancellation removes an unbounded logarithmic column.
It does not remove all logarithmic dependence: terms of the form
distance times log(distance) remain. Boundedness is not a bound
small enough to certify the returned scalar.

## BRT-T5. The chosen trial enters the complete source certificate

Retain OPC's mean-zero projection Q_0, positive block C, column
b=Q_0 T_L phi_0, and returned scalar s_L. Because
h_*=phi_2/sqrt(5) is mean-zero and in the operator domain,
\[
b-Ch_*=Q_0T_L(\phi_0-h_*)=r_*.
\]
GRR-T5 therefore gives
\[
\boxed{
q_L[g_*]-\frac{31}{2}\int_{-L}^{L}r_*(x)^2\,dx
\le s_L\le q_L[g_*].
}
\tag{13}
\]
Here q_L[g_*] is (9), and r_* is explicitly (10)--(12).
No unknown inverse occurs in the data of this sufficient certificate.

The exact error identity behind (13) is
\[
q_L[g_*]-s_L=\langle r_*,C^{-1}r_*\rangle.
\]
All unresolved mean-zero memory is retained through this error;
the trial does not replace the actual source by a two-mode model.

**Unproved task.** A nonnegative lower endpoint in (13) would
certify full source positivity at L=7/20. No such sign is proved
here. A negative lower endpoint would only mean this certificate
has failed; it would not prove that the actual source is negative.
There is no assertion that this particular quadratic trial will pass.

## BRT-T6. The specified trial has a strictly positive actual source value

At L=7/20, without a numerical source evaluation,
\[
\boxed{
q_L[g_*]>\frac{34733}{24000000}>\frac1{700}.
}
\tag{14}
\]
This is positivity for one specified trial, not for every source
profile. In particular, (14) alone is not a lower bound for s_L.

**Proof.** Let b_old be GRR's older MCW comparison, retaining
its actual prime term. MCW-T1 gives
\[
\begin{aligned}
\mathfrak b_{\rm old}[g]={}&
B_L N_g+J_L[g]+\int_I w_L g^2+\ell(g)^2\\
&-\frac3{25}m_0(g)^2+\frac1{25}m_0(g)m_2(g)
-a_2C_g(u_2),\\
w_L(x)={}&-\frac14\log(1-(x/L)^2).
\end{aligned}
\tag{15}
\]
OPC-T3's displayed rational estimates give
B_L>-206137/300000 at this endpoint.

For the trial (3), direct integration and its two Legendre
coefficients yield
\[
\begin{gathered}
N_{g_*}=\frac65,\quad
m_0^2=2L,\quad m_0m_2=\frac{2L^3}{5},\quad J_L[g_*]=\frac3{20},\\
\int_I w_Lg_*^2=\frac{47}{100}-\frac35\log2
>\frac{676}{12500}.
\end{gathered}
\tag{16}
\]
For the last identity, substitute x=Ls and use
\[
\begin{aligned}
\int_0^1\log(1-s^2)\,ds&=2\log2-2,\\
\int_0^1s^2\log(1-s^2)\,ds&=\frac23\log2-\frac89,\\
\int_0^1s^4\log(1-s^2)\,ds&=\frac25\log2-\frac{46}{75}.
\end{aligned}
\]
These follow by integration by parts with an endpoint cutoff
and polynomial division; the logarithmic endpoint terms cancel
before the cutoff is removed. Multiply their combination with
coefficients 1,-2,1 by -9/16. The strict bound in (16) uses
OPC's retained log2<1733/2500.

The trial is nonnegative. Since cosh(x/2)>=1+x^2/8,
\[
\ell(g_*)^2\ge m_0^2+\frac14m_0m_2
=2L+\frac{L^3}{10}.
\tag{17}
\]
The additional square of m_2/8 is nonnegative and may be dropped.
Combine (7), (15)--(17) at L=7/20:
\[
\begin{aligned}
\mathfrak b_{\rm old}[g_*]
>{}&-\frac{206137}{250000}+\frac3{20}
+\frac{676}{12500}+\frac7{10}+\frac{343}{80000}\\
&-\frac{21}{250}+\frac{343}{500000}
-\frac3{1000000}
=\frac{201}{400000}.
\end{aligned}
\tag{18}
\]

GRR-T3's exact comparison difference is
\[
\mathfrak b_*[g]-\mathfrak b_{\rm old}[g]
=d_LN_g-\frac{\Delta A}{2}m_0^2+\Delta Bm_0m_2,
\]
with
\[
d_L=\frac{14203}{7200000},\quad
\Delta A=\frac{23}{4800},\quad
\Delta B=\frac{29}{2400}.
\]
On this trial that difference is exactly
\[
\frac{14203}{6000000}-\frac{161}{96000}
+\frac{9947}{48000000}
=\frac{14357}{16000000}.
\tag{19}
\]
Finally GRR-T2 gives
q_L[g]>=b_*[g]+(91/2304000)N_g.
Using N_g=6/5, (18)--(19) give
\[
q_L[g_*]>
\frac{201}{400000}+\frac{14357}{16000000}
+\frac{91}{1920000}
=\frac{34733}{24000000}.
\]
The comparison with 1/700 follows from
34733*700=24313100>24000000. QED.

## BRT-T7. Higher residual modes admit a smaller error coefficient

At L=7/20 let P_2 be projection onto phi_2 within the complete
mean-zero space, and set
\[
\delta=\frac{93383}{1440000},\qquad
\delta_4=\delta+\frac7{24}=\frac{513383}{1440000}.
\]
Then the actual positive mean-zero block obeys the form inequality
\[
\boxed{C\succeq D:=\delta P_2+\delta_4(I-P_2)>0.}
\tag{20}
\]
This does not assert that C is diagonal or commutes with P_2.

**Proof.** Retain the Legendre energy rather than replacing all
of it by its least mean-zero eigenvalue. OPC-T2 absorbs the
actual prime correlation into the boundary weight. On mean-zero
profiles the mixed ordinary-moment terms vanish. The exact
GRR comparison gain is d_L times the norm squared. Consequently,
by the same scalar bounds as GRR-T5,
\[
q_L[h]\ge
\delta\|h\|^2+
\left(J_L[h]-\frac34\|h\|^2\right)
+\int_Iv_*h^2+\ell(h)^2 .
\]
The last two terms are nonnegative. MCW-T2's retained Legendre
inequality gives
J_L[h] >= (3/4)||P_2h||^2+(25/24)||(I-P_2)h||^2.
This proves (20) on the core and then on the closed form domain.
QED.

For r_* in (12), write eta=<r_*,phi_2> and
r_perp=r_*-eta phi_2. Inverse order gives the stronger
complete-source certificate
\[
\boxed{
s_L\ge q_L[g_*]
-\frac{1440000}{93383}\eta^2
-\frac{1440000}{513383}\|r_{\perp}\|^2.
}
\tag{21}
\]
The degree-two coefficient is below 31/2; the higher-mode
coefficient is below 20/7. Thus both improve (13).

No commutativity premise is needed. For example the variational
identity
\[
\langle r,C^{-1}r\rangle
=\sup_h(2\langle r,h\rangle-\langle h,Ch\rangle)
\le\sup_h(2\langle r,h\rangle-\langle h,Dh\rangle)
=\langle r,D^{-1}r\rangle
\]
proves the inverse comparison, with the first supremum understood
on the form domain and the second over the mean-zero L2 space.
All actual mixed couplings remain in C and in the source residual.

Combining (14) and (21), a sufficient explicit condition is
\[
\boxed{
\frac{1440000}{93383}\eta^2+
\frac{1440000}{513383}\|r_\perp\|^2
\le\frac{34733}{24000000}.
}
\tag{22}
\]
T8 below proves that the left side exceeds this threshold for the
fixed quadratic trial. Thus this particular sufficient test fails.
The exact returned scalar and the sign of (21) remain undetermined.

## BRT-T8. The stated sufficient threshold fails for this quadratic trial

Let
\[
\phi_4(x)=\frac3{\sqrt{2L}}P_4(x/L),\qquad
P_4(s)=\frac{35s^4-30s^2+3}{8},
\quad L=\frac7{20}.
\]
The actual residual (12) satisfies
\[
\boxed{\langle r_*,\phi_4\rangle<-\frac3{100}.}
\tag{23}
\]
Consequently the left side of the sufficient condition (22),
denoted E_bound, obeys
\[
\boxed{
E_{\rm bound}>\frac9{5000}>
\frac{34733}{24000000}.
}
\tag{24}
\]
Thus (22) cannot hold for the fixed trial g_*.

This conclusion is about the particular sufficient test formed
from the established lower trial margin and the inverse upper
bound. It does NOT show s_L<0, and it does NOT show that the
stronger expression q_L[g_*]-E_bound in (21) is negative.
A sharper trial-value bound or inverse estimate could change
that assessment. The exact error <r_*,C^{-1}r_*> is not being
lower-bounded by E_bound.

**Step 1: bound the Gamma remainder on the actual interval.**
GRR-T1 gives
\[
j(u)=\frac14-\frac u{48}-\frac{u^2}{32}+R_0(u),
\qquad R_0(u)=\frac{E(u)}{2\sinh u}\ge0.
\]
Its positive entire-series coefficients, for m>=2, are bounded
above respectively by 1/[24(2m-1)!] and 1/[16(2m-1)!].
Indeed the remaining parts of each displayed numerator in
GRR (7) are negative. Summing gives
\[
E(u)\le\left(\frac u{24}+\frac{u^2}{16}\right)(\sinh u-u).
\]
Using sinh(u)>=u and
sinh(u)-u<=u^3 cosh(u)/6, we obtain
\[
0\le R_0(u)\le
\left(\frac{u^3}{288}+\frac{u^4}{192}\right)\cosh u.
\]
For 0<=u<=7/10,
cosh(u)<=1/(1-u^2/2)<=200/151<4/3. Therefore
\[
\boxed{
0\le R_0(u)\le\frac{14063}{4320000}<\frac1{300}.
}
\tag{25}
\]
All bounds are on a finite interval; no polynomial approximation
is integrated to infinity.

**Step 2: compute the degree-four Gamma reading.**
The Gamma operator can equivalently be written on I as
\[
G_\Gamma g=
\left(F(2L)+\frac12\log2+\int_0^{2L}j(u)\,du\right)g
+J_{\rm op}g+w_Lg
-\frac12\int_I j(|x-y|)g(y)\,dy .
\tag{26}
\]
This follows by splitting the kernel at 2L and using the retained
logarithmic decomposition; J_op is the operator of J_L.

For g_*, both g_* and J_op g_* have degree at most two.
Their projections on phi_4 vanish. The exact logarithmic reading is
\[
\langle\phi_4,w_Lg_*\rangle
=\frac92\int_0^1P_4(s)(1-s^2)
       \left[-\frac14\log(1-s^2)\right]ds
=-\frac9{280}.
\tag{27}
\]
For verification, the two integrals with P_4 and s^2P_4,
respectively, are 1/40 and 9/280; their difference is -1/140.
They follow by the same elementary logarithmic integrals as
in T6, adding
integral_0^1 s^6 log(1-s^2) ds=(2/7)log2-352/735.

The constant and u^2 parts of j in (26) produce polynomials
of degree at most two and again vanish under phi_4.
For the u part, put V(x)=integral_I |x-y|g_*(y)dy.
Since V''=2g_*,
\[
V(x)=\text{constant}+A_Lx^2-\frac{A_L}{6L^2}x^4,
\qquad
\langle\phi_4,V\rangle=-\frac{2L^2}{105}.
\]
The u coefficient in the operator is +1/96, so its contribution
is -L^2/5040. Consequently
\[
\langle\phi_4,G_\Gamma g_*\rangle
=-\frac9{280}-\frac{L^2}{5040}
-\frac12\int_I\int_I
 \phi_4(x)R_0(|x-y|)g_*(y)\,dy\,dx .
\tag{28}
\]
Here ||g_*||_1=sqrt(2L) and ||phi_4||_1<=sqrt(2L).
Equation (25) bounds the absolute double-integral contribution
including its factor 1/2 by
\[
\frac{L}{300}=\frac7{6000}.
\tag{29}
\]

**Step 3: retain the pole and the actual prime.**
Orthogonality removes the constant and quadratic terms of cosh
from <phi_4,cosh(x/2)>. Its remaining Taylor term is bounded by
x^4 cosh(L/2)/384. Also p_*<1: for example
p_*<=sqrt(2L)cosh(L/2) and cosh(7/40)<11/10.
Using L<2/5, cosh(L/2)<6/5 and ||phi_4||_1<1 yields
\[
|p_*\langle\phi_4,\cosh(x/2)\rangle|<\frac1{12500}.
\tag{30}
\]

For the prime term use |P_4(s)|<=1 on [-1,1], which follows
directly from its endpoint and stationary values. Put
d_2=(2L-u_2)/L<1/50 as in T2. Either translated profile has
integral at most A_L L d_2^2 over its overlap with I:
the exact endpoint-cap integral is A_L L(d_2^2-d_2^3/3).
Hence
\[
\left|\frac{a_2}{2}
\langle\phi_4,g_*(\cdot-u_2)+g_*(\cdot+u_2)\rangle\right|
\le\frac{9a_2}{4}d_2^2<\frac9{20000}.
\tag{31}
\]
The constant Omega term has zero phi_4 projection.
The mean subtraction in (12) also has zero projection.
Combining (28)--(31) and dropping the negative -L^2/5040 gives
\[
\langle r_*,\phi_4\rangle
<-\frac9{280}+\frac7{6000}+\frac1{12500}+\frac9{20000}
=-\frac9{280}+\frac{509}{300000}<-\frac3{100}.
\]
The last comparison follows from
509/300000<3/1400=9/280-3/100. This proves (23).

**Step 4: assess the existing sufficient condition.**
The vector phi_4 is a unit vector orthogonal to both phi_0
and phi_2. Thus ||r_perp||^2>9/10000.
Since delta_4=513383/1440000<1/2,
\[
E_{\rm bound}
=\delta^{-1}\eta^2+\delta_4^{-1}\|r_\perp\|^2
>2\frac9{10000}=\frac9{5000}.
\]
Finally 9/5000=43200/24000000>34733/24000000.
This proves (24). In fact the lower certificate formed using
only the proved trial margin is strictly below
\[
\frac{34733}{24000000}-E_{\rm bound}
<-\frac{8467}{24000000}.
\]

**Gate accounting.** No new assumption has been added. The
original residual-positivity task is still open. The specified
quadratic trial and the current conservative margin/inverse
bounds do not satisfy its stated sufficient threshold. The
degree-four residual identifies a concrete obstruction in this
attempt; it does not create a new independent RH condition.

## BRT-T9. No nonzero polynomial trial cancels the complete residual

At the fixed first-prime window L=7/20, let p be any nonzero
real even polynomial on I, extended by zero. Then
\[
\boxed{Q_0T_Lp\ne0.}
\tag{32}
\]
More strongly, this residual has infinitely many nonzero even
Legendre coefficients. This statement concerns finite polynomial
trials for the retained source, not every possible observer.

**Proof.** In a sufficiently small one-sided neighborhood of x=L,
the exact operator in T4 has the form
\[
(T_Lp)(x)=-\frac14p(x)\log(L-x)+A_p(x),
\tag{33}
\]
where A_p is real analytic across x=L.

To verify this regularity, k(t)=1/(2t)+j(t), with j real analytic
on a neighborhood of [0,2L], implies
F(t)=-(1/2)log(t)+f(t), where f is analytic at zero.
The endpoint multiplier in T4 therefore contributes exactly
the logarithm displayed in (33); F(L+x) is analytic near x=L.

In the singular interior integral, p(x)-p(y) is divisible by
x-y. Splitting the integral at y=x leaves integrals of
polynomials with polynomial endpoints, hence a polynomial in x.
The remaining interior integrals have kernels j(x-y) and j(y-x)
on the two respective subintervals. These are analytic integrands
with analytic endpoints, and extend analytically near x=L.

The pole term and Omega term are analytic. Because
0<u_2<2L, on a sufficiently small neighborhood to the left of L
the first prime translate is p(x-u_2), strictly inside its support,
while the other translate is zero. Thus their sum also has an
analytic extension near L. This proves (33). The extension here
is of the interior formula; it does not assert smoothness of
the zero-extended trial across its support boundary.

A nonzero polynomial has a finite order m>=0 of vanishing at L:
p(L-d)=c d^m+O(d^(m+1)), c nonzero. The m-th derivative of
p(L-d)log d has a nonzero logarithmic divergence as d decreases
to zero. An analytic A_p cannot cancel it. Consequently T_Lp
cannot agree with a constant, or with any polynomial, on that
one-sided neighborhood.

ASM admits p to the operator domain. If Q_0T_Lp=0 in L2,
T_Lp would equal a constant almost everywhere, and then everywhere
on the interior neighborhood by continuity, a contradiction.
If Q_0T_Lp had only finitely many Legendre coefficients, it would
instead equal a polynomial almost everywhere, giving the same
contradiction. Evenness restricts these coefficients to even
degrees. QED.

For a polynomial p with <p,phi_0>=1, the exact return identity
therefore gives the strict, source-specific statement
\[
\boxed{q_L[p]-s_L
=\langle Q_0T_Lp,C^{-1}Q_0T_Lp\rangle>0.}
\tag{34}
\]
Strict positivity follows from C>=delta I and injectivity of its
positive inverse. This does not give a uniform positive lower
bound for the error over all polynomial trials. In particular,
it does not establish a minimum probability, prohibit convergence,
or imply that a finite rigorous certificate is impossible.

## BRT-T10. A quartic correction preserves endpoints but not exact closure

The endpoint-preserving mean-zero quartic direction is
\[
\psi=\phi_4-\frac3{\sqrt5}\phi_2
=\frac{21}{8\sqrt{2L}}(1-(x/L)^2)(1-5(x/L)^2).
\tag{35}
\]
Thus p_tau=g_*+tau psi has <p_tau,phi_0>=1 and vanishes at both
endpoints for every real tau. Its norm is
\[
\|p_\tau\|^2=\frac65+\frac65\tau+\frac{14}{5}\tau^2.
\]
These identities follow directly from the orthonormal Legendre
coefficients. T9 applies to every p_tau: no choice of this
coefficient, or of finitely many polynomial coefficients, cancels
the complete residual.

For clarity about what such a correction can accomplish, set
a=q_L[psi]>0 and beta=<psi,r_*>. Exact polarization gives
\[
q_L[p_\tau]=q_L[g_*]+2\tau\beta+\tau^2a,\qquad
Q_0T_Lp_\tau=r_*+\tau C\psi.
\tag{36}
\]
The positive mean-zero bound ensures a>0. The scalar
tau=-beta/a minimizes the trial value along this line and removes
the residual reading against psi. It does not remove the complete
residual or necessarily maximize a conservative lower certificate.
Indeed,
\[
q_L[p_\tau]-
\langle r_*+\tau C\psi,C^{-1}(r_*+\tau C\psi)\rangle=s_L
\tag{37}
\]
for every tau: both quadratic changes cancel exactly.

No coefficient beta, optimized lower bound, or new positivity
interval is asserted here. The useful conclusion is narrower:
finite polynomial exact-residual closure is ruled out, while
positivity certification by a rigorously bounded nonzero residual
remains possible.

**Gate accounting.** This closes the proposed route of making
all source memory vanish by finitely many polynomial corrections.
It neither closes the residual-positivity gate nor adds a new
independent hypothesis to RH. The required task remains a signed
bound on the same exact returned scalar; it does not require
that every residual coefficient vanish.

## BRT-T11. Retaining the next Legendre level strictly improves the same certificate

Keep the trial g_* and its complete residual unchanged. Use degree
indexing phi_n=sqrt((2n+1)/(2L)) P_n(x/L) for even n, and write
r_n=<r_*,phi_n>. Let P_4 project onto phi_4 and let
R_6=I-P_2-P_4 on the mean-zero space. Put
\[
\delta_6=\delta+\frac{19}{40}
=\frac{777383}{1440000},\qquad
D_6=\delta P_2+\delta_4P_4+\delta_6R_6.
\]
Then, as closed-form inequalities,
\[
C\succeq D_6\succeq D,\qquad
s_L\ge q_L[g_*]-E_6,
\quad
E_6=\frac{r_2^2}{\delta}
+\frac{r_4^2}{\delta_4}
+\frac{\|R_6r_*\|^2}{\delta_6}.
\tag{38}
\]
In particular, T12 below proves the quantitative improvement
\[
\boxed{
E_{\rm bound}-E_6>\frac1{50000},\qquad
q_L[g_*]-E_6>
q_L[g_*]-E_{\rm bound}+\frac1{50000}.
}
\tag{39}
\]
This improves the lower certificate for the same exact returned
scalar. It is not a new positivity theorem or a change of source.

**Coercivity proof.** Before T7 discards the higher Legendre energy,
it establishes
\[
q_L[h]\ge\delta\|h\|^2+
J_L[h]-\frac34\|h\|^2+\int_Iv_*h^2+\ell(h)^2.
\]
ASM's retained energy inequality is
\[
J_L[h]\ge\frac12\sum_{n=2,4,\ldots}H_n
|\langle h,\phi_n\rangle|^2 .
\]
Since H_2/2=3/4, H_4/2=25/24 and H_6/2=49/40,
and H_n increases, dropping the last two nonnegative terms
proves C>=D_6. The argument holds on the closed form domain.
The variational inverse comparison used in T7 gives (38);
no invariance of the mode subspaces under C is assumed.

The same argument also retains the entire diagonal comparison
with entries d_n=delta+H_n/2-3/4. For every finite even N>=4,
a finite outward expression is
\[
E^{(N)}
=\sum_{\substack{2\le n\le N\\n\ {\rm even}}}\frac{r_n^2}{d_n}
+\frac{\|r_*-\sum_{2\le n\le N,\ n\ {\rm even}}r_n\phi_n\|^2}
 {d_{N+2}},
\qquad s_L\ge q_L[g_*]-E^{(N)}.
\tag{40}
\]
These expressions decrease under refinement; E^(4)=E_6.
The inequality is derived using a bounded, finitely stepped
comparison operator, so no unproved operator-domain assertion
for a new trial or an unbounded inverse product is needed.
The full diagonal comparison has form domain given by the
weighted coefficient sum; source-form domination puts every
source-form vector in that domain. No convergence rate for
unweighted source residuals is claimed.

## BRT-T12. The endpoint logarithm forces a quantitative sixth-mode gain

At the fixed window,
\[
\boxed{r_6=\langle r_*,\phi_6\rangle<-\frac1{200}.}
\tag{41}
\]

**Exact logarithmic reading.** For even n>=2, elementary
Legendre integration gives
\[
I_n:=\int_0^1P_n(s)\left[-\frac14\log(1-s^2)\right]ds
=\frac1{2n(n+1)}.
\]
For completeness, multiply the Legendre equation
((1-s^2)P_n')'+n(n+1)P_n=0 by log(1-s^2)
and integrate with an endpoint cutoff. The boundary term tends
to zero; integration by parts reduces the logarithmic integral
to -2 integral_0^1 s P_n'(s) ds / [n(n+1)].
For even n>=2, integral_0^1 P_n=0 and P_n(1)=1, giving
integral_0^1 P_n log(1-s^2)=-2/[n(n+1)].
This proves the displayed I_n.

The two applications of the recurrence for sP_n give
\[
s^2P_6=\frac{56}{195}P_8+
\left(\frac{49}{195}+\frac{36}{143}\right)P_6+
\frac{30}{143}P_4.
\]
Consequently
\[
\begin{aligned}
\int_0^1(1-s^2)P_6(s)w(s)\,ds
&=\frac1{84}-\frac{56}{195\cdot144}
-\left(\frac{49}{195}+\frac{36}{143}\right)\frac1{84}
-\frac{30}{143\cdot40}\\
&=-\frac1{756},\\
\langle\phi_6,w_Lg_*\rangle&=-\frac{\sqrt{13}}{504}.
\end{aligned}
\tag{42}
\]
Here w(s)=-log(1-s^2)/4.

**Remaining source terms.** In (26), the scalar multiple of g_*
and J_op g_* have degree at most two. The constant and u^2
parts of j give polynomials of degree at most two. The u part
gives the polynomial V of degree four displayed in T8.
All of them are orthogonal to phi_6. Hence
\[
\langle\phi_6,G_\Gamma g_*\rangle
=-\frac{\sqrt{13}}{504}
-\frac12\int_I\int_I\phi_6(x)R_0(|x-y|)g_*(y)\,dy\,dx.
\]
Because ||phi_6||_1<=sqrt(2L), exactly the same remainder bound
(25) gives absolute error at most 7/6000.

The pole error is less than 1/12500 by the proof of (30):
phi_6 is orthogonal to the constant and quadratic terms,
||phi_6||_1<=sqrt(2L)<1, and p_*<1. No sharper Taylor
remainder is needed. The Omega term and mean subtraction vanish
in this reading.

The actual prime term is bounded by the same endpoint-cap
integration as (31), now with the normalization sqrt(13):
\[
\left|\frac{a_2}{2}
\langle\phi_6,g_*(\cdot-u_2)+g_*(\cdot+u_2)\rangle\right|
\le\frac{3a_2\sqrt{13}}4d_2^2
<\frac{3\sqrt{13}}{20000}<\frac3{5000}.
\]
The bound |P_6(s)|<=1 used here follows, for example, by
expanding the polynomial identity
P_6(s)=pi^{-1} integral_0^pi
(s+i sqrt(1-s^2) cos(theta))^6 dtheta.
The integrand's base has modulus at most one on [-1,1];
the identity itself follows by the finite binomial expansion.
No numerical evaluation is involved.

Since sqrt(13)>18/5,
\[
r_6<
-\frac1{140}+\frac7{6000}+\frac1{12500}+\frac3{5000}
=-\frac1{140}+\frac{277}{150000}
<-\frac1{200}.
\]
The last inequality uses 277/150000<3/1400.
This proves (41).

**Quantitative inverse-bound reduction.** Subtract the two
explicit comparison expressions:
\[
E_{\rm bound}-E_6
=\left(\frac1{\delta_4}-\frac1{\delta_6}\right)
 \|R_6r_*\|^2 .
\]
Now delta_6-delta_4=11/60, delta_4<3/8 and delta_6<9/16.
Thus
\[
\frac1{\delta_4}-\frac1{\delta_6}
>\frac{11/60}{27/128}=\frac{352}{405}>\frac45,
\qquad
\|R_6r_*\|^2\ge r_6^2>\frac1{40000}.
\]
Their product proves (39). The endpoint logarithm has supplied
a concrete nonzero residual reading and a strict quantified
improvement; residual memory has not been set to zero.

**Assessment of the full lower certificate.** The improvement
does not rescue the old conservative trial margin. T8 still gives
\[
E_6\ge\frac{r_4^2}{\delta_4}>
\frac9{10000\delta_4}>\frac9{5000}.
\]
Therefore the explicitly assessed conservative lower bound obeys
\[
\boxed{
\frac{34733}{24000000}-E_6<
-\frac{8467}{24000000}<0.
}
\tag{43}
\]
This is a failed sufficient lower certificate, not a negative
upper bound for s_L. Neither the sign of q_L[g_*]-E_6 nor the
sign of s_L is established. Even refining all degrees above four
in (40) leaves this same obstruction when paired only with the
trial margin (14), because the degree-four denominator is unchanged.
Further work must improve the degree-four inverse treatment,
the proved trial-value margin, or the admitted trial together with
its complete residual penalty. Merely refining degrees >=6 cannot
make that particular margin pass.

**Gate accounting.** Closed here: an explicit source-specific
reduction exceeding 1/50000 in the inverse-error upper bound,
and the assessment that this refinement still cannot pass with
the current conservative trial margin. Open: the first-prime
returned scalar's sign and all-width positivity. No new independent
assumption is added. These are written derivations, not independently
validated certificates; no simulation, scientific computation,
certificate replay, independent agent or workflow dispatch is used.

## BRT-T13. Boundary surplus improves the coupled second and fourth modes

The previous refinement left the degree-four denominator unchanged.
Here retain the boundary surplus already proved for the actual prime
source, instead of discarding it in T7 and T11. No new trial is used.

For degree-indexed normalized Legendre polynomials, define the
auxiliary odd functions phi_1 and phi_3 on I by the same formula
phi_n=sqrt((2n+1)/(2L)) P_n(x/L). They are used only to read
s h, where s=x/L and h is even; they are not added to the
even source domain. For mean-zero h let h_n=<h,phi_n>.
Since v_*(x)=s^2/(4-2s^2)>=s^2/4, Bessel's inequality gives
\[
\begin{aligned}
\int_I v_*h^2
&\ge\frac14\|s h\|^2\\
&\ge\frac14\left[
 \left(\frac2{\sqrt{15}}h_2\right)^2+
 \left(\frac3{\sqrt{35}}h_2+\frac4{\sqrt{63}}h_4\right)^2
 \right]\\
&=(h_2,h_4)M(h_2,h_4)^{\mathsf T},\\
M&=\begin{pmatrix}
11/84&1/(7\sqrt5)\\
1/(7\sqrt5)&4/63
\end{pmatrix}.
\end{aligned}
\tag{44}
\]
Indeed the normalized recurrence gives
s phi_1=phi_0/sqrt(3)+2 phi_2/sqrt(15) and
s phi_3=3 phi_2/sqrt(35)+4 phi_4/sqrt(63).
The phi_0 reading vanishes because h is mean-zero.
Higher even coefficients cannot affect either of these two
odd readings. Multiplication by s is bounded on L2, so the
argument requires no regularity or additional domain assumption.

Combining (44) with the retained T11 inequality yields
\[
\boxed{
C\succeq K:=A\oplus\delta_6 I_{R_6},\qquad
A=\begin{pmatrix}
\delta+11/84&1/(7\sqrt5)\\
1/(7\sqrt5)&\delta_4+4/63
\end{pmatrix}.
}
\tag{45}
\]
The direct sum describes the comparison K only, not C.
The actual source still couples the low modes to the tail.
The boundary surplus in (44) is what remains after OPC-T2
absorbs the actual signed prime term; that term has not been
counted twice or removed from the source.

There is a simple diagonal lower bound for this retained matrix:
\[
M\succeq\operatorname{diag}(1/30,1/60).
\tag{46}
\]
To check it, the diagonal entries of the difference are
41/420 and 59/1260, and its determinant is
\[
\frac{41\cdot59}{420\cdot1260}-\frac1{245}
=\frac{259}{529200}>0.
\]
Thus A and K are strictly positive.

Write a=delta+11/84, d=delta_4+4/63 and c=1/(7 sqrt(5)).
For the same complete residual r_* as before, the new error bound is
\[
\boxed{
\begin{aligned}
E_{\rm bdry}
&=\frac{d r_2^2-2c r_2r_4+a r_4^2}{ad-c^2}
 +\frac{\|R_6r_*\|^2}{\delta_6},\\
s_L&\ge q_L[g_*]-E_{\rm bdry}.
\end{aligned}
}
\tag{47}
\]
This is an explicit inverse of the proved comparison, not a
truncation or substitution for the actual source inverse.
The variational inverse-order argument in T7 proves (47)
from (45), without any commutativity premise.

**Quantitative gain for the actual residual.** Equation (46) implies
\[
E_{\rm bdry}\le
E_{\rm diag}:=
\frac{r_2^2}{\delta+1/30}
+\frac{r_4^2}{\delta_4+1/60}
+\frac{\|R_6r_*\|^2}{\delta_6}\le E_6.
\]
In fact, the known nonzero fourth reading forces a strict gain:
\[
\begin{aligned}
E_6-E_{\rm bdry}
&\ge E_6-E_{\rm diag}\\
&\ge
\frac{1/60}{\delta_4(\delta_4+1/60)}r_4^2
>\frac19\frac9{10000}=\frac1{10000}.
\end{aligned}
\tag{48}
\]
Here delta_4<3/8 and delta_4+1/60<2/5, so the product
in the denominator is less than 3/20. T8 supplies r_4^2>9/10000.
This estimate is independent of the sign and size of r_2.
Together with T12,
\[
\boxed{
E_{\rm bound}-E_{\rm bdry}>\frac3{25000}.
}
\tag{49}
\]
Thus the lower certificate improves by this amount relative
to (21), with the source, trial, and residual all fixed.

**Full-certificate assessment, including the mixed term.**
It would be incorrect to drop the mixed term in (47) when
seeking a lower bound on E_bdry. Instead complete its square:
\[
\frac{d r_2^2-2c r_2r_4+a r_4^2}{ad-c^2}
=\frac{r_4^2}{d}
+\frac{d}{ad-c^2}
 \left(r_2-\frac c d r_4\right)^2
\ge\frac{r_4^2}{d}.
\tag{50}
\]
Since delta_4<3/8 and 4/63<1/8, d<1/2. Therefore
\[
E_{\rm bdry}>\frac9{5000},\qquad
\boxed{
\frac{34733}{24000000}-E_{\rm bdry}
<-\frac{8467}{24000000}<0.
}
\tag{51}
\]
This proves that even the coupled boundary comparison cannot
pass when paired only with the currently proved trial margin.
It does not determine the sign of q_L[g_*]-E_bdry or s_L.

The practical distinction from T12 is that degree four has now
actually been strengthened, with a controlled mixed second-mode
term, and the failure assessment survives that strengthening.
The next missing input for this route is a sharper outward
trial-value and residual enclosure, or a stronger comparison
than (45); no sign closure is claimed from the gain alone.
A trial correction must still be assessed through its entire
residual-subtracted certificate.

**Gate accounting.** Closed here: a coupled boundary-surplus
inverse estimate with actual-residual improvement exceeding
1/10000 beyond T11, and a proof that the old trial-margin test
still fails for that estimate. Open: quantitative positivity
at L=7/20 and all-width positivity. No new independent assumption.
The proof uses finite written algebra, Bessel's inequality and
inverse form order; no scientific computation, simulation,
certificate replay, subagent or workflow dispatch is performed.

## Scope and continuation

The established outputs are an explicit boundary-selected correction,
its exact overlap polynomial, the strict prime penalty bound (7),
a bounded full source column, the positive trial bound (14), and the
stronger mode-resolved residual certificate (21), and the proved
failure of the specific threshold (22) for this trial. T9--T10 add a proof that finite polynomial corrections cannot make the exact residual vanish, with an explicit endpoint-preserving quartic family. T11--T12 then use the explicit endpoint logarithm to prove a sixth-mode residual below -1/200 and reduce the same inverse-error upper bound by more than 1/50000. The refined conservative certificate still fails by (43); this is not a positivity closure.
T13 retains the boundary surplus in a coupled second/fourth-mode comparison and reduces the error upper bound by a further amount greater than 1/10000. Its full mixed-term assessment still fails with the current conservative trial margin, as (51) records.
No new positive source interval is claimed. DWE's previously proved
unrestricted interval remains unchanged.

The trial is an ordinary source profile. Its quadratic shape is not
identified with EMK connection curvature, a physical field, or a
universal cut law. Its usefulness is justified by (6) and (10).

The source-specific remaining task is a sharper signed enclosure of
(21), a sharper inverse comparison, or a better admitted trial.
T8 does not determine the sign of (21); it rules out passing (22)
with the existing lower trial margin. Boundary
regularity and a small prime term cannot replace that calculation.
Even a first-prime certificate would leave the all-width RH task open.

Only this note and its README entry change. Validation consists of
written polynomial and logarithmic integration, rational comparisons,
source polarization, boundary limits, inverse form order, and
immutable-content readback. No scientific computation, numerical
source evaluation, independent agent, certificate replay or workflow
dispatch is performed. The commit uses [skip ci]. The PR is a draft
stacked on #137 and remains unmerged.
