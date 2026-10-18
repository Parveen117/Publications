# BRT-1 — Full-source positivity through the first-prime window

Monty · 18 October 2026.

Written continuation through BRT-T21, based on branch head
06dfad2975ed96884acbac820830daba2ee827bc.
The latest result is a complete residual-subtracted lower bound
at L=7/20. The explicit quartic trial (94) satisfies
\[
q_L[p]>\frac1{5000},\qquad
\langle Q_0T_Lp,C^{-1}Q_0T_Lp\rangle<\frac1{5184}.
\]
Consequently
\[
s_{7/20}>\frac{23}{3240000}>\frac1{150000},\qquad
q_L[g]\ge\frac1{260000}\|g\|^2
\quad(0<L\le7/20,\ g\in D_L).
\]
The full-family estimate has no mean-zero or cosh-neutrality
restriction. The quartic trial has nonzero endpoint values:
its logarithmic source column and infinitely many residual
coefficients are retained and bounded, not cancelled.

T1--T18 record the preceding quadratic-trial investigation,
failed sufficient tests, observer/transport audits and partial
improvements. Their statements that the first-prime gate was
open describe those historical stages. T19--T21 now settle
that gate by another fully assessed trial; the sign of the
particular nonpolynomial test (93) is not claimed.

All-width positivity and RH remain open. No new independent
assumption is introduced. These are written proofs in the
working source chain, not independently validated certificates.

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

## Observer-angle audit — two-sided ratios and quarter-turn signs

Monty requested a check of whether the pending sign is merely a
one-chart reading artifact: multiple observer readings, two-sided
ratios, and a coordinate sign flip at a quarter turn must be retained.
This audit examines the repository's written identities at
e19719811b9cdfb9a13ff0bd3a8dcc96d36d20a7. It does not replay
certificates or independently validate the entire RH programme.

### What is already in the native/source chain

- [IRT-1](papers/native-invariant-ratio-transport/THEOREM.md),
  sections 1--2: invariance is relative to specified transport.
  A co-transformed source, numerator and denominator preserve the
  ratio. Its declared denominator is positive and nonzero; the note
  does not license evaluating the same quotient at a zero denominator.
  A calibrated log-diagonal is not an arbitrary physical angle.
- [CHD-3](papers/native-cut-generated-history/CIRCULAR_OBSERVER_COMPLETION.md),
  sections 1 and 5: Pythagoras normalizes circular probes.
  Observer, form and reference metric transform together:
  p'=Sp, H'=S^(-*)HS^(-1), so p'^*H'p'=p^*Hp.
  Probe normalization does not assume positivity of the tested H.
- [QRI-1](RH_REFLECTION_RATIO_GATE.md), T1--T3:
  the raw thermodynamic response Hessian and the transported source
  metric need not have the same signature. Axis reflection changes
  the cross coefficient b to -b while leaving b^2/(ac) unchanged.
  For an actual source bank the normalized coupling criterion is
  retained explicitly. This already addresses quadrant signs; it
  is not a new positivity result of the present audit.
- [OQD-1](RH_OBSERVER_QUOTIENT_DESCENT.md), T4:
  the two-channel cut satisfies
  G_a(Jf,Jh)=B_a(f,h) exactly. Its retained diagonal carries the
  whole original polynomial pairing, not an unsigned half of it.
- [SPT-1](RH_SHUNYA_POTENTIAL.md), T1--T2:
  the invariant potential is (r^2+y(1-y)h^2)/2. Applying the actual
  source functional preserves the invariant identity, but does not
  make that functional positive. Its products at complex source
  points remain bilinear.
- [CWF-1](RH_CRITICAL_WAVE_FORMS.md), T1--T2:
  the compact-wave form is the specified positive-normalization
  limit of native polynomial pairings and has the full prime,
  Gamma and pole formula. This is a fixed-profile limit, not a
  claimed invertible frame change between all carriers.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md), T2, and ASM:
  the same source is exactly E_L[g]-c_L||g||^2, and this identity
  is retained on the closed source domain. In particular the
  subtraction is not the sign of one chosen coordinate reading.

Within this inspected chain, no step replaces a carried pairing
by the sign of one coordinate projection, or turns a bilinear
complex-source square into an absolute square. No missing-quadrant
sign correction was identified. This is a scoped identity audit,
not a claim that every source theorem has been fully re-proved.
The EMK geometry README also explicitly separates recognition
curvature, metric curvature and global sheet memory; none of
these is silently identified with the source form here.

### The current unresolved condition is already a two-sided ratio

At L=7/20, OPC-T6 gives alpha=q_L[phi_0]>17/200, and the
mean-zero block C is strictly positive. Define
\[
\beta=\langle b,C^{-1}b\rangle\ge0,\qquad
\chi_L=\frac{\beta}{\alpha}.
\]
Then the exact existing Schur identity says
\[
\boxed{
s_L=\alpha(1-\chi_L),\qquad
q_L\ge0\ \Longleftrightarrow\ \chi_L\le1.
}
\tag{52}
\]
This is QRI's normalized coupling ratio in the current scalar/
complete-complement split. It is an equivalent interpretation,
not a new positivity theorem. No bound chi_L<=1 is supplied here.

Reversing the constant observer orientation sends b to -b,
leaves alpha and beta unchanged, and therefore leaves chi_L and
the sign of s_L unchanged. An invertible coordinate change in the
complement, with its form carried by congruence, likewise preserves
the inverse contraction beta. Equivalently beta is the supremum
of 2<b,h>-q_L[h] over the complement; a bijective reparameterization
does not change that supremum. The source form domain and selected
observer functional must be carried as well.
Changing the selected physical observer/complement is different:
the numerical ratio need not be unchanged under that operation,
although a faithful coordinate presentation preserves source inertia.

If s_L is nonzero, the reciprocal retention ratio is
alpha/s_L=1/(1-chi_L). It is undefined when s_L=0.
A reciprocal-chart singularity is not permission to replace
s_L by its absolute value. The bounded ratio chi_L remains defined
there because alpha is strictly positive.

### A complete angular calculation, including 90 degrees

When beta>0 let
\[
u=\frac{C^{-1}b}{\sqrt{\beta}}.
\]
It is an admitted mean-zero operator-domain vector, with
q_L[u]=1 and <b,u>=sqrt(beta). The two-dimensional source form
on span{phi_0,u}, normalized by its positive diagonal entries,
is exactly
\[
\begin{pmatrix}1&\sqrt{\chi_L}\\
\sqrt{\chi_L}&1\end{pmatrix}.
\]
For the circular reference preparation
\[
g_\theta=\frac{\cos\theta}{\sqrt{\alpha}}\phi_0
+\sin\theta\,u,
\]
the positive reference energy alpha c^2+q_L[h] equals one
by Pythagoras, whereas the actual source is
\[
\boxed{q_L[g_\theta]=1+\sqrt{\chi_L}\sin(2\theta).}
\tag{53}
\]
Thus both quadrants and the mixed term are included. At pi/2,
the constant-coordinate reading is zero but q_L[g_(pi/2)]=1.
After crossing pi/2 the mixed term changes sign; this does not
by itself change the sign of the complete expression.

Over all angles its least value is 1-sqrt(chi_L).
Every remaining mean-zero vector decomposes into a multiple of u
and a C-orthogonal vector z; <b,z>=0 and q_L[z]>=0.
Thus the angular test in this coupled plane, together with the
positive remainder, yields exactly (52), not a further condition.
If beta=0, positivity of C^-1 implies b=0, and the source already
splits as alpha c^2+q_L[h]>0 for nonzero states.

The angle in (53) is defined using this stated positive reference
energy. It has not been identified with an EMK sheet angle or a
TVSP compass angle. Nor does a zero coordinate projection imply
absence of all observer readings: the other coordinate survives.
An observer-frame rotation carrying the state and form together
preserves the full reading; varying theta in (53) instead selects
different preparations in the fixed source.

### Another exact ratio retaining the source subtraction

For every nonzero admitted g, PTE gives c_L>0 and
\[
\mathcal R_L(g)=\frac{\mathcal E_L[g]}{c_L\|g\|^2},\qquad
q_L[g]=c_L\|g\|^2(\mathcal R_L(g)-1).
\tag{54}
\]
Both numerator and denominator are positive, but the required
condition is R_L(g)>=1 for every g, not just positivity of the
ratio. The ratio is invariant under nonzero scalar rescaling of g
and passive coordinate changes carrying both forms. It is not
asserted constant under arbitrary changes of physical preparation.

**Audit conclusion and gate accounting.** Quarter-turn coordinate
sign changes and two-sided invariant ratios are already compatible
with the retained source. They do not remove its signed comparison
against one. The proposed one-chart sign-artifact explanation was
not supported by the inspected matching identities. No new
independent gate or assumption is introduced. The same unresolved
task is chi_(7/20)<=1, equivalently s_(7/20)>=0; all-width positivity
remains open. This audit supplies no improved numerical bound and
is not labelled a new positivity theorem.

## Co-moving observer matching — MOR return and the compact-wave limit

This continues the angle audit by explicitly pulling CWF's observers
back through MOR's transport. The observer is a function with a source
pairing and normalization, not an isolated coordinate point. The result
is a matching identity, not a new positivity theorem or a new sign bound.

Use [MOR-T2](RH_MOVING_OBSERVER_RETURN.md),
[CWF-T1](RH_CRITICAL_WAVE_FORMS.md), and the same actual source.
Here q(a)=X'(a)/X(a) is the logarithmic derivative, distinct from q_L.
Set
\[
a=b^2=\rho a_0,\quad q_0=q(a_0)>0,\quad
z_\mu=\frac{a_0}{a_0+\mu},\quad
w^0_\mu=\frac{\nu_\mu}{q_0(a_0+\mu)}.
\]
The wave polynomial P_(b,g) is precisely CWF (4), for a real even
smooth compact profile g. Define
\[
d_\rho(z)=1+(\rho-1)z,\qquad
\zeta_\rho(z)=\frac{\rho z}{d_\rho(z)}.
\]
At every actual source point,
\[
d_\rho(z_\mu)=\frac{a+\mu}{a_0+\mu},\qquad
\zeta_\rho(z_\mu)=y_a(\mu).
\tag{55}
\]
All these values are nonzero. Since Re(z_mu)>0, Re(d_rho(z_mu))>0;
choose the square root continuously from rho=1 in that half-plane.

For a polynomial P at anchor a, its pulled-back reading is
\[
(\mathcal I_\rho P)(z)
=d_\rho(z)^{-1/2}P(\zeta_\rho(z)).
\tag{56}
\]
This is the inverse MOR formula on the actual source: from
1-r y_a=1/d_rho(z_mu), MOR gives
(O_rho p)(y_a)=sqrt(d_rho(z_mu)) p(z_mu).
Thus (56) carries the observer weight as well as its argument.

**Domain distinction.** Formula (56) is not generally a polynomial,
nor holomorphic on a neighborhood of the entire closed unit disk:
its possible singular point is -1/(rho-1). We do not invoke MOR's
polynomial statement on an inadmissible input. Instead all ensuing
identities are proved directly on source values with absolute
convergence. They require no extension of the polynomial operator
domain or a guessed logarithmic trial.

Define the positive error-control norm on source-value sequences by
\[
\|F\|_{\rm abs}^2=\sum_\mu |w^0_\mu|\,|F(z_\mu)|^2.
\]
On its finite-norm sequences the original bilinear source pairing is
\[
B_0(F,H)=\sum_\mu w^0_\mu F(z_\mu)H(z_\mu),\qquad
|B_0(F,H)|\le\|F\|_{\rm abs}\|H\|_{\rm abs}.
\tag{57}
\]
Absolute values in the norm control convergence only; they do not
replace the signed/complex weights or bilinear products in B_0.
Conjugation-compatible sequences give real pairings.

### Exact finite-anchor matching with the observer carried

Normalize the returned wave by
\[
F_{b,g}(z)=\sqrt{a q_0}\,
d_\rho(z)^{-1/2}P_{b,g}(\zeta_\rho(z)).
\]
Then, before taking any limit,
\[
\boxed{
B_0(F_{b,g},F_{b,h})
=a q(a)\,\mathcal L_a[P_{b,g}P_{b,h}].
}
\tag{58}
\]
Indeed termwise cancellation gives
\[
w^0_\mu(aq_0)d_\rho(z_\mu)^{-1}
=\frac{a\nu_\mu}{a+\mu}.
\]
These are exactly the CWF-normalized source coefficients.
For each fixed b, the polynomials are bounded on the actual support
|y_a|<1 and sum nu_mu/|a+mu| is finite, proving absolute convergence.
In particular
\[
\|F_{b,g}\|_{\rm abs}^2
=\sum_\mu\frac{a\nu_\mu}{|a+\mu|}
 |P_{b,g}(y_a(\mu))|^2<\infty.
\tag{59}
\]

### The co-moving limit exists with all source weights retained

For each fixed actual source point the limit is
\[
F_g(z_\mu)=\sqrt{q_0(a_0+\mu)}\,\widehat g(v_\mu),
\qquad \mu=v_\mu^2.
\tag{60}
\]
This follows from a/d_rho(z_mu)->a_0+mu and CWF's
P_(b,g)(y_a(mu))->hat(g)(v_mu). Root branches agree by
continuation from the positive real source case, or by their
arguments in the right half-plane.

The limit is not merely pointwise. CWF's absolute tail estimate
(8), applied to g twice, gives for b>=R>=2
\[
\sum_{\gamma>R}|w^0_\mu|\,|F_{b,g}(z_\mu)|^2
\le C_g\,\frac{\log R}{R^3}.
\tag{61}
\]
For the limit,
\[
|w^0_\mu|\,|F_g(z_\mu)|^2
=\nu_\mu|\widehat g(v_\mu)|^2,
\]
whose tail tends to zero by the smooth compact profile's Fourier
decay and the retained zero-count estimate. There are finitely
many source terms below each R; on those, pointwise convergence
is norm convergence. On the tail use
|F_b-F|^2<=2|F_b|^2+2|F|^2 and (61).
First increasing R and then b proves
\[
\|F_{b,g}-F_g\|_{\rm abs}\longrightarrow0.
\tag{62}
\]
Continuity of B_0 therefore gives the fully matched identity
\[
\boxed{
B_0(F_g,F_h)
=\sum_\mu\nu_\mu\widehat g(v_\mu)\widehat h(v_\mu)
=\mathscr Q_\xi(g,h)=q_L(g,h).
}
\tag{63}
\]
The last equality is CWF/PTE's arithmetic-source identity on the
smooth core, with polarization for g,h.

Thus the compact-window source is precisely the paired reading
of the limiting, explicitly co-moving observer at the reference
anchor. This is not a fixed-degree arcsine limit with the observer
left behind. In reference coordinates F_(b,g) itself varies with b;
it is not a fixed polynomial p simply being transported outward.
A fixed p carried by MOR keeps B_0(p,h) exactly, whereas CWF's
growing wave family converges to (63). These are compatible but
different families of preparations/readings.

### What this establishes for the current window and what it does not

The identities hold for every fixed smooth profile in every
finite L. The actual closed form q_L is the retained extension
from that core. Consequently positivity on that core and on D_L
are equivalent: the forward direction follows by form closure,
and the reverse direction by restriction. A strictly negative
closed-form vector would also have negative smooth approximants.
This is sufficient to match the sign question for the existing
BRT operator-domain trials. No absolute-source evaluation formula
for every nonsmooth BRT trial is asserted without a separate
continuity proof; (60)--(62) were proved for smooth compact inputs.

The comparison preserves every cross pairing of those profiles,
not just their norms or one coordinate reading. Under an additional
faithful frame change carrying the pairing, both sides of (63)
remain the same scalar. No source sign is created by the pullback.
A positive absolute error norm in (57) does not establish positivity
of B_0, just as a Pythagorean reference circle does not establish
positivity of the tested source form.

MOR moves the anchor a while carrying source coordinates and
observer weights. [DWE](RH_DILATION_WINDOW_EXTENSION.md) instead
changes support width with U_r g(x)=r^(-1/2)g(x/r); its actual
correlations become C_g(u/r), and its kernel and pole readings
change. Exact MOR covariance cannot replace those physical-source
dilation differences by zero. In particular it cannot transport a
proved small-L margin into L=7/20 without the additional estimates
already required by DWE.

**Gate accounting.** The specified MOR/CWF co-moving matching is
now explicit, including the nonpolynomial pullback domain and a
source-weighted convergence proof using the retained CWF estimate.
No frozen-observer sign defect is present in this construction.
This resolves that matching question, not the first-prime positivity
inequality or all-width positivity. No new independent assumption,
new numerical sign bound, EMK angle identification, state-dependent
cut law or sheet-gluing premise has been introduced. Verification is
written derivation review; no scientific computation, simulation,
certificate replay, subagent or workflow dispatch is performed.

## BRT-T14. A stronger trial margin and its complete-certificate assessment

Return to the quantitative source track with the same admitted g_*.
Retaining the exact profile in the Gamma remainder and recovering
discarded scalar tail surplus gives
\[
\boxed{
q_L[g_*]>m_{\rm new}:=\frac{28967}{16000000}
>\frac9{5000}.
}
\tag{64}
\]
The complete sufficient certificate using T13 nevertheless still
fails, as proved below. This is an improved trial lower bound,
not a positive lower bound for s_L.

**Step 1: retain the trial's exact remainder integral.**
Put f(u)=N_(g_*)-C_(g_*)(u) for 0<=u<=2L. T2 gives, with z=u/L,
\[
f(u)=\frac32z^2-\frac34z^3+\frac3{80}z^5\ge0.
\]
Its three needed moments at L=7/20 are
\[
X=\int_0^{2L}f(u)\,du=\frac{49}{100},\qquad
Z=\int_0^{2L}u f(u)\,du=\frac{231}{1000},\qquad
Y=\int_0^{2L}u^2 f(u)\,du=\frac{2401}{20000}.
\tag{65}
\]
For X and Y the retained GRR moment identities give
X=2LN-m_0^2/2 and Y=(8/3)L^3N-m_0m_2, using T6's
N=6/5, m_0^2=2L, m_0m_2=2L^3/5.
For Z, direct polynomial integration gives
\[
Z=L^2\left(6-\frac{24}{5}+\frac{24}{35}\right)
=\frac{66}{35}L^2=\frac{231}{1000}.
\]

GRR's exact comparison remainder at parameter 1/2 is
Delta_*(u)=(u-1/2)^2/48+R_0(u), with R_0(u)>=0.
Therefore
\[
\begin{aligned}
q_L[g_*]-\mathfrak b_*[g_*]
&\ge\frac1{48}\int_0^{2L}(u-1/2)^2 f(u)\,du\\
&=\frac1{48}\left(Y-Z+\frac X4\right)
=\frac1{48}\frac{231}{20000}
=\frac{77}{320000}.
\end{aligned}
\tag{66}
\]
This REPLACES T6's coarser 91/1920000 contribution; it is not
added to it. The kernel remainder is integrated only on its
declared finite interval, and the signed prime term remains in
the comparison exactly as before.

**Step 2: sharpen the tail lower bound without source evaluation.**
In OPC-T1, z=tanh(7/40)<(99/100)(7/40). Its bound on
-(1/2)log z used only -log(1-1/100)>1/100.
The positive logarithm series improves that step to
\[
-\log(1-1/100)>\frac1{100}+\frac1{20000}.
\]
Keeping every other OPC-T1 estimate unchanged therefore gives
\[
F(7/10)>
\frac{149}{100}+\frac{989}{8100000}+\frac1{40000}.
\tag{67}
\]
The earlier BRT-T6 bound used only F>149/100.
Since the coefficient of F in the trial comparison is N=6/5,
(67) strengthens (18) by
\[
\frac65\left(\frac{989}{8100000}+\frac1{40000}\right)
=\frac{989}{6750000}+\frac3{100000}.
\tag{68}
\]
No other Gamma or scalar bound is changed or counted twice.

**Step 3: combine with the unchanged complete comparison.**
T6's comparison difference (19) still contributes
14357/16000000. Equations (18), (66), and (68) give
\[
q_L[g_*]>
\frac{201}{400000}+\frac{14357}{16000000}
+\frac{77}{320000}+\frac{989}{6750000}+\frac3{100000}.
\tag{69}
\]
Here 989/6750000>7/50000, since 989>945.
Replacing that term by this smaller rational bound yields
\[
q_L[g_*]>
\frac{8040+14357+3850+2240+480}{16000000}
=\frac{28967}{16000000}.
\]
Also 9/5000=28800/16000000, proving (64).
Relative to the old proved margin,
\[
m_{\rm new}-\frac{34733}{24000000}
=\frac{17435}{48000000}>0.
\tag{70}
\]
Every step is written polynomial integration or a displayed
rational/series inequality, not an independent scientific computation.

**Step 4: assess the complete mixed inverse expression.**
T13's completion of the square proves
\[
E_{\rm bdry}\ge\frac{r_4^2}{d},\qquad
d=\delta_4+\frac4{63}.
\]
Use the slightly sharper rational denominator comparison
\[
\delta_4=\frac{513383}{1440000}<\frac9{25},\qquad
d<\frac9{25}+\frac4{63}
=\frac{667}{1575}<\frac37.
\]
Together with T8's r_4^2>9/10000 this gives
\[
\boxed{E_{\rm bdry}>\frac{21}{10000}.}
\tag{71}
\]
In particular, despite the increased trial margin, the corresponding
sufficient lower certificate still obeys
\[
\boxed{
m_{\rm new}-E_{\rm bdry}
<\frac{28967-33600}{16000000}
=-\frac{4633}{16000000}<0.
}
\tag{72}
\]
This is not an assertion that q_L[g_*]-E_bdry<0, nor that s_L<0:
the actual trial value may exceed its proved lower margin.
It decisively rules out using the specific improved margin (64)
with this comparison to claim positivity.

**Decision for continuation.** The next priority is a trial or
inverse treatment that addresses the fourth-mode residual, assessed
through the full lower certificate. Merely presenting (64) as
source positivity would be wrong. For the present comparison to
pass, its trial margin must exceed at least the necessary floor
21/10000 in (71), and must also cover the other retained terms.
No success is guaranteed by crossing that floor alone.
The present fixed-trial estimate is therefore not being declared
an indefinitely sufficient route based on successive small gains.

**Gate accounting.** Improved: the actual trial lower bound.
Assessed and failed: its pairing with the current complete coupled
inverse upper bound, by the explicit amount in (72). Unresolved:
s_(7/20)>=0 and all-width positivity. No new independent assumption
or new gate is introduced. No scientific computation, simulation,
certificate replay, subagent or workflow dispatch is performed.

## BRT-T15. Retained scalar surplus strengthens the same inverse comparison

The trial and complete residual remain g_* and r_*. T14's
sharper tail bound also applies to the mean-zero comparison:
the coefficient of F(2L) there is the norm squared. In addition,
GRR-T2's uniform remainder was not consumed in constructing
T13's K. Thus, on the source mean-zero form domain,
\[
C\succeq K+\sigma I,\qquad
\sigma=\frac{989}{8100000}+\frac1{40000}
       +\frac{91}{2304000}>\frac{11}{60000}.
\tag{73}
\]
The first two terms are the scalar tail surplus in T14; the last
is GRR's uniform remainder at L=7/20. The profile-specific
remainder in (66) is used for the trial value only. It is not
inserted as a uniform mean-zero estimate.

For the displayed rational comparison,
989/8100000>61/500000 and 91/2304000>39/1000000.
Consequently sigma>93/500000>11/60000.
Define the bounded positive comparison
\[
\widehat K=K+\frac1{6000}I.
\]
Then
\[
\boxed{C\succeq\widehat K+\frac1{60000}I.}
\tag{74}
\]
This retains a strict source surplus after strengthening K.

Put epsilon=1/6000, a_+=a+epsilon, d_+=d+epsilon and
kappa=delta_6+epsilon, using T13's a,d,c. The improved complete
inverse bound is
\[
\begin{aligned}
\widehat E
&=\langle r_*,\widehat K^{-1}r_*\rangle\\
&=\frac{d_+r_2^2-2cr_2r_4+a_+r_4^2}{a_+d_+-c^2}
 +\frac{\|R_6r_*\|^2}{\kappa},\\
s_L&\ge q_L[g_*]-\widehat E.
\end{aligned}
\tag{75}
\]
For this actual residual its improvement is quantitative:
\[
\boxed{E_{\rm bdry}-\widehat E>\frac1{2000000}.}
\tag{76}
\]

Indeed c<1/15, delta<1/15, and d<3/7. The two row sums
of the positive-entry matrix A are bounded above by
111/420 and 52/105, respectively, both less than 1/2.
Since delta_6>1/2, it follows that
||K||=delta_6<27/50 and ||widehat K||<11/20.
The shifted comparison operators commute with each other;
no commutation with C is asserted. Their exact inverse
difference and T8 give
\[
\begin{aligned}
E_{\rm bdry}-\widehat E
&=\epsilon\langle r_*,K^{-1}\widehat K^{-1}r_*\rangle\\
&\ge\frac{\epsilon}{\|K\|\,\|\widehat K\|}\|r_*\|^2
>\frac{1}{6000}\frac{1000}{297}\frac9{10000}
=\frac1{1980000}>\frac1{2000000}.
\end{aligned}
\]

This scalar improvement still does not rescue the proved trial
margin. Since d<667/1575 and
1/6000<8/1575, one has d_+<3/7. Completing the same mixed
square as in (50) proves
\[
\boxed{
\widehat E\ge\frac{r_4^2}{d_+}>\frac{21}{10000},
\qquad
m_{\rm new}-\widehat E<-\frac{4633}{16000000}.
}
\tag{77}
\]
Thus (76) is not being presented as a sufficient resolution.
The bound (77) concerns the specific proved margin and comparison,
not the sign of q_L[g_*]-widehat E or of s_L.

## BRT-T16. An admitted nonpolynomial correction improves the full certificate

Here a correction is selected by the complete lower certificate,
rather than by minimizing the trial value alone. Its domain is
proved before applying the source operator a second time.

### Domain of the source-adapted direction

Every real even function v that is Holder continuous of exponent
1/2 on the closed interval, then extended by zero, belongs to
D_L and to Dom(T_L). Its source column is ASM (16), interpreted
almost everywhere.

To see form-domain membership, a bounded Holder profile has
\[
\|v-T_uv\|_2^2=O(u)\qquad(u\downarrow0):
\]
the interior Holder difference gives O(u), and the endpoint
strips created by zero extension also give O(u). This is
integrable against k(u)=O(1/u). The assertion concerns the
declared closure domain, not merely a finite formal integral:
cut v off in strips of width e at the endpoints. For the removed
part v_e, the cutoff Lipschitz bound and the Holder bound give
\[
\|v_e\|_2^2=O(e),\qquad
\|v_e-T_uv_e\|_2^2\le C\min(u,e).
\]
For u<=e, the cutoff difference contributes O(u^2/e)<=O(u);
the internal profile difference and the endpoint-crossing strips
also contribute O(u). For u>=e use the support/norm bound.
Therefore its Gamma energy is O(e(1+|log e|)), tending to zero.
Each cutoff profile can then be approximated by even smooth
compact profiles. Mollification converges in L2; its uniform
Holder bound supplies an integrable dominating translation
estimate near u=0. This proves membership in D_L.

In ASM's column formula, the interior singular integral is
bounded by a constant times
integral_0^(2L) u^(-1/2) du. The boundary term is a bounded
profile times F(L-x)+F(L+x), which belongs to L2.
The pole and the finite prime translations are bounded in L2.
Thus the formula defines a vector in H_L. Polarization against
the smooth core and form-core density prove membership in
Dom(T_L), including identification of this vector as T_Lv.

The particular W and r_* in (10)--(12) have this Holder
regularity. Near either endpoint, write
F(t)=-(1/2)log t plus an analytic function. Since g_* vanishes
linearly, the only nonsmooth boundary products have the form
t log t times a smooth function, which is Holder of every
exponent below one. Also j is analytic through its removable
value at zero, so the remaining integral term in (10) is
Lipschitz. The zero-extended translated g_* terms are Lipschitz.
Mean subtraction preserves the Holder property.

Define the concrete mean-zero correction direction
\[
\boxed{
z=\widehat K^{-1}r_*
=z_2\phi_2+z_4\phi_4+\kappa^{-1}R_6r_*,
\quad
\binom{z_2}{z_4}
=\begin{pmatrix}a_+&c\\c&d_+\end{pmatrix}^{-1}
 \binom{r_2}{r_4}.
}
\tag{78}
\]
This is a scalar multiple of r_* plus a finite polynomial.
The preceding argument therefore proves
\[
z\in Q_0H_L\cap\operatorname{Dom}(T_L)
\subset\operatorname{Dom}(C),\qquad Cz=Q_0T_Lz.
\]
T9's infinitely many residual coefficients also show that z is
nonpolynomial. There is no assumed membership of an arbitrary
logarithmic ansatz. Equation (78) uses only the explicit residual
and the inverse of the proved bounded comparison, not C^{-1}.

### Optimize the complete lower certificate on this direction

Set
\[
\begin{gathered}
E=\langle r_*,z\rangle=\widehat E,\qquad
v=Cz-r_*,\\
B=q_L[z]-E=\langle z,v\rangle,\qquad
N=\langle v,\widehat K^{-1}v\rangle .
\end{gathered}
\tag{79}
\]
These quantities are finite by the domain proof. They are
specified by the actual source integrals; the signed prime,
the pole and all source couplings remain in Cz. They require
no inverse of the unknown C.

The retained strict surplus (74) gives the source-specific bound
\[
\boxed{
B\ge\frac1{60000}\|z\|^2>0.
}
\tag{80}
\]
Here r_4 is nonzero by T8, so r_* and z are nonzero. It follows
that v is nonzero, N>0, and
\[
\boxed{t_{\rm opt}=\frac{B}{B+N}\in(0,1),\qquad p_{\rm opt}=g_*-t_{\rm opt}z}
\tag{81}
\]
is a specified admitted nonpolynomial trial with constant
reading one. 
For arbitrary real t write p_t=g_*-tz and r_t=r_*-tCz.
Define its complete lower certificate by
\[
\mathcal L(t)=q_L[p_t]
-\langle r_t,\widehat K^{-1}r_t\rangle.
\]
Then
\[
\boxed{
\mathcal L(t)=q_L[g_*]-E+2tB-t^2(B+N).
}
\tag{82}
\]
In particular t_{\rm opt} is the unique maximizing step on this line and
\[
\boxed{
s_L\ge\mathcal L(t_{\rm opt})
=q_L[g_*]-\widehat E+\frac{B^2}{B+N}.
}
\tag{83}
\]
This is a strict improvement of the full lower certificate for
the same returned scalar, with its exact gain displayed.
In combination with (76), it improves the T13 certificate by
more than 1/2000000, plus the strictly positive gain in (83).

For a direct check, q_L[z]=E+B and Cz=r_*+v, so
\[
\begin{aligned}
q_L[p_t]&=q_L[g_*]-2tE+t^2(E+B),\\
\langle r_t,\widehat K^{-1}r_t\rangle
&=E-2t(E+B)+t^2(E+2B+N).
\end{aligned}
\]
Subtracting proves (82). This accounts for the entire residual,
including its mixed terms. It does not infer improvement merely
from a smaller q_L[p_t].

There is also a precise comparison with trial-value minimization.
Weighted Cauchy--Schwarz gives B^2<=EN. Therefore
\[
\boxed{
t_{\rm opt}=\frac{B}{B+N}\le\frac{E}{E+B}=t_{\rm trial},
}
\tag{84}
\]
where t_trial minimizes q_L[p_t]. Equality requires v to be a
scalar multiple of r_*; it is not assumed. This explains why
minimizing the trial value is generally a different operation
from maximizing a conservative lower certificate.

### Quantitative result and unresolved sign

Using the proved margin from T14, the new rigorous lower bound is
\[
\boxed{
s_L>m_{\rm new}-\widehat E+\frac{B^2}{B+N}.
}
\tag{85}
\]
The correction gain is proved strictly positive and is defined
by finite source quantities on a proved domain. T17 below now
settles the sign of the right-hand side of (85): the gain cannot rescue m_new, and that
specific lower certificate is strictly negative. It is retained
here as the assessed old-margin test, not as a still-unevaluated
condition. T18 supplies a stronger source-value input and states
the current unresolved certificate in (93).

This is a constructed certificate-improving step, not a new
positivity theorem or a claimed successful window. No assertion
is made that its complete residual vanishes. No residual
component is discarded or required to vanish.

**Gate accounting.** Established: the nonpolynomial direction's
source-operator domain, a quantitative inverse-bound improvement
greater than 1/2000000, and an explicit step maximizing and
strictly improving the complete lower certificate on that
direction. Assessed and failed: the strengthened inverse
comparison paired with m_new before the correction, by (77).
The later T17 assessment rules out sufficient gain for (85);
the stronger-input test in T18, the first-prime scalar's sign,
and all-width positivity remain unresolved. This is the same
original sign obligation, not a new independent RH assumption. No scientific computation,
numerical simulation, certificate replay, subagent or workflow
dispatch is performed.

## BRT-T17. Actual inverse memory rules out rescuing the old trial margin

T16 improved the complete certificate, but positivity requires
a sufficiently large gain. Here the actual inverse memory,
rather than just an inverse-comparison expression, has a
proved lower bound:
\[
\boxed{
\langle r_*,C^{-1}r_*\rangle>\frac1{500}.
}
\tag{86}
\]
Consequently the T16 gain cannot rescue T14's margin
m_new=28967/16000000. This is a restriction on that certificate's
input margin, not a negative-source result.

### A scalar enclosure with the cancelling terms kept together

Define, for T17 and the subsequent source estimates,
\[
\mathfrak d=-\frac{\Omega_0}{2}-F(7/10),\qquad
x=\frac7{40},\quad z_0=\tanh x.
\]
We will use the written rational enclosure
\[
\boxed{
\frac{11953}{10000}<\mathfrak d
<\frac{1435169}{1200000}<\frac{299}{250}.
}
\tag{87}
\]
In particular, no separate rounding of the cancelling pi/4
terms is needed.

Here are the scalar details. [FCP's](RH_FINITE_CUT_POSITIVITY.md) integral identity for the
Euler constant also supplies a lower bound:
\[
H_n-\log n-\frac1{2n}+\frac1{12n^2}-\frac1{120n^4}
<\gamma_E
<H_n-\log n-\frac1{2n}+\frac1{12n^2}.
\]
For the additional lower estimate, the positive entire series
\[
x\cosh x-(1+x^2/3-x^4/45)\sinh x
=\sum_{m\ge3}
\frac{16m(m-1)(m-2)(m+2)}{45(2m+1)!}x^{2m+1}>0
\]
proves coth x>1/x+x/3-x^3/45. Substitution of x=t/2 gives
\[
\frac1t-\frac1{e^t-1}<\frac12-\frac t{12}+\frac{t^3}{720}.
\]
Integrate this inequality in FCP's retained Laplace identity.
At n=4 the resulting two-sided estimate is
\[
\frac{60319}{30720}-2\log2<\gamma_E
<\frac{377}{192}-2\log2.
\]
Using F(7/10)=-(1/2)log z_0+pi/4-arctan z_0, the pi/4
terms cancel exactly and hence
\[
\frac{60319}{61440}+\frac12\log(2\pi z_0)+\arctan z_0
<\mathfrak d
<\frac{377}{384}+\frac12\log(2\pi z_0)+\arctan z_0.
\]

For the lower side, tanh x>x-x^3/3 gives
z_0>33257/192000. This elementary tanh inequality follows by
integrating (tanh x)'=1-tanh^2 x>1-x^2.
PRW's pi/4>7853/10000 then yields
\[
2\pi z_0>
\frac{261167221}{240000000}>\frac{136}{125}.
\]
The positive atanh series and the alternating arctan series give
\[
\frac12\log(2\pi z_0)>\frac{11}{261},\qquad
\arctan z_0>z_0-\frac{z_0^3}{3}
>x-\frac{2x^3}{3}=\frac{16457}{96000}.
\]
Therefore
\[
\mathfrak d>
\frac{60319}{61440}+\frac{11}{261}+\frac{16457}{96000}
=\frac{239}{200}+\frac{43729}{133632000}
>\frac{11953}{10000}.
\]

For the upper side, OPC gives z_0<693/4000 and pi<22/7,
so 2 pi z_0<1089/1000. For t=89/1000,
\[
\log(1+t)<t-t^2/2+t^3/3<\frac{853}{10000}.
\]
Also z_0>33257/192000>433/2500 gives
z_0^3/3>173/100000, while z_0^5/5<x^5/5<1/30000.
Thus
\[
\arctan z_0<
\frac{693}{4000}-\frac{173}{100000}+\frac1{30000}.
\]
Combining the upper estimates gives
\[
\mathfrak d<
\frac{377}{384}+\frac{853}{20000}
+\frac{693}{4000}-\frac{173}{100000}+\frac1{30000}
=\frac{1435169}{1200000}<\frac{299}{250}.
\]
This proves (87) by written rational and series comparisons.

### A positive second-mode residual reading

The same actual quadratic trial satisfies
\[
\boxed{r_2=\langle r_*,\phi_2\rangle>\frac1{125}.}
\tag{88}
\]
To control the Gamma remainder, integrate T8's pointwise
R_0 bound on the actual finite interval:
\[
J_R:=\int_0^{7/10}R_0(u)\,du
\le\frac{(7/10)^4}{864}+\frac{(7/10)^5}{720}
=\frac{55223}{108000000}<\frac1{1900}.
\]
Its positive bounded form is
integral_0^(2L) R_0(u)(N_g-C_g(u)) du, and its operator
mathcal R obeys 0<=mathcal R<=2 J_R I. Therefore
\[
|\langle\phi_2,\mathcal Rg_*\rangle|
\le2J_R\sqrt{6/5}<\frac{11}{9500}.
\]

For j_0(u)=1/4-u/48-u^2/32,
\[
J_0:=\int_0^{7/10}j_0(u)\,du=\frac{15967}{96000}.
\]
Projecting (26) and the scalar Omega term on phi_2 gives,
before the pole, prime and R_0 terms,
\[
\frac1{\sqrt5}\left(
\mathfrak d-J_0-\frac{61}{60}
+\frac{L^2}{112}+\frac{L^3}{48}\right).
\]
For verification, the logarithmic multiplier contributes
\[
\langle\phi_2,w_Lg_*\rangle
=\frac1{\sqrt5}\left(\frac{\log2}{2}-\frac4{15}\right);
\]
its log2 term cancels the scalar log2 term in (26).
The |x-y| convolution contributes L^2/(112 sqrt(5)),
and the (x-y)^2 convolution contributes L^3/(48 sqrt(5)).
These follow respectively from
<phi_2,V>=6L^2/(7 sqrt(5)) for T8's polynomial V and
<phi_2,x^2>=sqrt(2L) 2L^2/(3 sqrt(5)).

The pole reading is bounded below by L^3/(6 sqrt(5)).
Indeed ell(g_*)>=sqrt(2L), and
\[
\langle\phi_2,\cosh(x/2)\rangle
\ge\frac{\sqrt{2L}L^2}{12\sqrt5}.
\]
Every higher even Taylor term has nonnegative phi_2 integral,
since integral_-1^1 P_2(s)s^(2m) ds
=4m/((2m+1)(2m+3)) for m>=1.

The endpoint-cap argument in T8, now with normalization sqrt(5),
bounds the absolute prime contribution by
3 sqrt(5)/20000<27/80000. Consequently
\[
r_2>
\frac1{\sqrt5}\left(\mathfrak d-\frac{450761}{384000}\right)
-\frac{11}{9500}-\frac{27}{80000}.
\]
By (87), the parenthesis exceeds 171/8000; since
sqrt(5)<9/4, its quotient exceeds 19/2000.
Also
\[
\frac{11}{9500}+\frac{27}{80000}
=\frac{2273}{1520000}<\frac3{2000}.
\]
This proves r_2>16/2000=1/125.

### An upper energy bound for the existing quartic direction

Use T10's admitted mean-zero polynomial
\[
\psi=\phi_4-\frac3{\sqrt5}\phi_2,\qquad
A_\psi=q_L[\psi].
\]
Then
\[
\boxed{0<A_\psi<\frac45.}
\tag{89}
\]
Positivity follows from the already proved mean-zero coercivity.
For the upper bound, direct polynomial and logarithmic integration
gives
\[
\|\psi\|^2=\frac{14}{5},\qquad J_L[\psi]=\frac{287}{120},
\qquad
\int_I w_L\psi^2=\frac{2299}{1800}-\frac75\log2.
\]
For the last identity, use
psi=(21/(8 sqrt(2L)))(1-6s^2+5s^4), s=x/L.
The square has coefficients 1,-12,46,-60,25 in degrees
0,2,4,6,8. The logarithmic moment in degree eight needed in
addition to T6/T8 is
\[
\int_0^1s^8w(s)\,ds=\frac{563}{5670}-\frac{\log2}{18}.
\]
It follows by the same polynomial-division integration;
substitution of the five moments gives the displayed value.

The constant and u^2 convolution pieces of j_0 vanish on this
even mean-zero profile. The u piece is nonpositive:
\[
\int_I\int_I |x-y|\psi(x)\psi(y)\,dx\,dy
=-2\int_I\left(\int_{-L}^x\psi(y)\,dy\right)^2dx\le0.
\]
Combining the scalar term with the logarithmic multiplier
cancels log2 and bounds all non-remainder Gamma/scalar terms by
\[
\frac{14}{5}\left(-\mathfrak d+\frac{15967}{96000}\right)
+\frac{287}{120}+\frac{2299}{1800}
<
\frac{567787}{720000}<\frac{79}{100},
\]
where the weaker consequence mathfrak d>239/200 of (87)
already suffices.

The R_0 form contributes at most
2J_R ||psi||^2<14/4750<3/1000.
Since psi has mean zero and
cosh(L/2)-1<=49/3151<1/60,
\[
\ell(\psi)^2
\le\|\psi\|^2(2L)(1/60)^2
=\frac{49}{90000}<\frac1{1000}.
\]
Finally |psi(x)|<=7g_*(x) on I. T2 therefore gives the actual
signed prime term the upper bound
\[
-a_2 C_\psi(u_2)\le a_2|C_\psi(u_2)|
\le49a_2C_{g_*}(u_2)<\frac{147}{1000000}<\frac1{1000}.
\]
Altogether
A_psi<79/100+3/1000+1/1000+1/1000=159/200<4/5.
This proves (89) with the pole, prime and remainder retained.

### Actual inverse-memory lower bound and the optimized-gain ceiling

Let beta=<psi,r_*>=r_4-(3/sqrt(5))r_2. T8, (88), and
3/sqrt(5)>4/3 imply
\[
\beta<-\frac3{100}-\frac4{375}
=-\frac{61}{1500}<-\frac1{25}.
\]
Cauchy--Schwarz in the positive C form now gives
\[
\langle r_*,C^{-1}r_*\rangle
\ge\frac{\beta^2}{A_\psi}>
\frac{(1/25)^2}{4/5}=\frac1{500}.
\]
This proves (86). It is a lower bound for the actual inverse
memory; it is not obtained by treating an inverse upper
comparison as a lower estimate.

Write G=B^2/(B+N) for T16's correction gain. Since (83) is
a lower bound for the same exact scalar,
\[
q_L[g_*]-\widehat E+G
\le s_L=q_L[g_*]-\langle r_*,C^{-1}r_*\rangle .
\]
It follows that
\[
\boxed{
G<\widehat E-\frac1{500},\qquad
m_{\rm new}-\widehat E+G
<-\frac{3033}{16000000}<0.
}
\tag{90}
\]
Thus even the optimized nonpolynomial step cannot pass the
certificate formed by propagating the particular old margin
m_new. There is no need to enclose its nested source integrals
more tightly to decide that specific question.

More generally, changing the trial and propagating this same
base margin through exact quadratic differences cannot overcome
this obstruction, even with the exact residual inverse:
the best such lower value is
m_new-<r_*,C^{-1}r_*>. Directly deriving a stronger lower
trial value is a different input and is not excluded.

**Scope.** Neither s_L<0 nor q_L[g_*]-widehat E+G<0 follows.
The actual trial value was replaced by a strict lower margin
in (90). T18 below improves precisely that input.

## BRT-T18. Keeping the scalar cancellation yields a substantially stronger trial margin

The scalar enclosure (87) also improves the actual trial value:
\[
\boxed{
q_L[g_*]>
m_{\rm sharp}:=\frac{48823}{16000000}>\frac3{1000}.
}
\tag{91}
\]
Relative to T14 this is a gain of
\[
m_{\rm sharp}-m_{\rm new}=\frac{1241}{1000000}.
\tag{92}
\]
This is a stronger source-value lower bound for the same trial;
it is not itself a positive lower bound for s_L.

For completeness, derive it directly from the exact source.
Set f=N_(g_*)-C_(g_*) on (0,2L). T14's exact moments are
X=49/100, Z=231/1000 and Y=2401/20000.
The singular-kernel part has the exact value
\[
J_L[g_*]+\frac12\log2\,\|g_*\|^2+\int_Iw_Lg_*^2
=\frac3{20}+\frac{47}{100}=\frac{31}{50}.
\]
The two log2 terms cancel before any bound is taken. GRR's
exact expansion k=1/(2u)+1/4-u/48-u^2/32+R_0 and R_0>=0 give
\[
\begin{aligned}
q_L[g_*]
&\ge\ell(g_*)^2-\frac65\mathfrak d
+\frac{31}{50}+\frac X4-\frac Z{48}-\frac Y{32}
-a_2 C_{g_*}(u_2)\\
&>\frac{920463}{640000}-\frac65\mathfrak d
-\frac3{1000000}.
\end{aligned}
\]
The last line uses the same pole lower bound
ell(g_*)^2>=2L+L^3/10 and prime bound as T6.
The finite polynomial terms are used only on (0,2L);
the tail remains in mathfrak d. Substitution of (87) gives
\[
q_L[g_*]>
\frac{920463}{640000}
-\frac{1435169}{1000000}-\frac3{1000000}
=\frac{48823}{16000000}.
\]
Since 48823>48000, this proves (91).
The R_0 remainder was dropped as nonnegative in this derivation;
no earlier remainder contribution is added a second time.

The currently strengthened corrected certificate is therefore
\[
\boxed{
s_L>
\frac{48823}{16000000}-\widehat E+\frac{B^2}{B+N}.
}
\tag{93}
\]
Its sign remains undetermined. T17's obstruction applies to
m_new, not automatically to m_sharp. Conversely, the fact that
m_sharp exceeds 1/500 does not bound the actual inverse memory
from above, and does not establish (93)'s positivity.

**Continuation decision.** The old-margin correction test has
been resolved negatively by (90), so it should not be retried.
The source-value input has now been strengthened by (91).
The remaining quantitative work is to bound the complete
residual penalty after correction against this stronger input,
retaining both low-mode mixing and the logarithmic tail.
This remains the same first-prime sign problem, not a new
independent condition or an all-width conclusion.

**Gate accounting.** Established: a lower bound greater than
1/500 for the actual inverse memory; a rigorous ceiling on the
T16 gain that rules out rescuing the old margin; and the
stronger trial margin (91). Unresolved: the sign of (93),
first-prime positivity, and all-width positivity. No new
independent assumption is added. All estimates are written
source identities, polynomial/logarithmic integrations and
rational series comparisons. No scientific computation,
simulation, certificate replay, subagent or workflow dispatch
is performed.

## BRT-T19. An explicit quartic trial has a usable positive source margin

Fix L=7/20 and, in this continuation, set
\[
\boxed{
p(x)=\phi_0 f(x/L),\qquad
f(s)=1-\frac87P_2(s)+\frac15P_4(s)
=\frac{461}{280}-\frac{69}{28}s^2+\frac78s^4.
}
\tag{94}
\]
This is a particular polynomial, not an unevaluated quartic family.
Its coefficients are assessed below through the complete residual.
ASM admits its zero extension to the source operator domain.
It has constant reading one and
\[
p=\phi_0-\frac8{7\sqrt5}\phi_2+\frac1{15}\phi_4,\qquad
\|p\|^2=\frac{13954}{11025},\qquad
\|p-\phi_0\|^2=\frac{2929}{11025}.
\]
The polynomial f decreases as a function of s^2 on [0,1]:
its derivative there is at most -5/7. Hence
f(s)>=f(1)=2/35>0. Its endpoint values are nonzero.
The logarithmic source column is therefore retained, not
asserted bounded or removed.

The actual trial value satisfies
\[
\boxed{
q_L[p]>\frac{1013}{5000000}>\frac1{5000}.
}
\tag{95}
\]

### Pole, finite-kernel and boundary terms

Continue to use mathfrak d from (87) and
j_0(u)=1/4-u/48-u^2/32, with
J_0=integral_0^(2L) j_0=15967/96000.
For this unit-mean trial,
\[
m_0=\sqrt{2L},\qquad
m_2=\sqrt{2L}\,\frac{19L^2}{105}.
\]
Nonnegativity of p implies
\[
\ell(p)^2\ge(m_0+m_2/8)^2
=2L+\frac{19L^3}{210}+\frac{m_2^2}{64}.
\]
The last contribution obeys
\[
\frac{m_2^2}{64}
=\frac{6067327}{1128960000000}>\frac1{200000}.
\]

The singular-kernel part, including its boundary term, is exactly
\[
\begin{aligned}
\mathcal J_p
&:=J_L[p]+\frac12\log2\,\|p\|^2+\int_Iw_Lp^2\\
&=\frac12+\frac{1376}{3675}+\frac{1091}{141750}
-\frac4{21}+\frac1{100}-\frac4{245}
=\frac{194287}{283500}.
\end{aligned}
\tag{96}
\]
For a direct check, the diagonal entries of the multiplier
w_L+(log2)/2 on phi_0,phi_2,phi_4 are
1/2,41/60,1739/2520. Its off-diagonal entries are
sqrt(5)/12,3/40,3 sqrt(5)/28. Combine these with the
Legendre energies 0,3/4,25/24 and the coefficients in (94).
These entries follow from the same polynomial/logarithmic
integrals used above; T20 also gives the general off-diagonal
identity. No logarithm is rounded before its cancellation.

The constant convolution part of j_0 contributes -L/4.
Its quadratic convolution part contributes m_0m_2/32
=19L^3/1680. For the linear part put
\[
H(s)=\int_0^s f(t)\,dt
=\frac{461s-230s^3+49s^5}{280}.
\]
Direct integration gives
\[
\int_0^1H^2
=\frac1{78400}\left(
\frac{212521}{3}-\frac{212060}{5}
+\frac{98078}{7}-\frac{22540}{9}+\frac{2401}{11}
\right)<\frac{513}{1000}.
\]
For example the numerator in parentheses is less than
70841-42412+14012-2504+219=40156, and
40156/78400<513/1000.
The distance-kernel identity, using H(1)=1, is
\[
\int_0^{2L}u C_p(u)\,du
=L^2\left(1-\int_0^1H^2\right)>\frac{487L^2}{1000}.
\]
It follows either by integrating twice the |x-y| kernel,
or from its cumulative-integral formula. The linear convolution
contribution is therefore greater than
487L^2/48000=23863/19200000.

### Actual prime-2 term

Write d=(2L-u_2)/L<1/50 and F_c(t)=f(1-t) for 0<=t<=d.
Here F_c is a cap polynomial, distinct from the Gamma tail F.
Its exact expansion is
\[
F_c(t)=\frac2{35}+\frac{10}{7}t+\frac{39}{14}t^2
-\frac72t^3+\frac78t^4.
\]
On this cap,
\[
0<F_c(t)\le A+B t,\qquad A=\frac2{35},\quad B=\frac{11}{7}.
\]
Indeed the coefficient after factoring out t is at most
10/7+39/700+7/1000000<11/7.
Since the two endpoints are paired by the prime displacement,
\[
C_p(u_2)=\frac12\int_0^d F_c(t)F_c(d-t)\,dt.
\]
Using a_2<1/2 and the affine cap bound gives
\[
\begin{aligned}
0<a_2 C_p(u_2)
&<\frac14\left(\frac{A^2}{50}+\frac{AB}{2500}
+\frac{B^2}{750000}\right)\\
&=\frac{3841}{147000000}<\frac{27}{1000000}.
\end{aligned}
\tag{97}
\]
This bounds the actual negative prime term; it is not omitted.

### Complete lower source bound

The R_0 kernel remainder contributes a nonnegative quadratic
form, so it can be dropped for this lower estimate.
Combining the preceding identities gives
\[
\begin{aligned}
q_L[p]>{}&
(-\mathfrak d+J_0)\frac{13954}{11025}
+\frac{194287}{283500}
+\frac{168931}{240000}-\frac7{80}\\
&+\frac{6517}{13440000}
+\frac{23863}{19200000}
+\frac1{200000}-\frac{27}{1000000}.
\end{aligned}
\tag{98}
\]
All polynomial-kernel contributions are restricted to (0,2L);
the Gamma tail and scalar remain in mathfrak d.

For an outward rational check, (87) yields
\[
(-\mathfrak d+J_0)\frac{13954}{11025}
>-\frac{11494202834}{8820000000}
>-\frac{13031977}{10000000}.
\]
The other positive fractions in (98), in their displayed order,
are greater than
6853156/10000000, 7038791/10000000,
4848/10000000, and 12428/10000000.
The exact negative constant is -875000/10000000;
the last two terms are 50/10000000 and -270/10000000.
Thus (98) is strictly greater than
\[
\frac{-13031977+6853156+7038791-875000
      +4848+12428+50-270}{10000000}
=\frac{2026}{10000000}.
\]
This proves (95) entirely by the written source comparison.

## BRT-T20. The quartic trial's complete inverse residual is smaller than its margin

Let r=Q_0 T_Lp be the complete actual residual of (94).
We prove
\[
\boxed{
\langle r,C^{-1}r\rangle
\le\|r\|_{\mathcal D^{-1}}^2
<\left(\frac{829}{60000}\right)^2
<\frac1{5184}.
}
\tag{99}
\]
The bounded comparison here is
\[
\mathcal D=\frac3{20}P_2+\frac3{10}P_4+\frac12R_6,
\quad
\|v\|_{\mathcal D^{-1}}^2
=\frac{20}{3}v_2^2+\frac{10}{3}v_4^2+2\|R_6v\|^2.
\tag{100}
\]
This comparison is a lower bound for C, not a diagonalization
of the source. To check it, T13's low matrix A has
a>19/100, d>41/100 and c^2=1/245. The diagonal entries of
A-diag(3/20,3/10) exceed 1/25 and 11/100, and its
determinant exceeds 11/2500-1/245>0. Also delta_6>1/2.
Thus C>=K>=mathcal D, and the variational inverse-order
argument gives the first inequality of (99).

We bound four complete residual components in the norm (100).
There is no truncation of the logarithmic tail.

### Main column and its exact low-mode coefficients

Define r^(0) by using j_0 in place of j in (26), retaining
the actual scalar mathfrak d and the full logarithmic multiplier,
and using m_0(1+x^2/8) for the pole column. Omit only R_0 and
the prime in this intermediate definition. Their full columns,
and the pole difference, are restored below.

The exact second-mode coefficient is
\[
\begin{aligned}
\sqrt5\,r^{(0)}_2={}&
\frac{19}{42}+\frac{343}{48000}+\frac{343}{384000}
+\frac{491}{432000}
-\frac87\left(-\mathfrak d+\frac{51189}{32000}\right).
\end{aligned}
\tag{101}
\]
The four positive terms are respectively the logarithmic
off-diagonal readings, the approximate pole, the u^2
convolution, and the u convolution.
For example, if V(x)=integral_I |x-y|p(y)dy, then V''=2p and
\[
V(x)=\phi_0L^2\left(
\frac{461}{280}s^2-\frac{23}{56}s^4+\frac7{120}s^6
\right)+\text{constant}.
\]
This verifies the last term by direct projection.
Inserting the two endpoints of (87) into (101) gives
\[
-\frac{68323}{120960000}
<\sqrt5\,r^{(0)}_2
<\frac{124369}{604800000}.
\]
These bounds imply
\[
|r^{(0)}_2|<\frac3{10000}.
\]

Similarly,
\[
3r^{(0)}_4=
\frac9{40}-\frac{18}{49}
+\frac15\left(-\mathfrak d+\frac{15967}{96000}
+\frac{1091}{630}\right)-\frac{17}{176000}.
\tag{102}
\]
The parentheses lie between 702/1000 and 703/1000 by (87).
The resulting rational bounds give
-21/10000<3r^(0)_4<0, and hence
\[
|r^{(0)}_4|<\frac7{10000}.
\]

### The infinite logarithmic tail

For even n>m>=0,
\[
\int_0^1P_n(s)P_m(s)w(s)\,ds
=\frac1{2(n-m)(n+m+1)}.
\tag{103}
\]
To verify it, apply the Legendre equations to the Wronskian
(1-s^2)(P_n'P_m-P_m'P_n), integrate against w, and use
w'=s/(2(1-s^2)). The endpoint Wronskian times w tends to zero.
The remaining integral is
(1/2)integral_0^1 s(P_n'P_m-P_m'P_n) ds.
Its value is 1/2: integrate the first term by parts and use
orthogonality of P_n to the lower-degree even polynomials
P_m+sP_m' and sP_m'. This proves (103).

The exact logarithmic coefficient of T_Lp at even n>=6 is
\[
t_n=\frac{\sqrt{2n+1}}2\left(
\frac1{n(n+1)}-\frac{8/7}{(n-2)(n+3)}
+\frac{1/5}{(n-4)(n+5)}
\right).
\]
Writing t=n(n+1), this becomes
\[
t_n=
\frac{\sqrt{2n+1}(t^2-76t+2100)}
 {35t(t-6)(t-20)}.
\]
The numerator is (t-38)^2+656>0. For t>=42 it is smaller
than (t-6)(t-20), because 1980-50t<0. Consequently
\[
0<t_n<\frac{\sqrt{2n+1}}{35n(n+1)},\qquad
t_6=\frac{2\sqrt{13}}{3465}.
\]
For even n>=8 the decreasing sequence
1/n^2-1/(n+1)^2 can be paired with its preceding odd term.
Therefore
\[
\begin{aligned}
\sum_{n=6,8,\ldots}t_n^2
&<\frac{52}{12006225}
+\frac1{1225}\,\frac12\sum_{n=7}^\infty
 \left(\frac1{n^2}-\frac1{(n+1)^2}\right)\\
&=\frac{52}{12006225}+\frac1{120050}
<\frac{127}{10000000}.
\end{aligned}
\tag{104}
\]
For the last rational comparison the two fractions are less
than 434/100000000 and 834/100000000, respectively.
Thus the logarithmic-tail norm is less than 357/100000.

The u convolution in j_0 has one additional tail coefficient,
only in degree six, equal to
343 sqrt(13)/864864000<1/600000.
It follows from V above and
integral_0^1 P_6(s)s^6 ds=16/3003.
All other main-column terms have degree at most four. Hence
\[
\|R_6r^{(0)}\|<
\frac{357}{100000}+\frac1{600000}
=\frac{2143}{600000}.
\]
Combining this with (101)--(102) gives
\[
\|r^{(0)}\|_{\mathcal D^{-1}}^2
<
\frac{20}{3}\left(\frac3{10000}\right)^2
+\frac{10}{3}\left(\frac7{10000}\right)^2
+2\left(\frac{2143}{600000}\right)^2
=\frac{9988898}{360000000000}
<\left(\frac{53}{10000}\right)^2.
\tag{105}
\]

### Restoring the full Gamma remainder

The operator mathcal R of the finite R_0 form satisfies
||mathcal R||<=2J_R<1/950 by T17.
Also ||p||=sqrt(13954/11025)<8/7, and
||mathcal D^(-1/2)||=sqrt(20/3)<8/3.
Mean subtraction contracts L2. Thus
\[
\|Q_0\mathcal Rp\|_{\mathcal D^{-1}}
<\frac83\frac1{950}\frac87
=\frac{64}{19950}<\frac{13}{4000}.
\tag{106}
\]

### Restoring the actual prime column

Let Pp=-(a_2/2)(p(x-u_2)+p(x+u_2)) on I, with zero
extension as always. The two caps are disjoint. In addition
to the affine bound used for the form, the exact cap polynomial
satisfies
\[
F_c(t)\le\frac3{35}+\frac{39}{35000}
+\frac7{50000000}<\frac7{80}
\qquad(0<=t<=d<1/50).
\]
Since |P_2|,|P_4|<=1 on [-1,1], cap integration yields
\[
\|Pp\|^2<\frac{49}{5120000},\qquad
|\langle Pp,\phi_2\rangle|<\frac{7\sqrt5}{16000},
\qquad
|\langle Pp,\phi_4\rangle|<\frac{21}{16000}.
\]
For example the norm bound is
(a_2^2/4) integral_0^d F_c(t)^2 dt
<d(7/80)^2/16. The coefficient bound for degree n is
a_2 sqrt(2n+1) integral_0^d F_c(t)dt/2.

The complete weighted norm, including its tail, is bounded by
\[
\begin{aligned}
\|Q_0Pp\|_{\mathcal D^{-1}}^2
&\le2\|Pp\|^2+\frac{14}{3}
 |\langle Pp,\phi_2\rangle|^2
+\frac43|\langle Pp,\phi_4\rangle|^2\\
&<\frac{9947}{384000000}
<\left(\frac{51}{10000}\right)^2 .
\end{aligned}
\tag{107}
\]
The low-mode weights here are essential; the prime column is
not bounded solely by its small form value (97).

### Restoring the exact pole column

Write R_c(x)=cosh(x/2)-1-x^2/8.
The difference from the approximate pole, after mean subtraction,
is
\[
e_{\rm pole}=(\ell(p)-m_0)\,Q_0(x^2/8)+\ell(p)Q_0R_c.
\]
Nonnegativity of p, the exact m_2, and
cosh(L/2)<64/63 give
\[
0\le\ell(p)-m_0\le\frac{m_2}{8}\cosh(L/2)
<\frac3{1000},\qquad \ell(p)<1.
\]
Also
\[
\|Q_0(x^2/8)\|=\frac{\sqrt{2L}L^2}{12\sqrt5}<\frac1{192},
\quad
\|Q_0R_c\|\le\|R_c\|
\le\frac{\sqrt{2L}L^4\cosh(L/2)}{384}<\frac1{24192}.
\]
Consequently
\[
\|e_{\rm pole}\|_{\mathcal D^{-1}}
<\frac83\left(\frac1{64000}+\frac1{24192}\right)
=\frac1{24000}+\frac1{9072}<\frac1{6000}.
\tag{108}
\]

The actual residual is exactly
r=r^(0)+Q_0 mathcal Rp+Q_0Pp+e_pole.
The triangle inequality in the single positive norm (100),
using (105)--(108), now gives
\[
\|r\|_{\mathcal D^{-1}}
<
\frac{53}{10000}+\frac{13}{4000}
+\frac{51}{10000}+\frac1{6000}
=\frac{829}{60000}<\frac1{72}.
\]
Squaring proves (99). Every tail and every source term has been
restored before this final comparison.

## BRT-T21. First-prime positivity and the unrestricted window through 7/20

For the actual source at the endpoint,
\[
\boxed{
s_{7/20}>\frac1{5000}-\frac1{5184}
=\frac{23}{3240000}>\frac1{150000}>0.
}
\tag{109}
\]
Indeed the admitted trial (94) has unit constant reading,
so the exact identity is
s_L=q_L[p]-<r,C^{-1}r>. Apply (95) and (99).
This is a positive lower bound for the returned scalar,
not merely a positive trial value.

The corresponding full-family coercivity can be made explicit:
\[
\boxed{
0<L\le\frac7{20},\qquad g\in D_L
\quad\Longrightarrow\quad
q_L[g]\ge\frac1{260000}\|g\|^2.
}
\tag{110}
\]
There is no mean-zero or cosh-neutrality restriction in (110).

Here are the normalization details at the endpoint. Let
rho=C^(-1)b be the original mean-zero minimizing return.
Since r=b+C(p-phi_0),
\[
\rho=C^{-1}r-(p-\phi_0).
\]
The already proved delta>2/31 and (99) imply
\[
\|C^{-1}r\|^2
\le\delta^{-1}\langle r,C^{-1}r\rangle
<\frac{31}{10368}<\left(\frac2{25}\right)^2.
\]
Also ||p-phi_0||^2=2929/11025<(13/25)^2. Thus ||rho||<3/5.
For arbitrary g=c phi_0+h, h mean-zero, source square completion
gives
\[
q_L[g]=s_Lc^2+
\|C^{1/2}(h+c\rho)\|^2
\ge s_Lc^2+\delta\|h+c\rho\|^2.
\]
On the other hand,
\[
\|g\|^2
=c^2+\|h\|^2
\le\frac{43}{25}c^2+2\|h+c\rho\|^2.
\]
The comparisons
(43/25)/260000<1/150000 and
2/260000<2/31 prove (110) at L=7/20.

A smaller supported smooth profile has exactly the same source
value and norm after zero extension into the endpoint window.
Every prime displacement outside its own support contributes
zero correlation. The endpoint estimate therefore holds on
every smaller smooth core, and source-form closure extends it
to D_L. This is the usual support inclusion of the actual source;
no dilation law, change of observer, or extra prime assumption
is introduced.

**Scope and gate accounting.** The first-prime positivity gate
at L=7/20 is closed by the written estimates (95), (99), and
(109). The unrestricted positive interval extends to 7/20,
with the explicit margin in (110). Earlier stronger margins on
smaller intervals remain available. The particular nonpolynomial
test (93) need not be evaluated to obtain this result; its sign
is not being claimed.

All-width positivity, and therefore the RH task, remains open.
The positive scalar does not make the complete residual zero:
the polynomial obstruction in T9 still applies. The endpoint
logarithm, all higher even modes, the signed prime term and
the source's low/tail couplings were retained in the estimates.
No new independent assumption is added. These are written
proofs in the working source chain, not independently validated
certificates. No scientific computation, numerical simulation,
certificate replay, subagent or workflow dispatch is performed.


## Scope and continuation

The current quantitative result is (109): a strictly positive
returned scalar at L=7/20, obtained by subtracting an upper bound
for the complete inverse residual from a lower bound for the
actual source value. Equation (110) extends full-family coercivity
to every support window 0<L<=7/20. The prior DWE margin 1/1700
on L<=1/4+1/47600 remains a stronger local margin on its smaller
interval.

The original quadratic trial remains useful historical evidence.
T8, T13--T15 and T17 identify specific conservative certificates
that fail; they do not prove a negative source. T9's finite
polynomial exact-cancellation obstruction still applies to (94).
T16's admitted nonpolynomial correction and its exact gain remain
valid, while the sign of its particular stronger-margin test (93)
has not been assessed. That sign is no longer necessary for the
first-prime result, which follows from T19--T21.

The successful quartic trial is not the endpoint-preserving
one-parameter family of T10. Its nonzero endpoint value restores
a logarithmic source column; the explicit infinite-tail estimate
(104) shows why this is admissible and quantitatively manageable.
Every signed source term, including prime 2, is included in (99).
Small prime correlation alone is not used as a residual bound.

This result concerns the actual arithmetic form under its stated
domain and normalization. It does not identify the trial polynomial
or the positive comparison operator with native EMK curvature,
Riemann curvature, a physical field or a universal cut law.
A co-carried observer change preserves the full source pairing;
it does not supply the new numerical margins by itself.

The next research task is extension beyond 7/20 with quantitative
control of the actual signed source, retaining its memory and
any further active prime powers. No all-width control or RH
closure follows from the finite interval proved here.

Only this note and its README entry change. Validation consists
of written polynomial/logarithmic integration, rational comparisons,
source polarization, inverse form order and immutable-content
readback. No scientific computation, numerical simulation,
independent agent, certificate replay or workflow dispatch is
performed. The commit uses [skip ci]. The PR remains a draft
stacked on #137 and is not merged.

**Gate accounting.** Closed by the working written proof:
first-prime positivity at L=7/20 and unrestricted positivity
through that window. Open: all-width positivity and RH.
The particular test (93) remains unassessed. No new independent
assumption is added.
