# BRT-1 — Full-source positivity through L=3/4

Monty · 18 October 2026. Updated 23 October 2026 through BRT-T68.

The current working computer-assisted result is
\[
 \boxed{
 0<L\le\frac34
 \quad\Longrightarrow\quad
 q_L[g]\ge9.0\cdot10^{-18}\,\|g\|^2
 \quad(g\in D_L).
 }
\]
It covers the complete admitted real even family, with no neutral
moment or finite-mode restriction, and the three active prime
powers 2, 3, 4 at their actual weights. Earlier stronger margins on
smaller intervals remain valid; in particular q_L[g] >= 10^(-11)||g||^2
for 0<L<=11/20 (T63).

T56--T63 supplied the complement bound and the actual returned
matrix at L=11/20. T64--T68 now do the same at L=3/4 by two
independent complement routes (T53's three-prime graph comparison
and the new width-dependent boundary-energy comparison T64), after
T65 locates the returned-sign question at this width below
10^(-11). Both actual prime terms, the prime power 4, the
logarithmic boundary term and the complete infinite residual are
retained.

The [reproducible packet](proof_lab/rh_l075/README.md) contains the
standard-library verifier, the fixed trials, the configuration of
each run and the 220-digit evidence. The L=11/20 packet
[proof_lab/rh_l055](proof_lab/rh_l055/README.md) is unchanged and is
reproduced by the new verifier as a control. Historical sections
below retain their status at the time of writing; their statements
that 11/20 or 3/4 is unresolved are superseded by T60--T68.

Widths beyond 3/4 and all-width positivity/RH remain open.
No new independent mathematical assumption is introduced.
The result is a working computer-assisted proof, pending external
review; it is not a formal proof-assistant certification.

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


## BRT-T22. The active-prime comparison persists on an explicit width band

Write
\[
L_0=\frac7{20},\qquad
L_{\rm cap}=\frac{7001}{20000}=L_0+\frac1{20000},
\qquad h=L-L_0.
\tag{111}
\]
Throughout T22--T24, assume L_0<=L<=L_cap unless stated
otherwise. At each width use its own normalized Legendre
functions phi_n^(L), mean projection Q_0^(L), actual operator
T_L and mean-zero operator C_L. The comparison
\[
\mathcal D_L=\frac3{20}P_2^{(L)}
+\frac3{10}P_4^{(L)}+\frac12R_6^{(L)}
\]
has the same coefficients as (100).

The actual mean-zero source satisfies
\[
\boxed{C_L\succeq\mathcal D_L+\frac1{1250}I
\succeq\mathcal D_L.}
\tag{112}
\]
This is derived at each width from the actual source. It does
not presume that prime translation is continuous in the L2
operator norm, or that support dilation is a change of observer.

### The same boundary absorption covers the moving prime caps

Exactly prime 2 is active: 2L<1<log3, while log2<2L_0.
Also L<log2. The retained bound log2>6931/10000 gives
\[
99L\le\frac{693099}{20000}
<\frac{6931}{200}=50\frac{6931}{10000}<50\log2.
\]
Hence d_L=(2L-log2)/L<1/50, and the normalized prime
caps have |x|/L>49/50. OPC-T2's translated-pair estimate
and its actual boundary absorption consequently apply
throughout this band:
\[
a_2|C_g(\log2)|
\le\frac{a_2}{2}\int_{\mathcal C_L}|g|^2,\qquad
w_L-\frac{a_2}{2}\mathbf1_{\mathcal C_L}
\ge\frac{(x/L)^2}{4-2(x/L)^2}.
\]
In particular the signed prime is retained before it is bounded.

### Scalar drift leaves room in the mean-zero comparison

GRR-T2 with a=1/2 gives, on mean-zero h_0,
\[
q_L[h_0]\ge B_*(L)\|h_0\|^2+J_L[h_0]
+\int_Iv_Lh_0^2+\ell(h_0)^2,
\]
where
\[
B_*(L)=\frac{\Omega_0}{2}+F(2L)+\frac12\log2
+\frac{47L}{96}-\frac{5L^3}{36}.
\]
The omitted Gamma remainder is nonnegative. The last two
ordinary-moment terms in GRR's comparison vanish because
h_0 has mean zero.

The kernel k is decreasing. For t>=7/10,
\[
k(t)\le k(7/10)
<\frac{e^{7/20}}{7/5}
\le\frac{100}{91}<2.
\]
Here sinh(t)>t and e^z<=1/(1-z) for 0<=z<1 suffice.
Thus F(2L)>=F(2L_0)-4h. The polynomial part of B_* is
increasing on this band, since its derivative is
47/96-5L^2/12>0. Consequently
\[
B_*(L)+3/4\ge\delta-4h.
\]
Retain the higher Legendre energies and the two odd readings
of (x/L)h_0 exactly as in T11 and T13. This yields
C_L>=K_L-4hI, where K_L has T13's same low matrix A and
tail coefficient delta_6 in the width-L basis.

For a fully rational slack estimate,
\[
K_L\succeq\mathcal D_L+\frac1{1000}I.
\]
Indeed a>19/100, d>41/100 and c^2=1/245.
The two diagonal differences exceed 39/1000 and 109/1000,
and their determinant exceeds
\[
\frac{4251}{1000000}-\frac1{245}
=\frac{41495}{245000000}>0.
\]
Also delta_6>501/1000. Finally
1/1000-4h>=1/1000-4/20000=1/1250, proving (112).
The form-core argument and bounded boundary weight extend
this comparison to the complete mean-zero form domain.

## BRT-T23. The quartic source and full residual have uniform width estimates

Use the same shape as (94), now with the correct normalization
at each width:
\[
p_L(x)=\frac1{\sqrt{2L}}f(x/L),\qquad
f(s)=1-\frac87P_2(s)+\frac15P_4(s),\qquad
r_L=Q_0^{(L)}T_Lp_L.
\]
ASM admits each zero-extended polynomial to the corresponding
operator domain. Its mean reading is one and its norm squared
is N=13954/11025. We prove throughout (111)
\[
\boxed{
q_L[p_L]>\frac{1013}{5000000}-6h,\qquad
\langle r_L,C_L^{-1}r_L\rangle
<\left(\frac{829}{60000}\right)^2<\frac1{5184}.
}
\tag{113}
\]
Every residual component is bounded before the inverse
comparison is applied.

### Trial-value drift

Let mathfrak d(L)=-Omega_0/2-F(2L), and let
J_0(L)=integral_0^(2L)(1/4-u/48-u^2/32)du.
The scalar multiplier in T19 is -mathfrak d(L)+J_0(L).
Since j_0>0 on [0,2L_cap], J_0 is increasing, while
mathfrak d(L)-mathfrak d(L_0)<=4h. Hence that multiplier
loses at most 4h.

The singular boundary energy (96) is scale invariant for
this fixed normalized shape. The positive lower pole terms,
the u and u^2 convolution contributions in T19 are increasing
functions of L. The constant convolution contributes -L/4,
losing h/4. The same prime bound 27/1000000 applies because
d_L<1/50, and R_0 remains nonnegative.
Therefore the full lower estimate (98) at width L exceeds
its lower estimate at L_0 minus (4N+1/4)h.
Since N<4/3 and 4N+1/4<6, this proves the first part of
(113), including its strict endpoint inequality.

### The main column, including its infinite logarithmic tail

Use the width-L version of T20's main residual r_L^(0):
the exact scalar and logarithmic multiplier, j_0, and the
approximate pole column m_0(1+x^2/8).
Its exact second and fourth readings are
\[
\begin{aligned}
R_2(L):=\sqrt5\,r_{L,2}^{(0)}
={}&\frac{19}{42}+\frac{3L^3}{16}
+\frac{491L^2}{52920}\\
&-\frac87\left(-\mathfrak d(L)+J_0(L)+\frac{43}{30}\right),\\
R_4(L):=3r_{L,4}^{(0)}
={}&\frac9{40}-\frac{18}{49}
+\frac15\left(-\mathfrak d(L)+J_0(L)+\frac{1091}{630}\right)
-\frac{17L^2}{21560}.
\end{aligned}
\tag{114}
\]
They reduce to (101)--(102) at L_0.

Here 0<mathfrak d'(L)-J_0'(L)
=2(k(2L)-j_0(2L))<4: GRR gives k-j_0>0 and
j_0>0, while k(2L)<2. Differentiating (114), and using
L<3/8, therefore gives
\[
0<R_2'(L)<5,\qquad -1<R_4'(L)<0.
\]
The endpoint enclosures already established in T20 imply
\[
-\frac{57}{100000}<R_2(L_0)<\frac{21}{100000},\qquad
-\frac{41}{20000}<R_4(L_0)<0.
\]
For the fourth-mode lower bound, (102)'s parenthesis exceeds
702/1000; also 18/49<36735/100000 and
17/176000<1/10000. These give the stated -41/20000 bound
directly. It follows for 0<=h<=1/20000 that
\[
|R_2(L)|<\frac3{5000},\qquad
-\frac{21}{10000}<R_4(L)<0.
\]
Since sqrt(5)>2,
\[
|r_{L,2}^{(0)}|<\frac3{10000},\qquad
|r_{L,4}^{(0)}|<\frac7{10000}.
\]

The entire logarithmic coefficient sequence t_n in (103)--(104)
is unchanged under this normalized shape scaling. Its tail norm
is still less than 357/100000. The sole additional degree-six
coefficient from the u convolution is
\[
\frac{L^2\sqrt{13}}{308880}<\frac1{617760}<\frac1{600000},
\]
because L<351/1000, L^2<1/8 and sqrt(13)<4.
Thus the same complete main-column estimate (105) holds:
\[
\|r_L^{(0)}\|_{\mathcal D_L^{-1}}<\frac{53}{10000}.
\tag{115}
\]

### All three omitted columns remain within the same bounds

For the full Gamma remainder, the T8 series argument gives
R_0(u)<1/200 on 0<=u<=2L_cap. To check a uniform bound,
u^2<1/2 implies cosh(u)<4/3, and hence
\[
R_0(u)\le u^3/216+u^4/144
<1/512+9/4096=17/4096<1/200,
\]
where u<3/4 was used in the second comparison. Consequently
\[
J_R(L)\le\frac{55223}{108000000}+\frac h{100}
<\frac1{1950}+\frac1{2000000}<\frac1{1900}.
\]
The complete remainder operator still has norm <1/950.
The shape norm is unchanged, so (106) gives
\[
\|Q_0^{(L)}\mathcal R_Lp_L\|_{\mathcal D_L^{-1}}
<\frac{13}{4000}.
\]

The actual prime caps are disjoint and d_L<1/50. Their
normalized cap polynomial is exactly the same F_c as in
T19--T20. All norm and low-mode estimates in (107) therefore
apply without changing their constants:
\[
\|Q_0^{(L)}Pp_L\|_{\mathcal D_L^{-1}}<\frac{51}{10000}.
\]
This includes the entire prime-column tail. It does not
replace the prime column by its small quadratic-form reading.

For the pole correction, L^2<1/8 gives
cosh(L/2)<64/63 and sqrt(2L)<6/7. Hence
ell(p_L)<1, and the exact second moment gives
\[
0\le\ell(p_L)-m_0
\le\frac{m_2}{8}\cosh(L/2)
<\frac{19}{6615}<\frac3{1000}.
\]
The same elementary variance and Taylor estimates give
\[
\|Q_0^{(L)}(x^2/8)\|<1/192,\qquad
\|Q_0^{(L)}R_c\|<1/24192.
\]
Thus (108) remains valid:
\[
\|e_{{\rm pole},L}\|_{\mathcal D_L^{-1}}<\frac1{6000}.
\]

All four terms sum to the actual r_L. The triangle inequality
in this one fixed positive comparison norm gives
\[
\|r_L\|_{\mathcal D_L^{-1}}
<\frac{53}{10000}+\frac{13}{4000}
+\frac{51}{10000}+\frac1{6000}
=\frac{829}{60000}<\frac1{72}.
\tag{116}
\]
Use (112) and inverse form order to finish (113).
The endpoint logarithm, all even tail modes and the signed
prime column have survived the entire width comparison.

## BRT-T24. A positive full-family interval beyond 7/20

Combining (113) with the exact returned-scalar identity gives
the varying bound
\[
\boxed{
s_L>
\frac{1013}{5000000}-\frac1{5184}-6(L-L_0)
\quad(L_0\le L\le L_{\rm cap}).
}
\tag{117}
\]
In particular, put
\[
L_+=\frac{350001}{1000000}=\frac7{20}+\frac1{1000000}.
\]
For L_0<=L<=L_+,
\[
\boxed{
s_L>\frac{983}{5000000}-\frac1{5184}
=\frac{749}{202500000}>\frac1{300000}>0.
}
\tag{118}
\]
The corresponding unrestricted source estimate is
\[
\boxed{
0<L\le L_+,\quad g\in D_L
\quad\Longrightarrow\quad
q_L[g]\ge\frac1{520000}\|g\|^2.
}
\tag{119}
\]

For the norm conversion, (112) implies C_L>2I/31.
Together with (113), the same T21 argument gives
||C_L^(-1)r_L||<2/25 and
||p_L-phi_0^(L)||<13/25. Therefore
rho_L=C_L^(-1)b_L=C_L^(-1)r_L-(p_L-phi_0^(L))
has norm <3/5. For g=c phi_0^(L)+h_0 and z=h_0+c rho_L,
\[
q_L[g]=s_Lc^2+\|C_L^{1/2}z\|^2,\qquad
\|g\|^2\le(43/25)c^2+2\|z\|^2.
\]
The rational comparisons
(43/25)/520000<1/300000 and
2/520000<2/31 prove (119) on [L_0,L_+].
Support inclusion and source-form closure give every smaller
window. The stronger previous margin 1/260000 through L_0,
and the earlier DWE margin on its smaller interval, remain valid.

This is a small explicit extension with prime 2 already active.
It consumes a finite proved margin. It cannot be iterated
indefinitely without a new quantitative source estimate.
Indeed at the larger comparison-band endpoint L_cap, the
right side of (117) is
\[
-\frac{487}{5000000}-\frac1{5184}<0.
\tag{120}
\]
Thus this particular uniform lower bound does not prove
positivity at L_cap. Equation (120) is not an upper bound
for s_L and does not imply that the actual source is negative.

**Gate accounting.** The working written estimates extend
unrestricted positivity beyond 7/20 to L_+. The complete
residual is retained. Positivity farther out, including at
L_cap using a sharper estimate, and all-width positivity/RH
remain unresolved. No new independent assumption is added.
No independent certification, scientific computation, numerical
simulation, certificate replay, subagent or workflow dispatch
is performed.

## BRT-T25. The discarded Gamma remainder supplies a uniform positive margin

Keep the width band L_0=7/20<=L<=L_cap=7001/20000
and the actual source conventions of T22--T24. The finite
remainder operator used in T20 and T23 has quadratic form
\[
\langle g,\mathcal R_Lg\rangle
=\int_0^{2L}R_0(u)(\|g\|^2-C_g(u))\,du.
\]
It is already a bounded positive operator. Here its lower
bound is made quantitative:
\[
\boxed{
\langle g,\mathcal R_Lg\rangle
\ge c_R(L)\|g\|^2,\qquad
c_R(L)=\frac{7L^4}{8192}+\frac{31L^5}{4096}
>\frac1{20000}.
}
\tag{121}
\]
This is a contribution of the existing arithmetic kernel,
not an additional positive source assumption.

### Two positive kernel coefficients suffice

GRR-T1 gives R_0(u)=E(u)/(2 sinh u), with positive entire
coefficients of E from degree four onward. Its first two
terms are exactly
\[
E(u)\ge\frac7{5760}u^4+\frac5{768}u^5.
\]
On 0<u<=2L_cap, T23 gives cosh(u)<4/3, while
sinh(u)<=u cosh(u). Consequently
\[
\boxed{
R_0(u)\ge\frac7{15360}u^3+\frac5{2048}u^4.
}
\tag{122}
\]
The equality at u=0 is understood by continuity.

For L<=u<=2L, the two translated overlap caps are disjoint.
Thus |C_g(u)|<=||g||^2/2 and
||g||^2-C_g(u)>=||g||^2/2. Retaining only this displacement
range and then integrating (122) gives
\[
\begin{aligned}
\langle g,\mathcal R_Lg\rangle
&\ge\frac{\|g\|^2}{2}
\int_L^{2L}\left(\frac7{15360}u^3+\frac5{2048}u^4\right)\,du\\
&=\left(\frac{7L^4}{8192}+\frac{31L^5}{4096}\right)\|g\|^2.
\end{aligned}
\]
The coefficient increases with L. At L_0,
\[
c_R(L_0)=\frac{689087}{13107200000}
>\frac{655360}{13107200000}=\frac1{20000}.
\]
This proves (121), including on the full source form domain;
the remainder itself is bounded on L2.

For the fixed quartic shape, N=13954/11025>5/4.
Therefore the actual trial's previously discarded contribution
has the uniform bound
\[
\boxed{\langle p_L,\mathcal R_Lp_L\rangle>\frac1{16000}.}
\tag{123}
\]
T19 and T23 dropped this entire term in their trial-value
lower bounds. It can be added there once. Their complete
residual bounds already include the full remainder column;
those bounds are kept unchanged below.

## BRT-T26. Keeping the source cancellation reduces the width loss by a factor of 200

Define a particular lower comparison for the same normalized
quartic trial:
\[
\begin{aligned}
\mathcal B(L)={}&
N\bigl(-\mathfrak d(L)+J_0(L)\bigr)
+\frac{194287}{283500}+\frac{7L}{4}
+\frac{57L^3}{560}\\
&+\frac{361L^5}{352800}
+\frac{487L^2}{48000}-\frac{27}{1000000}.
\end{aligned}
\tag{124}
\]
Here mathfrak d(L)=-Omega_0/2-F(2L) and
J_0(L)=integral_0^(2L)j_0(u)du, exactly as in T23.
For every width in the stated band,
\[
q_L[p_L]\ge\mathcal B(L)
+\langle p_L,\mathcal R_Lp_L\rangle.
\]
This uses the same actual prime bound as (97), valid uniformly
because d_L<1/50. The logarithmic energy is (96).
The pole lower bound is
2L+19L^3/210+361L^5/352800; combining it with the
constant and quadratic convolution terms gives 7L/4
and 57L^3/560. The linear convolution lower bound is
487L^2/48000. Thus (124) retains exactly the terms already
used in T19, without rounding the positive width dependence
away. In particular mathcal B(L_0)>1013/5000000.

The exact kernel identity gives
\[
\frac{d}{dL}(-\mathfrak d(L)+J_0(L))
=-2(k(2L)-j_0(2L))
=-\frac1{2L}-2R_0(2L).
\]
Hence
\[
\begin{aligned}
\mathcal B'(L)={}&
-\frac{N}{2L}-2N R_0(2L)+\frac74
+\frac{171L^2}{560}\\
&+\frac{361L^4}{70560}+\frac{487L}{24000}.
\end{aligned}
\tag{125}
\]
This derivative is of the declared scalar lower comparison.
No derivative formula for the full source inverse, no operator
norm continuity of prime translation, and no sign of the
actual s_L derivative are asserted.

The rational bounds
\[
\frac54<N<\frac{19}{15}<\frac43,\qquad
R_0(2L)<\frac1{200},\qquad L\ge\frac7{20}
\]
now suffice. Retain the positive quadratic and linear terms
in (125), and drop only the positive L^4 term:
\[
\begin{aligned}
\mathcal B'(L)
&>-\frac{38}{21}-\frac1{75}+\frac74
+\frac{1197}{32000}+\frac{3409}{480000}\\
&=-\frac{23813}{840000}>-\frac3{100}.
\end{aligned}
\]
The bound on R_0 here is its already proved upper bound:
the negative endpoint-remainder term in (125) has not
been discarded.

Integrate on [L_0,L] and add the positive remainder (123).
The actual trial value therefore satisfies
\[
\boxed{
q_L[p_L]>
\frac{1013}{5000000}+\frac1{16000}
-\frac3{100}(L-L_0)
=\frac{2651}{10000000}-\frac3{100}(L-L_0).
}
\tag{126}
\]
The previous loss coefficient was 6; it is now 3/100.
This comparison improves the full returned-scalar lower bound
because the complete inverse-residual penalty from (113)
remains valid for the same actual p_L and r_L.

The improvement does not follow from changing sign conventions
or rotating the observer. It follows from retaining the positive
source contribution and the cancelling width derivatives in
the actual scalar and pole terms.

## BRT-T27. The whole former comparison band is now positive

Use the unchanged complete residual estimate (113) together
with (126). Since L-L_0<=1/20000,
\[
\begin{aligned}
s_L
&=q_L[p_L]-\langle r_L,C_L^{-1}r_L\rangle\\
&>\frac{2651}{10000000}
-\frac3{2000000}-\frac1{5184}\\
&=\frac{659}{2500000}-\frac1{5184}
>\frac7{100000}.
\end{aligned}
\tag{127}
\]
For the last rational comparison,
659/2500000-7/100000=121/625000, and
121*5184=627264>625000.
Thus the entire band [7/20,7001/20000], including the endpoint
where the old lower estimate (120) failed, has positive
returned scalar in the working written proof.

The full-family conclusion is
\[
\boxed{
0<L\le\frac{7001}{20000}
=\frac7{20}+\frac1{20000},\qquad g\in D_L
\quad\Longrightarrow\quad
q_L[g]\ge\frac1{25000}\|g\|^2.
}
\tag{128}
\]
There is no mean-zero or cosh-neutrality restriction.

For completeness, T24's normalization remains unchanged:
C_L>2I/31, ||rho_L||<3/5, and for
g=c phi_0^(L)+h_0, z=h_0+c rho_L,
\[
q_L[g]=s_Lc^2+\|C_L^{1/2}z\|^2,\qquad
\|g\|^2\le(43/25)c^2+2\|z\|^2.
\]
Now (43/25)/25000<7/100000 and
2/25000<2/31. These comparisons prove (128) on the band.
Zero extension of smaller smooth profiles into its endpoint,
followed by the retained source-form closure, proves it for
every smaller L. The previous BRT full-family margins on their smaller
intervals are superseded by this stronger bound. DWE's
1/1700 margin remains stronger on its much smaller interval.

The old bound (120) still correctly records failure of that
particular conservative estimate. Its failure no longer leaves
the source sign at L_cap open: (127) resolves it with a stronger
lower bound. The complete residual remains nonzero, with all
higher even modes, the Gamma remainder, the pole correction
and the signed prime-2 column included.

**Scope and gate accounting.** The formerly open larger-band
positivity gate at L_cap is closed by this working written
derivation. The unrestricted positive interval is (128).
Beyond that endpoint, the cap-width and residual constants
must be established again; none are extrapolated here.
All-width positivity and RH remain open. No new independent
assumption is added. No independent certification, scientific
computation, numerical simulation, certificate replay, subagent
or workflow dispatch is performed.

## BRT-T28. A sharper full Gamma-remainder bound on a wider band

For T28--T31 put
\[
L_0=\frac7{20},\qquad L_1=\frac{351}{1000},
\qquad 0\le h=L-L_0\le\frac1{1000}.
\]
Use the same admitted quartic shape p_L, with its width-L
normalization, and the actual residual r_L=Q_0^(L)T_Lp_L.
The comparisons D_L and K_L have the meanings in T22.
None of the old d_L<1/50 estimates is extrapolated here.

For 0<=u<=2L_1, the actual kernel remainder obeys
\[
\boxed{
0\le R_0(u)\le\frac5{6912}u^3+\frac{79}{23040}u^4
<\frac1{800}.
}
\tag{129}
\]
The strict final inequality includes u=0.

To derive it, separate the exact first two coefficients from
GRR-T1's positive entire series:
\[
E(u)=\frac7{5760}u^4+\frac5{768}u^5+E_{\ge6}(u).
\]
The coefficient bounds used in T8, now applied only for m>=3,
give
\[
0\le E_{\ge6}(u)
\le\left(\frac u{24}+\frac{u^2}{16}\right)
\left(\sinh u-u-\frac{u^3}{6}\right).
\]
Since sinh(u)>=u and the last bracket is at most
u^5 cosh(u)/120,
\[
R_0(u)\le\frac7{11520}u^3+\frac5{1536}u^4
+\frac{u^5\cosh u}{5760}+\frac{u^6\cosh u}{3840}.
\]
Here u<=351/500 has u^2<1/2, hence cosh(u)<4/3.
Bounding u^5<=u^3/2 and u^6<=u^4/2 yields the polynomial
in (129). All higher kernel coefficients have been bounded.

For its final pointwise bound, use u<71/100,
u^3<9/25, u^4<13/50, 5/6912<3/4000 and
79/23040<7/2000. The resulting upper bound is
59/50000<1/800.

Let J_R(L)=integral_0^(2L)R_0(u)du.
At L_0, integrating (129) gives
\[
J_R(L_0)\le
\frac{5(7/10)^4}{27648}+\frac{79(7/10)^5}{115200}
<\frac{44+116}{1000000}=\frac1{6250}.
\]
For the added displacement range, (129) gives
\[
\boxed{
J_R(L)<\frac1{6250}+\frac{2h}{800}
\le\frac{13}{80000}<\frac1{6000}.
}
\tag{130}
\]
Thus the complete finite-remainder operator has norm <1/3000.
Since ||p_L||<8/7 and ||D_L^(-1/2)||<8/3,
\[
\boxed{
\|Q_0^{(L)}\mathcal R_Lp_L\|_{\mathcal D_L^{-1}}
<\frac{64}{63000}<\frac{11}{10000}.
}
\tag{131}
\]
This replaces the previous 13/4000 remainder-column bound
by a smaller bound, on the larger band.

## BRT-T29. Wider actual prime caps retain boundary absorption and inverse control

Exactly prime 2 is active throughout [L_0,L_1]:
log2<2L_0, while 2L_1<1<log3 and L_1<log2.
Using log2>6931/10000,
\[
d_L=\frac{2L-\log2}{L}
<\frac{89}{3510}<\frac1{39}.
\]
The new cap estimate is d_L<1/39. The old bound d_L<1/50
is not used beyond its proved interval.

### Boundary absorption on the wider caps

Write t=(x/L)^2 and
\[
W(t)=-\frac14\log(1-t)-\frac{t}{4-2t}.
\]
Then
\[
W(0)=0,\qquad
W'(t)=\frac{t^2}{(1-t)(4-2t)^2}\ge0
\quad(0\le t<1).
\]
On a prime cap, |x|/L>38/39. At t_0=(38/39)^2,
1-t_0=77/1521<1/19 and t_0/(4-2t_0)=361/799<23/50.
Moreover
\[
\frac14\log19
=\log2+\frac14\log(19/16)
>\frac{6931}{10000}+\frac3{76}.
\]
Here log(1+y)>y/(1+y), with y=3/16, supplies the last
additional term. Therefore W(t_0)>1/4.
Monotonicity and a_2/2<1/4 prove, on the entire interval,
\[
w_L-\frac{a_2}{2}\mathbf1_{\mathcal C_L}
\ge v_L,\qquad
v_L(x)=\frac{(x/L)^2}{4-2(x/L)^2}.
\]
The translated-pair estimate a_2|C_g(log2)|
<= (a_2/2) integral_(mathcal C_L)|g|^2 is unchanged.
Thus the actual signed prime is again absorbed before
forming the mean-zero comparison.

The scalar-drift argument of T22 applies on this wider band:
k(2L)<2 and 47/96-5L^2/12>0 still hold.
Retaining the same Legendre and boundary readings gives
C_L>=K_L-4hI. There is enough unused matrix slack:
\[
K_L\succeq\mathcal D_L+\frac1{200}I.
\]
Indeed a=delta+11/84>2/31+11/84>39/200, d>41/100,
and c^2=1/245. The two low diagonal differences exceed
1/25 and 21/200, with determinant greater than
\[
\frac{21}{5000}-\frac1{245}
=\frac{145}{1225000}>0.
\]
Also delta_6>101/200. Since 4h<=1/250, this proves
\[
\boxed{
C_L\succeq\mathcal D_L+\frac1{1000}I.
}
\tag{132}
\]
The comparison holds on the complete mean-zero form domain.
No commutation or invariance of actual source modes is assumed.

### The actual prime form and its complete column

The cap polynomial remains
F_c(t)=2/35+(10/7)t+(39/14)t^2-(7/2)t^3+(7/8)t^4.
For 0<=t<=1/39 it satisfies
\[
0<F_c(t)\le A+B t<\frac1{10},\qquad
A=\frac2{35},\quad B=\frac{11}{7}.
\]
For the affine bound, after factoring out t, the positive
terms are at most 10/7+1/14+7/(8*39^3)<11/7.
Its maximum affine value is
A+B/39=133/1365<1/10.

The exact overlap integral and a_2<1/2 yield
\[
\begin{aligned}
a_2C_{p_L}(\log2)
&<\frac14\left(\frac{A^2}{39}
+\frac{AB}{39^2}+\frac{B^2}{6\cdot39^3}\right)\\
&<\frac{84+60+7}{4000000}<\frac{38}{1000000}.
\end{aligned}
\tag{133}
\]
The three terms inside the parentheses are respectively
4/47775, 22/372645 and 121/17439786; the displayed outward
bounds follow by integer comparison.

For the actual prime column
Pp_L=-(a_2/2)(p_L(x-log2)+p_L(x+log2)), the caps are
disjoint. Integrating their squared profile and using
|P_2|,|P_4|<=1 gives
\[
\|Pp_L\|^2<\frac1{62400},\qquad
|\langle Pp_L,\phi_2^{(L)}\rangle|<\frac{\sqrt5}{1560},
\qquad
|\langle Pp_L,\phi_4^{(L)}\rangle|<\frac3{1560}.
\]
As in (107), include the entire tail by using its weight 2,
then add the excess low-mode weights. This gives
\[
\begin{aligned}
\|Q_0^{(L)}Pp_L\|_{\mathcal D_L^{-1}}^2
&<\frac1{31200}+\frac{106}{7300800}
=\frac{17}{365040}\\
&<\frac{47}{1000000}
<\left(\frac{69}{10000}\right)^2.
\end{aligned}
\tag{134}
\]
The prime-column bound is separate from the prime-form
bound (133); neither substitutes for the other.

## BRT-T30. The complete residual remains controlled through 351/1000

The exact low-mode expressions R_2(L),R_4(L) in (114)
remain valid for the same quartic shape at every width.
Their scalar derivative now has the sharper estimate
\[
\mathfrak d'(L)-J_0'(L)
=\frac1{2L}+2R_0(2L)
<\frac{10}{7}+\frac1{400}.
\]
Differentiate (114). Using L<3/8 gives
\[
0<R_2'(L)<\frac74,\qquad -\frac3{10}<R_4'(L)<0.
\]
For explicit upper comparisons, the first derivative is less
than
81/1024+1473/211680+(8/7)(10/7+1/400)<7/4;
the magnitude of the second is less than
(1/5)(10/7+1/400)+51/86240<3/10.

The retained endpoint enclosures then give, for h<=1/1000,
\[
-\frac{57}{100000}<R_2(L)<\frac{49}{25000},\qquad
-\frac{47}{20000}<R_4(L)<0.
\]
Since sqrt(5)>11/5, it follows that
\[
|r_{L,2}^{(0)}|<\frac9{10000},\qquad
|r_{L,4}^{(0)}|<\frac8{10000}.
\]

The full logarithmic tail t_n is unchanged by normalized
shape scaling. Its norm is still below 357/100000.
The additional degree-six coefficient is
L^2 sqrt(13)/308880<1/600000, since L^2<1/8.
Thus ||R_6 r_L^(0)||<2143/600000, and
\[
\begin{aligned}
\|r_L^{(0)}\|_{\mathcal D_L^{-1}}^2
&<\frac{20}{3}\left(\frac9{10000}\right)^2
+\frac{10}{3}\left(\frac8{10000}\right)^2
+2\left(\frac{2143}{600000}\right)^2\\
&=\frac{11896898}{360000000000}
<\left(\frac{29}{5000}\right)^2.
\end{aligned}
\tag{135}
\]
This includes every higher even logarithmic coefficient.

For the pole correction, the premises used in T23 all remain
valid on the new band: L^2<1/8, sqrt(2L)<6/7,
cosh(L/2)<64/63, and the same positive shape and exact m_2.
They give ell(p_L)<1, ell(p_L)-m_0<3/1000,
||Q_0(x^2/8)||<1/192 and ||Q_0R_c||<1/24192.
Consequently the complete pole-correction norm is still
less than 1/6000.

Restore all four columns:
r_L=r_L^(0)+Q_0 mathcal R_Lp_L+Q_0Pp_L+e_(pole,L).
Equations (131), (134), (135) and the pole estimate give
\[
\begin{aligned}
\|r_L\|_{\mathcal D_L^{-1}}
&<\frac{29}{5000}+\frac{11}{10000}
+\frac{69}{10000}+\frac1{6000}
=\frac{419}{30000}.
\end{aligned}
\]
With (132) and inverse form order,
\[
\boxed{
\langle r_L,C_L^{-1}r_L\rangle
<\left(\frac{419}{30000}\right)^2<\frac1{5100}.
}
\tag{136}
\]
For the last comparison,
175561*5100=895361100<900000000.
The remainder column is smaller than before; the wider
prime caps and changed low readings are explicitly allowed
larger bounds. The complete source residual is retained.

## BRT-T31. Full-source positivity extends to 351/1000

The lower kernel estimate (122) also holds on this larger
band, since cosh(u)<4/3 still holds for u<=2L_1.
Its integrated coefficient c_R(L) increases with L.
Thus the actual trial remainder continues to exceed 1/16000.

Let mathcal B_1(L) be the explicit function (124) with its
constant prime allowance 27/1000000 replaced by 38/1000000.
Equation (133) now proves
q_L[p_L]>=mathcal B_1(L)+<p_L,mathcal R_Lp_L>.
In particular
\[
\mathcal B_1(L_0)>
\frac{1013}{5000000}-\frac{11}{1000000}.
\]
This is a newly justified lower comparison on the wider band,
not an extrapolation of the old prime allowance.

Its derivative is still (125), since only a constant allowance
has changed. Now (129) bounds R_0(2L) by 1/800. The same
positive source terms therefore give
\[
\begin{aligned}
\mathcal B_1'(L)
&>-\frac{38}{21}-\frac1{300}+\frac74
+\frac{1197}{32000}+\frac{3409}{480000}\\
&=-\frac{15413}{840000}>-\frac1{50}.
\end{aligned}
\]
The negative endpoint-remainder contribution is bounded,
not discarded. Integrating and adding the actual remainder
lower bound yields
\[
\boxed{
q_L[p_L]>
\frac{2541}{10000000}-\frac{L-L_0}{50}
\ge\frac{2341}{10000000}.
}
\tag{137}
\]
Combine (137) with the complete residual estimate (136):
\[
\boxed{
s_L>\frac{2341}{10000000}-\frac1{5100}
>\frac3{80000}
\quad(L_0\le L\le L_1).
}
\tag{138}
\]
For the final rational comparison,
2341/10000000-3/80000=983/5000000 and
983*5100=5013300>5000000.

The full-family estimate is
\[
\boxed{
0<L\le\frac{351}{1000},\quad g\in D_L
\quad\Longrightarrow\quad
q_L[g]\ge\frac1{50000}\|g\|^2.
}
\tag{139}
\]
It has no mean-zero or cosh-neutrality restriction.
Indeed C_L>2I/31 and (136) imply
\[
\|C_L^{-1}r_L\|^2<\frac{31}{10200}
<\left(\frac2{25}\right)^2.
\]
Since ||p_L-phi_0^(L)||<13/25, the original return rho_L
still has norm <3/5. Apply the same exact square completion
and norm comparison as T27:
\[
q_L[g]=s_Lc^2+\|C_L^{1/2}z\|^2,\qquad
\|g\|^2\le(43/25)c^2+2\|z\|^2.
\]
Now (43/25)/50000<3/80000 and
2/50000<2/31. This proves (139) on [L_0,L_1].
Support inclusion and the existing source-form closure give
every smaller width. The stronger prior margin 1/25000
through 7001/20000 remains available, as does DWE's stronger
margin on its smaller interval.

**Scope and gate accounting.** The working written result
extends unrestricted positivity through 351/1000. The wider
prime caps, mean-zero comparison, source-value lower bound
and complete inverse-residual upper bound have all been
established on that same band. No residual is made zero.
Positivity beyond this finite interval and all-width positivity/RH
remain open; the constants are not extrapolated further.
No new independent assumption is added. No independent
certification, scientific computation, numerical simulation,
certificate replay, subagent or workflow dispatch is performed.

## BRT-T32. Sharper scalar drift retains the inverse comparison through 44/125

For T32--T35 put
\[
L_0=\frac7{20},\qquad L_2=\frac{44}{125},
\qquad 0\le h=L-L_0\le\frac1{500}.
\]
Keep the same admitted normalized quartic p_L from T23,
the actual source T_L and residual r_L=Q_0^(L)T_Lp_L.
All comparisons below apply on this entire band.

Only prime 2 is active: log2<2L_0 and
2L_2=88/125<1<log3. Also L_2<log2, so the prime caps
remain disjoint. Their normalized width satisfies
\[
d_L=\frac{2L-\log2}{L}
<\frac{109}{3520}<\frac1{32}.
\]
These are newly justified caps, not an extrapolation of
the former 1/39 cap bound.

### Boundary absorption

Use the increasing difference W(t)=w(t)-t/(4-2t) from T29.
On a prime cap, |x|/L>31/32. At t_0=(31/32)^2,
\[
1-t_0=\frac{63}{1024},\qquad
\frac{t_0}{4-2t_0}=\frac{961}{2174}<\frac49.
\]
Since log(1+y)>y/(1+y) for y>0,
\[
W(t_0)>
\log2+\frac14\log(64/63)-\frac49
>\frac{6931}{10000}+\frac1{256}-\frac49
>\frac14.
\]
Thus, using a_2/2<1/4, the same actual-prime absorption
w_L-(a_2/2)1_(mathcal C_L)>=v_L holds on the wider caps.
Away from the caps it follows from W>=0. This retains the
signed prime before forming the mean-zero comparison.

### The scalar loss is less than 2h

The scalar in T22 is
\[
B_*(L)=\frac{\Omega_0}{2}+F(2L)+\frac12\log2
+\frac{47L}{96}-\frac{5L^3}{36}.
\]
Keep its positive polynomial derivative when bounding its loss.
Since k(2L)<100/91 and L^2<1/8,
\[
B_*'(L)
=-2k(2L)+\frac{47}{96}-\frac{5L^2}{12}
>-\frac{200}{91}+\frac7{16}>-2.
\]
Consequently B_*(L)+3/4>=delta-2h. Retaining T13's
higher Legendre and boundary readings gives
C_L>=K_L-2hI, with K_L expressed in the width-L basis.
The matrix slack K_L>=D_L+I/200 proved in T29 is
independent of L. Therefore
\[
\boxed{
C_L\succeq\mathcal D_L+\frac1{1000}I
\quad(L_0\le L\le L_2).
}
\tag{140}
\]
Indeed 1/200-2/500=1/1000. The same form-core argument
extends this bound to the complete mean-zero form domain.
No continuity in operator norm of prime translation is assumed.

### The full Gamma and pole upper bounds persist

For u<=2L_2, one still has u^2<1/2 and u<71/100.
Every premise of the coefficient proof of (129) therefore
holds on this larger displacement interval. In particular
R_0(u)<1/800 and
\[
J_R(L)<\frac1{6250}+\frac{2h}{800}
\le\frac{33}{200000}<\frac1{6000}.
\]
It follows exactly as in (131) that the complete remainder
column obeys
\[
\|Q_0^{(L)}\mathcal R_Lp_L\|_{\mathcal D_L^{-1}}
<\frac{11}{10000}.
\]
The pole proof also retains all its premises:
L^2<1/8, sqrt(2L)<6/7, cosh(L/2)<64/63, positivity
of the unchanged shape, and its exact second moment.
Hence the full pole-correction norm remains below 1/6000.
The lower kernel bound (122) remains valid as well.
These bounds have been checked on the enlarged band;
their old stated ranges are not silently extended.

## BRT-T33. Integrated cap bounds control the complete prime column

For the actual cap polynomial F_c and 0<=t<=1/32,
\[
0<F_c(t)\le A+Bt,\qquad A=\frac2{35},\quad B=\frac{11}{7}.
\]
To prove the affine bound, the coefficient after factoring
out t is at most
10/7+39/(14*32)+7/(8*32^3)<11/7.
Use the cap integrals themselves, instead of replacing
the whole cap by its maximum.

The actual prime form satisfies
\[
\begin{aligned}
a_2C_{p_L}(\log2)
&<\frac14\left(\frac{A^2}{32}
+\frac{AB}{32^2}+\frac{B^2}{6\cdot32^3}\right)\\
&=\frac14\left(\frac1{9800}+\frac{11}{125440}
+\frac{121}{9633792}\right)
<\frac{51}{1000000}.
\end{aligned}
\]
The three fractions inside the parentheses are respectively
less than 103/1000000, 88/1000000 and 13/1000000.

For the full prime column define the two upper integrals
\[
\begin{aligned}
I_1&=\int_0^{1/32}(A+Bt)\,dt
=\frac1{560}+\frac{11}{14336}<\frac{13}{5000},\\
I_2&=\int_0^{1/32}(A+Bt)^2\,dt
=\frac1{9800}+\frac{11}{125440}+\frac{121}{4816896}
<\frac{216}{1000000}.
\end{aligned}
\]
For I_1, bound the two terms by 9/5000 and 4/5000.
For I_2, bound them by 41/400000, 11/125000 and
51/2000000, whose sum is 216/1000000.

Disjointness of the caps and a_2<1/2 give
\[
\|Pp_L\|^2<\frac{I_2}{16},\qquad
|\langle Pp_L,\phi_n^{(L)}\rangle|
<\frac{\sqrt{2n+1}}4 I_1\quad(n=2,4).
\]
Using the tail weight 2 and adding the excess low-mode
weights, as in (107), yields
\[
\begin{aligned}
\|Q_0^{(L)}Pp_L\|_{\mathcal D_L^{-1}}^2
&<\frac{I_2}{8}+\frac{53}{24}I_1^2\\
&<\frac{27}{1000000}+\frac{8957}{600000000}
=\frac{25157}{600000000}\\
&<\frac{25350}{600000000}
=\left(\frac{13}{2000}\right)^2.
\end{aligned}
\]
Thus the two distinct estimates are
\[
\boxed{
a_2C_{p_L}(\log2)<\frac{51}{1000000},\qquad
\|Q_0^{(L)}Pp_L\|_{\mathcal D_L^{-1}}<\frac{13}{2000}.
}
\tag{141}
\]
The second includes every prime-column tail mode. Its bound
is smaller than (134)'s 69/10000 despite the wider caps,
because the affine profile was integrated before estimating
the norm. The small prime form has not been substituted
for the full column.

## BRT-T34. The actual quartic overlap restores more positive Gamma energy

For this same trial, a stronger lower bound than (123) is
\[
\boxed{
\langle p_L,\mathcal R_Lp_L\rangle>\frac9{100000}
\quad(L_0\le L\le L_2).
}
\tag{142}
\]
This replaces the former 1/16000 lower bound for the same
contribution; the two are not added together.

Write N=13954/11025>5/4. As before, for L<=u<=2L,
C_(p_L)(u)<=N/2. On the shorter range 3L/2<=u<=2L,
both arguments in the overlap have absolute normalized
coordinate at least 1/2, and the overlap length is at most L/2.
The shape f is positive and decreasing in |s|, with
\[
f(1/2)=\frac{4861}{4480}<\frac{11}{10}.
\]
It follows directly that
\[
0\le C_{p_L}(u)\le\frac{121}{400}
<\frac5{16}<\frac N4
\quad(3L/2\le u\le2L).
\]
Thus N-C is at least N/2 on [L,2L], with an additional
N/4 on [3L/2,2L]. Outside these ranges it remains
nonnegative. Integrate the lower kernel polynomial (122):
\[
\begin{aligned}
\langle p_L,\mathcal R_Lp_L\rangle
&\ge N c_R(L)
+\frac N4\int_{3L/2}^{2L}
\left(\frac7{15360}u^3+\frac5{2048}u^4\right)\,du\\
&=N c_R(L)
+N\left(\frac{245L^4}{786432}
+\frac{781L^5}{262144}\right).
\end{aligned}
\]
All three coefficients increase with L. At L_0,
\[
c_R(L_0)=\frac{689087}{13107200000}
>\frac{21}{400000},
\]
since 689087>688128=21*32768. Also, using N>5/4,
\[
\begin{aligned}
N\frac{245L_0^4}{786432}
&>\frac{2941225}{503316480000}>\frac{29}{5000000},\\
N\frac{781L_0^5}{262144}
&>\frac{65631335}{3355443200000}>\frac{39}{2000000}.
\end{aligned}
\]
Consequently the entire remainder contribution is greater than
\[
\frac{21}{320000}+\frac{29}{5000000}+\frac{39}{2000000}
=\frac{3637}{40000000}>\frac9{100000}.
\]
This proves (142) by integrating the retained actual-kernel
lower bound against a proved overlap loss. It neither
changes the source nor assumes that its residual vanishes.

## BRT-T35. The complete returned scalar is positive through 44/125

The exact main-column readings (114) still satisfy the
derivative bounds from T30:
0<R_2'(L)<7/4 and -3/10<R_4'(L)<0.
Their proof uses L<3/8 and R_0(2L)<1/800, both valid here.
For h<=1/500, the same initial enclosures give
\[
-\frac{57}{100000}<R_2(L)<\frac{371}{100000},
\qquad
-\frac{53}{20000}<R_4(L)<0.
\]
Since sqrt(5)>11/5,
\[
|r_{L,2}^{(0)}|<\frac{17}{10000},\qquad
|r_{L,4}^{(0)}|<\frac9{10000}.
\]
Every logarithmic tail coefficient remains in the unchanged
sequence (103)--(104). Its norm is below 357/100000;
the additional sixth coefficient is still below 1/600000
because L^2<1/8. Therefore
\[
\begin{aligned}
\|r_L^{(0)}\|_{\mathcal D_L^{-1}}^2
&<\frac{20}{3}\left(\frac{17}{10000}\right)^2
+\frac{10}{3}\left(\frac9{10000}\right)^2
+2\left(\frac{2143}{600000}\right)^2\\
&=\frac{17092898}{360000000000}
<\left(\frac{69}{10000}\right)^2.
\end{aligned}
\tag{143}
\]
Restore the full Gamma remainder, prime column and pole
correction, using T32 and (141):
\[
\|r_L\|_{\mathcal D_L^{-1}}
<\frac{69}{10000}+\frac{11}{10000}
+\frac{13}{2000}+\frac1{6000}
=\frac{11}{750}.
\]
Inverse form order from (140) gives the complete penalty bound
\[
\boxed{
\langle r_L,C_L^{-1}r_L\rangle<\frac{121}{562500}.
}
\tag{144}
\]

### Trial value with the same enlarged caps

Let mathcal B_2(L) denote (124) with its constant prime
allowance 27/1000000 replaced by 51/1000000.
Equation (141) justifies this lower comparison on the whole
new band. Its derivative is still (125); the argument of
T31 still gives mathcal B_2'(L)>-1/50, since the endpoint
remainder is bounded above by 1/800.

At L_0, mathcal B_2 exceeds 1013/5000000-24/1000000.
Add (142) once, then integrate the derivative bound:
\[
\boxed{
q_L[p_L]>
\frac{1013}{5000000}-\frac{24}{1000000}
+\frac9{100000}-\frac h{50}
\ge\frac{1143}{5000000}.
}
\tag{145}
\]
Every signed prime and endpoint-remainder loss is retained.
Subtract (144) from (145) in the exact returned-scalar identity:
\[
\boxed{
s_L>\frac{1143}{5000000}-\frac{121}{562500}
=\frac{607}{45000000}>\frac1{75000}
\quad(L_0\le L\le L_2).
}
\tag{146}
\]

### Full-family bound

The resulting unrestricted source estimate is
\[
\boxed{
0<L\le\frac{44}{125},\quad g\in D_L
\quad\Longrightarrow\quad
q_L[g]\ge\frac1{130000}\|g\|^2.
}
\tag{147}
\]
There is no mean-zero or cosh-neutrality restriction.

Indeed (140) implies C_L>2I/31, and (144) gives
\[
\|C_L^{-1}r_L\|^2
<\frac{3751}{1125000}<\left(\frac2{25}\right)^2.
\]
The original return
rho_L=C_L^(-1)r_L-(p_L-phi_0^(L)) therefore still has norm
below 3/5. For g=c phi_0^(L)+h_0 and z=h_0+c rho_L,
the same exact square completion gives
\[
q_L[g]=s_Lc^2+\|C_L^{1/2}z\|^2,\qquad
\|g\|^2\le(43/25)c^2+2\|z\|^2.
\]
Now (43/25)/130000<1/75000 and 2/130000<2/31,
proving (147) on [L_0,L_2]. Support inclusion and the retained
source-form closure give every smaller width.

The stronger margins 1/50000 through 351/1000 and
1/25000 through 7001/20000 remain valid, as does the
earlier DWE margin on its smaller interval. The new wider
interval is obtained by accepting a smaller uniform margin,
not by claiming that the source margin increases with width.

**Scope and gate accounting.** The working written estimates
close positivity through 44/125. The full source and nonzero
residual have been assessed on that same band: integrated
prime caps, a sharper mean-zero scalar drift, and the actual
large-shift overlap supply the improvement. Positivity beyond
this finite interval and all-width positivity/RH remain open.
No new independent assumption is added. No independent
certification, scientific computation, numerical simulation,
certificate replay, subagent or workflow dispatch is performed.

## BRT-T36. The actual prime coefficient permits a wider boundary comparison

For T36--T39 use
\[
L_0=\frac7{20},\qquad L_3=\frac{1761}{5000},
\qquad 0\le h=L-L_0\le\frac{11}{5000}.
\]
Keep the admitted quartic p_L, its unit constant reading,
the unchanged source T_L, and r_L=Q_0^(L)T_Lp_L.
No new trial or domain assumption is introduced.

Only prime 2 is active: log2<2L_0 and 2L_3<1<log3.
The caps remain disjoint since L_3<log2. OPC-T2 retains
6931/10000<log2<1733/2500. Since d_L=2-(log2)/L
increases with L,
\[
\boxed{
\frac{19}{1000}<\frac{17}{875}<d_L
<\frac{113}{3522}<\frac1{31}.
}
\tag{148}
\]
The lower bound uses L=L_0 and log2<1733/2500; the upper
bound uses L=L_3 and log2>6931/10000.

### A sharper arithmetic coefficient absorbs these caps

The elementary inequalities
sqrt(2)>140/99 and log2<7/10 give
\[
a_2=\frac{\log2}{\sqrt2}<\frac{99}{200}.
\]
The square-root inequality follows from 19600<19602.
Also log2>27/40 and sqrt(2)<3/2 give a_2>9/20;
this lower bound will be used for a signed residual reading.

Use the increasing function W(t)=w(t)-t/(4-2t) from T29.
On a prime cap, |x|/L>30/31. At t_0=(30/31)^2,
\[
1-t_0=\frac{61}{961},\qquad
\frac{t_0}{4-2t_0}=\frac{225}{511}<\frac{441}{1000}.
\]
Since log(1+y)<y,
\[
\begin{aligned}
W(t_0)
&=\log2-\frac14\log(976/961)-\frac{225}{511}\\
&>\frac{6931}{10000}-\frac{15}{3844}
-\frac{441}{1000}
>\frac{99}{400}>\frac{a_2}{2}.
\end{aligned}
\]
Thus w_L-(a_2/2)1_(mathcal C_L)>=v_L on the new caps,
and W>=0 handles their complement. The actual prime has
been absorbed using its proved coefficient, without assuming
that W exceeds 1/4 on these larger caps.

All premises of T32's scalar derivative bound remain valid:
k(2L)<100/91 and L^2<1/8. Consequently C_L>=K_L-2hI.
Together with K_L>=D_L+I/200, this proves
\[
\boxed{
C_L\succeq\mathcal D_L+\frac3{5000}I
\quad(L_0\le L\le L_3).
}
\tag{149}
\]
This comparison holds on the full mean-zero form domain.

### Remainder and pole bounds on the same band

For u<=2L_3, u^2<1/2 and u<71/100 still hold. The full
kernel proof of (129), not merely its old stated interval,
therefore gives R_0(u)<1/800 and
\[
J_R(L)<\frac1{6250}+\frac{2h}{800}
\le\frac{331}{2000000}<\frac1{6000}.
\]
Hence the complete Gamma-remainder column has D_L-inverse
norm below 11/10000. The pole proof also retains
L^2<1/8, sqrt(2L)<6/7 and the same positive shape and
moments; its complete correction norm is below 1/6000.

The lower kernel estimate (122) remains valid since
cosh(u)<4/3. Every step of T34's large-shift overlap proof
then applies on this band: its polynomial coefficients increase
with L, while the shape itself is unchanged. Thus
\[
\langle p_L,\mathcal R_Lp_L\rangle>\frac9{100000}
\]
on the enlarged band as well. No remainder is counted twice.

## BRT-T37. Signed prime readings are combined before taking the norm

Put d=d_L, A=2/35 and B=11/7. On 0<=t<=1/31,
the actual cap polynomial satisfies both bounds
\[
A+\frac{10}{7}t\le F_c(t)\le A+Bt.
\tag{150}
\]
For the lower bound, the remaining quadratic factor is
39/14-(7/2)t+(7/8)t^2>0 on this interval.
For the upper bound, use
10/7+39/(14*31)+7/(8*31^3)<11/7.

Write pi_n=<Pp_L,phi_n^(L)> for an even degree n.
The exact signed prime reading is
\[
\pi_n=-\frac{a_2\sqrt{2n+1}}2
\int_0^d F_c(t)P_n(1-d+t)\,dt.
\tag{151}
\]
The orientation and sign here come from the actual negative
translation column. Define V(d)=Ad+(5/7)d^2. Equation (150)
gives integral_0^d F_c>=V(d), and (148) gives
\[
V(d)>V(19/1000)=\frac{1881}{1400000}
>\frac{13}{10000}.
\]
For upper estimates,
\[
I_1(d):=Ad+\frac B2d^2
\le\frac2{1085}+\frac{11}{13454}<\frac{27}{10000}.
\]
The two endpoint fractions are respectively less than
185/100000 and 82/100000.

### The relevant low prime readings have proved negative signs

The polynomial identities
P_2'=3P_1, P_4'=7P_3+3P_1 and
P_6'=11P_5+7P_3+3P_1, together with |P_j|<=1,
give, for 30/31<=s<=1,
\[
P_2(s)\ge\frac{28}{31},\qquad
P_4(s)\ge\frac{21}{31},\qquad
P_6(s)\ge\frac{10}{31}.
\tag{152}
\]
These follow by integrating each derivative bound back from
P_n(1)=1. In particular pi_2, pi_4 and pi_6 are negative.
No common sign is asserted for all higher prime readings.

For degree two, (151) and sqrt(5)<9/4 give
\[
-\frac{243}{160000}<\pi_2<0.
\]
The exact main-column formulas (114) and their derivative
bounds from T30 give on this wider band
\[
-\frac{57}{100000}<R_2(L)<\frac{406}{100000},
\qquad
-\frac{271}{100000}<R_4(L)<0.
\]
Consequently, using sqrt(5)>11/5,
\[
-\frac{13}{50000}<r_{L,2}^{(0)}<\frac{37}{20000},
\qquad
-\frac{91}{100000}<r_{L,4}^{(0)}<0.
\]
Combining the signed second readings, rather than adding
their separate absolute upper bounds, yields
\[
|r_{L,2}^{(0)}+\pi_2|<\frac{37}{20000}.
\tag{153}
\]
Indeed 13/50000+243/160000=1423/800000
is less than 37/20000.

### A cap slope improves the fourth prime reading

On [30/31,1], P_4''>0 and
P_4'(30/31)=512550/59582>8. Therefore
\[
0<P_4(1-y)\le1-8y\quad(0<=y<=1/31).
\]
Inserting this and (150) in (151) bounds its integral by
\[
J(d)=Ad+\left(\frac B2-4A\right)d^2-\frac{4B}{3}d^3.
\]
This polynomial is increasing on [0,1/31], since
J'(d)=A+(B-8A)d-4Bd^2>=A-4B/31^2>0.
At the endpoint,
\[
J(1/31)
=I_1(1/31)-\frac8{33635}-\frac{44}{625611}
<\frac{27}{10000}-\frac{23+7}{100000}
=\frac{12}{5000}.
\]
Thus |pi_4|<(3/4)J(d)<9/5000. It follows that
\[
|r_{L,4}^{(0)}+\pi_4|
<\frac{271}{100000}<\frac{11}{4000}.
\tag{154}
\]

### The sixth logarithmic reading is partly cancelled by the actual prime

The main-column sixth reading is positive and explicit:
\[
r_{L,6}^{(0)}
=t_6+\frac{L^2\sqrt{13}}{308880},\qquad
t_6=\frac{2\sqrt{13}}{3465}.
\]
Using 18/5<sqrt(13)<181/50 and L^2<1/8,
\[
\frac1{500}<r_{L,6}^{(0)}
<\frac{209}{100000}+\frac1{600000}
<\frac{21}{10000}.
\]
Meanwhile a_2>9/20, sqrt(13)>18/5, (152) and V(d)>13/10000
give a strictly positive lower bound on its opposite prime
reading:
\[
|\pi_6|>
\frac{81}{100}\frac{10}{31}\frac{13}{10000}
>\frac3{10000}.
\]
The upper bound from (151) is |pi_6|<27/10000.
Therefore
\[
-\frac7{10000}<r_{L,6}^{(0)}+\pi_6<\frac9{5000},
\qquad
\boxed{|r_{L,6}^{(0)}+\pi_6|<\frac9{5000}.}
\tag{155}
\]
This is an actual signed arithmetic cancellation in this
specified component. The complete residual has not vanished;
the Gamma remainder and pole correction are restored in T38.

## BRT-T38. Projection accounting retains and controls the entire higher tail

Let R_8 project onto the even Legendre degrees n>=8.
The main-column polynomial terms have degree at most six.
The infinite logarithmic estimate in (104) thus gives directly
\[
\|R_8r_L^{(0)}\|^2<\frac1{120050},
\qquad
\|R_8r_L^{(0)}\|<\frac{29}{10000}.
\tag{156}
\]

For the prime tail, first subtract the already retained
degrees 0,2,4,6 from its exact L2 norm. Define
\[
\kappa_0=1,\quad\kappa_2=\frac{28}{31},\quad
\kappa_4=\frac{21}{31},\quad\kappa_6=\frac{10}{31},
\qquad
K=\sum_{n=0,2,4,6}(2n+1)\kappa_n^2
=\frac{10150}{961}.
\]
Thus 10<K<11. Write b=5/7 and
\[
U_2(d)=A^2d+ABd^2+\frac{B^2}{3}d^3,\qquad
U(d)=U_2(d)-K(Ad+bd^2)^2.
\]
Disjointness of the caps, the exact readings (151),
and their positive polynomial lower bounds (152) give
\[
\begin{aligned}
\|R_8Pp_L\|^2
&=\frac{a_2^2}{4}
\left[
\int_0^dF_c(t)^2dt
-\sum_{n=0,2,4,6}(2n+1)
\left(\int_0^dF_c(t)P_n(1-d+t)dt\right)^2
\right]\\
&\le\frac{a_2^2}{4}U(d).
\end{aligned}
\tag{157}
\]
This is norm accounting for the infinite complement, not
a truncation of it or an invariance assumption for T_L.

### The entire upper comparison is increasing on the cap range

Differentiate the declared polynomial U:
\[
U'(d)=A^2+2(AB-KA^2)d+(B^2-6KAb)d^2-4Kb^2d^3.
\]
Here AB-KA^2>0. Since K<11, the other coefficients give
\[
U'(d)\ge
A^2-\frac{11}{49}d^2-\frac{1100}{49}d^3
\ge A^2-\frac{1441}{1459759}>0.
\]
The last step uses A^2>3/1000 and
1441/1459759<1/1000. Thus U(0)=0 and U is positive
and increasing for 0<d<=1/31.

At that endpoint,
\[
U_2(1/31)
=\frac4{37975}+\frac{22}{235445}
+\frac{121}{4379277}
<\frac{106+94+28}{1000000}.
\]
Also V(1/31)=2/1085+5/6727>1/400: its terms exceed
18/10000 and 7/10000. Since K>10,
\[
U(d)\le U(1/31)
<\frac{228}{1000000}-\frac1{16000}
=\frac{331}{2000000}.
\]
Using a_2<1/2 in (157) therefore yields
\[
\boxed{
\|R_8Pp_L\|^2<\frac{331}{32000000}
<\left(\frac{33}{10000}\right)^2.
}
\tag{158}
\]
All higher prime modes survive in this upper bound.
No sign is required for their individual cross pairings.
Equations (156)--(158) imply
\[
\|R_8(r_L^{(0)}+Pp_L)\|<\frac{31}{5000}.
\]

### Combine orthogonal blocks, then restore the remaining source columns

Set v_L=r_L^(0)+Q_0^(L)Pp_L. Using (153)--(155) and
the complete higher-tail estimate,
\[
\begin{aligned}
\|v_L\|_{\mathcal D_L^{-1}}^2
&<\frac{20}{3}\left(\frac{37}{20000}\right)^2
+\frac{10}{3}\left(\frac{11}{4000}\right)^2
+2\left(\frac9{5000}\right)^2
+2\left(\frac{31}{5000}\right)^2\\
&=\frac{78831}{600000000}
<\left(\frac{23}{2000}\right)^2.
\end{aligned}
\tag{159}
\]
Only orthogonal projections of the comparison norm are used.
The actual operator is neither diagonalized nor assumed to
commute with these projections.

Restore the full Gamma and pole columns from T36:
\[
r_L=v_L+Q_0^{(L)}\mathcal R_Lp_L+e_{{\rm pole},L},
\qquad
\|r_L\|_{\mathcal D_L^{-1}}
<\frac{23}{2000}+\frac{11}{10000}+\frac1{6000}
=\frac{383}{30000}.
\]
By (149) and inverse form order, the complete penalty is
\[
\boxed{
\langle r_L,C_L^{-1}r_L\rangle
<\frac{146689}{900000000}<\frac1{6000}.
}
\tag{160}
\]
This improves (144)'s upper penalty on a larger interval.
The improvement uses proved signed readings and projection
accounting; it does not assert that conjugation, rotation
or a sign convention establishes positivity.

## BRT-T39. Both the positive interval and its frontier margin improve

The actual prime-form allowance must also be enlarged.
The same affine cap bound, now valid for d<1/31, gives
\[
\begin{aligned}
a_2C_{p_L}(\log2)
&<\frac14\left(\frac{A^2}{31}
+\frac{AB}{31^2}+\frac{B^2}{6\cdot31^3}\right)\\
&=\frac14\left(\frac4{37975}+\frac{22}{235445}
+\frac{121}{8758554}\right)
<\frac{54}{1000000}.
\end{aligned}
\tag{161}
\]
The three fractions are less than 106, 94 and 14,
respectively, in units of 1/1000000.

Let mathcal B_3(L) be (124) with its constant prime
allowance replaced by 54/1000000. Equation (161) justifies
it throughout the new band. Its derivative is still (125),
and T31's proof still gives mathcal B_3'(L)>-1/50,
because R_0(2L)<1/800. At L_0 its value exceeds
1013/5000000-27/1000000.

Add the actual remainder lower bound 9/100000 from T36
once and integrate the derivative estimate:
\[
\boxed{
q_L[p_L]>
\frac{1013}{5000000}-\frac{27}{1000000}
+\frac9{100000}-\frac h{50}
\ge\frac{277}{1250000}.
}
\tag{162}
\]
The complete inverse subtraction now yields
\[
\boxed{
s_L>
\frac{277}{1250000}-\frac{146689}{900000000}
=\frac{52751}{900000000}>\frac7{120000}
\quad(L_0\le L\le L_3).
}
\tag{163}
\]
For the final strict comparison, 7/120000=52500/900000000.

The full-family conclusion is
\[
\boxed{
0<L\le\frac{1761}{5000},\quad g\in D_L
\quad\Longrightarrow\quad
q_L[g]\ge\frac1{30000}\|g\|^2.
}
\tag{164}
\]
No mean-zero or cosh-neutrality restriction is imposed.

Indeed (149) gives C_L>2I/31, and the last bound in (160)
implies
\[
\|C_L^{-1}r_L\|^2<\frac{31}{12000}
<\left(\frac2{25}\right)^2.
\]
As before, ||p_L-phi_0^(L)||<13/25, so the original return
rho_L has norm below 3/5. For g=c phi_0^(L)+h_0 and
z=h_0+c rho_L, the exact source square and norm bound are
\[
q_L[g]=s_Lc^2+\|C_L^{1/2}z\|^2,\qquad
\|g\|^2\le(43/25)c^2+2\|z\|^2.
\]
The comparisons (43/25)/30000<7/120000 and
2/30000<2/31 prove (164) on the band. Support inclusion
and retained source-form closure give all smaller widths.

The new margin 1/30000 supersedes the former 1/130000
through 44/125 and 1/50000 through 351/1000.
The stronger 1/25000 through 7001/20000 remains available,
as does DWE's 1/1700 on its smaller interval.

**Scope and gate accounting.** The working written proof
now closes positivity through 1761/5000, with a stronger
uniform margin at the former frontier. Every infinite tail
and every actual source column is retained. Only the specified
low-mode cancellations were used; no sign was assigned to
all higher cross pairings. Positivity beyond this finite interval
and all-width positivity/RH remain unresolved.
No new independent assumption is added. No independent
certification, scientific computation, numerical simulation,
certificate replay, subagent or workflow dispatch is performed.


## BRT-T40. Multiplicative readings preserve ratios, not missing source signs

For a real zero-extended profile g, let T_u denote full-line
translation, T_ug(x)=g(x-u). For any finite list of positive
labels m_1,...,m_k define its observation Gram matrix by
\[
\boxed{
G_g(m_i,m_j)
=\langle T_{\log m_i}g,T_{\log m_j}g\rangle
=C_g(\log(m_i/m_j)),\qquad
G_g(Mm_i,Mm_j)=G_g(m_i,m_j)\quad(M>0).
}
\tag{165}
\]
It is positive semidefinite, since for real coefficients c_i,
\[
\sum_{i,j}c_ic_jG_g(m_i,m_j)
=\left\|\sum_i c_iT_{\log m_i}g\right\|^2\ge0.
\]
The common multiplier gives one common unitary translation.
This proves (165) without assuming positivity of the signed
arithmetic source. For example, the pair M,2M reads C_g(log 2);
a single large prime does not by itself contain that reading.
Common rescaling adds no independent correlations. Adding
new ratios can add constraints between existing readings.

The useful support pattern below has relative positions
\[
\boxed{
1,\quad \frac32,\quad 2,\quad 3;
\qquad
2,\quad 3,\quad 4,\quad 6
\text{ after multiplication by }2.
}
\tag{166}
\]
These are observation labels, not four active arithmetic
source terms. In the window used below only prime powers
2 and 3 contribute. The labels 4 and 6 do not insert new
source penalties, and the ratio 3/2 is not declared a prime.

This common-scale invariance is the elementary translation
Gram identity. It is not a newly proved positivity theorem.
It also does not identify translation of a full-line reading
bank with changing the physical support I or enlarging L.
The quantitative gain below uses the actual truncated prime
operators and their shared support, not Gram positivity alone.

## BRT-T41. The two active primes have an exact joint norm

Write
\[
u=\log2,\quad v=\log3,\quad d=v-u,\qquad
a=\frac{u}{\sqrt2},\quad b=\frac{v}{\sqrt3}.
\]
For this section suppose v/2<L<=u. Precisely prime powers
2 and 3 are active; at L=u the shift log 4=2L still has
zero overlap. Use PCM's positive-sign convention
\[
\begin{aligned}
(\mathcal P_2g)(x)&=\frac a2\bigl(g(x-u)+g(x+u)\bigr),\\
(\mathcal P_3g)(x)&=\frac b2\bigl(g(x-v)+g(x+v)\bigr),\\
\mathcal P_{23}&=\mathcal P_2+\mathcal P_3,\qquad
\langle g,\mathcal P_{23}g\rangle=aC_g(u)+bC_g(v).
\end{aligned}
\tag{167}
\]
All outside values are zero. The actual source subtracts
this operator. The sign convention does not assert that
the operator itself is positive semidefinite.

Put epsilon=2L-v>0 and A=(-L,L-v). Since 8<9,
\[
0<\epsilon\le2u-v<d,\qquad 3u<2v.
\]
The four translates A, A+d, A+u, A+v are disjoint up to
endpoints. Indeed their successive offsets differ by
d, 2u-v, d, all at least epsilon. For each x in A the
ordered path and its operator matrix are
\[
\boxed{
(x+u,\ x,\ x+v,\ x+d),\qquad
M_{23}=\frac12
\begin{pmatrix}
0&a&0&0\\
a&0&b&0\\
0&b&0&a\\
0&0&a&0
\end{pmatrix}.
}
\tag{168}
\]
The first and last edges have length u; the middle edge
has length v.

Here the path accounts for every interaction with these
four intervals. A v-edge must have its left endpoint
in A, so both endpoints are already included. Each
vertex can have at most one u-neighbor because 2u>=2L.
The two displayed u-edges give that neighbor for all
four vertices. Thus no further u-edge leaves the union.
Nor can a v-edge connect an inner vertex to a new point:
that would require 2L>2v-u, whereas 2L<=2u<2v-u.
This last observation also follows from the disjoint
description of all v-edge endpoints A and A+v.
On the complement there are no v-edges, only isolated
u-pairs and isolated vertices.

Translation preserves Lebesgue measure. Consequently
the operator is an orthogonal direct integral of the
four-by-four matrices (168), together with two-by-two
u-pair matrices and zero components. No smoothness or
pointwise sign of g is required for this L2 decomposition.

Path reversal splits (168) into symmetric and
antisymmetric two-dimensional sectors. On a symmetric
vector (c,e,e,c), the eigenvalue equations are
lambda c=ae/2 and lambda e=ac/2+be/2.
The four eigenvalues are
\[
\pm\frac{\sqrt{b^2+4a^2}+b}{4},
\qquad
\pm\frac{\sqrt{b^2+4a^2}-b}{4}.
\]
The largest absolute value exceeds a/2, the norm of an
isolated u-pair. Hence
\[
\boxed{
\|\mathcal P_{23}\|
=\Lambda_{23}:=\frac{b+\sqrt{b^2+4a^2}}4,\qquad
|aC_g(u)+bC_g(v)|\le\Lambda_{23}\|g\|^2.
}
\tag{169}
\]
The exact norm remains valid on the real even subspace.
For sharpness, choose a nonzero smooth compact function
chi on A symmetric under x -> -x-v. Multiply a positive
symmetric top eigenvector (c,e,e,c) of (168) by chi(x)
on the four path intervals, and set g=0 elsewhere.
Reflection reverses the path and preserves chi, so g is
even. It is a smooth compact test in I and satisfies
P_23 g=Lambda_23 g. Thus sharpness uses the established
source core, not an unproved domain for a new ansatz.

The separate disjoint-cap estimate would give (a+b)/2.
The exact improvement is
\[
\boxed{
\Delta(a,b)=\frac{a+b}{2}-\Lambda_{23}
=\frac{2a+b-\sqrt{b^2+4a^2}}4>\frac1{10}.
}
\tag{170}
\]
For the explicit rational threshold in (170), TPW-T3/T5
give log 2>27/40 and log 3>263/240>21/20.
Since sqrt 2<3/2 and sqrt 3<7/4,
a>9/20 and b>3/5. Each partial derivative of Delta is
strictly positive for a,b>0: they are
(2-4a/sqrt(b^2+4a^2))/4 and
(1-b/sqrt(b^2+4a^2))/4. Therefore
\[
\Delta(a,b)>
\Delta(9/20,3/5)
=\frac{3/2-\sqrt{117/100}}4>\frac1{10},
\]
where the last inequality follows from 117/100<121/100.
This compares the norm of the combined signed reading
with the sum of separate operator norm allowances.

We will also use an explicit rational upper bound.
The retained estimates a<1/2 and
log 3<11/10, sqrt 3>5/3 give b<33/50<2/3.
Monotonicity of the expression in (169) gives
\[
\boxed{
\Lambda_{23}<
\frac{2+\sqrt{13}}{12}<\frac{15}{32};
\qquad
\|t_2\mathcal P_2+t_3\mathcal P_3\|
\le\Lambda_{23}\quad(0\le t_2,t_3\le1).
}
\tag{171}
\]
Here sqrt 13<29/8. The coupling bound follows from the
same path matrices with a,b replaced by t_2a,t_3b;
their norm formula is nondecreasing in these nonnegative
weights. QED.

The strict lower endpoint L>v/2 matters. At L=v/2 the
v-caps are empty and the norm is a/2. Arbitrarily narrow
nonempty caps support the localized tests above, so no
continuity in operator norm across this support threshold
has been assumed. For L>u, additional chains and prime
powers occur; (169) is not extrapolated to that range.

## BRT-T42. The joint bound makes a specified polynomial tail positive

Now fix L_*=11/20. TPW's rational inequalities place
v/2<L_*<u. On the unchanged form domain D_(L_*), define
the diagnostic coupling family
\[
q_t[g]=\ell(g)^2+\frac{\Omega_0}{2}\|g\|^2
+\int_0^\infty k(z)(\|g\|^2-C_g(z))\,dz
-t_2aC_g(u)-t_3bC_g(v),\quad t\in[0,1]^2.
\tag{172}
\]
Its actual arithmetic endpoint is t=(1,1).
The changes in t are bounded self-adjoint perturbations;
the form and operator domains are common.

Use the degree labels of this note, rather than ASM's
index counting the even modes:
\[
\boxed{
V_5=\operatorname{span}\{\phi_0,\phi_2,\phi_4,\phi_6,\phi_8\},
\quad
\phi_{2j}(x)=\sqrt{\frac{4j+1}{2L_*}}P_{2j}(x/L_*),
\quad Q_5=I-P_5.
}
\tag{173}
\]
For h in Q_5 D_(L_*), its mean m_0 vanishes. GRR-T1/T2,
with its fixed comparison parameter 1/2, gives
\[
\begin{aligned}
q_t[h]\ge{}&
B_*(L_*)\|h\|^2+J_{L_*}[h]
+\int_Iw_{L_*}(x)|h(x)|^2\,dx+\ell(h)^2\\
&-\langle h,(t_2\mathcal P_2+t_3\mathcal P_3)h\rangle,\\
B_*(L)&=\frac{\Omega_0}{2}+F(2L)+\frac12\log2
+\frac{47L}{96}-\frac{5L^3}{36},\\
w_L(x)&=-\frac14\log(1-(x/L)^2)\ge0 .
\end{aligned}
\tag{174}
\]
Indeed GRR's polynomial is 47/192-5z^2/96.
Its two finite moment corrections are
-47m_0^2/384+5m_0m_2/96, which vanish here.
The remaining Gamma integral is nonnegative on the same
closed domain. The pole square and boundary multiplier
are displayed before any lower estimate discards them.

TPW-T3 supplies Omega_0/2>-2689/1000,
F(11/10)>59/50 and log 2>693/1000. At L_*,
\[
\frac{47L_*}{96}-\frac{5L_*^3}{36}
=\frac{14179}{57600}>\frac{49}{200}.
\]
Consequently
\[
\boxed{
B_*(L_*)>
-\frac{2689}{1000}+\frac{59}{50}
+\frac{693}{2000}+\frac{49}{200}
=-\frac{367}{400}.
}
\tag{175}
\]
ASM-T1's scale-invariant Legendre energy bound gives
J_(L_*)[h]>=H_10||h||^2/2 on this complement.
Since H_10=7381/2520, (171), (174) and (175) imply
\[
\boxed{
q_t[h]\ge\frac{3943}{50400}\|h\|^2
>\frac3{40}\|h\|^2
\quad(h\ne0,\ h\in Q_5D_{L_*},\ t\in[0,1]^2).
}
\tag{176}
\]
The first displayed coefficient follows from
7381/5040-367/400-15/32=3943/50400;
the proof in fact gives a strict first inequality for
nonzero h. The corresponding closed tail operators obey
C_5(t)>=3I/40 uniformly on the coupling square.

For a concrete comparison, using the same Gamma and
Legendre bounds but the separate allowances a<1/2,
b<2/3 would give only
\[
\boxed{
\frac{H_{10}}2-\frac{367}{400}-\frac7{12}
=-\frac{229}{6300}<0.
}
\tag{177}
\]
Replacing that allowance 7/12 by the proved joint
allowance 15/32 improves this particular lower bound
by 11/96 and makes it positive. Equation (177) records
failure of that specified estimate, not negativity of
the true tail. QED.

This is a quantitative improvement for a specified
polynomial cut. ASM already proved that some finite
cut has a positive tail at every fixed width. TPW
already proved positivity through L_* on its family
with the two neutral readings m_0=ell=0. Neither
earlier result is being claimed anew. Five is not a
minimal cut dimension or a universal observer count.
Its advantage here is that the five retained columns
are explicit admitted polynomial columns. The full
family still retains all five coordinates and their
coupling to the infinite tail.

## BRT-T43. The improved tail bound retains the complete returned memory

At L_* let U:R^5 -> V_5 synthesize the normalized modes
in (173), and let T(t) be the source operator of (172).
Define the actual source blocks
\[
\begin{gathered}
A(t)=U^*T(t)U,\qquad B(t)=Q_5T(t)U,\qquad
C(t)=C_5(t),\\
E(t)=U-C(t)^{-1}B(t),\qquad
\mathsf S(t)=A(t)-B(t)^*C(t)^{-1}B(t),\\
\boxed{0<C(t)^{-1}\preceq\frac{40}{3}I.}
\end{gathered}
\tag{178}
\]
ASM-T2 admits every column of U to the operator domain.
The restricted inverse maps into the tail operator
domain. Since the coupling columns are bounded, the
usual block domain argument also admits E(t) to the
full source domain: its pairing against a tail test
is represented in L2, and its retained pairing is
finite dimensional. Thus these source expressions do
not require a postulated logarithmic trial domain.

For any five-column operator-domain trial Z satisfying
P_5Z=U, put R=Q_5T(t)Z. Completing the same tail square
as ASM-T3/T5 gives
\[
\boxed{
\begin{aligned}
\mathsf S(t)&=Z^*T(t)Z-R^*C(t)^{-1}R,\\
Z^*T(t)Z-\frac{40}{3}R^*R
&\preceq\mathsf S(t)\preceq Z^*T(t)Z.
\end{aligned}
}
\tag{179}
\]
For completeness, write Z=U+Y with Q_5Y=Y.
Then R=B+CY. Expansion of
R^*C^{-1}R cancels exactly the mixed terms and Y^*CY
in Z^*TZ, leaving A-B^*C^{-1}B. Inequality (178)
proves the enclosure. No component of R is omitted.

The same bound controls a simultaneous insertion of
both primes. Write bold 0=(0,0), bold 1=(1,1),
V=P_23, E_0=E(bold 0), and D=Q_5VE_0. PCM-T4 gives
\[
\boxed{
\mathsf S(\mathbf1)=\mathsf S(\mathbf0)
-E_0^*\mathcal P_{23}E_0
-D^*C(\mathbf1)^{-1}D,
\quad
0\preceq D^*C(\mathbf1)^{-1}D
\preceq\frac{40}{3}D^*D.
}
\tag{180}
\]
One can also bound the latter matrix directly by
\[
\boxed{
\frac{40}{3}D^*D
\preceq\frac{40}{3}\Lambda_{23}^2 E_0^*E_0
\preceq\frac{375}{128}E_0^*E_0.
}
\tag{181}
\]
Here projection is contractive, (169) bounds V, and
(40/3)(15/32)^2=375/128. This last coarse enclosure is
optional; (180) retains the actual joint column D.
If D_i=Q_5P_iE_0 for i=2,3, then
D^*D=D_2^*D_2+D_3^*D_3+D_2^*D_3+D_3^*D_2.
The cross terms have not been silently replaced by
independent scalar costs. The direct term in (180)
need not have either fixed sign.

These are the existing exact source-memory identities
with a new explicit tail constant, not a new positivity
criterion inferred from coordinate invariance. They
show where the joint prime improvement enters the same
returned-source lower bound. Neither the entries of
S(bold 1) nor a positive lower matrix in (179) or (180)
have been established here. In particular (181) alone
does not establish their sign.

**Gate accounting.** The exact two-prime operator norm,
its improvement over separate cap costs, and positivity
of this specified infinite polynomial tail at L_*
are closed in the working written derivation.
The five-dimensional returned source sign at L_*
and all-width positivity/RH remain unresolved.
The unrestricted interval in (164) is unchanged.
No new independent assumption is added.


## BRT-T44. Group source translations by prime-power chains and shared ratios

Use the full source with every active prime power as in GRR-1.
For a prime p write tau=log p and let W=2L. Its complete
primitive-prime contribution, before subtraction from the
source, is
\[
\begin{aligned}
\mathcal P^{(p)}_L
&=\frac12\sum_{k\tau<W}\frac{\log p}{p^{k/2}}
\bigl(T_{k\tau}+T_{-k\tau}\bigr)\big|_I,\\
N_p&=\left\lceil\frac W\tau\right\rceil,\qquad
(M^{(p)}_n)_{ij}
=\begin{cases}
0,&i=j,\\
\dfrac{\log p}{2p^{|i-j|/2}},&i\ne j,
\end{cases}\quad 0\le i,j<n .
\end{aligned}
\tag{182}
\]
The translations use zero outside I. If N_p=1 the matrix
is zero and the prime is inactive. The strict support
cutoff is respected at W/tau an integer.

There is an exact chain decomposition and norm formula:
\[
\boxed{
\mathcal P^{(p)}_L
\simeq\int^\oplus M^{(p)}_{n(r)}\,dr,\qquad
\|\mathcal P^{(p)}_L\|=\|M^{(p)}_{N_p}\|.
}
\tag{183}
\]
To prove it, shift the interval to (0,W), and write each
point uniquely as r+j tau with 0<r<tau, j a nonnegative
integer, omitting empty fibres and endpoints of measure
zero. The fibre contains precisely the n(r) points still
in (0,W). A p^k translation connects indices differing
by k, with the coefficient in (182). Translation
preserves measure, so the L2 norm and quadratic form
split over these fibres. The largest occurring length
is ceil(W/tau), and occurs for a set of r of positive
measure, including when W/tau is an integer.

The matrices of smaller size are principal submatrices.
For each symmetric matrix with nonnegative entries,
|z^TMz|<=|z|^TM|z|. Thus its operator norm is its largest
nonnegative Rayleigh value, and these norms increase
under principal inclusion. This proves (183).
For N_p>=2 the top eigenvector can be chosen strictly
positive: take the absolute value of a maximizing
vector, and use the positive off-diagonal entries in
its eigenvalue equation. Multiply that eigenvector by a nonnegative smooth
cutoff supported in a subinterval of r on which the
chain length is N_p. This gives a smooth compact
source test attaining the norm. Reflection commutes
with the operator; adding its reflected test yields
a nonzero real even test with the same eigenvalue.
The domain statement therefore uses the existing
smooth core. The zero case N_p=1 is immediate. QED.

This is a useful classification of the source
translations at a declared width: powers of one
primitive prime share a chain. Distinct prime families
must then be joined where their observation points
coincide. A product or quotient of two bank labels
is a reading ratio; it is not automatically another
source coefficient.

For instance, the seven labels
\[
\boxed{
\mathcal B=\{6,8,9,12,16,18,24\},\qquad
\{6,12,24\},\quad\{6,18,9\},\quad\{24,8,16\}
}
\tag{184}
\]
form a triangle with ratios 2,2,4 and two attached
two-edge paths with ratios 3,2. The subgroups share
the labels 6 and 24. Multiplying every label by any
M>0 preserves all ratios by (165). M may itself be
a large prime; in that case all seven multiples are
composite observation labels and add no source terms.

The classification is by ratios, shared vertices
and allowed support, not a new partition of all
prime numbers into intrinsic classes. The next
active source after 2 and 3 is 4=2^2, whose weight
is log 2/2, not log 4/2. The explicit joint result
below retains that weight and the exact support.
No claim that all multi-prime support components
are finite at arbitrary widths follows from (183).

## BRT-T45. The 2,3,4 source block has an exact seven-point reduction

Keep u=log 2, v=log 3 and put
\[
\begin{gathered}
u<L\le v-\frac u2,\quad W=2L,\quad
\epsilon=W-2u,\quad w=2u-v,\quad h=2v-3u,\quad d=v-u=w+h,\\
a=\frac{u}{\sqrt2},\quad b=\frac{v}{\sqrt3},\quad
c=\frac u2,\qquad
\mathcal P_{234}=\mathcal P_2+\mathcal P_3+\mathcal P_4,\\
(\mathcal P_4g)(x)=\frac c2\bigl(g(x-2u)+g(x+2u)\bigr).
\end{gathered}
\tag{185}
\]
Here the subscript 2,3,4 labels an individual active
prime power. The full primitive-prime-2 operator
in (182) is P_2+P_4 on this range. Precisely 2,3,4
are active, since 4<e^W<=9/2<5.

The elementary inequalities 8<9 and 27<32 give
0<h<w, and 0<epsilon<=h. In the shifted interval
(0,W), the following disjoint components exhaust
all source translations, up to endpoints:

- For x in (0,epsilon), a seven-point component at
  x plus the offsets 0,u,2u,v,d,w,u+w.
- For x in (epsilon,w), a four-point component at
  x plus the offsets u,0,v,d, with the path (168).
- For x in (w+epsilon,d), an isolated u-pair
  x,x+u. This last family is empty if epsilon=h.

The seven-point offsets, in increasing order,
are 0,w,d,u,u+w,v,2u. Their multiplicative labels
are 1,4/3,3/2,2,8/3,3,4, which become (184)
after multiplication by 6. In the path-and-triangle
order 0,u,2u,v,d,w,u+w, their matrix is
\[
\boxed{
M_7=\frac12
\begin{pmatrix}
0&a&c&b&0&0&0\\
a&0&a&0&0&0&0\\
c&a&0&0&0&b&0\\
b&0&0&0&a&0&0\\
0&0&0&a&0&0&0\\
0&0&b&0&0&0&a\\
0&0&0&0&0&a&0
\end{pmatrix}.
}
\tag{186}
\]

Here is a direct support check. All 2u-edges have
left endpoint x in (0,epsilon), and are the long
edges of the displayed triangles. All v-edges have
left endpoint in (0,w+epsilon). Partitioning this
interval into (0,epsilon), (epsilon,w),
(w,w+epsilon) places each such edge in exactly one
seven-point or four-point component. At a seven-point
component, the u-edges are exactly
(0,u), (u,2u), (v,d), (w,w+u), with x added.
Any additional v-edge would require x>h or
x<epsilon-h, outside its base interval apart from
null endpoints. Any additional u-edge at an arm
would require x>w or x<epsilon-w, also impossible.
The only 2u-edge is the displayed one.

For x in (epsilon,w), the four-point component
has u-edges (0,u), (d,v) and v-edge (0,v).
The strict inequalities x>epsilon and x<w rule
out a second u-step at either end; the inequalities
epsilon<=h and w<d rule out additional v-edges.
No 2u-edge touches it. Finally x in (w+epsilon,d)
and x+u have only their mutual u-edge.

For completeness the intervals fill (0,W):
in increasing order the seven-point base intervals
start at 0,w,d,u,u+w,v,2u; the four-point intervals
fill (epsilon,w), (d+epsilon,u), (u+epsilon,u+w),
(v+epsilon,2u); the pair intervals fill
(w+epsilon,d) and (u+w+epsilon,v).
All omitted or shared endpoints have measure zero.
Thus no interaction crosses between the components.

The exact norm is
\[
\boxed{
\|\mathcal P_{234}\|=\lambda_7,\qquad
\lambda_7-\frac c2-\frac{a^2}{2\lambda_7}
-\frac{b^2\lambda_7}{4\lambda_7^2-a^2}=0,
\quad \lambda_7>\frac a2.
}
\tag{187}
\]
There is exactly one root in that range: the left
side tends to minus infinity at a/2 from above,
tends to plus infinity at infinity, and has derivative
\[
1+\frac{a^2}{2\lambda^2}
+\frac{b^2(4\lambda^2+a^2)}{(4\lambda^2-a^2)^2}>0.
\]
To identify it with the norm without a numerical
eigenvalue calculation, assign a reflection-symmetric
positive vector on (186) with values
(p,q,p,r,s,r,s), where
\[
p=1,\quad q=\frac a\lambda,\quad
r=\frac{2b\lambda}{4\lambda^2-a^2},\quad
s=\frac{ar}{2\lambda}.
\]
Its eigenvalue equations reduce exactly to (187).
For any positive vector z with Mz=lambda z and any
real y, applying
2|y_i y_j|<=(z_j/z_i)y_i^2+(z_i/z_j)y_j^2
to every edge gives |y^TMy|<=lambda||y||^2.
Hence this eigenvalue is the full norm of M_7.
The four-point matrix (168) and isolated pair are
principal submatrices, so their norms are no larger.

The norm is attained on real even smooth compact
tests. Choose a smooth function chi on (0,epsilon)
symmetric under x -> epsilon-x and multiply the
positive vector above by chi along the seven
translated intervals. In original coordinates
reflection reverses the two triangle endpoints
and exchanges the two arms; chi and the vector
are preserved. This supplies a nonzero even
source-core test. It proves (187) on the same
space used for the source, not only on unrestricted
point values.

A rational bound suffices for the next source estimate:
\[
\boxed{
\lambda_7<\frac35,\qquad
\left\|\sum_{n=2,3,4}t_n\mathcal P_n\right\|
\le\lambda_7<\frac35\quad(0\le t_n\le1).
}
\tag{188}
\]
Indeed a<1/2, c<7/20, and b<16/25.
The last estimate follows from log 3<11/10 and
sqrt 3>55/32, whose square is smaller than 3.
Use the positive vector
\[
z=(1,5/6,1,2/3,5/18,2/3,5/18).
\]
At a triangle endpoint the corresponding row of M_7z
is less than
\[
\frac7{40}+\frac5{24}+\frac{16}{75}
=\frac{179}{300}<\frac35.
\]
At the middle triangle vertex it is a<1/2=(3/5)(5/6).
At an inner arm vertex it is less than
8/25+5/72<2/5=(3/5)(2/3).
At an outer arm vertex it is less than
1/6=(3/5)(5/18).
The weighted edge inequality just used proves
||M_7||<3/5. Lowering any nonnegative edge weight
preserves this bound, as does passage to its
principal four-point and pair blocks. More precisely
the same inequality with the positive eigenvector
in (187) gives the uniform bound by lambda_7.

This is a genuine improvement over separate costs:
\[
\boxed{
\|\mathcal P_2\|+\|\mathcal P_3\|+\|\mathcal P_4\|
=\frac{3u}{4}+\frac b2,\qquad
\frac{3u}{4}+\frac b2-\lambda_7>\frac15.
}
\tag{189}
\]
On this range the u-operator has chains of length
at most three, and its three-point norm is a/sqrt 2=u/2.
The v- and 2u-operators each consist of disjoint
pairs with norms b/2 and c/2=u/4.
These maximal chains and pairs occur on intervals
of positive measure, so these are exact norms.
Since u>2/3 and b>3/5, their sum exceeds 4/5;
(188) proves the claimed saving.

Even grouping powers of 2 first still leaves a
strict loss if the prime-3 norm is then added:
\[
\|\mathcal P_2+\mathcal P_4\|
=\frac{c+\sqrt{c^2+8a^2}}4
=\frac{u(1+\sqrt{17})}{8},\qquad
\|\mathcal P_2+\mathcal P_4\|+\frac b2-\lambda_7
>\frac7{60}>\frac1{10}.
\tag{190}
\]
The three-point triangle gives the first equality
by its reflection-symmetric two-dimensional sector.
The lower estimate uses sqrt 17>4, u>2/3 and
b>3/5. Thus shared points between distinct
primitive-prime families yield additional gain
beyond grouping powers of one prime. QED.

The upper support boundary in (185) is essential.
Beyond it, extra v-edges connect the displayed
components; this seven-point decomposition cannot
be extrapolated. At its lower boundary epsilon=0
the seven-point fibres disappear. No continuity
of the operator norm across that threshold is used.

## BRT-T46. A positive eight-mode complement at L=3/4

Fix L_dagger=3/4. TPW-T3 gives
log(9/2)>451/300>3/2 and log 2<7/10<3/4.
Thus this endpoint lies strictly in (185).
Only source terms 2,3,4 are active. Let q_t be the
full Gamma/pole/scalar source minus
sum_(n=2,3,4) t_n a_n C_g(log n), with t in [0,1]^3.
Here a_4=log 2/2. The actual source has all t_n=1.
Independent diagnostic couplings are only a bounded
perturbation family; the primitive-prime insertion
family is the subfamily t_2=t_4.

At this width retain the eight normalized even
polynomials of degrees 0,2,...,14. Denote their span
by V_8, projection by P_8 and complement by Q_8.
ASM admits their source columns and preserves the
form domain under this cut. For h in Q_8D_(L_dagger),
m_0(h)=0. The same exact GRR comparison behind (174)
gives
\[
q_t[h]\ge
\left(B_*(L_\dagger)+\frac{H_{16}}2-\frac35\right)\|h\|^2,
\qquad
B_*(L_\dagger)=\frac{\Omega_0}{2}+F(3/2)
+\frac12\log2+\frac{79}{256}.
\tag{191}
\]
The polynomial contribution is exactly
47(3/4)/96-5(3/4)^3/36=79/256.
The pole square, boundary multiplier and positive
Gamma remainder have been dropped only after their
nonnegative signs and the vanishing m_0 corrections
are recorded, as in T42. Equation (188) treats all
active prime terms jointly.

A short written lower bound for the full Gamma tail
is sufficient. Expanding the positive geometric
series of k and integrating term by term gives
\[
F(3/2)=\sum_{j\ge0}
\frac{e^{-(2j+1/2)3/2}}{2j+1/2}
>2t+\frac25t^5,\qquad t=e^{-3/4}.
\]
The alternating exponential series gives
\[
t\ge
1-\frac34+\frac9{32}-\frac{27}{384}
+\frac{81}{6144}-\frac{243}{122880}
=\frac{58017}{122880}>\frac{47}{100}.
\]
Also (47/100)^2>11/50, so
(47/100)^5>5687/250000>1/50. Therefore
\[
\boxed{
F(3/2)>\frac{237}{250},\qquad
B_*(L_\dagger)>-\frac{2789}{2000}+\frac{79}{256}.
}
\tag{192}
\]
The second inequality uses the retained
Omega_0/2>-2689/1000 and log 2>693/1000.
This is an analytic tail estimate, not numerical
source evaluation.

Finally
\[
H_{16}=\frac{2436559}{720720}>\frac{169}{50}.
\]
The last comparison follows, for example, from
50(2436559)-169(720720)=26270>0.
Combining it with (191)--(192) gives
\[
\boxed{
q_t[h]\ge\frac{131}{32000}\|h\|^2
>\frac1{250}\|h\|^2
\quad(0\ne h\in Q_8D_{3/4},\ t\in[0,1]^3).
}
\tag{193}
\]
The rational difference used here is
\[
\frac{169}{100}-\frac{2789}{2000}
+\frac{79}{256}-\frac35=\frac{131}{32000}.
\]
The estimates actually make the first inequality
strict for nonzero h.

For comparison, the separate bounds
||P_2||<7/20, ||P_3||<8/25, ||P_4||<7/40
sum to 169/200. With exactly the same Gamma and
Legendre lower constants, they give only
\[
\frac{169}{100}-\frac{2789}{2000}
+\frac{79}{256}-\frac{169}{200}
=-\frac{7709}{32000}<0.
\tag{194}
\]
The joint allowance improves this specified lower
estimate by 49/200 and makes it positive.
Equation (194) is failure of that estimate, not
negativity of the true tail. QED.

No eight additional moment conditions have been
imposed on the full source family. They define the
complement in a decomposition of every admitted
profile; the retained eight coordinates and all
their mixed terms remain in the returned source.
Eight is not claimed to be the smallest possible
cut or an intrinsic observer count.

## BRT-T47. Shared prime groups retain their cross-memory terms

At L_dagger, let U_8 synthesize V_8 and define the
source blocks A_8(t), B_8(t), C_8(t) exactly as in
(178), with Q_8 in place of Q_5. Equation (193) gives
\[
\boxed{
C_8(t)\succeq\frac1{250}I,\qquad
0<C_8(t)^{-1}\preceq250I,\qquad
\mathsf S_8(t)=A_8(t)-B_8(t)^*C_8(t)^{-1}B_8(t).
}
\tag{195}
\]
For every admitted operator-domain column trial Z
with P_8Z=U_8, put R=Q_8T(t)Z. The existing exact
residual identity and the new quantitative bound are
\[
\boxed{
\mathsf S_8(t)=Z^*T(t)Z-R^*C_8(t)^{-1}R,\qquad
Z^*T(t)Z-250R^*R\preceq\mathsf S_8(t)
\preceq Z^*T(t)Z.
}
\tag{196}
\]
Every inverse-residual direction is included.
The polynomial trial Z=U_8 is admitted; improved
trials must be assessed by this complete lower
matrix, not by their source values alone.

To insert the three active powers as one group, put
E_0=U_8-C_8(0)^{-1}B_8(0) and
D_n=Q_8P_n E_0 for n=2,3,4. Set D=D_2+D_3+D_4.
PCM's unchanged insertion identity gives
\[
\begin{aligned}
\mathsf S_8(\mathbf1)
&=\mathsf S_8(\mathbf0)-E_0^*\mathcal P_{234}E_0
-D^*C_8(\mathbf1)^{-1}D,\\
0\preceq D^*C_8(\mathbf1)^{-1}D
&\preceq250D^*D
\preceq250\lambda_7^2E_0^*E_0
\preceq90E_0^*E_0.
\end{aligned}
\tag{197}
\]
The last coarse bound follows from lambda_7<3/5;
the preceding D^*D retains the actual column and
is preferable whenever its entries can be bounded.
In particular D^*D contains every D_i^*D_j term.
Grouping 2 with 4 first and then inserting 3 gives
the same exact endpoint only when all intermediate
returned-memory terms are carried. Simply adding
the group costs discards those cross terms.
No sign is assigned to the direct source term.

Equations (195)--(197) are ASM/PCM identities with
the new tail constant. The entries of S_8(1) and
the lower matrix in (196) have not been bounded
to establish their sign here. The scalar 90 is
only an inverse-memory allowance; it is not a
positivity margin or a successful full-source test.

**Gate accounting.** Prime-power chain grouping,
the exact seven-point source reduction, the
quantitative joint norm improvement, and the
specified positive eight-mode complement at 3/4
are closed in the working written derivation.
The returned source matrices at 11/20 and 3/4
remain unresolved. The unrestricted positive
interval (164) is unchanged, and all-width
positivity/RH remains open.
No new independent assumption is added.


## BRT-T48. Nested groups preserve one returned boundary form

Use the seven-point matrix M_7 from (186), with boundary
labels 6 and 24. The other five labels are internal.
Let rho>a/2 be a diagnostic spectral parameter and set
H(rho)=rho I-M_7. It is not assumed that H(rho) is positive.

The internal block D(rho) consists of the isolated middle
triangle label 12 and the two arm pairs (18,9), (8,16):
\[
D(\rho)=[\rho]\oplus
\begin{pmatrix}\rho&-a/2\\-a/2&\rho\end{pmatrix}
\oplus
\begin{pmatrix}\rho&-a/2\\-a/2&\rho\end{pmatrix}\succ0.
\]
Its exact Schur return to the two boundary labels is
\[
\boxed{
K(\rho)=
\begin{pmatrix}
\rho-\dfrac{a^2}{4\rho}
-\dfrac{b^2\rho}{4\rho^2-a^2}
&
-\dfrac c2-\dfrac{a^2}{4\rho}\\[3pt]
-\dfrac c2-\dfrac{a^2}{4\rho}
&
\rho-\dfrac{a^2}{4\rho}
-\dfrac{b^2\rho}{4\rho^2-a^2}
\end{pmatrix}.
}
\tag{198}
\]
The term a^2/(4rho) appears both on the diagonal and
off the diagonal because the two boundary points share
label 12. Each arm returns b^2 rho/(4rho^2-a^2) to
its attached endpoint. These terms are retained memory,
not independent edge costs.

More precisely, for a boundary vector v in R^2,
\[
\boxed{
\min_{y\in\mathbb R^5}(v,y)^TH(\rho)(v,y)
=v^TK(\rho)v.
}
\tag{199}
\]
This follows by completing the positive D(rho) square.
Eliminating the outer arm points first, then the inner
arm points, then label 12 gives the same minimum as
eliminating all five internal coordinates at once.
Indeed each procedure minimizes the same quadratic
function over the same internal variables; the
positive internal block makes every such conditional
minimum lawful and uniquely attained.

Thus the quantity preserved by these nested groups
is the returned boundary quadratic form at a fixed
rho, with the same boundary and source data.
Equation (199) is the finite support-graph version
of ASM-T4's coherent nested-cut identity. It is not
a new positivity theorem or an additional global
invariant supplied merely by renaming the groups.

The symmetric and antisymmetric boundary readings
diagonalize K. Its smaller eigenvalue is
\[
\begin{aligned}
f(\rho)&=\rho-\frac c2-\frac{a^2}{2\rho}
-\frac{b^2\rho}{4\rho^2-a^2},\\
f_{\rm anti}(\rho)&=
\rho+\frac c2-\frac{b^2\rho}{4\rho^2-a^2}
=f(\rho)+c+\frac{a^2}{2\rho},\\
\det H(\rho)&=
\rho(\rho^2-a^2/4)^2\det K(\rho).
\end{aligned}
\tag{200}
\]
The internal block is positive, so H(rho) is positive
exactly when K(rho) is positive, equivalently f(rho)>0.
This recovers the root equation (187) with its entire
nested return intact.

Common rescaling of observation labels preserves
this relative component and K; relabeling its two
boundary coordinates conjugates K by that permutation.
This statement concerns the same relative support
component. It does not physically shift an interval
without its boundary, change L, insert another prime,
or identify rho with a native curvature parameter.
Changing the retained boundary also changes the
returned matrix and requires its own elimination law.

## BRT-T49. The nested return sharpens the joint prime bound

On the same support band (185), the exact joint norm
satisfies the stronger bound
\[
\boxed{
\lambda_7<\frac{59}{100},\qquad
\left\|\sum_{n=2,3,4}t_n\mathcal P_n\right\|
\le\lambda_7<\frac{59}{100}
\quad(0\le t_n\le1).
}
\tag{201}
\]
This improves the former 3/5 allowance for the
same source and support; no new prime grouping is
assumed outside the seven-point band.

Here is a written rational estimate of the actual
returned boundary form. TPW-T3/T5 give
\[
\boxed{
a<\frac{123}{250},\quad a^2<\frac{243}{1000},\qquad
b<\frac{159}{250},\quad b^2<\frac{81}{200},\qquad
\frac c2<\frac{1733}{10000}.
}
\tag{202}
\]
For b, use log 3<11/10 and sqrt 3>173/100 to obtain
b<110/173<159/250. The square bound follows by
squaring this positive upper bound. The bound for
a is TPW's a/2<123/500; its displayed square
bound is weaker. Finally c/2=(log 2)/4.

Set rho=59/100. For the return through label 12,
\[
\frac{a^2}{2\rho}
<\frac{243}{1180}<\frac{103}{500}.
\]
For the complete return through either arm,
\[
4\rho^2-a^2>\frac{5747}{5000}>\frac87,\qquad
b^2\rho<\frac{4779}{20000}<\frac6{25},
\]
so
\[
\frac{b^2\rho}{4\rho^2-a^2}<\frac{21}{100}.
\]
Inserting these three complete returns in (200) gives
\[
f(59/100)>
\frac{59}{100}-\frac{1733}{10000}
-\frac{103}{500}-\frac{21}{100}
=\frac7{10000}>0.
\]
Therefore K(59/100) and H(59/100) are positive.
By (187), or by the norm characterization in T45,
lambda_7<59/100. Lowering nonnegative edge weights
preserves the norm upper bound, as already proved
for the diagnostic cube in T45. QED.

The positive K here bounds the arithmetic translation
operator. It is not the returned matrix of the full
Gamma/pole/scalar/prime source. That separate source
return is retained below.

## BRT-T50. Optimize the Gamma comparison together with its positive return

Let nu>0 be GRR's fixed kernel-comparison parameter.
It is independent of the prime coefficients and
diagnostic couplings. At any fixed finite L, let q_t
be the GRR source with every active prime term carrying
its diagnostic coupling. For mean-zero even g in the
established source domain, define
\[
\begin{aligned}
B_\nu(L)&=\frac{\Omega_0}{2}+F(2L)+\frac12\log2
+\frac L2-\frac{\nu L}{48}
-\frac{L^3}{12}-\frac{L^3}{36\nu},\\
\varepsilon(L,\nu)&=
\frac{(2L-\nu)^3-(L-\nu)^3}{576\nu}.
\end{aligned}
\]
GRR-T1/T2, including its retained positive remainder,
give
\[
\boxed{
\begin{aligned}
q_t[g]\ge{}&
\bigl(B_\nu(L)+\varepsilon(L,\nu)\bigr)\|g\|^2
+J_L[g]+\int_Iw_L|g|^2+\ell(g)^2\\
&-\left\langle g,\sum_n t_n\mathcal P_n g\right\rangle .
\end{aligned}
}
\tag{203}
\]
Every active prime power is included in the last sum.
Here the two polynomial moment corrections vanish
because m_0(g)=0. The positive surplus is the same
Gamma remainder from GRR, retained on L<=u<=2L
using N_g-C_g(u)>=N_g/2. It is not an extra source
term and is not counted a second time.

Adding that surplus before choosing nu gives
\[
\begin{aligned}
\mathcal G_\nu(L)&:=B_\nu(L)+\varepsilon(L,\nu)\\
&=\frac{\Omega_0}{2}+F(2L)+\frac12\log2
+\frac L2-\frac{L^3}{12}-\frac{L^2}{64}
-\frac L{64}\left(\nu+\frac{L^2}{\nu}\right),\\
\boxed{
\mathcal G_L(L)-\mathcal G_\nu(L)
=\frac{L(\nu-L)^2}{64\nu}\ge0.
}
\end{aligned}
\tag{204}
\]
Thus nu=L maximizes this particular uniform
comparison, with the displayed positive remainder
included. The maximizer is a fixed parameter for
each declared width and does not depend on g.
No nonlinear trial-dependent source is introduced.
This is optimality within the stated comparison
family, not optimality among all possible Gamma
estimates.

To verify (204), expand the numerator of epsilon
as 7L^3-9L^2 nu+3L nu^2. Its contributions combine
with -nu L/48 and -L^3/(36nu) to give
-L nu/64 and -L^3/(64nu), respectively.
The difference from nu=L is the displayed square.
QED.

At L=3/4 the width-dependent polynomial-plus-surplus
contribution is
\[
\frac L2-\frac{L^3}{12}-\frac{3L^2}{64}
=\frac{321}{1024}.
\]
The previous T46 contribution was 79/256=316/1024,
so the new uniform source allowance gains 5/1024.
For comparison, retaining GRR's surplus while still
choosing nu=1/2 gives 639/2048; choosing nu=L raises
it to 642/2048. The parameter improvement alone is
3/2048, and the rest of the gain restores the
previously dropped positive remainder.

## BRT-T51. The same full residual now has inverse allowance 53

Keep L_dagger=3/4, the same eight retained modes
of degrees 0,2,...,14, and the same actual source
family q_t from T46. No source, domain, trial
normalization or cut dimension changes.

For h in Q_8D_(3/4), apply (203)--(204) with
nu=L=3/4, retain J>=H_16||h||^2/2, and use (201).
The old scalar bounds (192), before their polynomial
addition, give
Omega_0/2+F(3/2)+(log 2)/2>-2789/2000.
Also H_16/2>169/100. Therefore
\[
\boxed{
q_t[h]\ge
\frac{2429}{128000}\|h\|^2
>\frac1{53}\|h\|^2
\quad(0\ne h\in Q_8D_{3/4},\ t\in[0,1]^3).
}
\tag{205}
\]
The rational coefficient is
\[
\frac{169}{100}-\frac{2789}{2000}
+\frac{321}{1024}-\frac{59}{100}
=\frac{2429}{128000}.
\]
It exceeds 1/53 because 2429(53)=128737>128000.
Relative to (193), the exact change in this lower
coefficient is
\[
\frac{2429}{128000}
=\frac{131}{32000}+\frac1{100}+\frac5{1024}.
\]
The two positive increments come from the sharper
complete prime return and the Gamma comparison
with its positive remainder, respectively.

The same tail operators consequently obey
\[
\boxed{
C_8(t)\succeq\frac{2429}{128000}I,\qquad
0<C_8(t)^{-1}\preceq\frac{128000}{2429}I
\preceq53I.
}
\tag{206}
\]
For the same admitted column trial Z, with
P_8Z=U_8 and its full R=Q_8T(t)Z, define
\[
\begin{aligned}
\mathcal L_{250}(Z,t)&=Z^*T(t)Z-250R^*R,\\
\mathcal L_{53}(Z,t)&=Z^*T(t)Z-53R^*R .
\end{aligned}
\]
The exact source return and improved lower matrix are
\[
\boxed{
\begin{aligned}
\mathcal L_{53}(Z,t)
&\preceq\mathsf S_8(t)
=Z^*T(t)Z-R^*C_8(t)^{-1}R,\\
\mathcal L_{53}(Z,t)-\mathcal L_{250}(Z,t)
&=197R^*R\succeq0.
\end{aligned}
}
\tag{207}
\]
Thus this continuation quantitatively improves the
same source-minus-complete-memory estimate. It
does not choose a correction merely because its
trial value falls. Every residual direction and
mixed prime-memory pairing remains.

For the concrete admitted polynomial trial Z=U_8,
the improvement is positive definite, not merely
semidefinite. To see this, put B_8(t)=Q_8T(t)U_8.
Every nonzero retained vector produces a nonzero
polynomial p of degree at most 14. The endpoint
argument in T9 extends to this source:
\[
\boxed{
T(t)p(x)=-\frac14p(x)\log(L-x)+A_{p,t}(x)
\quad\text{near }x=L,\qquad
\ker B_8(t)=\{0\}.
}
\tag{208}
\]
Here A_(p,t) is analytic across L as an extension
of the one-sided interior formula. The Gamma and
pole terms have the same local form as in T9.
For each active shift log 2, log 3, log 4<2L,
the inward translate of p is an interior polynomial
near L and the outward translate is zero there.
Multiplying these analytic contributions by t_n
does not alter the boundary logarithm.

A nonzero p has finite endpoint vanishing order,
so some derivative of p(x)log(L-x) still diverges
logarithmically. Thus T(t)p cannot be a polynomial
near L. If Q_8T(t)p=0, it would instead agree
almost everywhere with an element of V_8, hence
with a polynomial throughout that interior
neighborhood by continuity, a contradiction.
This proves injectivity and B_8(t)^*B_8(t)>0.

Consequently (207) gives the strict matrix improvement
197 B_8(t)^*B_8(t)>0 for that fixed polynomial trial.
No rational lower eigenvalue for this Gram matrix
has been established. The statement is for this
fixed finite cut; it does not assert an error floor
over all polynomial trials or increasing cut sizes.
Residual cancellation remains unnecessary.

The simultaneous prime insertion bound also improves.
With the same E_0 and D=D_2+D_3+D_4 as in (197),
\[
\boxed{
0\preceq D^*C_8(\mathbf1)^{-1}D
\preceq53D^*D
\preceq53\lambda_7^2E_0^*E_0
\preceq\frac{184493}{10000}E_0^*E_0
\preceq19E_0^*E_0.
}
\tag{209}
\]
This replaces the former coarse allowance 90.
The bound involving the actual D^*D retains all
cross terms and is preferable to the last scalar
allowance when estimating entries.

Neither (207)'s improved lower matrix nor the
exact S_8(1) has been shown positive semidefinite
here. A strict improvement of a lower matrix
does not by itself determine its sign. The
comparison parameter optimization and nested
return identity do not supply that missing sign.

**Gate accounting.** The explicit nested boundary
return, sharper joint norm, optimized Gamma
comparison with retained surplus, and reduction
of the full residual allowance from 250 to 53
are closed in the working written derivation.
The same returned source sign at 3/4 remains
open, as does the retained sign at 11/20.
The unrestricted interval (164) is unchanged;
all-width positivity/RH remains unresolved.
No new independent assumption is added.


## BRT-T52. Containment requires a coherent positivity transfer

For a finite set of primes S, the generated multiplicative
group is
\[
G(S)=\left\{\prod_{p\in S}p^{k_p}:k_p\in\mathbb Z\right\}
\subset\mathbb Q_{>0}^{\times}.
\tag{210}
\]
Multiplication adds exponents and inversion negates them,
so this is a group; S contained in T implies
G(S) contained in G(T). The set of primes itself is not
closed under multiplication. None of these containment
facts controls the sign of the arithmetic source.

Even two positive one-coordinate restrictions do not
force a positive joint form: x^2+y^2+4xy restricts to a
positive square on each axis but takes value -2 at
(1,-1). This is only a counterexample to the proposed
logical implication, not a negative test for q_L.

There is, however, a precise local-to-joint rule for the
actual nonnegative edge weights in (186). For any
symmetric zero-diagonal matrix M with nonnegative
entries and any positive vector w,
\[
\boxed{
\begin{aligned}
x^T(\rho I-M)x
={}&\sum_{i<j}M_{ij}w_iw_j
\left(\frac{x_i}{w_i}-\frac{x_j}{w_j}\right)^2\\
&+\sum_i\left(\rho-\frac{(Mw)_i}{w_i}\right)x_i^2.
\end{aligned}
}
\tag{211}
\]
Expanding each edge square gives its two diagonal
terms and -2M_ij x_i x_j, proving the identity.
If the final diagonal coefficients are nonnegative,
these smallest edge contributions assemble to a
positive semidefinite whole. Their common diagonal
budget is counted once, even at a shared point.
Separate positive block tests with incompatible
budgets do not supply this conclusion.

This rule is concrete for M=M_7 and rho=59/100.
In the order 6,12,24,18,9,8,16, use
\[
w=(1,q,1,r,s,r,s),\qquad
q=\frac a\rho,\quad
r=\frac{2b\rho}{4\rho^2-a^2},\quad
s=\frac{ar}{2\rho}.
\]
These entries are positive. At the five internal
vertices (Mw)_i=rho w_i. At the two boundary
vertices, rho-(Mw)_i/w_i=f(rho), where T49 proved
f(rho)>7/10000. Thus (211) becomes the exact
factorization
\[
\boxed{
x^T(\rho I-M_7)x
=\sum_{i<j}(M_7)_{ij}w_iw_j
\left(\frac{x_i}{w_i}-\frac{x_j}{w_j}\right)^2
+f(\rho)(x_6^2+x_{24}^2).
}
\tag{212}
\]
Here subscripts 6 and 24 are the observation labels.
If this expression vanishes, the endpoint readings
vanish and all connected edge ratios x_i/w_i are
equal. The seven-point graph is connected, so x=0.
Hence the comparison is positive definite.

The edge squares can be grouped into the triangle
and the two arms, or subdivided and regrouped;
their sum and shared boundary budget remain the
same. Principal subcomponents inherit this positive
comparison. Lowering any active nonnegative edge
weight leaves a nonnegative budget using the same w,
so the diagnostic coupling family is covered too.

This is a verified transfer from elementary edge
squares to this joint prime comparison. It explains
the valid part of the smallest-block intuition.
It does not prove positivity of the full source:
the comparison bounds the subtracted prime operator
by rho I, and the Gamma, scalar, pole and complete
returned memory still have to be included.
The identity itself is an explicit factorization
of the already established prime comparison, not
a new RH positivity theorem.

## BRT-T53. One controlled tail level propagates to every higher level

Keep L=3/4, the eight retained modes of degrees
0,2,...,14, and all diagnostic couplings t in [0,1]^3.
Use degree labels phi_(2j), and put
\[
\begin{gathered}
\delta=\frac{2429}{128000},\qquad
d_j=\delta+\frac{H_{2j}-H_{16}}2\quad(j\ge8),\\
\mathcal D\phi_{2j}=d_j\phi_{2j}
\quad\text{on }Q_8H_L.
\end{gathered}
\tag{213}
\]
The positive diagonal operator D is defined by its
closed quadratic form, with domain consisting of
the sequences for which sum d_j|h_j|^2 is finite.
It is a comparison operator, not an assertion that
the actual C_8(t) is diagonal.

Retaining ASM's full Legendre energy bound, rather
than replacing every level by H_16, gives
\[
\boxed{
q_t[h]\ge
\sum_{j\ge8}d_j|\langle h,\phi_{2j}\rangle|^2,
\qquad C_8(t)\succeq\mathcal D
\quad(h\in Q_8D_L).
}
\tag{214}
\]
Indeed T50--T51 and the full prime comparison give
\[
q_t[h]\ge
J_L[h]+\left(-\frac{3969}{2000}
+\frac{321}{1024}\right)\|h\|^2.
\]
ASM-T1 gives J_L[h]>=sum H_(2j)|h_j|^2/2.
The coefficient at degree 16 exceeds delta because
H_16/2>169/100; subtracting H_16/2 at every higher
level yields exactly (213)--(214).
This also proves that the source form domain is
contained in the displayed comparison form domain.
No equality of those domains is required.

The propagation rule is explicit:
\[
\boxed{
d_8=\delta>0,\qquad
d_{j+1}-d_j
=\frac1{2(2j+1)}+\frac1{2(2j+2)}>0.
}
\tag{215}
\]
Thus the first controlled tail level and a proved
positive increment control every higher level,
including the infinite complement. This implication
uses the same global source comparison, which has
already bounded all prime couplings. It is not
inferred from subgroup containment or separate
diagonal readings of the actual source.

For later use,
\[
\boxed{
d_8>\frac1{53},\qquad
d_9>\frac3{40},\qquad
d_j>\frac18\quad(j\ge10).
}
\tag{216}
\]
The first inequality is T51. For the second,
delta>3/160 and
(H_18-H_16)/2=35/612>9/160.
For the third,
(H_20-H_18)/2=39/760>1/20,
and 3/40+1/20=1/8. Monotonicity then covers all
j>=10. QED.

The modes here are comparison-energy levels,
not individual primes or intrinsic prime classes.
The transfer rule has a specified source and
width, and does not extend the support band.

## BRT-T54. Retain that hierarchy in the complete inverse residual

The form inequality (214) implies
\[
\boxed{
0<C_8(t)^{-1}\preceq\mathcal D^{-1},\qquad
\langle r,C_8(t)^{-1}r\rangle
\le\sum_{j\ge8}\frac{|\langle r,\phi_{2j}\rangle|^2}{d_j}.
}
\tag{217}
\]
No commutation of C_8(t) with the Legendre projections
is assumed. To prove the inverse comparison, use
\[
\langle r,C_8(t)^{-1}r\rangle
=\sup_{h\in Q_8D_L}\{2\langle r,h\rangle-q_t[h]\}.
\]
Replacing q_t by its smaller comparison form and
enlarging the supremum to that form's domain yields
<r,D^(-1)r>. The positive lower bound for D makes
this finite and proves (217).

For an admitted column trial Z with P_8Z=U_8, let
R=Q_8T(t)Z be its full residual as before. Write
r_16=phi_16^*R and r_18=phi_18^*R for its two row
readings, and
\[
R_{\ge20}
=\left(I-P_8-|\phi_{16}\rangle\langle\phi_{16}|
-|\phi_{18}\rangle\langle\phi_{18}|\right)R .
\]
Then
\[
\boxed{
\begin{aligned}
R^*C_8(t)^{-1}R
&\preceq\sum_{j\ge8}d_j^{-1}
(\phi_{2j}^*R)^*(\phi_{2j}^*R)\\
&\preceq
53\,r_{16}^*r_{16}
+\frac{40}{3}r_{18}^*r_{18}
+8\,R_{\ge20}^*R_{\ge20}
=:\mathcal E_{\rm hier}(R).
\end{aligned}
}
\tag{218}
\]
The infinite series is a finite-matrix norm-convergent
sum: every residual column is in L2, its squared
coefficient sum is finite, and d_j^(-1)<=delta^(-1).
Thus every higher residual mode is included.
The coarse three-level bound uses (216).

The complete lower source matrix is now
\[
\boxed{
\begin{aligned}
\mathcal L_{\rm hier}(Z,t)
&:=Z^*T(t)Z-\mathcal E_{\rm hier}(R)
\preceq\mathsf S_8(t),\\
\mathcal L_{\rm hier}(Z,t)-\mathcal L_{53}(Z,t)
&=\frac{119}{3}r_{18}^*r_{18}
+45R_{\ge20}^*R_{\ge20}\succeq0.
\end{aligned}
}
\tag{219}
\]
This follows from the unchanged exact identity
S_8=Z^*TZ-R^*C_8^(-1)R and the orthogonal decomposition
R^*R=r_16^*r_16+r_18^*r_18+R_(>=20)^*R_(>=20).
It is a quantitative improvement for the same
source, cut and trial, retaining its nonzero residual.

For the polynomial trial Z=U_8 this improvement is
positive definite. If R_(>=20)v vanished for a nonzero
retained vector v, then T(t)U_8v would be a polynomial
of degree at most 18. The endpoint logarithm in (208)
rules out agreement with any polynomial for this
nonzero polynomial input. Therefore R_(>=20) is
injective, and its finite Gram is positive definite.
No rational lower eigenvalue for that Gram is claimed.

The simultaneous prime insertion bound has the
same refinement: in (209), replace 53 D^*D by
E_hier(D), with D=D_2+D_3+D_4 unchanged.
Every row Gram of D and the full higher residual
Gram retain all mixed D_i^*D_j terms. Splitting into
orthogonal comparison levels does not assert that
those prime columns are independent. QED.

This hierarchy permits progressively sharper inverse
weights without adding retained source coordinates.
The actual source is not diagonalized and no
residual component is discarded. The entries of
the improved complete lower matrix have not yet
been bounded to establish its sign.

## BRT-T55. The same comparison fails one level below the controlled tail

It is useful to test the proposed downward shortcut
explicitly. With only seven retained polynomial
modes, degrees 0,2,...,12, the complement includes
degree 14. The diagonal comparison used above
has coefficient there
\[
\boxed{
\beta_{14}
=\frac{H_{14}}2-\frac{3969}{2000}
+\frac{321}{1024}
<-\frac{5251}{128000}<-\frac1{25}.
}
\tag{220}
\]
Indeed H_14/2=1171733/720720<163/100, and
163/100-3969/2000+321/1024=-5251/128000.
The admitted profile phi_14 therefore has a
negative reading in this particular lower
comparison. This comparison cannot certify
a positive seven-mode complement.

This does not prove q_t[phi_14]<0, does not
exclude a stronger seven-mode estimate, and
does not make eight a universal minimum number
of observers. It identifies the exact limit
of the current propagation argument.

Above degree 16, (215) propagates positivity to
every higher comparison level because the full
couplings have already been bounded. Below that
threshold, containment supplies no missing sign.
The eight retained coordinates and their complete
returned memory must still be assessed together.

**Gate accounting.** The actual prime comparison
has an explicit coherent edge-square transfer;
the complete infinite tail has a proved upward
propagation rule; and its inverse residual
allowances improve from uniform 53 to 53,40/3,8
on the stated levels. These are closed in the
working written derivation. The attempted
seven-mode shortcut fails at the explicit lower
coefficient (220). Returned source signs at
11/20 and 3/4 and all-width positivity/RH remain
open; the unrestricted interval (164) is unchanged.
No new independent assumption is added.

## BRT-T56. A newly active prime is small relative to boundary energy

23 October 2026. This continuation uses the same arithmetic source
and form domain. It retains a part of the Gamma boundary energy
that the uniform operator-norm estimate in T42 discarded.

For L<u<2L, put
\[
\begin{gathered}
E_u=(-L,L-u)\ \cup\ (u-L,L),\\
\mathcal W_L[g]=\int_I w_L(x)|g(x)|^2\,dx,\qquad
w_L(x)=-\frac14\log(1-(x/L)^2),\\
w_{\min}(L,u)=\frac14
 \log\frac{L^2}{u(2L-u)}>0 .
\end{gathered}
\tag{221}
\]
The two intervals are disjoint. For every admitted real profile,
\[
\boxed{
 |C_g(u)|\le\frac12\int_{E_u}|g|^2
 \le\frac{1}{2w_{\min}(L,u)}
              \int_{E_u}w_L|g|^2 .
}
\tag{222}
\]
Indeed write C_g(u) as the integral of g(x)g(x+u) on
(-L,L-u), apply 2|rs|<=r^2+s^2 and change variables in the
second square. On E_u one has |x|>=u-L, giving the displayed
minimum of w_L. This argument imposes no pointwise sign on g.

For a fixed positive translation weight a_u, (222) gives
\[
\boxed{
 |a_u C_g(u)|\le
 \rho_u(L)\mathcal W_L[g],\qquad
 \rho_u(L)=
 \frac{2a_u}{\log[L^2/(u(2L-u))]}
 \longrightarrow0
 \quad\left(L\downarrow\frac u2\right).
}
\tag{223}
\]
At L=u/2 the overlap is empty. The limiting statement is about
this relative form coefficient, with the weight w_L carried
with its window. It is not operator-norm continuity of the
truncated translation on L2. In particular it is compatible
with T41's operator-norm jump at prime activation.

The weighted integral is finite on D_L. The retained logarithmic
decomposition (GRR equation (4)) has nonnegative J_L and boundary
terms; its difference from the Gamma form on a finite interval
is a bounded quadratic form. The inequalities therefore extend
from the smooth core to the same closed form domain. No new
boundary condition or trial regularity is assumed.

Specialize now to L_*=11/20, v=log 3 and b=v/sqrt 3.
The positive atanh series gives the explicit lower bound
\[
\log3=\sum_{j=0}^{\infty}\frac{1}{(2j+1)4^j}
>
1+\frac1{12}+\frac1{80}+\frac1{448}
=\frac{7379}{6720}>\frac{549}{500}.
\tag{224}
\]
The last cross-product difference is
7379*500-549*6720=220>0. Together with the retained
log 3<11/10, this locates the nonempty prime-3 overlap.
On its two caps,
\[
\frac{|x|}{L_*}>\frac{274}{275},\qquad
1-(x/L_*)^2<\frac{549}{75625}<\frac1{128}.
\]
The last comparison follows from 549*128=70272<75625.
Since log 2>693/1000, already used in T42,
\[
 w_{L_*}(x)>
 \frac{\log128}{4}=\frac{7\log2}{4}>
 \frac{4851}{4000}>\frac65 .
\]
Also b<2/3 by T41. Consequently
\[
\boxed{
 |bC_g(v)|
 \le\frac{5}{18}\int_{E_v}w_{L_*}|g|^2
 \le\frac{5}{18}\mathcal W_{L_*}[g].
}
\tag{225}
\]
This spends at most 5/18 of the existing boundary energy.
It does not simultaneously spend the whole joint prime norm
from (171). The rest of the boundary energy remains available.

## BRT-T57. Positive tails after three, four or five retained modes

Keep the actual diagnostic family q_t in (172), t in [0,1]^2,
at exactly L_*=11/20. Since log 2>L_*, its two translation caps
are disjoint as well. With a=(log 2)/sqrt 2<1/2, (222)'s first
inequality gives |aC_g(log 2)|<=||g||^2/4.

For any mean-zero h in D_(L_*), the two moment terms in GRR's
polynomial comparison vanish. Retain (174)--(175), charge
prime 2 to the norm and prime 3 to (225). This yields
\[
\boxed{
 q_t[h]\ge
 J_{L_*}[h]-\frac{467}{400}\|h\|^2
 +\frac{13}{18}\mathcal W_{L_*}[h]+\ell(h)^2 .
}
\tag{226}
\]
Here 467/400=367/400+1/4. The nonnegative unused Gamma
remainder is dropped only once. There is no cancellation
assumption between the two prime readings.

For an integer N>=1 define
\[
 V_N=\operatorname{span}\{\phi_0,\phi_2,\ldots,\phi_{2N-2}\},
 \quad Q_N=I-P_N,\qquad
 d_N=\frac{H_{2N}}2-\frac{467}{400}.
\]
ASM's Legendre energy bound, in this note's degree notation,
then gives
\[
\boxed{
 q_t[h]\ge
 \sum_{j\ge N}\left(\frac{H_{2j}}2-\frac{467}{400}\right)
       |\langle h,\phi_{2j}\rangle|^2
 +\frac{13}{18}\mathcal W_{L_*}[h]+\ell(h)^2
 \ge d_N\|h\|^2
 \quad(h\in Q_ND_{L_*}).
}
\tag{227}
\]
The series inequality follows first for finite projections;
the nonnegative remainder energy and form closure give the
full statement. All higher even modes are retained.

In particular, the positive restricted source operators obey
\[
\boxed{
 C_3(t)\succeq\frac{23}{400}I,\qquad
 C_4(t)\succeq\frac{67}{350}I,\qquad
 C_5(t)\succeq\frac{1871}{6300}I
 \quad(t\in[0,1]^2).
}
\tag{228}
\]
The rational identities are
\[
\begin{aligned}
 \frac{H_6}{2}-\frac{467}{400}
 &=\frac{49}{40}-\frac{467}{400}=\frac{23}{400},\\
 \frac{H_8}{2}-\frac{467}{400}
 &=\frac{761}{560}-\frac{467}{400}=\frac{67}{350},\\
 \frac{H_{10}}{2}-\frac{467}{400}
 &=\frac{7381}{5040}-\frac{467}{400}=\frac{1871}{6300}.
\end{aligned}
\tag{229}
\]
Closedness, density and the self-adjoint restricted operators
are those of ASM-T1--T3; the polynomial cuts preserve the form
domain. The source and every actual prime weight are unchanged.

With only two retained modes this particular comparison gives
d_2=-151/1200. That is a failed lower estimate, not a negative
actual-source direction. Three is a sufficient retained count
for this estimate at this width, not a claimed universal or
minimal observer count.

## BRT-T58. The same five-column residual has a smaller complete cost

Use the five retained columns and notation A(t), B(t), C(t),
S(t), Z and R from (178)--(179). Equation (228) replaces the
earlier scalar inverse allowance by
\[
\boxed{
 0<C(t)^{-1}\preceq\frac{6300}{1871}I,\qquad
 Z^*T(t)Z-\frac{6300}{1871}R^*R
 \preceq\mathsf S(t)\preceq Z^*T(t)Z .
}
\tag{230}
\]
This is the existing exact identity
S=Z*TZ-R*C^(-1)R, with a stronger bound for the same C.
The trial must still satisfy P_5Z=U and have operator-domain
columns. The residual R=Q_5T(t)Z is complete.

Relative to (179), for the same admitted Z the lower matrix
improves by exactly
\[
\boxed{
 \left(\frac{40}{3}-\frac{6300}{1871}\right)R^*R
 =\frac{55940}{5613}R^*R\succeq0 .
}
\tag{231}
\]
The new scalar coefficient is 945/3742 of the former one.
Strict improvement in a given retained direction requires
its residual column to be nonzero. A smaller error allowance
does not establish that the improved lower matrix is positive.

There is also a stronger weighted version. Define the positive
diagonal comparison on Q_5H by
\[
 \mathcal D_5\phi_{2j}
 =\left(\frac{H_{2j}}2-\frac{467}{400}\right)\phi_{2j},
 \qquad j\ge5 ,
\]
with its maximal diagonal form domain. Equation (227) means
C(t)>=D_5 in form order, including the needed domain inclusion.
The inverse quadratic-form variational principle therefore gives
C(t)^(-1)<=D_5^(-1), without a commutation assumption. If
r_(2j)=<phi_(2j),R> is a row on the five retained coordinates,
\[
\boxed{
 R^*C(t)^{-1}R
 \preceq
 \sum_{j=5}^{\infty}
 \frac{r_{2j}^*r_{2j}}{H_{2j}/2-467/400}
 \preceq\frac{6300}{1871}R^*R .
}
\tag{232}
\]
The series converges in finite-matrix norm: its weights are
bounded by 6300/1871 and all five residual columns lie in L2.
Every component of a combined source is squared together in
each row, so mixed prime-memory terms remain.

For example, the exact simultaneous insertion identity (180),
with D=Q_5 P_23 E_0, now has the enclosure
\[
\boxed{
 \mathsf S(\mathbf1)\succeq
 \mathsf S(\mathbf0)-E_0^*\mathcal P_{23}E_0
 -\frac{6300}{1871}D^*D .
}
\tag{233}
\]
The direct signed term is not assigned an artificial sign,
and D_2*D_3+D_3*D_2 is not dropped. Equation (232) can sharpen
the same complete D further. These are quantitative uses of
the inherited inverse-order and Schur identities.

## BRT-T59. An exact three-dimensional sign target at this width

For N=3,4,5 use the actual source blocks A_N(t), B_N(t),
C_N(t) and synthesis U_N defined as in ASM-T3. All their
infinite inverses exist by (228). Put
\[
 \mathsf S_N(t)=A_N(t)-B_N(t)^*C_N(t)^{-1}B_N(t).
\]
Then, on the entire original test family,
\[
\boxed{
 q_t\ge0\text{ on }D_{L_*}
 \quad\Longleftrightarrow\quad
 \mathsf S_3(t)\succeq0
 \quad\Longleftrightarrow\quad
 \mathsf S_5(t)\succeq0 .
}
\tag{234}
\]
In particular the full source has at most three negative
directions at this fixed width. This is not a claim that any
negative direction exists or that none exists.

To exhibit the exact relation, split the five-coordinate
returned matrix into degrees (0,2,4) and (6,8):
\[
 \mathsf S_5(t)=
 \begin{pmatrix}H(t)&F(t)^*\\F(t)&D(t)\end{pmatrix}.
\]
Then
\[
\boxed{
 D(t)\succeq\frac{23}{400}I_2,\qquad
 \mathsf S_3(t)=H(t)-F(t)^*D(t)^{-1}F(t).
}
\tag{235}
\]
For the first assertion minimize q_t over Q_5 at fixed
degree-(6,8) coordinates b. Every candidate lies in Q_3D_L,
so (228) gives
q_t>=23(||b||^2+||h||^2)/400>=23||b||^2/400.
The attained minimum is b*D(t)b. For the second assertion,
minimize first over Q_5 and then over the remaining two
coordinates. ASM's coherent fibre-minimization identity
identifies this with minimization over all of Q_3. No memory
is omitted when the two intermediate coordinates are removed.

The three-coordinate carrier is the Legendre cut at L=11/20.
It is not the separate five-source eta-zero endpoint carrier
in MP, and it does not change that route's assumptions.

The next sign task is now explicit: establish S_3(1,1)>=0,
or establish a positive five-column lower matrix in (230)
or (232). Merely checking the retained block A_3, A_5 or H
does not do so. Neither actual returned sign has been
evaluated or proved in this continuation.


## BRT-T60. An exact source column with a bounded polynomial remainder

23 October 2026. Fix L=11/20 and the actual couplings t=(1,1).
Only prime powers 2 and 3 are active. Write u=log 2, v=log 3,
a=u/sqrt(2), b=v/sqrt(3), and let all translations use zero extension.
For this continuation, the unnormalized retained columns are
psi_j(x)=P_(2j)(x/L), j=0,...,4. Their squared norms are
n_j=2L/(4j+1). Thus the matrices below are congruent to, but do not
have the same numerical eigenvalues as, the earlier normalized matrices.

Put j(s)=k(s)-1/(2s), with j(0)=1/4, and let W denote
multiplication by w_L from (221). The actual operator on its
polynomial domain is
\[
\boxed{
\begin{aligned}
T&=J_{\rm op}+W+c_LI+|h\rangle\langle h|
       -\frac12K_j-\mathcal P_{23},\\
c_L&=-\frac12\{\gamma_E+\log(2\pi L)\},\qquad h(x)=\cosh(x/2),\\
(K_jf)(x)&=\int_{-L}^{L}j(|x-y|)f(y)\,dy,\qquad
J_{\rm op}P_n(x/L)=\frac{H_n}{2}P_n(x/L).
\end{aligned}
}
\tag{236}
\]
Here P_23 is the same signed-source subtraction operator
a(S_u+S_(-u))/2+b(S_v+S_(-v))/2; neither prime is suppressed.

**Proof.** The Gamma column is the one used in T4. Split its
integral into 1/(2s) and j(s). Since
F(s)=atanh(exp(-s/2))+atan(exp(-s/2)),
\[
F(s)+\int_0^s j(r)\,dr
 =\log2+\frac\pi4-\frac12\log s .
\tag{237}
\]
Both sides have derivative -1/(2s), and their constants agree
as s decreases to zero. Adding the two endpoint diagonal terms
and Omega_0/2 gives exactly c_L+w_L. The remaining integral is
-K_j/2 and ASM's singular Legendre operator is J_op. This proves
(236) with the original scalar constant retained.

For a finite and rigorously bounded evaluation, define
j_48(s)=sum_(m=0)^48 j_m s^m by the exact rational recurrence
\[
\begin{aligned}
f_n&=\frac12\left\{\frac1{2^nn!}
 -\sum_{q=1}^{\lfloor n/2\rfloor}
       \frac{2f_{n-2q}}{(2q+1)!}\right\},\qquad
j_m=f_{m+1},\\
(j_0,j_1,j_2,j_3,j_4)
 &=\left(\frac14,-\frac1{48},-\frac1{32},
                  \frac7{11520},\frac5{1536}\right).
\end{aligned}
\tag{238}
\]
Indeed sum f_n z^n=e^(z/2)/(2 sinh(z)/z). Also put
h_16(x)=sum_(r=0)^8 (x/2)^(2r)/(2r)!.
Then
\[
\boxed{
\begin{aligned}
\sup_{0\le s\le11/10}|j(s)-j_{48}(s)|
 &\le\epsilon_j:=\frac{100}{9}\left(\frac{11}{20}\right)^{49},\\
\sup_{|x|\le11/20}|h(x)-h_{16}(x)|
 &\le\epsilon_h:=\frac{2(11/40)^{18}}{18!}.
\end{aligned}
}
\tag{239}
\]
For the first estimate, on |z|<=2 one has
|sinh(z)/z-1|<=sinh(2)/2-1<5/6. For example, the first four
terms of sinh(2) sum to 1142/315 and its remaining tail is
at most (4/2835)/(1-1/25), giving sinh(2)<11/3.
Together with e<3, this gives |j(z)|<5 on |z|=2.
Cauchy's coefficient estimate and the geometric tail with
ratio 11/20 give (239). The cosh estimate follows from its
positive Taylor tail, whose successive-term ratio is below 1/2.

Let bar(T) be (236) with only j and h replaced by j_48 and h_16.
It keeps W, c_L and both actual prime translations. Schur's
kernel bound and the rank-one difference identity give
\[
\boxed{
 \|T-\overline T\|\le
 \epsilon:=L\epsilon_j+8L\epsilon_h
 <\frac{1159}{10^{15}}.
}
\tag{240}
\]
For the rank-one part use ||h||,||h_16||<=2 sqrt(2L) and
||h-h_16||<=sqrt(2L) epsilon_h. The difference is a bounded
operator on H_L. The unbounded logarithmic multiplier and
the source domain have not been approximated away. QED.

## BRT-T61. A fixed rational trial encloses the complete returned matrix

The packet [proof_lab/rh_l055](proof_lab/rh_l055/README.md) fixes
a rational 12-by-5 coefficient array c in trial.json, with common
denominator 10^14. It defines
\[
 Z_k(x)=\sum_{j=0}^{11}c_{jk}\psi_j(x),\qquad
 c_{jk}=\delta_{jk}\ (0\le j,k<5),\qquad
 P_5Z=(\psi_0,\ldots,\psi_4).
\tag{241}
\]
The trial is a finite choice of admitted columns, not a truncation
of the source space. Its origin was a numerical candidate search;
its coefficients are now fixed rationals. Candidate quadrature
and floating-point eigenvalues play no role in the certificate.

Here are the analytic integration identities used by the verifier.
Set y=(x+L)/(2L), so
P_n(2y-1)=sum_(r=0)^n (-1)^(n+r) binom(n,r) binom(n+r,r)y^r.
For nonnegative integers m,r,
\[
\int_0^1|y-s|^m s^r\,ds
=\sum_{i=0}^{m}
 \frac{(-1)^i\binom mi}{m-i+r+1}y^i
 +(1-(-1)^m)\frac{m!r!}{(m+r+1)!}y^{m+r+1}.
\tag{242}
\]
This follows by splitting the integral at s=y and integrating
the binomial expansion. Thus the nonsingular kernel columns
are exact polynomials before constant enclosures.

The boundary weight is -[log4+log y+log(1-y)]/4. With n=r+1,
the necessary moments, all over [0,1], are
\[
\begin{aligned}
\int y^r\,dy&=\frac1n,&
\int y^r\log y\,dy&=-\frac1{n^2},&
\int y^r\log(1-y)\,dy&=-\frac{H_n}{n},\\
\int y^r(\log y)^2\,dy&=\frac2{n^3},&
\int y^r(\log(1-y))^2\,dy&=\frac{H_n^2+H_n^{(2)}}n,\\
\int y^r\log y\log(1-y)\,dy
 &=\frac{H_n}{n^2}+\frac{H_n^{(2)}}n-\frac{\pi^2}{6n}.
\end{aligned}
\tag{243}
\]
These follow by differentiating the beta integral, or by
integrating the convergent logarithmic series. The necessary
weighted piecewise moments use the exact primitives
\[
\begin{aligned}
\int y^r\log y\,dy
 &=y^n\left(\frac{\log y}{n}-\frac1{n^2}\right),\\
\int y^r\log(1-y)\,dy
 &=\frac{(y^n-1)\log(1-y)-\sum_{k=1}^{n}y^k/k}{n}.
\end{aligned}
\tag{244}
\]
Endpoint limits are taken analytically. The five pieces have
endpoints 0, 1-v/(2L), 1-u/(2L), u/(2L), v/(2L), 1.
On each piece the translated polynomials are retained exactly.
All physical integrals carry the factor 2L.

Write bar(R)=Q_5 bar(T)Z. Its full Gram, not a sampled tail Gram, is
\[
\begin{aligned}
\overline A_{ij}&=\langle Z_i,\overline T Z_j\rangle,\qquad
B_{ki}=\langle\psi_k,\overline T Z_i\rangle,\\
\overline G_{ij}
 &=\langle\overline T Z_i,\overline T Z_j\rangle
        -\sum_{k=0}^{4}\frac{B_{ki}B_{kj}}{n_k}
   =(\overline R^*\overline R)_{ij},\\
t_z&=\operatorname{tr}(Z^*Z),\qquad
t_r=\operatorname{tr}\overline G,\qquad
\kappa=\frac{6300}{1871}.
\end{aligned}
\tag{245}
\]
In particular W^2, both prime columns, and every mixed product
are integrated. No projection of the remaining residual to
finitely many Legendre modes is made.

Let tilde(S)_5 be the actual returned matrix in the unnormalized
coordinates (241). T58 and (240) imply
\[
\boxed{
\begin{aligned}
\widetilde{\mathsf S}_5
 &\succeq\overline A-\kappa\overline G-\Delta I_5=:M,\\
\Delta&\ge
 \epsilon t_z+\kappa\left(2\epsilon\sqrt{t_zt_r}
                                  +\epsilon^2t_z\right).
\end{aligned}
}
\tag{246}
\]
The value of Delta in the packet is one fixed outward decimal
upper bound, hence an exact rational. To prove (246), use
||Z||<=sqrt(t_z), ||bar(R)||<=sqrt(t_r), and
||R-bar(R)||<=epsilon sqrt(t_z). The two Gram matrices differ
in operator norm by at most
2 epsilon sqrt(t_z t_r)+epsilon^2 t_z. Apply T58 to the actual
R=Q_5TZ and bound the remaining source-column difference.
This proves (246) for the full operator. QED.

## BRT-T62. Directed integration and an exact rational positive witness

The accompanying standard-library verifier evaluates (238)--(246)
with 220-digit directed Decimal intervals and exact Fraction
arithmetic. It produces the committed result-220.json. In particular,
\[
\boxed{
 \epsilon<1.159\,10^{-12},\qquad
 \Delta<2.485\,10^{-12},\qquad
 t_z<\frac{17}{10},\qquad t_r<\frac1{250}.
}
\tag{247}
\]
All displayed decimal endpoints in the evidence are exact decimal
rationals. They enclose constants; they are not empirical error bars.

For auditability, the exact rational witness is Q=10^(-18) times
the symmetric integer matrix
\[
\boxed{
\begin{pmatrix}
62347687989954241&61132422349207165&43318504512546330&-3720756156001538&11422107949888666\\
61132422349207165&60227209299108499&43564718753023856&-1353734009149340&13473160425239606\\
43318504512546330&43564718753023856&37571822221109323&17820509294298464&14390988184388402\\
-3720756156001538&-1353734009149340&17820509294298464&59738985484053935&11221692700169019\\
11422107949888666&13473160425239606&14390988184388402&11221692700169019&35613995099809425
\end{pmatrix}.
}
\tag{248}
\]
The exact entry-enclosure errors have maximum absolute row sum
less than 10^(-17). Therefore ||M-Q||<10^(-17), using symmetry.
Exact rational LDL^T elimination of
Q-(10^(-17)+10^(-10))I has five positive pivots. Their complete
fractions are in result-220.json; conservative lower bounds in
order are
\[
 \left(\frac6{100},\frac2{10000},\frac3{1000},
                   \frac1{10000},\frac5{100000}\right).
\tag{249}
\]
They can also be checked directly from the integer matrix (248).
The verifier separately performs interval LDL^T on M-10^(-10)I.
Both sign checks pass. Consequently
\[
\boxed{
 \widetilde{\mathsf S}_5\succeq10^{-10}I_5 .
}
\tag{250}
\]
This numerical margin refers to the unnormalized coefficient
coordinates. If D=diag(sqrt(n_j)), then tilde(S)_5=D S_5 D for
the earlier normalized matrix. Hence its sign agrees with S_5,
and T59 also gives the positive sign of the actual S_3 at this width.

**Constant and rounding enclosure.** Pi is enclosed by Machin's
identity with 110 alternating terms for each arctangent and a
first-omitted-term bound. Euler's constant is enclosed using
H_1024-log(1024)-1/(2*1024), the first 32 Bernoulli corrections,
and the absolute next correction as remainder. The positive-real
digamma expansion has remainder bounded by that first neglected
term; see [DLMF 5.11(i)--(ii)](https://dlmf.nist.gov/5.11).
Bernoulli numbers and all polynomial coefficients are rational.

Arithmetic endpoints use ROUND_FLOOR and ROUND_CEILING.
Logarithms and square roots use correctly rounded Decimal
operations expanded by one adjacent representable number on each
side; see the [Python 3.12 Decimal specification](https://docs.python.org/3.12/library/decimal.html).
Binary float input is rejected by the interval constructor.
There is no numerical quadrature in this verification.

Seven rejection controls cover an indefinite matrix with positive
diagonal, omitted residual debt, a zero pivot, a negative direction,
inflated uncertainty, an interval inverse containing zero, and
binary float input. The packet records all seven passing. These
checks test the verifier's decisions; they do not constitute
external mathematical review or a machine-checked formal proof.

This is a new deterministic computer-assisted certificate.
It does not replay any earlier certificate. Its mathematical
dependencies are the declared source/domain, T57's tail bound,
the exact residual identity, and the analytic bounds proved above.

## BRT-T63. The complete source is positive through L=11/20

For the original real even closed source family,
\[
\boxed{
 0<L\le\frac{11}{20}
 \quad\Longrightarrow\quad
 q_L[g]\ge10^{-11}\|g\|^2
 \quad\text{for every }g\in D_L.
}
\tag{251}
\]
There is no mean-zero, cosh-neutrality, polynomial, finite-mode,
or endpoint-cancellation restriction on g in this conclusion.

**Proof at the endpoint.** Keep L=11/20. Equations (240),
(245), (247) give, with the outward bounds in the packet,
\[
 \|Z\|<\frac43,\qquad
 \|R\|\le\sqrt{t_r}+\epsilon\sqrt{t_z}<\frac8{125},\qquad
 \|C_5^{-1}\|\le\kappa<\frac72 .
\tag{252}
\]
The exact minimizing lift E=Z-C_5^(-1)R therefore obeys
\[
 \|E\|<\frac43+\frac72\frac8{125}
       =\frac{584}{375}<\frac85,\qquad P_5E=(\psi_0,\ldots,\psi_4).
\tag{253}
\]
Every g in D_L has the unique form g=E alpha+z with z in Q_5D_L.
The exact completed-square identity and T57 yield
\[
\begin{aligned}
q_L[g]
 &=\alpha^*\widetilde{\mathsf S}_5\alpha
                    +\|C_5^{1/2}z\|^2\\
 &\ge10^{-10}|\alpha|^2+\frac{1871}{6300}\|z\|^2,\\
\|g\|^2&\le\frac{128}{25}|\alpha|^2+2\|z\|^2.
\end{aligned}
\tag{254}
\]
Since (128/25)10^(-11)<10^(-10) and
2*10^(-11)<1871/6300, this proves (251) at L=11/20.

For any smaller L, zero extension includes its even smooth core
in the endpoint source family. The source formula depends on g,
not on a chosen containing interval: prime terms with no overlap
are automatically zero. Apply the endpoint bound on that core,
then pass to D_L by the inherited form closure. This proves the
entire interval in (251). The earlier stronger margins on their
smaller intervals remain valid. QED.

The positive complement and the finite returned sign have now
both been supplied at 11/20. The result uses the actual prime
weights at t=(1,1); it does not assert positivity for every
diagnostic coupling merely because its complement is positive.
It is a bounded-support result, not all-width positivity or RH.

## BRT-T64. A width-dependent complement comparison with boundary-energy allocation

23 October 2026. Fix 0<L<1 and the actual source (1) with every
prime power n in N_L carrying its weight a_n=Lambda(n)/sqrt(n).
Write u_n=log n, and call n *cap-disjoint* at width L when
L<u_n<2L: its two overlap caps E_(u_n) in (221) are then disjoint.
For such n put
\[
w_n(L)=w_{\min}(L,u_n)=\frac14\log\frac{L^2}{u_n(2L-u_n)}>0,
\qquad
\rho_n(L)=\frac{a_n}{2w_n(L)} .
\tag{255}
\]
Let f_n in [0,1] be chosen for the cap-disjoint n with
\[
\sum_{n\ \text{cap-disjoint}} f_n\rho_n(L)\le1 .
\tag{256}
\]
Define the width constant and the complement constant
\[
\boxed{
\begin{aligned}
\mathcal G_L&=\frac{\Omega_0}2+F(2L)+\frac12\log2
+\frac L2-\frac{L^3}{12}-\frac{3L^2}{64},\\
d_{L,N}&=\mathcal G_L+\frac{H_{2N}}2
-\sum_{n\ \text{cap-disjoint}}(1-f_n)\frac{a_n}2
-\sum_{n\ \text{not cap-disjoint}}a_n .
\end{aligned}
}
\tag{257}
\]
Then for every N>=1 and every h in Q_N D_L,
\[
\boxed{
q_L[h]\ge
\sum_{j\ge N}\Bigl(d_{L,N}+\frac{H_{2j}-H_{2N}}2\Bigr)
|\langle h,\phi_{2j}\rangle|^2
+\ell(h)^2 ,
\qquad
C_N\succeq\mathcal D_{L,N},
}
\tag{258}
\]
where D_(L,N) is the positive diagonal comparison operator with
levels d_(L,j)=d_(L,N)+(H_(2j)-H_(2N))/2, j>=N, on Q_N. The
levels are strictly increasing in j, so the inverse comparison
(217)--(218) holds with these levels: for any admitted column
trial Z with P_N Z=U_N and complete residual R=Q_N T Z,
\[
R^*C_N^{-1}R\preceq
\sum_{j=N}^{N+m-1}\frac{r_{2j}^*r_{2j}}{d_{L,j}}
+\frac1{d_{L,N+m}}\,R_{\ge2(N+m)}^*R_{\ge2(N+m)}
\qquad(m\ge0).
\tag{259}
\]

**Proof.** Every h in Q_N D_L with N>=1 is mean-zero, so T50's
(203)--(204) with nu=L give
\[
q_L[h]\ge\mathcal G_L\|h\|^2+J_L[h]+\mathcal W_L[h]
+\ell(h)^2-\sum_na_nC_h(u_n).
\]
For a cap-disjoint n, (222) gives both
|C_h(u_n)|<=||h||^2/2 and
|C_h(u_n)|<=(2w_n)^(-1) integral over E_(u_n) of w_L|h|^2,
the latter at most W_L[h]/(2w_n). Splitting a_nC_h(u_n) into the
fractions f_n and 1-f_n and using (256),
\[
\sum_na_n|C_h(u_n)|\le\mathcal W_L[h]
+\sum_{\text{cap-disjoint}}(1-f_n)\frac{a_n}2\|h\|^2
+\sum_{\text{other}}a_n\|h\|^2 ,
\]
where the last sum uses only unitarity of each translation. The
boundary energy is therefore spent at most once. ASM (14) gives
J_L[h]>=sum_j H_(2j)|<h,phi_(2j)>|^2/2, and completeness on Q_N
gives (258). Monotonicity of the levels and the inverse-order
argument of T54 give (259). QED.

The greedy allocation that spends the boundary energy on the
cap-disjoint primes in decreasing order of w_n(L) maximizes the
saving in (257) among all allocations satisfying (256), since the
saving per unit of boundary energy spent on n is exactly w_n(L).
The verifier evaluates (255)--(257) with directed arithmetic and
records every quantity in `complement_report`.

Two instances. At L=5/8 both active primes 2 and 3 are
cap-disjoint. The outer cap belongs to prime 3, on which the whole
boundary energy is spent; prime 3 is then charged to the norm only
for the remaining fraction 0.327 of a_3/2, and prime 2 for
the full a_2/2 (its caps are too shallow to use). Then
d_(5/8,5)>0.1390. The packet `l0625_route` certifies the complete
source there with N=5 (T68 below). At L=3/4 the prime powers 3 and 4
are cap-disjoint and prime 2 is not (log 2<3/4). The greedy
allocation spends 0.545 of the boundary energy on prime
power 4 (fully) and the remainder on prime 3, which is charged to the
norm for the fraction 0.913; prime 2 is charged in full. Then
\[
d_{3/4,8}<0,\qquad d_{3/4,14}>0.1123,\qquad
d_{3/4,24}>0.3781 .
\tag{260}
\]
With eight retained modes this rule fails and T53's three-prime
graph comparison (213) is the available route; with fourteen
retained modes both routes are available and they treat the primes
by different arguments. That is the point of carrying both below.

## BRT-T65. The returned matrix at 3/4 is a 10^(-12)-scale sign question

Keep L=3/4 and t=(1,1,1). The packet `near_null` fixes three real
even profiles g_1, g_2, g_3 as rational combinations of
P_0,...,P_22 (denominator 10^12), with g_1 having the coefficient
vector
\[
(497403908704,\,-10^{12},\,805924272529,\,-420265991666,\,
142477108546,\,-26075763174,\,-1112135580,\,2264373717,\,
-759468193,\,181654611,\,-52514923,\,21006257)\cdot10^{-12}
\tag{261}
\]
on P_0, P_2, ..., P_22 of x/L. The verifier encloses
q[g_i]=<g_i,T g_i> two-sidedly with (240)'s operator error:
\[
\boxed{
\frac{q[g_1]}{\|g_1\|^2}\in[3.2049\cdot10^{-12},\,3.2050\cdot10^{-12}],\qquad
\frac{q[g_2]}{\|g_2\|^2}\in[1.8468\cdot10^{-8},\,1.8469\cdot10^{-8}],\qquad
\frac{q[g_3]}{\|g_3\|^2}\in[5.0527\cdot10^{-5},\,5.0528\cdot10^{-5}].
}
\tag{262}
\]
All three values are positive. Since alpha^*S_N alpha is the
minimum of q over the fibre P_N g=alpha, and g_1 has
|P_Ng_1|^2>=0.9977||g_1||^2 for every N>=4 (exact rational
check on (261)), the normalized returned matrix satisfies
\[
\boxed{
\lambda_{\min}(\mathsf S_N(\mathbf1))<3.3\cdot10^{-12}
\qquad(N\ge4),
}
\tag{263}
\]
and already S_3 has a direction below 3.4*10^(-12). The three
displayed quotients fall by four orders of magnitude each; the
24-mode candidate exploration behind the trials indicates a
smallest direction near 2*10^(-15), but only (262)--(263) are
certified.

This locates the difficulty that T46--T55 met. Their complete
inverse allowances 250, 53 and the hierarchy 53, 40/3, 8 multiply
residual Gram traces of order 10^(-1) to 10^(-3) uniformly; no such
bound can resolve a sign of size 10^(-12). A certificate at this
width must therefore produce a residual that is itself below
10^(-12) along the direction of (261). The fixed trials below do
that, not by a smaller allowance but because the minimizing tail
is nearly exact along that direction: on it, T g is nearly a
multiple of g, so the complete residual Q_N T Z nearly vanishes
there while remaining of order 10^(-1) elsewhere. T67 measures
exactly this.

## BRT-T66. Operator and analytic bounds at L=3/4

The actual operator on its polynomial domain is (236) with
\[
\mathcal P_{234}=\frac{a}{2}(S_u+S_{-u})+\frac b2(S_v+S_{-v})
+\frac{a_4}2(S_{w}+S_{-w}),\qquad
w=\log4,\quad a_4=\frac{\log2}2 ,
\tag{264}
\]
and c_L=-(gamma_E+log(3 pi/2))/2. The proof of (236) is unchanged:
(237) holds for every s>0 and the three translations are bounded.
The prime-2 caps overlap at this width; nothing in (236) uses cap
disjointness.

Replace j by its degree-200 Taylor polynomial j_200 from (238) and
h=cosh(x/2) by its twelve even Taylor terms h_24. Cauchy's bound
|j_m|<=5/2^m on |z|=2 and the geometric tail with ratio L give
\[
\boxed{
\sup_{0\le s\le3/2}|j-j_{200}|\le20\Bigl(\frac34\Bigr)^{201},\qquad
\sup_{|x|\le3/4}|h-h_{24}|\le\frac{2(3/8)^{24}}{24!},\qquad
\|T-\overline T\|\le\epsilon:=L\epsilon_j+8L\epsilon_h<1.158\cdot10^{-24} .
}
\tag{265}
\]
The Schur and rank-one arguments of T60 apply verbatim. In the
variable y=(x+L)/(2L) the translated polynomials switch on and off
at the seven-piece partition
\[
0<1-w'<1-v'<u'<1-u'<v'<w'<1,\qquad
u'=\frac u{2L},\ v'=\frac v{2L},\ w'=\frac w{2L},
\tag{266}
\]
whose order the verifier checks with directed arithmetic; the
identities (242)--(244) are then used piece by piece exactly as in
T61.

## BRT-T67. Fixed trials with the level-by-level complement charge

**T53 route.** The packet `t53_route` fixes a rational 24-by-8 array
(denominator 10^14) with identity on the first eight rows; the
columns Z_0,...,Z_7 are combinations of P_0,...,P_46 of x/L. The
complete residual Gram (245) is computed with N=8. The complement
inverse is charged by (218) extended through the hierarchy (213):
explicit levels j=8,...,23 with weights 1/d_j,
d_j=2429/128000+(H_(2j)-H_16)/2, and every level j>=24 with the
single weight 1/d_24<1.80. With Delta as in (246) using
kappa_max=128000/2429,
\[
\widetilde{\mathsf S}_8\succeq\overline A-\mathcal E_{\rm hier}(\overline R)-\Delta I_8=:M,
\qquad
\epsilon<1.158\cdot10^{-24},\quad \Delta<2.360\cdot10^{-23},\quad
t_z<2.401,\quad t_r<0.01214 .
\tag{267}
\]
Here t_r is the complete residual Gram trace. Its charged part,
tr E_hier(bar R)<0.02176, is what enters M; the first sixteen
levels contribute only the rational rounding of the trial.

The directed entry enclosure gives a rational symmetric matrix
Q=10^(-24)-grid rounding of M with ||M-Q||<10^(-22) (maximum
absolute row sum of the enclosure errors). Exact rational LDL^T of
Q-(10^(-22)+10^(-16))I has eight positive pivots, the smallest
exceeding 2.07*10^(-10). The interval LDL^T of M-10^(-16)I also passes.
Consequently
\[
\boxed{
\widetilde{\mathsf S}_8(\mathbf1)\succeq10^{-16}I_8
\qquad\text{(unnormalized coordinates, T53 complement).}
}
\tag{268}
\]
The near-null direction of T65 is retained inside these eight
coordinates, so (268) is a sign statement at the scale (263)
demands; it is not contradicted by (263), which bounds the same
smallest eigenvalue from above.

**Boundary route.** The packet `boundary_route` fixes a 24-by-14
array with identity on the first fourteen rows and uses T64's
levels d_(3/4,j), j>=14, with explicit levels 14,...,23 and tail
weight 1/d_(3/4,24)<2.65. With the same epsilon,
\[
\Delta<8.088\cdot10^{-24},\qquad t_z<2.585,\qquad t_r<0.02370,\qquad
\operatorname{tr}\mathcal E_{\rm hier}(\overline R)<0.06267,
\]
and the same two sign checks pass with the same shift:
\[
\boxed{
\widetilde{\mathsf S}_{14}(\mathbf1)\succeq10^{-16}I_{14}
\qquad\text{(unnormalized coordinates, T64 complement).}
}
\tag{269}
\]
The two routes share the operator, the trial construction and the
residual identity; they differ in the complement comparison, which
is the only place where the prime terms are bounded rather than
integrated. (268) relies on T41--T53's three-prime graph comparison,
(269) on T50, (222) and (14) alone.

## BRT-T68. The complete source is positive through L=3/4

For the original real even closed source family,
\[
\boxed{
 0<L\le\frac34
 \quad\Longrightarrow\quad
 q_L[g]\ge\mu\|g\|^2
 \quad\text{for every }g\in D_L,
 \qquad \mu=9.0\cdot10^{-18} ,
}
\tag{270}
\]
with no mean-zero, polynomial, finite-mode or endpoint restriction
on g. The two certificates give the two values
mu_T53=7.28*10^(-18) and mu_B=9.02*10^(-18); (270) states the larger.

**Proof at the endpoint.** Keep L=3/4 and the T53 route. Every g in
D_L is uniquely g=E alpha+z with E=Z-C_8^(-1)R and z in Q_8D_L, and
\[
q_L[g]=\alpha^*\widetilde{\mathsf S}_8\alpha+\langle z,C_8z\rangle
\ge10^{-16}|\alpha|^2+\frac{2429}{128000}\|z\|^2 .
\]
For the lift norm use C_8^(-1)<=D^(-1) in form order:
||C_8^(-1)R alpha||<=||C_8^(-1/2)|| ||C_8^(-1/2)R alpha||
<=sqrt(kappa_max) sqrt(alpha^*R^*C_8^(-1)R alpha), and
R^*C_8^(-1)R<=E_hier(R), whose trace is at most
tr E_hier(bar R)+kappa_max(2 epsilon sqrt(t_z t_r)+epsilon^2 t_z).
The packet records ||E||<2.63. Since
||g||^2<=2||E||^2|alpha|^2+2||z||^2, the endpoint bound follows with
mu_T53=min(10^(-16)/(2||E||^2), 2429/256000). The boundary route is
identical with d_(3/4,14) in place of 2429/128000 and gives mu_B.
For any smaller L the argument of T63 applies unchanged: the source
formula depends on g alone, prime terms with no overlap vanish
automatically, and zero extension places D_L inside the endpoint
family. QED.

The window result at L=5/8 in the packet `l0625_route` (N=5, T64
levels, margin 10^(-12), all-function bound 1.87*10^(-13)) is a
two-prime instance of the same machinery and is subsumed by (270)
except for its larger constant.

**Dependencies.** The declared source and domain; ASM's Legendre
energy bound (14); T50 (203)--(204) and T56 (222) for T64 and
the boundary route; T53 (213)--(216) and T54 (217)--(218) for the
T53 route; the analytic bounds (265); the exact residual identity;
and the supplied directed-arithmetic implementation. The packet is
reproducible evidence, not independent external certification.

**Closed in this continuation.** The actual returned sign at 3/4
by two complement routes, the all-function estimate (270), the
general-width complement comparison T64, and the quantitative
location (262)--(263) of the sign problem at this width.

**Unresolved.** Widths beyond 3/4: on (3/4, (log 5)/2) the same
three prime powers are active but T53's constants are specific to
3/4 and T64's caps must be re-evaluated; at 2L>log 5 the prime
power 5 activates. All-width positivity/RH remain open. No
percentage of an RH proof follows from this bounded-window advance.

## Scope and continuation

T64--T68 extend the unrestricted working positive interval from
0<L<=11/20 to 0<L<=3/4, with the explicit all-function margin
9*10^(-18) at the new endpoint. The 10^(-11) margin on (0,11/20], the
1/30000 margin on its smaller interval and all earlier stronger local
margins remain valid.

At the new endpoint two complement comparisons are available and
both are used: T53's three-prime graph comparison with eight retained
modes, and T64's width-dependent boundary-energy comparison with
fourteen retained modes. T65 shows that the returned matrix at this
width has a direction below 3.3*10^(-12), so the sign is decided not
by a smaller inverse allowance but by fixed trials whose complete
residual nearly vanishes along that direction. Rigorous operator
remainders, every level of the complement inverse, and all arithmetic
enclosure errors are subtracted before checking the sign; exact
rational LDL^T provides each finite sign witness, and T68 transfers
the sign to every admitted profile and every smaller support width.

The source/domain distinction from native curvature is unchanged.
The retained Legendre counts are not universal observer counts or
MP's separate five-label eta-zero carrier.

This continuation changes this note, its README entry, and adds
proof_lab/rh_l075/ (one verifier, five packet directories with their
configurations, fixed trials and evidence). The L=11/20 packet is
unchanged. The separate Yang--Mills note and earlier certificates are
unchanged. The PR is draft and unmerged, with [skip ci]. New
deterministic candidate exploration (truncated minimizing tails of the
24-mode operator matrix) and interval certification are explicitly
distinguished.

**Closed in the working derivation and accompanying computation:**
the actual returned sign at 3/4 by two routes, the all-function
estimate (270), the general-width comparison T64, and the enclosure
(262)--(263).

**Unresolved:** widths in (3/4, (log 5)/2) keep the same three prime
powers but need the comparison constants re-evaluated (T64 applies
with recomputed caps; T53's constants are specific to 3/4); beyond
(log 5)/2 the prime power 5 activates; and all-width positivity/RH.
No percentage of an RH proof follows from this bounded-window
advance.

**Dependencies:** the declared source/domain and inherited
Gamma/Legendre and Schur results, T50, T53--T54, T56 (222), the
analytic estimates above, and the supplied directed-arithmetic
implementation. The packet is reproducible evidence, not independent
external certification.
