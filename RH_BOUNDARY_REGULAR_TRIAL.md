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
complete certificate is NOT established. There is no new unrestricted positive interval,
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
The left side has not been bounded by this threshold in this note.
The positive trial reading in (14) is a completed part of the
calculation; the residual subtraction is the part still open.

## Scope and continuation

The established outputs are an explicit boundary-selected correction,
its exact overlap polynomial, the strict prime penalty bound (7),
a bounded full source column, the positive trial bound (14), and the
stronger mode-resolved residual certificate (21).
No new positive source interval is claimed. DWE's previously proved
unrestricted interval remains unchanged.

The trial is an ordinary source profile. Its quadratic shape is not
identified with EMK connection curvature, a physical field, or a
universal cut law. Its usefulness is justified by (6) and (10).

The source-specific remaining task is a signed enclosure of (21), or
a better admitted trial if that endpoint is negative. Boundary
regularity and a small prime term cannot replace that calculation.
Even a first-prime certificate would leave the all-width RH task open.

Only this note and its README entry change. Validation consists of
written polynomial and logarithmic integration, rational comparisons,
source polarization, boundary limits, inverse form order, and
immutable-content readback. No scientific computation, numerical
source evaluation, independent agent, certificate replay or workflow
dispatch is performed. The commit uses [skip ci]. The PR is a draft
stacked on #137 and remains unmerged.
