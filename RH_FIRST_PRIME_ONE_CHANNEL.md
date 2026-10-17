# OPC-1 — One retained source channel across the first prime

17 October 2026. Written continuation of MCW-1 at
749538e5e6c53bb815215bc73c9b3d4adf022345.

This step reaches an active-prime window while releasing one of
the two moment restrictions used in FCP. For the actual source,
\[
\boxed{
0<L\le\frac7{20},\quad g\in D_L,\quad \int g=0
\quad\Longrightarrow\quad
q_L[g]\ge\frac{18863}{300000}\|g\|_2^2+\ell(g)^2
\ge\frac1{16}\|g\|_2^2+\ell(g)^2.
}
\]
The cosh moment \(\ell(g)\) is arbitrary. Prime 2 is active at
the new endpoint, and its actual signed correlation is retained
before being absorbed by the logarithmic boundary term.

This also proves that the complement of the constant observer is
coercive throughout this window. Consequently the **full** source
has an exact one-scalar Schur reduction with all unresolved memory
retained. The sign of that scalar is not established here.
MCW's full-family positive interval remains \(L\le1/4\).

The actual pole can be returned exactly in this reduction. Its
denominator correction is nonnegative and less than \(1/4000\).
The remaining numerator and arithmetic memory still have to be
controlled; a small denominator correction alone does not prove RH.

## 1. Source, domain and the one mean-zero condition

Retain the full real even Hilbert space \(H_L=L^2_{\rm even}(-L,L)\),
zero extension outside \(I=(-L,L)\), and the source form domain
\(D_L\) of PTE/ASM/GCC. Write
\[
\begin{gathered}
N_g=\|g\|_2^2,\qquad C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,\\
\ell(g)=\int_Ig(x)\cosh(x/2)\,dx,\qquad
m_j(g)=\int_Ix^jg(x)\,dx,\\
\mathcal N_L=\{n\ge2:n\text{ is a prime power},\ \log n<2L\},
\qquad a_n=\Lambda(n)/\sqrt n .
\end{gathered}
\]
The actual form is
\[
\boxed{
q_L[g]=\ell(g)^2+\frac{\Omega_0}{2}N_g
+\int_0^\infty k(u)(N_g-C_g(u))\,du
-\sum_{n\in\mathcal N_L}a_nC_g(\log n),
}
\tag{1}
\]
with
\[
\begin{gathered}
\Omega_0=-\gamma_E-\pi/2-3\log2-\log\pi,\qquad
k(u)=e^{u/2}/(2\sinh u),\\
F(v)=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}),\\
w_L(x)=-\tfrac14\log(1-(x/L)^2),\qquad
J_L[g]=\tfrac18\int_I\int_I
 \frac{|g(x)-g(y)|^2}{|x-y|}\,dy\,dx .
\end{gathered}
\tag{2}
\]
Put
\[
\phi_0=(2L)^{-1/2},\qquad P_0g=\langle g,\phi_0\rangle\phi_0,
\qquad Q_0=I-P_0,\qquad \mathcal Z_L=Q_0D_L.
\tag{3}
\]
Thus \(\mathcal Z_L=\{g\in D_L:m_0(g)=0\}\).
It is larger than FCP's family with both \(m_0=0\) and \(\ell=0\).

ASM-T2 establishes that zero-extended polynomials, including
\(\phi_0\), belong to the source operator domain and form a core.
In particular P_0 and Q_0 preserve D_L. The actual operator
associated with (1) is denoted \(\mathcal T_L\).

MCW-T1 gives, for \(L\le9/16\),
\[
\begin{aligned}
q_L[g]\ge{}&
 B_LN_g+J_L[g]+\int_Iw_L|g|^2+\ell(g)^2\\
&-\frac3{25}m_0^2+\frac1{25}m_0m_2
-\sum_{n\in\mathcal N_L}a_nC_g(\log n),\\
B_L={}&\frac{\Omega_0}{2}+F(2L)+\frac12\log2
       +\frac{12L}{25}-\frac{8L^3}{75}.
\end{aligned}
\tag{4}
\]
The moment terms are retained for the full family. They vanish
on \(\mathcal Z_L\) because its mean is zero, not because its
cosh or quadratic moment is assumed zero.

## OPC-T1. The endpoint tail has a written rational lower bound

At \(L_*=7/20\),
\[
\boxed{F(7/10)>\frac{149}{100}.}
\tag{5}
\]

**Proof.** Let \(x=7/40\) and \(z=\tanh x\). Retain PRW's
polynomial bound
\[
\tanh x<x-x^3/3+2x^5/15 .
\]
At this x, \(x^2/3-2x^4/15>1/100\), so
\[
z<\frac{99}{100}x=\frac{693}{4000}.
\tag{6}
\]
The cosh series gives
\(\cosh x\le(1-x^2/2)^{-1}=3200/3151\), while
\(\sinh x>x\). Hence
\[
z>\frac{22057}{128000}>\frac{43}{250},\qquad z<x.
\tag{7}
\]

PRW/FCP retain
\[
\log2>\frac{6931}{10000},\quad
\log(5/4)>\frac{2231}{10000},\quad
\frac{\pi}{4}>\frac{7853}{10000}.
\]
Therefore
\[
\log5>\frac{16093}{10000},\qquad
\log(8/7)=2\operatorname{artanh}(1/15)>
\frac{1352}{10125}.
\]
Since \(40/7=5(8/7)\), (6) and \(-\log(1-t)>t\) yield
\[
-\frac12\log z>
\frac{16093}{20000}+\frac{676}{10125}+\frac1{200}.
\tag{8}
\]
The alternating arctangent bound gives
\[
\arctan z
<\frac{693}{4000}-\frac{169}{100000}+\frac1{30000}.
\tag{9}
\]
Here \((43/250)^3/3>169/100000\) and
\((7/40)^5/5<1/30000\).

Use the exact identity
\(F(7/10)=-\tfrac12\log z+\pi/4-\arctan z\).
Equations (8)--(9) give
\[
F(7/10)>
\frac{149}{100}+\frac{989}{8100000}>
\frac{149}{100}.
\]
All comparisons are finite rational inequalities or displayed
series bounds. No numerical source evaluation is used. QED.

## OPC-T2. The first prime is absorbed by the actual boundary surplus

At \(L_*=7/20\), the only active prime power is 2.
Put \(u=\log2\), \(a=(\log2)/\sqrt2\), and
\[
\mathcal C=\{x\in(-L_*,L_*):|x|>u-L_*\}.
\]
For every real L2 profile on this interval,
\[
\boxed{
a|C_g(u)|\le\frac a2\int_{\mathcal C}|g|^2.
}
\tag{10}
\]
For \(s=x/L_*\), define the bounded rational weight
\[
v_*(x)=\frac{s^2}{4-2s^2}.
\]
Then
\[
\boxed{
w_{L_*}(x)-\frac a2\mathbf1_{\mathcal C}(x)
\ge v_*(x)\ge0
\quad\text{almost everywhere}.
}
\tag{11}
\]

**Proof.** The retained bounds
\[
\frac{6931}{10000}<\log2<\frac{1733}{2500}<\frac7{10},
\qquad \log3>1
\]
show that exactly prime 2 is active and that \(L_*<u<2L_*\).
Each point has at most one neighbor at displacement u.
Integrating \(2|g(x)g(x-u)|\le|g(x)|^2+|g(x-u)|^2\)
over the overlap proves (10). The two endpoint caps have total
length \(2(2L_*-u)\); neither their size nor their sign is ignored.

Off the caps, MCW-T4 gives \(w_{L_*}\ge v_*\).
On the caps,
\[
\frac{|x|}{L_*}>
\frac{3431}{3500}>\frac{49}{50},
\qquad
1-(49/50)^2=\frac{99}{2500}<\frac1{25}.
\]
Thus
\[
w_{L_*}(x)>\frac14\log25>\frac34.
\]
The last bound follows from \(\log5>16093/10000>3/2\).
Also \(\log2<7/10\) and \(\sqrt2>7/5\) imply \(a/2<1/4\).
Consequently on the caps
\[
w_{L_*}-a/2>1/2\ge v_*.
\]
This proves (11). QED.

Combining (4), (10) and (11) gives the **full-family** lower bound
\[
\boxed{
q_{L_*}[g]\ge B_{L_*}N_g+J_{L_*}[g]
+\int_Iv_*|g|^2+\ell(g)^2
-\frac3{25}m_0^2+\frac1{25}m_0m_2.
}
\tag{12}
\]
The prime has been bounded by proved boundary energy. It has not
been removed from the source definition or declared positive.

## OPC-T3. One mean-zero condition suffices through the first prime

For all \(0<L\le7/20\) and all \(g\in\mathcal Z_L\),
\[
\boxed{
q_L[g]-\ell(g)^2
\ge\delta_*\|g\|_2^2,\qquad
\delta_*=\frac{18863}{300000}>\frac1{16}.
}
\tag{13}
\]
At the endpoint the stronger estimate retains the boundary weight:
\[
q_{L_*}[g]-\ell(g)^2
\ge\delta_*N_g+\int_Iv_*|g|^2 .
\tag{14}
\]

**Proof.** On the smooth mean-zero even core, the two residual
terms in (12) vanish. The Legendre energy has no constant mode,
so \(J_{L_*}[g]\ge3N_g/4\). Retain FCP's scalar bound
\(\Omega_0/2>-26871/10000\), OPC-T1 and the lower log2 bound:
\[
\begin{aligned}
B_{L_*}+\frac34
>{}&-\frac{26871}{10000}+\frac{149}{100}
+\frac{6931}{20000}+\frac{21}{125}
-\frac{343}{75000}+\frac34\\
={}&\frac{18863}{300000}>\frac1{16}.
\end{aligned}
\tag{15}
\]
This proves (14) on that core. No bound for \(\ell(g)\) is used.

A smaller supported test extends by zero to the endpoint without
changing its mean, norm, pole reading or actual source value.
This proves (13) on every smaller smooth mean-zero core.

For completeness, the smooth mean-zero core is dense in
\(\mathcal Z_L\). Approximate g in the source form norm by
smooth even g_n, then subtract \(m_0(g_n)\rho\), where rho is
one fixed smooth even test supported in I with integral one.
The correction tends to zero in form norm. Source form closure
and continuity of \(\ell\) extend (13) to \(\mathcal Z_L\).
The endpoint weight is bounded, so (14) extends as well. QED.

The advance is the removal of the cosh-neutrality condition in
this active-prime window. It does not remove the mean-zero
condition from the positivity statement (13).

## OPC-T4. The full source has exactly one retained sign condition

Fix \(0<L\le7/20\). Let C_L be the positive self-adjoint operator
of the restricted form \(q_L|_{\mathcal Z_L}\), acting on \(Q_0H_L\).
By (13),
\[
C_L\ge\delta_*I,\qquad
\|C_L^{-1}\|\le\delta_*^{-1}<16.
\tag{16}
\]
Define from the actual source
\[
\begin{gathered}
\alpha_L=q_L[\phi_0],\qquad
b_L=Q_0\mathcal T_L\phi_0,\qquad
\rho_L=C_L^{-1}b_L,\\
\boxed{s_L=\alpha_L-\langle b_L,C_L^{-1}b_L\rangle.}
\end{gathered}
\tag{17}
\]
For \(g=c\phi_0+h\), \(h\in\mathcal Z_L\),
\[
\boxed{
q_L[g]=s_Lc^2+
\|C_L^{1/2}(h+c\rho_L)\|_2^2.
}
\tag{18}
\]
In particular,
\[
\boxed{
q_L\ge0\text{ on }D_L\quad\Longleftrightarrow\quad s_L\ge0.
}
\tag{19}
\]

**Proof.** Domain preservation and closedness are retained from
ASM. The new estimate (13) supplies coercivity already for
N=1, with an explicit inverse bound in an active-prime window.
ASM-T2 gives \(\mathcal T_L\phi_0\in H_L\).
Polarization yields
\[
q_L[c\phi_0+h]=\alpha_Lc^2+
2c\langle b_L,h\rangle+\|C_L^{1/2}h\|^2.
\]
Complete the square using \(C_L\rho_L=b_L\) to obtain (18).
Taking \(h=-c\rho_L\) proves the reverse implication in (19). QED.

No infinite-dimensional information is lost in the scalar:
its returned-memory term contains the inverse on the complete
mean-zero space. The Schur identity itself is ASM-T3 specialized
to N=1; the new input is the first-prime coercivity in OPC-T3.

For clarity, the constant observer's actual column is explicit:
\[
\begin{aligned}
(\mathcal T_L\phi_0)(x)
={}&\ell(\phi_0)\cosh(x/2)+\frac{\Omega_0}{2}\phi_0\\
&+\frac{\phi_0}{2}\bigl(F(L-x)+F(L+x)\bigr)\\
&-\frac{\phi_0}{2}\sum_{n\in\mathcal N_L}a_n
 \bigl(\mathbf1_I(x-\log n)+\mathbf1_I(x+\log n)\bigr).
\end{aligned}
\tag{20}
\]
It is in L2 because its boundary singularities are logarithmic.
Equation (20), not a freely chosen positive matrix, defines b_L.

There are also precise spectral consequences:

- The negative index of the full form is zero if \(s_L\ge0\),
  and exactly one if \(s_L<0\).
- If \(s_L=0\), its kernel is exactly
  \(\operatorname{span}\{\phi_0-\rho_L\}\).
- If \(s_L>0\), a full-family lower bound is
  \[
  q_L[g]\ge
  \min\left\{\frac{s_L}{1+2\|\rho_L\|^2},\frac1{32}\right\}
  \|g\|^2.
  \tag{21}
  \]
- At most one eigenvalue of the source operator lies below \(1/16\).
  This follows from the codimension-one estimate and the retained
  compact resolvent.

For (21), put \(z=h+c\rho_L\), use
\(\|g\|^2\le(1+2\|\rho_L\|^2)c^2+2\|z\|^2\), and apply (18).
The index and kernel assertions follow from the invertible change
of variables \((c,h)\mapsto(c,h+c\rho_L)\).
These are statements about this finite-window source operator;
they do not count zeta zeros or assert an off-line zero.

## OPC-T5. The pole returns as an exact positive quotient

Separate only the explicitly known pole square:
\[
\mathfrak a_L[g]=q_L[g]-\ell(g)^2.
\]
Let A_L be the self-adjoint operator of
\(\mathfrak a_L|_{\mathcal Z_L}\). By (13),
\(A_L\ge\delta_*I\). Put
\[
\begin{gathered}
p=\ell(\phi_0),\qquad
u=Q_0\bigl(\cosh(x/2)\bigr),\qquad
\alpha_0=\mathfrak a_L[\phi_0],\\
v=Q_0\bigl(\mathcal T_L\phi_0-p\cosh(x/2)\bigr),\\
\mu=\langle v,A_L^{-1}v\rangle,\quad
\nu=\langle u,A_L^{-1}v\rangle,\quad
\kappa=\langle u,A_L^{-1}u\rangle .
\end{gathered}
\tag{22}
\]
These are actual-source quantities. They satisfy
\(\mu,\kappa\ge0\) and \(\nu^2\le\mu\kappa\).
The exact returned scalar is
\[
\boxed{
s_L=\alpha_0-\mu+\frac{(p-\nu)^2}{1+\kappa}.
}
\tag{23}
\]
Furthermore, throughout the stated interval,
\[
\boxed{0\le\kappa<\frac1{4000}.}
\tag{24}
\]

**Proof of the identity.** On the mean-zero space,
\[
C_L=A_L+u\otimes u,\qquad b_L=v+pu,\qquad
\alpha_L=\alpha_0+p^2.
\]
The rank-one inverse identity is
\[
C_L^{-1}
=A_L^{-1}-
\frac{(A_L^{-1}u)\otimes(A_L^{-1}u)}{1+\kappa}.
\tag{25}
\]
Its denominator is positive. Direct multiplication by
\(A_L+u\otimes u\) proves (25) on its operator domain.
Substitute in (17) and expand:
\[
\alpha_0+p^2-
\left(\mu+2p\nu+p^2\kappa-
\frac{(\nu+p\kappa)^2}{1+\kappa}\right)
=\alpha_0-\mu+\frac{(p-\nu)^2}{1+\kappa}.
\]
Cauchy--Schwarz for the positive inverse A_L^{-1} gives
\(\nu^2\le\mu\kappa\).

**Proof of the denominator bound.** The mean-zero cosh profile obeys
\[
\boxed{
\|u\|^2\le\frac{L^5}{360}\cosh^2(L/2).
}
\tag{26}
\]
Indeed, as a function of z=x^2, cosh(x/2) has derivative at most
\(\cosh(L/2)/8\) on \(0\le z\le L^2\).
The pairwise variance identity and this Lipschitz bound give
\[
\begin{aligned}
\|u\|^2
&=\frac1{4L}\int_I\int_I
  |\cosh(x/2)-\cosh(y/2)|^2\,dx\,dy\\
&\le\frac{\cosh^2(L/2)}{64}
  \int_I(x^2-L^2/3)^2\,dx
=\frac{L^5}{360}\cosh^2(L/2).
\end{aligned}
\]
Since \(A_L^{-1}\le16I\), it remains to bound the right side
at \(L_*=7/20\). The elementary comparisons
\[
L_*^5=\frac{16807}{3200000}<\frac{17}{3200},
\qquad
\cosh^2(7/40)\le(3200/3151)^2<\frac{21}{20}
\]
give
\[
\kappa<
\frac{16}{360}\frac{17}{3200}\frac{21}{20}
=\frac{119}{480000}<\frac1{4000}.
\]
Monotonicity of \(L^5\cosh^2(L/2)\) extends the bound to smaller L.
This proves (24). QED.

Thus the full sign condition can also be written
\[
\boxed{
(p-\nu)^2\ge(1+\kappa)(\mu-\alpha_0).
}
\tag{27}
\]
If the right side is nonpositive, the condition is automatic.
Otherwise its sign requires the actual arithmetic return
\(\mu,\nu\). A positive denominator, the pole square, or the small
bound (24) does not by itself establish (27).

This quotient is an exact response formula for the declared
source and cut. No probability interpretation or universal
physical constant is assigned to it.

## OPC-T6. The bare constant observer is positive at the new endpoint

At \(L_*=7/20\),
\[
\boxed{\alpha_{L_*}=q_{L_*}[\phi_0]>\frac{17}{200}.}
\tag{28}
\]
This does not determine the sign of s_(L_*).

**Proof.** For the normalized constant profile,
\[
C_{\phi_0}(t)=
\begin{cases}
1-t/(2L),&0\le t<2L,\\
0,&t\ge2L,
\end{cases}
\qquad
p^2=2L\left(\frac{\sinh(L/2)}{L/2}\right)^2.
\]
Consequently its exact value at the endpoint is
\[
\alpha_{L_*}
=p^2+\frac{\Omega_0}{2}+F(2L_*)
+\frac1{2L_*}\int_0^{2L_*}t\,k(t)\,dt
-a\left(1-\frac{u}{2L_*}\right).
\tag{29}
\]
The retained quadratic kernel bound and the positive sinh series
give
\[
\frac1{2L}\int_0^{2L}t\,k(t)\,dt
\ge\frac12+\frac{6L}{25}-\frac{2L^3}{25},
\qquad p^2\ge2L+\frac{L^3}{6}.
\]
Also \(a<1/2\) and \(u>6931/10000\), so
\[
a\left(1-\frac{u}{2L_*}\right)<\frac{69}{14000}.
\]
Using OPC-T1 and the retained scalar bound yields
\[
\begin{aligned}
\alpha_{L_*}>{}&
-\frac{26871}{10000}+\frac{149}{100}+\frac12
+\frac{21}{250}-\frac{343}{100000}
+\frac7{10}+\frac{343}{48000}-\frac{69}{14000}\\
={}&\frac{17}{200}+\frac{5773}{8400000}>\frac{17}{200}.
\end{aligned}
\]
The last step is a rational comparison. QED.

The constant reading and the full mean-zero block are both
positive. Their coupling still matters:
\[
s_{L_*}=\alpha_{L_*}-\langle b_{L_*},C_{L_*}^{-1}b_{L_*}\rangle .
\]
A positive bare observer value cannot replace this returned value.

## OPC-T7. The returned potential is invariant under a change of lift

The scalar has the variational meaning
\[
\boxed{
s_L=\min_{\substack{g\in D_L\\\langle g,\phi_0\rangle=1}}q_L[g],
\qquad g_*=\phi_0-\rho_L.
}
\tag{30}
\]
The minimizer is unique. Choosing a different representative
of the same unit-mean fibre leaves both s_L and g_* unchanged.

More explicitly, take \(h_0\in\operatorname{Dom}C_L\) and
\(e=\phi_0+h_0\). Then
\[
\begin{gathered}
\alpha'=\alpha_L+2\langle b_L,h_0\rangle+
              \langle C_Lh_0,h_0\rangle,\\
b'=b_L+C_Lh_0,\qquad
\rho'=C_L^{-1}b'=\rho_L+h_0,\\
\boxed{\alpha'-\langle b',C_L^{-1}b'\rangle=s_L,
\qquad e-\rho'=\phi_0-\rho_L.}
\end{gathered}
\tag{31}
\]

**Proof.** Fix c=1 in (18). Coercivity of C_L makes the minimum
unique at \(h=-\rho_L\), proving (30). Expanding the primed
quadratic expressions proves (31). QED.

This is a stated invariance under changing the lift while
preserving the mean functional and the actual source. A
different support width, a different observation functional,
or a different arithmetic source is not identified by (31).

## OPC-T8. One scalar has an outward residual enclosure

For any trial \(h\in\operatorname{Dom}C_L\), set
\[
r_h=b_L-C_Lh,\qquad
\widehat s_L(h)=
\alpha_L-2\langle b_L,h\rangle+\langle C_Lh,h\rangle
=q_L[\phi_0-h].
\]
Then
\[
\boxed{
\widehat s_L(h)-16\|r_h\|^2
\le s_L\le\widehat s_L(h).
}
\tag{32}
\]

**Proof.** Since \(r_h=C_L(\rho_L-h)\), (18) gives
\[
\widehat s_L(h)-s_L
=\langle r_h,C_L^{-1}r_h\rangle .
\]
Use \(0<C_L^{-1}\le16I\). QED.

A nonnegative lower endpoint in (32) certifies full positivity
at that L. A strictly positive endpoint gives a positive margin
through (21). A negative upper endpoint supplies an actual
negative source test. An interval containing zero does neither.

This is ASM's residual enclosure with a proved one-channel
inverse bound in the new window. No trial, residual norm, or
finite certificate is evaluated here. In particular, convergence
of unweighted operator residuals is not assumed merely from
form-core density; GCC retains the appropriate bounded-cut
certificate mechanism.

## What changed and what is still open

| Family or condition | Status after this note |
|---|---|
| Full even source, no neutral moments | MCW's \(L\le1/4\) positive interval is unchanged |
| Mean-zero even source, arbitrary cosh moment | Positive through \(L\le7/20\), with margin \(1/16\) even after subtracting the pole square |
| Mean-zero and cosh-zero source | FCP's larger \(L\le111/200\) interval is unchanged |
| Actual prime 2 at \(L=7/20\) | Its signed correlation is dominated by a proved boundary surplus |
| Complete source sign through \(7/20\) | Exactly one retained scalar \(s_L\), including all unresolved memory |
| Pole return denominator | \(1\le1+\kappa<4001/4000\) |
| Sign of \(s_{7/20}\) | Not established |
| All-width positivity, FWT growth condition, RH | Still open |

The next concrete task is to bound the actual scalar (17), or
the equivalent arithmetic-return inequality (27), with retained
errors. The full-family gate is not passed by testing only
\(\alpha_L\), by deleting its coupling, or by identifying a
positive mean-zero restriction with the full source.

The new result reduces the required observer aperture to one
through an active-prime window. It does not reduce the entire
infinite source to a freely specified scalar or prove a global
one-mode obstruction. No actual off-line zero, RH closure,
universal probability floor, or physical law is claimed.

## References and verification scope

- [MCW-1](RH_MIXED_CHANNEL_WINDOW.md): full source, quadratic moment residual, boundary weight and the separate full-family quarter window.
- [PRW-1](RH_POLE_RECOVERY_WINDOW.md): elementary constants, tanh estimate and the earlier pole-recovery comparison.
- [FCP-1](RH_FINITE_CUT_POSITIVITY.md): improved Gamma constants and the two-moment restricted interval.
- [TPW-1](RH_TWO_PRIME_WINDOW.md): polynomial-kernel identity and retained logarithm bounds.
- [BPC-1](RH_BOUNDARY_PRIME_CAPS.md): translated-pair estimates and boundary caps.
- [ASM-1](RH_ARITHMETIC_SEAM_REDUCTION.md): actual operator columns, domain-preserving polynomial cuts, exact Schur memory and residual enclosures.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): the closed actual source, energy representation and compact resolvent.
- [GCC-1](RH_GAMMA_CUT_CERTIFICATES.md): bounded Gamma cuts and finite-certificate completeness with errors retained.

The Schur complement, rank-one inverse, variance estimate and
Cauchy--Schwarz are standard tools. The new source-specific
content is the first-prime boundary absorption, the coercive
mean-zero interval without cosh neutrality, and the resulting
quantified one-channel and pole-return statements. No
literature-wide novelty or optimal endpoint is asserted.

Only this note and its README entry change. Verification consists
of written derivation review, document checks, immutable-head
content and PR metadata. No scientific computation, numerical xi
evaluation, independent agent, source replay, finite certificate
evaluation, or workflow dispatch is performed. Frozen sources,
engines, certificates and workflows remain unchanged.
The documentation commit includes [skip ci]. The PR is stacked
on draft #126 and remains unmerged.
