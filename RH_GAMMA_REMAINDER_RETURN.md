# GRR-1 — Global Gamma remainders and quantitative source return

17 October 2026. Written continuation of OPC-1 at
90b60ed0a775a8d10e2f0bd9fa8e901b040b15e2.

An exact positive remainder in the Gamma kernel improves the
complete arithmetic source. On the already established full
quarter window it gives
\[
\boxed{
q_L[g]\ge\frac1{850}\|g\|_2^2
\qquad(0<L\le1/4,\ g\in D_L).
}
\]
No neutral moment is imposed. This strengthens the guaranteed
margin; it does not enlarge the unrestricted support interval.

At the first-prime endpoint \(L_*=7/20\), the same remainder
gives a quantitative surplus above the older MCW comparison.
A sufficient condition for full source positivity becomes
\[
\boxed{\sigma_{\rm old}\ge-1/1600,}
\]
where \(\sigma_{\rm old}\) is the exact returned scalar of that
specified comparison, including its complete mean-zero inverse.
An explicit outward residual enclosure for this condition is
given below. Its premise is not established in this note.

The kernel estimate holds for every positive displacement.
It also removes the \(e^{-L}\) loss from the Legendre coefficient
in ASM's all-width positive-tail prescription. The sign of the
remaining finite Schur block is still required. Neither this
improvement nor the local margins prove RH.

## Source, domain and notation

Use the complete source and retained closed domain of
[MCW-1](RH_MIXED_CHANNEL_WINDOW.md) and
[OPC-1](RH_FIRST_PRIME_ONE_CHANNEL.md).
On \(I=(-L,L)\), functions are real and even, and are extended
by zero outside I. Write
\[
\begin{gathered}
N_g=\|g\|_2^2,\qquad
C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,\\
m_j(g)=\int_Ix^jg(x)\,dx,\qquad
\ell(g)=\int_Ig(x)\cosh(x/2)\,dx,\\
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
k(u)=\frac{e^{u/2}}{2\sinh u},\qquad
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}),\\
\Omega_0=-\gamma_E-\pi/2-3\log2-\log\pi,\qquad
\frac{\Omega_0}{2}>-\frac{26871}{10000}>-\frac{27}{10}.
\end{gathered}
\tag{2}
\]
The form domain \(D_L\) is the closure of real even smooth
compact tests in the Gamma form norm. Put
\[
J_L[g]=\frac18\int_I\int_I
\frac{|g(x)-g(y)|^2}{|x-y|}\,dy\,dx,\qquad
w_L(x)=-\frac14\log(1-(x/L)^2).
\tag{3}
\]
The retained logarithmic decomposition is
\[
\int_0^{2L}\frac{N_g-C_g(u)}{2u}\,du
=J_L[g]+\frac12\log2\,N_g+\int_Iw_L|g|^2.
\tag{4}
\]
Also \(C_g(u)=0\) for \(u\ge2L\) and
\(0\le N_g-C_g(u)\le2N_g\). Initially one can work on the smooth
core; the bounded differences below preserve the closed domain.

The positive parameter a introduced next is a kernel-comparison
parameter. It is not the native source anchor, a prime, a bank
depth, or the diagnostic coupling lambda.

## GRR-T1. A global kernel bound has an explicit positive remainder

For every \(u>0\),
\[
\boxed{
k(u)=\frac1{2u}+\frac14-\frac{u}{48}-\frac{u^2}{32}
+R_0(u),\qquad R_0(u)>0.
}
\tag{5}
\]
Consequently, for every fixed \(a>0\), define
\[
A_a=\frac14-\frac a{96},\qquad
B_a=\frac1{32}+\frac1{96a},\qquad
P_a(u)=A_a-B_au^2 .
\]
Then
\[
\boxed{
k(u)=\frac1{2u}+P_a(u)+\Delta_a(u),\qquad
\Delta_a(u)=\frac{(u-a)^2}{96a}+R_0(u)>0.
}
\tag{6}
\]

**Proof.** Set
\[
E(u)=e^{u/2}-
\left(1+\frac u2-\frac{u^2}{24}-\frac{u^3}{16}\right)
\frac{\sinh u}{u},\qquad
R_0(u)=\frac{E(u)}{2\sinh u}.
\]
The coefficients of E through degree three vanish. For \(m\ge2\),
the coefficients of \(u^{2m}\) and \(u^{2m+1}\), respectively, are
\[
\frac{(2m+1)/2^{2m}-1+(2m)(2m+1)/24}{(2m+1)!},
\qquad
\frac{2^{-(2m+1)}-1/2+(2m)(2m+1)/16}{(2m+1)!}.
\tag{7}
\]
Both are positive. In the even case, the numerator is \(7/48\)
at \(m=2\), and the quadratic term alone exceeds one for \(m\ge3\).
In the odd case, \((2m)(2m+1)/16\ge5/4>1/2\).
These are convergent entire series, so \(E(u)>0\) for \(u>0\).
This proves (5). Subtract \(P_a\) from the nonsingular part
of (5); the difference is
\((u-a)^2/(96a)\). This proves (6). QED.

In particular, the fixed choice \(a=1/2\) gives
\[
\boxed{P_*(u)=\frac{47}{192}-\frac5{96}u^2.}
\tag{8}
\]
It is not asserted that \(P_*\) pointwise dominates the older
polynomial \(P_{\rm old}(u)=6/25-u^2/25\). The complete forms
are compared explicitly in GRR-T3.

## GRR-T2. The remainder survives the complete source and its cut

For finite \(L>0\) and fixed \(a>0\), define the comparison form
\[
\begin{aligned}
\mathfrak b_{a,L}[g]={}&
\left(\frac{\Omega_0}{2}+F(2L)+\frac12\log2
+2A_aL-\frac83B_aL^3\right)N_g\\
&+J_L[g]+\int_Iw_L|g|^2+\ell(g)^2
-\frac{A_a}{2}m_0(g)^2+B_am_0(g)m_2(g)\\
&-\sum_{n\in\mathcal N_L}a_nC_g(\log n).
\end{aligned}
\tag{9}
\]
The exact identity is
\[
\boxed{
q_L[g]=\mathfrak b_{a,L}[g]
+\int_0^{2L}\Delta_a(u)(N_g-C_g(u))\,du.
}
\tag{10}
\]
Moreover,
\[
\boxed{
q_L[g]\ge\mathfrak b_{a,L}[g]+\varepsilon(L,a)N_g,\qquad
\varepsilon(L,a)=
\frac{(2L-a)^3-(L-a)^3}{576a}>0.
}
\tag{11}
\]

**Proof of the identity and domain statement.** Evenness gives
the exact moment relations
\[
\int_0^{2L}C_g(u)\,du=\frac12m_0^2,\qquad
\int_0^{2L}u^2C_g(u)\,du=m_0m_2 .
\]
Insert (6) on \(0<u<2L\), use (4), and keep the tail
\(\int_{2L}^\infty k(u)N_g\,du=F(2L)N_g\) unchanged.
This gives (9)--(10). The polynomial comparison is never
integrated to infinity.

Near zero \(R_0(u)=O(u^3)\), and on every finite interval
\(\Delta_a\) is integrable. The difference in (10) is a bounded
positive quadratic form, bounded above by
\(2N_g\int_0^{2L}\Delta_a(u)\,du\).
Thus \(\mathfrak b_{a,L}\) is a closed semibounded form on the
same \(D_L\), obtained by a bounded form perturbation of \(q_L\).
Identity (10) extends to that domain.

**Proof of the uniform surplus.** For \(L\le u\le2L\), the two
overlap caps \((u-L,L)\) and \((-L,L-u)\) are disjoint up to
endpoints. Applying \(2|rs|\le r^2+s^2\) to each translated pair
therefore gives
\[
|C_g(u)|\le\frac12N_g,\qquad N_g-C_g(u)\ge\frac12N_g.
\]
Retain only the square term in \(\Delta_a\) on that interval:
\[
\int_0^{2L}\Delta_a(u)(N_g-C_g(u))\,du
\ge\frac{N_g}{192a}\int_L^{2L}(u-a)^2du.
\]
This is (11), and the integral is strictly positive. QED.

For the fixed choice \(a=1/2\),
\[
\boxed{
\varepsilon(1/4,1/2)=\frac1{18432},\qquad
\varepsilon(7/20,1/2)=\frac{91}{2304000}.
}
\tag{12}
\]

There is also a precise profile-dependent interpretation.
For nonzero g, put \(f(u)=N_g-C_g(u)\) on \((0,2L)\) and
\[
X=\int f=2LN_g-\frac12m_0^2,\quad
Y=\int u^2f=\frac83L^3N_g-m_0m_2,\quad Z=\int uf .
\tag{13}
\]
The cap bound above implies \(X,Y>0\). Cauchy--Schwarz gives
\(Z\le\sqrt{XY}\). Let \(\mathcal H_L[g]\) denote (9) with
all \(A_a,B_a\) terms omitted. Equation (5) gives the identity
\[
\begin{aligned}
q_L[g]={}&\mathcal H_L[g]+\frac X4-\frac Y{32}
-\frac{\sqrt{XY}}{48}\\
&+\frac{\sqrt{XY}-Z}{48}+\int_0^{2L}R_0(u)f(u)\,du .
\end{aligned}
\tag{14}
\]
The first line is the best bound within the family (9) for
this particular g, attained at \(a_g=\sqrt{Y/X}\).
Indeed the a-dependent loss is \((aX+Y/a)/96\).

The normalized density \(f(u)/X\) is a mathematically specified
probability distribution on displacement for this fixed test.
Its contribution in the second line includes
\(X(\sqrt{\mathbb E[u^2]}-\mathbb E[u])/48\).
This does not assign a universal probability floor, diagnose
irrationality, or identify a physical probability law.
A g-dependent choice \(a_g\) also does not define a single
linear source operator. The Schur comparisons below use fixed a.

## GRR-T3. The fixed comparison has a strict two-channel gain

Let \(\mathfrak b_{{\rm old},L}\) be the right side of MCW (7),
with polynomial \(6/25-u^2/25\), keeping its pole, boundary,
moments and every active prime term. Let
\(\mathfrak b_{*,L}=\mathfrak b_{1/2,L}\).
For \(0<L\le9/16\), the former is MCW's proved lower comparison.

Set
\[
\Delta A=\frac{23}{4800},\qquad
\Delta B=\frac{29}{2400},\qquad
d_L=2L\Delta A-\frac83L^3\Delta B .
\]
Their exact difference, as forms, is
\[
\mathfrak b_{*,L}[g]-\mathfrak b_{{\rm old},L}[g]
=d_LN_g-\frac{\Delta A}{2}m_0^2+\Delta Bm_0m_2.
\tag{15}
\]
Use degree-indexed coordinates
\[
\phi_0=\frac1{\sqrt{2L}},\quad
\phi_2=\sqrt{\frac5{2L}}P_2(x/L),\quad
c_j=\langle g,\phi_j\rangle\ (j=0,2).
\]
The difference equals
\[
d_L\|g-c_0\phi_0-c_2\phi_2\|^2+
(c_0,c_2)M_L(c_0,c_2)^{\mathsf T},
\]
where
\[
\boxed{
M_L=
\begin{pmatrix}
L\Delta A-2L^3\Delta B&2L^3\Delta B/(3\sqrt5)\\
2L^3\Delta B/(3\sqrt5)&d_L
\end{pmatrix}.
}
\tag{16}
\]
Consequently,
\[
\boxed{
\begin{aligned}
\mathfrak b_{*,1/4}[g]&\ge
\mathfrak b_{{\rm old},1/4}[g]+\frac1{1250}N_g,\\
\mathfrak b_{*,7/20}[g]&\ge
\mathfrak b_{{\rm old},7/20}[g]+\frac3{5000}N_g .
\end{aligned}
}
\tag{17}
\]

**Proof.** The retained moment coordinates are
\[
m_0=\sqrt{2L}\,c_0,\qquad
m_2=L^2\sqrt{2L}\left(\frac{c_0}{3}+\frac{2c_2}{3\sqrt5}\right).
\]
Substitution in (15) proves (16).

At \(L=1/4\), the matrix entries are
\[
M_{00}=\frac{21}{25600},\quad
M_{22}=\frac{109}{57600},\quad
M_{02}=\frac{29}{230400\sqrt5}.
\]
The diagonal entries of \(M_L-I/1250\) are
\(13/640000\) and \(1573/1440000\), and its determinant is
\[
\frac{36\cdot13\cdot1573-125\cdot29^2}
 {36\cdot640000\cdot1440000}
=\frac{631039}{36\cdot640000\cdot1440000}>0.
\]
Also \(d_L=M_{22}>1/1250\).

At \(L=7/20\), the entries are
\[
M_{00}=\frac{2051}{3200000},\quad
M_{22}=\frac{14203}{7200000},\quad
M_{02}=\frac{9947}{28800000\sqrt5}.
\]
The diagonal entries of \(M_L-3I/5000\) are
\(131/3200000\) and \(9883/7200000\). Its determinant is
\[
\frac{180\cdot131\cdot9883-9947^2}{5(28800000)^2}
=\frac{134098331}{5(28800000)^2}>0.
\]
Again \(d_L=M_{22}>3/5000\).
The positive diagonal and determinant tests prove (17). QED.

The second estimate also has a scale-sensitive extension:
\[
\mathfrak b_{*,L}[g]-\mathfrak b_{{\rm old},L}[g]
\ge\frac{3L}{1750}N_g\qquad(0<L\le7/20).
\tag{18}
\]
Indeed \(M_L/L\) decreases in the matrix order as \(L^2\)
increases. Its derivative with respect to \(L^2\) is
\[
\Delta B
\begin{pmatrix}-2&2/(3\sqrt5)\\2/(3\sqrt5)&-8/3\end{pmatrix},
\]
which is negative definite; the unscaled determinant is \(236/45\).
The same monotonicity holds for \(d_L/L\). Apply the endpoint
bound (17) and multiply by L. A constant endpoint gain is not
asserted uniformly as \(L\) tends to zero.

## GRR-T4. The full quarter-window margin is more than tripled

For every \(0<L\le1/4\) and \(g\in D_L\),
\[
\boxed{
q_L[g]\ge
\left(\frac{667}{2000000}+\frac1{1250}+\frac1{18432}\right)N_g
\ge\frac1{850}N_g .
}
\tag{19}
\]
The last comparison is strict for nonzero g.

**Proof.** The proof of MCW-T5 bounds the old comparison itself:
\[
\mathfrak b_{{\rm old},1/4}[g]\ge
\frac{667}{2000000}N_g.
\]
To see that this is not an unsupported strengthening of MCW,
its proof starts with the right side of MCW (7), bounds the
retained pole and mixed moment, and then applies the boundary
weight and weighted Cauchy--Schwarz. Every subsequent lower
estimate therefore applies to this same comparison form.

Apply (11), (12) and the first line of (17) to obtain the first
bound in (19) at the endpoint. For the second, the first two
fractions sum to \(2267/2000000\), and
\[
\frac{2267}{2000000}-\frac1{850}
=-\frac{1461}{34000000},\qquad
34000000-1461\cdot18432=7070848>0.
\]
Thus the added \(1/18432\) is more than the remaining difference.

Zero extension preserves the actual source and norm of every
smaller supported smooth test. The endpoint bound consequently
holds on all smaller cores. Closure extends it to every \(D_L\)
in the stated interval. QED.

This result imposes no condition on \(m_0,m_2,\ell\).
The interval remains \(L\le1/4\); the first prime threshold
\((\log2)/2\) is outside it.

The compact formulation also improves. In PTE's notation,
\[
q_L[g]=\|G_L^{1/2}g\|^2-c_LN_g,\qquad K_L=c_LG_L^{-1}.
\]
No prime is active in this interval, and \(0<c_L=-\Omega_0/2<27/10\).
The spectral theorem applied to (19) gives
\[
\boxed{
G_L\ge(c_L+1/850)I,\quad
\|K_L\|\le\frac{2295}{2296},\quad
I-\tau(\lambda)K_L\ge\frac1{2296}I,
}
\tag{20}
\]
where \(\tau(\lambda)=\lambda/(1+\lambda)\), including
\(\tau(\infty)=1\). The ratio is
\((27/10)/(27/10+1/850)=2295/2296\).

GCC's finite-certificate completeness still applies on this
same interval. No specific finite certificate is evaluated.
PRW's stronger margin \(1/300\) on \(L\le1/5\) remains available.

## GRR-T5. The first-prime mean-zero inverse bound improves

Let
\(\mathcal Z_L=\{g\in D_L:m_0(g)=0\}\).
For every \(0<L\le7/20\),
\[
\boxed{
q_L[g]-\ell(g)^2
\ge\delta_{\rm new}N_g,\qquad
\delta_{\rm new}=\frac{93383}{1440000}>\frac2{31},
\quad g\in\mathcal Z_L.
}
\tag{21}
\]

**Proof.** At \(L_*=7/20\), OPC-T3's boundary and Legendre proof
gives
\[
\mathfrak b_{{\rm old},L_*}[g]-\ell(g)^2
\ge\frac{18863}{300000}N_g
\quad(g\in\mathcal Z_{L_*}).
\]
The proof uses this comparison before any positive source
remainder is discarded. When \(m_0=0\), (15) is exactly
\(d_{L_*}N_g\), and \(d_{L_*}=14203/7200000\).
Equation (10) therefore gives
\[
\frac{18863}{300000}+\frac{14203}{7200000}
=\frac{93383}{1440000}.
\]
The additional remainder in (11) is not needed for this bound.
The comparison with \(2/31\) follows from
\(93383\cdot31-2\cdot1440000=14873>0\).

Zero extension preserves the mean, pole, norm and actual source.
The mean-zero smooth core is dense by OPC-T3's fixed integral-one
correction argument. This transfers and closes the estimate
on each smaller interval. QED.

In OPC's notation, let \(C_L\) be the mean-zero operator of
\(q_L\), and \(A_L\) that of \(q_L-\ell^2\). Both satisfy
\[
\boxed{\|C_L^{-1}\|,\ \|A_L^{-1}\|<\frac{31}{2}.}
\tag{22}
\]
Their exact source Schur identities are unchanged. In particular,
\[
s_L=\alpha_0-\mu+\frac{(p-\nu)^2}{1+\kappa}
\]
still contains the full mean-zero return.

The denominator estimate also strengthens:
\[
\boxed{0\le\kappa<\frac1{4160}.}
\tag{23}
\]
Indeed OPC-T5 proves
\[
\|Q_0\cosh(x/2)\|^2\le
\frac{L^5}{360}\cosh^2(L/2)<\frac{119}{7680000}
\quad(0<L\le7/20).
\]
Combine this with (22):
\[
\kappa<\frac{3689}{15360000}<\frac1{4160},
\qquad3689\cdot4160=15346240<15360000.
\]
Here \(Q_0\) is the orthogonal projection onto the mean-zero
space; it is not the source form.

For an actual-source trial \(h\in\operatorname{Dom}C_L\),
OPC's outward enclosure becomes
\[
\boxed{
q_L[\phi_0-h]-\frac{31}{2}\|b_L-C_Lh\|^2
\le s_L\le q_L[\phi_0-h],
}
\tag{24}
\]
where \(b_L=Q_0\mathcal T_L\phi_0\) is the actual source column.
This follows from the same exact error
\(\langle b_L-C_Lh,C_L^{-1}(b_L-C_Lh)\rangle\).
No trial or sign of the actual returned scalar is evaluated here.

## GRR-T6. A signed comparison tolerance suffices for first-prime positivity

In this section fix \(L=L_*=7/20\).
Define the complete returned comparison scalars
\[
\sigma_{\rm old}=
\min_{\substack{g\in D_{L_*}\\\langle g,\phi_0\rangle=1}}
 \mathfrak b_{{\rm old},L_*}[g],
\qquad
\sigma_*=
\min_{\substack{g\in D_{L_*}\\\langle g,\phi_0\rangle=1}}
 \mathfrak b_{*,L_*}[g].
\tag{25}
\]
Both minima exist uniquely. Their mean-zero restrictions are
coercive by the old bound used in GRR-T5 and by (15).
Their constant columns belong to L2: these comparison forms
differ from the actual form by bounded forms, so ASM's
operator-domain assertion for \(\phi_0\) still applies.
Thus the usual exact Schur construction is valid.

OPC's actual returned scalar \(s_{L_*}\) obeys
\[
\boxed{
s_{L_*}\ge\sigma_*+\frac{91}{2304000}
\ge\sigma_{\rm old}+\eta_*,
\qquad
\eta_*=\frac3{5000}+\frac{91}{2304000}>\frac1{1600}.
}
\tag{26}
\]
In particular, the following implication is proved:
\[
\boxed{
\sigma_{\rm old}\ge-\frac1{1600}
\quad\Longrightarrow\quad
q_{L_*}[g]\ge\frac{167}{11520000}\|g\|_2^2
\quad(g\in D_{L_*}).
}
\tag{27}
\]
**The premise of (27) remains unproved.**

**Proof.** On the unit-mean fibre,
\(N_g\ge|\langle g,\phi_0\rangle|^2=1\).
Apply the uniform source gain (11)--(12) and the second line
of (17), then take the respective infima. This gives (26).

For the full norm bound in (27), write \(g=c\phi_0+h\),
\(h\in\mathcal Z_{L_*}\). The old comparison's exact Schur
square gives
\[
\mathfrak b_{{\rm old},L_*}[g]
=\sigma_{\rm old}c^2+
\|C_{\rm old}^{1/2}(h+c\rho_{\rm old})\|^2.
\]
Here \(C_{\rm old}\) is its complete mean-zero operator and
\(\rho_{\rm old}=C_{\rm old}^{-1}b_{\rm old}\).
Under the premise,
\(\mathfrak b_{{\rm old},L_*}[g]\ge-N_g/1600\).
Meanwhile (11) and (17) give the form inequality
\(q_{L_*}\ge\mathfrak b_{{\rm old},L_*}+\eta_*I\).
Finally,
\[
\eta_*-\frac1{1600}
=\frac{91}{2304000}-\frac1{40000}
=\frac{167}{11520000}>0.
\]
This proves (27). QED.

There is a concrete residual condition for its unproved premise.
Let \(\mathcal T_{\rm old}\) be the comparison's self-adjoint
operator, and put
\[
\alpha_{\rm old}=\mathfrak b_{{\rm old},L_*}[\phi_0],\qquad
b_{\rm old}=Q_0\mathcal T_{\rm old}\phi_0 .
\]
Then
\[
\sigma_{\rm old}=\alpha_{\rm old}
-\langle b_{\rm old},C_{\rm old}^{-1}b_{\rm old}\rangle,
\qquad \|C_{\rm old}^{-1}\|<16.
\]
The inverse bound follows from
\(18863/300000>1/16\), with the pole square retained.
For \(h\in\operatorname{Dom}C_{\rm old}\), set
\[
\widehat\sigma_{\rm old}(h)
=\mathfrak b_{{\rm old},L_*}[\phi_0-h],\qquad
r_{\rm old}(h)=b_{\rm old}-C_{\rm old}h.
\]
Completing the square gives
\[
\boxed{
\widehat\sigma_{\rm old}(h)-16\|r_{\rm old}(h)\|^2
\le\sigma_{\rm old}\le\widehat\sigma_{\rm old}(h).
}
\tag{28}
\]
Therefore the sufficient certificate is
\[
\boxed{
\widehat\sigma_{\rm old}(h)-16\|r_{\rm old}(h)\|^2
\ge-\frac1{1600}.
}
\tag{29}
\]
A proof of (29) for one admitted trial would give (27) and,
by zero extension, the same source margin for every smaller L.
No such trial certificate has been established in this note.

The allowed negative comparison tolerance in (29) comes from
a proved positive difference between the comparison and the
actual source. It is not permission to discard error or to
declare a negative actual returned scalar acceptable.
Equations (25)--(29) retain the whole mean-zero inverse.

## GRR-T7. Every width has a stronger positive-tail prescription

The global estimate (6) permits a fixed parameter chosen from
L alone:
\[
a_L=\frac{2L}{\sqrt3}.
\]
For any mean-zero profile, the bulk coefficient in (9) is
maximized within this family by that choice, and equals
\[
2A_{a_L}L-\frac83B_{a_L}L^3
=\frac L2-\frac{L^2}{12\sqrt3}-\frac{L^3}{12}.
\tag{30}
\]
Indeed for \(m_0=0\), (13) gives
\(X=2LN_g\), \(Y=8L^3N_g/3\), so \(\sqrt{Y/X}=2L/\sqrt3\).
This parameter does not depend on the profile.

Define
\[
\mathcal B_L=
\frac{\Omega_0}{2}+F(2L)+\frac12\log2
+\frac L2-\frac{L^2}{12\sqrt3}-\frac{L^3}{12},
\qquad S_L=\sum_{n\in\mathcal N_L}a_n .
\]
Use ASM's mode-count indexing
\[
\psi_j(x)=\sqrt{\frac{4j+1}{2L}}P_{2j}(x/L),\qquad
V_N=\operatorname{span}\{\psi_0,\ldots,\psi_{N-1}\},
\quad N\ge1 .
\]
Thus \(\psi_0=\phi_0\) and \(\psi_1=\phi_2\); the two index
conventions are explicit.

For every \(L>0\), and every \(h\in D_L\) orthogonal to \(V_N\),
\[
\boxed{
q_L[h]\ge d^{\rm new}_{L,N}\|h\|^2,\qquad
d^{\rm new}_{L,N}=
\mathcal B_L+\frac12H_{2N}-S_L,
\quad H_m=\sum_{j=1}^m\frac1j .
}
\tag{31}
\]
In particular, a sufficient finite aperture for margin at least
one is
\[
\boxed{
2N\ge\exp\left(2[1+(S_L-\mathcal B_L)_+]\right).
}
\tag{32}
\]

**Proof.** Since \(N\ge1\), orthogonality implies \(m_0(h)=0\).
Use (9)--(10) with the fixed \(a_L\). The ordinary moment
residual vanishes. Drop the nonnegative pole, boundary weight,
and remainder, and bound each prime correlation in absolute
value by \(\|h\|^2\). The retained Legendre energy satisfies
\(J_L[h]\ge H_{2N}\|h\|^2/2\).
This proves (31). Finally
\(H_{2N}\ge\log(2N+1)\), and substitution of (32) proves
\(d^{\rm new}_{L,N}\ge1\). QED.

There is no factor \(e^{-L}\) in front of \(H_{2N}/2\).
ASM's previous sufficient tail margin was
\[
d^{\rm ASM}_{L,N}
=F(L)-c_L+\frac{e^{-L}}2H_{2N},\qquad
c_L=S_L-\Omega_0/2 .
\]
Both are lower bounds; their maximum may be used. No claim
that the new expression dominates at every fixed \(L,N\) is needed.

For a comparison of explicit sufficient sizes without prime
asymptotics, the elementary estimate is
\[
S_L\le\sum_{2\le n<e^{2L}}\frac{\log n}{\sqrt n}
\le4Le^L .
\tag{33}
\]
This follows from \(\log n<2L\) and
\(\sum_{n=2}^{M}n^{-1/2}\le\int_1^M x^{-1/2}dx<2\sqrt M\).
The empty sum causes no exception. Using (2), \(F(2L)>0\),
and dropping other positive terms from \(\mathcal B_L\), a
fully elementary sufficient version of (32) is
\[
\boxed{
2N\ge
\exp\left(8Le^L+\frac{L^2}{6\sqrt3}
+\frac{L^3}{6}+\frac{37}{5}\right).
}
\tag{34}
\]
The same estimate in ASM's earlier prescription gave the
sufficient condition
\[
2N\ge\exp\left(8Le^{2L}+\frac{37}{5}e^L\right).
\]
Thus the displayed sufficient thresholds for \(\log(2N)\)
grow as \(O(Le^L)\) and \(O(Le^{2L})\), respectively, as
\(L\) tends to infinity. These are upper prescriptions, not
estimates of the actual minimum aperture.

The finite Schur block, its exact coupling, and all returned
memory must still be retained. A positive unresolved tail for
each width does not prove positivity of that finite block.
This result supplies no evaluated bank, full-source sign at
arbitrary width, or native RBP/UCB bank-depth extension.

## What is established and what remains

| Item | Result in this continuation |
|---|---|
| Gamma comparison range | Every \(u>0\), with an exact positive remainder |
| Full even source on \(L\le1/4\) | Margin at least \(1/850\), with the stronger exact sum (19) |
| Full unrestricted support interval | Still \(L\le1/4\) |
| Mean-zero source through \(7/20\) | Pole-subtracted margin \(93383/1440000>2/31\) |
| First-prime returned source | Strict quantitative surplus above a specified comparison |
| Sufficient endpoint certificate | The outward lower bound (29) need only reach \(-1/1600\) |
| Truth of that certificate | Not established |
| All-width positive unresolved tails | Improved explicit sufficient size (32)--(34) |
| Finite Schur signs at arbitrary width / RH | Still open |

The displacement distribution after (14) is defined from an
actual test and source correlation. No claim that irrationality
forces a universal nonzero error is used. Likewise, positive
kernel remainders do not remove the complete source's scalar
subtraction or signed prime correlations.

PRW's stronger smaller-window bound, FCP's separate two-moment
interval \(L\le111/200\), and native finite-bank ranges remain
available and unchanged. FWT's full-wave growth condition and
the all-width arithmetic positivity condition remain unproved.

## References and verification scope

- [OPC-1](RH_FIRST_PRIME_ONE_CHANNEL.md): first-prime boundary absorption, mean-zero coercivity, exact one-channel and pole returns, variance and residual bounds.
- [MCW-1](RH_MIXED_CHANNEL_WINDOW.md): the old polynomial comparison, mixed moments, and quarter-window proof.
- [TPW-1](RH_TWO_PRIME_WINDOW.md): polynomial-kernel moment identity and logarithmic boundary decomposition.
- [PRW-1](RH_POLE_RECOVERY_WINDOW.md) and [FCP-1](RH_FINITE_CUT_POSITIVITY.md): retained scalar constants, previous full and restricted windows.
- [ASM-1](RH_ARITHMETIC_SEAM_REDUCTION.md): polynomial domains, Legendre energy, positive tails and exact finite Schur memory.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): the actual closed source and compact energy representation.
- [GCC-1](RH_GAMMA_CUT_CERTIFICATES.md): finite-certificate completeness with errors retained.

Power-series comparison, weighted Cauchy--Schwarz, bounded form
perturbation and Schur completion are standard tools. The
source-specific content here is the global remainder identity,
its quantitative complete-form gains, the improved local margins
and the resulting signed certificate tolerance and tail-size
prescription. No literature-wide novelty or optimality is claimed.

Only this note and its README entry change. Verification consists
of written derivation review, document checks, immutable-head
content and PR metadata. No scientific computation, numerical xi
evaluation, independent agent, source replay, finite certificate
evaluation, or workflow dispatch is performed. Frozen sources,
engines, certificates and workflows remain unchanged.
The documentation commit includes [skip ci]. The PR is stacked
on draft #127 and remains unmerged.
