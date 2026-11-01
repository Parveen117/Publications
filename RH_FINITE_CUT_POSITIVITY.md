# FCP-1 — Finite-cut error and a wider positive prime window

17 October 2026. Written continuation at Publications commit
9a2ba9de1854942f5d18b8aa379d66c5ee5db9cc, following
[PCM-1](RH_PRIME_CUT_MEMORY.md) and
[ASM-1](RH_ARITHMETIC_SEAM_REDUCTION.md).

The finite prime record has a strictly positive discarded probability.
That fact does not prevent an exact positivity proof: the corresponding
source error can be bounded below an independently proved positive
margin. This note supplies both a general cut-and-memory enclosure
and an explicit application to the actual arithmetic source.

The concrete extension is
\[
\boxed{
0<L\le\frac{111}{200}
\quad\Longrightarrow\quad
q_L[g]\ge\frac{49}{245500}\|g\|_2^2
>\frac1{6000}\|g\|_2^2
\quad(0\ne g\in\mathcal W_L).
}
\]
Here \(\mathcal W_L\) is exactly the real even family with zero
constant and \(\cosh(x/2)\) moments used in FPW, BPC and TPW.
The endpoint increases from TPW's \(11/20\) to \(111/200\).
No new moment condition is imposed. TPW's stronger margin \(1/500\)
remains available on its earlier range.

At the new endpoint only primes 2 and 3 contribute. Their record cuts
at depth \(K=35\), retaining records \(0,\ldots,35\), have positive
discarded probabilities \(2^{-18}\) and \(3^{-18}\). Nevertheless their
combined source-error bound is less than \(1/12000\), so both the actual
source and its finite-record approximation are strictly positive on
this same family.

**Status: RH remains unproved.** The proved window is restricted to
\(\mathcal W_L\). The general finite-matrix enclosure below applies to
the full even source space, but it has not been evaluated to prove a
larger unrestricted window. All-width positivity remains the RH
obligation. The arguments are written analytic proofs; no numerical
source evaluation or scientific execution is used.

## 1. The retained source and the two different cuts

Let \(I=(-L,L)\). On the real Hilbert space \(H_L\) of even
\(L^2(I)\) functions, extend profiles by zero and retain
\[
\begin{gathered}
N_g=\|g\|_2^2,\qquad
T_ug(x)=g(x-u),\qquad C_g(u)=\langle g,T_ug\rangle,\\
\ell(g)=\int_I g(x)\cosh(x/2)\,dx,\qquad
a_n=\frac{\Lambda(n)}{\sqrt n},\\
\mathcal N_L=\{n\ge2:n\text{ is a prime power},\ \log n<2L\},\\
q_L[g]=\ell(g)^2+\frac{\Omega_0}{2}N_g+
\mathcal E_\Gamma[g]-\sum_{n\in\mathcal N_L}a_nC_g(\log n),\\
\Omega_0=-\gamma_E-\frac{\pi}{2}-3\log2-\log\pi,\qquad
\mathcal E_\Gamma[g]=\frac12\int_0^\infty k(u)\|g-T_ug\|_2^2\,du,\\
k(u)=\frac{e^{u/2}}{2\sinh u},\qquad
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}).
\end{gathered}
\tag{1}
\]
Use the form domain \(D_L\) and self-adjoint operator
\(\mathcal T_L\) constructed in PTE and ASM. The scalar
\(\Omega_0/2\), every active prime power and the pole square are retained.

The fixed restricted family is
\[
\mathcal W_L=
\left\{g\in C_c^\infty(I,\mathbb R):
g(-x)=g(x),\
\int g=0,\
\ell(g)=0\right\}.
\tag{2}
\]
On this family write
\(\mathscr A[g]=\Omega_0N_g/2+\mathcal E_\Gamma[g]\), so that
\(q_L[g]=\mathscr A[g]-\sum a_nC_g(\log n)\).

There are two distinct cuts in the following argument.

- A prime record cut truncates the geometric packet of PCM.
  Its omitted probability produces a bounded error in the source.
- A Legendre cut retains finitely many profile modes.
  Its unresolved source memory is the exact Schur correction of ASM.

Neither error is discarded, and neither is identified with the other.

## FCP-T1. A finite prime record gives a uniform source-error bound

For each active prime \(p<e^{2L}\), put
\[
r_p=p^{-1/2},\qquad
\kappa_p=\frac{1+r_p}{1-r_p},\qquad U_p=T_{\log p},
\]
where \(U_p\) is unitary on the full line. On the zero extension of g set
\[
\begin{aligned}
V_p&=\sqrt{1-r_p^2}\sum_{k=0}^\infty r_p^k U_p^k,\\
V_{p,K_p}&=\sqrt{1-r_p^2}\sum_{k=0}^{K_p}r_p^k U_p^k,\qquad
\epsilon_p=r_p^{K_p+1}.
\end{aligned}
\tag{3}
\]
The record weights are \(w_{p,k}=(1-r_p)r_p^k\), so
\[
\epsilon_p=\sum_{k>K_p}w_{p,k}>0
\]
at every finite depth. No renormalization after truncation is made.

Define the finite-record source by
\[
\widetilde q_{L,\mathbf K}[g]=
\ell(g)^2+\frac{\Omega_0}{2}N_g+\mathcal E_\Gamma[g]
+\frac12\sum_{p<e^{2L}}\log p
\left(N_g-\|V_{p,K_p}g\|_2^2\right).
\tag{4}
\]
All sums over p in this note are over primes. Then
\[
\boxed{
|\widetilde q_{L,\mathbf K}[g]-q_L[g]|
\le\eta_{L,\mathbf K}N_g,\qquad
\eta_{L,\mathbf K}=
\sum_{p<e^{2L}}\log p\,\kappa_p\epsilon_p .
}
\tag{5}
\]
The difference of the associated operators is bounded on \(H_L\),
with norm at most \(\eta_{L,\mathbf K}\). They have the same operator
domain and the same form domain.

**Proof.** The PTE/PCM resolvent identity gives
\[
\|V_pg\|^2=N_g+2\sum_{k\ge1}r_p^k C_g(k\log p).
\]
The correlations vanish when \(k\log p\ge2L\). Thus replacing each
\(V_{p,K_p}\) in (4) by \(V_p\) gives exactly (1), including all prime
powers.

The geometric norm bounds are
\[
\|V_p\|\le\sqrt{\kappa_p},\quad
\|V_{p,K_p}\|\le\sqrt{\kappa_p}(1-\epsilon_p),\quad
\|V_p-V_{p,K_p}\|\le\sqrt{\kappa_p}\epsilon_p.
\]
Consequently
\[
\left|\|V_pg\|^2-\|V_{p,K_p}g\|^2\right|
\le 2\kappa_p\epsilon_pN_g.
\]
Multiply by \(\log p/2\) and sum to obtain (5).
Compression of the bounded full-line differences to \(H_L\)
gives the bounded self-adjoint operator difference. The standard
bounded-perturbation domain identity gives the last assertion. QED.

For fixed L the prime set is finite and \(\eta_{L,\mathbf K}\)
tends to zero as all record depths increase. This is compatible
with \(\epsilon_p>0\) at every finite depth. It provides no
positive, depth-independent lower floor for the omitted probability.

## FCP-T2. Prime-cut error and unresolved source memory have one enclosure

Retain ASM's even orthonormal Legendre modes
\[
\phi_j(x)=\sqrt{\frac{4j+1}{2L}}P_{2j}(x/L),
\]
their synthesis \(U_N:\mathbb R^N\to H_L\), and
\(P_N=U_NU_N^*\), \(Q_N=I-P_N\). Put
\[
c_L=\sum_{n\in\mathcal N_L}a_n-\frac{\Omega_0}{2},\qquad
d=d_{L,N}=F(L)-c_L+\frac{e^{-L}}2H_{2N}.
\tag{6}
\]
Choose N with \(d>0\), as ASM-T1 permits for every finite L.
The tail operator C associated with the restriction of \(q_L\)
to \(Q_ND_L\) satisfies \(C\ge dI\).

Let
\[
A=U_N^*\mathcal T_LU_N,\qquad B=Q_N\mathcal T_LU_N,\qquad
S=A-B^*C^{-1}B
\tag{7}
\]
be the exact source-memory matrix. For any finite-column map
\[
Z:\mathbb R^N\longrightarrow\operatorname{Dom}(\mathcal T_L),
\qquad P_NZ=U_N,
\]
define
\[
G_Z=Z^*Z,\qquad
\widetilde A=Z^*\widetilde{\mathcal T}_{L,\mathbf K}Z,\qquad
\widetilde R=Q_N\widetilde{\mathcal T}_{L,\mathbf K}Z,\qquad
\eta=\eta_{L,\mathbf K}.
\tag{8}
\]
For every \(\theta>0\) one has the finite-matrix enclosure
\[
\boxed{
\begin{aligned}
S\succeq{}&
\widetilde A-
\left(\eta+\frac{(1+\theta^{-1})\eta^2}{d}\right)G_Z
-\frac{1+\theta}{d}\widetilde R^*\widetilde R,\\
S\preceq{}&\widetilde A+\eta G_Z .
\end{aligned}
}
\tag{9}
\]
If the left matrix in (9) is positive semidefinite, the actual
full even source \(q_L\) is nonnegative on \(D_L\).

**Proof.** Write
\(A_Z=Z^*\mathcal T_LZ\) and \(R=Q_N\mathcal T_LZ\).
Completing the positive tail square in
\(q_L[Za+h]\), \(h\in Q_ND_L\), gives the exact identity
\[
S=A_Z-R^*C^{-1}R.
\tag{10}
\]
Indeed \(P_N(Za+h)=U_Na\), and its infimum is the same
observation-fibre potential as in ASM-T4.

By (5),
\[
\widetilde A-\eta G_Z\preceq A_Z
\preceq\widetilde A+\eta G_Z.
\]
Also \(R=\widetilde R+
Q_N(\mathcal T_L-\widetilde{\mathcal T}_{L,\mathbf K})Z\).
For every a, the weighted two-term square inequality gives
\[
\|Ra\|^2\le
(1+\theta)\|\widetilde Ra\|^2+
(1+\theta^{-1})\eta^2\|Za\|^2.
\]
Use \(0<C^{-1}\preceq d^{-1}I\) in (10) to prove (9).
ASM's exact square decomposition gives
\(q_L\ge0\) if and only if \(S\succeq0\), proving the implication. QED.

There is a finite construction of useful trial columns.
If \(\eta<d\), choose \(M>N\) and minimize the finite-record form
over the intermediate Legendre modes \(N,\ldots,M-1\) for each
fixed first-N reading. Its intermediate block is at least
\((d-\eta)I\), hence invertible. The resulting Z has polynomial
columns, \(P_NZ=U_N\), and
\[
\widetilde R=Q_M\widetilde{\mathcal T}_{L,\mathbf K}Z
\]
by the finite normal equations. ASM-T2 ensures that these columns
are in the operator domain.

Formula (9) is a sufficient certificate, not an evaluated
certificate for a new unrestricted window. Its matrices are finite,
but their source integrals and residual Gram entries still require
rigorous outward bounds. In particular, a positive trial matrix
without both error terms does not establish the sign of S.
No residual-norm convergence rate is inferred from form-core density.

## FCP-T3. Improved elementary Gamma bounds

The bounds needed for the wider window are
\[
\boxed{
\frac{\Omega_0}{2}>-\frac{26871}{10000},\qquad
k(u)\ge\frac1{2u}+\frac6{25}-\frac{u^2}{25}
\quad(0<u\le9/8),\qquad
F(111/100)>\frac{21121}{18000}.
}
\tag{11}
\]

**The scalar Gamma bound.** First,
\[
\gamma_E<
H_n-\log n-\frac1{2n}+\frac1{12n^2}.
\tag{12}
\]
For a direct proof, telescope
\(\log(1+1/k)-1/(k+1)\) from \(k=n\) to infinity and use
the Laplace integral for \(1/x\). This gives
\[
H_n-\log n-\gamma_E
=\int_0^\infty e^{-nt}
\left(\frac1t-\frac1{e^t-1}\right)\,dt.
\]
The positive series
\[
(1+x^2/3)\sinh x-x\cosh x
=\sum_{j\ge2}\frac{4j(j-1)}{3(2j+1)!}x^{2j+1}>0
\]
implies
\(1/(e^t-1)<1/t-1/2+t/12\).
Integration proves (12).

At \(n=4\), (12) gives
\(\gamma_E<377/192-2\log2\).
Retain the bounds proved in TPW-T3,
\[
\frac{6931}{10000}<\log2<\frac{1733}{2500},\qquad
\log\pi<\frac{573}{500},\qquad \frac{\pi}{2}<\frac{11}{7}.
\tag{13}
\]
Then
\[
\begin{aligned}
\gamma_E+\pi/2+3\log2+\log\pi
&<\frac{377}{192}+\frac{11}{7}
+\frac{1733}{2500}+\frac{573}{500}\\
&=\frac{4514303}{840000}
<\frac{26871}{5000}.
\end{aligned}
\]
This proves the first bound in (11).

**The wider kernel interval.** Retain TPW's exact polynomial
identity, with
\(E_4(u)=1+u/2+u^2/8+u^3/48+u^4/384\) and
\(B(u)=1+12u/25-2u^3/25\):
\[
\begin{aligned}
&(1-u^2/20)E_4(u)
-B(u)(1-u^2/20+u^2/6)=uP_*(u),\\
P_*(u)={}&\frac1{50}-\frac u{24}
+\frac{119u^2}{6000}-\frac{7u^3}{1920}
+\frac{199u^4}{24000}-\frac{u^5}{7680}.
\end{aligned}
\tag{14}
\]
For \(0<u\le9/8\),
\[
\frac{119}{6000}-\frac7{1920}\frac98>\frac{157}{10000},
\qquad
\frac{199}{24000}-\frac1{7680}\frac98>\frac{81}{10000}.
\]
It follows that
\[
\begin{aligned}
P_*(u)
&\ge\frac1{50}-\frac u{24}
+\frac{157u^2}{10000}+\frac{81u^4}{10000}\\
&\ge\frac{1033u^2}{40000}-\frac u{24}
+\frac{431}{25600}>0.
\end{aligned}
\tag{15}
\]
The second line uses \((u^2-5/8)^2\ge0\).
Its final quadratic has positive leading coefficient and
discriminant \(-7007/2304000000\).

The series comparison
\[
\frac{\sinh u}{u}\le1+\frac{u^2}{6(1-u^2/20)},
\qquad e^{u/2}\ge E_4(u)
\]
is valid throughout this range. Both \(1-u^2/20\) and B(u)
are positive. Equations (14)--(15) therefore give
\(e^{u/2}\ge B(u)\sinh(u)/u\).
Dividing by \(2\sinh u\) proves the kernel bound in (11).

**The boundary tail.** TPW-T3 proved the unrounded bound
\[
F(11/10)>\frac{21241}{18000}.
\]
The kernel is decreasing, since
\((\log k)'=1/2-\coth u<0\).
Moreover
\[
e^{11/20}>
1+\frac{11}{20}+\frac{(11/20)^2}{2}
+\frac{(11/20)^3}{6}
=\frac{82991}{48000}>\frac{12}{7}.
\]
Writing \(r=e^{-11/20}<7/12\) gives
\[
k(11/10)=\frac{r}{1-r^4}
<\frac{12096}{18335}<\frac23.
\]
Consequently
\[
F(111/100)>
\frac{21241}{18000}-\frac1{100}\frac23
=\frac{21121}{18000}.
\]
This completes (11). QED.

## FCP-T4. A positive Gamma weight at the new endpoint

Set
\[
L_*=\frac{111}{200},\qquad
\alpha=\frac{49}{400},\qquad
d_*(x)=\alpha-\frac14\log(1-(x/L_*)^2).
\tag{16}
\]
For \(g\in\mathcal W_{L_*}\),
\[
\boxed{
\mathscr A[g]\ge\int_{-L_*}^{L_*}d_*(x)|g(x)|^2\,dx.
}
\tag{17}
\]
The weight is positive, even and log-convex.

**Proof of the source bound.** Retain the exact quadratic
correlation cancellation and Legendre estimate of TPW-T1:
\[
\begin{aligned}
\mathscr A[g]\ge{}&
\left(\frac{\Omega_0}{2}+F(2L)+\frac12\log2+
\frac{12L}{25}-\frac{8L^3}{75}+\frac{25}{24}\right)N_g\\
&+\int_I w_L(x)|g(x)|^2\,dx-\frac7{24}|b_2|^2,\\
w_L(x)={}&-\frac14\log(1-(x/L)^2),\\
|b_2|\le{}&\frac{\sqrt5}{96}L^2\cosh(L/2)\sqrt{N_g}.
\end{aligned}
\tag{18}
\]
The kernel hypothesis holds because \(2L_*=111/100<9/8\).
The cancellation uses exactly \(\int g=0\) and evenness; the
cosh moment controls \(b_2\) and does not set it to zero.

Using \(\sqrt5<3\), \(\cosh(L_*/2)<2\) and \(L_*<9/16\),
the lost coefficient is less than
\[
\frac{7L_*^4}{6144}
<\frac{45927}{402653184}<\frac1{8000}.
\]
Also
\[
\frac{12L_*}{25}-\frac{8L_*^3}{75}
=\frac{6204123}{25000000}.
\]
Substitution of (11) and (13) into (18) leaves the constant
\[
\begin{aligned}
&-\frac{26871}{10000}+\frac{25}{24}
+\frac{21121}{18000}+\frac{6931}{20000}
+\frac{6204123}{25000000}-\frac1{8000}\\
&\hspace{1cm}=\frac{49}{400}+\frac{1279}{28125000}
>\alpha .
\end{aligned}
\tag{19}
\]
This proves (17), retaining the negative Gamma scalar.

**Log-convexity.** For
\(D_\alpha(s)=\alpha-\frac14\log(1-s^2)\), \(|s|<1\),
put \(z=s^2\). Then
\[
D_\alpha''D_\alpha-(D_\alpha')^2
=\frac{2(1+z)D_\alpha-z}{4(1-z)^2}.
\]
Since \(D_\alpha\ge\alpha+z/4\) and \(\alpha\ge1/16\),
the numerator is at least
\[
2\alpha+(2\alpha-1/2)z+z^2/2
\ge\frac{4z^2-3z+1}{8}>0.
\]
Thus \(\log D_\alpha\) is convex. Scaling proves the claimed
properties of \(d_*\). In particular, evenness and convexity give
\[
d_*(x)d_*(x\pm u)\ge d_*(u/2)^2
\tag{20}
\]
whenever both arguments belong to I. QED.

## FCP-T5. The joint prime row remains strictly below one

At \(L_*=111/200\), the only active prime powers are 2 and 3.
Write \(u_n=\log n\) and define
\[
R(x)=\frac12\sum_{n=2,3}a_n
\left(
\frac{\mathbf1_I(x-u_n)}{\sqrt{d_*(x)d_*(x-u_n)}}+
\frac{\mathbf1_I(x+u_n)}{\sqrt{d_*(x)d_*(x+u_n)}}
\right).
\tag{21}
\]
Only terms with their shifted point in I are evaluated. Then
\[
\boxed{
\operatorname*{ess\,sup}_{x\in I}R(x)
\le\rho:=\frac{2451}{2455}<1.
}
\tag{22}
\]

**Active terms and coefficients.** From TPW, \(\log3<11/10\).
The positive atanh series also gives
\[
\log3>
1+\frac1{12}+\frac1{80}+\frac1{448}
=\frac{7379}{6720}>\frac{549}{500}.
\tag{23}
\]
Together with (13), these inequalities give
\[
L_*<u_2<u_3<2L_*<2\log2.
\]
Thus 2 and 3 are active, 4 and all larger prime powers are
outside the cutoff, and each shift has at most one neighbor
at a given point of I.

The elementary square bounds \(\sqrt2>140/99\) and
\(\sqrt3>12/7\) give
\[
\frac{a_2}{2}
<\frac{171567}{700000}<\frac{2451}{10000},
\qquad
\frac{a_3}{2}<\frac{77}{240}.
\tag{24}
\]

**Prime 2 on the entire window.** By (13),
\[
\frac{u_2}{2L_*}>\frac{78}{125},
\qquad \left(\frac{78}{125}\right)^2>\frac7{18}.
\]
For \(z=7/18\), the first six terms of the logarithmic series give
\[
-\log(1-z)>\sum_{j=1}^6\frac{z^j}{j}>\frac{123}{250}.
\tag{25}
\]
One entirely rational check of the last inequality is to bound the
six terms below, respectively, by
\[
\frac{3888}{10000},\quad
\frac{189}{2500},\quad
\frac{49}{2500},\quad
\frac{57}{10000},\quad
\frac{177}{100000},\quad
\frac{72}{125000}.
\]
Their sum is \(492046/1000000>123/250\).
Hence
\[
d_*(u_2/2)>
\frac{49}{400}+\frac{123}{1000}=\frac{491}{2000}.
\]
Equations (20), (24) and the one-neighbor property show that the
prime-2 part of R is everywhere at most
\[
\frac{2451/10000}{491/2000}=\rho.
\tag{26}
\]

**The prime-3 caps.** Let
\[
\mathcal C_3=\{x\in I:|x|>u_3-L_*\}.
\tag{27}
\]
These are exactly the points with a prime-3 neighbor, up to
endpoints of measure zero. For \(x\in\mathcal C_3\), (23) gives
\[
\frac{|x|}{L_*}>\frac{181}{185},
\qquad
1-\left(\frac{181}{185}\right)^2
=\frac{1464}{34225}<\frac1{23}.
\]
Since \(\log(1+t)<t\) for \(t>0\),
\[
\log23
=3\log3-\log(27/23)
>\frac{1647}{500}-\frac4{23}>\frac{78}{25}.
\tag{28}
\]
It follows that
\[
d_*(x)>
\frac{49}{400}+\frac{78}{100}>\frac9{10}
\quad(x\in\mathcal C_3).
\tag{29}
\]

If y is a prime-2 neighbor of such an x, \(u_2>L_*\) and
\(|x|<L_*\) give
\[
\frac{|y|}{L_*}\ge\frac{u_2}{L_*}-1>\frac{49}{200}.
\]
As \((49/200)^2>3/50\) and \(-\log(1-t)>t\),
\[
d_*(y)>
\frac{49}{400}+\frac3{200}=\frac{11}{80}.
\]
Therefore
\[
\sqrt{d_*(x)d_*(y)}>
\sqrt{\frac{99}{800}}>\frac{351}{1000},
\]
and the prime-2 contribution on \(\mathcal C_3\) is less than
\[
\frac{2451}{3510}=\frac{817}{1170}.
\tag{30}
\]

For the prime-3 pair use the midpoint estimate (20):
\[
\frac{u_3}{2L_*}>\frac{183}{185},\qquad
1-\left(\frac{183}{185}\right)^2
=\frac{736}{34225}<\frac1{46}.
\]
From (13) and (28), \(\log46>38131/10000\). Thus
\[
d_*(u_3/2)>
\frac{49}{400}+\frac{38131}{40000}
=\frac{43031}{40000}>\frac{43}{40}.
\]
Its contribution to R is less than
\[
\frac{77/240}{43/40}=\frac{77}{258}.
\tag{31}
\]
The joint bound on the caps is consequently
\[
R(x)<
\frac{817}{1170}+\frac{77}{258}
=\frac{25073}{25155}<\frac{2451}{2455}.
\tag{32}
\]
Outside \(\mathcal C_3\) the prime-3 contribution is zero and
(26) applies. This proves (22). QED.

The estimates in (30)--(32) are simultaneous at the same x.
The proof does not add the separate global maxima for the two primes.

## FCP-T6. Actual positivity extends to \(L=111/200\)

For every \(0<L\le111/200\) and \(g\in\mathcal W_L\),
\[
\boxed{
q_L[g]\ge\frac{49}{245500}\|g\|_2^2
\ge\frac1{6000}\|g\|_2^2 .
}
\tag{33}
\]
The second inequality is strict for nonzero g.

**Proof.** At \(L_*\), the weighted pair inequality of BPC-T3 is
\[
a|g(x)g(y)|
\le\frac{a}{2\sqrt{d_*(x)d_*(y)}}
\bigl(d_*(x)|g(x)|^2+d_*(y)|g(y)|^2\bigr).
\]
Integrate each translated pair once. Changing variables in its
second contribution gives exactly the two entries in (21), hence
\[
\sum_{n=2,3}a_n|C_g(u_n)|
\le\int_I R(x)d_*(x)|g(x)|^2\,dx.
\]
Combine (1), (17) and (22):
\[
q_{L_*}[g]\ge
\frac4{2455}\int_I d_*(x)|g(x)|^2\,dx
\ge\frac4{2455}\frac{49}{400}N_g
=\frac{49}{245500}N_g.
\]
No sign is assumed for the individual correlations.

For smaller L, the fixed moments give
\(\mathcal W_L\subset\mathcal W_{L_*}\).
Zero extension preserves the source form and norm. Any extra shift
with disjoint support has zero correlation, so the endpoint bound
applies unchanged. Finally \(49/245500>1/6000\). QED.

The restricted compact operator from FPW/BPC/TPW consequently obeys
\[
\boxed{
\|K_L^\circ\|\le\frac{24000}{24001}<1
\quad(0<L\le111/200).
}
\tag{34}
\]
Indeed the closed restricted energy has
\(G_L^\circ\ge(c_L+1/6000)I\), while
\[
0<c_L<
\frac{4902}{10000}+\frac{77}{120}
+\frac{26871}{10000}<4.
\]
With \(K_L^\circ=c_L(G_L^\circ)^{-1}\), this proves (34).
The extension from the smooth neutral core is by closure.

For the diagnostic parameter \(\tau(\lambda)=\lambda/(1+\lambda)\),
with \(\tau(\infty)=1\), the same consequence is
\[
I-\tau(\lambda)K_L^\circ\ge\frac1{24001}I
\quad(0\le\lambda\le\infty,\ 0<L\le111/200).
\tag{35}
\]
This parameter is a retained diagnostic coupling; it does not
identify the scalar zero with infinity.

## FCP-T7. Positive discarded probability coexists with exact positivity

At every \(L\le L_*\), use depth \(K_p=35\) for each active prime.
There are at most the primes 2 and 3. The bound in (5) satisfies
\[
\boxed{
\eta_{L,\mathbf K}<
6\,2^{-18}+8\,3^{-18}
<\frac7{131072}<\frac1{12000}.
}
\tag{36}
\]
At an empty active-prime set the left side is zero and the same
strict upper bounds hold. At the endpoint both probabilities
\(\epsilon_2=2^{-18}\) and \(\epsilon_3=3^{-18}\) are positive.

**Proof.** The exact factors are
\[
\kappa_2=3+2\sqrt2<6,\qquad
\kappa_3=2+\sqrt3<4.
\]
Use \(\log2<1\), \(\log3<2\) and
\(r_p^{K_p+1}=p^{-18}\) in (5). Then
\[
6\,2^{-18}+8\,3^{-18}<14\,2^{-18}=\frac7{131072},
\]
and \(7\cdot12000<131072\). This proves (36).
Combining (5), (33) and (36) gives
\[
\boxed{
\widetilde q_{L,\mathbf K}[g]
\ge(1/6000-\eta_{L,\mathbf K})N_g
\ge\frac1{12000}N_g,\qquad g\in\mathcal W_L.
}
\tag{37}
\]
For nonzero g the last inequality is strict. QED.

Here the exact source positivity was proved first in FCP-T6,
and (37) is its finite-record consequence. It is not an independent
numerical validation or a circular proof of (33).
For a later certification starting from finite-record data, (9)
provides the converse sufficient route with both errors retained.

The event-tail probability is not a scalar to append to every
arithmetic equality. Its source consequence is the weighted bound
(5). Density of the irrational numbers, by itself, does not choose a
probability measure or impose a nonzero probability floor.

## Results and the remaining RH problem

| Statement | Scope established here |
|---|---|
| Prime record error | Actual-source operator error at most \(\eta_{L,\mathbf K}\) for each finite L |
| Combined error enclosure | Equation (9), with both prime-record error and source-memory residual |
| Wider positive window | Same two-moment family, \(L\le111/200\), margin at least \(1/6000\) |
| Positive finite record | Depth \(K=35\), source error less than \(1/12000\), on that family |
| Restricted compact operator | Norm at most \(24000/24001\) in the new range |
| Full even source at the new width | Not proved positive by this note |
| All-width RH criterion | Still unproved |

The retained fixed-family criterion is
\[
\mathrm{RH}\quad\Longleftrightarrow\quad
q_L[g]\ge0\text{ for every }L>0,\ g\in\mathcal W_L,
\tag{38}
\]
as proved in FPW-T1. The present theorem extends one bounded
range in that criterion. It supplies no all-width bound.

The earlier stronger margins remain available: \(1/80\) for
\(L\le1/2\) from BPC and \(1/500\) for \(L\le11/20\) from TPW.
The established unrestricted range remains PTE's \(L\le1/32\).
No larger uniform native finite-bank range or new lower constant
for a first negative-bank scale is claimed.

The next all-width obstruction has a precise form. For the full
source, either establish the sign of every exact S in (7), or
prove a sufficient outward enclosure such as (9) at every required
width. Finite-record tails can be made small at any fixed width,
but this alone does not bound the source-memory residual or prove
that the retained matrix is positive. For the restricted route,
the joint Gamma--prime estimate must control every active prime
power as the support grows.

## Dependencies and verification scope

- [PCM-1](RH_PRIME_CUT_MEMORY.md): normalized geometric records,
  the coherent prime cut, exact scale and finite-record error.
- [ASM-1](RH_ARITHMETIC_SEAM_REDUCTION.md): explicit positive
  Legendre tail, source-column domain, exact Schur memory and
  finite trial residual.
- [TPW-1](RH_TWO_PRIME_WINDOW.md): quadratic correlation
  cancellation, exact kernel polynomial and retained rational bounds.
- [BPC-1](RH_BOUNDARY_PRIME_CAPS.md): weighted translated-pair
  comparison and boundary-cap method.
- [FPW-1](RH_FIRST_PRIME_WINDOW.md): fixed two-moment RH criterion,
  Legendre energy and control of the quadratic component.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): actual Gamma--prime source
  normalization and the closed energy construction.

The Schur complement, bounded perturbation estimate and weighted
square inequalities are standard. The work here is their explicit
assembly for the retained source and the stated rational extension
of its positive window; no literature-wide novelty claim is made.

Only this note and the README entry are changed. The document is
stacked on draft PR #114, with a [skip ci] documentation commit.
Verification is limited to written argument review, document
structure, repository links, exact immutable-head content and PR
metadata. No scientific computation, spectrum, source replay,
independent agent or workflow dispatch is performed. The PR is
left as a draft and is not merged.
