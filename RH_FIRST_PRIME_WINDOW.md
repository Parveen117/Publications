# FPW-1 — A first-prime window with fixed neutral moments

17 October 2026. Written continuation at Publications commit
8b77d145ca64475729f750bed1d051053c067bcd, following
[PTE-1](RH_PRIME_TRANSPORT_ENERGY.md).

The preceding unrestricted estimate stops before any prime translation
can meet the support. This note crosses the first prime threshold on
one fixed family of real even tests: their constant and cosh(x/2)
moments vanish. These two conditions remain the same at every width.

For this family the actual source form has the explicit lower bound
\[
\boxed{
0<L\le2/5
\quad\Longrightarrow\quad
\mathscr Q_\xi(g,g)\ge\frac1{250}\|g\|_2^2.
}
\]
At L=2/5 the prime 2 is active. Its correlation is bounded using
disjoint support caps. A Legendre energy estimate, the two moments,
and the exact Gamma boundary energy supply the positive margin.

The fixed restriction does not discard the global RH criterion:
positivity on these same tests at every finite width is equivalent
to RH. That statement is proved by adapting the retained finite
detector and approximating its rational profile while enforcing the
two moments exactly.

**Status: RH remains unproved.** The new bound is for the specified
restricted family, not every test at L<=2/5. The unrestricted
L<=1/32 result from PTE is retained. Positivity on the fixed family
for all L is still open. These are written analytic proofs, not
machine-checked certificates. No scientific computation, source
evaluation, independent agent or numerical replay is performed.

## Source, normalization and the fixed family

Use PTE's real even space
\[
\mathcal V_L=
\{g\in C_c^\infty((-L,L),\mathbb R):g(-x)=g(x)\},
\qquad
\widehat g(v)=\int_{\mathbb R}g(x)e^{ivx}\,dx.
\]
All functions are extended by zero. Write
\[
C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,\qquad
\ell(g)=\int_{\mathbb R}g(x)\cosh(x/2)\,dx.
\]
Define, for every L>0, exactly the same two neutral moments:
\[
\boxed{
\mathcal W_L=
\left\{g\in\mathcal V_L:
\int_{\mathbb R}g(x)\,dx=0,\quad \ell(g)=0\right\}.
}
\tag{1}
\]
These are two independent real linear conditions. Evenness makes
them equivalent to hat(g)(0)=hat(g)(i/2)=hat(g)(-i/2)=0.
The spaces W_L are infinite dimensional.

Retain the complete source identity from CWF-T2 and PTE-T2:
\[
\begin{aligned}
\mathscr Q_\xi(g,g)
&=\ell(g)^2+\mathscr A[g]
-\sum_{2\le n<e^{2L}}a_n C_g(\log n),\\
\mathscr A[g]
&=\frac1{4\pi}\int_{\mathbb R}|\widehat g(t)|^2\Omega(t)\,dt,\\
\Omega(t)
&=\Re\psi(1/4+it/2)-\log\pi,\qquad
a_n=\frac{\Lambda(n)}{\sqrt n}.
\end{aligned}
\tag{2}
\]
Here Lambda(p^j)=log p and it is zero otherwise. No RH assumption
enters this identity. In particular the Gamma term is not replaced
by a positive surrogate.

The exact energy identity and boundary decomposition are
\[
\begin{gathered}
\Omega_0=-\gamma_E-\frac{\pi}{2}-3\log2-\log\pi,\qquad
k(u)=\frac{e^{u/2}}{2\sinh u},\\
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}),\\
\mathscr A[g]=\frac{\Omega_0}{2}\|g\|_2^2+
\mathcal E_\Gamma[g],
\end{gathered}
\tag{3}
\]
\[
\begin{aligned}
\mathcal E_\Gamma[g]
={}&\frac14\int_{-L}^L\int_{-L}^L
k(|x-y|)|g(x)-g(y)|^2\,dy\,dx\\
&+\frac12\int_{-L}^L|g(x)|^2
\bigl(F(L-x)+F(L+x)\bigr)\,dx .
\end{aligned}
\tag{4}
\]
Since F is convex, the second line is at least F(L)||g||^2.
The singular double integral is finite for smooth g.
The proof and digamma references are retained in PTE-T2--T3.

The same two transform zeros occur in the related theorem of
Connes and Consani, *Weil positivity and Trace formula: the archimedean
place*, Theorem 1, for multiplicative support
[2^(-1/2),2^(1/2)]; see their
[author-hosted text](https://alainconnes.org/wp-content/uploads/Selecta.pdf).
The logarithmic coordinate has half-width (log 2)/2.
Their archimedean functional has twice the normalization of A here;
compare their Appendix B, equations (153)--(154).
That theorem supplies context for the constraints. Its positivity
estimate or numerical spectral calculations are not inputs to the
bounds below. The polynomial energy identity used below is proved
in full; no claim of a new classical spectral identity is made.

## FPW-T1. The fixed neutral family retains the global RH criterion

The following statements are equivalent:
\[
\boxed{
\mathrm{RH}
\quad\Longleftrightarrow\quad
\mathscr Q_\xi(g,g)\ge0
\text{ for every }L>0\text{ and }g\in\mathcal W_L .
}
\tag{5}
\]
This is an all-width equivalence. It does not assert that positivity
on W_L implies positivity on V_L at the same fixed width.

**Proof.** Retain CWF's source notation
\[
v_\mu=\gamma-i\delta,\quad
\mu=v_\mu^2,\quad
\delta=\beta-\tfrac12,\quad
a_0=3249/4,\quad y_\mu=\frac{a_0}{a_0+\mu}.
\]
Here beta+i gamma is a positive-height nontrivial zero, with
multiplicity nu_mu. The retained source bounds give Re(mu)>0,
|y_mu|<1, conjugation symmetry, and
\[
\sum_\mu\frac{\nu_\mu}{1+\gamma^2}<\infty.
\tag{6}
\]
For real even smooth compact g the exact source form is
\[
\mathscr Q_\xi(g,g)
=\sum_\mu\nu_\mu\widehat g(v_\mu)^2.
\tag{7}
\]
These are bilinear squares, not absolute squares.
Under RH all v_mu are real and hat(g) is real there, proving
the forward implication of (5).

For the converse suppose RH fails. Use the fixed even polynomial
\[
R(v)=v^2(v^2+1/4)
\]
and define weights
\[
\alpha_\mu
=\frac{\nu_\mu R(v_\mu)^2}{(a_0+\mu)^6}.
\tag{8}
\]
They are nonzero and conjugate-paired. They are absolutely summable:
at large height their absolute values are O(nu_mu gamma^(-4)).
Nonvanishing follows from Re(mu)>0, which excludes mu=0,-1/4.
The support y_mu is the same discrete set as in the retained
finite detector, with sole accumulation point 0.

For clarity, the needed detector argument works for these new
weights without a positivity assumption on them.
Choose a nonreal support point z. A real polynomial A annihilates
the finitely many other support points with modulus at least |z|,
while A(z) is nonzero. On y^m A(y) and y^(m+1) A(y), the selected
conjugate pair contributes
\[
\begin{pmatrix}
2\Re B_m&2\Re(B_m z)\\
2\Re(B_m z)&2\Re(B_m z^2)
\end{pmatrix},\qquad
B_m=\alpha_z z^{2m}A(z)^2 .
\]
Its determinant is
\[
-4|B_m|^2(\Im z)^2<0.
\tag{9}
\]
All remaining support has modulus at most s_0<|z|.
Absolute summability bounds each remainder entry by O(s_0^(2m)).
After division by |z|^(4m), the determinant error tends to zero,
and the negative main determinant has fixed nonzero magnitude.
Thus some real polynomial p satisfies
\[
\sum_\mu\alpha_\mu p(y_\mu)^2<0.
\tag{10}
\]
This is precisely the pair-isolation proof of
[the finite detector, D1--D2](RH_FINITE_OBSERVER_DETECTOR.md),
with the weights in (8).

Now put
\[
H_p(v)=\frac{p(a_0/(a_0+v^2))}{(a_0+v^2)^3},
\qquad F_p(v)=R(v)H_p(v).
\tag{11}
\]
Equation (10) says sum nu_mu F_p(v_mu)^2<0.
Both profiles are real even rational functions, holomorphic in
the strip |Im(v)|<=1/2.

We approximate F_p by transforms of tests in W_L with the two
moments exact at each stage. For m>=0 define
\[
k_m(x)=
\left.\frac{(-1)^m}{m!}\partial_A^m
\frac{e^{-\sqrt A|x|}}{2\sqrt A}\right|_{A=a_0},
\qquad
\widehat k_m(v)=(a_0+v^2)^{-m-1}.
\tag{12}
\]
If p(y)=sum_j p_j y^j, the inverse transform of H_p is
\[
h_p=\sum_j p_j a_0^j k_{j+2}.
\]
Each such kernel belongs to weighted W^(5,1), with weight
e^(|x|/2). Here is a direct regularity check. The elementary
kernel k_0 and its first weak derivative are integrable with
that weight because sqrt(a_0)=57/2. Also
k_m is the (m+1)-fold convolution of k_0 and
\[
D^2 k_m=a_0 k_m-k_{m-1}\quad(m\ge1).
\]
For m>=2 this gives
\[
D^4 k_m=a_0^2k_m-2a_0k_{m-1}+k_{m-2}.
\]
The right side and its first weak derivative are integrable
with the same weight. The intermediate derivatives follow
from the D^2 identity. This proves the asserted regularity.

Real even smooth cutoffs followed by real even mollification
give h_j in C_c^infinity(R) converging to h_p in that weighted
W^(5,1) norm. Define
\[
g_j=D^2(D^2-1/4)h_j.
\tag{13}
\]
The transform multiplier in (13) is exactly R(v).
Integration by parts against 1 and cosh(x/2) proves that
both moments in (1) vanish. Thus g_j belongs to some W_(L_j)
with finite L_j.

The g_j converge to D^2(D^2-1/4)h_p in weighted W^(1,1).
Their transforms converge pointwise to F_p on
|Im(v)|<=1/2. Integration by parts and the bounded weighted
norms give one constant C such that
\[
|\widehat g_j(v)|\le\frac{C}{1+|\Re v|}
\quad\text{throughout that strip, for all }j.
\tag{14}
\]
The same bound holds for F_p. By (6), dominated convergence
therefore applies to the bilinear squared source sums in (7).
Their limit is the strictly negative value in (10).
Some finite j has Q_xi(g_j,g_j)<0, contradicting the
right side of (5). This proves the converse. QED.

The construction only assumes a hypothetical off-line pair for
the converse. It neither finds such a pair nor exhibits an actual
negative test. The support width needed by a hypothetical witness
has no universal bound in this argument.

## FPW-T2. The Gamma kernel controls a Legendre energy

For 0<u<=1,
\[
\boxed{k(u)\ge\frac1{2u}+\frac15.}
\tag{15}
\]
For 0<L<=1/2 and real even mean-zero g in V_L, write
\[
\phi_n(x)=\sqrt{\frac{2n+1}{2L}}P_n(x/L),\qquad
c_n=\langle g,\phi_n\rangle_{L^2(-L,L)}.
\]
Then
\[
\boxed{
\mathscr A[g]\ge
\left(F(L)+\frac L5+\frac{25}{24}+\frac{\Omega_0}{2}\right)
\|g\|_2^2-\frac7{24}|c_2|^2 .
}
\tag{16}
\]
The c_n in this section are Legendre coefficients, not the
source moments c_n(a) of the earlier bank notes.

**Proof of the kernel bound.** Successive terms of the sinh
series after u^2/6 have ratio at most u^2/20. Therefore
\[
\frac{\sinh u}{u}
\le1+\frac{u^2}{6(1-u^2/20)}.
\]
The positive exponential series gives
e^(u/2)>=1+u/2+u^2/8+u^3/48. Exact multiplication yields
\[
\begin{aligned}
&(1-u^2/20)(1+u/2+u^2/8+u^3/48)\\
&\quad -(1+2u/5)(1-u^2/20+u^2/6)\\
&=u\left(\frac1{10}-\frac u{24}
-\frac{61u^2}{1200}-\frac{u^3}{160}-\frac{u^4}{960}\right)
\ge\frac u{4800}>0 .
\end{aligned}
\]
The last inequality uses u<=1. Dividing by 1-u^2/20 gives
e^(u/2)>=(1+2u/5)sinh(u)/u, proving (15).

**Proof of the polynomial energy bound.** On [-1,1] consider
\[
(Af)(x)=\int_{-1}^1\frac{f(x)-f(y)}{|x-y|}\,dy.
\]
For the Legendre polynomials,
\[
AP_n=2H_n P_n,\qquad
H_n=\sum_{j=1}^n\frac1j,\quad H_0=0.
\tag{17}
\]
To prove (17), on a monomial direct division and integration give
\[
A(x^n)=
\sum_{j=0}^{n-1}x^{n-1-j}
\frac{2x^{j+1}-1-(-1)^{j+1}}{j+1}.
\]
This is a polynomial of degree at most n with leading coefficient
2H_n. The operator is symmetric on polynomials, because
\[
\langle f,Ah\rangle
=\frac12\int_{-1}^1\int_{-1}^1
\frac{(f(x)-f(y))(h(x)-h(y))}{|x-y|}\,dy\,dx .
\]
It preserves each lower-degree polynomial space.
Consequently AP_n-2H_nP_n has degree less than n and is
orthogonal to every polynomial of degree less than n.
It is zero, proving (17).

Scaling to (-L,L), using finite Legendre projections and the
nonnegative energy of their remainders, gives
\[
\frac18\int_{-L}^L\int_{-L}^L
\frac{|g(x)-g(y)|^2}{|x-y|}\,dy\,dx
\ge\frac12\sum_{n\ge0}H_n|c_n|^2.
\tag{18}
\]
More explicitly, the energy cross term between a finite
projection and its L2-orthogonal remainder vanishes by (17).
This first proves the inequality for every finite partial sum;
passing to their increasing limit proves (18).
Completeness of the Legendre basis gives sum |c_n|^2=||g||^2.

Evenness removes odd n and mean zero removes n=0.
Since H_2/2=3/4 and H_4/2=25/24,
\[
\frac12\sum H_n|c_n|^2
\ge\frac{25}{24}\|g\|^2-\frac7{24}|c_2|^2.
\tag{19}
\]
The constant part of (15) also contributes
\[
\frac1{20}\int_{-L}^L\int_{-L}^L
|g(x)-g(y)|^2\,dy\,dx=\frac L5\|g\|^2
\tag{20}
\]
by mean zero. Since |x-y|<=2L<=1, apply (15) inside (4),
combine (18)--(20) with the boundary term F(L)||g||^2,
and add Omega_0||g||^2/2. This proves (16). QED.

## FPW-T3. The second moment condition controls the remaining low mode

For every g in W_L,
\[
\boxed{
|c_2|\le\frac{\sqrt5}{96}L^2\cosh(L/2)\|g\|_2.
}
\tag{21}
\]
In particular, for 0<L<=2/5,
\[
|c_2|\le\frac1{100}\|g\|_2.
\tag{22}
\]

**Proof.** Taylor's formula gives
\[
\cosh(x/2)=1+x^2/8+R_4(x),\qquad
|R_4(x)|\le\frac{\cosh(L/2)}{384}|x|^4.
\]
Both moments in (1) vanish, and
x^2=(L^2/3)(2P_2(x/L)+1). Hence
\[
\frac{L^2}{12}\langle g,P_2(x/L)\rangle
=-\langle g,R_4\rangle.
\]
Also
\[
\|R_4\|_{L^2(-L,L)}
\le\frac{\cosh(L/2)\sqrt{2L^9}}{1152}.
\]
Cauchy--Schwarz and the normalization of phi_2 prove (21). For L<=2/5, use sqrt5<3 and cosh(L/2)<2:
the resulting bound is (6/96)(4/25)||g||=||g||/100.
This proves (22). QED.

Thus cosh neutrality does not exactly remove the quadratic
Legendre mode. Its possible contribution is bounded explicitly.
From (16) and (22),
\[
\boxed{
\mathscr A[g]\ge
\left(F(L)+\frac L5+\frac{25}{24}
+\frac{\Omega_0}{2}-\frac7{240000}\right)\|g\|_2^2,
\quad g\in\mathcal W_L,\quad L\le2/5.
}
\tag{23}
\]

## FPW-T4. Explicit rational margins cross the first prime threshold

For every 0<L<=2/5 and every g in W_L,
\[
\boxed{
\mathscr A[g]\ge\frac{127}{500}\|g\|_2^2,\qquad
\mathscr Q_\xi(g,g)\ge\frac1{250}\|g\|_2^2.
}
\tag{24}
\]
At the endpoint L=2/5, the only active prime power is n=2.

**Proof of the Gamma margin.** We first establish the elementary
bounds
\[
\frac{\Omega_0}{2}>-\frac{27}{10},
\qquad F(2/5)>\frac{11}{6}.
\tag{25}
\]
All estimates below are inequalities of exact series and rational
numbers, with no decimal evaluation of the source or its spectrum.

For the first inequality, midpoint convexity of 1/x gives
gamma_E<H_n-log(n+1/2). At n=1,
\[
\gamma_E<1-\log(3/2)<119/200,
\]
because the positive atanh series gives
log(3/2)>2/5+2/375=152/375>81/200.
Also
\[
\pi/2<11/7,\qquad \log2<25/36,\qquad \log\pi<23/20.
\]
For the logarithm of 2, write
log2=2 sum_(j>=0) (1/3)^(2j+1)/(2j+1).
Bound all tail denominators for j>=1 below by 3 to obtain
log2<2/3+1/36=25/36.
For the logarithm of pi, the positive exponential series gives
\[
e^{23/20}>
\sum_{j=0}^5\frac{(23/20)^j}{j!}
=\frac{1211276443}{384000000}>\frac{22}{7}>\pi .
\]
Consequently
\[
\gamma_E+\pi/2+3\log2+\log\pi
<
\frac{119}{200}+\frac{11}{7}+\frac{25}{12}+\frac{23}{20}
=\frac{22679}{4200}<\frac{27}{5}.
\]
This proves the first inequality in (25).

For the second, the alternating exponential series gives
\[
e^{-1/5}\ge1-\frac15+\frac1{50}-\frac1{750}
=\frac{307}{375}.
\]
At r=307/375,
(1+r)/(1-r)=341/34>10 and (1-r)/(1+r)=34/341.
Monotonicity, the arctangent subtraction identity, and
arctan(t)<=t for t>=0 yield
\[
F(2/5)>
\frac12\log10+\frac{\pi}{4}-\frac{34}{341}.
\]
The classical Archimedean bounds give pi>223/71>157/50.
Further,
\[
\log2>\frac23+\frac2{81}+\frac2{1215}
=\frac{842}{1215},\qquad
\log(5/4)>\frac29 .
\]
Since 10=8(5/4), this gives log10>932/405>23/10.
It follows that
\[
F(2/5)>
\frac{23}{20}+\frac{157}{200}-\frac{34}{341}
=\frac{387}{200}-\frac{34}{341}>\frac{11}{6},
\]
where the last inequality is 61/600>34/341.
This proves (25).

The function F(L)+L/5 decreases for 0<L<=2/5,
since its derivative is -k(L)+1/5<0.
Therefore (23)--(25) give
\[
\begin{aligned}
\mathscr A[g]
&\ge
\left(\frac{11}{6}+\frac2{25}+\frac{25}{24}
-\frac{27}{10}-\frac7{240000}\right)\|g\|^2\\
&=\frac{61193}{240000}\|g\|^2
\ge\frac{127}{500}\|g\|^2.
\end{aligned}
\tag{26}
\]

**Proof of the prime bound.** The elementary bounds
2/3<log2<7/10 and log3>1 imply
\[
\log2<4/5<\log3,\qquad L\le2/5<\log2.
\tag{27}
\]
For example log3>1 follows from the positive atanh series
with parameter 1/2. Thus n<e^(2L) can include only n=2.
It is included when L>(log2)/2, in particular at L=2/5.

Put u=log2. If 2L<=u, C_g(u)=0.
Otherwise the two overlap caps
\[
I_-=[-L,L-u],\qquad I_+=[u-L,L]
\]
are disjoint because L<u. Changing variables in the correlation
and applying 2|ab|<=a^2+b^2 gives
\[
\begin{aligned}
|C_g(u)|
&=\left|\int_{-L}^{L-u}g(x)g(x+u)\,dx\right|\\
&\le\frac12\int_{I_-\cup I_+}|g(x)|^2\,dx
\le\frac12\|g\|^2.
\end{aligned}
\tag{28}
\]
Since sqrt2>7/5 and log2<7/10,
\[
\frac{\log2}{\sqrt2}|C_g(\log2)|
\le\frac{\log2}{2\sqrt2}\|g\|^2
<\frac14\|g\|^2.
\tag{29}
\]
For g in W_L the pole square ell(g)^2 is exactly zero.
Combining (2), (26) and (29) proves
\[
\mathscr Q_\xi(g,g)
\ge\left(\frac{127}{500}-\frac14\right)\|g\|^2
=\frac1{250}\|g\|^2 .
\]
This proves (24). QED.

The prime term is genuinely capable of contributing on W_(2/5).
For an explicit existence argument, let u=log2 and take a nonzero
real even smooth bump b of sufficiently small support around zero.
The even function b(x-u/2)+b(x+u/2) lies inside (-2/5,2/5)
and has correlation ||b||^2 at shift u. Correct its two moments
by two even smooth functions supported in |x|<u-2/5 with
independent moment vectors. Such functions exist because cosh(x/2)
is not constant there. These central corrections cannot interact
at shift u with any function supported in (-2/5,2/5).
The corrected function belongs to W_(2/5) and retains that positive
correlation. Thus the first-prime estimate is not merely an empty
prime sum or a family on which the prime term vanishes identically.

## FPW-T5. The restricted compact operator has a strict contraction bound

Use the Hilbert space and closed form domain H_L,D_L from PTE-T4,
complexifying the real form by Q(g_1+ig_2)=Q(g_1)+Q(g_2).
Define
\[
H_L^\circ=
\{g\in H_L:\langle g,1\rangle=\langle g,\cosh(x/2)\rangle=0\},
\qquad D_L^\circ=D_L\cap H_L^\circ.
\tag{30}
\]
Restrict the positive energy E_L of PTE to D_L^circ, and let
G_L^circ be its associated self-adjoint operator on H_L^circ.
Retain
\[
S_L=\sum_{2\le n<e^{2L}}a_n,\qquad
c_L=S_L-\Omega_0/2>0,\qquad
K_L^\circ=c_L(G_L^\circ)^{-1}.
\tag{31}
\]
Then K_L^circ is positive and compact, and
\[
\boxed{
\mathrm{RH}\quad\Longleftrightarrow\quad
\|K_L^\circ\|\le1\text{ for every }L>0.
}
\tag{32}
\]
The proven first-prime range has
\[
\boxed{
0<L\le2/5\quad\Longrightarrow\quad
\|K_L^\circ\|\le\frac{800}{801}<1 .
}
\tag{33}
\]

**Proof.** The two moments are bounded functionals on H_L.
Their common kernel is closed. To verify the needed smooth form
core, choose two fixed real even compact smooth functions with
independent moment vectors. They can be bumps near two different
absolute positions. Given a smooth approximation to any element
of D_L^circ, correct its two moment errors by a linear combination
of those fixed functions. The correction coefficients tend to zero,
because the moment matrix is invertible and the approximation
converges in L2. Hence the corrected approximations converge in
the form norm as well. The same argument with L2 approximation
proves density in H_L^circ.

The restricted closed form is therefore well defined, densely
defined and has W_L as a real form core. The bound G_L>=F(L)I
and compact embedding from PTE persist under this restriction.
Thus G_L^circ>=F(L)I and its inverse is compact.

The exact restricted source identity is still
\[
\mathscr Q_\xi[g]=
\|(G_L^\circ)^{1/2}g\|^2-c_L\|g\|^2.
\]
The spectral theorem proves Q>=0 on W_L if and only if
||K_L^circ||<=1. Equation (32) follows from T1.

For L<=2/5, T4 gives
G_L^circ>=(c_L+1/250)I.
Only n=2 can contribute to S_L, and
a_2=log2/sqrt2<1/2. Together with (25),
\[
0<c_L<\frac12+\frac{27}{10}=\frac{16}{5}.
\]
It follows that
\[
\|K_L^\circ\|
\le\frac{c_L}{c_L+1/250}\le\frac{800}{801}.
\]
This proves (33). QED.

Restriction is performed before taking the inverse.
In general (G_L^circ)^(-1) is not the orthogonal compression
of G_L^(-1). No bound for the full K_L at L=2/5 is inferred
from (33).

For the PTE diagnostic coupling tau(lambda)=lambda/(1+lambda),
define B_L^circ(lambda)=I-tau(lambda)K_L^circ and set tau(infinity)=1.
Equation (33) gives the explicit endpoint-uniform bound
\[
B_L^\circ(\lambda)\ge\frac1{801}I
\quad(0\le\lambda\le\infty,\ 0<L\le2/5).
\tag{34}
\]
Thus its negative count is zero throughout this parameter range.
For unrestricted L the all-width condition (32) remains to be
proved. Equality of ordinary Fredholm indices still gives no
substitute for that sign condition.

## What this continuation settles and what remains

| Reading | Established conclusion |
|---|---|
| Fixed family W_L | Exactly two neutral moments, independent of L |
| Global coverage | Positivity on W_L for every L is equivalent to RH |
| Gamma contribution, L<=2/5 | A[g]>=127/500 times ||g||^2 on W_L |
| First prime | n=2 can contribute; its adverse cost is less than ||g||^2/4 |
| Complete source, L<=2/5 | Q_xi[g]>=||g||^2/250 on W_L |
| Restricted operator | ||K_L^circ||<=800/801 in the same range |
| Unrestricted tests | PTE's L<=1/32 result is retained; no larger full-family bound here |
| All widths | The remaining positivity assertion is not proved |

The next arithmetic step must bound the complete prime correlations
as the support grows while retaining test coverage. The cap form
in (28) retains more information than its final norm bound and can
be compared with the exact position-dependent boundary term in (4).
That observation is a possible next estimate, not an established
all-prime domination theorem.

CWF's actual finite-bank return is for each fixed profile or finite
collection. The new boundary bound does not supply operator-norm
convergence for an entire growing bank. No enlargement of SQR's
uniform positive bank sizes, no new lower constant for
n_-(a)/sqrt(a), and no actual negative witness are claimed.

The constant 1/250 is a dimensionless analytic margin in the
specified normalization. No identification with light speed,
mass, a thermodynamic ratio, or an electromagnetic limit is made.
The arithmetic source, its subtraction term and both neutral
moments are stated explicitly; the remaining RH condition is not
hidden in a choice of a positive comparison source.

## Dependency and verification record

- [CWF-1](RH_CRITICAL_WAVE_FORMS.md): exact source sum, critical
  wave return, and the fully normalized Guinand-Weil identity.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): exact Gamma difference
  energy, boundary decomposition and compact closed-form setup.
- [Finite detector D1--D2](RH_FINITE_OBSERVER_DETECTOR.md): the
  conjugate-pair determinant and finite polynomial isolation.
- [Connes--Consani, Theorem 1 and Appendix B](https://alainconnes.org/wp-content/uploads/Selecta.pdf):
  related fixed transform-zero conditions and normalization context;
  their positivity or numerical estimates are not used in T2--T4.
- [NIST DLMF 5.9.16](https://dlmf.nist.gov/5.9.E16) and
  [5.4.19](https://dlmf.nist.gov/5.4.E19): retained digamma
  integral and quarter-argument value underlying (3)--(4).

All new inequalities are accompanied by analytic arguments.
Text, links and exact committed content are the repository checks.
No scientific test, source computation, numerical spectrum,
workflow run or independent replay is requested by this note.
