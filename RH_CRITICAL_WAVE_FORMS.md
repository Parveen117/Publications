# CWF-1 — Critical wave forms and finite prime windows

17 October 2026. Written continuation of the actual-source programme
at Publications commit 6f0f32a2dff30b2bd6d9bc317127593ebcea516f.

SQR proves positivity through one explicit square-root bank width.
BRR returns fixed reference banks at square-root depth and leaves
positivity at every scaled width as the remaining condition.

This continuation derives the actual form seen at that critical scale.
Smooth combinations of the native integer waves converge, after one
specified normalization, to a source quadratic form. Its arithmetic
expression contains the full Gamma term, the pole term and finitely
many prime powers in each finite window.

The resulting condition is a version of classical Weil positivity.
The derivation below identifies its exact relation to the native
moment banks and proves the converse using the retained finite-pair
detector. It also locates the first prime threshold and proves an
explicit positive margin in the small window inherited from SQR.

**Status: the RH proof is not complete.** Transport does not supply
the signs of all new directions. The full arithmetic inequality at
every window remains unproved. No actual negative witness is asserted.
These are written analytic derivations, not new machine-checked
source certificates.

## Actual source, wave and scale

Retain
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\qquad
a=b^2\ge a_0=3249/4.
\]
For a positive-height nontrivial zero beta+i gamma put
\[
v_\mu=\gamma-i\delta,\qquad \delta=\beta-\tfrac12,\qquad
\mu=v_\mu^2,\qquad y_a(\mu)=\frac a{a+\mu}.
\]
The multiplicity is nu_mu. The source and moment banks are
\[
\mathcal L_a[p]=\frac1{q(a)}
 \sum_\mu\frac{\nu_\mu}{a+\mu}p(y_a(\mu)),\qquad
c_n(a)=\mathcal L_a[y^n],\qquad
H_N(a)=(c_{i+j}(a))_{i,j<N}.
\tag{1}
\]
MOR retains q(a)>0, Re(mu)>0, |y_a|<1 and the absolute source
budget. The retained zero count gives
N(T)=O(T log(T+2)), hence sum nu_mu/(1+gamma^2)<infinity.

For the same native wave,
\[
\tau_k(y)=T_k(2y-1),\qquad
\omega_a(\mu)=\frac{b+i v_\mu}{b-i v_\mu},\qquad
\tau_k(y_a)=\frac{\omega_a^k+\omega_a^{-k}}2.
\]
CWR gives
\[
\max(|\omega_a|,|\omega_a|^{-1})
\le R_a=\frac{b+1/2}{b-1/2},\qquad \log R_a\le2/b.
\tag{2}
\]
No positive actual zero measure or RH assumption enters these bounds.

A window of half-width L means the real vector space
\[
\mathcal V_L=
\{g\in C_c^\infty((-L,L),\mathbb R):g(-t)=g(t)\}.
\]
Use the Fourier convention
\[
\widehat g(v)=\int_{\mathbb R}g(t)e^{ivt}\,dt.
\tag{3}
\]
The variable t is an auxiliary wave coordinate, not physical time.
L is a Fourier half-width. It is distinct from the scaled bank width
denoted L in BRR: a bank with N/sqrt(a) approaching L/2 has
N/sqrt(rho) approaching sqrt(a_0)L/2, where rho=a/a_0.

## CWF-T1. The actual critical wave limit exists with a uniform zero-tail bound

For g in V_L define h_b=2/b and the finite real polynomial
\[
\boxed{
P_{b,g}(y)=h_b\sum_{k\in\mathbb Z}g(kh_b)\tau_{|k|}(y).
}
\tag{4}
\]
It has degree at most floor(Lb/2). For g,h in V_L define
\[
\mathscr Q_\xi(g,h)
=\sum_\mu\nu_\mu\widehat g(v_\mu)\widehat h(v_\mu).
\tag{5}
\]
The sum is absolutely convergent and real. Then
\[
\boxed{
\lim_{b\to\infty}
b^2q(b^2)\,
\mathcal L_{b^2}[P_{b,g}P_{b,h}]
=\mathscr Q_\xi(g,h).
}
\tag{6}
\]
Thus every finite matrix of these boundary readings is a limit of
specified real congruences of the actual unshifted banks, multiplied
by the positive scalar b^2 q(b^2).

**Proof.** Evenness of g and (2) give the exact source identity
\[
P_{b,g}(y_a(\mu))
=h_b\sum_{k\in\mathbb Z}g(kh_b)\omega_a(\mu)^k.
\]
For each fixed v_mu,
log(omega_a)=2i v_mu/b+O_mu(b^(-3)).
On the support of the sum, |k|<=Lb/2. Consequently this is a
Riemann sum tending to the Fourier integral in (3).

The passage through the complete source sum needs a uniform bound.
For a finitely supported sequence A_k, two applications of
\[
(1-\omega)\sum_k A_k\omega^k
=\sum_k(A_k-A_{k-1})\omega^k
\]
give a second difference. For A_k=g(kh_b), that difference is a
double integral of g'' over two intervals of length h_b. Its
support extends the original interval by at most two mesh steps.
The weighted absolute sum is therefore at most C_g h_b^2:
the number of terms is O_g(h_b^(-1)), and
R_a^(Lb/2+2)<=exp(L+4/b) is uniformly bounded.
Here and below C_g is independent of b and of the source zero.

Since
\[
|1-\omega_a|=\frac{2|v_\mu|}{|b-i v_\mu|},\qquad
\frac{h_b}{|1-\omega_a|}
\le\frac1{\gamma}+\frac1b,
\]
we obtain
\[
\boxed{
|P_{b,g}(y_a(\mu))|
\le C_g(\gamma^{-1}+b^{-1})^2.
}
\tag{7}
\]
All actual gamma are positive. On any bounded height interval
there are only finitely many terms.

Also, since |delta|<1/2 and b>=57/2,
\[
\left|\frac{b^2}{b^2+\mu}\right|
\le\frac{2b^2}{b^2+\gamma^2}.
\]
For b>=R>=2, split the tail at gamma=b. The absolute contribution
of gamma>R to the left side of (6) is bounded by a constant
depending only on g,h times
\[
\sum_{\gamma>R}\frac{\nu_\mu}{\gamma^4}
+b^{-2}\sum_{\gamma>b}\frac{\nu_\mu}{\gamma^2}
=O_{g,h}\left(\frac{\log R}{R^3}\right).
\tag{8}
\]
The first sum bounds the R<gamma<=b part; adding its remaining
terms only enlarges the bound. Partial summation of the retained
zero count proves the last estimate, uniformly for b>=R.

For each fixed source term, b^2/(b^2+mu)->1. The uniform tail
bound and the finite-height Riemann-sum limit prove (6).
Repeated integration by parts in (3), on the strip |Im(v)|<=1/2,
gives rapid horizontal decay of its Fourier transform. That proves
absolute convergence of (5). Source conjugation makes it real.
QED.

The diagonal source sum in (5) is a sum of
nu_mu times hat(g)(v_mu)^2, not absolute squares. Only under RH
are all v_mu real, making those squares nonnegative.

For later comparison, the arcsine reference
dnu_*=dy/(pi sqrt(y(1-y))) gives the exact norm
\[
\int_0^1 P_{b,g}(y)^2\,d\nu_*(y)
=h_b^2\sum_{k\in\mathbb Z}g(kh_b)^2,\qquad
\lim_{b\to\infty}b\int P_{b,g}^2\,d\nu_*=2\|g\|_2^2.
\tag{9}
\]
Indeed the constant coefficient is h_b g(0), the positive wave
coefficients are 2h_b g(kh_b), and the arcsine squared norm
of each nonconstant tau_k is 1/2. The last limit is a real
Riemann sum.

Equation (6) concerns each fixed smooth profile or finite collection
of profiles. It does not assert operator-norm convergence for every
possible b-dependent vector in a growing Hankel bank.

## CWF-T2. Each finite window has an exact finite-prime arithmetic form

Let psi=Gamma'/Gamma and put
\[
\Omega(t)=\Re\psi(1/4+it/2)-\log\pi,\qquad
(g*g)(u)=\int_{\mathbb R}g(x)g(u-x)\,dx.
\]
For g in V_L,
\[
\boxed{
\begin{aligned}
\mathscr Q_\xi(g,g)
={}&\left(\int_{\mathbb R}g(t)\cosh(t/2)\,dt\right)^2\\
&+\frac1{4\pi}\int_{\mathbb R}
|\widehat g(t)|^2\Omega(t)\,dt\\
&-\sum_{\substack{n\ge2\\n<e^{2L}}}
\frac{\Lambda(n)}{\sqrt n}\,(g*g)(\log n).
\end{aligned}
}
\tag{10}
\]
The last sum is exact and finite. Lambda(n)=log(p) when n=p^j,
and is zero otherwise. All integrals converge absolutely.

**Proof.** Use the classical Guinand-Weil explicit formula with
H(v)=hat(g)(v)^2 and inverse transform K=g*g. With the Fourier
normalization (3), that formula is
\[
\sum_{\rho}H((\rho-1/2)/i)
=H(i/2)+H(-i/2)
+\frac1{2\pi}\int_{\mathbb R}H(t)\Omega(t)\,dt
-\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}
 [K(\log n)+K(-\log n)].
\tag{11}
\]
Its sum includes all nontrivial zeros with multiplicity, without RH.
The Fourier transform of a compact smooth function is entire and
rapidly decreasing in any fixed horizontal strip, as required here.

Both H and K are even. Pairing upper and lower zeros makes the
left side 2 Q_xi(g,g). The two pole terms coincide and equal
hat(g)(i/2)^2. On the real axis hat(g) is real, so H=|hat(g)|^2.
Finally K is supported strictly inside (-2L,2L). Dividing (11)
by two proves (10). QED.

Equation (11) is classical, not a new RH theorem. Its normalization
and contour derivation are given in Paul Garrett,
[Guinand's explicit formula, Theorem 0.1](https://www-users.cse.umn.edu/~garrett/m/mfms/notes_2019-20/04_Guinand_explicit_fml.pdf),
with the original references to Guinand (1947) and Weil (1952).
The same even-form normalization is displayed in the
[American Institute of Mathematics statement](https://aimath.org/WWN/rh/articles/html/75a/).
It comes from the same completed-zeta logarithmic derivative:
the Gamma factor supplies Omega, the poles at 0 and 1 supply the
two evaluations, and the full Euler series supplies the prime powers.
The new programme identification here is (6), with the normalization
and the full-source convergence proved above.

In particular,
\[
\boxed{
2L\le\log2\ \Longrightarrow\
\text{the prime sum in (10) is empty};
\qquad
p^j\text{ can enter only when }2L>j\log p.
}
\tag{12}
\]
The wave and arithmetic scales are consequently related by
\[
\boxed{
N_b(L)=\lfloor L\sqrt a/2\rfloor+1,\qquad
n<e^{2L}.
}
\tag{13}
\]
Equivalently a limiting wave width N/sqrt(a)=c corresponds to
prime-power reach n<exp(4c). This statement concerns the limiting
profiles in (4), not the absence of prime dependence from the
finite-anchor moments themselves.

A finite prime sum still acts on infinitely many test functions
in V_L. It is not automatically a finite-dimensional positivity
certificate.

## CWF-T3. The real even boundary tests already give an RH equivalence

The following statements are equivalent:
\[
\boxed{
\begin{aligned}
&\mathrm{RH};\\
&\mathscr Q_\xi(g,g)\ge0
 \quad\text{for every }L>0\text{ and every }g\in\mathcal V_L.
\end{aligned}
}
\tag{14}
\]
The second statement is precisely the full arithmetic inequality
in (10). The sufficiency is proved below with the retained native
finite-pair detector; it does not assume positivity of the complex
source weights.

**Proof.** Under RH every v_mu=gamma is real. The summands
nu_mu hat(g)(gamma)^2 are nonnegative, proving necessity.

For the converse, suppose there is a nonreal source point
z=a_0/(a_0+mu). Consider the auxiliary bilinear functional
\[
\mathscr M[p]=\sum_\mu
 \frac{\nu_\mu}{(a_0+\mu)^2}
 p\!\left(\frac{a_0}{a_0+\mu}\right).
\tag{15}
\]
Its weights are nonzero, conjugate-paired and absolutely summable.
Its support is the same discrete set with sole accumulation point
zero as in the retained finite detector.

The detector's pair-isolation proof applies to these weights too.
Explicitly, annihilate the finitely many other support points of
modulus at least |z| with a real polynomial A. On the two directions
y^m A(y), y^(m+1) A(y), the selected pair contributes
\[
\begin{pmatrix}
2\Re C_m&2\Re(C_m z)\\
2\Re(C_m z)&2\Re(C_m z^2)
\end{pmatrix},\qquad
C_m=\frac{\nu_z z^{2m+2}A(z)^2}{a_0^2}\ne0.
\]
Its determinant is -4|C_m|^2(Im(z))^2. All remaining support has
modulus at most some s_0<|z|, and its matrix entries are
O(s_0^(2m)). After division by |z|^(4m), the determinant error
tends to zero, while the selected determinant stays strictly
negative. Thus some real polynomial p has M[p^2]<0.
This uses exactly D1-D2's algebra and tail isolation, with the
new absolutely summable nonzero weights.

Put
\[
F_p(v)=\frac1{a_0+v^2}\,
p\!\left(\frac{a_0}{a_0+v^2}\right).
\]
Then
\[
\sum_\mu\nu_\mu F_p(v_\mu)^2=\mathscr M[p^2]<0.
\tag{16}
\]
This rational profile is the Fourier transform of a real even
exponentially decreasing function. In detail,
\[
k_j(t)=\left.\frac{(-1)^j}{j!}\,
\partial_A^j\frac{e^{-\sqrt A|t|}}{2\sqrt A}\right|_{A=a_0}
\quad\Longrightarrow\quad
\widehat k_j(v)=(a_0+v^2)^{-j-1}.
\tag{17}
\]
The j=0 identity follows by integrating the two exponentials on
the positive and negative half-lines; differentiating proves
the rest. Consequently g_p=sum_j p_j a_0^j k_j has transform F_p.

The function g_p and its weak first derivative are integrable
with weight exp(|t|/2), since sqrt(a_0)=57/2. Symmetric smooth
cutoffs followed by symmetric mollification produce real even
g_m in C_c^infinity(R) converging to g_p in the norm
\[
\|g\|_{1,1/2}+\|g'\|_{1,1/2},\qquad
\|g\|_{1,1/2}=\int_{\mathbb R}e^{|t|/2}|g(t)|\,dt.
\tag{18}
\]
For example, first cut off at |t|<=R with uniformly bounded
first derivative and let R increase; exponential decay controls
the removed part. Then let the mollifier radius decrease.
The exponential weight changes by a bounded factor under those
small translations, proving convergence in both terms.

These norms are uniformly bounded for g_m. Integration by parts
and the weighted L1 bound imply, uniformly on |Im(v)|<=1/2,
\[
|\widehat g_m(v)|+|F_p(v)|
\le\frac C{1+|\Re v|}.
\]
Weighted L1 convergence also gives pointwise convergence of the
transforms on this strip. Since sum nu_mu/(1+gamma^2) is finite,
dominated convergence in the source sum gives
\[
\mathscr Q_\xi(g_m,g_m)
\longrightarrow\sum_\mu\nu_\mu F_p(v_\mu)^2<0.
\]
Some smooth compact even g_m therefore has a negative reading,
contradicting the asserted positivity. This proves sufficiency.
QED.

This is an even real version of the classical Weil positivity
criterion, with an explicit adapter to the retained source.
For context, Connes and Consani,
[Weil positivity and Trace formula, the archimedean place](https://arxiv.org/abs/2006.13771),
Introduction and Appendix C, discuss compact-support positivity
and the fact that only finitely many primes enter each test.
Their paper does not supply positivity at all supports; no such
conclusion is imported here.

## CWF-T4. SQR supplies a positive boundary margin before the first prime

Retain kappa=e^(-9)/(2pi). For every g in V_(2kappa),
\[
\boxed{\mathscr Q_\xi(g,g)\ge\frac14\|g\|_2^2.}
\tag{19}
\]
The same holds in every smaller window.

**Proof.** Since g has compact support strictly inside
(-2kappa,2kappa), there is L'<2kappa containing its support.
For all sufficiently large b, the degree of P_(b,g) plus one
is at most floor(kappa b). SQR-T5 therefore applies:
\[
\mathcal L_{b^2}[P_{b,g}^2]
\ge\frac1{2\ell_b}\int P_{b,g}^2\,d\nu_*,
\qquad \ell_b=\log(b/(2\pi)).
\]
Multiply by b^2 q(b^2). SQR-T1 gives
b^2 q(b^2)=b ell_b/4+O(1), and (9) gives the limiting
reference norm. Taking the limit in (6) yields
(1/4)||g||_2^2, proving (19). QED.

This positive boundary margin belongs to the actual source.
Nevertheless its window is smaller than the first prime threshold:
\[
\boxed{2\kappa<\frac12\log2.}
\tag{20}
\]
For example kappa<1/3072 and log2>1/2 already prove this.
Thus the inherited boundary inequality lies in a window where the
finite prime sum in (10) vanishes.

This does not remove prime powers from SQR's finite-anchor proof:
that proof retained their full contribution. It identifies which
prime-translation terms survive in this particular critical limit.
In particular, the proven small window does not itself cover the
first nonzero prime term, let alone every finite prime window.

## CWF-T5. A finite window profile locates any possible hidden sign

Define the real spectral profile
\[
\lambda_\xi(L)=
\inf_{\substack{g\in\mathcal V_L\\g\ne0}}
\frac{\mathscr Q_\xi(g,g)}{\|g\|_2^2},\qquad L>0.
\tag{21}
\]
For every finite L this is a finite real number and is
nonincreasing as L increases. More explicitly,
\[
\boxed{
\lambda_\xi(L)\ge
\frac{\Omega(0)}2
-\sum_{\substack{n\ge2\\n<e^{2L}}}\frac{\Lambda(n)}{\sqrt n},
\qquad
\lambda_\xi(L)\ge\frac14\quad(0<L\le2\kappa).
}
\tag{22}
\]
Let
\[
L_*=\inf\{L>0:\lambda_\xi(L)<0\},\qquad \inf\varnothing=+\infty.
\]
Then
\[
\boxed{\mathrm{RH}\ \Longleftrightarrow\ L_*=+\infty.}
\tag{23}
\]
If RH fails, L_* is finite, L_*>=2kappa, and
\[
\boxed{
\kappa\le\liminf_{a\to\infty}\frac{n_-(a)}{\sqrt a}
\le\limsup_{a\to\infty}\frac{n_-(a)}{\sqrt a}
\le\frac{L_*}{2}<\infty.
}
\tag{24}
\]
No value of L_* or negative actual profile is evaluated here.

**Proof.** The digamma series used in SQR gives, for A>0,
\[
\Re\psi(A+iu)-\psi(A)
=\sum_{j=0}^\infty
\frac{u^2}{(j+A)((j+A)^2+u^2)}\ge0.
\]
Thus Omega(t)>=Omega(0). The same series, split at j of
order |u|, bounds its growth by O(log(2+|u|)).
Parseval gives integral |hat(g)|^2=2pi||g||_2^2, while
Cauchy--Schwarz gives |(g*g)(u)|<=||g||_2^2.
The pole square in (10) is nonnegative. These observations
prove the first bound in (22); T4 proves the second.
A nonzero smooth even test exists in every window and has
finite reading, giving a finite upper value for the infimum.
Inclusion V_(L1) subset V_(L2) for L1<L2 proves monotonicity.

The equivalence (23) is exactly T3: a negative reading exists
if and only if some window has negative infimum.
For any such reading g in V_L, (6) is negative for all sufficiently
large b. Its polynomial P_(b,g) has degree at most floor(Lb/2).
The scalar b^2 q(b^2) is positive, so
\[
n_-(b^2)\le\lfloor Lb/2\rfloor+1
\quad\text{for all sufficiently large }b.
\tag{25}
\]
Divide by b and take the upper limit. Taking negative windows
arbitrarily close to the infimum L_* proves the last bound in
(24). The first bound is SQR-T7. QED.

Neither convergence of n_-(a)/sqrt(a) nor equality with L_*/2
is claimed. Such conclusions would require control of all
changing bank vectors, stronger than the fixed-profile limit (6).

The profile (21) is a source property. It has not been identified
with a thermodynamic capacity, rest mass or propagation speed.

## CWF-T6. Enlarging a window requires a new coupled positivity inequality

For any finite real even test family g_1,...,g_m in V_L, its
exact boundary matrix is
\[
\boxed{
\begin{aligned}
\mathsf Q_{ij}
={}&\widehat g_i(i/2)\widehat g_j(i/2)\\
&+\frac1{4\pi}\int_{\mathbb R}
\widehat g_i(t)\widehat g_j(t)\Omega(t)\,dt\\
&-\sum_{\substack{n\ge2\\n<e^{2L}}}
\frac{\Lambda(n)}{\sqrt n}(g_i*g_j)(\log n).
\end{aligned}
}
\tag{26}
\]
Suppose a chosen old bank A is positive definite, and write
the enlarged matrix as [[A,C],[C^T,B]]. Then
\[
\boxed{
\mathsf Q\succeq0
\quad\Longleftrightarrow\quad
B-C^TA^{-1}C\succeq0.
}
\tag{27}
\]
The signs of the newly included prime translations are not
automatically positive updates.

**Proof.** Polarization of (10) gives (26).
Completing the square in the old coefficients gives
\[
x^TAx+2x^TCy+y^TBy
=(x+A^{-1}Cy)^TA(x+A^{-1}Cy)
+y^T(B-C^TA^{-1}C)y,
\]
which proves (27).

To see the sign issue for a single translation, fix u=log n>0.
Choose a nonzero real even smooth h supported in (-epsilon,epsilon),
where 2epsilon<u, and choose a window L>u+epsilon. Set
\[
g_\pm(t)=h(t)\pm[h(t-u)+h(t+u)].
\]
The separated convolution supports give exactly
\[
\boxed{(g_\pm*g_\pm)(u)=\pm2\|h\|_2^2.}
\tag{28}
\]
Indeed only the two products with centers 0 and u contribute at u;
all other center sums are at distance at least u.
Thus for Lambda(n)>0 the individual contribution
-Lambda(n)(g*g)(log n)/sqrt(n) takes both signs on the same
real even test space. QED.

Equation (28) is not a negative witness for the full form:
the Gamma term, pole term and all other active prime powers
are still present. It shows why adding prime contributions
one at a time is not a positive-semidefinite update rule.

The Schur condition is the same classical source sign test already
used in QRI and OQD, now with every entry supplied by a finite-prime
arithmetic expression. Neither observer transport nor old-bank
positivity proves its new residual sign.

## The exact remaining proof task

The all-width condition is now the explicit inequality
\[
\boxed{
\begin{aligned}
\sum_{\substack{n\ge2\\n<e^{2L}}}
\frac{\Lambda(n)}{\sqrt n}(g*g)(\log n)
\ \le\ &
\left(\int g(t)\cosh(t/2)\,dt\right)^2\\
&+\frac1{4\pi}\int|\widehat g(t)|^2\Omega(t)\,dt
\quad\text{for every }L>0,\ g\in\mathcal V_L.
\end{aligned}
}
\tag{29}
\]
Proving (29) would prove RH by T3. This continuation proves the
native adapter to that classical condition and the small-window
margin (19). It does not prove (29) at arbitrary width.

The outcome of the attempted propagation is therefore precise:
the signs of the inherited window are controlled, while widening
the observer introduces additional coupled arithmetic conditions.
A global argument controlling those conditions could close RH;
none follows from the current transport and congruence identities
alone. No theorem here rules out such a future argument.

| Established here | Limit of the conclusion |
|---|---|
| Exact critical-scale limit of actual wave banks | Fixed profiles; no operator-norm limit for all changing vectors |
| Exact arithmetic form with finitely many prime powers per window | Infinitely many observer functions still need control |
| Even real boundary positivity is equivalent to RH | Its universal premise remains unproved |
| Small-window margin at least 1/4 | The covered window is before the first prime threshold |
| A negative boundary profile would return a finite negative bank | No actual negative profile is exhibited |
| Finite arithmetic matrices and their exact Schur gate | New residual signs are not supplied automatically |

## Attribution and execution

Retained source lineage:
[finite detector](RH_FINITE_OBSERVER_DETECTOR.md),
[arithmetic source](RH_ARITHMETIC_GATE.md),
[MOR](RH_MOVING_OBSERVER_RETURN.md),
[CWR](RH_CUT_WAVE_RETURN.md),
[QRI](RH_REFLECTION_RATIO_GATE.md),
[OQD](RH_OBSERVER_QUOTIENT_DESCENT.md),
[SQR](RH_SQUARE_ROOT_BANKS.md), and
[BRR](RH_BOUNDARY_OBSERVER_RETURN.md).

The Guinand-Weil explicit formula and Weil positivity criterion are
classical results with the external attribution above. The actual
critical wave limit, its complete tail argument, the finite-detector
adapter, and the inherited boundary margin are derived here.
No literature-wide novelty or optimality claim is made.

Only this note and its README entry are added. No local scientific
computation, independent agent, numerical source evaluation or
scientific replay is performed. Existing data, engines, certificates
and workflows remain unchanged. The draft contains written proofs,
not new machine-checked certificates. The documentation commit
includes [skip ci]. It is stacked on PR #104 without a merge.
