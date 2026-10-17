# BPC-1 — Boundary energy and paired prime caps

17 October 2026. Written continuation at Publications commit
95ca1559eab67179b3890ed4cc0c51ec1d68cc3d, following
[FPW-1](RH_FIRST_PRIME_WINDOW.md).

The first-prime estimate bounded the Gamma boundary energy by one
constant and then charged the prime correlation against the full
norm. This loses the fact that a prime translation can only couple
specified parts of the support.

Keeping the position-dependent energy gives a stronger estimate.
The two points in a prime pair are compared together. An even
log-convex positive weight has its smallest product, among points
separated by a fixed translation, at the symmetric midpoint pair.
That product controls the correlation by weighted squares.

For exactly the fixed two-moment family of FPW, the resulting
written proof establishes
\[
\boxed{
0<L\le\frac12
\quad\Longrightarrow\quad
\mathscr Q_\xi(g,g)\ge\frac1{80}\|g\|_2^2.
}
\]
This extends the earlier L<=2/5, margin 1/250 result.
Prime 2 is active at the new endpoint. Prime 3 is not yet active.

Two general sufficient bounds are proved as well: a scalar bound
whose minimum is found at finitely many overlap boundaries, and a
weighted pair bound that retains more of the same geometry.

**Status: RH remains unproved.** The family still has the same two
fixed neutral moments. Its positivity for every L is RH-equivalent
by FPW-T1, but this note establishes only the stated finite range.
No larger unrestricted-family range or uniform finite-bank range
is asserted. These are analytic written proofs, not machine-checked
source certificates. No scientific computation, source evaluation,
independent agent or numerical replay is performed.

## Retained source and family

Let I=(-L,L), extend functions by zero, and retain
\[
\begin{aligned}
\mathcal W_L=\{g\in C_c^\infty(I,\mathbb R):\;&g(-x)=g(x),\\
&\int g(x)\,dx=0,\quad
\int g(x)\cosh(x/2)\,dx=0\}.
\end{aligned}
\tag{1}
\]
The same two conditions are used for every L. Write
\[
C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,\qquad
u_n=\log n,\qquad a_n=\frac{\Lambda(n)}{\sqrt n},
\]
and let
\[
\mathcal N_L=\{n=p^j:\ j\ge1,\ 2\le n<e^{2L}\}.
\]
The exact source form on W_L is
\[
\boxed{
\mathscr Q_\xi(g,g)=\mathscr A[g]
-\sum_{n\in\mathcal N_L}a_n C_g(u_n).
}
\tag{2}
\]
The pole square vanishes by the second moment in (1).
All prime powers in the support window are retained.

Use the PTE/FPW normalization
\[
\begin{gathered}
\mathscr A[g]=
\frac1{4\pi}\int_{\mathbb R}|\widehat g(t)|^2
\bigl(\Re\psi(1/4+it/2)-\log\pi\bigr)\,dt,\\
\widehat g(t)=\int_{\mathbb R}g(x)e^{itx}\,dx,\\
\Omega_0=-\gamma_E-\pi/2-3\log2-\log\pi,\\
k(u)=\frac{e^{u/2}}{2\sinh u},\qquad
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}).
\end{gathered}
\tag{3}
\]
The exact Gamma decomposition is
\[
\begin{aligned}
\mathscr A[g]
={}&\frac{\Omega_0}{2}\|g\|^2
+\frac14\int_I\int_I
k(|x-y|)|g(x)-g(y)|^2\,dy\,dx\\
&+\int_I\mathcal B_L(x)|g(x)|^2\,dx,\\
\mathcal B_L(x)
&=\frac12\bigl(F(L-x)+F(L+x)\bigr).
\end{aligned}
\tag{4}
\]
PTE-T2--T3 prove this identity from the digamma representation
[NIST DLMF 5.9.16](https://dlmf.nist.gov/5.9.E16).
The negative scalar Omega_0/2 remains in all lower bounds below.

FPW-T2--T3 give the two further inputs
\[
\frac18\int_I\int_I
\frac{|g(x)-g(y)|^2}{|x-y|}\,dy\,dx
\ge\frac{25}{24}\|g\|^2-\frac7{24}|b_2|^2,
\tag{5}
\]
\[
|b_2|\le
\frac{\sqrt5}{96}L^2\cosh(L/2)\|g\|,\qquad
b_2=\left\langle g,\sqrt{\frac5{2L}}P_2(x/L)\right\rangle.
\tag{6}
\]
Here b_2 is a Legendre coefficient. No source moment is renamed.

## BPC-T1. Prime translations have an explicit local norm cost

For x in I define the even step function
\[
\boxed{
V_L(x)=\frac12\sum_{n\in\mathcal N_L}a_n
\bigl(\mathbf1_I(x-u_n)+\mathbf1_I(x+u_n)\bigr).
}
\tag{7}
\]
For every real g supported in I,
\[
\boxed{
\sum_{n\in\mathcal N_L}a_n|C_g(u_n)|
\le\int_I V_L(x)|g(x)|^2\,dx.
}
\tag{8}
\]
At a positive radius 0<r<L, away from the finitely many endpoints,
\[
V_L(r)=\frac12\sum_{n\in\mathcal N_L}a_n
\bigl(\mathbf1_{\{r<L-u_n\}}+
      \mathbf1_{\{r>u_n-L\}}\bigr).
\tag{9}
\]

**Proof.** For u>0, let J_u=I intersect (I-u). Then
\[
|C_g(u)|
\le\int_{J_u}|g(x)g(x+u)|\,dx
\le\frac12\int_{J_u}
\bigl(|g(x)|^2+|g(x+u)|^2\bigr)\,dx.
\]
Changing variables in the second square gives (8), after summing.
The two indicators in (7) describe exactly those two portions
of the support. For 0<r<L their inequalities simplify to (9).
This proves the claims. QED.

The cases for one translation are consequently exact:

| Translation u | Local contribution to V_L(r), 0<r<L |
|---|---|
| 0<u<L | a for r<L-u, and a/2 for r>L-u |
| u=L | a/2 throughout, up to endpoints of measure zero |
| L<u<2L | zero for r<u-L, and a/2 for r>u-L |
| u>=2L | zero |

Thus a long translation charges only two edge caps.
For a shorter translation the overlap multiplicity can be two;
the disjoint-cap factor cannot then be used unchanged.

## BPC-T2. Boundary energy gives a finite sufficient scalar test

Suppose a real number sigma satisfies
\[
k(u)\ge\frac1{2u}+\sigma\qquad(0<u\le2L).
\tag{10}
\]
For every finite L some such number exists: k(u)-1/(2u)
extends continuously to u=0 with value 1/4 and has a finite
minimum on [0,2L]. In particular FPW-T2 proves that sigma=1/5
works whenever L<=1/2.

Define
\[
\begin{gathered}
w_L(x)=-\frac14\log\bigl(1-(x/L)^2\bigr),\quad x\in I,\\
\varepsilon_L=
\min\left(1,\frac{5L^4\cosh^2(L/2)}{9216}\right),\\
M(L,\sigma)=
\frac{\Omega_0}{2}+\frac{25}{24}
+F(2L)+\frac12\log2+2\sigma L-\frac7{24}\varepsilon_L .
\end{gathered}
\tag{11}
\]
Then
\[
\boxed{
\mathscr A[g]\ge
\int_I\bigl(M(L,\sigma)+w_L(x)\bigr)|g(x)|^2\,dx,
\qquad g\in\mathcal W_L .
}
\tag{12}
\]
Consequently, with
\[
\eta_L=\operatorname*{ess\,inf}_{0<r<L}
\bigl(w_L(r)-V_L(r)\bigr),
\]
the complete source satisfies
\[
\boxed{
\mathscr Q_\xi(g,g)\ge
\bigl(M(L,\sigma)+\eta_L\bigr)\|g\|^2 .
}
\tag{13}
\]

There is an exact finite expression for this particular scalar
lower bound. Put
\[
\mathcal D_L=\{0\}\cup
\{u_n-L:\ n\in\mathcal N_L,\ L<u_n<2L\}.
\]
Use right-hand limits for V_L at the elements of D_L.
Then
\[
\boxed{
\eta_L=\min_{r\in\mathcal D_L}
\bigl(w_L(r)-V_L(r+)\bigr).
}
\tag{14}
\]

**Proof of the energy bound.** Integrate (10) from v to 2L:
\[
F(v)\ge F(2L)+\frac12\log\frac{2L}{v}
+\sigma(2L-v),\qquad 0<v\le2L.
\]
Apply this to v=L-x and v=L+x and take the average. This gives
the position-dependent estimate
\[
\boxed{
\mathcal B_L(x)\ge
F(2L)+\frac12\log2+\sigma L+w_L(x).
}
\tag{15}
\]
The boundary weight tends to infinity at the two endpoints;
it is integrable against the smooth compact tests used here.

The interior part of (4), bounded using (10), is at least
the left side of (5) plus
\[
\frac{\sigma}{4}\int_I\int_I|g(x)-g(y)|^2\,dy\,dx
=\sigma L\|g\|^2.
\]
The identity uses mean zero and is valid even if sigma is negative.
Also (6) and orthonormality each bound |b_2|^2, so
|b_2|^2<=epsilon_L||g||^2.
Combining these facts with (4)--(5) and (15) proves (12).
Equation (8) then proves (13).

**Proof of the finite minimum.** The function w_L(r) increases
on [0,L), while V_L is constant between its finitely many jumps.
A minimum on each such interval is therefore its left endpoint
limit. By (9), a translation u_n<L makes V_L decrease at
r=L-u_n. Such a jump increases w_L-V_L and cannot create a
new minimum. A translation u_n>L can make w_L-V_L decrease
at r=u_n-L. These are exactly the additional points in D_L.
Simultaneous jumps are read using their combined right-hand limit.
A jump involving only the shorter translations is upward, and any
jump involving a longer one is included in D_L.
The initial limit r=0+ handles u_n=L as well.
Finally w_L tends to infinity as r tends to L, whereas V_L
is bounded. No missing endpoint can give a smaller value.
This proves (14). QED.

Equations (13)--(14) are sufficient tests for positivity.
They do not turn the full source form on W_L into a finite
matrix, nor give an equivalent finite criterion for RH.
A nonpositive lower bound here does not exhibit a negative
source test.

The scalar comparison still bounds each overlap point separately.
The next theorem retains the pairing and can be stronger.

## BPC-T3. Weighted pair domination uses the product of the two energies

Let d be a positive measurable function on I with inf_I d=m>0,
and suppose
\[
\mathscr A[g]\ge\int_I d(x)|g(x)|^2\,dx
\quad(g\in\mathcal W_L).
\tag{16}
\]
Assume the integrals on the right are finite for these tests.
Only evaluate the following summands when their shifted point
lies inside I:
\[
\begin{aligned}
R_{L,d}(x)=\frac12\sum_{n\in\mathcal N_L}a_n
\left(
\frac{\mathbf1_I(x-u_n)}{\sqrt{d(x)d(x-u_n)}}+
\frac{\mathbf1_I(x+u_n)}{\sqrt{d(x)d(x+u_n)}}
\right),\\
R_{L,d}=\operatorname*{ess\,sup}_{x\in I}R_{L,d}(x).
\end{aligned}
\tag{17}
\]
Then
\[
\boxed{
\mathscr Q_\xi(g,g)\ge
(1-R_{L,d})\int_I d(x)|g(x)|^2\,dx .
}
\tag{18}
\]
In particular if R_(L,d)<=1, then
\[
\mathscr Q_\xi(g,g)\ge
(1-R_{L,d})m\|g\|^2.
\tag{19}
\]

If d is even and log-convex, and every active translation
satisfies L<=u_n<2L, the simpler bound
\[
\boxed{
R_{L,d}\le
\sum_{n\in\mathcal N_L}\frac{a_n}{2d(u_n/2)}
}
\tag{20}
\]
holds. For the source here the translation condition is automatic
when L<=log2.

**Proof.** For a pair x,y=x+u inside I, weighted squares give
\[
a|g(x)g(y)|
\le
\frac{a}{2\sqrt{d(x)d(y)}}
\bigl(d(x)|g(x)|^2+d(y)|g(y)|^2\bigr).
\tag{21}
\]
Integrate over J_u, change variables in the second term, and sum
over n. The result is
\[
\sum_n a_n|C_g(u_n)|
\le\int_I R_{L,d}(x)d(x)|g(x)|^2\,dx.
\]
Equations (2) and (16) prove (18).
When 1-R_(L,d)>=0 the lower bound d>=m proves (19).
No such last inference is made with a negative coefficient.

For an even log-convex d, apply convexity of log d to x and -y.
Their midpoint is -u/2. Evenness gives
\[
\boxed{d(x)d(y)\ge d(u/2)^2.}
\tag{22}
\]
For u>=L the two caps J_u and J_u+u are disjoint, up to
a measure-zero endpoint when u=L. Thus a fixed x belongs
to at most one of the two caps for that translation.
Bound each denominator in (17) by d(u_n/2) and use that
multiplicity bound to obtain (20). QED.

This is a weighted square estimate with its overlap multiplicity
shown explicitly. The weight d must be supplied by the actual
Gamma lower bound (16); choosing a positive weight without
proving (16) would not establish source positivity.

## BPC-T4. An explicit Gamma weight at L=1/2 is even and log-convex

At L_*=1/2 define
\[
\boxed{
d_*(x)=\frac18-\frac14\log(1-4x^2),
\qquad |x|<1/2.
}
\tag{23}
\]
For g in W_(1/2),
\[
\boxed{
\mathscr A[g]\ge\int_{-1/2}^{1/2}d_*(x)|g(x)|^2\,dx .
}
\tag{24}
\]
The function d_* is even, positive and log-convex, with
inf d_*=1/8.

**Proof of the Gamma weight.** FPW's kernel bound gives
sigma=1/5 in (10). Its moment estimate (6), with
sqrt5<3 and cosh(1/4)<2, gives
\[
|b_2|\le\frac1{64}\|g\|.
\tag{25}
\]
The retained rational bounds include
Omega_0/2>-27/10 and log2>69/100.
It remains to prove the elementary estimate
\[
F(1)>\frac{31}{25}.
\tag{26}
\]

The alternating exponential series gives e^(-1/2)>=29/48.
For r=29/48 put q=(1-r)/(1+r)=19/77. Then
\[
\operatorname{artanh}(r)
=\log2+\frac12\log(77/76),\qquad
\arctan(r)=\frac{\pi}{4}-\arctan(q).
\]
The positive atanh series gives
\[
\log2>\frac{842}{1215}>\frac{693}{1000},\qquad
\frac12\log(77/76)>\frac1{153}>\frac3{500}.
\]
For 0<=q<=1/2, the alternating arctangent series gives
\[
\arctan q\le q-\frac{q^3}{3}+\frac{q^5}{5}
\le q-\frac{q^3}{4}.
\]
At q=19/77,
\[
q<\frac{247}{1000},\qquad q>\frac6{25},\qquad
\frac{q^3}{4}>\frac{54}{15625}>\frac3{1000}.
\]
Also pi/4>157/200 by the classical Archimedean lower bound.
Monotonicity of artanh and arctan now yields
\[
\begin{aligned}
F(1)
&\ge\operatorname{artanh}(29/48)+\arctan(29/48)\\
&>
\frac{693}{1000}+\frac3{500}+\frac{157}{200}
-\frac{247}{1000}+\frac3{1000}
=\frac{31}{25}.
\end{aligned}
\]
This proves (26) without evaluating a numerical spectrum.

Apply the proof of (12) using the sharper coefficient bound
(25). The constant part is greater than
\[
\begin{aligned}
&\frac{31}{25}+\frac{69}{200}+\frac15+\frac{25}{24}
-\frac{27}{10}-\frac7{98304}\\
&\hspace{2cm}=\frac{19}{150}-\frac7{98304}>\frac18.
\end{aligned}
\tag{27}
\]
For the last inequality, 19/150-1/8=1/600>7/98304.
The remaining position-dependent term is exactly w_(1/2)(x).
This proves (24).

**Proof of log-convexity.** Write
\[
D(t)=\frac18-\frac14\log(1-t^2),\quad |t|<1,\qquad
d_*(x)=D(2x).
\]
It is positive and even, and D(t)>=1/8+t^2/4. Direct derivatives
give
\[
D'(t)=\frac{t}{2(1-t^2)},\qquad
D''(t)=\frac{1+t^2}{2(1-t^2)^2}.
\]
Therefore
\[
D''(t)D(t)-D'(t)^2
\ge\frac{1-t^2+2t^4}{16(1-t^2)^2}>0.
\tag{28}
\]
This is exactly positivity of the numerator in (log D)''.
Scaling preserves log-convexity, proving the assertions. QED.

The full negative scalar from the Gamma term was included in
(27). The positive weight (23) is a consequence of the estimate,
not a replacement of that source term by assumption.

## BPC-T5. Pairing extends the positive window to 1/2 with margin 1/80

For exactly the family (1),
\[
\boxed{
0<L\le\frac12
\quad\Longrightarrow\quad
\mathscr Q_\xi(g,g)\ge\frac1{80}\|g\|^2
\quad\text{for all }g\in\mathcal W_L .
}
\tag{29}
\]

**Proof at the endpoint.** At L_*=1/2 the prime-power condition
is n<e. The elementary inequalities log2<1<log3 show that
only n=2 contributes. Put
\[
u=\log2,\qquad a=\frac{\log2}{\sqrt2}.
\]
Here 1/2<u<1, so its overlap caps are disjoint and T3 applies.
It remains to bound the midpoint weight d_*(u/2).

For 0<z<1 the positive atanh series gives
\[
-\frac14\log(1-z)
=\frac12\operatorname{artanh}\frac{z}{2-z}
>\frac{z}{2(2-z)}.
\]
Since u>69/100, take z=u^2>4761/10000 to obtain
\[
-\frac14\log(1-u^2)>
\frac{4761}{30478}>\frac2{13}.
\]
The last comparison is 61893>60956.
Consequently
\[
\boxed{d_*(u/2)>\frac18+\frac2{13}=\frac{29}{104}.}
\tag{30}
\]
FPW's elementary bounds log2<7/10 and sqrt2>7/5 give
a/2<1/4. Thus (20) and (30) imply
\[
R_{L_*,d_*}\le\frac{a}{2d_*(u/2)}
\le\frac{26}{29}<1 .
\tag{31}
\]
Equations (18), (23) and (24) yield
\[
\begin{aligned}
\mathscr Q_\xi(g,g)
&\ge\frac3{29}\int_{-1/2}^{1/2}d_*(x)|g(x)|^2\,dx\\
&\ge\frac3{232}\|g\|^2
\ge\frac1{80}\|g\|^2 .
\end{aligned}
\tag{32}
\]
The last comparison is 240>232.

**Transfer to smaller windows.** The two moments in (1) are fixed,
so W_L is a subspace of W_(1/2) whenever L<=1/2.
The actual source sum, hence Q_xi(g,g), is independent of which
larger window is used to describe the support. Any extra prime
translations at a larger cutoff have zero correlation with a
profile whose support is smaller. Norms are unchanged by the
zero extension. Applying the endpoint estimate proves (29).
QED.

More generally the restricted source bottom
\[
\lambda_\xi^\circ(L)=
\inf_{\substack{g\in\mathcal W_L\\\|g\|=1}}
\mathscr Q_\xi(g,g)
\tag{33}
\]
is nonincreasing in L by the same inclusion. Thus an endpoint
lower bound covers every smaller width. This statement concerns
the source form; no monotonicity of a differently normalized
compact operator is inferred from it.

The endpoint is beyond (log2)/2 and before (log3)/2.
Prime 2 remains capable of nonzero correlation on the family,
as shown by FPW's two-cap construction. The next prime threshold
is therefore still an outstanding step.

## BPC-T6. The restricted compact norm and coupling bounds improve

Use the restricted closed energy and compact operator from
FPW-T5:
\[
\mathscr Q_\xi[g]
=\|(G_L^\circ)^{1/2}g\|^2-c_L\|g\|^2,\qquad
K_L^\circ=c_L(G_L^\circ)^{-1},
\]
\[
c_L=\sum_{n\in\mathcal N_L}a_n-\Omega_0/2>0.
\tag{34}
\]
The same smooth neutral tests form the core; the moment constraints
are imposed before the inverse is taken.

For 0<L<=1/2, (29) extends by form closure and gives
G_L^circ>=(c_L+1/80)I. Only prime 2 can contribute, and the
retained bounds give 0<c_L<1/2+27/10=16/5. Hence
\[
\boxed{
\|K_L^\circ\|
\le\frac{c_L}{c_L+1/80}\le\frac{256}{257}<1 .
}
\tag{35}
\]
For tau(lambda)=lambda/(1+lambda), including tau(infinity)=1,
\[
\boxed{
I-\tau(\lambda)K_L^\circ\ge\frac1{257}I
\quad(0\le\lambda\le\infty,\ 0<L\le1/2).
}
\tag{36}
\]
This proves zero negative count throughout this restricted range.

The all-width equivalence remains
\[
\mathrm{RH}
\quad\Longleftrightarrow\quad
\mathscr Q_\xi\ge0\text{ on }\mathcal W_L\text{ for all }L
\quad\Longleftrightarrow\quad
\|K_L^\circ\|\le1\text{ for all }L .
\tag{37}
\]
It is the retained FPW-T1/T5 theorem, with no new moment conditions.
Equation (35) establishes only its stated range. It does not
bound the unrestricted K_L on all of V_L.

## Scope and the next unresolved estimate

| Result | FPW-1 | BPC-1 |
|---|---|---|
| Neutral moments | Constant and cosh(x/2) | The same two moments |
| Proven half-width | L<=2/5 | L<=1/2 |
| Stated source margin | 1/250 | 1/80 |
| Restricted compact norm | At most 800/801 | At most 256/257 |
| Active prime at endpoint | 2 | 2 |
| RH for all widths | Unproved | Unproved |

The pointwise bound (13) and weighted bound (18) preserve the
complete arithmetic source. For a new window, a valid kernel
estimate, an actual Gamma lower weight, and a bound on the full
active-prime sum are still required. Adding another prime adds
another positive cost to (17); the one-prime ratio (31) cannot
simply be reused.

The paired estimate explains why treating each cap by its worst
point lost information. If one member of a pair is near a support
edge, its boundary weight is large. Log-convexity controls the
product of the two weights and makes the symmetric pair the
least favorable one. The proof uses that relation while retaining
both fixed neutral moments.

The remaining task is to control the combined translated pairs
as further primes enter, with a sufficient Gamma weight on the
same family. The finite sufficient test and weighted comparison
provide explicit objects for that task; neither is a proof that
their bounds hold for every width.

CWF's return to the native banks is still for each fixed profile
or fixed finite collection. The larger positive boundary window
does not supply uniform convergence over an entire growing bank.
No new finite-bank positivity size or lower constant for
n_-(a)/sqrt(a) follows here.

## Dependencies and verification scope

- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): exact difference energy,
  boundary decomposition and closed-form compactness.
- [FPW-1](RH_FIRST_PRIME_WINDOW.md): fixed-family RH equivalence,
  Legendre estimate, moment control and elementary constant bounds.
- [CWF-1](RH_CRITICAL_WAVE_FORMS.md): actual source matching and
  the normalized Guinand-Weil form.
- [NIST DLMF 5.9.16](https://dlmf.nist.gov/5.9.E16): the retained
  primary digamma identity underlying the Gamma energy.

All new bounds above have written analytic proofs. The finite
endpoint reduction concerns a sufficient lower bound, and the
weighted comparison is proved from the square inequality (21).
Repository validation concerns text structure, links, exact committed
content and draft status. No scientific test, source computation,
numerical spectrum or independent replay is requested or performed.
