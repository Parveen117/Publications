# DWE-1 — A quantified dilation extension of the positive source window

Monty · 18 October 2026.

Written continuation at e324197d8d901be97f397cc65582360169bb3e8e.
This note uses the actual compact-test source and GRR's established
quarter-window margin. It proves the small unrestricted extension
\[
\boxed{
0<L\le L_+=\frac{11901}{47600}=\frac14+\frac1{47600}
\quad\Longrightarrow\quad
q_L[g]\ge\frac1{1700}\|g\|_2^2.
}
\]
No mean-zero or cosh-zero restriction is imposed. The new endpoint
is still below the first-prime threshold. This is a quantitative
continuity consequence of a previously positive margin, not a
global positivity argument or a proof of RH.

## Source, domain and retained input

Use the real even compact-test source of
[PTE](RH_PRIME_TRANSPORT_ENERGY.md),
[OPC](RH_FIRST_PRIME_ONE_CHANNEL.md) and
[GRR](RH_GAMMA_REMAINDER_RETURN.md):
\[
\begin{aligned}
q_L[g]={}&\ell(g)^2+\frac{\Omega_0}{2}N_g
+\mathcal E_\Gamma[g]
-\sum_{\log n<2L}\frac{\Lambda(n)}{\sqrt n}C_g(\log n),\\
N_g={}&\|g\|_2^2,\qquad
C_g(t)=\int_{\mathbb R}g(x)g(x-t)\,dx,\\
\ell(g)={}&\int_{\mathbb R}g(x)\cosh(x/2)\,dx,\\
\mathcal E_\Gamma[g]={}&
\int_0^\infty k(t)(N_g-C_g(t))\,dt,\qquad
k(t)=\frac{e^{-t/2}}{1-e^{-2t}},\\
\Omega_0={}&-\gamma_E-\frac\pi2-3\log2-\log\pi .
\end{aligned}
\tag{1}
\]
Functions are supported in (-L,L), extended by zero.
The domain D_L is the closure of the real even smooth compact
core in the retained source form norm.

GRR-T4 proves
\[
q_{1/4}[f]\ge\frac1{850}\|f\|_2^2
\quad(f\in D_{1/4}).
\tag{2}
\]
No all-width sign premise is used below.

The dilation parameter r is a ratio of support widths. It is
not the native anchor, lambda, derivative depth or physical time.

## DWE-T1. A global Gamma-energy dilation loss bound

Let r>=1 and define the norm-preserving dilation
\[
(U_rf)(x)=r^{-1/2}f(x/r).
\]
For every smooth compact real f,
\[
\boxed{
\mathcal E_\Gamma[U_rf]\ge
\mathcal E_\Gamma[f]-5(r-1)\|f\|_2^2 .
}
\tag{3}
\]

**Proof.** Direct substitution gives
\[
\|U_rf\|_2^2=N_f,\qquad C_{U_rf}(u)=C_f(u/r).
\]
Consequently
\[
\mathcal E_\Gamma[U_rf]-\mathcal E_\Gamma[f]
=\int_0^\infty(rk(rt)-k(t))(N_f-C_f(t))\,dt .
\tag{4}
\]
For t>0, integration on (0,2t) gives
\[
1-e^{-2rt}=r\int_0^{2t}e^{-rv}\,dv
\le r\int_0^{2t}e^{-v}\,dv=r(1-e^{-2t}).
\]
Thus
\[
rk(rt)\ge e^{-(r-1)t/2}k(t),\qquad
rk(rt)-k(t)\ge-\frac{r-1}{2}t k(t).
\tag{5}
\]
Cauchy--Schwarz gives |C_f(t)|<=N_f, hence
0<=N_f-C_f(t)<=2N_f. Multiplication of (5) by this
nonnegative factor and integration yields the loss
\[
(r-1)N_f\int_0^\infty t k(t)\,dt.
\]
The positive geometric series and monotone integration give
\[
\begin{aligned}
\int_0^\infty t k(t)\,dt
&=\sum_{n=0}^\infty\frac1{(2n+\tfrac12)^2}\\
&<4+\int_0^\infty\frac{dx}{(2x+\tfrac12)^2}=5.
\end{aligned}
\tag{6}
\]
The strict sum comparison follows because the summand's
continuous extension is strictly decreasing; each n>=1 term
is smaller than its integral over (n-1,n).

There is no subtraction of divergent integrals in (4):
k(t)=1/(2t)+O(1) near zero, so rk(rt)-k(t)=O(1);
both tails decay exponentially. Alternatively derive (4)
with finite cutoffs and then take their limits. This proves (3).
QED.

## DWE-T2. A uniform pole-square loss on the reference interval

For f supported in (-1/4,1/4) and 1<=r<=101/100,
\[
\boxed{
\ell(U_rf)^2\ge\ell(f)^2-2(r-1)\|f\|_2^2.
}
\tag{7}
\]

**Proof.** On I_0=(-1/4,1/4), set
\[
v_r(y)=\sqrt r\,\cosh(ry/2).
\]
Then ell(U_rf)=<f,v_r>. We verify
\[
\|v_r\|_{L^2(I_0)}\le1,\qquad
\|\partial_rv_r\|_{L^2(I_0)}\le1 .
\tag{8}
\]
Indeed |ry/2|<=101/800<1/7 and sqrt(r)<11/10.
For 0<=z<=1/7, the positive series give
\[
\cosh z\le\frac1{1-z^2/2}\le\frac{98}{97}<\frac65,
\qquad
\sinh z\le\frac{z}{1-z^2}\le\frac7{48}<\frac15.
\]
For the first inequality one may use (2n)!>=2^n;
for the second, discard the factorial denominators in the
positive odd series. Thus |v_r|<33/25 and
\[
\|v_r\|_2^2<\frac12(33/25)^2
=\frac{1089}{1250}<1.
\]
Moreover
\[
\partial_rv_r
=\frac{\cosh(ry/2)}{2\sqrt r}
+\sqrt r\,\frac y2\sinh(ry/2),
\]
so
\[
|\partial_rv_r|
\le\frac35+\frac{11}{10}\frac18\frac15
=\frac{251}{400}<1.
\]
Since I_0 has length 1/2, (8) follows. Integration of the
L2 derivative bound gives ||v_r-v_1||_2<=r-1. Therefore
\[
\begin{aligned}
|\ell(U_rf)^2-\ell(f)^2|
&=|\langle f,v_r-v_1\rangle
       \langle f,v_r+v_1\rangle|\\
&\le N_f\|v_r-v_1\|_2(\|v_r\|_2+\|v_1\|_2)
\le2(r-1)N_f.
\end{aligned}
\]
This proves (7). QED.

## DWE-T3. The full positive interval extends by 1/47600

For 1<=r<=101/100 and every real even smooth f supported
in (-1/4,1/4),
\[
\boxed{
q_{r/4}[U_rf]\ge
\left(\frac1{850}-7(r-1)\right)\|f\|_2^2.
}
\tag{9}
\]
Consequently
\[
\boxed{
q_L[g]\ge\frac1{1700}\|g\|_2^2
\quad\left(0<L\le\frac{11901}{47600},\ g\in D_L\right).
}
\tag{10}
\]

**Proof.** On the entire stated dilation range,
\[
2L=r/2\le101/200<6931/10000<\log2,
\]
using the retained elementary log2 bound in OPC-T2.
All prime-power correlation terms in (1) therefore vanish
by support, both before and after dilation. No active signed
prime term is dropped.

The constant term in (1) is unchanged because the norm is
preserved. Add (3) and (7), then use (2), to prove (9).

Choose
\[
r_+=1+\frac1{11900}<\frac{101}{100},\qquad
L_+=r_+/4=\frac{11901}{47600}.
\]
The remaining coefficient is
\[
\frac1{850}-\frac7{11900}
=\frac1{850}-\frac1{1700}=\frac1{1700}.
\]
Every real even smooth test supported in (-L_+,L_+) is
U_(r_+)f for a real even smooth test supported in (-1/4,1/4).
Thus the estimate covers the whole endpoint core, not just a
selected test family. It extends to D_(L_+) by the defining
form closure.

For smaller L, zero extension leaves the norm, pole, correlations
and actual source value unchanged. Apply the endpoint inequality
on the smaller core and then close in D_L. This proves (10).
QED.

The original stronger 1/850 margin remains valid on L<=1/4.
Equation (9) also gives a varying margin on the slightly larger
range where its displayed coefficient is positive. The endpoint
in (10) is selected to retain half the reference margin, not
asserted to be optimal.

## What this does and does not settle

| Statement | Status |
|---|---|
| Unrestricted real even source through L=1/4 | Retained margin 1/850 |
| Unrestricted source through L=1/4+1/47600 | Proved margin 1/1700 |
| Gamma-energy dilation loss | Global bound (3) |
| Complete-source dilation bound (9) | Restricted to the stated prime-free range |
| First-prime returned scalar at L=7/20 | Sign remains unproved |
| All-width returned-matrix positivity and RH | Remain unproved |

The quantitative advance is small. It uses the existing positive
margin as a finite budget; repeating a continuity argument does
not replenish that budget or prove unlimited positive extension.

The preceding chat also derived a boundary-weighted residual
comparison from OPC/GRR. No evaluated sign certificate for it
was obtained. That conditional comparison is not used to infer
(10), and this note does not claim it closes the first-prime gate.

## Sources and verification scope

- [GRR-1](RH_GAMMA_REMAINDER_RETURN.md), especially T4: the actual full-family reference margin.
- [OPC-1](RH_FIRST_PRIME_ONE_CHANNEL.md): source conventions, support thresholds and the still-open first-prime scalar.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): the source form and its closed domain.
- [CLA-1](RH_CLOSURE_AUDIT.md): conservation alone does not establish the missing source inequality.

The dilation, elementary kernel comparison, integral bound and
rank-one pole estimate are standard analytic tools applied here
to the retained source. No literature-wide novelty is asserted.

Verification is written derivation review plus exact document
readback and repository metadata checks. No numerical source
evaluation, scientific computation, independent agent, certificate
replay or workflow dispatch is performed. Only this note and its
README entry change. The commit uses [skip ci]; the PR is a
draft stacked on #136 and remains unmerged.
