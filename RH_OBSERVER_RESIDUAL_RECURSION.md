# ORR-1 — Exact observer recursion and inverse-residual curvature

18 October 2026 (Asia/Kolkata). Author: Monty.
Written continuation of DSG-1 at
b2914974b7eedaad64845f648317f7b98664f188.

## Result and scope

DSG gives the exact scalar sign test for adjoining one actual observer
to a fixed positive bank. This note determines what changes after
that observer has been accepted.

The complete filtered residual has an exact rational update.
The two mixed-observer kernels also have exact formulas in that
residual. Their Schur complements determine the loss of each
positivity margin at the next addition. Separate positive diagonal
readings do not establish positivity of their joint bank.

Taking the next observer to the newly accepted position gives
explicit second- and third-derivative conditions for the inverse
residual. On UCB's existing two-extra-position range these conditions
have strictly positive quantitative margins. In particular the inverse
residual is increasing, strictly concave, and has positive third
derivative throughout the observation interval.

These are additional actual-source identities and derivative bounds.
The finite bank-size range is unchanged. No unrestricted dense
positive extension, all-order derivative law, or RH proof is claimed.

## 1. Retained source and notation

Keep the actual source and the finite bank of DSG/RCD:
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,
\qquad A=\{t_1,\ldots,t_m\},\quad t_i>9/4,
\]
\[
L_{ij}=-q[t_i,t_j],\quad J_{ij}=(tq(t))[t_i,t_j],
\quad u_i=q(t_i),\qquad L,J\succ0.
\tag{1}
\]
Divided differences at coincidence mean their derivative limits.
SE48 supplies
\[
P(z)=J+zL,\quad v(z)=P(z)^{-1}u,\quad
q_A(z)=u^TP(z)^{-1}u
      =\sum_{j=1}^m\frac{r_j}{z+\eta_j},\quad r_j,\eta_j>0.
\]
Set
\[
d_A(z)=\prod_j(z+\eta_j),\quad
\omega_A(z)=\prod_i(z-t_i),\quad H(z)=H_A(z)=\omega_A(z)/d_A(z).
\]
RCD gives the retained analytic residual
\[
q(z)-q_A(z)=H(z)^2T(z),\qquad T(z)=\mathcal T_A(z)
 =\sum_\mu\frac{\nu_\mu h_A(\mu)^2}{z+\mu}.
\tag{2}
\]
Here \(\mu\) is the actual, possibly nonreal, source value and
\(h_A(\mu)=\prod_j(\mu-\eta_j)/\prod_i(\mu+t_i)\).
Its square is bilinear, not an absolute square.
The model poles \(\eta_j\) are not identified with actual xi zeros.

For positive \(z\), define the smooth vector
\[
b(z)=b_A(z),\qquad
b_i(z)=\frac{\prod_{k\ne i}(z-t_k)}{d_A(z)}.
\]
Away from A this is \(H(z)((z-t_i)^{-1})_i\).
Use the mixed coefficients
\[
K_L(x,y)=b(x)^TL^{-1}b(y),\qquad
K_J(x,y)=b(x)^TJ^{-1}b(y).
\tag{3}
\]
Their diagonal values agree with DSG's \(K_L,K_J\).
The two normalized scalar innovations are
\[
\mathfrak a(x)=-T'(x)-K_L(x,x)T(x)^2,
\]
\[
\mathfrak b(x)=T(x)+xT'(x)-x^2K_J(x,x)T(x)^2.
\tag{4}
\]
Thus the actual innovations are
\(\alpha_A(x)=H(x)^2\mathfrak a(x)\) and
\(\beta_A(x)=H(x)^2\mathfrak b(x)\).

Every derivative below is in the source position with the stated bank
fixed. In an update, x is fixed while z varies. An approach y to x
means a second, distinct candidate position tending to x; it does not
differentiate a moving bank. None of these derivatives is a physical
time, wave-coordinate or layer-count derivative.

## ORR-T1. The complete residual has an exact one-observer update

Let \(x\notin A\) be accepted, so \(\mathfrak a(x),\mathfrak b(x)>0\),
and set \(B=A\cup\{x\}\). DSG then gives \(T(x)>0\).
For real \(z>0\), define
\[
C_A(x;z)=b(x)^TP(z)^{-1}b(x),
\]
\[
\boxed{
\sigma_A(x;z)=T(x)+(x-z)T'(x)
             -(x-z)^2T(x)^2C_A(x;z).
}
\tag{5}
\]
Then
\[
\boxed{\sigma_A(x;z)\ge\mathfrak b(x)+z\mathfrak a(x)>0,}
\tag{6}
\]
and
\[
\boxed{
q_B(z)-q_A(z)=\frac{H(z)^2T(x)^2}{\sigma_A(x;z)}.
}
\tag{7}
\]
The model denominator and residual factors update by
\[
\boxed{
d_B(z)=d_A(z)\frac{\sigma_A(x;z)}{\mathfrak a(x)},\qquad
H_B(z)=\frac{\mathfrak a(x)(z-x)H(z)}{\sigma_A(x;z)},
}
\tag{8}
\]
\[
\boxed{
\mathcal T_B(z)=
\frac{\sigma_A(x;z)}
     {\mathfrak a(x)^2(z-x)^2}
\left[\sigma_A(x;z)T(z)-T(x)^2\right].
}
\tag{9}
\]
Equation (9) has its analytic removable value at z=x.
These formulas extend as meromorphic identities at other regular
arguments; (6) is asserted for real positive z.

**Proof.** Let \(R_x=H(x)^2T(x)\) and
\(d_x=((x-t_i)^{-1})_i\).
The actual new mixed columns are
\[
l=Lv(x)-R_xd_x,\qquad
j=Jv(x)+xR_xd_x.
\]
In SE49's block update put
\[
k(z)=j+zl=P(z)v(x)+(x-z)R_xd_x,
\]
\[
s(z)=q(x)+(x-z)q'(x)-k(z)^TP(z)^{-1}k(z),
\quad
e(z)=q(x)-k(z)^Tv(z).
\]
The rational residual-observer identity in RCD gives
\[
d_x^Tv(z)=\frac{H(z)/H(x)-1}{z-x};
\tag{10}
\]
its value at z=x is \(H'(x)/H(x)\).
Also
\[
v(x)^TP(z)v(x)=q_A(x)+(x-z)q_A'(x).
\]
Consequently the cross terms containing \(H'(x)\) cancel:
\[
s(z)=H(x)^2\sigma_A(x;z),\qquad
e(z)=H(x)H(z)T(x).
\tag{11}
\]
SE49 proves \(s(z)\ge\beta_A(x)+z\alpha_A(x)\) and
\(q_B-q_A=e^2/s\). Equations (6)--(7) follow.

For the determinants, block elimination gives
\[
\det P_B(z)=\det P(z)s(z),\qquad
\det L_B=\det L\,\alpha_A(x).
\]
Since \(d_A=\det P/\det L\), the first formula in (8) follows.
The relation \(\omega_B(z)=(z-x)\omega_A(z)\) gives the second.
Subtract (7) from (2), then divide by \(H_B^2\), to obtain (9).

At z=x, \(\sigma=T(x)\) and
\(\partial_z\sigma=-T'(x)\). Hence the bracket in (9) and its
first derivative both vanish there. Its double zero provides the
removable continuation; no value is assigned arbitrarily. QED.

This update retains the entire actual residual. Positivity of its
denominator proves a nonnegative model increment. It does not
prove that the updated residual is positive at every test position.

An exact value at the new fitted point follows from the same formula.
Where \(V=1/T\), put \(c_A(x)=C_A(x;x)>0\). Taylor expansion gives
\[
\boxed{
\mathcal T_B(x)=
-\frac{T(x)^4}{\mathfrak a(x)^2}
\left[\frac12V''(x)+c_A(x)\right].
}
\tag{12}
\]
Indeed the coefficient of \((z-x)^2\) in the bracket in (9) is
\[
\frac12T(x)T''(x)-T'(x)^2-T(x)^3c_A(x)
=-T(x)^3\left[\frac12V''(x)+c_A(x)\right].
\]
Thus positivity of this next filtered reading requires a strictly
negative inverse-residual curvature after the stated finite correction.
Acceptance of B alone does not establish that further positivity.

## ORR-T2. Mixed observer signs are divided-difference kernel signs

For candidates x,y outside A, let \(l(x),j(x)\) denote their actual
mixed columns with A, as in T1. Define the two conditional source
kernels before either candidate is added:
\[
C_L^A(x,y)=-q[x,y]-l(x)^TL^{-1}l(y),
\]
\[
C_J^A(x,y)=(zq(z))[x,y]-j(x)^TJ^{-1}j(y).
\]
Then
\[
\boxed{
C_L^A(x,y)=H(x)H(y)\,\mathcal A_A(x,y),\quad
\mathcal A_A(x,y)=-T[x,y]-T(x)T(y)K_L(x,y),
}
\tag{13}
\]
\[
\boxed{
C_J^A(x,y)=H(x)H(y)\,\mathcal B_A(x,y),\quad
\mathcal B_A(x,y)=(zT(z))[x,y]-xyT(x)T(y)K_J(x,y).
}
\tag{14}
\]
On the diagonal these give exactly \(\mathfrak a,\mathfrak b\).

For any finite set Y of distinct candidates outside A,
the complete union \(A\cup Y\) is positive if and only if both
matrices \((\mathcal A_A(x,y))_{x,y\in Y}\) and
\((\mathcal B_A(x,y))_{x,y\in Y}\) are positive definite.
This is a finite equivalence; their positivity is not presumed.

**Proof.** The model kernels equal \(v(x)^TLv(y)\) and
\(v(x)^TJv(y)\). Write \(R_z=H(z)^2T(z)\) and expand the
actual columns \(l(z)=Lv(z)-R_zd_z\),
\(j(z)=Jv(z)+zR_zd_z\).
For the L kernel the remaining terms are
\[
-R[x,y]+R_yv(x)^Td_y+R_xd_x^Tv(y)
           -R_xR_y d_x^TL^{-1}d_y.
\]
Equation (10) and its x,y interchange reduce the first three
terms to \(-H(x)H(y)T[x,y]\).
The last term is
\(-H(x)H(y)T(x)T(y)K_L(x,y)\).
This proves (13). The same expansion for J reduces its first
three terms to \(H(x)H(y)(zT(z))[x,y]\), and gives (14).
Coincident limits follow by analyticity.

For a finite Y these kernels are the two Schur-complement matrices
of the old positive blocks. Their congruence factors are the
invertible diagonal matrix with entries H(x). The usual Schur
criterion proves the claimed equivalence. QED.

Where \(T>0\), set
\[
V=1/T,\qquad F=-V/z.
\]
Define the real symmetric kernels
\[
\boxed{
D_L^A(x,y)=V[x,y]-K_L(x,y),\qquad
D_J^A(x,y)=F[x,y]-K_J(x,y).
}
\tag{15}
\]
Direct division of (13)--(14) gives
\[
C_L^A(x,y)=H(x)T(x)H(y)T(y)D_L^A(x,y),
\]
\[
C_J^A(x,y)=xH(x)T(x)yH(y)T(y)D_J^A(x,y).
\tag{16}
\]
For example,
\(-T[x,y]/(T(x)T(y))=V[x,y]\), while
\[
\frac{(zT(z))[x,y]}{xyT(x)T(y)}
=\frac{V(y)/y-V(x)/x}{x-y}=F[x,y].
\]
Thus the full mixed gate is positivity of two finite
divided-difference matrices after subtraction of the retained
finite Gram matrices.

DSG's scalar window is their diagonal test:
\[
D_L^A(x,x)=V'(x)-K_L(x,x)>0,\quad
D_J^A(x,x)=F'(x)-K_J(x,x)>0.
\tag{17}
\]
The matrices in (15), rather than just (17), retain the mixed
readings required for several new observers.

## ORR-T3. An accepted observer consumes an exact squared mixed margin

Suppose x is accepted and y is another distinct candidate.
Where T is positive, write
\[
p_L=D_L^A(x,x)>0,\qquad p_J=D_J^A(x,x)>0.
\]
The actual innovations relative to \(B=A\cup\{x\}\) are
\[
\boxed{
\alpha_B(y)=H(y)^2T(y)^2
\left[D_L^A(y,y)-\frac{D_L^A(x,y)^2}{p_L}\right],
}
\tag{18}
\]
\[
\boxed{
\beta_B(y)=y^2H(y)^2T(y)^2
\left[D_J^A(y,y)-\frac{D_J^A(x,y)^2}{p_J}\right].
}
\tag{19}
\]
Hence the remaining margin is the old diagonal margin minus the
displayed nonnegative mixed contribution.

**Proof.** Eliminate the old bank first, then eliminate x.
The first conditional block is (16). Its one-dimensional Schur
complement at x is exactly (18) or (19). This is also the Schur
complement obtained by eliminating the whole bank B at once. QED.

If x and y each separately pass the original scalar gate, define
\[
\rho_K(x,y)=
\frac{D_K^A(x,y)}
{\sqrt{D_K^A(x,x)D_K^A(y,y)}}\qquad(K=L,J).
\]
Then
\[
\boxed{
A\cup\{x,y\}\text{ is positive}
\quad\Longleftrightarrow\quad
|\rho_L(x,y)|<1\ \text{and}\ |\rho_J(x,y)|<1.
}
\tag{20}
\]
These are normalized mixed readings. Calling them correlations
does not assume a full-source probability law. That interpretation
would require the joint positivity being tested.

## ORR-T4. The close-observer test is a third-derivative inequality

Fix an accepted \(x\notin A\). In a neighborhood of x, \(T>0\).
For \(K=L,J\), use
\[
f_L=V,\quad M_L=L^{-1},\qquad
f_J=F,\quad M_J=J^{-1}.
\]
At x set
\[
p_K=f_K'-b^TM_Kb>0,\quad
r_K=\frac12 f_K''-b^TM_Kb',\quad
t_K=\frac16 f_K'''-(b')^TM_Kb',
\]
\[
\boxed{\kappa_K=t_K-r_K^2/p_K.}
\tag{21}
\]
Then, as \(y\to x\),
\[
D_K^A(y,y)-\frac{D_K^A(x,y)^2}{p_K}
=(y-x)^2\kappa_K+O((y-x)^3).
\tag{22}
\]
Consequently
\[
\boxed{
\lim_{y\to x}\frac{\alpha_B(y)}{(y-x)^2}
=H(x)^2T(x)^2\kappa_L,
}
\]
\[
\boxed{
\lim_{y\to x}\frac{\beta_B(y)}{(y-x)^2}
=x^2H(x)^2T(x)^2\kappa_J.
}
\tag{23}
\]

**Proof.** For the divided-difference kernel of an analytic f,
its value, first partial derivative, and mixed partial derivative
at (x,x) are \(f'(x),f''(x)/2,f'''(x)/6\).
Subtracting the Gram kernel gives respectively \(p_K,r_K,t_K\).
A Taylor expansion of a symmetric kernel D yields
\[
D(y,y)-D(x,y)^2/D(x,x)
=(y-x)^2\left[\partial_1\partial_2D
-(\partial_2D)^2/D\right]_{(x,x)}
+O((y-x)^3).
\]
The unmixed second derivatives cancel. This proves (22);
(18)--(19) prove (23). QED.

If both \(\kappa_L,\kappa_J>0\), every sufficiently close distinct
y is accepted after x. If either is negative, that channel rejects
every sufficiently close distinct y. A zero value requires higher
terms and is not decided by this test.

In particular the strict local condition in the L channel is
\[
\boxed{
\frac16V'''>(b')^TL^{-1}b'
+\frac{(\frac12V''-b^TL^{-1}b')^2}{V'-b^TL^{-1}b},
}
\tag{24}
\]
with the corresponding F,J inequality for the other channel.
This condition retains the mixed derivative record. A positive
scalar slope alone does not establish it.

## ORR-T5. The existing finite source range proves strict curvature margins

Use exactly UCB's comparison anchor and size budget:
\[
I_a=[a/2,3a/2],\quad
\ell_a=\log(\sqrt a/(2\pi)),\quad
N_*(a)=\left\lfloor\frac{e^{-4}\sqrt a}{4\pi}\right\rfloor\ge128,
\]
\[
M(a)=\lfloor N_*(a)/12\rfloor,\quad
\delta_k(a)=\frac6{\ell_a(441/4)^k},\quad
S_k(v)=\sum_{j=0}^{k-1}v^{2j}.
\]
Assume \(A\subset I_a\) has m distinct positions and
\[
\boxed{m+2\le M(a).}
\tag{25}
\]
This reserves two additional observations. Define, for every
\(x\in I_a\),
\[
\boxed{
E_A(x)=
\frac{q(a)\delta_{m+2}(a)d_A(x)^2}
 {a^{2m+3}S_{m+2}(x/a-1)T(x)^2}>0.
}
\tag{26}
\]
All functions in this formula are the retained actual residual
and its finite bank, with \(T>0\) already proved by DSG.
Then throughout \(I_a\), including old fitted positions,
\[
\boxed{
-\frac12V''(x)-c_A(x)\ge(a+x)E_A(x)>0,
}
\tag{27}
\]
\[
\boxed{
\frac16V'''-(b')^TL^{-1}b'
-\frac{(\frac12V''-b^TL^{-1}b')^2}
       {V'-b^TL^{-1}b}
=\kappa_L\ge E_A(x)>0,
}
\tag{28}
\]
\[
\boxed{
\frac16F'''-(b')^TJ^{-1}b'
-\frac{(\frac12F''-b^TJ^{-1}b')^2}
       {F'-b^TJ^{-1}b}
=\kappa_J\ge\frac{a}{x^2}E_A(x)>0.
}
\tag{29}
\]
Thus
\[
\boxed{V'>0,\qquad V''<0,\qquad V'''>0}
\tag{30}
\]
on this declared finite range. Also \(F'''>0\).
No sign pattern for every derivative is asserted.

**Proof.** First take \(x\in I_a\setminus A\).
UCB guarantees that \(B=A\cup\{x\}\) is positive.
Its size is m+1, and (25) leaves one further position available.
DSG-T4 therefore applies to B, including at its fitted point x:
\[
\mathfrak a_B(x)\ge\varepsilon_L^B(x),\qquad
\mathfrak b_B(x)\ge\varepsilon_J^B(x)=a\varepsilon_L^B(x),
\]
where
\[
\varepsilon_L^B(x)=
\frac{q(a)\delta_{m+2}(a)d_B(x)^2}
 {a^{2m+3}S_{m+2}(x/a-1)}.
\tag{31}
\]
From (8),
\[
d_B(x)=d_A(x)\frac{T(x)}{\mathfrak a(x)},\qquad
H_B'(x)=\frac{\mathfrak a(x)H(x)}{T(x)}.
\tag{32}
\]
DSG's confluent innovations give
\[
\lim_{y\to x}\frac{\alpha_B(y)}{(y-x)^2}
=H_B'(x)^2\mathfrak a_B(x),
\]
and the same equality with beta and \(\mathfrak b_B\).
Comparison with (23), then substitution of (31)--(32), gives
\[
\kappa_L\ge
\frac{\mathfrak a(x)^2}{T(x)^4}\varepsilon_L^B(x)=E_A(x),
\]
\[
\kappa_J\ge
\frac{\mathfrak a(x)^2}{x^2T(x)^4}\varepsilon_J^B(x)
=\frac a{x^2}E_A(x).
\]
This proves (28)--(29).

For (27), the residual sum identity for B is
\[
\mathcal T_B(x)
=x\mathfrak a_B(x)+\mathfrak b_B(x)
+\left(xK_L^B(x,x)+x^2K_J^B(x,x)\right)\mathcal T_B(x)^2.
\]
All terms in the last coefficient are nonnegative. Hence
\(\mathcal T_B(x)\ge(x+a)\varepsilon_L^B(x)\).
Insert (12) and use (31)--(32) again to obtain (27).

The functions \(T,V,F,b,c_A,E_A\) are smooth on the whole
interval, and the denominators in (28)--(29) have DSG's strict
positive margins there. Therefore these inequalities extend
continuously to \(x\in A\). Their limiting statements do not
make duplicated raw observer rows independent.

Finally DSG gives \(V'>0\); (27) gives \(V''<0\);
and the nonnegative terms on the right of (28), together with
\(E_A>0\), give \(V'''>0\). Equation (29) similarly gives
\(F'''>0\). QED.

The derivative bounds are full-source statements on the already
proved finite range. Their constants are explicit symbolic
expressions; no source value or derivative was numerically evaluated.

## Remaining source-identification problem

The scalar window now has an exact update, and the mixed signs
have an explicit local curvature test. Equations (18)--(19)
show precisely why acceptance at one stage does not propagate
for free: the next stage subtracts a squared mixed term.

Within (25), UCB controls that subtraction and the new curvature
margins are strictly proved. Repeating a limiting argument does
not reset the total-size budget. No inequality above has been
established uniformly after that budget is exhausted.

The retained sufficient route to RH remains a nested sequence
of positive actual banks whose positions tend to infinity and
satisfy
\[
\sum_n t_n^{-1/2}=\infty.
\]
NOE/DSG then identify the model limit with the actual source;
PWC supplies faithful-wave positivity. The new formulas determine
the residual and mixed quantities that such a construction must
control. They do not construct the required infinite sequence.

## Sources, attribution, and verification

- [DSG-1](RH_DENSE_OBSERVER_SLOPE_GATE.md), T2--T4:
  normalized scalar innovations, inverse-residual slope window,
  and finite confluent margins.
- [RCD-1](RH_RESIDUAL_COMPLETION.md), T1--T3:
  complete source residual, rational residual-observer identity,
  and fitted-point continuation.
- [UCB-1](RH_OBSERVER_BANK_UNION.md), T2--T3:
  actual union positivity and its unchanged finite size budget.
- [SE48](papers/native-prime-observation-closure/sources/SE48.md),
  T1: finite positive representation and interpolation.
- [SE49](papers/native-moving-zero-covariance/sources/SE49.md),
  T1: exact source Schur signs and response update.
- [NOE-1](RH_NESTED_OBSERVER_EXTENSION.md), T5--T7, and
  [PWC-1](RH_PROBABILITY_WAVE_CLOSURE.md), T3:
  model limits, the source-identification requirement, and the
  conditional faithful-wave conclusion.

Schur complementation, divided-difference kernels, rational
interpolation and Taylor limits are standard tools. The additions
here are their exact filtered-residual recursion for the retained
source and the quantitative curvature bounds on its declared range.
No literature-wide novelty is claimed.

Verification consists of the written algebra and analytic limits
above, followed by text and repository metadata checks. No scientific
computation, numerical source evaluation, independent agent,
certificate replay or workflow dispatch is performed.
Only this note and its README entry change. The commit uses
[skip ci]. Published as an unmerged draft stacked on #133.
