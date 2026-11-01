# PRW-1 — Pole recovery and the unrestricted positive window

17 October 2026. Written continuation at Publications commit
6ccef2211e492b86445c1d0a26ac2f8ca22cfefe, following
[GCC-1](RH_GAMMA_CUT_CERTIFICATES.md).

The next step is an actual source inequality. Keeping the pole
square and controlling the constant mode gives
\[
\boxed{
0<L\le\frac15
\quad\Longrightarrow\quad
q_L[g]\ge\frac{22673}{6250000}\|g\|_2^2
\ge\frac1{300}\|g\|_2^2
\quad(g\in D_L).
}
\]
The second inequality is strict for nonzero g.
Here D_L is the full even source form domain. No zero-mean or
cosh-neutrality condition is imposed. The established unrestricted
even range therefore increases from \(L\le1/32\) to \(L\le1/5\).
The earlier stronger margin remains available in its smaller range.

The proof uses the retained Gamma kernel, its boundary weight,
the ordinary even Legendre energy and the actual pole square.
It reduces the possible constant-mode loss to a rank-one scalar
inequality. A prime-aware version of that inequality is also given
for later windows.

The technical notes at the end record the requested relational
viewpoint and angular probability coordinates. Normalized cut
weights can be written as \(\cos^2\theta,\sin^2\theta\);
phase and source cross-terms remain explicit.

**Status: RH remains unproved.** The new unrestricted range is
still before the first prime threshold \(L=(\log2)/2\).
FCP's larger two-prime range remains a result on its two-moment
restricted family. No all-width positivity, new uniform native
finite-bank size or numerical certificate run is claimed.

## 1. Full even family and the actual source

For \(I=(-L,L)\), let H_L be the real Hilbert space of even
L2 functions on I, extended by zero. Retain
\[
\begin{gathered}
N_g=\|g\|_2^2,\qquad C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,\\
\ell(g)=\int_I g(x)\cosh(x/2)\,dx,\qquad
a_n=\frac{\Lambda(n)}{\sqrt n},\\
\mathcal N_L=\{n\ge2:n\text{ is a prime power},\ \log n<2L\},\\
q_L[g]=\ell(g)^2+\frac{\Omega_0}{2}N_g+
\mathcal E_\Gamma[g]-\sum_{n\in\mathcal N_L}a_nC_g(\log n),\\
\Omega_0=-\gamma_E-\frac{\pi}{2}-3\log2-\log\pi,\\
\mathcal E_\Gamma[g]
=\int_0^\infty k(u)(N_g-C_g(u))\,du,\qquad
k(u)=\frac{e^{u/2}}{2\sinh u},\\
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}).
\end{gathered}
\tag{1}
\]
D_L is PTE's closed Gamma form domain, equivalently the maximal
finite-energy domain identified in GCC-T1. Smooth real even
compact tests form its defining core.

The constant channel and the boundary weight are
\[
\phi_0(x)=\frac1{\sqrt{2L}},\qquad
a_0(g)=\langle g,\phi_0\rangle,\qquad
w_L(x)=-\frac14\log(1-(x/L)^2).
\tag{2}
\]
The symbol \(a_0(g)\) is a profile coefficient, not a prime weight.
The proof keeps this coefficient; it does not set it to zero.

## PRW-T1. The pole square compensates part of a rank-one constant-mode loss

For \(0<L\le1/2\) and \(0<\alpha<1\), define
\[
\begin{aligned}
e_{L,\alpha}
&=(\alpha^{-1}-1)\frac{L^5}{160}\cosh^2(L/2),\\
\gamma_{L,\alpha}
&=\frac{\Omega_0}{2}+F(2L)+\frac12\log2
+\frac{2L}{5}+\frac34-e_{L,\alpha},\\
\beta_{L,\alpha}
&=\frac34+\frac L5-2L(1-\alpha).
\end{aligned}
\tag{3}
\]
For every smooth real even g supported in I,
\[
\boxed{
q_L[g]\ge
\int_I(\gamma_{L,\alpha}+w_L(x))|g(x)|^2\,dx
-\beta_{L,\alpha}|a_0(g)|^2
-\sum_{n\in\mathcal N_L}a_nC_g(\log n).
}
\tag{4}
\]

**Proof.** FPW-T2 proves
\[
k(u)\ge\frac1{2u}+\frac15\qquad(0<u\le1).
\]
Apply TPW-T1's exact polynomial-kernel identity with
\(P(u)=1/5\). Its correlation residual is
\[
\mathcal R_P[g]
=\frac1{10}\left(\int g\right)^2
=\frac L5|a_0(g)|^2.
\]
Consequently
\[
\begin{aligned}
\frac{\Omega_0}{2}N_g+\mathcal E_\Gamma[g]
\ge{}&
\left(\frac{\Omega_0}{2}+F(2L)+\frac12\log2+
\frac{2L}{5}\right)N_g\\
&+J_L[g]+\int_Iw_L(x)|g(x)|^2\,dx
-\frac L5|a_0(g)|^2,
\end{aligned}
\tag{5}
\]
where
\[
J_L[g]=\frac18\int_I\int_I
\frac{|g(x)-g(y)|^2}{|x-y|}\,dy\,dx.
\]
The Legendre energy identity gives
\[
J_L[g]\ge\frac34\bigl(N_g-|a_0(g)|^2\bigr).
\tag{6}
\]
Indeed the constant mode has eigenvalue zero, odd modes are absent
by evenness, and every remaining even mode has
\(H_{2j}/2\ge H_2/2=3/4\).

To retain the pole square, write
\[
\ell(g)=\sqrt{2L}\,a_0(g)+
\int_I g(x)(\cosh(x/2)-1)\,dx .
\]
The weighted square inequality
\[
(A+B)^2\ge(1-\alpha)A^2-(\alpha^{-1}-1)B^2
\]
holds for real A,B. Taylor's bound gives
\[
|\cosh(x/2)-1|
\le\frac{x^2}{8}\cosh(L/2),\qquad
\|\cosh(x/2)-1\|_{L^2(I)}^2
\le\frac{L^5}{160}\cosh^2(L/2).
\]
Thus
\[
\boxed{
\ell(g)^2\ge
2L(1-\alpha)|a_0(g)|^2-e_{L,\alpha}N_g.
}
\tag{7}
\]
Combine (5)--(7) with the complete source (1) to obtain (4).
No moment is assumed to vanish. QED.

The pole is not replaced by a constant-mode square without error:
\(e_{L,\alpha}\) retains the difference between its cosh profile
and the constant profile. The parameter alpha controls this
estimate; it is not a probability or an observation angle.

## PRW-T2. A scalar rank-one gate includes all active prime translations

Suppose \(\gamma_{L,\alpha}>0\). Choose a positive weight v with
\[
0<v(x)\le\gamma_{L,\alpha}+w_L(x).
\]
For \(u_n=\log n\), define the weighted prime row
\[
R_v(x)=\frac12\sum_{n\in\mathcal N_L}a_n
\left(
\frac{\mathbf1_I(x-u_n)}{\sqrt{v(x)v(x-u_n)}}+
\frac{\mathbf1_I(x+u_n)}{\sqrt{v(x)v(x+u_n)}}
\right),
\tag{8}
\]
evaluating a term only when its shifted point lies in I.
Set
\[
d(x)=(1-R_v(x))v(x),\qquad
m=\operatorname*{ess\,inf}_{x\in I}d(x),\qquad
J_d=\frac1{2L}\int_I\frac{dx}{d(x)}.
\tag{9}
\]
If
\[
m>0,\qquad
\beta_+J_d<1,\qquad
\beta_+=\max\{\beta_{L,\alpha},0\},
\tag{10}
\]
then
\[
\boxed{
q_L[g]\ge(1-\beta_+J_d)m\|g\|^2
\quad(g\in D_L).
}
\tag{11}
\]
The stated weights are assumed to give finite integrals on the
smooth source core.

**Proof.** BPC's weighted pair inequality, which does not itself
require neutral moments, gives
\[
\sum_{n\in\mathcal N_L}a_n|C_g(u_n)|
\le\int_I R_v(x)v(x)|g(x)|^2\,dx.
\]
Equation (4) therefore implies
\[
q_L[g]\ge\int_I d(x)|g(x)|^2\,dx
-\beta_+|a_0(g)|^2.
\]
Weighted Cauchy--Schwarz gives
\[
|a_0(g)|^2
\le
\left(\int_I\frac{|\phi_0(x)|^2}{d(x)}\,dx\right)
\left(\int_I d(x)|g(x)|^2\,dx\right)
=J_d\int_I d(x)|g(x)|^2\,dx .
\]
Use (10) and \(d\ge m\) to obtain (11) on the core.
Its uniform lower bound extends to D_L by form closure. QED.

This is a sufficient source inequality for later prime windows.
The prime terms must actually satisfy (8)--(10); the rank-one
reduction alone does not establish that. If no prime is active,
\(R_v=0\) and the gate uses just \(J_v\).

## PRW-T3. Rational boundary constants at \(L=1/5\)

The retained constants and the new boundary estimate are
\[
\boxed{
\frac{\Omega_0}{2}>-\frac{26871}{10000},\qquad
\log2>\frac{6931}{10000},\qquad
F(2/5)>\frac{4597}{2500}.
}
\tag{12}
\]
The first two are proved in FCP-T3 and TPW-T3.
Here is an elementary proof of the third, including the additional
constants it uses.

**A lower arctangent constant.** The addition identity gives
\(\pi/4=\arctan(1/2)+\arctan(1/3)\).
Alternating series lower bounds give
\[
\begin{aligned}
\arctan(1/2)
&>
\frac12-\frac1{24}+\frac1{160}-\frac1{896}
+\frac1{4608}-\frac1{22528}
>\frac{4636}{10000},\\
\arctan(1/3)
&>
\frac13-\frac1{81}+\frac1{1215}-\frac1{15309}
>\frac{3217}{10000}.
\end{aligned}
\]
Hence
\[
\frac{\pi}{4}>\frac{7853}{10000}.
\tag{13}
\]
Also
\[
\log(5/4)=2\operatorname{artanh}(1/9)
>\frac29+\frac2{2187}=\frac{488}{2187}
>\frac{2231}{10000}.
\]
Together with the lower bound for log2, this yields
\[
\log10=3\log2+\log(5/4)>\frac{1439}{625}.
\tag{14}
\]

**A bounded hyperbolic parameter.** Let \(z=\tanh(1/10)\).
For
\[
p(x)=x-\frac{x^3}{3}+\frac{2x^5}{15},
\]
direct differentiation gives
\[
p'(x)-(1-p(x)^2)
=\frac{x^6}{225}(85-20x^2+4x^4)>0\quad(x>0).
\]
Since \((\tanh x)'=1-\tanh^2x\), subtracting the equations
and using an integrating factor proves \(p(x)>\tanh x\)
for x>0. Thus
\[
z<\frac{24917}{250000}
=\frac1{10}\left(1-\frac{83}{25000}\right).
\tag{15}
\]
The cosh series gives
\(\cosh(1/10)\le200/199<101/100\), while
\(\sinh(1/10)>1/10\). Therefore
\[
z>\frac{10}{101}>\frac{99}{1000}.
\tag{16}
\]

Write \(r=e^{-1/5}\). The relation
\(z=(1-r)/(1+r)\) gives the exact identity
\[
F(2/5)=-\frac12\log z+\frac{\pi}{4}-\arctan z.
\tag{17}
\]
By (15), \(-\log(1-t)>t\), and (14),
\[
-\frac12\log z>
\frac{1439}{1250}+\frac{83}{50000}.
\]
For \(0<z<1\), the alternating upper bound
\(\arctan z\le z-z^3/3+z^5/5\), together with (15)--(16),
therefore gives
\[
\begin{aligned}
F(2/5)>{}&
\frac{1439}{1250}+\frac{83}{50000}
+\frac{7853}{10000}
-\frac{24917}{250000}\\
&+\frac13\left(\frac{99}{1000}\right)^3
-\frac1{500000}\\
={}&\frac{4597}{2500}+\frac{13433}{10^9}
>\frac{4597}{2500}.
\end{aligned}
\tag{18}
\]
This proves (12). QED.

All these bounds are written series or rational inequalities.
They are not numerical source evaluations.

## PRW-T4. The full even source is positive through \(L=1/5\)

For every \(0<L\le1/5\),
\[
\boxed{
q_L[g]\ge\frac{22673}{6250000}\|g\|^2
\ge\frac1{300}\|g\|^2
\qquad(g\in D_L).
}
\tag{19}
\]
No neutral moment condition is imposed.

**The endpoint source bound.** Put \(L_*=1/5\) and
\(\alpha=1/100\). No prime is active because
\(2L_*=2/5<\log2\).
The pole-estimate loss in (3) satisfies
\[
e_{L_*,\alpha}
=\frac{99}{500000}\cosh^2(1/10)
<\frac{11979}{50000000}<\frac1{4000};
\tag{20}
\]
here \(\cosh(1/10)<11/10\).
Using (12) in (3),
\[
\begin{aligned}
\gamma_{L_*,\alpha}
>{}&-\frac{26871}{10000}+\frac{4597}{2500}
+\frac{6931}{20000}+\frac2{25}+\frac34-\frac1{4000}\\
={}&\frac{41}{125}=:\gamma .
\end{aligned}
\tag{21}
\]
The rank-one coefficient is exactly
\[
\beta_{L_*,\alpha}
=\frac34+\frac1{25}-\frac{99}{250}
=\frac{197}{500}=:\beta .
\tag{22}
\]
Since \(-\log(1-s^2)\ge s^2\), (4) now gives
\[
q_{L_*}[g]\ge
\int_I\left(\gamma+\frac{x^2}{4L_*^2}\right)|g(x)|^2\,dx
-\beta|a_0(g)|^2.
\tag{23}
\]

**The scalar rank-one integral.** For
\(v(x)=\gamma+x^2/(4L_*^2)\), its normalized inverse integral is
\[
J_\gamma
=\frac1{2L_*}\int_I\frac{dx}{v(x)}
=\int_0^1\frac{ds}{\gamma+s^2/4}
=\frac2{\sqrt\gamma}
\arctan\!\left(\frac1{2\sqrt\gamma}\right).
\tag{24}
\]
We claim
\[
J_\gamma<\frac{251}{100}.
\tag{25}
\]
Indeed \(\sqrt{41/125}>229/400\), so
\[
J_\gamma<
\frac{800}{229}\arctan(200/229).
\]
Put \(t=29/429\). Then
\(\arctan(200/229)=\pi/4-\arctan t\).
The retained bound \(\pi<22/7\) and the alternating series give
\[
\arctan t>t-\frac{t^3}{3}>\frac{337}{5000}.
\]
It follows that
\[
J_\gamma<
\frac{800}{229}
\left(\frac{11}{14}-\frac{337}{5000}\right)
=\frac{100564}{40075}<\frac{251}{100}.
\]
This proves (25).

**The complete source margin.** Weighted Cauchy--Schwarz in
(23), or the prime-free case of PRW-T2, gives
\[
\begin{aligned}
q_{L_*}[g]
&\ge(1-\beta J_\gamma)\int_Iv(x)|g(x)|^2\,dx\\
&\ge
\left(1-\frac{197}{500}\frac{251}{100}\right)
\frac{41}{125}\|g\|^2
=\frac{22673}{6250000}\|g\|^2.
\end{aligned}
\tag{26}
\]
The last rational coefficient exceeds \(1/300\).

For smaller L, every smooth even test in (-L,L) is also a
smooth even test in \((-L_*,L_*)\). Its pole square, Gamma energy
and norm are unchanged by zero extension. Extra disjoint prime
translations have zero correlation. Thus the endpoint inequality
applies to every smaller window. Form-core closure proves
(19) on all of D_L. QED.

The constant channel is retained in (23); its negative rank-one
coefficient is controlled by (24)--(25). No step discards that
channel, assumes the source is already positive or changes the
test family to the two-moment restricted one.

## PRW-T5. The unrestricted compact gate and finite-certificate consequence

On the full even source domain, retain PTE's positive energy
operator \(G_L\) and compact positive operator
\[
q_L[g]=\|G_L^{1/2}g\|^2-c_L\|g\|^2,\qquad
K_L=c_LG_L^{-1}.
\tag{27}
\]
For \(0<L\le1/5\), no prime is active and
\[
0<c_L=-\Omega_0/2<\frac{27}{10}.
\]
Equation (19) therefore gives
\[
\boxed{
G_L\ge(c_L+1/300)I,\qquad
\|K_L\|\le\frac{810}{811}<1.
}
\tag{28}
\]
For PTE's diagnostic parameter
\(\tau(\lambda)=\lambda/(1+\lambda)\), including
\(\tau(\infty)=1\), the consequence is
\[
I-\tau(\lambda)K_L\ge\frac1{811}I
\quad(0\le\lambda\le\infty,\ 0<L\le1/5).
\tag{29}
\]
The corresponding negative count is zero in this range.

**Proof.** The closed-form lower bound in (19) gives the first
inequality in (28). The spectral theorem then gives
\[
\|K_L\|\le\frac{c_L}{c_L+1/300}
\le\frac{27/10}{27/10+1/300}=\frac{810}{811}.
\]
Positivity of K_L and \(0\le\tau\le1\) prove (29). QED.

GCC-T5 now has its strict-positive hypothesis verified for the
unshifted full source at every \(L\le1/5\). Hence finite Gamma-cut,
prime-record and mode-cut certificates exist throughout this
enlarged unrestricted range. No particular finite matrix size,
Gamma cutoff or numerical implementation is evaluated here.
The direct analytic proof is (20)--(26).

The norm bound in (28) concerns the full even operator K_L.
It is distinct from FCP's bound for the restricted operator
\(K_L^\circ\) on profiles with two neutral moments.

## Technical notes — relations, numerical coordinates and angular cut probabilities

The working perspective requested for this continuation is to
take objects, relations, composition and invariance as primary.
Numerical coordinates are readings of that structure in a chosen
representation. Mathematics can organize relations before choosing
those coordinates. Algebraic laws, normalization, source matching
and proof obligations remain explicit constraints on the relations.

An angular coordinate is a useful exact example of this viewpoint.
The intrinsic data are the cut and its relations; theta and p are
alternative numerical coordinates on those data. This is the
precise sense of "mathematics beyond numbers" used here.
It changes how a normalized cut is read while preserving its
probability and source identities.

### A normalized two-channel cut

In a declared positive inner-product representation, let P be
an orthogonal projection and \(u\ne0\). Define
\[
p_{\rm seen}=\frac{\|Pu\|^2}{\|u\|^2},\qquad
p_{\rm unseen}=\frac{\|(I-P)u\|^2}{\|u\|^2}.
\]
Orthogonality gives \(p_{\rm seen}+p_{\rm unseen}=1\).
The angle
\[
\theta=\arccos\frac{\|Pu\|}{\|u\|}
\quad(0\le\theta\le\pi/2)
\]
therefore gives
\[
\boxed{
p_{\rm seen}=\cos^2\theta,\qquad
p_{\rm unseen}=\sin^2\theta .
}
\tag{30}
\]
This is an exact parametrization of normalized cut weights.
At the endpoint angles, probabilities can be exactly zero or one.
A strictly positive residual requires a separate property of the
specified source, such as the finite-depth geometric tail below.
A probability distribution over the angles would require its
own specified measure.

For the constant-mode cut of this note and nonzero g,
\(P_0g=a_0(g)\phi_0\), so
\[
\frac{|a_0(g)|^2}{N_g}=\cos^2\theta_g,\qquad
\frac{\|g-P_0g\|^2}{N_g}=\sin^2\theta_g .
\]
The Legendre estimate (6) becomes
\[
\boxed{\frac{J_L[g]}{N_g}\ge\frac34\sin^2\theta_g.}
\tag{31}
\]
The exact pole uses cosh rather than the constant profile;
its difference is still controlled by (7).
Thus (31) connects the angular reading to the actual proof,
without changing the source definition.

### Phase and cross-terms survive the probability reading

In a complexified two-dimensional channel space write
\[
u=\begin{pmatrix}\cos\theta\\e^{i\varphi}\sin\theta\end{pmatrix},
\qquad
M=\begin{pmatrix}A&B\\\overline B&D\end{pmatrix},
\quad A,D\in\mathbb R .
\]
Then
\[
\boxed{
u^*Mu=
A\cos^2\theta+D\sin^2\theta+
\sin(2\theta)\operatorname{Re}(Be^{i\varphi}).
}
\tag{32}
\]
The squared sine and cosine give the channel probabilities.
The double-angle term retains the coupling and relative phase.
Probabilities alone do not recover that phase or, in a real
representation, the missing relative sign.

For every angle and phase to give a nonnegative source reading,
the matrix must satisfy
\[
A\ge0,\qquad D\ge0,\qquad AD\ge|B|^2.
\tag{33}
\]
These are the two-dimensional source-positivity conditions.
The angular parametrization records them; normalization of the
probabilities does not establish them.

### Prime-record angles and the angle of the proved margin

PCM's geometric prime records have weights
\(w_{p,k}=(1-p^{-1/2})p^{-k/2}\).
For the direct record cut retaining \(0,\ldots,K\),
\[
\epsilon_{p,K}=\sum_{k>K}w_{p,k}=p^{-(K+1)/2}.
\]
Define
\[
\theta_{p,K}
=\arcsin\!\left(p^{-(K+1)/4}\right).
\]
Then
\[
\boxed{
\epsilon_{p,K}=\sin^2\theta_{p,K},\qquad
1-\epsilon_{p,K}=\cos^2\theta_{p,K}.
}
\tag{34}
\]
This describes retained versus discarded records. It does not
replace PCM's separate coherent prime-resolvent projection.
Every finite K gives a positive angle, while
\(\theta_{p,K}\to0\); no positive depth-independent angle floor
is imposed.

The rank-one coefficient in the endpoint proof has its own
angular reading. Having proved
\[
0<\rho:=\beta J_\gamma<\frac{49447}{50000}<1,
\]
define the certificate angle
\(\vartheta=\arccos\sqrt{\rho}\).
Equation (26) can then be read as
\[
\boxed{
q_{L_*}[g]\ge
\gamma\sin^2\vartheta\,N_g,\qquad
\sin^2\vartheta>\frac{553}{50000}.
}
\tag{35}
\]
This angle encodes the already established margin. It is a
deterministic certificate parameter; a physical probability law
is not being inferred from it. The real angle is available
because the actual source estimate has proved \(\rho<1\).

## Positive ranges and the remaining RH obligation

| Family | Half-width established | Available source margin |
|---|---|---|
| Full even family, earlier PTE range | \(L\le1/32\) | \(1/16\) |
| Full even family, new PRW range | \(L\le1/5\) | \(1/300\) |
| Constant- and cosh-neutral even family, FCP range | \(L\le111/200\) | \(1/6000\) |

TPW and BPC retain their stronger restricted margins on their
earlier ranges. The full even range and the larger restricted
range must not be identified.

The new endpoint is still prime-free. Extending the unrestricted
result into an interacting-prime window requires the complete
prime sum to satisfy a bound such as PRW-T2, or a different
actual-source certificate. The theorem does not claim that
(10) holds at every later width.

GCC's exact all-width, all-tolerance finite-certificate criterion
is unchanged. This note verifies a larger full-source range to
which its positive-certificate existence theorem applies.
The universal arithmetic sign remains unproved.

## Dependencies and verification scope

- [GCC-1](RH_GAMMA_CUT_CERTIFICATES.md): correct limiting
  source domain, spectral-bottom convergence and finite
  positive-certificate completeness.
- [FCP-1](RH_FINITE_CUT_POSITIVITY.md): sharper Gamma scalar,
  retained logarithm bounds and the restricted positive window.
- [TPW-1](RH_TWO_PRIME_WINDOW.md): polynomial-kernel identity
  with the nonzero constant-moment residual retained.
- [FPW-1](RH_FIRST_PRIME_WINDOW.md): kernel bound and
  Legendre energy identity.
- [BPC-1](RH_BOUNDARY_PRIME_CAPS.md): weighted translated-pair
  estimate for the complete active-prime sum.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): actual source,
  pole term, closed energy and unrestricted compact gate.
- [PCM-1](RH_PRIME_CUT_MEMORY.md): normalized prime records
  and the distinction between record and coherent-event cuts.

Weighted squares, rank-one domination, elementary series and
orthogonal projection probabilities are standard tools.
The work here is their explicit source-preserving use to prove
the stated unrestricted window and the accompanying angular
interpretation. No literature-wide novelty claim is made.

Only this note and its README entry are changed. The document
is stacked on draft PR #116 with a [skip ci] documentation commit.
Verification is written argument review, document structure,
repository links, immutable-head content and PR metadata.
No scientific computation, numerical source evaluation,
independent agent, replay or workflow dispatch is performed.
The PR remains an unmerged draft.
