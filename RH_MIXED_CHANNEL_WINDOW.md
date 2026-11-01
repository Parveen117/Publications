# MCW-1 — Mixed observer channels and a larger unrestricted positive window

17 October 2026. Written continuation of FWT-1 at
167eee544a2c5cfb5cef5b8032c82ac26a4b2e0d.

The next step is an actual source inequality, using the later
PRW/FCP estimates already available beyond CWF. Retaining and
bounding the constant/quadratic observer coupling gives
\[
\boxed{
0<L\le\frac14
\quad\Longrightarrow\quad
q_L[g]\ge\frac{667}{2000000}\|g\|_2^2
>\frac1{3000}\|g\|_2^2
\quad(0\ne g\in D_L).
}
\]
Here \(D_L\) is the full real even source form domain. No constant,
quadratic, or cosh moment is set to zero.

The unrestricted half-width increases from PRW's \(1/5\) to \(1/4\),
a factor \(5/4\). PRW's stronger margin on its earlier interval is
retained. The larger FCP interval \(L\le111/200\) still concerns
its two-moment restricted family.

The proof keeps the finite quadratic moment residual that vanished
in the restricted arguments. It supplies an explicit two-channel
comparison matrix with all active prime terms retained in a weighted
bound. At the new endpoint, the exact mixed-channel estimate and a
rational lower bound for the logarithmic boundary weight close the
source inequality.

The endpoint is still before the first prime threshold
\(L=(\log2)/2\). This is not an all-width positivity proof or an
RH closure. No scientific computation is used.

## 1. The complete source and the two observed channels

Use the actual source and closed form domain from CWF, PTE,
ASM, GCC, FCP and PRW. On \(I=(-L,L)\), functions are real
and even and are extended by zero. Write
\[
\begin{gathered}
N_g=\|g\|_2^2,\qquad
C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,\\
\ell(g)=\int_I g(x)\cosh(x/2)\,dx,\qquad
m_j(g)=\int_I x^j g(x)\,dx,\\
\mathcal N_L=
 \{n\ge2:n\text{ is a prime power},\ \log n<2L\},\qquad
a_n=\Lambda(n)/\sqrt n .
\end{gathered}
\]
The retained complete arithmetic form is
\[
\boxed{
q_L[g]=\ell(g)^2+\frac{\Omega_0}{2}N_g
+\int_0^\infty k(u)(N_g-C_g(u))\,du
-\sum_{n\in\mathcal N_L}a_n C_g(\log n),
}
\tag{1}
\]
where
\[
\begin{gathered}
\Omega_0=-\gamma_E-\pi/2-3\log2-\log\pi,\qquad
k(u)=\frac{e^{u/2}}{2\sinh u},\\
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}).
\end{gathered}
\tag{2}
\]
The domain \(D_L\) is the retained closed Gamma form domain.
Smooth real even compact tests form its defining core.
All initial inequalities below are proved on that core and
extended by closure when their uniform bound is established.

Define
\[
\phi_0(x)=\frac1{\sqrt{2L}},\qquad
\phi_2(x)=\sqrt{\frac5{2L}}P_2(x/L),\qquad
P_2(s)=\frac{3s^2-1}{2},
\]
\[
c_0=\langle g,\phi_0\rangle,\qquad
c_2=\langle g,\phi_2\rangle,\qquad
\mathbf c=(c_0,c_2)^{\mathsf T},\qquad
w_L(x)=-\frac14\log(1-(x/L)^2).
\tag{3}
\]
The indices 0 and 2 denote polynomial degrees. Neither c_0 nor
c_2 is a source anchor, prime coefficient, or assumed zero moment.

Retain the following already proved estimates:
\[
k(u)\ge\frac1{2u}+\frac6{25}-\frac{u^2}{25}
\quad(0<u\le9/8),
\tag{4}
\]
\[
\frac{\Omega_0}{2}>-\frac{26871}{10000},\qquad
\log2>\frac{6931}{10000},\qquad
\frac{\pi}{4}>\frac{7853}{10000}.
\tag{5}
\]
These are FCP-T3 and PRW-T3. Their written proofs use elementary
series and rational inequalities, not xi evaluations.

## MCW-T1. The quadratic Gamma improvement retains an explicit mixed channel

For \(0<L\le9/16\), put
\[
B_L=
\frac{\Omega_0}{2}+F(2L)+\frac12\log2
+\frac{12L}{25}-\frac{8L^3}{75}.
\tag{6}
\]
For every admitted smooth even g,
\[
\boxed{
q_L[g]\ge
B_LN_g+J_L[g]+\int_I w_L(x)|g(x)|^2dx
+\ell(g)^2-\frac3{25}m_0(g)^2
+\frac1{25}m_0(g)m_2(g)
-\sum_{n\in\mathcal N_L}a_n C_g(\log n),
}
\tag{7}
\]
where
\[
J_L[g]=\frac18\int_I\int_I
 \frac{|g(x)-g(y)|^2}{|x-y|}\,dy\,dx .
\]

**Proof.** Apply TPW-T1's general polynomial-kernel identity to
\(P(u)=6/25-u^2/25\). It applies throughout the stated range
because \(2L\le9/8\). Evenness gives \(m_1=0\), and the exact
correlation residual is
\[
\int_0^{2L}P(u)C_g(u)\,du
=\frac3{25}m_0^2-\frac1{25}m_0m_2.
\]
Also
\(\int_0^{2L}P(u)\,du=12L/25-8L^3/75\).
Insert these in the complete source (1). This proves (7). QED.

The two ordinary moments have the exact coordinates
\[
m_0=\sqrt{2L}\,c_0,\qquad
m_2=L^2\sqrt{2L}
 \left(\frac{c_0}{3}+\frac{2c_2}{3\sqrt5}\right).
\]
Consequently the retained mixed contribution is
\[
\boxed{
\frac1{25}m_0m_2
=\frac{L^3}{75}\,
\mathbf c^{\mathsf T}
\begin{pmatrix}2&2/\sqrt5\\2/\sqrt5&0\end{pmatrix}
\mathbf c.
}
\tag{8}
\]
Its two matrix eigenvalues, including the scale, are
\[
\frac{L^3}{75}(1+3/\sqrt5),\qquad
\frac{L^3}{75}(1-3/\sqrt5).
\]
Since \(c_0^2+c_2^2\le N_g\), this gives the uniform bound
\[
\boxed{
\frac1{25}m_0m_2\ge-\chi_LN_g,\qquad
\chi_L=\frac{L^3}{75}(3/\sqrt5-1)<\frac{L^3}{150}.
}
\tag{9}
\]
The last inequality uses \(\sqrt5>2\).

Thus the polynomial-kernel correction has a controlled negative
direction in its two-channel comparison. This is not a negative
direction of the complete arithmetic source. The pole, Gamma,
boundary, and prime terms remain in (7). No moment cancellation
has been imported from the restricted family.

## MCW-T2. An explicit two-channel gate retains the pole and prime couplings

For \(0<L\le9/16\) and \(0<\alpha<1\), retain PRW's pole-error estimate
\[
\ell(g)^2\ge
2L(1-\alpha)c_0^2-e_{L,\alpha}N_g,\qquad
e_{L,\alpha}=(\alpha^{-1}-1)
 \frac{L^5}{160}\cosh^2(L/2).
\tag{10}
\]
Its Taylor-bound proof is valid on this entire interval.
The even Legendre energy also gives
\[
J_L[g]\ge
\frac{25}{24}N_g-\frac{25}{24}c_0^2-\frac7{24}c_2^2.
\tag{11}
\]
This follows from the eigenvalues \(0,3/4,25/24,\ldots\) of
the even logarithmic energy: the constant and degree-2
coefficients are retained explicitly.

Define
\[
\Gamma_{L,\alpha}=B_L+\frac{25}{24}-e_{L,\alpha},
\qquad
\beta_{L,\alpha}=\frac{25}{24}+\frac{6L}{25}
                 -2L(1-\alpha),
\]
\[
\boxed{
T_{L,\alpha}=
\begin{pmatrix}
\beta_{L,\alpha}-2L^3/75&-2L^3/(75\sqrt5)\\
-2L^3/(75\sqrt5)&7/24
\end{pmatrix}.
}
\tag{12}
\]
Equations (7)--(11) give
\[
\boxed{
q_L[g]\ge
\int_I(\Gamma_{L,\alpha}+w_L)|g|^2dx
-\mathbf c^{\mathsf T}T_{L,\alpha}\mathbf c
-\sum_{n\in\mathcal N_L}a_nC_g(\log n).
}
\tag{13}
\]
This is a source-derived comparison; it is not an equality
with the complete source.

To retain the active prime terms, choose a positive even measurable
weight v such that
\[
0<v(x)\le\Gamma_{L,\alpha}+w_L(x).
\]
For \(u_n=\log n\), set
\[
R_v(x)=\frac12\sum_{n\in\mathcal N_L}a_n
\left(
\frac{\mathbf1_I(x-u_n)}{\sqrt{v(x)v(x-u_n)}}+
\frac{\mathbf1_I(x+u_n)}{\sqrt{v(x)v(x+u_n)}}
\right),
\]
\[
d(x)=(1-R_v(x))v(x),\qquad
m=\operatorname*{ess\,inf}_{I}d(x).
\tag{14}
\]
Only evaluate a shifted term when both of its arguments lie
inside I. Suppose \(m>0\), and the indicated weighted integrals
are finite on the source core. Define the two-by-two matrix
\[
\mathsf J_d=
\left(\int_I\frac{\phi_i(x)\phi_j(x)}{d(x)}dx\right)_{i,j\in\{0,2\}}.
\tag{15}
\]
It is positive definite. For \(0<\eta\le1\), the sufficient gate is
\[
\boxed{
\mathsf J_d^{-1}-T_{L,\alpha}
\succeq\eta\,\mathsf J_d^{-1}
\quad\Longrightarrow\quad
q_L[g]\ge\eta m\|g\|_2^2.
}
\tag{16}
\]

**Proof.** The retained weighted translation inequality gives
\[
\sum_n a_n|C_g(u_n)|
\le\int_I R_v(x)v(x)|g(x)|^2dx.
\]
It follows from (13) that
\[
q_L[g]\ge
\int_I d(x)|g(x)|^2dx-\mathbf c^{\mathsf T}T_{L,\alpha}\mathbf c.
\tag{17}
\]
This step bounds the actual signed correlations without assuming
that inserting each new prime is a positive update.

The weight d is even. The minimum of \(\int_I d|h|^2\)
over the real even weighted L2 space with
\(\langle h,\phi_i\rangle=c_i\), \(i=0,2\), is
\(\mathbf c^{\mathsf T}\mathsf J_d^{-1}\mathbf c\).
Indeed its minimizer is
\[
h_*(x)=d(x)^{-1}(\phi_0(x),\phi_2(x))
       \mathsf J_d^{-1}\mathbf c,
\]
and weighted orthogonality gives
\[
\int_I d|h|^2
=\mathbf c^{\mathsf T}\mathsf J_d^{-1}\mathbf c
+\int_I d|h-h_*|^2.
\tag{18}
\]
This is a variational identity in the comparison space; no
source-domain membership of h_* is required for the lower
bound applied to source profiles.

Under (16),
\[
\mathbf c^{\mathsf T}T_{L,\alpha}\mathbf c
\le(1-\eta)\mathbf c^{\mathsf T}\mathsf J_d^{-1}\mathbf c
\le(1-\eta)\int_I d|g|^2.
\]
Use (17) and \(d\ge m\) to obtain the asserted margin.
Its extension to D_L follows by source form closure. QED.

The non-strict matrix condition
\(\mathsf J_d^{-1}-T_{L,\alpha}\succeq0\) is necessary and
sufficient for nonnegativity of the real even weighted comparison
form in (17). It is only a sufficient condition for the actual source,
because (13) and the prime bound are lower estimates.
Neither this finite matrix nor the required positive weight is
declared to have the needed sign at every L.

There is a useful scalar consequence. Dropping only the
nonnegative higher Legendre excess gives
\(J_L[g]\ge(3/4)(N_g-c_0^2)\).
Combining this with (9) and (10), put
\[
\gamma^\flat_{L,\alpha}
=B_L+\frac34-e_{L,\alpha}-\chi_L,\qquad
\beta^\flat_{L,\alpha}
=\frac34+\frac{6L}{25}-2L(1-\alpha).
\]
Then
\[
\boxed{
q_L[g]\ge
\int_I(\gamma^\flat_{L,\alpha}+w_L)|g|^2dx
-\beta^\flat_{L,\alpha}c_0^2
-\sum_{n\in\mathcal N_L}a_n C_g(\log n).
}
\tag{19}
\]
This bound will suffice for the explicit unrestricted extension.

## MCW-T3. A rational boundary estimate at the new endpoint

The required new boundary constant is
\[
\boxed{F(1/2)>\frac{5111}{3000}.}
\tag{20}
\]

**Proof.** Put \(z=\tanh(1/8)\). The same exact identity used
in PRW gives
\[
F(1/2)=-\frac12\log z+\frac{\pi}{4}-\arctan z.
\tag{21}
\]
PRW's proved polynomial upper bound for tanh is
\[
\tanh x<x-\frac{x^3}{3}+\frac{2x^5}{15}\qquad(x>0).
\]
At \(x=1/8\), it yields
\[
z<\frac{30561}{245760}
=\frac18\left(1-\frac{53}{10240}\right)
<\frac{311}{2500}.
\tag{22}
\]
The elementary series bounds
\(\cosh(1/8)\le128/127\) and \(\sinh(1/8)>1/8\) also give
\[
z>\frac{127}{1024}>\frac{31}{250},\qquad z<\frac18.
\tag{23}
\]

Using (5), (22), and \(-\log(1-t)>t\),
\[
-\frac12\log z
>\frac32\log2+\frac{53}{20480}
>\frac{20793}{20000}+\frac{129}{50000}.
\]
The alternating arctangent upper bound and (23) give
\[
\arctan z\le z-z^3/3+z^5/5
<\frac{311}{2500}-\frac{63}{100000}+\frac1{160000}.
\]
Here \((31/250)^3/3>63/100000\), and
\(z^5/5<1/163840<1/160000\).
Insert these estimates and \(\pi/4>7853/10000\) in (21):
\[
F(1/2)>
\frac{1363003}{800000}>\frac{5111}{3000}.
\]
The last comparison is an integer cross-multiplication.
This proves (20). QED.

These constants are bounded by written series and fractions;
no numerical Gamma or xi evaluation is performed.

## MCW-T4. A rational boundary weight gives a controlled inverse integral

For \(0\le z<1\),
\[
-\log(1-z)\ge\frac{z}{1-z/2}.
\tag{24}
\]
Indeed the two sides agree at zero, and the derivative of their
difference is
\[
\frac1{1-z}-\frac1{(1-z/2)^2}
=\frac{z^2}{4(1-z)(1-z/2)^2}\ge0.
\]

At \(\gamma=23/100\), define
\[
v_L(x)=\gamma+
\frac{(x/L)^2}{4-2(x/L)^2}.
\]
Then \(v_L\le\gamma+w_L\), \(v_L\ge\gamma\), and the normalized
inverse integral is independent of L:
\[
J_\gamma=\frac1{2L}\int_I\frac{dx}{v_L(x)}
=\int_0^1
 \frac{4-2s^2}{4\gamma+(1-2\gamma)s^2}\,ds.
\]
It obeys the explicit rational bound
\[
\boxed{J_\gamma<\frac{317}{100}.}
\tag{25}
\]

**Proof of the constant.** Set \(r=27/46\), and write
\[
J_\gamma=-\frac{100}{27}+\frac{5000}{621}A,\qquad
A=\int_0^1\frac{ds}{1+rs^2}.
\]
A finite geometric expansion gives the strict upper bound
\[
A<
1-\frac r3+\frac{r^2}{5}-\frac{r^3}{7}
+\frac{r^4}{9}-\frac{r^5}{11}
+\frac{r^6}{13}-\frac{r^7}{15}+\frac{r^8}{17}.
\tag{26}
\]
The omitted remainder has negative sign: pointwise the displayed
even truncation exceeds \(1/(1+rs^2)\).

For completeness, the grouped fraction comparisons are
\[
\begin{aligned}
1-r/3=37/46&<804348/10^6,\\
r^2/5-r^3/7&<40016/10^6,\\
r^4/9-r^5/11&<6855/10^6,\\
r^6/13-r^7/15&<1546/10^6,\\
r^8/17&<829/10^6.
\end{aligned}
\]
Each follows by substituting \(r=27/46\) and multiplying positive
integer denominators. Their sum gives
\[
A<\frac{853594}{10^6}
<\frac{426857}{500000}.
\]
Finally
\[
-\frac{100}{27}
+\frac{5000}{621}\frac{426857}{500000}
=\frac{317}{100}.
\]
This proves (25). QED.

The rational weight keeps more of the actual logarithmic boundary
energy than the quadratic lower weight used in PRW. No new
source positivity assumption enters this inverse-integral bound.

## MCW-T5. Unrestricted actual positivity extends to \(L=1/4\)

For every \(0<L\le1/4\), the full even source satisfies
\[
\boxed{
q_L[g]\ge\frac{667}{2000000}\|g\|_2^2
\ge\frac1{3000}\|g\|_2^2,\qquad g\in D_L.
}
\tag{27}
\]
The second comparison is strict for nonzero g.

**Endpoint estimates.** Set \(L_*=1/4\), \(\alpha=1/100\).
There is no active prime term, because
\(2L_*=1/2<\log2\). The polynomial kernel estimate is valid
since \(2L_*<9/8\).

The pole loss is bounded by
\[
e_{L_*,\alpha}
=\frac{99}{163840}\cosh^2(1/8)
\le\frac{99}{161290}<\frac1{1600},
\]
using \(\cosh(1/8)\le128/127\).
Equation (9) gives
\(\chi_{L_*}<L_*^3/150=1/9600\).

Insert (5), (20), and these losses in the scalar coefficient
of (19):
\[
\begin{aligned}
\gamma^\flat_{L_*,\alpha}
>{}&-\frac{26871}{10000}+\frac{5111}{3000}
+\frac{6931}{20000}+\frac3{25}-\frac1{600}
+\frac34-\frac1{1600}-\frac1{9600}\\
={}&\frac{55373}{240000}>\frac{23}{100}=\gamma .
\end{aligned}
\tag{28}
\]
The remaining rank-one coefficient is exactly
\[
\beta^\flat_{L_*,\alpha}
=\frac34+\frac3{50}-\frac{99}{200}
=\frac{63}{200}=\beta .
\tag{29}
\]
No observed coefficient has been set to zero.

By (19), (24) and (28),
\[
q_{L_*}[g]\ge
\int_I v_{L_*}(x)|g(x)|^2dx-\beta c_0^2.
\]
Weighted Cauchy--Schwarz, with \(\phi_0=1/\sqrt{2L_*}\), gives
\[
c_0^2\le J_\gamma\int_I v_{L_*}(x)|g(x)|^2dx.
\]
Use (25), \(v_{L_*}\ge\gamma\), and
\(\beta J_\gamma<1\):
\[
\begin{aligned}
q_{L_*}[g]
&\ge(1-\beta J_\gamma)\int_I v_{L_*}|g|^2dx\\
&\ge
\left(1-\frac{63}{200}\frac{317}{100}\right)
\frac{23}{100}\,N_g
=\frac{667}{2000000}N_g.
\end{aligned}
\tag{30}
\]
Since \(667\cdot3000>2000000\), this exceeds \(N_g/3000\)
for nonzero g.

Every smooth even test supported in a smaller interval is also
a test at the endpoint. Zero extension leaves its actual source
form and norm unchanged; correlations at disjoint shifts vanish.
This transfers the bound to all smaller L. The retained defining
core and form closure extend it to D_L. QED.

This is a full-family source inequality. It does not use FCP's
assumptions \(m_0=0\) and \(\ell(g)=0\). Those restricted
results remain separate and retain their larger support interval.

## MCW-T6. The full compact gate and certificate range extend with the source bound

Retain PTE's full even positive energy representation
\[
q_L[g]=\|G_L^{1/2}g\|^2-c_L\|g\|^2,\qquad
K_L=c_LG_L^{-1}.
\]
For \(0<L\le1/4\), no prime is active, so
\(c_L=-\Omega_0/2\) and \(0<c_L<27/10\).
Equation (27) gives
\[
\boxed{
G_L\ge(c_L+1/3000)I,\qquad
\|K_L\|\le\frac{8100}{8101}<1.
}
\tag{31}
\]
For the retained diagnostic coupling
\(\tau(\lambda)=\lambda/(1+\lambda)\), with
\(\tau(\infty)=1\), it follows that
\[
\boxed{
I-\tau(\lambda)K_L\ge\frac1{8101}I
\quad(0\le\lambda\le\infty,\ 0<L\le1/4).
}
\tag{32}
\]

**Proof.** The closed-form source margin implies the lower
energy bound. The spectral theorem and \(c_L<27/10\) give
\[
\|K_L\|\le\frac{c_L}{c_L+1/3000}
\le\frac{27/10}{27/10+1/3000}=\frac{8100}{8101}.
\]
Positivity of K_L and \(0\le\tau\le1\) prove (32). QED.

GCC's finite-certificate completeness theorem now applies to
the strictly positive full source throughout this larger window.
It supplies existence of finite Gamma-cut and mode-cut certificates
with the required errors retained. No particular matrix size,
cutoff, or numerical certificate has been evaluated here.

The stronger PRW margin \(1/300\) and compact bound \(810/811\)
remain available on \(L\le1/5\). The new bounds cover a larger
interval; they do not replace the stronger earlier estimates
where those already apply.

## What advanced and what remains

| Item | Before this note | Established here |
|---|---|---|
| Full real even positive half-width | \(L\le1/5\), PRW | \(L\le1/4\) |
| Guaranteed margin across the enlarged full window | Not supplied by PRW | \(667/2000000>1/3000\) |
| Constant/quadratic kernel residual | Vanishes only on the declared mean-zero family, or must be retained | Exact two-channel matrix and explicit negative-part bound |
| Low-channel comparison with active primes | PRW's scalar gate | The two-channel gate (12)--(18), with all active prime correlations bounded explicitly |
| Full compact norm on the new window | Not supplied | At most \(8100/8101\) |
| Larger two-moment restricted window | \(L\le111/200\), FCP | Unchanged |
| Arbitrary-window positivity / RH | Open | Open |

The new endpoint is still below \((\log2)/2\). Thus this result
does not yet certify an unrestricted window with an active prime
translation. The prime-aware matrix in MCW-T2 is a sufficient
mechanism for testing such extensions, not their established sign.

The completed mixed-channel estimate is a concrete instance of
the Schur-sign work selected after FWT. It includes the ordinary
moment residual and the actual pole-profile error. Replacing the
full source by its positive Gamma terms, dropping the constant
channel, or imposing new neutral moments is not part of the proof.

CWF's all-width arithmetic inequality, BRR's unrestricted scaled
width condition, and FWT's full-wave subexponential bound remain
unproved. This finite support extension does not by itself give
any of them. RBP/UCB/MOP's native finite-bank size ranges are
unchanged. No global zero exclusion or actual off-line zero is
asserted.

## References and verification scope

- [PRW-1](RH_POLE_RECOVERY_WINDOW.md): retained full-source domain, pole-error estimate, prior unrestricted endpoint, and scalar prime gate.
- [FCP-1](RH_FINITE_CUT_POSITIVITY.md), T3: the wider quadratic Gamma bound and retained scalar constants; T6: the separate restricted range.
- [TPW-1](RH_TWO_PRIME_WINDOW.md), T1: exact polynomial-kernel moment residual and logarithmic boundary decomposition.
- [BPC-1](RH_BOUNDARY_PRIME_CAPS.md), T3: weighted simultaneous prime-translation domination.
- [ASM-1](RH_ARITHMETIC_SEAM_REDUCTION.md): domain-preserving source cuts, positive unresolved tails, and full Schur memory.
- [GCC-1](RH_GAMMA_CUT_CERTIFICATES.md): existence of finite positive certificates when a strict source margin is known.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): full positive energy and compact-gate representation.
- [CWF-1](RH_CRITICAL_WAVE_FORMS.md): the actual arithmetic form and all-width RH equivalence.
- [FWT-1](RH_FAITHFUL_WAVE_TRANSPORT.md): the faithful source wave and remaining global growth condition.

The two-channel variational calculation, Legendre energy, weighted
Cauchy--Schwarz and scalar series bounds use standard tools with the
retained actual source. No literature-wide novelty, optimal endpoint,
or machine-checked analytic proof is claimed.

Only this note and its README entry change. The new constants and
matrix identities are written derivations. No numerical source value,
derivative, wave, spectrum, probability, or finite certificate is
evaluated. No scientific computation, independent agent, source replay,
or workflow dispatch is performed. Frozen sources, certificates, engines
and workflows remain unchanged. The documentation commit uses [skip ci].
The PR is stacked on draft #125 and remains unmerged.
