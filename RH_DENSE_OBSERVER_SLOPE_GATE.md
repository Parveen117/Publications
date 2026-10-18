# DSG-1 — An exact residual-slope gate for dense observer extension

18 October 2026 (Asia/Kolkata). Author: Monty.
Written continuation of NOE-1 at
25e8418fcf2be2bebf430112305dc2341d53a89f.

## Result and scope

NOE proves that a positive finite actual-source bank can always be
extended sufficiently far away. Its displayed sufficient test also
requires the new position to be at least twice the previous maximum.
Repeated use of that restriction necessarily gives a summable
reciprocal-square-root sampling sequence. It cannot supply NOE's
sufficient source-identification density.

This note replaces that particular test with an exact scalar gate,
derived from the unchanged actual Schur innovations and RCD's complete
residual. The known interpolation zeros are divided out analytically.
The resulting gate is regular when an observer approaches an old
position and imposes no minimum spacing ratio.

On UCB's existing finite range, both normalized gate margins are
strictly positive, including their confluent limits. This proves that
the complete filtered residual decreases while its product with the
source position increases. It also supplies a two-sided enclosure of
the third derivative of the residual at each fitted point. No enlargement of UCB's total bank-size range is claimed.

The unrestricted dense positive extension remains open. The exact
gate and the coverage estimate below identify its two requirements;
they do not establish them at all orders. RH is not proved.

## 1. Retained source, bank and residual

Use the actual source
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q=\mathcal X'/\mathcal X,
\]
and a bank \(A=\{t_1,\ldots,t_m\}\), \(t_i>9/4\), with both actual
pencils strictly positive:
\[
L_{ij}=-q[t_i,t_j],\qquad J_{ij}=(tq(t))[t_i,t_j],
\qquad u_i=q(t_i).
\tag{1}
\]
Diagonal entries use derivatives. SE48 gives
\[
q_A(z)=u^T(J+zL)^{-1}u=\sum_{j=1}^m\frac{r_j}{z+\eta_j},
\qquad r_j,\eta_j>0.
\tag{2}
\]
Write
\[
R_A(z)=q(z)-q_A(z),\quad
\omega_A(z)=\prod_i(z-t_i),\quad
d_A(z)=\prod_j(z+\eta_j),\quad
H_A(z)=\omega_A(z)/d_A(z).
\]
RCD-T1 proves the complete-source identity
\[
\boxed{
R_A(z)=H_A(z)^2\,\mathcal T_A(z),\qquad
\mathcal T_A(z)=\sum_\mu
 \frac{\nu_\mu h_A(\mu)^2}{z+\mu},
\quad
h_A(\mu)=\frac{\prod_j(\mu-\eta_j)}
                   {\prod_i(\mu+t_i)}.
}
\tag{3}
\]
The actual \(\mu\) may be nonreal. The square \(h_A(\mu)^2\)
remains bilinear; it is not replaced by an absolute square.
The function \(\mathcal T_A\) is analytic at the positive fitted
positions. At those points its value is given by this retained
analytic continuation, not by assigning a value to \(0/0\).

Throughout, a prime denotes differentiation in the positive source
position \(x\), with the physical bank \(A\) fixed. It is not a
derivative of the wave in its auxiliary coordinate, a cut count, or
physical time. The symbol \(H_A\) here is the scalar rational
interpolation factor, not a thermodynamic Hessian.

## DSG-T1. The old far-step restriction cannot supply the sufficient density

Suppose successive new positions obey the side condition in NOE-T2,
\[
t_{n+1}\ge2\max_{j\le n}t_j.
\]
Then
\[
t_n\ge2^{n-1}t_1,\qquad
\boxed{\sum_{n\ge1}t_n^{-1/2}
 \le\frac{t_1^{-1/2}}{1-2^{-1/2}}<\infty.}
\tag{4}
\]
The same conclusion holds after any fixed finite initial bank.

**Proof.** Induction gives the first inequality; summing the geometric
majorant gives the second. QED.

Thus a sequence constructed using this side condition at every step
cannot satisfy NOE-T7's divergent-sum hypothesis. This is a restriction
of that sufficient test, not a proof that the actual admissible
positions must be so far apart. An exact test without the side
condition is derived next.

## DSG-T2. The exact Schur signs factor into smooth residual expressions

For real \(x>9/4\) outside \(A\), set
\[
d(x)=((x-t_i)^{-1})_{i=1}^m,\qquad
b_A(x)=H_A(x)d(x),\qquad P(x)=J+xL.
\]
The vector \(b_A\) extends smoothly through every old position:
\[
(b_A(x))_i=\frac{\prod_{k\ne i}(x-t_k)}{d_A(x)}.
\]
Define
\[
K_L(x)=b_A(x)^TL^{-1}b_A(x),\qquad
K_J(x)=b_A(x)^TJ^{-1}b_A(x).
\tag{5}
\]
Both are strictly positive for positive \(x\), including at the old
positions, because \(b_A(x)\ne0\).

Let \(\alpha_A(x),\beta_A(x)\) be the actual SE49 Schur innovations
when appending the actual value and slope at \(x\). Then
\[
\boxed{
\alpha_A(x)=H_A(x)^2\,\mathfrak a_A(x),\qquad
\mathfrak a_A=-\mathcal T_A'-K_L\mathcal T_A^2,
}
\tag{6}
\]
\[
\boxed{
\beta_A(x)=H_A(x)^2\,\mathfrak b_A(x),\qquad
\mathfrak b_A=\mathcal T_A+x\mathcal T_A'
                         -x^2K_J\mathcal T_A^2.
}
\tag{7}
\]
Consequently \(A\cup\{x\}\) is positive exactly when
\(\mathfrak a_A(x)>0\) and \(\mathfrak b_A(x)>0\).
No restriction such as \(x\ge2\max A\) occurs.

**Proof.** Let \(r_A=q_A(x)\) and \(v=P(x)^{-1}u\).
The finite model's mixed columns are
\[
(l_A)_i=\frac{u_i-r_A}{x-t_i},\qquad
(j_A)_i=\frac{x r_A-t_i u_i}{x-t_i}.
\]
Its spectral representation gives
\[
l_A=Lv,\quad j_A=Jv,\quad
l_A^TL^{-1}l_A=-q_A'(x),\quad
j_A^TJ^{-1}j_A=q_A(x)+xq_A'(x).
\tag{8}
\]
For example, put \(\Phi_{ji}=\sqrt{r_j}/(t_i+\eta_j)\).
Then \(L=\Phi^T\Phi\), \(J=\Phi^T\operatorname{diag}(\eta_j)\Phi\),
and \(\Phi\) is invertible. These identities give (8) directly.

With \(R=R_A(x)\), the actual mixed columns are
\[
l=l_A-Rd,\qquad j=j_A+xRd.
\]
Subtracting the model's identities (8) from the actual Schur formulas
therefore yields
\[
\alpha_A=-R'+2R\,d^Tv-R^2d^TL^{-1}d,
\]
\[
\beta_A=R+xR'-2xR\,d^Tv-x^2R^2d^TJ^{-1}d.
\tag{9}
\]

The cross term has the exact value
\[
\boxed{d(x)^Tv=\frac{H_A'(x)}{H_A(x)}.}
\tag{10}
\]
To verify it, use RCD's rational residual-observer identity
\[
1-(x+\mu)f(\mu)^Tv=H_A(x)h_A(\mu),
\qquad f_i(\mu)=(t_i+\mu)^{-1}.
\]
Differentiate in the auxiliary \(\mu\) at \(\mu=-x\).
The left side becomes \(d(x)^Tv\). On the right,
\(H_A(x)h_A(-x)=1\), and
\[
h_A'(-x)/h_A(-x)
=\sum_i\frac1{x-t_i}-\sum_j\frac1{x+\eta_j}
=H_A'(x)/H_A(x).
\]
This proves (10).

Finally insert \(R=H_A^2\mathcal T_A\) into (9). The
\(2H_AH_A'\mathcal T_A\) terms cancel exactly, while
\(H_A^2d^TL^{-1}d=K_L\) and similarly for \(J\).
This proves (6)--(7). QED.

These are the original source signs, with their known vanishing
interpolation factor removed. Dividing by \(H_A^2\) does not turn a
negative sign positive. Regularity of the normalized expressions
does not assert their positivity outside a proved range.

## DSG-T3. Positivity is an exact inverse-residual slope window

At a new position \(x\notin A\), both signs in (6)--(7) hold if and
only if \(\mathcal T_A(x)>0\) and, on its positive neighborhood,
the inverse residual
\[
V_A(x)=1/\mathcal T_A(x)
\]
obeys
\[
\boxed{
K_L(x)<V_A'(x)<\frac{V_A(x)}x-xK_J(x).
}
\tag{11}
\]
Equivalently, with \(g_A=xV_A'/V_A=-x\mathcal T_A'/\mathcal T_A\),
\[
\boxed{
xK_L\mathcal T_A<g_A<
1-x^2K_J\mathcal T_A.
}
\tag{12}
\]
In particular the admissible window requires
\[
V_A>xK_L+x^2K_J,\qquad 0<g_A<1.
\tag{13}
\]

**Proof.** From (6)--(7),
\[
x\mathfrak a_A+\mathfrak b_A
=\mathcal T_A-(xK_L+x^2K_J)\mathcal T_A^2.
\]
If both normalized signs are positive, this expression is positive.
It is nonpositive when \(\mathcal T_A\le0\), so
\(\mathcal T_A>0\) is necessary. Now divide each normalized
innovation by \(\mathcal T_A^2\):
\[
\mathfrak a_A/\mathcal T_A^2=V_A'-K_L,\qquad
\mathfrak b_A/\mathcal T_A^2=V_A-xV_A'-x^2K_J.
\]
This proves both directions of (11). Multiplication by
\(x/V_A=x\mathcal T_A\) gives (12), and (13) follows. QED.

For a fixed bank, on any interval where (11) holds, integration gives
\[
V_A(s)-V_A(r)>\int_r^s K_L(x)\,dx,
\]
\[
\boxed{
\frac{V_A(r)}r-\frac{V_A(s)}s>\int_r^s K_J(x)\,dx
\quad(r<s).
}
\tag{14}
\]
The second inequality is the integral of
\((V_A/x)'<-K_J\). The bank must remain fixed during these
integrations; changing the bank changes all these functions.

The gate is unchanged by any invertible change of coordinates on
the old finite bank. Indeed
\[
L\mapsto S^TLS,\quad J\mapsto S^TJS,\quad
u\mapsto S^Tu,\quad b_A\mapsto S^Tb_A
\]
leaves \(q_A,H_A,\mathcal T_A,K_L,K_J\) and the window unchanged.
The physical source positions are held fixed in this assertion.
No invariance under an arbitrary change of the source variable or
of the physical bank is asserted.

This removes dependence on the raw inverse-norm bounds used by
NOE's far-step sufficient gate. It does not remove the retained
source derivative that actually decides the sign.

## DSG-T4. The existing finite range proves both normalized margins and a third-derivative enclosure

Use UCB's comparison anchor \(a\), with
\[
I_a=[a/2,3a/2],\quad
N_*(a)=\left\lfloor\frac{e^{-4}\sqrt a}{4\pi}\right\rfloor\ge128,
\quad
M(a)=\lfloor N_*(a)/12\rfloor.
\]
Let \(A\subset I_a\) have \(m\ge1\) distinct positions and
\(m+1\le M(a)\). Put
\[
\ell_a=\log(\sqrt a/(2\pi)),\quad
\delta_{m+1}(a)=\frac6{\ell_a(441/4)^{m+1}},\quad
S_{m+1}(z)=\sum_{k=0}^m z^{2k}.
\]
For every \(x\in I_a\), INCLUDING the old positions, define
\[
\varepsilon_L(x)=
\frac{q(a)\delta_{m+1}(a)d_A(x)^2}
     {a^{2m+1}S_{m+1}(x/a-1)},\qquad
\varepsilon_J(x)=a\varepsilon_L(x).
\tag{15}
\]
Then
\[
\boxed{
\mathfrak a_A(x)\ge\varepsilon_L(x)>0,\qquad
\mathfrak b_A(x)\ge\varepsilon_J(x)>0.
}
\tag{16}
\]
Consequently the complete filtered source satisfies throughout \(I_a\)
\[
\boxed{
\mathcal T_A>0,\qquad
-\mathcal T_A'\ge K_L\mathcal T_A^2+\varepsilon_L>0,
}
\]
\[
\boxed{
\mathcal T_A+x\mathcal T_A'
\ge x^2K_J\mathcal T_A^2+\varepsilon_J>0.
}
\tag{17}
\]
Thus \(\mathcal T_A(x)\) strictly decreases, while
\(x\mathcal T_A(x)\) strictly increases, on this finite interval.
Both statements concern the full filtered actual source.

**Proof of the margins.** For \(x\notin A\), UCB-T3 gives
\[
\alpha_A(x)\ge\frac{q(a)}a\delta_{m+1}(a)
 \frac{\prod_i((x-t_i)/a)^2}{S_{m+1}(x/a-1)},
\]
with the same bound multiplied by \(a\) for \(\beta_A(x)\).
Divide by
\[
H_A(x)^2=\frac{a^{2m}\prod_i((x-t_i)/a)^2}{d_A(x)^2}.
\]
Equations (6)--(7) give exactly (15)--(16).
All normalized expressions are continuous at an old position, so
the bounds extend there. Positivity of \(\mathcal T_A\) follows
from the sum identity in T3, including at those limiting points.
Equations (6)--(7) then prove (17). QED.

For the inverse residual, the explicit strengthened window is
\[
\boxed{
K_L+\frac{\varepsilon_L}{\mathcal T_A^2}
\le V_A'
\le\frac{V_A}x-xK_J
       -\frac{\varepsilon_J}{x\mathcal T_A^2}.
}
\tag{18}
\]
No probability law for the entire actual xi source was inserted
in obtaining these inequalities.

**Confluent observation and the third source jet.** At \(t_i\in A\),
write
\[
h_i=H_A'(t_i)\ne0,\quad k_i=H_A''(t_i),\quad
T_i=\mathcal T_A(t_i).
\]
The regularized vector and coefficients have the exact values
\[
b_A(t_i)=h_i e_i,\quad
K_L(t_i)=h_i^2(L^{-1})_{ii},\quad
K_J(t_i)=h_i^2(J^{-1})_{ii}.
\tag{19}
\]
As a distinct new position approaches \(t_i\),
\[
\boxed{
\frac{\alpha_A(x)}{(x-t_i)^2}\to h_i^2\mathfrak a_A(t_i)>0,\qquad
\frac{\beta_A(x)}{(x-t_i)^2}\to h_i^2\mathfrak b_A(t_i)>0.
}
\tag{20}
\]
The raw innovations vanish at coincidence, as they must for duplicate
rows. The strict limits concern the divided-observer direction;
they do not make duplicate rows independent.

Taylor expansion of the retained identity \(R_A=H_A^2\mathcal T_A\)
gives
\[
T_i=\frac{R_A''(t_i)}{2h_i^2}>0,\qquad
R_A'''(t_i)=6h_i^2\mathcal T_A'(t_i)+6h_i k_i T_i.
\tag{21}
\]
Combining (17) with (21) supplies the explicit enclosure
\[
\boxed{
\begin{aligned}
6h_i k_iT_i+
6h_i^2\left[-\frac{T_i}{t_i}
+t_iK_J(t_i)T_i^2+\frac{\varepsilon_J(t_i)}{t_i}\right]
&\le R_A'''(t_i),\\
R_A'''(t_i)
&\le
6h_i k_iT_i-
6h_i^2\left[K_L(t_i)T_i^2+\varepsilon_L(t_i)\right].
\end{aligned}
}
\tag{22}
\]
This bounds the next derivative in terms of the retained second
residual derivative and finite model data. It does not assert that
\(R_A'''\) has a universal sign.

Equations (19)--(22) follow from the simple zeros of \(\omega_A\),
the smooth normalized gate, and the Taylor identity. All derivatives
are source-position derivatives of analytic functions. No termwise
differentiation of FWT's continuous wave is used.

RCD had proved \(R_A''(t_i)>0\) and \(\mathcal T_A>0\) on this
range. The extra conclusions here are the normalized two-channel
margins, the filtered derivative inequalities, and (22).
The size budget remains \(m+1\le M(a)\). A fixed finite budget
cannot be iterated indefinitely merely by making spacings smaller.

## DSG-T5. The exact gate and a retained coverage sum separate acceptance from recovery

For a positive bank containing a fixed reference \(a\ge3249/4\),
put \(b=\sqrt a\), \(s_i=\sqrt{t_i}\), and
\[
C_a=q(a)+\frac{41}{8a}+\frac1{8b}.
\]
At a real test \(w>1\) distinct from the sampled \(s_i\), define
\[
\mathcal I_A(w)=
2\sum_{t_i\in A}\log\left|\frac{w+s_i-2}{w-s_i}\right|.
\tag{23}
\]
Each summand is positive. The complete-source Schwarz estimate
retained in NOE-T7 gives
\[
\boxed{
|q(w^2)-q_A(w^2)|
\le C_a\frac{(w+b)^2}{w}\,e^{-\mathcal I_A(w)}.
}
\tag{24}
\]
At a sampled \(w=s_i\), interpolation already gives zero error;
one may read \(\mathcal I_A(w)=+\infty\).

When one new point \(x=s^2\) is accepted by the exact gate,
its contribution is
\[
\boxed{
\Delta\mathcal I(w;x)=
2\log\left|\frac{w+s-2}{w-s}\right|
=\frac{4(w-1)}s+O_w(s^{-2})
\quad(s\to\infty).
}
\tag{25}
\]
Thus, along positions tending to infinity and away from an exactly
sampled test, divergence of \(\mathcal I\) is equivalent to
divergence of the reciprocal-square-root sampling sum.
The exponential remainder bound (24) then tends to zero at every
fixed unsampled test.

**Proof.** Formula (24) is the retained product estimate with its
product written as an exponential. For real \(w>1,\ s>1\),
\((w+s-2)^2-(w-s)^2=4(w-1)(s-1)>0\), proving positivity.
Taking the logarithm of the new factor gives (25); expansion at
large \(s\) gives its asymptotic expression. The positive summands
are eventually comparable to \(1/s\), proving the divergence
equivalence. QED.

Accordingly a source-identifying construction can be specified by
two distinct requirements:

1. Each actual new position satisfies the scalar window (11),
   equivalently the two original source signs.
2. The accepted positions yield a divergent coverage sum, for example
   \(\sum_n t_n^{-1/2}=\infty\) with \(t_n\to\infty\).

With both requirements, (24) and PWC-T3 give actual-source recovery,
faithful-wave positivity and RH. Requirement 1 alone gives the
positive model limit of NOE, which need not recover the source.
Requirement 2 cannot be claimed for the geometrically separated
construction ruled out by T1.

A finite value of the right side of (24), or a nonvanishing Schwarz
product, is not a lower bound on actual error. These are sufficient
recovery estimates, not necessary conditions for every possible
source reconstruction.

## Remaining problem

The exact residual-slope test now permits investigation of close
positions and interlaced banks without NOE's imposed geometric
spacing. On UCB's established finite range its two margins are
strictly proved, including coincident limits and their next jet.
This is a concrete extension of the available source derivative
information.

What remains is an unlimited sequence whose actual normalized
margins stay positive while its coverage diverges. No bound here
controls those margins after UCB's total-size budget is exhausted.
No existence of that sequence, enlargement of the finite range,
or RH closure is claimed.

The inverse residual \(V_A\) is a defined scalar of the retained
source. Calling it a potential does not supply its derivative sign.
The inequalities above are the actual sign conditions and the
finite range where they have been proved.

## Sources, attribution, and verification

- [NOE-1](RH_NESTED_OBSERVER_EXTENSION.md), T2, T5--T7:
  far-step gate, positive model limit, non-recovery example,
  and retained density criterion.
- [RCD-1](RH_RESIDUAL_COMPLETION.md), T1, T3:
  complete residual factorization and positive second-derivative defect.
- [UCB-1](RH_OBSERVER_BANK_UNION.md), T2--T3, T6:
  finite actual source signs, spacing-normalized innovations,
  and the source remainder estimate.
- [SE48](papers/native-prime-observation-closure/sources/SE48.md),
  T1: finite positive spectral representation and interpolation.
- [SE49](papers/native-moving-zero-covariance/sources/SE49.md),
  T1--T3: exact source Schur tests, analytic remainder and density argument.
- [RLO-1](RH_RELATIONAL_OBSERVER_CONFLUENCE.md):
  divided observations and lawful confluent limits.
- [PWC-1](RH_PROBABILITY_WAVE_CLOSURE.md), T3:
  actual-source convergence implies faithful probability-wave closure.

Schur complements, rational interpolation and Taylor limits are
standard tools. The written additions are their exact residual-slope
form for this retained source, its finite actual margins, and the
third-derivative enclosure. No literature-wide novelty is claimed.

Validation consists of written algebra and analytic reasoning plus
text/metadata checks. No scientific computation, numerical source
evaluation, independent agent, certificate replay or workflow
dispatch is performed. Only this note and its README entry change.
The commit uses [skip ci]. Published as an unmerged draft stacked
on #132.
