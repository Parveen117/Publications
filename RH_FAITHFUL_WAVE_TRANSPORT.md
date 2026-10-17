# FWT-1 — Faithful continuous waves and the information left by a finite cut

17 October 2026. Written continuation of RCD-1 at
a744ea8bbb8be8955738c8b4733d584da21ec261.

A finite reading can match a source without recovering its full flow.
This note gives an exact continuous wave for the retained actual xi
source, derives it from carried native cut waves, and recovers the
source response by an exponential integral. No possible off-critical
pole is discarded.

CWR already proves the discrete wave-growth criterion. The continuous
criterion below is its source-bound counterpart, not a claimed new
general principle. The additional finite-cut result is an exact signed
wave error: its first exponentially weighted reading cancels, but the
next weighted moment is strictly positive on RCD's declared range.

The proposed phrase "irrational jump" is treated as an intuition about
information omitted by a cut. Arithmetic irrationality is not assumed
to cause that omission. A fully rational source/model example below
has the same nonzero cut error, with rational wave frequencies too.

No all-parameter wave bound, positive full-source probability model,
or RH proof is established. No scientific computation is used.

## Source and parameter conventions

Retain the actual source from CWR, CWF and SE39:
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q(z)=\frac{\mathcal X'(z)}{\mathcal X(z)}
=\sum_\mu\frac{\nu_\mu}{z+\mu}.
\tag{1}
\]
For each positive-height xi zero write
\[
\kappa_\mu=\gamma-i\delta,\qquad
\delta=\beta-\tfrac12,\qquad
\mu=\kappa_\mu^2,\qquad |\delta|<\tfrac12.
\]
The sum counts distinct source values with their positive integer
multiplicities \(\nu_\mu\). It is absolutely locally uniformly
convergent away from its poles. Conjugate source values and
multiplicities are paired. No positive actual source measure is
assumed.

For \(a\ge a_0=3249/4\), retain
\[
q(a)>0,\qquad
\frac1{q(a)}\sum_\mu\frac{\nu_\mu}{|a+\mu|}\le2,\qquad
\mathcal L_a[f]=\frac1{q(a)}
 \sum_\mu\frac{\nu_\mu}{a+\mu}
 f\!\left(\frac a{a+\mu}\right),
\tag{2}
\]
whenever the displayed series converges.

The real variable \(u\) below is an auxiliary continuous wave
coordinate, obtained from CWF's scaling of native wave steps. It is
not identified with spatial distance, physical time, the layer index,
LAT's logarithmic anchor coordinate, or LWP's bounded-generator
exponential parameter. No spatial distance metric is introduced.
The Laplace variable \(w\), source position \(z=w^2\), and anchor a
also keep their distinct roles.

## FWT-T1. A carried native wave has an exact continuous limit

Define
\[
\boxed{
\mathscr W_a(u)=\frac1{q(a)}
 \sum_\mu\frac{\nu_\mu}{a+\mu}\cos(\kappa_\mu u),
 \qquad u\in\mathbb R.
}
\tag{3}
\]
This is a continuous real even function, with
\[
\mathscr W_a(0)=1,\qquad
|\mathscr W_a(u)|\le2e^{D_\xi|u|},\qquad
D_\xi=\sup_\mu|\delta_\mu|\in[0,1/2].
\tag{4}
\]
The series converges absolutely, uniformly on bounded real
u-intervals.

It is the following limit of the existing native wave. For
\(B=b^2\ge a\), use
\[
y_B(\mu)=\frac B{B+\mu},\qquad
\tau_k(y)=T_k(2y-1),\qquad k_B(u)=\lfloor b|u|/2\rfloor.
\]
Then
\[
\boxed{
\mathscr W_a(u)=\lim_{B\to\infty}
\frac{q(B)}{q(a)}
\mathcal L_B\!\left[
 \frac{\tau_{k_B(u)}(y)}
      {1-(1-a/B)y}
\right].
}
\tag{5}
\]
The convergence is uniform on bounded real u-intervals. The
rational weight is precisely the squared constant-observer
transport factor of CWR, with reference anchor a.

**Proof.** For real u,
\[
|\cos((\gamma-i\delta)u)|
\le\cosh(\delta u)\le e^{D_\xi|u|}.
\]
The absolute source budget (2) gives uniform convergence and (4).
Conjugation symmetry makes the result real; cosine makes it even.

The expression before the limit in (5) equals exactly
\[
\frac1{q(a)}
\sum_\mu\frac{\nu_\mu}{a+\mu}\tau_{k_B(u)}(y_B(\mu)),
\]
because
\[
1-(1-a/B)y_B(\mu)=\frac{a+\mu}{B+\mu}.
\]
Thus the physical source weights stay at the fixed anchor a.

CWR gives
\[
\tau_k(y_B(\mu))
=\frac{\omega_B(\mu)^k+\omega_B(\mu)^{-k}}2,\qquad
\omega_B(\mu)=\frac{b+i\kappa_\mu}{b-i\kappa_\mu}.
\]
For each fixed source value,
\(\log\omega_B=2i\kappa_\mu/b+O_\mu(b^{-3})\).
Since \(2k_B(u)/b\to|u|\), this tends to
\(\cos(\kappa_\mu u)\), uniformly for u in a fixed bounded
interval. The floor contributes an error tending to zero there.

The bound
\[
|\tau_{k_B(u)}(y_B(\mu))|
\le\left(\frac{b+1/2}{b-1/2}\right)^{k_B(u)}
\le e^{|u|}
\]
holds for \(b\ge57/2\), using CWF's \(\log R_B\le2/b\).
Together with the fixed summable weights in (2), it makes the
source tail uniformly small. A finite-source truncation followed
by its uniform limit proves (5). QED.

This is a carried, growing-step limit. Dropping its rational weight
or taking a fixed wave step as B grows is a different limit and
can lose the signal, as CWR-T4 already demonstrates. Formula (5)
does not assert positive pairings for the transported rational
observer at every depth.

No termwise differentiation of (3) in u is assumed. The source
weight controls the continuous wave, not an unrestricted Taylor
series at \(u=0\).

## FWT-T2. Exponential readings recover the complete response and retain every pole

For \(\Re w>1/2\), define the ordinary one-sided Laplace reading
\[
\mathscr F_a(w)=\int_0^\infty e^{-wu}\mathscr W_a(u)\,du.
\]
Then
\[
\boxed{
\mathscr F_a(w)=
\frac{w}{w^2-a}\left(1-\frac{q(w^2)}{q(a)}\right).
}
\tag{6}
\]
At \(w=\sqrt a\), the apparent singularity is removable, with
value \(-\sqrt a\,q'(a)/q(a)\). Equivalently,
\[
\boxed{
q(w^2)=q(a)\left[
 1-\frac{w^2-a}{w}\mathscr F_a(w)\right].
}
\tag{7}
\]
Thus the complete wave together with its retained scalar \(q(a)\)
determines the actual meromorphic response by analytic continuation.
This is exact source recovery from the full function, not from
finitely many samples or intensities. It is a uniqueness statement;
no stable inversion bound for incomplete or noisy wave data is claimed.

**Proof.** Equation (4) bounds the absolute integrand by an
integrable exponential on \(\Re w>1/2\). The absolute source
budget also justifies interchanging the sum and integral.
For each actual source mode,
\[
\int_0^\infty e^{-wu}\cos(\kappa_\mu u)\,du
=\frac{w}{w^2+\mu}.
\]
Use
\[
\frac1{(a+\mu)(w^2+\mu)}
=\frac1{w^2-a}\left(\frac1{a+\mu}-\frac1{w^2+\mu}\right)
\]
to obtain (6). The derivative limit at \(w=\sqrt a\) and
rearrangement give the remaining claims. QED.

The right side of (6) gives a meromorphic continuation. At
either \(w_0=\pm i\kappa_\mu\), it has the exact residue
\[
\boxed{
\operatorname*{Res}_{w=w_0}\mathscr F_a(w)
=\frac{\nu_\mu}{2q(a)(a+\mu)}\ne0.
}
\tag{8}
\]
Indeed q has residue \(\nu_\mu\) at \(-\mu\), and the change
\(z=w^2\) contributes \(1/(2w_0)\). Its coefficient in (6)
is \(-w_0/[q(a)(w_0^2-a)]\), proving (8).
The actual zeros have \(\gamma>0\), so \(w_0\ne0\), and
\(a+\mu\ne0\). Distinct source values cannot give the same
\(w_0\); multiplicities add with the same coefficient.

In particular a hypothetical zero with \(\delta>0\) gives
a pole at \(w_0=\delta+i\gamma\) in the right half-plane.
The reflected partner supplies this case when the original
zero has \(\delta<0\). No such zero is asserted to exist.

The recovery identity holds whether or not RH is true. Thus exact
recoverability alone is not equivalent to RH. The additional global
boundedness or positivity requirement is made explicit in FWT-T3.

Standard Laplace-transform conventions and analyticity are
background tools; see [NIST DLMF 1.14(iii)](https://dlmf.nist.gov/1.14#iii).
The source identity and noncancellation above are derived explicitly.

## FWT-T3. Full-wave growth is the retained horizontal source defect

For every fixed \(a\ge a_0\), with \(\log0=-\infty\),
\[
\boxed{
\limsup_{u\to+\infty}
\frac{\log|\mathscr W_a(u)|}{u}=D_\xi.
}
\tag{9}
\]
Consequently the following are equivalent:
\[
\boxed{
\begin{aligned}
&\mathrm{RH};\\
&\mathscr W_a\text{ is bounded on }[0,\infty);\\
&\text{for every }\epsilon>0,\ 
 |\mathscr W_a(u)|\le C_\epsilon e^{\epsilon u}\ (u\ge0);\\
&|\mathscr W_a(u)|\le1\quad\text{for every }u\in\mathbb R;\\
&(\mathscr W_a(u_i-u_j))_{i,j=1}^n\succeq0
 \quad\text{for all finite real observation sets}.
\end{aligned}
}
\tag{10}
\]
Each assertion concerns the actual wave (3), with no fitted
source replacement.

**Proof of the growth law.** The upper bound in (4) gives the
upper inequality in (9). If the limsup were less than \(D_\xi\),
choose a real c strictly between them. Continuity and the
definition of limsup give
\(|\mathscr W_a(u)|\le C e^{cu}\) for every \(u\ge0\).

If \(c\ge0\), the definition of \(D_\xi>c\) and reflection
provide an actual zero with positive displacement \(\delta>c\).
If \(c<0\), any actual zero and its reflected partner provide
one with \(\delta\ge0>c\). The exponentially bounded wave has
a holomorphic Laplace transform throughout \(\Re w>c\).
Formula (6), initially valid on \(\Re w>1/2\), continues
meromorphically to this connected half-plane. Its nonzero
residue (8) at \(w=\delta+i\gamma\) contradicts holomorphy.
This proves (9), including the case \(D_\xi=0\).

Under RH, \(\mu=\gamma^2>0\) and the weights
\[
p_\mu(a)=\frac{\nu_\mu}{q(a)(a+\gamma^2)}
\]
are positive and sum to one. Thus (3) is the cosine transform
of a probability measure, obtained by splitting each \(p_\mu\)
equally between frequencies \(+\gamma\) and \(-\gamma\).
It follows directly that \(|\mathscr W_a|\le1\). Moreover
the displayed finite kernel form is a positive sum of
\(\left|\sum_j c_j e^{i\gamma u_j}\right|^2\) and its
negative-frequency counterpart, proving positivity.

Conversely, positivity of the two-position kernel at \(0,u\)
gives \(1-\mathscr W_a(u)^2\ge0\). This implies boundedness,
then the subexponential condition, and (9) gives \(D_\xi=0\).
The retained definition of \(D_\xi\) identifies this with RH.
QED.

If RH fails, some finite real u therefore satisfies
\[
\det\begin{pmatrix}1&\mathscr W_a(u)\\
                   \mathscr W_a(u)&1\end{pmatrix}<0.
\tag{11}
\]
No such u has been exhibited or evaluated. The statement concerns
arbitrary wave separation, not a two-number certificate obtainable
from the already covered finite range.

This is the continuous version of CWR-T3's discrete criterion,
with the same invariant \(D_\xi\). It does not prove the required
global bound. In particular, a local stationary value or stable
finite model cannot substitute for the universal assertions (10).

## FWT-T4. Anchor transport is an exact exponential convolution

For another anchor \(c\ge a_0\), define
\[
G_c(u)=\frac{e^{-\sqrt c\,|u|}}{2\sqrt c},\qquad
(G_c*f)(u)=\int_{\mathbb R}G_c(v)f(u-v)\,dv.
\]
There is no Fourier normalization factor in this convolution.
Then
\[
\boxed{
q(c)\mathscr W_c
=q(a)\left[\mathscr W_a+(a-c)G_c*\mathscr W_a\right].
}
\tag{12}
\]
The source response is unchanged; the observation weights have
been transported exactly.

**Proof.** Since \(\sqrt c>1/2\), (4) makes the convolution
absolutely convergent, including after the source sum is inserted.
The two half-line exponential integrals give
\[
G_c*\cos(\kappa_\mu\,\cdot)
=\frac{\cos(\kappa_\mu\,\cdot)}{c+\mu}.
\]
Apply this identity to (3), and use
\[
\frac1{a+\mu}\left(1+\frac{a-c}{c+\mu}\right)
=\frac1{c+\mu}.
\]
This proves (12). QED.

Denote the bracketed operator by
\(\mathsf T_{c\leftarrow a}=I+(a-c)G_c*\).
On the retained source-wave class its mode multiplier is
\((a+\mu)/(c+\mu)\). Thus
\[
\boxed{
\mathsf T_{d\leftarrow c}\mathsf T_{c\leftarrow a}
=\mathsf T_{d\leftarrow a},\qquad
\mathsf T_{a\leftarrow c}\mathsf T_{c\leftarrow a}=I.
}
\tag{13}
\]
The factors \(q(a)/q(c)\) telescope under normalized transport.
These composition statements are on the stated source class;
no assertion about arbitrary growing functions is needed.

Because \(\|G_c\|_1=1/c\), a bounded wave at a gives
\[
\|\mathscr W_c\|_\infty
\le\frac{q(a)}{q(c)}
 \left(1+\frac{|a-c|}{c}\right)\|\mathscr W_a\|_\infty.
\tag{14}
\]
Reversing the anchors proves that boundedness is independent
of the chosen anchor. This agrees with the exact growth
invariant (9). It does not establish boundedness at the first anchor.

The exponential kernel in (12) is a derived source reweighting.
Its appearance does not give a universal least probability
or identify the auxiliary coordinate with physical distance.

## FWT-T5. Every positive finite model has a stable, observable wave

Let \(A=\{t_1,\ldots,t_m\}\) have both actual finite pencils
\(L,J\succ0\), and retain SE48's
\[
q_A(z)=\sum_{j=1}^m\frac{r_j}{z+\eta_j},
\qquad r_j,\eta_j>0,
\tag{15}
\]
with distinct \(\eta_j\). For any positive anchor a set
\[
\pi_j^A(a)=\frac{r_j}{(a+\eta_j)q_A(a)},\qquad
\boxed{\mathscr W_{A,a}(u)=
 \sum_{j=1}^m\pi_j^A(a)\cos(\sqrt{\eta_j}\,u).}
\tag{16}
\]
The weights are positive and sum to one. This finite wave is
positive definite and bounded by one for every real u.

For an explicit realization in the positive spectral chart of
the acquired pairing, let
\[
\Gamma_A=\operatorname{diag}
 (\sqrt{\eta_1},\ldots,\sqrt{\eta_m},
  -\sqrt{\eta_1},\ldots,-\sqrt{\eta_m}),
\]
and let g have the corresponding entries
\(\sqrt{\pi_j^A(a)/2}\) on both halves. Then
\[
U_A(u)=e^{iu\Gamma_A},\qquad
U_A(u)^*U_A(u)=I,\qquad
\mathscr W_{A,a}(u)=g^*U_A(u)g.
\tag{17}
\]
This chart represents the retained finite positive model. It
does not supply an infinite positive actual-source space.

Every state of this \(2m\)-dimensional realization is observable
from its complex amplitude trajectory \(g^*U_A(u)x\) on any
open interval. Equivalently, for every \(T>0\),
\[
\boxed{
\mathsf O_{A,T}
=\int_0^T U_A(u)^*gg^*U_A(u)\,du\succ0.
}
\tag{18}
\]
Hence
\[
\int_0^T|g^*U_A(u)x|^2du
\ge\lambda_{\min}(\mathsf O_{A,T})\|x\|^2.
\tag{19}
\]

**Proof.** The finite positive spectral expansion proves (16)--(17).
A zero quadratic reading in (18) forces the continuous scalar
trajectory to vanish on the interval. Its distinct real
frequencies are \(\pm\sqrt{\eta_j}\), and every component of
g is nonzero. Linear independence of finitely many exponentials
therefore forces every component of x to vanish. One can see this
by differentiating at an interior point: the first \(2m\) equations
form an invertible frequency Vandermonde matrix. This proves (18);
the finite-dimensional minimum-eigenvalue inequality gives (19).
QED.

No lower bound uniform in m, T, or the model is claimed. Exact
amplitude observations are required: intensity alone can lose
phase. For example \(|e^{iu}|^2=|e^{2iu}|^2=1\) for every u,
although the trajectories have different frequencies.

Both finite stability and finite observability are thus available.
They establish completeness within that retained finite model,
not identification of its frequencies with every actual xi zero.
This distinction prevents the model's automatic stability from
being used as the missing premise of FWT-T3.

## FWT-T6. The RCD defect is a signed wave error with an exact positive moment

Now take \(a\in A\) with \(a\ge a_0\), so that SE48 gives
\(q_A(a)=q(a)\), \(q_A'(a)=q'(a)\). Define
\[
\mathscr E_{A,a}(u)=\mathscr W_a(u)-\mathscr W_{A,a}(u),
\qquad R_A(z)=q(z)-q_A(z).
\tag{20}
\]
The actual and finite versions of (6) give
\[
\boxed{
\int_0^\infty e^{-wu}\mathscr E_{A,a}(u)\,du
=-\frac{wR_A(w^2)}{q(a)(w^2-a)},\qquad \Re w>1/2.
}
\tag{21}
\]
Thus the actual residual has not been erased by passing to waves.

At the matched anchor, the exact first two damped readings are
\[
\boxed{
\int_0^\infty e^{-\sqrt a\,u}\mathscr E_{A,a}(u)\,du=0,
\qquad
\int_0^\infty u e^{-\sqrt a\,u}\mathscr E_{A,a}(u)\,du
=\frac{a}{q(a)}R_A''(a).
}
\tag{22}
\]
These are exponentially weighted moments of a continuous wave
error, not derivatives of the unsmoothed wave at \(u=0\).

**Proof.** Subtract (6) and its finite version to obtain (21).
Analyticity of q near a and the value/slope matching give
\[
R_A(w^2)=\tfrac12 R_A''(a)(w^2-a)^2
        +O((w^2-a)^3).
\]
The right side of (21) therefore vanishes at \(w=\sqrt a\)
and has derivative \(-aR_A''(a)/q(a)\) there.
Differentiation under the damped integral is justified by
\[
|\mathscr E_{A,a}(u)|\le2e^{u/2}+1.
\]
It yields (22). QED.

To obtain an actual strict sign and explicit bounds, retain
RCD's UCB extra-position range:
\[
A\subset[a/2,3a/2],\quad a\in A,\quad
N_*(a)=\left\lfloor\frac{e^{-4}\sqrt a}{4\pi}\right\rfloor
\ge128,\quad m+1\le\lfloor N_*(a)/12\rfloor.
\tag{23}
\]
Put
\[
b=\sqrt a,\quad \ell=\log(b/(2\pi)),\quad
\delta_{m+1}(a)=\frac6{\ell(441/4)^{m+1}},\quad
K_A=q(a)+\frac{41}{8a}+\frac1{8b},
\]
and denote the second integral in (22) by \(M_1(A,a)\).
Then RCD (17), (19), with fitted point and reference both a, give
\[
\boxed{
0<
\frac{4\delta_{m+1}(a)}{a}
 \prod_{\substack{t_k\in A\\t_k\ne a}}(t_k/a-1)^2
\le M_1(A,a)
\le
\frac{K_A b}{2q(a)(b-1)^2}
 \prod_{\substack{t_k\in A\\t_k\ne a}}
 \left|\frac{b-\sqrt{t_k}}{b+\sqrt{t_k}-2}\right|^2.
}
\tag{24}
\]
Empty products equal one.

Because the zeroth integral in (22) is zero but \(M_1>0\),
the continuous error takes both positive and negative values
on \(u>0\). If it had only one sign, its positive damping weight
would force it to vanish identically, contradicting \(M_1>0\).
Thus the finite-cut defect is a signed wave discrepancy.

This proves that the omitted information persists in the full
flow and becomes visible in a next weighted reading. It does
not turn \(\mathscr E_{A,a}\) into a probability density.
Its bounds depend on anchor, bank size and positions; no
positive universal error floor or least probability is obtained.

## FWT-T7. Nonzero cut error does not diagnose arithmetic irrationality

Consider the exact rational toy source
\[
q_{\mathrm{toy}}(z)=\frac{28}{z+1}+\frac{125}{z+4},
\qquad A=\{1\}.
\tag{25}
\]
This is not the xi source and is not a new source certificate.
At the sole observation, its finite pencils and probe are
\(L=12\), \(J=27\), and \(u=39\). Therefore
\[
q_{\mathrm{toy},A}(z)=\frac{u^2}{J+zL}
=\frac{507}{4z+9}.
\]
Direct polynomial subtraction gives
\[
\boxed{
q_{\mathrm{toy}}(z)-q_{\mathrm{toy},A}(z)
=\frac{105(z-1)^2}{(z+1)(z+4)(4z+9)}.
}
\tag{26}
\]
The source and model match value and slope at 1, but their
second-derivative difference is \(21/13>0\).

All source poles, residues, observation data, and the model pole
and residue are rational. Even the wave frequencies are rational:
the source frequencies are 1 and 2; the model frequency is \(3/2\).
The corresponding normalized waves are
\[
\boxed{
\mathscr W_{\mathrm{toy},1}(u)
=\frac{14}{39}\cos u+\frac{25}{39}\cos(2u),\qquad
\mathscr W_{\mathrm{toy},A,1}(u)=\cos(3u/2).
}
\tag{27}
\]
Both are stable positive finite-model correlations. Their
difference is nevertheless nonzero, with
\[
\int_0^\infty e^{-u}
 (\mathscr W_{\mathrm{toy},1}-\mathscr W_{\mathrm{toy},A,1})du=0,
\]
\[
\int_0^\infty u e^{-u}
 (\mathscr W_{\mathrm{toy},1}-\mathscr W_{\mathrm{toy},A,1})du
=\frac7{169}>0.
\tag{28}
\]
These follow from (22) with \(a=1\), \(q_{\mathrm{toy}}(1)=39\);
the finite rational source justifies that transform directly.
The xi-specific anchor bound is not invoked for this toy example.

The example is written exact algebra, not a numerical experiment.
It shows that a reduction from two source modes to one model mode
can create the observed defect without any irrational source or
model frequency. Such an error is therefore not, by itself,
evidence of an irrational-coordinate jump.

Conversely, \(q_{\mathrm{one}}(z)=1/(z+2)\) has the irrational
wave frequency \(\sqrt2\), yet its singleton pencil at any positive
observation t reconstructs it exactly: \(L=1/(t+2)^2\),
\(J=2/(t+2)^2\), and \(u=1/(t+2)\) give
\(u^2/(J+zL)=1/(z+2)\). Its residual is identically zero.
Within this response/cut construction, irrational frequency is
therefore neither necessary nor sufficient for a nonzero residual.

An irrational number can also have an exact finite specification,
for example the positive solution of \(x^2=2\). A nonterminating
decimal expansion does not alone impose a nonzero error floor.
Any proposed cut-resolution floor would need an additional
specified observation rule and its proof.

## What the flow formulation establishes

| Question | Result |
|---|---|
| Does the continuous wave come from native cut waves? | Yes, via the carried growing-step limit (5) |
| Is the actual source retained? | The Laplace recovery (7) and nonzero residues (8) retain it |
| What measures possible off-line instability? | The exact growth exponent is the same \(D_\xi\) as CWR |
| Does changing anchor erase it? | No; the exponential convolution composes and is invertible on the source class |
| Are the finite models stable and observable? | Yes, with respect to their own finite positive spectral realizations |
| Does finite-model completeness establish actual-source completeness? | No; the source/model difference remains (21) |
| Where is RCD's missing next layer in the flow? | In the strictly positive damped moment (22)--(24) |
| Is that error automatically a probability density? | No; the wave discrepancy has both signs |
| Does it require irrational coordinates or frequencies? | No; (25)--(28) give an entirely rational witness |
| Is RH now proved? | No; the all-parameter actual-wave bound remains unproved |

Under RH the normalized source wave is a probability characteristic
function on real frequencies. That conclusion is part of the
equivalence, not a premise available merely because a wave was
constructed. The actual complex source weights must remain in (3)
until the needed global property is proved.

The exponential in the integral and anchor kernel has a derived
transport/recovery role. The \(2\pi\) period of a real-frequency
phase is a rotation convention. Neither establishes a minimum
probability, a universal hidden-information quantum, or an
identification of scalar response curvature with native commutator
curvature or physical curvature.

The exact remaining target is to prove a global subexponential
bound for the actual \(\mathscr W_a\), equivalently any condition
in (10), from the acquired source relations. CWR already records
the discrete version of this obligation. This note provides the
continuous source bridge and exposes the finite-cut error within
it; it does not remove the all-parameter requirement.

UCB's bank range, MOP's dyadic range, BRR's all-width requirement,
and the established full/restricted support windows are unchanged.

## References, attribution and verification scope

- [RCD-1](RH_RESIDUAL_COMPLETION.md): complete residual factorization and strict second-derivative bounds.
- [CWR-1](RH_CUT_WAVE_RETURN.md), T1--T4: actual source budget, carried observer, native waves, uncancelled discrete poles and \(D_\xi\).
- [CWF-1](RH_CRITICAL_WAVE_FORMS.md), T1: continuous critical scaling of the native wave and the distinction between bilinear source forms and positive measures.
- [UCB-1](RH_OBSERVER_BANK_UNION.md), T6: actual finite response sign and complete-source Schwarz bound.
- [SE48](papers/native-prime-observation-closure/sources/SE48.md), theorem 1: positive finite realization with value and slope interpolation.
- [LWP-1](papers/native-layer-wave-propagation/THEOREM.md): the separate finite native wave generator, exact horizons and conditional infinite extension.
- [LAT-1](papers/native-layer-anchor-transport/THEOREM.md): exact native anchor transport and its normalization.
- [MOP-1](RH_MULTI_OBSERVER_PROBABILITY.md): retained model probabilities.
- [BRR-1](RH_BOUNDARY_OBSERVER_RETURN.md): the remaining unrestricted source-window condition.
- [NIST DLMF 1.14(iii)](https://dlmf.nist.gov/1.14#iii): standard Laplace-transform background.

The Laplace/cosine transform, spectral realization, finite observability,
and analytic pole argument are standard tools. No literature-wide
novelty or machine-checked analytic proof is claimed. The retained
discrete criterion is explicitly attributed to CWR rather than
presented as a newly solved global sign problem.

Only this note and its README entry change. All derivations are written;
no numerical xi value, derivative, spectrum, wave, or probability is
evaluated. The rational witness is checked by displayed exact algebra.
No scientific computation, independent agent, source replay, or workflow
dispatch is performed. Frozen sources, certificates, engines and
workflows remain unchanged. The commit uses [skip ci]. The PR is
stacked on draft #124 and remains unmerged.
