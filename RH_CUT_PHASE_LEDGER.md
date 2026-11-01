# CPL-1 — Native cut-phase records and infinite exponential flow

17 October 2026. Written continuation of ZMR-1 at
e3e6a41a8a7a6d8f49e637d4515eb8f4cd352c12.

The question is whether the rotation associated with each cut can be
retained as memory, and what that record shows after indefinitely
many repetitions. This note records the actual native two-reflection
step, rather than replacing its angle by a different continuous
generator. The record retains both the lifted angle and any
logarithmic amplitude change.

The resulting exponential equals the native wave at every integer
step. In a positive rotation sector the phase keeps accumulating:
there is no endpoint limit, although the mean reading is zero.
Its frequency-resolved memory and its mean square remain nonzero.
For the actual xi source the exact exponential growth rate is the
already established CWR defect, now read through this interpolation.
Its vanishing remains unproved.

These are analytic records and written proofs, not newly acquired
numerical phase data or an evaluated infinite list of zeros.

## Retained source and the meaning of a cut-step

Use the conventions of CWR and FWT:
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q(z)=\frac{\mathcal X'(z)}{\mathcal X(z)}
=\sum_\mu\frac{\nu_\mu}{z+\mu}.
\]
For a positive-height zero \(\beta+i\gamma\), put
\[
\delta=\beta-\tfrac12,\qquad
\kappa_\mu=\gamma-i\delta,\qquad
\mu=\kappa_\mu^2,\qquad |\delta|<\tfrac12.
\]
Distinct source values are counted with their multiplicities
\(\nu_\mu\). Fix \(a\ge a_0=3249/4\), and set
\[
b=\sqrt a,\qquad y_\mu=\frac a{a+\mu},\qquad
w_\mu(a)=\frac{\nu_\mu}{q(a)(a+\mu)}.
\tag{1}
\]
The retained source gives
\[
q(a)>0,\qquad \sum_\mu w_\mu(a)=1,\qquad
\sum_\mu|w_\mu(a)|\le2.
\tag{2}
\]
Weights and source values can be complex. No full positive measure
is assumed.

The native cut and reference reflection are
\[
P(y)=\begin{pmatrix}y&1-y\\y&1-y\end{pmatrix},\qquad
Z=\operatorname{diag}(1,-1),\qquad
g(y)=\operatorname{diag}(y,1-y).
\]
The associated two-reflection wave is
\[
\boxed{
W(y)=(2P(y)-I)Z
=\begin{pmatrix}2y-1&-2(1-y)\\2y&2y-1\end{pmatrix}.
}
\tag{3}
\]
Here \(P^2=P\), while \(W^TgW=g\) and \(\det W=1\).
The noninvertible projection \(P\) is not itself a rotation.
Throughout, a recorded cut-step means the existing operation (3),
with its reference reflection retained.

Cut count \(N\), real interpolation parameter \(t\), anchor \(a\),
and finite carrier depth are separate quantities. The parameter
\(t\) is not identified with physical time or the lambda layer index.

## CPL-T1. The one-step record has a unique chosen lift

The CWR wave multiplier is
\[
\omega_a(\mu)=\frac{b+i\kappa_\mu}{b-i\kappa_\mu}.
\]
Choose its logarithm in the upper half-plane:
\[
\boxed{
\ell_a(\mu)=\log\omega_a(\mu)
=\alpha_a(\mu)+i\theta_a(\mu),\qquad 0<\theta_a(\mu)<\pi.
}
\tag{4}
\]
Then
\[
\boxed{
\begin{aligned}
\alpha_a(\mu)
&=\frac12\log
 \frac{(b+\delta)^2+\gamma^2}{(b-\delta)^2+\gamma^2}\\
&=\operatorname{artanh}
 \frac{2b\delta}{a+\gamma^2+\delta^2},\\
\theta_a(\mu)
&=\arctan\frac{\gamma}{b+\delta}
 +\arctan\frac{\gamma}{b-\delta}.
\end{aligned}
}
\tag{5}
\]
In particular,
\[
\alpha_a(\mu)=0\ \Longleftrightarrow\ \delta=0.
\tag{6}
\]

**Proof.** Write the numerator as \(b+\delta+i\gamma\) and
the denominator as \(b-\delta-i\gamma\). Both real parts are
positive. Their arguments give (5), with the indicated range.
Equivalently
\[
\Im\omega_a(\mu)=
\frac{2b\gamma}{(b-\delta)^2+\gamma^2}>0,
\]
so the branch in (4) is well-defined. Taking the modulus gives
the logarithmic ratio; writing its numerator and denominator
as \(a+\gamma^2+\delta^2\pm2b\delta\) gives artanh.
The argument has absolute value less than one because
\((b-|\delta|)^2+\gamma^2>0\). This proves (6). QED.

For the conjugate source value, the positive-height choice gives
\[
\omega_a(\bar\mu)=\frac1{\overline{\omega_a(\mu)}},
\qquad
\ell_a(\bar\mu)=-\overline{\ell_a(\mu)}.
\tag{7}
\]
Thus conjugation reverses \(\alpha\) while retaining \(\theta\).
Reversing the wave exchanges \(\ell\) and \(-\ell\).

The interval in (4) is a specified lift convention. A single
endpoint matrix alone cannot reveal extra complete turns during
a physical cut: other interpolations can add integer multiples
of \(2\pi\) to eigenphases. We use the displayed branch and carry
it thereafter. No unobserved physical winding is inferred.

## CPL-T2. An exponential that returns the actual cut at each step

Put \(x=2y-1\) and retain LWP's skew part
\[
A(y)=\frac{W-W^{-1}}2=W-xI
=\begin{pmatrix}0&-2(1-y)\\2y&0\end{pmatrix}.
\tag{8}
\]
At each source value, let \(\ell=\ell_a(\mu)\) and
\[
s_\mu=\sinh\ell
=\frac{2ib\kappa_\mu}{a+\mu}\ne0.
\]
Then
\[
x=\cosh\ell,\qquad A^2=s_\mu^2 I=-4y(1-y)I.
\]
Define, mode by mode,
\[
\boxed{
B_a(\mu)=\frac{\ell}{\sinh\ell}A(y_\mu),\qquad
V_a(t;\mu)=e^{tB_a(\mu)}
=\cosh(t\ell)I+
 \frac{\sinh(t\ell)}{\sinh\ell}A(y_\mu).
}
\tag{9}
\]
This gives
\[
V_a(0;\mu)=I,\quad V_a(1;\mu)=W(y_\mu),\quad
V_a(N;\mu)=W(y_\mu)^N,
\]
\[
V_a(t+u;\mu)=V_a(t;\mu)V_a(u;\mu).
\tag{10}
\]

**Proof.** CWR gives \(x=(\omega+\omega^{-1})/2\);
direct subtraction gives the displayed formula for \(s_\mu\).
The source excludes \(\kappa_\mu=0\) and \(a+\mu=0\).
Thus \(B^2=\ell^2I\), and the even and odd exponential series
give (9). At \(t=1\) this is \(xI+A=W\). The group law is the
exponential law for the same fixed matrix. QED.

This differs from the earlier LWP/CEM exponential \(e^{tA}\).
To see the distinction exactly, take a positive real sector
\(0<y<1\), put \(\theta=\arccos(2y-1)\), and let
\[
S=\operatorname{diag}(\sqrt y,\sqrt{1-y}),\qquad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]
Then
\[
SWS^{-1}=e^{\theta R},\qquad
SAS^{-1}=\sin\theta\,R,\qquad
SBS^{-1}=\theta R.
\tag{11}
\]
The old exponential advances by \(t\sin\theta\); (9) advances
by \(t\theta\). Both are well-defined, with different specified
generators. LWP-T4 explicitly did not assert \(e^A=W\).

This also exposes the scalar orientation fold: the old balanced
return \(\cos(t\sin\theta)\) is unchanged by
\(\theta\mapsto\pi-\theta\), corresponding to \(y\mapsto1-y\).
The new balanced return \(\cos(t\theta)\), as a full function
of \(t\), distinguishes those angles except at their common
midpoint. ZMR repaired the old channel using reference motion;
the present construction instead retains the actual step angle.

No positive infinite Hilbert-space realization is assumed in
(9). It is a two-by-two source-mode identity. The coefficient
\(\ell/\sinh\ell\) need not stay bounded in the unweighted
matrix chart as source values approach its endpoint. The
convergent scalar source reading is proved separately below.

## CPL-T3. Repeated-cut memory accumulates instead of resetting

For each fixed source mode and anchor, retain the record
\[
\mathcal R_N(a,\mu)=
\bigl(a,\mu,\nu_\mu,w_\mu(a);
       N,\mathcal A_N,\Theta_N\bigr),
\quad
\mathcal A_0=\Theta_0=0,
\]
with the exact update
\[
\boxed{
\mathcal A_{N+1}=\mathcal A_N+\alpha_a(\mu),\qquad
\Theta_{N+1}=\Theta_N+\theta_a(\mu).
}
\tag{12}
\]
Therefore
\[
\mathcal A_N=N\alpha,\quad \Theta_N=N\theta,\qquad
\omega^N=e^{\mathcal A_N+i\Theta_N}.
\tag{13}
\]
If an endpoint angle is also stored, use
\[
\Theta_N=2\pi k_N+\varphi_N,\qquad
k_N\in\mathbb Z,\quad 0\le\varphi_N<2\pi.
\tag{14}
\]
Keeping only \(\varphi_N\) discards the completed turns.

The infinite-flow rates are exact:
\[
\boxed{
\frac{\mathcal A_N+i\Theta_N}{N}=\ell,\qquad
\frac1N\log|\omega^N|=\alpha,\qquad
\frac1N\log\rho_{\rm sp}(W^N)=|\alpha|
\quad(N\ge1).
}
\tag{15}
\]
Here \(\rho_{\rm sp}\) is spectral radius, not an assumed positive
source norm. The last equality follows from the two eigenvalues
\(\omega^N,\omega^{-N}\).

If \(\alpha=0\), the endpoint \(\omega^N\) has no limit as
\(N\to\infty\). Indeed its modulus is one; a hypothetical
limit \(L\) would satisfy \(L=\omega L\), forcing \(\omega=1\),
contrary to \(0<\theta<\pi\). Rational \(\theta/(2\pi)\) gives
a nonconstant cycle. The nonconvergence proof covers the other
angles too. The continuous endpoint \(e^{it\theta}\) likewise
does not converge.

If \(\alpha>0\), this eigenbranch grows exponentially; if
\(\alpha<0\), it decays but its reciprocal grows. Both branches
belong to the retained wave. Decay in one branch does not erase
the other branch or its recorded rate.

The single-step record determines the source position:
\[
\boxed{
y_\mu=\frac{2+e^\ell+e^{-\ell}}4,\qquad
\kappa_\mu=-ib\,\frac{e^\ell-1}{e^\ell+1}
=b\tan\frac{\theta-i\alpha}{2},\qquad
\mu=\kappa_\mu^2.
}
\tag{16}
\]
These follow by inverting the defining fractional transformation.
The denominator is nonzero because \(\omega\ne-1\).
Retaining the multiplicities then specifies the source product
response; retaining \(w_\mu(a)\) and \(q(a)\) fixes its normalized
and unnormalized readings. This is an exact identity for the
complete record, not a reconstruction algorithm from the currently
available finite source jet.

The update (12) applies to repetitions of the same step in the
same frame. For different, noncommuting cut matrices the record
must instead retain the ordered product
\(W_N\cdots W_2W_1\) and any frame transitions. Replacing that
product by the exponential of the sum of individual logarithms
requires a separate commutation argument.

## CPL-T4. The full recorded flow is entire and has the exact CWR growth rate

Define the balanced reading
\[
\boxed{
\mathscr P_a(t)=
\sum_\mu w_\mu(a)\cosh(t\ell_a(\mu)).
}
\tag{17}
\]
Each summand is the two-seed reading
\(\operatorname{tr}(g(y_\mu)V_a(t;\mu))\), since \(A\) has
zero diagonal and \(\operatorname{tr}g=1\).

Set
\[
\sigma_a=\sup_\mu|\alpha_a(\mu)|
=\log\mathfrak A(a),\qquad
R_a=\frac{\sqrt a+1/2}{\sqrt a-1/2},
\tag{18}
\]
where \(\mathfrak A(a)\) is CWR-T3's actual reciprocal-wave gain.
Then \(\sigma_a\le\log R_a\), and (17) is entire in complex \(t\).
For real \(t\) it is real and even, with
\[
\mathscr P_a(0)=1,\qquad
\mathscr P_a(N)=r_N(a),\qquad
|\mathscr P_a(t)|\le2e^{\sigma_a|t|}.
\tag{19}
\]
Its exact growth rate is
\[
\boxed{
\limsup_{t\to+\infty}
\frac{\log|\mathscr P_a(t)|}{t}=\sigma_a,
\qquad \log0=-\infty.
}
\tag{20}
\]

**Proof.** The bounds on \(\alpha\) and \(0<\theta<\pi\) give
a uniform bound on \(|\ell|\). For complex \(t\) in a fixed
compact set, \(|\cosh(t\ell)|\le e^{|t||\ell|}\). Together
with (2), this proves absolute local uniform convergence of
the entire series. Equation (7) and conjugate weights make
the result real on the real axis. Evenness and the value at
zero follow directly.

For real \(t\),
\[
|\cosh(t(\alpha+i\theta))|\le\cosh(t\alpha)
\le e^{|\alpha t|}.
\]
This proves the bound in (19). At integers, CWR's
\((\omega^N+\omega^{-N})/2=T_N(2y-1)\) proves equality with
its actual balanced readings \(r_N\).

The upper bound in (20) follows from (19). For the lower
bound, CWR-T3 proves
\(\limsup_{N\to\infty}|r_N(a)|^{1/N}=\mathfrak A(a)\).
Taking its logarithm and using the integer subsequence in
(19) gives the reverse inequality. QED.

Consequently, for any one fixed \(a\ge a_0\),
\[
\boxed{
\mathrm{RH}
\ \Longleftrightarrow\ \sigma_a=0
\ \Longleftrightarrow\
\mathscr P_a\hbox{ is bounded on }[0,\infty)
\ \Longleftrightarrow\
|\mathscr P_a(t)|\le1\quad(t\in\mathbb R).
}
\tag{21}
\]
Indeed \(\sigma_a=0\) is exactly (6) for all source modes.
Under RH, the weights are positive and sum to one, and
\(\mathscr P_a(t)=\sum_\mu w_\mu(a)\cos(t\theta_a(\mu))\),
which is bounded by one. Boundedness implies \(\sigma_a=0\)
by (20).

The bound at all real parameters, or equivalently zero gain
for every source mode, is not proved here. Equation (20)
inherits its global source detection from CWR's uncancelled
poles; the interpolation does not supply the missing sign.

## CPL-T5. Zero mean can coexist with persistent rotational memory

Assume RH for the full actual source in this theorem only.
Then all \(w_\mu(a)>0\), their sum is one, and the distinct
source modes have distinct angles
\[
\theta_\mu=2\arctan(\gamma_\mu/\sqrt a)\in(0,\pi).
\]
The long-flow readings satisfy
\[
\boxed{
\lim_{T\to\infty}\frac1T\int_0^T\mathscr P_a(t)\,dt=0,
}
\tag{22}
\]
but for every retained mode \(\mu\),
\[
\boxed{
\lim_{T\to\infty}\frac1T\int_0^T
 \mathscr P_a(t)e^{-i\theta_\mu t}\,dt
=\frac{w_\mu(a)}2>0.
}
\tag{23}
\]
Moreover
\[
\boxed{
\lim_{T\to\infty}\frac1T\int_0^T
 |\mathscr P_a(t)|^2\,dt
=\frac12\sum_\mu w_\mu(a)^2>0.
}
\tag{24}
\]
In particular \(\mathscr P_a(t)\) has no finite limit as
\(t\to+\infty\), even under RH.

**Proof.** Every \(\theta_\mu\) is nonzero. The mean of
\(e^{i\eta t}\) over \([0,T]\) tends to zero for real
\(\eta\ne0\), by direct integration. Expanding each cosine
into its two exponentials proves (22). When multiplied by
\(e^{-i\theta_\mu t}\), precisely the positive-frequency
half of the \(\mu\) term has zero frequency; all other
frequencies remain nonzero. This proves (23).

For (24), expand the absolutely summable double series.
Distinct positive frequencies have zero limiting product
mean, and each diagonal \(\cos^2(t\theta_\mu)\) has mean
\(1/2\). The weight sums \(\sum w_\mu=1\) and
\(\sum_{\mu,\nu}w_\mu w_\nu=1\) justify all interchanges
by domination. At least one positive weight exists, so the
sum of squares is strictly positive.

If \(\mathscr P_a(t)\to L\), the normalized integral of
\((\mathscr P_a(t)-L)e^{-i\theta_\mu t}\) would tend to zero,
as would the integral of \(L e^{-i\theta_\mu t}\).
This contradicts (23). QED.

The same three conclusions hold, without an RH assumption,
inside any retained positive finite model with distinct
support in \(0<y<1\): replace the source weights by that
model's positive normalized weights and put
\(\theta_j=\arccos(2y_j-1)\). They describe the finite model;
they do not assert that it contains every actual source mode.

Thus a zero long-time average does not mean that the cut
memory has vanished. The frequency-selected readings (23)
retain it, while (24) retains a positive mean square.
Calling that mean square physical energy would require an
additional physical identification, which is not made.

## CPL-T6. Infinitesimal per-step phase retains the source when its scale is carried

For each fixed source mode, as \(a\to\infty\),
\[
\boxed{
\frac{\sqrt a}{2}\ell_a(\mu)
\longrightarrow i\kappa_\mu=\delta+i\gamma.
}
\tag{25}
\]
In particular
\[
\frac{\sqrt a}{2}\alpha_a(\mu)\to\delta,\qquad
\frac{\sqrt a}{2}\theta_a(\mu)\to\gamma.
\tag{26}
\]
Although the unscaled one-step memory tends to zero, its
rescaled rate retains both horizontal displacement and height.

**Proof.** With \(b=\sqrt a\), for each fixed \(\kappa_\mu\)
the branch from (4) agrees, for large b, with
\[
\log\frac{1+i\kappa_\mu/b}{1-i\kappa_\mu/b}
=\frac{2i\kappa_\mu}{b}+O_\mu(b^{-3}).
\]
Multiplication by \(b/2\) proves (25)--(26). QED.

For fixed real \(u\), this gives the modewise scaling
\[
e^{(\sqrt a\,u/2)\ell_a(\mu)}\to e^{i\kappa_\mu u}.
\tag{27}
\]
Passing to a full source sum requires its weights to be
carried from a fixed reference anchor. FWT-T1 already proves
that summable carried limit. The pointwise expansion above
is not uniform over all source heights and is not used as
a substitute for FWT's tail argument.

CWR-T4 also already proves the exact global rate limit
\[
\boxed{
\lim_{a\to\infty}\frac{\sqrt a}{2}\sigma_a
=D_\xi,\qquad D_\xi=\sup_\mu|\delta_\mu|.
}
\tag{28}
\]
Taking the unscaled limit \(\sigma_a\to0\) cannot establish
RH: it holds even if \(D_\xi>0\). The retained scale in
(28), and the carried weights in the full wave limit, are
essential parts of the memory.

## Phase memory, curvature and the remaining proof condition

The record here is the lift of a specified transport phase.
A nonzero rotation along a path does not, by itself, prove
nonzero curvature. Repeating one fixed step gives one
parameter: in the connection convention
\(\nabla=d-B\,dt\), its curvature
\(d(-B\,dt)+(-B\,dt)\wedge(-B\,dt)\) is zero for constant B,
while its transport is still \(e^{tB}\).
A curvature interpretation needs independently specified
directions and a comparison of ordered closed-loop transports.
No such new two-direction law is assumed in this note.

Likewise the phase is not determined by the Jacobian
determinant alone: every matrix \(e^{\theta R}\) has
determinant one, for all angles \(\theta\). A determinant
can retain area change without identifying the rotation.
This is why (12) records the phase itself.

| Question | Exact result |
|---|---|
| What is retained after each repetition? | Log-amplitude and lifted phase, updated by (12) |
| Does the exponential reproduce the original cut-wave? | Yes, at every integer step by (10) |
| Does the record reset after a complete turn? | No; the turn count in (14) is retained |
| Is the infinite-flow endpoint a fixed constant? | Not even in the pure rotation case |
| Can an average vanish while memory survives? | Yes; (22)--(24), under the stated positive-source or finite-model premise |
| What detects possible actual-source growth? | The exact rate \(\sigma_a=\sup|\alpha_a|\), by (20) |
| Does an infinitesimal cut erase the source? | Its unscaled phase shrinks; the carried rate (25) retains it |
| What remains for RH? | Prove \(\sigma_a=0\), equivalently the actual global bound (21) |

The conclusions do not identify zero with a newly proved physical
potential, infer a universal speed, or deduce positive source
pairings from the existence of a phase record. The source windows,
finite certificates and unresolved all-depth positivity from the
preceding work are unchanged.

## Lineage and verification scope

- [CWR-1](RH_CUT_WAVE_RETURN.md), T1--T4: native two-reflection
  wave, exact source multipliers, uncancelled poles, discrete
  growth law and the scaled horizontal defect.
- [LWP-1](papers/native-layer-wave-propagation/THEOREM.md),
  T1 and T4: the skew generator and its separate exponential.
- [CEM-1](RH_CUT_EQUILIBRIUM_MEMORY.md): the earlier cut-return
  memory for that skew-generator flow.
- [ZMR-1](RH_ZERO_MEMORY_RECOVERY.md): recovery of its folded
  orientation channel using reference motion.
- [FWT-1](RH_FAITHFUL_WAVE_TRANSPORT.md), T1--T3: carried
  source-frequency limit and the continuous global-growth criterion.
- [HCM-2](papers/native-heat-observer-memory/PHASE_RECORD.md):
  the separate heat-response record and the warning that a
  finite phase label or seam count alone does not supply a
  continuous lifted angle. No identification with its heat law
  or the legacy UGD numeral product is asserted.

Matrix exponentials, logarithmic lifts and the elementary means
of distinct exponentials are standard tools. Their needed
calculations are displayed here. No literature-wide novelty,
physical interpretation or machine-checked analytic proof is
claimed. The RH criterion is explicitly inherited from the
retained source detector.

Only this note and its README entry change. No scientific
computation, numerical source evaluation, independent agent,
source replay or workflow dispatch is performed. The record
is specified analytically; no finite data are represented as
an acquired infinite spectrum. The commit includes [skip ci].
The change remains a draft stacked on PR #129 and is not merged.
