# PWC-1 — Probability tightness and faithful wave closure

18 October 2026 (Asia/Kolkata). Author: Monty Dabas.
Written continuation at Publications commit
18f9ba3ae536deff7fcba0fc6a15deccbb4a325c.

## Result and claim boundary

The positive finite responses already have probability waves. The missing
step is identification of their limit with the complete actual source.
This note gives that step under a precise source-convergence premise,
with an explicit finite wave-error bound.

At one fixed retained anchor, convergence of the normalized responses
on the positive ray forces tightness of the model probability weights,
convergence of every retained layer moment, and locally uniform
convergence to the faithful actual wave. Its positive-definite kernels
then follow by passage to the limit. By FWT-T3 this implies RH.

The source-convergence premise is not established here. It is itself
RH-equivalent when the models are the existing positive actual-source
interpolants. The presently proved moving-window budgets do not
establish it: their Schwarz products stay bounded away from zero at
each fixed test point. This last statement concerns that estimate,
not an impossibility theorem for the actual approximation error.

The new written contributions are the explicit probability-tail and
wave-error bounds, their connection to fixed-anchor normalization,
and this quantitative audit of the available moving-window estimate.
The normal-family/positive-measure closure principles are standard;
SE48 and SE49 already give conditional response-convergence criteria.
No literature-wide novelty, numerical certificate, or RH proof is claimed.

## 1. Source, finite models, and separate parameters

Keep the actual source and faithful wave of FWT:
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q(z)=\mathcal X'(z)/\mathcal X(z)
     =\sum_\mu\frac{\nu_\mu}{z+\mu},
\]
\[
\mu=(\gamma-i\delta)^2,\quad |\delta|<1/2,\quad
a\ge a_0=3249/4,\quad q(a)>0,
\]
\[
\mathscr W_a(u)=\frac1{q(a)}
 \sum_\mu\frac{\nu_\mu}{a+\mu}
 \cos((\gamma-i\delta)u).
\tag{1}
\]
The actual source can have nonreal \(\mu\); no positive actual spectral
measure is assumed. FWT supplies continuity, reality, evenness,
\(\mathscr W_a(0)=1\), and the bound \(2e^{|u|/2}\).

For any finite bank \(A\) whose actual pencils \(L_A,J_A\) are strictly
positive, SE48 supplies
\[
q_A(z)=\sum_j\frac{r_j^A}{z+\eta_j^A},\qquad r_j^A,\eta_j^A>0.
\tag{2}
\]
At the SAME fixed reference \(a\), whether or not \(a\in A\), define
\[
\pi_j^A=\frac{r_j^A}{(a+\eta_j^A)q_A(a)},\quad
y_j^A=\frac a{a+\eta_j^A},\quad
D_A(R)=\frac{q_A(R)}{q_A(a)},
\]
\[
c_k^A=\sum_j\pi_j^A(y_j^A)^k
 =\frac{(-a)^kq_A^{(k)}(a)}{k!\,q_A(a)},\qquad
\mathscr W_{A,a}(u)=\sum_j\pi_j^A\cos(\sqrt{\eta_j^A}\,u).
\tag{3}
\]
These weights sum to one. Each wave is bounded by one and is positive
definite, as proved in FWT-T5. The positive finite carrier is earned
from the actual finite pencils; it is not imposed on the actual
infinite source.

Here \(R\) is a positive spectral threshold/source test position,
\(k\) is layer-moment order, \(m\) below is an approximation degree,
and \(u\) is FWT's auxiliary wave coordinate. They are not a single
clock or a common layer index. In particular \(m\) is not bank size.
Model frequencies are not identified with actual xi zero ordinates.

## PWC-T1. A normalized response controls probability escaping to high frequency

For every positive finite model and \(R>a\),
\[
D_A(R)=\sum_j\pi_j^A\frac{a+\eta_j^A}{R+\eta_j^A},
\qquad
\boxed{\sum_{\eta_j^A>R}\pi_j^A
 \le\frac{2R}{a+R}D_A(R)\le2D_A(R).}
\tag{4}
\]
Thus a reading at \(R\) controls probability at wave frequencies
larger than \(\sqrt R\).

**Proof.** Substitute (3) into (2). For \(\eta>R\), the factor
\((a+\eta)/(R+\eta)\) is at least \((a+R)/(2R)\);
it is increasing in \(\eta\) because \(R>a\). Sum its positive
contribution to obtain (4). QED.

The actual response independently satisfies
\[
\boxed{\lim_{R\to\infty}q(R)=0.}
\tag{5}
\]
Indeed \(\Re\mu=\gamma^2-\delta^2>-1/4\). For \(R\ge a\), each
\(|R+\mu|\) is nondecreasing in \(R\), so
\[
|R+\mu|^{-1}\le |a+\mu|^{-1}.
\]
FWT's retained sum \(\sum_\mu\nu_\mu/|a+\mu|<\infty\) dominates
the series; every term tends to zero. Dominated convergence proves
(5), with no RH assumption.

Consequently, if a sequence of positive models has
\[
D_{A_n}(R)\longrightarrow q(R)/q(a)\quad\text{for every }R>a,
\tag{6}
\]
then
\[
\lim_{R\to\infty}\limsup_{n\to\infty}
 \sum_{\eta_j^{A_n}>R}\pi_j^{A_n}=0.
\tag{7}
\]
This is tightness of the probabilities at the fixed reference.
It is a consequence of the stated source-convergence premise,
not a consequence of probability normalization alone.

## PWC-T2. Finite layer errors and a tail reading bound the whole wave on an interval

For two positive finite models \(A,B\), set
\[
\Delta_m(A,B)=\max_{0\le k\le m}|c_k^A-c_k^B|,
\quad d=\frac a{a+R},\quad
L_{R,T}=\frac{aT^2}{2d^2}.
\]
For every \(R>a,\ T\ge0,\ m\ge1\),
\[
\boxed{
\sup_{|u|\le T}|\mathscr W_{A,a}(u)-\mathscr W_{B,a}(u)|
\le
4\bigl(D_A(R)+D_B(R)\bigr)
+\frac{L_{R,T}}{\sqrt m}
+3^m\Delta_m(A,B).
}
\tag{8}
\]
This estimate does not assume RH, nested banks, common atoms, or
positivity of an unacquired union. It only uses the two stated
positive finite models.

**Proof.** Push each probability to \(y=a/(a+\eta)\in(0,1)\).
Clip the wave frequency at \(\sqrt R\):
\[
h_u(y)=
\cos\!\left(u\sqrt{\min\{R,a(1-y)/y\}}\right)\quad(y>0),
\qquad h_u(0)=\cos(u\sqrt R).
\]
The difference from the original integrand is at most two, and
occurs only at \(\eta>R\). By (4) the resulting error is at most
\(4D_A(R)\) for \(A\), and at most \(4D_B(R)\) for \(B\).

The clipped function is constant on \([0,d]\). On \([d,1]\), put
\(v=a(1-y)/y\); differentiation and
\(|\sin(u\sqrt v)|\le |u|\sqrt v\) give
\[
|h_u'(y)|\le\frac{a u^2}{2y^2}\le L_{R,T}.
\]
The endpoint at \(y=1\) is read by continuity. Thus \(h_u\) is
globally Lipschitz with this constant, uniformly for \(|u|\le T\).

Its Bernstein polynomial
\[
B_mh_u(y)=\sum_{\ell=0}^m h_u(\ell/m)
 {m\choose\ell}y^\ell(1-y)^{m-\ell}
\]
has uniform error at most \(L_{R,T}/(2\sqrt m)\).
To see this directly, write the polynomial as the expectation of
\(h_u(K/m)\) for a binomial \(K\) with parameters \(m,y\);
Cauchy--Schwarz bounds the error by
\(L_{R,T}\sqrt{y(1-y)/m}\).
The two integrals therefore contribute at most \(L_{R,T}/\sqrt m\).

Expand \(B_mh_u(y)=\sum_{k=0}^m b_k(u)y^k\).
Since \(|h_u|\le1\),
\[
\sum_{k=0}^m|b_k(u)|
\le\sum_{\ell=0}^m {m\choose\ell}2^{m-\ell}=3^m.
\]
The difference of its two probability integrals is at most
\(3^m\Delta_m(A,B)\). Combining the three errors proves (8). QED.

The constants are explicit but not optimized. Increasing \(m\)
reduces the approximation term and amplifies moment error;
one cannot discard the \(3^m\) factor. A fixed finite list of
matched moments alone does not control the full unbounded-frequency
wave. The tail reading is the additional information needed by (8).

## PWC-T3. Source-response convergence forces faithful wave convergence and positivity

Fix \(a\ge a_0\). Suppose \(A_n\) is a sequence of actual finite
banks, each with \(L_{A_n},J_{A_n}\succ0\), for which (6) holds.
Then
\[
c_k^{A_n}\longrightarrow
c_k(a):=\frac{(-a)^kq^{(k)}(a)}{k!\,q(a)}
\quad\text{for every fixed }k,
\tag{9}
\]
and
\[
\boxed{\mathscr W_{A_n,a}\longrightarrow\mathscr W_a
\quad\text{uniformly on every bounded real }u\text{-interval}.}
\tag{10}
\]
In particular every finite translation kernel of the actual wave
is positive semidefinite, \(|\mathscr W_a(u)|\le1\), and RH follows
by FWT-T3. The premise (6) has NOT been proved for an unlimited
positive sequence in this note.

**Proof of the moment identification.** Write
\(\nu_n=\sum_j\pi_j^{A_n}\delta_{y_j^{A_n}}\) as a probability on
the compact interval \([0,1]\). Every subsequence has a weakly
convergent further subsequence. This is ordinary compactness of
probability measures on a compact interval; equivalently, select
convergent polynomial moments and extend by uniform polynomial
approximation.

For a subsequential limit \(\nu\), (6) gives, for every \(x>0\),
\[
\int_0^1\frac{d\nu(y)}{1+xy}
=\frac{q(a(1+x))}{q(a)}.
\tag{11}
\]
For \(|x|<1\), the left side has its uniformly convergent geometric
series. The right side is analytic near zero, since \(q(a)\ne0\)
and \(a\) is a regular source point. Equality for small positive
\(x\) identifies all coefficients, giving
\(\int y^k\,d\nu=c_k(a)\).
Probability measures on \([0,1]\) with all the same moments agree:
polynomials are uniformly dense in continuous functions.
Thus all subsequential limits are identical, and the full sequence
converges weakly, proving (9).

Moreover, letting \(R=a(1+x)\to\infty\) in (11), bounded
convergence on its left and (5) on its right give
\[
\boxed{\nu(\{0\})=0.}
\tag{12}
\]
The compactified endpoint \(y=0\), representing infinite frequency,
has acquired no probability atom. No lower probability floor is
introduced.

**Proof of wave convergence.** Given \(T\) and an error tolerance,
first choose \(R\) large using (5)--(6), then choose \(m\) large,
and finally choose both model indices large using (9).
Inequality (8) makes the model waves uniformly Cauchy on
\([-T,T]\). Their limits for different \(T\) agree, defining a
continuous real even function \(W\), with \(W(0)=1\) and
\(|W|\le1\). Every finite positive-definite kernel passes to the
limit.

For real \(w>\sqrt a\), the finite version of FWT-T2 is
\[
\int_0^\infty e^{-wu}\mathscr W_{A_n,a}(u)\,du
=\frac{w}{w^2-a}\bigl(1-D_{A_n}(w^2)\bigr).
\]
The common bound one justifies dominated convergence on the left;
(6) gives on the right the exact actual transform
\[
\int_0^\infty e^{-wu}W(u)\,du
=\frac{w}{w^2-a}\left(1-\frac{q(w^2)}{q(a)}\right)
=\int_0^\infty e^{-wu}\mathscr W_a(u)\,du.
\tag{13}
\]
Both functions are continuous and of exponential order.
Uniqueness of the Laplace transform implies \(W=\mathscr W_a\)
on the positive half-line, and evenness gives the whole line.
This proves (10) and the asserted positivity, then FWT-T3 gives
RH. No actual positive source measure was assumed in this
identification. QED.

Conversely, under RH, SE49's fixed integer bank
\(t_j=(j+3/2)^2\) has both positive actual pencils at every order;
its responses converge locally uniformly to \(q\) away from the
nonpositive real axis. Therefore (6) holds at any fixed \(a\).
Accordingly,
\[
\boxed{\mathrm{RH}\iff
\text{there exists a sequence of positive actual-source finite
models satisfying (6) at one fixed }a\ge a_0.}
\tag{14}
\]
This is a proved closure criterion, not a verified instance of
its universal premise.

Under that premise, letting the second model index in (8) tend
to infinity also gives the conditional quantitative bound
\[
\boxed{
\sup_{|u|\le T}|\mathscr W_{A_n,a}(u)-\mathscr W_a(u)|
\le4\left(D_{A_n}(R)+\frac{q(R)}{q(a)}\right)
+\frac{L_{R,T}}{\sqrt m}
+3^m\max_{k\le m}|c_k^{A_n}-c_k(a)|.
}
\tag{15}
\]
At fixed \(R,m\), its first two terms need not vanish. The order
of choices in the proof, or a justified diagonal schedule, matters.

## PWC-T4. Returning to a fixed reference retains a normalization residual

Let \(R_A(z)=q(z)-q_A(z)\). For real \(R>a\), exact subtraction gives
\[
\boxed{
D_A(R)-\frac{q(R)}{q(a)}
=\frac{q(R)R_A(a)-q(a)R_A(R)}
       {q(a)q_A(a)}.
}
\tag{16}
\]
If \(|R_A(a)|\le e_a<q(a)\) and \(|R_A(R)|\le e_R\), then
\[
\boxed{
\left|D_A(R)-\frac{q(R)}{q(a)}\right|
\le\frac{q(R)e_a+q(a)e_R}{q(a)(q(a)-e_a)}.
}
\tag{17}
\]
Here \(q(R)>0\) is the retained positive real source reading.
Both statements follow by substituting \(q_A=q-R_A\); the
lower bound on \(q_A(a)\) proves (17).

If \(a\in A\), the first residual vanishes by interpolation.
If the bank moves to larger and larger source positions, the
fixed \(a\) need not remain in it. Its normalization residual
must then be controlled too. LAT's exact transport within an
individual finite model does not establish either residual's
decay against the complete source.

These formulas give a concrete sufficient route to (6): obtain
a positive sequence and prove \(e_a\to0,\ e_R\to0\) at every
fixed \(R>a\). They do not supply those estimates themselves.

## PWC-T5. The proved moving-window size budget does not make its Schwarz product vanish at a fixed position

Use UCB's common anchor \(B=b^2\), distinguishing it from the
fixed reference \(a\). Assume UCB's \(N_*(B)\ge128\). Its positive-bank range is
\[
A\subset[B/2,3B/2],\qquad
|A|=n\le M(B)\le C_0 b,\qquad C_0=\frac{e^{-4}}{48\pi}.
\]
For a fixed real test \(w\ge2\), UCB's actual-source error
estimate contains
\[
P_A(w)=\prod_{t_i\in A}
 \left|\frac{w-\sqrt{t_i}}{w+\sqrt{t_i}-2}\right|^2.
\]
Whenever \(b/\sqrt2\ge4(w-1)\),
\[
\boxed{P_A(w)\ge
\exp\!\bigl[-8\sqrt2\,C_0(w-1)\bigr]>0.}
\tag{18}
\]

**Proof.** Put \(s_i=\sqrt{t_i}\ge b/\sqrt2>w\). Each unsquared
factor equals \(1-x_i\), where
\[
x_i=\frac{2(w-1)}{s_i+w-2}\le\frac{2\sqrt2(w-1)}b\le\frac12.
\]
For \(0\le x\le1/2\), \(\log(1-x)\ge-2x\), by differentiation
or integration of \(1/(1-x)\le2\). Therefore
\[
\log P_A(w)\ge-4\sum_i x_i
\ge-\frac{8\sqrt2(w-1)n}{b}
\ge-8\sqrt2 C_0(w-1),
\]
which proves (18). QED.

In UCB (28), with \(t_*\in A\), the full upper estimate is
\[
|R_A(w^2)|\le
\left(q(t_*)+\frac{41}{8t_*}+\frac1{8\sqrt{t_*}}\right)
\frac{(w+\sqrt{t_*})^2}{w}\,P_A(w).
\tag{19}
\]
The prefactor is at least \(\sqrt{t_*}/(8w)\), since \(q(t_*)>0\).
Together with (18), the numerical expression on the right of
(19) in fact tends to infinity as \(b\to\infty\) at this fixed
\(w\). It cannot establish decay to zero there.

A lower bound on an UPPER estimate is not a lower bound on the
actual error. This proves neither failure of actual convergence
nor failure of the framework. It identifies why the currently
available estimate, despite proving positivity for larger banks,
does not prove the fixed-reference premise (6).

SE49 uses a different arrangement: one fixed nested bank with
enough retained sampling density. Its Schwarz product vanishes,
but positivity of every prefix remains unproved. Combining these
two distinct results without their respective hypotheses would
silently omit exactly that unresolved step.

## PWC-T6. Wave error gives a finite kernel bound; signed residual is not probability

Suppose a positive finite-model wave obeys the independently
established estimate
\[
\sup_{|u|\le T}|\mathscr W_a(u)-\mathscr W_{A,a}(u)|\le\epsilon.
\]
For observer positions with \(|u_i-u_j|\le T\), and complex
coefficients \(z_i\),
\[
\boxed{
\sum_{i,j}\bar z_i z_j\mathscr W_a(u_i-u_j)
\ge-\epsilon\left(\sum_i|z_i|\right)^2.
}
\tag{20}
\]
Indeed the finite-model form is nonnegative and the error form
is bounded in absolute value by the right-hand error term.
For \(r\) observers this also gives the matrix lower bound
\(-r\epsilon I\). A positive error tolerance does not prove
exact positive semidefiniteness; a controlled zero-error limit
does, as in T3.

Probability, waves, and residuals therefore have linked but
different roles:

| Object | Established mathematical role |
|---|---|
| \(\pi_j^A\) | Nonnegative weights summing to one in each acquired positive model |
| \(\mathscr W_{A,a}(u)\) | A probability cosine transform; its value may be negative |
| Translation matrix \((\mathscr W_{A,a}(u_i-u_j))\) | Positive semidefinite correlation kernel |
| \(R_A=q-q_A\) | Complete source-response residual, retained in (16) |
| \(\mathscr W_a-\mathscr W_{A,a}\) | Signed wave discrepancy; FWT-T6 proves both signs on its declared finite range |
| Pole residue of \(q\) | Actual zero multiplicity, not a normalized probability |

No replacement of a bilinear source square by an absolute square
has occurred. No intensity measurement has replaced the signed
phase-sensitive wave. The probability interpretation is acquired
where positivity is proved, and the infinite interpretation requires
the source identification in T3.

## Remaining work and source attribution

The concrete outstanding task is to construct an unlimited family
of positive actual-source finite models and prove (6), or prove
another already established RH-equivalent actual-source sign
condition. T1--T3 remove a separate ambiguity about losing
probability at infinite frequency once (6) is supplied. They do
not make that RH-equivalent premise automatic.

- [FWT-1](RH_FAITHFUL_WAVE_TRANSPORT.md), T1--T3, T5--T6:
  faithful actual wave, source recovery, pole noncancellation,
  exact RH criterion, positive finite waves, signed discrepancy.
- [UCB-1](RH_OBSERVER_BANK_UNION.md), T2, T4, T6:
  finite positive-bank budget, compatible compressions and
  complete-source Schwarz remainder.
- [RCD-1](RH_RESIDUAL_COMPLETION.md), T1--T3:
  full residual factorization, retained off-line poles and
  next-layer defect.
- [SE48](papers/native-prime-observation-closure/sources/SE48.md),
  T1, T3: positive finite realization and the earlier
  conditional response-limit criterion.
- [SE49](papers/native-moving-zero-covariance/sources/SE49.md),
  T2--T3: finite source-error estimate and fixed integer-bank
  convergence under all-order positivity.
- [LAT-1](papers/native-layer-anchor-transport/THEOREM.md):
  exact positive-model anchor transport with retained weights.
- [NIST DLMF 1.14(iii)](https://dlmf.nist.gov/1.14#iii):
  standard Laplace-transform conventions and uniqueness.

Weak compactness of probabilities on a compact interval, polynomial
density, Bernstein approximation, and Laplace uniqueness are standard
analytic tools of the admitted representation, not newly claimed
native axioms. The quantitative estimates are proved explicitly above.

Validation is written algebra and analytic reasoning plus
text/metadata checks. No scientific computation, numerical source
evaluation, independent agent, certificate replay, or workflow
dispatch is performed. Only this note and its README entry change;
engines, frozen evidence and workflows are unchanged. The commit
uses [skip ci]. This continuation is published as an unmerged
draft stacked on #130.
