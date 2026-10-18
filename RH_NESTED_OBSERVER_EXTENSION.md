# NOE-1 — Infinite positive observer banks and source identification

18 October 2026 (Asia/Kolkata). Author: Monty.
Written continuation of PWC-1 at
d6ab4bf3293cfd45e40acacdcde42c1b5050bade.

## Result

An infinite nested bank of positive actual-source readings can be
constructed without assuming RH. Every already positive finite bank
admits every sufficiently distant additional real observation. An
explicit sufficient Schur gate is given below. Old observations,
their derivatives, and all mixed source pairings are retained.

For a nested unbounded bank containing a fixed reference, its positive
rational responses converge locally uniformly on the slit plane.
The normalized probability measures are tight, and their waves converge
locally uniformly to a positive-definite probability wave.

This completes existence of an infinite positive refinement and
convergence of its own model wave. It does not complete identification
with the actual source. In fact, this note constructs an actual-source
bank of that kind whose limiting response provably differs from the
actual response at a specified position. The construction uses exact
xi readings and UCB's source-deficit bound, not a toy substitute.

A source-identifying extension is still needed. The retained Schwarz
argument proves that a positive nested bank with
\(\sum_n t_n^{-1/2}=\infty\) suffices. The far-extension construction
does not prove that density condition and can be chosen to violate it.

These statements do not prove or disprove RH. They distinguish an
infinite positive carrier from a source-complete positive carrier.
No numerical source computation or new machine certificate is used.

## 1. Retained source and finite-bank identities

Keep
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q(z)=\mathcal X'(z)/\mathcal X(z),\qquad
F(w)=wq(w^2).
\tag{1}
\]
For distinct real source positions \(t_i>9/4\), write
\[
u_i=q(t_i),\quad T=\operatorname{diag}(t_i),\qquad
L_{ij}=-q[t_i,t_j],\quad J_{ij}=(tq(t))[t_i,t_j].
\tag{2}
\]
Diagonal divided differences are derivatives. A positive bank means
both actual matrices \(L,J\) are strictly positive definite.

SE48 then gives
\[
q_A(z)=u^T(J+zL)^{-1}u
      =\sum_{j=1}^{|A|}\frac{r_j}{z+\eta_j},
\qquad r_j,\eta_j>0,
\tag{3}
\]
with exact value and derivative interpolation at every point of \(A\).
The finite \(\eta_j\) are model spectral values, not certified xi zeros.

For a new \(x\notin A\), retain the actual \(r=q(x)\), \(d=q'(x)\),
and the mixed columns
\[
l_i=\frac{u_i-r}{x-t_i},\qquad
j_i=\frac{xr-t_i u_i}{x-t_i}.
\tag{4}
\]
SE49's exact extension conditions are
\[
\alpha_A(x)=-d-l^TL^{-1}l>0,\qquad
\beta_A(x)=r+xd-j^TJ^{-1}j>0.
\tag{5}
\]
For real \(z>0\), put
\[
P(z)=J+zL,\quad k(z)=j+zl,\quad
s(z)=r+(x-z)d-k(z)^TP(z)^{-1}k(z),
\]
\[
e(z)=r-k(z)^TP(z)^{-1}u.
\]
Under (5),
\[
\boxed{
q_{A\cup\{x\}}(z)-q_A(z)=\frac{e(z)^2}{s(z)}\ge0,\qquad
s(z)\ge\beta_A(x)+z\alpha_A(x)>0.
}
\tag{6}
\]
These are inherited identities. Positivity of separate banks is not
substituted for positivity of their mixed union.

## NOE-T1. The complete source supplies explicit distant-observer estimates

For \(w\ge64\), set \(\ell(w)=\log(w/(2\pi))\). Then
\[
\boxed{
\left|q(w^2)-\frac{\ell(w)}{4w}\right|\le\frac3{w^2},
}
\tag{7}
\]
\[
\boxed{
\left|-q'(w^2)-\frac{\ell(w)-1}{8w^3}\right|
 \le\frac6{w^4},
\qquad
\left|q(w^2)+w^2q'(w^2)-\frac{\ell(w)+1}{8w}\right|
 \le\frac6{w^2}.
}
\tag{8}
\]

**Proof.** SQR-T1 already proves, for the complete actual source,
\[
d(w):=F(w)-\tfrac14\log(w/(2\pi)),\qquad
|d(w)|\le2/|w|+3\,2^{-\Re w}\quad(\Re w\ge2).
\tag{9}
\]
The rational, Gamma, and full prime-power terms are included.
In particular \(|d(w)|\le3/w\) for real \(w\ge64\).

On the complex circle of radius \(w/2\) about that real point,
\(\Re z\ge w/2\) and \(|z|\ge w/2\). Cauchy's derivative estimate
therefore gives
\[
|d'(w)|\le\frac2w\left(\frac4w+3\,2^{-w/2}\right)
\le\frac9{w^2}.
\tag{10}
\]
The last bound uses \(6w\,2^{-w/2}\le1\) for \(w\ge64\);
the ratio decreases there and the inequality holds at 64.

Now use the exact chain rules
\[
q(w^2)=F(w)/w,\quad
q'(w^2)=\frac{wF'(w)-F(w)}{2w^3},\quad
q(w^2)+w^2q'(w^2)=\frac{F(w)+wF'(w)}{2w}.
\]
Substitution of (9)--(10) proves (7)--(8). QED.

In particular \(q(w^2)\sim\ell(w)/(4w)>0\) and tends to zero.
The derivative bounds were obtained from a complex analytic remainder,
not by differentiating an unspecified real-axis asymptotic error.
The estimates are written consequences of the retained source bound.

## NOE-T2. Every positive finite actual bank has a positive distant extension

Fix a positive bank \(A\) of size \(m\). Define finite constants
\[
t_{\max}=\max_{t_i\in A}t_i,\quad U=\|u\|_2,\quad V=\|Tu\|_2,
\quad C_L=\|L^{-1}\|_2,\quad C_J=\|J^{-1}\|_2.
\]
For \(w\ge64,\ w^2\ge2t_{\max}\), set
\[
r_+(w)=\frac{\ell(w)}{4w}+\frac3{w^2},\qquad
b_L(w)=\frac{2(U+\sqrt m\,r_+(w))}{w^2},
\]
\[
b_J(w)=2\sqrt m\,r_+(w)+\frac{2V}{w^2}.
\]
The following are sufficient explicit gates:
\[
\boxed{
G_L(w):=\frac{\ell(w)-1}{8w^3}-\frac6{w^4}
          -C_L b_L(w)^2>0,
}
\tag{11}
\]
\[
\boxed{
G_J(w):=\frac{\ell(w)+1}{8w}-\frac6{w^2}
          -C_J b_J(w)^2>0.
}
\tag{12}
\]
Both hold for every sufficiently large \(w\). Consequently there
exists \(X_A<\infty\) such that every real \(x>X_A\) extends \(A\)
to a positive actual bank.

**Proof.** Since \(x=w^2\ge2t_{\max}\), every \(x-t_i\ge x/2\).
The actual \(r\) is positive and at most \(r_+(w)\), by T1 and the
retained positive-real source bound. Hence (4) gives
\[
\|l\|_2\le b_L(w),\qquad \|j\|_2\le b_J(w).
\]
Use (8) and the inverse norms in (5). This proves
\(\alpha_A(w^2)\ge G_L(w)\) and \(\beta_A(w^2)\ge G_J(w)\).

As \(w\to\infty\), \(b_L=O_A(w^{-2})\) and
\(b_J=O_A(\ell(w)/w)\). Thus the subtracted terms in (11)
are smaller than its leading \(\ell(w)/(8w^3)\), and those in
(12) are smaller than its leading \(\ell(w)/(8w)\).
Both gates eventually hold for all \(w\), proving the claim. QED.

More precisely the actual innovations obey
\[
\alpha_A(w^2)=\frac{\ell(w)-1}{8w^3}+O_A(w^{-4}),
\qquad
\beta_A(w^2)=\frac{\ell(w)+1}{8w}
             +O_A(\ell(w)^2/w^2).
\tag{13}
\]
The constants and sufficient threshold depend on the existing bank,
including its inverse norms. No bound uniform in bank size is claimed.

The new position can therefore be required to be the square of an
integer above any prescribed threshold. This is an existence theorem
with an explicit sufficient inequality, not a numerically acquired
list of additional observations.

## NOE-T3. A sufficiently distant new observation changes fixed readings arbitrarily little

Fix a positive finite bank \(A\) and a compact real interval
\(K\subset(0,\infty)\). Write
\[
\omega_A(z)=\prod_{t_i\in A}(z-t_i),\qquad
d_A(z)=\prod_j(z+\eta_j).
\]
As \(x\to\infty\) along the positive extensions of T2,
\[
\boxed{
\frac{q_{A\cup\{x\}}(z)-q_A(z)}{q(x)}
\longrightarrow
2\left(\frac{\omega_A(z)}{d_A(z)}\right)^2
\quad\text{uniformly for }z\in K.
}
\tag{14}
\]
In particular
\[
0\le q_{A\cup\{x\}}(z)-q_A(z)
      =O_{A,K}(q(x))\longrightarrow0
\tag{15}
\]
uniformly on \(K\).

**Proof.** With \(r=q(x)\), direct expansion of the exact columns
(4), with the old bank fixed, gives
\[
l=O_A(x^{-1}),\quad j=r\mathbf1+O_A(x^{-1}),\quad
k(z)=r\mathbf1+O_{A,K}(x^{-1}).
\]
The inverse \(P(z)^{-1}\) is uniformly bounded for \(z\in K\).
Moreover \(x^{-1}/r\to0\), by T1. Therefore
\[
e(z)/r\longrightarrow 1-\mathbf1^TP(z)^{-1}u.
\]
The same source estimates give
\[
\frac{r+xq'(x)}r\longrightarrow\frac12,\qquad
q'(x)/r=O(x^{-1}).
\]
Since \(k^TP^{-1}k=O_{A,K}(r^2)\) and \(r\to0\),
\[
s(z)/r\longrightarrow\frac12
\]
uniformly on \(K\). Equation (6) now yields (14), with the
coefficient initially written as
\(2(1-\mathbf1^TP(z)^{-1}u)^2\).

RCD-T1, equation (9), identifies this coefficient. In its rational
identity
\[
1-(z+\mu)f(\mu)^TP(z)^{-1}u
=\frac{\omega_A(z)}{d_A(z)}h_A(\mu),
\]
let the auxiliary \(\mu\) tend to infinity. Then
\((z+\mu)f(\mu)\to\mathbf1\) and \(h_A(\mu)\to1\).
This proves the stated ratio and finishes the proof. QED.

At an old fitted position, \(\omega_A=0\), and the response increment
is in fact exactly zero by interpolation. Formula (14) retains this
memory of the original cut positions.

## NOE-T4. Infinite nested positivity can preserve all old readings and an arbitrarily small response budget

Suppose \(A_0\) is positive and contains a fixed \(a\ge a_0=3249/4\).
For every compact interval \(K\subset(0,\infty)\) and \(\varepsilon>0\),
there is a strictly nested sequence
\[
A_n=A_{n-1}\cup\{x_n\},\qquad x_n\to\infty,
\]
such that both actual pencils are positive at every stage, every
old value and derivative remains matched, and
\[
\boxed{
0\le q_{A_n}(z)-q_{A_0}(z)\le\varepsilon
\quad(z\in K,\ n\ge1).
}
\tag{16}
\]
The new square-root positions can be chosen to double at each step,
so that
\[
\boxed{\sum_{n\ge1}x_n^{-1/2}<\infty.}
\tag{17}
\]

**Proof.** At stage \(n\), T2 allows every sufficiently large next
position. T3 then allows that position to be chosen still larger so
that
\[
\sup_{z\in K}
\bigl(q_{A_n}(z)-q_{A_{n-1}}(z)\bigr)
<\varepsilon\,2^{-n}.
\]
Choose its square root to be an integer at least twice the preceding
largest square root. These finitely many requirements are compatible.
Induction proves both signs at every stage, with all mixed columns
taken from the actual source. Sum the nonnegative increments to
obtain (16). The doubling condition proves (17). QED.

Thus infinite nested positive actual banks exist unconditionally:
UCB supplies a positive initial bank containing a suitable reference.
The theorem does not assert positivity for every sequence of prescribed
positions or for the union of arbitrary separately positive banks.

## NOE-T5. Every such unbounded positive chain has a tight, convergent probability wave

More generally, let \(A_n\) be any nested sequence of positive actual
banks, containing a fixed \(a\ge a_0\), and containing retained positions
\(x_n\to\infty\). Then there is a holomorphic function
\(q_\infty\) on \(\Omega=\mathbb C\setminus(-\infty,0]\) such that
\[
\boxed{q_{A_n}\to q_\infty\quad\text{locally uniformly on }\Omega.}
\tag{18}
\]
For every retained position \(t\in\bigcup_n A_n\),
\[
q_\infty(t)=q(t),\qquad q_\infty'(t)=q'(t).
\tag{19}
\]
At the fixed reference, the positive probabilities
\[
\pi_j^{(n)}
=\frac{r_j^{(n)}}{(a+\eta_j^{(n)})q(a)}
\]
are tight as measures of \(\eta\in[0,\infty)\), and
\[
\boxed{
\mathscr W_{A_n,a}(u)\to W_\infty(u)
\quad\text{uniformly on every bounded real }u\text{-interval}.
}
\tag{20}
\]
The limiting wave is continuous, even, positive definite, equals
one at zero, and has modulus at most one. Its Laplace reading is
\[
\int_0^\infty e^{-wu}W_\infty(u)\,du
=\frac{w}{w^2-a}
 \left(1-\frac{q_\infty(w^2)}{q(a)}\right)
\quad(w>0,\ w\ne\sqrt a).
\tag{21}
\]
The omitted point is removable. These statements do not identify
\(q_\infty\) with \(q\).

**Proof of response convergence.** By (6), \(q_{A_n}(t)\) is
nondecreasing in \(n\) for every real \(t>0\). From the finite
positive representation and \(q_{A_n}(a)=q(a)\),
\[
0<q_{A_n}(t)\le q(a)\max\{1,a/t\}.
\tag{22}
\]
Thus a finite pointwise limit exists on the positive real axis.

For each compact \(K'\subset\Omega\), the finite constant
\[
C_{K'}=\sup_{z\in K',\,\eta\ge0}\frac{a+\eta}{|z+\eta|}
\]
gives \(|q_{A_n}(z)|\le C_{K'}q(a)\). The same analytic
compactness used in SE48 supplies locally uniformly convergent
subsequences. All subsequential limits agree on the positive real
axis by (22), hence agree everywhere on \(\Omega\) by the identity
theorem. This proves full convergence (18). Derivatives converge
locally too, giving (19).

The real functions \(q_{A_n}(t)\) are decreasing in \(t\).
Their limit is also decreasing and positive. By (19) and T1,
\(q_\infty(x_n)=q(x_n)\to0\). Consequently
\[
\boxed{q_\infty(R)\to0\quad(R\to\infty).}
\tag{23}
\]
PWC-T1 therefore bounds the probability tails uniformly in model
index, for \(R>a\):
\[
\sum_{\eta_j^{(n)}>R}\pi_j^{(n)}
\le\frac{2q_{A_n}(R)}{q(a)}
\le\frac{2q_\infty(R)}{q(a)}
\longrightarrow0.
\tag{24}
\]

For completeness, push the probabilities to
\(y=a/(a+\eta)\in[0,1]\). Their moments converge by (18) and its
derivatives at \(a\). Compactness and polynomial density then give
one limiting probability measure \(\nu_\infty\). Its atom at zero
vanishes, by (23) and the response identity
\[
\frac{q_\infty(z)}{q(a)}
=\int_0^1\frac{a\,d\nu_\infty(y)}{a+(z-a)y}.
\tag{25}
\]
Thus no probability is lost at infinite frequency. Push this measure
back to \(\eta=a(1-y)/y\), allowing a possible atom at \(\eta=0\).

Apply PWC-T2's finite wave-error bound, with (24) and convergence
of each fixed moment. First choose the spectral threshold, then
the polynomial degree, then the two model indices. The waves are
uniformly Cauchy on every bounded interval. Their limit is also the
cosine transform of the resulting probability on \(\eta\), and
inherits the stated kernel signs and normalization.

Finally, the finite Laplace identity passes to the limit by the
common bound one, proving (21). QED.

This proves model convergence and tightness without PWC's
actual-source convergence premise. The latter premise is replaced
here by nestedness, an unbounded retained bank, and convergence
to its own response \(q_\infty\). Source identification remains
a separate question.

## NOE-T6. An actual xi bank can have an infinite positive limit that misses a fixed source reading

Choose a reference \(a\) with UCB's
\[
N_*(a)=\left\lfloor\frac{e^{-4}\sqrt a}{4\pi}\right\rfloor\ge128,
\qquad
\ell_a=\log(\sqrt a/(2\pi)),\qquad
\delta_2(a)=\frac6{\ell_a(441/4)^2}.
\]
Take \(A_0=\{a\}\), \(z_*=3a/2\), and
\[
D=\frac{q(a)\delta_2(a)}2>0.
\tag{26}
\]
There exists an infinite nested bank of actual xi readings,
positive at every stage, whose limit from T5 satisfies
\[
\boxed{
q_\infty(a)=q(a),\quad q_\infty'(a)=q'(a),\qquad
q(z_*)-q_\infty(z_*)\ge D/2>0.
}
\tag{27}
\]
All other retained positions and derivatives are also matched
exactly. Its probability wave is tight and positive definite but
differs from FWT's faithful actual wave.

**Proof.** UCB-T2 makes \(A_0\) positive. UCB-T6 has enough budget
for one extra position and gives, at \(x=z_*/a-1=1/2\),
\[
q(z_*)-q_{A_0}(z_*)
\ge q(a)(1+3/2)\delta_2(a)
       \frac{(1/2)^2}{1+(1/2)^2}
=\frac{q(a)\delta_2(a)}2=D.
\tag{28}
\]
Apply T4 with \(K=[a,3a/2]\) and \(\varepsilon=D/2\).
T5 then supplies its limit and
\[
q_\infty(z_*)\le q_{A_0}(z_*)+D/2\le q(z_*)-D/2.
\]
The interpolation statements follow from (19). This proves (27).

For an exact wave witness, put \(w_*=\sqrt{z_*}\).
Subtract (21) from FWT's actual-source identity:
\[
\boxed{
\int_0^\infty e^{-w_*u}
 \bigl(\mathscr W_a(u)-W_\infty(u)\bigr)\,du
=-\frac{w_*[q(z_*)-q_\infty(z_*)]}
        {q(a)(z_*-a)}
\le-\frac{w_*D}{2q(a)(z_*-a)}<0.
}
\tag{29}
\]
Both integrals converge. Hence the waves cannot be identical. QED.

This is an actual-source counterexample to the proposed inference
"infinite positive, stable, tight, exact at every selected observer,
therefore complete source recovery." It is NOT an RH counterexample,
does not exhibit an off-critical zero, and is compatible with RH
being true. Even two positive waves can represent different sources.

No observation is fabricated: each bank uses the original actual
\(q,q'\), and its finite spectral measure is the SE48 interpolant.
The construction is a written existence proof; no new infinite
dataset or finite numerical certificate is claimed.

## NOE-T7. A positive extension with sufficient retained density does identify the source

Let \(A_n\) be nested positive actual banks containing the fixed
\(a\), with distinct retained positions \(t_j>9/4\) tending to
infinity. If
\[
\boxed{\sum_j\frac1{\sqrt{t_j}}=\infty,}
\tag{30}
\]
then
\[
q_\infty=q,\qquad W_\infty=\mathscr W_a,
\qquad \mathrm{RH}.
\tag{31}
\]
This extends the retained SE49 density argument from its specified
integer bank to the displayed class. The density condition is
sufficient; it is not asserted to be necessary for every source.

**Proof.** Put \(b=\sqrt a\), \(s_j=\sqrt{t_j}\), and
\[
C_a=q(a)+\frac{41}{8a}+\frac1{8b}.
\]
The analytic argument of UCB (28) applies to every positive bank
containing \(a\): it uses positive finite residues, the actual
complete-source bound, and exact value/derivative interpolation,
not the moving-window size gate. It gives
\[
\boxed{
|q(w^2)-q_{A_n}(w^2)|
\le C_a\frac{|w+b|^2}{|w|}
\prod_{t_j\in A_n}
 \left|\frac{w-s_j}{w+s_j-2}\right|^2
\quad(\Re w>1).
}
\tag{32}
\]
To recall the premises explicitly, \(wq_{A_n}(w^2)\) has positive
real part on \(\Re w>0\), with its value fixed at \(b\).
The half-plane Schwarz bound controls it by
\(q(a)|w+b|^2/\Re w\). SE49 controls the actual function by
\(41/8+|w|/4\). After division by \(C_a(w+b)^2\), the difference
is bounded by one on \(\Re w>1\); its double interpolation zeros
give the displayed product by successive Schwarz division.

Fix a real \(w>1\). A factor vanishes if \(w=s_j\), already
giving exact matching. Otherwise, after discarding finitely many
terms, \(s_j>w\) and
\[
\left|\frac{w-s_j}{w+s_j-2}\right|
=1-\frac{2(w-1)}{s_j+w-2}.
\]
The logarithm of its square is at most
\(-4(w-1)/(s_j+w-2)\). Assumption (30) makes the sum diverge
to minus infinity. Thus the product in (32) tends to zero.
Every fixed positive-ray response \(R>a\) therefore converges
to \(q(R)\), while normalization at \(a\) is exact.
PWC-T3 proves faithful wave convergence, actual kernel positivity,
and RH. QED.

The extension theorem T2 gives no estimate ensuring (30) through
an unlimited sequence of stages. Its thresholds depend on changing
inverse norms. T4 deliberately chooses a summable reciprocal-root
schedule; T6 proves that this can preserve a nonzero source defect.
One cannot replace (30) by mere unboundedness or by convergence
of the model wave.

## What is completed and what remains

| Requirement | Status after this note |
|---|---|
| Retain a fixed reference and every old reading while extending | Proved |
| Find a positive next actual observation after any positive finite bank | Proved, with explicit sufficient gates |
| Obtain an infinite nested positive bank | Proved as a written existence result |
| Converge to a holomorphic positive-model response | Proved |
| Prevent probability escaping to infinite frequency | Proved for every unbounded nested positive chain here |
| Converge to a stable positive probability wave | Proved |
| Identify that limit with the complete actual xi source | Not automatic; T6 proves an actual-source non-recovery example |
| Identify the source under positive extension plus (30) | Proved implication |
| Construct an unlimited positive extension satisfying (30), or another source-identifying condition | Open |
| Prove RH | Not completed |

The next concrete target is an estimate on the actual mixed source
innovations that permits sufficient retained sampling density.
This is more specific than assuming positivity of a limiting wave.
The source-identifying property cannot be inferred just by naming
a residue probability or by observing stable infinite flow.

## Dependencies, attribution, and verification

- [SQR-1](RH_SQUARE_ROOT_BANKS.md), T1: complete-source complex remainder.
- [SE48](papers/native-prime-observation-closure/sources/SE48.md),
  T1, T3: finite positive realization, interpolation, analytic compactness.
- [SE49](papers/native-moving-zero-covariance/sources/SE49.md),
  T1--T3: two actual Schur innovations, response update, source
  majorant, and the original fixed integer-bank density criterion.
- [UCB-1](RH_OBSERVER_BANK_UNION.md), T2, T6: positive initial
  source banks, strict actual deficit, and Schwarz remainder.
- [RCD-1](RH_RESIDUAL_COMPLETION.md), T1: exact residual-observer
  factorization used in the distant-update coefficient.
- [PWC-1](RH_PROBABILITY_WAVE_CLOSURE.md), T1--T3:
  probability tail, finite wave-error bound, and faithful limit criterion.
- [FWT-1](RH_FAITHFUL_WAVE_TRANSPORT.md), T2--T3, T5:
  actual Laplace recovery, RH criterion, and finite probability waves.

Cauchy's estimate, Schur elimination, normal-family compactness,
polynomial density, and the half-plane Schwarz argument are standard
analytic tools within the admitted representation. Their use is
explicit; no new native primitive or literature-wide novelty is claimed.
The new written deductions are the distant-observer gates, small-update
extension, infinite nested completion, and actual-source non-recovery
construction.

Validation uses written derivations and exact text/metadata checks only.
No scientific computation, numerical source evaluation, independent
agent, certificate replay, or workflow dispatch is performed.
Only this note and its README entry change. The commit uses [skip ci].
Published as an unmerged draft stacked on #131.
