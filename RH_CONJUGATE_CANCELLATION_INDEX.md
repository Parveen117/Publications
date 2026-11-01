# CCI-1 — Conjugate cancellation and the source negative index

17 October 2026. Written analytic continuation of CSS-1 at
`c14cba503632c327ba60df6c06a7a98d84b5a1b5`.

The proposal is to call a reading real when its imaginary contributions
cancel. This is a valid description of the resulting scalar. The actual
arithmetic source already has that property. The question for RH is what
happens to the source relations after this cancellation.

There is an exact product defect. Its simplest form is
\[
\frac{z^2+\bar z^2}{2}
-\left(\frac{z+\bar z}{2}\right)^2=-(\Im z)^2.
\]
Thus cancellation can give real readings while preserving a negative
real record of the nonreal input. This note derives the correctly
weighted identity for the actual source and proves a stronger
observer result:

> The supremum of the negative eigenvalue counts of the actual finite
> Hankel banks equals the number of distinct off-critical zero
> quartets. Every finite selection of such quartets, if present, has
> an explicit finite family of real polynomial observers whose
> complete source matrix is negative definite.

The old finite-detector result already isolates one hypothetical
pair, with an explicit sufficient degree. The extension here is
simultaneous isolation of arbitrarily many pairs, an explicit
negative-definite margin, and the exact count of negative directions.
No off-line arithmetic zero is exhibited. Proving that this count is
zero remains the RH task.

No scientific computation, source replay, independent agent or
workflow dispatch is used. Existing source packets and quantitative
positive ranges remain unchanged.

## 1. Keep the actual arithmetic source

Use the source and normalization of
[CSS-1](RH_CURVATURE_SQUARE_SYMMETRY.md) and
[OQD-1](RH_OBSERVER_QUOTIENT_DESCENT.md). For a fixed finite
\(a\ge a_0=3249/4\),
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,
\]
\[
\delta=\beta-\tfrac12,\quad \mu=(\gamma-i\delta)^2,\quad
u_\mu=\frac a{a+\mu},\quad
w_\mu=\frac{\nu_\mu}{q(a)(a+\mu)}
     =\frac{\nu_\mu u_\mu}{a q(a)}.
\tag{1}
\]
The support \(S_a=\{u_\mu\}\) is infinite and conjugation symmetric;
its only accumulation point is zero. Its values and weights are
nonzero, and
\[
q(a)>0,\quad |u_\mu|<1,\quad \sum_\mu|w_\mu|\le2,\quad
\sum_\mu w_\mu=1,\quad \Re w_\mu>0.
\tag{2}
\]
A real support point comes from \(\delta=0\), hence belongs to
\((0,1)\) and has a positive real weight. Indeed
\[
\Im u_\mu=\frac{2a\gamma\delta}{|a+\mu|^2},
\qquad \gamma>0.
\tag{3}
\]
All these properties are retained at every admitted anchor.

On the real polynomial carrier \(V=\mathbb R[y]\), define
\[
\mathcal L_a[p]=\sum_{u\in S_a}w_u p(u),\qquad
B_a(p,h)=\mathcal L_a[ph],\qquad
H_N(a)=(c_{i+j}(a))_{0\le i,j<N},
\quad c_j=\mathcal L_a[y^j].
\tag{4}
\]
Each value in (4) is real by conjugate cancellation. Products remain
bilinear; no absolute-square replacement of the source is made.

The real parts in (2) follow directly from
\(a+\Re\mu=a+\gamma^2-\delta^2>0\). For useful convergence control,
\[
\left|\frac{\Im w_\mu}{\Re w_\mu}\right|
=\frac{2\gamma|\delta|}{a+\gamma^2-\delta^2}
\le\frac{1}{2\sqrt{a-1/4}}.
\tag{5}
\]
The last inequality uses \(|\delta|<1/2\) and
\(a-1/4+\gamma^2\ge2\gamma\sqrt{a-1/4}\).

## CCI-T1. Weighted cancellation retains an exact product defect

Choose one nonreal support point \(z\) with \(\Im z>0\), and write
\[
w_z=\alpha_z+i\beta_z,\qquad \alpha_z>0,\qquad
\theta_z=\beta_z/\alpha_z.
\]
Normalize the conjugate pair's linear reading by its positive
total weight:
\[
\mathcal C_z[p]
=\frac{w_zp(z)+\bar w_zp(\bar z)}{2\alpha_z}
=\Re p(z)-\theta_z\Im p(z).
\tag{6}
\]
For every real polynomial \(p\), this is real and
\(\mathcal C_z[1]=1\). For any two real polynomials,
\[
\boxed{
\mathcal C_z[ph]-\mathcal C_z[p]\mathcal C_z[h]
=-(1+\theta_z^2)\Im p(z)\Im h(z).
}
\tag{7}
\]
In particular
\[
\boxed{
\mathcal C_z[p^2]-\mathcal C_z[p]^2
=-(1+\theta_z^2)(\Im p(z))^2.
}
\tag{8}
\]
Taking \(p(y)=y\) makes this strictly negative for every nonreal
pair.

**Proof.** Write \(p(z)=r+it\) and \(h(z)=s+iv\), with real
\(r,s,t,v\). The reading of their product is
\(rs-tv-\theta_z(rv+st)\).
The product of their readings is
\(rs-\theta_z(rv+st)+\theta_z^2tv\).
Subtract to obtain (7), and set \(h=p\) for (8). QED.

Equation (8) is a negative variance-like defect of a real, normalized
linear functional. It shows precisely why this pair functional is
not a positive probability expectation. Reality of its output does
not remove its nonreal source dependence.

For the actual full source, organize the support into real points
and conjugate pairs, indexed by \(g\). Set
\[
\pi_g=
\begin{cases}
w_u,&g=\{u\},\ u\in\mathbb R,\\
2\alpha_z,&g=\{z,\bar z\},
\end{cases}
\qquad
m_g(p)=
\begin{cases}
p(u),&g=\{u\},\\
\mathcal C_z[p],&g=\{z,\bar z\}.
\end{cases}
\]
Then \(\pi_g>0,\ \sum_g\pi_g=1\), and the original bilinear form
has the exact decomposition
\[
\boxed{
B_a(p,h)=P_a(p,h)-M_a(p,h),
}
\tag{9}
\]
\[
P_a(p,h)=\sum_g\pi_g m_g(p)m_g(h),
\]
\[
M_a(p,h)=
\sum_{\Im z>0}2\alpha_z(1+\theta_z^2)
                    \Im p(z)\Im h(z).
\tag{10}
\]
Both \(P_a\) and \(M_a\) are positive semidefinite real forms.
The negative term is retained source information, not a chosen sign
convention. The sums converge absolutely: polynomial values are
bounded on the unit disk, (5) bounds all \(\theta_z\), and the
positive group weights sum to one.

With \(\bar m(p)=\sum_g\pi_gm_g(p)=\mathcal L_a[p]\),
the corresponding full-source identity is
\[
\boxed{
\mathcal L_a[p^2]-\mathcal L_a[p]^2
=\sum_g\pi_g\bigl(m_g(p)-\bar m(p)\bigr)^2-M_a(p,p).
}
\tag{11}
\]
Consequently variation between group readings can mask the negative
within-pair term for a selected observer. The next theorem constructs
real polynomial observers that separate any finite collection of
hypothetical nonreal pairs.

The source \(B_a\) has not been replaced by \(P_a\). Such a replacement
would erase \(M_a\) and alter precisely the relations being tested.
Likewise, global multiplicativity
\(\mathcal L_a[ph]=\mathcal L_a[p]\mathcal L_a[h]\) is not the desired
condition: even a positive source with many real modes need not be
multiplicative.

## CCI-T2. Simultaneous finite isolation gives a negative-definite observer bank

Suppose \(z_1,\ldots,z_r\) are any \(r\ge1\) distinct actual support
points with positive imaginary part. They specify \(r\) different
conjugate pairs. Define
\[
s_*=\min_{1\le j\le r}|z_j|>0,\qquad \eta=s_*/2,\qquad
E=\{u\in S_a:|u|\ge\eta\},\qquad d=|E|.
\tag{12}
\]
The set \(E\) is finite, is closed under conjugation, and contains
all \(z_j,\bar z_j\). For each \(v\in E\), use its Lagrange polynomial
\[
L_v(y)=\prod_{\substack{u\in E\\u\ne v}}\frac{y-u}{v-u}.
\tag{13}
\]
All denominators are nonzero because \(E\) contains distinct
source points. For every integer \(n\ge0\), put
\[
\boxed{
p_{j,n}(y)=i\,y^n
 \left(\frac{L_{z_j}(y)}{z_j^n}
       -\frac{L_{\bar z_j}(y)}{\bar z_j^n}\right).
}
\tag{14}
\]
These are real polynomials of degree at most \(n+d-1\), with exact
interpolation readings
\[
p_{j,n}(z_k)=i\,\mathbf1_{j=k},\qquad
p_{j,n}(\bar z_k)=-i\,\mathbf1_{j=k},\qquad
p_{j,n}(u)=0
\quad(u\in E\setminus\{z_j,\bar z_j\}).
\tag{15}
\]
Complex evaluations in (15) do not make the observer coefficients
complex. The conjugation symmetry of \(E\) gives
\(\overline{L_z(\bar y)}=L_{\bar z}(y)\), which makes (14) a
real-coefficient polynomial.

Define the finite constants
\[
C_j=
\sum_{v\in\{z_j,\bar z_j\}}
 \prod_{\substack{u\in E\\u\ne v}}
       \frac{\eta+|u|}{|v-u|},
\qquad
S_C=\sum_{j=1}^r C_j^2,\qquad
\alpha_*=\min_{1\le j\le r}\Re w_{z_j}>0.
\tag{16}
\]
On the remaining support, where \(|y|<\eta\),
\[
|p_{j,n}(y)|\le C_j\,2^{-n}.
\tag{17}
\]
Let the actual \(r\)-observer bank be
\[
G_{r,n}(a)=
\bigl(B_a(p_{i,n},p_{j,n})\bigr)_{1\le i,j\le r}.
\]
Then
\[
\boxed{
G_{r,n}(a)=-2\operatorname{diag}(\Re w_{z_1},\ldots,\Re w_{z_r})
          +R_{r,n},
\qquad
\|R_{r,n}\|\le2S_C\,4^{-n}.
}
\tag{18}
\]
In particular, for
\[
\boxed{
n\ge n_0:=
\max\!\left(0,
\left\lceil\log_4\frac{2S_C}{\alpha_*}\right\rceil\right),
\qquad
G_{r,n}(a)\preceq-\alpha_* I_r\prec0.
}
\tag{19}
\]
The full monomial source bank \(H_{n+d}(a)\) therefore has at least
\(r\) negative eigenvalues.

**Proof.** Interpolation at the finite set \(E\) proves (15).
For \(|y|\le\eta\), each Lagrange numerator in (13) is bounded by
the corresponding product of \(\eta+|u|\). Both selected points
have modulus at least \(s_*=2\eta\). Equation (14) therefore
gives (17).

The pair \(z_j,\bar z_j\) contributes
\(w_{z_j}i^2+\bar w_{z_j}(-i)^2=-2\Re w_{z_j}\) to its own
diagonal reading. Different selected observers have zero cross
reading on \(E\), and all other points of \(E\) vanish.

For a real vector \(x=(x_1,\ldots,x_r)\), the remaining source obeys
\[
\begin{aligned}
|x^T R_{r,n}x|
&\le\sum_{u\notin E}|w_u|
       \left|\sum_jx_jp_{j,n}(u)\right|^2\\
&\le 2\,4^{-n}\left(\sum_j|x_j|C_j\right)^2\\
&\le2S_C\,4^{-n}\|x\|^2.
\end{aligned}
\tag{20}
\]
This uses absolute values only to bound the error in the unchanged
bilinear source. Conjugation symmetry makes \(R_{r,n}\) real
symmetric, so (20) proves its norm bound. Condition (19) makes
the error at most \(\alpha_*\), while the diagonal part is at
most \(-2\alpha_*I_r\).

The polynomials (14) are linearly independent, as their selected
evaluations in (15) show. Their coefficient matrix is a real
inclusion into the polynomials of degree below \(n+d\), and their
bank is its congruence compression of \(H_{n+d}\).
A negative-definite \(r\)-dimensional restriction requires at least
\(r\) negative eigenvalues of the larger real symmetric matrix.
QED.

The margin in (19) concerns coefficients in the specified observer
basis. It is not a lower bound of the same size for the negative
eigenvalues in the unnormalized monomial basis.

The bound depends on the selected hypothetical pairs, their weights,
and the entire finite outer set \(E\). No such off-line source data
are evaluated or asserted to exist. It supplies a conditional
finite-degree witness, not a universal fixed-size RH decision test.
The one-pair detector and its earlier explicit bound remain intact.

## CCI-T3. The negative index counts exactly the off-critical quartets

Define the negative index of the polynomial source by
\[
\operatorname{ind}_{-}B_a
=\sup\{\dim U:\ U\subset\mathbb R[y]\text{ is finite dimensional
 and }B_a|_U\text{ is negative definite}\}.
\]
It takes values in \(\{0,1,2,\ldots,\infty\}\). Let
\(n_-(H_N(a))\) count the strictly negative eigenvalues of the
finite real symmetric matrix, including their matrix multiplicities.
Every finite polynomial subspace has bounded degree, hence
\[
\operatorname{ind}_{-}B_a=\sup_N n_-(H_N(a)).
\tag{21}
\]

Let \(N_{\rm off}\) be the number of distinct nontrivial zeros with
\(\beta>1/2,\ \gamma>0\), possibly infinite. Each is the representative
of one off-critical quartet under conjugation and the functional
equation. Zero multiplicity is retained in \(w_\mu\); it does not
count additional distinct evaluation directions. Then
\[
\boxed{
\operatorname{ind}_{-}B_a
=\sup_N n_-(H_N(a))
=N_{\rm off}
\qquad\text{for every finite }a\ge a_0.
}
\tag{22}
\]

**Proof of the upper bound.** If \(N_{\rm off}=m<\infty\), consider
the real map
\[
T:\mathbb R[y]\longrightarrow\mathbb R^m,\qquad
Tp=(\Im p(z_1),\ldots,\Im p(z_m)).
\]
On its kernel the term \(M_a\) in (9) vanishes, so
\(B_a(p,p)=P_a(p,p)\ge0\). Any subspace of dimension greater than
\(m\) has a nonzero intersection with \(\ker T\), and therefore
cannot be negative definite. Thus the index is at most \(m\).
For \(m=0\), the original weights and support are already real
positive, and this conclusion is immediate.

**Proof of the lower bound.** For any finite selection of \(r\)
different nonreal pairs, T2 supplies a negative-definite polynomial
subspace of dimension \(r\). If there are finitely many \(m\)
pairs, select all of them. If there are infinitely many, \(r\)
is arbitrary and the index is infinite.

Finally, (3) identifies nonreal support pairs with off-critical
positive-height zero pairs. The maps from such zeros to \(\mu\)
and from \(\mu\) to \(u_\mu\) are injective: the square root
\(\gamma-i\delta\) has positive real part, and the fractional
source map is invertible. This proves the quartet count. QED.

For a finite \(N_{\rm off}=m\), the finite counts eventually
stabilize:
\[
n_-(H_N(a))=m\quad\text{for all sufficiently large }N.
\tag{23}
\]
Indeed the counts are nondecreasing under nested polynomial
subspaces, the upper bound is \(m\), and T2 attains \(m\) at
one finite degree. For infinitely many off-critical quartets the
finite counts are unbounded.

This is the negative index of a real quadratic form. It is not a
Fredholm index or an application of an index theorem with missing
operator-domain hypotheses.

## CCI-T4. The count is invariant under faithful observer and anchor changes

For an invertible real change of a finite observer basis,
\[
H\longmapsto R^THR,\qquad n_-(R^THR)=n_-(H).
\tag{24}
\]
To see this without a spectral completion, \(R\) bijects subspaces
with negative-definite restrictions in the two coordinate systems.
Thus the maximal such dimension is unchanged. A noninvertible
compression can lose directions; it is not covered by this claim.

Equation (22) also gives
\[
\boxed{
\operatorname{ind}_{-}B_a=\operatorname{ind}_{-}B_b
=N_{\rm off}\qquad(a,b\ge a_0).
}
\tag{25}
\]
Finite bank sizes or eigenvalues need not agree at different
anchors. The complete source count is unchanged because the same
arithmetic zero pairs remain distinct and nonreal under the
anchor maps.

OQD's inherited cut quotient is isomorphic to \((V,B_a)\) with its
original pairing. Its negative index is therefore also
\(N_{\rm off}\). In OQD's notation, if \(Jp=(p,p)\), then
\(PJ=J\) and
\(\mathcal G_a(Jp,Jh)=B_a(p,h)\).
Consequently the entire negative subspace constructed in T2 would
already lie on the retained cut. Deleting the complementary channel
would not eliminate it while preserving these readings.

There is a direct continuation of the previous curvature result.
CSS defines
\[
k(y)=4y(1-y),\quad t(y)=2y-1,\quad
d_j=\mathcal L_a[k^j],\quad
h_j=\frac{(2c_1-1)d_j-2a\,\partial_a d_j}{2j+1},
\]
\[
\mathcal J_N=
\begin{pmatrix}
(d_{i+j})&(h_{i+j})\\
(h_{i+j})&(d_{i+j}-d_{i+j+1})
\end{pmatrix}_{0\le i,j<N}.
\tag{26}
\]
Its polynomial basis is \(1,k,\ldots,k^{N-1},t,tk,\ldots,tk^{N-1}\).
CSS's invertible coefficient map gives
\(\mathcal J_N=S_N^T H_{2N}S_N\). Therefore
\[
\boxed{
\sup_N n_-(\mathcal J_N(a))=N_{\rm off}.
}
\tag{27}
\]
Even matrix sizes suffice for the supremum because they contain
every finite polynomial degree. The curvature square and its
properly normalized first derivative retain the same complete
negative-direction count. Keeping just the square can fold
different modes together, as CSS's exceptional-circle analysis
already showed.

## 2. Consequence for the proposed reality condition

The distinction is now exact:

| Proposed condition | Actual-source consequence |
|---|---|
| Every summed polynomial reading is real | Already true by conjugate cancellation |
| Every finite source bank is real symmetric | Already true; it may still have negative directions |
| A pair's cancellation preserves its polynomial product readings | Its defect (7) must vanish; a nonreal pair fails already for \(p=h=y\) |
| The complete source has no negative observer directions | By (22), \(N_{\rm off}=0\), exactly RH |
| A finite number of observer banks are positive | Excludes negative directions only in those finite carriers |

The third row is a condition on the resolved pair reading (6).
It does not assert that the full many-mode source should become
multiplicative. The fourth row concerns every polynomial observer,
not just a chosen set of real outputs.

Thus imaginary cancellation is compatible with a faithful definition
of a real reading. It cannot itself establish RH: the cancelled
component survives in real product relations. The construction in
T2 shows that, if off-critical pairs exist, their negative records
cannot all be hidden by the other source terms.

## 3. Retained positive progress and the next proof obligation

RBP already proves actual positivity through
\[
N\le n_*(a)=
\left\lfloor\frac{e^{-4}\sqrt a}{4\pi}\right\rfloor
\]
in its stated range. Accordingly every hypothetical witness from
T2 would have to satisfy
\[
n_0+|E|>n_*(a).
\tag{28}
\]
If its conditional source data instead satisfied the reverse
inequality, the selected configuration would contradict the
retained positive theorem. No new arithmetic configuration is
evaluated or excluded numerically here.

The next source-specific obligation is to prove
\(\operatorname{ind}_{-}B_a=0\), equivalently all-order positivity,
from the actual arithmetic source. This note identifies exactly
what must be eliminated and proves finite simultaneous
detectability. It does not establish that the count vanishes.

| Item | Status |
|---|---|
| Actual weighted cancellation and its product defect | Exact algebraic identity |
| Positive-minus-negative decomposition retaining the original source | Written proof, with convergence control |
| Simultaneous detection of any finite set of hypothetical off-line pairs | Written proof with finite degree and negative margin |
| Equality of source negative index and distinct off-critical quartet count | Written proof |
| Invariance under faithful observer coordinates and source anchor | Derived |
| Actual off-critical zero, negative source matrix, or numerical witness | None asserted or evaluated |
| All-order positivity / RH | Open |

The construction uses finite interpolation, elementary quadratic-form
inertia, and the retained actual-source estimates. These methods
retain their established status. No literature-wide novelty,
optimal witness degree or machine-checked analytic proof is claimed.
No Hilbert-space completion is needed as a primitive.

Only this note and its README entry are added. No predecessor
source, certificate, engine or workflow is changed. RBP's growing
relational range, PRW's unrestricted window, and FCP's restricted
window remain unchanged. The documentation commit uses
`[skip ci]`.

## References

- [CSS-1 — curvature-square symmetry](RH_CURVATURE_SQUARE_SYMMETRY.md), for the actual weight real parts, the folded square's ambiguity, and complete first-derivative reconstruction.
- [D1–D5 — finite observer detector](RH_FINITE_OBSERVER_DETECTOR.md), for the predecessor one-pair determinant, explicit isolation depth, and one-tower RH criterion.
- [OQD-1 — observer quotient descent](RH_OBSERVER_QUOTIENT_DESCENT.md), for source faithfulness and the cut's exact retained polynomial pairing.
- [CWR-1 — cut-wave return](RH_CUT_WAVE_RETURN.md), for the finite-anchor source estimates and observer transport.
- [SPT-1 — shunya potential](RH_SHUNYA_POTENTIAL.md), for the invariant two-channel potential and its bilinear source meaning.
- [RBP-1 — relational bank positivity](RH_RELATIONAL_BANK_POSITIVITY.md), for the unchanged finite positive range.
