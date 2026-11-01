# YM90 — Gauge-selected dynamic memory of the actual vacuum

YM89 bounded the actual product-link marginal and its compressed diagonal
generator, but left the size of the dynamic Schur memory open. This note
closes that memory estimate in the same proved small-theta window. It uses
internal gauge symmetry, the unretained-coordinate conditional gap, and
the actual mixed log-Hessian estimate from YM88. The bound is uniform in
the finite volume and the straight-path block factor.

This is a strong-coupling lattice result. It does not construct a
four-dimensional continuum theory, discharge a weak-coupling estimate,
or improve the numerical full-gap bound already obtained in YM88.
The analytic proof below is accompanied by exact finite controls; those
controls are not a formal machine verification of the analytic proof.

## 1. Statement and conventions

Use the YM62 native SU(2) Hamiltonian on the periodic spatial cubic lattice,

\[
 A=-\Delta-\theta\sum_p W_p,\qquad
 Ah=E_0h,\quad h>0,\quad \chi=\log h,\quad d\mu=h^2d\Phi.
\]

Normalize h in L2(Phi). The positive ground-transformed generator is
\(\mathcal L=-\Delta-2\nabla\chi\cdot\nabla\). Native link generators have
fundamental Casimir 3/4 and Ricci tensor (1/2)g. Take disjoint straight
b-link retained paths, b >= 2, on a fine torus of side Lb with L >= 4.
Write R for their fine links and Z for all unretained links. Internal
vertices of these paths have retained degree two and belong to only one
path. No arbitrary fine-coordinate pins are inserted.

Let pi_b be the ordered product on each path, w_b its actual mu marginal,
I_b f = f composed with pi_b, P = I_b I_b*, and Q = 1-P. In the orthogonal
coarse/hidden decomposition, write the closed forms as

\[
 \mathcal L=\begin{pmatrix}C_b&B^*\\B&D\end{pmatrix},\qquad
 C_b=b[-w_b^{-1}\operatorname{div}(w_b\nabla)],\quad B=Q\mathcal L I_b.
 \tag{1}
\]

Here B maps coarse functions to the hidden space. This is the transpose
convention to the finite matrix B in the YM64 runner. Define

\[
 \Sigma_b(s)=B^*(D+s)^{-1}B,\qquad K_b(t)=B^*e^{-tD}B.       \tag{2}
\]

All comparisons with C_b are quadratic-form comparisons on its energy
domain. In particular B is not asserted to be bounded on coarse L2.

**Theorem.** For the actual Wilson ground state, if 0 <= theta <= 1/400,
put

\[
 d=\tfrac12-80\theta,\qquad
 \sigma^2=\frac{100\theta^2}{d},\qquad
 \eta=\frac{100\theta^2}{d^2}\le\frac1{144}.
 \tag{3}
\]

Then the range of B is internal-gauge invariant, the restriction D_0 of D
to that invariant hidden sector has lower bound d, and

\[
 \boxed{\|Bf\|^2\le\sigma^2 C_b[f]},\qquad
 0\le\Sigma_b(s)\le\alpha(s)C_b,\quad
 \alpha(s)=\frac{\sigma^2}{d+s}\quad(s\ge0).                \tag{4}
\]

Moreover, for t >= 0 and integers n >= 0,

\[
 0\le K_b(t)\le\eta d e^{-dt}C_b,\qquad
 0\le M_n:=B^*D^{-(n+1)}B\le\eta d^{-n}C_b.                \tag{5}
\]

For s > 0 define R_b(s) = I_b*(L+s)^(-1)I_b, with L here denoting
\(\mathcal L\). The exact Schur resolvent and the diagonal approximation
obey

\[
 \begin{split}
 (C_b+s)^{-1}&\le R_b(s)\le[(1-\alpha(s))C_b+s]^{-1},\\
 0\le R_b(s)-(C_b+s)^{-1}
   &\le\frac{\alpha(s)}{1-\alpha(s)}(C_b+s)^{-1}
    \le\frac1{143}(C_b+s)^{-1}.                         \tag{6}
 \end{split}
\]

The same statements hold at s = 0 after restricting to centered coarse
functions. Thus 1/143 is a relative Loewner form error for the resolvent;
it is not an entrywise, trajectory, or semigroup error bound. An operator
exponential is not generally order preserving, so (6) alone does not
provide a corresponding ordering of reduced semigroups.

## 2. Gauge coordinates and the selected hidden gap

Let G be the product of SU(2) gauge groups at all retained-path internal
vertices, fixing the coarse endpoints and all other vertices. Its action
fixes every coarse product. Gauge invariance of h and of the native
product metric makes the action unitary on L2(mu), commuting with L.
The group average A_G is an orthogonal projection with
P A_G = A_G P = P. Consequently it commutes with Q and with the hidden
form D; B has range in the A_G-invariant subspace.

There is a global, Haar-preserving change of variables

\[
 (q_R,q_Z)\longleftrightarrow(Q_{\rm coarse},g,z').        \tag{7}
\]

Indeed, successive internal prefix gauges transform each retained path
to (1,...,1,Q_E). Apply the same gauges to its adjacent unretained links.
The internal prefixes and Q_E form Haar product coordinates on each
retained path. For fixed retained coordinates, left and right
translations of each unretained link preserve Haar. Hence the full
change has product Haar measure. Gauge invariance gives exactly

\[
 d\mu(\cdot\mid Q_{\rm coarse})=d\Phi_G(g)\,d\nu_Q(z'),
 \quad d\nu_Q(z')\propto h(\operatorname{section}(Q),z')^2d\Phi_Z(z').
 \tag{8}
\]

This is not a claim that the transformed fine metric is a product metric.
Only the probability disintegration is a product as written. In
particular, an invariant hidden v is a function V(Q,z') with
\(\int V(Q,z')d\nu_Q=0\).

Assume, for the moment, that every coordinate conditional nu_Q has
Poincare floor d > 0 in the native Z energy. Then

\[
 \|v\|^2\le d^{-1}\int\sum_{z\in Z}|\nabla_{z'}V|^2d\mu
          \le d^{-1}\mathcal E(v).                       \tag{9}
\]

For the second step, at fixed q_R the changes of q_Z are native
isometries; the other fine gradient energies are simply discarded.
This proves D_0 >= d. It does not identify D_0 with a conditional
Z generator and does not assume the full fine Hamiltonian gap.
At each fixed finite volume, smooth positive h and the product-map
coordinates make P bounded on H1; no uniform bound on that projection
is needed for this core statement. Smooth approximation followed by Q
and gauge averaging gives a smooth invariant hidden core, extending
(9) to the closed hidden form domain.

For completeness, the other internal gauge sectors also have a uniform
floor. A vertex generator is a signed sum of at most six incident native
link generators. Each link meets at most two internal vertices. Thus

\[
 \sum_{v\in G,a}|G_{v,a}u|^2\le12\sum_{e,a}|D_{e,a}u|^2,
 \quad\|u-A_Gu\|^2\le\tfrac43\sum_{v,a}\|G_{v,a}u\|^2
                  \le16\mathcal E(u).                  \tag{10}
\]

The middle inequality is product-group Haar Poincare with native gap
3/4, integrated against the invariant measure mu. Therefore all of D
has floor min(d,1/16). The stronger floor d is available in (4)--(5)
because B only couples to D_0. Physical Gauss-invariant states are
inside this invariant sector.

## 3. Conditional score variance bounds the dynamic coupling

For a fine e on path E, let O_e be the orthogonal adjoint map of the
ordered prefix before e, and set

\[
 S_E(q)=\sum_{e\in E}O_e(q_R)\nabla_e\chi(q).              \tag{11}
\]

The native product rule and Laplacian intertwining give

\[
 \mathcal L I_bf=I_b(-b\Delta f)-2\sum_E S_E\cdot I_b\nabla_E f.
 \tag{12}
\]

Comparison with (1), or marginal differentiation under Haar integration,
gives \(\mathbb E_\mu[S_E\mid Q]=b\nabla_E(\tfrac12\log w_b)\).
Consequently

\[
 Bf=-2\sum_E(S_E-\mathbb E[S_E\mid Q])\cdot I_b\nabla_E f.
 \tag{13}
\]

The vectors S_E are internal-gauge invariant: the prefix adjoint and the
fine gradient transform by inverse rotations. This also follows by
applying invariance to the carré du champ of chi and arbitrary product
pullbacks. It is consistent with the abstract range selection above.

Define the nonnegative rectangular cross-Hessian matrix

\[
 H_{ze}=\sup_q\|\operatorname{Hess}_{ze}\chi(q)\|,
 \qquad z\in Z,\ e\in R,\qquad \|H\|_{2\to2}\le a.       \tag{14}
\]

For fixed Q and coarse covectors u_E, a Z derivative of
\(\sum_E u_E\cdot S_E\) differentiates only chi: the prefix frames
depend on retained coordinates only. On the section, or in any rotated
fine frame, copying u_E to its b disjoint links has squared norm
\(b\sum_E|u_E|^2\). The block operator norm and (14) therefore imply

\[
 \sum_z\left|\nabla_{z'}\sum_Eu_E\cdot S_E\right|^2
       \le ba^2\sum_E|u_E|^2.                            \tag{15}
\]

There is no extra factor of three: H already uses the norm on each
three-component native tangent block. Apply conditional Poincare in
(8), put u_E = grad_E f(Q), and integrate Q. Since
\(C_b[f]=b\int\sum_E|\nabla_E f|^2w_b d\Phi\), this proves

\[
 \boxed{\|Bf\|^2\le\frac{4a^2}{d}C_b[f]}.               \tag{16}
\]

The factor b cancels against the exact coarse kinetic coefficient.
The formula first holds for smooth f and extends B continuously from
the coarse energy domain to hidden L2. The positive hidden lower bound
then makes the form Schur construction of YM64 applicable at s >= 0
(centered at zero). Spectral calculus in D_0 proves (4)--(5) with
sigma^2 = 4a^2/d. Inverse monotonicity in the exact form identity
\(R_b(s)=(C_b+s-\Sigma_b(s))^{-1}\) proves (6).

## 4. Insert the actual Wilson estimates

YM88 proves the following pointwise bounds for the actual chi, uniformly
in volume and the entire configuration:

\[
 \sup\|\operatorname{Hess}_{ee}\chi\|\le35\theta,
 \quad\sup_e\sum_{f\ne e}\sup\|\operatorname{Hess}_{ef}\chi\|
       \le5\theta.                                      \tag{17}
\]

The off-diagonal matrix is symmetric, so its column bound is also
5 theta. Every cross submatrix (14) thus has a <= 5 theta. Restricting
the same pointwise weighted Ricci bound to any set Z with R fixed gives

\[
 \operatorname{Ric}_Z-2\operatorname{Hess}_{ZZ}\chi
       \ge(\tfrac12-80\theta)I=dI.                       \tag{18}
\]

Compact weighted Bochner and the spectral theorem give the conditional
Poincare inequality used in (9) and (16), for every section value Q.
This uses the geometry in (17), not YM88's conclusion that the full
operator is gapped. It avoids assuming the conclusion when assembling
coarse and hidden estimates.

At theta = 1/400 the exact constants are

| Quantity | Bound |
| --- | --- |
| Cross-Hessian norm a | 1/80 |
| Conditional and selected-hidden floor d | 3/10 |
| Source energy factor sigma squared | 1/480 |
| Relative zero-frequency memory eta | 1/144 |
| Relative resolvent form error | 1/143 |

All these constants improve towards the free case, where B and the
memory vanish. They do not shrink with b under this estimate; the
exponentially small static score from YM89 alone is insufficient to
deduce exponentially small dynamic memory.

YM89 supplies a centered coarse floor c_b = b lambda_b with
lambda_b >= 227/480. The relative memory bound also gives the physical
two-block gap estimate

\[
 \operatorname{gap}_{\rm physical}(\mathcal L)
 \ge\frac{c_b+d-\sqrt{(c_b-d)^2+4\eta c_bd}}2
 \ge\frac{c_bd(1-\eta)}{c_b+d}.                           \tag{19}
\]

To check it, write the two diagonal energies as x squared and y squared.
The cross energy is bounded by 2 sqrt(eta)xy, using
\(\|D^{-1/2}Bf\|^2\le\eta C_b[f]\).
Subtracting lambda times the norm leaves a positive two-by-two form
when (1-lambda/c_b)(1-lambda/d) >= eta. Equation (19) follows. On the
unrestricted fine Hilbert space also take the minimum with 1/16 from
(10). The already proved direct YM88 floor d is better numerically;
the contribution here is control of the projected dynamics and a
noncircular way to assemble its blocks.

## 5. All-tail and clock consequences

For the expansion \(\Sigma_b(s)=\sum_{n\ge0}(-s)^nM_n\), (5) gives
energy-relative norm convergence for |s| < d, with remainder after n=N
at most

\[
 \eta\frac{(|s|/d)^{N+1}}{1-|s|/d}.                       \tag{20}
\]

For every real s >= 0, finite spectral division instead gives a signed
Loewner remainder:

\[
 0\le(-1)^{N+1}\left[\Sigma_b(s)-\sum_{n=0}^N(-s)^nM_n\right]
 \le\eta\frac{(s/d)^{N+1}}{1+s/d}C_b.                    \tag{21}
\]

The integrated time tail is bounded by eta exp(-dT) C_b. These are
full spectral-tail bounds, not certificates made by dropping a positive
tail. No commutation of different M_n is assumed or needed.

If a clock multiplier kappa > 0 replaces L by kappa L, then
C, B, D scale by kappa,
\(\Sigma_\kappa(s)=\kappa\Sigma(s/\kappa)\), and
\(K_\kappa(t)=\kappa^2K(\kappa t)\). The dimensionless eta is unchanged
and the selected hidden floor becomes kappa d. These identities do not
supply a physical lattice-spacing scaling or a continuum limit.

## 6. Nonzero memory despite an exactly Haar coarse marginal

An additional native example isolates what the theorem controls.
It is an engineered parent Hamiltonian, not the original cubic Wilson
model. Use five links with endpoints (0,1), (1,2), (0,3), (3,2), (0,2).
The first two paths have products Q=q_0q_1 and R=q_2q_3; z=q_4 is hidden.
Set

\[
 W=\operatorname{Re}(Qz^{-1}),\quad h_\tau\propto e^{\tau W},
 \quad V_\tau=\tfrac94\tau W-\tfrac34\tau^2(1-W^2).
 \tag{22}
\]

Since -Delta W = (9/4)W and Gamma(W,W) = (3/4)(1-W squared), h_tau is
the exact positive zero-energy ground state of -Delta-V_tau. Integrating
z makes the coarse (Q,R) marginal exactly product Haar for every tau.
For the physical coarse probe f=Re(Q R inverse), however,

\[
 C_b[f]=\tfrac34,\qquad
 Bf=-\tau\,\mathbf u\cdot\mathbf v,\quad
 u=Qz^{-1},\ v=QR^{-1}.                                  \tag{23}
\]

Bold symbols denote imaginary quaternion vectors. At the section point
q_0=q_2=q_3=z=1, q_1=(3/5,4/5,0,0), this source is -16 tau/25.
It is -8/125 at tau=1/10. It remains gauge invariant at all four vertices.

For 0 <= tau < 1, the hidden coordinate floor is at least (1-tau)/2;
only the two mixed entries linking z with q_0,q_1 are present, each at
most tau/4. Thus a squared <= tau squared/8. In particular at tau=1/10,
d=9/20, sigma squared=1/90, and eta=2/81. Exact Haar integration gives

\[
 \|Bf\|^2=\frac{\tau^2}{4}\frac{F(\tau)}{Z(\tau)}>0\quad(\tau>0),
 \quad Z(\tau)=\sum_{n=0}^\infty\frac{\tau^{2n}}{n!(n+1)!},
 \quad F(\tau)=\frac32\sum_{n=0}^\infty\frac{\tau^{2n}}{n!(n+2)!}.
 \tag{24}
\]

These formulas follow from the semicircle distribution of the scalar
coordinate of a Haar SU(2) quaternion. They integrate all representation
content of the exponential. The certificate encloses each positive
series with its partial sum plus a geometric bound on the entire tail;
it performs no finite-spin replacement.

Separately, the unchanged YM64 native control witnesses a nonzero
first-order source in the actual Wilson model, for a physical coarse
plaquette with b=3. Its coefficient at the pinned rational configuration
is -32/75. That finite-volume perturbative witness does not supply a
uniform Taylor remainder; the upper bound (4) uses the full actual
vacuum and needs no perturbation truncation.

## 7. Failure controls, reproduction, and remaining gate

The exact controls include noncommuting matrix moments, full inverse
comparisons, the selected gauge sector, native Gauss identities, the
block-copy norm, the entire positive-series tail in (24), clock scaling,
and the source-pin chain. One scalar counterexample guards a useful
distinction: C=1, D=100, B=5 gives Sigma(0)=1/4. With a valid but coarse
hidden floor d=1, this does not imply Sigma(1) <= (1/4)/(1+1): the actual
value is 25/101. The stronger source bound (16), rather than only a
zero-frequency memory bound, supplies the frequency decay in (4).

Run with Python 3.12 without optimization:

```sh
python -B papers/yang-mills-certified-benchmark/certificates/ym90_gauge_memory_bound.py --check
python -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym90_gauge_memory_bound.py -v
```

The workflow replays the unchanged YM89 static-block, YM88 pointwise,
and YM64 vacuum-adapted certificates, and the YM89 and YM64 tests. The
manifest pins 403 upstream files; the result records those, four new
local inputs, and its manifest, for 408 input hashes. Default and --check
are read-only. Only explicit --write regenerates this result and digest.

Conditional expectations and Mori-Zwanzig memory are established methods;
see, for example, T. Hudson and X. H. Li,
[Coarse-graining of overdamped Langevin dynamics via the Mori-Zwanzig formalism](https://arxiv.org/abs/1810.08175).
The development here is the gauge-selected bound (16) and its insertion
of the pinned actual-vacuum estimates, not a novelty claim for the
general projection formalism.

The next sufficient criterion is concrete: at the desired running
coupling, prove an actual unretained conditional floor d_Z > 0, an actual
cross-Hessian norm a with 2a < d_Z, and a surviving coarse centered floor.
Then (16)--(19) close that dynamic step. This criterion is sufficient,
not necessary; failure of these particular bounds proves no physical
gap failure. The present estimates do not establish those quantities in
weak coupling, do not iterate such a bound along a UV trajectory, and
do not establish nontrivial continuum observables or the required
four-dimensional continuum axioms. Those tasks remain open.
