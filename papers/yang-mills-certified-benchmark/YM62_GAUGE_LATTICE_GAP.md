# YM-62: native links, Gauss symmetry and a 3+1 lattice gap

Monty Dabas. 4 October 2026. **Python 3.12 only.**

This continuation constructs a locally gauge-invariant model on three
spatial lattice directions with continuous heat time. In the normalization
specified below, a sufficient strong-coupling window is

\[
 \boxed{|\theta|<1/8000.}                                      \tag{1}
\]

Every finite cubic box or periodic cube has a unique vacuum and an
interacting gap at least gamma(theta)>0, independent of its size. The
same bound holds in the Gauss-invariant sector. Subsequences of periodic
vacuum histories give nontrivial infinite-volume, reflection-positive
histories whose reconstructed vacuum Hamiltonian has the same lower gap.
The spatial lattice spacing stays fixed. This is not a continuum R4
Yang--Mills construction or the Clay mass-gap theorem.

**Model choice.** We retain the native compact quaternion carrier,
positive reference, ordered multiplication and derivative calculus.
The independent isotropic electric operator and closed plaquette
interaction below are new declared choices. They do not equal YM61's
correlated, rank-deficient chain operator. Its unresolved interacting
lower gap is not silently discharged by this result.

## Definition and lineage

Let Lambda be a cubic subset of Z3 with free spatial boundary, or a
periodic cube of side L>=4. Put one native unit quaternion q_e in Q on
each positively oriented link e=(x,x+mu). Reverse orientation means
q_(y,x)=q_(x,y)^(-1). There are three spatial directions; heat time is
the fourth coordinate. No fourth spatial direction is added.

Use the YM50 native completion, represented as L2(Q^E,Phi_E), where
Phi_E is the inherited product reference. On a link,
D_a f is differentiation along q -> exp(-t e_a/2)q. Define

\[
 \begin{split}
 L_E&=-\sum_{e\in E}\sum_{a=1}^3 D_{a,e}^{\,2},\\
 W_{x,\mu\nu}(q)&=\operatorname{Re}
   (q_{x,\mu}q_{x+\mu,\nu}q_{x+\nu,\mu}^{-1}q_{x,\nu}^{-1}),\quad\mu<\nu,\\
 V&=\sum_{p\in P}W_p,\qquad A_\Lambda(\theta)=L_E-\theta M_V.
 \end{split}                                                  \tag{2}
\]

There is **no factor 1/2** in L_E. Its one-link coordinate rate is 3/4,
and its product vacuum is 1. Thus b=|theta||P| bounds the potential
norm. Replacing -theta W by theta(1-W) only shifts the energy.

This electric-Casimir/ordered-plaquette Hamiltonian belongs to the
established Hamiltonian lattice-gauge construction; see
[Kogut and Susskind, Phys. Rev. D 11 (1975), 395](https://doi.org/10.1103/PhysRevD.11.395).
The work here supplies its native realization and an explicit conservative
gap budget from YM53's temporal-block argument. It is not a new claim
that strong-coupling lattice gauge theory is gapped. YM53-T1--T5 supplies
the heat, bridge and spectral arguments; the changed incidence and
refinement constants are derived here. YM60-T2 supplies the compactness
method, with the additional uniform-gap passage proved in T5 below.

Reflection-positive reconstruction uses the established positive-time
quotient construction; see
[Jorgensen, Neeb and Olafsson, SIGMA 12 (2016), 058, Remark 2.8(b)](https://sigma-journal.com/2016/058/sigma16-058.pdf).
Time reflection positivity alone does not establish spatial Euclidean
invariance, relativistic axioms or a time-zero Markov limit.

## YM62-T1: local gauge symmetry and the physical sector

For independent g_x in Q at the vertices, set

\[
 q_{xy}\longmapsto g_xq_{xy}g_y^{-1}.                         \tag{3}
\]

The ordered product around a closed face transforms by conjugation at
its base vertex. Its real part is invariant. An open-link overlap,
such as q_e dot q_f, generally is not invariant under (3). Permuting
the factors of the face generally spoils the cancellation. This is
where native noncommutativity is essential.

The electric Casimir commutes with both left and right link translations.
For a fixed colour tensor in a quadratic one-link operator, invariance
under every source-vertex rotation requires that tensor to be scalar:
R C R^T=C for every adjoint rotation R implies C=cI. This explains the
isotropic choice in (2); it is not an assertion about theories with
additional transforming tensor fields.

Let D_a^R differentiate q -> q exp(-t e_a/2). The local Gauss derivation is

\[
 G_{x,a}=\sum_{e\text{ outgoing at }x}D_{a,e}
          -\sum_{e\text{ incoming at }x}D^R_{a,e}.            \tag{4}
\]

On coefficient polynomials, G_(x,a)W_p=0 and [G_(x,a),L_E]=0.
The group action is unitary for Phi_E and commutes with the self-adjoint
operator A. Haar averaging the vertex action gives an orthogonal
projection P_G commuting with its semigroup. The physical Hilbert
space is the range of P_G, equivalently the invariant subspace; (4)
is its infinitesimal constraint on the smooth core. This avoids a
formal manipulation of unbounded Gauss operators outside their domains.

The finite-volume heat is positivity improving. Bounded real potential
perturbation preserves positivity improvement and compact resolvent.
Hence A has a unique strictly positive normalized vacuum h_Lambda.
Gauge transformations preserve positivity and normalization; uniqueness
therefore gives P_G h_Lambda=h_Lambda. Once a full-space gap is proved,
restriction to this invariant subspace cannot lower it.

## YM62-T2: the four-link block estimate

The combinatorial extension can be stated for bounded interactions
|W_p|<=1 involving at most q distinct link variables, with at most d
interactions containing any one link. A cubic three-dimensional lattice
has q=4 and d<=4. An interior link meets four plaquettes and has
d(q-1)<=12 other-link incidences. Incidences retain their multiplicity;
different faces containing the same link pair are not merged.

Put K_t=exp(-t L_E),
S_h=K_(h/2) exp(h theta V) K_(h/2), and
T_h=exp(h theta V/2) K_h exp(h theta V/2). Take 0<h<=1. YM53-T1
with C=I and beta=1 gives

\[
 19/20\le k_t(q,q')\le21/20\quad(t\ge9),\qquad
 \epsilon=4/5\le(19/21)^2.                                  \tag{5}
\]

This is an all-representation-content bound, from the convergent tail
sum_(n>=1)(n+1)^2 exp(-tn/2); it is not a spin cutoff. Let
r=ceil(9/h), ell=Jr, with J a positive integer, so 9<=hr<10.
The free bridge coupling for one changed temporal endpoint has
expected Hamming cost at most r/epsilon. For a conditional link block
with k<=ell interior time rows:

* Its log interaction tilt has magnitude at most d h |theta| k.
  Including normalization, its total-variation distance from the free
  bridge is at most 2d h |theta| k.
* Coupling each tilted law to its free law gives endpoint cost
  D(k)=r/epsilon+4d h |theta| k^2.
* Changing one other link at one time changes at most m shared face
  terms. Its conditional log-density difference is at most
  2m h |theta|; the normalized total-variation cost is at most
  4m h |theta| and Hamming cost at most 4m h |theta| k.

Here TV is half-L1. For example, a log tilt bounded by a has normalized
density ratio at least exp(-2a), hence TV<=1-exp(-2a)<=2a.
The changed-exterior estimate applies the same bound to the ratio of
the two unnormalized tilts. All bounds may exceed one without failing.
Admissible bridge supports and grouped endpoint changes are handled
exactly as in YM53-T2; no inadmissible intermediate bridge is inserted.

On N interior time rows use every clipped label starting at
2-ell,...,N in every link column. Retain duplicate clipped blocks.
Each coordinate is removed ell times, affects at most two temporal
endpoint labels and at most d(q-1)ell spatial labels with multiplicity.
Thus its total outgoing cost divided by ell is at most

\[
 \begin{split}
 \alpha&=\frac{2r}{\epsilon\ell}
       +[8d+4d(q-1)]h|\theta|\ell\\
 &\le\bar\alpha=\frac{2}{\epsilon J}
                  +40d(q+1)|\theta|J.
 \end{split}                                                  \tag{6}
\]

For the cubic plaquettes the fine coefficient is 80 and the coarse
coefficient is 800. Substitution d=2,q=2 recovers YM53's chain
coefficient 240, an independent normalization check.

Choose R>1 with R baralpha<1 and set gamma=log(R)/(10J).
Weights exp(gamma h min(u,N+1-u)) vary by at most R across a block
and its temporal endpoints. With M=|E|(N+ell-1), D=r/epsilon+16h|theta|ell^2,
coupled block replacement gives

\[
 W_{n+1}\le\left(1-\frac{\ell(1-R\bar\alpha)}M\right)W_n
           +\frac{2|E|\ell R D}{M}.                         \tag{7}
\]

This is the same marginal-preserving comparison as YM53-T3: charge
removal to each interior disagreement, charge the outgoing costs in
(6), and charge the two fixed time boundaries to the last term.
All spatial link variables are inside the strip. Free spatial faces
are simply absent. Iteration bounds the influence of temporal
boundary data on a middle-row observable by a finite prefactor times
exp(-gamma h times its boundary distance). The prefactor may depend
on volume and h; gamma does not.

For J=5,

\[
 \bar\alpha=\tfrac12+4000|\theta|,
 \quad R=\frac{1+\bar\alpha^{-1}}2,
 \quad \boxed{\gamma(\theta)=\frac1{50}
      \log\!\left(\frac{1+\bar\alpha^{-1}}2\right)>0}.       \tag{8}
\]

These choices work throughout (1). At theta=1/16000 they give
baralpha=3/4, R=7/6, margin 1/8, and gamma=log(7/6)/50>0.00308301.
The identity
J(baralpha-1)=(J-5)^2/10+(800|theta|-1/10)J^2 shows that at or above
1/8000 no integer J improves this sufficient budget. It does not show
that the true gap vanishes there.

## YM62-T3: full-space gap and continuous time

Let b=|theta||P|. The coarse row kernel satisfies

\[
 e^{-hrb}K_{hr}(x,y)\le T_h^r(x,y)\le e^{hrb}K_{hr}(x,y).     \tag{9}
\]

Together with (5), this gives a unique bounded strictly positive
vacuum h_h and top eigenvalue lambda_h. Conditional-strip comparison
from (7), followed by vacuum averaging at the endpoints, bounds
centered stationary row correlations. For every bounded f with
Phi(h_h^2 f)=0 and g=f h_h,

\[
 \langle g,(T_h/\lambda_h)^{2n}g\rangle
       \le C_g e^{-2n\gamma h}.                            \tag{10}
\]

By the spectral theorem, a spectral component above exp(-gamma h)
would contradict (10) as n tends to infinity. These g are dense in
the vacuum complement. Thus the complement ceiling is exp(-gamma h)
on the entire completed space. With X=K_(h/2) exp(h theta V/2),
T_h=X*X and S_h=XX* have the same positive nonzero spectrum; the
isolated vacuum and complement bound transfer to S_h. A zero spectral
component does not affect the estimate. There is no polynomial cutoff.

For time refinement a face is linear in each of its four link variables:

\[
 L_E W_p=3W_p,\qquad
 \Gamma(W_p)=\sum_{e,a}(D_{a,e}W_p)^2=1-W_p^2\le1.           \tag{11}
\]

The second identity holds on the unit carrier. Before imposing the
four unit-radius relations its right side is the product of the four
squared radii minus W_p^2, as checked exactly. For a convenient
conservative bound, Gamma(W_p)<=3 and K_s W_p=exp(-3s)W_p imply

\[
 K_s[(W_p-W_p(x))^2](x)\le2(1-e^{-3s})\le6s.
\]

Cauchy--Schwarz in the symmetric Markov kernel yields
||[K_s,M_(theta V)]||<=b sqrt(6s). The bounded-insertion packet proof
of YM44, used in YM53-T4, then gives for mesh delta and total time t

\[
 \left\|S_{h_n}\cdots S_{h_1}-e^{-tA_\Lambda}\right\|
       \le b t e^{bt}\sqrt{18\delta}.                      \tag{12}
\]

No Taylor expansion of the unbounded Casimir is used. This is a
fixed-volume norm limit; its prefactor is not uniform in volume.
The limiting semigroup has a unique positive vacuum by (9)'s
continuous-time Feynman--Kac counterpart. Its norm is at least exp(-bt).
Norm convergence transports the isolated vacuum projection and the
normalized complement ceiling: for h=t/n apply the discrete bound
to S_h^n, then let n tend to infinity. The isolation gate is
2 times the error in (12) < exp(-bt)(1-exp(-gamma t)), eventually true.

Writing E_Lambda=inf spec A, H_Lambda=A-E_Lambda and P_0=|h><h|,

\[
 \boxed{H_\Lambda\ge\gamma(I-P_0),\qquad
   \|(I-P_0)e^{-tH_\Lambda}(I-P_0)\|\le e^{-\gamma t}.}     \tag{13}
\]

T1 restricts (13) to the Gauss-invariant subspace. The inherited
ground-source cancellation gives, for smooth F,
<Fh,H Fh>=mu_Lambda(Gamma(F))>=gamma Var_(mu_Lambda)(F), where
d mu_Lambda=h^2 d Phi_E. The interacting weight is the actual vacuum,
not the free reference or a truncated perturbative exponential.

## YM62-T4: physical fluctuations remain nonzero

Use the periodic cubes, where |P|=|E|=3L^3. Translation and cubic
symmetry, and uniqueness of the vacuum, make <h,L_e h> the same on
every link. The constant trial vector has energy zero, so E_Lambda<=0
and

\[
 \langle h,L_E h\rangle=E_\Lambda+\theta\langle V\rangle
          \le |\theta||P|,\qquad \langle h,L_e h\rangle\le|\theta|.
                                                               \tag{14}
\]

For a set S of four face links let P_S project those links onto their
free constant vector, leaving all other links unrestricted. The free
one-link gap is 3/4; commuting product projections give

\[
 \varepsilon_S:=\|(1-P_S)h\|^2
       \le\frac43\sum_{e\in S}\langle h,L_e h\rangle
       \le\frac{16}3|\theta|<\frac1{1500}.                 \tag{15}
\]

The normalized vector v=P_S h/||P_S h|| factors as 1_S times a vector
on the remaining links. The two pure states of h and v have trace
distance 2 sqrt(epsilon_S): diagonalize their difference on their
two-dimensional span. Therefore for every bounded F supported in S,

\[
 |\mu_\Lambda(F)-\Phi_S(F)|\le2\|F\|_\infty\sqrt{\varepsilon_S}.
                                                               \tag{16}
\]

Since Phi(W_p)=0 and Phi(W_p^2)=1/4, (15)--(16) show

\[
 \operatorname{Var}_{\mu_\Lambda}(W_p)
 \ge\tfrac14-2\sqrt{\varepsilon_S}-4\varepsilon_S
 >\tfrac14-\tfrac1{19}-\tfrac1{375}>\boxed{19/100}.         \tag{17}
\]

Here sqrt(1/1500)<1/38. This is a bound on a genuine gauge-invariant
observable, uniform in volume throughout (1). It prevents the limit
below from becoming a one-dimensional vacuum space with no excitations.

## YM62-T5: nontrivial gapped infinite-volume histories

The ground transform h^(-1) exp(-tH) h is a stationary reversible
Markov semigroup on each finite Q^E with invariant mu_Lambda.
Its path laws are positive and reflection positive. For any real smooth
local F, the ground-source identity in T3 gives

\[
 \mathbb E[(F(t)-F(0))^2]\le2t\,\mu_\Lambda(\Gamma(F))
                       \le2t\|\Gamma(F)\|_\infty.           \tag{18}
\]

Only derivatives on F's finitely many links occur, so the bound is
independent of L. For a single face it is at most 2t by (11).

Embed a periodic cube configuration periodically in Z3, retaining its
three oriented link types. No convex averaging of distinct vacua is
used. Each embedded history law is exactly spatially stationary.
The compact countable product Q^(oriented Z3 links x rational times)
has a weakly convergent subsequence. Local finite-dimensional
distributions, stationarity, consistency, time reversal and reflection
positivity pass to this limit. Estimate (18) extends rational times
to real times in mean square, as in YM60-T2. A finite set of local
gauge transformations on any fixed observation window is represented
on sufficiently large tori. Gauge invariance therefore also passes.
The histories are invariant under time-independent local gauge
transformations, as appropriate to this Hamiltonian/temporal-gauge
description; arbitrary independent changes at each time are not asserted.

The gap requires more than positivity. Let G be any finite linear
combination of bounded continuous positive-time cylinder words. Its
finite-volume reflected vector v_(L,G) is obtained by the corresponding
ordered multiplications and heat semigroups applied to h_L. For the
reflected inner product and positive time translation tau_t,

\[
 \begin{split}
 \|G\|_{\mathrm{OS},L}^2&=\|v_{L,G}\|^2,\qquad
      \mathbb E_LG=\langle h_L,v_{L,G}\rangle,\\
 0\le\langle G,\tau_tG\rangle_{\mathrm{OS},L}-|\mathbb E_LG|^2
 &\le e^{-\gamma t}
       (\|G\|_{\mathrm{OS},L}^2-|\mathbb E_LG|^2).
 \end{split}                                                  \tag{19}
\]

This follows directly from (13) on v_(L,G). Its right-hand coefficient
is an actual reflected norm, not a possibly divergent mixing prefactor.
Every quantity in (19) is a finite cylinder expectation. Pass to the
single subsequence, then use mean-square continuity for real t.
The reflection-positive null quotient, completed over all positive-time
cylinders, carries the strongly continuous self-adjoint contraction
semigroup exp(-tH_infinity). Such vectors are dense by construction.
Thus (19), density and the spectral theorem imply

\[
 \boxed{H_\infty\ge\gamma(I-|\Omega\rangle\langle\Omega|).} \tag{20}
\]

The same argument on the algebra of local gauge-invariant readouts
gives the physical vacuum representation. At time zero, (17) gives
||(W_p-mu(W_p))Omega||^2>19/100, so it is nontrivial. The vacuum
kernel in each constructed representation consists only of Omega.
This does not prove that all subsequences give the same state, nor
identify the history reconstruction with a time-zero Markov generator.
It also makes no claim about other superselection sectors.

## YM62-T6: the response tower now starts with closed loops

Introduce independent real sources j_p in
A(j)=L_E-sum_p j_p W_p at fixed finite volume. The simple isolated free
vacuum has ordinary analytic bounded-perturbation response near j=0.
The convergence neighbourhood is not claimed uniform in volume.
Ordered multiplication and the product reference give

\[
 \Phi(W_p)=0,\quad \Phi(W_pW_q)=\tfrac14\delta_{pq},\quad
 L_E^{-1}W_p=\tfrac13W_p.                                   \tag{21}
\]

For p!=q there is a link occurring once in W_pW_q; integration over
its native reference kills the coefficient. Every face has four
distinct links and distinct faces have distinct supports in the
declared boxes and tori. The four-link product is itself Haar distributed
under the free reference, giving its second moment 1/4.

Consequently the first source equation and vacuum energy are

\[
 h^{(1)}=\tfrac13\sum_p j_pW_p,\qquad
 E(j)=-\tfrac1{12}\sum_p j_p^2+O_\Lambda(\|j\|^3),\qquad
 \left.\frac{\partial\langle W_p\rangle}{\partial j_q}\right|_0
       =\tfrac16\delta_{pq}.                               \tag{22}
\]

For uniform j_p=theta the quadratic coefficient per periodic spatial
cell is -1/4. This coefficient is volume independent; the error in
(22) is not asserted uniform. With P_p(j)=<W_p>=-partial_p E(j),
the exact finite-volume response obeys

\[
 \partial_q P_p=\partial_p P_q,\qquad
 \chi=-\nabla^2 E\succeq0,                                 \tag{23}
\]

where the Hessian is evaluated in a differentiable simple-vacuum
region. Mixed-derivative symmetry follows from the potential E;
positivity follows from ground-energy concavity. This is the precise
response-potential structure connecting the earlier Jacobian diagrams
to these gauge observables. Higher source derivatives must respect
closed-loop incidence and Gauss symmetry. A placement diagram alone
does not determine their values or the continuum dynamics.

As an independent established-theorem compatibility check, group the
three outgoing links at each vertex into one cell. Its free gap is
delta=3/4. For range B={0,e1,e2,e3}, set
h_x=(4/3)sum_(y in x+B)L_cell,y and
phi_x=-(16/3)theta sum_(mu<nu)W_(x,mu,nu).
Then sum h_x=(16/3)L_E, each h_x has unique product vacuum and gap
at least one, and ||phi_x||<=16|theta|. They are translation invariant
and diagonal in the product free spectral partitions; infinite-dimensional
cell Hilbert spaces are allowed. These verify the hypotheses of
[Yarotsky, arXiv:math-ph/0412040, Theorem 1](https://arxiv.org/abs/math-ph/0412040)
for sufficiently small coupling. That theorem supplies an independent
non-numeric perturbative window; it is not the source of the numerical
constant 1/8000 or an analytic-radius claim for all of (1).

## Verification and remaining work

The certificate checks exact ordered gauge covariance, the Gauss and
Casimir polynomial identities, a failing rank-one colour operator,
wrong face order and open-overlap counterexamples, plaquette moments,
finite lattice and clipped-block incidences, strict coupling gates,
outward time-refinement/rate bounds and the local-variance estimate.
It preserves every YM61 source pin and reuses the native arithmetic.
Finite controls support the written proofs; they are not a mechanical
formalization, a numerical computation of the interacting vacuum, or
independent expert review.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym62_gauge_lattice_gap.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym62_gauge_lattice_gap.py -v
```

The next physical gate is a spatial-spacing/coupling trajectory with
uniform control of the scaled gap and nontrivial continuum observables.
The weak-magnetic/electric-ratio bound above cannot be extrapolated
along an unproved weak-coupling continuum trajectory. Spatial rotational
symmetry restoration, the full continuum axioms, physical calibration
and extension to every compact simple group remain open. The
[official Clay statement](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf)
requires a nontrivial quantum theory on R4 with a positive gap; a fixed
spatial lattice with continuous time does not meet that requirement.
