# YM-44: full interacting time refinement with an explicit Smriti tail

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

**Result under the existing positive-functional heat-kernel adapter:** for
every fixed finite chain width, the full interacting time transfer converges
in operator norm as its step is refined. The result controls every source
and every content, rather than only a leading tiling or a compressed channel.
The new ingredient is a heat/bridge commutator estimate from two exact
character identities. The error is explicit.

This establishes the time-limit construction on the declared trajectory. It
does **not** establish a gap uniform in spatial width and time cutoff.
YM-43's three coarse cells do not automatically supply that stronger bound.

## Carrier, observer, source and prior results

Use the same declared single-site positive normalized functional, finite
product functional Phi, coefficient orthogonality/convolution dictionary,
null-seam quotient and positive-form completion as YM-9/12/19/43. Its SU(2)
reference integral remains an admitted adapter. It has not been identified
with a primitive-derived Phi_Sigma or with the NCG matter/gauge measure.

For m sites let K_t be the **full**, symmetric, positivity-preserving free
heat transfer, with K_t 1=1 and coefficient factors
Exp_Sigma(-t sum_i C_(j_i)), C_j=j(j+1). YM-9's exact coefficient semigroup
gives K_s K_t=K_(s+t). On the dense finite-content core, L multiplies each
coefficient by sum_i C_(j_i). Thus K_t has generator -L on that core.
The free family is strongly continuous on the completed source space by
finite-core approximation and contraction. L is unbounded.

Set

\[
v_i(x)=\tfrac12\chi_{1/2}(x_i x_{i+1}^{-1}),\qquad
V=\sum_{i=1}^{m-1}v_i,\qquad B=\theta M_V,\qquad
b=|\theta|(m-1).
\]

Here M means multiplication and theta is fixed. Since |v_i|<=1,
||B||<=b; equality holds in the compact-coordinate adapter. The chosen
trajectory is kappa(a)=theta a. The existing square-sourced time step is

\[
S_a=K_{a/2}\operatorname{Exp}_\Sigma(aB)K_{a/2}.
\]

The recognition target is the **entire** time-transfer operator on this
fixed carrier. An arbitrary bounded source/observer pair is therefore
controlled by the same operator error. The refinement index is a count;
t is a declared heat parameter, not a derived material clock. All products
retain their order, and the error is recorded as unresolved Smriti memory.

**Reuse rather than duplication.** YM-9 already treats the free semigroup
and a fixed-graph sandwich gap. YM-22 treats a leading tiling rate, and
YM-40 treats compressed channel rates. EMK-T2 already proves order-sensitive
transport in a nilpotent sector, explicitly excluding general refinement
limits. RKF T28 supplies the completion discipline; WC2 supplies bounded
ordered series and factorial tails. None is relabelled as a new theorem.
Here the unbounded full heat family receives a specific operator-norm
refinement estimate. No second shared operator engine is introduced.

The product-formula method has established lineage: H. F. Trotter,
[*On the product of semi-groups of operators*](https://doi.org/10.1090/S0002-9939-1959-0108732-6)
(1959). That theorem is attribution, not a premise replacing the proof below.
No priority or optimal convergence-rate claim is made.

## YM44-T1: a character identity controls the full commutator

For a single bridge v=v_i, the existing coefficient multiplication and
two-end heat action give

\[
K_s v=e^{-3s/2}v,\qquad
K_s(v^2)=\tfrac14+e^{-4s}(v^2-\tfrac14).
\]

**Proof.** The two endpoints each carry spin 1/2, hence total Casimir 3/2.
The exact fusion identity chi_(1/2)^2=chi_0+chi_1 gives
v^2=(1+chi_1)/4; the nonconstant term has total Casimir 4. This uses
the full free kernel, without truncating its action on other sources.

At a fixed x write z=v(x)^2 in [0,1]. The positive conditional square is

\[
d_s(z)=K_s[(v-v(x))^2](x)
=\frac{1-e^{-4s}}4+(1+e^{-4s}-2e^{-3s/2})z.
\]

It obeys 0<=d_s(z)<=3s for every s>=0. For the upper bound rewrite it as
(1-4z)(1-e^(-4s))/4+2z(1-e^(-3s/2)).
When z<=1/4, the inequality 1-e^(-r)<=r gives d_s(z)<=(1-z)s<=s.
When z>=1/4, the first term is nonpositive and d_s(z)<=3zs<=3s.
Positivity also follows purely algebraically: for u=e^(-s/2),
4 d_s(1)=(1-u)^2(5+10u+15u^2+12u^3+9u^4+6u^5+3u^6),
while 4 d_s(0)=1-u^8; d_s is affine in z.

For any source f, conditional cut-square Cauchy–Schwarz yields

\[
|[K_s,M_v]f(x)|^2
\leq d_s(v(x)^2)\,K_s(|f|^2)(x).
\]

Integrate with Phi and use its K_s invariance. Sum the bridge commutators
in their original order to obtain

\[
\boxed{\|[K_s,B]\|\leq b\sqrt{3s}.}
\]

This is an all-source operator bound. It does not assert that [L,B] is a
bounded operator. That distinction avoids an invalid unbounded-generator
Taylor expansion.

## YM44-T2: construct the full interacting continuation

Define ordered insertion packets

\[
D_k(t)=\int_{0<s_1<\cdots<s_k<t}
 K_{t-s_k}B K_{s_k-s_{k-1}}B\cdots B K_{s_1}\,ds_1\cdots ds_k,
\quad D_0(t)=K_t.
\]

The integrals mean limits of rational cut sums on each source. For a
fixed source their integrands are continuous by strong continuity of K
and boundedness of B. The bound

\[
\|D_k(t)\|\leq\frac{(bt)^k}{k!}
\]

follows from contraction of every K factor and the ordered simplex
volume t^k/k!, obtained by partitioning the cube into k! orderings.
Consequently

\[
U(t)=\sum_{k\geq0}D_k(t),\qquad \|U(t)\|\leq e^{bt}
\]

exists in operator norm. For N with bt/(N+2)<1, the omitted packet obeys

\[
\left\|\sum_{k>N}D_k(t)\right\|
\leq \frac{(bt)^{N+1}}{(N+1)!}\frac1{1-bt/(N+2)}.
\]

This is an explicit factorial Smriti tail, not a formal infinite series.

**Composition and dagger.** Split every ordered insertion list at a time
s. Absolute summability permits regrouping, and the free semigroup combines
the adjacent heat factors. Hence U(t+s)=U(t)U(s). Reversing the simplex
coordinates and the dagger of each factor gives U(t)^dagger=U(t).
U(t)=U(t/2)^dagger U(t/2) is therefore cut-square positive.

Separating the last insertion gives the exact integral identity

\[
U(t)=K_t+\int_0^t K_{t-s}B U(s)\,ds.
\]

It proves strong continuity and, on the finite-content core,

\[
\lim_{t\downarrow0}\frac{U(t)f-f}{t}=(-L+B)f.
\]

Thus the derived energy-form expression on that core is H=L-theta M_V.
The sign follows from the actual interaction weight; it is not selected
by an analogy to a classical Hamiltonian. A claim of norm differentiability
of U at zero, or of a new general unbounded-domain theorem, is not made.

In the compact-coordinate adapter the same insertion packets define a
kernel for t>0. Replacing B by B+bI makes every insertion nonnegative.
The simplex bounds then give, pointwise,

\[
e^{-bt}K_t(x,y)\leq U(t;x,y)\leq e^{bt}K_t(x,y).
\]

This also proves positivity preservation. At any t where the existing
heat-kernel coefficient estimate gives a strictly positive lower floor,
YM43-T4 constructs U(t)'s positive vacuum by the same geometric iteration.
The construction requires no new measure or spectral-existence premise.

## YM44-T3: an explicit norm error for arbitrary time partitions

For every h>0,

\[
\|U(h)-K_h e^{hB}\|\leq b h\sqrt{3h}\,e^{bh}.
\]

**Proof.** In a k-insertion packet, move each B past the heat factors to
its right, preserving all other factors. An induction comparing the packet
with K_h B^k bounds the difference by k b^k sqrt(3h): at the inductive
step the new commutator is [B,K_s] for one s<=h. Integrating and summing
uses sum_(k>=1) k (bh)^k/k! = bh e^(bh). No divergent heat-generator
power series has been used.

The bounded exponential series separately gives

\[
\|[e^{hB},K_{h/2}]\|
\leq h e^{bh}\|[B,K_{h/2}]\|.
\]

Together,

\[
\|S_h-U(h)\|
\leq(\sqrt3+\sqrt{3/2})b h^{3/2}e^{bh}
<3b h^{3/2}e^{bh}
\]

when b>0; for b=0 the error is exactly zero.

Let a positive partition have lengths h_1,...,h_n, total t and mesh
delta=max h_j. The later action multiplies on the left. A finite telescoping
of the two ordered products, with norms at most e^(b h_j), gives

\[
\boxed{
\|S_{h_n}\cdots S_{h_1}-U(t)\|
\leq 3b e^{bt}\sum_j h_j^{3/2}
\leq 3bt e^{bt}\sqrt{\delta}.
}
\]

For equal steps a=t/n this is 3bt e^(bt) sqrt(t/n). It proves convergence
in the **full operator norm**, independent of the chosen time partitions.
Both recognized and complementary projections inherit the same bound;
T28's common-carrier and tail obligations are explicit. Source-observer
readout error is at most this bound times ||source|| ||observer||.

For fixed m the result is uniform for t in a bounded parameter interval.
Its constants depend on m through b. No uniform-in-volume limit, physical
continuum or optimal sqrt(delta) rate is inferred from it.

## YM44-T4: the order residue survives a coarse subdivision

EMK-T2 already establishes that order is content. Here, for bounded
finite-cut generators A and B, direct ordered coefficients give

\[
e^{hA/2}e^{hB}e^{hA/2}-e^{h(A+B)}
=h^3\mathcal D+R_4(h),
\]
\[
\mathcal D=-\tfrac1{24}[A,[A,B]]-\tfrac1{12}[B,[A,B]].
\]

The coefficients of degrees 0,1,2 agree. At degree 3 the respective
coefficients of AAB, ABA, BAA, ABB, BAB, BBA in the difference are
-1/24, 1/12, -1/24, 1/12, -1/6, 1/12. Expanding the two nested
commutators proves the identity word by word.

If M>=||A||+||B||, the remainder has bound
2 sum_(k>=4) (|h|M)^k/k!, controlled by WC2's factorial tail.
Consequently one coarse symmetric step and two half steps differ first
by (3/4)h^3 D. A nonzero D obstructs exact subdivision.

For A=-L and B=theta M_V this identifies the algebraic order defect on
the finite-content core. The bounded-cut remainder is **not** promoted to
an all-content O(h^3) operator bound: L is unbounded. T3 supplies the lawful
global bound instead. The defect concerns these two specified transports;
it is not an identification with NCG gauge curvature or Riemann curvature.

The independent two-state positive fixture has
L=(3/4)[[1,-1],[-1,1]], B=theta diag(1,-1).
Its certificate separates S_1 from S_(1/2)^2 at theta=1/16, verifies the
ordered cubic coefficient, and compares refined products with an independent
closed-form exponential of -L+B. It is a control, not an SU(2) model substitute.

## YM44-T5: a gap can be transported only with its normalization

Let positive self-adjoint A have unit top vacuum h, value lambda>0 and
complementary operator bound q lambda, with 0<=q<1. Let positive
self-adjoint C have a unit top vacuum g and ||C-A||<=epsilon.
Assume a certified lower bound ell<=lambda and

\[
2\epsilon<\ell(1-q).
\]

Then, writing nu=||C||,

\[
\nu\geq\lambda-\epsilon,\qquad
\|Q_g C Q_g\|\leq q\lambda+\epsilon,\qquad
\boxed{\frac{\|Q_g C Q_g\|}{\nu}
\leq\frac{q\ell+\epsilon}{\ell-\epsilon}<1.}
\]

**Proof without importing a min-max theorem.** Testing C on h gives the
vacuum lower bound. Every vector orthogonal to h has C-quadratic value
at most c=q lambda+epsilon times its norm squared. If some v orthogonal
to g had Rayleigh quotient greater than c, then every nonzero vector in
span{g,v} would also have quotient greater than c: cross terms vanish,
and nu>c by the strict margin. That two-dimensional span contains a
nonzero vector orthogonal to h, a contradiction. Thus the complementary
quadratic ceiling holds. For a positive operator this bounds its norm:
Cauchy–Schwarz in its positive form gives |<u,Cv>|<=c||u||||v|| on that
complement. The ratio decreases as lambda increases; replace lambda
by ell to obtain the displayed certificate.

As a consequence, if a later proof supplies, for **every sufficiently
small** a=t/n at each m, the normalized bound

\[
\|Q_a(S_a/\lambda_a)Q_a\|\leq e^{-a\gamma}
\]

with gamma>0 independent of m and a, T3 transports it to U(t), at any
t with its positive vacuum constructed. Indeed use A=S_a^n,
q=e^(-t gamma), epsilon=T3's error, and
lambda=lambda_a^n=||S_a^n||. Norm convergence gives
lambda_a^n -> ||U(t)||>0, so the strict margin holds for sufficiently
fine a at each fixed m. Let epsilon go to zero in the displayed bound.

The resulting rate would be uniform in m even though the required
refinement depth can depend on m. This is a **conditional implication**,
not the missing small-a gap estimate itself. A constant per-step ratio
and a rate per unit-a must not be interchanged.

## YM44-T6: exact boundaries and refusal witnesses

1. **Interaction is not a free subdivision.** The positive two-state
   witness in T4 has S_1 != S_(1/2)^2. A coarse YM-43 bound therefore
   cannot simply be raised to a fractional power for finer steps.
2. **No hidden content cutoff.** For any a>0, the free coefficient
   sup_j |e^(-a C_j)-1| equals 1. Thus ||K_a-I||=1 on the full carrier,
   even though K_a tends strongly to I. Finite-content norm-Taylor
   bounds cannot be used uniformly over all contents.
3. **Volume cost is real in this estimate.** At aligned bridge arguments
   |V|=m-1, so the multiplication norm is b=|theta|(m-1). Replacing it
   by a width-independent |theta| is false.
4. **Scaling is a premise.** Keeping kappa nonzero and fixed while
   a tends to zero gives S_a -> e^(kappa M_V) strongly, rather than I.
   It cannot give the finite generator in T2. The declared linear
   trajectory is not an asymptotic-freedom derivation.
5. **Nonzero order curvature is not a gap certificate.** Tensor an
   interacting factor with a spectator transfer having eigenvalues
   1 and 1-epsilon. The active order defect remains nonzero, but the
   total normalized complementary ratio is at least 1-epsilon.
   Let epsilon decrease. The certificate verifies the exact tensor
   action, rather than treating noncommutation alone as sufficient.
6. **Normalization and strict margin matter.** Omitting epsilon from
   the vacuum denominator, allowing equality in the gap gate, omitting
   the cubic word residue, or discarding the factorial tail are refused.

## Status and reproduction

| Obligation | Status |
| --- | --- |
| Full interacting time-refinement limit, every fixed finite m, all sources/content | Proved under the declared heat/functional adapter, with explicit norm error |
| Generator on the finite-content core | -L+theta M_V derived |
| Arbitrary positive time partitions | Same limit and mesh-dependent bound proved |
| Coarse step equals its interacting subdivisions | False in general; exact coefficient and positive numerical witness |
| Gap transport when a uniform fine-step bound is supplied | Proved conditional gate |
| Uniform interacting gap as both m grows and a decreases | Open |
| Physical trajectory, native measure/NCG dictionary, 4D continuum, Clay/QG | Open |

The certificate checks the fusion and moment identities, all cubic words,
directed tails, partition budgets, independent two-state refinement
enclosures, cut-square gap controls and the listed refusal witnesses.
It binds sources, proof, code, tests and workflow. The infinite proof is
written, not mechanically formalized; finite fixtures are not its replacement.
Previous YM certificates and the canonical RKF engine are unchanged.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym44_time_refinement.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym44_time_refinement.py -v
~~~

--check is read-only and is the default. --write explicitly regenerates
this new certificate. The dedicated CI schedules Python 3.12 only.
