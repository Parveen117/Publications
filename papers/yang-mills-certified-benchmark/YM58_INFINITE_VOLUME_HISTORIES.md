# YM-58: exact local histories on the infinite correlated-turn chain

Monty Dabas. 4 October 2026. **Python 3.12 only.**

YM57 established a width-independent gap for a specified correlated
kinetic chain. This chapter constructs its infinite-volume heat and
stationary local histories. The key is exact support preservation:
native derivatives and their constant-coefficient squares do not make
a readout depend on another site's coordinates. A bond crossing the
readout's boundary still supplies noise at its inner endpoint.

Consequently, finite-volume heat histories stabilize **exactly** once
the chain contains the readout sites and one layer of neighbours. The
infinite-volume centered generator retains
\[
 \boxed{\operatorname{gap}(H_\infty)\ge\gamma_k
       :=\min\{1/630,k/420\}>0,\qquad k>0.}                 \tag{1}
\]
This is a theorem for YM57's declared model and heat clock. The finite
global-label turn approximation has a separate refinement issue: a
sufficient joint volume/time condition is \(W_N/n\to0\), where
\(W_N=N+k(N-1)\). Without a relative-rate condition, that particular
discrete approximation can freeze even though its heat limit has (1).

## Sources and scope

Use YM50's compact native quaternion carrier Q, invariant positive
reference Phi_Q, and polynomial/recognition completions; YM51's native
derivatives; YM52's exact square-energy spectrum; and YM57's correlated
turn source and uniform gap proof. Write A_i=D_{1,i}, Z_i=D_{2,i}.
Distinct sites commute, while the inherited native brackets at one
site remain noncommutative. No new one-site multiplication, spin engine
or curvature definition is introduced.

At finite open interval Lambda in the integer lattice, the operator is
\[
 H_\Lambda=-\tfrac12\sum_{i\in\Lambda}A_i^2
            -\tfrac{k}{2}\sum_{i,i+1\in\Lambda}(Z_i-Z_{i+1})^2.
                                                               \tag{2}
\]
Its normalized vacuum is the constant 1 and its reference is the product
Phi_Lambda. A local or cylinder readout depends on finitely many sites.
The construction concerns continuous local readouts, their uniform
completion, and their L2 completion in the product reference. The
history algebra uses finitely many heat times, then uniform completion.
It does not claim every measurable bounded history lies in that algebra.

These are written mathematical proofs with exact finite controls, not
proof-assistant formalization or independent external expert review.
Consistent Markov marginals and symmetric-semigroup reflection positivity
are established methods; their application and boundary identities here
are specific to this selected native model.

## YM58-T1: the correct marginal includes boundary noise

Let S be any finite subset of Z, not necessarily an interval. Let E(S)
be the nearest-neighbour edges with both ends in S and partial E(S) the
edges with exactly one end in S. On functions of S define
\[
 \begin{split}
 H_S^m={}&-\tfrac12\sum_{i\in S}A_i^2
       -\tfrac{k}{2}\sum_{\{i,j\}\in E(S)}(Z_i-Z_j)^2\\
       &-\tfrac{k}{2}\sum_{\substack{\{i,j\}\in\partial E(S)\\i\in S}}Z_i^2
       =-\tfrac12\sum_{i\in S}A_i^2-k\sum_{i\in S}Z_i^2
                           +k\sum_{\{i,j\}\in E(S)}Z_iZ_j .       \tag{3}
 \end{split}
\]
Each boundary edge is counted once, at its unique endpoint in S. Each
lattice vertex has two incident edges, which explains the coefficient
-k of Z_i squared in the expanded expression. The empty marginal is
the zero generator on constants.

To prove (3), apply (2) to a cylinder polynomial F supported on S. Every
derivative at an absent site kills F. On an edge meeting S once, its
mixed derivative and outer-end square vanish, leaving the inner-end
square. Edges disjoint from S vanish entirely. The remaining derivatives
preserve S and preserve the degree at every site. Therefore, if Lambda
contains S and all its nearest neighbours,
\[
 H_\Lambda F=H_S^m F.                                      \tag{4}
\]
For any finite S subset T the same calculation gives
\(H_T^m(F\otimes1)=H_S^mF\otimes1\). This is exact projective
consistency, including disconnected sets.

One must not replace H_S^m by the open-chain operator H_S. For a singleton,
\[
 H_{\{i\}}^m=-\tfrac12 A_i^2-kZ_i^2.                       \tag{5}
\]
It contains two boundary bonds. YM56's centered invariant w_i obeys
\(A_iw_i=0,\ Z_i^2w_i=-w_i\), hence \(H_{\{i\}}^m w_i=kw_i\).
Deleting both boundary bonds instead makes w_i stationary; retaining
only one gives rate k/2. These are different dynamics, not sign conventions.

## YM58-T2: exact stabilization and a consistent infinite-volume heat

Fix a finite S and a bound on the polynomial degree at each site.
This is a finite-dimensional invariant space for all generators in
(4). Iterating (4) gives equality of every power; convergent matrix
exponential series then give
\[
 e^{-tH_\Lambda}F=e^{-tH_S^m}F\qquad(t\ge0)                 \tag{6}
\]
for every Lambda containing the one-neighbour layer of S. This is
equality at every time, not a truncation estimate.

Choose any such finite Lambda. YM57-T2 makes its heat positive, unital,
symmetric and product-reference-preserving. Its restriction in (6)
has the same properties. Polynomial density and the supremum-norm
contraction extend (6) to every continuous readout on Q^S. Marginal
consistency follows first on polynomials and then by this same density.

Define Phi_infinity on cylinders by finite product integration and set
\[
 P_tF=e^{-tH_S^m}F\quad\hbox{for }F\hbox{ supported on }S.   \tag{7}
\]
Consistency makes this independent of the chosen support. It extends
to the uniform completion of continuous cylinders, with
\(P_t1=1\), positivity, \(P_tP_s=P_{t+s}\), invariant Phi_infinity,
and symmetry in its positive pairing. On each finite coefficient space
it is strongly continuous. Density and contraction extend strong
continuity to both the uniform and L2 completions.

The L2 semigroup has a nonnegative self-adjoint generator H_infinity.
Equivalently, take the nested invariant finite polynomial spaces,
quotient reference-null polynomials, and close their compatible
nonnegative symmetric blocks. Symmetry makes each invariant space
reducing; their union is dense and is a generator core. This also
specifies the closed form without an infinite formal matrix sum.

There is no propagation estimate hidden here: this particular kinetic
source preserves cylinder support exactly. Adding coordinate-dependent
drift or a potential/ground transform would require a new argument.

## YM58-T3: positive stationary Markov histories and time reflection

For ordered times t_0<=...<=t_r and continuous cylinder readouts F_j,
define the stationary history functional by
\[
 \Omega(F_0(t_0)\cdots F_r(t_r))=
 \Phi_\infty\!\left(F_0P_{t_1-t_0}
   \left(F_1P_{t_2-t_1}\left(\cdots P_{t_r-t_{r-1}}F_r\right)\right)\right).
                                                               \tag{8}
\]
Use linear extension, complex conjugation for the involution, and merge
insertions at equal times. Each nested multiplication or propagation
is supported on the union S of the insertion supports. Equation (6)
therefore proves that the finite-volume vacuum value of (8) is already
equal to its infinite-volume value once Lambda contains S's neighbour
layer. The limit is independent of how its two ends recede.

For completeness, finite-volume positivity does not follow merely
from the positivity of each displayed monomial value. On a compact
finite carrier, each positive unital transition evaluation is a
probability measure. Iterating these evaluations against Phi_Lambda
gives a positive joint functional on arbitrary finite-time functions.
Thus Omega(G*G)>=0 for all finite linear combinations of histories.
Deleting an identity insertion uses the semigroup law; integrating the
first or last time uses stationarity or unitality. These verify time
consistency. Positivity and normalization imply a uniform norm bound,
so Omega extends to the uniform history completion.

Symmetry of P_t gives time reversal. Let G be a finite linear
combination of histories at nonnegative times, and let V_G be its
conditional future value at time zero, computed with the same nested
formula including the initial heat from zero. Reflection reverses time
and conjugates coefficients. Applying symmetry and (8) across zero gives
\[
                    \Omega((\Theta G)G)=\Phi_\infty(|V_G|^2)\ge0.  \tag{9}
\]
A common finite support suffices for the whole linear combination;
completion preserves (9). This is reflection positivity in the chosen
heat time, not a four-dimensional Euclidean gauge-field assertion.

This model has exact row-Markov composition by (7), where a row is the
full native configuration. No further coarse-graining is asserted to
be Markov. It does not close
the distinct row-reconstruction obligation for the earlier scalar
potential/history models merely because they share a native carrier.

The histories are genuinely correlated across a neighbouring pair.
For k=1, let f=q_{0,0}q_{0,1} and g=q_{2,0}q_{2,1}. Native derivatives give
\[
 H_{\{0,1\}}^m(f+g)=f+g,\qquad
 H_{\{0,1\}}^m(f-g)=\tfrac12(f-g),\quad
 \Phi(f^2)=\Phi(g^2)=1/16,\quad\Phi(fg)=0.
\]
Consequently \(\Phi(fP_tf)=(e^{-t}+e^{-t/2})/32\). For sites 0 and 2,
the same product readout has correlation \(e^{-3t/4}/16\), since no
edge joins the supports. At t=4 log 2 these are respectively 5/512
and 1/128. More generally the marginal generator splits across connected
components of S; their entire histories factorize in the product
reference. This finite spatial dependence explains why the exact
construction is much simpler than a general interacting field model.

## YM58-T4: the volume-independent gap survives completion

Embed a centered polynomial F(S) into any finite interval containing
its neighbour layer. Product integration leaves its variance and norm
unchanged, while (4) leaves its energy unchanged. YM57-T4 then gives
\[
 \langle F,H_\infty F\rangle
      =\langle F,H_\Lambda F\rangle
      \ge\gamma_k\,\bigl(\Phi_\infty(|F|^2)-|\Phi_\infty(F)|^2\bigr).
                                                               \tag{10}
\]
Closure extends this inequality to the full form domain. Alternatively,
the finite-volume centered contraction transfers by (6) to cylinders,
and then directly by L2 density. Thus
\[
 \|P_t(F-\Phi_\infty F)\|_2
       \le e^{-\gamma_k t}\|F-\Phi_\infty F\|_2.           \tag{11}
\]
This proves (1) and that constants are the unique stationary L2 sector.
It does not assert uniqueness among every possible singular invariant
probability measure. For continuous local F,G, Cauchy--Schwarz gives
\[
 |\Phi_\infty(\overline F P_tG)
             -\overline{\Phi_\infty F}\Phi_\infty G|
 \le e^{-\gamma_k t}
      \sqrt{\operatorname{Var}_{\Phi_\infty}(F)
             \operatorname{Var}_{\Phi_\infty}(G)}.         \tag{12}
\]
The rational lower bound is conservative. The singleton marginal has
the sharper exact YM52 gap
\(\min\{(1+2k)/8,\min(1/2,k)\}\); a one-site rate is not asserted
to be the exact full-chain gap. Its centered w_i supplies the valid
full-chain upper witness gap(H_infinity)<=k.

## YM58-T5: a sufficient joint limit for the actual global-label turns

Exact volume stabilization applies to the heat semigroup. It does not
apply at finite step to YM57's globally normalized label average. In a
chain of N sites, let W_N=N+k(N-1), and let S_epsilon be the average of
the signed turns exp(+/-epsilon sqrt(W_N) X), with label weights 1/W_N
for A_i and k/W_N for each bond. Here epsilon=sqrt(t/n).

For a fixed support S, write
\[
 A_S=|S|+k\,|\{\hbox{lattice edges meeting }S\}|.
\]
Take a coefficient space with total polynomial degree at most d on S.
Every active label X has norm at most d/2 there. Indeed each linear
coordinate carries a native axis generator of norm 1/2, and on a
degree-d product it acts as the sum of the coordinate generators;
equivalently its unitary weights have magnitude at most d/2. The
invariant native pairing gives the same skew-generator norm bound.
For a bond use the degrees at its two ends; their sum is at most d.
Inactive labels act as identity. Hence
\(\|H_S^m\|\le A_Sd^2/8\).

For a skew matrix X, the paired exponential remainder after its
quadratic term has norm at most u^4||X||^4/24, by diagonalization and
the real cosine remainder. The heat remainder after I-tH/n has norm
at most t^2||H||^2/(2n^2). Their second-order terms agree. Both the
signed-turn average and the heat step are contractions, so telescoping
n steps proves, in the native L2 norm on this coefficient space,
\[
 \boxed{\left\|S_{\sqrt{t/n}}^n-P_t\right\|
 \le {t^2d^4\over n}
       \left({W_NA_S\over384}+{A_S^2\over128}\right).}     \tag{13}
\]
Assume throughout that the chain contains S's neighbour layer. Formula
(13) proves joint convergence for fixed t whenever N,n grow with
W_N/n tending to zero. The term involving A_S squared is volume-free.
The bound is sufficient, not a claimed necessary or optimal condition.

For a fixed polynomial history, multiplication by its finitely many
insertions stays within a finite coefficient space whose degree can
be bounded by the sum of their degrees. Apply (13) to each time gap
and telescope the nested products; bounded insertion multipliers give
the corresponding finite sum of errors. Finite-dimensional norm
equivalence also gives uniform convergence on this fixed support.
Uniform polynomial approximation and positive unital contraction extend
the result to continuous local histories. Each positive time gap may
use its own n_j, with W_N/n_j tending to zero. No unrestricted joint
limit of the global-label approximation is inferred from (6).

## YM58-T6: an exact obstruction to unrestricted joint refinement

At fixed step count n, there is an even simpler volume-first obstruction.
A global label meets S with probability A_S/W_N. Inactive labels fix
F(S), and every active turn is a supremum-norm isometry preserving S.
Contraction and telescoping give
\[
 \|S_{\sqrt{t/n}}^nF-F\|_\infty
                 \le {2nA_S\over W_N}\|F\|_\infty .       \tag{14a}
\]
Thus taking volume first freezes every local continuous readout. Taking
refinement first gives the nontrivial heat (7); the two iterated limits
of this global-label scheme disagree. A local update schedule would be
a different approximation and is not needed for the heat construction.

Choose one interior site's linear coordinate q_{0,i}, k=1 and t=1.
Every signed on-site or incident bond turn multiplies its average by
cos(sqrt(W_N/n)/2); all other labels leave it fixed. The active total
weight before division by W_N is a=3. Its n-step multiplier is therefore
\[
 r_{N,n}=\left[1-{3\over W_N}
              \left(1-\cos\!\left({\sqrt{W_N/n}\over2}\right)\right)
          \right]^n.                                    \tag{14}
\]
The correct heat multiplier is e^{-3/8}. If N=n^2, then W_N=2n^2-1.
For n>=2 the bracket is in [0,1], and Bernoulli's inequality gives
\[
             1-{6n\over2n^2-1}\le r_{n^2,n}\le1.
\]
Thus this volume/time sequence converges to the identity on that
nonconstant coordinate even though epsilon=1/sqrt(n) tends to zero.
The actual rotation amplitude epsilon sqrt(W_N) does not tend to zero.
At n=100 the multiplier is at least 19399/19999>96/100, whereas
e^{-3/8}<7/10: the latter follows already from the first three terms
of e^{3/8}>10/7. This is a strict rational separation, without floating
trigonometry. The boundary-noise omission and the fast-volume sequence
are both required negative controls.

The new result closes the infinite-volume local-history and L2-gap
step for the free correlated kinetic model. It also supplies a checked
joint-refinement gate for its specified records. Adding the old scalar
interaction, selecting a physical gauge action and clock, imposing the
relevant gauge constraints, removing spatial lattice spacing and
identifying four-dimensional Yang--Mills remain separate obligations.

## Reproduction and evidence

`certificates/ym58_infinite_volume_histories.py` reuses the existing
native polynomial and product-moment engine. Its exact rational checks
cover disconnected marginal consistency, open-volume stabilization,
boundary noise, powers and multi-time Taylor coefficients, symmetric
energy matrices, neighboring versus separated path correlations,
reflection Gram matrices, and the failed unrestricted global-label limit.
Finite series coefficients are labeled as such; they are not passed off
as positive heat operators or as the proof for all times.

The source manifest preserves the preceding 139 upstream files and pins
this chapter, implementation, tests and dedicated workflow. The written
proof supplies the infinite-volume and all-degree statements.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym58_infinite_volume_histories.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym58_infinite_volume_histories.py -v
```

Default and `--check` are read-only; `--write` regenerates only YM58's
own result and hash. Optimized Python and non-3.12 runtimes are refused.
