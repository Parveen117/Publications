# YM-45: a uniform fine-time gap from temporal blocks

Monty Dabas. Runtime: **Python 3.12 only**.

**Result within the declared full heat-kernel chain:** for every
\[
|\theta|<1/1680
\]
there is an explicit rate gamma>0, independent of finite spatial width m
and of 0<a<=1, such that the actual square-sourced transfer obeys
\[
\boxed{\|Q_a(S_a/\lambda_a)Q_a\|\leq
       \operatorname{Exp}_\Sigma(-\gamma a).}
\]
The bound covers every complementary source and all coefficient contents.
It also passes to YM-44's interacting time limit at every fixed finite m:
\[
\boxed{\|Q\,U(t)/\|U(t)\|\,Q\|\leq
       \operatorname{Exp}_\Sigma(-\gamma t),\qquad t>0.}
\]

The new step is to compare whole temporal blocks whose heat-parameter
length stays bounded as a decreases. The old single-slice overlap need
not stay useful. The interaction range above is a sufficient, conservative
window, not a phase boundary. In particular, the old theta=1/16 example
is outside this certificate.

## Carrier, premise order and scope

Use exactly YM-9/12/19/43/44's full positive symmetric SU(2) heat kernel,
normalized reference functional Phi and its finite products. The coefficient
semigroup, positive conditional normalization, product integration and
bounded measurable functions are the **admitted compact-coordinate adapter**.
The positive cut form, null-seam quotient, completion and cut-square
Cauchy–Schwarz are as in YM-43. No primitive-derived native measure is
asserted. The identification with Phi_Sigma and with the NCG action/measure
remains an explicit dictionary obligation.

For an open spatial chain of m sites set
\[
v_i(x)=\tfrac12\chi_{1/2}(x_ix_{i+1}^{-1}),\quad
V=\sum_{i=1}^{m-1}v_i,\quad B=\theta M_V,\quad b=|\theta|(m-1).
\]
Here |v_i|<=1 and the declared trajectory is kappa(a)=theta a. Write
\[
S_a=K_{a/2}e^{aB}K_{a/2},\qquad
T_a=e^{aB/2}K_a e^{aB/2}.
\]
K on a whole row is the product of the one-site kernels. T is the
convenient row ordering; T6 below transfers its bound to the existing S.
The observer is any bounded whole-row function, and ultimately any source
in the completed positive form. Conditional replacement is a proof device
for boundary memory, not a newly postulated physical event law.

All exponentials/logarithms denote the existing completed native cut-scalars.
The parameter a is a declared heat parameter, not a derived physical clock.
The model is the existing bounded-overlap chain, with at most two spatial
neighbours per site; it is not the full four-dimensional gauge lattice.
No infinite-volume state, AF trajectory, physical continuum or Clay theorem
is inferred merely from a rate uniform in finite widths.

**Reuse and lineage.** YM-19 supplies the full coefficient tail, YM-43 the
common-part coupling, vacuum iteration, even-moment extraction and ordering
intertwiner, and YM-44 the norm refinement and normalized gap-transport gate.
YM-33/35's proposed worldline/tiling expansion is not assumed convergent.
No second operator engine is introduced. The established block-comparison
method is credited to the Dobrushin–Shlosman lineage; see the authors' paper
[Rebeschini–van Handel, *Comparison Theorems for Gibbs Measures*](https://web.math.princeton.edu/~rvan/dobrushin130819.pdf).
The specific estimates are proved below rather than imported as a comparison
theorem. This is a scoped derivation and integration, not a priority claim.

## YM45-T1: coarse overlap controls a fine conditioned bridge

YM-19's unchanged all-content tail gives
\[
\mathcal S_s=\sum_{n\geq1}(n+1)^2e^{-s n(n+2)/4},\qquad
1-\mathcal S_s\leq K_s(x,y)\leq1+\mathcal S_s.
\]
The outward certificate verifies mathcal S_6<1/20; its summands decrease
in s. Consequently, for every s>=6,
\[
19/20\leq K_s(x,y)\leq21/20,\qquad
\epsilon:=4/5\leq(19/21)^2.
\]
For any common nonnegative future weight H with positive normalizer,
the normalized laws proportional to K_s(x,y)H(y) and K_s(x',y)H(y)
have a common part of mass at least epsilon. Both the kernel ratio and
the normalizer ratio are included. This is YM43-T1 applied before discarding
the future weight; H itself need not have a uniform lower bound.

Fix 0<a<=1 and r=ceil(6/a), so 6<=ar<7. Consider a one-column free
bridge with k interior fine-time sites, different left endpoints and a
common right endpoint z. Its skeleton at times r,2r,... uses the conditional
transition
\[
p_x(dy)\ \propto\
K_{ar}(x,y)K_{a(k+1-u-r)}(y,z)\,\Phi(dy)
\]
from skeleton time u, whenever u+r<=k. Couple each such transition
using its common part. While the two skeleton states differ, the next
pair agrees with probability at least epsilon. After they agree, use
identical future transitions.

Fill the fine path between successive skeleton states with the appropriate
free conditional bridge. If both endpoints agree, use an identical filling;
otherwise the segment contributes at most r mismatches. The final shortened
segment contributes at most r, conditional on a remaining mismatch.
Thus, including k<r,
\[
\boxed{\mathbb E\,d_H(X,\widetilde X)
 \leq r\sum_{j\geq0}(1-\epsilon)^j=r/\epsilon.}
\]
Here d_H counts differing interior coordinates. The unspent bound after
J skeleton trials is r(1-epsilon)^J/epsilon, an explicit geometric memory
tail. Reversing time proves the same bound for a changed right endpoint.
All bridges are finite normalized laws from the declared positive kernels.
No unconditional skeleton law is substituted for a conditioned one.

## YM45-T2: two different perturbation costs for an interacting block

Update k consecutive time sites in one spatial column, with all other
coordinates fixed. Relative to the free time bridge, the conditional tilt
is e^H with
\[
H=a\theta\sum_{u\ {\rm in\ block}}\sum_{j\ {\rm spatial\ neighbour}}
 v(x_u,y_{j,u}),\qquad |H|\leq2a|\theta|k.
\]
Its normalized density relative to the free law is at least
e^(-4a|theta|k). The half-L1 distance is at most
1-e^(-4a|theta|k)<=4a|theta|k. A common-part coupling therefore costs at
most 4a|theta|k^2 in expected Hamming distance.

For a change of **one temporal endpoint**, couple tilted to free, couple
the two free bridges by T1, and couple free to tilted. Gluing over the
shared intermediate marginals and using the Hamming triangle inequality gives
\[
D(k):=r/\epsilon+8a|\theta|k^2.
\]
For clarity, gluing means composing the conditional kernels of two joint
laws over their common marginal. It preserves each outer marginal exactly.
On finite cuts this is summation with the positive intermediate normalizer;
on the compact adapter the same lawful product-kernel integration applies.

For a change of **one external spatial coordinate at one time**, the
unnormalized conditional density ratio is in
[e^(-2a|theta|),e^(2a|theta|)]. Including its changed normalizer, the
half-L1 distance is at most 4a|theta|. Hence the block mismatch cost is
at most 4a|theta|k.

For several exterior changes, insert them one coordinate at a time and
glue the couplings. Costs add by the pointwise Hamming triangle inequality.
No independence of exterior differences is needed. In particular, the
temporal tilt cost must not be replaced by the smaller spatial cost.

## YM45-T3: clipped blocks give an exact removal and influence budget

Take a finite strip with fixed time rows 0 and L and interior rows
1,...,N, where N=L-1>=1. Choose an integer J>=1 and ell=Jr.
For every spatial column i and each start s=2-ell,...,N, use the label
\[
B_{i,s}=\{i\}\times([s,s+\ell-1]\cap[1,N]).
\]
All these blocks are nonempty. Clipped duplicates are distinct update labels.
There are M=m(N+ell-1) labels. Direct counting gives:

- Every interior coordinate lies in exactly ell blocks, so its old mismatch
  is removed exactly ell times in the sum over update labels.
- An interior coordinate is a temporal exterior endpoint of at most two
  labels, one on each side.
- It is a spatial exterior coordinate for at most 2ell labels.
- A fixed boundary coordinate at row 0 or L is a temporal endpoint of
  exactly ell labels, including when N<ell.

Let D=r/epsilon+8a|theta|ell^2. T2 therefore bounds the total new interior
mismatch attributable to one changed interior exterior coordinate by
\[
2D+8a|\theta|\ell^2
=2r/\epsilon+24a|\theta|\ell^2=\ell\alpha,
\]
where
\[
\alpha=\frac{2r}{\epsilon\ell}+24a|\theta|\ell
\leq\bar\alpha:=\frac2{\epsilon J}+168|\theta|J.
\]
The last inequality uses ar<7. It is this choice ell=J ceil(6/a) that
removes the bad fine-cutoff scaling. Fixing ell while a decreases would
not give the displayed uniform budget.

## YM45-T4: weighted boundary-memory comparison

Choose R>1 with R bar-alpha<1 and set
\[
\gamma=\frac{\operatorname{Log}_\Sigma R}{7J},\qquad
w_u=e^{\gamma a\min(u,L-u)}.
\]
The distance-to-boundary function changes by at most the time separation.
An influencing exterior coordinate and any site of its block have time
separation at most ell. Since a ell<7J, their weight ratio is at most R.
Thus T3's total outgoing influence becomes at most R ell bar-alpha times
the weight of its source coordinate.

Consider any two finite strip laws with arbitrary endpoint configurations.
Begin with any joint law with these marginals. Select a block label with
weight 1/M and replace that block by the coupled conditional laws of T2.
Each marginal is unchanged by this operation: integrate the normalized
conditional against its own exterior marginal. Denote expected weighted
Hamming distance after k proof updates by W_k. Summing the exact removals,
the weighted incoming costs and the at most 2m changed boundary coordinates,
\[
W_{k+1}\leq(1-d)W_k+F/M,\qquad
d=\frac{\ell(1-R\bar\alpha)}M>0,\quad
F=2m\ell RD.
\]
Also d<=1 because M>=ell. Consequently
\[
W_k\leq(1-d)^k W_0+\frac{F}{Md}
 \leq(1-d)^k W_0+\frac{2mRD}{1-R\bar\alpha}.
\]
W_0 is finite on this finite strip. The first term is the explicit
discarded geometric memory tail. A stationary joint distribution and an
infinite-volume comparison theorem are not required.

For any bounded real whole-row F at time n, one mismatch in that row
can change F by at most osc(F). Since every site of the row has weight
w_n, at each proof step
|E F-E' F|<=osc(F) W_k/w_n. Let k grow to obtain
\[
\boxed{
|E^{x,z}F-E^{x',z'}F|
\leq C_{m,a}\operatorname{osc}(F)
e^{-\gamma a\min(n,L-n)},\qquad
C_{m,a}=\frac{2mRD}{1-R\bar\alpha}.}
\]
This prefactor can depend on m and a. The exponential rate cannot.
The next extraction removes the prefactor from the operator bound.

## YM45-T5: construct a fine-step vacuum from a coarse power

No small-a lower floor from the single-slice coefficient-tail estimate
is assumed. For fixed m,a, positivity and the total potential budget in
r steps give the pointwise kernel inequalities
\[
e^{-arb}K_{ar}^{\otimes m}\leq T_a^r
\leq e^{arb}K_{ar}^{\otimes m}.
\]
T1 supplies strictly positive finite lower and upper bounds for this
coarse power. YM43-T4's normalized iteration, with its explicit geometric
tail, constructs bounded h>0, bounded away from zero, such that
T_a^r h=nu h.

Its normalized Doob map P_r f=T_a^r(hf)/(nu h) has a strictly positive
common-part overlap for any two starting rows. Hence osc(P_r f)<=zeta
osc(f) for some zeta<1, and every bounded fixed function is constant.
The function f=T_a h/h is bounded and positive: the heat map is Markov,
the potential is bounded and h has positive upper/lower bounds. Commutation
of T_a with its power gives P_r f=f. Therefore
\[
T_a h=\lambda_a h,\qquad \lambda_a>0,\qquad \lambda_a^r=\nu.
\]
The Doob map P_a f=T_a(hf)/(lambda_a h) preserves 1 and is symmetric
and contractive in the positive form with normalized weight h^2 Phi.
This follows by conditional cut-square Cauchy–Schwarz, exactly as in YM43-T4.
In particular lambda_a=||T_a||. No spectral-existence theorem is a premise.

## YM45-T6: all sources, the actual ordering, and the free factors

Finite P_a-path laws are consistent by finite integration. Conditioned
on both end rows, their h and lambda endpoint factors cancel, leaving
exactly the strip weight used in T4. Set L=2n and mix the conditional
expectations over the final row, separately for each initial row. Since
T4 bounds every pair of endpoints, every pair of such mixtures differs by
at most the same bound. Thus
\[
\operatorname{osc}(P_a^n F)\leq
C_{m,a}\operatorname{osc}(F)e^{-\gamma an}.
\]
This step does not interchange a strip-length limit with a width limit.

For a mean-zero bounded F in the invariant h^2 form, put
u_n=||P_a^n F||^2. Symmetry and cut-square Cauchy–Schwarz imply
u_n^2<=u_(n-1)u_(n+1). If u_0,u_1>0, this gives
u_n>=u_0(u_1/u_0)^n. On the other hand
\[
u_n=\langle F,P_a^{2n}F\rangle
\leq C_{m,a}\,\pi(|F|)\operatorname{osc}(F)e^{-2\gamma an}.
\]
Geometric growth would contradict this inequality if
u_1/u_0>e^(-2 gamma a). Zero cases are immediate. Therefore
||P_a F||<=e^(-gamma a)||F||. Density of the bounded core and contraction
extend it to the completed carrier; real and iota components give the
complex version. Multiplication by h transports it to T_a/lambda_a on
its vacuum-orthogonal source cut.

Now X=K_(a/2)e^(aB/2) gives T_a=X dagger X and S_a=X X dagger.
YM43-T6 constructs g_a=Xh/sqrt(lambda_a) and transports the complementary
bound to **S_a**, with the same top value and ratio. In particular,
\[
\|Q_a(S_a/\lambda_a)^nQ_a\|\leq e^{-\gamma an}
\]
for every integer n>=1, every finite m and every 0<a<=1.

The independent B-channels in YM-16 have free ratio e^(-3a/4).
The gate implies R<epsilon J/2, so
gamma=(log R)/(7J)<epsilon/14=2/35<3/4. Every nonvacuum block of the
orthogonal A/B tensor decomposition therefore obeys the same ceiling.
This is an arbitrary-source bound, not a selected-probe correlation claim.

## YM45-T7: an open interval and explicit certified cells

With epsilon=4/5, J=5 gives bar-alpha=1/2+840|theta|.
For every |theta|<1/1680 choose
\[
R=\tfrac12(1+\bar\alpha^{-1})>1.
\]
Then R bar-alpha=(1+bar-alpha)/2<1 and T4 gives a positive gamma.
This algebra proves the whole open interval; sampled cells are not used
to infer an unsampled continuum of parameters.

The following alternative rational choices give explicit rate enclosures:

| abs(theta) | J | R | 1-R bar-alpha | gamma, rounded down |
| --- | --- | --- | --- | --- |
| 1/8192 | 11 | 2 | 531/5632 | 0.009001911435 |
| 1/4096 | 8 | 4/3 | 7/48 | 0.005137179865 |
| 1/2048 | 6 | 21/20 | 117/2560 | 0.001161670575 |
| 1/1792 | 5 | 33/32 | 1/1024 | 0.000879190247 |

Signs of theta are immaterial to these absolute-value bounds.
For theta>=0 the identity
\[
J(\bar\alpha-1)=\frac{(J-5)^2}{10}
 +(168\theta-\tfrac1{10})J^2
\]
shows why this particular budget cannot pass at or above 1/1680, for any J.
Failure of this sufficient criterion is **not** a proof that the actual
gap vanishes there. The overlap floor, tilt comparison and weights are
deliberately conservative.

## YM45-T8: remove the time cutoff at fixed width

YM-44 has already constructed U(t) and shown
e^(-bt)K_t<=U(t)<=e^(bt)K_t pointwise. At t=6, T1 and YM43-T4
construct a bounded positive vacuum g for U(6), bounded away from zero.
Its Doob map has a strict common-part overlap. For any t>0,
U(t)g/g is a bounded positive fixed function for this map, by the
semigroup and its commutation. As in T5 it is constant, say rho(t)>0.
The symmetric Doob cut-square argument then gives
U(t)g=rho(t)g and ||U(t)||=rho(t). Thus the same positive top vacuum
is available at every t>0, without requiring a small-t coefficient floor.

For fixed m,t let a=t/n<=1 and A_n=S_a^n. T6 gives the complementary
ratio q=e^(-gamma t), independent of n and m. Testing S_a on 1 gives
lambda_a>=e^(-ab), hence ||A_n||>=ell_0:=e^(-bt)>0. YM44-T3 gives
\[
\delta_n=\|A_n-U(t)\|\leq3bt\,e^{bt}\sqrt{t/n}\longrightarrow0.
\]
For all sufficiently large n, 2 delta_n<ell_0(1-q). The already proved
normalized transport gate YM44-T5 therefore yields
\[
\frac{\|Q U(t)Q\|}{\|U(t)\|}
\leq \frac{q\ell_0+\delta_n}{\ell_0-\delta_n}.
\]
Let n grow to obtain the second boxed bound at the start of this chapter.
The necessary refinement depth may depend on m; gamma does not.
This is a uniform family of finite-width time-limit gap bounds, not an
exchange of time and infinite-volume limits. The latter state/observable
construction remains a separate task.

## Evidence, boundaries and reproduction

Eight general results above are written proofs under the stated adapter,
not mechanically formalized proofs. The exact rational/outward-interval
certificate checks:

- The unchanged full heat tail, four rate cells, fine-step budget samples,
  the interval gate and its endpoint identity.
- Eighteen conditioned skeleton couplings, twelve direct path marginals,
  twenty-five independent matrix-power identities, twenty tilted blocks
  and sixty single-coordinate spatial perturbations.
- Seventy-five clipped layouts, 910 interior incidence checks and 1,260
  weighted influence checks, including strips shorter than their blocks.
- Two exact joint block-update fixtures, 4,608 pair/label updates each,
  with both marginals preserved and weighted mismatch decreased.
- A coarse-power vacuum control, 24 fixed-width limit controls, and
  twelve refusal groups, including the future normalizer, missing temporal
  tilt costs, fixed fine-block length and omitted vacuum denominator.

The binary fixtures are independent controls for proof operations, not an
SU(2) replacement. The SU(2) input is the full coefficient tail and the
declared positive kernel, with no finite-content truncation. Source pins bind
the preceding results and these new proof/code/test/workflow inputs.

| Obligation | Status after YM-45 |
| --- | --- |
| All-source normalized gap, every finite m and 0<a<=1 | Proved for abs(theta)<1/1680 on declared kappa=theta a chain |
| All-content bound | Full kernel; no content cutoff introduced |
| Time-limit gap for all t>0 at each finite m | Proved using YM-44 with the same width-independent rate |
| Uniform positive per-step gap 1-q as a tends to zero | Not claimed; q(a)=exp(-gamma a) tends to one |
| Old theta=1/16 cell, arbitrary interactions, Wilson-grid E4D-C | Not closed by this small-bridge criterion |
| Infinite-volume state and order of limits | Open |
| Native Phi_Sigma/NCG dictionary, physical clock and AF trajectory | Open |
| Full 4D interacting continuum, Clay and quantum gravity | Not established |

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym45_temporal_blocks.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym45_temporal_blocks.py -v
~~~

Earlier certificate bytes and the canonical RKF operator engine are preserved.
