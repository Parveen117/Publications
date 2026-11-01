# YM-46: an infinite-volume local-observable state and its time map

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

**Result under the existing heat/positive-functional adapter:** for each
fixed 0<a<=1 and abs(theta)<1/1680, the interacting chain has a constructed,
boundary-independent infinite-volume positive functional on local space-time
observables. Its normalized row time map exists as a uniform limit on every
local observable and extends to the uniform completion. This is proved for
both the convenient T ordering and the actual square-sourced S ordering.
The completed positive forms retain YM-45's all-source gap,
\[
\|Q(P_a^\sigma)^nQ\|\leq e^{-\gamma an},
\qquad \sigma\in\{T,S\},\quad n\geq1,
\]
and the regulated history readouts have site- and bond-reflection positivity.

The new ingredient is a boundary estimate with spatial and temporal decay,
local-observable support factors and explicit Cauchy bounds. A uniform gap
alone would not have constructed the state. The present boundary constants
depend on a; this chapter does not exchange infinite-volume and a->0 limits.

## Declared source, readout and completion

Reuse YM-9/12/19/43/44/45's full SU(2) heat kernel K_t, finite reference
products Phi, positive normalized conditionals, coefficient semigroup and
bounded measurable functions. Native source, cut and positive-form completion
retain the same contract. The reference integral is an **admitted adapter**,
not a newly derived Phi_Sigma or NCG measure.

The spatial carrier is the countable chain of sites indexed by integers.
An interval I has the same open-chain interaction as before,
\[
V_I=\sum_{\{i,i+1\}\subset I}\tfrac12\chi_{1/2}(x_ix_{i+1}^{-1}),
\quad B_I=\theta M_{V_I},\quad \kappa(a)=\theta a,
\]
\[
T_{I,a}=e^{aB_I/2}K_a^{\otimes I}e^{aB_I/2},\qquad
S_{I,a}=K_{a/2}^{\otimes I}e^{aB_I}K_{a/2}^{\otimes I}.
\]
Finite normalized vacua and Doob time maps already exist by YM45-T5/6.
The gap rate is in declared heat-parameter units. Neither the site index
nor this parameter has acquired a derived physical metric or clock.

A local observable is a bounded function of finitely many site/time
coordinates. Its support is a specified finite set D, and osc(F) is the
supremum of |F(x)-F(y)|. The local-observable algebra and its uniform-norm
completion are the targets. A **state** here means a normalized positive
linear functional on this algebra/completion. No infinite product reference
measure, Kolmogorov extension theorem or spectral-existence theorem is a
construction premise. A separate countably additive path-measure representation
is not asserted.

YM-33's low-order spatial edge corrections, YM-35's local commutator
remainder and YM-37/42's finite-rail spatial transfer are retained at their
proved scopes. They are not relabelled as an all-content state construction.
YM-43/45's comparison operations, YM-44's fixed-width refinement and the
canonical RKF completion discipline are reused. The general comparison and
local-specification method has established Dobrushin–Shlosman/DLR lineage;
see [Rebeschini–van Handel, *Comparison Theorems for Gibbs Measures*](https://web.math.princeton.edu/~rvan/dobrushin130819.pdf).
The specific proof is below; no priority claim or external comparison theorem
is used in place of it.

## YM46-T1: a joint spatial/time gate on the same open window

Set epsilon=4/5, r=ceil(6/a), ell=Jr and eta=a|theta|.
Choose J a positive integer and R,rho>1 with
\[
\boxed{\beta=
R\left(\frac2{\epsilon J}+(112+56\rho)|\theta|J\right)<1.}
\]
Write gamma=Log_Sigma(R)/(7J), q=e^(-gamma a) and
D_0=r/epsilon+8 eta ell^2. The underlying YM-45 gate holds because
112+56rho>168, so its rate gamma applies at every finite width.

This stronger gate still covers the **whole** open interval
abs(theta)<1/1680. To see it, take J=5 and
alpha_0=1/2+840|theta|<1. For theta nonzero choose
\[
\rho=1+\frac{1-\alpha_0}{112|\theta|J},\quad
A=\frac2{\epsilon J}+(112+56\rho)|\theta|J
  =\frac{1+\alpha_0}{2}<1,\quad
R=\tfrac12(1+A^{-1}).
\]
Then rho,R>1 and beta=(1+A)/2<1. At theta=0 use rho=2, J=5,
A=1/2 and R=3/2. Thus no smaller interaction interval is silently imposed.

Four concrete choices retain the preceding four time-gap rates:

| abs(theta) | J | R | rho |
| --- | --- | --- | --- |
| 1/8192 | 11 | 2 | 3/2 |
| 1/4096 | 8 | 4/3 | 3/2 |
| 1/2048 | 6 | 21/20 | 5/4 |
| 1/1792 | 5 | 33/32 | 257/256 |

For example the second cell has 1-beta=7/96 and
gamma=Log(4/3)/56>=0.005137179865. Signs of theta do not change these bounds.

## YM46-T2: local boundary memory with explicit tails

Use a rectangle with spatial interval I=[l,h], interior time rows
1,...,N and fixed endpoint rows 0,N+1. At the two spatial sides allow
either the declared bond to a prescribed exterior history or an absent
bond (the open boundary). Comparing a present and absent bond costs no
more than comparing two arbitrary exterior values: the changed log-weight
at one time has absolute difference at most 2 eta.

For one target z=(i_0,u_0) assign the weight
\[
w_z(i,u)=\rho^{-|i-i_0|}q^{|u-u_0|}.
\]
For support D use w=sum_(z in D) w_z. Then w>=1 at every support site.
Use precisely YM-45's clipped block labels, including their duplicates.
Their removal count is ell at every interior coordinate.

A temporal exterior coordinate influences at most two labels, costing
D_0 per label. Its weight changes by at most q^(-ell)<=R within a block.
A spatial exterior coordinate influences at most 2ell labels, costing
4 eta ell per label; its weight changes by at most rho R. Consequently
the total weighted outgoing influence is bounded by
\[
R\{2r/\epsilon+(16+8\rho)\eta\ell^2\}\,w
\leq\ell\beta\,w.
\]
This explains the additional spatial weight in T1. Omitting rho would use
YM-45's temporal budget for a different assertion.

At either time boundary each coordinate meets ell clipped labels.
At either spatial boundary each time coordinate also meets ell labels.
The corresponding forcing costs, before division by the removal margin,
are bounded respectively by
\[
\ell R D_0\,w(i,u),\qquad
4\eta\ell^2\rho R\,w(i,u).
\]
These counts hold even if N<ell. Define
\[
H_\rho=\frac{\rho+1}{\rho-1},\qquad
H_q=\frac{1+q}{1-q},\qquad
C_s(a)=\frac{4\eta\ell\rho R H_q}{1-\beta},\qquad
C_t(a)=\frac{R D_0 H_\rho}{1-\beta}.
\]
The H factors are the exact sums of the two-sided geometric weights.

For a target z let d_l(z)=i_0-l+1 and d_h(z)=h-i_0+1 be its distances
to the exterior spatial columns. Its distances to the time endpoints
are u_0 and N+1-u_0. For any two lawful boundary data b,b',
\[
\boxed{
|\omega_\Lambda^b(F)-\omega_\Lambda^{b'}(F)|
\leq\operatorname{osc}(F)\sum_{z\in D}
\left[
C_s(a)(\rho^{-d_l(z)}+\rho^{-d_h(z)})
C_t(a)(q^{u_0}+q^{N+1-u_0})
\right].}\tag{1}
\]
Terms for unchanged boundary faces can be omitted.

**Proof of the limiting comparison operation.** Start any coupling of
the two exact finite laws and make the coupled conditional block updates
of YM45-T2. Each marginal is invariant under its own conditional replacement.
For M=|I|(N+ell-1) labels, weighted mismatch obeys
\[
W_{k+1}\leq(1-d)W_k+F_\partial/M,\qquad
d=\ell(1-\beta)/M\in(0,1].
\]
Thus W_k<=(1-d)^k W_0+F_partial/[ell(1-beta)].
The expectation difference is at most osc(F)W_k because w>=1 on D.
Sum the displayed boundary costs using the two geometric series, then let
k grow. The discarded term is the explicit geometric Smriti tail
(1-d)^k W_0. No stationary coupling is assumed.

## YM46-T3: construct the positive space-time functional and prove uniqueness

Fix a>0, theta and a finite-support F. Compare its readout in two
rectangles containing D. Condition the larger rectangle on the complement
of the smaller one. Nearest-neighbour factorization gives exactly one of
the smaller rectangle's lawful boundary specifications. The bound (1) is
uniform in that boundary, so it also bounds its mixture. Open boundaries
are included by the absent-bond comparison above.

As the distances of D to all four faces grow, (1) tends to zero. For
non-nested rectangles compare both with their common interior rectangle.
This proves the Cauchy property, independence of exhaustion, and independence
of prescribed time endpoints and prescribed or open spatial sides. Define
\[
\omega_a^T(F)=\lim_{\Lambda\uparrow\mathbb Z^2}\omega_\Lambda^b(F).
\]
Finite linearity, normalization, positivity and the bound
|omega_Lambda(F)|<=||F||_infinity pass to the limit. If F_n is a
uniformly Cauchy sequence of local observables, the last bound gives a
unique continuation to the uniform completion. This constructs the stated
positive functional rather than postulating an infinite product measure.

For an explicit modulus, if all spatial face distances are at least d and
time face distances at least n, the error is at most
\[
2|D|\operatorname{osc}(F)\{C_s(a)\rho^{-d}+C_t(a)q^n\}.
\tag{2}
\]
Increase d and n until each term is below half the requested tolerance.
The certificate gives outward examples. This is the native completion's
declared error, not a bare assertion that a limit exists.

Let Gamma_Lambda replace a finite region by its normalized local conditional
law, leaving exterior coordinates untouched. It sends a local F to a local
function, since only finitely many boundary coordinates enter. Finite
conditional integration gives omega_B Gamma_Lambda F=omega_B F once B
contains this support and Lambda. Passing to the limit yields
\[
\omega_a^T\Gamma_\Lambda=\omega_a^T.
\]
Conversely, for any normalized positive functional nu satisfying these
finite conditional identities, nu(F)=nu(Gamma_Lambda F). The pointwise
range of Gamma_Lambda F differs from any fixed-boundary readout by at most
(1). Positivity bounds the same difference after applying nu. Sending
Lambda outward gives nu(F)=omega_a^T(F). This is uniqueness **within the
stated local-specification class**, without a measure-representation premise.
Spatial/time translation and time reversal follow from that uniqueness and
the symmetries of the finite weights.

## YM46-T4: the normalized T time map is a uniform local-observable limit

Let P^T_(I,a) be the finite-volume Doob map of T, with invariant row
functional pi^T_(I,a). Extend it to a full spatial configuration by updating
I and leaving the exterior arguments of a function fixed. This auxiliary
extension is a positive unital uniform-norm contraction and obeys the finite
semigroup law. Denote it by hat-P_I.

For a row observable F supported in A subset I, compare hat-P_I^n F(x)
with hat-P_J^n F(x), J containing I, for the **same arbitrary initial row x**.
Use finite paths through a future row L>n. Condition the larger path on
exterior spatial histories and on its final row; the vacuum endpoint factors
cancel from the conditional strip. The common initial row causes no bottom
boundary forcing. Equation (1) leaves the spatial forcing and a top term
which vanishes as L grows. Its bounds do not depend on x or n. Therefore
\[
\boxed{
\|\widehat P_I^nF-\widehat P_J^nF\|_\infty
\leq\operatorname{osc}(F)C_s(a)
\sum_{i\in A}(\rho^{-d_l(i)}+\rho^{-d_h(i)}),\qquad n\geq1.}
\tag{3}
\]
In particular the limit P^T_(a,n)F is uniform in x, and the same error
holds for every integer n. A uniform limit of finite-support functions
belongs to the uniform completion, which is the meaning of quasilocal here.

Extend these contractions from the dense local core to its uniform completion.
For any F there, uniform boundedness and approximation imply
hat-P_I^n F -> P^T_(a,n)F. Now
\[
\widehat P_I^{n+k}F=\widehat P_I^n\widehat P_I^kF
\longrightarrow P^T_{a,n}P^T_{a,k}F:
\]
bound the inner-map replacement by its uniform error and use convergence
of the outer map on the fixed function P^T_(a,k)F. Thus
P^T_(a,n)=(P^T_a)^n. The identity at n=0, positivity and preservation of
1 follow directly. This proves composition, not merely convergence of
selected two-point correlations.

Let the time length tend to infinity first at fixed I in the finite
vacuum path law. Equation (1) identifies its local readouts with those of
omega_a^T as I grows. Hence pi^T_(I,a) converges on every local row F to
the row restriction pi^T_a of omega_a^T. To handle quasilocal F one may
extend each finite pi_I with any fixed exterior configuration and use uniform
approximation; the limit is independent of that exterior.
Passing the finite invariance and detailed-balance identities gives
\[
\pi_a^T P_a^T=\pi_a^T,\qquad
\pi_a^T(\overline F P_a^T G)=
\pi_a^T(\overline{P_a^T F}\,G).
\]

## YM46-T5: the actual square-sourced row ordering

Write H=K_(a/2), so K_a=H^2. Given adjacent T-layer configurations x,z,
insert at each site a midpoint y with normalized bridge density
\[
b(dy\mid x,z)=\frac{H(x,y)H(y,z)}{K_a(x,z)}\,\Phi(dy).
\tag{4}
\]
The denominator is strictly positive on the declared carrier, and
normalization is exactly the coefficient semigroup. For finitely many
midpoint observables use the finite product of these conditionals.
The resulting map J_a sends a local history F to a bounded local T-history
function, preserves 1 and positivity, and has norm at most one.

Define omega_a^S(F)=omega_a^T(J_a F). Consistent normalization in (4)
makes this definition independent of extra unused midpoint coordinates.
It is also the limit of the finite S vacuum path readouts. To verify the
ordering, the interleaved weight with fixed midpoint endpoints is
\[
\prod_j H(y_j,x_j)e^{aB_I(x_j)}H(x_j,y_{j+1}).
\]
Integrating x_j gives the actual S kernel. Integrating the interior y_j
instead gives full K_a bonds between the x_j, with H bonds only at the
two endpoints. This is an exact finite product identity. It does not set
T equal to S. At a single midpoint, integration of the finite stationary
T path gives the squared vacuum
g_I=H e^(aB_I/2)h_I/sqrt(lambda_I), exactly YM43-T6.

For a **uniform time-map limit** an endpoint check is also necessary;
convergence of stationary midpoint correlations alone would be insufficient.
Condition an S path on its initial/final y rows and integrate unobserved
midpoints. The x strip has spacing a internally and half-edges at its
two time boundaries. Internal block influences are still those of T2.
For a changed half-edge endpoint, the first skeleton segment can use r+1
sites, with elapsed heat parameter (r+1/2)a>=6. Subsequent segments have r
sites. The free mismatch cost is at most
\[
(r+1)+(1-\epsilon)r/\epsilon=1+r/\epsilon.
\]
If the whole block is shorter, its site count gives the same bound.
Thus only the fixed time-boundary cost changes from D_0 to D_0+1;
the interior contraction beta and side constant C_s are unchanged.
For a block with an interior endpoint, start the skeleton from that full
edge; its opposite half-edge remains a common future weight. The original
r/epsilon bound applies. This covers clipped blocks touching a time boundary.

A midpoint row observable with spatial support A becomes J_aF supported
on at most two adjacent x rows at those sites, with osc(J_aF)<=osc(F).
Applying the half-edge version of (1), using the same initial y row and
then removing the far future boundary, gives
\[
\boxed{
\|\widehat P_{I}^{S,n}F-\widehat P_{J}^{S,n}F\|_\infty
\leq2\operatorname{osc}(F)C_s(a)
\sum_{i\in A}(\rho^{-d_l(i)}+\rho^{-d_h(i)}).}
\tag{5}
\]
It is uniform in the initial row and in n>=1. The proof of T4 now constructs
the positive unital semigroup (P_a^S)^n on the row uniform completion.
Its invariant row functional pi_a^S is the row restriction of omega_a^S,
and it is symmetric in that positive form. The factor two and half-edge
normalizers are retained explicitly.

## YM46-T6: positive-form completion and the inherited all-source gap

For either ordering sigma, conditional cut-square Cauchy–Schwarz gives
|P_I F|^2<=P_I|F|^2 pointwise. Uniform convergence passes this inequality
to P_a^sigma. Invariance then gives
\[
\pi_a^\sigma(|P_a^\sigma F|^2)\leq\pi_a^\sigma(|F|^2).
\]
Take the null-seam quotient of the row algebra under this positive form
and complete it. The contraction descends and extends; finite detailed
balance gives self-daggerness. The finite transfers are manifest-square
positive, so pi_I(overline F P_I F)>=0. State and map convergence give
the same positivity on the completed infinite-volume carrier.

For a local F, YM-45 gives at every finite I
\[
\pi_I(|P_I^nF-\pi_I(F)|^2)
\leq e^{-2\gamma an}
\{\pi_I(|F|^2)-|\pi_I(F)|^2\}.
\]
Use (3) or (5) and state convergence to pass every term to the limit.
Then extend by positive-form density. With Q the complement of the unit
constant source, this proves the first displayed operator bound, on **all**
completed sources. In particular the fixed-source subspace of P consists
only of constants. The independent free B factors of YM-16 can be included
through their finite products and the same local completion; their rate
3/4 exceeds gamma and leaves the inherited ceiling unchanged.

## YM46-T7: regulated time correlations and reflection positivity

Cut-square Cauchy–Schwarz and T6 imply, for completed sources F,G,
\[
|\pi(\overline F P^nG)-\overline{\pi(F)}\pi(G)|
\leq e^{-\gamma an}
\sqrt{\operatorname{Var}_\pi(F)\operatorname{Var}_\pi(G)}.
\]
This is a normalized discrete-time statement on the constructed state.
It does not identify a transfer step with unitary physical real time.

For bounded local histories F in the future half, let Theta reflect the
time index and conjugate the value. In each finite stationary reversible
row chain, reflection through a time site gives
omega((Theta F)F)=pi(|f|^2)>=0, where f is the future conditional readout
given the row at that site. Reflection through a time bond gives
omega((Theta F)F)=pi(overline f P f)>=0, using the transfer's positive form.
The identities follow by conditioning the past/future path products at the
cut and using finite detailed balance. They hold for either ordering.

The history products have finite support, so T3/T5 permit taking the
volume limit of these nonnegative readouts. Both reflection inequalities
therefore survive at every fixed a. This extends the existing square-sourced
positivity to the regulated infinite-volume history functional. It is not
a claim that the remaining continuum reconstruction axioms have been proved.

## YM46-T8: the remaining limit is explicit

The gap rate gamma is independent of a, but the present spatial/state
error constants are not. Already
\[
C_t(a)\geq
\frac{6R H_\rho}{\epsilon(1-\beta)a}.
\]
For theta nonzero, eta ell>=6|theta|J and
H_q>=(gamma a)^(-1), so also
\[
C_s(a)\geq
\frac{24|\theta|J\rho R}{(1-\beta)\gamma a}.
\]
These are divergences of the proved **budgets**, not proofs of a failure
of the actual physical limit. They show precisely why a fixed spatial
truncation followed by a->0 is not controlled by (1)–(5).

YM-44 constructs continuous heat-parameter evolution at each finite width;
YM-46 constructs infinite-volume normalized discrete evolution at each
fixed a. Their two constructions are not yet identified. A cutoff-uniform
local boundary estimate, or another joint Cauchy argument, must connect
them before exchanging the limits or claiming a continuous-time
infinite-volume state/evolution. A conserved native process, the Phi_Sigma/
NCG dictionary, a physical trajectory, and a four-dimensional gauge theory
remain separate obligations.

## Evidence and reproduction

The eight general results are written proofs under the declared adapter.
The exact certificate checks four joint parameter cells, explicit local
Cauchy boxes, two-direction weighted incidence including short strips,
independent finite strip readouts with positive/negative interactions,
local specification consistency, normalized midpoint/path identities,
half-edge skeletons, finite positive time/reflection forms and refusal
witnesses. Binary fixtures test these operations; they do not replace SU(2)
or mechanically formalize the infinite theorem.

Source pins bind this proof, verifier, tests and workflow, plus the unchanged
upstream evidence. Previous certificate bytes and the canonical RKF engine
are preserved. The native measure dictionary is not certified by a new name
for the existing reference functional.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym46_infinite_volume.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym46_infinite_volume.py -v
~~~
