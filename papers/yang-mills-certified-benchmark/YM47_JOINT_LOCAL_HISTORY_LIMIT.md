# YM-47: joint volume and time refinement of local history readouts

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

**Result on the declared heat/positive-functional carrier:** for
abs(theta)<1/1680, the actual square-sourced chain has a unique limit of its
finite-volume vacuum history readouts as the interval exhausts the spatial
chain and a decreases to zero. No relation between those two rates is
required. Both iterated limits give this same positive history functional.
It retains the finite-volume correlation gap and reflection positivity.

The main improvement is to normalize **each** step before telescoping.
The resulting error is proportional to width times sqrt(a), with an
exponential depending on width times a, rather than width times the whole
observation time. Summable one-boundary extensions then construct the
continuous-parameter volume limit. Cropping very large intervals to a
logarithmic window proves the unrestricted joint statement.

This is a limit of local history readouts in the declared heat parameter.
An infinite-volume strongly continuous time-map/generator construction,
a derived physical clock, the native Phi_Sigma/NCG measure dictionary,
four-dimensional continuum Yang--Mills and Clay/QG remain separate.

## Source contract and one typographical clarification

Reuse the complete carrier, finite reference products, null-seam quotient,
completion and open-chain interaction of YM-44/45/46:
\[
B_I=\theta M_{V_I},\quad b_I=|\theta|(|I|-1),\quad
S_{I,a}=K_{a/2}^{\otimes I}e^{aB_I}K_{a/2}^{\otimes I}.
\]
The reference integral remains an admitted SU(2) adapter. The trajectory
kappa(a)=theta a is declared, not derived. All contents and sources of this
carrier are retained. No second shared RKF operator engine is introduced.

YM-44 constructs U_I(t). YM45-T8 supplies its common positive unit vacuum
h_I and its normalized all-source gap for every t>0. Write
\[
A_{I,a}=S_{I,a}/\lambda_{I,a},\qquad
V_I(t)=U_I(t)/\nu_I(t),\qquad \nu_I(t)=\|U_I(t)\|.
\]
The common vacuum implies nu_I(s+t)=nu_I(s)nu_I(t): evaluate the
semigroup identity on h_I. Thus V_I is a positive self-adjoint contraction
semigroup. Denote the positive unit vacuum of A_(I,a) by h_(I,a).

**Forward clarification of a frozen source:** the displayed equation (1)
in YM-46 omits the plus sign at the line break between its spatial and
temporal terms. Its proof adds the boundary forcing terms, and equation
(2) already uses their sum. The intended and used expression is
\[
\operatorname{osc}(F)\sum_{z\in D}
\big[C_s(a)(\rho^{-d_l(z)}+\rho^{-d_h(z)})
+C_t(a)(q^{u_z}+q^{N+1-u_z})\big].
\]
Unchanged faces contribute zero. This chapter uses only the spatial part
after removing the time endpoints. Prior proof/certificate bytes stay
frozen; the clarification is explicit rather than silently rewriting them.

Choose any YM46-T1 parameters J,R,rho with beta<1, rho,R>1, and let
\[
\gamma=\log(R)/(7J),\quad 0<g\leq\gamma,\quad d=g/(1+g),\qquad
C=\frac{28|\theta|J\rho R}{1-\beta}(1+2/g).
\]
Here g can be a certified rational lower bound. The gap at heat time at
least one is at least d, because 1-exp(-gamma)>=g/(1+g). Also
\[
\boxed{C_s(a)\leq C/a,\qquad 0<a\leq1.}\tag{1}
\]
Indeed a ell<=7J and
H_q=1+2/(exp(gamma a)-1)<=1+2/(gamma a)<=(1+2/g)/a.
These are upper bounds; YM-46's divergent budget was not evidence that
the actual local limit failed.

## YM47-T1: normalization before accumulation

For nonzero bounded operators X,Y, with lambda=||X|| and nu=||Y||,
\[
\|X/\lambda-Y/\nu\|\leq2\|X-Y\|/\lambda.
\]
Subtract using denominator lambda first, and use
|lambda-nu|<=||X-Y||. No equality of the vacuum normalizations is assumed.

YM-44 gives ||S_a-U(a)||<=3b a^(3/2) exp(ba), while
lambda_a>=exp(-ba). Hence
\[
\|A_{I,a}-V_I(a)\|\leq6b_Ia^{3/2}e^{2b_Ia}.
\]
Every normalized factor is a contraction. Telescoping at a fixed width
therefore gives, for every integer n>=0,
\[
\boxed{\|A_{I,a}^{\,n}-V_I(na)\|
\leq6b_I(na)\sqrt a\,e^{2b_Ia}.}\tag{2}
\]
This is an all-source operator bound. It does not require S_a to equal
its interacting subdivisions, and it does not commute any interaction
through a heat factor.

## YM47-T2: the vacuum also converges with an explicit bound

Let A,C be positive self-adjoint contractions with unit top vacua h,k
at value one. Suppose ||Q_h A Q_h||<=q<1 and ||A-C||<=epsilon.
Choose the phase so <h,k>>=0, as already holds for our positive vacua.
Writing w=Q_h k gives
\[
(I-A)w=Q_h(C-A)k.
\]
The geometric inverse on h's complement has norm at most 1/(1-q).
Consequently
\[
\|w\|\leq\epsilon/(1-q),\qquad
\|h-k\|^2=2(1-\langle h,k\rangle)\leq2\|w\|^2.
\]
In particular ||h-k||<=2 epsilon/(1-q). This proof does not infer
vacuum-vector convergence from a gap alone.

Apply (2) with n=ceil(1/a), so 1<=na<=2. Both finite vacua have already
been constructed, and YM-45 supplies q<=exp(-gamma). Thus
\[
\boxed{\|h_{I,a}-h_I\|
\leq24b_I\sqrt a\,e^{2b_Ia}/d.}\tag{3}
\]
For b_I=0 the two vacua coincide. Large bounds may always be replaced
by the trivial bound two; smallness is needed for a useful error, not
as an existence premise.

## YM47-T3: rounding observation times without norm continuity at zero

For a positive self-adjoint contraction semigroup V, and t,s>0,
\[
\boxed{\|V(t+s)-V(t)\|\leq\min(1,s/t).}\tag{4}
\]
Here is a proof on the existing completed source space. Put W=V(s) and
n=floor(t/s). The operators W^j(I-W) are positive and decrease with j.
Positivity follows by inserting V(js/2) on both sides of I-W.
Their consecutive differences are W^j(I-W)^2, also positive. Therefore
\[
(n+1)W^n(I-W)\leq\sum_{j=0}^nW^j(I-W)
=I-W^{n+1}\leq I.
\]
The positive-form Cauchy--Schwarz bound gives norm at most 1/(n+1).
Multiply by the remaining contraction V(t-ns). For n=0 use the
trivial bound one. This proves (4) without a spectral representation.

Observe a finite ordered list 0=t_0<t_1<...<t_r=T and put
Delta=min_j(t_j-t_(j-1)) when r>0. For r=0 set Delta=1.
Round each t_j to a multiple s_j(a) of a with s_0=0 and error <=a.
If a<=Delta/4, the rounded times remain distinct; each positive
separation changes by at most 2a. Equation (4) bounds each corresponding
operator change by 2a/(Delta-2a)<=4a/Delta.
The conclusion is independent of the chosen such rounding rule.

This does **not** assert ||V(t)-I||->0. The unbounded free carrier
still has ||K_t-I||=1 for every t>0.

## YM47-T4: finite-width history error

Take bounded local row observables F_0,...,F_r. Let
M=product_j ||F_j||_infinity, and let A be a fixed spatial interval
containing all their supports, with k=|A| and
p=sum_j |support(F_j)|. Constants may be removed, and the empty
observable is handled by normalization. Define
\[
\Gamma_I(F_0,\ldots,F_r;\mathbf t)
=\langle h_I,M_{F_0}V_I(t_1-t_0)M_{F_1}\cdots
V_I(t_r-t_{r-1})M_{F_r}h_I\rangle.
\]
Its discrete counterpart Gamma_(I,a) uses h_(I,a), A_(I,a) and the
rounded integer separations. Inserting the vacuum replacements costs
at most 2M times (3). Telescope the r contraction factors using (2)
and (4); their rounded total duration is at most T+2. For
0<a<=min(1,Delta/4),
\[
\boxed{|\Gamma_{I,a}-\Gamma_I|
\leq M\{B b_I\sqrt a\,e^{2b_Ia}+4ra/\Delta\},\qquad
B=48/d+6(T+2).}\tag{5}
\]
No pointwise lower bound on the vacuum or its Doob quotient is used.
Complex observables are allowed; all estimates use norms.

We also need a spatial comparison. For nested intervals I subset J,
whose changed faces are at distance at least n from A, YM46-T2/5 gives
for the actual S stationary history
\[
|\Gamma_{I,a}-\Gamma_{J,a}|
\leq4Mp\,C\,a^{-1}\sum_{\text{changed faces}}\rho^{-n}.
\tag{6}
\]
To see the factors, the product observable has oscillation at most 2M;
the normalized midpoint lift uses at most two T coordinates for each
of its p site/time coordinates. Its oscillation does not increase.
Condition outside I, apply the uniform face comparison, and remove
the temporal endpoints. An unchanged spatial face is omitted.
The constant is independent of the number of unobserved time rows.

## YM47-T5: a summable spatial extension error after time refinement

Choose u>1 with u^2<rho and use a_n=u^(-2n). Such a u exists for every
rho>1; for example u=2rho/(rho+1) satisfies the strict inequality.
We will use n large enough that a_n<=Delta/4.

Consider a one-site extension at one end, with changed-face distance
at least n from A and both widths at most k+2n+2. Equations (5)--(6)
and a triangle through this common mesh give
\[
|\Gamma_I-\Gamma_J|\leq M E_n,
\]
\[
E_n=2B|\theta|Z(k+2n+2)u^{-n}
+(8r/\Delta)u^{-2n}+4pC(u^2/\rho)^n,\quad
Z=\exp\{2|\theta|(k+2+2/(u^2-1))\}.
\tag{7}
\]
Indeed b_I,b_J<=|theta|(k+2n+2), and
n u^(-2n)<=1/(u^2-1) bounds their step exponent uniformly.
Each term in (7) is summable.

For a completely explicit error, put
\[
G(q,N)=q^N/(1-q),\quad
H(q,N)=q^N\left\{\frac{k+2+2N}{1-q}
+\frac{2q}{(1-q)^2}\right\}.
\]
Direct finite geometric sums and their vanishing remainders give
\[
\boxed{\mathcal E_N=\sum_{n\geq N}E_n
=2B|\theta|Z H(u^{-1},N)
+(8r/\Delta)G(u^{-2},N)+4pC G(u^2/\rho,N).}\tag{8}
\]
The certificate rounds exponentials and powers outward. It supplies
finite radii for prescribed tolerances, including the small spatial
decay of YM-46's near-edge cell.

## YM47-T6: continuous-parameter readouts have a volume limit

Let A_N extend A by N sites on each side. Passing from A_n to A_(n+1)
uses two one-site extensions satisfying T5. Thus Gamma_(A_N) is Cauchy,
with remaining error at most 2M Ecal_N. Define its limit Omega on the
specified product history.

An arbitrary interval I containing A has padding n_l,n_h and
N=min(n_l,n_h), taken large enough that u^(-2N)<=Delta/4.
First use A_N; then extend only its longer side to I.
At an extension labelled n>=N, the other padding is at most n and
the width is at most k+2n+2, so exactly the same T5 estimate applies.
The extra one-sided sum is at most M Ecal_N. Therefore
\[
\boxed{|\Gamma_I-\Omega|\leq3M\mathcal E_N.}\tag{9}
\]
This proves independence of spatial interval exhaustion, including
arbitrarily asymmetric intervals. It does not require a pre-existing
infinite-volume continuous-time state.

## YM47-T7: the unrestricted joint limit and both iterated limits

For 0<a<1 choose an integer D(a)>=1 with rho^(-D(a))<=a^2, taking the
least such integer. Then D(a)=O(log(1/a)). Crop an arbitrary I containing
A to J=I intersect A_(D(a)). Faces that were already closer remain
unchanged. Every changed face is at distance at least D(a).
Equation (6) therefore costs at most 8MpCa.

The cropped width is at most k+2D(a). Use (5) in J, then (9), and set
N=min(n_l,n_h,D(a)) and b_*=|theta|(k+2D(a)). Once N is large enough
for the tail in (8) and a<=Delta/4,
\[
\boxed{|\Gamma_{I,a}-\Omega|\leq M\left[
8pCa+B b_*\sqrt a\,e^{2b_*a}+4ra/\Delta
+3\mathcal E_N\right].}\tag{10}
\]
As a->0 and both spatial paddings tend to infinity, all four terms
vanish independently of their relative rates. The logarithmic crop
ensures b_* sqrt(a)->0 and b_*a->0. This is a joint local-readout
limit, not a claim of operator-norm convergence across different
volume Hilbert spaces.

At fixed I, (5) proves the time-first limit, followed by (9).
At fixed a, YM-46 constructs the volume-first stationary history
functional omega_a^S. Let I grow in the same comparison with the
fixed crop A_(D(a)); (10) then applies with N=D(a). Let a decrease.
Thus both iterated limits equal Omega.

The theorem concerns finite-volume **vacuum** histories, so the remote
time endpoints have already been removed. An arbitrary simultaneous
limit of finite rectangles with finite temporal padding is not asserted.
Spatial lattice spacing and coupling trajectory are not additional
cutoffs removed by this theorem.

## YM47-T8: positive functional, gap, reflection and exact remaining scope

The history algebra is the finite linear span of products of bounded
local row observables at finitely many distinct real heat times,
with equal times combined by multiplication. The theorem applies to
each such term. Finite positive history functionals give
|Gamma(F)|<=||F||_infinity, so the limit is well-defined, linear,
normalized and positive on this algebra, and extends to its uniform
completion. We do not silently enlarge this completion to **all**
bounded measurable functions on an infinite path space.

Adding an unused coordinate, inserting the constant one, changing the
spatial support interval, or choosing a different rounding convention
does not change the limit. Finite-volume continuous histories are
stationary and time-reversal invariant; these identities pass to Omega.
Their norms are continuous when observation times stay separated by a
fixed positive distance, by (4). Continuity at time collisions is a
separate issue.

For bounded local row F,G, the finite normalized gap and Cauchy--Schwarz
pass through the volume limit:
\[
|\Omega(\overline{F(0)}G(t))-\overline{\Omega(F)}\Omega(G)|
\leq e^{-\gamma t}
\sqrt{\operatorname{Var}_\Omega(F)\operatorname{Var}_\Omega(G)},\quad t>0.
\tag{11}
\]
Finite reflected histories factor as <f,f>, or <f,V_I(s)f> with s>=0.
Both are nonnegative. Finite linear combinations followed by (9) prove
reflection positivity of Omega about any real heat-time plane.

These statements do not identify products of limiting conditional time
maps. The theorem does not construct an infinite-volume strongly
continuous operator semigroup, its generator/domain, a unitary physical
time, or a countably additive path measure. In particular the fixed-a
maps from YM-46 are not promoted to a continuous-time map by a weak
correlation limit. That composition/completion obligation is a next step.

## Certification and reuse

YM-44's unnormalized finite-width estimate, YM-45's common vacua and gap,
and YM-46's face comparison are reused without modifying their evidence.
The new content is the normalized accumulation bound, controlled vacuum
replacement, separated-time rounding, summable continuous-parameter
spatial extensions and the cropping argument for the joint limit.
There is no claim of priority for general semigroup or thermodynamic-limit
methods.

The new certificate checks directed error budgets and tails, rational
vacuum perturbations, independent two-state normalized history readouts,
positive-contraction polynomial inequalities, asymmetric cropping
geometry, strict gates and refusal witnesses. These finite controls
support the written proof; they are not a mechanical formalization or
an SU(2) discretization proof. All four prior parameter cells are retained;
YM46-T1 supplies parameters for the entire abs(theta)<1/1680 window.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym47_joint_history.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym47_joint_history.py -v
~~~

Default/--check is read-only. --write explicitly regenerates only this
chapter's certificate. CI uses one Python 3.12 job.
