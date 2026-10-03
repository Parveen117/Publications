# YM-55: native curvature protocols retain their gap in a joint history limit

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

## Result and exact scope

Fix one site-indexed, time-independent native protocol profile on the
integer chain. Let each \(C_i\) be symmetric positive semidefinite, with
ordered eigenvalues \(0\leq a_i\leq b_i\leq c_i\). Assume
\[
 b_i\geq\beta>0,\qquad \operatorname{tr}C_i\leq M<\infty,
 \qquad |\theta|<\beta/2400.                              \tag{1}
\]
The actual square-sourced interacting chain has a **unique joint limit
of its local finite-volume vacuum history readouts**, as both spatial ends
recede and the heat step goes to zero. No relation between the rates of
these limits is required. Both iterated limits agree. The positive limit
functional has reflection positivity in heat time and exponential centered
time-correlation decay with an explicit positive rate. Rank two, differing
site orientations and spatially inhomogeneous profiles are allowed.

YM53 already supplies the finite-width interacting gap and norm time
refinement; YM54 supplies native curvature sources for its lower bound.
YM55 verifies the previously open joint-limit hypotheses for these
noncentral, possibly singular protocols. In particular, the global trace
ceiling controls refinement as the comparison window grows. Every finite
interval must use the restriction of the **same** infinite profile.

This is a stationary chain history limit, not a four-dimensional gauge
continuum or an identification of the rate with physical mass. Spatial
lattice spacing is not removed. Time-dependent interacting vacua,
state-dependent tensors and physical source/action/clock selection remain
separate. Extra Ideas [MP-2](https://github.com/Parveen117/extra-ideas/blob/bd3287e9dfef9953132e121534ea91969c1a10e8/meta-physics/mp2/THEOREM.md)
constructs a selected evolving *free* response; it is not substituted for
stationarity here.

## Native carrier and dimensionless bookkeeping clock

Reuse YM50's native reference Phi_Q and its completed compact coefficient
representation, YM51's symmetric turn construction, YM52's energy, and
YM53's interacting chain. No new Hilbert space, measure or operator engine
is postulated. Represented notation abbreviates the existing positive
recognition completion and native cut scalars.

Put
\[
 \widehat C_i=C_i/\beta,\quad\mu=\theta/\beta,\quad
 \lambda=|\mu|,\quad\Lambda=M/\beta,\quad a=\beta h,\quad s=\beta t.
\]
This is a change of units in the declared heat parameter, not physical
calibration. Necessarily \(\Lambda\geq2\). Below, time is s, \(0<a\leq1\),
and the second-eigenvalue floor is one. For a finite interval I define
\[
 L_I=\sum_{i\in I}L_{\widehat C_i},\quad
 K_I(s)=\bigotimes_{i\in I}E_s^{\widehat C_i},\quad
 W_I=\mu\sum_{\{i,i+1\}\subset I}M_{v_i},\quad
 v_i=x_i\cdot x_{i+1},\quad w_I=\lambda(|I|-1),
\]
\[
 S_{I,a}=K_I(a/2)e^{aW_I}K_I(a/2),\qquad
 T_{I,a}=e^{aW_I/2}K_I(a)e^{aW_I/2}.
\]
YM53 constructs the finite positive vacua, the norm limit \(U_I(s)\),
and \(V_I(s)=U_I(s)/\|U_I(s)\|\). Write
\(A_{I,a}=S_{I,a}/\|S_{I,a}\|\); let \(h_{I,a},h_I\) be their positive
unit vacua. A history readout is
\[
 \Gamma_I=\langle h_I,M_{F_0}V_I(s_1-s_0)M_{F_1}\cdots
 V_I(s_r-s_{r-1})M_{F_r}h_I\rangle,\quad
 0=s_0<s_1<\cdots<s_r=T.                                \tag{2}
\]
Its discrete version \(\Gamma_{I,a}\) rounds each \(s_j\) to a multiple of
a with error at most a and replaces V by powers of A. Let A be a support
interval of size k, \(p=\sum_j|\operatorname{supp}F_j|\), and
\(\mathcal M=\prod_j\|F_j\|_\infty\). Constants can be removed. Let Delta
be the smallest positive separation, or Delta=1 for one readout.
The target is YM47's local-history algebra and uniform completion, with
equal times combined by multiplication. It is not all bounded functions
on an infinite path space.

## YM55-T1: a spatial gate on the full YM53 open window

Set \(\epsilon=4/5\). Choose an integer \(J\geq1\) and \(R,\rho>1\) with
\[
 \boxed{\zeta=R\left\{\frac{2}{\epsilon J}
                 +(160+80\rho)\lambda J\right\}<1.}       \tag{3}
\]
Write \(\eta=1-\zeta\), \(g=\log R/(10J)\), \(q=e^{-ga}\).
This implies YM53's temporal gate and gives its dimensionless gap g;
the original-clock rate is \(\gamma=\beta g\).

The gate covers every \(\lambda<1/2400\). For positive lambda choose
\[
 J=5,\quad\alpha_0=\tfrac12+1200\lambda,\quad
 \rho=1+\frac{1-\alpha_0}{160\lambda J},\quad
 A_0=(1+\alpha_0)/2,\quad R=(1+A_0^{-1})/2.
\]
The braces in (3) equal \(A_0<1\), and \(\zeta=(1+A_0)/2<1\).
At lambda=0 use J=5, rho=2 and R=3/2. Thus the spatial comparison does
not shrink the old sufficient interval. At its endpoint the budget
abstains; it does not prove a gapless phase.

**Proof of the constants.** YM53-T1 supplies at every site a symmetric
nonnegative kernel with \(19/20\leq k_s^i\leq21/20\) for \(s\geq9\).
Let \(r_a=\lceil9/a\rceil\), \(\ell=Jr_a\), so \(9\leq ar_a<10\).
YM53-T2 gives temporal-endpoint cost
\[
 D_a=r_a/\epsilon+8a\lambda\ell^2,
\]
and spatial-coordinate cost \(4a\lambda\ell\). Use YM46's clipped
blocks, retaining duplicate labels, and weight
\(\rho^{-|i-i_0|}q^{|u-u_0|}\). Within a block the ratio is at most R;
to a spatial neighbour it is at most rho R. Each coordinate is removed
ell times, and its outgoing influence is at most
\[
 R\{2r_a/\epsilon+(16+8\rho)a\lambda\ell^2\}
 \leq\ell\zeta.                                         \tag{4}
\]
There are two spatial neighbours. Substitute \(a\ell<10J\) to get
(3). This uses neither centrality nor common diagonalization. The upper
trace bound is not used until T3.

## YM55-T2: admissible face comparison and the actual square ordering

Define
\[
 H_\rho=(\rho+1)/(\rho-1),\quad H_q=(1+q)/(1-q),\quad
 C_s(a)=\frac{4a\lambda\ell\rho R H_q}{\eta},\quad
 C_t(a)=\frac{R D_a H_\rho}{\eta}.                        \tag{5}
\]
For a bounded finite T-strip readout F on support D and any two
**admissible** boundary laws,
\[
 |\omega^b(F)-\omega^{b'}(F)|\leq\operatorname{osc}(F)
 \sum_{z=(i,u)\in D}\left[
 C_s(a)(\rho^{-d_l(z)}+\rho^{-d_h(z)})
 +C_t(a)(q^u+q^{N+1-u})\right].                          \tag{6}
\]
Unchanged faces are omitted. A spatial side can carry the specified bond
or no bond. The plus sign is explicit; YM47 already recorded that
typographical correction to YM46.

**Proof and singular-support repair.** The normalized skeleton overlap
after \(r_a\) steps is at least epsilon; its temporal mismatch cost is
\(r_a/\epsilon\). For a short bridge use its length. If both temporal
endpoints change, couple the two admissible laws directly as a group:
an intermediate pair can have a zero normalizer. Positive spatial tilts
preserve support. These are YM53-T2's rules with uniform constants.

Weight the block costs as in T1. With \(M_b\) labelled blocks, coupled
conditional replacement preserves both marginals and satisfies
\[
 W_{n+1}\leq(1-\ell\eta/M_b)W_n+F_\partial/M_b.
\]
Time-face forcing per coordinate is at most \(\ell R D_a w\), and
spatial-face forcing at most \(4a\lambda\ell^2\rho Rw\). Sum the geometric
weights and let n grow to prove (6). No stationary coupling or all-points
positivity of fine kernels is assumed.

For S ordering insert a midpoint using
\[
 \frac{k^i_{a/2}(x,y)k^i_{a/2}(y,z)}{k^i_a(x,z)}\,\Phi_Q(dy)       \tag{7}
\]
only for a positive denominator. A zero-denominator pair has zero
T-path weight and zero integrated numerator. Its assigned conditional
value is immaterial. The native semigroup normalizes admissible pairs.
Thus the midpoint map is positive and unital almost surely in each strip
law; no globally continuous quotient at null pairs is asserted.

The interleaved weight is exactly the product of
\(k_{a/2}(y_j,x_j)e^{aW_I(x_j)}k_{a/2}(x_j,y_{j+1})\).
Integrating x gives S; integrating internal y gives T full edges and
two half edges. This identity uses unnormalized products and remains
true at zero entries. A changed half-edge endpoint costs at most
\(1+r_a/\epsilon\): use \(r_a+1\) sites initially, of elapsed duration
\((r_a+1/2)a\geq9\), then \(r_a\)-site skeleton segments. Short blocks
again use grouped admissible comparisons. Only the fixed time-face cost
changes to \(R(D_a+1)H_\rho/\eta\); the side constant stays the same.

Remove the remote time endpoints using YM53's finite positive vacua.
Each midpoint readout uses at most two T-layer coordinates per original
site/time coordinate and does not increase oscillation. For nested
intervals \(I\subset J\) containing A, with changed faces at distance
at least n from A,
\[
 |\Gamma_{I,a}-\Gamma_{J,a}|
 \leq4\mathcal M p\,C_s(a)
             \sum_{\text{changed faces}}\rho^{-n}.       \tag{8}
\]
Conditioning larger strips gives admissible smaller-strip exteriors
almost surely; (6) therefore applies to their mixtures. Non-nested
intervals compare through their common interior. This constructs the
fixed-a volume limit of local vacuum histories independently of spatial
exhaustion. Admissible prescribed/open-side variants have the same limit.
We do not claim uniqueness for arbitrary states with arbitrary null-set
conditional definitions.

Finally \(H_q=1+2/(e^{ga}-1)\leq1+2/(ga)\). Since \(a\leq1\),
\[
 \boxed{C_s(a)\leq C_*/a,\qquad
 C_*=\frac{40\lambda J\rho R(1+2/g)}{\eta}.}               \tag{9}
\]
Certificates may replace g by a positive lower bound. This \(1/a\)
budget alone cannot justify the joint limit.

## YM55-T3: normalized anisotropic refinement and vacuum control

Choose any supplied \(c_*\geq2\sqrt{3\Lambda}\), and let \(d=g/(1+g)\).
YM53-T4 gives
\[
 \|S_{I,a}-U_I(a)\|\leq w_Ia^{3/2}e^{w_Ia}\sqrt{3\Lambda}.
\]
Both top eigenvalues are at least \(e^{-w_Ia}\), and their difference
is at most this norm error. Normalize first, then telescope contraction
factors. For all integers \(n\geq0\),
\[
 \boxed{\|A_{I,a}^{\,n}-V_I(na)\|
 \leq c_*w_I(na)\sqrt a\,e^{2w_Ia}.}                      \tag{10}
\]
There is no exponential in width times the full observation time.

Apply (10) at \(n=\lceil1/a\rceil\), so \(1\leq na\leq2\).
YM47-T2's elementary vacuum estimate
\(\|h-k\|\leq2\|A-C\|/(1-q_0)\), with the finite-width gap, yields
\[
 \|h_{I,a}-h_I\|\leq(4c_*/d)w_I\sqrt a\,e^{2w_Ia},        \tag{11}
\]
because \(1-e^{-g}\geq g/(1+g)=d\). No pointwise lower bound on the
vacuum or its Doob quotient is used. For \(w_I=0\) the vacua coincide.

YM47-T3 proves \(\|V(t+s)-V(t)\|\leq\min(1,s/t)\) for positive times
by finite positive powers. For \(a\leq\min(1,\Delta/4)\), rounding
changes a separation by at most 2a. Two vacuum replacements cost twice
(11), the rounded duration is at most T+2, and rounding costs \(4ra/\Delta\).
With \(B_*=c_*(8/d+T+2)\),
\[
 \boxed{|\Gamma_{I,a}-\Gamma_I|\leq\mathcal M
 \{B_*w_I\sqrt a\,e^{2w_Ia}+4ra/\Delta\}.}                 \tag{12}
\]
Lambda is uniform in I by (1). This supplies the upper budget missing
from a mere anisotropic gap estimate. At Lambda=3 and \(c_*=6\), (12)
recovers YM47's normalized coefficient, but the spatial gate remains
YM55's more conservative rank-two gate.

## YM55-T4: summable continuous-time spatial extension

Choose \(1<u\), \(u^2<\rho\); for example \(u=2\rho/(\rho+1)\).
Put \(a_n=u^{-2n}\) and take n large enough for \(a_n\leq\Delta/4\).
Compare a one-site extension whose changed-face distance is at least n
and whose two widths are at most k+2n+2. Apply (12) twice and (8) once
at \(a_n\). Since \(nu^{-2n}\leq1/(u^2-1)\), define
\[
 Z=\exp\{2\lambda(k+2+2/(u^2-1))\}.
\]
Then \(|\Gamma_I-\Gamma_J|\leq\mathcal M E_n\), where
\[
 E_n=2B_*\lambda Z(k+2n+2)u^{-n}
       +(8r/\Delta)u^{-2n}+4pC_*(u^2/\rho)^n.             \tag{13}
\]
Every term is summable. Define
\[
 G(q,N)=q^N/(1-q),\quad
 H(q,N)=q^N\left\{\frac{k+2+2N}{1-q}+\frac{2q}{(1-q)^2}\right\}.
\]
\[
 \boxed{\mathcal E_N=\sum_{n\geq N}E_n
 =2B_*\lambda ZH(u^{-1},N)
 +(8r/\Delta)G(u^{-2},N)+4pC_*G(u^2/\rho,N).}              \tag{14}
\]
These follow from finite geometric sums and vanishing remainders.
Let \(A_N\) pad A by N sites at each end. Two one-site extensions prove
\(\Gamma_{A_N}\) Cauchy, with tail \(2\mathcal M\mathcal E_N\); call its
limit Omega on this history. For arbitrary paddings \(n_l,n_h\), start
from \(A_N\), \(N=\min(n_l,n_h)\), and extend its longer side. At distance
n its width is still at most k+2n+2. The extra sum costs at most
\(\mathcal M\mathcal E_N\). Thus
\[
 \boxed{|\Gamma_I-\Omega|\leq3\mathcal M\mathcal E_N.}     \tag{15}
\]
No infinite product reference or exchange of operator spaces constructs
this limit: the displayed Cauchy modulus does.

## YM55-T5: unrestricted joint refinement and both iterated limits

Let D(a)>=1 be the least integer with \(\rho^{-D(a)}\leq a^2\).
It is \(O(\log(1/a))\). Crop I to \(J=I\cap A_{D(a)}\); closer existing
faces remain unchanged. Each changed face is at distance at least D(a).
By (8)--(9) the crop costs at most \(8\mathcal M p C_*a\).
The cropped width is at most k+2D(a). Put
\(N=\min(n_l,n_h,D(a))\) and \(w_*=\lambda(k+2D(a))\).
Once the separated-time and tail gates hold, (12) and (15) give
\[
 \boxed{|\Gamma_{I,a}-\Omega|\leq\mathcal M\left[
 8pC_*a+B_*w_*\sqrt a\,e^{2w_*a}
 +4ra/\Delta+3\mathcal E_N\right].}                      \tag{16}
\]
As a tends to zero and both paddings to infinity, N tends to infinity,
\(w_*\sqrt a\to0\), and \(w_*a\to0\). All terms vanish without a
constraint on relative rates. Equation (12) gives time-first at fixed I,
then (15) gives volume. At fixed a, T2 gives volume-first; let I grow
in (16) and then let a decrease. It converges to the same Omega with
N=D(a). Permitted rounding conventions and exhaustions give the same result.

These are finite-volume **vacuum** histories: remote temporal endpoints
were removed first. An unrestricted simultaneous finite-rectangle limit,
spatial-spacing limit, cutoff-dependent changing profile or operator-norm
convergence between different-volume carriers is not asserted.

## YM55-T6: positive history, reflection and surviving correlation rate

Finite positive history functionals are bounded by the uniform norm.
Their limits are therefore well-defined, linear, normalized and positive
on the local-history algebra, and extend to its uniform completion.
Unused coordinates and constant readouts do not affect Omega.
Time stationarity and reversal pass from each fixed-profile finite
vacuum history. No spatial translation or spatial reflection symmetry
is claimed for an arbitrary inhomogeneous profile.

For bounded local row readouts F,G, in the original heat clock,
\[
 \boxed{|\Omega(\overline{F(0)}G(t))-
       \overline{\Omega(F)}\Omega(G)|
 \leq e^{-\gamma t}
       \sqrt{\operatorname{Var}_\Omega(F)\operatorname{Var}_\Omega(G)},
 \quad\gamma=\frac{\beta\log R}{10J}>0.}                 \tag{17}
\]
**Proof.** Apply the finite completed-source gap and Cauchy--Schwarz to
centered \(Fh_I,Gh_I\). Variances and correlations are local history
readouts, so T4/T5 pass every term to the limit.

Reflection in a heat-time plane factors each finite future-history
quadratic form as \(\langle f,f\rangle\) or
\(\langle f,V_I(s)f\rangle\), \(s\geq0\). Both are nonnegative.
Taking finite linear combinations and the T4 limit proves reflection
positivity about every real heat-time plane. Anisotropy and varying site
orientations do not affect self-adjointness. Separated-time norm
continuity follows from the positive-power estimate in T3.

This chapter does not infer continuous-time row-Markov composition,
time-zero surjectivity, a physical unitary-time generator or a path
measure from correlation convergence. YM48/YM49's additional
reconstruction/domain obligations need auditing for this profile before
their stronger conclusions are transferred. The rate in (17) is not
renamed a four-dimensional Yang--Mills mass.

## YM55-T7: native curvature sources and certified anisotropic examples

Use YM54's selected response construction at every site: two shape
directions, its compact lift and four equally counted signed turns.
If \(|f_i|\geq f_0>0\) and \(\tau_i\leq T_0\), YM54 supplies both
\(b_i\geq4f_0^2/T_0\) and \(\operatorname{tr}C_i=\tau_i/4\leq T_0/4\).
Therefore choose
\[
 \boxed{\beta=4f_0^2/T_0,\quad M=T_0/4,\quad
 \Lambda=T_0^2/(16f_0^2),\quad
 |\theta|<f_0^2/(600T_0).}                               \tag{18}
\]
The same native curvature floor and response budget now control both
interaction mixing and joint-limit refinement. The profile must be
consistent across intervals and frozen in heat time. Fixed sitewise native
compact reframings are permitted, with their selections declared.

For \(C_i=O_i\operatorname{diag}(0,3/2,3/2)O_i^T\), with any fixed
sitewise native compact rotations, choose beta=3/2, M=3, theta=1/4096.
Then lambda=1/6144, Lambda=2, J=8, R=4/3, rho=3/2 give
\[
 \eta=7/72,\qquad\gamma=3\log(4/3)/160>0.00539403.
\]
The earlier finite-width example now has (16) with that rate.
Take u=9/8 and \(c_*=5\), since \(5^2\geq12\Lambda\).
Outward rational certificates give crop radii and remaining errors.
These conservative budgets do not predict experimental lengths or masses.

YM54's original curved thermo fixture also passes:
\[
 \beta=5/234,\quad M=65/648,\quad\Lambda=169/36,\quad
 \theta=1/262144,\quad J=6,\quad R=5/4,\quad\rho=3/2,
\]
with \(\eta=10249/98304\) and \(\gamma>0.000079467076\).
Its actual tensor and fixed compact rotations give an admissible infinite
profile. An arbitrary unconstrained tensor assignment does not thereby
acquire a polynomial thermodynamic source.

## Failure controls, provenance and replay

- A uniform finite-width gap alone does not prove state convergence;
  (6), (12) and (16) supply the required comparisons.
- Fine kernels can have zero normalizers. Exact controls check midpoint
  identities with actual zeros, arbitrary null-pair fillers and grouped
  endpoint comparisons; division by zero is refused.
- Bounds (1) do not permit a fixed site's source to alternate with the
  cutoff. Free linear readouts already distinguish two such profiles and
  give different subsequential correlations.
- Rank one fails the beta gate. The trace budget and curvature floor in
  (18) remain required. Refusal of this certificate does not establish
  nonexistence of another limit or a gapless phase.
- Rescaling \(C,\theta,\beta,M\) by k and \(t,h\) by \(1/k\) leaves
  a, lambda, Lambda, cutoffs and margins unchanged; gamma scales by k.
  Physical calibration and four-dimensional gauge dynamics remain open.

YM46/YM47's comparison, normalization and logarithmic-cropping methods
are reused, not claimed as new general techniques. The new content is
their audited anisotropic application with singular admissible supports,
uniform trace control, fixed-profile consistency and native curvature
sources. General block-comparison lineage: Rebeschini and van Handel,
[Comparison Theorems for Gibbs Measures](https://web.math.princeton.edu/~rvan/dobrushin130819.pdf).
The external paper supplies attribution, not a substitute proof premise.

Seven general results have written proofs and reproducible exact/outward
finite controls. These are not proof-assistant verification or independent
expert review. No native engine or prior certificate is rewritten.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym55_anisotropic_joint_limit.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym55_anisotropic_joint_limit.py -v
~~~

Default/--check is read-only; --write regenerates only this chapter's
result and expected digest. Pins bind the proof, code, tests and workflow.

