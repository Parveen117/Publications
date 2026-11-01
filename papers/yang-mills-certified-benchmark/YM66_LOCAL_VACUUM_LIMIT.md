# YM-66: boundary response, the local vacuum limit and its spectral tower

Monty Dabas. 4 October 2026. **Python 3.12 only.**

YM65 gives spatial decay at a fixed coupling but leaves the vacuum
history limit subsequential. This continuation extends the gap and
response estimates along independent plaquette-coupling paths. Integrating
the response makes distant changes of the Hamiltonian small in the
local state norm. Consequently all admissible exhaustions give the same
local vacuum, real-time correlations and heat histories. Every fixed
inverse-response moment now has a unique thermodynamic limit.

For each fixed finite set of coarse links, its actual vacuum density
also converges in L1, with an explicit spatial-buffer error. This is
control of a marginal density under distant Hamiltonian changes.
Pointwise conditional-density estimates, the full Schur memory kernel,
and the four-dimensional spatial continuum remain separate problems.

## Setup and relation to established methods

On a finite open cubic box or periodic cube of side at least four, use
the same full native link Hilbert spaces and electric operator as YM62.
Allow independent real face couplings:

\[
 A(j)=L_E-\sum_p j_p W_p,\qquad |j_p|\le\eta<1/8000.       \tag{1}
\]

A missing face means j_p=0. We retain every link's free Casimir when
deleting interactions. Write h_j for the positive normalized vacuum,
E_j for its energy, H_j=A(j)-E_j, and omega_j(B)=<h_j,Bh_j>.
The local algebra on a finite link set X is the entire B(H_X), with
its usual tensor inclusion in larger volumes. Gauge-invariant local
operators and bounded multiplication readouts are included.

Use the link graph and support distance from YM65, formed from **all**
geometrically allowed plaquettes, including faces whose coupling is
zero. Thus the distance is common along any path in (1). Put

\[
 \gamma_*=\gamma(\eta),\quad v_*=64\eta,\quad
 \rho_*=\frac{(\log2)\gamma_*}{2v_*+\gamma_*}>0.             \tag{2}
\]

The principle that local perturbations have local effects in a gapped
vacuum is established; see [De Roeck and Schuetz,
arXiv:1501.04571](https://arxiv.org/abs/1501.04571).
For infinite-dimensional sites and unbounded on-site Hamiltonians,
the strong-operator and infinite-volume dynamics framework is developed
by [Nachtergaele, Sims and Young, arXiv:1810.02428, Sections 3 and
5](https://arxiv.org/abs/1810.02428). Here the needed estimates are
derived from YM62 and YM65 with explicit native constants. No priority
claim is made for the general locality principle or compact moment
determinacy.

## YM66-T1: the whole bounded plaquette-coupling cube is uniformly gapped

Every Hamiltonian (1), including every straight segment between two
profiles in the cube, obeys

\[
 H_j\ge\gamma_*(I-|h_j\rangle\langle h_j|).                 \tag{3}
\]

To verify the extension, revisit YM62-T2 before specializing its
coefficients. For a temporal block on link e, its log interaction tilt
is bounded by hk sum_(p containing e)|j_p| <=4hk eta. An endpoint
change therefore has cost r/epsilon+16h eta k^2. Changing another link
at one time changes the log tilt by at most
2h sum_(p containing both links)|j_p|. Summing the shared-face
multiplicities, including repeated neighbours, gives exactly the
previous temporal contraction budget

\[
 \bar\alpha=\tfrac12+4000\eta<1,\qquad
 \gamma_*=\frac1{50}\log\frac{1+\bar\alpha^{-1}}2.          \tag{4}
\]

No step used equality or a common sign of the coefficients. The
finite-volume potential bound is now b=sum_p |j_p|. For time refinement,
each face still has electric rate three, so the triangle inequality
gives ||[K_s,M_(sum j_p W_p)]||<=b sqrt(6s). The fixed-volume norm
refinement and isolated-vacuum passage of YM62-T3 therefore apply with
this b. Positivity improvement and gauge invariance are unchanged.
This proves (3) on the full space and in its physical sector.

For propagation, the weighted link-incidence matrix has entries
sum_(p containing e,f)|j_p| and row sums at most 16 eta, including
the diagonal. The YM65 Gaussian-filter proof thus applies everywhere
in the cube with the common constants (2). In particular, for
self-adjoint B supported on X and a separated face p at distance d>=1,

\[
 |\partial_{j_p}\omega_j(B)|
 =|2\operatorname{Re}\langle Q_jBh_j,H_{j,\perp}^{-1}
                                      Q_jW_ph_j\rangle|
 \le\|B\|K_*(d),\quad
 K_*(d)=\min\left\{\frac2{\gamma_*},
        \frac{2(22\rho_*d+5)}{3\gamma_*}e^{-\rho_*d}\right\}. \tag{5}
\]

The derivative follows by differentiating the simple isolated
eigenvector, with <h_j,partial h_j>=0. Since the perturbation is
-W_p, its sign in (5) is positive before taking the absolute value.
The source response need not have a positive off-diagonal sign.
At overlapping support one may use 2||B||/gamma_*.

## YM66-T2: integrating the response controls distant changes

Let j,j' be profiles on a common finite geometric graph. Write
delta_p=j'_p-j_p and delta=sup_p |delta_p|<=2 eta. Suppose every
changed face has distance at least R>=1 from X. The straight segment
stays inside (1), so integrating (5) proves

\[
 |\omega_{j'}(B)-\omega_j(B)|
 \le\|B\|\sum_{p:\delta_p\ne0}|\delta_p|K_*(d(X,p)).       \tag{6}
\]

For an arbitrary finite link set X, the number of faces at distance
d is at most 3|X|(2d+5)^3. Choose a nearest link in X and a link on
the face. Each graph step moves the midpoint by at most one; allowing
the face's base-point offset puts that base inside a coordinate cube
of radius d+2 about the chosen link. There are three face orientations.
Counting separately for each link in X overcounts and is legitimate.
Wrapped representatives give the same upper bound on a torus.

Consequently, with a normalized observable, a convenient common bound is

\[
 \begin{split}
 \mathcal E_X(R;\delta)
 &=\min\left\{2,\frac{2|X|\delta}{\gamma_*}
     \sum_{d=R}^{\infty}(2d+5)^3(22\rho_*d+5)e^{-\rho_*d}\right\},\\
 \boxed{\|\varrho_{j',X}-\varrho_{j,X}\|_1
                       \le\mathcal E_X(R;\delta).}         \tag{7}
 \end{split}
\]

Here varrho_(j,X) is the reduced quantum density operator on H_X.
Trace-norm duality uses only self-adjoint norm-one observables; hence
(6) already gives the second line. No finite local Hilbert dimension
is assumed. The bound is O(|X| delta (1+R^4)exp(-rho_* R)).

For a reproducible tail calculation put q=exp(-rho_*) and
T_k(R,q)=sum_(d>=R) d^k q^d. Then

\[
 \begin{split}
 T_k(R,q)&=q^R\sum_{a=0}^k {k\choose a}R^{k-a}S_a^{(0)}(q),\\
 S_0^{(0)}&=(1-q)^{-1},\qquad
 S_a^{(0)}=\sum_{m\ge1}m^a q^m\quad(a\ge1).              \tag{8}
 \end{split}
\]

Use YM65's rational formulas for a=1,...,4. The coefficients multiplying
T_0,...,T_4 in (7) are respectively
625, 2750rho_*+750, 3300rho_*+300, 1320rho_*+40, 176rho_*.
In particular the R endpoint is included. Outward interval rounding
is used in the certificate.

For eta=delta=1/16000 and |X|=4, upward-rounded trace-norm bounds are:

| Buffer R in the link graph | Bound from (7) |
|---:|---:|
| 64 | 2 |
| 128 | 0.200366 |
| 192 | 0.000004121 |
| 256 | 0.00000000005490 |

These are conservative error budgets for distant coupling removal,
not numerical solutions of the vacuum.

## YM66-T3: all admissible exhaustions select one local vacuum

Fix an infinite profile on the elementary faces of Z3 with
sup_p |j_p|<=eta. An admissible finite-volume exhaustion agrees with
that profile on every fixed interior cube eventually, while its
remaining face couplings are arbitrary in [-eta,eta]. Open cubes
and centered periodic cubes, with any such distant couplings, qualify.
The observation region must recede from the geometric boundary;
no statement is made for a probe kept at an edge.

Choose an interior cube C containing a finite support X. In any larger
volume delete every face not entirely supported on the links E(C).
The Hamiltonian becomes A_C tensor identity plus free outside Casimirs.
Its vacuum factors as h_C tensor 1. The path to this decoupled
Hamiltonian stays in (1), and the changed coefficients have magnitude
at most eta. If their distance from X is at least R, (7) gives

\[
 \|\varrho_{\Lambda,X}-\varrho_{C,X}\|_1
                         \le\mathcal E_X(R;\eta).          \tag{9}
\]

Any two sufficiently large volumes can be compared to this same cube.
Thus their local density operators are trace-norm Cauchy. Trace-class
completeness supplies a positive trace-one limit varrho_(infinity,X).
Partial trace is contractive, so the limits are compatible as X varies.
They define one locally normal state omega_j on the quasi-local
observable algebra. Equation (9) also bounds the difference between
the interior cube and the infinite-volume limit.

This closes the ambiguity between admissible finite-volume vacuum
limits, for **all bounded local quantum observables**, not only the
plaquette mean. It is not a classification of every algebraic ground
state or of arbitrary imposed boundary states outside (1).
For a uniform profile, translation, lattice, and local gauge symmetries
pass from the periodic sequence to this unique limiting state.

## YM66-T4: real-time locality removes the history subsequence as well

For a real time t, bounded perturbation Duhamel comparison and YM65-T1
give, under the separation assumptions of T2,

\[
 \begin{split}
 \|\tau_t^{j'}(B)-\tau_t^j(B)\|
 &\le\|B\|\mathcal D_X(R,\delta,t),\\
 \mathcal D_X&=\min\left\{2,\,
 24|X|\delta\,J_{v_*}(|t|)
       \sum_{d=R}^{\infty}(2d+5)^3 2^{-d}\right\},\\
 J_v(t)&=\frac{e^{vt}-1}{v}-t\quad(v>0),\qquad J_0(t)=0.
                                                               \tag{10}
 \end{split}
\]

Indeed each commutator has S<=4 and bound
8||B||2^(-d)(exp(v_*u)-1), and integrating u from zero to |t|
gives J. The geometric shell count gives the factor twenty-four.
For eta=0 there is no interaction difference and the comparison is zero.
Strong integrals suffice for the unbounded on-link Casimirs; the norm
estimate concerns the bounded difference of the evolutions.

Comparing to the same interior cube proves a norm limit tau_t(B) for
each bounded local B, uniform on compact real-time intervals. Isometry
extends it to the quasi-local algebra, and the finite-volume composition
law passes to the limit. We do not assert pointwise norm continuity in
t on all B(H_X): unbounded on-site generators need not have that property.

For local B_0,...,B_k, the vacuum word

\[
 F_\Lambda(t_1,\ldots,t_k)=
 \langle h_\Lambda,B_0e^{it_1H_\Lambda}B_1\cdots
                                  e^{it_kH_\Lambda}B_kh_\Lambda\rangle
                                                               \tag{11}
\]

equals omega_Lambda(B_0 tau_(t_1)(B_1) ...
tau_(t_1+...+t_k)(B_k)). Approximate each evolved observable by
its evolution in a large fixed interior cube. Equation (10) controls
the norm errors uniformly in the outer volume, and T3 controls the
remaining local word. Thus (11) has one limit at every real time tuple.
Its magnitude is at most product_i ||B_i||.

Since H_Lambda>=0, the spectral theorem and the scalar Poisson identity give

\[
 e^{-uH_\Lambda}=\int_{\mathbb R}
       \frac{u}{\pi(u^2+t^2)}\,e^{itH_\Lambda}\,dt,\qquad u>0.
                                                               \tag{12}
\]

The Fourier transform of this probability kernel is exp(-u|E|);
nonnegative **vacuum-shifted** H is essential. Use (12) in every
positive-time factor in a heat word. The product kernel has mass one,
so the common norm bound and dominated convergence transfer the unique
limit of (11) to every heat word. Zero time increments simply combine
adjacent operators. For a finite time cutoff T,
the scalar omitted mass is at most 2u/(pi T), an explicit tail control.

For bounded continuous configuration readouts these words determine
the finite-dimensional history laws of YM62. All periodic subsequences
therefore coincide; the entire sequence converges to one history law.
Positivity, reflection positivity, the local time-continuity estimate
and YM62's reflected gap inequality pass to that law. Its reconstructed
physical vacuum representation retains gap at least gamma_* and, for
a uniform coupling, the plaquette variance greater than 19/100.
No independent spatial continuum or time-zero Markov identification
is supplied by this uniqueness argument.

## YM66-T5: fixed coarse vacuum marginals converge in L1

Let pi_b multiply each of finitely many disjoint fine-link paths into
a coarse link, as in YM63 and YM64. Fix a finite set Y of coarse links
and let X be the union of their **complete fine paths**. For a fixed
integer b this is a finite set, with |X|=b|Y| for the usual disjoint
straight paths.

The finite-volume configuration marginal on X has a density relative
to the product reference. Its classical total variation is bounded by
half the quantum trace distance in (7). Pushforward by pi_b contracts
total variation. Product reference pushes forward to coarse product
reference, by invariance of the native reference under multiplication.
Writing w_(Lambda,b,Y) for the actual coarse marginal density gives

\[
 \boxed{\|w_{\Lambda,b,Y}-w_{C,b,Y}\|_{L^1(\Phi_Y)}
                  \le\mathcal E_X(R;\eta).}               \tag{13}
\]

The densities are L1 Cauchy and converge to one probability density
w_(infinity,b,Y). Marginalization and any fixed finite composition of
product-link cuts commute with the limit. This is an actual coarse
vacuum marginal, with an error that depends on its fine support.
It does not assume that its density is a one-parameter plaquette ansatz.

There is also controlled conditioning on events of positive probability.
For probability measures mu,nu with TV(mu,nu)<=epsilon and an event C
with both probabilities at least p_0>0, elementary division yields

\[
 \operatorname{TV}(\mu(\,\cdot\,|C),\nu(\,\cdot\,|C))
                        \le\min\{1,2\epsilon/p_0\}.        \tag{14}
\]

For a fixed finite coarse event of positive limiting probability,
(13) therefore controls the corresponding conditional law at all large
volumes. The support X must include every path used by the conditioned
event as well as by the readout.

The factor 1/p_0 matters. YM65's rare-flag comparison has
TV=p/4 before conditioning and TV=1/4 after conditioning on a flag
of probability p. No bound uniform in vanishing p follows. Exact
continuous coarse configurations have probability zero and require
separate pointwise conditional control. L1 convergence does not establish
smooth positive density bounds or convergence of its logarithmic
derivatives. Therefore it does not yet control YM64's drift or full
Schur memory kernel in a spatial operator norm.

## YM66-T6: one infinite-volume inverse tower and its spectral content

For bounded local B,C, the unique heat-word limits and the common gap
bound justify dominated convergence in every fixed-order integral:

\[
 \begin{split}
 R_{\infty,n}(B,C)
 &=\lim_\Lambda\langle Q_\Lambda Bh_\Lambda,
               H_{\Lambda,\perp}^{-n}Q_\Lambda Ch_\Lambda\rangle\\
 &=\frac1{(n-1)!}\int_0^\infty
       u^{n-1}C_{\infty,BC}(u)\,du,\qquad n\ge1.          \tag{15}
 \end{split}
\]

Order zero is the unique equal-time covariance. YM65's spatial bounds
and its right-half-plane resolvent bounds pass to these limits. Thus
the response tower is attached to one state, rather than to an
unspecified history subsequence. Here n is an inverse spectral order;
it is not automatically an n-th nonlinear source derivative.

For a real smooth local configuration readout B with nonzero centered variance, let
sigma_B be its positive excitation spectral measure in the reconstructed
vacuum representation. Its support lies in [gamma_*,infinity). Then

\[
 m_n=R_{\infty,n}(B,B)=\int E^{-n}\,d\sigma_B(E),\quad n\ge0.
                                                               \tag{16}
\]

Push the measure to x=1/E in [0,1/gamma_*]. All its moments are m_n.
Polynomials are dense in continuous functions on this compact interval,
so the complete sequence, including m_0, determines this measure
uniquely. No finite list generally suffices. The native energy identity
gives (m_0-C_BB(u))/u<=||Gamma(B)||_infinity at finite volume.
Pass to the unique limit and let u decrease to zero to obtain
integral E d sigma_B(E)<=||Gamma(B)||_infinity. Hence the
mass with 0<x<=epsilon is at most epsilon ||Gamma(B)||_infinity.
There is no hidden atom at x=0 in this pushforward.

Positivity gives the block Hankel inequalities and the gap localizer:
for all finite scalar sequences c_i,

\[
 \sum_{i,j}\bar c_i c_j m_{i+j}\ge0,\qquad
 \sum_{i,j}\bar c_i c_j
       (\gamma_*^{-1}m_{i+j}-m_{i+j+1})\ge0.              \tag{17}
\]

In particular m_n m_(n+2)>=m_(n+1)^2. If E_B is the infimum of the
spectral support seen by B, compact-moment concentration gives

\[
 \gamma_*\le E_B\le\frac{m_{n+1}}{m_{n+2}}
                    \le\frac{m_n}{m_{n+1}},\qquad
 \lim_{n\to\infty}\frac{m_n}{m_{n+1}}=E_B.               \tag{18}
\]

To see the limit, the ratios m_(n+1)/m_n increase by log convexity
and are bounded by the largest supported x. Any smaller candidate
contradicts the positive mass arbitrarily close to that endpoint,
or equivalently the limit of m_n^(1/n). Inverting gives (18).

Thus the infinite tower determines the excitation spectrum **seen by
the chosen readout**. A readout can be blind to a lower excitation,
so E_B need not equal the full Hamiltonian gap. The certificate gives
an explicit two-atom reconstruction and an invisible lower-level
counterexample. It does not compute the interacting YM spectral
measure or a physical mass. A clock rescaling multiplies E_B and all
ratios in (18) by that clock factor; it leaves the spatial rates unchanged.

## Verification and remaining scope

The Python 3.12 controls check signed coupling profiles and path budgets,
native spatial cuts, shell and tail counts, a source-derivative sign,
real-time word identities, the need to shift the Poisson spectrum,
coarse pushforward contraction, rare conditioning, and positive moment
reconstruction. The finite spectral fixtures are labelled comparisons;
none is a finite-spin replacement for the actual link vacuum.
The analytic proofs above control the full Hilbert spaces. Arithmetic
certificates are reproducibility checks, not formal verification or review.

The completed step is uniqueness and quantitative local approximation
of the admissible lattice vacuum and its histories, together with L1
control of every fixed coarse marginal. Pointwise conditional-vacuum
control, spatial Schur-kernel estimates, a time-zero Markov identification,
and a nontrivial continuum trajectory remain open.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym66_local_vacuum_limit.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym66_local_vacuum_limit.py -v
~~~
