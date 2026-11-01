# YM-65: spatial decay of connected responses and the bulk plaquette limit

Monty Dabas. 4 October 2026. **Python 3.12 only.**

YM64 controls hidden inverses and nested blocks relative to coarse energy.
This continuation adds a spatial estimate for the actual interacting
gauge Hamiltonian. For separated bounded readouts, connected correlations
and every fixed inverse-response order decay exponentially in spatial
distance, with constants independent of volume. Applied to independent
plaquette sources, the estimate makes the susceptibility absolutely
summable. It also gives a unique bulk plaquette mean and a quantitative
finite-volume approximation bound for that mean.

The estimates hold at the full coupling theta inside YM62's window;
they are not truncations of a theta expansion. Their relation to the
remaining conditional-vacuum and Schur-kernel locality problem is stated
in T4 and T6. No four-dimensional continuum construction is asserted.

## Model, distance and established methods

Retain A_L=L_E-theta sum_p W_p on finite open cubic boxes or periodic
cubes of side L>=4. The link Hilbert spaces are the full native compact
completions. Write H=A_L-E_L>=0, with normalized vacuum h and
Q_h=I-|h><h|. YM62 gives

\[
 H\ge\gamma(\theta)Q_h,\quad
 \gamma(\theta)=\frac1{50}\log\!\frac{1+\alpha^{-1}}2,
 \quad\alpha=\frac12+4000|\theta|<1.                       \tag{1}
\]

All uniform-gap statements below use |theta|<1/8000. Propagation alone
does not require that restriction. One native unit of heat time is used
until the clock rescaling in T6.

The spatial graph has one vertex for each oriented positive lattice link.
Two distinct links are adjacent if an elementary plaquette contains both.
Let d(X,Y) be its shortest-path distance between link supports, and set
S=min(|X|,|Y|). Throughout the spatial estimates X,Y are nonempty and
disjoint, so d>=1. Norms are operator norms. Gauge-invariant Wilson
readouts are included; their supports are their complete fine boundaries.

The tools are established Lieb--Robinson propagation and Gaussian-filter
clustering. [Nachtergaele and Sims, arXiv:1410.8174,
Section 3](https://arxiv.org/abs/1410.8174) treats the strong-operator
calculus needed for unbounded on-site Hamiltonians. Their earlier
[arXiv:math-ph/0506030](https://arxiv.org/abs/math-ph/0506030) develops
the gap-to-clustering argument. Here these methods are applied with the
explicit native plaquette incidence and YM62 gap. We give the specialized
constants, the inverse-response integral, the blocked-resolvent meaning,
and the bulk approximation consequence. No priority claim is made for
the general propagation or clustering mechanism.

## YM65-T1: propagation depends on plaquette incidence, not a spin cutoff

Let N_(e,f) count plaquettes containing both e and f, including e=f.
It is a nonnegative symmetric matrix. Each plaquette has four links,
and each link belongs to at most four plaquettes. Consequently

\[
 \sum_f N_{ef}\le16.                                      \tag{2}
\]

On periodic cubes the sum is exactly sixteen: four diagonal incidences
and twelve off-diagonal incidences. The diagonal is part of this count.
For link midpoints in lattice units, every adjacent pair has L1 distance
at most one. The same holds for wrapped distances on a torus. Thus
d(X,Y) is at least the minimum midpoint L1 separation of their supports.

The unbounded part L_E is a sum of independent on-link Casimirs. Its
unitary evolution preserves every support. In the interaction picture,
each plaquette term is bounded, remains supported on its four links,
and has norm at most |theta|. The interaction-picture functions are
strongly continuous on the separable Hilbert space. Bounded strong
Dyson integrals and their norm estimates are legitimate at fixed volume;
no norm differentiability of the free on-link automorphism is assumed.

Removing terms internal to a support before differentiating its
commutator gives the usual boundary-plaquette integral inequality.
Iteration produces chains of overlapping plaquettes connecting X to Y.
A chain with n plaquettes selects an n-step walk in N by choosing a
link at each overlap and at the two ends. The walk sum overcounts the
chains. It vanishes for n<d, and (2) bounds it by |X|16^n. Swapping
the observables and reversing time replaces |X| by S. Therefore, for
tau_t(A)=exp(itH) A exp(-itH),

\[
 \|[\tau_t(A_X),B_Y]\|
 \le2S\|A\|\|B\|\sum_{n\ge d}
                      \frac{(32|\theta||t|)^n}{n!}
 \le\boxed{2S\|A\|\|B\|\,2^{-d}(e^{v|t|}-1)},\quad
 v=64|\theta|.                                            \tag{3}
\]

The bound can always be clipped at 2||A|| ||B||. The second inequality
uses 2^n>=2^d and n>=1. Constants are independent of volume and of
the unbounded Casimir spectrum. Restricting gauge-invariant operators
to the physical subspace preserves the bound: it was proved on the
tensor-product space before restriction, so no tensor factorization of
the constrained physical Hilbert space is presumed.

## YM65-T2: a space-and-heat bound for the full interacting vacuum

For bounded self-adjoint local A,B define
C_AB(u)=<A h,exp(-uH)Q_h B h>, u>=0. The same norm estimate holds for
complex readouts by replacing the first operator in the commutator by
its adjoint. Introduce

\[
 \rho=\frac{(\log2)\gamma}{2v+\gamma}>0,\qquad
 r=\rho d,\quad T=\frac{2r}{\gamma},\quad
 C_S(r)=1+\frac{2S}{3}+\frac1{3r}.                         \tag{4}
\]

Then

\[
 \boxed{|C_{AB}(u)|\le\|A\|\|B\| C_S(r)
       e^{-r-\gamma u^2/(2T)}\quad(0\le u\le T),}
                                                               \tag{5}
\]
and for every u>=0 the independent gap bound is
|C_AB(u)|<=||A|| ||B|| exp(-gamma u).

Here is a direct scalar-filter proof, retaining the normalization that
fixes the heat sign. Put f(t)=<h,[tau_t(A),B]h> and for a>0 define

\[
 J_a(u)=-\frac{e^{-au^2}}{2\pi i}
             \int_{\mathbb R}\frac{e^{-at^2}f(t)}{t+iu}\,dt.
                                                               \tag{6}
\]

At u=0 use the limit from u>0. Since the supports are disjoint, f(0)=0
and (3) bounds f(t)/t near zero, so this limit is controlled. For a
positive spectral energy E, the negative-frequency term exp(-itE)
in f maps to
exp(-uE) Phi_a(E-2au), where Phi_a(x)=erfc(-x/(2 sqrt(a)))/2.
The positive-frequency term maps to exp(uE)Phi_a(-E-2au).
These identities follow by expressing 1/(t+iu) as a one-sided Fourier
integral and completing the Gaussian square. Each of the two spectral
measures has total variation at most ||A|| ||B||. Vacuum poles cancel.

For E>=gamma and 2au<=gamma, erfc(x)<=exp(-x^2) for x>=0 yields

\[
 |C_{AB}(u)-J_a(u)|
       \le\|A\|\|B\|e^{-\gamma^2/(4a)-au^2}.              \tag{7}
\]

In particular, the two completed-square exponents are both
-E^2/(4a)-au^2. The negative Gaussian normalization in (6) is essential.
For the short-time part of (6), use (3), |t+iu|>=|t| and
integral_0^T (exp(vt)-1)/t dt <= exp(vT)-1. For its tail use the
universal commutator bound and
integral_T^infinity exp(-at^2)/t dt <= exp(-aT^2)/(2aT^2).
Together with (7) this gives

\[
 \frac{|C_{AB}(u)|}{\|A\|\|B\|}
 \le e^{-au^2}\left[
 e^{-\gamma^2/(4a)}+
 \frac{2S}{\pi}e^{-(\log2)d}(e^{vT}-1)+
 \frac{e^{-aT^2}}{\pi aT^2}\right].                       \tag{8}
\]

Choose a=gamma/(2T). Then aT^2=gamma T/2=r and
(log2)d-vT=r; the three exponents agree. Finally 1/pi<1/3 gives
(5). Only the earlier gap and bounded propagation were used. This is
a full-coupling estimate in the stated window, not a low-order vacuum
ansatz. At theta=0 the product vacuum gives C_AB(u)=0 exactly for
disjoint supports, a stronger separate statement.

For theta=1/16000, gamma=log(7/6)/50 and certified intervals give
rho>0.1928159848815. The following are upward-rounded bounds for two
plaquette readouts of norm at most one (S=4):

| Link-graph distance | Equal-time covariance | gamma times inverse order 1 | gamma^2 times inverse order 2 |
|---:|---:|---:|---:|
| 16 | 0.172609 | 1 | 1 |
| 32 | 0.007780 | 0.096011 | 0.592430 |
| 64 | 0.000016150 | 0.000398586 | 0.004918633 |
| 128 | 7.036e-11 | 3.473e-9 | 8.571e-8 |

The last two columns use T3. These are conservative upper bounds,
not measured correlations or predicted masses.

## YM65-T3: every fixed inverse order has spatial decay

Let R_n(A,B)=<Q_h A h,(H|_(h-perp))^(-n)Q_h B h> for n>=1.
Spectral calculus and the gap give the convergent integral

\[
 R_n(A,B)=\frac1{(n-1)!}\int_0^\infty u^{n-1}C_{AB}(u)\,du.
                                                               \tag{9}
\]

Use (5) up to T, dropping its additional Gaussian time decay, and
use the gap bound beyond T. Integration of the latter is an elementary
integer incomplete-gamma identity. The result is

\[
 \boxed{|R_n(A,B)|\le\frac{\|A\|\|B\|}{\gamma^n}
 \min\left\{1,
 C_S(r)e^{-r}\frac{(2r)^n}{n!}
       +e^{-2r}\sum_{k=0}^{n-1}\frac{(2r)^k}{k!}\right\}.   \tag{10}
\]

For each fixed n this decays exponentially times a polynomial in d.
The bound holds in arbitrary admissible volumes. It does not assert
one prefactor uniform over n, growing supports or changing couplings.
For n=0 use the equal-time bound (5), clipped at ||A|| ||B||.

The same integral argument applies to (H+s)^(-n) for Re s>=0 because
|exp(-su)|<=1. Thus it controls every fixed spectral derivative of
the connected resolvent in the closed right half-plane. Inside |s|<gamma
the usual inverse-response expansion is analytic on the centered space.
For the native continuous local readouts in YM62's subsequential history
limits, local heat correlations converge and their common gap bound
dominates (9). Hence the spatial estimate passes to those constructed
infinite-volume responses. No uniqueness of the whole limiting state
is inferred from this passage.

There is a direct physical source interpretation. Perturb the actual
finite-volume Hamiltonian by -jA-kB. Simplicity of its vacuum gives

\[
 -\partial_j\partial_k E(j,k)|_{0}=2\operatorname{Re}R_1(A,B).
                                                               \tag{11}
\]

Differentiate the eigenvalue equation with the vacuum derivative
orthogonal to h to obtain (11). This is the same mixed-potential
reciprocity used in YM61, now with a nonperturbative spatial bound for
the YM62 gauge model. The Hessian response matrix is symmetric and
positive semidefinite, though individual off-diagonal entries may be
negative. No electromagnetic field equation is being substituted here.

## YM65-T4: exact blocked resolvents inherit these readout bounds

Use YM64's actual vacuum measure, product-link isometry I_b and
spectral Schur operator F_b(s). For centered bounded coarse readouts
f,g, the exact identity is

\[
 \langle f,(F_b(s)+s)^{-1}g\rangle_{\mu_b}
   =\langle hI_bf,(H+s)^{-1}hI_bg\rangle_{\Phi}.           \tag{12}
\]

Apply (10) to the complete fine supports of I_bf and I_bg. For n=k+1,
it bounds the k-th spectral derivative of (12), including the factor k!.
At s=0 the inverse means the centered inverse. These estimates hold
through any finite sequence of blocks, since YM64 proves their exact
resolvent composition. They are nonperturbative locality estimates for
blocked **resolvent readouts**, rather than for a fitted scalar coupling.

A coarse link uses b disjoint fine links, so S in (10) is computed on
fine supports and can grow with b. No constant independent of the size
of a growing readout is claimed. These moments are also different from
YM64's Schur-correction moments M_(b,n)=T_b* D_b^(-n) T_b.
Taking an unbounded inverse of a spatially controlled resolvent does
not automatically give locality of its generator or its correction.

Ordinary clustering alone does not control every conditional measure.
For an explicit full-support probability comparison, put binary X and
a binary flag Z near the left site, and binary Y at distance d. Set
P(Z=1)=p=2^(-d). Conditional on Z=0, X,Y are independent fair signs;
conditional on Z=1, their probability is (1+XY/2)/4. Every bounded
left/right covariance has magnitude at most p||f|| ||g||/2, since the
total-variation distance from the product of the marginals is p/4.
But E[XY|Z=1]=1/2 at every d. Conditioning on the rare flag amplifies
the small unconditional correlation by 1/p.

This comparison is not the YM Hamiltonian. It isolates the missing
implication: even spatial clustering of all bounded readouts between
two regions does not, by itself, bound rare conditional densities.
The uniform locality of the actual coarse density and of the full
Schur memory kernel remains a separate obligation.

## YM65-T5: summable susceptibility and a controlled bulk plaquette mean

Let chi_(p,q)=2 Re R_1(W_p,W_q). Plaquettes have norm at most one.
At most sixteen plaquettes share a link with a fixed p, including p.
Use |chi_(p,q)|<=2/gamma for these distance-zero cases. For distance
d>=1, (10) with S=4 and exp(-2r)<=exp(-r) gives

\[
 |\chi_{pq}|\le\frac{2}{3\gamma}(22\rho d+5)e^{-\rho d}.
                                                               \tag{13}
\]

The number of q at link-graph distance d from p is at most
3(2d+5)^3. Indeed a shortest graph path moves link midpoints by at
most d in L1; accounting for the two plaquette offsets places the
base of q in a coordinate cube of radius d+2, with three orientations.
This deliberately loose shell count also bounds periodic cubes by
choosing wrapped representatives. Exponential decay is summable against
this polynomial count.

For q0=exp(-rho), put S_k(q0)=sum_(d>=1) d^k q0^d. Then a concrete
volume-independent susceptibility budget is

\[
 \begin{split}
 \sup_p\sum_q|\chi_{pq}|\le\mathcal B(\theta)
 =\frac2\gamma\big[16+176\rho S_4
 +(1320\rho+40)S_3+(3300\rho+300)S_2\\
 \hspace{35mm}+(2750\rho+750)S_1+625S_0\big].              \tag{14}
 \end{split}
\]

Here S_0=q0/(1-q0), S_1=q0/(1-q0)^2,
S_2=q0(1+q0)/(1-q0)^3, S_3=q0(1+4q0+q0^2)/(1-q0)^4 and
S_4=q0(1+11q0+11q0^2+q0^3)/(1-q0)^5. They follow by successive
q0 d/dq0 operations. At theta=1/16000 the certificate gives
B(theta)<3.005479e9. This large conservative constant establishes
finiteness; it is not an estimate of the actual susceptibility.

Now normalize the periodic vacuum energy by the number of plaquettes:
e_L(theta)=E_L(theta)/(3L^3). It is concave, and translation/cubic
symmetry gives -e_L'=mu_L(W_p). Summing (11) over both source labels
and using (14) proves

\[
 0\le-e_L''(\theta)\le\mathcal B(\theta).                 \tag{15}
\]

For a compact interval |theta|<=theta_*<1/8000 one may use throughout
the weaker uniform gap gamma(theta_*) and speed 64 theta_*. Formula
(14) then gives a single finite bound B_* for all volumes and all
theta in that interval.

For completeness the energy limit and its boundary error can be proved
without assuming an infinite-volume vacuum. Tile an open cube of side
m ell by m^3 open cubes of side ell. Keep the unused linking variables
with their free Casimirs; they have ground energy zero. Only crossing
plaquette potentials are removed. Their number is

\[
 3m\ell(m\ell-1)^2-m^3\,3\ell(\ell-1)^2
                       \le6(m\ell)^3/\ell.               \tag{16}
\]

Their total norm is at most |theta| times this count. The variational
principle bounds the difference of ground energies by that norm.
This comparison uses the full tensor-product Hilbert space. The ground
energies also equal those in the Gauss-invariant sector, since each
box has YM62's unique positive, gauge-invariant vacuum; the detached
free links have the invariant constant vacuum.
Arbitrary larger cubes have leftover strips of vanishing relative
volume for fixed ell; the same bounded-potential estimate handles
those strips. Letting first the large side and then ell grow shows
that E_open,L/(3L^3) has a limit e(theta), with error <=2|theta|/L.
Adding the wrapping links as free variables and then their at most
6L^2 wrapping plaquettes proves

\[
 \boxed{|e_L(\theta)-e(\theta)|\le4|\theta|/L.}            \tag{17}
\]

Thus there is no boundary-dependent choice of this bulk energy. On a
compact coupling interval, (15) makes the derivatives e_L' uniformly
bounded and Lipschitz, while (17) gives uniform convergence of e_L.
Any uniformly convergent subsequence of derivatives must integrate
to e. This proves uniqueness of the derivative limit and

\[
 \boxed{e\in C^{1,1}_{\rm loc}((-1/8000,1/8000)),\qquad
 m(\theta)=\lim_L\mu_L(W_p)=-e'(\theta).}                 \tag{18}
\]

The mean m is nondecreasing and locally Lipschitz. The open-lattice
center involution from YM63 makes the open energy even in theta;
(17) transfers that parity to e, so m is odd, including limits of
odd-sided periodic cubes. This rules out a discontinuous bulk
plaquette-mean jump in the certified interval. It does not prove
pointwise convergence of e_L'' or analyticity of e.

There is a quantitative approximation consequence. If
|theta|+delta<=theta_* and 0<eta<=delta, set epsilon_L=4 theta_*/L.
Comparing centered difference quotients and using the common Lipschitz
bound B_* gives

\[
 |\mu_L(W_p)-m(\theta)|
       \le\min\{2,\epsilon_L/\eta+B_*\eta\}.             \tag{19}
\]

Once eta=sqrt(epsilon_L/B_*) fits inside delta this is
2 sqrt(epsilon_L B_*). It controls one bulk vacuum readout, with a
very conservative explicit constant. It does not compute the whole
coarse vacuum density or establish uniqueness of every limiting state.

## YM65-T6: the spatial estimate gives an additional continuum constraint

Under H->kappa H, both gamma and v are multiplied by kappa, so rho
is unchanged. Heat cutoff T becomes T/kappa and R_n becomes
kappa^(-n) R_n. YM64's correction moments instead scale as
kappa^(1-n); these two towers are not interchangeable. The propagation
rate in (3) is a coupling- and clock-dependent bound, not a derivation
of the physical speed of light.

At spatial spacing a, let two supports have minimum physical midpoint
separation at least D>0. T1 gives d>=D/a. If theta(a) stays inside
one compact subinterval of the certified window, use a common rho_*>0.
For readout norms and support sizes growing at most polynomially in
1/a, (5) at u=0 gives

\[
 |\operatorname{Cov}(A_a,B_a)|
       \le\operatorname{poly}(1/a)e^{-\rho_*D/a}
                         \longrightarrow0.               \tag{20}
\]

This applies to bounded regulated readouts with any specified polynomial
normalization, including bounded Wilson loops whose boundary lengths
grow polynomially. The conclusion is independent of a time-clock
rescaling. It complements YM63's temporal scaling constraint: staying
uniformly inside this strong-coupling certificate while refining spacing
forces separated connected readouts in this class to decorrelate.

A nontrivial continuum route must therefore address its coupling
trajectory and observables; these estimates do not supply such a route.
They make no exclusion claim outside the controlled interval, or for
readouts outside the stated normalization class. The full conditional
coarse-vacuum and Schur-kernel locality estimates also remain open.

## Verification

The Python 3.12 certificate verifies native lattice incidence, two
independent chain/walk counts, midpoint distances, rational Gaussian
filter enclosures and signs, inverse-response/source-Hessian identities,
tiling counts, decay budgets and a rare-conditioning counterexample.
Finite spectral fixtures are labelled and do not stand in for the actual
YM vacuum. The proof controls the full on-link Hilbert spaces; the
finite arithmetic checks are not mechanical formalization or peer review.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym65_spatial_response_locality.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym65_spatial_response_locality.py -v
```
