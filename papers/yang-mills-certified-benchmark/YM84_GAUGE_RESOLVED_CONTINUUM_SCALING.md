# YM84 — Gauge-resolved response and the continuum scale

Monty Dabas. 5 October 2026. Python 3.12 evidence packet.

## What is established

The request to close the four-dimensional continuum problem makes it
necessary to check which response norms describe physical states. YM83's
full-space criterion is sufficient, but some large correlations can live
entirely outside the Gauss sector. This stage proves that this distinction
is quantitatively important, and constructs a sector-resolved replacement.

- In the **actual homogeneous lattice vacuum**, the complete response
  between two open paths around a plaquette is a central SU(2) convolution.
  Its representation channels are exactly normalized Wilson-character
  expectations. Weak-coupling concentration forces a growing family of
  almost-aligned channels and a lower bound on every finite-rank omitted
  response tail. All these open-path markers have zero Gauss projection.
- The conditional-projection gap argument can instead be run entirely in
  the physical Hilbert space, with explicitly defined physical response
  norms and integrated block energies. Their weak-coupling bounds in the
  full overlapping lattice remain to be proved.
- An **isolated original one-plaquette Hamiltonian**, used as a controlled
  comparison, exhibits the distinction exactly. Its physical two-block
  resampling gap is 2, even while the unrestricted resampling gap tends to
  zero. Its actual quantum physical gap obeys
  \((3/4)\sqrt\theta\le\Delta_\square\le9\sqrt\theta\) for \(\theta\ge1\).
- If an analogous local physical stiffness of order \(\sqrt\theta\) is
  established in the full lattice, a surviving finite-energy continuum
  state forces the sufficient correlation margin to close at least
  linearly in a. A constant margin would instead drive the certified
  physical spectrum to infinity.

This is a repaired continuum route, not a construction of four-dimensional
Yang--Mills or a proof of its positive mass gap. The new actual-vacuum
response statement is unconditional. The full-lattice physical block
estimates, nonzero continuum states and field-theoretic identification
remain open. Representation orthogonality, projection geometry and
oscillator comparison are established inputs, not new general theories.

## 1. Model and sectors

Keep the homogeneous periodic cubic model with side L>=4 (and even L
for any blocking that requires it):

\[
 A_\theta=L_E-\theta\sum_p W_p,\quad H=A_\theta-E_0,
 \quad d\mu=h^2d\Phi,\quad \mathcal L=h^{-1}Hh.
 \tag{1}
\]

h is the unique positive normalized ground state. The native link metric
has fundamental Casimir 3/4. Local gauge transformations preserve h, Haar
measure and \(\mu\). Let \(\Pi_G\) be averaging over the finite lattice
gauge group. Its range \(\mathcal H_G\subset L^2(\mu)\) is the physical
space. Constants represent the vacuum. No gauge fixing or replacement
of the actual measure by a chosen Gibbs density is made.

The physical calibration is unchanged:

\[
 \mathsf H_a=\kappa_a H_a,\quad
 \theta_a=c_B/(c_Eg_a^4),\quad
 \kappa_a=c_Eg_a^2/a,\quad
 \kappa_a\sqrt{\theta_a}=c/a,\quad c=\sqrt{c_Ec_B}>0.
 \tag{2}
\]

YM76 gives, for an elementary plaquette in this actual homogeneous model,

\[
 \delta:=\mu(1-W_p)\le\frac9{4\sqrt\theta},\qquad\theta\ge25,
 \tag{3}
\]

uniformly in total volume. The isolated one-plaquette comparison below
has its own vacuum and its own concentration estimate; (3) is not
silently transferred to it.

## 2. Theorem 84.1 — Actual open-path response and a growing omitted tail

Split the four edges of an elementary plaquette into two edge-disjoint
paths P,Q with distinct common endpoints s,t, oriented from s to t. Let
U and V be their ordered transports. Then \(W_p=\tfrac12\operatorname{Tr}
(UV^{-1})\). Their joint marginal in the actual vacuum has the form

\[
 d\nu(U,V)=q(UV^{-1})\,dU\,dV,
 \quad\int q(D)dD=1,
 \quad q(gDg^{-1})=q(D).
 \tag{4}
\]

q is smooth and strictly positive at each finite regulator; both
marginals are Haar. The centered conditional response operator is

\[
 Tv(U)=\int[q(UV^{-1})-1]v(V)dV.
 \tag{5}
\]

For n=1,2,... label the spin-n/2 representation by its dimension
\(d_n=n+1\), character \(\chi_n\), and normalized matrix entries
\(e_{n,ab}=\sqrt{d_n}\,D^{(n/2)}_{ab}\). Define

\[
 \lambda_n=\frac{\mu\chi_n(UV^{-1})}{n+1}.
 \tag{6}
\]

Then

\[
 Te_{n,ab}=\lambda_ne_{n,ab},\qquad
 c_{P,Q}=\sup_{n\ge1}|\lambda_n|,\qquad
 \Xi_{P,Q}=\|T\|_{\rm HS}^2
       =\sum_{n\ge1}(n+1)^2\lambda_n^2.
 \tag{7}
\]

In particular \(\lambda_1=\mu W_p=1-\delta\), but \(c_{P,Q}\)
cannot in general be identified with this one coefficient. All channels
must be retained or bounded. For every n,

\[
 0\le1-\lambda_n\le\frac{n(n+2)}3\,\delta.
 \tag{8}
\]

If \(N(N+2)\delta\le3/2\), set

\[
 D_N=\sum_{n=1}^N(n+1)^2
     =\frac{N(2N^2+9N+13)}6.
 \tag{9}
\]

There are D_N orthonormal response directions with eigenvalues at least
1/2. Every operator R of rank at most q_0 therefore satisfies

\[
 \boxed{\|T-R\|_{\rm HS}^2\ge\frac{(D_N-q_0)_+}{4}.}
 \tag{10}
\]

For \(\theta\ge256\), the choice
\(N=\lfloor\theta^{1/4}/2\rfloor\) is permitted by (3). Consequently

\[
 \boxed{\Xi_{P,Q}\ge\theta^{3/4}/768.}
 \tag{11}
\]

To keep the omitted Hilbert--Schmidt norm at most \(\varepsilon\), a
rank-q_0 approximation must have
\(q_0\ge D_N-4\varepsilon^2\ge\theta^{3/4}/192-4\varepsilon^2\).
This is a necessary rank bound for this open-path protocol, not a
statement that every useful physical observable needs that many modes.

**Proof.** Endpoint gauge transformations act as
\((U,V)\mapsto(g_sUg_t^{-1},g_sVg_t^{-1})\). Invariance under their
common right action makes the joint density a function of \(UV^{-1}\);
invariance under the common left action makes that function central.
The marginal is smooth and positive because disjoint path products can
replace one link on each path by a Haar-preserving smooth coordinate
change, after which the remaining variables are integrated out. Its
normalization and Haar marginals give (4).

Schur orthogonality and centrality imply
\(\int q(D)D^{(n/2)}(D^{-1})dD=\lambda_n I\), which proves (7).
For SU(2) the character is real and invariant under inversion, so the
coefficients are real. Smoothness at fixed regulator justifies the
Hilbert--Schmidt expansion. Equivalently it is Parseval for q-1.

If \(W_p=\cos r\), the normalized character is the average of
\(\cos(kr)\) for k=-n,-n+2,...,n. The integer inequality
\(1-\cos(kr)\le k^2(1-\cos r)\), obtained by summing k unit chords,
and \((n+1)^{-1}\sum k^2=n(n+2)/3\), prove (8).

On the D_N-dimensional span in (9), T has minimum singular value at
least 1/2. The restriction of R to that span has a kernel of dimension
at least D_N-q_0. Sum the squared norm of (T-R) on an orthonormal family
in that kernel. This proves (10), including for non-spectral rank-q_0
approximations.

Write \(s=\theta^{1/4}\ge4\). Then
\(N(N+2)\le s^2/4+s\le2s^2/3\), so (3) gives the required condition.
Also \(N\ge s/4\) and \(D_N\ge N^3/3\ge s^3/192\), proving (11).
QED.

### Scope of this actual-vacuum obstruction

For a two-block partition of *all* links containing P in one block and
Q in the other, the path response is a compression of the full block
response. The same rank obstruction follows by projecting its images
back to the path span. In particular its unrestricted maximal correlation
is at least \(1-\delta\).

For more than two blocks, YM83 uses an **exterior-conditioned** response.
The unconditioned calculation above does not bound each of those
conditional coefficients. No such substitution is made.

Crucially, for n>=1 every marker \(e_{n,ab}(U_P)\) or
\(e_{n,ab}(U_Q)\) obeys

\[
 \Pi_G e_{n,ab}(U_P)=\Pi_G e_{n,ab}(U_Q)=0.
 \tag{12}
\]

Average at endpoint s to see this. Thus the entire explicit large-response
tower is gauge noninvariant. It obstructs an unrestricted response-tail
certificate, but does **not** prove a small physical quantum gap.

## 3. Theorem 84.2 — A comparison entirely in the Gauss sector

Partition the links into disjoint nonempty blocks B_i. Use the actual
measure to define

\[
 P_i=\mathbb E_\mu[\cdot\mid X_{\bar B_i}],\quad Q_i=I-P_i,
 \quad P_{ij}=\mathbb E_\mu[\cdot\mid X_{\overline{B_i\cup B_j}}].
 \tag{13}
\]

All these projections commute with \(\Pi_G\). Define the *global
physical operator norm*

\[
 c^G_{ij}=\|(P_iP_j-P_{ij})|_{\mathcal H_G}\|,
 \qquad r_G=\max_i\sum_{j\ne i}c^G_{ij}.
 \tag{14}
\]

This is not a conditional covariance computed after discarding exterior
variables. It keeps the whole gauge-invariant Hilbert space.

Suppose numbers \(\rho_i^G>0\) satisfy the integrated inequalities

\[
 \mu|\nabla_{B_i}f|^2\ge\rho_i^G\|Q_i f\|_\mu^2
 \quad\text{for every smooth gauge-invariant f}.
 \tag{15}
\]

Then, if r_G<1,

\[
 \boxed{\operatorname{gap}H|_{\rm Gauss}
     \ge\rho_*^G(1-r_G),\qquad\rho_*^G=\min_i\rho_i^G.}
 \tag{16}
\]

YM83's ordinary conditional floors supply valid choices in (15), but
stronger physical floors may be possible. A block Hamiltonian's gap
with its *own* vacuum does not automatically supply (15).

**Proof.** A fixed gauge transformation acts separately on every link.
Invariance of \(\mu\) makes each conditional projection equivariant,
including P_ij, so their restrictions to \(\mathcal H_G\) are orthogonal
projections. Their common range is the range of P_ij. This follows from
strict positivity of the density: a function measurable off B_i and
off B_j is measurable off their union. The same identity remains true
after gauge restriction, since all three projections commute with
\(\Pi_G\).

Apply the two-projection decomposition to these restricted projections,
subtracting their common range. It gives

\[
 (Q_iQ_j+Q_jQ_i)|_{\mathcal H_G}
       \ge-c^G_{ij}(Q_i+Q_j)|_{\mathcal H_G}.
 \tag{17}
\]

For \(K_G=(\sum_iQ_i)|_{\mathcal H_G}\), summing proves
\(K_G^2\ge(1-r_G)K_G\). Its kernel is the constants, as in YM83.
Thus \(\sum_i\|Q_i f\|^2\ge(1-r_G)\operatorname{Var}_\mu f\).
Sum (15) and use the exact ground-state transform. Closure proves
(16) on the physical form domain. QED.

The unrestricted operator norm bounds the restricted norm, but its
large gauge-noninvariant channels need not occur in (14). Equation
(12) provides an actual example of precisely such channels. Gauge
projection alone does not bound the remaining norms in (14).

## 4. Theorem 84.3 — Exact one-loop separation and quantum stiffness

This section is a **separate isolated four-link cycle**, with its
unchanged original potential \(-\theta W\), rather than the overlapping
cubic lattice of (1). Denote its actual vacuum by \(h_\square\).

Split the cycle into two two-link paths. Every gauge-invariant function
is a class function of the loop holonomy. Conditional integration over
either path makes such a function constant: change that path product
to the relative loop holonomy, and use its central vacuum density.
Therefore on the physical space,

\[
 P_1=P_2=P_0,\qquad K_G=2(I-P_0),\qquad c^G_{12}=0.
 \tag{18}
\]

The physical resampling gap is exactly 2 at all finite theta. The full
resampling gap is \(1-\sup_{n\ge1}|\lambda_n|\), with the characters
of this loop's own vacuum, and will tend to zero by the estimate below.

For the actual quantum physical gap, for every \(\theta\ge1\),

\[
 \boxed{\frac34\sqrt\theta\le\Delta_\square(\theta)
       \le9\sqrt\theta.}
 \tag{19}
\]

**Radial proof.** Set W=cos r, 0<r<pi. Physical Haar has radial density
\((2/\pi)\sin^2r\,dr\), and the four-link electric operator on class
functions is \(-\partial_r^2-2\cot r\,\partial_r\). Multiplication
by \(\sqrt{2/\pi}\sin r\) is a unitary identification with ordinary
L2(0,pi); the closed physical form becomes the Dirichlet form for

\[
 B_\theta=-\partial_r^2+\theta(1-\cos r)
 \quad\text{on }H^1_0(0,\pi),
 \tag{20}
\]

up to the scalar shift -1-theta. Endpoint regularity gives Dirichlet
conditions; taking the form closure gives all of \(H^1_0\). Scalar
shifts do not change gaps. Let \(\epsilon_0<\epsilon_1\) be the first
two eigenvalues of (20).

On [0,pi], \(2r^2/\pi^2\le1-\cos r\le r^2/2\). Extension by zero
to the half-line and min--max comparison with its Dirichlet harmonic
oscillator give

\[
 \epsilon_1\ge7\sqrt{2\theta}/\pi.
 \tag{21}
\]

Indeed \(-\partial_r^2+\omega^2r^2\) on the half-line with a
Dirichlet endpoint selects the odd full-line Hermite functions, whose
first two energies are 3omega and 7omega.

For the upper ground energy, take the normalized sine on (0,b), extended
by zero, which is an admissible form-domain trial function. Its second
moment is \(b^2(1/3-1/(2\pi^2))\). Hence, with
\(A=1/6-1/(4\pi^2)>0\),

\[
 \epsilon_0\le\pi^2/b^2+\theta A b^2.
 \tag{22}
\]

Choose \(b^2=\pi/\sqrt{\theta A}\). This is <=pi^2 when theta>=1,
since \(\pi^2/6-1/4>1\). Consequently

\[
 \epsilon_0\le C_0\sqrt\theta,
 \quad C_0=\sqrt{2\pi^2/3-1},\quad
 \Delta_\square\ge(7\sqrt2/\pi-C_0)\sqrt\theta.
 \tag{23}
\]

The rational bounds \(\pi<22/7\), \(\sqrt2>140/99\) and
\(C_0<19/8\) give

\[
 7\sqrt2/\pi-C_0>
 3430/1089-19/8=6749/8712>3/4.
 \tag{24}
\]

For the upper gap, use the first two sine modes on (0,b). Their maximum
kinetic Rayleigh quotient is \(4\pi^2/b^2\), and the potential is at
most \(\theta b^2/2\). Choosing
\(b^2=2\sqrt2\pi/\sqrt\theta\le\pi^2\) for theta>=1 yields
\(\epsilon_1\le2\sqrt2\pi\sqrt\theta<9\sqrt\theta\).
Since \(\epsilon_0\ge0\), this completes (19). QED.

Equation (23) also proves for this loop's actual vacuum

\[
 \delta_\square:=\mu_\square(1-W)
       \le C_0/\sqrt\theta.
 \tag{25}
\]

The potential contribution is bounded by the total nonnegative energy
\(\epsilon_0\). Thus the unrestricted resampling gap is at most
\(\delta_\square\to0\), while (18) and (19) hold. This gives an
explicit example where the small unrestricted response margin comes
from excluded gauge channels and does not diagnose the physical gap.

For a physical function of r, each of the two equal paths contributes
exactly half of its total Dirichlet energy. Since Q_i removes its mean,
(19) supplies the sector bounds
\(\rho_i^G\ge(3/8)\sqrt\theta\) in this isolated model. Combining
them with the exact physical resampling gap 2 in (18), rather than
the general lower bound 1-r_G, recovers (19). This is a verified example
of the stronger *physical* local estimate sought in (15).

There is no claimed identity between \(h_\square\) and a conditional
vacuum of the full lattice. Overlapping plaquettes and conditional
exterior data remain part of the full-lattice problem.

## 5. Theorem 84.4 — The margin required for a finite physical scale

Return to the actual homogeneous lattice and the physical clock (2).
Suppose a sequence of physical block estimates (15) has
\(\rho_{*,a}^G\ge d\sqrt{\theta_a}\), with d>0 independent of a,
and \(\eta_a=1-r_{G,a}>0\). Then (16) gives

\[
 \operatorname{gap}\mathsf H_a|_{\rm Gauss}
       \ge\frac{cd}{a}\eta_a.
 \tag{26}
\]

If there are centered physical states u_a with
\(\|u_a\|^2\ge v>0\) and physical energy at most E<infinity, then

\[
 \boxed{0<\eta_a\le\frac{E}{cdv}\,a.}
 \tag{27}
\]

**Proof.** Rayleigh's principle applied to u_a gives an upper gap bound
E/v. Combine it with (26). QED.

Thus, for this improved-stiffness strategy, a uniform positive margin
is incompatible with a surviving bounded-energy state. If eta_a/a
diverges, (26) pushes every centered heat correlation of bounded norm
to zero at every positive physical time. A nonzero norm limit could
then not belong to a strongly continuous physical-time Hilbert-space
semigroup. This is the spectral escape mechanism isolated in YM73,
now expressed directly in the gauge-resolved block-response variables.

Conversely, if the actual estimates give
\(\eta_a\ge ma/(cd)\), then (26) certifies the lattice physical
lower gap m. Passing that bound to a nontrivial continuum theory still
requires convergent correlations on a dense physical state core,
nonzero local states, field identification, covariance and locality.
None of those inputs follows from the inequality alone.

The scale-matched target is therefore to establish, in the actual
overlapping vacuum, physical local stiffness and a controlled margin
of order a, together with state survival. Equation (27) is conditional
on those hypotheses; it is not a proof that the actual margin has that
order. Nor does a failed sufficient row estimate prove a vanishing gap.

The one-plaquette control itself has physical quantum gap between
\(3c/(4a)\) and \(9c/a\). It is a microscopic stiffness model, not a
finite-mass continuum theory. Simply repeating independent such cells
does not construct continuum Yang--Mills.

## 6. Verification and remaining work

The Python packet uses exact SU(2) Haar moments and character polynomials
to check every response channel for finite-degree positive comparison
densities. It checks Parseval and omitted tails independently, the native
character derivative constants, the rational oscillator margins, and the
physical clock algebra. A finite-group cycle verifies the gauge projection
and the contrast between full and physical resampling operators directly.
These comparison densities and finite groups are explicitly not the
homogeneous Yang--Mills vacuum. The actual-vacuum statements are the
written proofs using (3), symmetry and the exact ground-state transform.

The continuum problem is still open in this programme. What has changed
is the target: the actual large open-path response tower has been located
and shown to be absent after Gauss projection, and an exact model proves
that physical stiffness can grow while the unrestricted response test
degenerates. The remaining target (15), (14), (27) is explicitly in the
physical sector and at the physical scale.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym84_gauge_resolved_scaling.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym84_gauge_resolved_scaling.py -v
```

No independent external review or proof-assistant verification is claimed.
The manifest preserves the previous source and evidence chain. Finite
checks certify their declared arithmetic and bytes, not a continuum proof.

## References and lineage

- [YM76](YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS.md), the actual
  homogeneous vacuum concentration estimate and native normalization.
- [YM83](YM83_ACTUAL_VACUUM_CONDITIONAL_RESPONSE.md), the conditional
  projection criterion and finite-response omitted-tail obligation.
- [YM73](YM73_PHYSICAL_MASS_MEMORY_METRIC.md), physical-clock escape
  and the dense-core requirements for a continuum gap passage.
- [YM77](YM77_PHYSICAL_SMEARING_AND_CONTINUUM_SURVIVAL.md) through
  [YM81](YM81_LOCAL_THIRD_MOMENT_AND_VACUUM_RESPONSE.md), the open
  nonzero continuum-state survival problem.
- Q. Qin, *Analysis of two-component Gibbs samplers using the theory of
  two projections*, [arXiv:2201.12500v4](https://arxiv.org/html/2201.12500v4),
  for the established projection-spectrum machinery.
- A. Jaffe and E. Witten, *Quantum Yang--Mills Theory*,
  [official problem statement](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf),
  for the full four-dimensional existence and mass-gap target.
