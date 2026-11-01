# YM-74: all-coupling electric enclosures and actual weak-coupling control

Monty Dabas. 5 October 2026. **Python 3.12 only.**

This stage proves an independent hidden-energy bound and retention of
the **entire bounded-energy physical subspace** on any prescribed
finite-lattice coupling trajectory. It also gives an explicit,
volume-uniform weak-bare-coupling plaquette concentration estimate
for the actual interacting vacuum. These results use the unchanged
YM62 SU(2) Hamiltonian and do not assume a uniform interacting gap.

For positive theta, the actual vacuum satisfies

\[
 \boxed{1-\mu_\theta(W_p)
 \le \min\!\left\{1,\frac{3}{\sqrt{2\theta}}
                         +\frac{3}{16\theta}\right\}.}       \tag{1}
\]

The estimate holds for every elementary plaquette of every periodic
cubic lattice of side at least four. Here theta->infinity is the
weak-bare-coupling direction in the declared calibration.

A complete physical electric cutoff gives, at every finite coupling,
computable enclosures for the first two exact energies and their
gap, without first computing or guessing the vacuum. On the actual
periodic side-four lattice the complete cutoff at electric energy
three is 241 dimensional. Its exact algebra yields

\[
 \boxed{0.4904<\Delta_{4^3}(1/48)<3.07874}                 \tag{2}
\]

in bare energy units. The full link Hilbert space is untruncated
in this enclosure; the discarded part is bounded rigorously.
Equation (2) is a finite-volume interval, not a physical mass value
or a volume-uniform weak-coupling gap.

The hidden separation and slow-state coverage required by YM73 can
now be achieved by an explicit electric cutoff. The remaining
continuum problem is the uniform retained spectral cancellation
and the construction of nontrivial local four-dimensional fields.
This stage does not claim a completed continuum mass-gap proof.

## Actual model and the independent audit

Let Lambda be a periodic cubic spatial lattice of side L>=4, with
p=3L^3 plaquettes and the same number of positive oriented links.
Keep

\[
 A_\theta=L_E-\theta V,\qquad V=\sum_p W_p,\qquad
 E_0=\min\operatorname{spec}A_\theta,\qquad H=A_\theta-E_0 .
                                                               \tag{3}
\]

The gauge-invariant physical Hilbert space is the full invariant
subspace of the native compact link completion. Each edge spin j
has electric energy j(j+1); the fundamental value is 3/4.
The potential obeys ||V||<=p. Write b=|theta|p.

At every fixed finite volume and finite real theta, compact
ellipticity and positivity improvement give a unique strictly
positive full-space ground state. Gauge transformations commute
with A and preserve positivity, so that same ground state is
physical. Compact resolvent and uniqueness imply a positive
fixed-regulator physical gap. No uniform value is inferred from
this qualitative fact. The constant Haar trial gives E_0<=0,
and the potential norm gives E_0>=-b.

An independent audit of YM73 found its conditional gap bracket and
weighted-resolvent estimate correct. It also identified a scope
distinction used here: weighted-resolvent control is not a uniform
heat-semigroup convergence theorem. The latter needs a separate
argument because the weighted compression has value Z, rather
than I, at time zero. No heat-convergence claim is imported.

The Hamiltonian and Schur methods are established:
[Kogut and Susskind](https://doi.org/10.1103/PhysRevD.11.395) and
[Griesemer and Hasler](https://web.ma.utexas.edu/mp_arc/c/07/07-102.pdf).
The explicit native estimates and cutoff calculations below are
proved directly. No general variational or projection principle
is claimed as a new method. The continuum target remains the
[Jaffe--Witten statement](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf),
including local field construction and a finite positive gap on R4;
the current concrete gauge group is SU(2).

## YM74-T1: a genuine weak-coupling vacuum estimate

For one native link q=(q_0,q_1,q_2,q_3), take f_n(q)=(1+q_0)^n,
with integer n>=1. Normalize the product of these functions over
all links. Its density on one link is proportional to f_n^2 dPhi.
With y=(1+q_0)/2, that density is Beta(2n+3/2,3/2). Therefore

\[
 m_n=\mathbb E(q_0)=\frac{2n}{2n+3},\qquad
 e_n=\frac{\langle f_n,L_e f_n\rangle}{\|f_n\|^2}
       =\frac{3n^2}{4(4n+1)}.                            \tag{4}
\]

The second identity uses Gamma(q_0)=(1-q_0^2)/4 and the exact
Beta ratio
E[(1-q_0)/(1+q_0)]=3/(4n+1). The other quaternion components
have zero mean. Independence of the four distinct links of a
plaquette then gives its ordered real-trace mean m_n^4.
Noncommutativity has not been replaced by scalar multiplication
pointwise; factorization is used only for this product trial's
expectation.

For theta>=0, the variational principle gives

\[
 E_0\le p(e_n-\theta m_n^4),\qquad
 0\le E_0+\theta p\le p\,d_n(\theta),
\]
\[
 d_n(\theta)=\min\{\theta,\ e_n+\theta(1-m_n^4)\}.         \tag{5}
\]

The product trial is not gauge invariant. It is legitimate for
this **ground-energy** upper bound because the full-space ground
state and the physical ground state coincide. It cannot be used
as an arbitrary trial space for higher physical eigenvalues.
The Ritz spaces in T3 are genuinely gauge invariant.

In the actual vacuum, write
T_0=<h_vac,L_E h_vac> for its nonnegative electric energy. Then

\[
 E_0+\theta p=T_0+\theta\sum_p[1-\mu(W_p)].               \tag{6}
\]

Both terms on the right are nonnegative. The unique vacuum is
invariant under translations and cubic symmetries, so every
plaquette has the same expectation. Thus, for theta>0,

\[
 1-\mu(W_p)\le d_n(\theta)/\theta,\qquad
 T_0/p\le d_n(\theta).                                   \tag{7}
\]

These are bounds on the actual interacting vacuum, not on the
trial measure used to prove (5).

Since e_n<=3n/16 and 1-m_n^4<=6/n, choose
n=ceil(sqrt(32theta)). The bounds n>=sqrt(32theta) and
n<sqrt(32theta)+1, with a non-strict version at integer endpoints,
give

\[
 d_n(\theta)\le
 \min\{\theta,\;(3/\sqrt2)\sqrt\theta+3/16\}.             \tag{8}
\]

Equations (7)--(8) prove (1). For theta=0 the exact vacuum is
the free product state, and E_0=0.

For SU(2), the Frobenius identity
||U_p-I||_F^2=4(1-W_p) converts this into actual plaquette
holonomy concentration. In the declared convention
theta=c_B/(c_E g^4), c_E,c_B>0,

\[
 \boxed{\mu\|U_p-I\|_F^2
 \le \min\!\left\{4,\;
 6\sqrt2\sqrt{c_E/c_B}\,g^2+
 \frac34(c_E/c_B)g^4\right\}.}                           \tag{9}
\]

This controls the opposite coupling direction from YM62's
small-theta gap theorem. It is uniform in the finite periodic
volume. Markov's inequality also bounds the probability of
||U_p-I||_F>=u by the right side of (9) divided by u^2.
No estimate on the renormalized local continuum curvature,
obtained after division by powers of the lattice spacing, is
asserted by this concentration bound.

## YM74-T2: a cutoff that retains every bounded-energy state

Define the **known electric spectral projection**

\[
 P_R=1_{[0,R]}(L_E)\big|_{\mathcal H_{\rm phys}},\qquad Q_R=I-P_R .
                                                               \tag{10}
\]

This is a finite-rank, gauge-compatible projection with smooth
range. It reduces L_E and preserves its operator domain. It is
not a spectral projection of the unknown interacting Hamiltonian.
Let lambda_>(R) be the lowest physical electric energy strictly
above R; in particular lambda_>(R)>R.

For arbitrary real theta one may use delta=b. For positive theta,
use the sharper, completely explicit delta=p d_n(theta) from T1.
In either case E_0+b<=delta and the operator-form inequality is

\[
 L_E\le H+\delta I.                                    \tag{11}
\]

Indeed theta V<=bI and L_E=H+theta V+E_0.
The hidden compression therefore satisfies

\[
 \boxed{Q_RHQ_R\ge[\lambda_>(R)-\delta]Q_R.}             \tag{12}
\]

No interacting gap was assumed. If lambda_>(R)>delta, this is
an independent positive hidden bound at every coupling.

Let kappa>0 be the declared physical clock and M>=0 a physical
energy ceiling. Applying (11) to an arbitrary normalized vector
in ran 1_[0,M](kappa H) gives

\[
 \boxed{\|Q_R1_{[0,M]}(\kappa H)\|^2
 \le\min\left\{1,\frac{\delta+M/\kappa}{\lambda_>(R)}\right\}.}
                                                               \tag{13}
\]

This is a whole-subspace norm estimate, not a statement about a
selected list of loop correlations. It includes the vacuum, whose
electric tail is at most delta/lambda_>(R).

There is an explicit recipe that works on any sequence of finite
volumes and couplings. Choose a fixed physical reference energy
E_*>0 and t_a>=1 tending to infinity, and take

\[
 R_a\ge\max\{3,\;t_a(\delta_a+E_*/\kappa_a)\}.            \tag{14}
\]

Then the physical hidden floor is at least t_a E_*, and for
every fixed M,

\[
 \|Q_{R_a}1_{[0,M]}(\kappa_a H_a)\|^2
       \le \max(1,M/E_*)/t_a\longrightarrow0.             \tag{15}
\]

This discharges hidden separation and retention of all bounded
energies by enlarging a constructive electric space. Its dimension
can be enormous. It does not establish the existence of any
nonvacuum state with bounded physical energy: the bounded-energy
space could contain only the vacuum.

The exact vacuum need not be inserted into the computational
basis. For the raw cut (10), if (12) is positive, the Schur
matrix S_0=C-BD^(-1)B* of H has exactly one zero eigenvalue:
harmonic extension maps its kernel to the unique vacuum.
The matrix K=kappa Z^(-1/2)S_0 Z^(-1/2), Z=I+BD^(-2)B*,
has the same one-dimensional kernel. YM73's physical gap bracket
then uses the **first positive** eigenvalue of K. The upper bound
follows by the two-dimensional min--max principle; below the
lower bound, Schur inertia leaves exactly one negative direction.
Thus there is no circular requirement to know the vacuum vector
before choosing P_R.

## YM74-T3: finite physical pencils enclose the full gap

An alternative directly uses the unshifted A, avoiding even an
exact vacuum-energy input. For P=P_R write

\[
 A=\begin{pmatrix}C&B\\B^*&D_A\end{pmatrix},\quad
 C=PAP,\quad D_A\ge dI,\quad d=\lambda_>(R)-b,
 \quad \|BB^*\|\le q,\quad q=b^2\ \hbox{always suffices}. \tag{16}
\]

The off-diagonal part is -theta PVQ because the cutoff reduces
the electric operator. It is bounded, while the hidden operator
keeps its full infinite-dimensional domain.

Let c_0<=c_1 be the first two eigenvalues of C, including
multiplicities, and define

\[
 F(c,d,q)=\frac{c+d-\sqrt{(c-d)^2+4q}}2.                 \tag{17}
\]

For j=0,1,

\[
 F(c_j,d,q)\le E_j(A)\le c_j.                            \tag{18}
\]

**Proof.** The upper inequalities are Rayleigh--Ritz. For z<d,
the exact Schur complement satisfies
C-z-B(D_A-z)^(-1)B* >= C-[z+q/(d-z)]I.
For z<F(c_j,d,q), the expression in brackets is less than c_j,
so the Schur complement has at most j negative eigenvalues.
The hidden block is positive, and inertia gives E_j(A)>=z.
Take z upwards to the value in (17). Zero-coupling or equality
cases follow by the same limit or directly from the reducing
electric projection. QED.

In particular the exact bare physical gap has the enclosure

\[
 \boxed{
 \max\{0,F(c_1,d,q)-c_0\}
 \le\Delta_\Lambda(\theta)
 \le c_1-F(c_0,d,q).}                                  \tag{19}
\]

It is valid even when the lower endpoint is zero. Such an
inconclusive interval is not converted into a positive-gap claim.

For R>=3, the constant and one physical fundamental plaquette
are retained. Their trial space and the potential norm give
c_0<=0 and c_1<=3+b. If R>3+2b, then d>c_1 and

\[
 0\le c_j-F(c_j,d,q)\le\frac{b^2}{R-3-2b}.
\]

Consequently the **total width** of the physical-gap interval is
at most

\[
 \boxed{\frac{2\kappa b^2}{R-3-2b}.}                    \tag{20}
\]

For a desired width epsilon>0, the explicit choice

\[
 R\ge4+\left\lceil2b+\frac{2\kappa b^2}{\epsilon}\right\rceil
                                                               \tag{21}
\]

suffices. The factor two in (20)--(21) accounts for both ground
and excited energy errors. At theta=0 the physical gap is exactly
three before clock scaling.

At fixed volume and coupling, (20) tends to zero as R increases.
Since the exact finite-regulator physical gap is positive, the
enclosure eventually certifies a positive lower endpoint. This is
an all-coupling terminating procedure at each fixed regulator.
The dependence b=|theta|p is retained explicitly; it is not a
volume-independent gap theorem.

For certified finite-matrix eigenintervals c_j in [c_j^-,c_j^+],
replace the lower endpoint in (19) by
max(0,F(c_1^-,d_L,q_U)-c_0^+) and the upper by
c_1^+-F(c_0^-,d_L,q_U), where d_L is a lower hidden bound and
q_U an upper mixed-norm bound. F is increasing in c,d and
decreasing in q. Finite-eigenvalue interval widths add to the
truncation width and must be included in the final error budget.

## YM74-T4: the required matrices have finite native content

The space P_R is constructible without an interacting-vacuum
ansatz. A link representation n=2j has electric energy
n(n+2)/4 and matrix-coefficient dimension (n+1)^2.
Only tuples with

\[
 \sum_e n_e(n_e+2)\le4R                              \tag{22}
\]

can occur. Tensor products of the corresponding native harmonic
polynomials give a finite electric space. Haar averaging the
vertex gauge action projects it exactly to the physical space.
The native polynomial coefficients, electric matrix and Haar
Gram and plaquette matrix elements can be chosen rational;
null Gram directions are removed before forming the pencil.
This is an exact finite construction for rational theta, rather
than a numerical assumption of completeness.

For example, a crude finite rank bound is
[sum_(n=0)^nmax (n+1)^2]^(number of links), with
nmax=floor(sqrt(1+4R)-1). A sharper unprojected count is the
sum of coefficients through degree floor(4R) in
[sum_n (n+1)^2 x^(n(n+2))]^(number of links).
Neither count is claimed to describe feasible computational cost.

Multiplication by one plaquette changes its four spins by
j_e->j_e+/-1/2. Its largest energy increase is at most
sum_(four edges)j_e+3<=2sqrt(R)+3. Hence

\[
 VP_R=P_{R+2\sqrt R+3}VP_R .                            \tag{23}
\]

The rational larger shell R+2ceil(sqrt R)+3 suffices for exact
arithmetic. Thus the mixed-norm matrix

\[
 BB^*=\theta^2\{PV^2P-(PVP)^2\}                         \tag{24}
\]

is computable from a finite electric shell, with no hidden-spin
tail omitted from that identity. Its certified largest eigenvalue
can replace b^2 in (17)--(19). A finite shell for B does not mean
that the entire hidden inverse D^(-1) has finite support.

If theta is specified through an enclosing interval rather than
an exact rational, the potential norm gives
|E_j(theta)-E_j(theta')|<=p|theta-theta'|. Therefore the physical
gap changes by at most 2kappa p|theta-theta'| at a fixed clock.
That parameter error must also be included; a declared
asymptotic-freedom-shaped function is not treated as a derived
renormalization trajectory.

## YM74-T5: a complete 241-dimensional physical cutoff

Take the actual periodic cube of side four. It has 64 vertices,
192 links and p=192 elementary plaquettes. Its simple four-cycles
are precisely the 192 plaquettes and 48 winding cycles.

The full physical electric space at energy at most three consists
of the vacuum constant and the fundamental traces of these 240
cycles. To see completeness, decompose into edge representations.
The support of any nonconstant physical sector has no degree-one
vertex, so it contains a cycle. At energy <=3 only four
fundamental edges can occur. Their support must be a four-cycle;
degree-two intertwiner constraints leave its unique trace.
The traces phi_c=2W_c form an orthonormal set under product Haar.

The graph is simple and bipartite. A five-edge nonempty support
with minimum degree two is impossible: five vertices would force
a five-cycle, while a bipartite graph on at most four vertices
has at most four edges. On a four-cycle all spins must coincide;
after j=1/2 the next such energy is eight. Every other support
has at least six edges. Therefore the next electric energy is
at least 9/2, attained by a fundamental six-edge rectangular loop.

The retained C matrix has excited diagonal 3 and only the
vacuum-to-plaquette entries -theta/2. Winding entries to the
vacuum vanish. Its spectrum is

\[
 c_0=\frac{3-\sqrt{9+192\theta^2}}2,\quad
 3\ \hbox{with multiplicity }239,\quad
 c_+=\frac{3+\sqrt{9+192\theta^2}}2 .                    \tag{25}
\]

This is the actual projected operator, not a fitted star model.
The native selection rules and the full graph enumeration verify
it independently.

The mixed norm can also be evaluated exactly. Distinct cycles
have an edge absent from the other. Haar integration gives
Phi(W_c^2)=1/4, Phi(W_p^4)=1/8, and
Phi(W_c^2 W_d^2)=1/16 for distinct cycles. A nonzero product
requires even fundamental incidence on every edge. There is no
closed set of four distinct elementary plaquettes, and no two
plaquettes have the boundary of a pair involving a winding
four-cycle. In this finite geometry these statements are checked
by all 18,336 distinct unordered plaquette-pair parity masks,
as well as every four-cycle pair.

Using (24), in vacuum/plaquette/winding order the result is

\[
 BB^*=\theta^2\operatorname{diag}
   \left(0,\frac{191I+J}{4},48I\right),\qquad
 \boxed{\|BB^*\|=\frac{383}{4}\theta^2.}                 \tag{26}
\]

For clarity, the plaquette diagonal of PV^2P is 193/4,
and its off-diagonal entries are 1/2. Subtracting (PVP)^2
removes 1/4 from every plaquette-block entry. This leaves
diagonal 48 and off-diagonal 1/4, as in (26).

Substitute c_0 from (25), c_1=3,
d=9/2-192|theta| and q=383theta^2/4 into (19).
Directed arithmetic gives the following outward intervals:

| Coupling theta | Exact bare physical gap is enclosed by |
|---|---|
| 1/8000 | [2.9999992, 3.0000006] |
| 1/128 | [2.92452, 3.00293] |
| 1/64 | [1.48847, 3.01929] |
| 1/48 | [0.4904, 3.07874] |

At theta=1/40 this particular lower estimate is inconclusive.
It does not show that the actual finite-volume gap closes there.
Increasing the complete cutoff invokes the all-coupling procedure
in T3. The displayed finite cube is not extrapolated in volume
or to a continuum mass.

## YM74-T6: what now remains in the continuum limit

The hidden and observation conditions can be made simultaneous
with arbitrary spectral accuracy: choose R_a at least the
maximum of the thresholds (14) and (21), with t_a->infinity and
epsilon_a->0. This works on any prescribed sequence of finite
volumes, finite couplings and positive physical clocks.

Let g_R(a)=c_1(a,R_a)-c_0(a,R_a) be the bare Ritz gap. The
enclosures imply

\[
 |\Delta_{\rm phys}(a)-\kappa_a g_R(a)|
       \le\epsilon_a                                  \tag{27}
\]

when all arithmetic/parameter errors are included in that
budget. In particular, convergence of the physical lattice gaps
to a finite positive m is equivalent to convergence of these
controlled finite physical pencil gaps to m.

Equation (27) is a reduction with a quantified error, not a proof
that either limit exists or is positive. The retained dimension
and the volume both grow. The difficult cancellation between
the two lowest retained energies has not been evaluated uniformly
on a physical weak-coupling trajectory.

The actual weak-coupling estimate (1) proves that elementary
holonomies concentrate near identity. It does not control the
rescaled local curvature distributions after division by powers
of a, establish four-dimensional Euclidean covariance, or
produce a surviving nonvacuum local state of bounded physical
energy. Those remain requirements for a Yang--Mills continuum
construction and the passage theorem of YM73.

| Continuum requirement | Status after YM74 |
|---|---|
| Independent positive hidden floor at arbitrary finite coupling | Explicit electric-cutoff theorem |
| Retention of every bounded-energy physical state | Explicit operator-norm estimate and cutoff recipe |
| Truncation error along arbitrary finite-regulator sequences | Explicit accuracy recipe, with growing dimension |
| Actual weak-coupling plaquette concentration | Volume-uniform bound for the interacting vacuum |
| Positive finite limiting physical gap | Retained spectral limit remains unproved |
| Nontrivial local four-dimensional Yang--Mills fields | Construction remains unproved |

Earlier chain-model weak-coupling notes and fixed-graph positivity
results do not discharge these last two rows for the actual
three-spatial-dimensional gauge lattice. The present work keeps
that model distinction explicit. Neither the finite cube interval
nor the certificate's algebraic fixtures are a claimed solution
of the four-dimensional continuum problem.

## Reproduction and evidence

The accompanying Python 3.12 certificate checks the native trial
integrals, weak-coupling constants, cutoff and coverage recipes,
Schur eigenvalue enclosures, finite-shell arithmetic, and the
complete periodic-cube four-cycle/parity calculation. Exact
matrix examples independently test the spectral inequalities.
Upstream source pins are preserved.

These calculations accompany written proofs. They do not compute
the full interacting vacuum, verify a proof in a proof assistant,
or establish a continuum mass from a passing test suite.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym74_all_coupling_electric_enclosures.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym74_all_coupling_electric_enclosures.py -v
~~~
