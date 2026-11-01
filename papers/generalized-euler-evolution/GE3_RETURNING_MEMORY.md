# GE3: returning memory survives native clock refinement

3 October 2026. Python 3.12 only.

GE2 constructs native heat from rational clock-free records. Here an
additional observer sees only part of that native coefficient carrier.
Retaining its hidden component and resetting that component after every
step have different limits. In an exact quadratic native sector,

\[
c(\tau)=\tfrac12(e^{-\tau/2}+e^{-\tau}),\qquad
c_{\rm reset}(\tau)=e^{-3\tau/4}.                            \tag{1}
\]

The first law retains memory and has asymptotic rate 1/2; the second
describes a different intervention and has rate 3/4. A first derivative
alone therefore does not determine an observer's finite-time evolution.
The full native heat law is unchanged.

## Source order and scope

Use frozen GE2 at Publications commit
`179bcca307ddeefabf794218b5144a8f76b86df9`, and its pinned GE1/YM50–YM52
native coefficient calculus, positive pairings and scalar completion.
The [YM49 bounded memory equation](../yang-mills-certified-benchmark/YM49_TIME_ZERO_OBSERVABLES.md),
T24 typed cut chain law and RKF N09 already supply the block elimination
principle. We do not claim to discover that principle here. GE3 connects
it to GE2's intrinsic unequal-step refinement and supplies an explicit
nonclosing observer on actual native quadratic polynomials.

All block-generator calculations below are on a **fixed finite native
invariant coefficient space**, with its derived positive pairing.
They do not assume that blocks of an arbitrary unbounded generator are
defined. The concrete example is also invariant in the native Phi
function pairing. The observer is declared here; it is not identified
with the instantaneous row of the interacting infinite chain.

## GE3-T1: exact clock-free memory for a nonuniform word

Let W_j be admitted contractions on a finite native core, with a fixed
orthogonal observer cut P and Q=I-P. In P/Q coordinates write

\[
W_j=\begin{pmatrix}a_j&b_j\\c_j&d_j\end{pmatrix},\quad
z_{j+1}=W_jz_j,\quad x_j=Pz_j,\quad y_j=Qz_j.
\]

Before choosing any scalar clock, finite substitution gives

\[
x_{n+1}=a_nx_n+b_nd_{n-1}\cdots d_0y_0
 +\sum_{j=0}^{n-1}b_nd_{n-1}\cdots d_{j+1}c_jx_j.             \tag{2}
\]

An empty intervening product is the identity. The kernel depends on
both endpoints and ordered intervening arrows. Proof is induction in
y_{j+1}=c_jx_j+d_jy_j, followed by the x equation. This is N09's
elimination carried out for a nonuniform word. It retains T24's return
corner b_n c_{n-1}. Setting y=0 after every step instead gives only
x_{j+1}=a_jx_j and deletes the other terms in (2).

Even when every W_j is positive and self-adjoint, individual unequal-step
return products need not be positive: b and c can have different signs
at different steps. Positivity of the stationary kernel is a separate
consequence below, not an assumption on all clock-free memory words.

## GE3-T2: retained and reset limits from the same native arrows

Take GE2's symmetric rational record means S_j, fixed normalized shape
C and budgets h_j. Set tau=sum h_j and bound all raw turn norms by
epsilon. On degree d, let L=L_C and A=PLP restricted to ran P. GE2
proves 0<=L, ||L||<=d^2/4 and
||S_j-I+h_jL||<=B_d m_{4,j}, with B_d=d^2(d^2+2)/384.
Compression preserves the first-order bound and contraction. Consequently

\[
\|P S_n\cdots S_1P-Pe^{-\tau L}P\|_d
 \le K_d\epsilon^2\tau,
\]
\[
\|(PS_nP)\cdots(PS_1P)-e^{-\tau A}\|_{\operatorname{ran}P}
 \le K_d\epsilon^2\tau,
\quad K_d=d^2(2d^2+1)/96.                                   \tag{3}
\]

**Proof.** The first inequality compresses GE2-T4. For the second,
||A||<=d^2/4, and twice integrating the finite-core exponential gives
||e^{-hA}-I+hA||<=h^2d^4/32. Telescope the contraction products and
use m_4<=2h epsilon^2 and h<=epsilon^2/2, exactly as GE2 does.
No equal-step assumption is used. Thus both limits are established
from the same rational arrows, but they implement different observer
operations. These bounds assert a fixed-core result, not a new uniform
unbounded moving-cut theorem.

## GE3-T3: continuum memory and the second-order closure marker

In that finite core write

\[
L=\begin{pmatrix}A&B\\B^\dagger&D\end{pmatrix}\succeq0.
\]

The retained limit z(tau)=e^{-tau L}z_0 obeys the native scalar-chart
equation z'=-Lz. Solving only its hidden equation by native integration
and substituting gives

\[
x'=-Ax-B e^{-\tau D}y_0+
 \int_0^\tau K(\tau-s)x(s)\,ds,\qquad
K(r)=B e^{-rD}B^\dagger\succeq0.                            \tag{4}
\]

The parameter tau is GE2's accumulated turn budget. Equation (4) is the
continuum equation of the retained limit whose discrete histories obey
(2); no separate entrywise convergence of every discrete kernel for
arbitrary meshes is asserted. Finite ODE uniqueness proves equivalence
between the full retained evolution and (4) with its specified y_0.

If D>=delta I, delta>0, and ell=||B||, then ||K(r)||<=ell^2e^{-delta r}.
Omitting memory ages greater than R makes the right-hand side error at
most ell^2 e^{-delta R} sup_{s<=tau}||x(s)||/delta. The initial hidden
term is bounded by ell e^{-delta tau}||y_0||. These are kernel/source
bounds, not an unproved bound on a separately solved truncated process.
For z>0, finite block substitution also gives

\[
P(z+L)^{-1}P|_P=
 [z+A-B(z+D)^{-1}B^\dagger]^{-1}.                            \tag{5}
\]

Invertibility follows because z+L and z+D are strictly positive.
This is the finite continuous-clock counterpart of YM49's already
proved bounded discrete Schur identity.

Let F(tau)=Pe^{-tau L}P|_P. Self-adjointness gives the exact defect
F(2tau)-F(tau)^2=(Qe^{-tau L}P)^dagger(Qe^{-tau L}P). Expansion on
the finite core proves

\[
F'(0)=-A,\qquad
\lim_{\tau\to0}\frac{F(2\tau)-F(\tau)^2}{\tau^2}
 =BB^\dagger.                                               \tag{6}
\]

In particular, dividing the defect only by tau always gives zero,
including nonclosing cuts. F is a semigroup for all times iff B=0:
necessity follows from (6), and sufficiency from the block diagonal
exponential. Also ||[P,L]||^2=||B||^2, by squaring its off-diagonal
skew block matrix. This is a cut-coupling diagnostic; it is not a new
definition of gauge or connection curvature.

## GE3-T4: exact returning memory in native quadratic observables

On the unit native quaternion variable x=(x_0,x_1,x_2,x_3), put
rho=sum x_i^2 and choose C=diag(0,1/2,1/2). Define

\[
e=x_1^2+x_2^2-\rho/2,\qquad
f=x_0^2+x_1^2-\rho/2.
\]

Direct action of the unchanged native D_a gives
L_Ce=e/2 and L_Cf=f. Their native reference means vanish, and
Phi(e^2)=Phi(f^2)=1/12, Phi(ef)=0. Thus their span is an actual
invariant centered polynomial sector. Set u=e+f and v=e-f; they are
orthogonal with equal norm. The observer retains u and hides v. In
these coordinates the inherited operator is

\[
L=\begin{pmatrix}3/4&-1/4\\-1/4&3/4\end{pmatrix}.             \tag{7}
\]

For initial (x_0,y_0)=(1,0), its visible coefficient is (1) and its
hidden coefficient is (e^{-tau/2}-e^{-tau})/2. The exact memory law is

\[
x'=-\tfrac34x+\tfrac1{16}\int_0^\tau
 e^{-3(\tau-s)/4}x(s)\,ds,\qquad
K(r)=\tfrac1{16}e^{-3r/4}.                                  \tag{8}
\]

Both exponentials and the integral follow from native scalar completion
already earned in GE1; no classical heat operator is imported. The
visible closure defect is

\[
F(2\tau)-F(\tau)^2=
 \tfrac14(e^{-\tau/2}-e^{-\tau})^2>0\quad(\tau>0).             \tag{9}
\]

Its quadratic coefficient is 1/16, matching (6). The compression's
initial rate is 3/4 but its asymptotic rate is 1/2. The full sector's
positive gap is inherited despite the absence of a visible semigroup.
This observer has no visible state clock that can justify replacing its
memory by the autonomous scalar equation x'=-(3/4)x at fixed ledger tau.
One can reparameterize one scalar trajectory, but that does not give
closure for arbitrary hidden initial states.

There is an exact finite-arrow witness as well. For equally weighted
pairs of turns of size r along e_2 and e_3, GE2's rational quaternion
has scalar component c=(16-r^2)/(16+r^2). On e and f the record mean
has respective factors alpha=c^2 and beta=2c^2-1. On u,v it is
[[a,b],[b,a]], a=(alpha+beta)/2, b=(alpha-beta)/2. Its budget is
h=r^2/2 and b/h ->1/4. Retained and reset n-step visible factors are

\[
\tfrac12\left(\prod_j\alpha_j+\prod_j\beta_j\right),\qquad
\prod_j\tfrac12(\alpha_j+\beta_j).                           \tag{10}
\]

Applying (3) proves their different limits (1), even on unequal meshes.
At a single small step b^2=O(h^2), but the many possible return paths
produce the finite memory contribution. Small local defects do not
justify discarding all return paths.

## Evidence, lineage and remaining gate

GE3 adds four scoped written results and 13 new tests, with exact native-polynomial controls,
nonuniform retained/reset enclosures, word-memory recursion, hidden-source
and Schur checks, and negative controls for false closure/positivity.
`certificates/ge3_returning_memory.py --check` is read-only. Its pins
preserve GE1/GE2 and the earlier YM inputs. Finite tests support the
general written arguments; they are not formal proof-assistant or
independent expert certification.

Projection-memory equations have established lineage, including
Hazime Mori, [Transport, Collective Motion, and Brownian Motion](https://doi.org/10.1143/PTP.33.423),
Progress of Theoretical Physics 33 (1965), 423–455. The publisher abstract
was checked for attribution; no theorem from that paper is a premise.
Our new contribution is the GE2 rational-arrow/clock limit and the exact
native quadratic observer witness. General block elimination, Schur
complements and the YM49 closure criterion are credited predecessors.

GE3 does not derive independent stochastic noise, a heat-bath law,
physical seconds, a universal observer or a new mass gap. The actual
interacting row cut from YM49 has not been identified with this declared
quadratic cut. Its closure defect, unbounded-domain control for general
cuts, moving observers, four-dimensional Yang–Mills and Clay remain open.
