# Coupling constraints and conditional coframe closure

**Monty Dabas — research continuation v0.5, 25 September 2026.**

This step asks what selects the transport family left open by U17–U20, and
whether a fixed family constrains its geometric representation. It gives a
complete signature criterion for a declared class of component couplings,
an exact obstruction to one proposed minimality rule, and a coframe equation
within an explicitly supplied homogeneous adapter. No physical gravitational
field equation is an input. The selection of the primitive physical law remains
open; the positive result is conditional mathematical closure of this adapter.

## 1. Existing selection results and the representation boundary

The [four source pins](SELECTION_SOURCE_PINS.json) supplement the existing
EMK audit in [GRAVITY_SOURCE_PINS.json](GRAVITY_SOURCE_PINS.json).

- **P1, RKF 31 §§5–6:** minimum faithful observer rank and the variational
  tail formula concern a **supplied memory operator M**. They select an observer
  for its target, not the target, a physical signature or a gravitational law.
- **P2, lambda-geometry 16:** a minimizer within a fixed rectangle family can
  cease to minimize when its endpoint restrictions are relaxed. Restricted
  minimization must not be promoted to a canonical native normalization.
- **P3, information-invariance D1, Theorems I/J:** invariance alone already
  fails to select dynamics in that model. The positive power-law selection
  explicitly needs an additional scale-usage principle. We inherit this
  distinction; we do not claim to discover general nonselection here or use
  that paper's dimensionful energy dynamics as a premise.
- **P4, its corrected ledger:** the full emergent-gravity derivation was
  withdrawn. Naming a curved one-function adapter here does not reinstate it.

The new coupling graph below has **response components as vertices**. It is
not a graph of spacetime locations, nor the carrier-location graph of U17.
Its matrix realization is already a mathematical presentation with distinctions;
it is not an identification of a distinction-free uncut ground with R^n.
As in U18, real K maps are candidate invertible algebraic transports. Their
presence does not establish that they are the positive-energy dagger-unitary
flows of RKF 55. The multi-component assembly remains an explicit hypothesis.

**v0.6 clarification:** [U24–U26](NATIVE_GRADING_AND_SEAM_MEMORY.md) show that
the sign partition reconstructed below is not automatically RKF 48's native
grading. That local grading cannot be shared by overlapping zero-extended K/R
pairs on this vertex carrier. A larger edge carrier retains it, at an exact
target-relative memory cost; this does not invalidate U21's candidate matrix
calculation or promote it to a native physical assembly.

## 2. U21 — connected coupling graph fixes a form exactly when cycles balance

Let n >= 2 and let a connected simple undirected graph on {0,...,n-1} have
one sector label s_e in {+1,-1} per edge. For i < j define

\[
G_e=s_e E_{ij}+E_{ji},\qquad
B_e=(I+t_eG_e)(I-t_eG_e)^{-1}.
\]

s=+1 is the K pair; s=-1 is the R pair. Require t_e nonzero and both
I +/- t_e G_e invertible: for real parameters K excludes t=+/-1, while R
has no further restriction. Forms are real symmetric matrices on the declared
carrier; a positive metric or a preferred signature is not supplied.

**Proposition U21.** A nonzero symmetric H satisfying B_e^T H B_e=H for
every edge exists if and only if every graph cycle contains an even number
of K edges. When it exists, the entire invariant space is the line

\[
H=\lambda\operatorname{diag}(h_0,\ldots,h_{n-1}),\qquad
h_0=1,\quad h_j=-s_e h_i\ \text{on each edge}.
\]

Every h_i is +/-1; every nonzero member is nondegenerate. For lambda>0 its
signature is the number of positive and negative h_i, respectively. For
lambda<0 the counts swap. A cycle conflict forces H=0.

**Proof.** Multiplication by the invertible Cayley denominators gives

\[
(I-tG)^T(B^THB-H)(I-tG)=2t(G^TH+HG).
\]

Thus finite invariance is equivalent to G^T H+HG=0. For an edge (i,j),
the (i,i) equation kills H_ij. The (k,i) and (k,j) equations, with k distinct
from i,j, kill H_kj and H_ki. Every vertex is incident to an edge, so every
off-diagonal entry of H vanishes. The (i,j) equation is

\[
s_e H_{ii}+H_{jj}=0.
\]

Starting with H_00, propagate this equation along a spanning tree. A further
edge is consistent precisely when the product of the edge ratios -s_e around
its cycle is +1. This holds for all cycles precisely when each cycle has an
even number of K edges. If a cycle is inconsistent, it forces its diagonal
values to zero, and connectedness forces every diagonal to zero. Otherwise
the propagation is path-independent and determines the displayed line. ∎

This is elementary signed-graph consistency and linear algebra, not a claim
of a globally new graph theorem. It extends the star calculation of U18 to
this full connected coupling class. A balanced graph also supplies a component
partition: K edges cross between opposite signs and R edges join equal signs.
That algebraic partition has not been identified with a physical causal split.

**W18.** An all-K triangle has an odd K cycle and admits only H=0. A triangle
with K edges (0,1),(0,2) and R edge (1,2) instead admits diag(1,-1,-1).
Disconnected graphs are outside U21: they can leave several independent form
scales or unconstrained isolated components. The executable API rejects them
rather than returning the connected-graph classification.

## 3. U22 — exact compatibility plus minimum couplings does not select 1+3

Fix n=4, connected simple K/R coupling graphs, and the normalization H_00=1.
Consider this specific rule: require exact invariant-form compatibility, then
minimize the number of nonzero pair couplings. This is a defined combinatorial
objective, not an already established action principle for nature.

**Proposition U22.** This rule has minimizers with mutually different signatures.
Even restricting every edge to K does not uniquely select Lorentz type.

**Proof and W17.** Every connected four-vertex graph has at least three edges.
Every tree has three edges and no cycle constraint, so U21 supplies its
normalized nondegenerate invariant form. In particular:

| Couplings | Reconstructed diagonal | Signature | Edge count |
| --- | --- | --- | --- |
| K on (0,1),(0,2),(0,3) | (1,-1,-1,-1) | (1,3) | 3 |
| K on (0,1),(1,2),(2,3) | (1,-1,1,-1) | (2,2) | 3 |
| R on (0,1),(0,2),(0,3) | (1,1,1,1) | (4,0) | 3 |

All three attain exact compatibility and the global edge minimum within the
declared class. Their signatures cannot be equated by a real invertible frame
change, even allowing reversal of the overall sign. ∎

This rejects **this selection rule**. It does not prove that no further native
criterion can select a sector. Additional admissibility or observational data
must be stated and justified independently of the desired signature. Choosing
the star because it gives 1+3 would put the answer into the assumptions.
P3's prior nonselection result is the conceptual precedent; the explicit
signature/minimum-edge obstruction is the present programme's application.

## 4. U23 — a fixed homogeneous connection constrains the coframe

Now make a **separate post-cut geometric hypothesis**. Supply a smooth chart
(u,x^1,...,x^{n-1}), a connected interval in u, a fixed frame gauge, a nonzero
constant kappa, and either the K star (s=+1) or R star (s=-1). Put

\[
H_s=\operatorname{diag}(1,-s,\ldots,-s),\qquad
A_u=0,\quad A_i=\kappa(sE_{0i}+E_{i0}).
\]

Restrict the admissible coframes to the homogeneous class

\[
e^0=a(u)\,du,\qquad e^i=b(u)\,dx^i,
\qquad a(u)b(u)\ne0.
\]

Both a and b are smooth. This **ansatz, chart, rank and connection** are inputs;
homogeneity has not been derived from the primitive seam law. kappa is supplied,
and no physical unit calibration is claimed. The common scale of H is fixed
here by H_00=1, as in U22.

**Proposition U23.** The connection is H_s-compatible. Within this coframe
class it is torsion-free exactly when

\[
\boxed{b'(u)=\kappa a(u).}
\]

Consequently, given a and one initial value b(u_0),

\[
b(u)=b(u_0)+\kappa\int_{u_0}^{u}a(v)\,dv
\]

on any patch where b remains nonzero. The coordinate change d tau=a du gives
b(tau)=kappa tau+b_0 after choosing a coordinate origin. In this coordinate,

\[
g=d\tau^2-s\,b(\tau)^2\sum_i(dx^i)^2.
\]

For curvature convention R_{mu nu}=partial_mu Gamma_nu - partial_nu Gamma_mu
+ [Gamma_mu,Gamma_nu],

\[
R_{0i}=0,\qquad
R_{ij}=s\kappa^2(E_{ij}-E_{ji})=e^{-1}F_{ij}e\quad(i,j>0).
\]

The geometry is flat when n=2 and has nonzero curvature when n>=3. Both sectors
pass the same coframe closure equation; this condition does not choose K over R.

**Proof.** Each generator satisfies G_i^T H_s+H_s G_i=0. Also

\[
(de+A\wedge e)^0=0,\qquad
(de+A\wedge e)^i=(b'-\kappa a)\,du\wedge dx^i.
\]

Vanishing torsion is therefore equivalent to the boxed equation, whose integral
gives the solution. Since a is nonzero, tau is a valid local coordinate; it is
a representation coordinate, not a derived physical clock. In the tau gauge,
the constant connection has F_0i=0 and

\[
F_{ij}=\kappa^2[G_i,G_j]
=s\kappa^2(E_{ij}-E_{ji}).
\]

The coframe acts by the same scale b on every spatial component, so it commutes
with these spatial curvature matrices. U20 applies and gives R=e^{-1}Fe.
For an independent calculation, the metric's only nonzero Christoffel entries
are Gamma^0_ii=s kappa b and Gamma^i_0i=Gamma^i_i0=kappa/b, without summation
over i. Their derivatives and products give the displayed R_0i and R_ij. ∎

This genuinely constrains the coframe **within the declared class**: once a,
kappa and an initial value are supplied, b is no longer a freely supplied
function. It does not construct the chart, pick that class, select kappa or
establish a unique geometry over all coframes. The formula is standard
connection/coframe mathematics used as an explicit adapter, not a new
gravitational field equation. The input constants and their possible gauge
or calibration equivalences require a separate operational interpretation.

**W19, fixed-connection control.** Set a=1, kappa=1 and try b=1+u^2. Torsion
has coefficient 2u-1; it vanishes at u=1/2 but not on any open interval. Thus
one successful point check cannot certify patch-wide coframe closure. W15 in
v0.4 admitted arbitrary f by setting A=f'K dv, which changed the connection
with f. U23 instead fixes A first; these are different admissibility problems.

**W20, alternative sector.** With n=3, a=1, kappa=1 and b=3+u, both K and R
are torsion-free on b!=0. At u=0 their metrics are diag(1,-9,-9) and
diag(1,9,9), respectively. Their nonzero spatial curvature matrices have
opposite signs. Conditional coframe closure has not removed sector nonselection.

## 5. Exact evidence and the remaining selection problem

The implementation uses exact rational arithmetic and retains all earlier
evidence. `selection_certificate.py` compares graph sign propagation against
the full symmetric invariant-form nullspace of finite Cayley maps. That oracle
does not assume diagonal forms or receive the expected signature.

| Declared domain | Coverage |
| --- | --- |
| Every labelled simple graph with absent/K/R edges, n=2,3,4 | 759 graphs |
| Connected graphs compared against full invariant-form equations | 646 |
| Balanced / inconsistent connected cases | 322 / 324 |
| Disconnected cases explicitly outside U21 | 113 |
| Supplied homogeneous chart points, n=2..5, both sectors | 48 |
| Connection versus independently metric-derived curvature matrices | 240 |
| Additional named regressions / package total | 13 / 63 |

For the chart checks, b=3+kappa u, kappa in {-1,1/2}, and u in {-1,0,2}.
These are finite checks of the identities, not a numerical proof of the
smooth theorem. A separate nonzero-second-derivative control checks that the
generic Christoffel oracle can detect curvature outside the linear profile.
The source-bound certificate additionally binds this manuscript, its four-file
source manifest, the implementation, evidence generator and regression module.

The next physical selection problem is now sharper: identify native data that
justify a coupling graph and sector **before** choosing the desired signature,
and supply a target on which rival choices predict different observations.
Then derive, enlarge or reject the homogeneous adapter from those data. There
is no derivation yet of a universal gravitational coupling, a quantum readout
law, a physical spacetime dimension, or a full quantum-gravity theory. These
remain proof obligations, not consequences of passing the finite certificate.
