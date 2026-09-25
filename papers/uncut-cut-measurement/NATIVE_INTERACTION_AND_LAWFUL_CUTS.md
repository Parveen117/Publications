# Native-compatible interaction, lawful cuts and finite response reconstruction

**Monty Dabas — research continuation v0.7, 25 September 2026.**

U24–U26 retained local native grading on a direct sum of edge carriers, but
separate-edge steps commuted. This continuation supplies an explicit graded
interaction between those carriers, proves which real-linear cuts preserve
its entire declared action family, and constructs a minimal finite transcript
from one scalar readout. A genuine information-losing quotient preserves the
nonzero interaction. Neither spacetime nor a classical physical field equation
is assumed in this module. Its interaction graph and strengths are declared
candidate data; their physical selection remains a separate obligation.

## 1. What the repository already supplies

The [four-file source audit](INTERACTION_SOURCE_PINS.json) establishes:

- **I1, RKF 51:** anti-self-dagger generators give energy-preserving native
  Cayley steps. On the real scalar subcarrier, the dagger is transpose, so
  real antisymmetric generators satisfy the rule. Its step parameter h equals
  2q in our notation below. The parameter labels a protocol, not a physical
  clock derived by this construction.
- **I2, RKF 55:** R is a licensed flow generator; real K alone is not. Dagger
  grading and the EMK cut grading are different. Both generators below are real
  anti-self-dagger, although one is cut-even and the other cut-odd. We do not
  identify their commutator with I1's distinct iota-turn exchange law.
- **I3, RKF 48:** the local EMK cut is K, with R cut-odd. Nonzero even/odd
  commutators are already part of the native algebra.
- **I4, EMK-T1, T1–T3 and claim boundary:** J=K and cross-channel commutator
  terms already appear in the master-tensor capsule. Its component list and
  physical channel choices remain declarations. The existence of a mixed
  commutator is therefore not itself a new result here.

Target-relative observer rank, exact quotient descent and the native grading
discipline are inherited from the earlier source pins. This module specifies
and solves one further representation problem under those rules; it does not
derive a universal EMK tensor from an uncut primitive.

## 2. Declared carrier and interaction family

Retain m>=1 two-component edge carriers and write

\[
W=\mathbb R^m\otimes\mathbb R^2,\qquad J=I_m\otimes K,
\qquad K=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]

The factor R^m labels retained edge identities from the preceding model.
Supply a simple interaction graph Gamma on these labels. Its vertices are
**edge carriers**, not spacetime locations or the original endpoint vertices.
The graph can be disconnected. It is not claimed to be forced by the native
primitive or to define physical locality.

For each label a, and for each link {a,b} with a<b, define

\[
G_a=E_{aa}\otimes R,\qquad
B_{ab}=(-E_{ab}+E_{ba})\otimes I_2.
\]

G_a rotates within an edge carrier; B_ab mixes the two retained carriers
without discarding their identities. Include one independently addressed G_a
step for every a and one B_ab step for every link. All their parameters are
known and nonzero. For any generator A in this list use

\[
T_q(A)=(I+qA)(I-qA)^{-1},\qquad q\in\mathbb R\setminus\{0\}.
\]

The mathematical state, observer and subspace claims below are **real-linear**
and finite-dimensional. They do not quantify over every representation of the
full native scalar algebra, nor over noisy, nonlinear or probabilistic cuts.

## 3. U27 — genuine interaction with one common native grading

**Proposition U27.** Every displayed generator is real anti-self-dagger. Its
Cayley step is orthogonal and preserves the declared quadratic recognition
form x^T x. The common cut satisfies

\[
JB_{ab}J=B_{ab},\quad JG_aJ=-G_a,
\]

and the finite steps satisfy

\[
JT_q(B_{ab})J=T_q(B_{ab}),\quad
JT_q(G_a)J=T_q(G_a)^{-1}.
\]

For every link {a,b},

\[
\boxed{[B_{ab},G_a]=(E_{ab}+E_{ba})\otimes R\ne0.}
\]

For nonzero finite parameters the corresponding two Cayley steps do not
commute, and their group commutator is not the identity.

**Proof.** Both generators have transpose equal to their negative, so I-qA
is invertible: for a real vector v,

\[
\|(I-qA)v\|^2=\|v\|^2+q^2\|Av\|^2.
\]

The commuting polynomials I+qA and I-qA give T_q(A)^T T_q(A)=I and
T_q(A)^(-1)=T_(-q)(A). The grading identities follow from K^2=I and KRK=-R.
The bracket follows by multiplying E_aa with -E_ab+E_ba on both sides.

On the two linked carriers write U=T_h(B_ab) and V=T_t(G_a). With

\[
c_h=\frac{1-h^2}{1+h^2},\quad s_h=\frac{2h}{1+h^2},\quad
Q_t=T_t(R),
\]

their active blocks are

\[
U=\begin{pmatrix}c_hI&-s_hI\\s_hI&c_hI\end{pmatrix},
\qquad V=\begin{pmatrix}Q_t&0\\0&I\end{pmatrix}.
\]

Hence the two off-diagonal blocks of UV-VU are s_h(Q_t-I). Both s_h and
Q_t-I are nonzero, and Q_t-I is invertible. Thus UV!=VU, equivalently
UVU^(-1)V^(-1)!=I. No compression is used in this calculation. ∎

**Represented quarter-turn.** The same family commutes with
Z=I_m tensor R, where Z^2=-I and Z^T=-Z. It therefore preserves a represented
complex structure as well as the quadratic form. J reverses Z, so Z is not
central in the entire EMK algebra. This calculation does not identify the
operator Z with the native central scalar iota, nor supply a Born/readout law
or a physical quantum state interpretation.

**W27, an exact order test.** Take m=2, one mixing link, h=t=1/2, input
x=(0,0,0,1), and read out its first coordinate after the protocol. Then

\[
\ell UVx=0,\qquad \ell VUx=16/25.
\]

Rightmost operators act first: R followed by mixing gives zero; mixing followed
by R gives 16/25. This is a calibrated model prediction of actual transport
noncommutativity, not observed physical data or a claim of metric curvature.

## 4. U28 — exact cut classification and minimum linear memory

Let Gamma have connected components Gamma_1,...,Gamma_s. Write
W_c=R^(Gamma_c) tensor R^2. Include all the separately addressed controls
specified in section 2, with nonzero parameters.

**Proposition U28.** The common invariant real subspaces of this complete
family are exactly the direct sums of whole component carriers W_c. Therefore:

1. A surjective real-linear cut C admits exact descended steps for every
   declared control if and only if ker C is a direct sum of whole W_c.
   Any such quotient also inherits J, since J preserves each W_c.
2. For any nonzero initial real-linear readout C, let O be the union of
   components on which C is not identically zero. The smallest closed linear
   observer retaining C has rank 2|O|. Its minimum additional number of
   independent real scalar channels is

   \[
   \boxed{2|O|-\operatorname{rank}C.}
   \]

   Here |O| counts retained edge labels, not original graph vertices. A connected
   interaction graph makes every nonzero exact closed observer faithful on W.
3. The complete space of compatible symmetric bilinear forms is

   \[
   \boxed{H=\bigoplus_{c=1}^{s}\lambda_c I_{2|\Gamma_c|}.}
   \]

   In particular, on a connected graph every nonzero compatible form is definite.

**Proof of the subspace classification.** A finite-dimensional subspace
invariant under U=T_q(A) is invariant under A. Indeed, U+I is invertible and

\[
A=q^{-1}(U-I)(U+I)^{-1}.
\]

U+I maps the subspace injectively into itself, hence onto itself; its inverse
preserves it. The converse follows in the same way from I-qA. Thus invariance
under the finite controls is equivalent to invariance under their generators.

Let N be invariant and suppose x in N has a nonzero component at label a.
Since G_a^2=-E_aa tensor I, N contains that isolated component x_a. The real
vectors x_a and G_a x_a span its entire two-dimensional carrier: a nonzero real
vector and its R rotation are linearly independent. For every mixing link
{a,b}, B_ab maps this plane bijectively to the b plane. Connectedness propagates
it through the entire component. Every component on which any vector of N is
nonzero is therefore fully contained in N. Conversely, sums of whole W_c are
visibly invariant. This proves the classification.

**Proof of descent and minimality.** The exact quotient condition is the
inherited kernel-invariance gate. The largest invariant subspace contained in
ker C is precisely the sum of components on which C vanishes identically.
Any closed observer retaining C must quotient by a subspace of that kernel,
so its rank is at least 2|O|. Projection onto the coordinates in O attains the
bound and contains C in its row space. Extending independent rows of C to a
basis of this observer gives the stated minimum repair count.

**Proof of the form classification.** Finite Cayley invariance of a symmetric
H is equivalent to A^T H+HA=0. Applying this to every G_a forces all off-diagonal
two-component blocks of H to vanish, since R is invertible. Each diagonal
symmetric block commuting with R is lambda_a I_2. The B_ab equation then
forces lambda_a=lambda_b along each link. This gives the stated form space,
and each such form satisfies every equation. ∎

The classification needs the full independently addressed catalogue and real
scalar type. Restricting the controls, changing the target, or changing the
representation changes the conclusion. It is not a universal prohibition on
information-losing cuts or on other physical geometries.

**W28, a proper information-losing interacting quotient.** Take m=3 and only
the mixing link (0,1). The component {2} is separate. Let C keep the first four
coordinates and discard the last two. Every local R step and the mixing step
descends, as does J. The (0,1) interaction remains noncommuting on the quotient:

\[
[\bar U,\bar V]C=C[U,V]\ne0.
\]

The source target form diag(I_4,0_2) descends exactly to I_4, in agreement with
U19. The full source form I_6 does not descend through this lossy cut; it is a
different target. Thus both transport and the declared bilinear target pass
the exact cut gates while two source coordinates genuinely remain lost.

**W29, a forbidden partial cut.** With two linked carriers, keeping only one
carrier fails mixing descent. Its two visible scalars need two additional
scalars for complete closure. The lost coordinates belonged to the same
interacting component and can later return to the readout.

**W30, cancellation is not an unobserved component.** With two unlinked
carriers, C=(1,0,-1,0) vanishes on x=(1,0,1,0). A local R step at label 0
with t=1/2 changes the output to -2/5. Both components are individually seen
by C; their cancellation at one input does not make either a blind component.
The smallest closed observer has rank four, even though the initial readout
has rank one and the mixing graph has no link.

**W31, form boundary.** The connected two-carrier example reconstructs the
single form line lambda I_4. It does not recover the Lorentz-type 1+3 form of
the earlier real-K candidate family. This is a proved property of the selected
energy-preserving representation, not a physical signature-selection theorem.
A geometric target, if proposed, must have its own explicit adapter.

## 5. U29 — a minimal response transcript from one scalar readout

Assume Gamma is connected. Pick a root r and use the single scalar readout
ell extracting the first coordinate at that root. Fix a rooted spanning tree.
For each label a, take the ordered list of mixing links along the tree path
from r to a, and define the **row product**

\[
w_a=\ell\,U_{e_1}\cdots U_{e_d},\qquad
w_r=\ell.
\]

Include the two scalar readings

\[
y_a=w_a x,\qquad z_a=w_a V_a x,
\quad V_a=T_{t_a}(G_a).
\]

The chronological action on x is the reverse of the written row-product
order. In the second experiment V_a acts first. The implementation records
the chronological word explicitly and verifies it separately.

**Proposition U29.** These 2m scalar readings reconstruct every x in W by a
triangular decoder. Their transcript matrix has rank 2m, the minimum possible
number of readings in a fixed real scalar linear transcript recovering W.
Each protocol has at most d(r,a)+1 steps, hence at most m steps, where d is
distance in the chosen spanning tree. The choice of a breadth-first tree gives
shortest-path distances; no global optimality claim for other experiment costs
is made.

**Proof.** Write x_a=(p_a,q_a). Every w_a has support only on the first
coordinates of labels on its tree path. Its coefficient alpha_a at a is the
product of the signed, nonzero mixing coefficients s_h along that path, with
alpha_r=1. In increasing tree depth, the reading y_a therefore determines
p_a after subtracting the already recovered ancestor terms.

Only the a term changes under V_a. Put
c_a=(1-t_a^2)/(1+t_a^2), s_a=2t_a/(1+t_a^2). Then

\[
z_a-y_a=\alpha_a\big((c_a-1)p_a-s_aq_a\big),
\qquad
q_a=\frac{(c_a-1)\alpha_ap_a-(z_a-y_a)}{s_a\alpha_a}.
\]

Every denominator is nonzero. This reconstructs all 2m real coordinates and
proves full transcript rank. A smaller number of scalar linear readings has
rank below 2m and cannot be injective. The word-length bound follows directly
from the tree paths. ∎

This is an exact finite **model transcript**. Physical use requires a justified
preparation/repetition interface, known control parameters, and access to the
signed readout. It does not assume one can copy an unknown quantum state or
jointly measure incompatible observables. Noise and conditioning are separate
obligations; small nonzero parameters can make exact algebraic decoding poorly
conditioned in an experiment.

**W32, missing-control guards.** For two linked channels and the first scalar
readout, mixing alone sees only the two first coordinates, giving rank two.
Local R controls without mixing see only the observed channel, also rank two.
The complete connected family sees all four coordinates. The reconstruction
therefore depends on the declared controls, not just the names of the carriers.

## 6. Evidence, lineage and next physical target

The certificate compares the graph classification to an independent matrix
row-space closure that knows no graph. Full symmetric-form nullspaces are
computed independently on the smaller stated domain. Transcript decoders are
checked on every basis vector of each tested carrier, and chronological action
is checked against the predicted scalar readings on a separate signed probe.

| Declared finite domain | Coverage |
| --- | ---: |
| Every simple interaction graph on 1–4 retained channels | 75 |
| Connected graphs in that domain | 44 |
| Common-grading and orthogonality checks for finite steps | 490 |
| Preserved represented quarter-turn checks | 490 |
| Graph formula versus matrix readout-closure comparisons | 360 |
| Exact component-projection intertwiners | 700 |
| Minimal transcripts, every root of every connected graph | 167 |
| Exact decoder reconstructions on basis vectors | 1,298 |
| Chronological word checks on a signed probe | 1,298 |
| Full invariant-form nullspaces, all graphs on 1–3 channels | 11 |
| Additional regressions / package total | 16 / 93 |

The exhaustive graph campaign uses t=h=1/2. Additional regressions include
negative parameters, parameters equal to one, other roots, independent Cayley
inversion, and a symmetric-mixing rejection control. Earlier evidence remains
in the source-bound certificate. Written general proofs, these finite checks,
source provenance, formal verification and physical validation remain distinct.

Cayley unitarity, cross-channel commutators, quotient descent and linear
observability have explicit predecessors; local proposition numbers do not
claim priority for them. The concrete advance is a declared interacting
family with one native cut grading, classified exact cuts, a genuine lossy
noncommutative quotient and an explicit minimal response protocol.

The next physical target is to specify **which observable of this native
interaction is proposed to be gravitational**, justify its selection, and
construct its geometric or non-geometric readout. The selected family supplies
a preserved positive form and represented quarter-turn; it does not yet give
a Born law, a Lorentzian metric, a physical clock, a spacetime dimension or
quantum-gravity closure. Those require additional typed constructions rather
than relabelling the native loop commutator as physical curvature.
