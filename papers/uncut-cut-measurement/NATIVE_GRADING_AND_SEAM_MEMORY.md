# Native grading, gluing and minimum seam memory

**Monty Dabas — research continuation v0.6, 25 September 2026.**

The next selection step exposes a prior representation obligation: the native
EMK grading of a two-component carrier must survive its proposed assembly and
cut. U24 proves that the usual overlapping zero-extended K/R pairs cannot all
retain their native grading under one involution. U25 constructs a lawful
grading observer on a larger declared carrier and derives its exact minimum
linear memory. U26 shows why preserving that grading alone does not preserve
every separately addressed edge response, and gives explicit separating
transcripts. This step uses no spacetime or classical physical field equation.

## 1. Source audit: three structures that must stay distinct

The [four-file audit](GRADING_SOURCE_PINS.json) pins these native results:

- **N1, RKF 41:** every declared involution gives a unique even/odd grading.
  The theorem does not say that independently graded pair representations share
  the same involution after they are assembled on overlapping coordinates.
- **N2, RKF 48:** on the native two-component EMK carrier,

  \[
  K=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
  R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad J_E=K.
  \]

  K is J_E-even and R is J_E-odd. Its negative control explicitly rejects
  diag(1,-1) as a replacement yielding that same grading.
- **N3, RKF 55:** dagger grading is different again. R is a lawful native
  flow generator; real K alone is not the positive-energy dagger-unitary
  generator. K-sector Cayley maps below remain algebraic comparison candidates.
- **N4, RKF 31:** minimum observer rank is relative to a supplied memory
  target. Spectral II's minimum linear repair principle and RSC exact descent
  are also inherited; this note computes them for a specific assembly.

**W21.** The metric-sign involution H=diag(1,-1) makes K odd, since HKH=-K.
The native involution K makes K even, since KKK=K. This distinction is not
removed by saying both operators are involutions. A simultaneous change of
basis preserves the distinction in their action on the represented generators.

U18/U21's invariant-form calculations remain valid for their declared
candidate matrices. This extension closes an ambiguity in their interpretation:
those calculations did **not** derive a common native EMK grading for the
multi-pair assembly. A form's signature partition is not that derivation.

## 2. U24 — obstruction to a common grading on overlapping pair coordinates

Let a finite simple undirected graph have vertices {0,...,n-1}. For every edge
{i,j}, i<j, put K_ij=E_ij+E_ji and R_ij=-E_ij+E_ji, zero on all other
coordinates. These are the same pair embeddings used in the earlier candidate
construction; their squares on the full carrier are pair projections or their
negatives, not +/-I on all n components.

**Proposition U24.** There exists a single real involution J satisfying

\[
JK_{ij}J=K_{ij},\qquad JR_{ij}J=-R_{ij}
\]

for every edge if and only if the graph is a matching together with possible
isolated vertices. In particular, a connected graph with n>=3 admits no such
involution. More strongly, for a connected graph with n>=3 the full linear
system JK_ij=K_ij J and JR_ij=-R_ij J has only J=0, even before imposing J^2=I
or self-adjointness.

**Proof.** Commutation with K_ij implies commutation with
P_ij=K_ij^2, so J preserves the pair plane and its coordinate complement. On
that plane a matrix commuting with K has the form aI+bK. Anticommutation with R
forces a=0. Hence J e_i=b_ij e_j and J e_j=b_ij e_i. If edges {i,j} and {i,k}
overlap with j!=k, the two expressions for J e_i force both coefficients to
zero. The zero then propagates along every edge of that connected component.
On a connected nontrivial component with overlap, J vanishes identically.

An involution cannot vanish on such a component. Conversely, on disjoint pairs
choose J to swap the paired basis vectors, and on isolated coordinates choose
J=+I. This real symmetric involution has precisely the required grading. ∎

The restriction concerns **this embedding and these simultaneous pure-grading
requirements**. It does not forbid larger native representations, a different
assembly, or a chosen global grading that mixes the formerly pure pair sectors.
The one-involution claim must be proved for whichever alternative is proposed.

## 3. A larger carrier that retains the local native grading

From now on let the graph be connected, n>=2, with m edges. Introduce a separate
two-component carrier for each edge:

\[
W=\bigoplus_{e\in E}\mathbb R^2_e,\qquad \dim W=2m.
\]

For e={i,j}, write its two basis vectors p_(e,i),p_(e,j). Coordinates retain
both the endpoint and the edge identity. Let J_0 swap these two vectors on
every edge. Thus J_0 is the direct sum of the native K involutions. Local
K_e and R_e now have disjoint support and obey

\[
J_0K_eJ_0=K_e,\qquad J_0R_eJ_0=-R_e.
\]

Define the **declared equal-weight endpoint readout** C:W -> R^n by

\[
Cp_{(e,i)}=e_i.
\]

C sums all copies of a vertex and forgets which incident edge carried them.
Connectedness gives rank C=n. This direct-sum presentation and readout are
explicit model data, not a derivation of the uncut primitive or a physical
many-body interaction. C is a surjective linear readout, not an involution.

**W22.** On edges (0,1),(0,2), use port order (01:0,01:1,02:0,02:2). Then

\[
z=(1,0,-1,0),\qquad Cz=0,\qquad CJ_0z=(0,1,-1).
\]

The zero state and z have the same present shadow but different shadows after
grading. Thus J_0 does not descend through C. More generally, every overlapping
pair of edges gives this witness, so C admits exact grading descent only in
the connected two-vertex case. This is a concrete failure of the inherited
kernel-invariance gate, not a failure of the local native EMK algebra.

## 4. U25 — exact minimum memory for the grading target

The target here is to retain C and all its continuations under the single
involution J_0. Since J_0^2=I, the complete target is

\[
D=\begin{pmatrix}C\\CJ_0\end{pmatrix}.
\]

Its codomain is **im D**, not all of R^(2n). On that image, J_0 descends by
swapping the two displayed output blocks: D J_0=S D with S(y,z)=(z,y).
The swap preserves im D.

**Proposition U25.** Among full-row-rank linear observers A retaining the
initial readout (C=L A) and admitting an exact grading intertwiner
A J_0=J_A A, the minimum observer rank is rank D. The minimum number of
additional independent scalar channels beyond C is

\[
\boxed{
r_J=\operatorname{rank}D-n
=\begin{cases}
n-2,&\text{graph bipartite},\\
n-1,&\text{graph non-bipartite}.
\end{cases}}
\]

Equivalently, in the supplied Euclidean port presentation define

\[
Q=I-C^T(CC^T)^{-1}C,\qquad M_J=CJ_0Q.
\]

Then Q projects onto ker C and r_J=rank M_J: the additional memory is exactly
the dimension of initially hidden information returned to the visible readout
by J_0. This is a linear channel count, not a bit count, thermodynamic entropy,
or a derived action law for nature.

**Proof of minimality.** For any admissible A,
CJ_0=LJ_A A. Hence D factors through A, and rank A>=rank D. Conversely, choose
a row basis of D containing the n independent rows of C. The resulting A has
rank D and contains the initial readout. Its row space is J_0-invariant, so an
exact J_A exists; J_A^2=I follows from A's full row rank and J_0^2=I. Thus the
bound is attained with rank D-n additional rows. The original row space of C
is orthogonal to ker C. Projecting the additional rows CJ_0 onto that kernel
gives M_J, so rank D=n+rank M_J.

**Proof of the graph formula.** A left-null vector (alpha,beta) of D satisfies,
for every edge {i,j},

\[
\alpha_i+\beta_j=0,\qquad \alpha_j+\beta_i=0.
\]

Writing a=alpha+beta and b=alpha-beta gives a_i=-a_j and b_i=b_j on each
edge. Connectedness leaves one constant degree of freedom for b. The alternating
assignment a has one degree of freedom on a bipartite graph; an odd cycle
forces a=0 on a non-bipartite connected graph. The left nullity is therefore
two or one, respectively. Thus rank D=2n-2 or 2n-1, proving the formula. ∎

The generic target-factorization/rank principle is inherited from the Spectral
and native observer results. The graph formula computes that principle for
this particular endpoint assembly. It is elementary linear algebra, not a
claim of global mathematical priority.

**W23.** A three-vertex tree has source rank 4, initial readout rank 3, and needs
one additional scalar channel. Any tree with n vertices has 2m=2n-2, so its
minimal grading observer is faithful on the whole source. For a four-vertex
star or path it has rank 6 and needs two extra channels. The repair does not
select the star over the path.

## 5. U26 — the admitted experiment family determines the memory requirement

Now enlarge the target to include a separately addressed step on **each** edge.
For each edge choose a known nonzero real t_e and the local R-sector step

\[
U_e=(I+t_eR_e)(I-t_eR_e)^{-1}.
\]

U_e acts as the identity on all other edge carriers. R is a licensed native
flow sector by N3, although physical implementation of these particular
controls and readouts still needs its own adapter. Define the model transcript

\[
\mathcal T x=(Cx,(CU_ex)_{e\in E}).
\]

All components refer to the same input x. Treating this stack as an actual
experiment requires justified preparation/repetition and access to the signed
readouts; neither simultaneous quantum measurability nor arbitrary copying of
an unknown state is asserted here.

**Proposition U26.** The transcript map T is injective. Therefore any linear
observer retaining C and descending all the U_e must be faithful on W. The
minimum additional linear channel count is

\[
\boxed{r_{\rm edge}=2m-n.}
\]

The result also holds algebraically if any edge step is replaced by its K-sector
Cayley candidate with nonzero t_e!=+/-1; that substitution does not upgrade K
to a licensed native dagger-unitary flow.

**Proof.** If T x=0, then Cx=0 and C(U_e-I)x=0 for every edge. On edge e={i,j},
C embeds its two coordinates as the distinct visible coordinates i,j. The
two-dimensional matrix U_e-I is invertible:

\[
U_e-I=2t_e R(I-t_eR)^{-1}
\]

on that edge. Hence the two coordinates x_e vanish. Doing this for every edge
gives x=0. More explicitly, the difference CU_ex-Cx, restricted to coordinates
i,j, reconstructs x_e by that two-dimensional inverse.

If C=L A and A U_e=V_e A, then every transcript component factors through A.
Injectivity of T forces rank A=2m. This is attainable by retaining the full
source, extending the independent rows of C to a basis with 2m-n additional
rows. The K candidate argument uses the same invertibility identity. ∎

This exposes a practical target distinction: preserving the grading operation
does not in general preserve the full declared action catalogue.

| Declared assembly | Source rank | Initial C rank | Minimum grading repair | Minimum addressed-edge repair |
| --- | ---: | ---: | ---: | ---: |
| Three-vertex tree | 4 | 3 | 1 | 1 |
| Four-vertex star or path | 6 | 4 | 2 | 2 |
| Triangle | 6 | 3 | 2 | 3 |
| Four-cycle | 8 | 4 | 2 | 4 |

## 6. Separating transcripts and the curvature guard

**W24/W25.** For triangle edges (0,1),(0,2),(1,2), order the endpoint copies as
(01:0,01:1,02:0,02:2,12:1,12:2). Take

\[
z=(1,-1,-1,1,1,-1).
\]

Both Cz=0 and CJ_0z=0, so z is invisible to the entire grading-only target.
An addressed (0,1) step with t=1/2 distinguishes it from zero. For the two
declared algebraic candidates the signed outputs are

\[
CU_{01}^{K}z=(-2/3,\,2/3,\,0),\qquad
CU_{01}^{R}z=(2/5,\,6/5,\,0).
\]

Their maximum coordinate separation is 16/15. If the input, step parameter and
readout calibration are fixed as stated, this is a model-level test that can
reject one candidate. It is not an observed gravitational measurement or a
physical derivation selecting a sector. The R branch alone already establishes
the extra information requirement without using a real K flow.

**W26.** The larger direct-sum carrier has a further consequence: steps on
distinct edges commute. If a surjective readout C genuinely descended two such
steps, their descended steps would also commute, since

\[
[V_e,V_f]C=C[U_e,U_f]=0.
\]

Surjectivity would give [V_e,V_f]=0. On the three-vertex star, however, the
averaging section S=C^T(CC^T)^(-1) gives compressed steps CU_eS and CU_fS with
a nonzero commutator; both exact descent tests reject them. This is an EMK-edge
instance of U19/W13's existing false-commutator control.

Thus retaining local native grading by a direct sum has **not** supplied the
interacting noncommutative transport of the earlier overlapping star. A physical
bridge needs a native inter-edge transport and its lawful readout. Compression
cannot be used to manufacture that missing interaction or call its artifact
gravitational curvature.

## 7. Evidence and the next native construction

The implementation uses exact rational arithmetic. The graph-rank test builds
C and J_0 from raw endpoint coordinates and compares rational elimination with
the combinatorial formula; it does not feed the expected rank to the solver.
The common-grader oracle solves every matrix entry with no assumed symmetry,
diagonal form or involution. Local response inverses are separately reconstructed.

| Domain | Coverage |
| --- | ---: |
| All labelled simple graphs on 2–5 vertices | 1,098 |
| Connected graphs checked for the rank formula | 771 |
| Bipartite / non-bipartite connected cases | 218 / 553 |
| Connected graphs on 2–4 vertices: full common-grader equations | 43 |
| Minimal grading observers and exact intertwiners constructed | 43 |
| Independent cut-return memory-rank checks | 43 |
| Separately addressed R families checked for full source recovery | 43 |
| Local K/R grading identities | 308 |
| Exact endpoint reconstructions from response differences | 154 |
| Additional regressions / package total | 14 / 77 |

The written proofs cover their stated finite-dimensional hypotheses; these
enumerations provide finite computational evidence, not proof-assistant
verification or physical validation. The certificate binds this manuscript,
its source pins, implementation and tests, alongside all earlier modules.

The next construction must provide **native interaction between the retained
edge carriers**, with its common grading specified and its experimental target
declared. It must then prove a lawful cut to any proposed metric carrier and
compare predictions on shared readouts. Native grading repair, physical sector
selection, spacetime emergence and quantum-gravity closure are distinct tasks.
This step resolves the grading/memory problem for the declared assembly and
identifies exactly what its direct-sum repair still cannot provide.

**v0.7 continuation:** [U27–U29](NATIVE_INTERACTION_AND_LAWFUL_CUTS.md) now
construct a declared native-compatible inter-edge mixing family, classify its
exact real-linear cuts, and give a proper lossy quotient preserving nonzero
interaction. The coupling data and physical gravitational interpretation
remain separate selection obligations.
