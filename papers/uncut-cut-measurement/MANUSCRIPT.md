# Physics from Cuts: Uncut Ground, Measurement and Lawful Continuation

**Monty Dabas**  
Working manuscript v0.2 — 25 September 2026

## Abstract

We propose that an account of physics can begin with uncut ground and admissible
distinctions, with spacetime reserved for a derived representation. The proposal
applies to its own mathematical language: every formal description declares
distinctions and therefore has a scope of recognition. We develop a finite
operational layer in which this limitation can be stated and tested without
identifying a finite state set with the uncut. A cut retains an observation;
a target specifies what must be recovered. We prove exact recovery and minimum
memory criteria, a refinement bound, and a constructive theorem for the coarsest
refinement on which a declared clock-free transition closes. Counterexamples
separate agreement from completeness, pairwise compatibility from global
realizability, and continuation consistency from physical symmetry. The
recognition-seam composition law connects the finite discussion to the existing
native repository. The v0.2 companion extends closure to an entire declared
family of arrows, proves its exact memory requirement, and constructs shortest
distinguishing continuations. The resulting programme seeks laws compatible
across cuts; the physical derivation of quantum and classical descriptions
remains open.

## 1. The starting point includes its own limitation

The conceptual starting point is **uncut ground**. This phrase does not already
specify a set of individually labelled physical states, a manifold, a metric,
a probability distribution or a clock. Those would add distinctions.

The second commitment is to investigate **admissible distinction and
continuation**. A wholly unspecified primitive cannot select a particular law;
the admissibility rule must acquire mathematical content. Whether continuation
is primitive at the ground level or becomes meaningful only through a first
cut is an open foundational choice. The present results require only declared
transitions within a presented sector and do not settle that choice.

The reflexive principle is:

> A mathematical description must disclose the distinctions it retains and
> the distinctions it cannot resolve for its stated target.

This is a modelling commitment. It is not a theorem that all mathematics is
physically incomplete, nor a proof that a conceptual distinction causes a
physical disturbance. A descriptive cut and a physical intervention need
separate operational identifications.

The native source hierarchy [S1] proceeds from pre-distinction ground through
cut, trans-cut arrows and recognition before completion and operator
representations. Our finite layer is placed after an admitted presentation.
It does not reverse that hierarchy by making a finite set the absolute ground.

## 2. A finite presentation and its measurement record

Fix a nonempty finite set X of **presented possibilities**. Its labels and
equality relation are declared mathematical structure. X may itself be a coarse
presentation. No completeness claim about nature follows from enumerating X.

A cut is a function c: X -> A, with A taken to be its actual image. The
observation is a = c(x). The fibre

$$
F_c(a)=\{x\in X:c(x)=a\}
$$

is the collection of presented possibilities left indistinguishable. Define
x ~_c y iff c(x)=c(y). This is observational equality, not an assertion that
the possibilities or their continuations are identical.

A target is a predeclared function g: X -> V. Its ambiguity at a is

$$
G_c(a)=\{g(x):x\in F_c(a)\}.
$$

The mathematical measurement record is

$$
\mathcal R=(\mathrm{id}(X,c,g),a,G_c(a),\text{provenance}).
$$

Raw observations are retained. Model revisions append a new interpretation
with new provenance; they do not overwrite the observation. In an empirical
application, the true fibre or g may be unknown. A declared model may supply
an enclosure for G_c(a), conditional on its coverage assumptions. Enumerating
a model does not prove those assumptions. Instrument error and intervention
effects require their own model; this v0.1 concerns exact deterministic cuts.

The unresolved seam at a is represented here by the distinctions in F_c(a)
needed for g or for continuation. This is an operational use of the word
**seam**, not an identification of the primitive seam with a spatial boundary.

## 3. Recovery and the limit of correction

### Proposition U1 — Exact target descent

There is a function f: A -> V with g=f composed with c if and only if g is
constant on each fibre of c. When it exists, f is unique on A=c(X).

**Proof.** If g=f composed with c, equal c-values give equal g-values.
Conversely, assign f(a) the common g-value on the nonempty fibre F_c(a).
Constancy makes this well-defined, and surjectivity onto A gives uniqueness. □

**Correction consequence.** If c(x)=c(y) but g(x) differs from g(y), any
deterministic algorithm using only the readout, fixed metadata and the same
model must give the same output for both possibilities. It therefore cannot
recover g correctly for both. Metadata that differs between x and y supplies
additional observation and must be counted as a further channel.

An exact claim in this situation is the set G_c(a), not an arbitrarily chosen
member disguised as a unique correction. Non-singleton G_c(a) does not imply
that information has been destroyed in nature; it diagnoses this observation.
Conversely, a cut can be lossless for a particular g even when it loses other
distinctions. Correction need not be nonzero at every measurement.

**Witness W1.** Let X={(0,0),(0,1),(1,0),(1,1)}, c(u,v)=u and g(u,v)=v.
Every fibre has both target values. No function of u alone recovers v.
For the different target g'(u,v)=u, no added memory is required.

## 4. A finite cost of the cut

A repair channel m: X -> M is target-faithful when g factors through the joint
readout (c,m). The model permits arbitrary deterministic channels m. Restricting
the physically admissible channels can increase the minimum or prevent repair.

Define the worst-fibre target multiplicity

$$
K_g(c)=\max_{a\in A}|G_c(a)|.
$$

### Proposition U2 — Minimum repair alphabet

The least possible size of M for a target-faithful repair is exactly K_g(c).
Consequently a fixed-length binary encoding needs

$$
b_g(c)=\lceil\log_2 K_g(c)\rceil
$$

bits, where K=1 means zero added bits. This is a counting cost, not entropy,
energy, work or a thermodynamic conversion law.

**Proof.** Choose a fibre containing K_g(c) distinct target values and one
representative for each. Their memory labels must be different; otherwise
two equal joint readouts would have different targets, contradicting U1.
Thus |M| >= K_g(c). For each fibre, enumerate its distinct target values by
labels from {0,...,K_g(c)-1}. Assign x the label of g(x) in its own fibre.
The decoder knows a and hence its fibre-specific enumeration; it recovers g.
The largest fibre uses all K labels. A binary word of length b has 2^b labels,
giving the displayed minimum. □

The construction proves existence of a distinguishing channel. It does not
make g observable for free: a proposed instrument must actually implement m.
The linear target-relative observer result in [S3] measures the analogous
obstruction by rank on a blind subspace. U2 is a finite counting analogue;
rank and finite alphabet cardinality are not interchanged.

## 5. Refinement and accumulated loss

A cut d refines c when c=r composed with d for a well-defined readout map r.
Equivalently, equal d-values imply equal c-values. Refinement is a declared
comparison of descriptions, not elapsed physical time.

### Proposition U3 — Refinement and repair bounds

If d refines c, then for every b in d(X),

$$
G_d(b)\subseteq G_c(r(b)),\qquad K_g(d)\le K_g(c).
$$

Let K_d(c) count the largest number of distinct d-values within a c-fibre.
Then

$$
K_g(c)\le K_d(c)K_g(d).
$$

**Proof.** Each d-fibre is contained in the indicated c-fibre, which gives
inclusion and monotonicity. A c-fibre splits into at most K_d(c) d-fibres.
Each contributes at most K_g(d) distinct target values. Counting the union
gives the product upper bound. □

The bound need not be equality. On the four possibilities (u,v), let c be
constant, d=u and g=v. Then K_g(c)=2 while K_d(c)K_g(d)=2*2=4. Target values
can repeat across the intermediate fibres. Blind addition or multiplication
of correction costs can overcount the same missing distinction.

## 6. Compatible descriptions without an absolute viewpoint

For several cuts c_i of the same declared X, a tuple of readouts (a_i) is
globally realizable exactly when

$$
\bigcap_i F_{c_i}(a_i)\ne\varnothing.
$$

This is consistency relative to the presentation X. For independently supplied
presentations, a joint comparison carrier and its coverage must first be
justified; a universal carrier is not supplied by the word "uncut".

**Witness W2 — Pairwise agreement can fail globally.** Take
X={000,011,101,110}, with three cuts reading the three coordinates. Proposed
readouts (1,1,1) are realizable for every pair of cuts: 110, 101 and 011 supply
the three pairwise witnesses. No member of X realizes the whole tuple.
The overlaps therefore do not justify inventing a globally compatible state.

### Proposition U4 — Common target descent

Given a finite family of cuts on X, form an undirected graph joining x and y
whenever some cut identifies them. A target g factors through every cut if
and only if g is constant on every connected component of this graph.

**Proof.** Factorization through each cut makes g equal at the endpoints of
every graph edge by U1, hence along paths. Conversely, each fibre lies within
one component, so componentwise constancy gives every factorization by U1. □

**Witness W3 — Common scalar invariants may be too weak.** On all four pairs
(u,v), take cuts u and v. The graph is connected, so a target accessible
separately through both must be constant. Nevertheless the joint readout (u,v)
distinguishes every presented possibility.

Thus a physical law need not be a scalar readable through every individual
cut. A stronger candidate is a family of compatible relations or transition
laws, with declared comparison maps. No universal nontrivial invariant follows
from the mere existence of multiple descriptions.

## 7. Lawful continuation without an external clock

Declare a deterministic transition T: X -> X. The arrow supplies successor
structure; no duration, external time parameter, distance or spacetime metric
is assumed. Calling an arrow clock-free does not derive the arrow's law.

### Proposition U5 — Closed observed dynamics

There is a unique observed transition t: A -> A satisfying

$$
c\circ T=t\circ c
$$

if and only if c(x)=c(y) implies c(Tx)=c(Ty).

**Proof.** Apply U1 to target g=c composed with T. □

When this condition fails, a readout a defines only the successor set
{c(Tx):x in F_c(a)}. Choosing one successor without additional information
is a new assumption. A cut supplies an initial observed condition, which
need not specify a unique full presented condition.

If d refines c through r and T descends to both cuts, their observed laws
t_d and t_c obey r composed with t_d = t_c composed with r. To prove this,
evaluate on b=d(x): both sides equal c(Tx). Every b has such a representative.
This commuting relation is a precise notion of law compatible across cuts.

## 8. Constructing the least memory needed for continuation

Define the finite observation histories

$$
h_n(x)=(c(x),c(Tx),\ldots,c(T^n x)),\qquad E_n=\ker h_n.
$$

The integer n counts composition depth. It is not a spacetime coordinate.
Computing these histories uses the declared model T; it does not give an
experimental observer advance access to future data.

### Proposition U6 — Coarsest dynamically closed refinement

Let N=|X| and q=|c(X)|. The sequence E_n stabilizes at an index n <= N-q.
At its first stable index, h_n refines c and T descends to h_n. Moreover,
every cut d refining c on which T descends also refines h_n. Hence h_n is
the coarsest such refinement, up to renaming its labels.

**Proof.** E_(n+1) refines E_n. Each strict change increases the number of
classes by at least one, starting at q and never exceeding N. There can
be at most N-q strict changes. At equality E_(n+1)=E_n, equal histories
through n have equal histories through n+1; removing their first coordinates
shows h_n(Tx)=h_n(Ty). U5 gives descent on h_n. This stability also prevents
any later splitting, by induction under T.

For minimality, let d refine c and let T descend to d. If d(x)=d(y), repeated
use of its observed transition gives d(T^k x)=d(T^k y) for every k>=0.
Since c factors through d, all corresponding c-values agree. Thus equal
d-values imply equal h_n-values; U1 supplies h_n=s composed with d. □

Applying U2 with target h_n gives the minimum memory alphabet needed in
addition to c for any dynamically closed refinement. For the lower bound,
every such refinement recovers h_n by U6 and must separate its values inside
each c-fibre. For attainment, the U2 construction labels the h_n-values within
each fibre. Since c already factors through h_n, this joint readout and h_n
factor through each other and have the same partition, hence the same closure.
This is a representation minimum under the declared finite transition, not
proof that nature selects a minimum-cost physical trajectory.

**Witness W4 — Delayed memory matters.** On X={0,1,2,3}, set
c=(0,0,0,1) and T=(1,2,3,3). At depth 0 there are two classes, at depth 1
three, and at depth 2 four. Depth 2 is stable and meets the bound N-q=2.
The initial readout 0 conceals three distinct eventual histories, so three
memory labels, or two fixed-length bits, are necessary and sufficient.

For infinitely many possibilities, the finite termination argument does not
apply. A finite family of deterministic arrows is treated in
[U8–U12](FAMILY_CONTINUATION.md). Nondeterminism, noisy observations, restricted
word admissibility and operationally accessible memory require further work.

## 9. Connection to the native seam composition identity

After a linear representation has been declared, let P_x be recognized
projections, Q_x=I-P_x, and U_gamma event transports with lawful composition.
The inherited definitions [S2] are

$$
R_\gamma=P_yU_\gamma P_x,\quad
J_\gamma=P_yU_\gamma Q_x,\quad
C_\gamma=Q_yU_\gamma P_x.
$$

### Proposition U7 — Represented composition residue (inherited)

For composable arrows gamma:x->y and delta:y->z,

$$
R_{\delta\gamma}-R_\delta R_\gamma=J_\delta C_\gamma.
$$

**Proof.** Insert I=P_y+Q_y into
P_z U_delta U_gamma P_x. The two resulting terms are exactly the two terms
on the right of R_delta_gamma=R_delta R_gamma+J_delta C_gamma. □

The return of an intermediate memory channel can therefore change the visible
composite. In the scalar 2-by-2 block case, U=V=[[0,1],[1,0]] gives visible
one-step corners zero but the visible corner of VU equal to one. Removing
the memory term predicts the wrong composite.

The cut involution P-Q squares to I and is reversible in this representation;
the observation P alone can be many-to-one. They are distinct operations.
The finite cut model and the linear representation share an observability
question, but an adapter identifying their memories must be constructed for
each application. Hilbert structure is not silently inserted as a primitive
of sections 1–8.

## 10. Alignment, symmetry and the quantum/classical target

Compatibility under composition means that successive descriptions fit their
declared comparison rules. A symmetry additionally needs specified invertible
transformations preserving the relevant laws (or specified covariant relabelling
rules). Transitivity alone supplies no such transformation. The three-state
chain T=(0,0,1) has associative iterated continuation but only the identity
permutation commuting with T: its unique fixed point and the successive
preimage structure distinguish every state.

The proposed physical identification is that quantum and classical descriptions
arise as different presentations of one native continuation structure.
"Before" and "after" a cut initially denote levels of representation. No
physical clock order is implied. A cut by itself does not establish classical
equations, and retaining memory by itself does not establish quantum theory.
The probability law, interference, observable algebra, empirical adapter and
classical regime need derivations. Existing quantum witness capsules are
comparisons, not hidden quantum axioms of this package.

Similarly, a spacetime adapter would have to construct and validate whichever
notions of locality, order, dimension, metric and physical clocks it uses.
Neither a finite arrow count nor K_g(c) is automatically proper time or energy.
RH and Yang–Mills papers are not premises for the finite propositions here.
No open endpoint in those programmes is treated as solved.

## 11. Reflexive correction as a versioned programme

A correction is itself a declared map within a presentation. It therefore
inherits a target and scope rather than escaping the cut principle. Further
observations may justify a refinement; new hypotheses may instead revise X or
T. The latter is a model change, not automatically an information-preserving
refinement. Keeping those operations distinct makes disagreement informative.

Certification is also scoped: a written proof establishes an implication from
its hypotheses; exact enumeration checks declared finite instances; hashes
identify the checked source bytes. None establishes that the chosen physical
presentation is exhaustive. The first test of this programme is whether it
can predict, before additional observations, which missing distinctions matter
and what additional observation will repair them.

## 12. Sources and provenance

Only primary sources from this repository programme are used. Exact commits,
paths and content digests appear in [SOURCE_PINS.json](SOURCE_PINS.json).

- **[S1]** Recognition-Kernel-Framework, `theorum/27_rsc_primitive_to_completion_hierarchy.md`:
  pre-distinction ground, recovered identity and the placement of completion.
- **[S2]** Recognition-Kernel-Framework, `theorum/24_clock_free_recognition_seam_cut_calculus.md`:
  declared arrows, cut corners and the composition-residue identity.
- **[S3]** Recognition-Kernel-Framework, `theorum/31_cut_variational_minimal_observer_theorem.md`:
  target-relative blindness and minimum linear repair. Its continuum and
  spectral statements are not premises of U1–U6.

All three are pinned to correction commit
`927cdb6ca98221c0b4285da953a2c8b689fb202c`. U1–U6 have self-contained proofs
using the finite definitions above. U7 restates and proves the inherited
composition identity. The exact certificate checks finite instances and
counterexamples; it does not convert the physical proposals into theorems.

## 13. Family continuation extension (v0.2)

[FAMILY_CONTINUATION.md](FAMILY_CONTINUATION.md) contains U8–U12 with complete
proofs. It derives a common closed observer and its exact minimum memory for
all words in a declared family of finite deterministic arrows. It also proves
compatibility under joint observation and invariance under adjoining an
already generated composite, and supplies a shortest distinguishing word
when two candidates have different future observations.

The new W5 example separates per-arrow closure from closure under mixed words:
the sequence A then B requires distinctions absent from both unmixed histories.
W6 prevents the fixed-target refinement bound from being misapplied when the
future-prediction target itself changes. All of these statements remain
relative to a declared presentation and operationally available arrows.
