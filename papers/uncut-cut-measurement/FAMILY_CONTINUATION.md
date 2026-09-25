# Lawful Continuation Under a Family of Cuts and Arrows

**Monty Dabas — v0.2 extension, 25 September 2026.**

This extends U1–U7 in [MANUSCRIPT.md](MANUSCRIPT.md). All carriers below are
declared finite presentations, not an identification of the uncut ground with
a set. The results added to this draft concern exact deterministic arrows. No external clock,
spacetime metric, quantum probability or classical field equation is a premise.

## 1. The composition question

An observer may be sufficient when one transition is repeatedly used and
insufficient when different admissible transitions are composed. The task is
to retain exactly the distinctions needed by **every admitted continuation**.
The family of arrows must be fixed before claiming that the observer closes.

Let X be nonempty and finite, c: X -> A an initial cut, and

$$
\mathcal F=(T_a:X\to X)_{a\in I}
$$

a finite labelled family of total deterministic maps. The family may be empty.
Every finite word in these generators is admitted in this version. State- or
history-dependent restrictions on admissible words are a further problem.
The arrow label is known when executing or predicting an experiment; hidden
choices require a different, possibly set-valued, observer law.

A word w=(a_1,...,a_k) is written in execution order:

$$
T_w=T_{a_k}\circ\cdots\circ T_{a_1},\qquad T_{\varnothing}=\mathrm{id}_X.
$$

The length of w counts compositions and does not assign a physical duration.
Define continuation equivalence by

$$
x\equiv_{\mathcal F,c}y
\quad\Longleftrightarrow\quad
c(T_wx)=c(T_wy)\text{ for every finite word }w.
$$

The empty word makes this relation finer than equality of initial readouts.
Its classes are the distinctions relevant to the whole declared continuation
family. Denote the class map by H_F(c), with arbitrary class labels.

## 2. Constructive family closure

Start with p_0=c up to renaming its labels and recursively form the partition
whose signature at x is

$$
p_{n+1}(x)=\bigl(p_n(x),(p_n(T_ax))_{a\in I}\bigr).
$$

Only equality of these tuples is used. Implementations canonically rename the
resulting classes; numerical label values have no physical meaning.

### Proposition U8 — Coarsest closure for an admitted family

Let N=|X| and q=|c(X)|. The recursion stabilizes at a depth n <= N-q. At
stability its partition equals H_F(c). Every T_a induces a unique transition
on these classes. Any cut d refining c on which all T_a descend also refines
H_F(c). Thus H_F(c) is the coarsest common dynamically closed refinement.

**Proof.** Inductively, p_n(x)=p_n(y) exactly when all observations along
words of length at most n agree. For n+1, the displayed signature compares
the initial histories and every one-generator prefix followed by a word of
length at most n. These are precisely the words through depth n+1.

Each strict refinement increases the number of classes. Starting with q
classes, there can be at most N-q strict increases. If p_(n+1) and p_n have
the same partition, equality of p_n implies equality of p_n after every
T_a. U5 therefore gives the induced transitions. By induction on word length,
all subsequent readouts agree within every stable class. This proves equality
with continuation equivalence and prevents later splitting.

Finally, if d refines c and all arrows descend to d, equal d-readouts remain
equal after each word. Hence all c-readouts after words agree, so equal
d-values imply equal H_F(c)-values. U1 gives the factorization through d. □

The construction is relative to the declared X and F. An empty family gives
H_empty(c)=c. A one-arrow family has the partition already proved in U6.
No assumption says that all possibilities or physically allowed arrows have
been listed merely because an enumeration terminates.

## 3. Exact memory required by the family

Let F_c(a) denote an initial readout fibre as in U1. Define

$$
K_{\mathcal F}(c)=
\max_{a\in c(X)}\left|H_{\mathcal F}(c)(F_c(a))\right|.
$$

This counts continuation-distinct classes hidden behind the same initial
readout. A repair is a deterministic channel m on X, and the joint observer
d=(c,m) must admit an induced map for every arrow in F.

### Proposition U9 — Minimum family memory and added-arrow monotonicity

With arbitrary deterministic repair channels allowed, the minimum repair
alphabet has size K_F(c). The minimum fixed binary length is
ceil(log2 K_F(c)). If F is a subfamily of G on the same X, then

$$
H_{\mathcal G}(c)\text{ refines }H_{\mathcal F}(c),
\qquad K_{\mathcal F}(c)\le K_{\mathcal G}(c).
$$

**Proof.** Any closed joint observer must recover H_F(c) by U8. Within an
initial c-fibre it must separate all the continuation classes, so U2 gives
the lower bound K_F(c). For attainment, label the continuation classes within
each initial fibre and reuse labels across fibres. The resulting joint
observer and H_F(c) factor through each other: the initial cut and this
constructed memory are functions of the continuation class, and the joint
readout decodes that class. The partitions are identical, so every arrow
descends. The binary length follows by counting labels.

Every F-word is a G-word. Equality after every G-word therefore implies
equality after every F-word. A finer partition has at least as many classes
inside each fixed c-fibre, yielding the inequality. □

**Acquisition scope.** An implementable repair channel must distinguish the
initial possibilities before information needed to do so becomes unavailable.
Naming m from a known model proves existence of a function, not an instrument
that can infer m from c alone. Once acquired, the sufficient observer state
can be updated from the known arrow label without reacquiring a full history.
Restricting available instruments can increase the cost or prevent repair.

This is the minimum memory needed to predict the declared family. It is not
a physical energy cost or a theorem selecting nature's minimum-action path.

## 4. Composition and the choice of generators

### Proposition U10 — Word descent and generator-presentation invariance

Let h=H_F(c), and let t_a be the induced map with h T_a=t_a h. Every finite
word descends by composition in the same execution order. If G is obtained
from F by adjoining identities, duplicate arrows, or maps already expressible
as F-words, then H_G(c) and H_F(c) have the same partition and memory cost.
More generally this holds whenever two families generate the same collection
of endomaps under finite composition and identity.

**Proof.** Repeated substitution of h T_a=t_a h gives
h T_w=t_w h. Therefore any observer closed under F is closed under every
map generated by F. Adjoining those maps changes none of the closed observers
refining c. Their coarsest member is unchanged by U8. The argument is symmetric
for two families generating the same collection. U9's class count is then
unchanged. □

Adding a composite as a named generator can shorten its word length and the
number of refinement rounds. The final partition and memory budget, rather
than generator-dependent depth, are the invariant objects here. Generator
order in a list changes no partition; changing execution order in a word can
change its observation.

### Counterexample W5 — Separately sufficient memories fail under mixing

Take X={0,1,2,3}, with the following completely declared data:

| x | c(x) | A(x) | B(x) |
| --- | --- | --- | --- |
| 0 | 0 | 0 | 0 |
| 1 | 0 | 2 | 0 |
| 2 | 0 | 2 | 3 |
| 3 | 1 | 3 | 3 |

Repeated A preserves c, so H_A(c) has classes {0,1,2}, {3}. Repeated B
gives H_B(c) with classes {0,1}, {2}, {3}. Their joint observation still
identifies 0 and 1. Neither a string of A's nor a string of B's can
distinguish these two candidates in c.

But the word **A followed by B** does distinguish them:

$$
c(B(A(0)))=0,\qquad c(B(A(1)))=1.
$$

The common family closure is the full four-class partition, reached with
class counts 2 -> 3 -> 4. Behind initial readout 0 there are now three
continuation classes. Thus

$$
K_{\{A\}}(c)=1,\quad K_{\{B\}}(c)=2,
\quad K_{\{A,B\}}(c)=3.
$$

The separate memory banks would yield two labels; the family requires three,
or two fixed binary bits. In particular there is no general bound of the form
K_F_union_G(c) <= K_F(c) K_G(c). This does not contradict U3, which bounds
a fixed target across a refinement. Here admitting mixed words changes the
target to be predicted.

The reversed word B followed by A gives readout 0 from both initial candidates.
The order effect is an exact observational statement in this finite model.
No probability, interference, physical curvature or quantum identification
follows just from this noncommutation.

## 5. Consistency between different initial cuts

For cuts c and d on the same X, their joint cut c join d records (c(x),d(x)).
Fix the **same** admissible family F throughout the next proposition.

### Proposition U11 — Closure of observations and joint cuts

H_F is extensive, monotone under refinement and idempotent as an operation
on partitions. It also preserves joint cuts:

$$
H_{\mathcal F}(c\mathbin{\vee}d)
\simeq H_{\mathcal F}(c)\mathbin{\vee}H_{\mathcal F}(d),
$$

where equivalence means equality of partitions, up to renaming labels. If
d refines c, the unique map r from H_F(d)-labels to H_F(c)-labels intertwines
every induced transition.

**Proof.** Extensivity and idempotence follow from U8. If d refines c, equality
of all d-readouts along words implies equality of the corresponding c-readouts,
which proves monotonicity. The joint of two F-closed partitions is F-closed:
equal joint labels have equal labels in each component after every arrow.
It refines c join d. Conversely, every F-closed refinement of c join d must
refine both H_F(c) and H_F(d) by U8. This gives the displayed equality.

For the intertwiner, evaluate on a label H_F(d)(x). Applying r after the
induced d-transition gives H_F(c)(T_a x), as does applying the induced
c-transition after r. Every d-label has a representative, so the maps agree. □

This proposition uses one common family. W5 joins closures computed for
different families, which need not already be closed under their union.

**Counterexample W6 — More initial detail can demand more continuation memory.**
On X={0,1,2}, let c be constant, d=(0,0,1), and T=(0,2,2). The constant
observation is closed with K_T(c)=1, while the more informative d requires
the full three-class closure and K_T(d)=2. A more detailed present observation
can ask a harder future-prediction question. Consequently U11's monotonicity
of partitions must not be misread as a general decrease of K_F with initial
refinement. U3's fixed-target information bound remains valid.

## 6. A concrete experiment distinguishing a hidden pair

### Proposition U12 — Shortest distinguishing continuation

If H_F(c)(x) differs from H_F(c)(y), there is a distinguishing word w with
length at most N-q. For candidates already distinguished by c, the empty word
suffices. For equivalent candidates no finite distinguishing word exists.

A breadth-first search on ordered pairs (u,v) with transitions
(u,v) -> (T_a u,T_a v) finds a shortest distinguishing word when one exists
and terminates after visiting at most N^2 pairs.

**Proof.** U8 identifies continuation classes by observations through a stable
depth n <= N-q. Distinct classes therefore differ on at least one word through
that depth. Conversely, equal classes agree on all words by definition.

Each vertex reached in the pair graph records exactly the endpoints of the
same word applied to both candidates. Breadth-first search explores shortest
paths in nondecreasing edge count. The first reached vertex with unequal
c-labels supplies a shortest distinguishing word. Revisiting a pair cannot
improve future reachability and is unnecessary. There are at most N^2 ordered
pairs, so exhaustion without unequal labels proves no word distinguishes. □

In W5, the shortest word for 0 and 1 is (A,B). This yields a model-level
prediction before a new measurement: apply A, then B, then read c. A physical
application must first justify that the two candidate preparations, arrows and
readout can be realized with a controlled error model. A mathematical witness
is not empirical validation of that adapter.

## 7. What this contributes to the uncut programme

The first version tracked distinctions lost by a single observation and a
single repeated continuation. This extension derives how much memory is needed
for an entire declared family and its compositions. It also supplies a finite
diagnostic sequence when a proposed identification loses a future distinction.

This makes one aspect of seam alignment precise: the retained observation and
memory must support compatible continuation across all declared arrows and
their composites. It does not infer a symmetry from transitivity or identify
the primitive seam with the resulting finite equivalence classes.

The next substantive question is whether different physically admissible cuts
and interventions admit a common operational presentation, including noise and
state disturbance. That question cannot be settled by assuming an exhaustive
X in advance. The quantum/classical and spacetime identifications remain
derivation targets, with no physical framework imported as an axiom here.

## 8. Evidence and provenance

U8–U12 are self-contained proofs from the finite definitions and U1–U6.
They are elementary closure and finite-observation results; no priority claim
is made for these general mathematical methods. The named W5/W6 calculations
are explicit controls for how the methods are applied in this programme.

The v0.3 [lineage audit](NOVELTY_AND_LINEAGE.md) records the exact Spectral I
predecessor of U1, the Spectral II minimum-repair analogue and established
finite-state refinement methods. Partial arrows require the distinct contracts
in [ADMISSIBLE_CONTINUATION.md](ADMISSIBLE_CONTINUATION.md); this module's total
arrow results do not authorize execution of a forbidden operation.

`family_certificate.py` exhaustively checks all ordered pairs of total arrows
and all cut partitions on carriers of size one through three. It compares
closure with direct enumeration of words and all candidate closed partitions,
and compares the diagnostic search with brute-force shortest words. Additional
four-state fixtures include W5, the U6 sharp-depth chain, and all identity
transitions on each initial partition. These additional fixtures are not
described as an exhaustive four-state two-arrow campaign.

`CERTIFICATE.json` binds this document, the implementations and tests. Its
status is PASS_FINITE_CHECKS. The general proofs, finite evidence, operational
instrument requirements and physical interpretation retain distinct scopes.
