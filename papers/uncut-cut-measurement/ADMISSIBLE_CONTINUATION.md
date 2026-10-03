# Admissible continuation and the limits of a shadow quotient

**Monty Dabas — v0.3 extension, 25 September 2026.**

This module extends the finite presentation of U1–U12 to partial arrows. It
does not supply a physical law deciding which arrows exist. The native RSC
requirement to recover the middle identity before composition already appears
in sources S1 and S6. Here we specify one finite operational representation of
an admissibility gate and expose a failure of a tempting quotient definition.
See [novelty and lineage](NOVELTY_AND_LINEAGE.md) before interpreting the
proposition numbering as a claim of discovery.

## 1. Two different observation contracts

Let X be a nonempty finite presented carrier, |X| = N. Fix a cut c:X → A,
q = |c(X)|, and a finite labelled family of partial deterministic arrows
T_a:E_a → X, E_a ⊆ X. A word is admitted at x exactly when each successive
arrow is defined at the current state. Its length is a composition index,
not an assumed physical duration. The empty word is admitted everywhere.

For typed arrows one may use a disjoint union of the declared object carriers
and restrict each domain to its source object. Any type distinctions that the
observer must report must also be included in c. History-dependent permission
requires an explicitly justified state augmentation; it is not automatically
represented by this Markovian finite model.

We keep two contracts separate:

- **Domain-aware observation:** the observer must retain whether each labelled
  action is enabled, as well as the readout after every admitted continuation.
  Operationally using availability as data requires a justified interface
  that reports it without performing the forbidden operation. The theorem
  supplies no such instrument and assigns no physical cost to that query.
- **Common-domain experiments:** to compare two candidates, every applied
  action must be admitted at both. Only actual readouts distinguish them;
  unobserved availability is not an extra measurement.

An undefined action is neither a zero-valued measurement nor a permitted
failed test. If a real instrument returns a refusal record, that interaction
must be modelled as its own permitted observation, with any disturbance stated.

## 2. Proposition U13 — Exact partial descent

A partial observed arrow t_a:c(E_a) → c(X) reproduces both the enabled domain
and the next readout, independently of the representative, precisely when:

1. c(x)=c(y) implies x∈E_a iff y∈E_a;
2. c(x)=c(y) and x,y∈E_a imply c(T_a x)=c(T_a y).

The first condition says E_a is a union of complete c-fibres. Under both
conditions set t_a(c(x))=c(T_a x) for x∈E_a. This is well-defined, and its
domain membership agrees exactly with that of T_a. Conversely, exact domain
agreement implies condition 1 and a well-defined next readout implies 2.
The observed arrow is unique on its declared domain. **Proof complete.**

If every arrow satisfies the gate, every word descends with exactly its true
domain. Prove this by induction on word length: domain agreement licenses the
next arrow at both levels, and next-readout agreement preserves the induction
relation. Equal shadows alone cannot license a middle composition.

This is quotient factorization with an explicit domain condition, not a new
derivation of the full native recovered-identity structure.

## 3. Proposition U14 — Closure and its minimum memory

Introduce a fresh formal symbol ⊥ outside A, and define a mathematical record

\[
o_w(x)=\begin{cases}c(T_wx),&w\text{ admitted at }x,\\
\bot,&w\text{ not admitted at }x.\end{cases}
\]

The symbol denotes a domain fact in the model; it is not a result obtained by
executing an inadmissible word. Define x≡_D y iff o_w(x)=o_w(y) for every word.
This is an equivalence relation because it is equality of complete records.
Let H_D(c) be its partition.

Start P_0=c up to renaming, and repeatedly distinguish states by

\[
P_{k+1}(x)=\left(P_k(x),\big(s_{a,k}(x)\big)_a\right),\qquad
s_{a,k}(x)=\begin{cases}(1,P_k(T_ax)),&x\in E_a,\\(0),&x\notin E_a.
\end{cases}
\]

All tuples are equality signatures, not additional primitive physical objects.
The sequence stabilizes after at most N−q strict refinement rounds. Its stable
partition is H_D(c), the coarsest refinement of c satisfying U13 for all arrows.

**Proof.** Induction shows that P_k equality is exactly equality of the records
o_w for |w|≤k. Prepending one action either reports matching undefined domains
or compares the successor records. A fixed partition satisfies U13, hence
preserves all words by induction. Every strict step increases the class count
by at least one; there can be at most N−q such steps. If d refines c and satisfies
U13, then equality under d preserves every record o_w. Thus d refines H_D(c),
establishing coarseness as well as identification with the stable partition.

If arbitrary deterministic repair maps m:X → M are allowed, the minimum
alphabet size for a closed augmentation (c,m) is

\[
K_D(c)=\max_{a\in c(X)}\#\{H_D(c)(x):c(x)=a\}.
\]

Indeed, any closed augmentation refines H_D(c), so it needs at least this many
labels within each initial fibre. Assigning distinct local labels to its
H_D(c)-classes, reusing labels between different c-fibres, makes (c,m) have
**exactly** the H_D(c) partition. It is therefore closed and attains the bound.
Merely refining a closed partition does not itself guarantee closure; the
upper-bound construction uses equality of partitions. **Proof complete.**

The abstract cost is a worst-case label count; fixed binary storage needs
ceil(log2 K_D(c)) bits. Spectral II's scalar probe count and this alphabet size
are different quantities. A restricted instrument catalogue may not realize
the constructed map at all. No energy, entropy production or minimum action
law follows from this counting argument.

### Mathematical totalization and its operational limit

Adjoin one absorbing state ∗, give it a fresh output different from every real
readout, and send every undefined transition to ∗. All real outputs are tagged
to prevent a collision even when c itself uses zero, None or tuple values.
Two original states have the same totalized future-output record exactly when
they have the same partial record o_w: an undefined word reaches ∗, while an
admitted word reaches a real tagged readout. Hence restricting the total
family closure to X gives H_D(c).

This justifies a mathematical implementation comparison, not an experimental
permission. Reading the sink as an observed outcome without a refusal/availability
interface would change the observation contract.

## 4. Proposition U15 — Diagnostics under the stated contract

Search ordered candidate pairs breadth-first. In common-domain mode traverse
only arrows enabled at both states; stop at differing actual readouts. This
returns a shortest commonly admitted distinguishing word, or reports none.
There are at most N² pairs. Any shortest path to differing readouts has no
repeated pair, so its length is at most N²−1.

In domain-aware mode also allow a terminal availability query at a pair where
one candidate enables an arrow and the other does not. Return one of:

- `output`: a word admitted at both candidates with different final readouts;
- `availability`: a common admitted prefix and the next arrow whose enabled
  flags differ. The final arrow is queried, **not executed** at the forbidden
  candidate.

Length counts the common executed prefix plus one position for an availability
query, or just the output word length. Breadth-first search must queue such
queries at their proper next depth; returning one immediately could skip a
shorter output witness already in the queue.

**Proof.** In safe mode the pair graph has exactly the common admissible
transitions, so paths and executable comparison words coincide. Standard
breadth-first enumeration visits them in nondecreasing length. For domain-aware
mode use the tagged totalization graph, stopping as soon as exactly one
coordinate reaches the sink and reporting it as a query. Both coordinates
have a common valid prefix before that first mismatch. If a shortest record
distinction were undefined for both, it would not distinguish; if only one
were undefined earlier, the earlier mismatch would be a shorter witness.
Thus the two reported types exhaust possible shortest record distinctions.
U14 proves that distinct H_D(c)-classes already differ at depth ≤N−q, so a
domain-aware witness exists within that bound. Same-class candidates have none.
The totalized search visits at most (N+1)² pairs. **Proof complete.**

These are model predictions for two candidates. If a larger set of candidates
remains possible, an executable experiment must be admitted at every remaining
candidate. Pairwise diagnostics alone do not solve that identification problem.

## 5. Proposition U16 — Common-domain agreement need not form a quotient

Define x∼_S y when all words admitted at both have equal final c-readouts.
This relation is reflexive and symmetric, but need not be transitive.

**W7.** Let X={0,1,2,3}, c=(0,0,0,1), with one arrow

\[
T=(0,\mathrm{undefined},3,3).
\]

For (0,1) and (1,2), only the empty word is admitted at both, and both initial
readouts are zero. Thus 0∼_S1 and 1∼_S2. But the word T is admitted at 0 and 2
and gives readouts 0 and 1. Thus 0 is not ∼_S-related to 2. In particular there
is no map h whose fibre equality is exactly ∼_S in this example: equality of
h-values would force the missing transitive relation. **Proof complete.**

Taking transitive closure would erase the valid 0-versus-2 experiment; requiring
domain equality instead would introduce availability distinctions. Either is
a changed contract and must be declared. This does not prohibit all useful
observers or quotients; it prohibits this proposed exact fibre characterization.

For a fixed pair, deterministic adaptive choice from observed input/readout
history cannot bypass absence of a safe word. Until the first different readout,
both candidates have the same history, so the policy chooses the same action.
If every chosen action is allowed at both, that branch is a common admitted
word. A first distinguishing readout would therefore provide precisely a word
excluded by x∼_S y. Conversely, a safe distinguishing word gives a fixed policy.
This statement assumes no hidden-state access or additional availability data;
it is not a general theorem about all physical adaptive instruments.

## 6. Further separating controls

**W8 — availability is not execution.** With c=(0,0) and T=(0,undefined), no
safe output word separates the states. An availability query separates them
immediately, and the domain-aware closure needs two memory labels. Replacing T
by total identity needs only one. Thus restricting an arrow's domain can
increase domain-aware memory: U9's monotonicity for adding labelled arrows
must not be reinterpreted as monotonicity under domain restriction. An undefined
value replaced by a self-loop or by measured zero would conceal this distinction.

**W9 — delayed availability.** With c=(0,0,0) and T=(1,2,undefined), refinement
has class counts 1,2,3 and depth 2=N−q. States 0 and 1 share a valid first step;
after it, their next enabled flags differ. The shortest domain-aware witness
is one applied T followed by an availability query for T. All actual readouts
are zero, so no pair has a safe output distinction. The domain-aware observer
nevertheless requires three memory labels within the initial fibre.

## 7. What has and has not advanced

The added value in this programme is a precise distinction between two
observation contracts, an executable composition gate, and controls that reject
false diagnostic experiments and an invalid quotient. Partial-transition
minimization itself is established mathematics; these are self-contained
finite formulations, not priority claims. The RSC source already required
recovered middle identity; this module does not rediscover that principle.

An empirical adapter must now specify which domains are justified, what
reports availability, whether that report disturbs the state, and which repair
channels the instrument can realize. A model flag is not an observed seam
coordinate. The finite presentation still is not identified with uncut ground.
No quantum/classical transition, spacetime construction or complete quantum
gravity follows yet.

## 8. Evidence

`partial_certificate.py` checks all **20,646** cut/ordered-two-partial-arrow
cases on carriers of size one through three. Independent pairwise quotient
gates, candidate partitions, direct word signatures and backward pair-distance
relaxation cross-check the implementation. Coverage includes **369,944**
pair/mode diagnostic checks, **30,392** exclusions of smaller memory maps,
**432,376** short-word descent checks and **20,646** tagged-totalization
comparisons. W7 is an additional named four-state control, not exhaustive
four-state coverage. The package has **37 regression tests** overall.

Written proofs, exact finite enumeration, source hashes, operational
realizability and physical validation retain separate statuses in CLAIMS.md.
