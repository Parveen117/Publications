# Claim ledger — v0.4

The word **proof** below refers to the written argument under the manuscript's
declared assumptions. **PASS_FINITE_CHECKS** refers only to the finite domains
recorded in `CERTIFICATE.json`. Neither status denotes proof-assistant checking,
independent review or physical validation.

See [NOVELTY_AND_LINEAGE.md](NOVELTY_AND_LINEAGE.md): numbering is local to
this draft. U1 is inherited exactly from Spectral I Theorem 18.1; minimum repair
has a Spectral II linear predecessor; quotient refinement is established
mathematics. Written proofs here do not imply priority.

| ID | Statement | Assumptions | Status and evidence |
| --- | --- | --- | --- |
| P0 | Uncut ground precedes formal distinctions; the framework is itself a cut | Proposed interpretation and admissibility programme | PROPOSED; manuscript §§1–2 |
| U1 | Exact target recovery iff target is fibre-constant; readout-only correction cannot separate a merged pair | Nonempty finite presented X, deterministic c and g | INHERITED_STATEMENT from Spectral I Theorem 18.1; self-contained proof and exhaustive decoder comparison through four states |
| U2 | Minimum repair alphabet is K_g(c); minimum fixed binary length is ceil(log2 K) | Arbitrary deterministic repair channels allowed; cost is worst-case label count | WRITTEN_PROOF; constructive decoder and exhaustive exclusion of smaller alphabets through four states |
| U3 | Refinement reduces ambiguity and K_g(c) <= K_d(c) K_g(d) | c factors through d on the same X | WRITTEN_PROOF; exhaustive partition triples through four states; strict inequality control |
| U4 | A common target factors through every cut iff constant on generated components | Finite family on a declared joint presentation | WRITTEN_PROOF; exhaustive pairs of cuts and target equality patterns through four states |
| U5 | Observed deterministic continuation exists iff next readouts are fibre-constant; lawful refinements intertwine | Declared total T on X; comparison maps well-defined | WRITTEN_PROOF; every endomap and cut through four states, plus every applicable refinement |
| U6 | History refinement terminates by depth N-q and is the coarsest dynamically closed refinement | N presented states, q initial readouts; one deterministic T; unrestricted access to model signatures for the construction | WRITTEN_PROOF; exhaustive comparison against all closed refining partitions through four states |
| U7 | R_delta_gamma - R_delta R_gamma = J_delta C_gamma | Declared composable linear transports; complementary cut projections | INHERITED_WRITTEN_PROOF from RSC; all scalar 2x2 block pairs with entries -1,0,1 checked |
| U8 | Coarsest observer closed under every arrow in a declared family; termination by N-q | Finite nonempty X; finite family of total deterministic arrows; all words admitted | WRITTEN_PROOF in FAMILY_CONTINUATION.md; exhaustive pairs of arrows and cuts through three states; 17 extra four-state fixtures |
| U9 | Minimum family repair alphabet equals maximum continuation-class count within an initial fibre; adding arrows cannot reduce it | Fixed X and initial cut; arbitrary deterministic repair channels permitted | WRITTEN_PROOF; smaller memory alphabets exhaustively rejected on checked cases |
| U10 | All words descend; adjoining generated composites or identities preserves closure and memory cost | Same generated collection of endomaps; known arrow labels | WRITTEN_PROOF; generator and observed-composition checks; W5 separates invariant memory from generator-dependent depth |
| U11 | Family closure is extensive, monotone, idempotent and preserves joint cuts; refinement maps intertwine | Same arrow family and carrier for all compared cuts | WRITTEN_PROOF; exhaustive cut pairs for two-arrow families through three states |
| U12 | A shortest distinguishing word exists within depth N-q for inequivalent candidates; pair search terminates | Finite deterministic family; all words admissible and labels known | WRITTEN_PROOF; 33,206 candidate-pair checks against direct shortest-word enumeration |
| U13 | A partial arrow descends with exact domain iff enabledness and successor readouts are constant on initial fibres; words then descend | Declared finite partial deterministic family, exact domain preservation | WRITTEN_PROOF in ADMISSIBLE_CONTINUATION.md; independent pairwise quotient gates and 432,376 short-word descent checks |
| U14 | Coarsest domain-aware closure stabilizes within N-q; exact repair alphabet is the maximum count of its classes per initial fibre | Fixed partial family; arbitrary deterministic repair maps; enabled-domain facts included in target | WRITTEN_PROOF; 20,646 exhaustive cases, all closed refining partitions compared, 30,392 smaller memory maps rejected |
| U15 | Shortest common-domain output experiment, or shortest domain-aware output/query witness, with separate bounds | Two candidate states; availability witnesses require a justified observation interface | WRITTEN_PROOF; 369,944 pair/mode checks against backward distance relaxation; all prefixes validated |
| U16 | Common-domain agreement need not be transitive or equal any observer's fibre relation; deterministic adaptation cannot separate an agreeing pair using only safe readouts | Partial family; deterministic policies use observed history, without hidden-state or availability access | WRITTEN_PROOF and W7 exact counterexample; adaptive statement is a general written argument, not exhaustive enumeration of policies |
| U17 | Compatible bilinear form fields correspond exactly to root forms fixed by all rooted loop comparisons | Connected graph; rank-n carriers; invertible typed linear transports | WRITTEN_PROOF in GRAVITY_BEFORE_CURVATURE.md; 48 integer single-arrow cases, 1,296 candidate forms, 27 edge checks, 12 frame changes |
| U18 | Complete chosen K-family selects a Lorentz-type form line; complete R-family selects a definite form line; their simultaneous invariant form is zero | Declared overlapping (0,i) EMK embeddings; nonzero nonsingular Cayley parameters; n is supplied | WRITTEN_PROOF; 24 K families in dimensions 2–5, 20 sector comparisons, 100 independent Cayley inversion checks; no physical signature-selection claim |
| U19 | A full bilinear form descends through a linear cut iff its kernel lies in the form radical; transport requires an exact intertwiner | Surjective linear cut; symmetric form; stated target is full bilinear pairing | WRITTEN_PROOF; eight cut-form checks and W13's rejected false commutator |
| U20 | Metric-compatible transport with a rank-matched invertible coframe and zero torsion induces Levi-Civita curvature by conjugation | Supplied smooth chart, internal form, connection and coframe; torsion equation holds | WRITTEN_PROOF of standard geometric adapter; 24 exact chart points checked by independent Christoffel/connection routes; no chart-emergence claim |
| W1 | One visible bit does not recover the other; cost depends on the target | Four presented pairs | EXACT_COUNTEREXAMPLE; named regression |
| W2 | Pairwise realizability need not give global realizability | Four even-parity triples | EXACT_COUNTEREXAMPLE; all pair intersections and empty triple intersection checked |
| W3 | Common scalar targets can be constant while joint readout is injective | Two coordinate cuts on four pairs | EXACT_COUNTEREXAMPLE; named regression and exhaustive common-target checks |
| W4 | One-step information can miss distinctions needed later | Four-state chain with declared cut | EXACT_WITNESS; depth 2, three repair labels, two bits |
| W5 | Joining separate single-arrow sufficient observers can miss a distinction revealed by a mixed word | Four-state model with declared A, B and c | EXACT_COUNTEREXAMPLE; individual memory sizes 1 and 2; family size 3; shortest separating word (A,B) |
| W6 | A finer initial observation can require more future-prediction memory | Three-state model; same T, different initial cuts | EXACT_COUNTEREXAMPLE; memory size rises from 1 to 2; does not contradict U3's fixed-target bound |
| W7 | Pairwise agreement on all commonly admitted words is not transitive | Four-state partial-arrow example | EXACT_COUNTEREXAMPLE; 0 agrees with 1 and 1 with 2, while word T separates 0 and 2 |
| W8 | Availability distinction need not be an executable output distinction; restricting a domain can increase memory | Constant two-state cut; total identity versus a partial restriction | EXACT_COUNTEREXAMPLE; output search returns none, availability query separates, memory rises 1 to 2 |
| W9 | Future permission differences can require memory even with constant actual readouts | Three-state partial chain; domain-aware target | EXACT_WITNESS; class counts 1,2,3, depth 2, memory 3; safe output distinction absent |
| W10 | A selected K family fixes a form but permits nontrivial metric-preserving loop transport | Declared EMK candidate maps, rank 3/4 examples | EXACT_WITNESS; invariant form reconstructed without supplying it to solver |
| W11 | The algebra permits distinct signature choices; the combined transport families have no nonzero common form | Same rank, different declared K/R families | EXACT_COUNTEREXAMPLE to signature selection by algebra alone |
| W12 | Missing couplings leave metric reconstruction ambiguous | One K pair in a four-component carrier | EXACT_WITNESS; invariant symmetric solution dimension four |
| W13 | Inadmissible compression can invent noncommutativity from a commuting native pair | Three-component cycle and square, lossy two-component projection | EXACT_COUNTEREXAMPLE; both exact transport descent gates reject |
| W14 | Metric compatibility alone does not identify fibre and metric curvature | Flat coframe and A=uK dv | EXACT_COUNTEREXAMPLE; nonzero torsion, fibre curvature nonzero, metric curvature zero |
| W15 | A supplied nonconstant coframe realizes a curved bridge; a constant one realizes a flat bridge | e=(du,f dv), A=f'K dv; f nonzero on the chart patch | EXACT_WITNESS and written identities; f remains an input |
| W16 | A genuine memory-losing cut can faithfully carry a nondegenerate geometric target | Declared memory-extended K family; source form radical equals the discarded memory direction | WRITTEN_COROLLARY and exact n=2 example; memory-sensitive targets remain lost |
| NG0 | A pregeometric recognition-transport law represents gravity before curvature | Need a selected native law, universal probe coupling and falsifiable gravitational predictions | PROPOSED; U17–U20 construct conditional representation bridges, not the physical identification |
| S0 | Associative continuation alone implies no nontrivial symmetry | Three-state chain | EXACT_COUNTEREXAMPLE; every permutation tested for commutation |
| Q0 | Quantum and classical physics are representations of one uncut law | Need native probability, interference, observables, dynamics and physical adapters | OPEN; no quantum axiom imported into U1–U6 or U8–U20 |
| ST0 | Spacetime is a derived optional representation | Need construction and physical interpretation of order, locality, dimension and metric | OPEN; v0.4 derives an invariant form only for a selected transport sector and uses a supplied smooth chart for curvature |
| A0 | Nature minimizes a cut cost along its continuation | Need an admissible path class, an independently defined action and a selection theorem or empirical test | OPEN; U2/U6/U9 minimize representation memory, not nature's action |
| CI0 | Cut ambiguity equals physical curvature or a universal information tensor | Need a typed adapter, invariant and scope-specific proof | OPEN; U7 alone is insufficient |
| G0 | Complete quantum gravity; RH or Yang–Mills closure | Not premises of this package | NOT CLAIMED |

## Certification coverage

The certificate enumerates all finite cuts/targets up to renaming labels for
1 <= |X| <= 4, all total endomaps on those carriers, and the specifically
declared represented block family. The extension separately enumerates all
ordered pairs of arrows and cuts through three states (3,678 cases), plus
17 declared four-state fixtures. It does not claim exhaustive two-arrow
coverage on four-state carriers. Hashes bind sources and evidence. This
coverage does not imply that a finite carrier models every physical possibility.

The v0.3 module adds all cuts and ordered pairs of partial endomaps through
three states (20,646 cases), 369,944 pair/mode comparisons, and named W7 on
four states. It does not exhaust all four-state partial families or certify
an availability-reporting instrument. Total regression count: 37.

The v0.4 module uses finite-dimensional linear carriers, whose state sets are
not finite. It adds the declared rational cases in GRAVITY_BEFORE_CURVATURE.md,
eight source pins and 13 regressions, bringing the package total to **50**.
Its 1,296 candidate-form comparisons exhaust only the stated integer grid
for 48 single-arrow matrices; they are not an exhaustive search of all forms
or geometries. U20's general proof and its 24 sampled chart checks are distinct.

The finite results are elementary self-contained mathematical consequences.
Their role is to make the proposed programme precise and refutable; no priority
claim is made for quotient factorization, finite refinement or memory counting.
