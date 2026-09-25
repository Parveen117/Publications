# Claim ledger — v0.1

The word **proof** below refers to the written argument under the manuscript's
declared assumptions. **PASS_FINITE_CHECKS** refers only to the finite domains
recorded in `CERTIFICATE.json`. Neither status denotes proof-assistant checking,
independent review or physical validation.

| ID | Statement | Assumptions | Status and evidence |
| --- | --- | --- | --- |
| P0 | Uncut ground precedes formal distinctions; the framework is itself a cut | Proposed interpretation and admissibility programme | PROPOSED; manuscript §§1–2 |
| U1 | Exact target recovery iff target is fibre-constant; readout-only correction cannot separate a merged pair | Nonempty finite presented X, deterministic c and g | WRITTEN_PROOF; exhaustive decoder comparison on all equality patterns through four states |
| U2 | Minimum repair alphabet is K_g(c); minimum fixed binary length is ceil(log2 K) | Arbitrary deterministic repair channels allowed; cost is worst-case label count | WRITTEN_PROOF; constructive decoder and exhaustive exclusion of smaller alphabets through four states |
| U3 | Refinement reduces ambiguity and K_g(c) <= K_d(c) K_g(d) | c factors through d on the same X | WRITTEN_PROOF; exhaustive partition triples through four states; strict inequality control |
| U4 | A common target factors through every cut iff constant on generated components | Finite family on a declared joint presentation | WRITTEN_PROOF; exhaustive pairs of cuts and target equality patterns through four states |
| U5 | Observed deterministic continuation exists iff next readouts are fibre-constant; lawful refinements intertwine | Declared total T on X; comparison maps well-defined | WRITTEN_PROOF; every endomap and cut through four states, plus every applicable refinement |
| U6 | History refinement terminates by depth N-q and is the coarsest dynamically closed refinement | N presented states, q initial readouts; one deterministic T; unrestricted access to model signatures for the construction | WRITTEN_PROOF; exhaustive comparison against all closed refining partitions through four states |
| U7 | R_delta_gamma - R_delta R_gamma = J_delta C_gamma | Declared composable linear transports; complementary cut projections | INHERITED_WRITTEN_PROOF from RSC; all scalar 2x2 block pairs with entries -1,0,1 checked |
| W1 | One visible bit does not recover the other; cost depends on the target | Four presented pairs | EXACT_COUNTEREXAMPLE; named regression |
| W2 | Pairwise realizability need not give global realizability | Four even-parity triples | EXACT_COUNTEREXAMPLE; all pair intersections and empty triple intersection checked |
| W3 | Common scalar targets can be constant while joint readout is injective | Two coordinate cuts on four pairs | EXACT_COUNTEREXAMPLE; named regression and exhaustive common-target checks |
| W4 | One-step information can miss distinctions needed later | Four-state chain with declared cut | EXACT_WITNESS; depth 2, three repair labels, two bits |
| S0 | Associative continuation alone implies no nontrivial symmetry | Three-state chain | EXACT_COUNTEREXAMPLE; every permutation tested for commutation |
| Q0 | Quantum and classical physics are representations of one uncut law | Need native probability, interference, observables, dynamics and physical adapters | OPEN; no quantum axiom imported into U1–U6 |
| ST0 | Spacetime is a derived optional representation | Need construction, hypotheses and physical interpretation of order, locality, dimension and metric | OPEN IN THIS PACKAGE; no new audit or promotion of other spacetime claims |
| A0 | Nature minimizes a cut cost along its continuation | Need an admissible path class, an independently defined action and a selection theorem or empirical test | OPEN; U2/U6 minimize representation memory, not nature's action |
| CI0 | Cut ambiguity equals physical curvature or a universal information tensor | Need a typed adapter, invariant and scope-specific proof | OPEN; U7 alone is insufficient |
| G0 | Complete quantum gravity; RH or Yang–Mills closure | Not premises of this package | NOT CLAIMED |

## Certification coverage

The certificate enumerates all finite cuts/targets up to renaming labels for
1 <= |X| <= 4, all total endomaps on those carriers, and the specifically
declared represented block family. Hashes bind sources and evidence. This
coverage does not imply that a finite carrier models every physical possibility.

The finite results are elementary self-contained mathematical consequences.
Their role is to make the proposed programme precise and refutable; no priority
claim is made for quotient factorization, finite refinement or memory counting.
