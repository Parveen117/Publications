# Claim ledger — v1.4 gravity/phase translator

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
| U21 | Connected K/R component couplings admit a nonzero invariant form iff every cycle has even K parity; when admitted, the form line and signature are fixed by sign propagation | Declared connected simple coupling graph, one sector per edge, nonzero invertible Cayley steps; dimension supplied | WRITTEN_PROOF in SECTOR_AND_COFRAME_SELECTION.md; all 646 connected labelled graphs through four components compared to full symmetric invariant-form equations |
| U22 | Exact compatibility plus minimum edge count does not select 1+3; even all-K trees can yield 2+2 | Four components, connected K/R graphs, normalized form; stated combinatorial objective | WRITTEN_PROOF and W17; application of established nonselection distinction, not a no-go for every possible native selector |
| U23 | Fixed homogeneous star connection gives zero torsion iff b'=kappa a, fixing b from a and one initial value; explicit LC curvature follows | Supplied chart, constant nonzero kappa, K or R star, coframe class e=(a du,b dx^i), ab nonzero | WRITTEN_PROOF of conditional adapter; 48 rational chart points, 240 independent curvature comparisons; both sectors pass |
| U24 | Zero-extended pair embeddings share a native K-even/R-odd involution iff the graph is a matching plus isolated vertices | Finite simple graph, usual real K/R embeddings, one common grading required | WRITTEN_PROOF in NATIVE_GRADING_AND_SEAM_MEMORY.md; full linear equations on 43 connected graphs through four vertices, plus disjoint-pair positive control |
| U25 | Minimum additional linear memory for endpoint readout and grading closure is rank(CJQ), equal to n-2 for connected bipartite graphs and n-1 otherwise | Declared equal-weight endpoint-summing cut on 2m edge coordinates, local direct-sum grading J, unrestricted linear repair | WRITTEN_PROOF; rank formula on 771 connected graphs through five vertices; 43 explicit minimal observers, grading intertwiners and cut-memory rank comparisons |
| U26 | Initial readout plus separately addressed nonzero edge R steps observes the whole edge carrier; minimum repair is 2m-n | Declared connected edge model; initial/readout transcripts refer to same input; all edge controls included | WRITTEN_PROOF; 43 addressed R families and 154 endpoint response inverses; actual measurement accessibility remains open |
| U27 | Real antisymmetric inter-edge mixing and local R steps obey one EMK cut grading and have a genuine nonzero commutator; the family preserves a represented quarter-turn | Declared real m-channel carrier, interaction graph and nonzero Cayley parameters; real dagger is transpose | WRITTEN_PROOF in NATIVE_INTERACTION_AND_LAWFUL_CUTS.md; 490 grading/orthogonality and quarter-turn checks, direct ordered-response witness |
| U28 | Common invariant real subspaces are sums of interaction components; minimum closed observer rank is 2 times the number of labels in components seen by C; invariant symmetric forms are scalar identities per component | All individually addressed local R and mixing controls included, nonzero parameters, real-linear finite carrier | WRITTEN_PROOF; 360 graph-versus-matrix closure comparisons, 700 projection intertwiners and 11 full invariant-form solves; proper lossy interacting quotient |
| U29 | Two scalar readings per channel along a rooted tree reconstruct the connected carrier with a minimal fixed linear transcript | One first-coordinate readout, calibrated nonzero local/mixing controls, stated same-input transcript contract | WRITTEN_PROOF; 167 full-rank transcripts, 1,298 basis-vector decoder checks and chronological-word probes; physical preparation/readout not established |
| U30 | Calibrated native responses define a Hermitian PSD tensor with exact grading conjugation, transport invariance and frame covariance | Real paired carrier, positive metric, represented quarter-turn; no physical probability law | WRITTEN_PROOF in NATIVE_RESPONSE_TENSOR_AND_CUT_LEDGER.md; direct complex Gram oracle and 490 transport checks |
| U31 | Full tensor splits under every response catalogue iff the orthogonal cut commutes with Z; minimum paired repair is rank([C;CZ])-rank(C) | Real output-coordinate cut; quarter-turn Z squared equals -I; ambient metric transported with frame | WRITTEN_PROOF; 340 coordinate cuts, 14 general repairs, cross-seam and indefinite-gap controls |
| U32 | Signed native responses obey inherited total covariance and tower ledgers; exact response and squared-norm recovery have the same input-kernel gate | Supplied finite ensemble or stated whole-carrier linear target; separate cut contracts | WRITTEN_PROOF; 405 ledgers, 1,620 tower checks, 90 exact scale/quotient checks; no gravitational identification |
| U33 | Complete order response determines phase t and mixing sine but has precisely s versus 1/s blindness; one direct reading repairs it | Oriented two-channel native pair, finite nonzero t,s, fixed calibration | WRITTEN_PROOF in PROBE_CALIBRATION_AND_NATIVE_IDENTIFICATION.md; 100 inversions, 4,950 pair comparisons and 40 nontrivial reciprocal collisions |
| U34 | Equal augmented calibrated triples identify the same U,V pair and generated words across probes | Declared four-dimensional pair family; independent two-sided frames | WRITTEN_PROOF; 300 frame checks, 400 held-out basis predictions, post-hoc fit and hidden-channel controls |
| U35 | Separated-denominator interval inversion encloses every compatible model; disjoint data/parameter boxes reject a shared pair | Independently justified finite error bounds; conditional model class; overlap is not acceptance | WRITTEN_PROOF; 2,700 enclosure and ratio-bound checks, small-signal and large-mixing margin controls |
| U36 | Declared quadratic source cost has a unique stationary response, native affine Cayley continuation, cut-flux ledger and exact stationary elimination | Positive grounded graph, linear source, independently aligned native coordinates; conservative flow does not imply settling | WRITTEN_PROOF in SOURCE_RESPONSE_AND_GRAVITY_TESTS.md; 86 weighted graphs, 332 source solutions, 2,204 cut flux and 2,204 stationary elimination checks |
| U37 | Source-probe cross energy gives inverse-square finite-difference response iff calibrated shell capacity is quadratic in independently assigned radius | Shell-constant sector, quadratic cost, positive displacement calibration; source/probe physical adapters additional | CONDITIONAL_WRITTEN_PROOF; 96 radial profiles and 176 pair-energy force checks; other exponents remain allowed |
| U38 | A nonlinear mismatch cost can duplicate inverse-square distance dependence while changing source-strength scaling; fitted radius can fake the exponent | Declared p>1 cost family, weak-probe first variation, fixed independent calibration | WRITTEN_PROOF and negative controls; 108 nonlinear scaling fixtures, 54 bounded-error ratio checks; physical source/distance comparison not performed |
| U39 | Existing native graph/generator rules cannot select a unique shell-capacity exponent | Current U27 arbitrary supplied connected graph, unit positive link weights, common degree four and 49 carriers | WRITTEN_COUNTERMODEL in SHELL_CAPACITY_SELECTION.md; capacity 4,12,20 versus 4,6,6; equal rank/weight controls |
| U40 | Supplied d-direction nearest-neighbour graph has exact combinatorial shell capacity with leading degree d-1; d=3 yields 12n²+12n+6 | Explicit three-direction integer graph, graph distance, unit weights; physical ruler separately calibrated | CONDITIONAL_WRITTEN_PROOF; 16 lattice profiles and 40 independent enumerated edge cuts; no primitive dimension selection |
| U41 | Aggregate inverse capacity need not equal local or shell-mean response; shell equitability is sufficient for radial stationarity | Declared rooted grounded network and U36's positive source cost | WRITTEN_PROOF and exact rational counterexample: cubic second mean drop 13/378 versus 1/30; tree control radial |
| TF-1 | Stationary elimination of hidden graph coordinates constructs a positive quadratic equilibrium fundamental relation on two visible coordinates, with a Schur complement Hessian | U36 supplied grounded graph, positive quadratic cost and two selected visible coordinates | CONDITIONAL_WRITTEN_PROOF in THERMODYNAMIC_FOUNDATION_BRIDGE.md; source-bound exact rational graph witness and hidden-basis control; graph and cost law remain inputs |
| TF-2 | If the two visible coordinates are physically calibrated as entropy and volume and the projected cost as energy, Maxwell/contact and response identities follow | TF-1 plus positive chart and independent entropy/volume/energy adapter | CONDITIONAL_WRITTEN_PROOF; source-bound exact rational witness gives T=1, P=8, Cp/Cv=Ks/Kt=35/26 and a planted Maxwell failure; physical adapter not derived |
| QB-1 | U36's paired real Cayley flow around the common stationary source is a complex unitary Cayley evolution with the **same** grounded graph stiffness | U27 represented quarter-turn and norm, U36 supplied grounded graph, affine source, native step parameter | CONDITIONAL_WRITTEN_PROOF in QUANTUM_CLASSICAL_BRIDGE.md; source-bound exact rational defining-equation and unitarity checks; no physical clock or hbar claimed |
| QB-2 | Native sitewise quadratic norm shares give positive normalized squared-amplitude weights, unique within the declared local quadratic/rotation/permutation-invariant readout class | Supplied orthogonal detector splitting, positive norm, local homogeneous quadratic and detector-covariance assumptions | CONDITIONAL_WRITTEN_PROOF and exact normalization witness; physical event-frequency rule and detector choice not derived |
| QB-3 | A declared phase-erasing cut converts one unitary step to a doubly stochastic classical update; the off-diagonal seam-memory residue exactly measures the failure of uncut diagonal closure | QB-1 and QB-2 with declared site-basis dephasing; repeated dephasing for classical protocol | CONDITIONAL_WRITTEN_PROOF; same graph yields coherent 16/25, phase-erased 8/25 and exact residue 8/25; same-diagonal different-future control |
| GP-1 | Complete **phase-referenced** coherent Cayley matrix of a supplied grounded graph reconstructs its native Green response `H^{-1}`, and hence every U36 source and cross cost | QB-1/U36 with positive graph stiffness, known nonzero native step, fixed source units and an absolute phase reference | CONDITIONAL_WRITTEN_PROOF in GRAVITY_PHASE_TRANSLATOR.md; exact rational coherent inverse and step-rescaling control; ordinary quantum state tomography does not supply that global phase |
| GP-2 | Even the complete one-step quantum density channel need not identify the same native Green/source response | Supplied two-site connected grounded graphs, known native step and same site detectors; no external absolute-phase reference | WRITTEN_COUNTERMODEL: Cayley matrices differ by global phase `(4-3i)/5`, so all quantum state channels and classical transitions agree, yet cross costs are `-1/3` versus `-5/8`; exact rational controls |
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
| W17 | K-star, K-path and R-star tie at the coupling minimum but give three different signatures | Four response components and three edges per graph | EXACT_COUNTEREXAMPLE; signatures (1,3), (2,2), (4,0) |
| W18 | Odd K cycle destroys the invariant form; a balanced mixed triangle admits one | Three components; all-K versus two-K/one-R triangle | EXACT_COUNTEREXAMPLE and exact positive control |
| W19 | A pointwise torsion zero need not give a patch solution for a fixed connection | a=kappa=1, b=1+u^2 | EXACT_COUNTEREXAMPLE; torsion coefficient 2u-1 vanishes at one point only |
| W20 | Both K and R permit homogeneous coframe closure | Same n=3, kappa=1 and b=3+u; different declared sectors | EXACT_COUNTEREXAMPLE to selection by the torsion gate alone |
| W21 | Metric-sign involution and the local native EMK cut grade K differently | Same two-component generator matrices | EXACT_COUNTEREXAMPLE to identifying the involutions while keeping the generators fixed |
| W22 | A forgotten difference between two endpoint copies returns to visible readout after grading | Three-vertex star, state (1,0,-1,0) | EXACT_COUNTEREXAMPLE to bare gluing as an exact grading quotient |
| W23 | Three-vertex tree needs one extra grading channel; four-vertex trees need two | Declared endpoint readout and grading target | EXACT_WITNESS; repaired observer and intertwiner constructed |
| W24 | Grading-sufficient memory can still miss an addressed edge response | Triangle, grading-blind state (1,-1,-1,1,1,-1) | EXACT_COUNTEREXAMPLE; grading observer rank five, addressed-response rank six |
| W25 | Declared K/R candidate steps give different signed transcripts on the same input | Triangle witness, addressed edge (0,1), known parameter 1/2 and calibration | EXACT_MODEL_PREDICTION; maximum coordinate separation 16/15; no physical realization claimed |
| W26 | Compressing commuting native edge steps can invent a commutator | Direct-sum R steps, endpoint-summing cut and averaging section | EXACT_COUNTEREXAMPLE; specialized U19/W13 control, both descent gates reject |
| W27 | Actual step order changes the readout before compression | Two linked carriers, parameters 1/2, input (0,0,0,1) | EXACT_MODEL_PREDICTION; outputs 0 and 16/25, nonidentity group commutator |
| W28 | A genuine lossy quotient preserves grading, a declared form and nonzero interaction | Three carriers, mixing link (0,1), cut discards isolated carrier 2 | EXACT_WITNESS; source dimension six, target four, all step intertwiners pass |
| W29 | A partial cut inside an interacting component fails exact descent | Two linked carriers, cut keeps only one | EXACT_COUNTEREXAMPLE; full closure needs two additional scalars |
| W30 | Aggregate cancellation is not a blind interaction component | C=(1,0,-1,0), two unlinked individually controlled carriers | EXACT_COUNTEREXAMPLE; zero initial output becomes -2/5 after a local R step |
| W31 | Connected energy-preserving interacting family has only a definite invariant form line | Full stated family on two carriers | EXACT_CONTROL and U28 form proof; no Lorentzian metric derived |
| W32 | Dropping mixing or local R controls changes response observability | One scalar readout on two carriers | EXACT_COUNTEREXAMPLE; rank two for either reduced catalogue, rank four for the full connected family |
| W33 | A real half-channel cut loses an imaginary seam and can give an indefinite full-minus-visible tensor | One paired channel, identity response catalogue | EXACT_COUNTEREXAMPLE; determinant -1, exactly one extra scalar repairs it |
| W34 | Paired tensor closure need not close the dynamics | Two interacting channels, one retained | EXACT_COUNTEREXAMPLE; Z descends, mixing does not |
| W35 | Nonzero two-channel order response has a positive scalar Gram form of rank four | Nonzero finite Cayley parameters; explicit formula for kappa | WRITTEN_IDENTITY and exact controls; at parameters 1/2, kappa=64/125 |
| W36 | A proper six-to-four quotient preserves the signed response | A linked pair and one isolated channel | EXACT_WITNESS; quotient response intertwiner |
| W37 | Zero average does not mean zero response; the full Gram tensor loses an overall response sign | Opposite signed preparations/responses | EXACT_COUNTEREXAMPLES; nonzero covariance and equal tensors of L and -L |
| W38 | Noncommuting SLDs may have zero antisymmetric state average | Declared qubit comparison at rho=I/2, SLDs sigma_x/sigma_y | EXACT_COUNTEREXAMPLE in corrected QTH-1 certificate and regression; no quantum axiom used to prove U30–U32 |
| W39 | Different native mixing laws can have identical full signed response and tensor | t=1/2, s=1/2 versus s=2 | EXACT_COUNTEREXAMPLE; direct c readings are 3/5 and -3/5 |
| W40 | The augmented comparison rejects reciprocal probe disagreement | Fixed signed calibration and the W39 family | EXACT_CONTROL; SHARED_PAIR_REJECTED |
| W41 | Post-hoc output decoding can force order-response agreement | Invertible native response Lp and fitted D=L0 Lp inverse | EXACT_COUNTEREXAMPLE; not the admissible independently fixed two-sided frame |
| W42 | A hidden third native channel can differ while all five raw readings agree | Selected pair plus separately controlled unobserved channel | EXACT_COUNTEREXAMPLE; extra R readings 3/5 versus 12/13 |
| W43 | A direct readout separates reciprocal laws with a finite error budget | W39 gap 6/5, explicit interval radii | EXACT_CONTROL; radius 1/10 rejects, radius 3/5 leaves the implemented comparison unresolved |
| W44 | Exact inverse has no uniform denominator margin over all nonzero parameters | Weak phase and very large mixing examples | EXACT_CONTROL; INSUFFICIENT_RESOLUTION rather than false acceptance |
| W45 | Quadratic shell stiffness at uniform spacing gives exact inverse-square response | C_n=n^2; unit source/probe and independently fixed spacing | CONDITIONAL_EXACT_WITNESS; responses -1, -1/4, -1/9, -1/16 |
| W46 | Same native generator rules admit different response exponents | C_n=n^alpha with alpha=0,1,2,3 | EXACT_NONSELECTION_CONTROL; no primitive choice of alpha=2 |
| W47 | Inverse-square distance dependence can coexist with nonlinear source response | Cubic cost, quartic stiffness growth, weak probe | EXACT_CONTROL; source times four gives response times two |
| W48 | Onsite pinning breaks the source-free conserved-flux law | Positive added diagonal stiffness | EXACT_CONTROL; outward currents 27/61, 20/61, 18/61 |
| W49 | Erasing hidden coordinates without the effective cost changes source response | Chain stiffnesses 1,4,9 | EXACT_CONTROL; source effective stiffness 36/49 rather than 1 |
| W50 | Fitting the ruler to the response can manufacture a distance exponent | A 1/n response relabelled by r=sqrt(n) | EXACT_CONTROL; independent calibration required |
| W51 | Conservative source continuation does not imply relaxation | Nonstationary affine Cayley orbit | WRITTEN_NONCONVERGENCE_ARGUMENT; exact excess-cost preservation checks |
| W52 | Same 49 carriers, degree four and equal coupling weights yield different rooted capacity profiles | Periodic square versus two-jump cyclic graph | EXACT_COUNTERMODEL; 4,12,20 versus 4,6,6 |
| W53 | A three-direction unit graph has quadratic leading shell count, with linear and constant finite corrections | L1 nearest-neighbour graph | EXACT_IDENTITY; 12n²+12n+6, not pure 12n² at every graph radius |
| W54 | Aggregate capacity counting alone does not determine the radial field | Grounded cubic ball, source at origin | EXACT_COUNTEREXAMPLE; unequal shell-two potentials and 13/378 versus 1/30 |
| W55 | A quadratic upper capacity bound does not force quadratic saturation | One-dimensional rooted path, unit edges | EXACT_CONTROL; constant capacity 2 |
| NP0 | Selected probes share one calibrated native U,V pair | Predeclared family, coordinates, preparations and readout; includes augmented direct response | MODEL AGREEMENT GATE; U34 proves its sufficiency in that family; physical universality and gravitational identity remain OPEN |
| NG0 | A pregeometric recognition-transport law represents gravity before curvature | Need a selected native law, universal probe coupling and falsifiable gravitational predictions | PROPOSED; U36-U38 add an explicit source candidate and conditional inverse-square sector with a test contract; physical selection and validation remain open |
| S0 | Associative continuation alone implies no nontrivial symmetry | Three-state chain | EXACT_COUNTEREXAMPLE; every permutation tested for commutation |
| Q0 | Quantum and classical physics are representations of one uncut law | Need physically selected native law, detector frequencies, interference, observables, calibrated dynamics and physical adapters | PARTLY CONSTRUCTED MATHEMATICALLY by QB-1–QB-3 on one declared graph; physical Born frequency law and universal quantum/classical identification OPEN; no quantum axiom imported into U1–U6 or U8–U35 |
| ST0 | Spacetime is a derived optional representation | Need construction and physical interpretation of order, locality, dimension and metric | OPEN; coframe constraint remains conditional; v0.6 repairs native grading on a declared larger carrier but does not construct a spacetime |
| A0 | Nature minimizes a cut cost along its continuation | Need an admissible path class, an independently defined action and a selection theorem or empirical test | OPEN; U36 supplies a declared stationary cost; W51 proves its conservative continuation does not force relaxation |
| CI0 | Cut ambiguity equals physical curvature or a universal information tensor | Need a typed adapter, invariant and scope-specific proof | OPEN; U7 alone is insufficient |
| G0 | Complete quantum gravity; RH or Yang–Mills closure | Not premises of this package | NOT CLAIMED |
| SR-1 | A unitary sector return has nonnegative target-relative squared defect, while holonomy operator norm is blind to every unitary phase | Declared cut carrier, unitary loop continuation and target map | CONDITIONAL_WRITTEN_IDENTITY; exact +1/-1 checks; not an entropy or an uncut observable |
| SR-2 | Local flatness can coexist with global -1 return, and principal endpoint +1 can hide integer phase winding | 3x3 torus sign transport; separate abelian path lift | EXACT_COUNTEREXAMPLES; global winding control consumes the prior RKF lift theorem; no physical metric or thermodynamic equivalence |
| SH-1 | A scalar -1 has distinct flat spin and conical vector realizations; on a cone with deficit pi the vector return is -I and spin return is ±i | Declared spin double cover and spin structure; supplied 2D cone metric and deficit range | CONDITIONAL_NONSELECTION; exact quarter-turn and flat monodromy controls; no identification of native sign with spin connection |
| SH-2 | Shrinking contractible-loop holonomy yields connection curvature, which becomes metric curvature under U20's coframe and torsion gates | Supplied smooth connection, oriented loop areas, nondegenerate coframe and compatible torsion-free connection | CONDITIONAL_STANDARD_GEOMETRIC_BRIDGE; exact f=1+u² rectangle coefficient; no physical gravitational source law |
| SH-3 | Dimensionless sign/angle alone cannot select a numerical dimensionful hbar | Sign-only or angle-only return data, no physical action calibration | DIMENSIONAL_NONSELECTION; richer primitive not excluded |
| AI-1 | A declared acceleration/phase adapter maps the native stationary response to a source-configuration double-difference phase | Fixed grounded graph; calibrated source, signed gradient/trajectory rows, physical clock, laser phase response and error budget | CONDITIONAL_WRITTEN_INTERFACE and exact kernel/linearity controls; physical mapping not selected by the primitive |
| AI-R | Nominal finite-cylinder Newtonian reference gives approximately 0.5304493433 rad; published total is 0.547870 rad | External reference G, nominal geometry, fixed point-cloud trajectories, instantaneous pulses; supports/recoil/cloud average omitted | RETROSPECTIVE_NUMERICAL_CONTROL; approximately 3.18% shortfall; NOT a native prediction, precision reproduction or held-out test |

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

The v0.5 module adds 13 regressions, bringing the total to **63**. It enumerates
759 simple labelled absent/K/R graphs on two through four vertices, compares
all 646 connected cases against full invariant-form equations, and excludes
113 disconnected cases from U21's scope. It checks 48 supplied homogeneous
chart points and 240 curvature matrices by independent metric/connection routes.
These domains are finite; they do not exhaust arbitrary transport families or
coframes. The source-bound certificate now binds 27 files. A source hash is
provenance, not independent verification of every claim in the source.

The v0.6 module adds 14 regressions for a total of **77** and binds 32 source
files. Its 1,098 simple graph enumerations include 771 connected cases for the
grading-memory rank formula. The 43 connected graphs through four vertices
add full common-grader linear equations, constructed minimal observers,
cut-memory ranks and complete addressed-R response ranks, including 154 exact
endpoint inverses. No exhaustive enumeration of physical instruments or all
native many-body assemblies is asserted. General repair principles are inherited;
these calculations specialize them to the stated edge/readout model.

The v0.7 module adds 16 regressions for a total of **93** and binds **37** source
files. All 75 simple interaction graphs on one through four channel labels
are checked at phase/mixing parameters 1/2. Evidence includes 490 orthogonality/
grading checks, 490 quarter-turn commutations, 360 readout-closure comparisons,
700 exact component projections, and 167 minimal transcripts with 1,298
basis-vector reconstructions and signed-probe chronological checks. Full
invariant-form nullspaces are separately solved on all 11 graphs through three
channels. The campaign does not enumerate all parameters, instruments or
possible native interaction laws.

The finite results are elementary self-contained mathematical consequences.
Their role is to make the proposed programme precise and refutable; no priority
claim is made for quotient factorization, finite refinement or memory counting.

## v0.8 response and source correction coverage

The native package adds 18 regressions for **111** total. Exact coverage adds
490 tensor transport identities, all 340 real coordinate cuts through eight
dimensions (30 preserve complete coordinate pairs), 14 noncoordinate repairs,
90 signed-response scale/quotient comparisons, 405 finite-ensemble covariance
ledgers and 1,620 nested ledgers. Pair-distance and iterative-closure oracles
are independent of the conditional-mean and one-step repair constructions.
There are **49 bound sources**, including six corrected CID/QTH surface files.

CID-1's preexisting certificate matches its pinned digest read-only. Its
20 tests and the corrected QTH-1's 30 tests pass separately from the native
suite. The QTH correction replaces an invalid general converse with a checked
maximally mixed counterexample and clarifies the difference identity.
The ten external source pins record both current research and unmerged
correction commits. Written proofs, finite PASS, source hashes, formal proofs,
external review and empirical validation remain separate.

## v0.9 calibration and identification coverage

The package adds 19 regressions for **130** total and binds **54 source files**.
The new finite campaign uses ten signed nonzero rational parameters in each
Cayley slot, giving 100 pairs, 4,950 pairwise order-matrix comparisons and
exactly 40 distinct reciprocal collisions. The augmented decoder separates
all distinct pairs. Its recovered operators predict 400 vector responses on
the held-out e2/e4 preparations, and 300 known frames preserve the protocol.
For each parameter pair, all 27 endpoint/centre perturbations at reported error
radius 1/10000 are checked: 2,700 enclosures and 2,700 phase ratio bounds.

The source audit pins six complete files (ID2 is used for its identification
section). There is no exhaustive coverage of instruments, unknown native
graphs, quantum measurements or physical probe species. NP0 is an internal
agreement gate, not a proof of gravitational universality. Non-overlap can
reject; overlap at finite error is unresolved. Exact certificates and written
proofs are separate from formal verification, external review and experiment.

## v1.0 conditional source-response coverage

The package adds 22 regressions for **152** total and binds **60 source files**.
The campaign considers every connected labelled graph on 2, 3 or 4 vertices,
each with unit and a fixed rational nonuniform weight pattern: 86 networks.
It checks 332 source solutions against edge-current sums, 2,204 cut-flux
ledgers, 2,204 stationary eliminations including hidden-source constants,
332 cost-completion identities and 86 source-coupled native affine laws.
Independent radial formulas match 96 matrix solutions; 176 source-probe
cross-energy differences match the claimed signed responses. There are 108
nonlinear distance/source controls and 54 bounded-error ratio fixtures.

Seven immutable repository sources are pinned. Two published experiments are
benchmark references only: their raw data and apparatus models have not been
ingested. The inverse-square sector assumes quadratic cost and calibrated
quadratic shell capacity. It is not a primitive-only gravity derivation or
an empirical PASS. The proposed physical protocol has no acquired data yet.
