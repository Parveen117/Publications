# Physics from Cuts: Uncut Ground, Measurement and Lawful Continuation

**Monty Dabas — working research programme, version 0.8, 25 September 2026.**

The starting proposal is reflexive: definitions create distinctions, so the
mathematical framework used to describe the uncut is itself a cut. A measurement
must therefore retain its declared observational context and unresolved
information. Spacetime is a possible derived representation of the programme;
it is not an input to the finite results here.

Read the [novelty assessment](NOVELTY_AND_LINEAGE.md),
[foundational manuscript](MANUSCRIPT.md),
[family-continuation extension](FAMILY_CONTINUATION.md),
[admissibility extension](ADMISSIBLE_CONTINUATION.md),
[gravity-before-curvature bridge](GRAVITY_BEFORE_CURVATURE.md),
[sector and coframe constraints](SECTOR_AND_COFRAME_SELECTION.md),
[native grading and seam memory](NATIVE_GRADING_AND_SEAM_MEMORY.md),
[native interaction and lawful cuts](NATIVE_INTERACTION_AND_LAWFUL_CUTS.md),
[native response tensor and cut ledger](NATIVE_RESPONSE_TENSOR_AND_CUT_LEDGER.md),
[information/thermodynamics source audit](RESPONSE_SOURCE_AUDIT.md),
and [claim ledger](CLAIMS.md).
The [research programme](RESEARCH_PROGRAMME.md) specifies the next derivations.

**Lineage:** U1 already appears exactly in Spectral I Theorem 18.1; minimum
repair has a linear counterpart in Spectral II Theorem 7.2; recovered middle
identity is already required by RSC 27/30. Partition refinement is established
mathematics. This package adds explicit finite formulations, operational
distinctions and executable controls; proposition numbering is not a priority claim.

## Foundational results U1–U7

Within an explicitly declared finite presentation, the manuscript proves:

1. A target can be recovered from a cut exactly when it is constant on every
   readout fibre. A deterministic correction cannot create missing distinctions.
2. The minimum extra memory alphabet is the maximum number of distinct target
   values hidden behind one readout. This is a target-relative cost.
3. Refinement reduces target ambiguity; sequential repair obeys a product bound,
   which need not be an equality.
4. A clock-free transition descends to observed dynamics exactly when equivalent
   observations have equivalent next observations.
5. Repeated distinction by future readouts terminates on a finite presentation
   and constructs the coarsest dynamically closed refinement of a cut.
6. Common scalar invariants may be trivial even when several cuts jointly
   distinguish all presented states. Pairwise consistency need not give a global
   state compatible with all observations.

The inherited RSC cut-corner composition identity is restated and checked
separately as a downstream linear representation. These are elementary finite
derivations and a connection to the existing native framework, not a claim of
priority for the underlying mathematical facts.

## Added in v0.2: memory for all admitted total continuations

[U8–U12](FAMILY_CONTINUATION.md) extend the construction to a declared finite
family of arrows and every word they generate:

- The coarsest observer closed under the whole family exists constructively.
- Its minimum additional memory is exactly the largest number of distinct
  continuation classes hidden behind a single initial readout.
- Naming an already available composite as a new arrow changes neither that
  partition nor its memory cost. Refinement depth can still change.
- Closing two initial observations jointly agrees with joining their closures,
  provided the admitted arrow family is the same.
- A finite pair search finds a shortest experiment distinguishing two candidates,
  or proves that no word in the declared family distinguishes them.

**Separating example:** an observer sufficient for repeated A and an observer
sufficient for repeated B can still miss a distinction revealed by A followed
by B. In the four-state W5 example the family needs **three memory labels**,
although the separate requirements are one and two. The shortest diagnostic
word is (A,B). W6 shows why extra present detail can also increase the future
prediction task rather than universally reducing its memory cost.

The v0.2 evidence comprises **3,678 exhaustive cut/two-arrow-family
cases through three states**, 17 additional four-state fixtures, and 33,206
candidate-pair diagnostic checks. The original v0.1 exhaustive domains are
retained separately in the certificate. These counts are finite coverage,
not a replacement for the written proofs.

## Added in v0.3: partial arrows and valid experiments

[U13–U16](ADMISSIBLE_CONTINUATION.md) separate an observer that must retain
enabled-domain information from one restricted to actual experiments admitted
at both candidates. Exact partial descent requires both whole-fibre domain
agreement and next-readout agreement. Closure again gives a coarsest sufficient
observer and an exact minimum repair alphabet under unrestricted repair maps.

**W7** shows that agreement on all commonly admitted experiments need not be
transitive, so it cannot always define observer fibres. **W8** separates an
availability query from executing a forbidden action, and shows that restricting
a domain can increase domain-aware memory. **W9** gives delayed availability
distinctions even when every actual readout is zero.

The v0.3 evidence includes **37 regression tests**, plus **20,646 exhaustive
cut/two-partial-arrow cases through three states**, **369,944 pair/mode diagnostic
checks**, and the retained v0.1/v0.2 evidence. A model's availability flag is not
automatically an observable; operational use needs a justified interface.

## Added in v0.4: transport before metric and curvature

[U17–U20](GRAVITY_BEFORE_CURVATURE.md) reconstruct compatible bilinear forms
from declared invertible transports. A chosen EMK K-sector family determines
diag(1,-1,...,-1) up to scale; four supplied components give a Lorentz-type
1+3 signature. The R-sector instead determines a definite form, so the native
algebra alone does not select a physical signature or a dimension.

An exact cut must preserve the geometric target and intertwine the transports.
A memory-extended example passes both gates: the cut discards a target-blind
memory coordinate and recovers a nondegenerate quotient form.
A counterexample shows that inadmissible compression can invent noncommutativity
from a commuting native pair. With a separately supplied smooth coframe and a
metric-compatible, torsion-free connection, the connection curvature becomes
metric curvature by an explicit intertwiner. A valid curved chart is checked
by two independent routes; a torsion control rejects a false identification.

The v0.4 evidence: **50 regression tests**, 24 K-family reconstructions,
1,296 integer-grid candidate-form comparisons, 24 supplied chart points and
the retained earlier certificates. The new source audit pins eight files.
The K maps are candidate algebraic transports, not asserted to be RKF 55's
dagger-unitary flows. Selecting the gravitational sector, coframe, dynamics
and shared quantum readout remains the next research problem.

## Added in v0.5: coupling constraints and conditional coframe closure

[U21–U23](SECTOR_AND_COFRAME_SELECTION.md) extend the chosen star to connected
component-coupling graphs. A nonzero invariant form exists precisely when every
cycle has an even number of K edges; the coupling pattern then fixes its
signature up to an overall sign. Exact compatibility plus the minimum number
of couplings still admits 1+3, 2+2 and definite signatures on four components.
This explicitly rejects that proposed sector-selection rule.

There is a positive coframe constraint: for a supplied constant star connection
and homogeneous coframe e=(a du,b dx^i), zero torsion is equivalent to
**b'=kappa a**. Given a and an initial value, b is determined on a nondegenerate
patch. Both K and R sectors pass, so this conditional adapter does not select
the physical sector, chart or homogeneity. No gravitational field equation is
imported. The existing information-invariance nonselection result is credited.

The v0.5 evidence: **63 regression tests**, all **646 connected K/R coupling
graphs on two through four labelled components** compared with a full
invariant-form solver, and **240 independent curvature-matrix comparisons**
at 48 supplied chart points. All earlier evidence is retained. The certificate
binds 27 source files; that provenance manifest pins four selection sources.

## Added in v0.6: native grading and the exact memory cost of gluing

[U24–U26](NATIVE_GRADING_AND_SEAM_MEMORY.md) test whether the local native EMK
grading survives the multi-pair assembly. Its cut J_E=K is distinct from a
metric's signature partition. Overlapping zero-extended K/R pairs on a connected
carrier with at least three components admit no common involution keeping each
K even and each R odd. The earlier candidate invariant-form results remain
valid; their native-grading identification was not established.

A direct sum retaining both endpoint and edge identities does carry the local
native grading. Summing endpoint copies into one vertex readout loses it.
The exact minimum repair is **n-2 extra linear scalar channels for a connected
bipartite graph, or n-1 otherwise**, when the target is grading closure.
If the target includes separately addressed R-sector steps on all m edges,
the minimum becomes **2m-n**. On a triangle those costs are two and three.
An explicit hidden state produces distinct K/R candidate transcripts, making
the model-level comparison concrete; physical instrument realization is open.

The v0.6 evidence: **77 regressions**, rank checks on **771 connected graphs**
through five vertices, **43** complete small-graph grading/observer/response
checks and **154 exact endpoint reconstructions**. Earlier campaigns remain
bound in its certificate, covering **32 source files**. The direct-sum
repair has commuting separate-edge transports; native interaction between
edges must be supplied before recovering a noncommutative gravity candidate.

## Added in v0.7: genuine graded interaction and a lawful lossy quotient

[U27–U29](NATIVE_INTERACTION_AND_LAWFUL_CUTS.md) construct a declared
antisymmetric mixing family between retained edge carriers. Together with
local R controls it preserves one native cut grading, a positive form and a
represented quarter-turn, while producing actual nonzero transport commutators.
An explicit order test gives scalar outputs **0 and 16/25** for the same input.

Exact real-linear cuts can discard whole interaction components. A proper
**6-to-4 dimensional quotient** loses one unobserved component while preserving
grading, the declared form target and nonzero interaction. Within a connected
component, every nonzero exact closed readout needs its full relevant memory.
A fixed transcript of **2m scalar readings** reconstructs m connected
two-component carriers from one scalar readout; the linear count is minimal.

The v0.7 evidence: **93 regressions**, all **75 interaction graphs on one through
four channels**, **360 independent readout-closure comparisons**, **700 exact
projection intertwiners**, and **167 minimal transcripts**. Their decoders pass
1,298 basis-vector reconstructions. That version binds **37 source files**
and retains every earlier campaign. Interaction graph and strengths remain
declared; physical gravitational/quantum identification is still open. The
full invariant bilinear form of a connected family is definite, so a Lorentzian
geometric target requires a separate adapter.

## Added in v0.8: native response tensor and two cut ledgers

[U30–U32](NATIVE_RESPONSE_TENSOR_AND_CUT_LEDGER.md) construct
Q(Y)=Y^T Y-iY^T ZY from calibrated native response columns, with Z the
represented quarter-turn. The tensor preserves both squared response size and
an oriented pairing. It is covariant under a change of units only when its
metric is transported too.

A real coordinate cut always splits the real Gram form. Its **full complex
ledger closes for every response catalogue exactly when the cut projector
commutes with Z**. Otherwise an explicit imaginary cross-seam term remains;
the minimum paired repair has rank([C;CZ])-rank(C) extra scalar channels.
A half-channel cut needs one extra scalar. Preserving this tensor structure
is weaker than closing the full interaction dynamics.

For the signed order response L=UV-VU, CID-1's inherited total covariance gives
a separate exact statistical information ledger and a nested-cut tower law.
A two-channel response has L^T L=kappa I_4 with explicit positive kappa, so
its exact squared norm needs the same input faithfulness as its signed vector.
Zero mean, zero skew pairing and zero operator response are kept distinct.

Current native evidence: **111 regressions**, **490 tensor transport checks**,
**340 coordinate cuts**, **14 general readout repairs**, **90 exact response
scale/quotient checks**, **405 covariance ledgers** and **1,620 nested ledgers**.
The certificate binds **49 source files**, including the QTH correction.
The adjacent CID/QTH suite has **50 passing tests**; CID-1's original certificate
also passes read-only regeneration.

The ten-source audit confirms the relevant existing work and uses the corrected
thermodynamic paper at its pinned PR #4 commit. QTH-1 receives a targeted
correction: noncommuting SLDs can have zero state-averaged commutator.
Neither its declared quantum model nor a Lorentzian propagation law is
imported into the native construction. Gravitational identification is open.

## Reproduce

Python 3.11 or 3.12; standard library only. From the repository root:

```bash
python -m unittest discover -s papers/uncut-cut-measurement/tests -v
python papers/uncut-cut-measurement/certificate.py --check
```

The second command is read-only. It regenerates exact exhaustive checks in
memory and compares both the committed certificate bytes and SHA-256 pin.
It also binds the local manuscript, implementation, tests, source manifest,
claim ledger and research programme. An intentional revision uses
`certificate.py --write`, followed by review of the changed source and pin.

`CERTIFICATE.json` records **PASS_FINITE_CHECKS** with the enumerated domains.
The written proofs cover the hypotheses stated in the manuscript; Python
enumeration is not formal verification of those proofs. Source pins identify
the inspected repository versions, not an independent certification of them.

## Scope of the physical programme

Uncut ground is a proposed primitive, not identified with the finite carrier.
The carrier and its functions are declared mathematical representations.
No spacetime metric, external clock, probability law, quantum state space,
Hamiltonian or classical field equation is assumed by U1–U19. U17–U19 introduce
declared finite-dimensional linear carriers. U20 is a conditional geometric
adapter on a supplied smooth chart with a coframe and connection. U21/U22 are
algebraic coupling results; U23 constrains a supplied homogeneous coframe class
for a fixed connection on a supplied chart.
U24–U26 are clock-free finite-dimensional grading and observer results with
declared algebraic controls; no spacetime or physical field equation is used.
U27–U29 add native-compatible real interaction, classified exact cuts and
model-level response reconstruction under a declared control catalogue.
U30–U32 construct a native response tensor, its cut seam, and a classical
finite-ensemble covariance ledger with supplied positive weights. These weights
are not a physical probability law derived from the uncut primitive.
No physical quantum/classical correspondence, spacetime reconstruction, minimum
action law, universal curvature-information identity or quantum-gravity closure
has yet been derived by this package.

No merge of [correction PR #4](https://github.com/Parveen117/Publications/pull/4)
is required to run this package. The v0.8 source audit cites that correction
edition at its immutable commit, with its open-PR status recorded.
It uses pinned native sources and does not treat archived publication claims as
new premises. See [SOURCE_PINS.json](SOURCE_PINS.json).
