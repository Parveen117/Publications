# Physics from Cuts: Uncut Ground, Measurement and Lawful Continuation

**Monty Dabas — working research programme, version 0.5, 25 September 2026.**

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
[sector and coframe constraints](SECTOR_AND_COFRAME_SELECTION.md), and [claim ledger](CLAIMS.md).
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

Current evidence: **63 regression tests**, all **646 connected K/R coupling
graphs on two through four labelled components** compared with a full
invariant-form solver, and **240 independent curvature-matrix comparisons**
at 48 supplied chart points. All earlier evidence is retained. The certificate
binds 27 source files; the new provenance manifest pins four additional sources.

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
No physical quantum/classical correspondence, spacetime reconstruction, minimum
action law, universal curvature-information identity or quantum-gravity closure
has yet been derived by this package.

This addition is independent of pending [correction PR #4](https://github.com/Parveen117/Publications/pull/4).
It uses pinned native sources and does not treat archived publication claims as
new premises. See [SOURCE_PINS.json](SOURCE_PINS.json).
