# Quantum gravity: canonical research record

Updated 16 October 2026. New smooth observable algebra and quantum orientation mixing: **GOA-1**. Critical-origin and two-scale construction: **COS-1**. Complementary clock and two-axis physical completion: **CCS-1**. Complete test-dual kernel and conditional measure uniqueness: **LSC-1**. Regular physical sector and geometric reduction: **LPS-1/2**. Local constraint algebra: **LQA-1**. Homogeneous comparison: **HQS-1**. Completed finite quadratic-observer verification: **CHD-1/2/3**. The full local physical quantum-gravity construction and physical FTL remain open. See the [closure ledger](RESEARCH_CLOSURES.md).

## Working agreement

The research belongs in this repository: derivations, code, source pins, certificates, rejected stronger claims and the next exact obligation. Chat reports summarize repository work; standalone documents do not substitute for it. Continue autonomously from the strongest verified endpoint, reuse previous results, and avoid parallel duplicate research. The goal is a rigorous quantum-gravity construction, not a growing count of loosely connected examples.

Preserve frozen predecessor packets. Distinguish an algebraic proof, a written analytic theorem, a selected model and a physical identification. No certificate count by itself closes a continuum or physical claim. Use Python 3.12 for orchestration and the existing native arithmetic/rewrite engines for the scientific fixtures.

## Current result

[GOA-1](papers/native-global-observable-algebra/THEOREM.md) classifies the
complete globally smooth classical Dirac value algebra of the selected
local source. On the two axes its values are c+P a(Q,P) and c+K b(Phi,K),
with arbitrary smooth a,b and one common critical value c. A periodic
mean-zero primitive of k-K gives a smooth full-local Hamiltonian lift with
exact canonical zero-mode brackets. The classification concerns physical
values; it neither identifies the full local constraint ideal nor quantizes
all arbitrary smooth functions.

Smooth classical observable flows preserve the branch and nonzero momentum
sign along finite regular trajectories. In the quantum theory the minimal
momentum/dilation algebra has five central blocks: P+,P-,K+,K-,origin.
This is not the complete observable algebra implied by the smooth source.
The source-allowed observable P[1+(1/2)cos(Q_sharp)], with the inherited
full-line canonical Q and symmetric reduced ordering, is self-adjoint on
the physical momentum domain. Its exact momentum-shift action couples the
two P orientations. A corresponding independent K observable couples the
K orientations; no local q/phi interchange symmetry is assumed.

Two normalized physical packets with identical momentum probabilities and
identical restrictions to the entire minimal affine algebra have new
observable readings -59/200 and -47/200: phase contrast -3/50. A positive
momentum packet has small-observable-parameter transition coefficient
1/768 into negative momentum. These are exact mathematical observer results,
not calibrated gravity predictions or physical time dynamics. The original
vacuum clock Hamiltonian stays zero.

Adding these two actual self-adjoint observable groups gives the exact
double-commutant completion B(H_P) direct_sum B(H_K) direct_sum C_origin,
with center C^3. Relative orientation phase is now measurable inside either
regular branch. P-to-K and regular-to-origin transitions are absent in this
constructed algebra; a universal no-coupling theorem for other quantizations
is not claimed. All previous physical norms and the two-scale distinction
are preserved.

[COS-1](papers/native-critical-origin-sector/THEOREM.md) resolves the
critical classical quotient and gives an explicit two-scale quantum option.
For P=K=0, r'=k has a periodic primitive and
Q_0=(1/2)log[(integral exp(q+r))(integral exp(q-r))] is gauge invariant.
The complete critical orbit labels are (Q_0,Phi); their restricted canonical
two-form is zero, but the original gauge flows do not identify them.
Globally continuous Dirac observables take one common value on all those
critical orbits. This is proved by exact gauge-equivalent regular sequences
and exposes the non-Hausdorff orbit quotient, without deleting its labels.

The same Gaussian average has the exact leading asymptotic
B_epsilon=log(1/epsilon) conjugate(e_0(f))e_0(g)+O(1), where e_0=(rho f)(0,0).
Division by the logarithm gives a positive one-dimensional origin quotient.
The maximal smooth-test domain with a finite unscaled average is e_0(f)=0;
its completion remains H_P direct_sum H_K. A filtered direct sum with an
additional origin line retains both rates using separate quotient maps.
No single scalar normalization retains both nonzero norms. No continuous
positive all-test form can retain the old regular weights and all scalar
momentum-adjoint rules; a positive subtraction construction explicitly
fails momentum symmetry. Relative origin normalization remains a choice.

There is also a source-derived observable algebra across the strata.
The Q_0 expression extends smoothly as Q_* with a chosen primitive and
coincides with LPS's Q on P!=0. D_P=P Q_* and D_K=K integral phi are global
weak Dirac observables, independent of the auxiliary convention on shell.
Their self-adjoint quantum realizations are iP partial_P and iK partial_K
on the active regular summands and zero on the origin. Logarithmic momentum
coordinates give explicit H^1 domains and affine commutators. These preserve
both regular sectors and give a meaningful observable contract for the
filtered option, without assigning a universal length to the origin line.

[CCS-1](papers/native-complementary-clock-sector/THEOREM.md) completes the
complementary P=0,K!=0 branch for the selected vacuum 1+1 representation.
The actual smooth periodic constraints force p=0 and phi'=0 whenever P=0.
For P!=0 they force both A currents nowhere zero and K=0, so the LPS patch
restriction is automatic on smooth solutions (not on off-shell data).

On P=0 the remaining q,k profiles transform as two real affine connections.
An explicit global circle-diffeomorphism straightening reduces them to
q=Q,k=K. For K!=0 geometry Q is a clock; Q=tau leaves the canonical pair
(Phi,K), reduced one-form K dPhi and zero reduced Hamiltonian. The scalar
is constant in space and in this clock; L(tau)=L_ref exp(tau) is the chosen
clock reading, not a new propagation law.

The universal oscillator dressing is exp[-sum U_n V_n/(2n)] per chirality,
independent of K. It solves every unchanged local constraint as a test-dual
state. Contraction followed by the declared PK average gives the positive
physical space H_K=L^2(dK/|K|), with self-adjoint Phi=i[partial_K-1/(2K)]
on the transported full-line H^1 domain. This supplies the physical norm
and a quotient observable intertwiner, separately from kinematic positivity.

One common scalar average on tests excluding the joint origin completes to
H_P direct_sum H_K. The full-line matching and relative branch normalization
remain selected quantization choices. At P=K=0 the same positive average
diverges logarithmically on nonzero-origin tests. A separate positive
rank-one origin form exists, but it is not the divergent average's value
and supports neither canonical clock. COS-1 now derives that form from a
logarithmically normalized limit and specifies the filtered option; a
single unchanged all-test extension remains obstructed.

[LSC-1](papers/native-local-sector-classification/THEOREM.md) closes the
remaining oscillator-sector classification for the fixed LQA representation.
An explicit operator D_N obeys Im<v,D_N v> >= N||v||^2 on every positive
oscillator level. Thus no extra primary exists and the vacuum descendants
span the entire Fock module for every real zero mode, without importing a
Kac determinant. This is an all-level written proof, not an extrapolation
from finite Gram checks.

On compact smooth zero modes times finite oscillator tests, every continuous
solution bra is F composed with a smooth descendant reduction rho, where
F=U(a_0)delta(b_0)+delta(a_0)V(b_0). The two axis contributions overlap only
by a transferable multiple of delta(a_0)delta(b_0). On P!=0 the pairing has
oscillator multiplicity one and all remaining freedom is a distribution
in P times delta(K). CCS-1 now supplies the complementary open axis with a physical norm,
clock and reduced observables; the joint origin remains separate.

LSC also separates measure freedom from representation freedom. Positive
rigging forms compatible with momentum multipliers are positive measures.
Symmetry of the fixed HQS Q=i[partial_P-1/(2P)] forces c_+/|P| and c_-/|P|
on the two signs. Its specified full-line H^1 matching forces c_+=c_-;
normalization then gives the LPS/HQS measure. This is uniqueness under the
inherited observable contract, not selection by the constraints alone.


[LPS-1](papers/native-local-physical-sector/THEOREM.md) constructs a nonzero
generalized all-mode solution family for LQA's unchanged reversed-vacuum
constraints. Writing each chiral constraint as A_m(u)-A_-m(w), with auxiliary
positive c=7 factors, gives a descendant pairing on the test-space dual.
The auxiliary intercepts cancel; the local constraint central charge stays
zero. The paired states have divergent Fock norm, but a declared oscillator
contraction followed by HQS regular-clock averaging gives a positive physical
quotient L^2(dP/|P|). An explicit surjective unitary intertwines HQS physical
states, norms and the specified self-adjoint observable domains. LSC-1 now
proves this oscillator pairing exhaustive on the regular
test dual and gives conditional measure uniqueness; CCS-1 supplies the
complementary open axis. COS-1 supplies a normalized critical-origin
quotient under a different averaging scale; uniqueness across quantizations
and a faithful quantization of stratum-only origin labels remain open.

[LPS-2](papers/native-local-physical-sector/LOCAL_GEOMETRY.md) then derives
a regular reduction of the full classical local constraints. On the patch
where both A_+=p+phi' and A_-=p-phi' are nowhere zero,
C=q-(1/2)log|A_+A_-| is spatially constant and weakly invariant under both
local gauge families. Independent density-straightening transformations
produce a homogeneous slice with Q=C+log|P|, K=0 and reduced one-form P dQ.
The relational length L_ref exp(Q) has precisely the positive reduced
quantization carried by LPS-1. Raw circumference is generally not invariant
under normal gauge transformations. This geometric identification uses
quantization after regular reduction; pointwise quantum inverse currents,
field exponentials and the uniqueness of Dirac quantization remain open.


[LQA-1](papers/native-local-quantum-anomaly/THEOREM.md) derives the chiral
current representation of the actual periodic V=0, matter-free LNC source.
The original null improvement does not cancel the c=2 anomaly in the
conventional oscillator vacuum. Equal left/right central charges cancel
in the DD and HH brackets but add in DH. No scalar zero-mode shift repairs
the unchanged strong all-mode constraints in that representation.

A different, explicitly selected vacuum reverses the negative-sign
oscillator ordering. Its two current contributions are +1 and -1, so the
original local constraints close without a central term. This representation
has a positive kinematic Fock inner product, a common dense invariant
algebraic domain, and the original homogeneous PK=0 zero-mode equation.
The null improvement is retained, not dropped. This adapts standard
free-current/Virasoro methods; it is not discovery of a new anomaly mechanism.

The physical distinction is essential: no nonzero strong common kernel
exists inside the finite-polynomial test domain. LPS-1 supplies its regular
generalized physical sector, norm, observable domains and HQS intertwiner.
LSC-1 now classifies the complete continuous test-dual kernel in this
representation; CCS-1 physically completes both open axes under the
declared averaging and observable contracts.
The mode grading is unbounded in both directions
and is not identified with physical energy. A separate linear-counterterm
family closes iff r*s=gamma=-1/12 (hbar=1), but shifts the homogeneous
equation to PK-pi/3=0; it is documented and not adopted. The reversed-vacuum
representation is the one used in LPS, LSC and CCS. All predecessor packets and
the finalized CHD manuscript remain frozen.

[HQS-1](papers/native-homogeneous-quantum-sector/THEOREM.md) quantizes the
regular homogeneous V=0, matter-free sector inherited from LNC/REC. The
canonical variables q=log(L/L_ref), p=L pi give C=p k. Product ordering
has a self-adjoint joint-spectral constraint and a complete unitary gauge
group. Averaging on the declared p!=0 test sector derives the positive
physical measure dp/|p|. An explicit unitary map to reduced quantization
provides self-adjoint length and dilation operators, positive length, and the exact
affine relation [L,P]=iL on a common core.

The physical momentum-space log-length operator includes the necessary
weight term i[partial_p-1/(2p)] and a selected full-line matching domain.
Ordinary half-line momentum instead has deficiency indices (1,0), so it
cannot be adopted as a self-adjoint observable in that scalar representation.
Both clock orientations are retained; the singular p=0 constraint sector
and alternative gluing choices remain separate. The reduced V=0 geometry
is static in phi, consistent with REC, not a propagating gravitational model.

The completed CHD observer theorem now acts on an explicit normalizable
physical packet span. Its first and second momentum compressions differ
from a closed finite dynamics by (1/12)I; the full momentum closure is
infinite. General self-adjointness, averaging and domain statements have
written proofs, with finite native arithmetic/rewrite certificates scoped
separately. This is a selected homogeneous quantum-sector closure, not
quantization of the full local constraint family.

**Manuscript closure, 16 October 2026.** The integrated eight-page paper
[Cut-Generated Information Transport and Complete Observer Verification](papers/native-cut-generated-history/manuscript/paper.pdf)
is finalized as version 1.0, with editable LaTeX, bibliography and a provenance
manifest. It closes the finite quadratic-observer verification study; it does
not close physical quantum gravity, the local quantum-constraint target or
the FTL state-law/arrival-time obligations. This manuscript supersedes the
earlier chapter as the current integrated write-up for this mathematical
scope. The historical chapter and all certified source proofs are preserved.
No journal submission or publication is implied.

[CHD-3](papers/native-cut-generated-history/CIRCULAR_OBSERVER_COMPLETION.md) closes the finite quadratic-observer verification project. A calibrated n-channel complex carrier has a minimum complete family of n^2 fixed scalar quadratic readings for an unrestricted Hermitian observable. The probes are axes and two-coordinate circular positions (3/5,4/5), (3/5,4iota/5); an explicit rational inverse recovers every matrix entry. Real-symmetric forms require n(n+1)/2 readings; known trace supplies one datum. These counts are not bit lengths or universal pure-state tomography bounds.

The same complete family certifies T*G_out T=G_in if and only if every probe defect vanishes, and positive metrics give the left decoder G_in^(-1)T*G_out. Applying CHD's coherent history W transports the complete observable certificate on its known image corner. The inherited thermo source reconstructs from (2,5,122/25,98/25); fixed classical branch readout loses the off-diagonal entry, detected by a -24/25 circular reading defect. Common-frame covariance and a deterministic reconstruction error bound are proved. Norm preservation alone neither recovers a nonlinear erased record nor forces a flat transport loop.

This settles a representation and verification project within declared finite quadratic/linear contracts. It rederives standard polarization and informational completeness in the native exact engines, with prior source results reused. It does not prove the observation contract applies to every physical process, select gravitational dynamics, or upgrade GSA's conditional signalling to physical FTL. General proofs, finite certificates and physical interpretation remain distinguished.

(Philosophical speculation: gravity imagined as a conscious presence awaiting observation, with information travelling through it faster than light. This is a speculative metaphor, not a derived or experimentally established claim of this research.)

[GSA-1](papers/native-gravity-signalling-audit/THEOREM.md) makes faster-than-light gravitational information transfer an active operational research target. The test is a sender-controlled difference in unconditional receiver statistics before a matched earliest-light-arrival benchmark. No speed is inferred from phase retention, entanglement or a dimensionless operation count alone.

The exact candidate is an explicit steering extension of the old EGR/QGC nonlinear relaxation. A single joint source and sender Z/X measurements produce QGC's equal-density ensembles. If each remotely prepared pure branch physically follows that nonlinear law without a delivered preparation label, the receiver probabilities at s*=(235/74)log(256/81) are 4/5 and 9/10, a contrast of 1/10. The initial probability contrast rate is 888/29375. This additional branch-sensitive physical premise conflicts with the admitted density-only mixture-equivalence contract; it is not established by the source theorem and is not adopted silently as the consistent quantum model.

Local trace-preserving operations, including coherent CHD histories, give zero remote contrast. The admitted MGH joint operation gives reference contrast 16/25 and optimized effect contrast 4/5, but its registers have no supplied separation/duration. GCP's existing local shared-cone theorem remains unchanged. The next signalling gate is to derive a physical branch/state or spatial coupling law from the native source, then supply localized observables and a calibrated light benchmark. GSA provides a conditional mechanism and exact controls, not physical FTL evidence or a completed spacetime theory.

[CHD-2](papers/native-cut-generated-history/BINARY_INVARIANT_SYNTHESIS.md) closes the declared finite positive-quadratic synthesis problem. A rank-r form has a minimal r-channel linear amplitude record with positive rational weights; a binary Pythagorean tree needs r-1 merges and depth ceil(log2 r). The record is faithful on the quotient by the form's kernel, and every form-preserving linear map descends to a weighted isometry. Additional cut operations can increase the minimum; CRC's established closure gives the exact enlargement, with a source example 1 -> 2. Channels retain phase/sign, unlike a total squared amount; these counts are not finite bit lengths for arbitrary real values.

The direct CHD-to-source compatibility test is now explicit: the same branch-local metric H is conserved by a projective binary cut exactly when [P,H]=0. CHD's full coherent observable W H W* remains recoverable even when that test fails. On the inherited thermo H, phase-related source values 122/25 and 74/25 both become 98/25 under the fixed branch-metric readout; coherent transport preserves them. Standard factorization, CRC closure and field-degree results are attributed. General root descriptors and binary storage are not identified with generation using square roots alone.

[CHD-1](papers/native-cut-generated-history/README.md) treats the uncut level as a declared scalar observable algebra, with no independent observable distinction. Linear unital transports act trivially there. A nontrivial projection supplies a binary algebra; successive specified projective cut interactions construct an isometric history `W`, with `W*W=I`, without an intervening source Hamiltonian, physical time or exponential. The entire represented input operator algebra is faithfully transported into the image corner and decoded by `W*`. This is an explicit complete coherent-record construction, reusing standard measurement dilation; it is not yet a full UGD numeral representation.

Its exact information boundary matters: a classical history of all outcomes and conditional states can lose cross-branch phases. On the inherited CLD fixture, different input phases give identical classical histories but distinct recoverable coherent records. Matched two-cut orders obey `Delta*Delta=4[B,A]*[B,A]`, evaluated as `(576/625)I`. The actual A-B-A projective protocol has final A-return probability `337/625`; a repeated A cut does not force change. Protocol choice and the projective implementation remain declared. The cut order supplies an algebraic obstruction, not a newly identified spacetime curvature tensor.

The resulting next interface is to derive or calibrate a cut protocol on the selected local matter/geometry sector and show that its retained maps preserve the physical constraints and observable algebra. Zero observable distinction is not geometric zero dimension, and does not assert the nonexistence of an inaccessible source. MZC and HCM remain available as sector/recovery tools. The pre-existing LNC/REC quantum representation gate below is not discharged by the new encoding.

## Earlier geometric and constraint results

[SRT-1](papers/native-singular-response-transport/README.md) retains the full axis-response record rather than infer it from a diagonal invariant. For a declared analytic chart `A(t)=diag(t^w_i)`, `C=A^-T H A^-1` extends analytically through zero exactly when each `H_ij` vanishes to order at least `w_i+w_j`. Its coefficient at that order gives `C_ij(0)`; higher boundary derivatives require correspondingly higher jets. A concrete two-axis family recovers cross-response from the first jet and second-axis response from the second jet. The collapsed point alone cannot do so.

The packet reuses NT's response element, CRC's kernel criterion and the existing tower contracts. Its new contribution is the explicit singular-chart reconstruction criterion and minimal boundary jet depths. It includes diagonal ambiguity, positive-response loss, wrong scaling order and omitted frame-term controls. The signed TVSP/EM compass map uses a partial Legendre transform to preserve the potential differential. Transport of electromagnetic field equations additionally requires a differential-chain intertwiner; the current two-axis adapter is not yet that full field construction.

These are analytic response-chart results. The weights have not been derived from cut dynamics; grade is not identified with quantum level, and chart rank loss is not identified with a physical curvature singularity. The latest step provides a concrete candidate record for testing against the local cut/matter model, with those identifications still open.

[REC-1](papers/native-relational-event-clock/README.md) reduces LNC's selected homogeneous, matter-free sector using the scalar phi as an event-clock, on `L>0`, `pi!=0`. Its Hamiltonian is `H_red=L V(phi)/pi`; `K=L/pi` and `M=pi^2+2 W(phi)` are weak Dirac observables for the single homogeneous constraint. Explicit complete observables at a clock reading have canonical bracket one. This does not reduce the full inhomogeneous constraint family.

The equations derive a normalized diagonal comparison. Two explicit potentials give the same endpoint readings and invariant but different interior readings. A differentiable pi record identifies `V=-(pi^2)'/2`; endpoint data identify an integral only. This answers the event-comparison proposal without inferring a unique law from a 45-degree drawing. IRT's general calibration result is reused. The `V=0` baseline has constant length/momentum and a running scalar clock.

The new periodic independent-coframe boundary choice is explicit. A varying proper length is incompatible with additionally imposing a smooth global single-angle cut of fixed winding; local cut charts alone do not resolve that obstruction. This global gluing obligation remains open. The packet also distinguishes algebraic closure from measurement completeness using a closed noncommutative matrix algebra with nonzero state uncertainty.

[LNC-1](papers/native-local-normal-closure/README.md) closes the three local classical secondary-constraint brackets in a selected 1+1 canonical coframe-plus-scalar model. Its geometry density is `h_g=e*pi*p_phi+(e'/e)*phi'-phi''+e^2 V(phi)`, with spatial density `d_g=-e*pi'+p_phi*phi'`. Combined with LCF's massive matter generator, it gives `{D,D}=D`, `{D,H}=H`, `{H,H}=D` with the displayed densitized smearings. The proof holds for arbitrary smooth potential and admissible compact smearings, not just the finite fixtures.

The geometry now has normal evolution `dot e|N=N e p_phi`. The existing matter source `-NJ` survives in coframe momentum evolution, and the local circle-cut variation remains `-(NJ)'/2`. The full canonical bracket supplies the geometry direction missing from LCF's matter-only normal action. Smooth solutions preserve the constraint surface under the stated boundary contract.

The new canonical pair, scalar and action ansatz are selected structure of a known dilaton-type construction; uniqueness from native cut information is not asserted. An exact control gives two closed theories with the same initial constraint data but accelerations `0` and `-1`. Symmetry/closure does not select the physical potential. Classical commuting matter fields, not physical quantized fermions, are used. Positivity, quantum anomalies and local quantum-field completion remain open.

All global-mode AGD/IRT/ZCR results remain frozen in their original scope. A physical 3+1 theory is not obtained by renaming this local 1+1 construction.

## Lineage to reuse

| Stage | Established content | Boundary to preserve |
|---|---|---|
| [CID-D2](papers/native-cut-transport-feedback/README.md) | Native cut and transport feedback from a selected conservative action | Finite carrier; not quantum gravitational constraints |
| [NJC-1](papers/native-jacobian-cut-certification/README.md) | Native Jacobian/cut/connection identities and retained free connection | A cut does not fix all ambient connection data |
| [EGR-1](papers/native-exponential-gravity-response/README.md) / [EGR-2](papers/native-nonlinear-gravity-continuation/README.md) | Regge stationary response and local nonlinear branch | Selected action, boundary and local hypotheses |
| [CRC-1](papers/native-cut-record-completeness/README.md) | Operation-dependent record completeness and source-law controls | A short invariant list need not close an arbitrary operation |
| [QGC-1](papers/native-quantum-gravity-consistency/README.md) | Ensemble obstruction and consistent quantum adapters | Old nonlinear relaxation is not an autonomous quantum channel |
| [N4H-1](papers/native-four-dimensional-hinge/README.md) | Actual four-dimensional Euclidean hinge metrics and lifted curvature | Lorentzian quantum dynamics not obtained by relabelling |
| [MGH-1](papers/native-matter-geometry-holonomy/README.md) | Independent metric register, geometric phase and exact quantum witness | Phase interaction and metric mixing selected |
| [CCL-1](papers/native-common-cone-limit/README.md) | Shared-refinement cone, local mixing and a conditional curved symbol | Shared scaling and physical dictionary explicit |
| [GCP-1](papers/native-geometry-controlled-propagation/README.md) | Geometry controls derivative, connection and curvature; reciprocal coupling and exact reduction obstruction | Geometry mixing and profiles selected |
| [AGD-1](papers/native-action-geometry-dynamics/README.md) | Inherited cut action fixes an infinite geometry generator under an explicit quantization contract; joint existence, source law and finite-section bound | Global circle geometry; local gravitational constraints, physical action selection and lapse profile remain open |
| [ZCR-1](papers/native-zero-reference-recovery/README.md) | Exact packet moments, finite observation-sector error and variable-profile flat-metric recovery | Joint preparation/action-scale limit; matter remains quantum; general classical gravitational field limit open |
| [IRT-1](papers/native-invariant-ratio-transport/README.md) | Coupled reflection-weight ratio, classification of the common geometry invariant algebra, and calibrated rate/diagonal theorem | Specified invariant algebra is not a complete state record or local gravitational constraint algebra |
| [CLD-1](papers/native-cut-local-dynamics/README.md) | Changing-cut memory, finite covariant closure/conditioning budget, and AGD moving-reflection criterion | Cut-relative locality is not spatial locality; infinite memory norm estimates and local gravitational constraints are not proved |
| [RCG-1](papers/native-relational-cut-space/README.md) | Based cut geometry, rank-one overlap/cycle reconstruction, common-transport residual and its exact geometry-change law | Standard reconstruction method credited; operational geometry does not select physical spacetime or an interacting gravitational action |
| [LCF-1](papers/native-local-coframe-covariance/README.md) | Cut-gradient coframe density, cancellation of mixed local residual, covariant matter action and local cut-source variation | Normal geometry transport, full constraints and quantum local-field completion remain open; regular spatial patch required |
| [LNC-1](papers/native-local-normal-closure/README.md) | Selected canonical local geometry action closes all three classical 1+1 brackets with massive matter; source and constraint preservation retained | Extra scalar/action inputs explicit; physical reduction, quantum representation/anomalies, positivity and 3+1 extension open |
| [REC-1](papers/native-relational-event-clock/README.md) | Homogeneous scalar-clock reduction, sector complete observables, derived comparison invariant and potential reconstruction from differentiated readings | Homogeneous matter-free branch only; clock failure and global winding explicit; full local quantum construction open |
| [SRT-1](papers/native-singular-response-transport/README.md) | Analytic singular-response extension iff graded divisibility holds; minimal boundary jet record, signed thermo/EM adapter and finite differential-chain controls | Weights and frame declared; relation to actual cut dynamics, full Maxwell fields and physical quantum gravity open |
| [MZC-1/2/3](papers/native-moving-zero-covariance/README.md) | Moving balance covariance, circular admissibility and maximal common equilibrium via the joint cut/sector record; applied to pinned RH definitions | Declared transport; no new zeta data or universal innovation bound; RH and physical QG remain open |
| [HCM-1/2/3](papers/native-heat-observer-memory/README.md) | Source heat memory, EMK phase records, coprime gain recovery and calibrated rational correction | Selected observer and finite invariant coefficient block; full interacting row closure, RH and QG completion not established |

The earlier classical MG/SM action results are pinned and included as source snapshots in CCL-1. Reuse their existing classical coframe/source variation instead of proposing it as new. They do not supply the missing positive interacting quantum theory.

## Next proof target

**Derive physical implementation or dynamics for the now-explicit observer
family, or extend the source beyond the vacuum baseline.** GOA-1 classifies
the globally smooth classical observable values and constructs a nonlinear
self-adjoint quantum probe with exact phase contrast. Reuse its three-block
completed quantum algebra and the COS two-scale rigging distinction.
A physical operation implementing that probe, a source-selected time law,
or P/K/origin coupling requires additional derivation. Arbitrary off-diagonal
Hilbert-space matrices do not establish such a coupling. Nonzero potential,
quantized matter and direct kinematic relational fields remain separate
construction targets.

CHD supplies a cut-generated record map instead of assuming pre-cut evolution. To integrate it here, specify cuts and their implementation from this model and prove constraint/observable compatibility on the retained image. Complete coherent transfer alone does not choose the cuts or a gravitational coupling law. No prime-phase extraction is required for this interface.

SRT supplies a response-record tool for this target. Its immediate integration gate is to derive the relevant axis weights and comparison frame from actual local cut/matter solutions, then test whether the resulting jet record obeys the local differential/constraint intertwiners. It does not replace the quantum representation obligation with a label change or a finite chart example.

LNC resolves LCF's normal-bracket obligation within its explicit classical extension. REC supplies a homogeneous classical reduction as a comparison target; it does not supply the full physical sector. The next work must:

1. Reuse LSC's all-level cyclicity and full dual-kernel theorem. Do not
   repeat a search for extra oscillator primaries in this fixed real-momentum
   representation; the theorem rules them out.
2. Reuse COS's two-scale obstruction and GOA's smooth value algebra.
   Orientation coupling is now constructed under the specified reduced
   ordering; inter-branch coupling still needs a source representative,
   physical domain and compatibility with the rigging/adjoint contract.
3. Develop physical implementation or dynamics for GOA's nonlinear
   observable family, using the LPS/CCS local reductions.
   Keep the smooth-field hypotheses and chosen quantum matching visible;
   reduced gauge-slice observables do not yet construct kinematic field operators.
4. Treat potential, physical matter, ordering/anomalies and global cut gluing
   as new obligations rather than consequences of the vacuum classification.

Physical 3+1 interactions and empirical calibration remain downstream obligations. The model-choice ledger must distinguish selected assumptions, derived identities and independent predictions.

## Replay status

The eight integrated predecessor packets (NJC-1, angular calibration, EGR-1, EGR-2, CRC-1, QGC-1, N4H-1, MGH-1) remain frozen. CCL-1 has 69 exact fixture checks and 3 controls. GCP-1 has 115 exact checks, 6 native word proofs/replays and 9 controls. AGD-1 adds 92 exact checks, 8 native word proofs/replays and 8 controls. ZCR-1 adds 169 exact checks, 11 native word proofs/replays and 13 controls. IRT-1 adds 103 exact checks, 9 native word proofs/replays and 16 controls. CLD-1 adds 167 exact checks, 15 native word proofs/replays and 14 controls. RCG-1 adds 108 exact checks, 9 native word proofs/replays and 10 controls. LCF-1 adds 95 exact checks, 8 native word proofs/replays and 17 controls. LNC-1 adds 66 exact checks, 4 native word proofs/replays and 14 controls. REC-1 adds 115 exact checks, 4 native word proofs/replays and 13 controls. SRT-1 adds 143 exact checks, 5 native word proofs/replays and 21 controls. MZC-1/2/3 adds 194 exact checks, 3 native word proofs/replays and 22 controls for moving-zero covariance, sector admissibility and common equilibrium. HCM-1/2/3 adds 265 exact checks, 3 native proofs/replays and 28 controls on the inherited YM54 heat source, phase-record adapter and common-gain recovery. CHD-1/2/3 adds 594 exact checks, 12 native word proofs/replays and 89 controls for cut-generated coherent histories, minimal binary synthesis and complete circular-observer verification, including exact reconstruction, minimum probes and transport defects. GSA-1 adds 150 exact checks, 3 native proofs/replays and 7 controls for the conditional steering mechanism and local/joint-operation signal audit. HQS-1 adds 147 exact checks, 10 native word proofs/replays and 14 controls for the selected homogeneous quantum sector and physical observer interface. LQA-1 adds 407 exact checks, 9 native word proofs/replays and 25 controls for the two local current representations, central terms, counterterm classification and positive kinematic alternative. LPS-1/2 adds finite descendant-pairing and regular local-reduction witnesses; see its certificate for exact counts. LSC-1 adds the all-level primary/dual classification and conditional measure proof, with finite native witnesses recorded in its certificate. CCS-1 adds 688 finite exact checks, 5 native word proofs/replays and 20 controls for the complementary clock, universal dressing, two-axis completion and origin obstruction. COS-1 adds 275 finite exact checks, 3 native word proofs/replays and 27 controls for the critical invariant, averaging scales, domain/adjoint obstructions and global dilation algebra. GOA-1 adds 192 finite exact checks, 3 native word proofs/replays and 14 controls for smooth observable closure, periodic lifts, phase contrast, orientation mixing and the algebra centers. The repository runner replays all 30 stages read-only; CI separately replays the source CID-D2 action packet. Analytic proofs and their assumptions are documented beside the finite certificates.

```sh
python3 scripts/replay_quantum_geometry.py
```
