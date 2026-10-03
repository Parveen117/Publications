# Yang–Mills certified programme

## Current continuation: YM-54

[YM-54's response/protocol proof](YM54_RESPONSE_CURVATURE_PROTOCOL_BRIDGE.md)
connects NT's existing native curvature to an actual compact turn process.
Remove the trace response, factor the positive H and express its two
shape words as Y_i=p_i K+q_i L. The shape Gram G has det G=16 f^2.
The compact lift iota Y_i supplies four equally counted signed turns;
their second moment is the three-by-three tensor C=VV^T/2, while
G=2V^TV. Consequently

\[
\beta_{\rm protocol}=\lambda_{\min}(G)/4\geq4 f^2/\tau,
\qquad |\theta|<f_0^2/(600 T)
\]

is a sufficient YM53 interaction window when |f_i|>=f0>0 and tau_i<=T.
The original curved thermo fixture gives beta_floor=5/234 and window
abs(theta)<1/112320. At theta=1/262144, J=6, R=5/4 the inherited
interacting rate exceeds 0.000079467076 in the declared heat clock.

[The certificate](certificates/YM54_RESULT.json) and
[source pins](certificates/YM54_SOURCE_PINS.json) bind seven written results,
48 exact response jets, 16 independent coefficient generator replays,
36 outward finite-turn refinement enclosures and the original potential's
derivative replay. Twenty new and 171 related tests pass on **Python 3.12**.
All 110 upstream pinned sources remain unchanged. The general proofs are
written, not mechanically formalized or externally expert-certified.

The source jet, two directions, equal signed counts, compact lift and
heat clock remain explicit selections. Same curvature with growing
shape trace can lose a uniform floor; inserting idle records changes
the clock without changing response curvature. This chapter constructs
a concrete source for YM53's hypotheses, not physical selection of the
protocol or a physical mass. State-dependent diffusion, anisotropic
joint/volume limits, actual row closure and 4D/Clay remain open.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym54_response_protocol.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym54_response_protocol.py -v
~~~

## Previous continuation: YM-53

[YM-53's anisotropic interaction proof](YM53_ANISOTROPIC_INTERACTING_GAP.md)
extends the declared open chain to site-dependent positive semidefinite
native protocol tensors. Their ordered eigenvalues obey b_i>=beta>0;
a_i may be zero and the orientations need not agree. The sufficient window

\[
|\theta|<\beta/2400
\]

provides an explicit full-source gap uniform in every finite width and
0<h<=1/beta. An all-content coarse kernel bound replaces isotropy in the
temporal-block argument. Short conditional bridges with zero normalizer
are excluded lawfully, with grouped endpoint comparisons where required.
The actual square-sourced transfer, fixed-width norm time limit and
interacting ground-source derivative energy carry the same rate.

For C=diag(0,3/2,3/2), theta=1/4096, J=8 and R=4/3,
gamma=3 log(4/3)/160>0.00539403. This is a rigorous conservative lower
bound in the declared chain, not a physical mass or a phase boundary.

[The certificate](certificates/YM53_RESULT.json) and
[source pins](certificates/YM53_SOURCE_PINS.json) bind seven written results,
18 new tests, 108 finite spin inequalities, exact infinite-tail algebra,
eight inhomogeneous two-site polynomial controls, twelve norm-gap
transport cells and independent support/ground-energy fixtures. YM45's
coupling controls are reused. All 151 related tests pass on **Python 3.12**.
The 103 upstream pinned files are unchanged. General proofs are written,
not mechanically formalized or externally expert-certified.

The anisotropic infinite-volume/joint-cutoff construction and actual row
closure remain open, as do physical clock/protocol/state/action selection,
NCG quantum measure and 4D/AF/Clay. The old isotropic results retain their
original scope and certificates.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym53_anisotropic_interaction.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym53_anisotropic_interaction.py -v
~~~

## Previous continuation: YM-52

[YM-52's energy and observer proof](YM52_ENERGY_NOISE_AND_OBSERVER_GAP.md)
uses YM51's general symmetric protocol tensor C. It derives

\[
\Phi(f^2)=\Phi((E_t f)^2)+\Phi\!\left(E_t(f^2)-(E_t f)^2\right),
\qquad
\partial_t\Phi((E_t f)^2)=-2\mathcal E_C(E_t f).
\]

Thus discarded record variance accounts exactly for the mean's lost
squared distinction. A defined relative entropy satisfies
H'(rho_t)=-Phi(Gamma_C(rho_t)/rho_t) for bounded positive densities.
This is not a physical heat/work or energy-unit identification.

The weighted native frame brackets derive the tensor cof(C). With
0<=a<=b<=c its exact optimal centered relaxation rates are

\[
\gamma_{\rm full}=\min\{(a+b+c)/4,a+b\},\qquad
\gamma_{\rm even}=a+b.
\]

Rank two is enough for positive relaxation. The old YM51 family has
exact full rate min(3/4,2 eta). The classical positive-definite formula
is credited to Lauret; finite native spin and polynomial energy controls
support the written all-content proof, including singular tensors.

For fixed trace s, the full gap is maximal for every c<=3s/4.
The even-observer gap is maximal only at C=(s/3)I. A deficit epsilon
below its maximum bounds the tensor's operator-norm anisotropy by
2 epsilon. This conditional selection criterion does not fix the
physical observer objective or clock budget.

[The certificate](certificates/YM52_RESULT.json) and
[source pins](certificates/YM52_SOURCE_PINS.json) bind six written results,
16 new tests, independent word-variance controls, 108 finite spin gap
inequalities, 40 polynomial energy matrices and six rational entropy
enclosures. Current verification uses Python 3.12 only.
The fixed isotropic interacting chain proofs remain unchanged.
YM53 supplies the declared anisotropic chain extension. Physical selection
and 4D/Clay remain open.

## Previous continuation: YM-51

[YM-51's audit and written proof](YM51_TOOL_DEPENDENCY_AND_HEAT_SELECTION.md)
checks the tools used by the current YM route. The
[dependency ledger](certificates/YM51_DEPENDENCY_LEDGER.json) records 29
nodes, nine explicit selection assumptions and three open extensions.
It preserves the existing native algebra, typed tensors, connection
curvature, observer/memory corrections, response metrics and completion
results. It does not recommend redefining every tool.

The decisive witness is a family of symmetric native turn protocols with

\[
C_\eta=\operatorname{diag}(3-2\eta,\eta,\eta),\qquad
L_C=-\sum_{a,b}C_{ab}D_aD_b,\qquad 0<\eta\leq1.
\]

Every member has the same Phi_Q, native derivative brackets and decay
3/4 on all linear coefficients. But the same nonzero centered quadratic
source has decay 2*eta. Every axis is sampled at each eta>0; nevertheless
no positive decay rate is uniform as eta tends to zero. Curvature,
reference invariance and linear probes do not select the heat law or gap.
This does not change YM45's gap for its fixed isotropic chain.

For protocol moments m2 and m4, the native refinement error on degree d is
d^4*t^2*(m4+3*m2^2)/(96*n), recovering YM50's bound in its special case.
Covariance transforms C along with the frame. An additional fixed-law
invariance under four native conjugations forces C=cI. Six quadratic
readouts recover unrestricted symmetric C; isotropy still leaves the
clock scale c, and clock rescaling changes the relative interaction theta/c.
These internal symmetries do not derive local gauge Ward laws or physical
spatial isotropy.

[The certificate](certificates/YM51_RESULT.json) binds six written results,
25 protocol factorizations and 25 all-source finite energy controls,
70 reference checks, twelve bracket identities, sixty covariant actions,
24 linear-blindness checks, six quadratic witnesses, eighteen tensor
recovery checks, fourteen independent labelled counts, 72 outward
refinement enclosures and twelve clock-rescaling controls.
Seventeen new and 117 related tests pass on **Python 3.12 only**.
The test suite also rejects nine invalid dependency-ledger mutations.
All 88 pinned upstream files and the canonical engines are preserved.

**Scope:** the audit covers the current YM43--50 tool route. Fixed-law
isotropy, physical clock/state/action selection, general RH E5C/E6,
the actual chain's row closure, NCG quantum measure and 4D/AF/Clay/QG
remain open. Alternative free protocols are not assigned the old
interacting gap. Written proofs, exact controls and external expert
certification remain distinct; neither formalization nor expert review
is claimed.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym51_heat_selection.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym51_heat_selection.py -v
~~~

## Previous continuation: YM-50

[YM-50's written proof](YM50_NATIVE_REFERENCE_HEAT_BRIDGE.md) constructs
the compact reference and heat law from specified native counted records.
It uses the already derived EMK quaternion sector and canonical
stationary-path Phi_Sigma. On finite diagonal record kernels,

\[
\varphi_N(p)=\frac{\Phi_\Sigma(D_{p,N})}{N12^{N-1}},
\qquad
|\varphi_N(p)-\Phi_Q(p)|\leq\frac{2B_p}{N}.
\]

Finite coefficient energy and a solved Poisson equation derive positivity,
invariance, uniqueness and exact moments. Only then does comparison show
Phi_Q=Phi_ref on coefficient polynomials and their uniform completion.
The resulting recognition completion is isometric to the old coefficient
carrier, with bounded observables intertwined. Hilbert space is this
positive completion's representation.

A six-turn native protocol also constructs the full coefficient heat law:

\[
\|S_{\sqrt{6t/n}}^n-\operatorname{Exp}_\Sigma(-tL)\|_d
\leq\frac{3d^4t^2}{8n}.
\]

The existing bounded face weights, finite-chain transfers and their
declared limiting coefficient histories agree under the representation.
[The certificate](certificates/YM50_RESULT.json) binds eight written results,
24 independent word-transfer checks, 18 padded trace checks, two complex
pairings, five exact stationary spaces, 20 telescoping cases, 160 tail
checks, 70 independently projected moments, 495 earlier-reference matches,
140 invariance checks, 36 multiplier controls, 70 native Casimir identities,
18 outward heat refinements and twelve refusal groups.
Seventeen new and 100 related tests pass on **Python 3.12 only**.

**Scope:** the sector, record protocol, independent reference copies,
isotropic heat clock, chain and kappa(a)=theta*a remain specified choices.
The [origin ledger](certificates/YM50_ORIGIN_LEDGER.json) records them.
Raw endpoint records are not recognition-Cauchy; distinct paths may have
the same observed endpoint. General RH E5C/E6, physical selection, the NCG
quantum measure, the actual YM-49 row defect, full bounded-history density,
4D, AF, Clay and QG remain open. The proof is written, not mechanically
formalized or externally expert-certified. All 78 pinned upstream files
and the shared RKF/RH engines are unchanged.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym50_native_reference.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym50_native_reference.py -v
~~~

Earlier admitted-reference statements retain their historical scope.
YM-50 supplies the above compact bridge; the general integration and
physical identification obligations remain separate.

## Previous continuation: YM-49

[YM-49's written proof](YM49_TIME_ZERO_OBSERVABLES.md) constructs bounded
time-zero observables on YM-48's reflected carrier. A uniform energy
multiplier bound controls the time-zero collision. It proves the norm
bound, dagger/product laws and ordered coefficient-history readouts:

\[
\|M_0(F)\|\leq\|F\|_\infty,\qquad
M_0(F)M_0(G)=M_0(FG),\qquad M_0(F)e=JF.
\]

The instantaneous-row compression has an exact positive closure test:

\[
C(2t)-C(t)^2=L(t)^\dagger L(t),\qquad
L(t)=(I-JJ^\dagger)T(t)J.
\]

All these defects vanish exactly when the coefficient-history sector is
its instantaneous row. Otherwise the retained complement obeys RKF's
exact memory equation, now with gap-controlled tails and a lawful
bounded Schur inverse. A same-chain source proves
||[T(t),M_0(x_0)]e|| >= (1-exp(-gamma t))/2 > 0.
This is an observable/time ordering result, not an identification with
native gauge curvature.

[The certificate](certificates/YM49_RESULT.json) binds eight written
results, 48 independent time-zero path readouts, nine complex products,
three dagger and three norm checks, three all-source form bounds,
18 weighted product-energy controls, 35 collision/Cauchy controls,
nine mixed-time defects, twelve positive memory-coefficient bounds,
36 exact memory recursions, 108 memory-tail checks, three Schur inverses,
twelve Schur-tail checks and twelve refusal groups.
Seventeen new and 83 related tests pass on **Python 3.12 only**.

**Scope:** the interacting chain's row defect is still unevaluated.
The certificate's closing/nonclosing finite observers are controls,
not evidence deciding that infinite-chain question. Coefficient-history
density in the larger bounded-history completion is also separate.
The admitted abs(theta)<1/1680 heat/reference-functional carrier is
unchanged; native measure/NCG, physical gauge observables, clock,
relativistic/spatial continuum, AF, Clay and QG remain open.
The infinite proof is written, not mechanically formalized.
All 60 inherited source hashes and the canonical RKF engine are preserved.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym49_time_zero.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym49_time_zero.py -v
~~~

## Previous continuation: YM-48

[YM-48's written proof](YM48_REFLECTED_TIME_GENERATOR.md) constructs
continuous heat-time action on the reflected positive-form completion of
YM-47's histories. Literal translation descends through the null space,
is contractive and composes, and extends strongly continuously.
Its positive self-adjoint semigroup has a nonnegative self-adjoint generator:

\[
R_\lambda=\int_0^\infty e^{-\lambda t}T(t)\,dt,\qquad
\mathcal D(H)=\operatorname{Ran}R_\lambda,\qquad
HR_\lambda=I-\lambda R_\lambda,\qquad H\geq\gamma Q.
\]

The resolvent has explicit time-cut and quadrature tails. Its range on
finite histories is a graph core. A local coefficient energy identity
also constructs an isometric time-zero row embedding. The centered
fundamental character at heat time 1/4 has reflected norm squared at
least **5/32**, proving a nonzero excitation in this interacting carrier.
The full prior abs(theta)<1/1680 window is retained.

[The certificate](certificates/YM48_RESULT.json) binds nine written results:
48 independent reflected path identities, 30 mixed-source gap/composition
checks, four translated null controls, eight two-sided inverse checks,
six resolvent identities, 24 generator-domain/gap checks, twelve outward
quadratures, 35 coefficient identities, 108 integration-by-parts checks,
twelve weighted energy controls and twelve refusal groups.
Fourteen focused tests use Python 3.12 only.

**Scope:** the carrier is the reflected history completion. YM-49 now
supplies bounded time-zero coefficient observable action and a closure
diagnostic. The row embedding is not yet proved onto or invariant, and arbitrary future
multiplication does not preserve reflection-null histories. Native
measure/NCG identification, physical gauge-observable selection, material
time and the full 4D continuum remain open. The nonzero coefficient source
is not a physical particle or Clay nontriviality certificate. The general
proof is written, not mechanically formalized; the canonical RKF engine
and all prior evidence remain unchanged.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym48_reflected_time.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym48_reflected_time.py -v
~~~

## Previous continuation: YM-47

[YM-47's written proof](YM47_JOINT_LOCAL_HISTORY_LIMIT.md) constructs the
joint limit of local vacuum history readouts as the spatial interval grows
and the time step decreases, with no restriction on their relative rates.
It applies to the actual square-sourced S chain throughout
abs(theta)<1/1680. Both iterated limits give the same positive functional
on the stated history algebra and its uniform completion.

Normalizing each step first improves the finite-width error to

\[
\|(S_{I,a}/\lambda_{I,a})^n-V_I(na)\|
\leq6b_I(na)\sqrt a\,e^{2b_Ia}.
\]

An explicit vacuum estimate, separated-time rounding, summable one-face
extensions and a logarithmic crop give the joint local Cauchy bound.
The correlation gap and reflection positivity pass to the limit.
The [certificate](certificates/YM47_RESULT.json) retains four parameter
cells and checks 16 spatial budgets, six rational vacuum comparisons,
eight independent three-time readouts, 17 positive-contraction bounds,
27 off-grid time roundings, 54 exact tail identities, 490 asymmetric
crops, 1,120 one-face width checks and ten refusal groups.
Twelve focused tests use Python 3.12.

**Scope:** this is a positive history functional in a declared heat
parameter. YM-48 now constructs its reflected continuous-time action and
generator, with time-zero continuity for local coefficient row sources.
Full instantaneous-row identification, arbitrary bounded time-collision
claims, a countably additive path measure, physical clock, native measure/NCG
dictionary and 4D continuum remain open.
Written proofs and finite exact controls are distinct from mechanical
formalization. The prior evidence is frozen; YM-47 explicitly clarifies
the missing plus sign in YM-46's displayed boundary equation (1).

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym47_joint_history.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym47_joint_history.py -v
~~~

## Previous continuation: YM-46

[YM-46's written proof](YM46_INFINITE_VOLUME_STATE.md) constructs the
infinite-volume positive local-observable state at each fixed 0<a<=1
on the declared full heat chain, for abs(theta)<1/1680. Explicit spatial
and temporal boundary tails prove independence of boundaries/exhaustion
and uniqueness within the stated local conditional specification.

The normalized time map converges uniformly on every local row observable,
with its discrete semigroup law proved. A normalized heat half-step bridge
constructs the actual square-sourced S ordering as well as the convenient
T ordering. The completed positive forms retain YM-45's all-source rate:

\[
\|Q(P_a^\sigma)^nQ\|\leq e^{-\gamma an},\qquad \sigma=T,S.
\]

Regulated site- and bond-reflection positivity also pass to the limit.
The [certificate](certificates/YM46_RESULT.json) binds eight written results:
four joint rate cells and explicit Cauchy boxes; 27 weighted layouts,
216 interior and 306 boundary controls; 24 independent strip readouts;
16 local-specification conditionals; 84 normalized bridge/entry checks,
40 path identities and twelve half-edge skeletons; twelve reflection
identities, 36 mixed-source powers and twelve refusal groups.
Twelve new tests and 40 related tests pass on Python 3.12.

**Scope:** the boundary-error constants depend on a. YM-47 now identifies
this fixed-step volume construction and YM-44's time construction for
local vacuum history readouts. The state is a positive functional
on local observables and their uniform completion; a primitive-derived
measure, physical real-time evolution, native NCG dictionary and full
4D continuum are not claimed. Finite fixtures control the written proof;
they are not SU(2) replacements or mechanical formalization.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym46_infinite_volume.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym46_infinite_volume.py -v
~~~

**Runtime policy: Python 3.12 only.** Earlier evidence bytes are preserved.

## Previous continuation: YM-45

[YM-45's written proof](YM45_TEMPORAL_BLOCK_UNIFORM_GAP.md) establishes
the missing fine-time gap in an explicit small-bridge window:

\[
|\theta|<1/1680,\quad\kappa(a)=\theta a
\quad\Longrightarrow\quad
\|Q_a(S_a/\lambda_a)Q_a\|\leq e^{-\gamma a}
\]

for every finite spatial width m and every 0<a<=1, with gamma>0 independent
of both. All sources and full coefficient contents are covered. For example,
at abs(theta)=1/4096, gamma=Log(4/3)/56>=0.005137179865.
The same rate passes to YM-44's normalized time-limit operator for all t>0
at each finite width.

Eight written results supply conditioned bridge coupling, bounded tilt
costs, clipped-block counting, weighted memory comparison, fine-step vacuum
construction, all-source extraction, the open parameter interval and
time-limit transport. The [certificate](certificates/YM45_RESULT.json)
records 75 clipped layouts, exact conditioned path/marginal controls,
two joint block updates with preserved marginals, four rate cells,
24 time-limit controls and twelve refusal groups. Ten targeted tests pass.
These are written proofs plus exact finite evidence, not a formal-assistant
verification of an infinite theorem.

**Scope:** the reference functional and linear trajectory remain admitted.
The result does not cover the old theta=1/16 example, arbitrary interactions,
the original truncated Wilson family or the full 4D gauge lattice.
YM-46 supplies the infinite-volume state/time map at fixed a, and YM-47
joins the local-history limits. YM-48 supplies the reflected time
action/generator; full row identification, native measure/NCG identification
and the physical continuum remain open.
The small per-step gap
vanishes as a tends to zero; the positive rate per unit-a stays bounded below.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym45_temporal_blocks.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym45_temporal_blocks.py -v
~~~

**Runtime policy: Python 3.12 only.** All prior certificate bytes are preserved.

## Previous continuation: YM-44 (3 October 2026)

[YM-44's written proof](YM44_INTERACTING_TIME_REFINEMENT.md) constructs the
full interacting time limit on the declared trajectory kappa(a)=theta a.
For every fixed finite width m, every source and all contents,

\[
\|S_{h_n}\cdots S_{h_1}-U(t)\|
\leq 3bt\,e^{bt}\sqrt{\max_j h_j},
\qquad b=|\theta|(m-1),\quad t=\sum_j h_j.
\]

The heat/bridge commutator bound follows from two exact character identities.
The generator on the finite-content core is -L+theta M_V. Ordered insertion
tails, a cubic order-residue identity and a correctly normalized conditional
gap-transport gate make the refinement assumptions and remaining gap explicit.
YM-9's fixed-graph result, YM-22's leading tiling rate and EMK-T2's nilpotent
transport are reused at their original scopes.

The [certificate](certificates/YM44_RESULT.json) records 136 character-moment
enclosures, 30 ordered-word identities, 60 refinement cases, 48 independent
exponential-entry comparisons, unequal partitions and ten refusal groups.
Ten targeted tests pass on Python 3.12. General proofs are written, not
mechanically formalized; finite two-state controls are not an SU(2) substitute.

**At the YM-44 stage:** the gap uniform in width and time cutoff was missing.
YM-45 now supplies it only in its stated small-bridge window. A physical
trajectory, the native measure/NCG dictionary and the interacting continuum
remain open. The refinement error depends on width; YM-43's three coarse
cells do not by themselves furnish the fine-step gap premise.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym44_time_refinement.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym44_time_refinement.py -v
~~~

**Runtime policy: Python 3.12 only.** Previous certificate bytes are preserved.

## Previous continuation: YM-43 (3 October 2026)

[YM-43's proof](YM43_NATIVE_TIME_TRANSFER_BOUND.md) derives the missing
time-direction estimate on the existing YM-19 coarse heat-kernel chain.
The criterion `2 c_s + c_t (q + 1/q) < 1` gives a bound on **every**
vacuum-orthogonal source, uniform in all finite chain widths and all
transfer powers, within the declared positive-functional carrier.

| a | kappa | normalized time-operator ratio <= |
| --- | --- | --- |
| 6 | 1/16 | 1/2 |
| 8 | 1/8 | 1/4 |
| 12 | 1/8 | 1/64 |

Seven written results replace the cited comparison/extraction step with
explicit positive-sum and cut-square arguments. The certificate verifies
the parameter hypotheses, 2,964 strip-barrier checks, 64 common-part cases,
192 path/transfer identities, 360 whole-row marginal comparisons, finite
transfer controls and eight refusal groups. This is not a formal-assistant
proof. The measure dictionary, original YM-42 Wilson-grid operator problem,
cutoff uniformity, interacting continuum and Clay problem remain open.

```bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym43_native_time_transfer.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym43_native_time_transfer.py -v
```

**Runtime policy: Python 3.12 only.** YM-42's workflow also now uses this
single version; its runtime/proof input hashes are rebound with every
mathematical output unchanged. Previous dual-runtime evidence is historical.

## Previous continuation: YM-42 (2 October 2026)

[YM-42's proof and scope correction](YM42_SILENCE_TIME_LADDER.md) continue
the already-named t=3 step on YM-37's 52-dimensional even invariant carrier.
A general rational criterion gives a structural silence ball and a geometric
Smriti tail for every spatial step. At kappa=1/8, 1/4, 1/2 the certified
three-rail contraction ceilings are 0.007132026236, 0.028765231931 and
0.118984522771. These are declared truncated rational instances.

The exact counterexample in YM42-T5 proves that spatial contraction alone,
even uniform in rail count, cannot establish a uniform time gap. The operator
estimate on its original truncated Wilson grid and the three-rail content
release remain open; YM-43 supplies a different, scoped heat-kernel route. T75's
two-rail release retains its declared higher-content Gram-uniformity.
The broad historical 'content gone / t only' shorthand is superseded by
this scoped ledger; previous finite certificate bytes are unchanged.

Reproduce without rewriting evidence:

```bash
python papers/yang-mills-certified-benchmark/certificates/ym42_silence_time_ladder.py --check
python -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym42_silence_time_ladder.py -v
```

The dedicated `ym42-silence-time-ladder.yml` workflow checks Python 3.12 only.
The original v1 paper and historical capsules keep their versioned scopes;
the current continuation and correction are in the linked written proof.

## YM-1: historical starting point

First rung of resuming the MP Yang–Mills adapter with the certified RNKE
machinery: the reduced finite-cutoff gap of the SU(2) one-holonomy Wilson
benchmark, previously a float, is now an exact two-sided rational enclosure.

## Certified statement

For the normalized SU(2) Wilson convolution operator on class functions
(source: `MP/adapters/yang_mills/cutoff_gap_benchmark.tex`), with exact
character spectrum `lambda_j(beta) = I_{2j+1}(beta)/I_1(beta)`:

```text
Delta_red(a=1, beta=2) = -log( I_2(2) / I_1(2) )
  in [0.83672330623158891305006780105454768321... ]  (width < 1e-30, exact ℚ)
```

Verdict engine: directed rational interval arithmetic only. No floats in any
verdict. Controls: (C1) `lambda_1` strictly inside (0,1) hence gap strictly
positive; (C2) tampering `I_2 -> I_3` produces a separated bracket; (C3)
width budget.

## Claim boundary (fail-closed)

- CERTIFIED: the reduced finite-cutoff benchmark number only.
- NOT CERTIFIED / OPEN: full-lattice gap, theta-graph interacting transfer,
  ultraviolet/infinite-volume uniformity, OS continuum reconstruction, the
  Clay existence and mass-gap predicate. Adapter verdict remains `hold`.

## Ledger

| Item | Status |
|---|---|
| **YM-51 tool audit and native heat-selection boundary** | PROVED for specified symmetric turn protocols: same curvature/reference and linear decay do not select quadratic decay; fixed-law invariance forces C=cI; six-channel recovery and explicit refinement bound; physical selection OPEN |
| **YM-52 energy, records and observer relaxation** | PROVED on the compact free protocol carrier: energy/noise and bounded-density entropy balances, cofactor bracket tensor, exact full/even rates including rank two; fixed-budget even-objective isotropy selector; physical objective/units and anisotropic interacting transfer OPEN |
| **YM-50 counted native reference, positive completion and symmetric-turn heat bridge** | PROVED for the specified quaternion/readout protocol; coefficient reference, observables, declared transfers and closed coefficient-history sector intertwined; general UGD, physical selection, NCG quantum measure and 4D OPEN |
| **YM-49 bounded time-zero coefficient observable representation, ordered readouts, positive row-closure defect and gap-controlled memory; nonzero time/observable commutator** | PROVED on the declared heat-functional carrier; YM-50 supplies its compact reference/heat bridge; actual chain row defect, full bounded-history density, physical observable/DICT and 4D continuum OPEN |
| **YM-48 reflected continuous-time action, self-adjoint generator/domain/core and gap; coefficient-row embedding and nonzero source** | PROVED on the declared small-bridge heat-functional carrier; YM-49 supplies time-zero observable action; instantaneous-row/physical identification, DICT and 4D continuum OPEN |
| **YM-47 joint volume/time limit of local vacuum histories at arbitrary relative cutoff rates; both iterated limits coincide** | PROVED under the declared abs(theta)<1/1680 heat-functional adapter; positive history functional, correlation gap and reflection positivity retained; YM-48 constructs the reflected time action/generator; full row identification, DICT and 4D continuum OPEN |
| **YM-46 infinite-volume positive local state and normalized T/S time maps at fixed a, with inherited gap and regulated reflection positivity** | PROVED under the declared small-bridge heat-functional adapter; YM-47 joins local-history limits; exact controls PASS; DICT and 4D continuum OPEN |
| **YM-45 temporal-block gap, uniform in finite width and fine time step; same rate in every fixed-width time limit** | PROVED for abs(theta)<1/1680 on declared kappa=theta a full heat chain; YM-46 supplies fixed-a volume construction, YM-47 the joint local-history limit; YM-48 constructs the reflected time action/generator; full row identification, DICT and 4D continuum OPEN |
| **YM-44 full interacting time refinement at every fixed finite width, with an all-content operator-norm tail, order-residue control and conditional normalized gap transport** | PROVED under declared heat/functional adapter and linear trajectory; YM-45 supplies its gap premise in a small-bridge window; DICT OPEN |
| **YM-43 local overlap to full time-operator bound; all finite chain widths and transfer powers in three coarse heat-kernel cells** | PROVED under declared positive-functional carrier; exact hypotheses/controls PASS; cutoff/DICT OPEN |
| **YM-42 structural silence ball, exact t=2/3 continuation, all-spatial-step Cauchy tail, and space/time non-implication witness** | PASS (scoped rational grid; its time operator gap OPEN) |
| YM-1 reduced gap enclosure (a=1, beta=2) | PASS (pinned) |
| YM-2 interacting theta-graph gap, certified small coupling (kappa < Delta_red/6) | PASS (pinned) |
| YM-3 first-order crossing direction: rank-one, A<->B transported, slope lambda_half/4 | PASS (pinned) |
| YM-4 symmetry-protection theorem [S,T_kappa]=0 + finite-kappa certified lower bounds | PASS (pinned) |
| YM-5 certified TWO-SIDED gap at finite kappa, beyond the sandwich threshold (Kill-Lemma realization) | PASS (pinned) |
| YM-6 seam-integer dock: EXACT eigenvalue counts on the native 5x5, certified to kappa=1/2 (D03 congruence) | PASS (pinned) |
| YM-7 V7 carrier: kappa=7/10 unlocked (exact count), certified 7-eigenvalue crossing curves | PASS (pinned) |
| YM-8 CAPSTONE: theta-graph gap is a THEOREM at every coupling (certified kernel floors + pinned Jentzsch anchor); remaining Millennium content = uniformity in the cutoff | PASS (pinned) |
| YM-9 FIRST UNIFORMITY: exact heat-kernel refinement family; free gap exactly cutoff-independent; interacting uniform bound 3/8 along a declared trajectory | PASS (pinned) |
| YM-10 blindness ledger (SPECTRAL-1/2 framework): exact multiplicity law, exact probe counts, composition ledger, YM-9 T4 correction | PASS (pinned) |
| YM-11 three gate verdicts on the toy carrier: gauge CLOSED; volume/IR SPLIT (free closed, sandwich route proven insufficient, critical volume exact); universality SPLIT (counting closed, metric open) | PASS (pinned) |
| YM-12 the T03 move: positivity SQUARE-SOURCED (classical anchor demoted to simplicity-only); volume+cutoff-uniform interacting gap 3/8 on the factorized chain; governance capsule | PASS (pinned) |
| YM-13 verification-overlap beta audit: beta 4/5 -> 3/5 via independent witnesses (CF Bessel, compound-interest exp, det/trace inertia); three lineage overlaps retained with refusal notes; RH remains lowest at 1/5 | PASS (pinned) |
| YM-14 interaction-overlap dock, first carrier: bridged pair — overlap seam rank-one and swap-transported (same theta integral), certified k=1 to kappa=1/2, overlap price kappa*lam/4 on one line; kappa=1 carrier limit recorded | PASS (pinned) |
| YM-15 1D block-transfer dock opened: exact closed form of the chain compression for EVERY m (vacuum f0^{m-1}, A-line f0^{m-1}·KMS_m(r), r=I2/I1; SU(2) convolution lemma, machine-checked as a formal polynomial identity m=2..8, no truncation remainder); m-UNIFORM ratio bracket λ(1−r)/(1+r)<ρ_m<λ(1+r)/(1−r), certified positive normalised gap for all m at κ≤1; κ=2 fails closed. A-line carrier ONLY — complement not bounded | PASS (pinned) |
| YM-16 chain dock: exact B-factorisation ⇒ gap(S_m) ≤ Δ_red for ALL m,κ (m-uniform upper bound, exact eigenvector); vacuum volume bracket f0^{m−1} ≤ λ₁ ≤ e^{κ(m−1)}; certified two-sided gap (ratio ≤ 9/10) of S_m for m ≤ m*(κ) — m*=13 (κ=1/8), 7 (1/4), 4 (1/2) — with m* computed exactly and matched by the dock; m*+1 refused. Sup route provably non-uniform: death volume has coordinates | PASS (pinned) |
| YM-17 interleaving seam: m_κ = B_odd·B_even; each dressed half-chain is a tensor product of the two-site chain at 2κ ⇒ product spectrum ⇒ gap ratio EXACTLY m-uniform (pair ratio certified ≤ 9/10 at κ≤1/2; bracket [1,3] honest at κ=1); vacuum bracket narrowed to rate κ² (per-bridge upper end 0.125→0.030 at κ=1/8); whole volume problem reduced to ONE inequality: λ₁(T) vs σ₁(X)σ₁(Y) (vacuum tracking of interleaved halves) | PASS (pinned) |
| YM-18 CALIBRATION + DESIGN: enlarged exact carrier W'_m (vacuum + sites + bridge-pairs, dim 2m) gives certified Rayleigh lower bounds on λ₁(T_A,m), m=2..6; per-bridge vacuum rate φ_m^lo is m-independent to ~1e-5 (0.002177 @κ=1/8, 0.008698 @1/4, 0.03461 @1/2) — vacuum tracking behaves as a clean per-bridge rate; bracket width = looseness of the upper end. Cluster dock DESIGNED (monomers = odd pairs at 2κ, activities = even bridges, KP on a line) — no claim | PASS (calibration, pinned) |
| YM-19 Dobrushin dock — **historical capsule DEMOTED to ANCHORED; scoped replacement in [YM-43](YM43_NATIVE_TIME_TRANSFER_BOUND.md)** (governance rule: no classical borrowing; T3 rests on a cited Dobrushin/Föllmer anchor, same fate as YM-8 Jentzsch). T1 (coefficient bound) and T2 (α<1 arithmetic) are native and stand; the m-uniform GAP claim is WITHDRAWN until derived natively. Was: on the bounded-overlap chain — Dobrushin coefficient C_ij ≤ 1−1/δ² (elementary, proved), δ_s = e^{2κ} exact, δ_t ≤ (1+S_a)/(1−S_a) certified; α = 2(1−1/δ_t²)+2(1−1/δ_s²) < 1 ⇒ (pinned anchor: Dobrushin/Föllmer/Künsch covariance decay; extraction lemma proved) λ₂/λ₁ ≤ α for EVERY m. Certified gap/step ≥ 0.263 (a=6,κ=1/16), ≥ 0.145 (a=8,κ=1/8, ON YM-9's trajectory θ=1/64), ≥ 0.235 (a=12,κ=1/8); refused at (6,1/8),(4,1/16),(10,1/4). Coupling ceiling of the route exact: κ < log2/4. STANDING CORRECTION to YM-17 T3 (vacuum tracking insufficient; correlation-decay route replaces it). NOT uniform in the cutoff (δ_t→∞ as a→0) | ANCHORED (pinned, not a native theorem) |
| YM-20 NATIVE ORIGIN AUDIT: machine-checked origin ledger for YM-1..19 in the framework's vocabulary (carrier = COORDINATE_SHADOW throughout; counting = NATIVE_DERIVED throughout; anchors: YM-8 and YM-19 demoted imports, YM-17 SV-Weyl as technique). CAYLEY FORM certified: YM-15 ceiling (1+r)/(1−r) and YM-19 δ_t are Exp_Σ(A_Σ(y)) (F00-G Thm 5.1 pinned) with seam coordinates r, S_a; native odd-series log agrees with YM-1 route to 1e-20; uniform gap in native form = −Log_Σ(λ) − A_Σ(r). Bessel coefficients = native factorial series (F00-E Lemma 2.1 shape); Haar identification = the shadow. NAMED YM-21: bridge-local recognition-energy contraction under the seam-involution flow (F00-E Thm 6.2) | PASS (pinned) |
| **YM-F1 THE CHAIN AS A RECOGNITION FABRIC (native carrier)**: T1 SU(2) = unit sphere of the EMK block over C_Σ — span{I, R, ιK, ιRK} is Hamilton's quaternion algebra, det = Δ∥+Δ⊥ becomes the quaternion norm (twisted seam channels), χ½ = 2×identity-sector coefficient; T2 chain = ladder fabric, bridge faces = plaquette holonomies H_l = A_i⁻¹A_{i+1}; T3 Recognition-Stokes EXACT on the chain (interior rungs cancel, E_int = 0, orientation tamper bites); T4 action = sum of face residues ρ_l = 1−Tr H_l/2 ∈ [0,2], weight = Exp_Σ(κ(m−1))·Exp_Σ(−κΣρ). Source: vault fabric appendix + EMK-1 + F00-E. Remaining shadow: Haar ↔ Φ_Σ (declared). Native form of YM-21 posed on the fabric | PASS (pinned) |
| YM-21 THE TILING LAW (native, on the fabric): STANDING CORRECTION to YM-F1's posed question (boundary residue cannot control face residue — witness g, g⁻¹, certified); the right object is the tiling. T2 Stokes sub-telescoping ⇒ YM-15 entry = r^{#faces crossed}: SPATIAL decay along the chain exactly geometric, rate −Log r, every m. T3 leading-order TIME two-point = (λ·KMS_m(r))^t, bracketed λ^t((1∓r)/(1±r))^t for every m ⇒ YM-15/20's uniform rate −Log λ − A_Σ(r) DERIVED as a tiling count (t=1..6, m=2..8 inside bracket). T4 remainder named natively: higher-content/branching tilings, expansion parameter = face-coefficient ladder r_j (certified decreasing: r_½ ≫ r_1 ≫ r_{3/2}); uniform convergence = native cluster statement | PASS (pinned) |
| **YM-22 THE TILING RATE SURVIVES THE CUTOFF**: on the heat-kernel family along κ=θa, the time-face log is EXACTLY 3a/4 (no transcendental), so the leading tiling rate per unit time γ(a) = 3/4 − A_Σ(r(θa))/a. Certified over SIX DECADES a=1…1e-6 at θ=1/16: **23/32 ≤ γ(a) ≤ 3/4 for every a** — the limit 3/4 − θ/2 is a uniform LOWER bound in the cutoff (approached from above); KMS bracket holds with λ_a for m=2..8 at every a ⇒ leading-order rate uniform in volume AND cutoff simultaneously. Beats YM-9's sandwich 3/4 − 6θ by factor 12 in the θ-coefficient; trajectory ceiling θ < 3/2 (θ=2 refused). Higher face coefficients per unit time r_1/a, r_{3/2}/a → 0: in the cutoff limit only j=½ branching survives. SCOPE: leading tiling order | PASS (pinned) |
| **YM-23 THE WEAK-COUPLING TURN**: T1 tiling route has an EXACT weak-coupling boundary — ladder collapses (r_½,r_1,r_{3/2} = 0.970, 0.922, 0.859 at κ=50: no small parameter), leading rate sign change at r_c(a) = tanh(3a/8) in native Cayley form, κ_c(a) bracketed (1.574 @a=1, 0.759 @1/2, 0.376 @1/4, 0.188 @1/8), and on a DECLARED AF-shaped trajectory κ = ¼log(1/a) the route is left at a = 1/8 — the continuum lies outside the tiling route. T2 native linearization: ρ = 2s/(1+s) in the stereographic chart; odd(H₁H₂) = p₁+p₂+½[p₁,p₂] EXACTLY — the commutator is the only non-additive term = EMK-1 T5 / RST-1 T4 curvature on the fabric. T3 EXACT NEGATIVE: the abelianized fabric's softest spatial stiffness ≤ 12κ/(m(m+1)) (Rayleigh identity certified m=2..40) — gapless in volume. T4 native statement: a weak-coupling gap lives entirely in the commutator sector (non-healable curvature) — the mass gap is non-perturbative, said in the framework's words | PASS (pinned) |
| YM-24 ABELIAN SUBFABRIC + NON-ABELIAN RESIDUE (exact algebra, no gap claim): T1 commuting faces = F00-E's circular Euler orbit, compose by the native addition law, Stokes additive in the Euler coordinate (the rotor chain), m=2..8 on rational circle points; T2 non-abelian residue = Gram defect |p×q|² = |p|²|q|² − ⟨p,q⟩² (zero iff one subfabric); T3 exact second-order Stokes odd₂ = Σp_i + Σ_{i<j} p_i×p_j, even₂ = 1 − Σ⟨p_i,p_j⟩, remainder degree ≥ 3 by exact interpolation; T4 YM-23's soft mode lies in an EXACT abelian subfabric, E_na = |Σ_{i<j}p_i×p_j|² vanishes there, positive generically, rotation-invariant; order tamper = 2 p_i×p_j (EMK-T2: order is content). Consequence: a weak-coupling gap must come from transitions BETWEEN abelian subfabrics (direction changes) | PASS (pinned) |
| YM-25 THE DIRECTION-CHANGE SECTOR: two exact negatives — the time kernel is direction-blind (bi-invariant: one scalar per spin-j block, d_j² coefficients, refinement-invariant) and the strong-coupling carriers are direction-blind (class functions, rotation-invariant: structural SPECTRAL-1 blindness) — and one exact positive: direction change between adjacent faces IS the 6j recoupling (MP gold/01 seed), squared overlaps between (12)3 and 1(23) schemes exactly [[1/4,3/4],[3/4,1/4]], unitary, no √ evaluated (RST-2 discipline), singlet tamper breaks unitarity. REDUCTION: a weak-coupling gap can live only in the intertwiner transfer | PASS (pinned) |
| **YM-26 DOCK TO theorum/28 (Recognition-Complete Finite-to-Infinite Cut Theorem)** — owner's correction: the framework's own convergence machinery. Dictionary: P_n = A-line carrier (YM-15), Q_n = complement, floor λ₁ ≥ f₀^{m−1}, threshold μ = νf₀^{m−1}, e_n = ‖QSQ‖/μ, N₊(B_n−I) = Haynsworth inertia; outward certificate e_n < 1 ∧ inertia count = 1 ⇒ λ₂ < νλ₁. T1: LAWFUL (e_n present, not a shadow certificate) for every m ≤ m*(κ) = 13, 7, 4 at κ = 1/8, 1/4, 1/2. T2: e_n(m) certified increasing, crosses 1 exactly at m*+1 — theorum/28's hypothesis 3 (memory-channel Cauchy bound uniform in m) is NOT YET BUILT; full hypothesis ledger 1–7 written in the theorem's vocabulary. T3: via the tiling law, item 3 reduces to a Cauchy bound uniform in m for the multi-insertion content-½ memory channel (the j=½ branching sector) | PASS (pinned) |
| **YM-27 DOCK TO RH-FRAMEWORK T01** (owner: use the RH journey's proved results, cite them): T1 the sup route IS T01-E4C (bounded native simple multiplier, U_Σ(w) = e^κ exactly; YM-16's e_n reproduced to the last digit) ⇒ theorum/28 hypothesis 3 = a memory-channel bound beating T01-E4C by a factor uniform in m — measured slack per face 56.5 / 27.8 / 13.5 at κ = 1/8, 1/4, 1/2 (YM-18 calibration vs E4C price). T2 the abelian subfabric is a T01-E5A rational UGD packet: wrap law = rotor-chain composition, winding = seam memory μ; (g,g⁻¹) witness vs full-turn fabric have the same visible boundary but μ = 0 vs 1 — the native form of EMK-1 winding. T3 the Haar shadow of YM-F1 is exactly RH's OPEN T01-E5C/E6 — YM inherits RH's item, invents none; aliasing structure of class functions on E1 grids certified (exact for n > 2j). T4 T01-B/C recognition energy = chain vacuum form, floor f₀^{m−1} | PASS (pinned) |
| **YM-28 T01-E4D ON THE CHAIN** (RH-Framework PR #19 applied): T1 the mean law is EXACT with zero deviation on single insertions — I_rec(W·χ½(A_i)) = f₀(2κ)^{m−1} for every m (formal identity, χ₁ branch dies at the open end); T2 two insertions: deviation = r₁(2κ)^{|i−j|}, a geometric correlator; T3 |S| ≤ 4, m ≤ 7: deviation ≤ (1+2r₁)^{|S|−1} — PER INSERTION, not per face (m does not enter); T4 E4D per-face price log(√f₀(2κ)/f₀(κ)) = 0.001949 / 0.007742 / 0.03016 accounts for 90 / 89 / 87 % of YM-18's measured true rate — the slack of YM-27 is now DERIVED (E4C/E4D = 63 / 31 / 16); T5 conditional projection m*_E4D ≥ 400 / 200 / 50 vs sup-route 13 / 7 / 4 IF E4D-C held on superpositions | PASS (pinned) |
| YM-29 CLOSED-FORM VACUUM FLOOR FOR EVERY m: native face-independence lemma (bridge faces independent contents, W a product over faces ⇒ face-product integrals factorise; fusion-rule certified: ∫wχ½ = 2f_½, ∫wχ½² = f₀+3f₁); T0^{1/2}χ½(H_ℓ) = λχ½(H_ℓ); ansatz ψ = 1 + cΣχ½(H_ℓ) gives λ₁ ≥ f₀^{m−1}·Q_m(c*) in closed form for ALL m, Q_m > 1 always, ~4λ²r²(m−1) for m ≫ 1/(4λ²r²) (Q₂₀₀₀ = 2.58 / 6.99 / 24.1). Honest: linear floor moves the dock's m* by 0 — the floor at the TRUE rate (0.000228/face margin, YM-28) needs the product ansatz, whose T-expectation is the one-time-step ladder = 6j recoupling contraction | PASS (pinned) |
| **YM-30 THE 2-ROW RECOUPLING ENGINE — BUILT AND VALIDATED**: exact surd arithmetic (Q adjoined √), Clebsch–Gordan by Racah's closed form (no floats), 3j, three-matrix-coefficient vertex integral, ε-removal of conjugates, ladder (one-time-step fabric) evaluated as a transfer over cut indices with label sums folded in. V1 every evaluation rational (surds cancel — the built-in consistency check); V2 one plaquette = 1/d_a²; V3 single-row limit recovers the YM-15 chain; V4 two plaquettes: Σ_c d_c N = 1/4 and the split d₀N₀ : d₁N₁ = 1/4 : 3/4 = YM-25's recoupling squares, now from the engine; V5 rail-swap and reversal symmetry. FIRST USE (calibration, κ=1/8, content ≤ 1, point coefficients): the W^{1/2}-dressed vacuum for T'' = W^{1/2}T0W^{1/2} beats the old floor f₀^{m−1} by a GEOMETRIC factor 1.000732 per face (m=2..6) — floor rate 0.002685/face > YM-18's bound 0.002177 > E4D excitation rate 0.001949: margin 0.000736/face; single-face excitation (c ≠ 0) does not help | PASS (pinned) |
| **YM-31 OUTWARD m-UNIFORM DRESSED VACUUM FLOOR — theorum/28 ITEM 6 AT THE TRUE RATE**: trial φ = W^{1/2}·(W_pt/W) with W_pt a product of rational truncated faces (no error by definition); numerator = ladder ⟨W_pt, T0 W_pt⟩ (YM-30 engine), rung truncation is a rigorous LOWER bound (T0 = Π_i Σ_c λ_c Π_c: a sum of nonnegative projection terms, monotone in λ_c); denominator ∫W_pt²/W = [Σ d_jd_kd_l f_j^pt f_k^pt (−1)^{2l} f_l N_jkl]^{m−1} exact with outward l-tail; NATIVE LOG-CONVEXITY Z_m² ≤ Z_{m−1}Z_{m+1} (Z_m = ⟨g, S^{m−2}g⟩, S = C^{1/2}M_K C^{1/2} symmetric PSD, cut-square inequality) ⇒ ratio nondecreasing ⇒ **λ₁(S_m) ≥ FLOOR_{m₀}·(ρ_lo/I_hi)^{m−m₀} for every m ≥ 6**, rung tail at m₀ bounded (5.9e-9). Per-face floor rate **0.002685 / 0.010722 / 0.042607** at κ = 1/8, 1/4, 1/2 vs old floor 0.001952 / 0.007802 / 0.031089 and E4D excitation rate 0.001949 / 0.007742 / 0.030161 — margin **0.000736 / 0.002980 / 0.012446 per face**, uniformly in m | PASS (pinned) |
| YM-32 DRESSED SINGLE-INSERTION ENERGY, m-UNIFORM (engine extended to rung insertions via CG reduction of 4-valent vertices; validated: j=0 reproduces YM-30, one-rail insertion vanishes, reversal identity, mixed-rung tamper smaller): E_m(p) = ⟨W_ptχ½(A_p), T0 W_ptχ½(A_p)⟩/⟨W_pt,T0W_pt⟩ exact for m=2..7, all p — bulk value m-independent to < 1e-6, ends differ; E_bulk = 0.43368 / 0.43532 / 0.44171 at κ = 1/8, 1/4, 1/2, i.e. λ·(1 + 1.27e-3 / 5.05e-3 / 1.98e-2), inside YM-15's bracket (I first wrote 'below λ'; the engine corrected me). Dressed denominator factorises (rung insertion free) ⇒ the J=½ sector's dressed excitation/vacuum ratio is exactly E_bulk, uniformly in m — a certified LOWER bound on that sector's top, not an upper bound on λ₂ | PASS (pinned) |
| YM-33 t-ROW FABRIC ENGINE + EXACT TIME TWO-POINT FUNCTION (owner's framing: time carries no infinity, only the flow λ^τ = Exp_Σ(−τ·gap); expand the space coupling only): general cut transfer over t rails with rungs contracted up the column (interior vertices 5-leg: two rungs, two faces, insertion — two CG fusions); validated: t=2 reproduces YM-32 exactly, t=1 the chain, 3-row fabrics rational with exact time reflection. Dressed sequence a_τ = ⟨φ, T''_pt^τ φ⟩ and vacuum z_τ exact for τ=1,2, m=5..8: ratio sequence ρ₁ = E_bulk (YM-32), ρ₂ ≥ ρ₁, and the SECOND-ORDER connected correction δ₂ = ρ₂ − ρ₁ is a BULK CONSTANT — 0.00076998 / 0.00305278 / 0.01178845 at κ = 1/8, 1/4, 1/2 — with boundary influence decaying geometrically: moving an end from distance 2→3 shifts δ₂ by 1.3e-9 (either end, identically), 3→4 by 1.6e-12, ratio certified in [r²/2, 2r²]. The framework's own strong-coupling series for the gap, order by order as exact rationals, m-uniform with exponentially localised edges. Lower-bound side (Rayleigh data) — item 3's upper bound still open | PASS (pinned) |
| YM-34 E4D-C RESTATED (native Laplace principle): the iterated multiplier WITHOUT smoothing, q_n = ∫wⁿ/∫wⁿ⁻¹, is strictly increasing to sup w = e^{2κ} (certified n ≤ 60) — so E4D-C as written in the RH ledger (deviation control of the multiplier iterated on one state) is unprovable, and the exponential of YM-16/26 is this growth. WITH the time kernel between multiplications the ratio is nondecreasing (log-convexity) and converges to λ₁(K^{1/2}WK^{1/2}) = f₀(1+δ), δ = 0.003582 / 0.014032 / 0.051911 (κ=1/8,1/4,1/2) vs sup price 0.274 / 0.598 / 1.405. E4D-C is restated for the alternating product (K^{1/2}M_wK^{1/2})ⁿ: λ₂/λ₁ ≤ 1−δ_gap uniform in m — the vacuum side of this is YM-31; only the excitation side is open | PASS (pinned) |
| Remaining: theorum/28 item 3 — the memory channel's OPERATOR norm (not basis-state norms) uniform in m, via the same engine (E4D-C quadratic form on content-½); with YM-31's floor the dock's e_n becomes (‖QTQ‖/FLOOR) and the sup numerator is the last exponential; (interval coefficients, content tail, Perron bound on the column transfer → m-uniform floor at rate ≥ 0.00268) and the E4D-C quadratic form on content-½ via the same engine; (exact rational ladder evaluations) — closes item 6 at the true rate AND item 3 (E4D-C quadratic form on content-½) in one build; E4D-C on SUPERPOSITIONS of content-½ basis states (a quadratic-form statement on the multi-insertion channel) = theorum/28 item 3 (content-½ multi-insertion memory channel, Cauchy-uniform in m); (weak) YM-27 the intertwiner (spin-network) transfer on the chain fabric — recoupling matrices as weak-coupling tiling weights, floor uniform in m? (declared: chain ~ O(4)-type rotor model in 1+1D, continuum gap at the frontier for every method); (strong) FULL-ORDER tiling convergence uniform in m and a (j=½ branching); (weak) the COMMUTATOR-RESIDUE TRANSFER on the fabric — does the non-abelian sector alone carry a gap; AF trajectory derivation; (convergence of the face-coefficient tiling expansion) — the cluster statement in fabric language; then (YM-43 replaces the cited comparison at three coarse heat-kernel cells; YM-44 constructs full time refinement at fixed width; YM-45 proves a width/fine-time uniform gap for abs(theta)<1/1680 on declared kappa=theta a using temporal blocks); YM-46 constructs the infinite-volume local state and normalized time map at each fixed a in that window; YM-47 joins the local vacuum history limits at arbitrary relative cutoff rates; YM-48 constructs the reflected continuous-time generator and a nonzero coefficient source; uniformity outside the window and full instantaneous-row/observable identification remain; AF trajectory, tightness, OS, non-triviality, metric universality, 2D lattice, Clay | OPEN |
| (superseded as the main route) cluster dock YM-18 proper: certified Kotecký–Preiss radius for the interleaving polymer gas with two-site pairs as monomers; per-bridge/local control in the content basis — certified strong-coupling cluster radius), 2D lattice, AF trajectory, tightness, OS reconstruction, non-triviality, metric universality, Clay predicate | OPEN |
| Continuum existence / mass gap | OPEN |

## Reproduce

New in YM-2: for beta=2 the interacting theta-graph transfer
`T_kappa = M^{1/2}(K x K)M^{1/2}` has a certified positive spectral gap for
every coupling `0 <= kappa < kappa_0 = Delta_red/6 = 0.13945388...`, via the
min-max sandwich `m_- T_0 <= T_kappa <= m_+ T_0` with
`m_+/m_- = e^{6 kappa}` from `|Tr U| <= 2`. At `kappa = 1/8` the reduced gap
is certified `>= 0.08672330623...`. Fail-closed above threshold.

```bash
python papers/yang-mills-certified-benchmark/certificates/ym1_certified_gap.py
python papers/yang-mills-certified-benchmark/certificates/ym2_theta_interacting_gap.py
python -m pytest papers/yang-mills-certified-benchmark/tests -v
```

New in YM-3 (exact, zero coupling): the free top-excited eigenspace is exactly
two-dimensional; the exact theta integral `Int chi12(A)chi12(B)chi12(AB^-1) = 1/2`
splits it at first order with derivatives `+-lambda_half/4`; the vacuum derivative
is exactly 0. Hence `r'(0) = lambda_half/4` (certified `0.10828185668...`), the
first-order CROSSING DIRECTION is the single symmetric line
`chi12(A)+chi12(B)` — rank-one and transported by the graph symmetry — and
YM-2's sandwich slope is slack by the exact factor 24.

Single-command reproduction (regenerates and pins all three):

```bash
python papers/yang-mills-certified-benchmark/certificates/ym12_square_sourced.py
```

YM-12 transplants the RH journey's decisive step (T03: re-source
positivity as a square with declared defects) to Yang-Mills:

- **T1** manifest square: `T_kappa(a) = S*S` with
  `S = [K_{a/2} x K_{a/2}] m_{kappa/2}` — both factors are the program's
  own halved-parameter objects (YM-9 semigroup + YM-5 doubling).
  Positivity is now **square-sourced**; YM-8's Jentzsch anchor is demoted
  to top-eigenvalue *simplicity* only. The compressed defect
  `A - N*G^{-1}N` is exactly a Gram of `(I-P)S phi` — a source-restriction
  object, the T03 shape.
- **T2** on the theta chain `L_m` (vertex-glued, no shared faces) the
  transfer tensor-factorizes exactly, so `lambda_2/lambda_1` is
  `m`-independent, and with YM-9: `Delta(L_m, a, a/16) >= 3/8` for
  **every** `m` and **every** `a` — an interacting gap uniform in volume
  and cutoff simultaneously. This closes YM-11 gate 2's interacting half
  on the factorized family and localizes the remaining volume obstruction
  to **interaction overlap** (shared faces): the chain (no sharing) is
  uniform, `B_n` (complete sharing) kills the sup route, and a physical
  lattice's bounded sharing sits strictly between the two certified
  extremes.
- **T3** governance capsule (RH 07-20 analog): claim/consumption chain for
  YM-1..12, standing corrections (YM-9 T4 amendment, YM-8 anchor
  demotion), and a do-not-reopen list — machine-checked for completeness.

YM-11 takes three dependency-map gates and gives each an exact verdict on
the toy carrier — the generalized theta ("n-banana") graphs `B_n`, two
vertices joined by `n` edges, so `b_1 = n-1` holonomies and
`F(n) = n(n-1)/2` plaquettes; `B_3` is the theta graph of YM-1..10.
**These are toy-carrier verdicts, not Clay-sense closures**, and two of
the three are partly negative.

- **Gauge — CLOSED.** Tree gauge-fixing on a connected graph is exact (no
  Gribov obstruction), leaving `b_1 = E-V+1` holonomies and a residual
  diagonal `G` acting by conjugation, so gauge-invariant states are
  exactly `L^2(G^{b_1})^Ad`. For `B_3` that **is** the carrier every
  capsule has used — now certified rather than assumed, with Wilson loops
  spanning it and YM-10's multiplicity law counting them.
- **Volume / IR — SPLIT.** The *free* reduced gap is `C_(1/2) = 3/4` for
  every `n` and every `a`: volume- **and** cutoff-uniform, closed. The
  *interacting* sandwich gives `Delta >= C_(1/2) - 2 F(n) theta`, which
  degrades quadratically in `n`: at `theta = 1/16` the bounds are `3/8`
  (n=3), exactly `0` (n=4), `-1/2` (n=5). The route is **proven
  insufficient** for volume-uniformity, with the critical volume computed
  exactly. What a replacement must do is named: control the interaction
  per-plaquette (dock route, YM-6) rather than by a global sup (envelope
  route, YM-5) — the same lesson one level up.
- **Universality — SPLIT.** For the whole regulator class
  `K = sum_j d_j c_j chi_j` with `c_j > 0`, the **counting layer**
  (content grading, multiplicity law, carrier dimensions, blindness
  ledger, probe counts) depends only on SU(2) representation theory and
  is therefore regulator-independent — closed, and it means YM-10's
  ledger is kinematics, not a choice of action. The **metric layer** (the
  gap value) genuinely differs between regulators (witnessed) and remains
  the real universality gate.

YM-10 applies the source-restriction framework of the companion
manuscripts *When Spectra Forget Order* and *Stable Recovery Beyond
Spectral Blindness* to our own compressions. Every dock computes on a
compressed carrier `V`; in that language `V` is a source restriction, so
the honest question is not whether the complement bound is tight but
whether the carrier is **blind** to the target.

- **T1** exact multiplicity law `m(j1,j2) = min(2j1,2j2)+1` on the
  Ad-invariant sector — it reproduces `dim V5 = 5` and `dim V7 = 7`
  independently, and re-derives YM-6's `(1/2,1/2)` two-dimensionality from
  representation theory rather than from the Gram matrix.
- **T2 correction to YM-9 T4** (self-audit): YM-9's uniform count counted
  *contents*; `k_Sigma` counts *eigenvalues with multiplicity*. At `s = 2`
  that is 4 versus 5. Both are `a`-independent, so YM-9's uniformity
  conclusion is unaffected — only the arithmetic. YM-9 is amended to
  report both.
- **T3/T4** blindness ledger and exact probe count
  `rank(Pi_s | ker E) = b(V,s)`. Result worth noting: `b(V5, s) = 0` for
  `s <= 2` — the YM-6 carrier is **exactly blindness-free in its own
  threshold window**, so that choice was correct, not lucky. At `s = 3`,
  V5 needs 6 probes, V7 needs 4, V9 needs 0 more than V7 minus 4.
- **T5** composition ledger: blindness is nonincreasing along
  `V5 ⊂ V7 ⊂ V9` with nonnegative stage increments.
- **T6** the whole ledger depends only on contents and Casimir levels, so
  it is exactly `a`-independent — one ledger covers the entire YM-9
  refinement family.

Honest remainder: the ledger is exact for the **free** target; the
interacting faces mix the content grading, so `b(V,s)` is a lower bound
there and YM-6/7's declared truncation remainder is the extra obligation.

YM-9 — the program's first statement that is uniform in the cutoff rather
than at a fixed cutoff. Switching the link action from Wilson to the
HEAT KERNEL `K_a(g) = sum_j d_j e^{-a C_j} chi_j(g)` makes refinement
EXACT (`K_a * K_b = K_{a+b}`; n links of spacing a/n compose to one link
of spacing a with no discretization error), and the cutoff dependence then
cancels identically:

- **T2** free reduced gap `= C_(1/2) = 3/4` EXACTLY for every `a > 0`;
- **T3** along the declared trajectory `kappa(a) = theta*a`,
  `Delta(a, kappa(a)) >= 3/4 - 6*theta` for EVERY `a > 0` — at
  `theta = 1/16` a cutoff-independent gap of exactly **3/8**; fail-closed
  at `theta >= 1/8`;
- **T4** the seam count `k(a, e^{-as})` is exactly `a`-independent
  (combinatorial), so the dock's threshold grammar transports across the
  whole family at once;
- **T5** Wilson contrast, computed: at fixed `beta` the reduced gap grows
  like `1/a` and diverges — no cutoff-independent gap without a running
  coupling.

Honest remainder (in-certificate): the trajectory is DECLARED not derived
(the physical one is fixed by asymptotic freedom, `beta ~ log(1/a)`, not
modelled); the graph is FIXED (no infinite-volume growth); the heat-kernel
action is a CHOICE (universality gate open); `Delta` is a reduced graph
gap, not a reconstructed physical mass gap. Net movement: one gate from
"untouched" to "touched on a toy carrier".

Capstone (YM-8): the interacting theta-graph transfer has a strictly
positive spectral gap at EVERY coupling — certified pointwise kernel
floors (`k >= e^{-3 kappa} [c0^{-1} e^{-beta}]^2 > 0`, exact enclosures)
plus the pinned classical Jentzsch/Krein-Rutman anchor (simple Perron
eigenvalue; CIRC-1 discipline: cited, never rederived). Quantified window
`kappa in [0, 7/10]` carries exact counts and enclosures (YM-2..7). The
honest remainder is named exactly: nothing here is uniform in the lattice;
a gap at every fixed cutoff is compatible with the gap closing in the
continuum limit — the fixed-regulator lesson. The Millennium content is
precisely the open uniformity/OS/vacuum gates of the MP dependency map.

New in YM-7: carrier enlarged to V7 (adds `chi_1(A), chi_1(B)`; Gram stays
block-diagonal, all seven exact T0 eigenvectors), dropping the complement
top from `lambda_1` to `lambda_1*lambda_half` (certified ordering) — which
UNLOCKS the previously refused `kappa = 7/10` cell: exact seam count
`k_Sigma = 1` at `(7/10, 13/10)`. Also ships certified brackets for all
seven compressed eigenvalue curves (inertia bisection; every step an exact
LDL count), published with honest kappa-drift — adopting the RH-line's
post-pause lesson that limited-denominator "exact constants" are never
promoted from bisection.

New in YM-6 (framework-native): the gap becomes an exact threshold COUNT
`k_Sigma(kappa,mu) = N(spec(T_kappa) > mu)` via a Haynsworth congruence —
the D03 "infinite ko finite" move — on the native 5x5 carrier
`{1, chi12(A), chi12(B), chi12(A)chi12(B), chi12(AB^-1)}` (all exact T0
eigenvectors; exact Gram has a single 1/2 overlap). Certified exact counts:
`k=1` at (kappa,mu) = (1/8,3/5), (1/4,3/5), (1/2,1) — so
`lambda_2 <= mu < lambda_1` EXACTLY, past YM-5's reach — and constant
count across the window = no seam spectral flow (thermodynamics/02
instance). The (7/10,5/4) cell is honestly REFUSED (bracket [1,5]).

New in YM-5: certified TWO-SIDED spectral gap of the interacting theta
transfer at finite coupling — including BEYOND YM-2's sandwich threshold
`kappa_0 = 0.1394...`. Complement control uses the exact doubling identity
`m_kappa^2 = m_{2 kappa}`: `|PMQ|^2 <= lammax(P m_{2k} P - (P m_k P)^2)`,
with `|QSQ| <= lambda_half^2 e^{3k}` and a Weyl bound. Certified rows
(ratio upper bound, gap lower bound): kappa=1/8: 0.5003, gap>=0.6924;
kappa=1/4: 0.5610, gap>=0.5781; kappa=3/10: 0.5867, gap>=0.5332.
Fail-closed at kappa=1. Lineage: realizes the old lambda-manuscripts'
Kill-Lemma skeleton unconditionally (see LINEAGE.md); no claim imported.

New in YM-4: (T1) the graph swap S:(A,B)->(B,A) commutes with T_kappa for
EVERY coupling — `Tr(BA^-1) = Tr(AB^-1)` on SU(2) — so YM-3's rank-one
crossing line `chi12(A)+chi12(B)` is symmetry-protected at all kappa: a
theorem, not a first-order accident. (T2) first certified finite-kappa
spectral data from the exact character expansion
`exp[(k/2)TrU] = sum d_j f_j chi_j`, `f_j = 2 I_{2j+1}(k)/k`, with exact
Clebsch-Gordan ring arithmetic, the exact pairing tensor
`Int chi_p(A)chi_q(B)chi_r(AB^-1) = delta_{p=q=r}/d_p`, and a certified
truncation remainder. Sector-resolved lower bounds on the grid show the
swap-even branch rising and the swap-odd branch falling with coupling.

CI runs exactly that one command plus pin diff and pytest.

## YM-35 — T54′ energy version + commutator remainder (E4D-C with the RKF operator ladder)

Framework inputs: theorum/54 (constant-μ chain is exactly its uncovered case — `sum_mu_diverges`), theorum/50 A4 (polarization visibility = YM-28's superposition gap), theorum/53 (energy, not mass), theorum/51 + 61 (commutator calculus). On the fabric, with the YM-30 engine as the only evaluator:

- **T1 locality** — `[K^{1/2}, B_l]` vanishes on states with no content at sites l, l+1; non-adjacent face pairs do not mix (site content exactly ½). Commutator support is `{l−1, l+1}` — bridge-local, every m.
- **T2 exact adjacent commutator energy** — `E_Σ([K^{1/2},B_l] B_{l+1}) = λ²·(3/4)·(1−√λ₁)²`, κ- and m-independent (site weights 1/4 : 3/4 = YM-25 recoupling squares, from the engine).
- **T3 T54′ ratio theorem (commuting model, exact m ≤ 8)** — vacuum and excitation carry the same per-face factor; ratio exactly m-independent; `Σμ<∞` not needed for a ratio.
- **T4 second-order remainder** — relative size ρ = μ²·(3/4)(1−√λ₁)² = 2.9e-4 / 1.2e-3 / 4.6e-3 at κ = 1/8, 1/4, 1/2 (a = 1); kill criterion ρ < 1/10 passed.

**Verdict:** E4D-C's excitation bound opens at the first non-commutative order — the first obstruction is an exact local number, not a sup. **Not done:** the full face expansion (all words, contents j ≥ 1) — the worldline cluster problem (YM-33). E4D-C OPEN. Pin `EXPECTED_YM35.sha256`.

## YM-36 — T54′ on the fabric at content ½ (face-word carrier)

The alternating product `K^{1/2} M_w K^{1/2}` compressed to ALL superpositions of content-½ face insertions (2^{m−1} orthonormal words), exact rational matrices m = 2..6 at μ = 1/32, 1/16, 1/8 (declared rational-square kernel √λ½ = 11/16, √λ₁ = 3/8; YM-33 engine with per-column face labels). Excitation/vacuum ratio ρ_m increases with geometrically decaying increments (consecutive-increment ratios in [2/5, 3/5] certified); under declared q ≤ 2/3 the m-uniform limit is bracketed, e.g. ρ∞ ∈ [0.2225014, 0.2225544] at μ = 1/32. The top eigenvector is a genuine superposition (single-insertion weight 0.98 / 0.92 / 0.73), yet exceeds the one-face ratio by only 0.24 % / 0.93 % / 3.5 % — superpositions open no gap-closing channel at content ½. **Boundary:** compressed-operator statement (interlacing: lower bounds), contents j ≥ 1 outside the carrier (first effect = YM-35 commutator, ~μ²), m ≤ 6. E4D-C OPEN; it holds at content ½ on this carrier in compression with a geometric tail. Pin `EXPECTED_YM36.sha256`.

## YM-37 — The space transfer: m-uniformity of the dressed time two-point ratio as a theorem

**(a) refused with witness.** theorum/53 T4 bounds cut-square *weights* (`d_k ≤ S_kk`); a block with every diagonal ≤ μ can still have λmax > μ (2×2 witness, exact). So "complement λmax from T53 weights" is not available; item 3's upper half stays OPEN.

**What is a theorem.** Every YM-33/36 quantity is `⟨boundary, M τ^{m−1} boundary⟩` for one column map τ = B⁻¹M (native quaternion carrier: faces/rungs are class functions of products, so only dot products appear; Haar = rational S³ moments; cut space = harmonic degree ≤ 2 per rail). τ commutes with simultaneous conjugation and with the centre of SU(2) acting on all rails; both boundary vectors are even, so everything lives in the even invariant sector: dim 8 (t = 2), 52 (t = 3). Validated **exactly** against `ym33.fabric_partition`. σ₁, σ₂ of the symmetric pencil (M, B) bracketed by exact inertia (native T53 symmetric elimination; the earlier Bareiss narration was incorrect); q := σ₂/σ₁ certified:

| κ | q (t = 2) | q/r² | q (t = 3) | q/r² | ρ₁^∞ (two-sided) | certified err, interior p ≥ 3 |
|---|---|---|---|---|---|---|
| 1/8 | 0.0014578 | 1.495 | 0.0021668 | 2.222 | 0.43367710465900 | 9.3e-7 |
| 1/4 | 0.0057991 | 1.492 | 0.0086057 | 2.215 | [0.4353154056078, 0.4353154056137] | 1.5e-5 |
| 1/2 | 0.0226981 | 1.483 | 0.0334686 | 2.187 | [0.4417107170, 0.4417107229] | 2.6e-4 |

so `|ρ₁(m,p) − ρ₁^∞| ≤ err(m,p)` for **every** m and interior p, with err explicit in σ₁, σ₂, ‖N‖_B, ‖e₀‖_B — YM-33's "bulk constant with geometric edges" is now certified with its rate (both rates have the r² scale; the t=3 values lie above the earlier 2r² guessed ceiling), and ρ₁^∞ = N₁₁/σ₁ is a single-column Rayleigh quotient (YM-32's E_bulk reproduced to every digit). **Boundary:** fixed time order only; no complement bound; no m-uniform upper bound on the true gap — item 3's m-dependence is now a one-column pencil question, its upper half unchanged. Runtime ~16 min (own CI step). Pin `EXPECTED_YM37.sha256`. **Standing correction (same day):** T3's for-all-m route uses the finite spectral split — classical import; status ANCHORED pending YM-38. **Resolved by YM-38 (same day): the for-all-m statement now rests on the native route; see the YM-38 entry.**

## YM-38 — Native m-uniformity (YM-37's spectral import removed)

The YM-37 T3 statement re-proved with the framework's own machinery, nothing else: deflation with the explicit vector w = x₁₂ (defect exactly B-orthogonal, `r_w'Bw = 0`), complement ceiling `‖Πτy‖_B ≤ σ''‖y‖_B` from **one T53 elimination sign check** on the restricted doubled pencil (no eigenvectors anywhere), recognized-channel Cauchy with declared geometric tail and outward certificate **β = L < 1** (theorum/28 §3–5, §9), window products handled by theorum/54 §2's `∏(1+x) ≤ 1+2s` bound. Result per κ, two regimes:

`|ρ₁(m,p) − ρ_c| ≤ A·(L^{j−p₀} + L^{j_c−p₀})` for **all** m, p with `j = min(p−1, m−p) ≥ p₀`, where ρ_c = exact ρ₁(26,13):

| κ | L (= β) | A (p₀ = 2, all interior p) | check rows |
|---|---|---|---|
| 1/8 | 0.0014578 | 1.93e-8 | 8/8 inside |
| 1/4 | 0.0057991 | 6.34e-7 | 8/8 inside |
| 1/2 | 0.0226981 | 2.17e-5 | 8/8 inside |

The coarse constant beats the anchored YM-37 one (1.9e-8 vs 9.3e-7 at κ = 1/8); ρ_c lies inside YM-37's bracket at every κ. **YM-37's standing correction is RESOLVED** — its T3 now rests on this capsule; also corrected there: the pinned dim-52 run used the native elimination, not Bareiss (LINEAGE). Certificates refused three drafts (missing reference-tail term; backwards extrapolation below p₀; display-precision comparison) — all recorded. Pin `EXPECTED_YM38.sha256`.

## YM-39 — Content ladder: centre superselection (exact) + certified channel ordering

Session item (b). **T1:** the centre of SU(2) (total flip on every rail) grades the fabric exactly — faces and rungs are centre-blind class functions of products, χ_j picks (−1)^{2j} — so every mixed half-integer/integer fabric quantity vanishes **identically**: the mixed insertion matrix is the zero matrix (exact rationals) and the YM-33 engine's mixed partitions are 0; a same-parity pair (½ with 3/2) is nonzero, so the vanishing is parity, not accident. Superselection at every coupling, every m, every time order. **T2/T3:** per-channel bulk ratios by the YM-38 native route (same τ, same L, channel-specific N), with the m-uniform band on both sides of every comparison:

| κ | ρ(½) | ρ(1) | ρ(3/2) | ordering with margins |
|---|---|---|---|---|
| 1/8 | 0.43367710 | 0.13403748 | 0.00016958 | ✓ |
| 1/4 | 0.43531541 | 0.13491222 | 0.00067672 | ✓ |
| 1/2 | 0.44171072 | 0.13837705 | 0.00268149 | ✓ |

ρ(1) sits at the declared λ₁ = 9/64 scale (anchor check). So on this carrier the J = ½ excitation is certified **lightest, uniformly in m**, and no content ≥ 1 channel can mix into it — YM-36's boundary concern closed at the carrier level. **Boundary:** truncated faces/rung kernel (contents ≤ 1), so 3/2 is a composite row (data); compressed statements; E4D-C OPEN. Named next: the same ladder on the YM-22 heat-kernel trajectory, quantifying the collapse onto content ½ as a → 0. Pin `EXPECTED_YM39.sha256`.

## YM-40 — Content ladder on the heat-kernel trajectory: the Casimir pinch

YM-39's ladder run along YM-22's family (heat-kernel time kernel, Wilson faces at κ = a/16), a ∈ {1, 1/4, 1/16, 1/64}, everything per-a by the YM-37/38/39 native machinery (superselection exact at every a; every γ bracket carries the m-uniform band, so it holds for every m and interior p; native exp = YM-13 compound interest, native log = F00-G odd series, with a cross-route re-enclosure control):

| a | γ_½ | γ₁ | Δγ | ε = \|Δγ − 5/4\| | r₁/r_½ |
|---|---|---|---|---|---|
| 1 | 0.7497494 | 1.9993628 | 1.2496134 | 6.4e-4 | 0.010416 |
| 1/4 | 0.7499914 | 2.0001512 | 1.2501597 | 1.6e-4 | 0.002604 |
| 1/16 | 0.7499999 | 2.0000692 | 1.2500693 | 7.0e-5 | 0.000651 |
| 1/64 | 0.750000005 | 2.0000196 | 1.2500196 | 2.0e-5 | 0.000163 |

**Certified:** Δγ(a) ≥ 1 at every a, and the two-sided pinch |Δγ(a) − 5/4| ≤ ε(a) with ε strictly decreasing — the content-1 channel's excess rate lands on the **free Casimir gap C₁ − C_½ = 5/4** in the cutoff direction, while the face ladder collapses (r₁/r_½ ↓, YM-22 tie). Session item (b) is now closed in both directions: fixed-a ordering with margins (YM-39) and cutoff-direction Casimir pinch (here). Build notes: the draft's "approach from below" was refused by the brackets (Δγ crosses 5/4; γ_½ crosses 3/4 — signed dressing), and the first run violated the standing round-before-log rule (12 min → 4 min). **Boundary:** order-1 dressed two-point, not the operator gap; declared trajectory; E4D-C OPEN. Pin `EXPECTED_YM40.sha256`.

## YM-41 — The tools dock: E4D-C restated with the operator-branch tools

The RKF operator-branch tools consumed by pin and re-verified on this repo's own pencil: **silence-channel contraction** β_sil < 1 at every κ (0.00256 / 0.01029 / 0.04173 — identical to RKF's digits, both roundings coincide at 1e-9), the **ladder law** f_{c+½}/f_c ≤ κ/(2(2c+2)) with the sharper-claim discrimination (tight to 0.07%), and the **seam-flow meter** (sector cascade on the column pencil, every bracket's count route agreeing with the native det route). Selective release (stationarity, seam-count retention, u+e < 1 for the untruncated column) consumed by pin — RKF theorum/75 is the derivation of record. **Current scope correction:** YM-37/38 control spatial convergence on the declared column; T75 supplies a two-rail release conditional on its declared higher-content Gram budget. The time-transfer operator upper bound remains open. [YM-42](YM42_SILENCE_TIME_LADDER.md) now supplies the finite t=3 silence contraction, with an explicit non-implication witness showing why spatial contraction alone is insufficient. Neither t-uniformity nor a three-rail content release is claimed. Pin `EXPECTED_YM41.sha256`.

## NG-1 — The native gap statement, its teeth, and the dictionary declared open

A governance capsule with **no new numbers**. Diagnosis behind it: every wall this program has hit (cutoff uniformity, OS, infinite volume, Haar↔Φ_Σ) sits exactly where a native object meets a classically-stated goal — progress has been fast inside and zero at that boundary. So the statement being pursued is now written in the framework's own ontology, and the translation is named as a separate open problem, explicitly rather than implicitly.

**NG.** There is γ\* > 0 and a declared trajectory κ(a) such that every member of the recognition fabric family has recognition gap γ(a, m, t) ≥ γ\* for all a, m, t — each instance certified by native verdicts (elimination sign patterns, cut-tail mass, recognition energy, seam counts), every infinite direction taken through theorum/28's outward certificate, never a completion or a limit axiom. **NG is not the Clay statement and is never claimed to be.**

**Teeth** (a statement with no refutation condition is a definition, not a theorem) — five ways to refute NG with a single certified counterexample: volume (R1, a route already died this way in YM-11), cutoff (R2, live — YM-23 shows the tiling route exits at a = 1/8), sector cascade (R3, the seam-flow meter is the instrument), a channel lighter than J = ½ (R4), and any verdict that survives its own tamper (R5).

**Ledger re-scored against NG:** O1 single-member gap DELIVERED; O2 volume DELIVERED on the column and, by YM-45, for the small-bridge full heat chain; the original Wilson-family operator upper bound remains OPEN; O3 content CONDITIONAL on T75's two-rail/tail contract on that route; **O4 time DELIVERED on the declared full heat chain for abs(theta)<1/1680 and 0<a<=1 by YM-45; the generic NG/Wilson claim remains OPEN**; **O5 cutoff PARTIAL: YM-44 constructs full time refinement at fixed width, and YM-45 supplies the same width-independent gap rate on its small-bridge window; YM-46 constructs the fixed-a infinite-volume state/time map; YM-47 constructs the joint local-history limit; YM-48 constructs the reflected generator and coefficient-row embedding; full row/observable identification, physical trajectory and AF remain OPEN**; O6 evidence standard IN FORCE; **DICT (translation to Clay) OPEN, separate problem** (RH T01-E5C/E6). Sharpened by today's finding: theorum/75's grading was built classically and then regraded natively with identical rationals — one classical-looking ingredient proved replaceable; whether the measure itself is, is exactly DICT. Pin `EXPECTED_NG1.sha256`.
