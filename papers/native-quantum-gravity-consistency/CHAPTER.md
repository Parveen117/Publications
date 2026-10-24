# Quantum-consistent cut transport in a finite information–curvature model

Research chapter, 14 October 2026. Target venue requested: Nature Physics.
Current status: exact finite results and conditional predictions; a complete
physical theory of quantum gravity has not been established.

## Abstract

We examine the quantum consistency of a native cut/operator construction
coupled to a finite four-dimensional Regge geometry. Earlier results give
two exact linear curvature-response modes, a local nonlinear stationary
branch and a family of source laws with the same first response but
different higher response and stability. We show that the selected
state-dependent cut relaxation cannot be promoted to an affine quantum
evolution: two preparations of the same mixed state have unequal averaged
initial velocities. A fixed exponential filter admits an explicit joint
unitary realization when success, failure and source amounts are retained.
The two response observables also require added noise in an unbiased
commuting joint readout; an exact finite construction saturates the bound.
Finally, two initial volume slopes identify the nonlinear source coefficient
and the unknown source scale, while a third slope provides a held-out
prediction. These results make the operational and empirical contracts
precise. They do not yet construct independent quantum geometry or its
continuum gravitational dynamics.

## 1. What the construction already supplies

The native programme starts from composable cuts, transport, retained
records and invariant readouts. A cut is an idempotent operator, not
automatically a physical qubit; the trace/probability interpretation is
an additional operational contract. The precise mathematical content
should be distinguished from an analogy about information being bound
by gravity.

The existing work has the following dependencies:

| Packet | Reused result | Scope that remains attached |
|---|---|---|
| QTH-1 / curvature-information duality | Information and curvature structures in declared finite statistical models | Does not identify every information metric with physical spacetime curvature |
| CID-D2 / native cut-transport feedback | Selected coupled cut/transport action and conserved quantity | Choice of action and operation parameter is declared |
| CID-D3 / native coupling reconstruction | Identifiability, degenerate protocols and held-out observations | Coefficients are recoverable only under sufficient independent response data |
| NJC-1 / native Jacobian-cut certification | Cut-compatible transport and independent normal-connection freedom | A cut and its derivative record do not determine all curvature |
| N4H-1 / four-dimensional hinge | Finite geometric gluing, lifted curvature and selected phase contract | Euclidean discrete construction; no quantum path measure or continuum limit follows |
| EGR-1 / exponential gravity response | Rank-two response with stiffness `diag(10,47)` and a selected cut relaxation | Response coefficients are fixture dependent; gradient flow is not asserted to be fundamental quantum evolution |
| EGR-2 / nonlinear gravity continuation | Local analytic stationary branches and exact response jets | Local statement away from declared degeneracies, not global singularity resolution |
| CRC-1 / cut-record completeness | Minimal records for specified operations and source-law ambiguity | Completeness depends on the next operation and target reading |
| QGC-1 / this chapter | Ensemble obstruction, outcome-preserving extension, readout noise and source-law identification | Finite quantum operational model and conditional predictions |

This continuation reuses those results. It does not rebrand established
quantum instruments or joint-measurement limitations as newly discovered
physical laws.

## 2. Ensemble consistency of the selected relaxation

The old selected pure-cut relaxation has

\[
 F(P)=[[P,H(P)],P],\qquad
 H(P)=-\left(\frac{r_z}{10}Z+\frac{r_x}{47}X\right).
\]

It preserves a pure cut's idempotence. That fact alone is insufficient
for quantum dynamics of an unknown source. Randomize between the two
`Z` cuts with weights `4/5,1/5`, or between the tilted cuts with directions
`(3/5,4/5)` and `(3/5,-4/5)` with equal weights. Both preparations give
exactly `rho=diag(4/5,1/5)`. Their mean initial velocities are respectively

\[
 0\quad\hbox{and}\quad(888/29375)Z.
\]

Thus no preparation-independent affine quantum operation can reproduce
all these trajectories. The old certificates prove the stated gradient
identities; they do not certify the stronger physical interpretation.
The inconsistency appears only when that interpretation is added.

Calling the preparation label hidden information does not by itself
remove the problem. A theory that retains such a label must specify its
physical status, its evolution, its observable effects and compatibility
with the quantum statistics it intends to recover. Alternatively, a
known external preparation record can legitimately control a protocol,
but then that enlarged record is part of the input.

## 3. What the exponential can consistently do

Fix a cut `Q` and write `Qbar=I-Q`. For `a^2+b^2=1`, the two operators
`Ks=Q+a Qbar`, `Kf=b Qbar` form a complete instrument. Their dilation

\[
 W=Q\otimes I+\bar Q\otimes
       \begin{pmatrix}a&-b\\b&a\end{pmatrix}
\]

is unitary in the native algebra. With `a=exp(-s)`, the successful filter
is `exp(-s Qbar)`, and its nonselective reduction is an affine dephasing
semigroup. The exponential supplies a continuous family of operations;
it does not determine their interpretation as gravitation.

For the displayed rational source and `a=3/5`, the successful direction
becomes `(91/109,60/109)`, but its probability is only `109/125`. The full
source has direction `(3/5,12/25)`. This gives a precise version of the
requirement to retain hidden scale and outcome information: use the
unnormalized branch sources, then add them. A larger normalized alignment
is not a larger unconditional source.

The full joint operation preserves information before any outcome
register is discarded. A reduced or postselected record generally does
not. This extension does not keep the prohibited state-dependent
relaxation as an autonomous unknown-state law.

## 4. A quantitative quantum restriction on the response record

At fixed source amount, the mean response implies source observables
`A1=chi Z/10`, `A2=chi X/47`. They obey

\[
 [A_1,A_2]=\iota\chi^2Y/235.
\]

Any commuting joint pointer unbiased on all source states therefore
requires added second-moment noise with

\[
 \|\mathsf N_1\|\|\mathsf N_2\|\ge\chi^4/220900.
\]

The proof is the compression identity
`[A1,A2]=L2^dagger L1-L1^dagger L2`, with
`N_i=Li^dagger Li`. A four-outcome exact native dilation saturates the
bound. This identifies what a classical readout must lose or broaden
when it attempts to encode both response modes simultaneously.

This is a restriction on the declared readout of a quantum source.
It is not yet a quantum state space for geometry, a universal length
uncertainty, or evidence that an independently propagating metric has
been quantized.

## 5. Turn an arbitrary source coefficient into an identifiable model

CRC-1 found that a quadratic coefficient `kappa` changes the nonlinear
source law while leaving the first curvature response unchanged. Merely
choosing a stable member does not derive the physical value. For the
same fixed geometry and source, QGC-1 instead derives an inverse protocol.

Let `v_o,v_cap,v_in` be initial log-volume slopes of the three cell
classes, measured against one common source control. Two slopes identify
`kappa` and the unknown source scale. The remaining slope must satisfy

\[
 v_{\rm in}=\frac{640}{141}v_o
       -\frac{1860931}{621932529}v_{\rm cap}.
\]

The relation is independent of those two fitted parameters. It and the
equality of slopes within each symmetric class are genuine predictions
of this declared finite family. A discrepancy can reject the family even
when the two calibration channels are fitted perfectly.

The volume slopes are derivatives at zero source, not measurements at
an arbitrary finite control setting. A usable experimental claim must
include the finite-range bias, parameter uncertainty and a physical
measurement corresponding to each registered cell volume. Those items
are not supplied by exact arithmetic alone.

## 6. Audit of the requested claim “quantum gravity completely solved”

| Necessary part of that physical claim | Evidence currently present | What has not been established |
|---|---|---|
| A specified source law | A finite family and a two-channel inverse protocol | A physical principle or observation selecting its law and coefficients |
| Consistent quantum operations | Fixed finite instruments, positive amplitudes and a complete outcome ledger | One joint quantum dynamics of matter and independent geometry reproducing the earlier nonlinear response |
| Gravitational degrees of freedom | A finite Euclidean Regge geometry and response tensors | Lorentzian causal dynamics, the physical constraint structure and relational observables |
| Recovery of known gravity | Local stationary response of one complex | A controlled regime recovering general relativity, including curved classical solutions and the appropriate weak-field limit |
| Control across scales | Exact local identities at one discretization | Refinement/continuum control or a physically justified fundamental discretization with controlled macroscopic predictions |
| Physical units and matter coupling | Dimensionless source, cut and geometry variables | A derived relation to measured stress-energy, length, action scales and gravitational coupling |
| Distinguishing predictions | Instrument statistics, a joint-readout noise bound and a held-out finite-geometric slope | A calibrated gravitational observable that discriminates this theory from alternatives |
| Singularities | Local admissibility and concentration relations in predecessor packets | Global evolution, accessible singular sets or a proof of singularity resolution |

These are substantive mathematical and physical gaps, not missing
formatting sections. In particular, deleting a primitive time coordinate
does not remove the need to define physical correlations or causal
structure; a relational clock can replace a coordinate only after its
role in the dynamics is specified.

Likewise, the cut spectrum or amount alone is not yet a stress-energy
tensor. Two normalized pure cuts have the same spectrum but can have
different directions relative to the geometry/readout. A statement that
“more invariant information bends more curvature” needs a specified
invariant, source map, geometric observable and monotonicity domain.
The model-dependent response does not prove that statement universally.

## 7. Submission position and the next decisive result

The present defensible contribution is a finite consistency and
identifiability analysis. A title, abstract or cover letter claiming
“complete solution of quantum gravity” would exceed the results.
The declared finite alternatives also reproduce established quantum
information structures, so their presence alone is not a claim of broad
new physical novelty.

[Nature Physics's publication criteria](https://www.nature.com/nphys/editorial-policies/peer-review)
ask for a consequential advance supported by strong evidence. This does
not mean every theoretical paper must already have experimental detection.
It means the evidence must support the specific claim being made. The
current certificates support the finite theorems in this chapter; they
do not support the larger physical claim. No editorial acceptance or
submission readiness has been certified, and no manuscript has been sent.

The next decisive derivation is a single joint matter–geometry model:
identify its independent geometry algebra, specify a state-independent
quantum operation or constrained amplitude rule, derive the earlier
source response as a controlled conditional or semiclassical limit,
and calculate one calibrated relational gravitational observable. Its
parameters must be fixed separately from that observable. A fixed
instrument is a consistency building block for this task, not the task's
completion.

Until such a result exists, this chapter closes the finite claims stated
above and records a precise boundary. It does not declare the physical
quantum-gravity programme finished.

## Reproducibility and references

The proof and calculation details are in [THEOREM.md](THEOREM.md), with
the prospective tests in [PREDICTIONS.md](PREDICTIONS.md). The native
certificate uses the original engines and pinned predecessor results.
Replay instructions are in [README.md](README.md). General analytic and
operator-norm arguments are written proofs, not proof-assistant outputs.

Relevant prior literature is cited for context, without a claim that the
present formulas were experimentally observed:

- Simon, Buzek and Gisin, [The no-signaling condition and quantum dynamics](https://arxiv.org/abs/quant-ph/0102125).
- Busch and Heinosaari, [Approximate joint measurements of qubit observables](https://arxiv.org/abs/0706.1415).
- Bose et al., [A Spin Entanglement Witness for Quantum Gravity](https://arxiv.org/abs/1707.06050).
- Marletto and Vedral, [Gravitationally-induced entanglement between two massive particles is sufficient evidence of quantum effects in gravity](https://arxiv.org/abs/1707.06036).

