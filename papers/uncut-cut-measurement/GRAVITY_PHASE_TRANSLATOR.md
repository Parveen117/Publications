# Gravity direction: recovering native source response from coherent phase

Status: **exact conditional finite translation**, with a written proof and rational negative controls. This continues QB-1–QB-3 toward U36's gravity-before-curvature source law. It does not equate that source law to measured gravity, select a physical metric, or derive `hbar`.

## GP-1. Coherent phase determines the native stationary response

Take the **same supplied** positive grounded-graph stiffness `H` in two readings:

1. U36's stationary source law `phi_*=H^{-1}b`, and its source/probe cross cost `E_cross(b,p)=-b^T H^{-1}p`;
2. QB-1's paired native continuation, represented by `U_h=(I+i h H)(I-i h H)^{-1}` for a known nonzero native step `h`.

No manifold or physical clock enters this construction. Since `H` is positive definite, `U_h-I` and `U_h+I` are invertible. Exact Cayley inversion gives

`i h H = (U_h-I)(U_h+I)^{-1}`,

`H^{-1} = i h (U_h+I)(U_h-I)^{-1}`.

Thus the **complete matrix including its global phase reference** `U_h`, together with its calibrated step `h` and source `b`, determines the native stationary field and every source/probe bilinear cross cost. Conversely, `H^{-1}` and `h` determine `U_h`. This is a two-way algebraic translator *within the supplied U36 model*. Ordinary quantum state readouts identify a unitary channel only up to a global phase; that operational distinction is load-bearing, as GP-2 shows. This is not an empirical identification of gravity.

**Proof.** Multiply `U_h(I-i hH)=I+i hH` (the two factors commute), rearrange to `(U_h-I)=i hH(U_h+I)`, and invert. `H>0`, `h!=0` imply eigenvalues of `U_h` are `(1+i h lambda)/(1-i h lambda)` with `lambda>0`, so none is `1` or `-1`. The source and cross-cost statements follow by substitution into U36. □

The scale qualification is essential: if `h` is not calibrated, `U_h` specifies only `hH`. Replacing `(H,h)` by `(kH,h/k)`, `k>0`, leaves `U_h` fixed while changing `H^{-1}` by `1/k`. A known coherent matrix without a known step and source-energy calibration therefore does not give a dimensionful force or `hbar`.

## GP-2. Even a quantum state channel can lose the source response

QB-3's classical protocol records only `P_ij=|(U_h)_ij|²` at every native step. This **does not** determine `H^{-1}`, even with known `h=1`, connected grounded graphs, positive link and grounding weights, and the same two site detectors. The explicit example below proves something stronger: the full density-matrix evolution `rho -> U_h rho U_h†` also fails to identify it unless the global phase is fixed by an independent reference.

Here are two exact grounded two-site graphs; `c` is their link stiffness and `g` each site's grounding. Both have `H=[[g+c,-c],[-c,g+c]]`:

| Graph | `g` | `c` | Eigenvalues of `H` | Off-diagonal `H^{-1}_{12}` |
|---|---:|---:|---|---:|
| A | `1` | `1` | `1,3` | `1/3` |
| B | `1/2` | `5/12` | `1/2,4/3` | `5/8` |

At `h=1`, coherent matrices are

`U_A=[[(-2+4i)/5,(2+i)/5],[(2+i)/5,(-2+4i)/5]]`,

`U_B=[[(4+22i)/25,(11-2i)/25],[(11-2i)/25,(4+22i)/25]]`.

In fact `U_B=((4-3i)/5) U_A`, with `|(4-3i)/5|=1`. Therefore both graphs yield the **same entire quantum state channel** for *every* input density matrix and every number of applications. They also yield the same one-step classical update

`P_A=P_B=[[4/5,1/5],[1/5,4/5]]`,

hence the same `P^n` for *every* repeated phase-erased experiment. Yet for unit sources/probes at distinct sites their stationary cross costs are `-1/3` and `-5/8`. **Neither the classical probabilities nor ordinary quantum state tomography** can select their native source coupling. A *phase-referenced* coherent matrix distinguishes them and the inversion above recovers each graph. This qualification prevents an unobservable global phase from being presented as a physical result.

With the shared real source `b=(1,2)`, A responds with `phi_A=(4/3,5/3)` and B with `phi_B=(21/8,27/8)`. The graph-current predictions differ as well. These are candidate source readings, not observations of gravitational mass or force.

**Proof.** Positive weights make both grounded graph matrices positive definite. Diagonalize each in the symmetric/antisymmetric site basis: the eigenvalue pairs above give the displayed Cayley entries by direct rational arithmetic. Multiplication shows the common scalar `(4-3i)/5`; it cancels with its conjugate in `U rho U†` for every density matrix. Direct inversion of `H` gives the two distinct off-diagonal entries. □

## What has and has not been derived

GP-1 is a precise **phase-referenced** quantum-form ↔ native static-source-response identity on *one declared law*. GP-2 is a stronger no-go theorem: even all density-matrix evolution data at one native step can coincide for two different grounded source laws. QB-3's seam memory remains relevant to future coherent readings, but **relative phase memory alone does not fix the global phase reference needed for GP-1**. This identifies a separate calibration obligation in the gravity direction.

These results keep the spacetime-first assumption out of the native start. To call the result physical gravity, independently fix an operational phase reference or other source-response calibration, select the native graph and law, identify universal probe coupling and measured energy/distance, recover at least a gravitational observable across held-out settings, and then test the spacetime curvature adapter. The hbar action scale is neither a prerequisite for the finite algebra nor derivable by renaming `h` a physical time. Physical quantum measurement likewise requires an operational detector adapter.

Run the [source-bound exact certificate](certificates/gravity_phase_translator.py) to check both Cayley matrices, their exact global-phase relation and equal quantum channels, full phase-referenced inversion, matching classical `P`, unequal Green responses, sources, the step-scale ambiguity and controls.

Parents: [U36 source/cross-cost law](SOURCE_RESPONSE_AND_GRAVITY_TESTS.md), [TF-1/TF-2 thermodynamic foundation](THERMODYNAMIC_FOUNDATION_BRIDGE.md), [QB-1–QB-3 quantum/classical seam](QUANTUM_CLASSICAL_BRIDGE.md).
