# MGH-1 — Reproducible matter–geometry quantum demo

This is an exact finite simulation and a prospective measurement protocol.
No hardware job or gravitational experiment was performed for this packet.

## 1. What the two registers represent

| Register | Label 0 | Label 1 |
|---|---|---|
| Matter occupancy | Matter outside the chosen observation cell | Matter inside that cell |
| Geometry configuration | Six-cell metric `q0`, selected hinge deficit zero | Six-cell metric `q1`, selected hinge phase `(-7+24 i)/25` |

The particle has a fixed mass ratio `mu=m/m_*=1` in the main demo. Its
local mass operator is the occupancy cut. The two geometry labels include
the whole finite metric region and its boundary, whose lengths differ.
The basis geometry observables are

\[
 \widehat\kappa=\operatorname{diag}(0,32/25),\qquad
 \widehat V=\operatorname{diag}(16/75,61/375).
\]

These diagonal values come from the actual simplex metrics and geometric
transport calculation. Their coherent quantum encoding and mixing are
declared parts of the model.

## 2. Prepare and apply the fixed operation

Prepare matter in the pure cut `(I+X)/2` and geometry in `Q0`.
Let `B=((1+i)/2)(I+R)`, with `R=[[0,-1],[1,0]]`. In basis order
`00,01,10,11`, take

\[
 D=\operatorname{diag}(1,1,1,z),\qquad
 z=(-7+24\iota)/25,\qquad
 W=(I\otimes B)D(I\otimes B^\dagger).
\]

For conventional right-to-left application, perform `B^dagger` on the
geometry register, then `D`, then `B`. This entire sequence is one
specified unitary and is independent of the unknown input state.

For the witness basis below, finally apply the local diagonal phases

\[
 L_m=\operatorname{diag}(1,(3-4\iota)/5),\qquad
 L_g=\operatorname{diag}(1,\iota).
\]

These phase calibrations commute with occupancy, compact-defect and
volume readouts. They do not change entanglement.

The resulting exact density matrix is

\[
 \rho=\frac12
 \begin{pmatrix}
 1&0&3/5&4/5\\
 0&0&0&0\\
 3/5&0&9/25&12/25\\
 4/5&0&12/25&16/25
 \end{pmatrix}.
\]

## 3. Predictions before fitting any output

| Observable | Prediction |
|---|---:|
| Matter occupancy probability | `1/2`, conserved |
| Local mass first moment | `1/2`, conserved |
| Local mass second moment | `1/2`, conserved |
| Joint purity | `1` |
| Either reduced-state purity | `17/25` |
| Entanglement linear entropy | `8/25` |
| Curved geometry probability | `8/25` |
| Flat geometry probability | `17/25` |
| Mean compact hinge defect | `256/625` |
| Mean four-volume | `616/3125` |

For another pure matter preparation with occupancy probability `p`, the
same unit-mass operation predicts

\[
 S_L=\frac{32}{25}p(1-p),\qquad
 \Pr(q_1)=\frac{16}{25}p.
\]

Changing the dimensionless mass ratio replaces the phase by
`exp(i mu delta_h)`. Integer masses use powers of the certified phase.
Fractional masses require the retained geometric angle lift. The
certificate includes masses `0,1/2,1,2,3` and several pure preparations.

## 4. Detect quantum correlations with five settings

Measure the following joint correlators on independently prepared copies
of the phase-calibrated state:

| Setting | Exact expectation | Coefficient in the witness |
|---|---:|---:|
| `X_m X_g` | `4/5` | `-1/5` |
| `Y_m Y_g` | `-4/5` | `1/4` |
| `Z_m Z_g` | `16/25` | `-1/5` |
| `X_m Z_g` | `3/5` | `-3/20` |
| `Z_m X_g` | `-12/25` | `3/20` |

Add `1/4` to the weighted sum. The predicted result is

\[
 \boxed{\langle W_{\rm ent}\rangle=-2/5.}
\]

All separable states give a nonnegative value under the declared trusted
measurement operators. THEOREM.md proves this by partial transpose and
an explicit positive rank-one projector. This is an entanglement witness,
not a Bell test or a device-independent conclusion.

For `rho_v=v rho+(1-v)I4/4`, the exact prediction is

\[
 \langle W_{\rm ent}\rangle_v=(5-13v)/20.
\]

It is negative for `v>5/13`. A real implementation must estimate the
correlators and their uncertainty, including preparation and measurement
calibration. Demonstrate that an upper uncertainty bound on the witness
lies below zero. The native exact value is a prediction, not a measured
significance or an experimental confidence interval.

## 5. Required controls

1. **Classical source mixture.** Replace the pure matter cut by `I/2`.
   Its mass moments and final geometry marginal are identical to those
   above, but its joint output has an explicit separable decomposition.
   Geometric population or marginal mixedness alone cannot demonstrate
   entanglement.
2. **Source present with certainty.** Prepare occupancy `p=1`. Geometry
   responds, but the joint state remains a product. This separates mean
   source response from source–geometry entanglement.
3. **Zero source.** Set mass ratio zero or occupancy `p=0`. The operation
   returns the initial flat geometry and creates no entanglement.
4. **Equivalent ensembles.** Use the two decompositions of
   `diag(4/5,1/5)` from QGC-1. Their joint outputs agree exactly, checking
   that the earlier nonlinear-ensemble obstruction is avoided.
5. **Omitted moment information.** The four-level source pair in
   THEOREM.md has identical first two energy moments but different finite
   responses. A two-moment description needs the projector closure
   assumption or another independently verified closure rule.
6. **Omitted winding information.** The old flat/excess-angle pair has
   identical compact holonomy but different half-mass phases. Retain the
   geometric lift when the coupling does not descend to the compact loop.

## 6. What a demonstration establishes

An implementation on encoded quantum registers can test this exact
geometric-phase circuit, the predicted joint state and its entanglement.
It would show that the declared finite model is implementable. Inferring
quantum gravitational physics additionally requires physical access to
the proposed metric configurations and an independently justified
matter–geometry coupling. The mass normalization and configuration
mixing are choices in this packet, and no Newton/Einstein limit is inferred
from the gate identities.

The controlled phase and witness are established quantum techniques.
The exact simplex realization and the conditional relation
`S_L=p(1-p) kappa_h` supply the concrete connection to this programme's
geometric curvature result.
