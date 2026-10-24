# MGH-1 — A finite matter–geometry holonomy model

14 October 2026. All quantities are dimensionless unless a reference unit
is displayed. The native rational cut field, trace, dagger, matrix operations
and word compiler are reused without modification.

## 1. Local matter and its observational invariants

Let `Pi_o` be a matter localization cut: `Pi_o^2=Pi_o=Pi_o^dagger`.
In a fixed rest-mass sector define

\[
 M_o=\mu\Pi_o,\qquad \mu=m/m_*,\qquad
 p=\operatorname{Tr}(\rho_m\Pi_o),\qquad
 I_n=\operatorname{Tr}(\rho_m M_o^n)=\mu^n p.                  \tag{1}
\]

Here `rho_m` is normalized. `m_*` is a declared reference/calibration
scale. The two matter labels mean outside/inside the selected observation
cell; the particle's rest mass need not itself be in superposition.
`Pi_o` selects the interaction location. Projectively measuring it first
would be a different preparation and can remove the required coherence.

For a general positive local energy operator one can instead use
`M_o=Pi_o H_m Pi_o/E_*`. This specifies a compressed energy observable;
it does not assert that compression retains every moment of the original
uncompressed Hamiltonian. The general result below applies directly to
the chosen finite Hermitian `M_o`.

Under a common matter-frame change, transform both `rho_m` and `M_o` by
conjugation. Every `I_n` and `p` is invariant. For an unnormalized source
`N rho_m`, retain `N`; its unnormalized moments are `N I_n`. Normalize
only when probabilities or entropies are requested.

The trace probability rule and tensor composition are the explicit finite
quantum adapter already used in QGC-1. A cut algebra alone does not select
that physical interpretation or a value of the mass-to-phase coupling.

## 2. Actual four-dimensional metrics fix the phase

Reuse the N4H-1 hinge star with six cells

\[
 \sigma_i=[A,B,C,D_i,D_{i+1}],\quad i\bmod6.
\]

In the local frame of each cell take `A=0`, `B=e1`, `C=e2`, `D_i=e3`,
`D_(i+1)=c_i e3+s_i e4`, where

\[
 c_i=\frac{1-t_i^2}{1+t_i^2},\quad
 s_i=\frac{2t_i}{1+t_i^2},\quad t_i>0.                         \tag{2}
\]

Every cell has positive volume `s_i/24`. Repeated edge lengths agree.
The shared-tetrahedron Gram matrices agree, and N4H-1's unique unfolding
isometry is

\[
 T_{i+1,i}=F(E^\dagger E)^{-1}E^\dagger-n_tn_s^\dagger.
\]

In these frames it equals `P_h+c_i Q_h-s_i R_h`, with
`P_h=diag(1,1,0,0)`, `Q_h=I-P_h`, and `R_h e3=e4`, `R_h e4=-e3`.
Thus the actual loop transport is

\[
 H_h=P_h+\cos\delta_h Q_h+\sin\delta_h R_h,\quad
 \delta_h=2\pi-\sum_i\theta_i,\quad\theta_i\in(0,\pi).
                                                               \tag{3}
\]

Take two metrics on this same triangulation:

| Configuration | Cayley parameters |
|---|---|
| `q0` | `(1/2,1/2,3/4,1/2,1/2,3/4)` |
| `q1` | `(1/2,1/2,2/11,1/2,1/2,2/11)` |

Let `alpha` have cosine `3/5` and sine `4/5`, and `beta` cosine `7/25`
and sine `24/25`. Both lie in `(0,pi/2)`. Their exact phase products give
`2 alpha+beta=pi`. The angle with parameter `2/11` has cosine `117/125`
and sine `44/125`; it is `beta-alpha`, which is positive because
`cos beta < cos alpha`. Consequently

\[
 \delta(q_0)=0,\qquad \delta(q_1)=2\alpha,\qquad
 z=e^{\iota\delta(q_1)}=-\frac7{25}+\iota\frac{24}{25}.
                                                               \tag{4}
\]

The compact hinge defect and total four-volumes, at unit reference length,
are

\[
 \kappa_h=\tfrac14\|H_h-I\|_F^2=1-\cos\delta_h=\frac{32}{25},
 \qquad V_0=\frac{16}{75},\quad V_1=\frac{61}{375}.             \tag{5}
\]

The hinge area is `1/2` in both configurations. The two radial boundary
edges corresponding to the changed parameters differ. The configuration
register therefore includes the complete region and its boundary metric;
these are not claimed to obey one externally fixed boundary condition.
Moving the two changed parameters continuously between their positive
endpoints supplies a nondegenerate metric interpolation.

The character in (4) is extracted from the geometric transport, not
inserted as an unrelated matrix entry. In the native complex completion,

\[
 E_+=\tfrac12(Q_h-\iota R_h),\qquad
 H_h E_+=z E_+,\qquad \operatorname{Tr}(E_+H_h)=z.              \tag{6}
\]

Independent cell-frame transformations conjugate the loop and the normal
eigen-cut at its base point. The trace in (6) is unchanged. Common spatial
scaling also leaves the phase and defect unchanged; coordinate volumes
scale with the fourth power of length.

## 3. Specify the quantum configuration model

Introduce orthogonal configuration cuts `Q0,Q1` representing the two
complete metrics. This is a declared finite quantum configuration algebra.
Its diagonal geometric observables are

\[
 \widehat\kappa=\kappa_h Q_1,\qquad
 \widehat V=V_0Q_0+V_1Q_1,\qquad Z_g=Q_0-Q_1.
                                                               \tag{7}
\]

Metric gluing supplies these geometric values and the phase `z`.
Coherent superposition and transitions between the labels require a
quantum operation prescription. Choose the native mixing

\[
 B=\frac{1+\iota}{2}(I+R_g),\quad
 R_g=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
 BB^\dagger=I,\quad BZ_gB^\dagger=X_g.
                                                               \tag{8}
\]

Here `X_g` is the off-diagonal configuration reflection. This mixing is
an active operation on metric labels, distinct from a coordinate-frame
change on a four-simplex. Static metric data alone do not fix it.

Select the mass-weighted holonomy phase interaction

\[
 D_\delta=\exp(\iota\delta M_o\otimes Q_1),\qquad
 W_\delta=(I\otimes B)D_\delta(I\otimes B^\dagger)
          =\exp(\iota\delta M_o\otimes Q_X),\quad
 Q_X=\tfrac12(I-X_g).                                         \tag{9}
\]

The coefficient multiplying `M_o delta` is set to one as a model
calibration. A physical theory must determine or independently measure
this coefficient/reference scale. Equation (9) is not inferred uniquely
from the existence of a holonomy.

**Theorem 1 (fixed joint operation and conserved source moments).**
`W_delta` is unitary, its generator is independent of the input state,
and it commutes with `M_o^n tensor I` for every nonnegative integer `n`.
It preserves the joint spectrum, trace and every source moment `I_n`.
Equivalent density-operator ensembles have identical joint outputs.

**Proof.** The generator `M_o tensor Q_X` is Hermitian. Its exponential
is unitary, and it commutes with every displayed source observable.
Cyclicity gives moment conservation. Conjugation by a fixed unitary is
linear in the density operator, proving the ensemble claim. In the
projector-mass case (9) can also be verified directly as a sum of
orthogonal source blocks. QED.

This supplies a consistent new finite operation; it does not retain the
state-dependent gradient trajectories excluded by QGC-1.

### 3.1 Recover the proposed coupling with the correct half-angle

For `M_o=mu Pi_o`, equation (9) factors as

\[
 W_\delta=\exp(\iota\delta M_o/2)\otimes I\;
          \exp(-\iota(\delta/2)M_o\otimes X_g).               \tag{10}
\]

The first factor is a local matter phase and cannot create entanglement.
The second is exactly the proposed chat coupling, with
`lambda=delta/2`. For the new geometry, its unit-mass half-angle phase
is `(3+4 iota)/5`. This relation fixes the geometric factor of two.

No primitive clock is needed to define the operation. Interpreting a
sequence of these operations as physical time evolution would require
an additional dynamical/clock prescription.

## 4. Exact information–curvature relation

Initially take a pure matter cut `rho_m` and the flat configuration cut
`Q0`, with no matter–geometry correlation. Define the matter characteristic
function for the chosen coupling by

\[
 \chi(\delta)=\operatorname{Tr}(\rho_m e^{\iota\delta M_o}).
                                                               \tag{11}
\]

**Theorem 2 (general finite-source response).** After (9), the geometry
state in the `(q0,q1)` basis is

\[
 \rho_g=\frac12
 \begin{pmatrix}
 1+\Re\chi&\iota\Im\chi\\
 -\iota\Im\chi&1-\Re\chi
 \end{pmatrix}.
\]

Hence

\[
 \boxed{S_L=1-\operatorname{Tr}\rho_g^2
       =1-\operatorname{Tr}\rho_m'^2
       =\frac{1-|\chi|^2}{2},\qquad
       \Pr(q_1)=\frac{1-\Re\chi}{2}.}                         \tag{12}
\]

**Proof.** In an eigen-cut of `M_o` with eigenvalue `m_j`, the geometry
amplitudes are `a_j=(1+e^{iota delta m_j})/2` and
`b_j=(1-e^{iota delta m_j})/2`. Tracing the matter register gives the
displayed matrix. Its squared trace is `(1+|chi|^2)/2`. The full output
is pure, so the two reductions have the same nonzero squared singular
values and therefore the same purity. QED.

The geometry matrix formula remains valid for mixed matter. Its linear
entropy is an entanglement measure here only with the stated globally
pure input. Section 7 gives a separable mixed-source control with the
same geometry marginal.

For the local mass cut, `chi=1-p+p exp(iota mu delta)`. Thus

\[
 S_L=p(1-p)[1-\cos(\mu\delta)],\qquad
 \Pr(q_1)=\tfrac p2[1-\cos(\mu\delta)].                        \tag{13}
\]

**Corollary (unit-mass information–hinge relation).** For `mu=1` and
the flat reference used above,

\[
 \boxed{S_L=p(1-p)\kappa_h,\qquad
 \langle\widehat\kappa\rangle=\tfrac p2\kappa_h^2.}           \tag{14}
\]

Here `kappa_h` is the defect of the curved basis geometry; the second
quantity is the defect expectation in the evolved configuration state.
They have different meanings. For a nonflat reference, the same phase
argument involves the relative deficit, under a specified comparison of
the hinge normals.

The finite-state expectation in (14) is an output probability response
of the selected unitary. It is not an Einstein/Regge field equation or a
derivation of the earlier EGR stationary branch.

### 4.1 Connect to the native configuration-loop curvature

Using `G=M_o tensor X_g` and `J=I tensor Z_g`, the native infinitesimal
configuration-loop commutator is

\[
 \Omega=[-\iota G,-\iota J],\qquad
 \Omega^\dagger\Omega=4M_o^2\otimes I.
\]

Subtracting the mean source only in this diagnostic, not in the
physical evolution generator, gives the connected squared strength

\[
 \mathcal C=4(I_2-I_1^2).
\]

For the projector source, with `lambda=delta/2`,

\[
 \boxed{S_L=\frac{\lambda^2}{2}\mathcal C
       \left[\frac{\sin(\mu\lambda)}{\mu\lambda}\right]^2.}    \tag{15}
\]

The ratio has its continuous value one at zero. Equation (15) is the
chat relation, now linked through (14) to the explicit geometric hinge.
`Omega` is curvature of two configuration operations; `kappa_h` is the
compact defect of a piecewise Euclidean metric. The equations relate
them in this selected model without identifying the two objects.

## 5. Exact quantum demo and geometric observables

Take `mu=1`, matter cut `(I+X_m)/2` and initial geometry `Q0`.
The joint output remains pure. Local diagonal phase calibrations,
which commute with the occupancy and geometric-defect readings, put
its amplitude array in the form

\[
 \Psi=\frac1{\sqrt2}
 \begin{pmatrix}1&0\\3/5&4/5\end{pmatrix}.
                                                               \tag{16}
\]

Rows label matter occupancy, columns metric configurations. The
certificate stores the density matrix, whose entries are all native
rational cuts; it does not numerically approximate `sqrt(2)`.
Direct partial traces yield

\[
 \operatorname{Tr}\rho_m'^2=\operatorname{Tr}\rho_g^2=17/25,
 \quad S_L=8/25,\quad\Pr(q_1)=8/25,
\]
\[
 \langle\widehat\kappa\rangle=256/625,\qquad
 \langle\widehat V\rangle=616/3125.                            \tag{17}
\]

The volume expectation follows from the same metric-label probability;
it is not an additional independent fitted observable. The quantum
witness below tests coherence that the diagonal geometric readings
alone cannot identify.

## 6. A five-setting entanglement witness

In the phase-calibrated basis of (16), put
`v=(-1,-3,3,-1)^T`, `P_v=v v^dagger/20`, and
`W_ent=(P_v)^(T_g)`, with partial transpose in the geometry register.

**Theorem 3 (operational witness).** Every separable state has
`Tr(W_ent rho)>=0`, while (16) has

\[
 \operatorname{Tr}(W_{\rm ent}\rho)=-2/5.                      \tag{18}
\]

**Proof.** For a product positive state, its geometry partial transpose
is still positive. Convex sums preserve positivity. Therefore its trace
against the positive rank-one `P_v` is nonnegative. For (16), exact
multiplication gives `rho^(T_g) v=-(2/5)v`, proving (18). QED.

With the native Pauli reflections and `Y=iota R`,

\[
 W_{\rm ent}=\tfrac14\left[
 I\otimes I+Y\otimes Y-\tfrac45(X\otimes X+Z\otimes Z)
 -\tfrac35X\otimes Z+\tfrac35Z\otimes X\right].                \tag{19}
\]

The five correlators `XX,YY,ZZ,XZ,ZX` are respectively
`4/5,-4/5,16/25,3/5,-12/25`. For the specified white-noise family
`rho_v=v rho+(1-v)I4/4`,

\[
 \langle W_{\rm ent}\rangle_v=(5-13v)/20,
\]

so this witness detects entanglement for `v>5/13`. The threshold concerns
this noise family and this witness. It is not an experimental accuracy
claim or a universal noise threshold.

## 7. Records that cannot be discarded

**Source coherence.** Replace the pure matter cut in section 5 by `I/2`.
Its occupancy probability and all moments of the projector mass remain
the same. The resulting geometry marginal is also identical. Nevertheless
the joint output is explicitly

\[
 \tfrac12\Pi_0\otimes Q_0
 +\tfrac12\Pi_1\otimes|g_1\rangle\langle g_1|,
\]

and is separable. The exact certificate constructs this decomposition.
The source invariant alone cannot certify entanglement; preparation
coherence or a joint witness must be retained.

**Higher energy moments.** For a general bounded finite `M_o`, (12) uses
its whole characteristic function. Expanding it gives

\[
 S_L=\frac{\delta^2}{2}(I_2-I_1^2)
 -\frac{\delta^4}{24}(I_4-4I_1I_3+3I_2^2)+O(\delta^6).        \tag{20}
\]

This follows by expanding the even characteristic function of the
difference of two independent draws from the spectral weights of `M_o`.
The projector relation `M_o^2=mu M_o` closes this moment tower exactly,
which explains why (13) needs only `mu,p` in the pure preparation sector.

For an exact counterexample to a universal two-moment closure, take
`M_o=diag(0,1,2,3)` and pure density operators `a a^dagger/30`,
`b b^dagger/30`, with `a=(1,5,0,2)` and `b=(2,4,3,1)`.
Both have `I1=37/30`, `I2=61/30`, and variance `461/900`.
At the geometric phase (4), their linear entropies differ:

\[
 S_L(a)=848872/3515625,\qquad S_L(b)=1463272/3515625.
\]

Their fourth-order brackets in (20) differ by `16/25`. This is an
explicit instance of the programme's existing operation-dependent
record-completeness principle.

**Curvature lift.** N4H-1 already gives two metrics with deficits `0`
and `-2 pi` but the same compact holonomy. A fractional mass-weighted
phase, here `mu=1/2`, distinguishes them as `+1` and `-1` and yields
`S_L=1/2` in the corresponding coherent configuration protocol.
Therefore the full geometric lift must be retained for arbitrary real
mass weights. Integer phase characters descend to compact holonomy.

**Phase wrapping.** Equation (13) is not globally monotone in mass.
For fixed positive deficit it increases on `0<mu delta<pi` and then
turns. The certificate includes a mass-2/mass-3 counterexample to an
unrestricted monotonicity claim. This separates an increasing unwrapped
phase from its periodic compact reading.

## 8. What is derived, selected, and still physical work

| Item | Status |
|---|---|
| Two positive four-dimensional metrics and their frame-covariant holonomies | Derived using N4H-1's geometric contract and exact native arithmetic |
| Compact curvature, lifted angle and geometric volumes | Derived; the two boundary metrics are included in the register |
| Coherent quantum metric labels and mixing `B` | Declared finite quantum configuration model |
| Mass-weighted phase and its normalization | Selected interaction prescription |
| State-independent joint operation, mass moments, entropy identities and witness | Derived under those premises |
| Half-angle realization of the preceding chat coupling | Derived exactly from the displayed geometric phase |
| Physical quantum metric preparation, Newton/Einstein limit and measured coupling | Open |
| Derivation of the EGR nonlinear stationary response from this joint dynamics | Open |

The geometry register now has independently specified metric observables,
and its state changes under a joint unitary with matter. This closes a
finite modelling gap left by the abstract two-generator chat example.
The existence of this model does not select it as nature's gravitational
dynamics. A quantum-device implementation would demonstrate the encoded
model; it would not on its own demonstrate quantum gravitational physics.

## 9. Proof status and references

The certificate rechecks all shared-edge and whole-face metric conditions,
actual transport products, frame changes, exact quantum states, source
moments, polynomial identities and negative controls. Symbolic witnesses
are serialized and replayed by the original native compiler. A separate
read-only run recomputes the certificate and checksum. Angle branches,
analytic exponentials, positivity and the arbitrary-source theorem have
written proofs here; no general-purpose proof-assistant verification is
claimed.

Controlled phases, partial-trace purity and entanglement witnesses are
established quantum methods. Relevant context is
[I-Sheng Yang, The Entanglement Timescale](https://arxiv.org/abs/1707.05792)
and [Asher Peres, Separability Criterion for Density Matrices](https://arxiv.org/abs/quant-ph/9604005).
The geometric realization uses the existing N4H-1 Regge adapter. The
source-record controls reuse CRC-1, and ensemble consistency reuses the
contract of QGC-1. These general methods are not claimed as new discoveries.
