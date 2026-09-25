# Gravity before curvature: a conditional transport-to-geometry bridge

**Monty Dabas — v0.4 research extension, 25 September 2026.**

## 1. Research question and the result obtained here

The proposed order is native law, admissible cut, geometric representation,
then curvature. A metric or gravitational field equation must not be silently
installed at the first stage. This module establishes a conditional part of
that programme: a specified family of invertible recognition transports can
determine a compatible bilinear form, including its signature, without giving
that form to the reconstruction. A separate geometric adapter identifies its
connection curvature with metric curvature when the required conditions hold.

For an explicitly chosen family assembled from the EMK K block, the only
invariant symmetric forms are multiples of diag(1,-1,...,-1). Four carrier
components therefore give a Lorentz-type 1+3 signature, up to overall sign.
The dimension four is an input, not a derived dimension of spacetime. Choosing
the R block instead gives definite forms. The common native algebra alone
therefore does not select the transport sector or the physical signature.

This is a **conditional mathematical bridge**, not a derivation of a physical
gravity law from an entirely unspecified uncut ground. The unknown gravitational
law, its universal coupling to probes, and its quantum representation remain
research targets. The constructions below make their missing obligations precise.

## 2. Native proposal and explicitly added representation data

**NG0 — proposed pregeometric gravitational carrier.** Investigate a common
law of admissible recognition transport, with recovered middle identity and
Smriti memory, whose representations govern the relative responses of probes.
The operational target is a predicted change in probe readouts under a declared
continuation protocol, including loop return. Calling that law gravitational
requires showing that the same native law couples to the relevant probes and
matches gravitational observations. General continuation, memory loss or a
nonzero loop residue alone does not establish that identification.

A mathematical presentation already makes distinctions. Neither the graph
below nor a vector space is identified with Śūnya or literal uncut ground.
They are candidate representations in which NG0 can be investigated.

Unlike U1–U16's finite sets, this module uses finite-dimensional vector spaces
over the rationals, with real extensions when discussing signature. Their
state sets are not finite. The additional data are:

- typed recognition objects x, each with a declared rank-n response carrier V_x;
- lawful linear invertible maps U_e:V_x → V_y for declared arrows e:x→y;
- composition only at matching recovered objects, as required by S1/S6;
- a specified response target and cuts, with exact intertwining checked;
- an unknown symmetric form H_x to be reconstructed if the transport permits it.

No manifold, metric, physical clock, action functional, quantum state space or
gravitational field equation is an input to U17–U19. Rank, linearity,
invertibility and the choice of arrows **are** new candidate assumptions.
They are not consequences of having named an uncut primitive.

The higher-rank star assembly of overlapping EMK planes is an additional
representation choice. It is not automatically one global copy of the original
two-dimensional grading, and its mixed commutator is not automatically the
source capsule's odd-memory observable.

## 3. Source audit: where this differs from existing papers

The [pinned source audit](GRAVITY_SOURCE_PINS.json) identifies exact files.

| Source | What is already present | What it does not supply for this route |
| --- | --- | --- |
| RSC 24/27/30, S1/S2/S6 | Clock-free transport, cut memory, recovered middle identity | A physically selected gravitational transport law |
| EMK-2, G3; RKF 48, G1 | K²=I, R²=-I, exact rational even-sector composition and cut-graded commutators | A unique spacetime signature or gravitational interpretation |
| RKF 55, G2 | Native dagger-compatible flow generators and EMK loop-residue bindings | Permission to call real K alone a dagger-unitary generator; the source explicitly distinguishes it from iota K |
| EMK-G1, G4 | Curvature of a declared warped metric; separation of recognition and metric curvature | Derivation of that metric from the primitive carrier |
| Thermodynamic propagation section, G5 | Metric reconstruction by polarization once a scalar Lorentzian principal symbol is assumed | Selection of that signature from unstructured native data |
| Thermodynamic action and non-selection, G6/G7 | Field equations in a declared action class; constitutive freedom made explicit | A uniquely selected action, coupling or potential without the extra contract |
| CID-1, G8 | A metric-independent loop obstruction and metric-dependent information decomposition in its declared response algebra | A universal identification of that obstruction with spacetime curvature |

The thermodynamic sources here are the corrected working edition on public
PR #4, not a claim that those changes are already merged. They are comparative
sources only; none of their spacetime or action assumptions is a premise of
U17–U19. All general methods below are elementary linear algebra or standard
connection geometry. The contribution is their explicit assembly, exact checks
and failure controls in this native research programme, not a priority claim
for invariant forms, quotient forms or the Levi-Civita characterization.

## 4. Proposition U17 — Recovering compatible forms from transport

Consider a connected finite graph of rank-n carriers and invertible edge maps.
The compatibility equations are

\[
U_e^T H_y U_e=H_x \quad(e:x\to y).
\]

Choose a root r and a spanning tree. Let T_x:V_r→V_x be its composite transport,
with T_r=I. An inverse used to traverse a tree edge backwards is an algebraic
frame comparison; it is not automatically an executable physical intervention.
For each edge define the rooted loop comparison

\[
L_e=T_y^{-1}U_eT_x.
\]

Then all compatible symmetric form fields are exactly

\[
H_x=T_x^{-T}H_rT_x^{-1},\qquad
L_e^T H_rL_e=H_r\quad\text{for every edge}.
\]

**Proof.** The tree compatibility equations successively force the expression
for H_x. Substituting it in an arbitrary edge equation and multiplying by T_x
and its transpose gives precisely the loop equation. Reversing these steps
constructs a compatible field from any root solution. Congruence by invertible
T_x preserves rank and real signature. Thus a nondegenerate root solution
gives one everywhere; if no such solution exists there is no compatible
nondegenerate metric field for this transport. **Proof complete.**

Solve the root equations as a homogeneous linear system on n(n+1)/2 unknown
symmetric entries. The solver takes the transports and n, not a proposed H.
An invariant solution space can have dimension zero, one or more. A
one-dimensional space generated by a nonsingular matrix determines a form
up to nonzero scale. Overall sign, normalization and physical units still
require a convention or a separate calibration.

Frame changes U_e↦S_yU_eS_x^{-1} send H_x to S_x^{-T}H_xS_x^{-1}; a chosen
coordinate frame does not change existence or signature. A loop L=2I forces
4H=H and hence H=0, showing that compatible nondegenerate forms need not exist.

## 5. Proposition U18 — Conditional signature from an EMK family

Let n≥2. Choose a distinguished component 0 and, for each i=1,...,n−1, embed
the native K block in the (0,i) plane. As a full n×n generator K_i exchanges
those two components and is zero elsewhere. For rational t_i≠0,±1 define

\[
B_i=(I+t_iK_i)(I-t_iK_i)^{-1}.
\]

The (0,i) block is

\[
\frac1{1-t_i^2}
\begin{pmatrix}1+t_i^2&2t_i\\2t_i&1+t_i^2\end{pmatrix},
\]

and the complement is identity. These are real invertible EMK-algebra
candidate maps. They are **not** claimed to be the native dagger-unitary flows
of RKF 55; that source requires its own generator typing.

The common invariant symmetric forms are exactly

\[
\boxed{H=\lambda\operatorname{diag}(1,-1,\ldots,-1).}
\]

In particular, every nonzero solution is nondegenerate and indefinite.

**Proof.** For any G with I−tG invertible and t≠0, Cayley invariance is equivalent
to the generator equation:

\[
C_t(G)^THC_t(G)=H
\iff (I+tG)^TH(I+tG)=(I-tG)^TH(I-tG)
\iff G^TH+HG=0.
\]

For G=K_i the entries imply H_0i=0, H_ii=−H_00, and H_ij=H_0j=0 whenever
j is outside the selected pair. Imposing every i leaves exactly the stated
diagonal form. Direct substitution proves sufficiency. **Proof complete.**

This is stronger than putting a Lorentzian form into the solver and checking
that the chosen matrices preserve it: the equations determine the full
solution space. But the family itself encodes a choice. Its distinguished
component, coupling pattern and K-sector selection were declared, not derived.

### Alternative sector and the positive-unitary obstruction

Embed the native R block instead, with R_i e_0=e_i and R_i e_i=−e_0, and use
the same nonzero rational Cayley parameters. The generator equation is now
HR_i−R_iH=0. Entry comparison gives H_ii=H_00 and all off-diagonal entries zero.
The common forms are therefore exactly λI. Requiring both complete K and R
families on the same carrier and same form forces H=0.

More generally, if a transport family preserves a positive-definite form E,
then E is already in its invariant symmetric solution space. If that space
is one-dimensional, its nonsingular members are definite, up to sign; they
cannot consist only of Lorentz-type forms. Thus a positive-unitary recognition
carrier must not be silently equated with the sought spacetime metric carrier.
A different representation or an additional identification is required.

These facts do not identify the R family with quantum physics or the K family
with gravity. They show exactly why a sector-selection theorem is needed.

**Controls W10–W12.** The full four-component K family selects one form line;
the R family selects a different one; the combined family has none. Keeping
only the (0,1) K pair in dimension four leaves four independent invariant
symmetric forms, so missing couplings destroy uniqueness. A mixed K-family
commutator loop in dimension three is nonidentity while preserving H: metric
compatibility does not force every loop to close.

## 6. Proposition U19 — What a cut may carry without inventing geometry

Let C:V→W be a surjective linear cut and H a symmetric bilinear form on V.
There is a unique h on W satisfying

\[
H=C^T h C
\]

if and only if ker C is contained in rad H={v:H(v,z)=0 for every z}. The
descended h is nondegenerate exactly when ker C=rad H.

**Proof.** If H pulls back from h, Ck=0 implies H(k,z)=0. Conversely set
h(Cv,Cw)=H(v,w). Changing either representative by ker C does not change the
value, so h is well-defined and unique by surjectivity. Its radical is the
image of rad H, which is zero exactly when rad H=ker C. **Proof complete.**

Consequently, a noninjective cut cannot carry an entire nondegenerate H by
exact pullback. A target-restricted form or a repaired observation can still
be useful; that is a different declared target, not an exception to the theorem.

### W16 — a lawful lossy cut with a nondegenerate geometric target

There is a positive example as well as the obstruction. Extend a candidate
n-component response by one memory coordinate m. For the K-family B_i of U18
declare invertible lifts and a separate memory update

\[
\widetilde B_i=\begin{pmatrix}B_i&0\\r_i&2\end{pmatrix},\qquad
M=\operatorname{diag}(I_n,2),\qquad C=(I_n\;0),
\]

where each r_i is a chosen row. Invariance under M forces all memory-row and
memory-column entries of a symmetric source form to vanish. Invariance under
the lifted B_i then reduces exactly to U18, so the common source forms are
λ diag(1,-1,...,-1,0). For λ≠0 their radical is precisely ker C.
Consequently C discards a real memory distinction while the descended form
is nondegenerate and C tilde(B_i)=B_i C, CM=C. Thus metric-target faithfulness
does not require recovering the whole presented state.

The n=2 example uses t=1/2 and r=(1,0), and is checked exactly. Its quotient
form is the internal H used in W15's geometric adapter. This links a genuine
lossy cut to the curvature example without claiming that the three-component
source is literal uncut ground. The memory update, target and lifted family
remain candidate inputs; another target sensitive to m would still be lost.

A transport U descends exactly when ker C is U-invariant, equivalently when
there is V with VC=CU. A right inverse S of C gives the candidate V=CUS, but
the full equation VC=CU must still be checked. If A and B commute and both
descend, their observed maps commute: apply their commutator to the surjective
C and use the intertwining equations.

**W13 — false curvature from inadmissible compression.** Let

\[
A=\begin{pmatrix}0&0&1\\1&0&0\\0&1&0\end{pmatrix},\quad B=A^2,\quad
C=\begin{pmatrix}1&0&0\\0&1&0\end{pmatrix},\quad S=C^T.
\]

Then AB=BA=I, but

\[
[CAS,CBS]=\operatorname{diag}(-1,1)\ne0.
\]

Both exact descent tests fail. This apparent noncommutativity is caused by
discarding and then resetting memory; it is not evidence that a faithfully
represented commuting native pair has acquired gravitational curvature. This
is the linear version of the admissibility discipline in U13, and is also
consistent with the inherited RSC cut-corner composition residue U7.

## 7. Proposition U20 — When transport curvature becomes metric curvature

U17–U19 do not construct a smooth spacetime. For this next, explicitly
**post-cut representation**, supply:

1. an n-dimensional smooth local chart M and a rank-n response bundle;
2. a constant nondegenerate internal form H in a chosen frame, and a smooth
   connection A satisfying A^TH+HA=0;
3. an invertible bundle map e:TM→V, represented by a coframe, with
   de+A∧e=0.

The coframe is the missing type bridge: it identifies chart tangent responses
with the internal response carrier. It is additional data, not a notation
that follows from solving the invariant-form equations.

Define

\[
g=e^THe,\qquad \nabla_Xv=e^{-1}\bigl(d_X(ev)+A(X)ev\bigr),\qquad
F=dA+A\wedge A.
\]

Then ∇ is the Levi-Civita connection of g and

\[
\boxed{R^g(X,Y)=e^{-1}F(X,Y)e.}
\]

**Proof.** Differentiating g(v,w)=H(ev,ew) and using A^TH+HA=0 proves metric
compatibility. Directly computing ∇_X Y−∇_Y X−[X,Y] gives
e^{-1}(de+A∧e)(X,Y), so the declared torsion condition makes it zero.
The unique metric-compatible torsion-free connection is Levi-Civita: the
difference of any two such connections is symmetric in its two vector inputs
by torsion freedom, while lowering its output with g makes it antisymmetric
in the last two slots by metric compatibility; cycling these identities
forces the difference to vanish. Finally applying the commutator of the
covariant derivatives to ev gives (dA+A∧A)(ev), yielding the displayed
curvature relation. **Proof complete.**

This is the standard geometric identification with every type bridge exposed.
It is not a new gravity field equation. A finite loop matrix does not by itself
prove local curvature: a smooth connection and infinitesimal-loop identification
are still required; global monodromy can persist with zero local curvature.

### W14 — metric compatibility alone is insufficient

Take H=diag(1,-1), e=(du,dv), and A=uK dv on a two-dimensional chart.
The connection preserves H and F=K du∧dv is nonzero. But e gives the constant
flat metric du²−dv², whose Levi-Civita curvature is zero. Here
de+A∧e=(0,-u du∧dv), so the torsion gate fails. This supplies a concrete reason
for the recognition/metric-curvature separation already present in EMK-G1.

### W15 — a valid curved adapter, computed by two routes

On any patch where f(u)≠0 set

\[
e^0=du,\quad e^1=f(u)dv,\quad A=f'(u)K\,dv,\quad H=\operatorname{diag}(1,-1).
\]

Directly, de+A∧e=0 and A preserves H. Thus

\[
g=du^2-f(u)^2dv^2,\qquad F=f''(u)K\,du\wedge dv,
\]

and in the coordinate frame

\[
R^g_{uv}=\begin{pmatrix}0&f f''\\ f''/f&0\end{pmatrix}.
\]

Computing from the metric independently gives Γ^u_vv=f f' and
Γ^v_uv=Γ^v_vu=f'/f and exactly the same curvature. Taking f=1+u² produces
nonzero curvature; f=1 produces zero curvature with the same internal H.
The entire calculation uses exact rational arithmetic at rational chart
points for polynomial f. The written differential identities hold on the
stated smooth patches, not merely at the sampled points.

The warp f is a **supplied adapter**. The native algebra and reconstructed H
have not selected it. Both choices obey the bridge, so the bridge does not
determine the physical gravitational field. This identifies the next missing
selection problem rather than hiding it in a metric ansatz.

## v0.6 native-grading clarification

[U24](NATIVE_GRADING_AND_SEAM_MEMORY.md) now proves that overlapping
zero-extended K/R pairs cannot all retain RKF 48's local grading under one
common involution on a connected carrier with at least three components.
The invariant-form statements in this module remain valid for their declared
candidate transports. Their metric-sign partition is not the native EMK cut.
The larger grading-preserving assembly and its exact readout-memory costs are
given in U25/U26; a native inter-edge interaction remains to be constructed.

## 8. Evidence, claim boundary and next research target

`gravity_model.py` reconstructs invariant forms by rational row reduction,
checks exact cuts and handles rational EMK candidate maps. The certificate
compares these with direct form invariance, generic Cayley inversion,
changes of frame, typed graph edges and independent metric curvature.
Coverage includes 24 K-family reconstructions in dimensions two through five,
20 R-sector/combinations comparisons, 100 two-route Cayley checks, 12 changes
of frame, 48 integer single-arrow cases with 1,296 candidate-form comparisons,
27 graph-edge checks, eight cut-form controls, one lawful memory-cut bridge
and 24 supplied warped-chart
points. There are 50 regression tests across the full package.
The exact finite coverage is recorded in CERTIFICATE.json; none of these
grids is an exhaustive check over real geometries. General proofs above are
written arguments and are not proof-assistant verification.

Completed here: a conditional metric-signature reconstruction, an exact loss
gate for geometric targets, and a conditional identification of connection and
metric curvature, with valid and invalid adapters. This goes upstream of the
existing principal-symbol construction by solving for H from a chosen algebraic
transport family, but it does not subsume the thermodynamic analytic results.

The next substantive target is to derive a **transport-sector and adapter
selection law** from declared native seam data. It must explain why the K
family, its rank and coupling pattern, and a particular coframe/connection
are selected for gravitational probes. It must predict a response that can
reject the candidate, keep lawful quantum response on an appropriately typed
carrier, and establish any shared quantum/geometric observational target.
No Einstein equation, Born rule, spacetime continuum, matter coupling,
quantum-gravity closure, RH or Yang–Mills closure is claimed in this module.
