# RCG-1 — Relational cut space and its transport defect

## Question and declared construction

The proposed intuition is that measurement belongs to an observation cut and need not close within its retained sector, while something shared between measurements may constitute space. CLD-1 establishes the first part. The second part needs a definition and a preservation law: nonclosure alone does not logically imply a conserved space.

We define a **finite operational cut geometry** from relationships among a labelled family of cuts brought to a common comparison carrier. Its dimensionless distances and cycle readings survive common reversible transport. We give a constructive test for common transport, and exact examples where distance-only records or actual measurement channels lose geometry. This supplies a candidate relational meaning of space; it does not derive physical spacetime, its dimension, or gravitational field equations.

The finite positive matrix representation, dagger, trace and unitary arrows are the existing native representation contract. They are stated inputs, not a claim that a physical Hilbert space or spatial metric has been derived from measurement alone. All operation parameters below are uncalibrated.

Lineage: reuse [CRC-1](../native-cut-record-completeness/THEOREM.md)'s operational quotient and record-sufficiency criterion; [QG-R1](../native-relative-curvature/THEOREM.md)'s common-base comparison and loop information; [CID-D1](../native-loop-concentration/CUT_INFORMATION_DYNAMICS.md)'s cut metric and curvature tensor; and [CLD-1](../native-cut-local-dynamics/THEOREM.md)'s moving-cut defect. No new priority is claimed for projective-unitary invariants. The finite rank-one reconstruction below is an explicit specialization of the known Gram/cycle method of Chien and Waldron, *A characterisation of projective unitary equivalence of finite frames*, [arXiv:1312.5393](https://arxiv.org/abs/1312.5393). We give the argument here instead of importing an unproved reconstruction step. The contribution to this research chain is the common-transport residual, its relation to CLD, the falsification fixtures, and the resulting observation-space contract.

## 1. Comparisons require transport

Let `P_j` be rank-r orthogonal cuts (`r>0`) on finite equal-dimensional carriers, and let declared unitary links `L_j` bring them to one base carrier. Set

\[
 \widehat P_j=L_jP_jL_j^*,\qquad
 \kappa_{ij}=\frac1r\operatorname{Tr}(\widehat P_i\widehat P_j),\qquad
 d_{ij}^2=1-\kappa_{ij}
 =\frac1{2r}\|\widehat P_i-\widehat P_j\|_{HS}^2. \tag{1}
\]

This follows by expanding the square and using `Tr(P_i)=r`. Thus `0<=kappa<=1`; `d` is a pseudometric on labels and a metric after identifying identical based projectors. Triangle inequality is inherited from the norm, not from the squared distance. This normalization requires equal ranks.

Independent frame changes have
`P_j -> G_j P_j G_j^*`, `L_j -> G_0 L_j G_j^*`. All based cuts then conjugate by the same `G_0`. Hence (1) and all ordered cycle readings

\[
 B_{j_1\ldots j_m}=\frac1r\operatorname{Tr}
 (\widehat P_{j_1}\cdots\widehat P_{j_m}) \tag{2}
\]

are invariant. A local frame change without transforming its comparison link changes the experiment.

Different comparison paths are not automatically equivalent. If `L_j'=h L_j` with loop transport `h`, the based cut is unchanged exactly when `[h,widehat P_j]=0`. Otherwise the path is part of the record. This reuses QG-R1's based-loop requirement; the measurement labels alone do not supply a path-independent space.

For a smooth family `P(x)` compared along a reference connection, write `X_a=D_a P`. The small-displacement expansion of (1) gives

\[
 d^2(P(x),P(x+dx))
 =\frac1{2r}\operatorname{Tr}(X_aX_b)dx^a dx^b+o(|dx|^2)
 =\frac1r g^{cut}_{ab}dx^a dx^b+o(|dx|^2). \tag{3}
\]

The second equality uses the offdiagonal cut identity `PX_aP=QX_aQ=0`: the two diagonal blocks of the symmetrized product have equal trace. This is exactly CID-D1's already established sensitivity metric, normalized by rank, with the second cut parallel-compared to the first. Its null directions are directions invisible to the cut. A quotient manifold needs regularity/constant rank beyond this local expansion. Coordinates `x` here are controls, not presumed spatial coordinates; (3) does not replace CID-D1's antisymmetric curvature information.

## 2. Distances do not in general determine the oriented cut geometry

In a two-mode carrier take normalized columns

\[
 v_0=(1,0)^T,\quad v_1=(3/5,4/5)^T,\quad
 v_2^{\pm}=(3/5,\ \pm4i/5)^T,\qquad P_j=v_jv_j^*.
\]

The labelled triples with either sign have the same pairwise overlap record,

\[
 \kappa_{01}=\kappa_{02}=9/25,\qquad \kappa_{12}=337/625,
\]

but

\[
 B_{012}^{\pm}=\frac{81\pm144i}{625}. \tag{4}
\]

A common unitary preserves the full complex cyclic trace, so no common unitary maps one labelled triple to the other. Complex conjugation does map them, but that antiunitary operation is outside the declared unitary equivalence. Distance-only geometry cannot distinguish these two orientations.

This also specifies a measurement limitation: sequential rank-one projection success probabilities contain products of overlap magnitudes and do not by themselves expose the complex phase in (4). A phase-sensitive comparison/interference protocol is additional access. We do not claim ordinary cut outcome counts already measure every complex cycle invariant.

## 3. A finite complete record in the rank-one case

For normalized rank-one cuts, let `G_ij=<v_i,v_j>`. Form the nonzero-overlap graph: edge `ij` exists iff `G_ij!=0`. For a connected graph choose a rooted spanning tree. The following data determine the labelled cut family up to common unitary on its span:

1. all pairwise magnitudes `|G_ij|^2`, including the zero pattern;
2. for each non-tree edge, the oriented cycle product of Gram entries along that edge and the unique return path in the tree.

**Proof.** Choose phases of the `v_j` recursively from the root so every tree-edge inner product is positive real. Their values are the positive square roots of the known magnitudes. Each fundamental cycle then has just one unknown non-tree-edge Gram entry (or its conjugate, according to orientation). Divide by the nonzero known tree product to recover it. Nonedges are zero. Thus the full Gram matrix in this gauge is reconstructed. If two column families have the same Gram matrix, the map `sum a_j v_j -> sum a_j w_j` is well defined and isometric, since squared norms are the same Gram quadratic form. It is onto their spans and extends to a unitary on equal-dimensional ambient carriers. Conversely all listed data are unchanged by phases and common unitaries. This proves completeness.

Disconnected components have mutually orthogonal spans; perform the construction on each component. The Gram rank fixes the occupied span dimension; it says nothing about unused ambient dimensions unless they were specified. For a complete graph rooted at label 0, the fundamental cycle record is just the triples `B_0ij`. For sparse graphs longer cycles may be essential. This theorem is rank-one: no assertion is made that scalar triples classify arbitrary higher-rank projector families.

On an interval of C1 cuts with constant listed data, there exists a C1 common unitary transporting the entire family. To see this, choose C1 normalized representatives locally, perform the tree gauge above, and obtain a constant Gram matrix. Select a fixed independent column subfamily `V_B(s)` and set `E(s)=V_B(s) G_B^{-1/2}`, where the positive inverse square root is constant. Then `E^*E=I`. With `A=E^*E'`,

\[
 K=E'E^*-EE'^*-EAE^* \tag{5}
\]

is anti-Hermitian and `KE=E'`. Solving `U'=KU` extends the span isometry to the full carrier. Every other column has constant coefficients in this basis because its Gram pairings are constant. Hence `P_j(s)=U(s)P_j(0)U(s)^*`. The local constructions extend along an interval; endpoint equality of invariants alone does not identify the intervening transport history.

This is a reconstruction theorem for relational geometry, not a derivation of a physical source generator. Different extensions on the unused complement are invisible to the cuts.

## 4. Common-transport reconstruction and an exact residual

Now allow arbitrary equal positive ranks. At a given operation parameter let `X_j=P_j'` be Hermitian cut tangents. Define the real linear map from anti-Hermitian generators

\[
 \mathcal L(K)=([K,P_1],\ldots,[K,P_m]),
\]

and the squared best-fit residual

\[
 \boxed{\varepsilon^2=\min_{K^*=-K}
 \sum_j\|X_j-[K,P_j]\|_{HS}^2.} \tag{6}
\]

The domain is the full anti-Hermitian matrix algebra for this theorem. Restricting it to a smaller admitted generator algebra answers a different operational question and can increase the residual.

With real Hilbert–Schmidt inner products, the adjoint is
`L^*(Y)=sum_j[Y_j,P_j]`. Therefore the normal equation is

\[
 \boxed{\sum_j[[K,P_j],P_j]=\sum_j[X_j,P_j].} \tag{7}
\]

Proof: expand `Re Tr([K,P_j]^*Y_j)` and cycle the trace. The left operator is `L^*L`, positive semidefinite, with kernel

\[
 \ker\mathcal L=\{K:K^*=-K,\ [K,P_j]=0\text{ for every }j\}. \tag{8}
\]

Finite-dimensional orthogonal projection onto `range L` gives a minimizer and (7). Minimizers differ exactly by (8), so all have the same residual tuple `R_j=X_j-[K,P_j]`. Thus `epsilon=0` iff one common infinitesimal transport explains all observed cut tangents. On an interval, a continuous solution `K(s)` of the zero-residual equations integrates to common transport; pointwise solvability alone is not silently promoted to a global regularity theorem. Section 3 gives a constructive regular lift under its rank-one hypotheses.

Each tangent separately can always be written `[K_j,P_j]` using the inherited minimal `K_j=[X_j,P_j]`. The nontrivial condition in (6) is **one and the same generator for the entire family**.

### Direct connection to CLD

For a specified physical source generator `H`, CLD defines
`D_j=P_j'+(i/eta)[H,P_j]`. Substitute `X_j=D_j` in (6). A zero residual means the cuts share an additional relative motion `K`, so their total generator is `K-iH/eta`. It does not mean `D_j=0` individually. In particular, all cuts can rotate rigidly relative to a stationary source, preserving their mutual geometry while their measured source weights change. This separates invariant relational space from conservation of individual measurement readings.

### Geometry-change budget

For raw derivatives `X_j=P_j'`, substitute `X_j=[K,P_j]+R_j` into the derivative of (1). Commutator terms cancel under the trace, giving

\[
 \boxed{\kappa_{ij}'=\frac1r\operatorname{Tr}(R_iP_j+P_iR_j),}
 \qquad
 |\kappa_{ij}'|\le\frac{\|R_i\|_{HS}+\|R_j\|_{HS}}{\sqrt r}
 \le\sqrt{2/r}\,\varepsilon\quad(i\ne j). \tag{9}
\]

The same cancellation holds when fitting CLD's `D_j`, since the known common `H` is also a commutator. For an ordered cycle,

\[
 B_{j_1\ldots j_m}'=\frac1r\sum_{a=1}^m
 \operatorname{Tr}(P_{j_1}\cdots R_{j_a}\cdots P_{j_m}). \tag{10}
\]

Thus the residual is an exact source of relational change, measured independently of a common moving frame. For a parameter-dependent frame `G(s)`, its anti-Hermitian `G'G^*` is absorbed in the minimizing generator, while each residual conjugates by `G`. Consequently `epsilon` is frame invariant for the full admitted anti-Hermitian algebra. If that algebra is restricted, transform the admitted class consistently as well.

## 5. An infinitesimal distance-blind deformation

Use the real columns `v_0=(1,0)`, `v_1=(3/5,4/5)`, `v_2=(5/13,12/13)`. Hold the first two cuts fixed and change only the phase of the second component of `v_2`:

\[
 v_2(s)=(5/13,(12/13)e^{is}),\qquad
 X_2(0)=\begin{pmatrix}0&-60i/169\\60i/169&0\end{pmatrix}.
\]

All pairwise overlap derivatives vanish at zero: real symmetric projectors have zero trace pairing with this imaginary offdiagonal tangent. Nevertheless

\[
 \left.\frac{d}{ds}B_{012}\right|_0=\frac{720i}{4225}\ne0. \tag{11}
\]

No common generator produces these three tangents. A generator fixing both `P_0` and the nonorthogonal distinct `P_1` in two dimensions is a scalar multiple of `iI`, which cannot move `P_2`. Therefore (6) has strictly positive residual. The native certificate obtains `epsilon^2=3600/25621` and verifies (7), residual orthogonality and (9)–(10). This value uses the full two-mode anti-Hermitian algebra, with its scalar `iI` kernel fixed to zero, and the unweighted sum in (6).

This is a local failure of a distance-derivative-only test, in addition to Section 2's global mirror ambiguity. It shows why a loop record matters when evolving the candidate space.

## 6. Actual measurement need not preserve the geometry

For a projective nonselective measurement with orthogonal outcomes `{Q_a}`, the effect/state pinching map in this self-adjoint case is
`E(A)=sum_a Q_a A Q_a`. It is an orthogonal projection in Hilbert–Schmidt operator space. For `Delta=P_i-P_j`,

\[
 \|\Delta\|_{HS}^2
 =\|E(\Delta)\|_{HS}^2+\|(I-E)\Delta\|_{HS}^2. \tag{12}
\]

The proof is that diagonal and offdiagonal blocks are orthogonal. Equality of before/after distances holds exactly when the difference already lies in the measurement's block-diagonal algebra. `E(P_j)` need not be a projector, so the after-distance here is the same normalized Hilbert–Schmidt comparison of the resulting operators, not a claim that the map remains inside the pure-cut manifold.

For `P_1,P_2^+` from Section 2, pinching in the `P_0,I-P_0` basis makes both outputs `diag(9/25,16/25)`. Their before squared distance `288/625` becomes zero. Measurement openness or cut dependence therefore cannot imply that every measurement channel preserves relational geometry. The invariant claim has a specified reversible transport domain; information-erasing interventions are a different operation class.

## 7. Research result and the next gravitational obligation

The defensible replacement for “space is automatically invariant because measurement is not closed” is:

> A declared family of observation cuts has a dimensionless relational geometry. Common reversible transport preserves it; pair distances plus sufficient cycle data reconstruct it in the finite rank-one case. A nonzero common-transport residual detects deformation of that geometry, while a measurement channel can erase parts of it.

The metric and antisymmetric tensor already existed in CID-D1. We have supplied a family-level reconstruction/transport criterion, identified distance-blind changes, and separated common frame motion, individual source readings and measurement disturbance. No external spacetime coordinate is needed for these finite statements.

Next, an actual gravitational extension must specify which operational directions become spatial, choose a local interacting action and derive its constraints. The cut relations do not yet select spatial dimension, Lorentzian causal structure, physical length scale or the source law for (6). Existing AGD/GCP matter–metric coupling remains a selected model; the residual `g R N f'` in its naive spatial constraint algebra is still open. This packet changes the observation contract available for that construction, not the status of the physical quantum-gravity completion claim.

## Certification scope

`verify.cjs` uses the unchanged native rational and rewrite engines. It reconstructs an exact common generator, verifies a nonzero least-squares residual, distinguishes mirrored cut triples and a distance-blind tangent, checks independent frame/link covariance, path dependence and pinching loss, and replays algebraic witnesses. Python 3.12 performs I/O only. General rank-one reconstruction, smooth lifting and norm inequalities have written proofs; they are not proof-assistant-certified by the finite fixtures. Prior source packets are frozen and pinned.
