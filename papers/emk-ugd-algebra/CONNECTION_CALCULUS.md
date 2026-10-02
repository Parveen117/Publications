# EMK-C1 — Native connection and curvature calculus

Author: Monty Dabas. Development edition: 2 October 2026.

This chapter continues the existing EMK/UGD algebra and Morphic derivation programme. Its mathematical subject is **noncommutative differential calculus on an admitted native algebra**, with a possibly noncommuting direction frame. General proofs, exact finite execution and physical realization have separate scopes. The methods are standard Lie-algebra/connection identities; no external priority is claimed.

## Sources and scope

The native scalar/algebra contract is RKF's [N01](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md). [EMK-1/EMK-2](README.md) supply the relations and cut grading; the [Morphic Algebra A-2 source](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/sources/rkf_reference/theorum/morphic_algebra/main.tex) supplies the derivation programme. [T55](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/55_generalized_euler_emk_dock_theorem.md) binds generalized Euler transport to EMK generators and distinguishes cut grading from dagger grading.

We work in a specified associative, unital algebra over the established characteristic-zero cut field. This is a concrete algebraic sector of the broader typed Morphic programme. It does not prove every historical overlay, category, spectral or universality axiom. The executable sector uses exact rational cut scalars and the unchanged canonical RKF symbolic compiler.

## C1. Direction frame and derived EMK derivatives

Let \(\mathcal A\) be that algebra, with derivations \(\delta_i\) satisfying
\[
[\delta_i,\delta_j]=\sum_l c_{ij}^{\,l}\delta_l.
\]
The supplied constants \(c_{ij}^{\,l}\) are central, antisymmetric and obey the Lie Jacobi identity. They are constant under every \(\delta_i\). This is a declared direction contract, without a physical clock or spacetime coordinate.

For the native EMK sector derive
\[
G_1=R,\quad G_2=K,\quad G_3=RK,\qquad \delta_i(a)=[G_i,a].
\]
Associativity gives
\[
[G_i,ab]=[G_i,a]b+a[G_i,b],\quad
[\delta_i,\delta_j]=\operatorname{ad}_{[G_i,G_j]}.
\]
The primitive relations give exactly
\[
[G_1,G_2]=2G_3,\quad [G_2,G_3]=-2G_1,\quad [G_3,G_1]=2G_2.
\]
Thus the derivatives and their frame constants are derived in this admitted EMK algebra. Three directions here label algebraic generators; they do not select three spatial dimensions.

**Proof.** Expand the two commutators and cancel the middle products. The derivative-commutator identity follows by the same expansion. The three displayed brackets follow from \(R^2=-I\), \(K^2=I\), \(KR=-RK\). Their Jacobi identities follow directly; for the three distinct generators each double bracket vanishes. Repeated-index cases cancel by antisymmetry. ∎

## C2. Connection and its correctly typed curvature

On the regular right module \(\mathcal A\), choose connection coefficients \(A_i\in\mathcal A\) and define
\[
\nabla_i x=\delta_i x+A_i x.
\]
This is a right-module connection: Leibniz gives
\(\nabla_i(xa)=(\nabla_i x)a+x\delta_i a\).
Left multiplication by the coefficients is therefore compatible with the stated module convention; no commutation of \(A_i\) with \(a\) is assumed.
These coefficients are supplied connection data. Their physical selection is an additional problem. The derived connection curvature is
\[
\mathcal F_{ij}:=[\nabla_i,\nabla_j]-\sum_l c_{ij}^{\,l}\nabla_l
                =L_{F_{ij}},
\]
where \(L_a(x)=ax\) and
\[
\boxed{F_{ij}=\delta_i A_j-\delta_j A_i+[A_i,A_j]
                    -\sum_l c_{ij}^{\,l}A_l.}
\]

**Proof.** Apply the commutator to \(x\), use Leibniz, and cancel \(A_j\delta_i x\) and \(A_i\delta_j x\). The remaining derivative of \(x\) is \([\delta_i,\delta_j]x\); the frame subtraction cancels it. What remains is the displayed left multiplication. Since \(L_a(1)=a\), this regular-module curvature vanishes exactly when \(F_{ij}=0\). ∎

For the inner EMK frame, put \(H_i=G_i+A_i\). The same formula becomes
\[
F_{ij}=[H_i,H_j]-\sum_l c_{ij}^{\,l}H_l,\qquad
\nabla_i x=H_i x-xG_i.
\]
This is a finite exact expression in native words.

The raw order defect \([\delta_i,\delta_j]\) and connection curvature relative to the declared frame are different targets. At \(A_i=0\), the EMK frame has nonzero order defect while \(F_{ij}=0\). This does not erase an earlier native commutator target or prove that observation makes it flat. It identifies precisely what is being measured by this connection calculus.

## C3. Gauge transport and curvature covariance

For any two-sided invertible \(u\in\mathcal A\), define
\[
A'_i=u^{-1}A_i u+u^{-1}\delta_i(u).
\]
Then
\[
\nabla'_i=L_{u^{-1}}\nabla_i L_u,\qquad F'_{ij}=u^{-1}F_{ij}u.
\]
In particular, a pure gauge \(A'_i=u^{-1}\delta_i(u)\) obtained from \(A_i=0\) is flat for this target.

**Proof.** Leibniz applied to \(\delta_i(ux)\) proves the first identity. Conjugating the curvature-operator definition proves the second: the central frame constants commute with the conjugation maps. Evaluate on \(1\) to recover the coefficient identity. Flatness is preserved because conjugation is invertible. ∎

For the inner EMK frame \(H'_i=u^{-1}H_i u\). An arbitrary conjugation of \(A_i\) alone omits the derivative term and generally fails this covariance contract.

## C4. Bianchi identity with a noncommuting frame

Define the covariant derivative of a coefficient by
\[
D_i b=\delta_i b+[A_i,b].
\]
The exact Bianchi identity is
\[
\boxed{\sum_{\mathrm{cyc}(i,j,k)}
 \left(D_iF_{jk}+\sum_l c_{jk}^{\,l}F_{il}\right)=0.}
\]

**Proof.** Operator commutators obey Jacobi by associativity. Substitute
\([\nabla_j,\nabla_k]=\sum_l c_{jk}^{\,l}\nabla_l+L_{F_{jk}}\)
into its cyclic Jacobi sum. Leibniz gives
\([\nabla_i,L_b]=L_{D_i b}\).
The terms containing \(\nabla_m\) vanish by the Jacobi identity of the supplied frame constants. The surviving left multiplication is the displayed cyclic sum, and evaluation on \(1\) proves the result. ∎

For \(G=(R,K,RK)\), the frame contribution happens to vanish for the three distinct directions because it involves \(F_{ii}=0\). This simplification is specific to that frame. It is not permission to omit the frame term generally.

An exact counterexample uses the same native algebra with
\(G=(RK,K+R,I)\), whose only nonzero independent bracket is
\([G_1,G_2]=2G_2\).
For \(A=(R,K,R)\),
\[
\sum_{\mathrm{cyc}}D_iF_{jk}=-8RK,\qquad
\sum_{\mathrm{cyc},l}c_{jk}^{\,l}F_{il}=8RK.
\]
The full identity closes; the incomplete one does not. The inner action of the central third frame direction is zero, which is allowed: this example does not assert a faithful spatial-frame representation.

## C5. Mixed channels are part of curvature

If \(A_i=\sum_\alpha A_i^\alpha\), then
\[
F_{ij}=\sum_\alpha\left(\delta_i A_j^\alpha-\delta_j A_i^\alpha
+[A_i^\alpha,A_j^\alpha]-\sum_l c_{ij}^{\,l}A_l^\alpha\right)
+\sum_{\alpha\ne\beta}[A_i^\alpha,A_j^\beta].
\]

**Proof.** Distribute the derivatives, frame sum and commutator. Partition the double commutator sum into equal and unequal channel labels. ∎

This extends the existing EMK-T1 mixed-channel identity to a specified derivation/connection contract. It does not derive the master tensor's eight-channel declarations, tolerances or physical couplings. Dropping mixed channels can change the curvature verdict even when each separately inspected channel is flat.

## C6. Observation of curvature requires descent and faithfulness

Let \(E:V\to Y\) be a supplied linear observer on an admitted module with connection operators \(\nabla_i\). Descent exists exactly when
\[
\nabla_i(\ker E)\subseteq\ker E
\]
for each \(i\). On \(\operatorname{ran}E\) there is then a unique \(\bar\nabla_i\) with \(E\nabla_i=\bar\nabla_i E\), and
\[
E\mathcal F_{ij}=\bar{\mathcal F}_{ij}E.
\]

**Proof.** Define \(\bar\nabla_i(Ex)=E\nabla_i x\). This is well defined exactly under the displayed kernel condition. Compose the intertwining identity twice and subtract the central frame term to obtain curvature descent. ∎

Vanishing observed curvature is not automatically vanishing native curvature. For this particular target the exact additional requirement is
\[
\ker E\cap \operatorname{ran}\mathcal F_{ij}=\{0\}
\]
for every pair whose flatness is being tested. Under that hypothesis \(E\mathcal F_{ij}=0\) implies \(\mathcal F_{ij}=0\). The converse is immediate. This is the existing observer-faithfulness mechanism, specialized here; its minimum repair is consumed from RKF N03, rather than reimplemented.

On the regular right module \(\mathcal A\), the left-curvature formula above applies. For another supplied module, its representation and possible kernel must be stated separately. A state-observer kernel is not automatically a two-sided algebra ideal. An exact two-copy regular-module test observes a flat second copy while the first copy retains nonzero curvature; descent alone therefore cannot certify native flatness.

## C7. Classical comparison and next native obligations

Under a later admitted smooth-coordinate representation with commuting derivatives \(\delta_i=\partial_i\), the constants are zero and the formula reduces to
\[
F_{ij}=\partial_i A_j-\partial_j A_i+[A_i,A_j].
\]
Thus ordinary gauge curvature is a special representation of the stated algebraic contract. Recovering Riemann curvature additionally needs a tangent-module interpretation and a metric-compatible, torsion-free connection; those structures are not selected by this chapter.

The next native foundation targets are: connections on typed multi-object modules; admissible tensor products and covariant differentiation; a separately proved metric/dual pairing; and bounded or unbounded completion/domain control for flows. Existing UGD phase/scale/seam digits and carry ledgers remain richer data than a single scalar cut-field coordinate. None is replaced by this finite EMK example.

Physical source selection, spacetime metric, noise preparation, clock/length calibration and field equations belong to downstream adapters in [extra-ideas](https://github.com/Parveen117/extra-ideas). This chapter establishes identities those adapters must preserve.

## Reproduction and evidence

Run from Publications with the pinned RKF checkout supplied separately:

~~~bash
node papers/emk-ugd-algebra/certificates/connection_calculus.cjs \
  --rkf-root ../Recognition-Kernel-Framework --check
~~~

The [source pins](certificates/EMKC1_SOURCE_PINS.json) bind this chapter, the new adapter and the canonical compiler files. [EMKC1_RESULT.json](certificates/EMKC1_RESULT.json) contains exact finite native-word checks, replayable proof witnesses and negative controls. The general statements above have written proofs; finite execution checks their admitted EMK instances and implementation. Formal proof-assistant verification, full Morphic-manuscript certification and empirical physics are not claimed.
