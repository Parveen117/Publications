# IRT-1 — Dimensionless invariant ratios and their transport

15 October 2026. Continues AGD-1 and the operational record programme.

**Result.** Invariance is always stated relative to specified transport; an invariant is not classified here as stable or unstable. Equal relative rates preserve a ratio. A 45-degree line represents this equality in a specifically calibrated logarithmic chart, and constant-rate exponential transport additionally requires a semigroup law. For the actual AGD quantum geometry, a reflection symmetry supplies exactly conserved dimensionless sector-weight ratios even with the spatially varying matter coupling. We classify the bounded geometry observables commuting separately with the inherited geometry generator and lapse: they are precisely the span of identity and reflection. This is a complete record for that stated invariant algebra, not a complete quantum state or gravitational theory.

## 1. Contract, reuse and the word invariant

Use the unchanged native dagger algebra, positive pairing and trace. Quantum states/readings use QGC-1's explicit statistical adapter. Infinite circle geometry and the joint self-adjoint operator are those of AGD-1. NAC-1 supplies the calibrated-angle convention; CRC-1 supplies operation-dependent record completeness. These premises and earlier proofs are reused.

For a positive source `rho`, Hermitian reading `A`, positive nonzero reference reading `B`, and `Tr(B rho)>0`, define

\[
\mathcal I_{A:B}(\rho)=\frac{\operatorname{Tr}(A\rho)}{\operatorname{Tr}(B\rho)}.
                                                               \tag{1}
\]

If numerator and denominator share the same units, (1) is dimensionless. Otherwise use independently specified nonzero calibration readings and form `(n/n_ref)/(d/d_ref)`. Being dimensionless alone does not prove conservation. The definition of an invariant under transport `Phi_s` is

\[
\mathcal I_{A:B}(\Phi_s\rho)=\mathcal I_{A:B}(\rho)          \tag{2}
\]

on its declared domain. No ratio is evaluated at a zero denominator. Frame covariance is a different statement: conjugating `rho,A,B` by the same unitary leaves (1) unchanged by trace cyclicity. Holding a physical observation fixed while transporting the state requires (2), not just that change-of-frame identity.

## 2. Equal relative rates and the precise diagonal condition

For positive differentiable readings `n(s),d(s)`, direct differentiation gives

\[
\boxed{\frac{d}{ds}\log\frac{n}{d}
=\frac{n'}n-\frac{d'}d.}                                  \tag{3}
\]

Thus the ratio is invariant on an interval exactly when their relative rates coincide. In that case there is a common positive multiplier `q(s)` with
`n(s)=q(s)n(0)`, `d(s)=q(s)d(0)`.

In the logarithmic coordinates

\[
X=\log[n(s)/n(0)],\qquad Y=\log[d(s)/d(0)],
\]

equation (3) gives `X=Y`. With equal axis scales and the Euclidean chart metric `dX²+dY²`, this is the diagonal line at 45 degrees to the axes. Its nonzero tangent can point forward or backward; if `q'=0` there is no tangent direction to assign. Under independent rescalings `X_new=alpha X`, `Y_new=beta Y`, its plotted slope is `beta/alpha`. If the metric is transformed too, intrinsic angles are unchanged; redrawing both axes with a new Euclidean metric is a different calibration. Dimensionless coordinates do not prohibit such rescaling.

On the raw reading plane, `n/d=r` is instead the ray `n=r d`, whose slope depends on `r`. The invariant statement is (3), not an uncalibrated universal geometric angle.

**IRT-T1 (when the common multiplier is exponential).** Suppose in addition `q(0)=1`, `q(s+t)=q(s)q(t)` for `s,t>=0`, and `q` is positive and continuous. Then `q(s)=exp(kappa s)` for a real constant `kappa`. Conversely this exponential has all these properties.

**Proof.** `h=log q` is continuous and additive. Integer multiplication and subdivision give `h(r)=r h(1)` for nonnegative rational `r`; continuity extends this to all nonnegative real arguments. Exponentiation proves the assertion. QED.

Common variable-rate multipliers also preserve the ratio and the log diagonal. For example `q(s)=1/(1+s)` has relative rate `-1/(1+s)` and does not obey the semigroup law: `q(2)=1/3` while `q(1)²=1/4`. Thus a diagonal alone does not prove constant-rate exponential evolution. Conversely conserved readings under unitary transport may be constant, with the log plot fixed at the origin; their invariance requires no nonzero 45-degree trajectory.

## 3. An all-state ratio criterion for a linear quantum operation

Let `L` be a real-linear Hermiticity-preserving generator on a finite matrix carrier, with dual `L*` in the trace pairing. Assume `A,B` are linearly independent Hermitian observables, `B>=0`, `B!=0`.

**IRT-T2 (common dual rate).** The infinitesimal ratio derivative vanishes for every positive definite initial source if and only if a single real scalar `kappa` satisfies

\[
\boxed{L^*A=\kappa A,\qquad L^*B=\kappa B.}              \tag{4}
\]

Then their readings acquire the same `exp(kappa s)` factor under `exp(sL)`, so their ratio is invariant wherever it is defined. If `A=cB`, the ratio is identically `c` and no such condition is needed; this degenerate case is excluded from the theorem.

**Proof.** Write real linear forms on the full Hermitian vector space as
`a(rho)=Tr(A rho)`, `b(rho)=Tr(B rho)`, `u(rho)=Tr((L*A)rho)` and `v(rho)=Tr((L*B)rho)`. The derivative numerator is the homogeneous polynomial `u b-a v`. Vanishing on the open positive cone forces it to vanish everywhere. On `ker(a)`, points where `b!=0` are dense because `a,b` are independent. Hence `u` vanishes on all of `ker(a)` and must equal `kappa a`. Substitution gives `a(kappa b-v)=0` identically, so `v=kappa b`. Nondegeneracy of the trace pairing yields (4). The reverse direction is substitution; exponentiation of the two dual eigenrelations gives the finite transport statement. QED.

This theorem concerns all initial sources for one specified generator. A ratio that happens to be conserved on one trajectory need not satisfy (4) on the whole carrier.

For finite unitary transport, `L*A=(iota/eta)[H,A]`. Equation (4) forces `kappa=0`: multiply its first relation by `A`, take the trace, and use cyclicity to obtain `kappa Tr(A²)=0`. Independence ensures `A!=0`. Thus, in this finite all-state setting,

\[
\boxed{\mathcal I_{A:B}\text{ is invariant}
\quad\Longleftrightarrow\quad [H,A]=[H,B]=0.}             \tag{5}
\]

This does not say that invariant ratios distinguish all quantum states. It identifies the observables needed for this specified universal conservation statement.

## 4. A surviving invariant in the coupled AGD model

On the periodic geometry carrier let

\[
(S\phi)(\theta)=\phi(-\theta),\qquad S|n\rangle=|-n\rangle,
\quad \Pi_\pm=(I\pm S)/2.
\]

This is **reflection parity**, not AGD's half-turn translation `T_pi|n>=(-1)^n|n>`. The latter was correctly shown to fail under the lapse coupling. Reflection obeys

\[
S^\dagger=S,\quad S^2=I,\quad S U S=U^\dagger,
\quad S L S=-L.
\]

Consequently, for the inherited operators,

\[
H_g=\eta^2 L^2+8-(U^2+U^{\dagger2}),\qquad
V(x)=\left(1+\frac{af(x)}2\right)I-\frac{af(x)}4(U+U^\dagger),
\]

one has `[S,H_g]=[S,V(x)]=0`. The spatial derivative, direction matrices, lapse derivative and mass term also commute with `S`.

**IRT-T3 (exact transported geometry ratio).** In the full variable-profile AGD evolution, for every joint density operator,

\[
w_\pm(s)=\operatorname{Tr}[(I_{\rm matter}\otimes\Pi_\pm)\rho(s)]
=w_\pm(0).
\]

Therefore

\[
\boxed{\mathcal I_{+:-}=w_+/w_-}                          \tag{6}
\]

is an exact dimensionless invariant whenever `w_->0`. The normalized alternative `w_+/(w_++w_-)` is defined for every nonzero source.

**Proof.** Reflection preserves AGD's smooth core and commutes with its generator there. Unitary conjugation by `S` therefore maps the essentially self-adjoint operator to itself, and its unique closure and spectral exponential commute with `S`. Both bounded projections reduce the group. Trace cyclicity establishes the claim for trace-class densities without imposing finite geometry support. QED.

For the normalized initial geometry vector
`(|0>+|1>-|-1>)/sqrt(3)`, the weights are `1/3,2/3`, and the invariant is exactly `1/2`. This holds for any initial matter factor and all admitted AGD operation parameters. It does not require the geometry vector itself to be stationary.

These are conserved sector weights, not local gravitational gauge constraints. Both sectors still contain nonconstant lapse and quantum geometry dispersion.

## 5. Complete geometry record for the specified common commutant

**IRT-T4 (classification).** For `eta>0` and a lapse profile nonconstant in geometry at at least one spatial point, the bounded geometry operators commuting with every `V(x)` and strongly commuting with `H_g` are exactly

\[
\boxed{\mathcal A_{\rm inv}=\{\alpha I+\beta S:\alpha,\beta\in\mathbb C\}.} \tag{7}
\]

Hermitian readings have real coefficients. Their expectation values are determined by the two weights `w_+,w_-`; for normalized sources only one independent weight is required.

**Proof.** Commutation with a nonconstant `V(x)` is equivalent to commutation with the bounded multiplication operator `cos(theta)`. It then implies commutation with `cos(2theta)` as well. Since `H_g` has the same domain as `L²` and differs from `eta²L²` by a bounded polynomial in `cos(theta)`, the hypotheses imply that the bounded operator and its adjoint preserve that domain and commute with `L²`.

Use the orthonormal basis `c_0=|0>`, `c_n=(|n>+|-n>)/sqrt(2)`, `s_n=(|n>-|-n>)/sqrt(2)` for `n>=1`. The `L²` eigenvalue zero is simple, so `A c_0=alpha c_0`. Multiplication by cosine sends `c_0` to `c_1/sqrt(2)`. For `n>=2` it sends `c_n` to `(c_(n-1)+c_(n+1))/2`, and sends `c_1` to `c_0/sqrt(2)+c_2/2`. Commutation therefore propagates `A c_n=alpha c_n` to every `n`. Applying the same argument to `A†` shows that the complementary sine subspace is invariant too.

Within that subspace each `L²` eigenspace is one-dimensional, hence `A s_n=beta_n s_n`. The cosine links `s_1` to `s_2/2` and `s_n` to `(s_(n-1)+s_(n+1))/2`, forcing all `beta_n` equal. Thus `A=alpha Pi_++beta_1 Pi_-`, equivalent to (7). Conversely both identity and reflection commute with the two generators. QED.

Equation (7) classifies the **common geometry-only commutant** of the two specified generator ingredients. It is not a classification of every invariant of the full interacting theory: functions of the full joint Hamiltonian, for example, have a different scope. Nor does it classify observables conserved only at a single specially tuned matter mode.

The finite rank fixtures classify the commutant of compressed `L²` and `cos(theta)`. Truncating `H_g` itself can introduce accidental symmetries: at `eta=1` in the three-mode section `n=-1,0,1`, the compressed `H_g` commutes with compressed cosine, although the infinite operators do not. A five-mode section already removes that accidental commutation. The infinite classification is proved by the basis recurrence above, not inferred from the small matrix ranks.

The invariant record is not the complete state. The states `|0>` and `(|1>+|-1>)/sqrt(2)` both lie wholly in the reflection-even sector. In the homogeneous lapse control `V=13I/10-3(U+U†)/20`, both have mean `13/10`, but their lapse variances are respectively

\[
\boxed{9/200\quad\hbox{and}\quad27/400.}                  \tag{8}
\]

The common invariant algebra cannot distinguish these states; a second-moment geometry reading can. This is an application of CRC's record-completeness rule, not a contradiction in invariance.

## 6. Common attenuation: a quantum example of the log diagonal

To distinguish conservation of a ratio from conservation of its amount, introduce a separate explicit finite adapter with logical cut `P_L=|0><0|+|1><1|` and erasure flag `|e>`. For `0<q<=1` take

\[
K_0=\sqrt q\,P_L+|e\rangle\langle e|,
\quad K_j=\sqrt{1-q}|e\rangle\langle j|\quad(j=0,1).
\]

Their dagger squares sum to identity; the resulting channel `Phi_q(rho)=sum K_j rho K_j†` is positive under tensor extension and trace preserving. Its block action is

\[
\Phi_q\begin{pmatrix}\rho_L&b\\b^\dagger&r_e\end{pmatrix}
=\begin{pmatrix}q\rho_L&\sqrt q\,b\\\sqrt q\,b^\dagger&r_e+(1-q)\operatorname{Tr}\rho_L\end{pmatrix}.
                                                               \tag{9}
\]

Composition gives `Phi_q Phi_r=Phi_(qr)`. Every logical observable `A=P_L A P_L` has `Phi_q* A=qA`, including the denominator `P_L`. Thus the entire conditional logical state is preserved wherever its survival weight is nonzero, while the erasure record changes. This common-attenuation model is specified here, not inferred as AGD's closed gravitational dynamics.

With `q(s)=exp(-gamma s)`, all positive logical readings and their survival denominator lie on the equal-scaled log diagonal. Reparameterizing the same family by `q(s)=1/(1+s)` preserves their ratios but removes the constant-rate semigroup in `s`.

An exact native fixture uses logical populations `4/5,1/5` with retained coherence and `q=9/25`. Then the first population is `36/125`, logical survival is `9/25`, their ratio stays `4/5`, and the erasure weight is `16/25`. A ratio alone does not recover that survival weight. Common attenuation transports the retained conditional information while storing a separate outcome record.

## 7. Scope, certificates and next obligation

The packet supplies an explicit invariant in the actual coupled AGD model and classifies its stated common geometry observable algebra. It also specifies exactly when the 45-degree representation and a constant exponential rate apply. No invariant is relabelled unstable; a proposed quantity either satisfies its specified conservation law or does not.

The native certificate verifies exact finite commutants, reflection and broken-half-turn controls, rational unitary transport, conserved weights, distinct state variances, ratio derivatives, quantum-channel identities and the chart/semigroup counterexamples. The infinite common-commutant classification and the general linear-form argument are written proofs, not consequences of checking finitely many examples. The usual commutant viewpoint is established mathematical infrastructure; see, for context, [Moudgalya and Motrunich, arXiv:2209.03370](https://arxiv.org/abs/2209.03370). No general novelty claim for symmetry conservation or the exponential functional equation is made.

The next gravity obligation remains to derive local coframe/source constraints and determine which of these invariant records survive that extension. This packet does not identify a physical gravitational rate, impose a universal physical 45-degree angle, or complete quantum gravity.
