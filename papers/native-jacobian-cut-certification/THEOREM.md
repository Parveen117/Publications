# NJC-1 — Certified core of the Jacobian and exponential manuscripts

Development date: 14 October 2026. This packet audits eight supplied files,
including two late-added exponential manuscripts. It certifies the statements
below under their stated contracts, rather than certifying every assertion in
the source documents. See [CLAIM_AUDIT.md](CLAIM_AUDIT.md) for rejected and
unestablished assertions.

## Contract and lineage

Arithmetic and algebraic proof witnesses use the unchanged native cut field
and Pāṇinian compiler. The executable coefficients are exact rationals and
their native quarter-turn extension. No SymPy, NumPy, floating-point physics,
or replacement arithmetic engine is used by the new certificate.

Three meanings must stay distinct: the real-linear reflection `J`, the
differential Jacobian `D Lambda`, and a selected evolution generator `G`.
Giving them similar names does not identify their types or dynamics.

The real two-component realization of `J(z)=iota conjugate(z)` is

\[
J=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad RJ=-JR.
\]

`R` is multiplication by the original plane's quarter-turn. A central cut
scalar acting on a subsequently complexified two-component module is a
different operation. No primitive Hilbert-space postulate is needed for these
finite modules, but scalar extension and a positive pairing are declared
structures, not consequences of the word `J`.

Analytic limits below use the existing completed scalar interpretation and,
where applicable, a bounded finite-dimensional generator. Word certificates
prove identities in declared associative quotients; the written analytic
proofs are not a Lean/Coq formalization of completion, geometry or physics.

The connection construction is **reused**, not a new discovery: CID-D1 D2e in
[CUT_INFORMATION_DYNAMICS.md](../native-loop-concentration/CUT_INFORMATION_DYNAMICS.md)
already proves the minimal cut-compatible correction and its curvature.
[EMK-C1](../emk-ugd-algebra/CONNECTION_CALCULUS.md) supplies connection covariance
and the frame-corrected curvature convention. The recovered RKF T42 supplies
the response tower; the recovered MP bridge supplies the quadratic-generator
trichotomy. This packet corrects the uploaded formulations, adds native
witnesses, and records their precise use in the current normal-plane problem.

## T1. Reflection, seam and sine residual

For a self-adjoint involution `J` over characteristic zero, put

\[
P_\pm=\tfrac12(I\pm J).
\]

Then these are complementary orthogonal projectors, and

\[
Jz-z=-2P_-z,\qquad
\|Jz-z\|^2=4\|P_-z\|^2.
\]

In the two-component realization the fixed set is the **whole real line**
`{a(1+iota): a real}`, including zero. For
`z=r exp(iota(pi/4+delta))`, with `r>=0`, the completed angular representation
gives exactly

\[
Jz-z=-2\iota r e^{\iota\pi/4}\sin\delta.
\]

**Proof.** Expand the projectors using `J^2=I`; self-adjointness gives
orthogonality. The residual identity follows by subtraction. Conjugate the
angular expression and subtract; the two exponentials differ by
`-2 iota sin(delta)`. This identity is exact at all angles, not just near zero.
The normalized defect is `sin^2(delta)` when `r>0`.

The numerical angle `pi/4` describes the chosen coordinate presentation. The
intrinsic statement is the fixed projector and its complementary residual.
It neither selects a unit radius nor establishes a measurement-lag law.

## T2. Correct tensor grading and comparison seam

For a module `V=V_+ direct-sum V_-`, opposite eigenspaces intersect only at
zero: `v=Jv=-Jv` implies `2v=0`. Consequently equality of an element of `V_+`
with an element of `V_-` cannot define a nontrivial diagonal seam.

There are two consistent constructions:

* A single state decomposes as `z=P_+z+P_-z`, with seam `P_-z=0`.
* Observer and measured are two copies of a common fibre with a specified
  comparison isomorphism `T`. Their comparison seam is `m=T e`; a distance
  such as `||m-T e||` additionally requires a compatible pairing.

For an `n`-slot tensor use the projectors on **each** slot, and the dual
projectors on covector slots. The `2^n` products of slot projectors are
pairwise disjoint and sum to identity. This gives the unique full parity
decomposition. A three-way grouping into all-plus, all-minus and mixed terms
is then defined by grouping these projectors.

**Proof.** Distribute `(P_++P_-)` over the slots. Idempotence and disjointness
follow slot by slot. This also gives an induction in `n`. The native
certificate explicitly replays the two-slot algebra. Ordinary vector and
covector frame laws follow from preserving evaluation; a phase factor is
admitted only after specifying its actual representation.

For real matrices, `C(A)=JAJ` and `T(A)=A^T` commute when `J^T=J` and `J^2=I`.
Their joint sectors are

\[
A^{\sigma,\tau}=\tfrac14
(A+\sigma JAJ+\tau A^T+\sigma\tau JA^TJ).
\]

They reconstruct `A` and have the stated cut and transpose signs. The
orthogonality condition is material: `J=[[1,1],[0,-1]]` is an involution for
which these two operations need not commute. The recovered T42 finite
fixture uses a diagonal orthogonal cut and passes; its unrestricted wording
needs this hypothesis.

## T3. Cut-compatible connection and the freedom it leaves

Work in a declared flat reference trivialization with commuting derivatives.
Let `J(x)^2=I`, `J(x)^dagger=J(x)` be a differentiable field, and let

\[
D_i=\partial_i+A_i,\qquad D_iJ=\partial_iJ+[A_i,J].
\]

All cut-compatible coefficients are exactly

\[
\boxed{A_i=\tfrac12J\partial_iJ+B_i,\qquad [B_i,J]=0.}
\]

For a compatible positive pairing require `B_i^dagger=-B_i`. The first term
is already anti-dagger. Setting `B=0` is the declared minimal correction rule.
For that choice,

\[
\boxed{F_{ij}=\tfrac14[\partial_iJ,\partial_jJ]
       =[\partial_iP_+,\partial_jP_+],\qquad [F_{ij},J]=0.}
\]

**Proof.** Differentiate `J^2=I` to obtain `J J_i+J_i J=0`. Thus
`[J J_i/2,J]=-J_i`. Conversely the odd component of any compatible `A_i` is
`(A_i-JA_iJ)/2=J J_i/2`, leaving an arbitrary commuting component. The
derivative part of the curvature is `[J_i,J_j]/2`; the coefficient
commutator is `-[J_i,J_j]/4`. Their sum is the displayed result. A product of
two cut-odd derivatives is cut-even. These identities have universal native
rewrite witnesses.

In contrast, the scalar form in the uploaded tensor draft vanishes:

\[
\operatorname{tr}(J\partial_iJ)=0.
\]

Indeed cyclicity equates `tr(J J_i)` and `tr(J_i J)`, while the differentiated
involution makes their sum zero. Taking the trace has discarded the
off-diagonal transport. It cannot supply the draft's claimed nonzero
connection term.

For a non-flat reference connection use CID-D1 directly. With
`X_i=D_i^0 P`, the corrected connection and curvature are

\[
A_i^\Sigma=A_i^0+[P,X_i],\qquad
F^\Sigma_{ij}=PF^0_{ij}P+(I-P)F^0_{ij}(I-P)+[X_i,X_j].
\]

A frame change must transform the reference connection too. One may not
conjugate `J` by a position-dependent frame and retain the same bare `d`
without its connection term. None of these formulas selects a torsion-free
tangent connection or a gravitational action by itself.

## T4. Quadratic-generator and semigroup bridge

For a bounded generator `G` with `G^2=Delta I` and real central `Delta`,

\[
e^{tG}=\begin{cases}
\cos(\omega t)I+\sin(\omega t)G/\omega,&\Delta=-\omega^2<0,\\
I+tG,&\Delta=0,\\
\cosh(\alpha t)I+\sinh(\alpha t)G/\alpha,&\Delta=\alpha^2>0.
\end{cases}
\]

**Proof.** Induction gives `G^(2n)=Delta^n I` and
`G^(2n+1)=Delta^n G`; split the norm-convergent exponential series. The native
certificate supplies the quadratic rewrite and exact coefficient replays
through order ten. The induction and bounded-series argument establish the
general formula; ten coefficients alone would not do so.

For `Delta>0`, `P_±=(I±G/alpha)/2` give

\[
e^{tG}=e^{\alpha t}P_++e^{-\alpha t}P_-.
\]

These are algebraic projectors; orthogonality requires additional adjoint
conditions. Stable-sector observation loses `P_+V`. A linear target `Pi`
factors through `P_-` exactly when `Pi P_+=0`, as follows by decomposing every
state into its two sectors. Both sectors reconstruct the full state.
The first-order stable rate is `sqrt(Delta)`, whereas the generator
`G^2=Delta I` has rate `Delta`. A threshold written `Delta=rho+1/4` is a
coordinate calibration unless a further normalization fixes it.

For exact finite transport, when `1-t^2 Delta` is invertible,

\[
(I-tG)^{-1}(I+tG)
=\frac{(1+t^2\Delta)I+2tG}{1-t^2\Delta}.
\]

Multiplication verifies the inverse and the formula. Circular spectral
behaviour alone does not imply positive-norm unitarity; a nonorthogonal
similarity of a rotation provides a counterexample.

## T5. Relaxation, real-linear exponentials and forced series

For the original real-linear reflection,

\[
e^{tJ}=e^tP_++e^{-t}P_-.
\]

It has inverse `e^(-tJ)` at every finite parameter value. For nonzero retained
component its norm grows; distance to the fixed **line** does decay as
`e^(-t)||P_-z||`. Those two statements are compatible. It does not converge
to a finite fixed point for generic initial data.

A bounded relaxation preserving the retained component is

\[
\boxed{T(t)=e^{-\lambda tP_-}=P_++e^{-\lambda t}P_-,\quad\lambda>0.}
\]

Its generator, norm contraction, composition law and limit `T(t)->P_+`
follow immediately from the two projectors. It remains invertible for finite
`t` on the finite carrier; projection and information loss occur in the
limit. Identifying the parameter with physical time or the discarded norm
with entropy production requires a physical adapter.

The late exponential draft also moves the original imaginary unit through
an anti-linear `J`. In the real two-component realization,

\[
RJ=-JR,\qquad (RJ)^2=+I.
\]

Therefore exponentiating **left multiplication by the plane's quarter-turn
followed by J** is hyperbolic, not the stated circular formula. If instead
the real module is complexified and `J` is extended complex-linearly, a
central scalar `iota` commutes with `J`, and

\[
(\iota J)^2=-I,\qquad e^{\tau\iota J}=\cos\tau I+\iota\sin\tau J.
\]

This latter transport is unitary for a self-adjoint `J` and the declared
positive pairing. It is a different, explicitly typed representation.

For `T_a=P_++aP_-`, the forced recurrence `z_(n+1)=T_a z_n+b_(n+1)` has

\[
P_+z_n=P_+z_0+\sum_{k=1}^nP_+b_k,\qquad
P_-z_n=a^nP_-z_0+\sum_{k=1}^na^{n-k}P_-b_k.
\]

**Proof.** Project the recurrence and induct. Both components must satisfy
the desired convergence condition. A vanishing transverse component alone
does not imply convergence: `b_k=1+iota`, `z_0=0` gives
`z_n=n(1+iota)` entirely on the seam.

## T6. Correct time-dependent memory propagator

Assume the draft's constitutive law `M(t)=R(1-exp(-lambda*t))/lambda`, with
`R,lambda>0`, and interpret scalar subtraction as a multiple of identity.
Then `L(t)=J-(lambda/R)M(t)I=J-m(t)I`, `m(t)=1-exp(-lambda*t)`. Coefficients
commute at distinct parameters, so the initial-value propagator is

\[
\boxed{U(t)=\exp\left(tJ-
\left[t-\frac{1-e^{-\lambda t}}{\lambda}\right]I\right).}
\]

**Proof.** Differentiate the displayed expression; its exponent derivative
is `L(t)` and commutes with it. The initial value is identity. This yields
`U(t)->exp(1/lambda)P_+`. The limit depends on the initial retained component,
so it is not a uniquely selected maximum-entropy state.

The pointwise expression `V(t)=exp(t L(t))` instead satisfies

\[
V'(t)=[J-m(t)I-tm'(t)I]V(t),
\]

not the source's `[J-2m(t)I]V(t)`. Their coefficient difference is
`m-tm'=1-(1+x)exp(-x)=x^2/2+O(x^3)`, `x=lambda*t`, which is nonzero.
The integral of the assumed positive cost is valid, but calling it entropy
or memory is a constitutive choice requiring calibration.

## T7. Response tower with explicit hypotheses

Specify connections on both the response bundle and tangent bundle, together
with compatible cuts. Then `L_1=Lambda`, `L_(n+1)=nabla L_n` has type
`E tensor (T*M)^(n-1)`. Each induced involution squares to identity, and its
two parity projectors give the cut grading. In a flat affine chart with
constant linear cut, differentiating `Lambda(J_x x)=J_E Lambda(x)` gives

\[
D^k\Lambda(J_xx)[J_xh_1,\ldots,J_xh_k]
=J_E D^k\Lambda(x)[h_1,\ldots,h_k].
\]

For a nonlinear base involution use the covariant, cut-compatible statement;
ordinary higher derivatives acquire higher derivatives of the base map.
Polynomial termination requires the flat affine derivative adapter. The
tower exponential is a formal graded generating series unless a completion
and convergence topology have been specified.

For linearized observation `A_N` and linear target `Pi`, the actual blind
target is `B_N=Pi(ker A_N)`. Stacking measurements shrinks these images, and

\[
\dim B_N=\operatorname{rank}\binom{A_N}{\Pi}
          -\operatorname{rank}A_N.
\]

**Proof.** Apply rank-nullity to `Pi` restricted to `ker A_N`. Full recognition
means `ker A_N subset ker Pi`; strict target repair means strict decrease of
the blind image. Intersecting the observation kernel with an unspecified
`support(Pi)` is not an equivalent criterion. For `A=(1,-1)` and `Pi=(1,0)`,
the blind vector `(1,1)` has a nonzero target even though the kernel misses
the orthogonal complement of `ker Pi`.

The recovered T42's ten legacy tests pass and its old payload hash matches.
They are historical Fraction-based fixtures. The new native polynomial
certificate rederives its value, Jacobian and Hessian covariance and exact
quadratic termination. Neither test suite certifies a physical TVSP
constitutive map that has not been supplied.

## T8. Explicit use in the four-dimensional normal-plane problem

Take `J0=diag(1,1,-1,-1)`, so `P_+` retains a two-dimensional face and `P_-`
is its normal plane. Let `X=E31-E13` and `Y=E41-E14` (indices 1 through 4),
and realize a smooth field by conjugating `J0` with the corresponding local
orthogonal transports. At the origin,

\[
J_x=[X,J_0],\qquad J_y=[Y,J_0],\qquad
F_{xy}=\tfrac14[J_x,J_y]=E_{34}-E_{43}.
\]

Its exact squared native matrix norm is `2`; it acts in the normal plane.
The certificate checks the coefficient derivatives independently, including
the mixed second derivative of `J`, rather than testing only a frozen
connection commutator. This is a finite geometric example, not a measured
curvature or a prediction of Newton's constant.

There is also a precise information limitation. Keep `J=J0` constant. Both
`B=0` and `B_x=0, B_y=xC`, with `C=E43-E34`, preserve the identical cut and
its entire derivative tower. The first connection is flat and the second
has `F_xy=C`. Thus **even the full cut-Jacobian tower cannot determine an
independently allowed normal-plane connection**. A transport/holonomy
measurement, action principle or additional condition fixing `B` is needed.
This specializes the existing EMK-C1 separation of response and connection
data to the current gravity construction.

The MP trichotomy gives the normal rotation and boost branches once the
normal generator and metric signature are declared. The positive normal
rotation is unitary; a Lorentz boost preserves its indefinite metric and
must not be relabeled a positive-norm quantum unitary.

## What the certificate does not assert

The primitive reflection alone does not select a gravitational action,
stress tensor, Einstein constraint algebra, physical clock, universal
coupling, thermodynamic entropy, integer topological charge or a solution of
RH/Yang–Mills/Navier–Stokes. Rewriting a proposed identification as a
definition would not prove it. The useful results here are the corrected
cut, transport, generator and observation statements above.
