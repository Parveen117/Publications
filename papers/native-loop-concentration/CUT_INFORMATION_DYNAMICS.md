# CID-D1 — A covariant cut–information–curvature equation

14 October 2026. Companion to QG-R2, derived on the finite native carrier.

The existing [CID-1 / QTH-1](../curvature-information-duality/README.md)
already connect information and curvature. CID-1 identifies one obstruction
across information metrics in its response algebra. QTH-1 exhibits the
symmetric information and antisymmetric curvature parts of one tensor in
declared rational state families. Here an evolving native cut produces a
tensor, its induced connection curvature, and a single covariant transport
equation for both parts. No universal physical field law is presumed.

## D1. The evolving cut and its generator

Work on a finite native positive coefficient carrier. Let
`P^2=P=P^dagger` be a differentiable observation cut and `K=2P-1` its
involution. Let an anti-dagger reference connection `A_a` define

\[
D_a Z=\partial_a Z+[A_a,Z],\qquad
F_{ab}=\partial_aA_b-\partial_bA_a+[A_a,A_b].
\]

The control derivatives commute in this presentation; they are not assumed
spacetime coordinates. For a noncommuting direction frame, subtract the
EMK-C1 frame term in curvature and include its corresponding term in
commuted derivatives. The displayed equations use the commuting case.
Put `X_a=D_aP`. Differentiating `P^2=P` gives

\[
PX_a+X_aP=X_a,\qquad PX_aP=0.
\tag{D1a}
\]

Thus cut change occupies the off-diagonal cut blocks. Along an admitted
control parameter `lambda`, a relative anti-dagger generator `V_lambda`
gives the evolution

\[
\boxed{D_\lambda P=[V_\lambda,P].} \tag{D1b}
\]

For any given cut derivative `X_lambda`, its unique off-diagonal generator is

\[
V_\lambda^{\rm min}=[X_\lambda,P]
=\tfrac14[D_\lambda K,K],\qquad
[[X_\lambda,P],P]=X_\lambda.
\]

Any additional generator commuting with `P` is invisible to that cut
derivative. In a reference gauge `A_lambda=0`, a **constant** generator
has the native exponential solution

\[
\boxed{P(\lambda)=\operatorname{Exp}_\Sigma(\lambda V)
P(0)\operatorname{Exp}_\Sigma(-\lambda V).} \tag{D1c}
\]

For a varying generator use ordered transport solving
`partial_lambda U=(V_lambda-A_lambda)U`, then `P=UP_0U^dagger`.
The finite-core native exponential and its completion are inherited from
GE1; the ordered-transport requirement is already EMK-T2. An exponential
of the cut itself commutes with that cut and leaves it unchanged under
conjugation. The exponent driving cut motion is its transport generator.
Neither an arbitrary generator nor a unique physical lambda-clock is
selected by this identity.

## D2. One cut tensor, information and curvature

Define an operator-valued tensor on the retained cut and its scalar reading:

\[
\mathsf T_{ab}=PX_aX_bP,\qquad
\boxed{\mathcal Q_{ab}=\operatorname{Tr}(PX_aX_bP)
=g^{\rm cut}_{ab}+\iota\,\Omega^{\rm cut}_{ab}.} \tag{D2a}
\]

Here

\[
g^{\rm cut}_{ab}=\tfrac12\operatorname{Tr}\big(P\{X_a,X_b\}P\big),
\quad
\iota\Omega^{\rm cut}_{ab}
=\tfrac12\operatorname{Tr}\big(P[X_a,X_b]P\big).
\]

The first is real symmetric and positive semidefinite; the second is real
antisymmetric after the displayed `iota` factor. Indeed, for native
coefficients `z_a`, the quadratic form of `Q` is the native squared norm
of `sum z_a X_aP`. Its null directions are precisely the directions that
do not change the retained cut in this norm. This is a derived cut
sensitivity metric. A general identification with a statistical Fisher
metric requires a state/measurement representation.

Under a native unit frame change, transform both
`P'=G^dagger P G` and `A'_a=G^dagger A_aG+G^dagger partial_aG`.
EMK-C1 then gives `X'_a=G^dagger X_aG` and
`T'_ab=G^dagger T_abG`; cyclic trace makes `Q_ab` invariant, including
for variable frames. Changing the cut while holding its reference fixed
is different from merely changing the frame of the same experiment.

There is a precise link to QTH-1 in the trace-one rank-one cut sector:
use `rho=P` and the valid score choice `L_a=2X_a`. D1a gives
`D_a rho=(rho L_a+L_a rho)/2`, and therefore

\[
\operatorname{Tr}(\rho L_aL_b)=4\mathcal Q_{ab}. \tag{D2b}
\]

Thus its information and antisymmetric halves are the corresponding
QTH-type tensor with an explicit factor of four. This is a pure-cut
representation result, not an identification for every mixed state or a
claim that all previous information definitions are identical.

The induced connection on retained sections is `nabla^P=P nabla`.
Its curvature is

\[
\boxed{F^P_{ab}=PF_{ab}P+P[X_a,X_b]P.} \tag{D2c}
\]

Consequently the information-curvature equation is

\[
\boxed{\mathcal Q_{ab}
=g^{\rm cut}_{ab}
+\tfrac12\operatorname{Tr}(F^P_{ab}-PF_{ab}P).} \tag{D2d}
\]

The final term is anti-dagger/pure-turn, equal to `iota Omega_ab`.
This identifies the **cut-induced curvature increment**, with the
ambient curvature retained explicitly. It does not identify an arbitrary
metric curvature with arbitrary recognition curvature; EMK-G1's separation
continues to hold.

**Proof.** Apply `P nabla_a P nabla_b` to a retained section and
antisymmetrize. The connection curvature contributes `PF_abP`; the two
derivatives of `P` contribute `P[X_a,X_b]P`, using D1a. Symmetric and
antisymmetric decomposition of `PX_aX_bP` proves D2d. All operations
are native composition, dagger, cut and coefficient trace.

An equivalent full-carrier connection that parallel-transports the cut is

\[
A^\Sigma_a=A_a+[P,X_a],\qquad D^\Sigma_aP=0,
\]
\[
\boxed{F^\Sigma_{ab}
=PF_{ab}P+(1-P)F_{ab}(1-P)+[X_a,X_b].} \tag{D2e}
\]

It is the unique correction with no cut-diagonal component that makes the
cut parallel. This is a stated minimal correction rule; extra cut-diagonal
transport remains possible. To prove the curvature formula, differentiate
`[P,X]`: its exterior covariant derivative is
`2 X wedge X+[P,[F,P]]`, while `[P,X] wedge [P,X]=-X wedge X`.
Finally `F+[P,[F,P]]=PFP+(1-P)F(1-P)`. Restricting D2e to `P` recovers D2c.

## D3. The covariant information–curvature evolution equation

Let `X_lambda=D_lambda P=[V_lambda,P]` and

\[
Y_{\lambda a}=D_aX_\lambda+[F_{\lambda a},P].
\]

EMK-C1's commuted covariant derivatives give
`D_lambda X_a=Y_(lambda a)`. Define the derivative of a retained operator
by `mathscr D_lambda T=P(D_lambda T)P`. Then the requested joint
evolution law is

\[
\boxed{\begin{aligned}
\mathscr D_\lambda\mathsf T_{ab}
 &=P\big(Y_{\lambda a}X_b+X_aY_{\lambda b}\big)P,\\
\partial_\lambda\big(g^{\rm cut}_{ab}+\iota\Omega^{\rm cut}_{ab}\big)
 &=\operatorname{Tr}\!\left[
P\big((D_a[V_\lambda,P]+[F_{\lambda a},P])D_bP
+D_aP(D_b[V_\lambda,P]+[F_{\lambda b},P])\big)P\right].
\end{aligned}} \tag{D3}
\]

Its symmetric part evolves the information metric; its antisymmetric
part evolves the cut-induced curvature. Their source is the same cut
transport `V_lambda`, with the reference curvature correction explicit.
The scalar equation uses the ordinary derivative because the trace of
the connection commutator vanishes.

**Proof.** Differentiate `PX_aX_bP` by Leibniz, obtaining four terms.
After sandwiching by `P`, the two terms differentiating the outer cuts
vanish: `PX_lambdaP=0`, and `X_aX_b` commutes with `P` by D1a.
The remaining two terms are the first line. Substitute the commuted
derivative identity and D1b. The trace of the discarded off-diagonal
terms is also zero, giving the scalar line. The native compiler certifies
the corner cancellation and the free noncommuting derivative identity.
No claim that D3 determines `V_lambda` from the initial cut is made:
it is a structural unification and transport equation, not a selected
closed physical source law.

## D4. Exact two-mode calibration and a changing-curvature example

At the reference cut `P_0=diag(1,0)`, any Hermitian cut tangent has form
`X_a=[[0,z_a^dagger],[z_a,0]]`. Thus `Q_ab=z_a^dagger z_b` and, for two
control directions,

\[
\boxed{\det(g^{\rm cut})=(\Omega^{\rm cut}_{12})^2.} \tag{D4a}
\]

Writing `z_1=x+iota y`, `z_2=u+iota v` proves this as
`(x^2+y^2)(u^2+v^2)-(xu+yv)^2=(xv-yu)^2`.
Native unit frame covariance transports it to every rank-one two-mode
cut. Its equality is special to this carrier; in larger carriers the
corresponding two-direction bound follows from positivity, without
asserting equality.

With flat ambient reference, the minimal cut connection has
`A^Sigma=[P,dP]`. In this two-mode sector its relative-connection Gram
matrix from QG-R1 equals `g^cut`, and
`F^Sigma_12=2 iota Omega_12(2P-1)`. Hence wherever the Gram determinant
is positive,

\[
\boxed{\frac{\langle F^\Sigma_{12},F^\Sigma_{12}\rangle}
{4\det g^{\rm cut}}=1.} \tag{D4b}
\]

This is QG-R1's regular small-cell return-density calibration for the
selected minimal cut connection. It is not a universal value for arbitrary
connections or an assertion that bounded density proves regularity.
For the direct proof at `P_0`, `[P_0,X_a]` has off-diagonal entries
`z_a^dagger,-z_a`. Its native anti-dagger pairing with `[P_0,X_b]` is
`Re(z_a^dagger z_b)`, while `[X_1,X_2]=2 iota Omega_12(2P_0-1)`.
Covariance transports both identities to an admitted unit frame.

An explicit family shows information and curvature changing together.
Put `a=1+lambda`, and take
`U=Exp_Sigma(a x R) Exp_Sigma(y J)`, `P=UP_0U^dagger`, `A=0`.
At `x=y=0`, let `S=[[0,1],[1,0]]`; then

\[
X_x=aS,\quad X_y=\iota R,\quad
\mathcal Q=\begin{pmatrix}a^2&\iota a\\-\iota a&1\end{pmatrix},
\quad F^P_{xy}=2\iota aP_0,
\]
\[
\partial_\lambda g_{xx}=2a,\quad
\partial_\lambda\Omega_{xy}=1,\quad
\partial_\lambda F^P_{xy}=2\iota P_0.
\]

Here `X_lambda=0` at the base but `partial_x X_lambda=S`; D3 retains
that derivative and gives exactly these answers. The cut at one point
can remain fixed while its response to neighbouring comparisons changes.
The full curvature includes derivative terms: the frozen bracket of the
minimal connection alone has the opposite sign in this example.

Rational Cayley charts with the same first jets give exact finite
certification without numerical exponentials. The executable also checks
generic rational moving cuts, the Gram/orientation identity, frame
invariance, the score normalization, and the minimal-connection calibration.

## Lineage and scope

- **CID-1 and QTH-1:** existing information ledger and tensor duality,
  with their finite-model scope preserved; D2b specifies the new cut map.
- **EMK-C1:** connection covariance, curvature and commuted derivatives.
- **EMK-G1/G2/G3:** independence of arbitrary metric, recognition and
  global return quantities; no universal identification is substituted.
- **GE1 and EMK-T2:** native finite-core exponential and ordered transport.
- **QG-R1 and QG-R2:** loop readout, concentration admission and retained
  records; scalar loop defect is not silently renamed as this cut tensor.

The identities are derived within finite native cut algebra and smooth
admitted control charts. Executable algebra has native witnesses; general
differentiation and exponential/limit statements have written proofs using
the cited native completion. A physical choice of cut, reference and
generator, and a gravitational geometry/source interpretation, remain
required to turn the structural equation into a predictive field law.
