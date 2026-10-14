# CID-D2 — Coupled cut and transport from information–curvature feedback

14 October 2026. Native finite-carrier continuation of
[CID-D1](../native-loop-concentration/CUT_INFORMATION_DYNAMICS.md).

The state is now the pair `(P,V)`: an observation cut and its instantaneous
transport generator. The transport record `U` is retained as well. The cut
changes the information/curvature functional; its variation changes `V`;
`V` changes the cut. A conservative action is selected explicitly and its
equations are derived. This closes a mathematical model, without claiming
that the action or its dimensionless weights are forced by nature.

The new results are the feedback force, its conserved joint quantity, an
exact solution with **both P and V moving**, and a return distinction:
`P` and `V` can both return while the accumulated transport `U` does not.

## 1. Native state and the selected action

Use the finite positive coefficient carrier, composition, dagger and trace
already admitted in CID-D1. No primitive Hilbert space, physical time,
spacetime metric, or external numerical algebra engine is added. All
parameters and operators in the selected model are dimensionless.

Let

\[
P^2=P=P^\dagger,\qquad V^\dagger=-V,\qquad PV+VP=V.
\tag{1}
\]

The last condition fixes the cut-diagonal generator freedom by selecting
the minimal off-diagonal transport. Set `X_lambda=D_lambda P=[V,P]`.
Then `V=[X_lambda,P]` and `||V||²=||X_lambda||²`, where
`||Z||²=Tr(Z^dagger Z)`. The operator `V` is not the full transport record:
in the gauge `A_lambda=0`, `U'=VU`, `U(0)=1` and `P=UP(0)U^dagger`.

Choose finitely many admitted anti-dagger control operators `B_a`, with
`D_lambda B_a=0`. They specify comparison directions, not spatial axes.
Define

\[
X_a=[B_a,P],\quad W_{ab}=[X_a,X_b],\quad
C_{ab}=PW_{ab}P,\quad
\mathcal Q_{ab}=\operatorname{Tr}(PX_aX_bP)=g_{ab}+\iota\Omega_{ab}.
\tag{2}
\]

This is an exact sector of CID-D1: take `P` independent of the control
coordinates and `A_a=B_a`. Its ambient curvature is `[B_a,B_b]`; the
cut-induced increment is `C_ab=F^P_ab-P[B_a,B_b]P`. Accordingly `C` in
this paper is **the induced increment**, not the entire ambient curvature.
The mixed curvature is zero since the controls are covariantly constant
along lambda. Noncommuting control operators are allowed.

For nonnegative dimensionless `alpha,beta,w_a,h_ab`, select

\[
\begin{aligned}
\Phi(P)&=\Phi_I(P)+\Phi_C(P),\\
\Phi_I(P)&=\frac{\alpha}{2}\sum_a w_a\|[B_a,P]\|^2
          =\alpha\sum_a w_a g_{aa},\\
\Phi_C(P)&=\frac{\beta}{2}\sum_{a<b}h_{ab}\|C_{ab}\|^2,\\
\mathcal S[P]&=\int\left(\frac12\|D_\lambda P\|^2-\Phi(P)\right)d\lambda.
\end{aligned}\tag{3}
\]

The cut constraint and native positive pairing determine the allowed
variations. **The quadratic kinetic term, the functional (3), controls,
weights, and initial conditions are selected model data.** The following
equations are consequences of this selection, not a uniqueness theorem
for all possible dynamics. `lambda` orders operations; it is not identified
with a physical clock. This action fixes a parameterization and is not
asserted to be invariant under arbitrary changes of lambda.

## 2. An explicit information–curvature force

Write the first variation as `delta Phi=Tr(H delta P)`, with `H=H^dagger`.
The following ambient gradient is valid at a cut:

\[
\begin{aligned}
H_I&=-\alpha\sum_a w_a[B_a,X_a],\\
H_C&=\beta\sum_{a<b}h_{ab}\left(
 -(W_{ab}C_{ab}+C_{ab}W_{ab})
 +[B_a,[X_b,C_{ab}]]+[B_b,[C_{ab},X_a]]\right),\\
H&=H_I+H_C.
\end{aligned}\tag{4}
\]

Only the cut-off-diagonal part of `H` affects the force. Thus gradients
that differ by an operator commuting with `P` give the same cut dynamics.

**Proof.** For the information term,
`delta X_a=[B_a,delta P]` and
`Tr(X_a[B_a,delta P])=-Tr([B_a,X_a]delta P)`.
For one curvature term use `C^dagger=-C` and

\[
\delta C=\delta P\,WP+PW\,\delta P
 +P\big([\delta X_a,X_b]+[X_a,\delta X_b]\big)P.
\]

Its variation is `-beta h Tr(C delta C)`. Since `PC=CP=C`, moving
the factors cyclically gives the first term of `H_C`; moving each
commutator off `delta X` and then off `delta P` gives its other two terms.
This proves (4) in every finite admitted carrier. No derivatives of a
numerically estimated energy are used.

## 3. The coupled equation and conserved total

Stationarity of (3), with cut-tangent endpoint-fixed variations, gives

\[
\boxed{\begin{aligned}
D_\lambda P&=[V,P],\\
D_\lambda V&=[P,H_I(P)+H_C(P)],\\
D_\lambda U&=VU.
\end{aligned}}\tag{5}
\]

In the third equation `U` maps the fixed initial carrier into the current
one; in a variable frame its ordinary equation is
`partial_lambda U=(V-A_lambda)U`. The first two equations constitute the
closed finite model. The third retains its accumulated transport.

**Variational proof.** Write `delta P=[eta,P]` with anti-dagger `eta`.
Integration by parts, trace cyclicity and fixed endpoints give
`[P,D_lambda²P+H]=0`. Differentiating
`V=[D_lambda P,P]` gives `D_lambda V=[D_lambda²P,P]=[P,H]`.
Conversely (5), with (1), gives that stationary equation. This uses
only the inherited scalar differentiation/integration and native products.

**Constraint proof.** Differentiating `P²-P` gives zero on (5).
Dagger conditions are preserved. Differentiating `PV+VP-V` gives
`[V²,P]`, which is zero when `PV+VP=V`. Thus no constraint-restoring
projection is needed by the continuous equation. The rank of `P` is
preserved by unitary conjugation. The unbounded-carrier analogue is not
asserted here.

**Conserved quantity.** If the controls and weights are covariantly
constant along lambda, then

\[
\boxed{\begin{aligned}
\mathcal E
 &=\frac12\|V\|^2+\Phi_I(P)+\Phi_C(P)\\
 &=g_{\lambda\lambda}
   +\alpha\sum_a w_a g_{aa}
   +\frac{\beta}{2}\sum_{a<b}h_{ab}\|C_{ab}\|^2,
\qquad \partial_\lambda\mathcal E=0.
\end{aligned}}\tag{6}
\]

Indeed `partial_lambda(||V||²/2)=-Tr(V[P,H])` and
`partial_lambda Phi=Tr(H[V,P])`; their sum is zero by trace cyclicity.
The identity `||V||²/2=g_lambda,lambda` follows from the off-diagonal
cut blocks. If the controls or weights have explicit lambda-dependence,
their explicit variation is a source term in the energy balance and
conservation is not asserted.

In the rank-one two-mode sector, `C_ab=2 iota Omega_ab P`, so (6) becomes

\[
\boxed{\mathcal E=g_{\lambda\lambda}
 +\alpha\sum_a w_a g_{aa}
 +2\beta\sum_{a<b}h_{ab}\Omega_{ab}^{,2}.}\tag{7}
\]

This is a conserved combination of **transport information, directional
information, and squared induced curvature**. It is not a claim that
each summand, Shannon entropy, or every information statistic is constant.
The continuous system is reversible under `(lambda,V)->(-lambda,-V)`.

## 4. Information and curvature now feed back into their own evolution

CID-D1's transport identity is reused with
`Y_lambda,a=[B_a,[V,P]]`:

\[
\boxed{\partial_\lambda(g_{ab}+\iota\Omega_{ab})=
\operatorname{Tr}\left(P\left([B_a,[V,P]]X_b
 +X_a[B_b,[V,P]]\right)P\right).}\tag{8}
\]

Equations (4), (5) and (8) form the requested feedback: `Q` gives the
selected scalar functional; its variation changes `V`; the updated `V`
changes `P` and hence `Q`. `V` no longer needs to be supplied as a whole
external path. Its **initial value** is still required, just as the
initial cut is required.

The second-order cut equation makes the two contributions explicit:

\[
D_\lambda^2P=[[P,H],P]+[V,[V,P]].\tag{9}
\]

The first term is the selected information–curvature force. The second
is the normal acceleration required to stay on the cut constraint.
Freezing `V` deletes the first term and fails whenever `[P,H]` is nonzero.

Under a lambda-dependent unit frame `G`, transform `P,V,B_a,H` by
`Z'=G^dagger ZG` and
`A'_lambda=G^dagger A_lambda G+G^dagger partial_lambda G`. Then (4)–(9)
are covariant and their scalar readings are invariant. `U` transforms
as `U'=G(lambda)^dagger U G(0)`. A zero connection in one frame is not
silently held zero in another variable frame.

## 5. Exact two-strand solution with an evolving generator

Use the inherited native matrices

\[
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
J=\begin{pmatrix}0&\iota\\\iota&0\end{pmatrix},\quad
L=RJ=\begin{pmatrix}-\iota&0\\0&\iota\end{pmatrix}.
\]

Set `B=(R,J,L)`, `alpha=2`, `beta=1`, `w=(1,1,2)` and
`h_RJ=1,h_RL=h_JL=2`. For a general rank-one cut

\[
P=\frac12\begin{pmatrix}1+z&x-\iota y\\x+\iota y&1-z\end{pmatrix},
\qquad x^2+y^2+z^2=1,
\]

native multiplication gives

\[
\Phi_I=6-2z^2,\quad\Phi_C=4-2z^2,\quad\Phi=10-4z^2,
\quad [P,H]=-4[P,[L,[L,P]]].\tag{10}
\]

Both the information and the curvature terms contribute half the force.
With `W=2L`, the following solves the full selected system:

\[
\boxed{P(\lambda)=e^{\lambda W}P_*e^{-\lambda W},\qquad
V(\lambda)=[[W,P(\lambda)],P(\lambda)].}\tag{11}
\]

Here exponentials mean the inherited native `Exp_Sigma`. The first
equation follows from cut tangency. Differentiation gives
`V'=[W,V]=-4[P,[L,[L,P]]]=[P,H]`, proving the second equation.
Except at special latitudes, **V itself is changing**, even though the
auxiliary rotating-frame operator `W` used to write the solution is constant.

Choose the exact cut

\[
P_*=\begin{pmatrix}4/5&2/5\\2/5&1/5\end{pmatrix},\qquad z=3/5.
\]

Then

\[
P_{12}(\lambda)=\frac25 e^{-4\iota\lambda},\qquad
V_{12}(\lambda)=\frac{24\iota}{25}e^{-4\iota\lambda},\qquad
V_{11}=-\frac{32\iota}{25}=-V_{22}.\tag{12}
\]

The cut and generator records therefore have the same rotating phase and
a fixed quarter-turn offset. Appending the unwrapped operation coordinate
lambda to their real/turn components produces two helical records. This
is an exact interpretation of the two-strand intuition in the selected
solution; it neither imposes a three-dimensional physical helix nor proves
that every coupled solution is helical.

The values are

\[
g_{\lambda\lambda}=64/25,\quad\Phi_I=132/25,\quad
\Phi_C=82/25,\quad\mathcal E=278/25.\tag{13}
\]

Although these totals are constant on this orbit, individual directional
information readings and individual curvature components change. The
certificate checks this directly through (8), including nonzero derivatives
of both kinds. The general force can also exchange the summands; a separate
rational state in the certificate has nonzero `partial_lambda Phi` with
the opposite kinetic derivative.

## 6. A returning cut and generator can retain transport memory

The ordered transport solving `U'=V(lambda)U` for (11) is exactly

\[
U(\lambda)=e^{\lambda W}e^{\lambda(V_*-W)}.\tag{14}
\]

This is verified by differentiating the two factors in order.
`V_*-W` commutes with `P_*`, so (14) transports the same cut as (11).
For (12), `V_*-W=(6 iota/5)(2P_*-1)`. The first cut/generator period is
`T=pi/2`, yet

\[
\boxed{P(T)=P_*,\quad V(T)=V_*,\quad
U(T)=-\exp\!\left(\frac{3\pi\iota}{5}(2P_*-1)\right)\ne1.}\tag{15}
\]

Let `q=Exp_Sigma(2 pi iota/5)`. The return operator is

\[
U(T)=q^{-1}P_*+q(1-P_*),\quad U(T)^5=1,
\quad U(T)^n\ne1\ (n=1,2,3,4).\tag{16}
\]

The native certificate proves this using the exact fifth-root polynomial,
not decimal sine/cosine or a sampled period. Thus the transported carrier
returns after five visible cycles in this example. The phase is a full
carrier transport record; a measurement using only the cut is blind to it.
No sheet-charge increment is inferred from this finite phase return.

EMK-G3 already separates visible return from lifted memory. The new step
here is a noncommuting, dynamically generated example in which **both the
cut and its current generator return before the accumulated transport**.
Resetting the transport at every visible cycle would erase that record.

## 7. Certification and limits

The unchanged native `F`, `Cut`, `Expr` and matrix operations execute all
algebra. Python 3.12 only orchestrates source integrity and the recursive
CID-D1/QG-R2/QG-R1/phase/six-packet replay. The new certificate contains:

- dimension-independent word identities for constraints, acceleration,
  cyclic-trace variation and energy cancellation;
- a native polynomial presentation of the rank-one sphere, proving the
  two force contributions, orbit and conserved values symbolically;
- exact rational cases checking the explicit gradient against direct
  variation, scalar covariance, tensor derivatives and energy exchange;
- exact fifth-root transport words, with nonreturn controls.

The variational calculus, continuous ordered-exponential solution and
helical interpretation have written proofs. Finite exact checks are not
advertised as formal proofs of every analytic assertion. No numerical ODE
simulation or floating-point fit supplies any claimed identity.

The remaining physical question is why a particular action, comparison
operators, weights and clock/measurement adapter describe a real system.
This paper provides a closed reversible information–curvature model and
an exact feedback example; it does not establish Einstein's equations,
quantum gravity, a singularity theorem, or universal gravitational constants.

For external context, projector-based information metrics and curvature
already belong to established projector geometry; they are not claimed
as a new general discovery here. See Mitscherling, Avdoshkin and Moore,
[Phys. Rev. B 112, 085104 (2025)](https://doi.org/10.1103/qscv-qxqt).
That reference supplies context, not an imported proof or execution engine.
