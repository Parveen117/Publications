# CID-D3 — Recover the information–curvature coupling ratio from a cut

14 October 2026. Continuation of
[CID-D2](../native-cut-transport-feedback/THEOREM.md).

CID-D2 selected a conservative action and derived coupled cut/transport
dynamics. This paper determines which of its coefficients an observation
can identify. The symmetric helical protocol has an exact ambiguity:
**even its complete cut, generator and transport histories cannot separate
the information and curvature couplings.** A calibrated change in one
control operator removes that ambiguity. At an initially stationary cut,
two scalar readings of its covariant acceleration recover the coupling
ratio without a physical clock or a known speed of the operation parameter.

This is an inverse theorem within the selected action family. It does not
derive a universal action or infer a measured constant from synthetic data.

## R1. The source family and an exact ambiguity

Keep the finite native positive carrier, anti-dagger controls `B_a`,
Hermitian cut `P`, minimal anti-dagger generator `V`, and retained transport
`U` of CID-D2. Write its action functional with unknown nonnegative
coefficients:

\[
\Phi(P)=\alpha I(P)+\beta C(P),\qquad
I=\frac12\sum_a w_a\|[B_a,P]\|^2,\quad
C=\frac12\sum_{a<b}h_{ab}\|P[[B_a,P],[B_b,P]]P\|^2.
\tag{1}
\]

The controls, their relative calibration, weights and action family are
known. The constants `alpha,beta` are the coefficients to identify; they
are not the electromagnetic fine-structure constant or a claimed
gravitational constant. Define the unit-coefficient gradients and forces

\[
\delta I=\operatorname{Tr}(H_I^0\delta P),\quad
\delta C=\operatorname{Tr}(H_C^0\delta P),\qquad
T_I=[P,H_I^0],\quad T_C=[P,H_C^0].
\]

Their explicit formulas are inherited from CID-D2, with respectively
`(alpha,beta)=(1,0)` and `(0,1)`. Thus

\[
D_\lambda P=[V,P],\qquad
D_\lambda V=\alpha T_I(P)+\beta T_C(P).\tag{2}
\]

For the symmetric protocol `B=(R,J,L)`, `w=(1,1,2)` and
`h_RJ=1,h_RL=h_JL=2`, every rank-one two-mode cut

\[
P=\frac12\begin{pmatrix}1+z&x-\iota y\\x+\iota y&1-z\end{pmatrix},
\quad x^2+y^2+z^2=1
\tag{3}
\]

satisfies

\[
I=3-z^2,\qquad C=4-2z^2=2I-2,\qquad T_C=2T_I.
\tag{4}
\]

Consequently only `kappa=alpha+2 beta` occurs in (2). Any two coefficient
pairs with the same kappa define the same vector field on `(P,V)`, and
therefore the same history wherever the smooth finite flow exists, for
the same initial data. Since `U'=VU`, the full transport is also the same.
Adding its holonomy record does not repair this particular ambiguity.
This is equality of the full equations, not an inference from finitely
many samples of the helical orbit.

For example `(alpha,beta)=(2,1),(3,1/2),(4,0),(0,2)` all have kappa 4.
They reproduce CID-D2's same cut/generator orbit and five-cycle transport
return. Their potentials differ by additive constants. An independently
calibrated absolute potential reading would be additional observation
data; it is not part of the cut/transport protocol considered here.

## R2. Native force readings and coefficient recovery

For anti-dagger operators use the full native positive pairing

\[
\langle A,B\rangle=-\operatorname{ReTr}(AB).
\]

This has no extra factor of one half. At a declared cut form

\[
G=\begin{pmatrix}
\langle T_I,T_I\rangle&\langle T_I,T_C\rangle\\
\langle T_I,T_C\rangle&\langle T_C,T_C\rangle
\end{pmatrix},\qquad \Delta=G_{11}G_{22}-G_{12}^2.
\tag{5}
\]

Let `A` be the observed generator acceleration in a calibrated affine
operation parameter and set `b_i=<T_i,A>`. If `Delta>0`, the force
directions are independent. Provided `A` is in their span, the unique
coefficients are

\[
\widehat\alpha=\frac{G_{22}b_1-G_{12}b_2}{\Delta},\qquad
\widehat\beta=\frac{G_{11}b_2-G_{12}b_1}{\Delta}.
\tag{6}
\]

The span condition is checked by the native residual
`A-alpha_hat T_I-beta_hat T_C`. It must be zero. In a rank-one two-mode
cut the real horizontal anti-dagger space has dimension two, so two
independent forces span it. In that sector a single two-reading fit alone
cannot validate the chosen action family; independent preparations are
needed for a meaningful prediction test.

**Proof.** Pair `A=alpha T_I+beta T_C` with both forces and eliminate
the two unknowns. The two exact numerator identities are
`G22 b1-G12 b2=Delta alpha` and
`G11 b2-G12 b1=Delta beta`. Positivity gives independence when
`Delta>0`. This is native finite Gram elimination already used in GE4;
the specific force signatures and their ambiguity are the new target.

## R3. Initial rest removes operation-speed calibration

Prepare `D_lambda P=0`, equivalently the minimal `V=0`, at the
observation. Let `lambda=f(tau)` be any twice-differentiable local
reparameterization with `f'(tau_0) != 0`, and use the transformed
reference connection in covariant derivatives. Then

\[
D_\tau^2P=f'^2D_\lambda^2P+f''D_\lambda P
\quad\Longrightarrow\quad
\mathcal A:=[D_\tau^2P,P]
=f'^2(\alpha T_I+\beta T_C)
\quad\text{at initial rest}.\tag{7}
\]

To see the last equality, CID-D2 gives
`D_lambda²P=[[P,H],P]` at rest and its cut identity recovers
`[[[P,H],P],P]=[P,H]`.

Apply (6) to `A=mathcal A`. The recovered pair is
`(alpha_hat,beta_hat)=f'^2(alpha,beta)`. The positive unknown multiplier
cancels in

\[
\boxed{\rho=\frac{\beta}{\alpha}
=\frac{G_{11}b_2-G_{12}b_1}{G_{22}b_1-G_{12}b_2},
\qquad b_i=\langle T_i,[D_\tau^2P,P]\rangle.}\tag{8}
\]

For finite rho require `alpha>0`; the denominator in (8) is then
positive. The pure-curvature branch `alpha=0,beta>0` is recorded as
`[alpha_hat:beta_hat]=[0:1]`, not evaluated by division. A zero response
contains no nonzero coupling ratio. Negative recovered coefficients
are incompatible with the declared nonnegative action family at rest.

The quantities in (5), (8) are unit-frame invariant: transform the cut,
controls and covariant acceleration together, and use trace cyclicity.
The ratio is dimensionless. No estimate of `f'`, `f''` or physical seconds
is required. Initial rest and a regular parameter are load-bearing:
away from rest, (7) has the additional term `f'' V`, which can rotate
the apparent force and corrupt the inferred ratio.

The first cut derivative is zero in every initial-rest experiment in this
family. It cannot separate the coefficients. The second covariant cut
derivative supplies the sufficient observation under `Delta>0`.
Reconstructing this derivative from noisy finite measurements is a
separate experimental task; the exact theorem does not provide a noise
model or declare an arbitrarily small determinant usable.

## R4. A calibrated control deformation that removes the ambiguity

Use the same weights, and change one control physically within the
declared protocol:

\[
B=(sR,J,L),\qquad s\ne0,\qquad u=s^2.
\tag{9}
\]

This is a relative control-amplitude change with the same known weights
and couplings. It is not a unit-frame transformation or a relabeling of
control units; such a relabeling would require transforming the weights
and would not manufacture new observations.

For (3), native multiplication now gives

\[
I=u+3-x^2-u y^2-2z^2,\qquad
C=4u x^2+4y^2+2u z^2.\tag{10}
\]

The force Gram determinant factors exactly:

\[
\boxed{\Delta=256(u-1)^2(u-2)^2x^2y^2z^2.}\tag{11}
\]

Thus `u` must avoid **both** 1 and 2, and all three cut coordinates
must be nonzero. Asymmetry alone is insufficient: at `u=2`, the other
exact dependence `C=4I-8` gives `T_C=4T_I`.

**Derivation of (11).** For checking the native two-mode matrix products,
write their real coefficient triples as `n=(x,y,z)`,
`d=(x,u y,2z)`, `e=(4u x,4y,2u z)`. The two anti-dagger force triples
are proportional to `-2 n cross d` and `2 n cross e`; the native trace
pairing gives a factor 8 on their squared vector pairings. Their Gram
determinant is `64 (n dot (d cross e))²` because `n dot n=1`.
The triple product is `-2(u-1)(u-2)xyz`, proving (11). These are
coefficient identities of the native matrices, not additional physical
three-space primitives. The certificate also computes the gradients
from (1), verifies these force forms, and proves (11) in the native
polynomial quotient `x²+y²+z²=1`.

## R5. Exact reconstruction example

Choose `s=2` and prepare at rest

\[
(x,y,z)=(2/3,2/3,1/3),\qquad
P_* =\begin{pmatrix}2/3&(1-\iota)/3\\(1+\iota)/3&1/3\end{pmatrix}.
\]

Then

\[
G=\frac1{81}\begin{pmatrix}1312&5120\\5120&20992\end{pmatrix},
\qquad\Delta=\frac{16384}{81}>0.\tag{12}
\]

Generate a synthetic exact response with the CID-D2 example's
`alpha=2,beta=1`, and let its uncalibrated parameter multiplier be any
`c=f'^2>0`. The two readings are

\[
b=\frac{c}{81}\begin{pmatrix}7744\\31232\end{pmatrix}.
\]

Equation (8) returns **rho=1/2** for every such c. This is a reconstruction
fixture, not a measurement of nature. The different formerly indistinguishable
pair `(3,1/2)` returns rho `1/6` under this deformed protocol.
The pure-information and pure-curvature endpoints are separated as well.

In particular the previously chosen ratio is now **testable**: it is no
longer protected by the symmetric protocol's ambiguity. Its truth in an
actual system requires its measured initial cut response.

## R6. A prediction check on independent preparations

After reconstructing a finite rho on one preparation, retain it and
prepare another allowed cut at rest. Compute
`F_rho=T_I+rho T_C` using its calibrated controls. For a nonzero response
`A` and nonzero `F_rho`, use the dimensionless diagnostic

\[
\mathcal D_\rho=
1-\frac{\langle A,F_\rho\rangle^2}
{\|A\|^2\|F_\rho\|^2},\qquad
\langle A,F_\rho\rangle>0.\tag{13}
\]

The model predicts `D_rho=0` and the displayed positive orientation, even
if the second experiment has a different regular operation speed. More
sharply, if `A=c(alpha T_I+beta T_C)`, native expansion gives

\[
\boxed{\|A\|^2\|F_\rho\|^2-\langle A,F_\rho\rangle^2
=c^2(\alpha\rho-\beta)^2\Delta.}\tag{14}
\]

At an independent preparation with `Delta>0` and `alpha>0`, this
vanishes exactly at the correct ratio. A negative collinear response
also makes the squared diagnostic zero, so the orientation test cannot
be dropped. Zero denominators and a blind preparation are not accepted
as successful tests.

The certificate learns rho on (12), predicts responses at other rational
cuts and control amplitudes without refitting, and rejects a wrong ratio.
These remain synthetic exact consistency tests. Passing one or several
such tests on real data would constrain the proposed action family; it
would not by itself prove a universal gravitational law.

## Provenance and claim coverage

CID-D1 supplies the information/induced-curvature tensor. CID-D2 supplies
the explicit gradient, dynamics, force extraction, and conserved total.
GE4 supplies the established native Gram-reconstruction discipline and
its observability boundary. We do not rebrand general inverse-problem or
Gram-inversion theory as new. The specific contributions are (4), the
initial-rest parameter cancellation (7)–(8), the factorized protocol
criterion (11), and the predictive residual (14) for the coupled model.

All executable mathematics uses the unchanged RKF `F`, `Cut`, `Expr`
and matrix engine. Python 3.12 only orchestrates pinned source integrity
and recursively replays CID-D2, CID-D1/QG-R2, QG-R1, phase certification
and the six original packets. Native algebra, exact rational cases and
negative controls are separated from the written chain-rule and
flow-equivalence arguments.

This continuation also corrects the literal comma in CID-D2 equation (7)'s
printed squared-curvature exponent to `Omega_ab^2`. Its mathematics, code,
39 certificates, 150 checks and 11 controls are unchanged; its source pin
and deterministic report are regenerated for the corrected document.
