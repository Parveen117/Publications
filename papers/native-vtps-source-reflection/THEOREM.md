# VSR-1 — Actual five-reading source in VTPS coordinates

17 October 2026. Continue the unchanged TLM/TET source from Publications
commit 0c98adf69682f6fee48be85dcbcb7667d77dea2e. This packet realizes the
user's two-sided tower proposal on the actual constitutive example whose
visible reading at the cut is 5. It does not assume that a numeral determines
a geometry or that every constitutive source has a reflection symmetry.

## VSR-T1. Exact source reflection about the calibrated cut

Retain TLM's source and its deformation parameter:
\[
u_\mu(s,v)=10s-10v+s^2+sv+\tfrac52v^2+\tfrac\mu2s^2v,\quad
S=1+s,\quad V=1+v,\quad T=u_s,\quad P=-u_v.
\]
Put
\[
a=2+\mu v,\quad b=1+\mu s,\quad
H_\mu=\begin{pmatrix}a&b\\b&5\end{pmatrix},\quad
\Delta=5a-b^2.
\]
For |\mu|<=5, |s|,|v|<1/20, TLM gives a>=7/4 and
Delta>=115/16>0. The involution is
\[
(\mu,s,v)\longmapsto(-\mu,-s,-v).
\]
The exact identities are
\[
u_{-\mu}(s,v)=u_\mu(-s,-v)+20s-20v,\qquad
H_{-\mu}(s,v)=H_\mu(-s,-v),
\]
\[
T_{-\mu}(s,v)=20-T_\mu(-s,-v),\qquad
P_{-\mu}(s,v)=20-P_\mu(-s,-v).
\tag{1}
\]

**Proof.** Substitute (-s,-v) into the cubic polynomial. The quadratic
part is even, the cubic term changes sign, and the displayed affine term
restores the calibrated linear part. Differentiate once and twice. The
regular source domain is preserved by the paired parameter/coordinate
reflection. QED.

Here mu is the inherited signed deformation control. It may label the two
directions of this example, but it is not identified with any individual
thermodynamic lambda response. In particular, mu=0 is a flat source in this
family while its six thermodynamic response values at the cut are nonzero.
The derivative depth N remains nonnegative.

## VSR-T2. Actual VTPS chart and Maxwell response matrix

Fix the order **VTPS=(V,T,P,S)**. Use controls (v,theta), theta=T-10, with
responses (P,S). The exact inverse chart is
\[
\theta=v+(2+\mu v)s,\qquad
s=\frac{\theta-v}{a},\qquad a=2+\mu v.
\tag{2}
\]
All quantities b and Delta below are evaluated at this s. The chart domain
is the image of the stable source patch in T1. The reflection becomes
(mu,v,theta)->(-mu,-v,-theta).

The same source's Helmholtz potential is
\[
F_\mu(v,\theta)=u_\mu-(10+\theta)(1+s)
=-10-\theta-10v+\tfrac52v^2-\frac{(\theta-v)^2}{2a}.
\tag{3}
\]
Consequently dF=-P dv-S dtheta and the VTPS response matrix is
\[
\boxed{\mathsf M_\mu=
\frac{\partial(P,S)}{\partial(v,\theta)}
=\frac1a\begin{pmatrix}-\Delta&-b\\-b&1\end{pmatrix}.}
\tag{4}
\]
Its off-diagonal equality is precisely
(∂P/∂T)_V=(∂S/∂V)_T. Its determinant is -5/a<0.

**Proof.** Equation (2) follows from T=10+2s+v+mu s v. Substitute it into
u-TS for (3). Equivalently ds=(-b/a)dv+(1/a)dtheta and
dP=-b ds-5dv, which immediately gives (4) and its determinant.
The same result follows from -D²F. QED.

This is a Legendre chart of the actual source, not just a permutation of
the original Hessian. The response matrix M is indefinite and must not be
used as the acquired positive metric. A compass ordering alone does not
supply the source reflection (1).

Define the centered potential
\[
\Phi_\mu=F_\mu+10+\theta+10v.
\]
Then
\[
\Phi_{-\mu}(y)=\Phi_\mu(-y),\quad
\binom{P_{-\mu}(y)-10}{S_{-\mu}(y)-1}
=-\binom{P_\mu(-y)-10}{S_\mu(-y)-1},\quad
\mathsf M_{-\mu}(y)=\mathsf M_\mu(-y).
\tag{5}
\]
These follow directly from (3), including the affine calibrations.

## VSR-T3. Diagonal positive metric with retained nonzero curvature

Let x=(s,v), y=(v,theta). Its coordinate Jacobian and positive pullback
metric are
\[
Q=\frac{\partial x}{\partial y}
=\begin{pmatrix}-b/a&1/a\\1&0\end{pmatrix},\qquad
\boxed{g_\mu=Q^T H_\mu Q
=\begin{pmatrix}\Delta/a&0\\0&1/a\end{pmatrix}.}
\tag{6}
\]
In particular det Q=-1/a: the coordinate orientation changes.

The old connection A_alpha=(1/2)H^(-1)∂_alpha H is the Levi-Civita
connection of H in the original Hessian chart. Its coordinate transform is
\[
\boxed{\Gamma_i=
Q^{-1}\Bigl(\sum_\alpha A_\alpha Q^\alpha{}_i\Bigr)Q
+Q^{-1}\partial_i Q.}
\tag{7}
\]
It is also obtained directly from g by the metric Levi-Civita formula.
The curvature transforms as
\[
\mathcal F^y_{v\theta}
=-\frac1a Q^{-1}\mathcal F^x_{sv}Q,\qquad
\mathcal F^x_{sv}=\frac{\mu^2}{4\Delta^2}
\begin{pmatrix}5b&25\\-5a&-5b\end{pmatrix}.
\tag{8}
\]
The retained convention gives the scalar
\[
\boxed{\kappa_\mu=\frac{(g\mathcal F^y_{v\theta})_{12}}{\det g}
=\frac{5\mu^2}{4\Delta^2}.}
\tag{9}
\]

**Proof.** Multiply HQ and Q^THQ for (6). The connection product rule
applied to a vector w_x=Qw_y gives (7); uniqueness of the torsion-free,
metric-compatible connection identifies the direct metric formula.
Curvature is an endomorphism-valued two-form, giving the determinant
factor and conjugation in (8). Substitute the inherited NT/TLM formula.
For clarity, (8) simplifies to
\[
\mathcal F^y_{v\theta}
=\frac{5\mu^2}{4a\Delta^2}
\begin{pmatrix}0&1\\-\Delta&0\end{pmatrix};
\]
multiplying by (6) and dividing by det g=Delta/a² gives (9). QED.

Thus a diagonal metric can retain nonzero curvature. Reselecting
(1/2)g^(-1)dg in the new chart would give a diagonal, flat connection:
each diagonal entry is one half an exact logarithmic derivative.
It omits the second term in (7), and in this example does not reproduce
the source connection. Dropping Q^(-1)dQ is also wrong even when the
other term is retained. These failures are tested for nonzero mu.

## VSR-T4. Both sides of the actual five-reading matrix

At the cut y=0,
\[
g_0=\operatorname{diag}(9/2,1/2),\quad
\mathsf M_0=\begin{pmatrix}-9/2&-1/2\\-1/2&1/2\end{pmatrix},
\]
\[
g_v(0)=\frac\mu4\operatorname{diag}(3,-1),\qquad
g_\theta(0)=\frac\mu2\operatorname{diag}(-1,0).
\tag{10}
\]
In the first-layer ordering (1,v,theta), the single matrix record is
\[
\mathcal G_\mu=
\begin{pmatrix}g_0&0&0\\g_v(0)&g_0&0\\g_\theta(0)&0&g_0\end{pmatrix}.
\tag{11}
\]
With D_1=diag(I_2,-I_2,-I_2),
\[
\boxed{\mathcal G_{-\mu}=D_1\mathcal G_\mu D_1,\qquad D_1^2=I.}
\tag{12}
\]
The original TLM six-by-six H record obeys the same relation.

More generally (5)-(6) give, for every multiindex alpha,
\[
\partial^\alpha g_{-\mu}(0)=(-1)^{|\alpha|}
\partial^\alpha g_\mu(0),\qquad
\partial^\alpha\mathsf M_{-\mu}(0)=(-1)^{|\alpha|}
\partial^\alpha\mathsf M_\mu(0).
\]
For any finite jet depth N, let D_N act by (-1)^|alpha| on each
monomial block. Then
\[
\boxed{\operatorname{Reg}_N(g_{-\mu})
=D_N\operatorname{Reg}_N(g_\mu)D_N,}
\tag{13}
\]
and likewise for M. The characteristic polynomials agree. If R_-=D R_+ D,
x_-=D x_+, and O_-=O_+D, then O_-R_-x_-=O_+R_+x_+.

**Proof.** Differentiate the source reflection at the fixed point.
A regular matrix block is a coefficient with degree difference
|alpha|-|beta|. Its parity equals |alpha|+|beta| modulo two, exactly the
left/right D factors. The observer equality follows from D²=I. QED.

This is the user's same-result rule with its actual source and observer
transport supplied: even layers agree and odd layers retain orientation.
The regular representation remains TLM's algebra encoding, not a faithful
positive self-adjoint quantum representation.

The actual mechanical reading is r=-lambda_S=Vc=5(1+v). Therefore
\[
r_\mu(0)=5,\qquad
\kappa_{+5}(0)=\kappa_{-5}(0)=125/324.
\tag{14}
\]
At reflected off-cut points, r_mu(v,theta)=5(1+v) and
r_-mu(-v,-theta)=5(1-v). Their calibrated readings r/V both equal 5.
Uncarried raw readings are not generally equal. The certificate retains
three exact off-cut mirror pairs to make this boundary visible.

## VSR-T5. Connection and transport across the reflection

The exact metric identity g_-mu(y)=g_mu(-y) says that y->-y is an isometry
between the paired source patches. The Levi-Civita connection satisfies
Gamma_-mu,i(y)=-Gamma_mu,i(-y). For a path gamma in the negative-parameter
patch, the transport equation therefore identifies its matrix transport
with that on the reflected path -gamma in the positive-parameter patch
in the corresponding coordinate frames. Co-carried observer pairings agree.

**Proof.** Insert the metric reflection in the Levi-Civita formula; its
first derivatives acquire one minus sign. In the path equation, the
derivative of -gamma supplies the same minus sign. Uniqueness gives the
transport identity. The differential of the isometry is -I, whose
endpoint conjugations cancel in the transport matrix. QED.

Reflecting a path about the cut is not generally traversal of the same
path backwards; no unsupported U(-parameter)=U(parameter)^(-1) statement
for arbitrary varying connections is inferred.

## Scope and continuation

Written proofs establish the all-finite-depth and path statements.
GitHub native fixtures verify actual source differentiation, both connection
constructions, curvature, graded matrices, observer pairing and off-cut
controls. Passing predecessor execution is reused only after immutable-run
and source-byte checks. No local scientific computation or audit agent is used.

The source-specific two-sided tower and VTPS realization are now explicit.
Extension to another constitutive source requires its own cut-equivariance,
and physical interpretation still requires observable preparation and units.
