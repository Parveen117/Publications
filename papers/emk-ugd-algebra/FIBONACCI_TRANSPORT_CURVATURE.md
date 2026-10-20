# EMK-F1 — Fibonacci transport, quadratic memory and native curvature

Author: Monty. Written development: 21 October 2026.

## Result and scope

The admitted two-channel EMK algebra contains the generalized Fibonacci
transfer family. Its alternating quadratic invariant and reflection /
exponential factorization follow exactly. Coupling its traceless part
to the existing native connection calculus gives an explicit curvature
family, including both flat and nonflat cases at the Fibonacci value.

Consequently, the primitive algebra and this connection contract do not
select the golden ratio as a universal physical curvature. This note
supplies a concrete algebraic model and a parameter-selection test, not
a spacetime field equation, empirical prediction, or RH result.

The recurrence and connection coefficients below are explicitly chosen
families. They have not been derived from the lambda tower. The proofs
establish consequences of those choices. No priority is claimed for the
classical Fibonacci matrix, its quadratic form, or exponential identities.

Sources, inspected at Publications commit
b7323b72928a6300bf915c38482822c79812960a:

- [EMK-1/2 algebra](README.md): primitive relations and cut grading.
- [EMK-C1](CONNECTION_CALCULUS.md): inner derivations, full curvature,
  gauge covariance and observer faithfulness.

We use the real 2-by-2 representation of the admitted algebra. Square
roots and exponentials use this specified real extension; they are not
asserted to live in the rational-only executable certificate sector.

## F1. Transfer recurrence and its parameter

Take
\[
K=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad L=RK.
\]
Here L denotes an algebra element, not the support width of the separate
RH programme. These matrices satisfy
\[
K^2=L^2=I,\quad R^2=-I,\quad KR=-RK.
\]

For a real central constant a>0, choose
\[
Q_a=\frac a2(I+K)+L
=\begin{pmatrix}a&1\\1&0\end{pmatrix}.
\]
Then
\[
\binom{x_{n+1}}{x_n}=Q_a\binom{x_n}{x_{n-1}}
\quad\Longleftrightarrow\quad x_{n+1}=a x_n+x_{n-1}.
\]
Its characteristic identity is
\[
Q_a^2=aQ_a+I,\qquad \det Q_a=-1.
\]
Its two eigenvalues are
\[
\rho_a=\frac{a+\sqrt{a^2+4}}2>1,\qquad -\rho_a^{-1}.
\]

For strictly positive initial readings x_0,x_1, the dominant coefficient
in
\[
x_n=c_+\rho_a^n+c_-(-\rho_a^{-1})^n
\]
is \(c_+=(x_1+\rho_a^{-1}x_0)/(\rho_a+\rho_a^{-1})>0\).
Thus \(x_{n+1}/x_n\to\rho_a\). Equivalently the ratio map is
\(r\mapsto a+1/r\), whose derivative at the positive fixed point is
\(-\rho_a^{-2}\). This is asymptotic alternating contraction of ratio
errors, not contraction of the full growing state.

At a=1 the recurrence is Fibonacci and rho_a is the golden ratio phi.
At a=2 the same primitive algebra gives rho_a=1+sqrt(2).
Therefore the primitive relations alone do not force a=1.

## F2. Alternating quadratic memory

Define
\[
J_a=\begin{pmatrix}1&-a/2\\-a/2&-1\end{pmatrix},\qquad
\mathcal I_a(v)=v^TJ_av=x^2-a xy-y^2.
\]
For v=(x,y), direct substitution gives
\[
\mathcal I_a(Q_av)=(ax+y)^2-a(ax+y)x-x^2
=-\mathcal I_a(v).
\]
Hence
\[
Q_a^TJ_aQ_a=-J_a,\qquad
(Q_a^2)^TJ_aQ_a^2=J_a.
\]
The signed reading alternates at each step, its absolute value is
conserved, and the two-step map preserves the signed form. This is an
anti-isometry followed by an isometry, not loss of state information:
Q_a is invertible. It is not a universal measurement or cut law.

Since det J_a=-(1+a^2/4)<0, the form has signature (1,1).
Its two null slopes x/y are rho_a and -rho_a^{-1}.
Thus the asymptotic ratio direction is a null direction of this
particular conserved indefinite form. A null quadratic reading need
not mean a zero vector or an absent state.

A real linear change of variables puts this form into X^2-Y^2.
That is a Lorentz-type algebraic structure; identifying either
variable with a physical clock or length requires additional data.

## F3. Exact exponential factorization

Set
\[
B_a=Q_a-\frac a2I=\frac a2K+L,\qquad
d_a=1+\frac{a^2}{4},\qquad N_a=B_a/\sqrt{d_a}.
\]
The cross term vanishes because KL=-LK, giving
\(B_a^2=d_a I\) and \(N_a^2=I\).

For
\[
\eta_a=\operatorname{arsinh}(a/2)=\log\rho_a
\]
we have sinh eta_a=a/2 and cosh eta_a=sqrt(d_a). Thus
\[
\boxed{Q_a=N_a e^{\eta_aN_a},\qquad
Q_a^2=e^{2\eta_aN_a}.}
\]
Indeed \(e^{tN_a}=\cosh(t)I+\sinh(t)N_a\);
multiply by N_a and use N_a^2=I. Squaring then proves the second identity.

Also \(N_a^TJ_a+J_aN_a=0\), as direct multiplication verifies.
It follows by differentiation that
\[
(e^{tN_a})^TJ_ae^{tN_a}=J_a
\]
for every real t. At a=1 the two-step parameter is 2 log(phi).

A one-step Q_a cannot itself be the exponential of a real 2-by-2
matrix: det(exp X)=exp(tr X)>0, whereas det Q_a=-1.
The involution factor is therefore essential in this real representation.
The two-step map does admit the displayed real exponential.

This constructs continuous algebraic interpolation for a specified
two-step map. It does not calibrate physical time or derive an
electromagnetic wave.

## F4. Complete coupling to the native connection

Use precisely EMK-C1's fixed inner frame
\[
G_1=R,\quad G_2=K,\quad G_3=L,\qquad \delta_i=\operatorname{ad}_{G_i}.
\]
The nonzero independent constants are
\(c_{12}^{3}=2,\ c_{23}^{1}=-2,\ c_{31}^{2}=2\).
On the regular right module, write
\[
\nabla_i z=H_i z-zG_i,\qquad A_i=H_i-G_i.
\]
EMK-C1 gives
\[
F_{ij}=[H_i,H_j]-\sum_l c_{ij}^{\,l}H_l.
\]
The frame term is retained. The labels denote native algebraic
directions, not established spacetime coordinates.

For a fixed real central coupling t>0, choose
\[
H_1=R,\qquad H_2=tB_a,\qquad H_3=tRB_a.
\]
Both a and t are constant under the inner derivations. No variable
parameter field or omitted derivative term is being assumed.

Because RB_a=-B_aR and B_a^2=d_aI,
\[
[R,B_a]=2RB_a,\quad
[B_a,RB_a]=-2d_aR,\quad
[RB_a,R]=2B_a.
\]
For example \(B_aRB_a=-d_aR\) and \(RB_aB_a=d_aR\), proving
the middle identity. Consequently
\[
\boxed{F_{12}=0,\quad F_{31}=0,\quad
F_{23}=2(1-t^2d_a)R.}
\]

All three independent curvature components have now been calculated,
not merely the raw commutator of two generators. In particular:

- With t=1, \(F_{23}=-(a^2/2)R\).
- With \(t=d_a^{-1/2}\), all curvature components vanish for every a>0.
- At a=1, the first choice gives \(F_{23}=-R/2\); the second is flat.
- At a=2, the normalized choice is also flat, despite a non-golden ratio.

These are different specified connections on the same fixed frame.
Normalization changes the supplied connection; it is not a gauge
operation erasing the curvature of the t=1 connection.

The Bianchi identity is consistent directly: its only potentially
nonzero cyclic derivative is \(D_1F_{23}=[H_1,F_{23}]=0\);
the other components vanish and the frame terms contain F_ii=0.

For the fixed 23 direction pair, a conjugation-invariant scalar is
\[
\mathcal K_{23}:=-\frac12\operatorname{tr}(F_{23}^2)
=4(1-t^2d_a)^2.
\]
This uses R^2=-I and tr I=2. It is invariant under simultaneous
internal gauge conjugation, not under arbitrary changes of the direction
frame. It is dimensionless here, not a Ricci scalar or physical curvature
in inverse-length-squared units. It also forgets the sign of the
curvature coefficient.

## F5. Why the flat case is a genuine native gauge example

In the normalized case choose theta with
\[
\cos\theta=\frac{a}{2\sqrt{d_a}},\qquad
\sin\theta=\frac1{\sqrt{d_a}},
\]
and put \(u=e^{-\theta R/2}\). Rotation conjugation gives
\[
u^{-1}Ru=R,\quad
u^{-1}Ku=\cos\theta K+\sin\theta L=N_a,\quad
u^{-1}Lu=RN_a.
\]
Thus the normalized tuple is \(H_i=u^{-1}G_i u\).
Its connection is \(A_i=u^{-1}\delta_i(u)\), a pure gauge of the zero
connection, exactly as in EMK-C1. Its flatness is independently
confirmed by F4.

This does NOT make Q_a and Q_b gauge equivalent for arbitrary a,b.
Their traces are a and b, which conjugation preserves. The normalized
direction tuples are pure gauge; their associated discrete transfer
maps still carry different eigenvalues and step parameters.

## F6. Consequence for a physics programme

Within this model the two channel recurrence, the invariant quadratic
form, and the connection curvature are related but distinct data.
The same Q_a can be embedded in different chosen connections through t.
Therefore its golden-ratio special case does not by itself determine
curvature.

A physically usable next result must independently select:

1. The recurrence/coupling rule, including why a=1 if Fibonacci is claimed.
2. The connection normalization t or a dynamical law fixing it.
3. Which state readings, clocks and lengths realize the algebraic channels.

These are unresolved model-selection and physical-interpretation tasks,
not newly proved universal laws. Setting a curvature value and solving
for t would be parameter fitting unless that value was independently
justified. Likewise, setting a=1 in the recurrence is not a derivation
of the golden ratio from the lambda tower.

A useful subsequent mathematical comparison is to identify which of
the already specified native transport rules actually selects Q_a,
then ask whether its admitted coupling determines t. No connection
between this model and the actual RH signed source has been established.

## Validation and gate accounting

Validation in this development is written algebra: direct matrix
substitution, eigenvalue solution, exact factorization, three complete
native-curvature components, Bianchi consistency and pure-gauge
comparison. No numerical experiment, scientific computation,
certificate replay, formal proof-assistant check, or workflow dispatch
was performed. Existing certificate verdicts do not certify this note.

Closed calculation: the proposed family now has an explicit invariant,
exponential two-step transport, and complete native curvature.

Resolved limitation: primitive EMK algebra plus the displayed coupling
does not uniquely select Fibonacci or a universal nonzero curvature;
explicit alternative parameters and flat examples satisfy the same
contracts.

Open: derivation of the selected family from existing native dynamics,
physical parameter selection, and a spacetime realization. No RH,
Yang-Mills continuum, or quantum-gravity proof gate closes here.
