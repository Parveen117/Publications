# CEM-1 — Native cut equilibrium and source memory

17 October 2026. Written continuation at Publications commit
c8e8da9b42f997b857bf6820a3e575603e51e22a, following
[TPW-1](RH_TWO_PRIME_WINDOW.md) and the native cut algebra of
[OQD-1](RH_OBSERVER_QUOTIENT_DESCENT.md) and
[RFT-1](RH_REFLECTION_TRANSFER.md).

A cut can produce a stationary effective evolution even when
the ambient generator has no nonzero stationary polynomial state.
For the retained native diagram this statement has an exact
realization, without changing its source or assuming positivity.

Take the antisymmetric part of the native wave,
A=(W-W^(-1))/2, and its specified exponential continuation U(t)=e^(tA).
Then PAP=0 although AJf is nonzero for every nonzero cut profile f.
Repeated cutting freezes the retained channel:
\[
\bigl[P\,U(t/n)\,P\bigr]^n\longrightarrow P.
\]
The information lost from that limit has an exact rescaling.
If R_(a,n)(t;f,h) is its source return reading, then
\[
\boxed{
-\frac n4\,\partial_t^2 R_{a,n}(0;f,h)
=\mathcal M_a(f,h):=\mathcal L_a[y(1-y)fh]
\quad\hbox{for every integer }n\ge1.
}
\]
This same form is the source pairing of the transverse velocities,
divided by four. It is invariant when the cut and profiles are
transported together. It has zero radical. Its positivity on
all real polynomial squares is equivalent to RH.

Thus the frozen limit alone cannot prove RH: it holds regardless
of the sign of the source memory. The rescaled response retains
the sign question that the unscaled limit discards.

**Status: RH remains unproved.** The results below are written
analytic proofs in the stated source algebra. The Zeno mechanism
and general compression ideas are classical; no general novelty
or priority claim is made for them. The contribution here is the
explicit native generator, its exact cut-depth rescaling, and
the identification of the retained source positivity condition.
No physical mass, consciousness variable, thermodynamic law,
or Yang--Mills conclusion is inferred.

## Retained source, cut and parameter conventions

Fix any finite a>=a_0=3249/4. Retain the actual source of OQD/RFT:
\[
\begin{gathered}
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\\
\mu=(\gamma-i\delta)^2,\qquad \delta=\beta-\tfrac12,\qquad
u_\mu=\frac a{a+\mu},\\
w_\mu=\frac{\nu_\mu}{q(a)(a+\mu)},\qquad
\mathcal L_a[f]=\sum_\mu w_\mu f(u_\mu),\qquad
B_a(f,h)=\mathcal L_a[fh].
\end{gathered}
\tag{1}
\]
The sum is over distinct source values, with multiplicity nu_mu.
The retained unconditional bounds and properties are
\[
q(a)>0,\quad 0<|u_\mu|<1,\quad u_\mu\ne1,\quad
\sum_\mu|w_\mu|\le2,\quad \mathcal L_a[1]=1.
\tag{2}
\]
The infinite support is conjugation symmetric, has only 0 as an
accumulation point, and all weights are nonzero. OQD-T2 proves
rad B_a={0}. Neither the weights nor B_a are assumed positive.

Set V=R[y], E=V direct-sum V, and
\[
\begin{gathered}
g(y)=\operatorname{diag}(y,1-y),\qquad
P(y)=\begin{pmatrix}y&1-y\\y&1-y\end{pmatrix},\\
Z=\operatorname{diag}(1,-1),\qquad
W=(2P-I)Z
=\begin{pmatrix}2y-1&-2(1-y)\\2y&2y-1\end{pmatrix},\\
\mathcal G_a(x,z)=\mathcal L_a[x^Tg z],\qquad
Jf=(f,f).
\end{gathered}
\tag{3}
\]
These are bilinear source products, even when u_mu is nonreal.
They are not replaced by squared absolute values. Pointwise,
P^2=P, P^Tg=gP, and W^TgW=g.

The following parameters are kept distinct.

| Symbol | Meaning |
|---|---|
| a | Fixed finite actual-source anchor |
| t | Real parameter of the exponential continuation specified below |
| n | Number of repeated cuts during a fixed total t |
| rho=e^t | Optional positive multiplicative scale for that same continuation |
| Native lambda coordinates | The existing coordinates retain their definitions; no identification with t, n or rho is asserted |

In particular n tending to infinity is an increasing number
of cuts, not a substitution for every earlier lambda-infinity
limit. The operator U(t) specified here need not equal W^t.

For clarity, the intuition called "zero potential" is formalized
here as a zero compressed generator PAP. It is not the zero
vector: the zero vector is fixed by every linear evolution.
A zero generator inside a cut can conceal a nonzero transverse
velocity.

## CEM-T1. The native generator has zero cut motion and nonzero transverse motion

Put v(y)=y(1-y) and define
\[
\boxed{
A:=\frac{W-W^{-1}}2
=\begin{pmatrix}0&-2(1-y)\\2y&0\end{pmatrix}.
}
\tag{4}
\]
For k_h=((1-y)h,-yh), the exact identities are
\[
\boxed{
A^Tg+gA=0,\quad A^2=-4vI,\quad
AJf=-2k_f,\quad Ak_h=2vJh,\quad PAP=0.
}
\tag{5}
\]
Moreover ker A={0} on the polynomial carrier E.

**Proof.** The determinant of W is one, so its explicit inverse
changes the signs of its off-diagonal entries. This gives (4).
Direct multiplication gives (5). Equivalently use the polynomial
change of basis
\[
x=Jr+k_h,\qquad
r=yx_1+(1-y)x_2,\qquad h=x_1-x_2 .
\]
Its change-of-basis determinant is -1, so it is invertible
over R[y]. In the (r,h) coordinates,
\[
P=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
g_{\mathrm{cut}}=\begin{pmatrix}1&0\\0&v\end{pmatrix},\qquad
A=\begin{pmatrix}0&2v\\-2&0\end{pmatrix}.
\tag{6}
\]
These matrices prove all the displayed relations and show that
A maps the retained channel into the complementary channel.

If Ax=0, then A^2x=-4vx=0. Since v is a nonzero polynomial
and R[y] has no zero divisors, x=0. In particular AJf is nonzero
whenever f is nonzero. Nevertheless its cut projection vanishes.
QED.

Consequently exp(t PAP) is the identity on ran P, although no
nonzero polynomial vector is fixed by exp(tA) for every t.
The latter assertion follows by differentiating at t=0.
Fixedness at one isolated time is a different question.

## CEM-T2. Exponential transport preserves the moving cut and the source pairing

Define entire functions of t and y by their power series:
\[
\begin{aligned}
c_t(v)&=\sum_{m=0}^\infty
\frac{(-4v)^m t^{2m}}{(2m)!},\\
s_t(v)&=\sum_{m=0}^\infty
\frac{(-4v)^m t^{2m+1}}{(2m+1)!}.
\end{aligned}
\tag{7}
\]
No square-root branch is needed. Equation (5) gives
\[
\boxed{
U(t):=e^{tA}=c_t(v)I+s_t(v)A,\quad
U(-t)=U(t)^{-1},\quad U(t)^TgU(t)=g.
}
\tag{8}
\]
For a transported cut and embedding,
\[
P_t=U(t)PU(-t),\qquad J_t=U(t)J,
\]
one has
\[
\boxed{
P_t^2=P_t,\quad P_tJ_t=J_t,\quad
\mathcal G_a(J_t f,J_t h)=B_a(f,h),\quad
\frac{dP_t}{dt}=[A,P_t].
}
\tag{9}
\]

**Proof.** Reduce the even and odd powers of A using A^2=-4vI.
The exponential group law gives the inverse. Differentiating
U(t)^TgU(t) and using A^Tg+gA=0 shows that this matrix is
constant in t; its value at zero is g. Conjugation proves the
projection and range identities, and differentiating it proves
the commutator equation. Since J^TgJ=1, the pairing in (9)
equals L_a[fh].

The series converge uniformly for t in compact complex sets
and y in the closed unit disk. The source sum therefore extends
to these functions by (2). It can be differentiated in t locally
uniformly. These statements follow from the factorial denominators,
the bound |v(y)|<=2 on that disk, and sum |w_mu|<=2.
No positivity or Hilbert completion is required. QED.

Thus the preserved source reading concerns a jointly transported
frame and profile. The projection P itself is not a finite
matrix exponential: det P=0, whereas a finite matrix exponential
is invertible. Its moving conjugate in (9) is the precise
exponential-transport statement supplied here.

Under the optional convention rho=e^t, U(log rho) is a
multiplicative-scale group. This convention does not give limits
at rho=0 or infinity. For example, at the illustrative algebraic
fiber y=1/2,
\[
A=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
U(t)=\begin{pmatrix}\cos t&-\sin t\\\sin t&\cos t\end{pmatrix}.
\tag{10}
\]
The moving cut rotates, and does not converge as t tends to
either infinity. This is a counterexample to deriving endpoint
equality from the cut identities alone. The illustrative fiber
is not asserted to be an actual source point.

## CEM-T3. Repeated native cuts freeze the channel unconditionally

On the fixed cut define
\[
Z_n(t)=\bigl[P\,U(t/n)\,P\bigr]^n,\qquad n\ge1.
\]
Then
\[
\boxed{
P U(t)P=c_t(v)P,\quad
Z_n(t)=c_{t/n}(v)^nP,\quad
Z_n(t)\longrightarrow P .
}
\tag{11}
\]
The convergence is uniform for bounded t and |y|<=1.
It also holds after taking any fixed source pairing with
polynomial profiles.

**Proof.** Equation (5) gives PAP=0, so (8) gives the first
identity. The scalar c commutes with P, and P^2=P gives the
second. Uniformly on bounded t and the closed unit disk,
\[
c_{t/n}(v)=1-\frac{2vt^2}{n^2}+O(n^{-4}).
\tag{12}
\]
For sufficiently large n this lies uniformly in a fixed disk
around 1 on which the analytic logarithm is defined. Hence
\[
n\log c_{t/n}(v)=-\frac{2vt^2}{n}+O(n^{-3}),
\]
which tends uniformly to zero. Exponentiation proves (11).
The absolute source budget in (2) allows passage to any fixed
polynomial pairing. QED.

This is a bounded analytic realization of the classical
repeated-projection or Zeno mechanism. A physical probability
interpretation would require a positive carrier; the algebraic
limit above does not assume one.

At every n, Z_n(0)=P. At infinite cut count the entire evolution
on ran P becomes the identity. This freezing happens whether
or not the actual source form is positive. It therefore cannot
by itself supply the missing RH sign condition.

## CEM-T4. Rescaled curvature retains exactly the same source memory at every depth

For real polynomial f,h put
\[
R_{a,n}(t;f,h)
=\mathcal G_a(Jf,Z_n(t)Jh)
=\mathcal L_a[c_{t/n}(v)^n f h]
\]
and define
\[
\boxed{
\mathcal M_a(f,h)=\mathcal L_a[vfh]
=\frac14\,\mathcal G_a(AJf,AJh).
}
\tag{13}
\]
Then, for every integer n>=1,
\[
\boxed{
-\frac n4\,\partial_t^2 R_{a,n}(0;f,h)
=\mathcal M_a(f,h).
}
\tag{14}
\]
For each fixed real t!=0 there is also the exact asymptotic limit
\[
\boxed{
\lim_{n\to\infty}
\frac{n}{2t^2}\bigl(B_a(f,h)-R_{a,n}(t;f,h)\bigr)
=\mathcal M_a(f,h).
}
\tag{15}
\]
Moreover
\[
\boxed{
\mathcal M_a(f,h)
=\frac14\mathcal G_a(AJ_t f,AJ_t h),\qquad
\operatorname{rad}\mathcal M_a=\{0\}.
}
\tag{16}
\]

**Proof.** Since k_f^Tg k_h=vfh, AJf=-2k_f proves (13).
In (7), c_0=1, c'_0=0 and c''_0=-4v. Differentiating
c_(t/n)^n twice at zero gives n(-4v/n^2)=-4v/n.
The locally uniform source differentiation established in T2
proves (14).

The logarithm estimate in T3 gives, uniformly on the source disk,
\[
c_{t/n}(v)^n=1-\frac{2vt^2}{n}+O(n^{-2}).
\]
Subtract from one, multiply by n/(2t^2), and apply the bounded
source functional to prove (15).

Since A commutes with U(t), equation (8) preserves the pairing
of AJf and AJh under U(t). This proves the first part of (16).
For the radical, if M_a(f,h)=0 for every polynomial h, then
B_a(vf,h)=0 for every h. OQD-T2 gives vf=0. The polynomial
v is nonzero, so f=0. QED.

The unscaled limiting reading is R_(a,infinity)(t;f,h)=B_a(f,h),
independent of t. Its second derivative is zero. Equation (14)
shows that multiplying the prelimit second derivative by n
retains M_a exactly at every depth. No noncommutation of
unscaled differentiation and convergence is asserted: the
additional scaling is essential.

"Memory" here means this specified bilinear response coefficient.
It does not mean a subjective memory or an unmeasured physical
quantity. Even a zero self-reading M_a(f,f) would not imply that
all cross-readings M_a(f,h) vanish; only the latter defines the
radical in (16).

## CEM-T5. Positivity of the retained memory is exactly the all-order RH condition

Write c_j(a)=L_a[y^j], as in the earlier notes. The memory
matrix on the monomial bank of size N is
\[
\boxed{
D_N(a)=
\bigl(c_{i+j+1}(a)-c_{i+j+2}(a)\bigr)_{0\le i,j<N}.
}
\tag{17}
\]
For every fixed finite a>=a_0,
\[
\boxed{
\mathrm{RH}
\ \Longleftrightarrow\
\mathcal M_a(f,f)\ge0\quad\hbox{for every }f\in\mathbb R[y]
\ \Longleftrightarrow\
D_N(a)\succeq0\quad\hbox{for every }N.
}
\tag{18}
\]

**Proof.** Under RH the actual u_mu are in (0,1) and w_mu>0.
Hence v(u_mu)>0, and
\[
\mathcal M_a(f,f)=
\sum_\mu w_\mu v(u_\mu) f(u_\mu)^2\ge0.
\]
It is in fact strictly positive for any nonzero polynomial:
the support is infinite and a nonzero polynomial has finitely
many roots. Expanding a polynomial in monomials proves the
equivalence with all matrices (17).

For the converse retain the finite-detector argument of OQD-T5,
but with the exact modified weights
\[
\widetilde w_\mu=w_\mu v(u_\mu).
\tag{19}
\]
These weights are conjugation symmetric, nonzero, and have
sum |wtilde_mu|<=4. They have the identical discrete support.
They need not sum to one; the detector does not require that.

If RH fails, that support contains a nonreal point z.
Choose a real polynomial p annihilating every other source
point of modulus at least s=|z| except z and its conjugate.
There are only finitely many such points. The remaining support
has maximum modulus s_0<s, while p(z)!=0.

For observers y^m p(y) and y^(m+1)p(y), their memory Gram
matrix is
\[
\begin{pmatrix}
2\Re A_m&2\Re(A_m z)\\
2\Re(A_m z)&2\Re(A_m z^2)
\end{pmatrix}
+E_m,\qquad
A_m=\widetilde w_z z^{2m}p(z)^2,\qquad
(E_m)_{ij}=O(s_0^{2m}).
\tag{20}
\]
The main determinant is -4|A_m|^2(Im z)^2. Since
|A_m|=|\widetilde w_z p(z)^2|s^(2m) with a nonzero prefactor,
division of the full determinant by |A_m|^2 gives the limit
-4(Im z)^2<0. Thus some finite memory compression is
indefinite, contradicting positivity on every polynomial.
This proves the converse in (18). QED.

This argument identifies an equivalent test family, not an
actual negative witness or a proof of its universal positivity.
The vanishing of v at the accumulation point 0 does not remove
any actual source point, because 0<|u_mu|<1 and u_mu!=1.
It may change quantitative detector sizes; no new size bound
is claimed.

For any single fixed n, the matrix of rescaled second derivatives
in (14) on the monomial bank is exactly D_N(a). Increasing the
cut count changes none of these rescaled finite readings.
The missing positivity condition remains all N.

## CEM-T6. The exact sign balance identifies what the cut discards

The scalar identity following from (7) is
\[
c_t(v)^2+4v\,s_t(v)^2=1.
\]
Let P_perp=I-P. For real polynomial f,h,
\[
\begin{aligned}
\mathcal G_a(PU(t)Jf,PU(t)Jh)
 &=\mathcal L_a[c_t(v)^2fh],\\
\mathcal G_a(P_\perp U(t)Jf,P_\perp U(t)Jh)
 &=4\mathcal L_a[v\,s_t(v)^2fh].
\end{aligned}
\tag{21}
\]
Their sum is B_a(f,h). Their cross-pairings vanish.

**Proof.** Equation (8) gives U(t)Jf=c_t Jf-2s_t k_f.
OQD's cut decomposition makes J and k orthogonal for G_a.
Their pairings are respectively L_a[fh] and L_a[vfh].
This proves (21), including its sum and the vanishing cross terms.
QED.

The complementary term in (21) starts as
4t^2 M_a(f,h)+O(t^4). Thus the curvature memory is exactly the
leading transverse source response.

Before the condition in (18) is established, this is a bilinear
sign balance, not a sum of nonnegative probabilities or energies.
The all-depth memory positivity cannot be inferred by writing
M_a=G_a(AJ.,AJ.)/4, because the carried G_a is precisely the
possibly indefinite source form. Under RH it becomes a positive
orthogonal decomposition. No positive Hilbert interpretation of
the native source has been inserted into the converse argument.

## CEM-T7. A single cut can hide a channel that the first time response detects

For x=Jr+k_h, define the complete polynomial observation families
\[
\begin{aligned}
O_0(f;x)&=\mathcal G_a(Jf,x),\\
O_1(f;x)&=\left.\frac{d}{dt}
\mathcal G_a(Jf,U(t)x)\right|_{t=0}.
\end{aligned}
\]
Then
\[
\boxed{
O_0(f;x)=B_a(f,r),\qquad
O_1(f;x)=2\mathcal M_a(f,h).
}
\tag{22}
\]
Consequently, if both readings vanish for every polynomial f,
then x=0. The complete value and first-response families together
separate all polynomial two-channel states, without assuming RH.

**Proof.** OQD's orthogonal decomposition gives the first formula.
Equation (5) gives Ax=-2k_r+2vJh. Pairing this with Jf gives
the second. The zero radicals of B_a and M_a then force r=h=0.
QED.

At the vector level the same identity is Px=Jr and PAx=2vJh.
A single instantaneous cut discards h, while the derivative
contains it. On the polynomial carrier this map is injective
because v is not a zero divisor.

This is an exact observability statement using every polynomial
test, not a finite-data or stable reconstruction guarantee.
It does not claim that a black-hole exterior can implement these
readings. In particular v(u_mu) tends to zero along the source
support, so formal division by v is not uniformly bounded on
the support. Injectivity alone supplies no noise bound.

## CEM-T8. The cut algebra alone cannot select a universal propagation speed

For every real alpha>0, retain the same source, P, g and J,
but set
\[
A_\alpha=\alpha A,\qquad U_\alpha(t)=U(\alpha t).
\]
All the transport, cut-freezing and zero-radical statements persist.
For the corresponding repeated-cut source return,
\[
\boxed{
P A_\alpha P=0,\quad
A_\alpha^Tg+gA_\alpha=0,\quad
-\frac{n}{4\alpha^2}
\partial_t^2R^{(\alpha)}_{a,n}(0;f,h)=\mathcal M_a(f,h).
}
\tag{23}
\]

**Proof.** The first two identities are homogeneous in A.
The replacement t by alpha t in T2--T4 proves the remaining
statements and multiplies every second time derivative by
alpha^2. The cut and its source pairing have not changed. QED.

Therefore these cut relations allow an arbitrary rescaling of
the transport rate. They cannot by themselves determine a
physical universal speed c. A spatial metric, a physical clock
and a causal propagation law have not been specified here.
Fixing such additional structure could constrain a physical
model, but that would require further premises and derivation.

There is a limited mathematical information statement:
U(t) is invertible, whereas P identifies distinct states.
Thus an instantaneous restricted reading can lose distinctions
that the full evolution preserves. T7 describes their recovery
in this particular algebra. This is not an identification of
an event horizon with P, nor a model of an evaporating black hole.

In general relativity the light cones determine local causal
propagation, and physical information cannot locally outrun
light; see the Max Planck Institute's
[causal-structure explanation](https://www.einstein-online.info/en/explandict/causal/).
NASA's [black-hole anatomy](https://science.nasa.gov/universe/black-holes/anatomy/)
distinguishes the event horizon from the predicted singularity
and notes that the singularity may mark the limits of the
classical description. Neither establishes a literal zero
reservoir that stores information. These physical statements
do not enter the algebraic proofs above.

## What this development establishes

| Statement | Scope |
|---|---|
| PAP=0, but ker A={0} on polynomial states | Unconditional native algebra |
| Jointly transported cut and profile preserve B_a | Unconditional source identity |
| Repeated cuts converge to a frozen channel | Unconditional locally uniform limit |
| Rescaled second-order response is M_a at every cut count | Exact unconditional identity |
| No nonzero polynomial is invisible to all M_a readings | Unconditional zero radical |
| Complete value and first-response readings separate both channels | Unconditional polynomial observability |
| The same cut permits arbitrary time rescaling | No speed selection from the cut algebra alone |
| M_a is positive on every polynomial square | Equivalent to RH; unproved |
| Physical mass, consciousness or universal endpoint identification | No such conclusion is supplied |

The most useful invariant supplied by this construction is the
rescaled response form in (14), together with its covariant
version (16). It is not necessary to treat the apparent fixed
profile of the repeatedly cut evolution as an ambient fixed point.

This gives a precise way to test the intuition that cutting can
hide motion: keep the second-order transverse response rather
than reading only the frozen limit. In this native source the
response retains the full all-order sign problem. Any next RH
step must establish positivity of the actual matrices (17), or
an equivalent complete source condition; repeated cutting alone
does not do so.

The independently established TPW estimate on the fixed smooth
neutral family remains Q_xi[g]>=||g||^2/500 for L<=11/20.
CEM supplies no larger support window, uniform finite-bank
positive range, or endpoint equality. The source anchor a, the
support width L and the new repetition index n remain distinct.

## Attribution, dependencies and execution

- [OQD-1](RH_OBSERVER_QUOTIENT_DESCENT.md), especially T2, T4
  and T5: faithful actual source, native cut decomposition,
  and the finite detector at every retained anchor.
- [RFT-1](RH_REFLECTION_TRANSFER.md): native W and P, faithful
  folded transfer, and the separation of positivity from formal
  adjoint products.
- [Finite observer detector](RH_FINITE_OBSERVER_DETECTOR.md):
  the two-observer nonreal-source determinant used in (20).
- [TPW-1](RH_TWO_PRIME_WINDOW.md): the preceding independent
  two-prime support estimate and its unchanged scope.
- P. Facchi, V. Gorini, G. Marmo, S. Pascazio and E. C. G.
  Sudarshan, [Quantum Zeno Dynamics](https://arxiv.org/abs/quant-ph/0004040),
  Phys. Lett. A 275 (2000), 12: classical repeated-projection
  dynamics. The elementary analytic case used here is proved
  directly in T3.

The term memory is a declared name for (13); the general Zeno
effect is not being renamed as a new discovery. The source
specialization and all stated identities are supplied with
proofs. Repository verification covers text structure, internal
links, immutable-head content, intended files and draft state.
No scientific computation, numerical spectrum, source evaluation,
independent agent, replay or workflow dispatch is performed.
