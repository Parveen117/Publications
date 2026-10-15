# MZC-1 — Moving balance sets and the RH sign obligation

15 October 2026. Research continuation for Monty Dabas.

## 1. Source diagnosis and scope

The pinned RKF T23 calls the cut a primitive object, realized by projectors
P,Q and the involution P-Q. A zero reading, target faithfulness and positivity
are additional conditions on derived maps. Object and condition are different
roles, not mutually exclusive interpretations. SE26 already separates the
comparison diagonal from a zero scalar, spatial diagonal and critical line.
This packet applies the QG record-completeness lesson to that existing work;
it does not claim to discover covariance or Schur complements.

The mathematical adapter below is finite-dimensional over R or C unless
explicitly stated otherwise. Transport, its generator and its parameter are
declared data. The parameter is not assumed to be physical time. A cut alone
does not select a dynamics. Analytic infinite-dimensional RH assertions remain
subject to the source's domain, completion, arithmetic and tail hypotheses.

## 2. Balance-set covariance

Let C0:X->Y be linear and U(t):X->X invertible. With a fixed output frame set

    x(t)=U(t)x0,   C(t)=C0 U(t)^(-1).

Then C(t)x(t)=C0x0 and ker C(t)=U(t)ker C0.

Proof: cancel U^(-1)U in the reading. For the kernel, C(t)y=0 iff
U^(-1)y belongs to ker C0. Both implications use invertibility.
If the output frame also moves by an invertible V(t), the correct law is
C(t)=V(t)C0U(t)^(-1); the reading becomes V(t)C0x0 and the same kernel law holds.

For differentiable U with U'=G(t)U, one obtains C'=-CG and
x'=Gx, hence (Cx)'=0. A constant G gives U=exp(tG). A variable generator
requires its fundamental solution; an ordinary exponential of the integral
is not generally valid without appropriate commutation.

A fixed observer instead reads C0U(t)x0, which need not vanish. The numerical
zero vector always satisfies U(t)0=0. What can move is a nontrivial balance
set and its nonzero members. For a reference point r0 and an affine frame
F(t)x=U(t)x+b(t), set r(t)=U(t)r0+b(t). The comparison
C0U(t)^(-1)(y-r(t)) has zero set F(t)(r0+ker C0).
This is a coordinate reference, not a new zero of a fixed analytic function.

## 3. Operators, forms and the retained metric

An endomorphism transforms by similarity D(t)=U D0 U^(-1); its kernel moves
by U. If D0 is self-adjoint and U is not unitary, D(t) need not be self-adjoint
in the original Euclidean metric. Retain

    g(t)=U^(-*) U^(-1),
    H(t)=U^(-*) H0 U^(-1).

Then x(t)^*H(t)x(t)=x0^*H0x0. Consequently H(t)>=0 iff H0>=0.
Every negative witness transports to a negative witness. The dimensions of
positive and negative subspaces are unchanged because U gives a bijection
between them. If H0=D0=D0*, g(t)D(t)=H(t), so D(t) is self-adjoint in g(t).
This identifies the missing metric record when an operator is transported by
nonunitary similarity. Dropping it changes the meaning of positivity.

For H0(t) with genuine intrinsic evolution, define the displayed form
Htilde=U^(-*)H0(t)U^(-1) and Gamma=U'U^(-1). Differentiation yields

    Htilde' + Gamma* Htilde + Htilde Gamma = U^(-*) H0' U^(-1).

Thus the covariant derivative separates frame motion from source evolution.
Pure frame motion has zero right-hand side; it cannot remove a negative sector.

## 4. Why a moving zero record is insufficient

Take Hplus=diag(0,1), Hminus=diag(0,-1). They have exactly the same kernel,
and so do all their congruent transports under the same U. Yet one is
nonnegative and the other has a negative direction. Even the complete
trajectory of the zero set does not determine the sign of a form.

SE29 supplies a second obstruction directly in the RH adapter:

    Sb=N^(1/2)(I-T)N^(1/2),
    Delta=(I-T)(I+T)^(-1),   K=T(I+T)^(-2),
    Delta^2+4K=I.

Under its stated bounded positive-source hypotheses, Sb>=0 iff Delta>=0
iff T<=I. Setting Delta=0 identically would impose T=I, a different and
stronger requirement than nonnegativity. For gains 9/16 and 16/9,
K=144/625 in both cases, but Delta=+7/25 and -7/25 respectively.
Curvature alone loses the sign. Selecting a positive square root without
excluding the actual negative spectral sector would assume the desired result.

SE26 also retains a common-mode warning: its Gamma contrast annihilates q
while the pole reading on q equals 45607/433125. This packet cites that source
calculation; it does not independently recompute its integral. A zero contrast
cannot stand in for the full source, common mode, mixed readings and tails.

## 5. Innovation signs survive lawful changes of the acquired frame

This reuses SE49 Theorem 4, equation (49.12), which already proves the
real upper-triangular LDL-pivot scaling law. The block proof below permits
a complex coefficient and any invertible change within the old subspace;
its use here is the connection to moving-cut records, not a new RH criterion.

Let a Hermitian acquired form have an invertible positive old block A:

    M = [[A,b],[b*,d]],       sigma=d-b* A^(-1)b.

Change the old basis and the new representative without changing the old
subspace: W=[[U,v],[0,a]], with U invertible and scalar a nonzero.
The transformed form is M'=W*MW. Its old block is A'=U*AU and
cross block b'=U*(Av+ba). Its new diagonal is
d'=v*Av+v*ba+conj(a)b*v+|a|^2 d. Using
(U*AU)^(-1)=U^(-1)A^(-1)U^(-*) gives

    b'* A'^(-1)b' = v*Av+v*ba+conj(a)b*v+|a|^2 b*A^(-1)b,
    sigma' = |a|^2 sigma.

This proves sign invariance, including zero. Equivalently, the exact
completion of the square is
(e+A^(-1)b z)*A(e+A^(-1)b z)+sigma|z|^2.
The result uses a filtration-preserving basis change. Arbitrary mixing of
old and new subspaces changes the conditional question and is not this theorem.
The native symbolic witness checks the real scalar cancellation; the written
argument proves the general complex block statement.

Apply this separately to SE49's L and J forms. Their actual next-node margins
are alpha=-d-ell^T L^(-1)ell and beta=r+xd-j^T J^(-1)j. Lawful acquired-frame
changes multiply each margin by a strictly positive scale. They cannot turn
an actual failed extension into a successful one. Choosing a different physical
node changes the acquired form, so this is not a prediction of unknown readings.
It can improve numerical enclosure: SE49's original interval boxes containing
zero were inconclusive, not negative witnesses, and its lawful basis change
made positivity certifiable. Sign invariance does not forbid that improvement.

## 6. Precisely what remains for RH

The pinned SE49 source reports an outward certificate on V24 (48 positive
innovations) and leaves universal positivity open. We have not rerun that RH
certificate here. Its proposed all-node criterion requires actual alpha_n>0
and beta_n>0 for every admitted extension, together with its analytic bridge.
The present packet verifies coordinate covariance and exact control examples;
it neither proves that universal condition nor supplies new zeta evaluations.

The next substantive RH target is a uniform, source-derived lower bound for
the two actual Schur innovations, or an equivalent exclusion of the negative
source sector. A moving frame must carry the common mode, signed form/metric,
cross readings and source/tail identification. Exponential transport preserves
an established bound; it does not create the missing bound.

The proof system replays specified finite algebraic identities and controls.
It is not a general proof-assistant certification of the analytic RH bridge,
and this packet is not a solution of RH or quantum gravity.
