# Native cut transport, energy descent and an internal response clock

Monty Dabas · CT-1–CT-7 · 30 September 2026

## What is supplied, and what is derived

This continues [CF-1–CF-8](../coherence-first-thermodynamics/THEOREM.md).
Its native carrier is the canonical cut-field EMK algebra, with
R²=-1, K²=1, KR=-RK, R†=-R and K†=K. Put L=RK.
Its recognition target retains the complete distinction D=sK+vL, both
aperture channels, the ordered transport and any lifted winding record.
The shared RKF operator engine is imported unchanged at commit
`3cc5a33b05c16d59c90994ddda69dedc0d392424`.

We keep the explicitly admitted CF thermodynamic chart x=(s,v), r=|x|,
and energy W=κr^ν+εs²/2, with κ>0, 1<ν<2, ε≥0.
The new process hypothesis is **uniform attenuation of the two reference-
subtracted conjugate responses**:

\[
p=\nabla W=(T-T_0,-P+P_0),\qquad p\longmapsto a p,
\quad 0<a\le1. \tag{1}
\]

Here a is a positive central cut scalar. Equation (1) is a selected lawful
protocol, not a theorem that coherence uniquely forces that protocol. It
acts in response coordinates, not by identifying native numeral multiplication
with transport. We derive its unique state path, energy balance, native
transport and internal additive clock. Physical seconds, spatial distance,
inertial mass and a universal gravitational law are not supplied.

## CT-1. A unitary cut rotation cannot dissolve the contrast amplitude

D²=r²1 and N=D/r is an involution on r>0. For every dagger-unitary U,

\[
(UDU^\dagger)^2=UD^2U^\dagger=r^2 1. \tag{2}
\]

Therefore changing orientation by such conjugation cannot change r. This is
a native multiplication result, without a Hilbert-space premise. It says
nothing against amplitude loss in a recognized channel after compression:
that is a different operation and requires the missing-channel ledger.

For a differentiable distinction path, let ρ=r'/r and ω=φ', where
N=cosφ K+sinφ L in the admitted coefficient chart. Since
[R,K]=2L and [R,L]=-2K,

\[
D'=\rho D+[\omega R/2,D]
   =G D+D G^\dagger,\qquad
G=\tfrac12\rho 1+\tfrac12\omega R. \tag{3}
\]

The coefficients ρ,ω are uniquely determined within this scalar-plus-R
ansatz by the specified nonzero path. The scalar part changes amplitude;
the dagger-skew part transports its cut. Equation (3) does not assert
uniqueness among all possible native generators.

## CT-2. Exact finite transport and the record that endpoint cuts forget

For two nonzero states on a path with a continuous angle lift, define

\[
U_{21}=\cos(\Delta\phi/2)1+\sin(\Delta\phi/2)R,
\quad B_{21}=\sqrt{r_2/r_1}\,U_{21}. \tag{4}
\]

Sine and cosine here denote the completed factorial flow of R, subsequently
expressed in the coefficient chart; R is not replaced by a scalar i.
Native multiplication gives

\[
U_{21}^\dagger U_{21}=1,\quad
D_2=B_{21}D_1B_{21}^\dagger,\quad
B_{21}^\dagger B_{21}=(r_2/r_1)1,
\quad P_2=U_{21}P_1U_{21}^\dagger. \tag{5}
\]

With consistent lifted angles B32 B21=B31. The operator is invertible for
positive endpoint radii. A path approaching the centre has B†B→0; its
limiting map loses invertibility. This is dissolution of this distinction
record, not a proof that the underlying whole is the zero state.

One full turn of D returns its aperture but gives U=-1; two full turns give
U=1. Thus even U retains only winding parity. Full winding is a separate
integer/path record whenever the recognition target asks for it. In the
exact witness two successive U=R steps return K with total U=-1. Resetting
the endpoint lift to +1 erases this difference. This is the stated algebraic
lift; physical spin statistics and a spacetime spin connection are not implied.

For a change of cut with angle Δφ, P1 P2 P1=cos²(Δφ/2)P1 and
P1 Q2 P1=sin²(Δφ/2)P1. These are algebraic compression coefficients, not
Born probabilities. With an intervening arrow A,

\[
P_2 A\psi=P_2 A P_1\psi+P_2 A Q_1\psi. \tag{6}
\]

For the exactly matched arrow U21 the cross term vanishes, because
P2 U21=U21 P1. For a mismatched arrow it need not vanish. The verifier tests
both: aligned transport and the original PRQRP=-P returning-memory witness.

## CT-3. Existence, uniqueness and a semigroup from response attenuation

The map F(x)=∇W(x) is a continuous bijection of the entire coefficient plane,
smooth with a smooth inverse away from the centre. To prove this, first note
that |x|^ν is strictly convex for ν>1; its addition to εs²/2 remains strictly
convex. For each finite p, W(x)-p·x is coercive since κr^ν-|p|r→∞. It
therefore attains a unique minimum; its stationarity equation is F(x)=p.
This proves surjectivity and uniqueness. On r>0, the positive definite CF
Hessian proves local smooth invertibility. Inverse continuity at zero follows
from p·x=κνr^ν+εs²≥κνr^ν, hence
|p|≥κνr^(ν-1). The minimum for p=0 is x=0.
These analytic statements are on the already admitted finite thermodynamic
chart; they are not a new primitive construction of that chart.

Consequently (1) defines an exact state arrow

\[
\Phi_a(x)=F^{-1}(aF(x)),\qquad
\Phi_b\circ\Phi_a=\Phi_{ab},\qquad\Phi_1=\mathrm{id}. \tag{7}
\]

For x≠0 and a>0 its state is never the centre. As a↓0 it tends to the
centre. Keeping a>0 makes each arrow invertible on the plane; the inverse
uses the amplification factor 1/a. Restricting allowed forward arrows to
0<a≤1 declares an orientation, not a fundamental irreversibility theorem.

These are global statements for the mathematical energy family. Physical
thermodynamic interpretation is restricted to the CF neighbourhood with
positive S,V,T,P. Choosing a sufficiently small invariant disk is possible:
CT-5 proves that r decreases, and all conjugate deviations tend to zero.

## CT-4. A native additive clock and a derived differential law

Use F00-G's positive radial logarithm and set

\[
\sigma(a)=-\operatorname{Log}_\Sigma(a),\qquad
\sigma(ab)=\sigma(a)+\sigma(b). \tag{8}
\]

This is an internal Morphic clock on the selected attenuation monoid. Any
continuous additive monotone clock on this monoid is cσ with c>0: compose
it with a=ExpΣ(-u), use additivity first for nonnegative rational u and then
continuity. Thus normalization remains free, and no SI second is derived.
For a known nonzero initial response, a can be recovered from any nonzero
component ratio p_i/p_i(0); the ratio is the same in both channels.

In this native clock p(σ)=ExpΣ(-σ)p(0). Differentiating F(x)=p gives

\[
\boxed{\frac{dx}{d\sigma}=-H^{-1}\nabla W,\qquad
\frac{dW}{d\sigma}=-p^T H^{-1}p<0\quad(x\ne0).} \tag{9}
\]

These equations follow from (1), rather than declaring a spacetime action.
The solution exists uniquely for every finite σ≥0, by (7), and approaches
the centre only as σ→∞. The Hessian's singularity there is not a finite-σ
breakdown along this protocol.

For fixed x, (9) is also the unique minimizer over tangent increments ξ of
J_x(ξ)=p·ξ+ξ^T H ξ/2: completing the square gives
J_x(ξ)=J_x(-H^-1p)+(ξ+H^-1p)^T H(ξ+H^-1p)/2.
This is a derived local variational characterization of the chosen protocol,
not a uniquely derived fundamental action. Multiplying W and H by the same
positive constant leaves the state law invariant while rescaling energy.

## CT-5. Inward motion, angular drift, and the native generator

Let e=(c,d)=(s/r,v/r), δ=ν-2, h=κνr^(ν-2), and
Δ=detH=h²(ν-1)+εh(1+δd²)>0. Direct inversion of the actual CF Hessian gives

\[
\boxed{\rho=\frac{1}{r}\frac{dr}{d\sigma}
=-\frac{h(h+\epsilon)}{\Delta}<0,\qquad
\omega=\frac{d\phi}{d\sigma}
=\frac{h\epsilon(\nu-2)cd}{\Delta}.} \tag{10}
\]

For example, the component derivatives are
s'=-hs[h+ε+εδd²]/Δ and v'=-hv[h+ε-εδc²]/Δ.
Taking s s'+v v' and s v'-v s' proves (10). Substituting these coefficients
in (3) provides the native generator for this constitutive trajectory.
Both inward movement and cut rotation are now consequences of a single
supplied response arrow, rather than unrelated prescribed paths.

For sv≠0, the conserved response ratio is

\[
\frac{p_s}{p_v}=\frac{(h+\epsilon)s}{hv}=\mathrm{constant}. \tag{11}
\]

Angular drift is nonzero off the axes when ε>0, and vanishes in the isotropic
control. The contrast angle φ is not the CF-5 principal-response angle θ.
For the latter set b=(Hss-Hvv,2Hsv). Along (9), b'=∂s b s'+∂v b v';
where |b|²>0, θ'=(b1 b2'-b2 b1')/(2|b|²). At isotropic Hessian points the
principal seam is undefined; (9) and the contrast cut can remain regular.
These calculations therefore do not substitute CF-5's fixed-ray derivative
for the derivative along a moving trajectory.

## CT-6. Exact energy-density scaling and a convergent clock ledger

When ε=0, (7) has the closed form

\[
r(a)=a^{1/(\nu-1)}r_0,\quad
\frac{W(a)}{W_0}=a^{\nu/(\nu-1)},\quad
\frac{\det H(a)}{\det H_0}=a^{2(\nu-2)/(\nu-1)}. \tag{12}
\]

At ν=3/2, κ=2/3, r0=1, a=1/2, the exact values are

| Observable | Before | After |
|---|---:|---:|
| Radius r | 1 | 1/4 |
| Formation energy W | 2/3 | 1/12 |
| Hessian determinant | 1/2 | 2 |
| Radial native transport factor sqrt(r/r0) | 1 | 1/2 |

A nonradial exact endpoint also exists: keep ν=3/2 and κ=2/3, but take
ε=7/5. The attenuation a=3/16 takes x0=(3/5,4/5) exactly to
x1=(1/20,3/80). Their radii are 1 and 1/16, and their directions are
(3/5,4/5) and (4/5,3/5). Direct substitution gives p1=(3/16)p0;
CT-3 establishes that this is the unique endpoint. Formation energy falls
from 689/750 to 73/6000. Thus finite inward motion and angular drift are
verified together without fitting an inverse solver.

These are model predictions for the supplied attenuation protocols, not measured
constants. The internal interval is LogΣ(2). It can be certified without
floating-point logarithms using F00-G's Cayley series:

\[
\operatorname{Log}_\Sigma2=2\sum_{j\ge0}\frac{(1/3)^{2j+1}}{2j+1}.
\]

After terms j=0,...,N, the positive tail is bounded by
2(1/3)^(2N+3)/[(2N+3)(1-1/9)]. The executable controls retain this bound
and the matching body/tail refinement ledger. They do not call a finite
polynomial the exact logarithm. This is the same correction-transfer grammar
as NC-4, now applied to the native internal clock.

For a finite chain, define released formation energy E_j=W_j-W_(j+1)>0.
The sum telescopes to W_initial-W_final exactly. A ledger recording this
amount restores an energy balance; identifying it with an actual reservoir's
heat or radiation requires a physical interaction model. No environment
has been experimentally identified by this bookkeeping identity.

## CT-7. What has closed and what is still a physical selection problem

For the supplied energy and response attenuation law, the mathematical bridge
now closes: unique state arrows, global finite-clock existence, an additive
internal clock, strict energy descent, native amplitude/rotation transport,
and exact memory obligations are derived together.

The law itself is not uniquely selected by coherence. For example, a direct
contrast attenuation x→ax is another lawful semigroup on the same chart.
In the isotropic case its response scales by a^(ν-1), while with anisotropy
the two response contributions scale differently. One cannot silently infer
(1) merely from the existence of positive native scalars.

For an external reading t and differentiable σ(t),

\[
\dot x=\dot\sigma X(x),\qquad
\ddot x=\dot\sigma^2(DX)X+\ddot\sigma X,
\quad X=-H^{-1}p. \tag{13}
\]

Changing that clock changes apparent acceleration without changing the
ordered intrinsic trajectory. A flat clock cannot resolve the transition.
No physical space adapter, force law, Einstein equation, c, G, or mass-energy
calibration follows just by renaming (13). The next discriminating experiment
is to identify an actual process with equal fractional attenuation of both
conjugate deviations, retain its records and measure its external clock.
The protocol predicts (11), (12), and the direction-dependent drift (10).
A different attenuation pattern rejects this protocol for that process;
it does not by itself refute every possible recognition law.

General written proofs, exact finite controls, source pins, proof-assistant
verification, independent review and experimental evidence remain distinct.
