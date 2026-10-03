# Critical native return channels and a third-probe prediction

Monty Dabas · CR-1–CR-5 · 30 September 2026

This continues AS-1–AS-5 and NI-1–NI-4 using the canonical R2 native return
family. It derives an inverse-response pole and its exact source/observer
selection rule, controls its coefficient with retained aperture uncertainty,
and predicts a third probe from the previously identified parameters.
**No numerical fine-structure constant is derived.**

## Object and interpretation contract

Let L=KR, L²=1, L†=L in the canonical EMK algebra. The completed boundary
return F_b(z)=1+x_b(z)L obeys

\[
bz x_b(z)^2+x_b(z)-(b+1)z=0,\quad b>0,\ z>0,\ x_b(z)>0.
\tag{1}
\]

The depth-paired periodic cells are ((b+1)z,bz). Existence and boundary-tail
independence are those of the unchanged canonical R2 theorem. The coefficient
x_b increases strictly, with x_b(1)=1. Put P±=(1±L)/2.

F_b is a boundary return operator (an inverse corner of the original depth
system). We now examine its **own algebraic inverse**, where it exists. This
is not the same observable as F_b and is not automatically a propagator.
A singular inverse of a bounded return is not by itself a diverging measured
signal, a massless particle or a physical instability. Such statements need
a source experiment or an action that actually uses this operator.

## CR-1. Exact channel decomposition and simple pole

Native multiplication gives complementary dagger-symmetric idempotents P±,
with P+P−=0 and P++P−=1. Therefore

\[
F_b=(1+x_b)P_+ +(1-x_b)P_-,\qquad
G_b=F_b^{-1}=\frac{P_+}{1+x_b}+\frac{P_-}{1-x_b},\quad z\ne1.
\tag{2}
\]

For positive z, x_b>0, so the only possible zero denominator is x_b=1;
substitution into (1) shows this occurs precisely at z=1. At that point
F_b(1)=2P+, its complementary channel vanishes, and a two-sided inverse
does not exist. The inverse on the retained P+ corner alone is P+/2;
it cannot be promoted to an inverse on the full algebra.

Set t=1-z and approach the cut from 0<z<1. Then 0<x_b<1. Rearrangement
of (1), retaining scalar commutation only, gives the exact identity

\[
\boxed{(1-x_b)[1+bz(1+x_b)]=1-z=t.} \tag{3}
\]

Thus, without exchanging a depth limit and a derivative,

\[
\boxed{\lim_{t\downarrow0}tG_b(1-t)=(1+2b)P_-.} \tag{4}
\]

The plus-channel contribution vanishes after multiplication by t; the
minus-channel coefficient follows directly from (3) and x_b→1.
This proves a first-order pole in the declared parameter t. In terms of
z-1 the residue has the opposite sign. AS-1's susceptibility s=1/(1+2b)
is the reciprocal of this normalized pole strength, not an independent
constant. Exactly closed cuts with different b still have different residues.

The next finite term, using AS-1's first derivative in (3), is

\[
G_b(1-t)=\frac{1+2b}{t}P_-
+\frac12P_+-\frac{b(3+4b)}{1+2b}P_-+O(t).
\tag{5}
\]

All statements concern this native algebraic inverse. No spacetime dimension,
physical clock, electromagnetic field or source normalization is inferred.

## CR-2. The pole is visible only to the matching source and observer

On an admitted native module, take a source ψ and a linear observation ℓ.
Equation (4) implies

\[
\lim_{t\downarrow0}t\,\ell(G_b\psi)
=(1+2b)\ell(P_-\psi). \tag{6}
\]

The pole contributes to this observation exactly when ℓ(P−ψ) is nonzero.
For a source entirely in the P+ channel it is absent, even though G_b has
an operator pole. The native R arrow exchanges the channels since RL=-LR:
P−R=RP+. It can therefore move a retained source into the critical channel.
An observer restricted to P+ would still erase that contribution. Source,
transition, target and observer are all necessary parts of a coupling claim.

These are module and normal-form statements. A source-source positive
quadratic residue can be defined only after a compatible pairing and source
normalization are admitted. It must also retain AS-3's hidden-source dressing.
An algebraic residue cannot be relabelled an electron charge by dropping these
contracts or the complementary-channel memory.

## CR-3. Stable certification close to the singular channel

Define the normalized minus-channel response

\[
D_b(t)=\frac{t}{1-x_b(1-t)}=1+b(1-t)[1+x_b(1-t)]. \tag{7}
\]

This exact form avoids division by the small difference 1-x_b. If the native
finite-aperture solver gives x_b∈[l,u], then

\[
D_b(t)\in[1+bz(1+l),\ 1+bz(1+u)],\quad z=1-t,
\tag{8}
\]

with interval width bz(u-l). The direct inverse coefficient would require
u<1 and has width (u-l)/[(1-u)(1-l)]. If a finite interval crosses one,
it cannot certify a finite inverse bound through that formula. This refusal
is a precision issue, not proof that the completed operator is singular away
from z=1. Formula (8) remains a lawful enclosure of the normalized quantity.

There is also an explicit limiting error. Since (3) gives 0<1-x_b≤t,

\[
0\le(1+2b)-D_b(t)
=b[2t+(1-t)(1-x_b)]\le3bt. \tag{9}
\]

If m=(l+u)/2 and epsilon=(u-l)/2, the finite certificate is

\[
\left|(1+2b)-[1+bz(1+m)]\right|
\le3bt+bz\,\epsilon. \tag{10}
\]

It explicitly separates distance from the critical point and retained-depth
uncertainty. Finite cutoffs are never called exact infinite returns. The
verifier checks these bounds using the unchanged native aperture solver,
and replays both inverse products in the native algebra at rational witnesses.

## CR-4. A third probe without a newly fitted parameter

NI's fixed calibration fixture has normalized readings x(1+δ)=6/5 and
x(1-δ)=2/5. Its unique solution is b=5/7, δ=3/4; those readings and their
certificate were committed before this packet. Choose the specified third
relative probe z3=1+δ/2=11/8. The model now predicts, with no extra fitted
parameter, the positive root of

\[
\boxed{55x_3^2+56x_3-132=0,
\qquad x_3=\frac{-28+2\sqrt{2011}}{55}.} \tag{11}
\]

The code produces an outward rational enclosure through the original R2
solver and verifies that the polynomial changes sign across it. This is a
held-out *model calculation*: no experimental third reading has been taken.
It is not an experimentally validated fundamental constant. The relative
probe factor 1/2 is the declared test setting; it does not add a fit parameter.

A future actual experiment could calibrate b and δ from the two readings,
freeze them, and compare this third response. Its probe must act as the same
common multiplier of native cell products; an arbitrary control dial is not
that adapter. A discrepant third reading would reject this periodic-response
model for that process, even if the first two readings were fitted exactly.

In the same fixture, the algebraic inverse-pole strength is 1+2b=17/7.
This is predicted from the identified model parameter. It is not alpha.

## CR-5. The remaining photon and normalization bridge

If one *additionally* supplies a spectral adapter t=eta k², eta>0, equation
(4) gives a 1/k² inverse coefficient (1+2b)/eta. This is an algebraic
substitution, not a derivation that t is physical momentum squared. A quadratic
cost with kernel lambda F_b would multiply the inverse coefficient by 1/lambda.
Under a separately established physical photon, charged-source and action
matching, the corresponding scalar candidate would therefore take the form

\[
\alpha_{\rm candidate}
=\frac{|q_-|^2(1+2b)}{4\pi\lambda\eta}. \tag{12}
\]

Here q− denotes the matched and source-dressed coupling to the critical
channel. Every identification in this sentence is additional: transverse
photon degrees of freedom, physical momentum, quantum action scale and the
electron source are not supplied by (1)–(11). There is no basis for setting
q−=lambda=eta=1 by fiat and calling the resulting number measured alpha.
The dimensional types must be fixed in a common dimensionless-action
normalization before even making that comparison.

The new result is the exact critical-channel residue, its observer selection
rule, a finite error budget and a third-probe model prediction. The native
source law still has to select b rather than infer it from readings, and the
physical matching must fix q−, lambda and eta together. **Alpha is NOT DERIVED.**
Written proofs, finite controls, physical experiment and independent review
remain separate evidence. The canonical engine and earlier packets are unchanged.
