# QG-R1 — Native relative curvature without a clock

Development 1.0. Monty research programme, 14 October 2026.

The next construction is a **dimensionless curvature contrast per native
comparison area**, computed from source and reference transports on the same
four-vertex cell. It is exactly invariant under a different native unitary
frame at every vertex. Its numerator and denominator use the same operator
data; no physical time, length calibration or spacetime metric is required.
The extension R7 gives a bounded return-defect numerator and a precise
concentration condition on shrinking native cells. Under stated comparison
admission, divergence obstructs a regular connection completion.

This continues [Phase Memory and Angular Current](../native-phase-observation-certification/README.md).
That stage supplied phase currents and a frozen-generator loop coefficient.
Here the loop uses all four possibly different edges, and the reference is
inside the construction. The general connection formula is consumed from
[EMK-C1](../emk-ugd-algebra/CONNECTION_CALCULUS.md), not proposed again.

## Native inputs and observation contract

Use the unchanged native cut field and its faithful two-mode representation:

\[
R^2=J^2=-1,\qquad RJ=-JR,\qquad J=\iota K,\qquad L:=RJ.
\]

Consequently `L^2=-1`, and the three words are anti-dagger, with
`[R,J]=2L`, `[J,L]=2R`, `[L,R]=2J`. These are internal operator directions,
not three assumed spatial axes. A native unit transport has the form

\[
Q=s1+xR+yJ+zL,\qquad s^2+x^2+y^2+z^2=1
\]

with real native coefficients. The cyclic native scalar part is
`sc=Tr/2` in this representation. On anti-dagger words define

\[
\langle X,Y\rangle=-\operatorname{sc}(XY).
\]

For `X=xR+yJ+zL`, its square-length is `x^2+y^2+z^2`.
This pairing is derived within the selected native sector; identifying it
with a physical metric would require an additional representation.

The input is a specified or independently reconstructed native link on each
edge, with a matched reference link. The preceding scalar channels do not
automatically reveal arbitrary additional central phases or all multi-axis
links. The real unit-transport sector, anchors and reconstruction/calibration
contract must therefore be supplied. Arbitrary detector gain on the four
scalar readings below is not an admitted symmetry.

Label the square vertices `0,1,2,3` in order, with base vertex `0`.
Use positive edges `x:0->1`, `y+:1->2`, `x+:3->2`, `y:0->3`.
Let their source links be `U_e` and reference links be `V_e`. Define

\[
H_U=U_y^\dagger U_{x+}^\dagger U_{y+}U_x,
\quad H_V=V_y^\dagger V_{x+}^\dagger V_{y+}V_x,
\]
\[
P_x=V_x^\dagger U_x,\quad P_y=V_y^\dagger U_y,
\quad P_\square=H_V^\dagger H_U.
\]

All three are endomorphisms at the same base vertex. Order is essential.

## R1. Exact cancellation of a moving native frame

Allow an independent native unit frame `G_v` at every vertex. Both the
source and reference transform by

\[
U_e'=G_{t(e)}^\dagger U_eG_{s(e)},\qquad
V_e'=G_{t(e)}^\dagger V_eG_{s(e)}.
\]

Then each relative object transforms only at the base:

\[
P_x'=G_0^\dagger P_xG_0,\quad P_y'=G_0^\dagger P_yG_0,
\quad P_\square'=G_0^\dagger P_\square G_0.
\]

**Proof.** Substitute the edge laws. Every adjacent target/source frame
pair cancels by `G_v G_v^dagger=1`. For an edge,
`(V_e')^dagger U_e'=G_s^dagger V_e^dagger U_e G_s`.
The four-edge product telescopes to base conjugation, and the same
cancellation applies to `H_V^dagger H_U`. These are native word identities.

Changing the physical reference links is a change of comparison data.
R1 concerns a simultaneous frame change, not invariance under replacing
the reference experiment by a different one.

## R2. The finite curvature/area quotient

For a unit word away from its antipode define the rational native readout

\[
\mathcal D(Q)=2(1-Q)(1+Q)^{-1}.
\]

It is anti-dagger and has no scalar part. Direct native multiplication gives

\[
\mathcal D(Q)=-\frac{2(xR+yJ+zL)}{1+s},\qquad
\langle\mathcal D(Q),\mathcal D(Q)\rangle
=\frac{4(1-s)}{1+s}.
\]

The domain is `s != -1`, equivalently `Q != -1` in this unit sector.
Set `D_x=D(P_x)`, `D_y=D(P_y)`, `K_square=D(P_square)` and

\[
g_{ij}=\langle D_i,D_j\rangle,\quad
\mathcal A^2:=g_{xx}g_{yy}-g_{xy}^2>0.
\]

The required two native comparison directions must be independent. The
Gram determinant is a sum of three squared two-by-two coefficient minors.
Define

\[
\boxed{\displaystyle
\mathfrak C_\square=
\frac{\langle K_\square,K_\square\rangle}{4\mathcal A^2}.}
\tag{R2}
\]

This is nonnegative, exactly invariant under R1, and has no clock or external
length unit. The factor four fixes the calibration used in R5 below.

**Proof.** Conjugation commutes with products, sums and inverses, hence
with `D`. The cyclic scalar part preserves the three pairings and the
numerator. The native coefficient identity proves positivity of the Gram
determinant on independent directions. Division gives the claim.

The exact same value results if labels `h,k != 0` are assigned afterwards,
using `D_x/h`, `D_y/k` and `K_square/(h*k)`: numerator and denominator both
acquire the factor `1/(h^2*k^2)`. Those labels can be discarded. Changing
the actual links or refining the cell generally changes the finite value.

### Four scalar readings suffice

Let

\[
a=\operatorname{sc}(P_x),\quad b=\operatorname{sc}(P_y),\quad
c=\operatorname{sc}(P_x^\dagger P_y),\quad
d=\operatorname{sc}(P_\square),
\]
\[
\Delta=(1-a^2)(1-b^2)-(c-ab)^2.
\]

On the admitted domain `-1<a,b<1`, `-1<d<=1`, and `Delta>0`,

\[
\boxed{\displaystyle
\mathfrak C_\square=
\frac{(1-d)(1+a)^2(1+b)^2}{16(1+d)\Delta}.}
\tag{R2s}
\]

**Proof.** Write `P_x=a+v`, `P_y=b+w` with anti-dagger vector parts.
Native multiplication gives `c=ab+<v,w>` and the unit identities give
`<v,v>=1-a^2`, `<w,w>=1-b^2`. Thus
`A^2=16*Delta/((1+a)^2*(1+b)^2)`. Substitute the displayed norm formula
for `K_square`. This proves the four-reading reconstruction without roots.

The cross-reading `c` is necessary for this unrestricted interface: rotate
the second comparison direction from `J` to `3R/5+4J/5`, holding its length,
the first edge, and the loop conjugacy class fixed. Then `a,b,d` agree but
the area-square is multiplied by `16/25`; the quotient is multiplied by
`25/16`. The executable constructs the complete two sets of links realizing
this comparison. This is a specific three-reading insufficiency result,
not a lower bound for every possible instrument.

## R3. Variable fields: the full mixed coefficient

For the optional smooth two-control representation, use the existing
connection convention `nabla=d+A`; a forward link has first jet
`1-h*A_i`. The control labels need not be time or spatial coordinates.
With real central `eps^2=eta^2=0`, the variable-edge word is

\[
\begin{aligned}
H_A={}&(1+\eta A_y)
(1+\varepsilon A_x+\varepsilon\eta\,\partial_y A_x)\\
&\quad(1-\eta A_y-\varepsilon\eta\,\partial_x A_y)
(1-\varepsilon A_x)\\
={}&1-\varepsilon\eta F_A,\qquad
F_A=\partial_x A_y-\partial_y A_x+[A_x,A_y].
\end{aligned}
\tag{R3}
\]

**Proof.** Expand the four factors with native noncommutative multiplication.
The pure first-order terms cancel, the derivative terms retain the signs
shown, and the mixed product is `-[A_x,A_y]`. The compiler certifies this
as one general noncommuting-word identity. Freezing derivatives recovers
the prior AC4 coefficient. Its two plus-oriented elementary generators
give the same frozen quadratic sign; the convention above fixes the signs
of the new derivative terms.

Let `B` be the reference connection, `C=A-B`, and `F_rel=F_A-F_B`. Then

\[
\boxed{F_{\rm rel}=D_B C+C\wedge C},\qquad
\boxed{P_\square=1-\varepsilon\eta F_{\rm rel},\quad
K_\square=\varepsilon\eta F_{\rm rel}}.
\tag{R3r}
\]

In components the first formula means
`partial_x C_y-partial_y C_x+[B_x,C_y]+[C_x,B_y]+[C_x,C_y]`.
It follows by subtracting the existing two curvature formulas. The second
and third follow by multiplication and the finite nilpotent inverse.
Relative curvature zero does not imply that either connection is flat.

## R4. The missing third-channel derivative is measurable in principle

In a flat-reference gauge `B=0`, expand the general anti-dagger connection
in the already generated three-word sector:

\[
C=R\alpha+J\beta+L\gamma.
\]

No new primitive generator is introduced: `L=RJ`. The complete curvature is

\[
\boxed{\begin{aligned}
F_R&=d\alpha+2\beta\wedge\gamma,\\
F_J&=d\beta+2\gamma\wedge\alpha,\\
F_L&=d\gamma+2\alpha\wedge\beta.
\end{aligned}}
\tag{R4}
\]

**Proof.** Substitute in `dC+C wedge C`; use the three native brackets.
All coefficients are real central forms. The executable certifies every
matrix entry of the component formula.

An oriented comparator `L` gives `F_L=<L,F>`. Therefore, in this fixed
reference frame,

\[
\boxed{\partial_x\gamma_y-\partial_y\gamma_x
=\langle L,F_{xy}\rangle
-2(\alpha_x\beta_y-\alpha_y\beta_x).}
\tag{R4r}
\]

This recovers the missing **curl of the third connection coefficient** if
the oriented curvature and the two stated coefficient directions are
available. It does not uniquely recover `gamma` itself: adding an exact
form preserves this equation. A cotransformed comparator makes the
projection frame-consistent. The four scalar readings R2s alone provide
only a magnitude quotient and cannot supply this signed projection.

At a base point with `C_x=R`, `C_y=J`, two continuations illustrate why this
matters. A fixed two-axis continuation has `F=2L`. A continuation with
`partial_y C_x=2L`, `partial_x C_y=0` has `F=0` at the same point.
The instantaneous coefficient pair and its Gram matrix agree; the missing
derivative changes the curvature. Exact endpoint-generated links
`U_e=T_t^dagger*T_s` realize flat loops while their two base comparison
directions can remain independent. Thus the frozen bracket is insufficient.

## R5. Continuum comparison and the two-current calibration

Assume C2 source/reference coefficient fields on an admitted two-control
chart. Generate each forward edge by the native Cayley rule
`U_i(x;h)=cayley(A_i(x),-h/2)` and likewise for `B`; use its inverse for
the reversed edge. For `h,k` tending to zero, on `det g^0>0`,

\[
\frac{D_x}{h}\to C_x,\quad \frac{D_y}{k}\to C_y,\quad
\frac{K_\square}{hk}\to F_{\rm rel},
\]
\[
\boxed{\mathfrak C_\square\to\mathfrak C_0
=\frac{\langle F_{\rm rel},F_{\rm rel}\rangle}
{4\det g^0},\qquad g^0_{ij}=\langle C_i,C_j\rangle.}
\tag{R5}
\]

**Proof.** The Cayley first derivative gives the edge limit. The loop
function equals the identity on both axes `h=0` and `k=0`. Its mixed
derivative at the origin is R3. Two native scalar integrations therefore
express `(H-I)/(hk)` as the average of that continuous mixed derivative
over the rectangle; it tends to `-F`. The same argument and the regular
inverse at `P_square=I` give the relative-loop limit. Rescale the Gram
entries as in R2 and divide by its nonzero limiting determinant.

The C2 limit argument is a written native-completion proof; the executable
certifies its exact mixed jet and finite rational sequences, not arbitrary
analytic quantifiers. A nonsingular change of the two smooth control
coordinates sends `g^0` to `M^T*g^0*M` and `F_12` to `det(M)*F_12`.
The quotient is unchanged. This smooth-coordinate statement is distinct
from changing the finite lattice itself.

Pointwise Cayley sampling here is a chosen discretization. An exact finite
frame change acts on the existing links by R1. Regenerating links only from
the transformed connection value at the start of an edge generally agrees
to first order, not at every finite step. For example `G(x)=cayley(R,x+x^2)`
has `G^dagger dG=2R dx` at zero, but its exact edge is
`cayley(R,-h-h^2)`, whereas that pointwise rule gives `cayley(R,-h)`.
The executable rejects their identification.

For the fixed-axis, flat-reference subcase `C=R*alpha+J*beta`, put
`j=alpha_x*beta_y-alpha_y*beta_x != 0`,
`u=(d alpha)_xy`, `v=(d beta)_xy`. Native pairing gives

\[
\boxed{\mathfrak C_0=1+\frac{u^2+v^2}{4j^2}.}
\tag{R5c}
\]

For actual smooth one-axis angular currents `alpha=d theta_Q`,
`beta=d theta_G`, their local closure makes `u=v=0` and `C_0=1`.
Nonzero `u,v` instead describe a nonintegrable transport-coefficient
protocol or a failure of that angular-current contract. They cannot be
assigned to an everywhere smooth scalar phase by relabelling.

This value is a calibration of the stated transport family, not a universal
curvature constant or a value of Newton's constant. The full three-channel
case permits derivative cancellation and `C_0=0` with nonzero comparison
area. Recognition curvature and metric curvature retain the separation
already proved in EMK-G1/G2.

## R6. Exact finite benchmarks

With identity reference and constant links
`U_x=cayley(R,-q)`, `U_y=cayley(J,-p)`, the inherited loop formula gives
`chi=16*q^2*p^2/((1+q^2)^2*(1+p^2)^2)`. Since
`D_x=2qR`, `D_y=2pJ`, native substitution yields

\[
\boxed{\mathfrak C_\square=
\frac{1}{(1+q^2)^2(1+p^2)^2-16q^2p^2}.}
\tag{R6}
\]

For `q=1/2,p=1/3` the exact result is **324/481**. It is unchanged by
independent frames at the four vertices. As `q,p` tend to zero it tends
to `1`; its finite value is not asserted to equal the limit. At `q=p=1`,
the loop is `-I`: the group element is regular but this Cayley readout is
outside its chart. The verifier rejects it and never calls it a physical
singularity. Dependent comparison directions are likewise outside R2.

A variable-field example `A_x=(1+y)R`, `A_y=(1+x)J`, `B=0` has at the
origin `F=-R+J+2L`, hence `C_0=3/2`. Exact finite Cayley rectangles test
that sign and target. A second example from `theta_Q=x+xy`, `theta_G=y`
has cancelling mixed scalar derivatives and target `C_0=1`.

## R7. Return-defect concentration and a regularity obstruction

The concentration proposal needs a function that remains finite when only
the loop Cayley chart fails. Define the **return-defect information density**
operationally by

\[
\chi_{\rm ret}(P)=\frac{1-\operatorname{sc}P}{2}
=\frac14\operatorname{sc}\big[(1-P)^\dagger(1-P)\big],
\qquad 0\leq\chi_{\rm ret}\leq1,
\]
\[
\boxed{\mathcal I_{\rm ret}(\square)
=\frac{\chi_{\rm ret}(P_\square)}{\mathcal A^2}.}
\tag{R7a}
\]

The denominator is the **squared native comparison area** from R2. This is
a dimensionless measure of loop mismatch relative to comparison area, not
an assumed Shannon entropy or information amount. A physical information
interpretation requires a specified encoding and measurement contract.
No time variable enters this definition.

**Finite theorem.** On the admitted edge charts and positive area, this
density is invariant under all R1 vertex-frame changes, including when
`P_square=-1`. On the smaller loop Cayley chart it obeys

\[
\boxed{\mathcal I_{\rm ret}=(1-\chi_{\rm ret})\mathfrak C_\square,\qquad
\mathfrak C_\square=\frac{\chi_{\rm ret}}
{\mathcal A^2(1-\chi_{\rm ret})}.}
\tag{R7b}
\]

The same four scalar readings recover it:

\[
\boxed{\mathcal I_{\rm ret}
=\frac{(1-d)(1+a)^2(1+b)^2}{32\Delta},\qquad
\Delta=(1-a^2)(1-b^2)-(c-ab)^2>0.}
\tag{R7c}
\]

Here `-1<a,b<1` as in R2, but **`-1<=d<=1`** is admitted. Edge antipodes
and zero comparison area are still excluded. The directional shape record
is

\[
\sigma_\square=\frac{\mathcal A^2}{g_{xx}g_{yy}}
=\frac{\Delta}{(1-a^2)(1-b^2)}\in(0,1].
\tag{R7d}
\]

**Proof.** Native unitarity gives
`(1-P)^dagger(1-P)=2(1-sc P)1`. Scalar cyclicity and R1 give frame
invariance. R2 gives `||K||^2=4 chi_ret/(1-chi_ret)`; divide by
`4 A^2` to obtain R7b. Substitute the already proved scalar Gram
determinant to obtain R7c and R7d. At `P=-1`, the new numerator is
exactly `1`, so the density is `1/A^2`, finite at any fixed positive area.
Thus the old loop-chart pole by itself is not a concentration singularity.
For the R6 benchmark, `chi_ret=144/625`, `A^2=4/9`, and
`I_ret=324/625`, whereas the Cayley quotient was `324/481`.

**Concentration condition.** Let `square_n` be a declared family of cells
shrinking towards a point on a path, with positive native areas tending
to zero and with admitted edge charts. Then

\[
\boxed{\limsup_{n\to\infty}\mathcal I_{\rm ret}(\square_n)=+\infty
\quad\Longleftrightarrow\quad
\chi_{\rm ret}(P_{\square_n})\ne O(\mathcal A_n^2).}
\tag{R7e}
\]

This equivalence is the definition of asymptotic boundedness after dividing
by the positive area. It says more than "a large defect": even a defect
tending to zero can fail to decay as fast as the squared native area.
Cells attached to a path, their reference, and their comparison directions
are required data. A scalar assigned to a bare one-dimensional path does
not supply these loops or their transverse comparison area. A result along
one declared cell family does not locate every possible singular point.

**Regularity obstruction.** Suppose cells are generated by R5's Cayley
sampling rule from C2 source/reference connections on one neighborhood of
a point `p`, with bases `p_n -> p`, control increments `h_n,k_n -> 0`, and
`det g^0(p)>0`. All cells and their corners stay in that neighborhood.
Then

\[
\boxed{\mathcal I_{\rm ret}(\square_n)\longrightarrow
\frac{\langle F_{\rm rel}(p),F_{\rm rel}(p)\rangle}
{4\det g^0(p)}<\infty.}
\tag{R7f}
\]

Consequently R7e excludes a C2 realization of this kind. If the measured
normalized Gram matrices

\[
G_n=\begin{pmatrix}
g_{xx}/h_n^2 & g_{xy}/(h_nk_n)\\
g_{xy}/(h_nk_n) & g_{yy}/k_n^2
\end{pmatrix}
\]

remain bounded and uniformly positive definite, a hypothetical C2
realization cannot avoid the obstruction by allowing its relative
comparison directions to collapse at `p`. The increments are optional
control labels for this completion theorem, not a primitive physical clock.
The finite observable R7a does not use them. The shape record `sigma` is
a useful direction-collapse diagnostic, but by itself does not establish
all calibration or smooth-realization assumptions.

**Proof.** R5 applies uniformly on a sufficiently small compact
neighborhood: the mixed loop derivative is continuous jointly in base
point and increments. Hence `K/(hk) -> F_rel(p)`,
`A^2/(h^2 k^2) -> det g^0(p)>0`, and `P -> 1`, so `chi_ret -> 0`.
Use R7b and R5 to obtain R7f. Under a hypothetical C2 realization,
`G_n -> g^0(p)`; uniform positivity makes its determinant positive.
Contraposition proves the obstruction. This is a written completion
argument, not a conclusion inferred from finitely many test cells.

**Exact families separating the possibilities.** Put `t_n=2^(-n)`,
`n>=1`, and use a flat reference. With the native convention
`cayley(X,u)=(1-uX)^(-1)(1+uX)`, choose

\[
P_x=\operatorname{cayley}(R,-t/2),\quad
P_y=\operatorname{cayley}(J,-t/2).
\]

Then `D_x=tR`, `D_y=tJ`, `A^2=t^4`, `sigma=1`, and the normalized Gram
matrix for `h=k=t` is the identity. Any target loop `P` is realized by
the full cell `U_x=P_x`, `U_y=P_y`, `U_y+=1`,
`U_x+=P_x P^dagger P_y^dagger`. This multiplication gives `H_U=P`
exactly, so the examples are complete cell protocols.

| Family | Target loop | Return defect | Density | Result |
| --- | --- | --- | --- | --- |
| Bounded calibration | `cayley(L,-t^2)` | `t^4/(1+t^4)` | `1/(1+t^4)` | Tends to 1 |
| Concentration | `cayley(L,-t)` | `t^2/(1+t^2)` | `1/[t^2(1+t^2)]` | Diverges with `sigma=1` |

For the concentrating family, all four links tend to the identity and the
loop stays away from the Cayley pole. Nevertheless

\[
\mathcal I_{\rm ret}\geq\frac{4}{5t^2}
=\frac45\,4^n\longrightarrow\infty,\qquad
\frac{K_\square}{t^2}=\frac{2L}{t}.
\]

The bound holds for every `n>=1` because `1+t^2<=5/4`. The first six
densities run from `16/5` to `16777216/4097`. These cells therefore
cannot all arise by R5 sampling a common C2 pair across the limiting point
with the stated calibration. They are an exact mathematical construction,
not observed gravitational data. A common smooth realization of the
bounded calibration family is not asserted merely from its bounded density.

A separate control exposes false concentration from collapsing directions.
Replace the `J` direction in `P_y` by the unit word

\[
B_t=\frac{1-t^2}{1+t^2}R+\frac{2t}{1+t^2}J
=a_tR+b_tJ,
\]

and retain `P=cayley(L,-t^2)`. Then

\[
\sigma=b_t^2\to0,\quad
\mathcal I_{\rm ret}=\frac1{(1+t^4)b_t^2}\to\infty,
\quad \frac{K_\square}{t^2}=2L.
\]

The curvature estimate is bounded here; the relative comparison area
collapses. This is why density divergence alone cannot establish curvature
blow-up or a physical singularity. At a fixed independent unit-area cell
with `P=-1`, conversely, `I_ret=1` while the old Cayley readout is undefined.
The verifier exercises both controls and independent vertex frames.

**Relationship to the earlier information paper.** The audited
[information-invariance v2](../information-invariance/paper/information_invariance_v2.tex)
defines a different object, the signed loop functional
`I(gamma)=integral_gamma omega`. Its `chi` is a microscopic statistical
weight ratio, not `chi_ret`. No equality with R7 is assumed or proved.
Its [ledger](../information-invariance/LEDGER.md), including the result
that invariance alone selects no dynamics and the withdrawn general
emergent-Einstein derivation, remains unchanged. R7 contributes a specific
native observable and a conditional continuation obstruction. It neither
revives that derivation nor proves that quantum gravity is uniquely a
condition on one scalar information function. Bounded `I_ret` also has
no proved converse guaranteeing regularity.

## Lineage, status and next physical step

| Existing source | Reused here | New obligation handled here |
| --- | --- | --- |
| Native lattice / phase-word packets | Unit words, loop defect, native engine and phase convention | Exact source/reference lattice readout and four-scalar quotient |
| PM/AC certification | Anchors, current identities, fixed-word sign, limits of finite endpoint data | Four varying edges and the third-channel derivative |
| EMK-C1 C2/C3/C5 | Full connection curvature, covariance, pure-gauge flatness, mixed channels | Match those formulas to the relative finite-cell observable |
| Native thermodynamic curvature, NT-1 | Cyclic native scalar part and finite pairing | Use it on the anti-dagger relative-transport sector |
| YM54 T1/T2 | Gram/commutator-area method | A different relative-link quotient; no reuse of YM heat/gap conclusions |
| EMK-G1/G2 | Independence of metric and recognition curvature | Preserve that distinction while constructing an operator comparison area |
| Information-invariance v2 and ledger | Distinction between loop invariance, microscopic reversibility and unselected dynamics | A separately defined return-defect density and its conditional regularity obstruction |

Lattice gauge invariance and relational observables have substantial prior
literature: Wilson, *Confinement of quarks*, Phys. Rev. D 10, 2445 (1974),
[DOI](https://doi.org/10.1103/PhysRevD.10.2445), and Rovelli, *Partial
observables*, Phys. Rev. D 65, 124013 (2002),
[arXiv:gr-qc/0110035](https://arxiv.org/abs/gr-qc/0110035).
These are contextual references, not additional proof engines or claims
that those papers prove this quotient. No broad priority claim is made.

R1–R7 are scoped mathematical results with written proofs. The executable
coverage is native rewrite witnesses plus exact finite tests and rejection
controls, detailed in the README and certificate. No external Hilbert or
spacetime structure is used as a primitive.

For the quantum-gravity programme this supplies a clock-free, frame-invariant
curvature observable, its variable-field interpretation, a precise
missing-channel reconstruction equation and a concentration obstruction
for regular native connection realizations. A gravitational identification
still requires a physical source/reference realization and a relation of
this connection to physical geometry and dynamics. A quantum state/measure
and a continuum gravitational law are not selected by this quotient.
The next focused task is to supply a native transport/source law that fixes
the third-channel response and can be compared with an independently known
gravitational observable, while retaining the reference and area records.
