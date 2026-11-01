# QG-R2 — Many centres, one observation base

Native loop fusion, retained concentration and the missing orientation record.
Development 1.0, 14 October 2026.

The next step after [QG-R1](../native-relative-curvature/THEOREM.md) is to
compare several loop experiments at a single observation base. This gives
an exact composition rule, but a combined loop can hide nonzero individual
defects. A labelled record retains them. Three loops already expose a new
measurement obligation: pairwise scalars miss an oriented triple reading.

## Motivation and inherited scope

The suggested "Kennedy theorem" is interpreted as the Aronhold–Kennedy
three-centre theorem. In planar rigid-body kinematics its three relative
instantaneous centres are collinear. See West Virginia University's
[MAE 342, Analysis of Mechanism Velocity, slide 10](https://community.wvu.edu/~bpbettig/MAE342/Lecture_9_velocity_analysis.pdf).
That theorem concerns relative motion; it does not transport density or
prove a singularity. It motivates the common-reference question here and
is not a premise of the native derivation. No claim of a mathematical
generalization of that kinematic theorem is made.

Use QG-R1's unchanged native sector
`R^2=J^2=L^2=-1`, `RJ=L`, `JL=R`, `LR=J`, with reversed products negative.
Write a unit word `Q=a1+v_R R+v_J J+v_L L`, so
`a^2+v_R^2+v_J^2+v_L^2=1`. Native scalar part is `sc=Tr/2`.
The return defect is the already defined observable

\[
\chi(Q)=\frac{1-\operatorname{sc}Q}{2}
=\frac18 E(1-Q)\in[0,1],
\qquad E(X)=\sum_{r,s}|X_{rs}|_{\rm cut}^2.
\]

This operational defect is not identified with an entropy or the earlier
information-invariance paper's signed loop functional. The positive native
coefficient norm is derived from the cut field; no primitive Hilbert space,
physical clock or spacetime metric is added.

## F1. Reversible base transport and reference-correct fusion

At site `i` retain source loop `H_i^U`, reference loop `H_i^V`, and a
declared native unit connector `T_i:o->i`. Their common-base versions are

\[
u_i=T_i^\dagger H_i^U T_i,\quad
v_i=T_i^\dagger H_i^V T_i,\quad p_i=v_i^\dagger u_i.
\]

For independent unit frames `G_o,G_i`, transform local loops by
`H_i'=G_i^dagger H_i G_i` and connectors by
`T_i'=G_i^dagger T_i G_o`. Then `u_i,v_i,p_i` all transform by the same
base conjugation `X'=G_o^dagger X G_o`. Scalar parts, defect values and
transported comparison Gram areas are unchanged. Keeping `T_i` makes
rebasing reversible: `H_i=T_i(T_i^dagger H_i T_i)T_i^dagger`.

To concatenate source and reference experiments in order, define

\[
U=u_m\cdots u_1,\quad V=v_m\cdots v_1,\quad
W_{i-1}=v_{i-1}\cdots v_1,\quad W_0=1,
\]
\[
Q_i=W_{i-1}^\dagger p_iW_{i-1}.
\]

Then the **exact relative fusion law** is

\[
\boxed{P=V^\dagger U=Q_m\cdots Q_1.} \tag{F1}
\]

All `Q_i` transform by common base conjugation and `chi(Q_i)=chi(p_i)`.
The reference prefixes cannot generally be omitted.

**Proof.** Base-frame cancellation follows by direct substitution and
`G_i G_i^dagger=1`. For two loops,
`V^dagger U=v_1^dagger p_2 v_1 p_1`. More generally, writing
`u_i=v_i p_i` in the full product inserts `W_{i-1}^dagger p_i W_{i-1}`
successively; adjacent prefixes cancel. This is an induction on `m`.
The native compiler certifies the general rebasing word and the two- and
three-loop prefix identities. The finite-m induction is written here.

Different connectors or a different reference are different protocol data.
They can change cross comparisons and the fused word even while preserving
every individual defect. Moving the observation base in this manner does
not move matter, reduce a physical region, or identify distinct sites.
Also, `P` is a concatenated loop protocol; identifying it with the boundary
of a geometrically tiled region requires separate shared-edge and connector
compatibility. No such identification is assumed.

## F2. Three-loop fusion requires an orientation record

For three dressed words `Q_i=a_i1+V_i` write their real native vector
coefficients as `v_i`. Define six pairwise readings and one oriented reading:

\[
a_i=\operatorname{sc}Q_i,\qquad
c_{ij}=\operatorname{sc}(Q_i^\dagger Q_j),\quad
g_{ij}=c_{ij}-a_i a_j=\langle V_i,V_j\rangle\ (i\ne j),
\]
\[
g_{ii}=1-a_i^2,\qquad
\tau=\operatorname{sc}(V_3V_2V_1)=\det[v_1\ v_2\ v_3].
\]

The determinant uses the inherited orientation `(R,J,L)`; the scalar
expression for `tau` makes its common-frame invariance explicit. Then

\[
\boxed{\operatorname{sc}(Q_3Q_2Q_1)
=a_1a_2a_3-a_3g_{12}-a_2g_{13}-a_1g_{23}+\tau,}
\tag{F2a}
\]
\[
\boxed{\tau^2=\det(g_{ij}).} \tag{F2b}
\]

Thus three scalar parts and three pairwise comparisons determine the
magnitude of the missing orientation but generally leave its sign open.
An oriented triple reading supplies it without choosing a square-root
branch. If the Gram determinant is zero, `tau=0` and that ambiguity closes.
These are reconstruction statements for actual admitted unit words,
not an assertion that arbitrary seven scalars are realizable.

**Proof.** Native multiplication gives
`(a,v)(b,w)=(ab-v dot w, a w+b v+v cross w)` in this representation.
Multiply in the stated order and take the scalar coefficient. The term
from three vectors is `+det[v_1 v_2 v_3]` for `Q_3 Q_2 Q_1`.
Expansion of the three-by-three Gram determinant gives F2b.
Both identities and the orientation reversal are native polynomial
certificates over arbitrary real central coefficients.

**Exact blindness witness.** The triples `(R,J,L)` and `(R,J,-L)` have
identical `a_i=0` and `c_ij=0`. The first has `tau=1` and product
`LJR=1`; the second has `tau=-1` and product `-LJR=-1`.
Their fused defects are respectively `0` and `1`, although all six
pairwise scalars agree. Each individual defect is `1/2` in both cases.
Reversing the order also changes the result. Therefore the sign record
is necessary on this fibre; pairwise density data cannot determine fusion.

This is a specific orientation repair for QG-R1's native loop data. The
more general principle that a quotient can lose orientation and require
a retained bit already appears in RST-1 T2/F4; it is reused, not renamed
as a new general principle.

## F3. Exact positive bound on coherent fusion

For any fixed finite number `m>=1` of unit words, set

\[
S=\sum_{i=1}^m\chi(Q_i),\qquad P=Q_m\cdots Q_1.
\]

Then

\[
\boxed{0\leq\chi(P)\leq\min\{1,mS\}.} \tag{F3}
\]

**Proof from native squares.** Define
`D_i=Q_m ... Q_(i+1)(1-Q_i)`. Direct telescoping gives
`1-P=sum D_i`. A unit left factor preserves `E`, so `E(D_i)=8 chi(Q_i)`.
For each native real and turn coefficient, expand

\[
m\sum_i E(D_i)-E\left(\sum_iD_i\right)
=\sum_{i<j}E(D_i-D_j)\geq0.
\]

Dividing by eight proves the bound; `chi(P)<=1` is inherited. This uses
only a sum of native coefficient squares. The compiler certifies the
three-term scalar identity and noncommuting telescope; their arbitrary
finite-m extension follows by counting each square and cross term.

There is **no strictly positive universal lower bound** proportional to
`S`: choose `Q_3=(Q_2Q_1)^dagger` with nontrivial `Q_1,Q_2`.
Then `P=1`, while `S` can be positive. The upper factor `m` cannot be
uniformly reduced: for three aligned words `q=cayley(R,-t)`,

\[
\chi(q^3)=\frac{t^2(3-t^2)^2}{(1+t^2)^3},\quad
\frac{\chi(q^3)}{3\sum_{i=1}^3\chi(q)}
=\frac{(3-t^2)^2}{9(1+t^2)^2}\longrightarrow1.
\]

For general fixed `m`, the native first-order product gives the same
sharpness limit using `q^m=1-2mtR+O(t^2)` and the native norm definition.
This limiting sharpness argument is written, not inferred from samples.

## F4. Retained records merge without cancellation

Keep the labelled tuple `(Q_1,...,Q_m)`; a direct block record is
`Q_record=diag(Q_1,...,Q_m)`. Then

\[
\frac18 E(1_{2m}-Q_{\rm record})=S,\qquad
S=0\quad\Longleftrightarrow\quad Q_i=1\text{ for every }i.
\tag{F4a}
\]

This is faithful for detecting any nonzero return defect. The scalar `S`
does not reconstruct the tuple; preserving full operator information still
requires the labelled record, reference data and connectors. Its record
capacity grows with the number of retained loops.

For each local cell retain its positive squared native area `b_i=A_i^2`
from QG-R1 and its density `I_i=chi(Q_i)/b_i`. Merge the **numerator and
budget**, not just density values:

\[
(s,b)\oplus(s',b')=(s+s',b+b'),\quad
B=\sum_i b_i,\quad
\overline{\mathcal I}=\frac{S}{B}
=\sum_i\frac{b_i}{B}\mathcal I_i.
\tag{F4b}
\]

The merge is associative and independent of grouping. Positivity gives

\[
\boxed{\min_i\mathcal I_i\leq\overline{\mathcal I}
\leq\max_i\mathcal I_i.} \tag{F4c}
\]

Thus a uniform bound on all child densities survives repeated merging.
Bringing their records to one observation base alone cannot violate that
bound. This is a statement about a finite family or a uniformly bounded
family, not a bound inferred from separately finite entries.

The budget `B` is the sum of squared comparison areas. It is an additive
record budget; it is **not asserted to be the squared geometric area of
a union of cells**. The distinction is essential under refinement.

**Proof.** Block norm additivity and positive native squares prove F4a.
Addition of numerator/budget pairs proves associativity. Divide by `B>0`
to obtain positive weights summing to one; subtract each proposed bound
and use nonnegative terms to prove F4c. No external probability model or
arithmetic engine is introduced.

## F5. Concentration at a chosen receiving cell

Supply an actual receiving comparison pair `X_*,Y_*` with positive
native Gram determinant `b_*`, using QG-R1's observation contract. Define

\[
\mathcal I_{\rm ret,*}=\frac{S}{b_*},\qquad
\mathcal I_{\rm fused,*}=\frac{\chi(P)}{b_*},\qquad
\kappa=\frac{B}{b_*}.
\]

These satisfy the exact accounting identities

\[
\boxed{\mathcal I_{\rm ret,*}=\kappa\overline{\mathcal I},\qquad
0\leq\mathcal I_{\rm fused,*}\leq m\mathcal I_{\rm ret,*}.}
\tag{F5a}
\]

When `S>0`, record the measured composition factor
`c=chi(P)/S` in `[0,m]`. Then

\[
\boxed{\mathcal I_{\rm fused,*}=c\,\mathcal I_{\rm ret,*}.} \tag{F5b}
\]

For fixed `m`, if `c>=c_0>0` along a refinement family, retained and
fused density are bounded together and their unbounded limsups coincide.
Without that condition, coherent cancellation can hide retained divergence.
Neither `c` nor `kappa` is fixed by invariance alone; both belong in the
observation record. Choosing a smaller denominator on paper does not
establish physical focusing: `b_*` must come from the receiving experiment.

**Conditional regularity result.** Suppose a fixed finite family of child
cells admits the C2 source/reference realizations and nondegenerate
calibrations of QG-R1 R7 near its limiting sites. If `kappa` stays bounded,
then both receiving densities stay bounded. Therefore divergence of
`I_ret,*`, under that fixed-count and bounded-budget-ratio admission,
excludes the simultaneous regular realizations.

**Proof.** QG-R1 R7 bounds every `I_i`. A finite number gives a uniform
bound. F4 bounds their weighted mean; bounded `kappa` and F3/F5a finish
the proof. With finitely many indices, an unbounded weighted mean implies
an unbounded subsequence in at least one child density. It need not make
all children singular or localize all their limiting points at the receiver.
This remains a native continuation obstruction, not a gravitational
singularity theorem.

For loops of a common small area increment `epsilon`, after admitted
transport and reference dressing, write
`Q_i=1-epsilon F_i+o(epsilon)`. Then

\[
P=1-\epsilon\sum_iF_i+o(\epsilon),\quad
\frac{S}{\epsilon^2}\to\frac14\sum_i\langle F_i,F_i\rangle,
\quad
\frac{\chi(P)}{\epsilon^2}\to
\frac14\left\langle\sum_iF_i,\sum_iF_i\right\rangle.
\tag{F5c}
\]

The first-order product is a native word identity; the norm limits follow
from `chi=E(1-Q)/8` and `E(F)=2<F,F>` for anti-dagger words.
Retained squares and the square of the combined curvature differ by
cross terms. Keeping this distinction is the continuum form of the record
requirement, with a fixed finite number of loops and the stated remainders.

## F6. Three concentrating cells with an exactly flat combined loop

Take `t=2^(-n)`, `n>=1`, identity reference and identity connectors. Set

\[
Q_1=\operatorname{cayley}(R,-t),\quad
Q_2=\operatorname{cayley}(J,-t),\quad
Q_3=(Q_2Q_1)^\dagger.
\]

All three are nontrivial and tend to the identity; their combined loop
is exactly `P=Q_3Q_2Q_1=1` at every refinement. Choose independent
comparisons `D_x=tR,D_y=tJ` for each child and the receiver, giving
`b_i=b_*=t^4`, shape `1` and `kappa=3`.
Every target `Q_i` is realized by a full four-edge cell using QG-R1's
construction `U_x=P_x`, `U_y=P_y`, `U_y+=1`,
`U_x+=P_x Q_i^dagger P_y^dagger`, with
`P_x=cayley(R,-t/2)`, `P_y=cayley(J,-t/2)`.

Exact native multiplication gives

\[
S=\frac{2t^2(2+t^2)}{(1+t^2)^2},\qquad
\boxed{\mathcal I_{\rm ret,*}
=\frac{2(2+t^2)}{t^2(1+t^2)^2}\longrightarrow\infty,\qquad
\mathcal I_{\rm fused,*}=0.} \tag{F6}
\]

For every `n>=1`, `I_ret,* >= 64/(25t^2)`, since `t^2<=1/4`.
The budget-normalized density is `I_ret,*/3` and also diverges.
The first six steps are exact rational certificate cases; the displayed
formula and bound prove the all-refinements claim. This is an explicit
constructed family, not physical gravitational data or a common smooth
connection asserted to exist across the limiting point.

As a separate control, use `m_n=2^n` copies of the regular loop
`cayley(L,-t^2)` with the same individual budget `t^4`. Each child and the
budget-normalized mean have density `1/(1+t^4)<=1`. Yet the total divided
by one receiving budget is `m_n/(1+t^4)`, which diverges because
`kappa=m_n` grows. This records increasing multiplicity, not a proved
local curvature singularity. It cannot pass F5's admission unchanged.

## Source audit and certification boundary

| Source | Reused content | Addition here |
| --- | --- | --- |
| QG-R1 R1/R2/R7 | Native frames, positive area, bounded return defect and C2 obstruction | Many-site base transport, reference-prefix fusion and retained/fused distinction |
| RST-1 T2/F4 and EMK-TOP-1 T3/T6 | Lost orientation, global cancellation and retained sector records | An explicit three-loop orientation fibre and positive native concentration ledger |
| GE3 returning memory | Keeping versus deleting hidden records changes observation | Context only; its heat generator and dynamics are not imported |
| Kennedy three-centre theorem | Verified kinematic motivation | No density or singularity claim is attributed to it |

F1–F6 have written proofs with explicit inputs. The unchanged native
compiler certifies the free-word identities, three-vector polynomial
identities and three-term positive-square identity, with serialized replay.
Exact native matrix cases independently exercise reconstruction, arbitrary
site/base frames, reference prefixes, cancellation, budget accounting and
invalid inputs. General finite-m induction, positivity implications and
analytic limiting arguments are written deductions, not machine-formalized
quantifiers. The predecessor and all inherited source packets replay.

The result provides a disciplined mechanism for collecting operator
information at one observation base while preserving the records needed
to diagnose concentration. A physical source/transport law selecting the
connectors, receiving geometry, record multiplicity and composition factor
is still required. Neither a universal gravitational source law nor a
mechanism that physically concentrates arbitrary information into a
singular point is derived here.
