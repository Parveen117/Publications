# CRC-1 — Cut records, operation closure and hidden source-law information

14 October 2026. Native cut/operator development, with dimensionless
coefficients and no primitive Hilbert-space premise.

The motivating claim is that an invariant carrier must retain the
information needed by its next operation. This packet makes that claim
target dependent and constructive. It identifies missing cut direction,
distinguishes count from rotation history, and applies the same reasoning
to an otherwise hidden quadratic source coupling in the gravity model.

The general observer-kernel, retention and minimal-realization ideas
already exist in the programme. They are explicitly reused below.
The additions are the particular cut-readout closure, its shortest
displayed distinguishing protocols, and the nonlinear source-law witness.

## 1. What it means for an invariant record to be complete

Specify an admitted source set X, an encoding E:X->C, allowed operations,
and the readings to be predicted. Let R_w(x) be the reading after a
specified finite operation word w, including the empty word.

**Factorization criterion.** There is a well-defined prediction F_w on
the encoded record with R_w=F_w composed with E for every allowed w if
and only if

\[
\boxed{E(x)=E(y)\ \Longrightarrow\ R_w(x)=R_w(y)
\quad\text{for every admitted }x,y,w.}
\tag{1}
\]

Proof: factorization implies equal outputs on each fibre of E. Conversely,
define F_w(E(x))=R_w(x); (1) makes the definition independent of which
source in that fibre is chosen. This is an elementary quotient argument,
not a new general theory of sufficiency.

The equivalence x~y iff all R_w agree is the coarsest operational record
for this question. A faithful record of the full source may contain more
information than that quotient needs. Enlarging the numerical value of
an invariant does not refine its fibres.

For instance a/b cannot distinguish (2,1) from (20,10). An operation
reading a+b distinguishes them. If b!=0, retaining r=a/b and s=b
reconstructs (a,b)=(rs,s). Scale is either retained elsewhere in the full
carrier or absent from the normalized record; transport does not recreate
a distinction that its input and operation law do not contain.

An exact pair violating (1) certifies incompleteness under the stated
model. An unexplained experimental error alone does not: calibration,
numerical error, noise or an incorrect operation law must also be considered.

## 2. Reuse the native operation-closure engine

On a finite real coefficient space V, let the admitted fixed operations
be linear matrices T_g and the initial record be A0 v. Define

\[
W_0=\operatorname{row}(A_0),\qquad
W_{k+1}=W_k+\sum_g W_kT_g.
\tag{2}
\]

The existing native `rowClosure` routine computes this ascending chain.
At equality W_(k+1)=W_k it has stabilized, since multiplication by every
generator preserves the space. There are at most dim(V)-rank(A0) strict
rank increases. Its final space W_* is exactly the span of all rows
A0 T_w for finite allowed words.

**Closure and minimality.** A full-row-rank matrix A containing the seed
readings is an autonomous sufficient linear record precisely when there
are matrices C_g satisfying

\[
\boxed{AT_g=C_gA\quad\text{for every generator }g.}
\tag{3}
\]

Equivalently ker(A) is invariant under every T_g. Then all future records
are computed by multiplying the C_g, without choosing an underlying source.
Every such A must contain W_*. A basis of W_* attains the minimum dimension.

Proof: (3) implies invariance of the kernel. Conversely, if ker(A) is
invariant, define C_g(Av)=AT_gv, which is well defined by the same fibre
argument as (1). Induction shows any closed row space containing W0
contains each Wk. At stabilization the computed space itself is closed.
The native `descendedAction` routine tests and constructs C_g.

If an allegedly sufficient record A misses a future row b, then b is
not in row(A). Linear algebra supplies d with Ad=0 but bd!=0. This is
a distinguishing direction; adding b repairs that particular failure
with one more independent row. The certificate includes actual cut-source
pairs, rather than only abstract coefficient differences.

This minimality concerns real linear records on the declared source
space. It does not compare arbitrary nonlinear encodings or establish
an infinite-dimensional closure. For a positive source cone with nonempty
interior, missing linear directions also give admitted positive witnesses
by perturbing a sufficiently positive scalar matrix in both directions.
Fixed trace, fixed rank or other prior information must be accounted for
before interpreting the number of independent unknowns.

## 3. How much information the current cut sector needs

Use the same native matrices as EGR-1/EGR-2:

\[
Z=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
X=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad Y=\iota R.
\tag{4}
\]

Here iota is the existing central cut-field scalar; Y is self-adjoint.
This is the scalar extension already typed in NJC-1, not an identification
of an anti-linear plane operation with a central scalar.

For a self-adjoint native carrier rho, define

\[
N=\operatorname{Tr}\rho,\quad z=\operatorname{Tr}(Z\rho),\quad
x=\operatorname{Tr}(X\rho),\quad y=\operatorname{Tr}(Y\rho).
\]

Trace orthogonality gives exact reconstruction

\[
\boxed{\rho=\tfrac12(NI+zZ+xX+yY).}
\tag{5}
\]

For the real cut source rho=NP of the gravity packets, y=0 and
z=N r1, x=N r2. A positive rank-one source obeys z^2+x^2=N^2 in that
sector. Thus three polynomial record coordinates are not three independent
continuous parameters on that rank-one surface. N is a source weight;
it is not identified with entropy or every notion of quantum information.

The coordinates are relative invariants: transform rho and its calibrated
readout operators Z,X,Y together under the same native unit frame and
the traces do not change. Transforming only rho is a different source
orientation relative to a fixed observer.

For real-plane transport U=cI+sR with c^2+s^2=1,

\[
\begin{pmatrix}N'\\z'\\x'\\y'\end{pmatrix}=
\begin{pmatrix}
1&0&0&0\\
0&c^2-s^2&-2cs&0\\
0&2cs&c^2-s^2&0\\
0&0&0&1
\end{pmatrix}
\begin{pmatrix}N\\z\\x\\y\end{pmatrix}.
\tag{6}
\]

If the receiving cut is Q=(I+aZ+bX)/2, a^2+b^2=1, the exact reading is

\[
\boxed{\operatorname{Tr}(QU\rho U^\dagger)
=\tfrac12\{N+[a(c^2-s^2)+2bcs]z
+[b(c^2-s^2)-2acs]x\}.}
\tag{7}
\]

Therefore (N,z,x) is sufficient for every real-plane rotation and real
receiving cut. Even an unknown y is invisible to this entire operation
and readout class. It must not be inferred from those readings.

Take the seed record (N,q_Z), where q_Z=Tr((I+Z)rho/2), and the exact
transport U=3I/5+4R/5. Native row closure has ranks **2 -> 3**, with final
basis (N,z,x). The source cuts

\[
P_\pm=\tfrac12(I+\tfrac35 Z\pm\tfrac45 X)
\tag{8}
\]

both have N=1 and initial q_Z=4/5. After this same U their readings are
respectively 4/125 and 4/5. One rotation exposes the omitted x component.
No change of source amount or spectrum is needed.

### How many cut settings?

Let q_plus=Tr((I+Z)rho/2), q_minus=Tr((I-Z)rho/2), and
q_x=Tr((I+X)rho/2). Then

\[
\boxed{N=q_{plus}+q_{minus},\quad z=q_{plus}-q_{minus},\quad
x=2q_x-N.}
\tag{9}
\]

Two calibrated cut orientations, Z and X, supply these three scalar
readings when both outcomes of the Z split are retained. If N is already
known, q_plus and q_x suffice for the real source. These are exact trace
or ensemble-readout statements. They do not claim state reconstruction
from a single quantum specimen or ignore measurement disturbance.

## 4. When a fourth channel becomes necessary

Allow the full self-adjoint carrier (5), with unknown y, and add the
native unit transport

\[
V=\tfrac35I+\tfrac45\iota X.
\tag{10}
\]

It mixes the Y component into the Z reading. The two positive rank-one
cuts rho_plus=(I+Y)/2 and rho_minus=(I-Y)/2 have the identical real record
(N,z,x)=(1,0,0). No real-plane continuation in section 3 separates them.
After V their Z-cut readings become

\[
\boxed{q_{Z,+}=1/50,\qquad q_{Z,-}=49/50.}
\tag{11}
\]

Closure for the same seed under both U and V has ranks **2 -> 4** and
retains (N,z,x,y). A Y cut supplies y=2q_y-N. Thus three calibrated cut
orientations Z,X,Y give the four scalar coordinates when N is unknown
and both outcomes of the first split are retained.

The fourth requirement depends on admitting an unknown Y component in
the source class. If y=0 is known initially, transport creates no new
unknown initial parameter: (N,z,x) already reconstructs that source, and
the generated fourth coordinate can be calculated. A current-state
record that forgets this generated coordinate need not remain closed
under the enlarged operation class.

This is the precise sense in which the next allowed operation decides
which previously hidden information matters. It is an operational
distinction, not a claim that a larger numerical invariant stores more
information. The closure engine emits an operation word for each added
independent reading and verifies its action exactly.

## 5. Cut direction, reflection products and rotation

For a real rank-one cut,

\[
P_\phi=e^{\phi R}P_0e^{-\phi R},\quad
J_\phi=2P_\phi-I=\cos(2\phi)Z+\sin(2\phi)X.
\tag{12}
\]

The native relations Z^2=X^2=I, XZ=R=-ZX give

\[
\boxed{J_\beta J_\alpha=e^{2(\beta-\alpha)R}.}
\tag{13}
\]

Proof: expand the product. Its identity coefficient is
cos(2beta)cos(2alpha)+sin(2beta)sin(2alpha); its R coefficient is
sin(2beta)cos(2alpha)-cos(2beta)sin(2alpha). The completed scalar addition
identities give (13). The direct conjugation transporting P_alpha to
P_beta has angle beta-alpha. The factor two in (13) belongs to the
product of reflections. Multiplying measurement projectors instead is
not generally unitary and is not this operation.

The cut orientation is modulo pi; the rotation matrix is modulo 2pi.
Rank and spectrum do not fix either relative orientation. The executable
checks native reflection identities and rational relative-rotation products,
including reversal of order. This instantiates the standard reflection
algebra; it is not a new general theorem about rotations.

## 6. Count and history are necessary only for the target that uses them

The empty word and the word ZZ both have endpoint I, with operation
counts 0 and 2. The words ZZ and XX have the same endpoint and the same
count, but different intermediate cuts. These are exact witnesses:

* endpoint alone does not recover count;
* endpoint plus count does not recover the full word;
* count alone does not specify relative angles or order.

For an order-sensitive reading, take J1=Z and J2=3Z/5+4X/5. Then
J2 J1=U, J1 J2=U^dagger, and both have trace 6/5 and count two. Starting
from rho=(I+X)/2, their final Z-cut readings are 1/50 and 49/50.
Thus count and unsigned rotation data can miss an observable distinction.

Conversely, if two histories yield the same final carrier rho_final,
then every identical subsequent fixed transport V and receiver Q gives
the same Tr(Q V rho_final V^dagger). Count cannot secretly alter this
formula. It becomes relevant only if the requested target explicitly
reads count or if the operation law includes a memory variable.

For a declared continuous real rotation path with fixed initial angle,
an endpoint U determines its lifted angle only modulo 2pi. Keeping a
real angle lift retains the accumulated angle; four positive quarter-turns
return U to I while their lift changes by 2pi. A cut returns already
after a half-turn. Discrete operation count supplies neither a continuous
interpolation nor its winding without a path contract. Even an angle lift
does not reconstruct the entire path's parameterization.

The existing NM1 and N4H-1 results already establish endpoint/lift loss,
and QG-R2 establishes retained versus fused loop information. Here these
are used to decide whether a particular next-reading target needs a
ledger. There is no universal instruction to append all history to every
source carrier.

## 7. Gravity: same source, a hidden quadratic coupling, different response

Now apply the completeness criterion to EGR-2's finite metric problem.
Retain its geometry, boundary conditions, source cut r=(3/5,4/5), amount,
angle lifts, matrices A,S,G,K, full third derivative T and positive gamma.
Write h=q-q0, epsilon=alpha*N/gamma, and B10=[S G]. Define

\[
\Gamma=\text{last eight rows of }B_{10}^{-1},\qquad
\Gamma G=I_8,\quad\Gamma S=0.
\tag{14}
\]

The inherited null-lift matrix is

\[
M(r)=\begin{pmatrix}-10r_1-2r_2/5&2r_2\\2r_2&-10r_2\end{pmatrix}
\otimes I_4.
\tag{15}
\]

To test identifiability, explicitly compare the source-law family

\[
\boxed{F_\kappa(h)=r^TAh+
\frac\kappa2(\Gamma h)^TM(r)(\Gamma h),\qquad
\mathcal S_\kappa/\gamma=\mathcal R(q)-\epsilon F_\kappa(h).}
\tag{16}
\]

Kappa is a declared dimensionless coupling parameter for this comparison.
It is not inferred from the source invariants and is not being selected
as a new physical law. Kappa=0 is exactly EGR-2's original action. All
members have F_kappa(0)=0 and the identical first derivative

\[
F'_\kappa(0)=A^Tr,
\quad F''_\kappa(0)=\kappa\Gamma^TM(r)\Gamma.
\tag{17}
\]

Thus the full source operator NP, its scale, direction and spectrum,
the source gradient, and the first-order curvature readout are identical.
Their source-to-geometry laws differ at quadratic order. That difference
is model information, not another component already implied by NP.

### Exact compatibility and local continuation

Let z0=K^{-1}r, d0=S z0, b=G^T T[d0,d0], and let eta0 be the EGR-2
solution b/2+M eta0=0. For an expansion h=epsilon d+epsilon^2 e+...,
the first equation still gives d=d0+G eta. The next equation is

\[
H e+\tfrac12T[d,d]-F''_\kappa(0)d=0.
\tag{18}
\]

EGR-2's pure-null cubic identity reduces its eight compatibility
conditions to

\[
\boxed{\tfrac12b+(1-\kappa)M\eta=0,
\quad\eta_\kappa=\frac{\eta_0}{1-\kappa}\quad(\kappa\ne1).}
\tag{19}
\]

The source term F''_kappa d0 vanishes because Gamma S=0. The readout's
second coefficient is

\[
z_{2,\kappa}=K^{-1}S^T
\left[F''_\kappa(0)d_\kappa-\tfrac12T[d_\kappa,d_\kappa]\right].
\tag{20}
\]

The certificate checks (18) in all ten edge components for e=S z2.
As before, the actual branch can have an additional G component in e;
it does not affect A e=z2.

For kappa!=1 and r1*r2!=0, the same desingularized analytic equations
as EGR-2 have Jacobian blocks K and (1-kappa)M. Both are invertible.
The analytic implicit-function theorem therefore supplies an actual local
branch, with all cells nondegenerate for sufficiently small source.
This argument retains the prior positivity and angle-branch hypotheses.

At kappa=1, b is nonzero for the specified r, while (19)'s eta coefficient
vanishes. No branch differentiable at epsilon=0 with h=O(epsilon) can
satisfy that compatibility equation. This does not exclude nonanalytic
scalings, a remote solution or a branch based at a different flat point.

### Source-law ambiguity changes stability

The eight-dimensional metric Hessian block becomes

\[
\boxed{G^T\left[\mathcal R''(q(\epsilon))
-\epsilon F''_\kappa\right]G
=\epsilon(1-\kappa)M+O(\epsilon^2).}
\tag{21}
\]

The top two-dimensional block is K+O(epsilon), and the cross block is
O(epsilon). The Schur-complement inertia argument from EGR-2 applies.
At the fixed positive source direction (3/5,4/5), M consists of four
negative-definite blocks of determinant 48. Therefore:

| Source-law parameter | Same first readout coefficient | Small-positive-source metric Hessian |
|---|---|---|
| kappa=0 | (3/50,4/235) | 2 positive, 8 negative: saddle |
| kappa=2 | (3/50,4/235) | 10 positive: strict local minimum |

The first area-weighted curvature derivative is also unchanged, since
the predecessor proves it factors through A d. The exact second
readout coefficients are different:

\[
z_{2,0}=\binom{521098203253/22485570048000}{336288173/11711234400},
\quad
z_{2,2}=\binom{2628357123901/22485570048000}{2771105167/11711234400}.
\tag{22}
\]

The initial logarithmic cell-volume rates likewise distinguish the laws:

| Cell class | kappa=0 | kappa=2 |
|---|---:|---:|
| Four outer cells | 163103/848256 | 1729253/4241280 |
| One cap cell | 4410869/1060320 | -4410869/1060320 |
| Four inner cells | 912209/1060320 | 1975471/1060320 |

The local minimum at kappa=2 contracts the cap and expands the other
cells. Stability does not itself imply contraction everywhere. This
comparison does not modify the EGR-2 complementary-cut result: that result
remains true for its specified kappa=0 action and negative source direction.

## 8. The resulting record rule

The record needed for a prediction follows its target:

| Target | Required data in the stated example |
|---|---|
| Real-plane next cut reading | N,z,x and calibrated transport/receiver |
| Phase-sensitive next reading with unknown Y component | N,z,x,y and the expanded protocol |
| Final carrier after a known history | Its actual endpoint transport; a full word is optional for this target |
| Number or intermediate sequence of cuts | A corresponding operation ledger |
| Accumulated continuous rotation | A declared path and appropriate angle lift |
| EGR-2 nonlinear metric response | Source, geometry, and the specified nonlinear coupling law |

For the metric response through the next order, a relevant additional
law datum is F''(q0); its projection G^T F''(q0)G controls the leading
null-direction lift. Other components of that derivative can also affect
the response. In the constructed family, the missing distinction is the
single parameter kappa. Higher response orders can require higher
derivatives of the law. The finite carrier dimension does not bound
the information needed to specify an arbitrary nonlinear law.

There is no implication that any large invariant is universally complete,
that every residual is a hidden physical variable, or that count alone
creates curvature. A distinguishing exact protocol proves a specific
omission; closure under the stated protocols certifies its repair.

## 9. Lineage and certification boundary

| Existing source | Reused result | Addition in CRC-1 |
|---|---|---|
| Native operator engine | `rowClosure`, `descendedAction`, exact matrix algebra | Cut-transfer adapter and explicit operation witnesses |
| NJC-1 T1/T7 | Cut reflection and kernel-target criterion | Concrete direction/phase repairs |
| NM1 T2/T3 | Endpoint fibres, continuation blindness, count versus endpoint | Target-dependent decision on which ledger is needed |
| GE4 | Minimal observable realization depends on an admitted source class | Fixed-transport cut records, not a new moment reconstruction theorem |
| QG-R2 | Retained records versus cancelling fused readings | Distinguish retained path targets from endpoint targets |
| EGR-1/EGR-2 | Geometry, full T tensor, response and local continuation | Same-source quadratic-law ambiguity and exact stability comparison |
| N4H-1 | Lifted curvature can differ at the same compact endpoint | Preserve lift data when the requested gravitational target uses it |

Finite-dimensional observability and state reconstruction have established
literature. For relevant primary context see D. D'Alessandro,
[On Quantum State Observability and Measurement](https://arxiv.org/abs/quant-ph/0307127),
J. Phys. A 36 (2003) 9721-9735. Its source class uses finite quantum
control dynamics and indistinguishable-state classes. Here the elementary
proofs and concrete coefficients are given within the declared native
carrier; no general novelty claim or outside arithmetic engine is used.

The executable verifies source pins, native algebra identities, real
closure ranks, descended updates, reconstruction and frame covariance,
count/order examples, the inherited geometric tensor contractions,
all ten second-order metric equations and stability inputs. Symbolic
identities carry serialized Paaninian witnesses and independent replay.
The certificate is regenerated in a separate read-only check process.

General all-word induction, arbitrary real-angle identities, quotient
factorization, analytic existence and sufficiently-small-source inertia
have written proofs. Exact finite tests do not replace those arguments
or provide a generic proof-assistant formalization. The work selects no
physical source law, proves no global singular collapse, and adds no
continuum quantum measure. All earlier packets and engines remain unchanged.
