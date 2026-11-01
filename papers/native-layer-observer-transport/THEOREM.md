# LOT-1 — Weighted observer transport and the resolved cut boundary

17 October 2026. Continues LCW-1, ILC-1 and LDC-1 on the same
completed xi source. The physical interpretation remains open. This
stage changes the carried channel pairing explicitly, obtains an
exact nested projection update, and identifies its compression with
the already certified LCW continuation position. It does not remove
LDC's obstruction for its different, fixed natural inclusion.

Retain c_0=1 and the source forms
H_N=(c_(i+j)), A_N=(c_(i+j+1)), B_N=H_N-A_N. All finite statements
assume A and B are strictly positive at every carrier used. The
natural polynomial inclusion E_(M,N) appends zero coefficients.
Transpose formulas below use real source coordinates; vector
quadratic forms extend with conjugation to complex states.

## LOT-T1. An exact channel pairing without numerical square roots

On pairs of degree-below-N polynomials define

\[
G_N=\operatorname{diag}(A_N,B_N),\qquad
Y_N=H_N^{-1}A_N,\qquad Z_N=\operatorname{diag}(I,-I),
\]
\[
\boxed{P_N=\begin{pmatrix}Y_N&I-Y_N\\Y_N&I-Y_N\end{pmatrix}.}
\tag{1}
\]

Then P_N is the G_N-orthogonal projection onto the diagonal pairs
(x,x). The operator Z_N is a G_N-self-adjoint involution, and
W_N=(2P_N-I)Z_N is G_N-unitary. All matrices in (1) are rational
expressions in the source forms.

For M>=N put V_(M,N)=diag(E_(M,N),E_(M,N)) and
V_(M,N)^sharp=G_N^(-1)V_(M,N)^T G_M. These maps satisfy

\[
\boxed{V_{M,N}^{\sharp}V_{M,N}=I,\quad
Z_MV_{M,N}=V_{M,N}Z_N,\quad
V_{L,M}V_{M,N}=V_{L,N}.}
\tag{2}
\]

**Proof.** The diagonal map D_N x=(x,x) is an isometry from the
H_N pairing to the G_N pairing because A_N+B_N=H_N. Its adjoint
is D_N^sharp=(Y_N,I-Y_N), so P_N=D_N D_N^sharp is the asserted
orthogonal projection. The two reflections are self-adjoint
involutions; their product is unitary. The source Hankel blocks
give E^T A_M E=A_N and E^T B_M E=B_N. These prove the isometry
in (2); the other identities follow by block multiplication. QED.

This is the ILC cut in another explicitly specified finite pairing.
Let the positive roots be taken in the acquired H_N pairing and set

\[
\mathcal R_N=\operatorname{diag}(\sqrt{Y_N},\sqrt{I-Y_N}),\qquad
\mathcal H_N=\operatorname{diag}(H_N,H_N).
\]

At finite strict depth, Rcal_N is invertible,
Rcal_N^T Hcal_N Rcal_N=G_N, and

\[
\mathcal R_NP_N\mathcal R_N^{-1}
=\begin{pmatrix}Y_N&\sqrt{Y_N(I-Y_N)}\\
\sqrt{Y_N(I-Y_N)}&I-Y_N\end{pmatrix}=\Pi_N.
\tag{3}
\]

Consequently J_(M,N)=Rcal_M V_(M,N) Rcal_N^(-1) is an isometry
between the finite ILC channel pairings. It composes across depths
and intertwines the fixed channel reflection Z. These facts require
no numerical square-root implementation: (1)-(2) give the exact
weighted-coordinate execution, while the positive roots establish
the mathematical equivalence (3).

## LOT-T2. The deeper cut adds one orthogonal boundary projection

For the adjacent inclusion V=V_(N+1,N), retain the LDC monic vector
and its positive squared norm

\[
v_N=\binom{-H_N^{-1}\boldsymbol h}{1},\qquad
\boldsymbol h=(c_N,\ldots,c_{2N-1})^T,\qquad
h_N=c_{2N}-\boldsymbol h^T H_N^{-1}\boldsymbol h>0.
\]

Write w_N=(v_N,v_N). Then

\[
Q_N=\frac{w_Nw_N^T G_{N+1}}{h_N},\qquad
\boxed{P_{N+1}=VP_NV^{\sharp}+Q_N.}
\tag{4}
\]

Here Q_N is a rank-one orthogonal projection, orthogonal to the
embedded old projection. In particular VP_NV^sharp<=P_(N+1).

**Proof.** The old diagonal polynomial space is contained in the
new diagonal polynomial space. Their induced pairing is precisely
the H pairing. LDC-T1 shows that v_N is orthogonal to all old
polynomials and has squared H_(N+1) norm h_N. Hence w_N is the
one missing diagonal direction, has squared G_(N+1) norm h_N,
and its projection is Q_N. Orthogonal decomposition of these
nested diagonal spaces proves (4). QED.

This preserves the old positive cut sector exactly:
P_(N+1)VP_N=VP_N. It does not intertwine the whole cut on all
old states. In fact P_(N+1)V-VP_N=Q_NV, and this term is nonzero.

For any old observable O_N, define the co-carried observable
O_(N+1)=VO_NV^sharp. For every old state z,

\[
\boxed{\langle Vz,O_{N+1}Vz\rangle_{G_{N+1}}
=\langle z,O_Nz\rangle_{G_N}.}
\tag{5}
\]

The extension in (5) is zero on the orthogonal complement of the
carried old state space; further observable action there can be
specified independently. In particular (5) does not identify the
newly measured full cut P_(N+1) with the co-carried old cut.

## LOT-T3. The new readout is exactly the LCW even position

Let ell=e_(N-1), d=(ell,-ell), and set

\[
\boxed{K_N=h_NG_N^{-1}dd^T
=h_N\begin{pmatrix}
A_N^{-1}\ell\ell^T&-A_N^{-1}\ell\ell^T\\
-B_N^{-1}\ell\ell^T&B_N^{-1}\ell\ell^T
\end{pmatrix}.}
\tag{6}
\]

The compression of the deeper cut is

\[
\boxed{C_N:=V^{\sharp}P_{N+1}V=P_N+K_N,\qquad
P_NK_N=K_NP_N=0.}
\tag{7}
\]

K_N is a nonzero positive semidefinite rank-one operator in G_N.
Its sole nonzero eigenvalue is

\[
p_N=h_N\ell^T(A_N^{-1}+B_N^{-1})\ell
=\frac{h_N}{\mathcal W_N}\in(0,1),\qquad
\mathcal W_N=U_N-L_N.
\tag{8}
\]

Here L_N,U_N are LCW-T1's sharp endpoints for the even moment
c_(2N). Thus p_N is exactly the LCW even continuation position,
not an additional source parameter. In particular

\[
\boxed{K_N^2=p_NK_N,\qquad
C_N-C_N^2=(1-p_N)K_N.}
\tag{9}
\]

For z=(x,t), the additional cut reading is

\[
\boxed{\langle Vz,P_{N+1}Vz\rangle_{G_{N+1}}
-\langle z,P_Nz\rangle_{G_N}
=h_N|\ell^T(x-t)|^2.}
\tag{10}
\]

**Proof.** Orthogonality of v_N gives
v_N^T H_(N+1)E=0. The multiplication identity in LDC-T1 gives
v_N^T A_(N+1)E=h_N ell^T, and consequently
v_N^T B_(N+1)E=-h_N ell^T. Thus
w_N^T G_(N+1)V=h_N d^T. Compressing Q_N in (4) yields (6)-(7).
Also G_NK_N=h_Ndd^T, which proves its positivity and rank.
The two rows of P_N coincide, so d^TP_N=0. Conversely
Y_NA_N^(-1)=(I-Y_N)B_N^(-1)=H_N^(-1), giving P_NG_N^(-1)d=0.
These prove both annihilation identities in (7).

The only nonzero eigenvalue of h_NG_N^(-1)dd^T is
h_Nd^TG_N^(-1)d. LCW-T1 gives
Wcal_N^(-1)=ell^T(A_N^(-1)+B_N^(-1))ell. Since
h_N=c_(2N)-L_N and strict next forms place c_(2N) strictly
between L_N and U_N, (8) follows. Rank-one multiplication now
gives (9); evaluating the compressed quadratic form gives (10).
QED.

Thus the compressed cut has an old eigenvalue-one sector, one new
eigenvalue p_N in the old negative sector, and its remaining
eigenvalues zero. It is not a projection under strict growth.
The normalized position in (8) is the classical scalar canonical
moment at this even index. The contribution here is its exact
identification with this source cut compression, derived from the
retained LCW and LDC identities.

## LOT-T4. Wave compression loss and wave transport error differ

Use W_N=(2P_N-I)Z_N and put T_N=V^sharp W_(N+1)V. Because V
intertwines Z, (7) implies

\[
T_N=(2C_N-I)Z_N=W_N+2K_NZ_N,
\]
\[
\boxed{I-T_N^{\sharp}T_N
=4(1-p_N)Z_NK_NZ_N.}
\tag{11}
\]

The compressed wave has one nonzero squared-norm-loss eigenvalue
4p_N(1-p_N), lying in (0,1]. The exact difference between carrying
the old wave and applying the new wave is a different operator:

\[
\mathcal D_N=W_{N+1}V-VW_N=2Q_NVZ_N,
\qquad
\boxed{\mathcal D_N^{\sharp}\mathcal D_N=4Z_NK_NZ_N.}
\tag{12}
\]

Its sole nonzero squared singular value is 4p_N. Thus natural
weighted wave intertwining fails at every strict adjacent extension,
even though the old positive cut sector and the co-carried readings
in (5) are retained exactly.

**Proof.** Both P and Z are self-adjoint in their specified metrics.
Therefore T_N^sharp T_N=Z_N(2C_N-I)^2Z_N; apply (9) to obtain
(11). Subtract the two wave actions, use Z_(N+1)V=VZ_N, and
then (4) to obtain Dcal_N=2Q_NVZ_N. Since Q_N^2=Q_N, its squared
norm operator is 4Z_NV^sharp Q_NVZ_N=4Z_NK_NZ_N. The spectra
follow from (8) and the unitary nature of Z_N. QED.

Compression loss alone cannot measure transport error: as p_N
approaches one the first expression tends to zero while the
second does not. No endpoint p_N=0 or p_N=1 is asserted for a
strict accepted source extension.

## LOT-T5. One added carrier resolves the cut action on old pairs

For any M>=N+1 with the intervening positive source forms, the
weighted cut action on an old pair z=(x,t) is the diagonal
polynomial pair

\[
P_MV_{M,N}z=(b,b),\qquad b(y)=yx(y)+(1-y)t(y),
\tag{13}
\]

interpreted in the M-dimensional polynomial basis. The degree of
b is at most N. Consequently

\[
\boxed{P_MV_{M,N}
=V_{M,N+1}P_{N+1}V_{N+1,N},\qquad
V_{M,N}^{\sharp}P_MV_{M,N}=P_N+K_N.}
\tag{14}
\]

**Proof.** Multiplication by y sends an old degree-below-N
polynomial into degree at most N. At every M>=N+1, the compressed
source multiplication Y_M therefore acts on that old polynomial
exactly. Substitute into (1) to obtain (13). This polynomial pair
already belongs to the (N+1)-carrier, proving the first identity
in (14). Adjoint composition and the source isometry give the
second. QED.

Equation (14) stabilizes this specified action on old polynomial
pairs. It does not identify every wave iterate, every source
derivative, or every observer with a finite-depth reading. The
state in (13) may occupy the new boundary direction even when
its compression back to the old carrier is unchanged by further
extensions.

## LOT-T6. Conditional weighted completion and endpoint scope

Assume now the additional all-order actual-source positivity
premise of ILC-T2. It supplies the unique compact moment measure
nu on [0,1]. Work intrinsically in the weighted Hilbert space

\[
\mathcal G_\infty=L^2(y\,d\nu)\oplus L^2((1-y)\,d\nu).
\tag{15}
\]

Each weighted carrier is the corresponding polynomial-pair
subspace. Define

\[
\boxed{P(f,g)=(b,b),\qquad b=yf+(1-y)g.}
\tag{16}
\]

This is the orthogonal projection onto the closed diagonal copy
of L^2(nu). The finite cuts, extended by zero off their weighted
carriers, are increasing orthogonal projections onto diagonal
polynomials. They converge strongly to P. On any fixed polynomial
pair of degree below N, the action is already exact at every
M>=N+1, as in (13).

**Proof.** Weighted equivalence classes make (16) well defined,
including at y=0 and y=1: the coefficient of any undefined channel
there vanishes. Pointwise Cauchy--Schwarz gives
|yf+(1-y)g|^2<=y|f|^2+(1-y)|g|^2, so b belongs to L^2(nu).
The diagonal embedding D:L^2(nu)->Gcal_infinity is an isometry,
and its adjoint is D^sharp(f,g)=yf+(1-y)g. Thus P=DD^sharp is
the asserted orthogonal projection. Polynomials are dense in
L^2(nu): continuous functions are dense for a finite Borel measure
on the compact interval, and uniform polynomial approximation
then gives the L2 approximation. Therefore the nested diagonal
polynomial ranges have closure equal to the full diagonal range.
Orthogonal projection onto nested dense ranges converges strongly:
approximate a vector in that range by a polynomial and use the
contractivity of every projection; vectors perpendicular to the
range are annihilated at every depth. This proves the limit.
The last assertion is the polynomial identity (13). QED.

No uniform invertibility of the finite Rcal_N is assumed in this
limit. Multiplication by diag(sqrt(y),sqrt(1-y)) maps (15)
isometrically into L^2(nu) direct_sum L^2(nu), but endpoint atoms
can prevent surjectivity onto that entire unweighted pair space.
In particular a first-channel value at y=0, or a second-channel
value at y=1, is absent from the relevant weighted space. The
intrinsic construction (15)-(16) does not invert those zero
weights or assume away their endpoint equivalence classes.

## Relation to the retained results and claim boundary

LDC compares the ILC cuts under E direct_sum E in the fixed
H_N direct_sum H_N pairing. LOT uses E direct_sum E in the
different, explicitly carried A_N direct_sum B_N pairing. Its
physical-coordinate inclusion is J=Rcal_M V Rcal_N^(-1), which
is generally different from LDC's natural inclusion. Hence (4)
does not contradict the LDC obstruction. LOT also retains its
own explicit whole-cut and wave defects, (7), (9), and (12).

The pinned LDC/LAA source jets at s=29,31,33 contain the data for
N=1,...,9 adjacent extensions. No new source derivative or matrix
square-root approximation is required for the rational identities,
the canonical positions, or their wave-defect scalars. The
infinite completion in LOT-T6 still depends on the unproved
all-order actual-source premise; finite nested projections do not
supply the missing positive source extensions.

Orthogonal projections, weighted isometries, compression defects
and canonical moments are classical tools. The normalized moment
coordinate is discussed, including its matrix version, in Holger
Dette and Jan Nagel, *Matrix measures, random moments and Gaussian
ensembles*, Section 3, <https://arxiv.org/abs/0904.3847>. No
priority claim is made for those general tools. The source-specific
relations above are proved directly from ILC-T3, LCW-T1/T2 and
LDC-T1 rather than importing an assumed positive xi-zero measure.
