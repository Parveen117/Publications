# CHD-3 — Complete circular-observer verification

16 October 2026. This closes the finite quadratic-observation verification
project begun in CHD-1/2. The input is a declared observation contract; the
output is an explicit, minimum-size separating record, an exact inverse,
and a necessary-and-sufficient transport certificate. No new physical law
is assumed in the proof.

## 1. What is being verified

Work on a finite n-channel carrier over Q+iota Q, iota^2=-1, with dagger
and a calibrated coordinate frame. The same arguments apply over C. A
quadratic observable is a Hermitian matrix H and its reading at a prepared
probe p is q_H(p)=p*Hp. Hermiticity, quadratic dependence, the carrier
dimension and calibration are hypotheses, not conclusions of finite data.
The candidate H, or the defect being tested, is not assumed invariant.

Observers may access distinct preparations or reproducible expectation
values. These are not simultaneous, noiseless observations of one unknown
quantum specimen. The individual probe projectors below need not sum to I;
each can be a separate binary cut measurement. No SIC-POVM is asserted.

Let a=3/5, b=4/5. The observer positions are

    p_i=e_i,
    p_ij^R=a e_i+b e_j,
    p_ij^I=a e_i+iota b e_j,              i<j.             (1)

All obey p*p=1 because a^2+b^2=1. The two-coordinate probes lie on unit
circles in the real planes spanned by (e_i,e_j) or (e_i,iota e_j). Each
P_p=pp* is a binary cut: P_p*=P_p, P_p^2=P_p. No claim says that every
possible observation geometry must be circular. Here Pythagoras specifies
the admissible probe geometry; the different observer positions provide
the distinguishing information.

## 2. CHD-T6 — Exact reconstruction and minimum observer count

For every Hermitian H, the n^2 readings (1) uniquely determine H. Write
d_i=q_H(e_i), r_ij=q_H(p_ij^R), s_ij=q_H(p_ij^I). Then

    H_ii = d_i,
    Re H_ij = (25 r_ij-9 d_i-16 d_j)/24,
    Im H_ij = (9 d_i+16 d_j-25 s_ij)/24,
    H_ji = conjugate(H_ij).                               (2)

Proof. Put H_ij=u+iota v. Direct expansion of p*Hp gives

    25 r_ij=9 d_i+16 d_j+24 u,
    25 s_ij=9 d_i+16 d_j-24 v.                            (3)

Solving these scalar equations yields (2). Conversely, (2) constructs a
Hermitian H for every real reading vector, and substituting it into (3)
recovers that vector. Thus the observation map is an isomorphism of real
vector spaces (of rational vector spaces in the native exact field).

Minimality. Hermitian n-by-n matrices have n diagonal real parameters and
two real parameters for each unordered pair: n^2 parameters. Any fixed
family of m scalar quadratic readings is a real-linear map from this
n^2-dimensional space to R^m. For m<n^2 it has a nonzero kernel, so it
cannot distinguish all H. Formula (2) attains the bound. In particular,
removing any one of our probes leaves an explicit invisible Hermitian
direction: reconstruct the reading vector which is one at that probe and
zero elsewhere. The executable constructs all these dual witnesses.

This lower bound concerns fixed scalar quadratic readings for an unrestricted
unknown Hermitian observable. It also holds on all positive definite forms:
if D is an invisible nonzero Hermitian direction, cI and cI+epsilon D are
both positive for sufficiently small epsilon and have identical readings.
It is not a lower bound for single-state pure-state tomography, adaptive
queries, or arbitrary nonlinear encodings. With known trace, one axis
reading can instead be supplied by trace minus the others; n^2-1 independent
scalar data suffice for arbitrary trace-fixed Hermitian matrices. If the
matrix is known real symmetric, the imaginary probes are unnecessary and
the exact minimum is n(n+1)/2.

This is a native rational specialization of standard polarization and
informational completeness, not a claim to have discovered tomography.

## 3. CHD-T7 — A finite certificate for all-input transport

Let T be a specified linear map from n input channels to m output channels,
and let G_in,G_out be specified Hermitian forms. Define

    D_T=T* G_out T-G_in.

Then the following are equivalent:

    q_(G_out)(T p)=q_(G_in)(p) for every probe in (1);
    D_T=0;
    q_(G_out)(T x)=q_(G_in)(x) for every input x.           (4)

Proof. The finite differences are exactly q_(D_T)(p). Equation (2) says
they all vanish if and only if the Hermitian D_T is zero. Matrix expansion
gives the final equivalence. This is a separating-family argument; it does
not assume that the tested T preserves the form. The same test certifies
T* B T=A for any declared input/output observables A,B, not just metrics.

If G_in and G_out are positive definite, (4) makes T injective and supplies
the explicit left decoder

    L=G_in^(-1) T* G_out,       L T=I.                    (5)

The full output amplitude record is recoverable. Equality of one scalar
norm reading alone does not reconstruct that record. For semidefinite
forms CHD-T4/T5 apply on the quotient by the invisible kernel, and CRC's
closure enlarges the record if later allowed cuts expose that kernel.
Linearity matters: x -> ||x|| e_1 preserves every Euclidean norm but loses
direction and phase. It satisfies the scalar invariant while falling
outside (4)'s linear-map inference of a decoder.

Two certified transports compose: if T1*G1 T1=G0 and T2*G2 T2=G1, then
(T2 T1)*G2(T2 T1)=G0. Thus the local certificates concatenate without
assuming return to the original point or a physical time parameter.

## 4. CHD-T8 — Complete verification survives coherent cut histories

Reuse the already proved CHD history isometry W*W=I. Set J(H)=W H W*.
At each transported observer position Wp,

    q_(J(H))(Wp)=q_H(p),       W*J(H)W=H.                (6)

Because (1) separates all input H, these transported readings recover H
by (2), and therefore recover the entire image-corner observable J(H).
They do not determine arbitrary operators on the unused output complement.
The source dimension n, not the larger history carrier dimension, determines
the n^2 count for this known image corner.

More generally an output B agrees with H on every source preparation iff
its n^2 transported readings agree, equivalently W* B W=H. This certifies
only the compression: off-image blocks are outside the target. Replacing
J(H) with the same H on every classical branch is a different contract.
CHD-2 gives W*(I_label tensor H)W=H-[P,[P,H]], so the observer certificate
passes exactly when [P,H]=0. The circle probes expose the lost cross terms.

For density or Gram records R=xx*, (2) recovers R exactly. It determines a
nonzero x up to one common unit phase: equality xx*=yy* gives equal
one-dimensional ranges, so y=cx and |c|=1. Quadratic observations cannot
recover that global phase. CHD's coherent amplitude decoder (5) assumes
the amplitude record itself is retained; it is not inferred from one
classical radius or from an unknown single-copy measurement.

## 5. CHD-T9 — Frame covariance and finite-error verification

For a declared invertible frame change S, set

    p'=S p,      H'=S^(-*) H S^(-1),
    G_frame'=S^(-*) S^(-1).

Then p'*G_frame'p'=1 and p'*H'p'=p*Hp. These are the same circular
positions relative to the transported reference metric; a nonunitary S
need not preserve their Euclidean coordinate norm. Inverting the readings
in the reference coordinates and applying the stated congruence recovers
H'. The observer, the form and calibration transform together. Choosing an
unknown new calibration independently for each reading is not licensed.

If each measured reading differs from its exact value by at most epsilon,
the linear reconstruction has diagonal error at most epsilon, and real
and imaginary off-diagonal errors at most 25 epsilon/12. Consequently

    ||H_hat-H||_op <= C_n epsilon,
    C_n=1+25(n-1)/6.                                     (7)

Proof. Apply the triangle inequality to (2). Each complex off-diagonal
entry has magnitude at most the sum of those two component bounds,
25 epsilon/6. The maximum absolute row sum is at most C_n epsilon;
for a Hermitian matrix its operator norm is bounded by that row sum
(for example apply the eigenvector equation at a maximum component).
For a tested transport defect the same bound gives
|x*D_T x| <= C_n epsilon ||x||^2 if all exact probe defects are bounded
by epsilon. Measurement uncertainty and calibration errors must be included
in epsilon. Finite-shot estimates do not prove an exactly zero defect.

## 6. Exact source demonstrations and rejected alternatives

The inherited YM54 form H=[[2,1],[1,5]] has readings

    (d_1,d_2,r_12,s_12)=(2,5,122/25,98/25).

Equation (2) reconstructs its off-diagonal value 1+iota 0. After the
axis cut, the fixed branch metric reads (2,5,98/25,98/25); its reconstructed
defect is [[0,-1],[-1,0]]. Its circular real probe fails by -24/25 even
though both axis readings pass. The coherent J(H) passes all four readings.
The verifier repeats this on a two-cut coherent history, including complex
Hermitian basis elements, not just this source matrix.

The positive complex fixture [[2,1+2iota],[1-2iota,5]] has readings
(2,5,122/25,2), reconstructing both components exactly. Real-only probes
cannot distinguish states (3/5,+4iota/5) and (3/5,-4iota/5), but the
imaginary circular probe distinguishes their Gram records. A common global
phase does not change any quadratic reading, as required.

Norm preservation is not flatness or a gravitational law. The exact unitary
maps U=[[3/5,-4/5],[4/5,3/5]] and V=diag(1,iota) give a nontrivial loop
K=V*U*VU. Both I and K preserve every norm; K changes a cut observable.
The full observer family detects K*P_0K-P_0. This is an algebraic loop
comparison, not a derived spacetime curvature or a faster-than-light signal.

## 7. Completion and proof scope

Closed project: finite quadratic-observer information verification, including
an explicit inverse, the exact minimum fixed probe count, a finite iff
certificate for all-input linear invariant transport, decoder, coherent
history reuse, common-frame covariance and a deterministic error bound.
CHD-1 supplies the retained history; CHD-2 its minimum weighted amplitude
record; CHD-3 supplies a complete observer certificate for its observables.

The hypotheses specify what verification means. They do not certify their
own physical applicability. Purely algebraic transport is verified within
this declared contract without assuming a universal absolute origin.
The n^2 count concerns independent scalar data, not bits or experiments
needed to estimate arbitrary real numbers exactly.

General results above are proved by written algebra, dimension arguments
and inequalities. The unchanged native engines separately certify rational
finite fixtures, basis completeness, every omitted-probe witness, symbolic
identities and serialized replay. This is not general proof-assistant
verification. A gravitational action, physical branch law, local quantum
constraints and experimental speed comparison remain separate projects.
