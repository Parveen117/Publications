# CHD-2 — Minimal binary synthesis of positive quadratic records

15 October 2026. This finishes a finite-class representation question raised
by CHD-1: how many independent scalar channels are needed to compose a
Pythagorean record, and when do they remain sufficient after the next cut?

## 1. The closed target and its contracts

Let H be a positive semidefinite Hermitian d-by-d matrix over the native exact
field Q+iota Q, with rank r. Define q_H(x)=x*Hx. It is an invariant under a
specified source map T precisely when T*HT=H; positivity alone does not make
a quadratic reading conserved by arbitrary operations.
The Pythagorean summands need not be individually invariant when a transport
mixes channels; it is their specified total that is conserved.

We count **linear scalar amplitude channels**, retaining their signs/phases,
with positive scalar weights. A channel cannot conceal an arbitrary vector.
Binary synthesis means a tree of two-input additions of the squared weighted
amplitudes, equivalently Pythagorean combination of their nonnegative lengths.
This is not the number of classical bits encoding an arbitrary real amplitude.

**The finite-class result.** The minimum number of such channels realizing
q_H is r. They can be grouped by a binary tree with r-1 merges and maximum
depth ceil(log_2 r), when r>=1. The full amplitude record is faithful on the
quotient by ker H. Every H-preserving linear operation descends exactly to a
weighted isometry of that r-channel record. For additional non-preserving
cut operations, CRC's already established closure determines the larger
minimum record. Thus this target is settled without claiming every possible
invariant or every physical process belongs to the class.

## 2. Constructive weighted Pythagoras without evaluating roots

**CHD-T4.** There are a full-row-rank r-by-d matrix L over Q+iota Q and
positive rational weights delta_1,...,delta_r such that

    H=L* D L,  D=diag(delta_j),
    q_H(x)=sum_(j=1)^r delta_j |(Lx)_j|^2.                 (9)

Here and below the zero form has r=0 and an empty record.

Proof by exact completion of the square. For a nonzero positive semidefinite
residual H_k, choose a diagonal pivot a=(H_k)_{pp}>0 and set the row
ell=(H_k)_{p,:}/a. Then

    H_k = a ell*ell + H_(k+1).

The residual has zero p-th row/column. After moving p to the front, its
remaining block is the Schur complement C-b*b/a and is positive semidefinite:
substitute the first coordinate -b y/a into the original positive form.
Congruence with an invertible triangular matrix shows that the rank falls
by exactly one. A positive semidefinite matrix with every diagonal entry zero
is zero: positivity of each two-coordinate restriction forces every
off-diagonal entry to vanish. Thus the algorithm ends after exactly r pivots.
The pivot rows are independent and all arithmetic is native rational/complex
arithmetic. This is standard square completion/LDL factorization, reused
constructively, not a new factorization theorem.

Minimality: any realization H=sum_(j=1)^m delta_j ell_j*ell_j has rank at most
m. The construction attains m=r. Equivalently, over the positive square-root
extension one can write E=D^(1/2)L and q_H(x)=||Ex||^2. The executable uses
(D,L) and never approximates or silently adds a root to the native scalar type.

For a fixed input, let s_j=delta_j |(Lx)_j|^2. A binary tree combines the s_j
by exact addition. Its root is q_H(x); interpreting a node as a length yields
sqrt(s_left+s_right), the Pythagorean merge. A tree with r leaves has r-1
binary internal nodes, and a depth-k binary tree has at most 2^k leaves.
Recursively split the leaf set as evenly as possible to attain depth
ceil(log_2 r). For r=1 the depth and merge count are zero; for r=0 no tree
is needed. These are algebraic synthesis costs, not a physical measurement
depth or a claim that each subtree sum contains all its input information.

## 3. Which information the minimal record actually retains

Put y=Lx. Positivity in (9) gives ker L=ker H, so y is a faithful coordinate
on the quotient V/ker H. For r=d it reconstructs x; for r<d it does not recover
the missing kernel. Keeping only s_j, or just their sum, loses sign/phase
even when r=d. A nonlinear single channel sqrt(q_H(x)) preserves the scalar
value, but is outside the linear-channel minimum and is not a complete source
encoding. This qualification prevents a false universal bit-count theorem.

**CHD-T5 (invariant transport on the minimum record).** Suppose T*HT=H. For
x in ker H, (Tx)*H(Tx)=0, hence Tx in ker H. Therefore there is a unique
r-by-r operator U with

    L T=U L,       U* D U=D.                              (10)

Construct S=L*(L L*)^(-1), so LS=I, and set U=LTS. Kernel invariance shows
x-SLx lies in ker L, whence LT(x-SLx)=0 and
LT=UL. For every y choose x with Lx=y; then
y*U*DUy=(Tx)*H(Tx)=x*Hx=y*Dy, proving (10). Since D>0, U is injective and,
in finite equal dimensions, invertible. The original T may still discard
kernel information. Compositions satisfy U_(T2 T1)=U_T2 U_T1.

The normalized mathematical record D^(1/2)y has unitary transport
D^(1/2) U D^(-1/2). This is a derived representation of an already specified
H-preserving operation, not the selection of a physical evolution or a
clock. The construction uses the same quotient/descent principle as CRC.

For a finite list of other operations T_g, reuse CRC's exact theorem:

    R_* = span{ rows(L T_w) : finite operation words w },
    m_* = dim R_*.                                       (11)

A basis A_* of R_* is the smallest autonomous linear record containing y;
it obeys A_* T_g=U_g A_*. The native rowClosure and descendedAction routines
already compute it. If all T_g preserve H then m_*=r. Otherwise r<=m_*<=d,
and the full state is required when m_*=d. No new minimality claim is made
for the existing closure theorem; the addition is its explicit use with
the weighted Pythagorean synthesis and cut-history interface.

## 4. Complete coherent transport versus conserving a fixed branch metric

CHD-1 transports every represented operator by J(H)=W H W*, with exact
recovery and expectation preservation on the image. For a single cut P,
Q=I-P, W x=|0>Qx+|1>Px, a different question is whether the *same* H on each
final branch conserves the quadratic amount. Its pullback is

    W*(I_label tensor H)W=P H P+Q H Q
       =H-[P,[P,H]].                                    (12)

Consequently the fixed branch metric is conserved for every input iff

    [P,H]=0.                                            (13)

Proof: expand Q=I-P. If (12) equals H, sandwich the missing off-diagonal
part by P and Q to obtain PHQ=QHP=0; this is equivalent to commutation.
Conversely commutation makes both cross blocks vanish. No positivity
assumption is needed for this Hermitian operator identity.

Thus an isometric coherent cut lift always transports the complete observable
J(H), but does not automatically conserve the branch-local H. The transformed
observable generally contains cross-history blocks. Retaining them is exactly
the information condition in CHD-1. Copying H to every branch is a different
physical readout contract and must pass (13). This is the concrete compatibility
gate relevant to later geometry/matter sources.

## 5. Exact inherited source and minimum-record examples

Use the published YM54 thermo fixture already evaluated by HCM:

    H=[[2,1],[1,5]],
    L=[[1,1/2],[0,1]],  D=diag(2,9/2),
    q_H(a,b)=2|a+b/2|^2+(9/2)|b|^2.                      (14)

Rank is two: two scalar linear channels are necessary and sufficient, and
one binary merge evaluates the total. No decimal approximation of sqrt(2)
or sqrt(9/2) is required. The existing HCM factor B=[[2,1],[0,3]] obeys
B*B=2H and gives the equivalent representation (||Bx||^2)/2.

At x_+=(3/5,4/5) and x_-=(3/5,-4/5), the values are respectively 122/25
and 74/25. The axis cut P=diag(1,0), if read with the fixed H on each branch,
gives 98/25 for both. The omitted interference contributes +24/25 and
-24/25. In contrast, W*J(H)W=H preserves both source values exactly.

For a rank-one initial metric H0=diag(1,0), L0=[1,0], allow the already used
cut B0=u u*, u=(3/5,4/5). Native closure grows 1 -> 2. The source e2 and zero
both have initial L0-reading zero, but after B0 their squared H0-readings
are 144/625 and zero. The original single amplitude is not complete for that
next cut; two amplitudes are necessary and sufficient. This is a concrete
answer to the requested composition count, with the operation class retained.

The verifier also constructs an exact H-preserving rational transport and
its descended weighted isometry, plus a rank-deficient example whose kernel
cannot be recovered. These are mathematical witnesses, not physical time laws.

## 6. Root values: finite symbolic transport and Pythagorean generation differ

The positive sqrt(2) has the exact finite descriptor (t^2-2, 1<t<3/2).
It is also the Pythagorean length of (1,1). This descriptor can be transported
as finite binary data; its positional binary expansion does not terminate.
There is no need to throw away its irrational information or round it to a
rational value when the intended operations accept the exact descriptor.

The same distinction matters for general roots. Starting from rationals,
arithmetic and finitely many square-root/Pythagorean steps stay inside a
tower of extensions of degree at most two per step. Its degree over Q is a
power of two. By the degree tower law, every contained element has degree
dividing that power. The real cube root of 2 has degree three: t^3-2 has
no rational root (the only candidates are +-1,+-2), and a reducible cubic
over a field has a linear factor. Therefore it cannot be generated by those
square-root steps from rational inputs.

Nevertheless (t^3-2, 1<t<2) is a finite exact descriptor of that real root.
The verifier checks the native rational companion matrix C with C^3=2I
and independent I,C,C^2. This certifies the polynomial representation; it
does not make the matrix itself a chosen real root or a UGD numeral. Root
selection is supplied by the isolating interval. The field-degree argument
is a written standard theorem application, not a finite search proof.

Thus finite binary symbolic storage can handle a wider root class than
finite Pythagorean arithmetic. Unlimited arbitrary real inputs require an
additional computability/representation contract. Nothing here proves every
invariant has a finite binary expansion or fixes all UGD operations.

## 7. Completion verdict

Closed: exact weighted binary synthesis and minimum scalar-channel count for
every finite positive quadratic form in the stated field; exact descended
transport for H-preserving operations; and a finite algorithm for the minimum
autonomous linear record under additional specified cuts, reusing CRC.
The fixed branch-metric conservation gate is the exact commutator test (13).

Still separate: selecting physical cuts and gravitational coupling, deriving
a full UGD arithmetic representation, local quantum constraints/domains and
physical 3+1 gravity. This closes a mathematical representation/transport
problem, not quantum gravity. It does not assume every invariant is quadratic.

All generic statements have explicit written proofs. Native rational matrix
checks and serialized polynomial word witnesses certify the displayed finite
inputs and identities; they are not a general proof-assistant verification.
