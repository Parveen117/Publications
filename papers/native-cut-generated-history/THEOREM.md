# CHD-1 — From cut distinctions to recoverable operation histories

## 1. Contract and the change from CLD

T23 constructs objects from a cut P, its complement Q and an event lift.
CLD already constructs dynamics after an observation cut, but supplies arrows
T_j as independent data. Here the operation word consists solely of declared
binary cuts. No evolution of an uncut source, external clock, generator or
exponential is used to construct the word's record.

Work on an admitted finite positive carrier H, with its adjoint, identity and
scalar field, as in the pinned native operator representation. This carrier is
a representation contract, not a derivation of physical spacetime from nothing.
All numerical witnesses use Q+iota Q. Written statements extend to finite
complex carriers. No historical UGD arithmetic is modified.

There are three distinct inputs: (i) which observables are available, (ii) which
cut is applied next, and (iii) how that cut is implemented. A projection alone
does not fix a physical state update. We explicitly select the minimal
projective implementation E_1=P, E_0=I-P, without additional within-branch
rotations. This is the standard projective measurement/dilation contract, also
used in QGC; its elementary identities are rederived below.

## 2. Uncut observability, binary distinction and the absence of a clock

Define the uncut observable algebra A_empty=C I. Its number of independent
nonconstant linear coordinates is dim_C(A_empty/C I)=0. This use of zero does
not mean a zero-dimensional carrier or topological space.

**CHD-T1.** Every complex-linear unital map T on A_empty is the identity.
Indeed T(aI)=aT(I)=aI. Every normalized linear state has the same restriction
phi(aI)=a. Thus there is neither distinguishability nor detectable evolution
in this observable algebra. No conclusion about inaccessible underlying change
follows. An operation-order index is not required for this proof.

For an orthogonal cut P with P neither 0 nor I, set Q=I-P. Then

    P^2=P, Q^2=Q, PQ=QP=0, P+Q=I,
    A_P={aP+bQ} is isomorphic to C direct-sum C.

Proof: all polynomial words in P reduce to span{I,P}; the two nonzero
orthogonal projections are linearly independent. The isomorphism sends aP+bQ
to (a,b), with componentwise product and conjugation. The chosen outcome names
0 and 1 are labels for this distinction; this does not derive all binary logic
or arithmetic from an otherwise unspecified cut.

A repeated identical ideal cut is idempotent, not a compulsory clock: PP=P
and QP=PQ=0. A continuous family of automorphisms of A_P starting at the
identity is constant, since its only automorphisms are identity and the swap.
Thus a binary distinction enables nonconstant readings but does not by itself
force their change. A finite family of cuts generates a larger algebra by
composition; its dimension can be computed by the existing native closure.

## 3. Construct the history before discussing a flow

For a fixed ordered list of cuts P_1,...,P_n, let E_{j,1}=P_j and
E_{j,0}=I-P_j. For a binary word h=(h_1,...,h_n), define

    K_h=E_{n,h_n} ... E_{1,h_1},
    W_n x = sum_h |h> tensor K_h x.                         (1)

The record carrier is the orthogonal direct sum of one H for each word.
Labels name branches; W_n x retains their amplitudes and relative phases as
well as the final carrier. It is not merely a list of observed bits.
For n=0 the empty word has K_empty=I. This bookkeeping identity does not make
all source observables accessible before a cut: the available algebra there
remains the stipulated A_empty.

**CHD-T2 (complete coherent history).** For every finite cut word,

    sum_h K_h* K_h=I,       W_n* W_n=I.                    (2)

Proof by induction. For each prefix h, append b in {0,1}; the contribution
of its two extensions is

    sum_b K_h* E_{n+1,b}* E_{n+1,b} K_h
       = K_h*(P_{n+1}+Q_{n+1})K_h=K_h*K_h.

Sum over prefixes and use the induction hypothesis. The argument also allows
P_{n+1}=P_{n+1,h} to depend on the prefix, provided each branch uses a declared
orthogonal partition. Such controlled branching defines a coherent operation;
actually reading and classicalizing the controls changes the retained record.

Let V_{n+1} send |h> tensor y to sum_b |hb> tensor E_{n+1,h,b}y. Then

    V_{n+1}*V_{n+1}=I,      W_{n+1}=V_{n+1}W_n.           (3)

These are composable isometric arrows between increasing record carriers.
They are the operation-generated dynamics of this protocol. They do not
require or derive a one-parameter flow on a fixed carrier. The cut choices
remain specified protocol data, not an emergent gravitational field equation.

## 4. Complete transfer preserves operations, not only one numerical reading

Put R_n=W_n W_n*. This is the projection onto the coherent image. For any
input operator O define J_n(O)=W_n O W_n*. Then

    J_n(AB)=J_n(A)J_n(B),   J_n(A*)=J_n(A)*,
    J_n(I)=R_n,             W_n* J_n(O) W_n=O.             (4)

Proof: insert W_n*W_n=I. In particular J_n is injective and is a unital
*-isomorphism onto the corner R_n End(record) R_n, whose identity is R_n,
not the identity of the entire larger carrier. For a source vector or density,

    W_n*(W_n x)=x,    W_n*(W_n rho W_n*)W_n=rho.           (5)

The record extension intertwines these operations exactly:

    J_{n+1}(O)=V_{n+1}J_n(O)V_{n+1}*.

No separate pre-cut evolution F_s is supplied. What is transported is the
source's represented information through successive cut interactions.
On vector records the inverse has norm one: the reconstruction error from an
additive record error e is W_n*e, of norm at most ||e||. This is a finite
isometric encoding, not the inverse of HCM's contracting heat law.

This theorem applies to a UGD record only after a faithful representation of
that record, its equality and the intended operations has been supplied.
HCM's missing full UGD homomorphism and independent seam ledger are not solved
by renaming W_n x as a UGD numeral. Nor is an arbitrary real residual turned
into a finite classical data string.

## 5. Where classical measurement loses information

The density of the coherent record has off-diagonal history blocks
K_h rho K_k*. Reading histories as classical outcomes removes h != k blocks:

    M_n(rho)=direct_sum_h K_h rho K_h*.                    (6)

Discarding the label as well gives sum_h K_h rho K_h*. For a single cut,
this is P rho P+Q rho Q, and the missing part is P rho Q+Q rho P.
The missing operator is zero exactly when [rho,P]=0: multiplication by P
and Q isolates each off-diagonal block. This criterion concerns equality of
the original source with its dephased version; it is not a universal recovery
criterion for arbitrary restricted coding families.

**CHD-T3 (phase blindness persists after the first classical cut).** If
E_{1,b}rho E_{1,b}=E_{1,b}sigma E_{1,b} for b=0,1, then M_n(rho)=M_n(sigma)
for every later specified cut tree beginning with that first cut. Each final
branch sandwiches exactly one of those same first-cut blocks, proving the
claim. No subsequent processing of identical records distinguishes them.

Thus keeping every classical outcome and every final conditional carrier is
still weaker than keeping the full coherent history. Equation (5) requires
access to the coherences as well; it is not a decoder of outcome bits alone.
This is the standard measurement-dilation distinction, applied here to the
native complete-record question, not a newly discovered no-loss principle.

## 6. An order obstruction made solely from two cuts

Take cuts A and B. Match labels (a,b) in both orders and define

    Delta x = direct_sum_(a,b) (B_b A_a-A_a B_b)x,
    C=[B,A].

Each of the four commutators is C or -C. Consequently

    Delta*Delta=4 C*C,      ||Delta x||^2=4||Cx||^2.       (7)

This proves order independence for every input iff [B,A]=0. It is an exact
algebraic order obstruction, not automatically spacetime curvature. A detector
can still be blind to it: outcome probabilities need not determine Delta.
No cut schedule or universal force law is selected by (7).

## 7. Reuse the existing CLD fixture, without an intervening source flow

Use CLD's A=diag(1,0), u=(3/5,4/5)^T, B=uu* and input e_1. All intervening
source arrows are identity. Then

    B=[[9/25,12/25],[12/25,16/25]],
    [B,A]=[[0,-12/25],[12/25,0]],
    Delta*Delta=(576/625)I.                               (8)

The algebra dimensions generated by I, by A, and by A,B are respectively
1,2,4. The last is the full two-by-two matrix algebra; this is a representation
dimension, not an inferred spatial dimension.

For the actual projective A-B-A protocol, histories with the first A outcome
1 have probabilities

    111:81/625, 110:144/625, 101:256/625, 100:144/625.

The total is one; the final A=1 probability is 337/625. CLD already gives
the 111 value and distinguishes intervention from mere decomposition. Here
we retain all branches and explicitly certify the coherent decoder and the
order identity. The final A=0 probability 288/625 measures a possible change
of reading produced by the specified cut interactions, without a free source
Hamiltonian. Repeating A alone leaves the A reading unchanged.

For x_+=(3/5,4/5)^T and x_-=(3/5,-4/5)^T, the two pure densities are different,
yet their first A branch densities and every subsequent classical history
record are identical. Their complete coherent W_n encodings remain different
and decode exactly. The verifier also checks imaginary relative phases and
a genuinely complex cut to avoid relying only on real-matrix examples.

## 8. What this advances and what it does not

The chain now has a precise cut-first order: stipulated scalar observability;
binary cut algebra; composable cut records; exact coherent recovery; and a
measurable-information boundary when classical histories are extracted.
It supplies an explicit operation-preserving complete-transfer construction
without assuming uncut evolution. Standard dilation mathematics is credited;
the contribution here is its integration, source fixture, exact order budget
and audit of the proposed interpretation.

It does not prove that all changes in nature are caused by measurement, that
an observer creates an otherwise nonexistent carrier, or that lack of an
observable implies lack of hidden information. An exponential additionally
needs a suitable fixed-carrier composition and regularity contract; (3) alone
does not provide it. The pre-existing LNC/REC physics results remain intact.

For quantum gravity, the next concrete gate is to connect these cut-generated
arrows to the selected geometry/matter constraints and geometric observables,
with a derived or calibrated protocol law. This packet supplies no local
constraint quantization, physical 3+1 completion or gravitational prediction.
It also makes no prime generator or RH claim.

## 9. Verification boundary

The unchanged native matrix engine checks the source fixture, adaptive and
fixed history isometries, exact operator transport, recovery, phase loss and
order discrepancy. The unchanged native word engine proves projector and
finite-history identities by audited rewriting and serialized witness replay.
The arbitrary-depth induction and observability theorem are written proofs;
their general quantifiers are not certified by enumerating finite cases.
