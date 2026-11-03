# Native response tensor and the cut ledger

Monty Dabas. Working research programme v0.8, 25 September 2026.

## 1. What the repository audit changes

The information/curvature programme already exists. CID-1 gives total
covariance and a metric-independent residue in a declared response algebra;
RKF 32 constructs cut covariance from an independently specified event lift;
EMK-T1 contains mixed-channel commutators; QTH-1 computes an SLD Gram tensor
for declared qubit states. These are predecessors, not results first found here.
The [ten-source audit](RESPONSE_SOURCE_AUDIT.md) and
[immutable pins](RESPONSE_SOURCE_PINS.json) state exactly what is reused.

The new task is to construct a response tensor from U27's native represented
transports, identify its cut compatibility, and attach the inherited covariance
ledger to an explicit signed response. No spacetime, density operator, Born
rule, thermodynamic equation of state or classical field equation is an input.
The finite carrier is still a presentation with distinctions, not uncut ground.

## 2. Declared datum and operational target

Use U27's real carrier W=R^(2m), represented quarter-turn
Z=I_m tensor R, native grading J=I_m tensor K and positive form H=I.
Thus Z^2=-I, Z^T=-Z, JZJ=-Z. Each selected transport U is orthogonal and
commutes with Z. Z is a represented operator, not identified with the native
central scalar iota. Complex notation below is algebraic bookkeeping.

For two calibrated controls U,V define the **signed order response**

    L = UV - VU,                 r(x)=Lx.

The first term applies V then U; the second applies U then V. To read a
component one prepares the same input under each word and subtracts the signed
outputs. The preparation/repetition and signed calibration contract remains
explicit; it is not yet a physical instrument. U29 supplies reconstruction
within that model contract. Graph, strengths and preparation labels are inputs.

One may also use the group-loop displacement UVU^-1V^-1-I. The identity

    UVU^-1V^-1-I = (UV-VU) U^-1 V^-1

relates it to L without assuming a chart. They have the same rank, but have
different input conventions and cannot be silently substituted on a fixed x.

Choose calibrated response columns Y=[y_1,...,y_k], for example y_a=Lx_a.
They are constructed before any tensor norm or desired positivity is evaluated,
as required by RKF 32's non-circularity rule.

## U30. Native Hermitian response tensor and representation covariance

Define

    Q(Y)=G(Y)+i A(Y),
    G=Y^T Y,                  A=-Y^T ZY.

Then G is real symmetric, A is real antisymmetric, and Q is Hermitian positive
semidefinite. Its real part measures squared response size; A is an oriented
pairing of responses. These descriptions do not identify either part with a
physical information metric or with connection curvature.

For every common orthogonal transport U commuting with Z,

    Q(UY)=Q(Y),               Q(JY)=conjugate(Q(Y)).

For a real change of response labels B, Q(YB)=B^T Q(Y) B. For an arbitrary
invertible carrier frame F, set H'=F^-T F^-1, Z'=FZF^-1 and Y'=FY. Then

    Y'^T H'Y' - i Y'^T H'Z'Y' = Q(Y).

Keeping a numerical Euclidean metric fixed under a nonorthogonal change of
units does not obey this identity. This is the unit-typing obligation in the
thermodynamic paper, not a new physical covariance postulate.

**Proof.** Identify each real pair (p,q) with p+i q. Since multiplication by
i corresponds to R, the complex Gram product of the resulting columns is
exactly Y^T Y-iY^T ZY. For any complex coefficient vector a, a*Q a is the
squared norm of the combined complex response, hence nonnegative. Transpose
gives the stated symmetries. Orthogonality, UZ=ZU and JZJ=-Z prove the two
transport identities. The other formulas follow by substitution. No quantum
interpretation is used. In particular, |Q_ab|^2 <= G_aa G_bb is a Gram bound.

## U31. Exact real-cut ledger, its seam term and minimum repair

Let C act on the **output response carrier**. In this section a cut is a
real linear readout; let P be the orthogonal projector onto row(C), Q0=I-P.
Use the ambient metric to define the visible and hidden tensors from PY and
Q0Y. Then

    Q(Y)=Q(PY)+Q(Q0Y)+i S_C(Y),
    S_C(Y)=-Y^T(P Z Q0+Q0 Z P)Y.

The real symmetric ledger always adds with no cross term. The full complex
ledger adds with no cross term **for every response catalogue Y** if and only if

    [P,Z]=0.

Equivalently, ker(C) is Z-invariant; equivalently, Z has an exact descended
action on im(C). A particular small Y may fail to detect a nonzero seam term;
the universal criterion must not be inferred from one zero reading.

The coarsest linear refinement of C supporting that action has row space

    row(D)=row(C)+row(CZ),
    minimum extra independent real channels = rank([C;CZ])-rank(C).

For the zero readout the minimum is zero. Every nonzero completed row space
has even real dimension. This is the prior observer-closure construction
specialized to Z^2=-I, not a new general minimum-memory theorem.

**Proof.** Expand Y=PY+Q0Y in the Gram formula. Orthogonality kills the real
cross terms, leaving the displayed imaginary term. If it vanishes for every
Y, take Y=I: PZQ0+Q0ZP=0. Multiplication on the left by P and right by Q0
gives PZQ0=0; its transpose gives Q0ZP=0. Thus Z is block diagonal and
commutes with P. The converse follows immediately. For a skew-adjoint Z,
invariance of ker(C) also implies invariance of its orthogonal complement,
so the equivalences follow from the exact descent gate.

The sum row(C)+row(CZ) is invariant under right multiplication by Z since
Z^2=-I. Any invariant row space containing row(C) must contain row(CZ),
which proves minimality and the rank formula. A Z-invariant real space is a
complex vector space under Z, hence has even real dimension. For arbitrary
surjective coordinates C, the induced metric on its image is (CC^T)^-1;
using a numerical identity metric there without calibration changes the target.

**W33: a cut can preserve squared size accounting while losing the pairing.**
Take m=1, Y=I_2, C=(1,0). Then

    Q(Y) = [[1,i],[-i,1]],
    Q(PY) = diag(1,0),         Q(Q0Y)=diag(0,1).

The missing imaginary off-diagonal entries are the seam term. In fact
det(Q(Y)-Q(PY))=-1: this real-coordinate projection is not a monotone quantum
information channel. One extra scalar channel repairs it. It supplies no
counterexample to the quantum data-processing inequality, whose channel and
metric hypotheses have not been installed here.

**W34: Z closure is weaker than full dynamical closure.** In m=2, keeping
both coordinates of channel 0 commutes with Z and gives a valid tensor ledger.
It does not descend mixing with channel 1. U28/U29 require the entire connected
component for the full action catalogue. Tensor repair and dynamical repair
answer different target questions.

## U32. Signed response information ledger and target-faithful descent

Now a different cut is used: let finite preparations x_a have independently
declared positive weights p_a summing to one, signed responses r_a=Lx_a, and
observation labels c_a (for instance c_a=C_in x_a). The weights are classical
ensemble bookkeeping supplied for this calculation, not probabilities derived
from the uncut primitive or a Born rule.

Let Sigma be the covariance of r, mu_c its conditional means, p_c its block
weights, and Sigma_c its conditional covariances. Then

    Sigma = Cov_c(mu_c) + sum_c p_c Sigma_c.

Both terms are positive semidefinite. The discarded term is zero exactly when
the signed response is constant on each observed fibre of the finite ensemble.
For nested observation partitions fine -> coarse, the discarded covariance
obeys the exact additive tower ledger

    Delta_coarse = Delta_fine
                 + Cov_within_coarse(E[r | fine]).

All of this is CID-1 total covariance applied to a newly constructed native
response. It is not the real-coordinate tensor ledger of U31, nor the QTH-1
quantum measurement gap.

**Proof.** Write r-mu=(r-mu_c)+(mu_c-mu). The conditional mean of the first
summand is zero, so both cross expectations vanish. This proves the identity.
Each covariance is a sum of positive weighted outer products. If the total
discard vanishes, its trace is a sum of nonnegative p_a||r_a-mu_c||^2;
strict positivity of the weights forces every term to vanish. Apply the same
decomposition within each coarse block to get the tower identity. Notice that
the vanishing of conditional cross expectations is a statistical operation;
it is not the orthogonality/complex-invariance condition of U31.

On the whole real input carrier, exact recovery of Lx through a linear cut
C_in is equivalent to ker(C_in) subset ker(L). The same condition is necessary
and sufficient to recover the scalar ||Lx||^2 for all x. The minimum extra
linear channels are rank([C_in;L])-rank(C_in).

**Proof.** The vector claim is the inherited factorization theorem. If the
squared norm factors through C_in, compare n in its kernel with 0: ||Ln||^2=0,
so Ln=0. Conversely vector recovery recovers its norm. Stacking the rows gives
the smallest target-containing row space, proving the rank formula. A finite
ensemble test alone does not prove this whole-carrier condition.

If C_in also descends U and V, then

    C_in L = (U_bar V_bar - V_bar U_bar) C_in.

The same identity holds for group-loop displacement with the corresponding
quotient products. Thus the **observable signed order target**, not a renamed
positive norm, survives a lawful quotient. This finite square is the algebraic
counterpart of the thermodynamic transport square. To turn it into smooth
connection curvature requires the paper's bundle, connection and small-loop
hypotheses. To turn that into metric curvature also requires U20's coframe and
torsion conditions. A physical gravitational interpretation remains open.

## 3. Explicit responses and rejection controls

For a mixing link (a,b) with nonzero parameter s and local rotation on a with
nonzero parameter t, the response has the exact scale law

    L^T L = kappa P_ab,
    kappa = 16 t^2 s^2 / ((1+t^2)(1+s^2)^2) > 0,

where P_ab is the projector onto those two real-pair channels. In their
four-dimensional block, write U=[[cI,-bI],[bI,cI]] and V=diag(T,I).
Direct multiplication gives L=[[0,b(T-I)],[b(T-I),0]]. Since T is the
local R Cayley rotation, (T-I)^T(T-I)=4t^2/(1+t^2) I and
b=2s/(1+s^2), proving the formula. Outside that block L is zero. This proves
rank(L)=4 for all finite nonzero parameters, not just the tested grid.

**W35.** In U27's m=2 example with both Cayley parameters 1/2, r(e_4)'s first
coordinate is -16/25. The four-column map L has full real rank four. Hence
every exact whole-carrier predictor of its signed response or squared norm
needs four independent real input channels. This requirement is not evaded
by replacing the response with its positive tensor.

**W36.** In m=3 with only interaction link (0,1), L ignores the last channel.
The 6-to-4 whole-component quotient preserves the signed response and tensor
on the retained response carrier. It is a proper information-losing example,
not a faithful six-dimensional relabelling.

**W37.** Preparations x and -x with equal weights have mean response zero,
but covariance (Lx)(Lx)^T, nonzero when Lx is nonzero. A zero average is not
absence of response. Separately, L and -L have the same complete Gram tensor
Q(L), while opposite signed order readings distinguish them. The tensor does
not replace oriented calibration or the raw response record.

**W38.** For the declared QTH-1 comparison, at rho=I/2 the SLDs sigma_x and
sigma_y have nonzero commutator 2i sigma_z but zero state average. The original
nonzero witness remains valid; its converse cannot be generalized. The QTH
source, regression and certificate are corrected in this PR.

## 4. What this achieves and what it leaves to select

The retained native interaction now has a common positive/oriented response
tensor, a precise seam term for incompatible real cuts, a minimal paired
repair, and a statistical information ledger tied to an explicit signed target.
Inherited Gram identities, total covariance and observer-rank arguments are
credited. There is no claim of priority for those mathematical methods.

This supplies a candidate response object for a gravity-before-curvature
programme; it does not establish that the object is gravitational. No Born
law, quantum metric, Lorentzian metric, equivalence principle or force law is
derived. The next discriminating task is to specify a preparation/readout law
and a physical response criterion that separates this candidate from rival
native transport families. The smooth thermo adapter can then be tested on
that same target without inserting its desired geometric conclusion as input.
