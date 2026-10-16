# HCM-1 — Exact observer memory in the published YM-54 heat protocol

15 October 2026. Application of MZC/CRC to the existing free Yang–Mills benchmark.

## 1. Source and scope

YM54's published thermo fixture fixes H=[[2,1],[1,5]], response derivatives
H_1=[[0,1],[1,0]], H_2=[[1,0],[0,0]], and factor B=[[2,1],[0,3]], B*B=2H.
Its trace-free factor-conjugated response words give protocol vectors
v1=(1/9,1/3,0), v2=(2/9,-1/6,0), with equal weights. Thus

    C=(v1 v1^T+v2 v2^T)/2=diag(5/162,5/72,0).

The certificate reconstructs these values with the unchanged native engine
and compares them to YM54_RESULT.json. This tensor and its heat generator are
inherited source data, not a newly fitted dynamics. The mixed observer below
is a declared choice; no physical observer selection is claimed.

The target is one exact three-dimensional invariant coefficient block of the
free single-site heat process. It is not YM49's entire instantaneous interacting
row or a solution of the infinite-chain closure gate. The general compression
identity and memory-elimination method already appear in YM49 and RKF N09;
they are reused, not claimed as new general theorems.

## 2. Exact block inside the native coefficient carrier

On YM50's unit quaternion x=(x0,x1,x2,x3), take one column of Ad_x:

    r1=x0^2+x1^2-x2^2-x3^2,
    r2=2(x1 x2+x0 x3),
    r3=2(x1 x3-x0 x2).

These are even, centered coefficient functions. The reference pairing gives
Phi(r_i r_j)=delta_ij/3, so their common Gram scale cancels from projections.
For D_a f(x)=d/du f(exp(-u e_a/2)x)|_0, differentiation of Ad_x shows
D_a r=-J_a r, with J_a the three-dimensional cross-product matrix.
Since J_a J_b=e_b e_a^T-delta_ab I, symmetric C yields

    L_C r = Lambda r,     Lambda=tr(C)I-C.

This is the existing spin-one/adjoint block, derived directly rather than
assumed from a spectral fit. The native replay also differentiates the explicit
four-variable quadratic forms. For the actual fixture,

    Lambda=diag(lambda1,lambda2,lambda3)
          =diag(5/72,5/162,65/648).

For f_z=sum z_i r_i the coefficient evolution is z(t)=exp(-t Lambda)z(0).
This is heat evolution, not unitary dynamics or physical Lorentzian time.
Its diagonal spectral sectors are invariant and its derivative energy is
nonnegative. That does not make every mixed observation channel autonomous.

## 3. A source-specific failure of one-channel closure

Choose orthonormal coefficient directions

    u=(3/5,4/5,0),   v=(-4/5,3/5,0),   w=(0,0,1).

The fixed readout is q(t)=u^T z(t). The projector P=u u^T is an observation
projector; it is not an allowed spectral-sector projector of Lambda.
Native rowClosure(u^T,[-Lambda]) has rank 2. Its kernel is span{w}; retaining
u^T and v^T is the minimal autonomous linear record for this chosen channel.
Thus a three-mode source does not require all three modes for this target,
but one mode is insufficient. This is CRC's target-relative minimality.

Start from z(0)=v. The cut reading is initially zero, but

    q(t)=(12/25)(exp(-5t/162)-exp(-5t/72)) > 0  for t>0,
    q'(0)=1/54.

The common fixed/co-moving zero carrier is span{w}, not all of ker u^T.
An algebraic moving readout u^T exp(t Lambda) would keep the initial reading
zero, but is a backward-heat compensation with norm growing exponentially.
It is not a unitary frame change or a stable inverse on the whole completion.

For a unit initial source in direction u, define its compressed propagator

    a(t)=u^T exp(-t Lambda)u
        =(9/25)exp(-5t/72)+(16/25)exp(-5t/162).

The YM49-style two-time defect is exactly

    a(2t)-a(t)^2
       =(144/625)(exp(-5t/72)-exp(-5t/162))^2 > 0  for t>0.

Proof: substitute the two exponentials and expand their weighted variance.
The native symbolic certificate proves this polynomial identity in two
independent scalar variables. Since lambda1 differs from lambda2, the defect
is strictly positive for t>0. Thus a(t) cannot be a scalar semigroup. Moreover

    lim_(t->0) [a(2t)-a(t)^2]/t^2 = 1/2916.

For 0<t<72/5, the integral of exp(-s t) between lambda2 and lambda1 gives

    (t^2/2916)(1-5t/72)^2 <= a(2t)-a(t)^2 <= t^2/2916.

At one declared heat tick the explicit lower bound is 4489/15116544>0,
and the upper bound is 1/2916. No floating-point exponentiation or data fitting
is needed to certify nonclosure. The exact Taylor enclosures in the executable
are additional checks of these bounds, not their general proof.

## 4. Exact retained-memory repair and a quantitative tail

In the basis (u,v), the positive heat matrix is

    Hobs=[[29/648,-1/54],[-1/54,1/18]].

Set h(t)=v^T z(t). Direct matrix multiplication gives the closed two-record law

    q'=-(29/648)q+(1/54)h,
    h'=(1/54)q-(1/18)h.

Eliminate h by variation of constants:

    q'(t)=-(29/648)q(t)+(exp(-t/18)/54)h(0)
          + integral_0^t [exp(-(t-s)/18)/2916] q(s) ds.

Thus the exact continuous memory kernel of this source/observer pair is

    K(t)=exp(-t/18)/2916.

This is a new evaluated instance of the already established memory mechanism,
not a new general memory formalism. Keeping one hidden channel closes the law
exactly. Dropping it creates a wrong Markovian model even though the full heat
process is positive and has a relaxation gap.

If only memory ages below T are retained, then for t>=T the omitted convolution
is bounded by

    exp(-T/18)/162 * sup_(0<=s<=t-T)|q(s)|.

The omitted hidden-initial-data term is bounded separately by
exp(-t/18)|h(0)|/54. Integrating the positive kernel proves the tail bound.
These constants use the declared heat clock and selected observer. They are
not physical masses, energy units or a new interacting-chain gap.

## 5. Interpretation and next boundary

This application has moved from selected abstract matrix controls to the
published YM54 response-derived free generator. It supplies: an exact
nonclosure value, its shortest leakage witness, the minimal retained record,
the continuous memory kernel and an explicit truncation-error budget.

It also separates three questions: the heat tensor is positive; its spectral
sectors are invariant; an observer mixing distinct rates may still need memory.
Calling a mixed observer an admissible sector does not make it spectral or
remove its coupling. At isotropy or along a single eigenchannel, the defect
vanishes, as the negative controls require.

The next interacting application must preserve the actual ground-source
pairing, bounded evolution, observation embedding and spatial/tail limits.
This free block calculation does not evaluate YM49's full interacting defect,
establish continuous-time row-Markov closure or solve physical 4D Yang–Mills.
No claim about RH source positivity or completed quantum gravity follows.
