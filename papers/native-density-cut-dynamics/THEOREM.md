# DCD-1 — Density and coherent cuts generate response transport

16 October 2026. Frozen parent: `59c68caac773a6c086c16d3c5c565cd7275cb266`.

## Source contract and result

Keep QRT's selected LMD sector H=L²(mu_h), 0<h<1, its positive density
measure M, X_i=M(w_i), and its two disjoint smooth nonnegative windows.
No representation, orientation sector, native engine or original current
is replaced. Write T=M(S1) for total density mass; T is not a clock.

For a normalized exponential-filter state psi_theta, write

    Y_i=X_i-m_i, E_i=Y_i psi_theta, H=E*E,
    Pi=|psi_theta><psi_theta|, Q=E H^(-1)E*, R=I-Pi-Q.

Pi and Q are mutually orthogonal projections. Along a prescribed smooth
control path theta(t), define the self-adjoint multiplication operator

    D(t)=dot(theta)(t).X

by joint spectral calculus. D is a density combination, not the inherited
Cameron--Martin momentum D(k). We retain the latter notation only with
its argument k. The path parameter is a control label, not an independently
derived physical time.

The main exact identity is

    boxed: B=[Q',Q]=[Pi+Q/2,D],     K'=BK.                    (1)

This supplies a density-plus-coherent-cut implementation of QRT transport.
The coherent phase controls are explicitly declared. They do not follow
merely from the existence of a spectral measurement. We prove that the
old configuration-preserving density/current flow class cannot implement
the same small-loop state and response action. Smooth bounded-mass
preparations retain QRT's nonzero curvature and close its current-domain
gap. The final theorem identifies precisely why density-response geometry
alone does not reconstruct the original material metric.

## DCD-T1. Smooth preparations with nonzero curvature and current domains

The joint spectral law of (X_1,X_2) under mu_h is equivalent to Lebesgue
measure on the positive quadrant. Consequently every nonzero normal
vector state has an absolutely continuous joint X law. If its second
moments are finite, its covariance is strictly positive. No nonzero
polynomial in X can vanish almost surely in such a state.

There exist smooth bounded-mass preparations psi_* with each of the
three QRT cells carrying probability 1/3, such that

    1/8 < kappa(0) < 3/8,                                  (2)
    psi_theta in Dom R_(a_h)[v] intersect Dom D(k)           (3)

for every smooth real v,k and every real theta. All theta derivatives
have these domains and depend continuously on theta in each fixed
current graph norm. Every density monomial domain is also included.

**Proof of the spectral statement.** QRT's independent Cameron--Martin
rescalings imply that the law nu of (log X_1,log X_2) is equivalent to
each translate of itself. For every Borel A, Tonelli gives

    integral_R² nu(A-t) dt=Leb(A).

If Leb(A)=0, some translate has zero nu measure, and translation
equivalence gives nu(A)=0. Conversely nu(A)=0 implies nu(A-t)=0 for
every fixed t and hence Leb(A)=0. No simultaneous samplewise action for
all t is required. Exponentiation proves the claim for X. A vector-state
law is absolutely continuous with respect to this law. Nonzero polynomial
zero sets have Lebesgue measure zero, by the elementary one-variable root
bound and induction using Fubini. In particular a nonzero affine variance
cannot vanish.

**Preparation.** Use delta=1/10000 and the three cells B_j centered at
(1,1), (2,1), (1,2), with coordinate half-width delta. Choose nonnegative
smooth nonzero bumps chi_j compactly supported inside these cells. Their
positive sets have positive spectral probability. Since T is finite
almost surely, there is a common finite L such that every bump is nonzero
on a positive-probability set with T<=L. Choose a smooth function chi_0
on the real line equal to one on [0,L], zero above 2L, and between zero
and one on the positive axis. Set

    f_j=chi_0(T)chi_j(X), p_j=E_mu |f_j|²>0,
    psi_*=(1/sqrt(3)) sum_j f_j/sqrt(p_j),
    psi_theta=Z(theta)^(-1/2)exp(-theta.X/2)psi_*.           (4)

The supports are disjoint, so psi_* is normalized and gives each cell
weight 1/3. The X values remain unbinned. The covariance and third-moment
error bounds of QRT-T4 depend only on these cell widths and weights,
not on the within-cell law. They apply unchanged and prove (2).
All real theta are allowed; the support has T<=2L and bounded X.
The spectral statement proves H(theta)>0 everywhere, with uniform
positive lower eigenvalues on compact theta patches.

**Current domains.** Each vector in (4), each E_i and each theta
derivative is a bounded smooth function F(T,X) times the vacuum, with
a common bounded-mass support. For the corrected material current
R_(a_h)[v], a_h=1+h/2, LMD proves that T is invariant and that conjugation
transports window observables by w_i composed with g_v(t). On T<=2L,
the first two parameter derivatives of these window observables are
bounded by 2L times the corresponding uniform window derivative bounds.
Smooth composition therefore makes the conjugated multiplier M_F
differentiable in operator norm, with a bounded derivative.

For a Cameron--Martin profile k, the conjugated density involves
M(e^(t k)w_i) and T_t=M(e^(t k)). With K_k=||k||_infinity,

    exp(-|t|K_k)T <= T_t <= exp(|t|K_k)T.

For |t|<=t_0 all transported F vanish when T>2L exp(t_0 K_k).
On the remaining common bounded-mass region, the first two derivatives
M(k^r e^(t k)w_i) and M(k^r e^(t k)), r=1,2, are uniformly bounded.
This proves norm differentiability even though the mass cutoff moves.

If a multiplier has this norm-differentiable conjugation, differentiating
its product with a vector in a Stone generator's domain proves domain
preservation and the bounded commutator formula. LMD already places the
vacuum in Dom R_(a_h)[v] and Dom D(k). Applying the product rule proves
(3). The same uniform estimates with theta derivatives give graph-norm
continuity and local uniform graph bounds. This uses actual Stone
domains, not a formal replacement of sharp indicators by smooth notation.

## DCD-T2. The source generator and its uniquely required cut weights

The commutator in (1), initially on Dom D, extends to a bounded
anti-self-adjoint finite-rank operator. It is exactly QRT's B. Moreover

    B psi_theta=psi_theta',
    K(t)psi_0=psi_theta(t),
    K(t)E(0)=E(t)U_TVSP(t).                                (5)

For nonzero velocity dot(theta), the two-cut ansatz
[a Pi+b Q,D] equals B precisely at a=1,b=1/2. A scalar multiple of I
can of course be added without changing a commutator. This uniqueness
is within the stated ansatz, not among all possible control Hamiltonians.

Put v=dot(theta) and q=E v. Then

    ||B||=||Q'||=(1/2)||(I-Q)DQ||,
    ||B|| >= (1/2)sqrt(v^T H v).                            (6)

On the three-cell support, the density combination D has range width
at most (1+2delta)(|v_1|+|v_2|), so

    ||i h B|| <= h(1+2delta)(|v_1|+|v_2|)/4.                (7)

The self-adjoint control generator i h B has domain all of H. Equation
(7) is a bound in the declared canonical units and control parameter;
it does not assign a laboratory energy or clock rate.

**Proof.** Write d=<D>. QRT's differentiation formulas give

    psi'=-(D-d)psi/2=-q/2,
    E'=psi(Hv)^T-(D-d)E/2,
    RD Pi=Pi DR=0.

Differentiating Q and using Pi,Q,R as three orthogonal blocks yields

    Q'=(Pi DQ+QD Pi-RDQ-QDR)/2,
    [Q',Q]=(Pi DQ-QD Pi+QDR-RDQ)/2=[Pi+Q/2,D].              (8)

The products D Pi and DQ are bounded finite rank, because their ranges
are spanned by vectors in Dom D. Their adjoints define the bounded
extensions of Pi D and QD. This proves the domain meaning of (8).
Acting on psi gives B psi=-q/2, and uniqueness of the bounded evolution
proves the first transport in (5). QRT proves the response intertwiner.

The Pi-to-Q block of [a Pi+b Q,D] has coefficient a-b and its Q-to-R
block has coefficient b. The first is nonzero because q=Ev!=0. The
map RD restricted to Ran Q is injective: otherwise
(v.X)(c.(X-m)) would equal an affine polynomial almost surely for
some c!=0. Its quadratic leading part is nonzero, contradicting T1.
Thus both block coefficients are forced, giving b=1/2 and a=1.
The off-diagonal block between Q and I-Q consequently has rank two,
and B has rank four for v!=0 in the continuum preparation.

In a saturated three-atom comparison R=0, so only a-b=1/2 is detected.
In particular the incorrect W=Pi/2 transports the base state correctly
but misses the higher response block in the continuum. Generic four-
and five-atom certificates expose this missing block.

The off-diagonal block of B is (Pi-R)DQ/2. Pi-R is isometric on
Q's orthogonal complement, proving the norm equality in (6).
Act on psi for its lower bound. Subtract the midpoint of D's spectral
range on the fixed support; a scalar does not affect the commutator.
The centered operator has norm at most half the range width, proving
(7). No extrapolation from a finite atom count establishes these bounds.

## DCD-T3. Preservation of the original current domains

For the preparations in T1, K(t) and K(t)^(-1) preserve each of

    Dom R_(a_h)[v],   Dom D(k),
    Dom(X_1^alpha_1 X_2^alpha_2),                           (9)

for every fixed smooth real v,k and nonnegative integer multiindex
alpha. They also preserve any intersection of these domains.

**Proof.** Fix one of the closed operators A in (9). The explicit
finite-rank expression for B uses psi,E and their first path derivatives.
T1 puts all of them in Dom A, continuously in its graph norm and with
uniform bounds on compact smooth path segments. Hence AB(t) extends
to a uniformly bounded, norm-continuous map H->H. In

    K(t,r)f=f+integral_r^t B(s)K(s,r)f ds,

the integral is in Dom A by graph-norm approximation and closedness.
Thus K maps Dom A into itself. The backward two-parameter evolution
has the same left-generator argument and gives equality of domains.
Apply this to every specified A. No norm convergence of pulses alone
would prove this unbounded-domain result; it follows from these graph
bounds. We do not identify K with a material-current flow.

## DCD-T4. Density and coherent-cut pulse synthesis with an error bound

Let W=Pi+Q/2. Assume coherent phase control of Pi and Q, together with
the original density spectral phases. The cut phase is explicitly

    exp(-i delta_p W)
      =I+(exp(-i delta_p)-1)Pi+(exp(-i delta_p/2)-1)Q.       (10)

It is a coherent unitary, not a selective or nonselective measurement.
Freeze the controls at one point and use, with rightmost action first,

    G_delta_p=exp(-i delta_p W) exp(i delta_p D)
              exp(i delta_p W) exp(-i delta_p D).           (11)

Let P_B be the joint spectral projection of the fixed three-cell union.
It reduces D and W=P_B W P_B. Put D_B=P_B D P_B, extended by zero.
The product (11) is exactly unchanged if D is replaced by D_B:
outside P_B the cut phases are identity and the density phases cancel.
Thus the implementation uses the original full density phases, while
the proof uses bounded operators.

For w=||W||=1, d=||D_B|| and b=||[W,D_B]||,

    ||G_delta_p-exp(delta_p² B)||
       <= |delta_p|³(w+d)b + delta_p^4 b²/2.                (12)

Consequently for tau>=0,

    ||G_sqrt(tau/n)^n-exp(tau B)||
       <= (w+d)b tau^(3/2)/sqrt(n)+b²tau²/(2n).             (13)

Negative transport uses the inverse sequence. These are operator-norm
estimates on the continuum Hilbert space, not finite-matrix evidence
for a limit. Scalar recentering of D can reduce d without changing
either the full four-pulse product or B.

For a varying C² control path on [0,t_f], choose m intervals of width
Delta_t=t_f/m, freeze W,D at their left endpoints, and take n cycles
per interval with delta_p=sqrt(Delta_t/n). Later intervals act on the
left. If w,d,b bound the frozen norms and B is Lipschitz with constant
L_B, the resulting product Khat obeys

    ||K(t_f)-Khat|| <= L_B t_f Delta_t/2
       +(w+d)b t_f sqrt(Delta_t/n)
       +(b²/2)t_f Delta_t/n.                              (14)

Thus even n=1 converges as the mesh tends to zero. For piecewise C²
paths, use a partition containing every breakpoint, put Delta_t equal
to its maximum interval length, and use delta_p,j=sqrt(ell_j/n) on
an interval of length ell_j. Bound (14) holds with uniform segmentwise
bounds and this maximum mesh. For merely continuous B,
replace the first term by t_f omega_B(Delta_t), its uniform modulus
of continuity bound.

**Proof.** Put A=-iW, C=iD_B, so [A,C]=B, and
C_delta_p=exp(delta_p A)C exp(-delta_p A). For delta_p>=0 the exact
unitary Duhamel identity gives

    G_delta_p-I=integral_0^delta_p integral_0^delta_p
       exp(s C_delta_p) exp(r A) B exp(-r A) exp(-s C) dr ds.

The integrand differs from B by at most 2sdb+2rwb. Integrating gives
||G-I-delta_p²B||<=(w+d)b delta_p³. The integral Taylor remainder
for the anti-self-adjoint B bounds ||exp(delta_p²B)-I-delta_p²B||
by delta_p^4 b²/2. Negative delta_p follows by changing A,C signs.
Unitary telescoping proves (13). Unitary Duhamel comparison between
the moving generator and its frozen version contributes at most
L_B t_f Delta_t/2; summing the cycle errors gives (14).

Each pulse preserves the density monomial domains: density phases
commute with their spectral multipliers, and cut phases differ from
identity by vectors supported in P_B. The current-domain claim for
the limiting K is independently established in T3. An arbitrary full
density phase is not here claimed to preserve every current domain.
The construction assumes the coherent cut phases; it does not derive
a laboratory actuator or its strength and duration calibration.

## DCD-T5. Why the old density/current flow class cannot supply the loop

Let A_conf=L-infinity(mu_h) be the full configuration multiplication
algebra on H. A unitary normalizes it when V A_conf V*=A_conf. Actual
LMD density phase unitaries, Cameron--Martin displacement unitaries,
and material-flow half-density unitaries normalize A_conf. Every finite
product does so, as does any strong limit of these products whose limit
is unitary.

Every nonidentity QRT lift K(t) is excluded from this class. More
strongly, for T1's actual smooth preparations, sufficiently small
nonzero curvature loops admit no normalizer V with both

    V psi_0=psi_0,           V E(0)=E(0)U_loop.              (15)

Thus this obstruction persists even if V is allowed to differ from K
away from the base and response states. It excludes the stated flow
class, not arbitrary higher-order momentum Hamiltonians, associative
operator sums, ancillary controls or all possible quantum dynamics.

**Full-space proof.** The Gaussian configuration probability space is
nonatomic. A compact multiplication operator on such a space is zero:
if its modulus exceeds epsilon on a positive-measure set, an infinite
orthonormal sequence supported in that set contradicts compactness.
Also this full multiplication algebra is maximal abelian. For completeness,
on a probability space a commuting bounded operator is multiplication
by its action on 1, first on bounded functions and then by density;
boundedness makes that multiplier essentially bounded.

If V-I is compact and V normalizes A_conf, then
V M_f V*-M_f is both compact and multiplication for every bounded f.
It is zero. Thus V commutes with A_conf and is itself multiplication;
compactness of V-I then forces V=I. On the other hand,

    K(t)-I=integral_0^t B(s)K(s) ds

is compact, because the integrand is norm-continuous and finite rank.
Hence nonidentity K cannot normalize A_conf. Normalizer closure under
unitary strong limits follows because both V_n and V_n* converge strongly,
and A_conf is strongly closed. LMD's half-density form proves the stated
normalizer property of its original flows without any formal sum-domain
assumption.

**Response-action proof.** T1's cell bound gives a nonzero third central
moment: its u_111 differs from -2/27 by at most 6delta. Whiten the centered
X variables using the oriented response frame; their covariance becomes
I and their symmetric third moment tensor remains nonzero. If (15) held,
the normalizer's configuration automorphism would rotate these variables
by the oriented loop angle Phi while preserving their probability law.
To justify this with possibly unbounded variables, use their bounded
spectral truncations and the identity V psi_0=psi_0; normalizing conjugation
acts on the affiliated multiplication operators and preserves the vector
state. The response equality then identifies their transformed variables
almost everywhere in that state.

The symmetric cubic tensor would be fixed by that rotation. Over the
complex numbers the action on homogeneous cubics has eigenvalues
exp(i Phi), exp(-i Phi), exp(3i Phi), exp(-3i Phi). Equivalently,

    det(I-Sym³ Rot(Phi))
       =4(1-cos(Phi))²(2cos(Phi)+1)².                      (16)

For 0<|Phi|<2pi/3 there is no nonzero fixed cubic tensor. This contradicts
the nonzero third moment. T1's positive curvature and QRT's loop phase
formula supply arbitrarily small loops in that range. This proves (15)'s
impossibility for the actual smooth continuum states, not merely for a
finite categorical fixture.

## DCD-T6. What density-response geometry cannot determine

Multiply a state by the configuration phase

    V_d(s)=exp[-i<s,d_*>/h],                               (17)

with smooth real d_*. Every density spectral law, filter normalizer,
covariance H, selected TVSP curvature and loop phase is unchanged.
The prepared base and response vectors and their coherent cuts are
conjugated by this fixed unitary. On the corresponding Stone domains,

    V_d* D(k) V_d = D(k)+integral d_* k.                    (18)

Thus conjugate momentum information can change while the entire
density-response metric remains identical.

The original classical material metric has exactly this missing data.
In the inherited inverse chart take homogeneous j=j_0>0 and d=d_0.
With vartheta=Y^(-1)(d_0/j_0), m=j_0/a(vartheta) is constant and
P=integral m=m on the unit circle. Therefore

    r=vartheta,    ds_material²=exp(2vartheta) dz².          (19)

Varying d_0 leaves the complete homogeneous density j fixed but changes
the material circumference and metric. Consequently no density-only
response data can determine that material metric on the full inherited
classical chart. LMD has not constructed a pointwise quantum r; (19)
is a classical nonidentifiability theorem, not a new quantum composite.

**Proof.** (17) commutes with all density multiplication operators, so
the equalities of spectral laws and filter metrics are exact. Conjugating
the Cameron--Martin group, or differentiating its phase on its Stone
domain, gives (18) with the displayed positive sign. Substitute the
homogeneous profiles into the inherited inverse map for (19). Since
Y is strictly increasing, different d_0/j_0 give different r. A metric
depending on both s and d cannot be inferred from a family insensitive
to d's phase information.

## Closure and next obligation

The canonical QRT transport now has an explicit density/coherent-cut
commutator, a convergent pulse implementation with a uniform error bound,
and preparations preserving the original current domains. The need to
leave the configuration-normalizing flow class is proved even at the
level of small-loop state and response action. These conclusions are
stronger than assigning a name to i h B.

The extra coherent phase controls remain an operational assumption.
The original configuration-normalizing flows cannot replace them in the
declared task. Material geometry needs conjugate information in addition
to density response; this is now an explicit nonidentifiability result.
The next physical construction must specify a coherent source interaction
or a conjugate-sensitive response family and compare its actual domains
and geometric readouts. No universal spacetime equation, mass/energy law
or 3+1 quantum-gravity closure follows from this packet.
