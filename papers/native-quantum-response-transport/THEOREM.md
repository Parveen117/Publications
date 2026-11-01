# QRT-1 — Quantum density responses realize TVSP transport

16 October 2026. Parent: `9a0621637b2cc12fb57eb6fd36e56daef8cf5f34`.

## Source contract and result

Work in one of LMD's existing orientation sectors
H_epsilon^log=L²(mu_h), 0<h<1. Its positive random measure M, physical
density operators J_h(w)=epsilon M(w), Cameron--Martin displacement law
and spectral domains are retained. The opposite orientation is a separate
copy; the two signs are not the two response channels below.

Choose two nonzero smooth nonnegative calibration windows w_1,w_2 with
disjoint compact supports in separated material arcs, and set

    X_i=epsilon J_h(w_i)=M(w_i),    i=1,2.                       (1)

These are positive, strongly commuting self-adjoint multiplication
operators on the existing physical Hilbert space. The windows, state
preparation and observer-control path are declared inputs. No extra
two-state probe, new configuration measure or gravitational field is added.

The result supplies the intertwiner absent in TET: actual density-response
vectors in this Hilbert space induce exactly A_i=(1/2)H^(-1)partial_i H.
Their parallel transport is the restriction of a full-space unitary with
an everywhere-defined bounded instantaneous generator. A prepared state
using the actual unbinned X_i has rigorously nonzero response curvature.
The geometry is on observer-control space. Identification with the old
material metric exp(2r) dz² or spacetime curvature is not inferred.

The covariance Hessian of an exponential family and projection parallel
transport are established constructions. Their source-specific realization,
domain statements and curvature bounds are proved directly below. The
source map credits those tools; no adiabatic approximation is assumed.

## QRT-T1. Actual source states and a strictly positive response Hessian

For theta=(theta_1,theta_2) in the positive quadrant define

    Z(theta)=<1,exp(-theta_1 X_1-theta_2 X_2)1>,
    psi_theta=Z(theta)^(-1/2) exp(-theta.X/2)1,
    u(theta)=log Z(theta),
    m_i=<psi_theta,X_i psi_theta>,
    H_ij=partial_i partial_j u=<(X_i-m_i)(X_j-m_j)>_theta.       (2)

Then Z>0, the normalized vectors are smooth in Hilbert norm and in every
joint polynomial domain of X_1,X_2, and H is strictly positive definite.
The joint spectral law of (X_1,X_2) under mu_h has full support relative
to (0,infinity)². On compact positive theta patches all derivatives and
polynomial graph norms needed below are uniformly bounded.

**Proof.** LMD full support of M gives 0<X_i<infinity almost surely.
Bounded positive spectral filtering gives 0<Z<=1. For every multiindex
alpha and c>0, x^alpha exp(-c(x_1+x_2)) is bounded on the positive
quadrant. Differentiation under the integral and in L² therefore holds
to every finite order on compact positive patches. It gives

    partial_i u=-m_i,
    partial_i m_j=-H_ij,
    partial_i H_jk=-<Y_iY_jY_k>_theta,  Y_i=X_i-m_i,
    partial_i psi_theta=-(1/2)Y_i psi_theta.                     (3)

The last identity also follows by differentiating the normalization.

Choose smooth Cameron--Martin profiles k_i equal to one on supp(w_i)
and zero on the other window support. LMD's displacement law gives
X_i(s+t k_i)=exp(t)X_i(s), with the other X unchanged. Quasi-invariance
implies equivalence of the original joint law and each independently
rescaled law. Thus its support is invariant under all positive independent
rescalings. It contains a point with both coordinates positive, so it
contains every such point. In particular every nonempty open rectangle
in the quadrant has positive probability.

If a^T H a=0, then a.X is almost surely constant under the tilted law,
which is equivalent to mu_h. Independent rescaling would preserve that
affine relation. Subtracting the original and rescaled relations forces
a_i(exp(t)-1)X_i=0, so a_1=a_2=0. This proves strict positivity.
Continuity and compactness give a uniform positive lower eigenvalue and
positive lower bound for Z on any compact patch.

More generally replace 1 in (2) by a normalized base vector psi_*. The
same formulas hold on any open parameter region where the required
exponential moments exist, provided the resulting covariance is positive.
In particular a base vector supported in a bounded joint spectral region
admits every real theta and every polynomial domain. QRT-T4 constructs
such a vector from the actual continuum density, not from an atomic
replacement law.

## QRT-T2. The exact TVSP connection is induced by quantum response vectors

Define the response frame E_theta:C² -> H_epsilon^log by its two columns

    E_i(theta)=Y_i psi_theta.

Then

    E^*E=H,
    (E^*partial_i E)_kj=(1/2)partial_i H_kj,
    A_i=H^(-1)E^*partial_i E=(1/2)H^(-1)partial_i H.              (4)

Thus the orthogonal projection of Hilbert-space differentiation onto
the actual response plane induces exactly NT/TET's selected connection.
Its curvature and oriented scalar marker are

    F_12=-[H^(-1)partial_1H,H^(-1)partial_2H]/4,
    Delta=det H,                 kappa=(H F_12)_12/Delta.        (5)

All coefficients are determined by quantum expectation values of the actual density channels.
First derivatives of H are third central moments, so curvature uses these
moments rather than just the means or covariance eigenvalues.

**Proof.** The Gram identity is the covariance definition. By (3),

    partial_i E_j=H_ij psi_theta-(1/2)Y_iY_j psi_theta.

Pairing this with E_k annihilates the first term and gives minus one half
of the third central moment. Equation (3) proves (4). The source coordinates
commute and H is a smooth Hessian, so NT-3/4 or direct inverse differentiation
gives (5). No arbitrary matrix has been inserted into a quantum state.

For explicit comparison with TLM's six normalized responses, translate
theta to local coordinates s,v, choose positive affine S,V on the patch,
and use the response potential

    u_cal(s,v)=u(theta(s,v))+b_T s-b_P v.                        (6)

The constants are declared calibrations. They do not change H, (4) or
(5). Set T=partial_s u_cal, P=-partial_v u_cal and use NT-1's actual
six response formulas. On a compact patch b_T can be chosen so T>0;
b_P>=0 gives P>0 here. The resulting layers are generated by this
quantum response source. Equation (6) is a control-response potential,
not a derivation of Kelvin temperature or physical internal energy.

The raw state derivative metric is <partial_i psi,partial_j psi>=H_ij/4.
We use the response/Fisher normalization H=E^*E throughout. Rescaling
the metric by a constant leaves its connection and loop transport unchanged
but rescales its scalar curvature. In particular numerical curvature
values below belong to H, not to H/4.

## QRT-T3. Full quantum unitary lift and domains

The operator

    Q(theta)=E(theta) H(theta)^(-1) E(theta)^*                   (7)

is an orthogonal rank-two projection in the existing Hilbert space. For
a piecewise smooth path gamma within a regular parameter patch, put

    B(t)=[Q'(t),Q(t)],
    K'(t)=B(t)K(t),                  K(0)=I.

There is a unique unitary K(t) on the full Hilbert space. B(t) is bounded,
finite rank, norm-continuous on each smooth segment and anti-self-adjoint.
The instantaneous self-adjoint generator i h B(t) has domain the entire
Hilbert space. The path parameter is not thereby identified with time.

Let U solve TET's U'=-A_gamma U, U(0)=I. Then the exact intertwiner is

    K(t)Q(0)K(t)^*=Q(t),
    boxed: K(t)E(0)=E(gamma(t))U(t).                             (8)

For an H(0)-orthogonal tangent cut P_0, its quantum realization

    Pi_0=E(0)P_0H(0)^(-1)E(0)^*

is an orthogonal cut on the response subspace, and

    K Pi_0 K^*=E(t)[U P_0 U^(-1)]H(t)^(-1)E(t)^*.               (9)

Rank-one choices therefore realize EGR's positive-cut transport on actual
states of the existing sector, not on an attached formal jet module.

**Proof.** The Gram identity proves Q²=Q=Q^*. Differentiating Q²=Q
gives Q'Q+QQ'=Q' and QQ'Q=0. Consequently B^*=-B and [B,Q]=Q'.
The norm-convergent ordered integral series for a bounded continuous B
constructs K. Differentiating K^*K gives I; the reverse evolution proves
surjectivity. Differentiating K^*QK proves the first equality in (8).

For W=E U, equation (4) gives

    W'=(I-Q)E'U.

Differentiate QE=E to obtain Q'E=(I-Q)E'. Since QQ'Q=0 and QE=E,
BW=Q'EU=(I-Q)E'U. Thus W and K E(0) solve the same initial-value
problem, proving (8). Metric preservation U^*H(t)U=H(0) and (8) prove
(9), including self-adjointness of the displayed cuts.

For the vacuum family of QRT-T1, K and its inverse preserve each domain

    D_alpha={f:X_1^(alpha_1)X_2^(alpha_2) f in H},                (10)

and hence their intersection. Indeed E,E' and Q' are finite sums of
vectors polynomial in X times psi_theta. The intersection is dense,
since it contains vectors supported in bounded joint spectral regions.
On a compact positive patch,
X^alpha B(t) is uniformly bounded and norm-continuous as an operator
into H. In

    K(t,r)f=f+integral_r^t B(s)K(s,r)f ds,

the integral belongs to D_alpha: approximate it in the graph norm,
using the bound on X^alpha B and closedness of X^alpha. Therefore K maps
D_alpha into itself. The backward two-parameter evolution has the same
left-generator argument and proves equality of the transported domains.
This is not an inference from finite matrix unitarity.

For the compact-cell family of QRT-T4 there is a simpler proof. If P_B is
its bounded-region joint spectral projection, then B=P_B B P_B and K is
identity on its complement. Each X_i is bounded on P_B H. Thus every
domain (10) is preserved for all compact paths in its real control plane.

These statements do not identify K with LMD's material-flow group T_g or
prove preservation of every conjugate-current domain Dom R_a[v]. The
cell-indicator preparations in particular are not claimed to be in those
current domains. The unitary is the source-response projection transport;
a dynamical implementation by the original current would need a further
source Hamiltonian and domain comparison.

## QRT-T4. Nonzero curvature for actual unbinned density states

Fix delta=1/10000. Let B_0,B_1,B_2 be disjoint open squares of coordinate
half-width delta centered respectively at

    c_0=(1,1),       c_1=(2,1),       c_2=(1,2).

Their spectral probabilities p_j=mu_h((X_1,X_2) in B_j) are positive by
QRT-T1. The following is an actual normalized state in H_epsilon^log:

    psi_delta=(1/sqrt(3)) sum_(j=0)^2 1_Bj(X)/sqrt(p_j).         (11)

Filter this fixed vector as in (2). All real controls are allowed because
X_i are bounded on its support. At theta=0, its covariance is positive
and its NT curvature obeys the explicit bound

                         1/8 < kappa_delta(0) < 3/8.           (12)

The X_i have not been replaced by their cell labels. The bound is about
the original unbinned density observables in the prepared state (11).
No numerical value for p_j or vacuum density-response curvature is needed.

**Proof.** The state gives each cell probability exactly 1/3. Couple the
random vector X to its cell center C. Then |X_i-C_i|<=delta, and the
centered variables Y=X-E X and Z=C-E C obey |Y_i-Z_i|<=2delta.
At the categorical center law,

    H_c=(1/9)[[2,-1],[-1,2]],       Delta_c=1/27,
    (u_111,u_112,u_122,u_222)_c=(-2,1,1,-2)/27.                  (13)

The sign is from the negative exponential in (2). Since |Z_i|<=2/3,
we have |Y_i|<=1. Product expansion consequently gives

    |H_ij-(H_c)_ij|<=4delta,
    |u_ijk-(u_c)_ijk|<=6delta.                                 (14)

For this delta all covariance entries have absolute value <=1/4 and
all third derivatives absolute value <=1/12. Write a=H_11,b=H_12,c=H_22.
Expansion of the two-by-two determinant yields

    |Delta-1/27|<=4delta,           a>=2/9-4delta>0.

In two dimensions, direct expansion of (5) gives

    D_3=det [[a,b,c],[u_111,u_112,u_122],[u_112,u_122,u_222]],
    kappa=-D_3/(4Delta²).                                       (15)

For (13), D_(3,c)=-1/729. Each of the six products in this determinant
changes by at most

    (4delta)/144+2(1/4)(6delta)/12=5delta/18.

Thus |D_3+1/729|<=5delta/3. The preceding positive determinant lower
bound ensures positive definiteness of H. Substitution into (15) gives

    (1/729-5delta/3) / [4(1/27+4delta)²]
       <= kappa_delta(0)
       <= (1/729+5delta/3) / [4(1/27-4delta)²].                  (16)

At delta=1/10000, the lower endpoint is greater than 1/8 and the upper
less than 3/8, by exact rational arithmetic. This proves (12) for every
h in LMD's range, without an atomic approximation claim for the state.
Continuity supplies a neighborhood of zero with positive curvature.

Every finite real tilt is mutually absolutely continuous with the
theta=0 spectral law of this compact-cell state. A zero affine variance
at any theta would therefore give zero variance at zero. Thus H remains
strictly positive throughout the real control plane, with a uniform
positive lower eigenvalue along each compact path.

The equal-weight center comparison has kappa=1/4 at zero. It supplies the algebraic
comparison and certificate fixtures, not a claim that the true measure
places atoms at those centers. The preparations have not been proved to
be obtainable with a specified laboratory apparatus or success rate.
Their Hilbert-space and spectral-filter existence is what is established.

## QRT-T5. Fixed-observer quantum return probabilities

Choose the positive oriented factor S(theta) with S^*S=H and define
F(theta)=E(theta)S(theta)^(-1). Its two columns are orthonormal actual
quantum response states. For a closed rectangle inside a regular chart,
TET's rank-two phase identity and (8) give

    K_loop F(0)=F(0) exp(R Phi),
    R=[[0,-1],[1,0]],
    Phi=integral_Sigma kappa sqrt(Delta) dtheta_1 dtheta_2.      (17)

Prepare eta_1=F(0)e_1 and measure the fixed rank-one projection
Pi_2=|eta_2><eta_2|, eta_2=F(0)e_2, after this transport. The exact
transition probability is

    <K_loop eta_1,Pi_2 K_loop eta_1>=sin²(Phi).                  (18)

Pi_2 is bounded on the entire physical Hilbert space, and the response
states belong to the density polynomial domains already established.
In QRT-T4's positive-curvature neighborhood, sufficiently small nonzero
rectangles have nonzero probability. For a square of side epsilon,

    probability=kappa(0)² Delta(0) epsilon^4+O(epsilon^5).      (19)

The center comparison law gives coefficient 1/432; the true unbinned
state has its own strictly positive coefficient, bounded using (16).

The probability sin²(Phi) alone loses the orientation sign. A second
fixed detector at 45 degrees in the response plane resolves that sign
for sufficiently small loops. With eta_+=(eta_1+eta_2)/sqrt(2) and
Pi_+=|eta_+><eta_+|,

    p_+=<K_loop eta_1,Pi_+ K_loop eta_1>
       =1/2+(1/2)sin(2Phi),
    p_+-1/2=kappa(0)sqrt(Delta(0))epsilon²+O(epsilon³).          (20)

For |Phi|<pi/4 its displacement from 1/2 has the sign of Phi. The
diagonal refers to the orthonormal response frame, not an arbitrary
coordinate diagonal. Neither detector alone determines a globally
unwrapped phase. For exact comparison, cos(Phi)=3/5 and sin(Phi)=+/-4/5
give the same p_2=16/25 but give p_+=49/50 and 1/50 respectively.

**Proof.** F is an isometry from C² into the existing Hilbert space.
The induced connection is the same source connection as TET, so its
oriented orthonormal-frame return is exp(R Phi). Apply (8). The off-diagonal
entry of this real rotation is sin(Phi), proving (18). Smoothness and
the area expansion of (17) prove (19). Nonzero curvature and sufficiently
small area exclude both Phi=0 and nonzero multiples of pi.
The inner product with eta_+ is (cos(Phi)+sin(Phi))/sqrt(2), which
proves (20), including its small-area expansion and orientation claim.

Transporting the detector as K Pi_2 K^* instead preserves its initial
reading. Equation (18) compares to the retained base detector. It is also
not a Berry phase of the one-dimensional base-state family psi_theta:
the latter satisfies <psi_theta,partial_i psi_theta>=0 in the fixed real
filter gauge. The nonzero holonomy is of its rank-two response plane.

## What has closed and what remains

The TET connection now has an explicit realization by physical response
vectors, a full-space quantum unitary, preserved density domains and a
fixed-observer spectral probability in the existing LMD representation.
Nonzero curvature is established for an actual unbinned density state,
not only for the finite exact fixtures.

This closes the declared quantum density-response transport bridge. The
control metric is a covariance/Fisher metric; no equation identifies it
with the original material-slice metric exp(2r) dz² or with spacetime
curvature. LMD's unresolved inverse r,m construction and 3+1 constraints
remain unchanged. No source mass/energy coupling or physical time law
has been added by naming the bounded generator a Hamiltonian.

The next dynamical obligation is to show whether this response-plane
transport is induced by a specified original density/current coupling,
on its actual domains, or to prove the precise obstruction. A comparison
to material geometry must retain the conjugate data needed to reconstruct
r,m. A claimed universal quantum-gravity conclusion cannot be inferred
from the response-control curvature alone.
