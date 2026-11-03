# A source-coupled native candidate and the inverse-square challenge

Monty Dabas — research continuation v1.0, 25 September 2026.

**Result.** An explicit source-coupled continuation can be built from U27's
native generators. Its stationary response follows an exact cut-flux ledger.
In a declared radial sector, quadratic shell capacity gives an inverse-square
finite-difference interaction response. The same native algebra also admits
other exponents, and a nonlinear cost can reproduce the distance exponent
while failing source proportionality. This gives a two-variable gravity test.
It is a conditional physical candidate, not a derivation of gravity from the
uncut primitive alone or an experimental validation.

## 1. Dependency audit: what is added

The [source manifest](SOURCE_RESPONSE_PINS.json) pins the earlier native
interaction, tensor and calibration results, the corrected thermo action and
non-selection results, and the represented graph/length lineage. The thermo
field equations already exist, conditional on a declared spacetime action.
They cannot supply a spacetime-free derivation merely by being cited. No
Einstein, Newton, Poisson or Born equation is used as a premise below.

This extension **does add constitutive hypotheses**. Avoiding a classical
field equation does not make those hypotheses consequences of the primitive.

| Item | Status in this candidate |
|---|---|
| Retained two-component carriers, local R generators and mixing brackets | Inherited from U27 |
| A finite connectivity graph and positive link stiffnesses | Declared; not physical space or inferred from U27 |
| An independently aligned response coordinate in each carrier | Declared calibrated representation, following U34 |
| Additive quadratic link mismatch cost and a linear source term | New candidate constitutive law |
| A grounded outer boundary, and a stationary preparation | Declared boundary/preparation conditions |
| Constant field-shift symmetry of the source-free bulk cost | Candidate symmetry; excludes an onsite pinning term |
| Shell-constant response and quadratic shell-capacity growth | Additional conditions for the inverse-square sector |
| Physical distance, source mass, work/force and probe calibration | Required empirical adapters; not supplied by a change of notation |

The word graph here means connectivity between retained carriers. RKF 35's
operator graph `(f,Cf)` is a different construction. The morphic geometry
graph reduction likewise does not select our weights or a three-dimensional
physical space. The native length result fixes length only inside its declared
sector and up to scale; it does not identify laboratory separation with a
shell index. All such identifications must be fixed independently of the
force data.

## 2. U36 — source-cost continuation and the exact stationary cut

Let a finite undirected graph have link stiffnesses c_uv>0. Ground a nonempty
boundary B; every connected component must meet B. At each vertex retain the
real two-component carrier of U27. First consider an aligned scalar field
phi, with phi|B=0, and supplied sources b on free vertices I. Define

    E(phi;b) = (1/2) sum_{links uv} c_uv (phi_u-phi_v)^2 - sum_{v in I} b_v phi_v.

Edge locality, an exactly quadratic cost, additivity and invariance under a
common scalar shift restrict each link's symmetric quadratic form to a
multiple of (phi_u-phi_v)^2. Positivity fixes the multiple's sign. Quadratic
degree and locality themselves are hypotheses, not selected by this argument.

Write E=(1/2)phi^T H phi-b^T phi on the free vertices. Then H is positive
definite and the unique minimizer is phi*=H^-1 b. Its edge current

    j_uv=c_uv(phi_u-phi_v)

satisfies, for every subset S of free vertices,

    sum_{u in S, v outside S} j_uv = sum_{u in S} b_u.                 (36.1)

Orientation is outward from S. This is a finite incidence identity and a
stationary response law; it does not presuppose a spatial surface integral.

**Proof.** The quadratic form is a sum of positive squared differences. A
zero-cost field is constant on each component and zero at its boundary, so
it vanishes. Thus H is positive definite. Differentiating E at each free
coordinate gives the sum of outgoing currents equal to b_v. Summing over S
cancels internal edges and proves (36.1). Completing the square gives

    E(phi;b)-E(phi*;b) = (1/2)(phi-phi*)^T H (phi-phi*) >= 0.

### An explicit native continuation law

For both coordinates put K_H=H tensor I_2, Z=I tensor R and G=H tensor R.
Then G^T=-G and JGJ=-G for J=I tensor K. Moreover

    G = sum_a H_aa G_a + sum_{a<b} H_ab [B_ab,G_a],

using exactly the local and bracket generators in U27. Hence the homogeneous
Cayley map T_h=(I+hG)(I-hG)^-1 is licensed by the same anti-self-dagger rule.
For a declared two-coordinate source b, define the affine continuation by

    (I-hG)x_next = (I+hG)x - 2h Zb,       h != 0.                    (36.2)

Its unique fixed point is x*=K_H^-1 b. It equals

    x_next = x* + T_h(x-x*).

Since K_H commutes with G, this map preserves
`(1/2)x^T K_H x-b^T x` and the squared distance from x* in the K_H metric.
The fixed point is exactly the stationary cost solution, not a separately
inserted inverse-square field. Selecting the stiffness matrix is still a
constitutive choice. The protocol parameter h is not a derived physical time.

**W51: stationarity is not healing.** A nonstationary initial state has positive
constant excess cost under (36.2), so it cannot converge to x*. A reservoir,
preparation or dissipative law is needed to realize settling. The conservative
native rule by itself does not prove minimum-cost relaxation in nature.

### Eliminating hidden field coordinates

Partition free coordinates into visible V and hidden M. Completing the square
over hidden coordinates gives the exact stationary effective cost

    H_eff = H_VV - H_VM H_MM^-1 H_MV,
    b_eff = b_V - H_VM H_MM^-1 b_M,
    E_eff(v) = (1/2)v^T H_eff v - b_eff^T v
               - (1/2)b_M^T H_MM^-1 b_M.                           (36.3)

H_eff is positive definite, its minimizer is the visible part of H^-1 b, and
its minimum cost equals the full minimum. This follows directly by completing
the hidden square. Keeping just H_VV erases a physical response in this
candidate. **W49:** a grounded chain with successive stiffnesses 1,4,9 has
effective source stiffness 36/49, not 1.

This is elimination at stationarity. It is not an exact quotient of every
time-dependent native word, so it does not bypass U28's closure obstruction.
In a more general graph it can create effective couplings between visible
vertices that had no direct link before the cut.

## 3. U37 — shell response and precisely when inverse-square follows

Consider a path of shell coordinates phi_0,...,phi_N, with phi_N=0. The link
from n-1 to n has aggregate stiffness C_n>0. Put source Q at 0 and no other
source. This is itself a declared graph model. It is also an exact reduction
of a larger network when its solution is constant on each shell and all links
cross consecutive shells: C_n is then the sum of their stiffnesses. A symmetry
that is transitive on each shell and preserves the source, boundary and weights
is a sufficient condition, by uniqueness. Shell constancy is not automatic.

For every shell link the exact solution is

    C_n (phi_(n-1)-phi_n)=Q,
    phi_k=Q sum_{n=k+1}^N 1/C_n.                                  (37.1)

Equation (36.1) gives the first identity; summation from the grounded boundary
gives the second. The outer boundary affects phi but not the drops for fixed
Q and C. No continuum, physical area or angular coordinate has been assumed.

To isolate the source-probe interaction from the probe's own boundary energy,
use the polarization of the stationary cost:

    E_cross(b,p)=E_min(b+p)-E_min(b)-E_min(p) = -b^T H^-1 p.         (37.2)

For source Q at 0 and a probe charge q at k, this is -q phi_k. Here q is a
candidate coupling, not inertial mass by definition. Across a separately
calibrated displacement ell_n>0, the inward signed finite-difference response is

    F_n = -[E_cross(n)-E_cross(n-1)]/ell_n
        = -q Q/(ell_n C_n) = -q Q/A_n,      A_n := ell_n C_n.       (37.3)

The minus sign means lower interaction cost nearer the source for qQ>0.
It gives attraction only under the physical work/force adapter. Opposite signs
repel. The dimensionful multiplier between native cost and measured energy
is an additional calibrated factor, omitted from the native units here.

Call A_n the **effective shell capacity**. It is not defined to be a Euclidean
area. For nonzero fixed q,Q, equation (37.3) is inverse-square in independently
assigned sample radii r_n exactly when

    A_n = A_0 (r_n/r_0)^2.                                       (37.4)

This is an iff statement about this candidate, not a proof of (37.4). Under
the constitutive identifications Q=k_s M, q=k_p m, the magnitude becomes
`G_eff M m/r_n^2`, with `G_eff=k_s k_p r_0^2/A_0` in native energy units.
The mass identifications, universal probe ratio and G_eff are not derived.
In a physical continuum reading one must additionally show that the sampled
finite-difference response approaches the measured local force. No finite
graph is silently called a smooth three-dimensional space.

**W45:** ell=1, C_n=n^2, q=Q=1 gives inward responses
`-1, -1/4, -1/9, -1/16`. Source Q->4Q multiplies every value by four.
This is an exact synthetic example satisfying the added growth hypothesis.

**W46: exponent non-selection.** With the same two-component native carrier
rules and the same quadratic cost, C_n=n^alpha, ell=1 gives `|F_n|=qQ/n^alpha`
for alpha=0,1,2,3. Exponential capacities are also allowed. Thus the native
grading/Cayley rules do not select inverse-square response. Total flux alone
also does not guarantee equal response on individual links of a nonsymmetric
shell; the certificate includes such a counterexample.

**Stability.** If A_n=A_0(r_n/r_0)^2(1+delta_n), |delta_n|<=epsilon<1,
then relative to the ideal magnitude F_n^0,

    F_n/F_n^0 = 1/(1+delta_n),
    |F_n/F_n^0-1| <= epsilon/(1-epsilon).

This transfers an independently justified capacity bound to a response bound;
it supplies no universal size or sign for a physical correction.

## 4. U38 — a distance law needs a source law and a fixed ruler

Changing the constitutive cost to the convex candidate

    E_p(phi;Q) = (1/p) sum_n C_n |phi_(n-1)-phi_n|^p - Q phi_0,
    p>1, Q>0,

gives the unique positive stationary drops

    C_n (phi_(n-1)-phi_n)^(p-1)=Q,
    phi_(n-1)-phi_n=(Q/C_n)^(1/(p-1)).                             (38.1)

**Proof.** With a grounded endpoint, link differences are independent
coordinates and phi_0 is their sum. The cost is a sum of strictly convex
one-variable terms `C_n |d_n|^p/p-Q d_n`. Differentiation gives (38.1).
The small-probe cross response is the first variation with respect to probe
source, `-q phi_k`, by stationarity; finite-probe superposition is asserted
only for p=2. Thus at fixed spacing the small-probe response is proportional
to `Q^(1/(p-1)) n^(-alpha/(p-1))` when C_n is proportional to n^alpha.

**W47:** p=3 and C_n=n^4 again give a 1/n^2 distance profile. But multiplying
the source by four multiplies the response by two, whereas the quadratic
candidate multiplies it by four. The two-variable test distinguishes them
without fitting a new amplitude at each source setting. No p=3 native
continuation generator is claimed; it is a competing positive mismatch-cost
model excluded only after adopting or testing the quadratic constitutive law.

**W48:** adding `(mu/2) sum phi_v^2`, mu>0, retains positivity but changes the
stationary matrix to H+mu I. Source-free vertices then absorb `mu phi_v` in
the cut ledger; radial outward flux need not be constant. Excluding this term
requires the candidate's bulk field-shift symmetry. The repo's graph-norm
term I+C* C must not be silently replaced by C* C to obtain unscreened gravity.

**W50:** without an independently calibrated ruler, a 1/n profile becomes
1/r^2 by the fitted assignment r=sqrt(n). Defining distance from the response
and then claiming an inverse-square prediction is circular. Likewise defining
Q from the observed force makes source proportionality vacuous.

The proposed observational discriminants are therefore simultaneous:

| Controlled change, others fixed | Quadratic cost and quadratic capacity | Cubic cost and quartic capacity |
|---|---|---|
| Independently calibrated distance r -> 2r | Response magnitude /4 | Response magnitude /4 |
| Independently calibrated source charge Q -> 4Q | Response magnitude times 4 | Small-probe response magnitude times 2 |
| Change probe material at equal independently measured inertial mass | Equality requires a common q/m ratio | Also not supplied by the cost law |

These ratios need an apparatus model for source extent, backgrounds, gain,
boundary and displacement uncertainty. They are hypotheses to test, not new
measured facts. [The experimental contract](GRAVITY_EXPERIMENT_CONTRACT.md)
specifies the rejection gate and current evidence status.

## 5. What the inverse-square challenge establishes here

The candidate law is now explicit, with a source term, stationary solution,
interaction cost, native affine continuation and sign-sensitive response.
It supplies a conditional inverse-square sector and controls that reject
several tempting shortcuts. Selecting the quadratic cost, its unscreened
symmetry, the physical source/probe couplings, the independent ruler and
quadratic capacity growth from primitive seam data remains unfinished.

Accordingly, the challenge **derive physical inverse-square gravity from the
existing native primitive without extra hypotheses is not solved in v1.0**.
Its exact missing hypotheses and a test that can discriminate candidate laws
are now stated. The familiar quadratic-network and elimination mathematics
is not claimed as a new general theorem; the native generator construction,
dependency audit and coupled source/distance controls are the present work.

Quantum probability, light propagation, relativistic radiation, empirical
gravity validation and quantum-gravity completion are not consequences of
this static sector. Their absence is not repaired by calling the source cost
curvature. The formalization is a representation after a cut; it is not an
identification of literal uncut ground with a finite graph.
