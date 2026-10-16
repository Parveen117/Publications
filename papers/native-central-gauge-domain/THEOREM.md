# CGD-1 — Central gauge reduction, a common source core and a test-space obstruction

16 October 2026. Frozen parent: IGQ-1 at
fa98e69a2c7b685187a2c5df3dd843193d1dbe11. The MRF action, IGQ inverse
embedding measure, off-shell canonical correction and selected strongly
invariant probe representative are fixed. No new matter or interaction is
introduced. Python is I/O only; continuum results below are written
theorems, not consequences of the finite native certificates.

There are two advances. First, an actual central subgroup admits positive
averaging compatible with the entire remaining gauge action. After this
reduction the source constraints have a dense common invariant core and
their all-smearing algebra closes without a central term. Second, the
natural test space used for this first average has *zero* full continuous
invariant dual. This is a theorem about that topology. It shows why a
stronger test topology in the shape variables is required for further
generalized physical states; it does not prove their nonexistence in
other riggings.

## 1. The central subgroup and the correct quotient

IGQ uses

    G=(Diff_tilde_+^infinity(S1) x Diff_tilde_+^infinity(S1))/Z_diag,
    chi_epsilon(x)=alpha_epsilon+I_sigma_epsilon(x),
    beta=alpha_+ mod 1, delta=alpha_--alpha_+,
    nu=d beta rho(delta)d delta dmu_0(sigma_+)dmu_0(sigma_-),
    rho(delta)=(2 pi)^(-1/2)exp(-delta^2/2).              (1)

The shapes are in H^3_0, with real Fourier covariance (1+n^2)^(-4).
The gauge action is chi_epsilon -> g_epsilon composed with chi_epsilon.
Its unitary implementation is the positive square-root Radon--Nikodym
representation. The actual source moment maps are

    T_+=F_+'Pi_+, T_-=-F_-'Pi_-, F_epsilon=chi_epsilon^-1. (2)

Let z be the class of (x->x+1,id). Every degree-one lift commutes with
integer translations, so Z_c=<z> is a central copy of Z. The projection
onto the two circle maps has precisely this kernel:

    1 -> Z_c -> G -> G_c -> 1,
    G_c=Diff_+^infinity(S1) x Diff_+^infinity(S1).        (3)

The element z acts by beta->beta and delta->delta-1. It is the endpoint
of the plus rotation flow, hence already belongs to the connected source
gauge group. Imposing its trivial action is a consequence of full gauge
invariance, not a new periodic physical-clock assumption or an extra
identification of the unreduced source configurations.

Flatten the Gaussian relative coordinate by J psi=sqrt(rho) psi. This
is a unitary from H_g=L^2(nu) to

    H_flat=L^2(d beta d delta dmu_0 dmu_0).              (4)

Set delta=d+n, 0<=d<1, n in Z. Then

    H_flat = ell^2(Z; H_c),
    H_c=L^2(d beta d d dmu_0 dmu_0)
       =L^2(d alpha_+ d alpha_- dmu_0 dmu_0).           (5)

The angles alpha_+=beta, alpha_-=beta+d are now taken modulo one.
Each Haar factor is normalized. The generator z acts as the bilateral
shift (S v)_n=v_(n+1). No normalizable fixed vector is inserted.

## 2. A positive average with a full remaining gauge action

Define the concrete dense Frechet test space

    Phi_Z = J^-1 s(Z;H_c),
    p_r(v)^2=sum_n (1+n^2)^r ||v_n||_Hc^2, r=0,1,...,
    P_Z v=sum_n v_n.                                   (6)

These are rapid weighted L2 tests in the integer relative coordinate,
with the Hilbert topology in the other variables. They need not be
differentiable in delta or in the shapes. All statements about this
space initially use finite gauge transformations, not an unjustified
infinitesimal domain.

**CGD-T1 (central periodization and residual unitaries).** Phi_Z is
invariant under every finite source gauge transformation. The series
P_Z is continuous and onto H_c, and

    sum_k <f,U_(z^k)h>_Hg
        = <P_Z Jf,P_Z Jh>_Hc.                          (7)

The series on the left is absolutely convergent. The form is positive;
its null quotient completes to H_c. There is a genuine unitary
representation U^c of G_c such that

    P_Z J U_g = U^c_(projection g) P_Z J.              (8)

Proof. Cauchy--Schwarz against (1+n^2)^(-r) proves convergence and
continuity in (6) for integer r>=1. A continuous right inverse E:H_c->
s(Z;H_c) is (Ev)_0=v, (Ev)_n=0 for n!=0. The double sum in (7) is
absolutely summable, since sum_n ||f_n|| and sum_n ||h_n|| are finite.
Reindexing its terms gives (7). Surjectivity and the norm formula give
the completed quotient.

For descent, the measure in (4) is exactly the lift of the measure in
(5) with counting measure in each Z_c fiber. A degree-one gauge lift
acts on those fibers by an integer translation. Its Radon--Nikodym
derivative for (4) is therefore the pullback of the derivative for (5);
it is independent of n. The shift integer is a measurable function of
the base angles. For each fixed g it is uniformly bounded: the smooth
periodic functions g_epsilon(x)-x are bounded, so |delta_new-delta| is
bounded. This gives p_r(JU_gJ^-1 v)<=C_(g,r)p_r(v). Thus Phi_Z is
invariant. Reindexing the fiber sum gives (8), first on finite sequences
and then by continuity. The induced measure on (5) is the product of
the two Haar-angle/inverse-shape measures. Its Radon--Nikodym unitaries
are exactly U^c, independently of the chosen lifts.

This succeeds because Z_c is normal. IGQ's real rotation subgroup was
not normal, so its earlier average could not provide (8). Moreover (7)
is covariant under the remaining gauge group, not already annihilated
by all of it. Covariance of a positive form is not the full Dirac
constraint condition in either argument.

Generalized Z_c-invariant histories are constant sequences in the flat
description, or rho(delta)^(-1/2) times a periodic base wave in the
original description. Their physical norm for this first reduction is
the H_c norm. Nonzero such histories have infinite original kinematic norm.

## 3. Exact coinvariants and a full-dual obstruction for this topology

**CGD-T2 (complete continuous invariant kernel on Phi_Z).**

    ker P_Z = (I-S)s(Z;H_c),
    s(Z;H_c)/(I-S)s(Z;H_c) is topologically H_c.         (9)

Every continuous Z_c-invariant linear functional on Phi_Z is uniquely
ell(f)=<v,P_Z Jf> for some v in H_c. The only continuous functional on
Phi_Z invariant under the *whole* G is zero. In particular no nonzero
separately continuous full Dirac rigging form exists on this Phi_Z,
positive or otherwise.

Proof of (9). If sum_n v_n=0, set

    h_n=sum_(k>=n) v_k = -sum_(k<n) v_k.                (10)

The first expression controls the positive tail and the second the
negative tail. Weighted Cauchy--Schwarz gives

    p_r(h)<=C_r p_(r+2)(v).                            (11)

For example at large positive n, the bound for ||h_n||^2 uses
sum_(k>=n)(1+k^2)^(-r-2)=O(n^(-2r-3)); multiplication by
(1+n^2)^r leaves a summable O(n^-3). The negative tail is the same.
Thus h is a rapid sequence and (I-S)h=v. Conversely summing a shift
difference gives zero. Together with the continuous right inverse E,
this proves the topological quotient. An invariant functional kills
I-S, so factors through P_Z. Its restriction to E(H_c) is Hilbert
continuous; Riesz gives precisely the claimed v.

There are no nonzero G_c-invariant vectors in H_c. If v were such a
vector, |v|^2 times the probability measure in (5) would be a finite
G_c-invariant measure. The alpha_+=chi_+(0) marginal would then be a
finite measure on S1 invariant under every smooth increasing circle
diffeomorphism. Rotation invariance makes it c times Haar measure;
a non-rigid diffeomorphism does not preserve Haar measure, so c=0.
The total mass ||v||^2 is zero. This argument also applies to an
arbitrary spectator Hilbert space by taking its squared norm.

Using (8), a full G-invariant functional in (9) would require just such
an invariant v and therefore vanishes. For a separately continuous
Dirac form eta(f)[h], fix f: its h functional is continuous and full
gauge invariant, hence zero. This proves the form statement without
assuming positivity.

The quantifiers matter. This is not a theorem excluding all generalized
physical states, all test topologies or all quantizations. The topology
in (6) permits distributions in the integer relative coordinate, but
does not permit unbounded distributional dependence on the residual
shape/angle Hilbert variables after periodization. A nonzero further
rigging must refine that residual test topology or use a different
justified reduction. Simply increasing the decay powers r in (6) cannot
help: all powers were already included in (9).

## 4. A finite differentiability extension supplies a common source core

The central quotient removes the need to apply a theorem about circle
diffeomorphism representations to a nontrivial covering group. Work on
H_c with the actual descended representation, not an assumed equivalent
representation of G on H_g. Let

    H_one=L^2(S1_alpha x H^3_0,d alpha dmu_0),
    H_c=H_one tensor H_one.                            (12)

Each circle factor acts only on its own copy.

**CGD-T3 (finite-C6 extension and a dense common invariant core).**
U^c extends to a strongly continuous unitary representation of
Diff_+^6(S1) x Diff_+^6(S1). Let H_one^infinity denote the smooth vectors
for its restriction to Diff_+^infinity(S1), in the usual smooth orbit-map
sense. Then

    D_c=H_one^infinity tensor_alg H_one^infinity         (13)

is dense in H_c, invariant under the smooth gauge group and under all
its derived generators. Every finite mixed product of the smoothly
smeared chiral generators is defined on D_c. Each real smeared chiral
generator, and each real H[f],D[f], is essentially self-adjoint there;
its closure is the Stone/joint-spectral operator of the finite action.

Proof of the extension. In IGQ's Ramer argument, C-infinity regularity
of a fixed g is stronger than needed. For g in Diff_+^6,
log g' is C5 and g''/g' is C4. The displacement

    u_g(alpha,sigma)=P_0 log g'(alpha+I_sigma)

is continuous from H^3 to H^4, with continuous Cameron--Martin
derivative. The coefficient (g''/g')(chi) is H^4; multiplication by it
is bounded on H^4. IGQ's Hilbert--Schmidt estimate O(1/n) for normalized
H^4 Fourier directions is unchanged. The finite alpha derivative is
finite rank. The inverse diffeomorphism is also C6, so its derivative
gives the bounded inverse of I+Du. These verify the same Ramer hypotheses
for C6. The action depends continuously on (g,alpha,sigma) in the C6
and H3 topologies. The finite generating-partition construction of
Radon--Nikodym derivatives is jointly Borel. Hence the representation
is a Borel homomorphism into U(H_one) with its strong topology.

For completeness, a Borel homomorphism between Polish groups is
continuous here. For a small target neighborhood V with V^-1 V inside
a prescribed neighborhood W, second countability gives a countable
cover by translates of V. A corresponding inverse image A is
nonmeagre and has the Baire property. It is comeagre in a nonempty
open set up to a meagre set; comparing that set with a sufficiently
small translate proves that A^-1 A contains an identity neighborhood.
Its image lies in V^-1 V. This proves continuity at the identity and
then everywhere. The C6 diffeomorphism group is Polish (an open
separable Banach-manifold group with continuous inversion), as is the
strong unitary group on this separable Hilbert space.

Density in (13) uses a stated external theorem, not finite-dimensional
Garding smoothing with a fictitious Haar measure. Shimomura's published
2001 smooth-vector theorem says that, for a compact smooth manifold,
a unitary representation of the identity component of its smooth
diffeomorphism group that extends continuously to some finite-Ck
diffeomorphism group has dense smooth vectors. Apply it to M=S1 and
k=6 in each factor. SOURCE_MAP.md records the checked published
statement and the access limit; the earlier 1999 conjectural discussion
is not substituted for the 2001 theorem.

Smooth orbit maps remain smooth after a fixed group translation or a
Lie derivative. Thus (13) is invariant under group elements and derived
generators, and all finite products exist. To check the core assertion,
let A be any of the real Stone generators and restrict it to D_c.
This is a dense symmetric restriction invariant under exp(-itA).
A vector in a deficiency space would have matrix coefficients against
exp(-itA)v, v in D_c, satisfying a scalar differential equation with
an exponentially growing solution in one time direction. Unitarity
bounds the coefficients for all real t, so they vanish. Density forces
the deficiency vector to be zero. Both deficiency spaces vanish;
the restriction is essentially self-adjoint and its closure is A.
The same proof covers H[f],D[f], whose one-parameter groups are the
commuting products of the two chiral groups.

This is a common core after the proved central gauge reduction. It is
not a claim that arbitrary cylinders in H_g already lie in every
unreduced generator domain, or that the entire smooth-vector space is
automatically a suitable rigging space for nonzero physical states.

## 5. The actual all-smearing source algebra on that core

Let Q_epsilon(f) generate left composition of chi_epsilon by the flow
of the real smooth periodic vector field f, with
U^c_(g_f(t))=exp(-it Q_epsilon(f)). Retain IGQ's source signs

    T_epsilon[f]=-epsilon Q_epsilon(f),
    H[f]=T_+[f]+T_-[f], D[f]=T_+[f]-T_-[f],
    b(f,g)=f g'-f' g.                                  (14)

**CGD-T4 (common-core quantum source closure).** On D_c,

    [T_epsilon[f],T_epsilon[g]]=i epsilon T_epsilon[b(f,g)],
    [T_+[f],T_-[g]]=0,
    [H[f],H[g]]=i D[b(f,g)],
    [D[f],H[g]]=i H[b(f,g)],
    [D[f],D[g]]=i D[b(f,g)].                            (15)

There is no central term. These hold for arbitrary smooth smearings,
not only a tested Fourier range. At circle length one, the complex
mode combinations

    L_n^epsilon=(epsilon/(2 pi)) T_epsilon[exp(2 pi i n x)]

satisfy [L_m^epsilon,L_n^epsilon]=(m-n)L_(m+n)^epsilon
on the same complex core, with opposite chiralities commuting. Their
formal adjoint pairing on the core is L_n <-> L_(-n). Real smearings
have the self-adjoint closures proved above. No bounded-below physical
energy interpretation is assigned to mode grading.

Proof. Differentiate the genuine finite group law on smooth orbit maps,
where both mixed derivatives are now defined. For composition of
diffeomorphisms the group Lie bracket identified with Eulerian vector
fields is -b(f,g). Equivalently, the configuration fundamental fields
X_f(chi)=f composed with chi obey [X_f,X_g]=X_b, and the pullback
differential representation has the opposite sign. Consequently
[Q(f),Q(g)]=-i Q(b). Substitution of (14) gives the first line of (15);
the tensor factors give the second. Addition and subtraction give the
H,D relations. The exact finite group law has no phase, so its mixed
derivative supplies no scalar cocycle. This uses CGD-T3's common core;
IGQ's mere list of separate Stone domains would not justify the step.
The Fourier formula follows by inserting b=2 pi i(n-m)exp(2 pi i(m+n)x).

The constraints quantize the exact moment maps (2) of the frozen MRF
source. This result does not construct the separate nonlinear original
field operators q,p,phi,k,rho or an inverse-p ordering. It does not
transfer LQA's oscillator proof to the reference fields or establish
unitary equivalence with that older representation.

## 6. Physical observables and the remaining full pairing

**CGD-T5 (spectator observables and an exact remaining criterion).**
Tensoring (7) with MPQ's physical Hilbert space yields precisely
H_c tensor H_red after the central gauge reduction. MPQ's spectral
observables, their maximal domains and the declared clock/pointer pulse
pass through unchanged on the second factor. In particular the selected
source probe G_w^* descends to I tensor J(w). On D_c tensor_alg D_red,
the source algebra (15) and the specified MPQ calculation algebra act
on independent factors. Their unitary/spectral actions strongly commute.

Proof. Apply P_Z tensor I to simple tensors, then quotient its positive
form. Surjectivity of P_Z and completion give the claimed Hilbert product.
The quotient intertwines every bounded spectral function of an MPQ
observable. Its spectral measure therefore transports the full unbounded
operator domain as well. The same argument applies to the pointer and
clock factors already proved in MPQ/RPM.

Full source reduction still needs a stronger dense test space Psi_c in
H_c on which the residual group/generators act, a nonzero generalized
invariant kernel and a positive pairing. If its physical quotient is
one-dimensional, the tensor-product argument of IGQ-T6 recovers H_red
and its measurement instrument. If more channels survive, they must
remain in the physical space. CGD-T2 rules out doing this continuously
on the unrefined Hilbert residual topology of Phi_Z; it does not decide
the answer on Psi_c. Dense smooth vectors solve the operator-product
problem but do not by themselves choose the physical test-dual pairing.

The next target is therefore the *residual* shape/angle rigging, using
the actual descended representation and common core. A proposed test
topology must support nonzero invariant distributions while preserving
the gauge action and physical adjoints. The earlier partial averages,+a repeated classical quotient, or an abstract Hilbert-space isomorphism
cannot replace that derivation. Nonlinear quantum geometry, other
orderings, measure selection, critical-current patches, local apparatus
and physical 3+1 identification remain separate obligations. All frozen
predecessors are unchanged.
