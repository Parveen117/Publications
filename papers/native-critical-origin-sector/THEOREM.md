# COS-1 — Critical origin, observable separation and two averaging scales

16 October 2026. Frozen predecessor: CCS-1 at
34b5ed1772ba383eddf76531ce1d3b3b28345af5. Retain the smooth real periodic,
period-one, V=0 matter-free LNC source and the LQA reversed representation.
All zero-mode units and the dP dK normalization are those of CCS. No
counterterm, new potential or matter dynamics is introduced.

## 1. The critical classical orbit space has two labels

CCS-T1 implies p=0 and phi=Phi on P=0. On the critical stratum K=0,
choose any smooth periodic primitive r of k: r'=k. Set

    E_+=exp(q+r), E_-=exp(q-r), Z_+=integral E_+, Z_-=integral E_-,
    J_0=Z_+ Z_-, Q_0=(1/2)log J_0, L_0=L_ref sqrt(J_0).       (1)

The integration is over one period. Z_+,Z_->0. The additive ambiguity
r -> r+c rescales the two Z's reciprocally and leaves J_0 unchanged.

**COS-T1 (critical local quotient).** The connected local gauge orbits on
P=K=0 are classified by (Q_0,Phi) in R^2. Each orbit has exactly one
homogeneous representative q=Q_0,k=0,p=0,phi=Phi. The pullback of the
canonical one-form and two-form to this representative family is zero.
This does not make distinct values of Q_0 or Phi gauge equivalent.

Proof. The inherited plus chiral flow gives

    delta_+ q=(f c_+ +f')/2, delta_+ k=(f c_+ +f')'/2,
    c_+=q'+k, c_-=q'-k.

Thus delta_+ r=(f c_+ +f')/2+alpha, where alpha is a spatial constant
depending on the chosen primitive convention. Consequently

    delta_+ E_+=(f E_+)'+alpha E_+, delta_+ E_-=-alpha E_-.

For the minus chiral flow the corresponding expressions are

    delta_- E_+=alpha E_+, delta_- E_-=-(f E_-)'-alpha E_-.

Periodicity gives delta J_0=0 in both cases. Phi is invariant because
p=phi'=0. This proves invariance under all connected smooth gauge flows.
CCS's zero-mean affine straightening applies independently to c_+,c_-
and sends both to zero, hence k=q'=0. On that homogeneous slice J_0=e^(2q),
so q=Q_0. Any two fields with the same Q_0,Phi therefore reach the same
representative, while unequal labels cannot be gauge related. Every real
pair occurs. Finally integral(p dq+k dphi)=0 on this family, proving the
two-form statement. The kernel of a form restricted to a singular stratum
must not be substituted for the original Hamiltonian gauge distribution.

In particular L_0 is a positive critical-stratum geometric invariant.
It is generally different from the raw circumference L_ref integral exp(q).
For example E_+=(1+u cos(2 pi x))^2, E_-=1, |u|<1 gives

    (L_raw/L_ref)^2=1, (L_0/L_ref)^2=1+u^2/2.                  (2)

The formulas use q=log(1+u cos(2 pi x)), r=q, k=r' and phi constant,
so they are actual smooth critical constraint solutions.

## 2. Which origin information extends to the whole constraint surface?

On homogeneous configurations the actual local source restricts to

    Theta=P dQ+K dPhi, C=PK,
    delta Q=N K, delta Phi=N P.                              (3)

These constant-lapse flows are a subset of the full local gauge flows.
All homogeneous tuples with PK=0 solve every local constraint. Their
ordinary parameter topology embeds continuously in the smooth field
topology. Thus a continuous gauge-invariant observable of the full local
constraint surface restricts to a continuous invariant function here.

**COS-T2 (failure of continuous separation at the origin).** Every such
observable takes the same value on all critical orbits (Q_0,Phi). Hence
the orbit quotient by the original connected gauge group is not Hausdorff.
The separate critical invariants Q_0,Phi,L_0 cannot be extended as
continuous invariant functions to the full constraint surface.

Proof. For every nonzero P, (Q,P,Phi_1,0) and (Q,P,Phi_2,0) lie on the
same orbit: use the constant-lapse flow parameter (Phi_2-Phi_1)/P.
Let P tend to zero. Continuity forces equality of the observable at
(Q,0,Phi_1,0) and (Q,0,Phi_2,0). Similarly the K!=0 flows show equality
at (Q_1,0,Phi,0) and (Q_2,0,Phi,0). Combining these two equalities makes
the origin value independent of both labels. The flow parameter diverging
in these sequences is allowed; each pair at nonzero momentum is exactly
gauge equivalent. It does not create a finite gauge flow between the
distinct origin representatives of T1.

More explicitly, any saturated open neighborhoods of two origin orbits
with the same Q and different Phi intersect: small equal nonzero P
homogeneous representatives in the two neighborhoods are gauge related.
These distinct quotient points cannot have disjoint neighborhoods. This
proves non-Hausdorffness directly, without confusing failure of scalar
separation with equality of the underlying orbits. The argument with K
is identical. Local gauge transformations can only enlarge these orbit
relations, and T1 still guarantees the distinct origin labels.

This theorem specifies an observable contract, not deletion of classical
information by fiat. An observer restricted to globally continuous Dirac
observables cannot distinguish the critical labels. An observer algebra
defined only on the critical stratum can distinguish them. Quantizing the
latter requires an additional representation/measure choice; its zero
restricted two-form supplies no canonical pair or preferred measure.

## 3. Exact logarithmic normalization of the origin average

Let Phi_all be LSC's compact smooth zero-mode tests tensor the finite
oscillator domain, with its stated inductive-limit topology. Write
a_f=rho f, e_0(f)=a_f(0,0). The inherited reduction is smooth, continuous,
surjective on scalar vacuum tests and satisfies rho(L_m^epsilon f)=PK g
for smooth compact g, with an irrelevant fixed normalization absorbed.
Define the Gaussian-regulated scalar average

    B_epsilon(f,g)=integral dP dK D_epsilon(PK) conjugate(a_f) a_g,
    D_epsilon(y)=exp[-y^2/(4 epsilon)]/sqrt(4 pi epsilon).       (4)

Equivalently this inserts exp(-epsilon s^2) into CCS's ds/(2 pi) average.
Take 0<epsilon<1 in the fixed units; changing a fixed reference inside the
logarithm below changes only a bounded term.

**COS-T3 (normalized origin limit).** For every f,g in Phi_all,

    B_epsilon(f,g)=log(1/epsilon) conjugate(e_0(f))e_0(g)+O(1),
    B_origin(f,g)=lim B_epsilon(f,g)/log(1/epsilon)
                =conjugate(e_0(f))e_0(g).                     (5)

The resulting null quotient completion is C with its positive usual norm.
It is the CCS origin form with normalization one, now derived from a
specified renormalization of the same regulator, rather than simply added.

Proof of the coefficient and bounded remainder. For h in C_c^infinity(R^2)
write h=h(0,0)chi(P)chi(K)+P h_1+K h_2, with smooth compact h_i and
chi=1 near zero. Such a compact decomposition follows from the radial
Hadamard formula near zero, multiplied by a cutoff, and division by
P^2+K^2 on the compact part away from zero. If the support is in [-R,R]^2,

    integral |P| D_epsilon(PK) |h_1| dP dK
      <= ||h_1||_infinity integral_(-R)^R dP <= 2R||h_1||_infinity,

because integrating the Gaussian over the whole K line gives 1/|P|
for P!=0. The K h_2 term has the same uniform bound.

For a square [-d,d]^2, direct integration gives

    integral_square D_epsilon(PK) dP dK
      =2 integral_0^(d^2/(2 sqrt(epsilon))) erf(t) dt/t
      =log(1/epsilon)+O(1).                                 (6)

Here erf(t)=2/sqrt(pi) integral_0^t exp(-u^2)du. The integral of erf(t)/t
is finite at zero, and integral_1^infinity (1-erf(t))/t is finite by the
Gaussian tail bound. Thus the large-upper-limit coefficient is exactly
two; after substituting d^2/(2 sqrt(epsilon)) it is one in log(1/epsilon).
Choose 0<=chi<=1, equal to one on a smaller square and supported in a
larger square. The positive square bounds give the same coefficient
with O(1) error for chi(P)chi(K). Applying the preceding decomposition
to h=conjugate(a_f)a_g proves (5), including off-diagonal pairs.

Constraint annihilation follows from e_0(L_m f)=0. Vacuum tests make
e_0 surjective. Its kernel is exactly the null space, proving the
one-dimensional completion and continuity. Smooth momentum multipliers
z(P,K) act on this line by z(0,0), with the correct complex adjoints;
P=K=0 as operators. No canonical conjugate with commutator iI exists on
this line. The functional e_0 includes CCS's full oscillator dressing;
it is not a bare vacuum evaluation on nonvacuum oscillator tests.

The regulator normalization is a quantization prescription. No assertion
of regulator independence beyond fixed changes of the logarithmic reference
is needed. The result is compatible with the continuous-observer collapse
in T2, but is not a faithful quantization of the two-label orbit space T1.

## 4. The largest smooth domain with a finite unscaled average

Define the continuous invariant test subspace

    Phi_van={f in Phi_all:e_0(f)=0}.                          (7)

It contains all punctured-plane tests used by CCS and is dense in the
kinematic Hilbert space, although it is a closed codimension-one subspace
of Phi_all's finer test topology. Do not confuse those two topologies.

**COS-T4 (maximal finite Gaussian-average domain).** A smooth test has a
finite limit of B_epsilon(f,f) iff it belongs to Phi_van. On this space
the unscaled limit is

    B_fin(f,g)=integral dP/|P| conjugate(a_f(P,0))a_g(P,0)
             +integral dK/|K| conjugate(a_f(0,K))a_g(0,K).      (8)

Both integrals converge at zero. The null quotient completion remains
exactly H_P direct_sum H_K with the CCS norms and matching conventions;
no extra regular states or origin atom appear in this unscaled completion.

Proof. Necessity follows from the positive leading term in T3. For
sufficiency, a_f(0,0)=a_g(0,0)=0 implies, on a small square,

    |conjugate(a_f)a_g| <= C(P^2+K^2).

Uniformly in epsilon,

    integral_[-d,d]^2 (P^2+K^2)D_epsilon(PK) dP dK <= 2d^2.  (9)

For the P^2 term integrate K over the whole real line, leaving
integral_(-d)^d |P| dP=d^2; interchange P,K for the other term.
The proposed axis limit has the same O(d^2) bound near zero. Outside
this square a smooth cutoff and CCS-T5's noncritical delta limit give
(8). Let d decrease to zero to finish the uniform near/far argument.

Restrictions of a_f to either axis vanish linearly at zero, so their
weighted norms are finite. Conversely the punctured vacuum tests already
have dense independent range in the CCS direct sum, proving the completion.
Null space consists exactly of tests whose two axis restrictions vanish.
Constraint invariance follows from rho L=PK g. Momentum-multiplier adjoints
are retained, as (8) is integration against the same positive axis weights.
The branch canonical operators keep their CCS domains; they need not map
every smooth amplitude in (7) into their operator domains. In particular,
a quotient completion is not an assertion that every test is a vector in
the domain of every unbounded observable. For example a smooth axis
amplitude equal to P near zero has finite norm in (8), but
chi/sqrt(|P|)=sign(P)sqrt(|P|) has a derivative outside L^2 near zero;
it is not in the inherited canonical Q domain.

## 5. Keeping both rates requires a filtered construction

The inclusions 0 subset Phi_van subset Phi_all give two different positive
forms: B_fin on Phi_van and B_origin on Phi_all/Phi_van. Their completed
associated pieces may be assembled as a declared Hilbert direct sum

    H_filtered=(H_P direct_sum H_K) direct_sum C.              (10)

This does not assert a single rigging map from Phi_all with all the
unchanged regular observables. Each piece has its own quotient map and
normalization. A relative positive factor for the origin summand remains
free when these two separately normalized scales are assembled.

**COS-T5 (no single scalar normalization or compatible subtraction).**
No scalar positive normalization r(epsilon) of (4) has finite nonzero
limits on both a nonnull regular test and a nonzero-origin test. Nor can
a continuous positive form on all of Phi_all agree with (8) on Phi_van and retain
the adjoint rules for all smooth scalar momentum multipliers.

For the first statement choose f in Phi_van with B_fin(f,f)=b>0.
A finite nonzero normalized limit implies r(epsilon)->c/b>0. For g with
e_0(g)!=0, T3 then forces divergence. Conversely a finite nonzero limit
on g forces r(epsilon)log(1/epsilon) to a nonzero finite constant, and
hence r(epsilon)->0, killing every regular limit.

For the second statement it suffices to restrict a putative continuous
positive form to scalar vacuum tests. Multiplier compatibility makes this
restriction a positive Radon measure, by the local positive-distribution
argument in LSC-T5. Agreement on punctured tests forces its measure there
to be dP/|P| and dK/|K| on the two axes. Those weights have infinite mass
in every neighborhood of the origin, contradicting local finiteness.
Adding any nonnegative atom cannot repair this contradiction.

The tempting subtraction construction shows exactly what would be lost.
Choose a smooth real compact scalar s with s(0,0)=1 and define

    a_f^s=a_f-e_0(f)s,
    B_s(f,g)=B_fin(a_f^s,a_g^s)+c conjugate(e_0(f))e_0(g), c>0. (11)

Here B_fin on scalar amplitudes means (8). This is positive, finite,
constraint-annihilating and agrees with B_fin on Phi_van, but violates
momentum symmetry. On the scalar tests s,Ps,P^2s,

    B_s(Ps,Ps)=integral |P| |s(P,0)|^2 dP >0,
    B_s(s,P^2s)=0.                                          (12)

These would have to agree if multiplication by P were symmetric. A smooth
compact multiplier equal to P on the test supports gives the same control.
Thus subtraction plus an atom is not an extension with the old observable
contract, even though its norm is positive.

## 6. The precise observable content of the filtered option

There is a concrete source-derived algebra across the strata, beyond the
momentum multipliers. Define on general smooth fields, with a fixed spatial
base point, r(x)=integral_0^x k, even when K!=0, and

    Q_*= (1/2)log[(integral exp(q+r))(integral exp(q-r))],
    Phi_av=integral phi, D_P=P Q_*, D_K=K Phi_av.              (13)

These are smooth real functionals; r need not be periodic off K=0.

**COS-T6 (global weak Dirac dilation observables).** On the full smooth
constraint surface, D_P,D_K are well-defined gauge invariants, independent
of the primitive/base-point convention in their physical values. Their
restrictions to the open branches are P Q and K Phi, respectively, and
both vanish on the critical stratum. They satisfy

    {D_P,P}=P, {D_K,K}=K,
    {D_P,K}={D_K,P}=0, {D_P,D_K} approximately 0.             (14)

Here Q is exactly the LPS reduced log-length, not a new coordinate assigned
to the first branch. In the filtered Hilbert option they have self-adjoint
quantizations iP partial_P and iK partial_K on their active summands,
acting as zero elsewhere, with the affine commutators iP and iK.

Proof of the classical identification. On the P!=0 branch K=0 and r is
periodic. The on-shell current ODEs and their mean P give exactly

    A_+=P exp(q+r)/Z_+, A_-=P exp(q-r)/Z_-.

Consequently LPS's Q=q-(1/2)log|A_+A_-|+log|P| equals Q_*.
Thus D_P is the regular invariant P Q there and is zero on P=0. P itself
is invariant under every local generator. The product formula proves weak
invariance also at P=0 without differentiating a discontinuous extension.
On K!=0, CCS gives p=0,phi=Phi, so D_K=K Phi; on K=0 it is zero.
K is invariant under all generators at V=0. These facts prove global weak
invariance, including the origin, directly on the actual local source.
Primitive/base-point changes cancel in Q_* on K=0, and P Q_*=0 on its
complementary constraint branch. There is therefore no physical base-point
dependence. Q_* by itself is not claimed to be a global Dirac observable.

The canonical P flow translates q by a constant, shifting Q_* by the same
constant; the K flow translates phi. This gives the first four brackets.
Q_* depends only on q,k, so {Q_*,K}=0, and {P,Phi_av}=0. Therefore
{D_P,D_K}=PK {Q_*,Phi_av}, which vanishes on the full smooth constraint
surface by CCS-T1. The coefficient is smooth. This also verifies that the
bracket statement does not infer a symplectic form from the critical
stratum's degenerate pullback.

For quantization, on H_P set

    (U_P chi)_sigma(y)=chi(sigma exp(y)), sigma=+1,-1.

This is unitary from L^2(dP/|P|) to L^2(dy) direct_sum L^2(dy).
The maximal domain U_P^-1[H^1(R) direct_sum H^1(R)] makes
hat D_P=U_P^-1(i partial_y)U_P self-adjoint. Its unitary dilation group is
chi(P)->chi(exp(t)P); t is a transformation parameter, not physical time.
It coincides on compact regular tests with
(P hat Q+hat Q P)/2=iP partial_P using the inherited weighted Q. The
operator just defined is the self-adjoint closure of that test expression,
not a claim that the unbounded anticommutator has that domain automatically.
The K construction is identical. Take their direct sums with the zero
operators on the other sectors. On the common invariant core of compact
smooth amplitudes on the two open axes plus the origin line, direct
differentiation gives [hat D_P,hat P]=i hat P, [hat D_K,hat K]=i hat K and
all cross commutators zero. The logarithmic coordinate sends zero momentum
to an infinite endpoint, so these generators require no finite-point
matching across zero. Their zero origin action is fixed by (13).

On (10), bounded smooth momentum functions have the explicit representation

    pi(z)(chi_P,chi_K,c)
       =(z(P,0)chi_P,z(0,K)chi_K,z(0,0)c).                   (15)

It is a unital star representation when constants are included. The real
unbounded momenta have the usual maximal multiplication domains on their
active regular summands and act by zero elsewhere. Their joint zero
eigenspace is exactly the added origin line: multiplication by nonzero
P or K has no zero eigenvector in its regular L^2 summand. The joint
spectral projector onto this line is therefore a well-defined bounded
observable, even though it is not a continuous scalar function of P,K.

One may choose the block-diagonal algebra generated by (15), this projector
and the dilation unitary groups just constructed. The already constructed
canonical branch observables can additionally act on their own summands.
It makes the summands superselected by that choice. No assertion rules
out all bounded mixing maps under a different observable contract. Any
bounded intertwiner from a regular summand into the origin that intertwines
all momentum multipliers is zero: its adjoint would give a regular vector
with spectral support at the single zero point, which has zero measure.

The origin line carries only scalar operators; it does not carry Q_0,Phi
as independently measurable labels. The stratum-only classical lengths
in T1 therefore cannot be assigned a universal origin eigenvalue from
this derivation. The regular log-length and scalar operators are retained
on their respective branches, not silently extended as canonical pairs
on the origin. A physical law coupling branches, a stratum-only observable
quantization, or interacting quantum fields would be additional work.

## 7. Scope, attribution and next obligation

Closed here: the critical classical orbit labels; their failure of global
continuous separation; the exact leading Gaussian logarithm and positive
origin quotient; the maximal finite smooth-test domain; the global weak
Dirac dilation algebra with self-adjoint domains; and the obstruction
to keeping both averaging rates in one unchanged rigging/observable contract.
The filtered direct sum is explicit and preserves both predecessor regular
Hilbert spaces, but its assembly is a declared two-scale construction.

Renormalized group averaging and divergence-related sectors are established
methods. Louko and Molgado studied their relation to singular reduced phase
spaces in a different constrained model (gr-qc/0505097); Louko reviews
convergence, positivity and superselection in gr-qc/0512076. Neither paper
is used as proof of the LNC formulas above. All source-specific coefficients,
orbit arguments and domain statements are derived here. No novelty claim
is made for renormalization or non-Hausdorff constrained quotients generally.

Next: select a physical observable contract across these strata and derive
any coupling or direct quantum relational-field action from the source.
Full local interacting quantum gravity, nonzero potential, quantized matter,
global cut gluing and physical 3+1 gravity remain open. The logarithmic
origin state alone supplies no dynamics, energy spectrum or experimental
prediction. Native finite Fourier, rational and word certificates are
witnesses for the displayed algebraic ingredients, not proofs of the
Gaussian asymptotics, topology or operator domains. Those are written proofs.
