# LPS-1 — A positive regular physical sector of the local constraints

16 October 2026. Predecessor: LQA-1 at 99a6547. We keep its periodic,
matter-free, V=0 source, reversed vacuum, shared zero modes and null
derivative improvement. This constructs a selected generalized solution
sector of **all** its Fourier constraints. It does not classify their full
kernel or derive a unique physical quantization. All predecessor packets
remain frozen. Units are hbar=1.

## 1. Auxiliary positive-energy factors; unchanged local constraints

For either chirality put w_n=v_-n, including w_0=v_0. Both u and w now
have [j_m,j_n]=2m delta_(m+n,0), positive Fock norm, j_n annihilating
the vacuum for n>0, and j_n*=j_-n. For a real zero-mode value t define

    T_m(j) = (1/4) sum_r :j_(m-r) j_r:,
    A_m(j) = T_m(j) + i m j_m/2 + (1/4) delta_m0.

**LPS-T1.** On the finite oscillator domain,

    [A_m,A_n]=(m-n)A_(m+n)+(7/12)(m^3-m)delta_(m+n,0),
    A_m*=A_-m,       A_0 = (t^2+1)/4 + sum_(n>0) n N_n,
    L_m^rev = A_m(u) - A_-m(w).                         (1)

Proof. LQA's positive-current calculation gives the T central term
(m^3-m)/12 and [T_m,j_n]=-n j_(m+n). The linear-term commutator adds
m^3/2. Writing the result in terms of A subtracts m/2, giving the
displayed central term. Adjointness includes conjugation of i. Substituting
w_n=v_-n in the difference cancels the two added 1/4 constants exactly
and reproduces LQA (18), including i m(u_m+v_m)/2. Thus no intercept
has been added to LQA or HQS. The auxiliary c=7 factors have equal and
opposite contributions to the constraint algebra, whose c is still zero.
Their positive A_0 grading is not a newly identified physical energy.

Take the regular branch

    P=2a_0 !=0,       K=2 pi b_0=0.

Then u_0=w_0=P/2 and both factors have

    h(P)=P^2/16+1/4 >0.                                (2)

The two original chiralities remain distinct and share P,K. We will solve
each of their u/w pairing constraints, not glue the two chiralities to
one another or introduce two independent zero-mode momenta.

## 2. Generalized solutions without an irreducibility assumption

At fixed real t, let V_t be the cyclic Virasoro submodule generated from
the normalized Fock vacuum by the A_-n, n>0. Its level-l subspace is
finite dimensional. A_m preserves V_t by commuting through a finite
descendant word. By the adjoint relation it also preserves its graded
orthogonal complement in the finite Fock space. We need neither the claim
that V_t exhausts Fock space nor a Kac determinant or irreducibility theorem.

The Gram matrix of ordered descendant words is real and depends only on
c=7 and h: commute positive modes to the right with (1). Its null vectors
are actual zero vectors, because the inherited Fock norm is positive.
Quotient by this null space at each level, and choose a real orthonormal
basis in the descendant-word coordinates. Apply the same real coefficients
in the two copies, obtaining e_(l,alpha)^u and e_(l,alpha)^w. A fixed
ordered Gram-Schmidt procedure, skipping zero residuals, suffices. The
Gram entries are polynomials in t; this choice is Borel measurable even
if ranks change. No continuity across a hypothetical rank change is needed.

**LPS-T2 (all-mode distributional kernel).** The formal paired vector

    I_t = sum_(l>=0) sum_alpha e_(l,alpha)^u tensor e_(l,alpha)^w   (3)

is a well-defined functional on finite oscillator polynomials and obeys

    (A_m(u)-A_-m(w)) I_t=0       for every integer m.     (4)

Here the ket pairs anti-linearly against test vectors using the Fock
inner product. The bra R_t(f)=<I_t,f> is linear. Both conventions solve
the constraints because the family is closed under adjoints.

Proof. Only finitely many levels pair with a finite polynomial. In the
chosen real descendant basis every A_m has real matrix entries: the
commutation relations, c,h and the basis changes are real. Adjointness
therefore says its matrix transpose is the matrix of A_-m. Equating the
coefficient of e_a tensor e_b on the two sides of (4) proves cancellation.
Tests in a graded orthogonal complement give zero on both sides by its
invariance. This proves the identity on the full finite polynomial test
space, not only on descendant tests. At level l the norm of the partial
paired vector is sqrt(dim V_(t,l)), bounded by sqrt(p(l)), where p(l) is
the number of integer partitions. Hence its action is uniformly bounded
on each finite oscillator subspace as t varies. In particular it is a
continuous functional for the locally convex direct-sum test topology.

Equivalently, if a descendant Gram matrix G_l is invertible,

    I_(t,l)=sum_(lambda,mu) (G_l^-1)_(lambda,mu)
                     d_lambda^u tensor d_mu^w.          (5)

If it is singular, use the real Gram quotient above. The all-level proof
does not extrapolate invertibility from the finite fixture calculations.
The tensor in (3) is independent of simultaneous real orthogonal changes
of basis. Arbitrary complex phases in only one copy would not preserve it.

The vacuum coefficient is one, so the solution is nonzero. It is not a
Fock vector: ||A_-1^l Omega||^2=l! product_(r=0)^(l-1)(2h+r)>0,
as follows inductively from the sl(2) subalgebra. Thus every level has at
least one basis vector and ||sum_(l<=M) I_(t,l)||^2 >= M+1. This is a
proof of divergence for these solutions, not a no-go for every possible
Fock-space solution. It also explains why the finite-polynomial no-kernel
theorem in LQA is entirely compatible with (3).

An explicit first nontrivial term is

    I_t = 1 + [(t-i)/(2(t+i))] U_1 V_1 + higher levels.   (6)

Indeed A_-1 Omega=(t-i)U_1/2 and its norm is 2h=(t^2+1)/2.
The phase in (6) is required by the original improvement. A naive pairing
with coefficient 1/2 fails; imposing a_n I=0 instead would force P=0
and would miss the regular HQS sector.

This is the established Ishibashi/contragredient pairing mechanism,
rederived here for LQA's specific currents and improvement. It does not
make this nonrational module into a rational CFT, impose Cardy sewing
conditions, or identify a physical spacetime boundary. For background see
Behrend, Pearce, Petkova and Zuber, hep-th/9908036v3, section 2.1.

## 3. Test topology, rigging prescription and physical norm

Use the dense invariant subspace of LQA's kinematic domain

    Phi = C_c^infinity((R\{0})_P x R_K) tensor_alg F_fin,

where F_fin contains both original chiralities and both oscillator species.
Give it the inductive-limit topology over compact zero-mode supports and
finite oscillator subspaces, with the usual smooth test-function topology.
The constant Jacobian between (a_0,b_0) and (P,K) is absorbed into a
unitary kinematic normalization. Each L_m^epsilon and each finite product
maps Phi continuously into itself. This follows directly from LQA's finite
action argument: a fixed finite set of input modes gives a finite set of
output modes, with coefficients polynomial in P,K.

At K=0 define the product of the two independently paired solutions

    I(P)=I_(P/2)^+ tensor I_(P/2)^-,
    (R f)(P)=<I(P), f(P,0)>_osc.                       (7)

This is a Borel function with compact support away from P=0. The level
bounds above make R continuous into L^2(dP/|P|) on each test-space step.
Theorem T2 gives R L_m^epsilon=0, including m=0. Vacuum tests give
R[chi(P)b(K)Omega]=chi(P) if b(0)=1, so R has dense range.

We explicitly select the vacuum-cyclic pairing (3), its vacuum coefficient
one, and HQS's regular-clock averaging prescription. These are quantization
inputs. The local constraint algebra alone does not uniquely select this
pairing sector, its momentum weight, or the full-line observable matching.

**LPS-T3 (positive physical completion).** For these inputs the form

    B(f,g)=integral_(P!=0) dP/|P| (Rf)(P)^* (Rg)(P)     (8)

is a positive constraint-annihilating rigging form. Its null quotient
completion is exactly

    H_phys^reg = L^2(R\{0},dP/|P|).                    (9)

Derivation of the measure, independently of the divergent Fock norm:
define the finite-test contraction

    (T f)(P,K)=<I(P), f(P,K)>_osc.

This auxiliary map uses the on-shell paired tensor even off shell; it is
not claimed to intertwine the off-shell constraints or to be a bounded
map on completed kinematic Fock space. It has compact support, bounded
P dependence and arbitrarily many uniformly bounded K derivatives on
each test-space step. Therefore the integral

    integral ds/(2 pi) integral dP dK
                  (Tf)^* exp(i s P K) (Tg)             (10)

exists: repeated integration by parts in K bounds the inner integral by
any inverse power of |s| since |P| is bounded away from zero. Fourier
inversion, or the Gaussian approximate identity with y=PK, reduces (10)
to (8). The absolute Jacobian is 1/|P|. This is a specified oscillator
contraction followed by the inherited one-constraint average, **not** a
Haar average over the full local gauge group or a derivation that such a
group representation exists.

Positivity now follows from (8), and its null space is exactly ker R.
Dense range, established by vacuum tests, proves (9). Define eta(g)[f]
=B(f,g), an anti-linear functional of f; it is continuous on Phi.
For every mode, B(f,L_m g)=B(L_m f,g)=0. Every chi in (9) likewise defines
a generalized solution eta_chi[f]=integral (Rf)^* chi dP/|P| by
Cauchy-Schwarz. This extension is injective by dense range. Positivity
is thus a property of a constructed physical quotient, not inherited by
declaring the infinite tensor (3) to have finite Fock norm.

In distributions the solution has K-support delta(K) and oscillator
coefficients supplied by I(P), with the weight dP/|P|. This is a
nontrivial local-constraint dressing of the zero-mode packet, not the
oscillator vacuum alone. No statement that these are all generalized
solutions is made; other cyclic modules, intertwiners, shifted zero-mode
supports and the singular P=0 branch need classification.

## 4. Actual HQS intertwiner and physical observables

**LPS-T4.** Let R_H f=f(P,0) be HQS's restriction map. On its test space
the kinematic injection iota f=f tensor Omega obeys

    R iota = R_H,
    B_LPS(iota f,iota g)=B_HQS(f,g).                    (11)

It therefore descends to an isometry of the physical null quotients.
Its range contains all compactly supported smooth amplitudes away from
zero, so it extends to a surjective unitary U:H_HQS -> H_phys^reg,
which is the identity on momentum amplitudes. The **physical** solution
corresponding to iota f is eta(R_H f), including the entire paired tensor;
iota itself is not a map into a kinematic strong kernel. Equation (11),
the all-mode identity T2 and the observable intertwiners below are the
claimed controlled recovery. This proves equivalence of the selected
physical sectors, not equivalence of the full local theory to a homogeneous
truncation or a classical global gauge reduction.

Transport HQS's specified reduced observable representation along U. With
J chi=chi/sqrt(|P|), J is unitary onto L^2(R,dP). Define

    P_phys chi=P chi,      Dom P_phys={chi:P chi in H_phys^reg},
    Q_phys=J^-1 (i partial_P) J,
    Dom Q_phys=J^-1 H^1(R),
    Q_phys chi=i[chi'-chi/(2P)]  away from zero.         (12)

Full-line H^1 matching through zero is a selected domain, as in HQS.
Writing A=F^-1 J, with the Fourier convention of HQS,

    L_phys=A^-1 [L_ref exp(q)] A,
    Dom L_phys={chi:exp(q) A chi in L^2(dq)}.           (13)

These are self-adjoint by unitary transport of real multiplication or
the full-line momentum operator. L_phys is positive with spectrum
[0,infinity) and has no zero eigenvector. On A^-1 C_c^infinity(R_q),

    [Q_phys,P_phys]=i I,       [L_phys,P_phys]=i L_phys.  (14)

This is a common invariant core for the displayed products; the domain
statements and density follow from HQS's proved reduced construction.
For each such observable O, the map eta_chi -> eta_(O chi) is well-defined
on chi in Dom O and preserves the all-mode distributional solution space.
Its physical adjoint relation follows from (8), not from a naive derivative
of oscillator coefficients. Moreover U O_HQS=O_LPS U, including domains.
This is a declared reduced observable algebra on the constructed local
solution sector. It does not yet identify L_phys with a relational
quantization of integral exp(q(x)) dx in the original local field theory.
The companion LOCAL_GEOMETRY.md now proves this identification at the
classical regular-reduction level and quantizes the resulting observable;
a direct quantum-field construction remains open.

Physical spectral projections of P, bounded finite-rank preparations and
their projective cuts consequently exist. For example HQS's two disjoint
packets chi_+=sqrt(|P|)1_[1,2], chi_-=sqrt(|P|)1_[-2,-1] lift to
normalizable physical states eta_(chi_+/-) of the local constraints.
For chi=(3 chi_+ +4 chi_-)/5, the sign probabilities remain 9/25,16/25,
<P>=-21/50 and <P^2>=7/3. Their measurability is now on this physical
solution space. These packets need not belong to Dom L_phys. No physical
propagation, time Hamiltonian, minimum length or empirical prediction is
inferred from them.

## 5. Closure and verification boundary

Closed under the stated choices: a nonzero distributional all-mode
solution family on a justified test dual; a positive rigging quotient;
self-adjoint reduced observables with domains; and an explicit unitary
intertwiner of physical states, norms and observable algebras with HQS.
The stronger equation-only comparison of LQA has been replaced by (11).

Open: uniqueness/completeness of the local kernel and rigging prescription;
uniquely selecting this sector from the local model; constructing quantum
relational field observables directly before reduction (the companion
LOCAL_GEOMETRY.md supplies the classical regular reduction and length);
integration of constraint smearings to a gauge group; singular-clock
sectors, potential, quantized matter, global gluing and physical 3+1 gravity.

The all-level pairing, topology, measure and self-adjointness arguments
above are general written proofs. The native certificate separately checks
finite descendant Gram matrices, inverse-Gram pairings, gluing against
finite tests, improvement-sensitive phases, auxiliary intercept cancellation,
norm growth and wrong pairings. No finite cutoff state is called an exact
all-mode state; no Gram checks are promoted to an all-level determinant
proof. Python 3.12 is I/O only; native exact arithmetic and Paninian word
rewriting perform the scientific fixture calculations.
