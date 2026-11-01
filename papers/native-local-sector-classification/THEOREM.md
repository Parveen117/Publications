# LSC-1 — Complete test-dual kernel and observable-selected physical measure

16 October 2026. Frozen predecessor: LPS-1/2 at 064f593. Keep LQA's
periodic V=0 matter-free source, reversed oscillator vacuum, real shared
zero modes and original null derivative improvement. All results concern
this representation and the test space specified below. They are not a
classification of other representations, interactions or physical 3+1 gravity.

## 1. A positive-level obstruction to additional primary vectors

Use LPS's auxiliary positive oscillator j_n, [j_m,j_n]=2m delta_(m+n,0),
j_n*=j_-n, j_0=t real, and its unchanged operators

    T_m=(1/4) sum_r :j_(m-r)j_r:,
    A_m=T_m+i m j_m/2+(1/4)delta_m0.

The full finite Fock space at oscillator level N is denoted F_(t,N).
Its dimension is the number p(N) of integer partitions of N. On this
finite-dimensional space define

    D_N=sum_(n=1)^N n j_-n A_n,
    E_N=(1/2)sum_(n=1)^N n^2 j_-n j_n.

**LSC-T1 (all-level primary obstruction).**

    D_N=H_N+i E_N,             H_N*=H_N,
    E_N=sum_(n=1)^N n^3 N_n >= N I,                  (1)
    ||D_N v|| >= N ||v||       for N>=1.             (2)

Here the adjoint uses the positive Fock inner product, including its
monomial weights. In particular no nonzero positive-level vector is
annihilated by all A_n, n>0.

Proof. For positive n, normal ordering gives

    T_n=(t/2)j_n+(1/2)sum_(k>0)j_-k j_(n+k)
                      +(1/4)sum_(a=1)^(n-1)j_a j_(n-a).

Multiply by n j_-n and sum. Symmetrizing n,k in the first cubic sum yields

    H_N=(t/2)sum_(n>0)n j_-n j_n
        +(1/4)sum_(a,b>0)(a+b)
           [j_-a j_-b j_(a+b)+j_-(a+b) j_a j_b].      (3)

Every term acts finitely on level N; in the cubic sums only a+b<=N
can contribute. The two cubic words are adjoints, and the first term is
self-adjoint because t is real. The improvement supplies exactly i E_N.
Since j_-n j_n=2n N_n, E_N has the value sum n^3 N_n on a monomial.
At level N, sum n N_n=N and n^3>=n, proving (1). Consequently

    N||v||^2 <= Im<v,D_N v> <= |<v,D_N v>| <= ||v|| ||D_N v||.

This proves (2). If every positive A_n killed v, D_N v would vanish,
contradicting (2). No determinant formula, irreducibility theorem, finite
cutoff extrapolation or change of improvement is used in this proof.

The hypotheses matter. If the improvement is deleted, at t=0 the level-one
vector j_-1 Omega is killed by every positive T_n. If real t is replaced
by t=-i, the same vector is killed by every positive A_n; the adjoint and
Hermiticity assumptions then fail. These are actual controls, not extensions
of the theorem to nonreal momenta or another quantization.

## 2. Vacuum descendants exhaust the full oscillator representation

**LSC-T2 (cyclicity and graded irreducibility).** For every real t, the
ordered descendants A_-lambda Omega, lambda a partition of N, form a
basis of F_(t,N) at every N. Their Gram matrix is strictly positive and
invertible. The finite Fock module has no proper nonzero invariant
submodule under all A_m.

Proof by induction on N. Level zero is the vacuum. Suppose lower levels
are spanned by descendants. If their descendant span at level N were
proper, take a nonzero orthogonal vector v. For every n>0 and every
z in level N-n,

    <z,A_n v>=<A_-n z,v>=0,

since A_-n z is a level-N descendant. For n>N the image vanishes by
level. Thus v is a forbidden positive-level primary by T1. The span is
therefore the full level. There are exactly p(N) ordered descendants,
matching the monomial dimension, so they form a basis. Negative words
can be ordered using [A_-m,A_-n]=(n-m)A_-(m+n); the commutator
shortens the word and introduces no central term. Positivity of
Fock space makes their Gram matrix strictly positive.

For the invariant-submodule assertion, a nonzero finite vector has finitely
many A_0 eigencomponents. Polynomials in A_0 isolate these components.
Choose the smallest level present in the submodule. Positive A_n kill
that component, so T1 forces its level to be zero. It contains the vacuum
and therefore all descendants. This is an algebraic assertion; no claim
of essential self-adjointness of every infinite smearing is needed.

Descendant coefficients are polynomials in t with complex rational
coefficients. Their square change-of-basis matrices are invertible for
all real t, so their inverses are smooth on the entire real t-axis and
bounded, with derivatives of every finite order, on compact intervals.
This removes the rank-change possibility left open in LPS-1. It also
justifies smooth zero-mode changes of basis on each finite test-space step.

## 3. Classification at fixed real zero modes

Write a=a_0, b=b_0, u_0=a+b and w_0=a-b. The auxiliary weights are

    h_u=((a+b)^2+1)/4,   h_w=((a-b)^2+1)/4,
    Delta=h_u-h_w=ab=PK/(4 pi).                        (4)

For one original chirality the constraint is L_m=A_m(u)-A_-m(w).
Use the algebraic dual of the finite two-factor Fock module; it contains
all formal levelwise coefficient arrays, not only finite polynomials.

**LSC-T3 (fiber kernel).** Its simultaneous constraint-annihilating dual
space is zero if ab!=0 and is one dimensional if ab=0. At ab=0 it is
precisely the LPS descendant pairing, extended to u_0=-w_0 when a=0.
For both original chiralities the kernel remains one dimensional, the
product of their pairings, with a,b shared.

Proof. Let S be a linear solution bra: S(L_m f)=0 for every m,f.
Moving each negative generator on the u factor to the w factor gives

    S(A_-n x tensor y)=S(x tensor A_n y).             (5)

If a positive-level descendant remains on the w side with vacuum on the
u side, applying L_n to Omega tensor z shows its pairing is zero. Therefore
for descendants d_lambda,d_mu,

    S(d_lambda^u tensor d_mu^w)
       = delta_(|lambda|,|mu|) G_(lambda,mu)(h_w) s,  (6)
    s=S(Omega tensor Omega).

The Gram entry comes from commuting the positive word in (5) through the
right descendant to its vacuum. It is real and polynomial in h_w, with
auxiliary c=7 fixed. The L_0 equation on the vacuum gives Delta s=0.
T2 says these descendants span every test, so there are no other free
coefficients. If Delta=0, the equal-weight Gram construction in LPS proves
existence with s=1, even when t_u=-t_w. Repeating (5) in the second
chirality leaves only the common four-factor vacuum coefficient and gives
the product of the two Gram entries. Thus there is no additional chiral
degeneracy or family at integer nonzero Delta.

L_0 alone would allow grade offsets when Delta is an integer. The other
modes rule them out; equality of zero-mode grading is not a substitute
for the all-mode argument (5).

## 4. Full distributional classification, including the axes' intersection

Choose the dense invariant test space

    Phi_all=C_c^infinity(R_a x R_b) tensor_alg F_fin,

with the inductive-limit topology over compact zero-mode supports and
finite oscillator subspaces. F_fin includes both original chiralities.
This enlarges LPS's regular test space to allow the axes and their
intersection. Each finite product of L_m^epsilon acts continuously on it.
The constraints have no zero-mode derivatives: a,b act by multiplication.
Here Phi_all' denotes the continuous **linear** dual (solution bras).
Conjugation gives the equivalent anti-dual convention for solution kets.

Expand a test in the smooth full descendant bases supplied by T2, in all
four oscillator factors. Define a smooth scalar reduction rho by replacing
each two-factor descendant pair with the right Gram entry in (6), or zero
when its levels differ, and multiplying the results for the two chiralities.
The remaining coefficient functions are unchanged. Thus

    rho:Phi_all -> C_c^infinity(R^2),   rho(f Omega_all)=f.       (7)

This map is continuous on every finite test-space step. It agrees with
LPS's pairing contraction on b=0, and with the analogous equal-weight
pairing on a=0. It is a computational off-shell extension, not an off-shell
constraint intertwiner.

**LSC-T4 (complete continuous test-dual kernel).**

    {S in Phi_all': S L_m^epsilon=0 for all m,epsilon}
       = {F composed with rho: F in D'(R^2), ab F=0}.            (8)

Moreover every such scalar distribution has the form

    F=U(a) tensor delta(b)+delta(a) tensor V(b),                  (9)

where U,V are arbitrary one-variable distributions. The only ambiguity
in (9) is (U,V) -> (U+c delta,V-c delta).

Proof of (8). Define F(f)=S(f Omega_all). The vacuum L_0 identity gives
F(ab f)=0. Moving descendant words as in (5), now with smooth coefficient
functions, is legitimate because these functions commute with the modes.
It gives S=F composed with rho on the full test space; thus the map is
injective as well as necessary. Conversely, rho(L_m f) vanishes on both
axes by the equal-weight pairing proof. A smooth function vanishing on
a=0 and b=0 is ab times a smooth function, by applying the one-variable
Hadamard integral formula twice. The quotient is compactly supported when
the original function is: outside its support the quotient vanishes on
the dense set ab!=0 and hence by continuity everywhere. Therefore
rho(L_m f)=ab g with g in C_c^infinity, and ab F=0 proves sufficiency.

For (9), choose compact smooth cutoffs chi_a,chi_b equal to one near zero.
Every compact smooth f differs from

    f(a,0)chi_b(b)+chi_a(a)f(0,b)-f(0,0)chi_a(a)chi_b(b)        (10)

by ab times a compact smooth function. Applying F gives a sum of two
one-variable distributions on the axes, after absorbing the overlap term
into either one. If two such decompositions coincide, their difference
is supported at the origin in each coordinate. Independence of the
derivatives of delta shows that only opposite multiples of delta(a)delta(b)
can cancel between the two summands. This proves the stated ambiguity.

For example delta'(a)delta(b) and delta(a)delta''(b) are allowed scalar
coefficients; delta'(a)delta'(b) is not, since ab times it equals
delta(a)delta(b). These expressions specify the **vacuum coefficient F**;
the full solution is F composed with rho, whose zero-mode derivatives also
act on the oscillator dressing. They are not bare delta products tensored
with an unchanged vacuum. They have no asserted positive physical norm.

In regular coordinates P=2a!=0, K=2 pi b, (8) reduces to

    F=U(P) tensor delta(K),                            (11)

with arbitrary distribution U on R\{0}. No transverse delta'(K), shifted
K shell, extra oscillator module or independent chiral solution remains.
Thus LPS-1's pairing exhausts the oscillator dependence of the regular
continuous test-dual kernel. Its normalizable physical amplitudes form
the Hilbert subspace selected by its rigging norm; not every distribution
U is a normalizable state. The a=0 branch and its intersection have now
been classified algebraically, but not supplied with a physical norm,
clock, geometry or complete observable representation.

There is no nonzero ordinary kinematic Hilbert kernel: even L_0^+ alone
acts as real multiplication by ab+N_u-N_w in its orthogonal monomial
decomposition. Each zero set has two-dimensional Lebesgue measure zero,
and there are countably many grade pairs. Its maximal multiplication
realization therefore has zero kernel. This does not remove the
distributional solutions just classified or their physical Hilbert norm.

## 5. Which part of the physical measure is fixed?

On the regular branch let Rf=rho f(P,0). T2 makes R smooth, continuous and
surjective onto C_c^infinity(R\{0}); vacuum tests give a continuous local
right inverse. Consider continuous positive sesquilinear rigging forms
which annihilate the constraints in either slot. The classification above
implies they factor through R: a zero Rf is annihilated by every solution
functional, and hence lies in their null pairing with every test.

Additional physical input must be explicit. Require the full scalar
momentum-multiplier algebra to have its declared adjoints:

    B_amp(f,z g)=B_amp(conjugate(z) f,g)                (12)

for smooth compactly supported complex z. This condition implements the
usual interpretation of momentum functions as observables; it is stronger
than the local constraint algebra alone.

**LSC-T5 (conditional uniqueness of the regular physical measure).**
Under (12), a continuous positive amplitude form is integration against
a positive Radon measure mu on R\{0}. If in addition the **fixed inherited**
observable expression

    Q_0=i[partial_P-1/(2P)]                            (13)

is symmetric on compact smooth amplitudes, then

    dmu=c_+ dP/P   on P>0,
    dmu=c_- dP/|P| on P<0,       c_+,c_- >=0.          (14)

For a nondegenerate representation of both signs, c_+,c_->0. Requiring
Q_0 on HQS's already selected full-line matching domain

    chi/sqrt(|P|) in H^1(R)                            (15)

forces c_+=c_-. An overall normalization then recovers exactly dP/|P|.

Proof of the measure statement. Locally choose a real cutoff zeta equal
to one on the supports of f,g. Multiplier compatibility gives
B_amp(f,g)=B_amp(zeta,conjugate(f)g). This defines a compatible local
distribution ell(h), independent of enlarging the cutoff, which is positive
on h=|f|^2. It is positive on every nonnegative smooth test: add an arbitrarily
small positive cutoff square on its compact support, take a smooth square
root there, and let the added coefficient decrease to zero. A positive
distribution is a measure: on a fixed compact support, comparison with a
nonnegative cutoff bounds |ell(h)| by a constant times the sup norm of h;
extension to continuous tests and the Riesz representation theorem give a
positive Radon measure. These local measures agree on overlaps by (12).
This standard positive-distribution argument is also recorded in Patrick
Gerard's 2019 distribution notes, Proposition 1.9.11 and Remark 1.9.12;
see SOURCE_MAP.md.

Symmetry of (13), with arbitrary compact smooth f,g, gives

    integral [(conjugate(f)g)'-(conjugate(f)g)/P] dmu=0.

Every compact smooth h can be written conjugate(f)g using a cutoff f=1
on its support. Thus distributionally mu'=-mu/P. Multiplication by the
smooth positive |P| on either half-line gives (|P|mu)'=0. A distribution
with zero derivative on a connected interval is constant (test functions
of integral zero are derivatives of compact smooth functions). Positivity
makes that constant nonnegative, proving (14) without assuming a density.

Write chi=sqrt(|P|)psi and eta=sqrt(|P|)varphi. On the matching domain (15)
the integration-by-parts boundary form is

    B(Q_0 chi,eta)-B(chi,Q_0 eta)
          = i(c_+-c_-) conjugate(psi(0))varphi(0).     (16)

Full-line H^1 allows an arbitrary common trace at zero; hence symmetry
forces equality of the constants. Conversely equal positive constants
give exactly HQS's unitary full-line derivative realization, self-adjoint
on (15). Its positive relational-length functional calculus and domains
are then unchanged from LPS-1/2.

This is a compatibility theorem, not a derivation of the observable or
the matching domain from constraints. Without (13), every positive Radon
measure gives a constraint-annihilating form; for example dP and dP/|P|
are distinct in the fixed amplitude convention. Without the matching
condition, the two constants in (14) need not agree. If a different
weight-dependent expression for Q is adopted, further equivalent
representations exist. A newly required physical law has not been hidden
inside a claim that the algebra alone picks a unique inner product.

## 6. Research endpoint and finite verification contract

Closed for the stated representation/test topology: all-level absence of
additional primaries; full Fock cyclicity; fixed-fiber all-mode kernel;
the complete continuous zero-mode-distribution kernel, including its axis
decomposition; and uniqueness up to normalization of the regular physical
measure under the inherited observable and matching contract.

This upgrades LPS's regular pairing from a selected oscillator subsector
to the full oscillator dependence of the regular test-dual solution space.
It does not assert uniqueness across other quantum representations. The
singular axis has an algebraic solution description, not a finished physical
quantization. Constraints alone still leave physical measure choices.

Next: construct relational field observables and decide the physical status
of the singular-axis solutions. Potential, quantized matter, global cut
gluing and physical 3+1 gravity retain their separate obligations.

The native verifier computes D_N, its weighted Hermitian/positive split,
finite Gram changes of basis, forbidden-primary controls, unequal-weight
pairing defects, distribution-jet multiplication and measure/boundary
identities. The all-level proofs and distribution/measure arguments above
are written proofs; finite certificates do not formalize the general
theorems. Scientific arithmetic/rewrite work uses the frozen native JS
engines; Python 3.12 performs only I/O, hashing and orchestration.
