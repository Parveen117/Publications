# LMD-1 — Logarithmic material density and local flow

16 October 2026. Frozen parent: MRM-1 at
e274b8e9890e9ee087d29ba2fcb87b76e5beb715.

MRM found a source total charge but obstructed nonconstant canonical
material flows in MPQ's smooth Gaussian representation. This packet
constructs a different, explicitly declared reduced representation.
It preserves the classical separating chart as the starting algebra,
but changes its quantum configuration measure and defines the density
by Wick renormalization. No equivalence with the old MPQ quantum
representation or its pointwise inverse geometry is asserted.

The finite-flow construction uses the Feldman--Hajek equivalence
criterion [FH]. The density construction is the elementary L2 part
of Gaussian multiplicative chaos [GMC], with its regularization
identification. These established methods are credited; the result
here is their source-specific integration and exact scope.
References and hypotheses are in SOURCE_MAP.md.

Throughout h=hbar is in the chosen dimensionless canonical units,
0<h<1, and the material circle has length one. This is an explicit
L2 construction range, not a derived value of Planck's constant.
The two orientation sectors are copies with sign epsilon=+1,-1.

## 1. Which covariance can remove the leading local-current divergence?

MRM's classical source current is unchanged:

    rho=-d'+d s', R[v]=integral d(v'+v s'), {d,s}=delta.   (1)

Consider paired Fourier covariances lambda_n=c |n|^(-2p), n!=0,
with c>0 and a positive zero-mode variance. Put sigma_n=sqrt(lambda_n).
The Weyl cutoff's second-chaos coefficient, divided by 2 pi, is

    S_mn/(2 pi)=(i/2)v_(m-n)
       [n sigma_n/sigma_m-m sigma_m/sigma_n].            (2)

**LMD-T1 (necessary exponent for this vacuum-domain criterion).**
For fixed k!=0, m=n+k and n -> +infinity, the bracket in (2) tends
to (2p-1)k. If p!=1/2 and v is nonconstant, its squared
Hilbert--Schmidt norm grows at least linearly in the cutoff, exactly
as in MRM-T5. The unique power cancelling this leading same-sign
divergence is p=1/2.

Proof. The ratio sigma_n/sigma_(n+k)=(1+k/n)^p
is 1+pk/n+O(n^-2). Its inverse gives the opposite coefficient,
leaving (2p-1)k in the difference. A nonzero v_k supplies order N
nonvanishing entries. MRM's orthogonal-chaos argument applies.
For p=1/2 the same-sign bracket is exactly zero, not only asymptotically.

This is a criterion for the stated Weyl vacuum domain. It alone
does not prove measure equivalence or exclude every alternative
domain for other exponents. Those distinctions from MRM are retained.

## 2. An alternative Gaussian representation with actual finite flows

Select instead

    lambda_0=1, lambda_n=1/(2|n|), n!=0,
    C(z,y)=K(z-y)=1-log(2|sin(pi(z-y))|), z!=y.           (3)

The real sine/cosine basis is orthonormal. Let mu_h be the Gaussian
distribution with independent coefficients of variance h lambda_j.
It is a countably additive probability measure on E=H^(-delta)(S1)
for any fixed delta>0, since sum_n n^(-1-2delta)<infinity.
Its Cameron--Martin space is H_C=H^(1/2), including constants.
Set H_epsilon^log=L2(mu_h), H_red^log=H_+^log direct_sum H_-^log.

For smooth real f, S(f)=<s,f> has its maximal multiplication domain.
For smooth real k the usual Cameron--Martin unitaries define
self-adjoint D(k), with

    D(k)=i h partial_k-(i/2) beta_k,
    beta_k=<C^{-1}k,s>, [D(k),S(f)]=i h integral k f.    (4)

Their Stone domains specify these operators. This positive norm is a
selected reduced quantization, not a consequence of anomaly cancellation.

For a real constant a and smooth increasing circle diffeomorphism g put

    F_g^a(s)=s composed with g+a log g'.                 (5)

Pullback of a distribution is defined by duality:
<s composed with g,f>=<s,(f/g') composed with g^{-1}>.
Our convention has F_g F_b=F_(b composed with g).

**LMD-T2 (quasi-invariance and local Stone generators).**
(F_g^a)_*mu_h is equivalent to mu_h for every smooth g and every real a.
The positive half-density operators

    rho_g=d[(F_g^a)^(-1)_*mu_h]/dmu_h,
    T_g^a psi(s)=sqrt(rho_g(s)) psi(F_g^a(s))            (6)

form a genuine, strongly continuous unitary representation, with
T_g T_b=T_(g composed with b). For the flow g_v(t) of +v, define

    T_(g_v(t))^a=exp[-it R_a[v]/h],
    Dom R_a[v]={psi:lim_(t->0) i h(T_(g_v(t))^a psi-psi)/t
                              exists in L2(mu_h)}.     (7)

Each real smooth v has a self-adjoint R_a[v]. This is a construction
of the nonconstant generators, not just a formal commutator.

Proof. The covariance change of s composed with g is the kernel

    K(g(z)-g(y))-K(z-y)
       =-log | sin(pi(g(z)-g(y)))/sin(pi(z-y)) |.         (8)

The positive ratio extends smoothly across the circle diagonal with
value g'(z); elsewhere it is smooth and nonzero. Thus (8) is a
smooth kernel on the compact two-torus. Whitening it on both sides
by C^(-1/2), an operator of order 1/2, leaves a smoothing, hence
Hilbert--Schmidt, operator. Composition with g and g^{-1} preserves
H^(1/2): use its Gagliardo norm, bounded Jacobians and two-sided
bounds on chord distances. The Cameron--Martin spaces therefore
agree. The shift a log g' is smooth and lies in that space.

These are precisely [FH]'s range, mean and Hilbert--Schmidt
conditions. To use its Hilbert-space statement, take E=H^(-delta);
the Fourier weights used to identify E with l2 cancel in the
whitened covariance comparison. The covariance is trace class on E.
Thus the two actual Gaussian measures are equivalent.

The positive Radon--Nikodym chain rule proves the exact group law
and norm preservation in (6), without any projective phase. For
g -> id in C^infinity, (8) tends to zero in all smooth kernel norms,
the whitened change tends to zero in Hilbert--Schmidt norm, and the
mean tends to zero in H_C. The Gaussian Hellinger formula then gives
sqrt(rho_g) -> 1 in L2: its logarithmic covariance factors are
O(sum_j (r_j-1)^2) for relative eigenvalues r_j near 1, and the mean
term tends to zero. On bounded continuous cylinder tests F_g s -> s;
dominated convergence proves strong continuity, then density gives
it on H. Stone's theorem supplies (7).

Rigid translations have rho_g=1 and log g'=0. Their total-charge
operator has the same angular expression and spectrum 2 pi h Z as
MRM-T2, now in a different Gaussian representation.

## 3. A positive renormalized density and its spectral domain

Write s_N=P_Ns and V_N=h(1+sum_(n=1)^N 1/n). Define positive measures

    M_N(s,dz)=exp[s_N(z)-V_N/2] dz.                     (9)

The subtraction is half the complete variance, including the zero
mode. Omitting only the zero-mode term multiplies the limit by
exp(h/2); omitting the divergent nonzero-mode subtraction does not
give this normalized limit.

**LMD-T3 (L2 density measure and observable).** For 0<h<1, M_N
converges weakly almost surely to a finite positive random measure M.
For every bounded measurable w, M_N(w) converges in L2 to M(w), with

    E M(w)=integral w,
    E[M(w)M(v)]=integral integral w(z)v(y)e^{h K(z-y)} dz dy. (10)

The measure is nonzero almost surely, has full support and no atoms.
For real bounded w,

    J_h(w)=epsilon M(w),
    Dom J_h(w)={psi:M(w)psi in L2(mu_h)}                 (11)

is self-adjoint. The vacuum belongs to this domain. Smoothing the
field with a smooth nonnegative local approximate identity gives
the same measure when its own full variance is subtracted.

Proof. Independent Fourier shells make M_N(w) a positive martingale
for w>=0. Summation by parts gives

    K_N(r)<=C+log(1/dist(r,Z)), r not in Z,

uniformly in N: split the harmonic cosine sum at 1/dist(r,Z),
and bound the rest by the geometric-partial-sum estimate.
Therefore e^{h K_N(r)} is bounded by an integrable constant times
dist(r,Z)^(-h). Gaussian exponentiation gives (10) at finite N;
dominated convergence gives its limit and L2 boundedness. Martingale
convergence on a countable dense set of continuous tests gives a
finite positive weak limit; the same L2 bounds extend to bounded w.

For any nonempty arc I, E M(I)=|I|>0. The event M(I)=0 is unaffected
by any finite change of Fourier coefficients, which multiplies the
limiting measure by a strictly positive continuous function. It is
a tail event and has probability 0 or 1; its positive mean and L2
integrability exclude 1. A countable arc basis proves full support.
Moreover E M(I)^2<=C_h |I|^(2-h). For dyadic intervals the sum of
these bounds is O(2^(-n(1-h))). Markov and Borel--Cantelli imply
the largest interval mass tends to zero, proving absence of atoms.

The bounded-test joint second-moment argument also identifies
Fourier and mollifier limits; their cross covariances converge off
the diagonal to K and have the same integrable logarithmic bound.
This is the L2 version of the regularization uniqueness in [GMC].
For a diffeomorphically transported mollifier the same estimates
hold because g and g^{-1} are smooth and bi-Lipschitz. The multiplication
spectral theorem now proves (11).

No pointwise density function exp(s(z)) has been constructed. The
field is a distribution and the output is a measure. The old
pointwise inverse r,m from MPQ/RMG is not supplied by (9).

## 4. The exact density law fixes a quantum correction

Use g^*M for measure pullback: (g^*M)(B)=M(g(B)). Its smooth-density
version is M(g(z))g'(z) dz.

**LMD-T4 (renormalization weight and its necessary correction).**
For every fixed smooth g and real a, almost surely,

    M(F_g^a s,dz)=(g'(z))^(a-1-h/2) g^*M(s,dz).         (12)

Consequently the density law M(F_g s)=g^*M(s) holds for all g
in the family (5) exactly when

    a=a_h=1+h/2.                                       (13)

The corresponding full spectral-domain covariance is

    (T_g^a)^* J_h(w) T_g^a
       =J_h((w composed with g)(g')^(1+h/2-a)).          (14)

For a=1 the old classical configuration flow is implementable,
but its Wick density has the extra factor (g')^(h/2) in (14).
For a=a_h the factor disappears. This is an explicit quantum
modification of the local current, not retention of every old
observable identity by changing notation.

Proof. In a mollifier cutoff, the diagonal value of the smooth
covariance difference (8) gives

    Var[(s composed with g)_eta(z)]-Var[s_eta(z)]
                                   -> -h log g'(z).

Subtracting the fixed reference variance, rather than the transported
variance, therefore contributes (g')^(-h/2).
The intrinsic Wick exponential of s composed with g is the
pullback M with its ordinary coordinate Jacobian divided out:
(g')^(-1) g^*M. This follows by the regularization identification
in T3 and change of variables; smoothing before or after g has
the same intrinsic limiting measure. The deterministic mean in
(5) contributes (g')^a. Multiplication proves (12).
Quasi-invariance ensures that the same measurable M is defined at
F_gs; the identity is for each fixed g, and also its continuous
unitary action, not an unproved simultaneous samplewise statement
for an uncountable group.

Equation (14) follows by conjugating multiplication in (6) and
applying (12) to g^{-1}. If its weight were one for every g
while a!=a_h, choose g with nonconstant derivative; full support
of M detects a nonzero multiplier on an open arc. This proves
necessity within the declared affine-flow and Wick-density family.

For smooth k the Cameron--Martin displacement gives
M(s+t k)=e^{t k}M(s). Thus the inherited smeared affine law
[D(k),J_h(w)]=i h J_h(kw) survives on its differentiable domains.
The new correction is in the material-flow generator.

## 5. Local current action, domains and source feedback

Let A_v f=v f', b_v=v', and in a real Fourier basis write
B_v=C^(-1/2) A_v C^(1/2), S_v=(B_v+B_v^*)/2.
The unbounded skew part is allowed; S_v is Hilbert--Schmidt.
The linear functional b_C=C^(-1/2)b_v is in l2.

**LMD-T5 (finite vacuum action and coherent domain).** R_a[v]
acts on a finite polynomial cylinder in smooth smeared fields by
the half-density derivative

    R_a[v]=i h[X_v^a+(1/2)div_(mu_h) X_v^a],
    X_v^a s=v s'+a v'.                                 (15)

The divergence is defined by the L2 limit of the centered
finite Gaussian expressions, not an unrenormalized infinite trace.
For the vacuum, with standard Gaussian coordinates x,

    R_a[v]1=-(i sqrt(h)a/2)b_C.x
                     -(i h/2):x^T S_v x:,
    ||R_a[v]1||^2=(h a^2/4)||b_C||^2
                              +(h^2/2)||S_v||_HS^2.     (16)

For a real v with complex Fourier coefficients v_n,

    ||S_v||_HS^2/(2 pi)^2
       =sum_(n>=1) [(n^3-n)/3+n/2] |v_n|^2,
    ||b_C||^2/(2 pi)^2=4 sum_(n>=1) n^3 |v_n|^2.        (17)

All terms are finite for smooth v. Weyl cutoff actions converge
in L2 on smooth finite polynomials and on the coherent states below.
This is an action convergence statement; a general cutoff strong
resolvent theorem is not asserted here.

Proof. Formula (2) is exactly zero for same-sign nonzero indices.
For opposite signs m=-q,n=p, p,q>0, its squared modulus after
division by (2 pi)^2 is pq |v_(p+q)|^2. Summing both signs and
using sum_(p=1)^(n-1) p(n-p)=(n^3-n)/6 gives the first part of
(17). The zero-mode rows and columns contribute n|v_n|^2/2.
Whitening b gives the second formula. The first and second Gaussian
chaos terms in (16) are orthogonal, with the centered quadratic
norm equal to twice the Hilbert--Schmidt norm.

These are also the derivatives of the actual Gaussian half-density
unitary (6). Indeed the relative covariance and mean curves are
smooth, respectively in the Hilbert--Schmidt and Cameron--Martin
norms, by (8). Near identity their square-root densities are L2
differentiable: expand the finite Gaussian density, center its
quadratic part, and pass to the limit using those two norms.
The higher expansion remainders are controlled by Gaussian
quadratic exponential moments on a sufficiently small parameter
interval. This gives (16). A polynomial cylinder's pullback derivative is
a finite polynomial in smooth Gaussian linear functionals, so the
product rule gives (15). The same Gaussian bounds apply to its
higher finite moments. The covariance whitening in (17) makes
the cutoff quadratic terms Cauchy, and the cylinder derivative
tails converge in their Gaussian variances.

For smooth s_*,d_* use the MPQ-form normalized coherent state,
now with the new C:

    psi_*=exp[beta_(s_*)/(2h)-||s_*||_C^2/(4h)]
                            exp[-i<d_*,s>/h].           (18)

It is in Dom R_a[v], and

    <R_a[v]>=integral d_*(v s_*'+a v'),
    Var R_a[v]
      =h C((v d_*)',(v d_*)')
          +(h/4)||v s_*'+a v'||_C^2
          +(h^2/2)||S_v||_HS^2.                         (19)

To verify this, set s=s_*+sqrt(h) C^(1/2)x under
|psi_*|^2 mu_h. The real fluctuation is the linear Gaussian
<d_*,v(s-s_*)'>. The imaginary fluctuation is
-(sqrt(h)/2)<C^(-1/2)(v s_*'+a v'),x>
-(h/2):x^T S_v x:. Their squared norms give (19);
the linear/quadratic cross term has odd Gaussian degree.
Smooth shifts and the Hilbert--Schmidt bound justify differentiation
in L2 and the cutoff action convergence.

On this calculation domain, changing a gives the exact identity

    R_(a_h)[v]=R_1[v]+(h/2)D(v').                       (20)

The closed realization is always specified by its Stone group,
rather than by an unproved sum-domain assertion.

For a=a_h, M(S1) is invariant under (5). Its bounded spectral
projections therefore reduce every T_g. On the subspace
M(S1)<=L, each J_h(w) is bounded by L||w||_infinity and

    (i/h)[R_(a_h)[v],J_h(w)]=J_h(v w')                  (21)

on Dom R_(a_h)[v] in that subspace. To justify the domain and
commutator, differentiate (14): its right side is differentiable
in bounded-operator norm there, since w is smooth. These mass
subspaces increase to the whole Hilbert space, giving a dense
domain prescription. This avoids assuming all powers of GMC
observables lie in a common L2 algebra.

For the unchanged density pulse U=exp[-i kappa J_h(w) tensor Pi/h],

    U^*(R_(a_h)[v] tensor I)U
       =R_(a_h)[v] tensor I-kappa J_h(v w') tensor Pi.   (22)

Use the preceding mass subspaces and bounded pointer-momentum
subspaces to prove the differentiated identity, then take their
dense union. The self-adjoint realization of the right side has
domain U^*Dom(R_(a_h)[v] tensor I). This is exactly MRF's smeared
reference-momentum balance. For a=1 an additional
(h/2)J_h(v' w) appears in that balance; it cannot be discarded.

An ideal spectral pointer and free-clock quotient can now be
constructed for each R_(a_h)[v] by MRM-T3's general self-adjoint
recipe. This is a selected reduced apparatus, with the classical
source lift from MRM-T1. It does not establish an independently
ordered unreduced local composite.

## 6. Controlled classical recovery and retained limits

**LMD-T6 (classical smeared density and current).** In the coherent
state (18), for bounded real w,

    E J_h(w)=epsilon integral w e^{s_*},
    Var J_h(w)=integral integral w(z)w(y)e^{s_*(z)+s_*(y)}
                                  [e^{h K(z-y)}-1] dz dy. (23)

As h -> 0 these spectral readings concentrate in probability at
the classical density integral. Equation (19), with a_h -> 1,
likewise gives concentration at the original R[v] of (1).

Proof. Cameron--Martin translation multiplies the Wick measure by
e^{s_*}; its phase d_* changes no configuration law. Apply (10).
For 0<h<=h0<1, the integrand is dominated by a constant times
1+dist(z-y,Z)^(-h0), which is integrable. Dominated convergence
gives variance -> 0; Chebyshev gives concentration. In (19) the
smooth-center coefficients remain bounded as a_h -> 1, while
the displayed powers of h vanish.

This is an alternative continuum reduced quantization with actual
local-flow generators, a positive smeared density, their physical
domains and a classical comparison. It is not a unitary continuation
of MPQ: the old measure is supported on H^3 functions, while the
new Gaussian series is not even L2 almost surely (its positive
Fourier-square sum diverges). Their configuration measure classes
are singular. A norm-preserving scalar-space isomorphism does not
intertwine those configuration spectral algebras.

The old selected gauge functional ell can still be tensored with
this declared new reduced Hilbert space to give the positive form
conj(ell(psi)) ell(phi)<u,v> and its quotient H_red^log. That elementary
factorization does not import Fable's unmerged completeness claims,
nor prove an equivalence to frozen RGR's MPQ physical observable
contract. Reconstruction of the original pointwise r,m,rho and
comparison with an unreduced nonlinear ordering remain open.

## 7. A concrete three-spatial-dimensional diagnostic

The one-dimensional cancellation in (8) must not be extrapolated
to arbitrary spatial dimension. A particularly explicit test
already fails on a three-torus.

**LMD-T7 (3D scalar logarithmic shear obstruction).** On T^3 choose
an isotropic scalar Gaussian Fourier covariance
lambda_k=|k|^(-3), k!=0, with any positive constant-mode variance.
Let the volume-preserving torus diffeomorphism be g(x)=A x modulo 1,

    A = [ [1,-2,-2], [0,1,0], [0,0,1] ], det A=1.      (24)

The laws of s and s composed with g are mutually singular. Adding
a log-Jacobian term changes nothing for this g.

Proof. Fourier pullback replaces the coefficient at k by the
old coefficient at A^(-T)k. For k=n e_1,

    A^(-T)e_1=(1,2,2), |A^(-T)e_1|=3,
    lambda_(A^(-T)n e_1)/lambda_(n e_1)=1/27.           (25)

The normalized real cosine coefficients along this ray are
independent centered Gaussians. Their empirical mean square
converges almost surely to 1 under the original normalized law,
and to 1/27 under its transformed law. These disjoint full-measure
sets prove singularity, exactly as a path-statistic argument,
not an extrapolation from a finite matrix. Since log det A=0,
no coefficient of a log-Jacobian shift repairs it.

This is a genuine three-spatial-dimensional test of the proposed
scalar representation, not construction or refutation of 3+1
gravity. It does not exclude tensor-valued/dynamical metric
measures, a different scalar measure class or a smaller symmetry
contract. The current research must supply an actual 3+1 source
and its gravitational constraints before claiming that target.

## 8. Evidence and next obligation

The analytic results use the explicit Gaussian and Wick contracts,
[FH], Stone's theorem, standard Gaussian moment/martingale facts
and the stated [GMC] regularization identification. Native tests
verify covariance weights, exact same/opposite-sign cancellations,
the Hilbert--Schmidt sum, the quantum density-weight coefficient,
coherent variance coefficients, feedback signs and the 3D shear.

No finite check proves quasi-invariance, the continuum measure,
operator domains or the three-dimensional singularity theorem.
The construction is a source application of established methods,
not a claim to have discovered Gaussian multiplicative chaos.

Next: specify a geometric observable recoverable from this
measure-valued density and distributional conjugate, or derive
an explicit tensor/source representation for the 3+1 bridge.
The old pointwise inverse cannot be reused without that work.
Fable's separate core-completeness result and its remaining
density audit are recorded in FABLE_REVIEW.md.
