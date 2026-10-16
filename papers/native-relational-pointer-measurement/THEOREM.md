# RPM-1 — Relational pointer measurement of the physical phase

16 October 2026. Frozen parent: GOA-1 at
2c4a2bd07ea1a44a9a6b6a839b0e4a17ae14eb07. Units have hbar=1.
This packet constructs an **ideal apparatus extension of the already reduced
physical theory**. It does not derive an apparatus from the local LNC action.

## 1. Inherited system and explicit new inputs

Use the LPS/CCS/COS filtered physical option

    H = H_P direct_sum H_K direct_sum C_0,
    H_P=L^2(R,dP/|P|), H_K=L^2(R,dK/|K|).

The separate regular and origin rigging maps and their selected positive
relative normalization are unchanged. The generalized all-mode constraint
solutions and their physical norms are inherited, not rederived from a
finite Fock truncation. The original vacuum reduced Hamiltonian remains zero.

Let F be GOA-T4's self-adjoint P-branch operator, extended by zero on H_K
and C_0. In the flat P representation J chi=chi/sqrt(|P|),

    (F psi)(P)=P psi(P)+(1/8)[(2P-1)psi(P-1)+(2P+1)psi(P+1)],       (1)
    Dom F={chi: P chi in H_P} direct_sum H_K direct_sum C_0.

It quantizes the smooth weak Dirac observable P[1+(1/2)cos Q_sharp]
with the inherited full-line matching and symmetric reduced ordering.
These choices and their self-adjointness proof are frozen GOA inputs.

Add a pointer H_A=L^2(R,dX), X its multiplication operator and
Pi=-i partial_X on H^1(R). Add a separate ideal clock L^2(R,dT),
Pi_T=-i partial_T. It has the whole real momentum spectrum: no positivity
of clock energy or identification with the inherited mode grading is made.
Specify a real g in C_c^infinity(R), choose T_0 before its support and put

    G(T)=integral_(T_0)^T g(s) ds, kappa=integral_R g(s) ds.

There is no apparatus free Hamiltonian in this model. The clock, switching
profile, pointer preparation, coupling strength and readout are new inputs.
T is not identified with the old scalar or geometry clocks. Spatial locality,
physical apparatus construction and source-selected pulse dynamics are open.

## 2. A self-adjoint clock constraint and a positive physical rigging

Put H_SA=H tensor H_A. The commuting spectral measures of F tensor I and
I tensor Pi define the self-adjoint operator A with spectral function lambda p:

    Dom A={v: integral lambda^2 p^2 d<v,E_F(lambda) tensor E_Pi(p)v><infinity},
    U_s=exp(-i s A), W(T)=U_(G(T)).                              (2)

This joint spectral definition is the closure of the tensor product on
spectrally bounded tensor vectors; it is not a formal product of arbitrary
unbounded operators. Multiplication by W(T) is a unitary mathcal W on
K=L^2(R_T;H_SA).

**RPM-T1 (clock constraint and rigging).** Define

    C=mathcal W Pi_T mathcal W^*,
    Dom C=mathcal W H^1(R_T;H_SA).                              (3)

It is self-adjoint. On its dense core mathcal W(C_c^infinity(R) tensor D),
where D=union_n Ran 1_[-n,n](A), its differential expression is

    C=-i partial_T+g(T) A.                                    (4)

For the test space Phi=mathcal W(S(R) tensor_alg D), with the transported
Schwartz/inductive-limit test topology, the group-average form exists and is

    B(f,h)=integral_R ds <f,exp(i s C)h>_K
          =<Rf,Rh>_SA,    Rf=integral_R W(T)^* f(T) dT.          (5)

Consequently the quotient Phi/ker R has a positive definite inner product
and completion H_SA. Its generalized solutions and physical norm are

    Psi_v(T)=W(T)v,  C Psi_v=0,
    <Psi_v,Psi_w>_phys=<v,w>_SA
                         =<Psi_v(tau),Psi_w(tau)>_SA.           (6)

Proof. Self-adjointness and the core in (3) follow by unitary conjugation
of Pi_T. D is dense and every vector in it lies in all powers of Dom A.
On this core W'=-i g A W; differentiating W^* gives precisely the plus sign
in (4). Formula (3), rather than a presumed domain for a sum, defines the
full operator even where the two separate terms of (4) need not exist.

For f=mathcal W a and h=mathcal W b, exp(i s Pi_T)b(T)=b(T+s).
Their finite sums of Schwartz functions with Hilbert coefficients make the
double integral absolutely integrable. Fubini and u=T+s give (5). With the
unitary Fourier transform this is 2 pi <hat a(0),hat b(0)>;
the normalization in (5) is fixed, not an unspecified delta-function factor.
R maps onto D by choosing a Schwartz scalar of integral one. Hence its
null quotient is isometric to dense D and completes to H_SA. Integration
of a derivative proves R C f=0. Pairing a test f with Psi_v gives <Rf,v>,
so these are continuous generalized solutions. Transforming by W^* reduces
their equation to partial_T v(T)=0. In the class of H_SA-valued weak
distribution solutions the general solution is a constant Hilbert vector;
the physical completion selects exactly (6). No claim about every possible
larger algebraic dual is needed. Unitarity gives the last equality of (6).
There is no nonzero solution in K: a nonzero constant-norm history has
infinite integral over T. This does not make its physical norm infinite.

Thus positivity here is obtained from the new constraint's rigging form,
using the previously derived gravitational physical norm as an input. It
is not inferred from positivity of an unreduced kinematic Fock space.

## 3. Relational observables, domains and the pointer shift

**RPM-T2 (relational readout).** Every bounded instantaneous SA observable O
at clock value tau has physical representative B_tau=W(tau)^* O W(tau).
Its decomposable history lift

    mathcal O_tau(T)=W(T) B_tau W(T)^*                         (7)

strongly commutes with C: after mathcal W conjugation it is the constant
operator B_tau, which commutes with all clock translations. An unbounded
self-adjoint O has the corresponding self-adjoint physical conjugate on
W(tau)^* Dom O; its spectral projections have the bounded lifts (7).
This is an actual Dirac readout for the added clock constraint, not a
kinematic projection onto the nonnormalizable event T=tau.

After the pulse, the relational pointer is

    X_out=U_kappa^*(I tensor X)U_kappa,
    X_out=X+kappa F,        F_out=F.                          (8)

The domain of X_out is U_kappa^* Dom(I tensor X). Equation (8) initially
holds on the dense core given by finite sums v tensor phi with v in a bounded
spectral subspace of F and phi in S(R_X); its closure is the displayed
self-adjoint conjugate. One must not replace that domain by an unproved
domain intersection. The identity also holds on Dom X intersect Dom(F tensor I),
which suffices for the moment calculations below.

Proof. In the joint F spectral/X coordinate representation,

    (U_kappa z)(lambda,X)=z(lambda,X-kappa lambda).            (9)

Conjugating multiplication by X gives (8). Cut off the F spectrum to
[-n,n] and approximate the pointer in its multiplication graph norm by
Schwartz functions. Translation by the bounded shift kappa lambda preserves
this graph norm equivalently on each such subspace. These cutoffs commute
with U_kappa and X and converge in the conjugate graph norm, proving the
core assertion. The tensor F spectral projections commute with A and U,
which proves F_out=F including its spectral calculus. These arguments also
prove (7) without using a merely formal commutator of unbounded operators.

## 4. Normalized instrument and explicit disturbance

Choose the normalized real Gaussian pointer

    phi_sigma(X)=(2 pi sigma^2)^(-1/4) exp[-X^2/(4 sigma^2)], sigma>0.

**RPM-T3 (physical instrument).** The premeasurement isometry and outcome
operators on the original physical H are

    V_kappa psi=U_kappa(psi tensor phi_sigma),
    (V_kappa psi)(x)=K_x psi,   K_x=phi_sigma(x-kappa F),
    integral_R K_x^* K_x dx=I.                               (10)

For each Borel set B,

    E(B)=integral_B K_x^* K_x dx,
    I_B(rho)=integral_B K_x rho K_x^* dx                       (11)

is a normalized positive-operator-valued measure and a completely positive
instrument on trace-class states, with probability Tr[rho E(B)]. If that
probability is nonzero, I_B(rho)/Tr[rho E(B)] is a physical conditional state.
For a normalized psi with spectral measure mu_psi of F,

    p_psi(x)=integral |phi_sigma(x-kappa lambda)|^2 dmu_psi(lambda). (12)

For psi in Dom F, its readout has

    mean X=kappa <F>,
    Var X=sigma^2+kappa^2 (||F psi||^2-<F>^2).                (13)

Proof. Bounded functional calculus defines every K_x even for states outside
Dom F. Formula (9) gives (10). Tonelli's theorem and normalized translates
of |phi_sigma|^2 prove the weak integral identity and V^*V=I. The isometric
dilation proves complete positivity and countable additivity; the integrals
in (11) converge in trace norm for positive trace-class states and then by
linearity. Spectral calculus gives (12). Integrating the first two Gaussian
moments gives (13). Here the second spectral moment is ||F psi||^2: it does
not require applying F twice to an arbitrary Dom F vector.

Discarding the pointer gives a trace-preserving channel. In the spectral
representation of F, its multiplier on the (lambda,mu) density kernel is

    D(lambda,mu)=exp[-kappa^2(lambda-mu)^2/(8 sigma^2)].        (14)

This follows by completing the square in integral phi_sigma(x-kappa lambda)
phi_sigma(x-kappa mu) dx. Kernel notation can equivalently be read as the
double spectral integral; trace-class operators are Hilbert-Schmidt, so
there is no assumption of normalizable F eigenvectors. The spectral
probabilities of F are unchanged, but its distinct spectral coherences are
attenuated. This is nondemolition of F, not nondisturbance of all observables.

On H_K and the origin, K_x=phi_sigma(x)I. All K_x commute with the three
branch projections. Thus this instrument and its dilation preserve the
inherited physical sectors and their selected norms, with no inter-branch
coupling or new identification of their two rigging scales. Preservation is
on the completed physical quotient; no strong kinematic L_m commutator or
invariance of the original smooth test domain is asserted.

## 5. Exact phase readout with a finite-variance discrimination bound

Use the frozen GOA normalized flat-momentum packets

    psi_+=2 1_[1/4,1/2], psi_-=2 1_[-3/4,-1/2],
    Psi_s=(3/5)psi_+ + s(4/5)psi_-, s in {+1,-1}.             (15)

The physical packet is J^-1 Psi_s=sqrt(|P|)Psi_s in H_P. These are in
Dom F; their jump discontinuities do not place them in Dom Q, nor is that
required here. Their restrictions to GOA's entire minimal affine algebra
coincide. Their F means are -59/200 and -47/200 respectively.

**RPM-T4 (readout statistics).** The exact new second moments and variances are

| State | ||F Psi_s||^2 | Var F | Mean pointer |
|---|---|---|---|
| s=+1 | 1233/3200 | 23863/80000 | -59 kappa/200 |
| s=-1 | 1169/3200 | 24807/80000 | -47 kappa/200 |

Therefore the pointer mean contrast is -3 kappa/50. For kappa!=0 the two
normalized outcome distributions differ. This is a conditional prediction
of the selected apparatus model, not an independently calibrated gravity
prediction.

Proof. Apply (1) before integrating. Each input interval has length 1/4;
its output is supported on that interval and its two translates by one.
Exact integration of the linear output polynomials gives

    <psi_+,F psi_+>=3/8, <psi_-,F psi_->=-5/8,
    <psi_-,F psi_+>=-1/32,
    ||F psi_+||^2=25/128, ||F psi_-||^2=61/128,
    <F psi_+,F psi_->=1/96.                                 (16)

For the last identity only the two central intervals overlap: their
integrands are P(2P-1)/2 on [1/4,1/2] and P(2P+1)/2 on [-3/4,-1/2].
Integration gives -1/96 and 1/48, hence 1/96. Squaring the coefficients in
(15) and including twice the signed cross term proves the table. Subtracting
the square of each mean gives its variance. Formula (13) proves the readout.

More explicitly its characteristic function is

    integral exp(itx)p_s(x)dx
       =exp(-sigma^2 t^2/2)<Psi_s,exp(i kappa t F)Psi_s>.      (17)

Its derivative difference at t=0 is -3 i kappa/50. Thus Gaussian smearing
has not erased this distinction. Different finite first moments alone also
prove the distributions differ; no pointwise eigenstate assumption is used.

For N independently prepared identical runs, classify the known pair by
which of its two means is closer to the sample mean (with the inequality
reversed if kappa<0). Chebyshev gives, separately under either hypothesis,

    Pr(error|s) <= min(1, (10000/(9N))
                         [sigma^2/kappa^2+Var_s(F)]).         (18)

The distance to the decision midpoint is 3|kappa|/100. Independence divides
the sample-mean variance by N; these two facts prove (18). At kappa=sigma=1,
the larger bound is 104807/(72N), so N=30000 gives a bound below 1/20.
This conservative existence bound assumes exact known preparations and
calibration; it is neither an optimal sample complexity nor laboratory
performance. At kappa=0 both readouts are the same Gaussian and the bound
with kappa in the denominator is inapplicable.

## 6. Reuse of CHD: coherent records versus classical readout

**RPM-T5 (physical record interface).** For any pointer Borel bin B let
P_B=I tensor 1_B(X). The binary coherent record map on H_SA is

    C_B z=|1> tensor P_B z+|0> tensor (I-P_B)z,
    W_B=C_B V_kappa,       W_B^* W_B=I_H.                     (19)

For every bounded physical system operator O, set J_B(O)=W_B O W_B^*.
Exactly as in CHD-T2/T3,

    J_B(O_1 O_2)=J_B(O_1)J_B(O_2), J_B(I)=W_B W_B^*,
    W_B^* J_B(O) W_B=O.                                    (20)

Proof. P_B and I-P_B are an orthogonal partition, so C_B^* C_B=I;
(10) then proves (19). Inserting W_B^*W_B=I proves (20). The identity of
the image algebra is the range projection, not the whole record identity.
For an unbounded self-adjoint O the corresponding operator on Ran W_B has
domain W_B Dom O. These formulas encode operations faithfully on the
coherent image; they do not assert independent experimental access to all
encoded operations. Restriction to GOA's three-block algebra respects each
block and leaves the origin normalization choice intact.

The pointer bin probability is <psi,E(B)psi>. Reading the bit and tracing
out the continuous pointer gives precisely I_B and I_(R\B) of (11), not
necessarily the square-root/Lueders instrument of the two effects. The
identity follows by tracing the pointer coordinate on the diagonal inside
each bin. Keeping the bit alone or classicalizing its branches does not
retain the full isometry inverse in (20). This is the existing CHD coherent
record distinction applied to a specified physical observable and its
continuous pointer instrument, not a new general measurement principle.

## 7. Closure and remaining source obligation

Closed under the stated inputs: a self-adjoint added clock constraint;
its nonzero generalized solutions and positive group-averaged norm;
Dirac pointer observables with domains; a normalized physical instrument;
exact phase-dependent readings, variances and a statistical bound; and the
CHD coherent-record interface. The clock norm proof is analytic. The native
certificates check finite exact identities, not continuum self-adjointness.

Still open: derive a localized apparatus/matter coupling and its backreaction
from an extended LNC action, recheck all local first-class constraints and
quantum anomaly/physical norms, or prove a controlled reduction to this
apparatus model. A weak classical observable does not by itself justify
adding its coupling to unreduced quantum constraints. The ideal clock and
global pulse above do not settle that issue. Source-selected dynamics,
inter-branch/origin transitions, general field operators, physical 3+1
recovery and empirical calibration remain open. No inherited proof changes.
