# HQS-1 — A positive quantum sector of the inherited homogeneous geometry

16 October 2026. This is a quantization of the specified regular homogeneous
zero-potential, matter-free LNC/REC sector. It is not a quantization of all
local LNC constraints or of physical 3+1 gravity. General operator-domain
and averaging statements have written proofs below; finite native checks
have their separate, narrower scope.

## 1. Source action and decisions that must remain visible

REC-1 already derives the homogeneous action

    S = integral [pi L_dot + k phi_dot - N C] ds,
    C = L pi k + L^2 V(phi),       L>0.

Keep its periodic independent-coframe boundary choice, homogeneous lapse,
zero shift and zero matter. Select the diagnostic V=0 already identified
in the research plan. This does not select V=0 as nature's potential.
Write q=log(L/L_ref), p=L pi with a fixed positive reference L_ref. Then

    pi dL = p dq,       C=p k,       {q,p}={phi,k}=1.       (1)

The reference sets a coordinate scale, not a measured fundamental length.
All quantum formulas below use a declared dimensionless action unit hbar=1;
no physical value of Planck's constant has been derived.

The regular classical clock branch has p!=0, phi_dot=Np, hence k=0.
Its reduced action is integral p dq and its phi-Hamiltonian is zero. Both
signs of p are retained in the quantum construction. We do not demand that
every geometry observable preserve a single sign branch. The singular
p=0 constraint sector, including solutions with k!=0, is excluded by the
regular-clock test-space choice. It is not declared nonexistent.

Quantization choices: the canonical log-length representation; product
ordering C_hat=p_hat k_hat; regular-sector group averaging; and a full-line
Fourier representation for the reduced geometry observables, including
its specified matching through p=0. These are explicit inputs. Positivity,
self-adjointness, the averaged measure and consistency statements are
derived for those inputs, not asserted unique among all quantizations.
The global fixed-winding single-angle restriction from REC remains separate.

## 2. HQS-T1 — Self-adjoint constraint and unitary gauge group

Let H_kin=L^2(R_q x R_phi,dq dphi). Use unitary Fourier transform with
kernel exp[-i(p q + kappa phi)] and its usual Plancherel normalization.
The standard self-adjoint momenta are p_hat=-i partial_q and
k_hat=-i partial_phi. In Fourier coordinates they multiply by p,kappa.
Define

    (C_hat f)^(p,kappa)=p kappa f^(p,kappa),
    Dom C_hat={f: integral p^2 kappa^2 |f^|^2 < infinity}.  (2)

This maximal real multiplication operator is self-adjoint. Explicitly,
the adjoint test against compactly supported functions forces its adjoint
to be multiplication by the same real function on exactly the domain (2).
Equivalently (C_hat +/- i)^(-1) is bounded multiplication by
1/(p kappa +/- i), giving both full resolvent ranges. Smooth compactly
supported Fourier functions are a core: truncate in (p,kappa), then
approximate smoothly on compact sets, in the graph norm of multiplication.

The product is unambiguous because the two momenta strongly commute: their
spectral projections multiply in different coordinates. Thus it is the
self-adjoint joint-spectral product, not an arbitrary product on an
unspecified domain. It agrees with -partial_q partial_phi on the common
smooth core. The constraint group is

    (U(s)f)^(p,kappa)=exp(i s p kappa) f^(p,kappa).         (3)

Unit modulus and dominated convergence prove unitarity and strong
continuity. Its generator is (2). This is gauge-parameter transport, not
a positive energy Hamiltonian or a physical faster-than-light evolution.
A single constraint has a trivial self-commutator; this supplies no result
about central terms in the inhomogeneous constraint algebra.

## 3. HQS-T2 — Positive physical inner product on the regular clock sector

Choose the dense invariant test space in Fourier coordinates

    Phi=C_c^infinity((R\{0})_p x R_kappa).

Its elements have compact p-support bounded away from zero. It is dense
in H_kin, since the omitted line has Lebesgue measure zero, and is a core
for (2) by the same truncation argument. It deliberately excludes the
distributional p=0 solution sector from this averaging prescription.

Define the sesquilinear average, anti-linear in f, by

    B(f,g)=integral_R ds/(2 pi_math) <f,U(s)g>_kin.

For this test space the integral is well-defined: repeated integration by
parts in kappa makes the integrand decay in s faster than any fixed
inverse power, uniformly on p-support away from zero. Using Fourier
inversion in kappa (or a Gaussian-regularized delta limit) gives

    B(f,g)=integral_(p!=0) dp/|p| f^(p,0)^* g^(p,0).     (4)

For the regularized derivation, the s-integral kernel is
(4 pi_math epsilon)^(-1/2) exp[-(p kappa)^2/(4 epsilon)].
Set y=p kappa at fixed p. The absolute Jacobian is 1/|p|, not 1/p.
The resulting approximate identity converges uniformly under the compact
support and smoothness bounds. This proves (4) without an uncontrolled
exchange at p=0.

The form is positive semidefinite and

    B(f,f)=0 iff f^(p,0)=0 almost everywhere.

Quotient by this null space and complete. The resulting physical Hilbert
space is

    H_phys = L^2(R\{0}, dp/|p|),
    eta(f) represented by chi_f(p)=f^(p,0).              (5)

Every compactly supported smooth chi away from zero is attained: take
f^(p,kappa)=chi(p)b(kappa), with b smooth compact and b(0)=1. Those chi
are dense in (5). Hence (5) is precisely the completion, with no residual
null vectors. Also B(f,C_hat g)=0 and eta(U(s)g)=eta(g), because evaluating
at kappa=0 kills the constraint factor and the gauge phase. The generalized
solution functional has density delta(kappa) chi_g(p)/|p| on test inputs.

This is a positive physical inner product for the chosen sector. It does
not assert a positive spectrum for C_hat; (2) has both signs, as expected
for this constraint. The averaging construction is a source-specific use
of established refined algebraic quantization, explicitly credited.

## 4. HQS-T3 — Equivalence with reduced quantization and positive geometry

The map

    (J chi)(p)=chi(p)/sqrt(|p|)

is unitary from (5) to L^2(R,dp), where the value at zero is irrelevant.
Let F_q be the unitary Fourier transform on R and define

    A = F_q^(-1) J : H_phys -> L^2(R,dq).                (6)

This is an explicit unitary equivalence with quantization of REC's regular
V=0 reduced pair (q,p). Physical observables are transported from the
reduced representation through A; we do not assert that a naive
kinematical length operator automatically intertwines the rigging map.

On L^2(dq) define

    Q psi=q psi,                Dom Q={psi:q psi in L^2},
    P psi=-i psi',              Dom P=H^1(R),
    L psi=L_ref exp(q) psi,     Dom L={psi:exp(q)psi in L^2}. (7)

Q and L are real maximal multiplication operators and self-adjoint by the
same argument as (2); P is the inverse Fourier image of real multiplication
by p, hence self-adjoint. L is positive, has spectrum [0,infinity), and has
no zero eigenvector since exp(q)>0 almost everywhere. There is no positive
minimum length inferred from this: zero is the spectral infimum.

Compactly supported smooth functions in q form a common invariant domain
for L,P and their displayed algebraic products. Direct differentiation gives

    [L,P]=i L,      [Q,P]=i I.                            (8)

The complete unitary dilation group D(s)=exp(-i s P) acts by
D(s)psi(q)=psi(q-s), and

    D(s)^* L D(s)=exp(s) L.                              (9)

This follows by substitution on the multiplication domain, which is
mapped into itself with the appropriate scalar factor. It describes a
geometry symmetry, not time evolution under the zero reduced Hamiltonian.
Conjugating (7)-(9) by A gives self-adjoint positive geometry on H_phys.

An important weight correction is visible in the physical p coordinates:

    Q_phys chi = i[chi' - chi/(2p)]             (p!=0),
    P_phys chi = p chi.                                  (10)

The domain of Q_phys is J^(-1)H^1(R). Thus chi/sqrt(|p|) has the standard
full-line H^1 matching across zero. This matching is part of the declared
reduced representation, not selected by the constraint alone. The naive
i partial_p is not symmetric for dp/|p|: integration by parts gives the
missing logarithmic weight term. Ignoring it would invalidate the
observable construction. For L_phys one uses the unitary functional
calculus A^(-1)L A, not a formal complex shift in momentum.

The classical q,L are weak regular-sector Dirac observables:
{q,C}=k=C/p, {L,C}=Lk, while {p,C}=0. The kinematical differential
expressions satisfy [C,Q]=-i k and [C,L]=-i L k on a common suitable
core. They do not strongly commute with C on all kinematical states.
Preservation of a regular kappa=0 solution alone supplies neither physical
symmetry nor the correct rigging-map intertwiner. Equation (6), including
the weight correction and domain, specifies the physical observable action.

Both p signs are included. The chosen full geometry algebra is not claimed
to preserve a superselected p>0 or p<0 subspace. A sharp clock orientation
restriction would require a separate observable/domain analysis.
At a relational clock reading tau the reduced V=0 evolution is exactly
the identity, so L and P are tau independent. This reproduces REC's
stationary geometry with a running scalar clock. It is not a matter or
gravitational-wave dynamics result.

## 5. HQS-T4 — Why ordinary half-line momentum cannot be assumed

The unitary map from L^2(dq) to L^2((0,infinity),dL) is

    (B psi)(L)=L^(-1/2) psi(log(L/L_ref)).                (11)

It sends length to multiplication by L and dilation to

    B P B^(-1)=-i(L partial_L + 1/2).

The latter is self-adjoint on the transported domain. By contrast the
ordinary canonical momentum Pi_0=-i partial_L on C_c^infinity(0,infinity)
is symmetric but has deficiency indices (1,0). Indeed its adjoint is
the distributional derivative on the absolutely continuous L^2 functions
with L^2 derivative. Solving (Pi_0^*-i)u=0 gives u=c exp(-L), an L^2
solution; solving (Pi_0^*+i)u=0 gives u=c exp(L), not in L^2 unless c=0.
The deficiency-index criterion therefore forbids a self-adjoint extension
of this operator on this scalar half-line Hilbert space.

This does not prohibit alternative enlarged theories or different physical
observables. It establishes the domain obstruction for this proposed
ordinary momentum and motivates the valid affine pair (L,P), where
P represents L pi rather than pi. The construction is a known affine
quantization method applied to the inherited sector, not a new general
quantization theorem or a singularity-resolution claim.

Nor can an exact finite positive matrix carrier represent (8) with both
operators Hermitian and L strictly positive: trace[L,P]=0, but trace(iL)
is nonzero. Finite CHD matrix tests are compressed observable protocols;
they are not a finite-dimensional realization of the full affine algebra.

## 6. HQS-T5 — The completed observer theorem applies on this physical space

Let E:C^n -> H_phys be any specified isometric finite preparation map.
For any Hermitian n-by-n H, the finite-rank physical observable

    H_phys,E = E H E^*

is bounded and self-adjoint, and its readings at Ep are p^*Hp. CHD-3's
n^2 circular probes therefore reconstruct it on the known image corner.
For any declared finite cut tree W, E_record=(I_labels tensor E)W is an
isometry into the corresponding physical record space and has the same
coherent decoder. These are statements about actual normalizable physical
states of (5), not merely kinematical matrix notation.

To realize a finite cut P on the whole physical carrier, use the projection
E P E^* and its full complement I_phys-E P E^*. On prepared inputs Ex these
act as E P x and E(I-P)x. Thus the finite history is the restriction of
genuine orthogonal physical cuts, not an incomplete partition of H_phys.

One explicit preparation family uses disjoint momentum intervals in the
unitarily equivalent reduced momentum space L^2(dp). The two unit vectors
are a_+(p)=1_[1,2](p) and a_-(p)=1_[-2,-1](p). Their physical versions are
chi_+/-=sqrt(|p|)a_+/-. They are orthonormal by (4). The full constraint
is solved in the averaged sense; its singular p=0 sector is not used.

For c=(3/5,4/5), the state has unit physical norm and sign probabilities
9/25 and 16/25. Its exact dilation-momentum moments are

    <P> = -21/50,       <P^2> = 7/3,
    Var(P)=16177/7500.                                  (12)

The packet verifier integrates these polynomials using exact native
rational antiderivatives. Box packets belong to the physical Hilbert
completion and the domains of momentum powers, not necessarily the length
operator domain; no finite length moment for them is asserted. Length's
domain and commutator proof use the separate common q-space core above.

The finite packet span is not invariant under P: its compressed first
moment is diag(3/2,-3/2), but its compressed second moment is (7/3)I,
not the square (9/4)I of that first compression. The difference is (1/12)I.
In fact {p^j 1_[1,2]: j>=0} is linearly independent, since a polynomial
vanishing almost everywhere on an interval is identically zero. Thus
momentum operation closure of this preparation record is infinite. CHD's
finite observer interface certifies its declared bounded corner observables,
not every future continuum geometry operation from two amplitudes alone.

The finite source observable [[2,1],[1,5]] can be embedded on this packet
span, yielding exactly the inherited CHD circular reconstruction and
coherent history test. It is a bounded diagnostic observable, not asserted
to equal length, curvature or a matter coupling. The cuts and preparations
remain specified operations; deriving them from an interacting gravitational
model is still separate.

## 7. Certification and closure verdict

Closed for the stated regular homogeneous V=0 matter-free quantization:

- a self-adjoint inherited constraint and complete unitary gauge group;
- an explicitly derived positive group-averaged physical inner product;
- unitary equivalence to the reduced geometry representation;
- self-adjoint positive length and self-adjoint dilation, with domains;
- an explicit ordinary-momentum domain obstruction and a valid replacement;
- a concrete physical Hilbert-space interface for the completed cut/observer
  certificate, including normalizable packet examples.

The continuum arguments use standard Fourier analysis, real multiplication
operators, group averaging and the deficiency-index theorem, with their
applications proved here. Native exact computations separately certify
source canonical identities, measure Jacobians, packet moments, finite
history/observer identities, affine word identities and rejected stronger
claims. They do not formalize the general analytic proofs in a proof assistant.

Open: the singular clock sector and alternative global matching choices;
nonzero V and matter; quantization of the inhomogeneous local constraints,
central terms and domains; the independent global cut-gluing contract;
physical 3+1 recovery and calibrated predictions. The vanishing single
constraint self-commutator cannot close any of these obligations.

This is a positive quantum completion of one explicitly selected geometry
sector. It advances the quantum-gravity construction beyond a classical
homogeneous reduction and beyond an abstract observer theorem. It is not a
complete physical theory of quantum gravity, nor evidence for FTL.
