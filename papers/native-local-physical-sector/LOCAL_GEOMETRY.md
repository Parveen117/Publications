# LPS-2 — Local regular reduction and relational length

This follows LPS-1 without changing the frozen LNC, REC, HQS or LQA packets.
It supplies the classical geometric interpretation missing from an abstract
transport of HQS observables. It does not define inverse quantum currents
or a pointwise exponential of a quantum field.

## 1. The regular patch and exact constraint solutions

Keep the period-one, V=0, matter-free source. Fix a length reference L_ref
and use q=log(e/L_ref). This is only a constant shift of LQA's q=log e;
the constraints depend on q' and are unchanged. Set

    A_+=p+phi', A_-=p-phi', B_+=k+q', B_-=k-q',
    T_+=A_+ B_+/2-A_+'/2,
    T_-=A_- B_-/2+A_-'/2.

Select the open classical patch of smooth periodic real fields for which
both A_+ and A_- are everywhere nonzero. Each has a constant sign on the
circle. Periodicity gives integral A_+=integral A_-=P=integral p, so
their signs agree and P!=0. This is a genuine restriction: P!=0 alone
does not imply pointwise nonvanishing of these currents.

**LPS-T5 (regular local invariant).** On this patch the constraints are
equivalent to

    B_+=A_+'/A_+,       B_-=-A_-'/A_-,
    q'=(log|A_+A_-|)'/2,
    k=(log|A_+/A_-|)'/2.                              (1)

Consequently

    C=q-(1/2)log|A_+A_-|                              (2)

is independent of x, and K=integral k=0. Conversely, choose any two smooth
nowhere-zero periodic densities with the same integral P and same sign,
choose C and phi_0, and put

    p=(A_++A_-)/2,  phi'=(A_+-A_-)/2,
    q=C+(1/2)log|A_+A_-|,
    k=(1/2)(log|A_+/A_-|)'.                           (3)

The zero integral of phi' supplies a periodic phi with arbitrary additive
constant. Substitution gives both original constraints identically zero.
This is a full parametrization of their smooth solutions on this patch.

Proof: divide T_+=T_-=0 by the nonzero currents, add/subtract, and integrate
the derivatives of periodic log absolute values. No Fourier truncation
or assumption that all fields were homogeneous is used.

## 2. Local gauge invariance and a global slice inside the patch

Let delta_epsilon be the Hamiltonian flow of T_epsilon[f] with a smooth
real periodic f. The canonical brackets give

    delta_+ A_+=(f A_+)',   delta_+ A_-=0,
    delta_- A_-=-(f A_-)',  delta_- A_+=0,
    delta_+ q=(f B_+ + f')/2,
    delta_- q=(f B_- - f')/2.                         (4)

**LPS-T6.** The on-shell C and P are invariant under both complete families
of local infinitesimal gauge transformations. They have reduced bracket
{C,P}=1.

Proof. Integrating delta A shows delta P=0. Substitution of (1) into (4)
gives, for the plus flow,

    delta_+ C = (f A_+'/A_+ + f')/2
                -(f A_+)'/(2A_+) =0,

and the minus flow cancels with the opposite sign in the same way.
For an off-shell functional representative use the spatial average of
(2). Since {q(x),P}=1 and {A_epsilon(x),P}=0, its bracket with P is
exactly one, even before restriction. Its value on shell is C. Thus this
is weak Dirac invariance under the actual local generators, not merely
invariance under homogeneous lapse.

Every solution in this patch is gauge-equivalent, under the orientation-
preserving circle reparametrizations generated independently by the two
families in (4), to A_+=A_-=P. To see existence explicitly, the maps

    y_epsilon(x)= (1/P) integral_0^x A_epsilon(z) dz    (5)

are degree-one lifts of smooth circle diffeomorphisms: their derivatives
are strictly positive and y(x+1)=y(x)+1. A density becomes the constant
P in this coordinate. The inverse or direct flow is chosen according to
the sign convention in (4); both belong to the connected orientation-
preserving diffeomorphism group. Such maps are reachable by a smooth
time-dependent periodic vector field, for example by interpolating their
lifts with the identity, whose derivative remains positive. Along the
flow the density remains nowhere zero. The two chiral flows commute on
currents and can straighten them independently. Equations (3) and the
invariance of C then determine the corresponding smooth q,k,phi fields;
the remaining additive phi constant changes by a finite smooth flow.

On the resulting slice, p=P, phi'=k=q'=0 and

    q=Q=C+log|P|.                                    (6)

Residual constant smearings keep the currents constant and shift phi_0
by P(f_++f_-)/2. Since P!=0, phi_0 can be fixed to any scalar-clock reading
tau. The other constant combination is a spatial translation acting
trivially on these homogeneous fields. Any two slices with the same C,P
and chosen phi_0 coincide; gauge transformations cannot identify different
C,P by their invariance. Thus (Q,P), with both P signs and P!=0, parameterize
the quotient by these connected gauge flows on this patch. This is not
a claim about gauge transformations outside the patch or a quantum
representation of the entire gauge group.

Pull back the original canonical one-form to this gauge slice:

    integral [p dq+k dphi] dx = P dQ.                 (7)

It gives the reduced canonical bracket {Q,P}=1. The regular reduction of
the full local constraints on this patch therefore yields precisely the
HQS geometry pair, with a vanishing V=0 clock Hamiltonian.

## 3. The geometric length and its quantization

**LPS-T7 (relational local length).** The proper circumference evaluated
on the simultaneous current-straightening and scalar-clock gauge slice is

    L_rel = L_ref exp(Q) = L_ref |P| exp(C).           (8)

It is a positive weak Dirac observable of the original local model on
this patch, and {L_rel,P}=L_rel. It coincides with the original proper
circumference integral e dx on the slice, not on every gauge representative.

Indeed the unevaluated circumference is

    L_raw=L_ref exp(C) integral sqrt(|A_+A_-|) dx
          <= L_ref exp(C)|P| = L_rel.                (9)

Cauchy-Schwarz proves the inequality because both density integrals have
absolute value |P|. Equality holds exactly when A_+=A_- on this patch:
proportionality and their equal integrals force equality. This equality
can also occur on nonuniform representatives related by a common spatial
reparametrization. Under a normal deformation delta q=N k,

    delta_N L_raw = integral N e k dx,               (10)

which is generally nonzero. For an explicit smooth witness take
A_+=P(1+epsilon sin(2 pi x)), A_-=P, 0<|epsilon|<1, and use (3).
Then k is not identically zero. Choosing the smooth lapse N=k makes
(10) strictly positive. A spatial shift, in contrast, changes e by a
total derivative and leaves the circumference invariant. Thus replacing
L_rel by the raw local integral would fail the full Dirac test.

Quantize the reduced pair (Q,P) with the ordering and full-line matching
already specified in HQS, and use the LPS-1 unitary U. This defines the
positive self-adjoint quantization of the classical relational length (8)
on the local distributional solution sector:

    L_rel,phys = A^-1 [L_ref exp(q)] A,
    Dom L_rel,phys={chi:exp(q) A chi in L^2(dq)}.       (11)

Its physical commutator with P and spectral properties are exactly LPS-T4.
This is quantization **after** the demonstrated regular local reduction,
followed by an explicit all-constraint solution-space embedding. It provides
a classical local geometric meaning for that reduced operator. It is not
a construction of exp(q(x)), log(A_+(x)A_-(x)) or inverse currents as
operators on the kinematic Fock space. No assertion that quantum states
are localized on the classical nonvanishing-current patch is used.

The quantum matching through P=0 is still an extra HQS domain choice:
the classical patch has two components and does not select that matching.
Nor does this reduction prove that LPS's selected distributional sector
exhausts Dirac quantization. Comparing quantization before and after
reduction beyond the explicit embedded sector remains a separate issue.

## 4. New endpoint

LPS-1 supplies actual generalized local solutions and a positive physical
Hilbert space. LPS-2 identifies a regular classical patch of the full local
model whose reduction produces their chosen geometry observables. The
selected sector now has a justified relational length, not just an
abstract observable renamed as length.

Next: investigate the remaining local solution modules/rigging ambiguity,
the nonvanishing-current patch boundary and quantum relational field
operators. Matter, potentials, global cut gluing and physical 3+1 gravity
remain open. Finite native jet identities below certify source substitutions
and gauge cancellations; the global diffeomorphism, positivity and operator
arguments above are written proofs, not claims of machine formalization.
