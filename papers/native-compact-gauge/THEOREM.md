# Native compact gauge completion

Author: Monty Dabas. 2 October 2026. Results NCG-1–NCG-8.

The existing EMK connection calculus already supplies curvature and Bianchi
identities. The massive-matter construction SM supplies a real phase multiplier
and a U(1) action, while the Yang–Mills benchmark supplies SU(2) separately.
This chapter constructs a bridge: the maximal internal symmetry preserving the
declared native norm, Clifford action and massive multiplier is U(n); its first
non-Abelian member has an explicit two-EMK-factor realization. Connection,
positive quadratic readout, source and constraint identities then follow on
that carrier. These are mathematical results under the contracts below, not a
selection of nature's gauge group, action, number of copies or physical units.

The quaternion/SU(2) identification, unitary commutant, Yang–Mills variational
identities and Noether mechanism are established mathematics. The contribution
here is their explicit integration with the source-bound native carrier,
massive multiplier and transport conventions, with exact executable controls.
No external priority claim is made. Applications and remaining physics gates
are maintained in extra-ideas, not silently added to the old YM gap certificate.

## Contract and notation

Use the associative real EMK algebra with R²=-1, K²=1, KR=-RK,
R†=-R, K†=K. In this chapter **L=RK**. SM/NP instead write L=KR;
their spin matrices are consumed with that opposite L sign explicitly.
Different numbered EMK factors commute. The coefficient-of-1 functional is sc.
No Hilbert space, Haar measure, physical time or spacetime metric is primitive.
Real scalar completion is admitted for the continuous groups; finite controls
use exact rational scalars and the unchanged RKF rewriting engine.

For the module theorem, admit the SM real four-component coefficient module S:
its kinetic and massive Clifford operators generate End_R(S)=M₄(R).
The replicated carrier is V=W⊗S, W=R^(2n), with the derived coefficient norm
and a fixed normalized, invertible skew multiplier C on W:
C†=-C, C²=-1. Copies, their mutual orthogonality, and this fixed multiplier are
declared data. A proper one-sided seam is not represented faithfully in these
finite coefficient modules. The general UGD container is not identified with
this associative EMK sector; UGD-A0 remains in force.

Smooth base directions, a coframe, a metric and a local action class enter only
where explicitly admitted below. Mathematical proofs, finite executions,
proof-assistant evidence and physical identification have separate statuses.

## NCG-1. Positive native pairing and the one-factor limit

For x=a+bK+cR+dL,

\[
\operatorname{sc}(x^\dagger x)=a^2+b^2+c^2+d^2,
\qquad \operatorname{sc}(xy)=\operatorname{sc}(yx).
\]

On a product of independent EMK factors, the same statement is the sum of
squared coefficients in the tensor-word basis. The group U†U=1 preserves this
pairing by both left multiplication and conjugation. In one real EMK factor
the skew space is exactly R·R, so its connected unitary group is only SO(2).
The three noncommuting generators R,K,L are not three compact gauge generators.

**Proof.** Multiplication gives b†b=1 for each basis word b and scalar part zero
for b†c with b≠c. Tensoring these orthogonal identities proves the coefficient
formula. Checking the sixteen basis products proves cyclicity and bilinearity
extends it. Thus sc((Ux)†Ux)=sc(x†x), and cyclicity proves conjugation invariance.
Only R is skew in the one-factor basis. The representation
R=[[0,-1],[1,0]], K=diag(1,-1), L=[[0,1],[1,0]] is faithful, preserves dagger,
and identifies its full unitary group with O(2), whose identity component is
SO(2). The representation checks the identification after the native identities.
For contrast, -sc(K²)=-1, while sc(K†K)=1. A positive coefficient norm cannot be
replaced by minus the square trace on an arbitrary noncompact generator. ∎

## NCG-2. The action-compatible internal symmetry is U(n)

Restrict internal frame changes to maps that commute with every Clifford
operator, preserve the positive coefficient norm, and preserve the fixed
massive multiplier C. Their full group is

\[
\mathcal G_C=\{g\in O(2n):g^T Cg=C\}
           =\{g\in O(2n):[g,C]=0\}\simeq U(n).
\tag{1}
\]

Equivalently these are the norm-preserving Clifford-commuting maps preserving
the SM bilinear H=-C⊗Γ⁰. Their Lie algebra consists of X†=-X, [X,C]=0.
Two real copies of S give U(1); four give U(2), the first non-Abelian group
in this fixed massive-module class. They have respectively 8 and 16 real
coefficient components. This is not a minimality statement for all possible
modules, Grassmann variables, auxiliary fields or alternative actions.

**Proof.** Since the Clifford words span M₄(R), commuting with all of them
forces an endomorphism of W⊗S to be g⊗1. This can be checked by commuting with
the sixteen matrix units of End(S): the S blocks become scalar multiples of
the identity. Norm preservation gives g^Tg=1. Preservation of H is exactly
g^TCg=C, and orthogonality turns this into commutation with C. Choose an
orthonormal basis (v₁,Cv₁,…,v_n,Cv_n), obtained inductively because C is
orthogonal and skew. C then acts as multiplication by the derived quarter-turn
scalar. Real orthogonal maps commuting with it are precisely complex-linear
unitaries, proving (1). A real invertible skew multiplier requires even
multiplicity (its determinant changes by (-1)^r under transpose); multiplicity
two is Abelian, and four is sufficient by NCG-3. Normalizing a general skew
multiplier to C²=-1 is a stated restriction; unequal frequency blocks can
reduce the symmetry to a product of smaller unitary groups. ∎

## NCG-3. Two native factors give the quaternion sector and U(2)

In two independent EMK factors put

\[
C=R_1,\qquad e_1=R_2,\quad e_2=R_1K_2,\quad
e_3=R_1R_2K_2,\qquad t_a=e_a/2,\quad t_0=C/2.
\tag{2}
\]

Then C²=e_a²=-1, all four generators are skew, [C,e_a]=0, and

\[
e_ae_b=-\delta_{ab}+\epsilon_{abc}e_c,\qquad
[t_a,t_b]=\epsilon_{abc}t_c,\quad [t_0,t_a]=0.
\tag{3}
\]

The full skew centralizer of C in the sixteen-dimensional tensor-word algebra
is span_R{C,e₁,e₂,e₃}. The unit elements in span{1,e₁,e₂,e₃} form SU(2);
including central C gives U(2)=(U(1)×SU(2))/{(1,1),(-1,-1)}.

**Proof.** The defining EMK relations give e₁e₂=e₃, e₂e₃=e₁, e₃e₁=e₂;
reversing each order changes its sign. They also give every square, dagger and
commutation assertion. In the sixteen-word basis, first impose [X,R₁]=0,
leaving first-factor coefficients only in span{1,R₁}; imposing X†=-X leaves
exactly the four stated words. For u=a₀+Σa_ae_a, u†u=Σ_(a=0)^3 a_a².
The faithful two-by-two cut-scalar adapter

\[
u\longmapsto
\begin{pmatrix}
a_0+\iota a_1&a_2+\iota a_3\\
-a_2+\iota a_3&a_0-\iota a_1
\end{pmatrix},\qquad \iota^2=-1,
\tag{4}
\]

preserves products and dagger; its determinant is Σa_a². Every determinant-one
unitary two-by-two matrix has this form, proving SU(2). C becomes the central
quarter-turn on the complex two-dimensional module. A unitary matrix's
determinant has a square root, so dividing by it leaves SU(2); the two choices
give the stated two-element kernel. In particular the surviving U(1) is not
silently discarded or confused with one of the non-Abelian directions. ∎

## NCG-4. Native non-Abelian curvature, balance and commuting sectors

Admit derivations δ_i with [δ_i,δ_j]=c_ij^k δ_k as in EMK-C1 and connection
coefficients A_i=Σ_(a=0)^3 A_i^a t_a with real scalar functions. On the module
∇_i=δ_i+A_i,

\[
F_{ij}=\delta_iA_j-\delta_jA_i+[A_i,A_j]-c_{ij}^kA_k,
\tag{5}
\]
\[
F^a_{ij}=\delta_iA_j^a-\delta_jA_i^a-c_{ij}^kA_k^a
          +\epsilon_{bca}A_i^bA_j^c\quad (a=1,2,3).
\tag{6}
\]

The central component has no bracket term. For g†g=1 and [g,C]=0,
A'_i=g⁻¹A_i g+g⁻¹δ_i g, F'_(ij)=g⁻¹F_(ij)g. EMK-C1's Bianchi identity applies.
If A_i and Ω_i act on independent internal and Clifford factors, the curvature
of δ_i+A_i⊗1+1⊗Ω_i is F^A_(ij)⊗1+1⊗F^Ω_(ij); there is no mixed commutator.
For a prior adjoint-valued memory one-form M, the balanced curvature is
F(A-M)=F(A)-D_A M+M∧M, as in NT-5. Gauge covariance does not set this to zero.

**Proof.** Expand [∇_i,∇_j]-c_ij^k∇_k and use (3); this proves (5)–(6) and
identifies curvature as the native transport-order defect. Substitution of
the transformed coefficients cancels derivative cross terms. The Jacobi
identity gives Bianchi, including the direction-frame terms. Independent
tensor factors commute, proving the split. Expanding A-M gives the memory
formula. For A_x=t₁,A_y=t₂ on a commuting base, every coefficient derivative
vanishes but F_xy=t₃≠0. For A=g⁻¹dg the derivative terms instead cancel A∧A.
Neither all off-diagonal entries nor all noncommuting coefficients certify
curvature without this full calculation. With ∂_s w=-A(∂_s)w and later path
segments multiplying on the left, the positive x,y,-x,-y square is
U=1-h²F_xy+O(h³), retaining NT's sign convention. ∎

## NCG-5. Positive invariant readout and its remaining coefficients

On the real skew centralizer define B(X,Y)=-4 sc(XY). Then B(t_A,t_B)=δ_AB
for A,B=0,1,2,3. It is positive and invariant under the group (1).
Every invariant symmetric real bilinear form on u(2) is

\[
b_0 X^0Y^0+b_1\sum_{a=1}^3X^aY^a;
\tag{7}
\]

it is positive definite exactly when b₀,b₁>0. On the su(2) sector alone there
is one free positive coefficient. Thus invariance fixes equal color weights
and forbids central/color mixing, but does not fix either coupling magnitude.

For a small unitary based-loop defect D=U-1=-h²F+O(h³),

\[
4h^{-4}\operatorname{sc}(D^\dagger D)=B(F,F)+O(h).
\tag{8}
\]

This is a positive loop record, distinct from LR's oriented gravitational
pairing of complementary planes. It requires a declared loop/readout protocol.

**Proof.** Equation (3) and the four scalar parts give orthonormality.
Cyclicity proves B([Z,X],Y)+B(X,[Z,Y])=0. For a general symmetric 4×4
coefficient matrix, invariance under the three t_a forces its central/color
vector to vanish (it is fixed by every three-dimensional rotation) and its
color block to be scalar (commuting with rotations in all coordinate planes).
This proves (7) without assuming a chosen Euclidean matrix chart is invariant.
Expansion of D†D proves (8). The real one-factor noncompact algebra is an
essential control: conjugating R by g=(5+3K)/4 changes sc(R†R)=1 to 257/32.
Thus the positive norm is not invariant under arbitrary invertible native
changes of frame. Restricting to the admitted dagger-unitary group is load-bearing. ∎

## NCG-6. Variational curvature law and Gauss propagation

Admit a fixed smooth Lorentzian metric q, w=√(-det q), commuting coordinates,
compactly supported variations and the su(2) quadratic action

\[
S_A=-\frac1{4g_*^2}\int w B(F_{\mu\nu},F^{\mu\nu})\,d^4x,
\qquad g_*^2>0.
\tag{9}
\]

Let a gauge-invariant module action have variation
δS_m=∫δψ^T E_ψ-∫w B(j^ν,δA_ν), and gauge change
δψ=-εψ, δA_ν=D_νε. Then

\[
D_\mu(wF^{\mu\nu})=g_*^2wj^\nu,\qquad
\bigl[D_\nu(wj^\nu)\bigr]^a=E_\psi^T t_a\psi.
\tag{10}
\]

Put E_A^ν=D_μ(wF^{μν})/g_*²-wj^ν. On the matter equation,
D_νE_A^ν=0. Consequently the temporal equation
G=D_i(wF^{i0})/g_*²-wj⁰=0 propagates: if all spatial E_A^i=0,
D₀G=0. Zero initial Gauss defect remains zero. This is a classical differential
identity, not a quantum constraint-algebra or continuum-existence theorem.

**Proof.** δF_μν=D_μδA_ν-D_νδA_μ. Antisymmetry, invariance of B and integration
by parts give δS_A=∫B(D_μ(wF^{μν})/g_*²,δA_ν). Adding the matter variation
gives the first equation. The stated infinitesimal gauge change and integration
by parts give the second coefficient by coefficient, off shell. Finally
D_νD_μ(wF^{μν})=½[F_νμ,wF^{μν}]=0: the contraction pairs each color-commutator
with a symmetric coefficient in the two antisymmetric index pairs. This holds
for variable q and density w as well. It proves the divergence and Gauss
identities. The action, metric, locality, derivative order and coefficient in
(9) are admitted inputs. Positivity of the internal pairing is not by itself
reflection positivity of a quantum measure. ∎

## NCG-7. A massive module with an actual non-Abelian source

Use four copies of SM's spin module. With the internal C,t_a from (2), set
γ^a=1₄⊗Γ^a, H=-C⊗Γ⁰, \(\widehat t_a=t_a\otimes1_4\),
Ω_μ=ω_abμ γ^aγ^b/4 and D_μ=∂_μ+Ω_μ+A_μ^a \(\widehat t_a\).
The independent spin connection and coframe are the SM geometric adapter.
The commuting real module action

\[
S_m=\frac Z2\int w\left(\psi^T H\gamma^\mu D_\mu\psi
                         -m\psi^T H\psi\right)d^4x
\tag{11}
\]

is invariant under local internal SU(2) (or the full U(2) extension) and the
existing local spin action. Its current is

\[
j_a^\mu=-\frac Z2\psi^TH\gamma^\mu\widehat t_a\psi.
\tag{12}
\]

It obeys (10), including the non-Abelian color-precession term. As in SM, the
actual independent-connection matter equation contains the torsion trace:
(γ^μD_μ-m+½T^ν_(μν)γ^μ)ψ=0. This construction supplies a classical source;
it does not give fermion statistics or a positive quantum matter Hamiltonian.

**Proof.** C commutes with all internal generators and all Clifford matrices.
H is symmetric and Hγ^a is skew. Each \(\widehat t_a\) is skew and commutes
with H and γ^a, so \(\widehat t_a^T H+H\widehat t_a=0\). The spin bilinear
identities of SM hold on every replicated block. These facts prove local
invariance when connections transform with their inhomogeneous terms.
Differentiating (11) with respect to each A_μ^a gives (12). Varying ψ and
integrating its derivative gives
H[γ^μ∂_μ+(2w)⁻¹∂_μ(wγ^μ)+½{γ^μ,Ω_μ+A_μ}-m]ψ;
the SM tetrad identity gives the displayed torsion trace. Gauge invariance
then proves (10) off shell. Internal and spin variations commute by NCG-4.
The multi-copy axial spin density may be inserted in SM's linear torsion
equation, but SM's special eight-component quartic Fierz simplifications
must not be copied to this sixteen-component module without a new proof. ∎

## NCG-8. Complete coefficient variations and the Abelian-gradient trap

Reuse CP's stable response patch
U=10s-v+s²+sv+5v²/2+s²v/2,
Δ=5(2+v)-(1+s)², Q=s, P=1/6-1/(2√Δ), P_v=5/(4Δ^(3/2))>0.
Admit the coefficient extraction

\[
A=\sum_{a=1}^3\sum_{\mu=0}^3 t_a P_{a\mu}\,dQ_{a\mu}.
\tag{13}
\]

Twelve independent stable channels suffice locally near A=0 to realize every
su(2) connection coefficient and every compactly supported coefficient
variation: set Q_aμ=x^μ, P_aμ=A_μ^a and invert the CP chart. At the reference
point the variation matrix is (5/108)I₁₂, including at F=0. All twelve gauge
Euler equations are therefore tested. Sixteen channels give the U(2) version.
Independent coframe channels preserve MG's full metric variations; sharing
channels requires a new joint-rank proof. No minimal channel count is asserted.

**Proof.** The local inverse is
v=[(1+Q)²+1/(4(P-1/6)²)]/5-2, with P<1/6. Around Q=P=0 it stays inside
the same strict stability, temperature and pressure inequalities as CP.
With Q fixed, δA_μ^a=P_v,aμ δv_aμ independently for every pair. Strict
positivity proves surjectivity, also at A=F=0. Equation (13) uses the selected
Darboux part of each thermo form, a constitutive extraction rather than a
claim that all exact terms can be discarded as independent non-Abelian gauges.
Indeed A=t₁ dx+t₂ dy has each coefficient individually exact but F_xy=t₃.
Deleting both as 'pure gradients' would falsely flatten it. A lawful gauge
transformation uses g⁻¹Ag+g⁻¹dg, transporting all components together. Thus
the Abelian CP reasoning needs this explicit extraction/variation interface
before its components are assembled into the non-Abelian connection. ∎

## Evidence and boundary

The verifier checks the complete finite linear centralizer and invariant-form
problems, native multiplication/word replays, independent rational matrix
identities, variable connection jets, action derivatives and negative controls.
The eight general arguments above are written proofs; finite fixture checks
do not establish arbitrary smooth-function identities by sampling. No new Lean
formalization, external review, empirical result or whole-engine rerun is claimed.

Neither n=2 nor SU(2) instead of U(2) is selected by primitive laws. SU(3), the
Standard Model's representations, chirality, anomaly cancellation, spin-statistics,
coupling constants and quantum gravity require further work. The YM benchmark's
time-direction operator bound E4D-C, cutoff trajectory and Clay dictionary remain
open. A common finite algebra is not an intertwiner between different regulated
theories. Its measure, transfer operator and continuum limits need separate proofs.

## Primary and native sources

- Yang and Mills, *Conservation of Isotopic Spin and Isotopic Gauge Invariance*,
  Physical Review 96, 191 (1954), https://doi.org/10.1103/PhysRev.96.191.
- Jaffe and Witten, official Yang–Mills existence and mass-gap problem:
  https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf.
- EMK-C1 and NT-1/2/5 in Publications at d98c2644ab1512024861e12d5468f32028903cdf.
- SM-1/2/3, CP-5/8, MG and LR on the research physics branch at
  f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56; that branch is not silently merged here.
- Canonical RKF engine at 3cc5a33b05c16d59c90994ddda69dedc0d392424.

Exact file pins and the source-to-result map accompany this chapter.
