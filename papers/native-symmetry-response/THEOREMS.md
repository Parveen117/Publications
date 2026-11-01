# NSR1: electromagnetic response curvature and symmetry-selected higher response

Monty Dabas. 5 October 2026.

## Constitutive and native contract

Use a local, reversible, nondispersive electromagnetic energy density
U(D,B). Conduction, hysteresis, frequency dispersion, and spatial nonlocality
are outside this adapter. In full vector notation its conjugates are
E=∂U/∂D and H_mag=∂U/∂B. The constitutive energy and electromagnetic Maxwell
equations are admitted physical input; neither follows from the native
multiplication table alone. The Hamiltonian convention is [1], equations
(1.8)–(1.9).

For the two-channel sector use positive fixed scales U_*, D_*, B_* and
q=(x,y)=(D_y/D_*,B_z/B_*), u=U/U_*. Set

\[
e=u_x,\quad h=u_y,\qquad
H=D^2u=\begin{pmatrix}a&b\\b&c\end{pmatrix},\quad
\Delta=ac-b^2,\qquad a>0,\quad\Delta>0.                 \tag{1}
\]

H denotes a response matrix throughout, not the magnetic field. Assume u is
C4 on this stable affine chart unless stronger smoothness is explicitly used.
Physical electric and magnetic field scales are E_*=U_*/D_* and
H_*=U_*/B_*. A raw Hessian under an arbitrary nonlinear coordinate replacement
must not be substituted for the transported response metric.

The inherited native algebra is K²=1, R²=-1, KR=-RK, L=RK. NT-1/3/4 give

\[
\mathbf H=\frac{a+c}{2}1+\frac{a-c}{2}K+bL,\qquad
A_i=\frac12\mathbf H^{-1}\partial_i\mathbf H,\qquad
F_{ij}=-\frac14[X_i,X_j],\quad X_i=\mathbf H^{-1}\partial_i\mathbf H.
                                                               \tag{2}
\]

The connection is selected by compatibility with the positive response
pairing and the absence of its H-skew part in the declared trivialization.
Its faithful two-dimensional representation is A_i=H⁻¹∂_iH/2. Since H is
a Hessian, this is also the Levi-Civita connection of the response metric.
For an arbitrary positive matrix family the selected connection still makes
sense, but the Levi-Civita identification needs the Hessian condition.

This response connection acts over the **space of bias values**. It differs
from the U(1) connection over spacetime. Noncommuting response matrices do not
change the electromagnetic gauge group.

## NSR1-T1. Constrained electric and magnetic response

The incremental constitutive equations are

\[
de=a\,dx+b\,dy,\qquad dh=b\,dx+c\,dy.                 \tag{3}
\]

The normalized constrained susceptibilities and mixed responses are

| Response | Constraint | Value |
| --- | --- | --- |
| ε_y = (∂x/∂e)_y | fixed magnetic flux | 1/a |
| ε_h = (∂x/∂e)_h | fixed magnetic field | c/Δ |
| μ_x = (∂y/∂h)_x | fixed electric displacement | 1/c |
| μ_e = (∂y/∂h)_e | fixed electric field | a/Δ |
| (∂x/∂h)_e | fixed electric field | -b/Δ |
| (∂y/∂e)_h | fixed magnetic field | -b/Δ |

Define r=b/√(ac), so |r|<1. Then

\[
\boxed{\frac{\epsilon_h}{\epsilon_y}
       =\frac{\mu_e}{\mu_x}
       =\frac{ac}{\Delta}=\frac1{1-r^2}\geq1.}          \tag{4}
\]

The corresponding relaxed stiffnesses are a-b²/c and c-b²/a.
For physical susceptibilities multiply the two ε values by D_*²/U_* and
the two μ values by B_*²/U_*. Mixed responses carry D_*B_*/U_*.
The gain (4) is independent of these units.

**Proof.** Hold y or x fixed in (3), or invert H. For fixed h eliminate
dy=-b dx/c; for fixed e eliminate dx=-b dy/a. Division by the constrained
values gives (4). Stability is precisely the strict Schur condition.
The equality of mixed derivatives is potential integrability, not a general
Onsager assertion about dissipative transport at magnetic bias. ∎

**Use.** Constraint-dependent permittivity and permeability must be compared
under a specified preparation. The common gain is a consistency diagnostic
for this conservative two-channel model. These Schur and inverse-response
identities are established mathematics, not newly discovered material laws.

## NSR1-T2. Curvature from coupling and impedance

Set s=√(ac)>0 and Z=√(a/c)>0. Thus

\[
H=s\begin{pmatrix}Z&r\\r&Z^{-1}\end{pmatrix}.
                                                               \tag{5}
\]

In NT-4's oriented orthonormal frame the connection is Rα, where

\[
\alpha=\frac{db-(b/a)\,da}{2\sqrt\Delta}
       =\frac{dr-r\,d\log Z}{2\sqrt{1-r^2}}.
                                                               \tag{6}
\]

Its curvature coefficient two-form is

\[
\boxed{f=d\alpha
       =-\frac{dr\wedge d\log Z}{2(1-r^2)^{3/2}}.}      \tag{7}
\]

It vanishes at a point exactly when the two covectors dr and d log Z are
linearly dependent there. Constant r or constant Z on a patch is sufficient
for flatness. Merely having r=0 at one point is not sufficient, since dr
need not vanish there. Pointwise reciprocity of H does not impose (7)=0.

For H=D²u the scalar Gaussian response curvature is

\[
\kappa=-\frac{f_{xy}}{\sqrt\Delta}
       =-\frac{1}{4\Delta^2}
          \det\begin{pmatrix}
          a&b&c\\a_x&b_x&c_x\\a_y&b_y&c_y
          \end{pmatrix}.                               \tag{8}
\]

Thus second and third derivatives of u determine the curvature; the fourth
derivatives cancel. A common positive multiplier H↦w(q)H leaves F and f
unchanged, while the normalized readout -f_xy/√Δ scales by 1/w. The rescaled
matrix need not itself be a Hessian, so its readout is not automatically a
Gaussian curvature.

**Proof.** Insert a=sZ and b=sr into the NT-4 expression: db-(b/a)da
=s(dr-r d log Z), and √Δ=s√(1-r²), proving (6). The derivative of a
function of r times dr vanishes, while
d[r/(2√(1-r²))]=dr/[2(1-r²)^(3/2)], proving (7). Equivalence of its zeros
uses |r|<1. NT-4 gives f_xy=-κ√Δ. Substitute the Hessian third-derivative
symmetries into (2), multiply the represented F by H, and take (HF)_12/Δ;
the result is (8). Independently this follows by differentiating the full
Christoffel formula. Under multiplication by w, X_i gains the central term
∂_i log w, which cancels from every commutator. Also r and Z are unchanged. ∎

**Use.** The formula isolates a response-shape diagnostic that common
stiffness changes cannot imitate. For parallel transport defined by
dw/dt=-A(γ_dot)w, the oriented orthonormal return angle is -∮α=-∫f on
a contractible smooth patch. This is the declared response-transport
protocol. Maxwell propagation has not thereby been shown to implement it.

## NSR1-T3. Wave measurements reconstruct the response geometry

Assume the (D_y,B_z) sector is dynamically closed, all fields depend only
on spatial coordinate ξ and time t, and there are no sources. The omitted
field components must not couple linearly into this sector. This condition
holds for the Born–Infeld sector in T5; it is not automatic for a general
anisotropic medium. Put v_*=U_*/(D_*B_*), which has velocity units. Maxwell's
equations reduce to

\[
\partial_t q+v_*\begin{pmatrix}0&1\\1&0\end{pmatrix}
                  \partial_\xi(\nabla u)=0.            \tag{9}
\]

About a constant bias, the signed speeds normalized by v_* are

\[
\boxed{v_\pm=b\pm\sqrt{ac}=s(r\pm1).}                 \tag{10}
\]

Their incremental electric-to-magnetic field ratios are ±Z in the
normalized field units. In physical units the impedance magnitude is
(B_*/D_*)Z. Therefore

\[
s=\frac{v_+-v_-}{2},\quad
b=\frac{v_++v_-}{2},\quad
r=\frac{v_++v_-}{v_+-v_-},\quad a=sZ,\quad c=s/Z.       \tag{11}
\]

Measurements of these quantities across two independent bias directions
determine f by (7). A measurement at just one bias does not.

**Proof.** The linearized matrix in (9) is [[b,c],[a,b]]. Its characteristic
polynomial is (v-b)²-ac. Since Δ>0, v_+>0>v_-. Its eigenvectors satisfy
δy/δx=±√(a/c). Substitution into (3) gives δe/δh=±Z, using the nonzero
eigenvalues. Solving for the three entries of H proves (11). ∎

For every admissible constant bias, frequencies obey ω=v_*v_±k in this
linearized sector. They approach zero with k. Consequently a nonzero
response-curvature diagnostic here does **not** create a photon mass gap.
The two signs label opposite propagation directions, not two polarization
indices; their difference is not a birefringence claim.

## NSR1-T4. Symmetry selects the tower without fixing its coefficients

Let a Cn potential u on a finite-dimensional real response space satisfy
u(Gq)=u(q) for a specified constant invertible linear symmetry G. Then

\[
D^nu(Gq)[Gv_1,\ldots,Gv_n]=D^nu(q)[v_1,\ldots,v_n].    \tag{12}
\]

For the Hessian and a continuous generator G(t)=exp(tT),

\[
H(Gq)=G^{-T}H(q)G^{-1},\qquad
\boxed{(Tq)\cdot\nabla H+T^TH+HT=0.}                  \tag{13}
\]

The native response connection and its induced tensor connections can be
transported through this same representation; the negative action on each
covariant slot is NT-8. Symmetry constrains all orders, but (12) supplies
neither a unique u nor numerical material coefficients.

For a smooth even potential u(-q)=u(q), every odd derivative at zero
vanishes. In particular ∂_iH(0)=0 and F(0)=0. If u is C6 near zero,

\[
F_{ij}(q)=-\frac14
 [H_0^{-1}C_i(q),H_0^{-1}C_j(q)]+O(|q|^4),\quad
(C_i(q))_{ab}=\sum_k u_{abik}(0)q_k.                  \tag{14}
\]

Thus fourth-order energy response can first appear as a curvature signal
quadratic in bias even when the zero-bias signal is forced to vanish.
Independent electric and magnetic parity symmetries give further selection
rules. For u(-x,y)=u(x,y)=u(x,-y), the Taylor monomial x^m y^n can appear
only when both m and n are even. In a parity- and time-reversal-invariant
electromagnetic preparation these flips follow from D being polar and
time-even and B axial and time-odd. Fixed symmetry-breaking material order
or external parameters must also be transformed before applying the rule.

Full rotational symmetry also permits curved response. In an abstract
two-channel sector consider

\[
u(q)=\tfrac12|q|^2+\tfrac\beta4|q|^4,\qquad\beta>0.
\]

It is SO(2)-invariant and globally strictly convex. At q=(t,0),

\[
\boxed{\kappa(t,0)=
-\frac{2\beta^2t^2}{(1+3\beta t^2)^2(1+\beta t^2)^2}.} \tag{15}
\]

**Proof.** Differentiate the invariance identity n times by the chain rule
to obtain (12). Differentiating the Hessian identity in t gives (13).
Parity removes odd Taylor tensors; inserting ∂_iH=C_i(q)+O(|q|³) and
H⁻¹=H_0⁻¹+O(|q|²) into (2) proves (14). Independent sign changes prove
the two parity selection rules. The radial model has
H=(1+β|q|²)I+2βqq^T. At (t,0) its eigenvalues are 1+3βt² and 1+βt²,
and its third-derivative rows in (8) are (6βt,0,2βt) and (0,2βt,0).
Their determinant is 8β²t², proving (15). ∎

Equation (15) supplies a general rotation-symmetry example. Rotation of
these two abstract coordinates is not identified with electromagnetic
duality of the transverse pair (D_y,B_z). Full electromagnetic duality
mixes vector D and B component by component and requires all its channels.
The established nonlinear duality literature is [2].

## NSR1-T5. Exact nonlinear electromagnetic witness

Take the established Born–Infeld Hamiltonian [1], equation (3.42), with
its positive scale set to one and its constant vacuum value subtracted:

\[
u_{\rm BI}(\mathbf D,\mathbf B)
=\sqrt{1+|\mathbf D|^2+|\mathbf B|^2
                      +|\mathbf D\times\mathbf B|^2}-1.
\]

In the closed transverse sector D=(0,x,0), B=(0,0,y),

\[
u(x,y)=\sqrt{(1+x^2)(1+y^2)}-1.                       \tag{16}
\]

Put Q=(1+x²)(1+y²). Direct differentiation gives

\[
a=\frac{\sqrt{1+y^2}}{(1+x^2)^{3/2}},\quad
b=\frac{xy}{\sqrt Q},\quad
c=\frac{\sqrt{1+x^2}}{(1+y^2)^{3/2}},\quad
\Delta=\frac{1-x^2y^2}{Q}.
\]

The positive-response patch is |xy|<1. On it,

\[
r=xy,\quad Z=\frac{1+y^2}{1+x^2},\quad
v_\pm=\frac{xy\pm1}{\sqrt Q},                         \tag{17}
\]

\[
\boxed{\kappa_{\rm BI}(x,y)=
\frac{x^2+y^2+2x^2y^2}{\sqrt Q\,(1-x^2y^2)^2}.}         \tag{18}
\]

This is positive except at zero. The signed speed magnitudes are at most
one: Q-(xy±1)²=(x∓y)²≥0. This response curvature therefore coexists with
bounded wave speeds in an established nonlinear electromagnetic model.

At the exact bias x=y=1/2,

\[
H_0=\begin{pmatrix}4/5&1/5\\1/5&4/5\end{pmatrix},\quad
\Delta=3/5,\quad r=1/4,\quad Z=1,\quad
v_+=1,\quad v_-=-3/5,\quad \kappa=128/225.             \tag{19}
\]

A quadratic potential with this same H_0, chosen also to match the energy
and gradient at that bias, has identical point susceptibilities, speeds,
and impedance but zero response curvature everywhere. In particular Z=1
at the point does not mean dZ=0 there.

**Proof.** Differentiate (16) twice. Its full-vector conjugates stay along
the indicated y and z axes when all fields depend only on ξ, proving
closure under the curl equations. Equation (17) follows from T2/T3, while
dr=y dx+x dy and d log Z=-2x dx/(1+x²)+2y dy/(1+y²). Substitution into
(7) gives

\[
f=-\frac{x^2/(1+x^2)+y^2/(1+y^2)}{(1-x^2y^2)^{3/2}}
           dx\wedge dy.
\]

Divide by -√Δ to obtain (18); evaluation gives (19). A quadratic potential
has constant Hessian, so all its X_i vanish. ∎

Born–Infeld's absence of birefringence is an established property [1].
The distinction proved here is that response curvature probes bias
derivatives beyond point linear propagation; it is not a new birefringence
prediction. At |xy|=1 the selected positive response metric degenerates.
Its curvature divergence is not asserted to be a singularity of the full
Born–Infeld field theory. No measured vacuum Born–Infeld scale is inferred.

## NSR1-T6. Internal-mode release can create curved effective response

Let U(q,z) be a smooth joint potential and let z=z_*(q) solve U_z=0 with
C=U_zz positive and invertible. The implicit function theorem gives a local
smooth equilibrium branch. For U_eff(q)=U(q,z_*(q)),

\[
\boxed{H_{\rm eff}=U_{qq}-U_{qz}C^{-1}U_{zq}.}         \tag{20}
\]

If the full Hessian is positive, so is H_eff. Its varying Schur correction
can yield nonzero (7) even when the response at fixed z is a constant matrix.

For an exact local example admit

\[
U(x,y,z)=\tfrac12(x^2+y^2+z^2)-zxy.
\]

Then z_*=xy and U_eff=(x²+y²)/2-x²y²/2. At x=y=1/2, z_*=1/4,

\[
H_{\rm clamped}=\begin{pmatrix}1&-1/4\\-1/4&1\end{pmatrix},\quad
H_{\rm eff}=\begin{pmatrix}3/4&-1/2\\-1/2&3/4\end{pmatrix},\quad
\boxed{\kappa_{\rm clamped}=0,\qquad\kappa_{\rm eff}=128/25.} \tag{21}
\]

Both response matrices and the full three-variable Hessian are positive
at this point, hence in a neighborhood. No global stability of this cubic
joint potential is asserted. Under electric or magnetic sign reversal z
must reverse too for the joint interaction to be invariant; z is an
admitted internal mode with those parities, not automatically ordinary strain.

**Proof.** Differentiating U_z(q,z_*(q))=0 yields
∂z_*/∂q=-C⁻¹U_zq. The chain rule and U_z=0 give (20). Completing the
quadratic form proves the Schur positivity claim. In the example the
clamped Hessian is independent of x,y. The effective Hessian is
[[1-y²,-2xy],[-2xy,1-x²]], whose determinant at the specified point is
5/16. The determinant numerator in (8) is -4(x²+y²), so its curvature
there is (x²+y²)/Δ²=128/25. The full Hessian is positive because C=1
and this Schur complement is positive. ∎

**Use.** Compare clamped and quasistatically relaxed measurements to locate
an internal contribution to nonlinear response. Dynamic elimination needs
an evolution law and generally a frequency-dependent memory kernel. The
static identity (20) does not replace the Yang–Mills Schur-memory problem.

## What is established and what is next

Six written results extend the native thermo-response construction into a
specified electromagnetic constitutive setting and general symmetry sectors.
The new synthesis is the explicit coupling/impedance formula, its wave
reconstruction, and these comparisons within the same native connection.
Standard Hessian geometry, symmetry differentiation, Schur complements,
and Born–Infeld dynamics are acknowledged inputs, not priority claims.

An experiment could estimate the diagnostic from a two-bias response map.
Calling its holonomy an actual optical phase would additionally require a
proof or calibrated apparatus that implements the selected transport.
Lossy/retarded media require an extension retaining frequency, preparation,
and hidden state. No electromagnetic mass, constant, or Yang–Mills continuum
gap has been derived here.

## Sources

[1] L. Mezincescu, J. G. Russo, P. K. Townsend, *Hamiltonian birefringence and
Born–Infeld limits*, JHEP 02 (2024) 186,
https://arxiv.org/html/2311.04278v3. Equations (1.8)–(1.9) fix the Maxwell
and constitutive conventions; (3.42) supplies the known Born–Infeld model.

[2] G. W. Gibbons, D. A. Rasheed, *Electric-Magnetic Duality Rotations in
Non-Linear Electrodynamics* (1995),
https://arxiv.org/abs/hep-th/9506035. Used for established nonlinear
electromagnetic duality context, not as a source for the new native readout.

[3] [Native thermodynamic curvature, NT-1/3/4/8](../native-thermodynamic-curvature/THEOREM.md),
and [EMK-C1 connection calculus](../emk-ugd-algebra/CONNECTION_CALCULUS.md).
These supply the native response carrier, selected connection, orientation,
Hessian representation, and tensor-slot action reused above.
