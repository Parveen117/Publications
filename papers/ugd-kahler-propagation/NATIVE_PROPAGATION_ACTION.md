# A native variational multiplier for a selected propagation law

Author: Monty Dabas. Development edition: 29 September 2026.

## Inputs and meaning of “native”

We use the real EMK coefficient relations

\[
R^2=-I,\quad K^2=I,\quad KR=-RK,\quad L=KR,\qquad
R^\dagger=-R,\quad K^\dagger=K,\quad L^\dagger=L.
\]

Admit two commuting tensor factors and the explicit real representation

\[
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
K=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
L=\begin{pmatrix}0&-1\\-1&0\end{pmatrix}.
\]

These matrices intertwine dagger with transpose. The canonical operator-engine
note uses K_can=[[0,1],[1,0]]; the present choice is K=K_can R and L=−K_can,
inside the same coefficient algebra. We are not silently changing its relations.

This is an admitted finite coefficient realization. It neither represents nor
eliminates the independent proper seam: finite-dimensional TS=I forces ST=I.
The continuum chart (t,x¹,x²,x³), a real four-component commuting field ψ,
its smoothness, and a local linear propagation ansatz are additional inputs.
The algebra alone does not select these, a physical clock, or an energy unit.
Hereafter † on matrices in this realization means transpose; a formal
differential adjoint also includes integration by parts.

## NP-1. Clifford closure and the positive pairing in this representation

Define

\[
\Gamma_0=R\otimes I,\quad
\Gamma_1=K\otimes K,\quad
\Gamma_2=K\otimes L,\quad
\Gamma_3=L\otimes I.
\]

The defining products give

\[
\{\Gamma_a,\Gamma_b\}=2\eta_{ab}I_4,\qquad
\eta=\operatorname{diag}(-1,1,1,1),
\]

with Γ₀†=−Γ₀ and Γᵢ†=Γᵢ. Set

\[
A_i=\Gamma_0\Gamma_i.
\]

Then Aᵢ†=Aᵢ, {Aᵢ,Aⱼ}=2δᵢⱼI, and tr Aᵢ=0. For constant coefficients,

\[
\mathcal D=\Gamma_0\partial_t+\Gamma_i\partial_i,
\quad \mathcal D^2=(-\partial_t^2+\Delta)I,
\quad \mathcal D\psi=0\iff B_0\psi=0,
\quad B_0=\partial_t-A_i\partial_i.
\]

The Euclidean component pairing makes the time coefficient I positive and
each Aᵢ symmetric: this is a symmetric hyperbolic system. That pairing is
explicitly compatible with the admitted dagger, not an axiom about a primitive
Hilbert space. On the coefficient algebra the functional τ=tr/4 gives
τ(X†X)>0 for X≠0. In the 16-element tensor basis
{I,R,K,L}⊗{I,R,K,L}, τ(E_a†E_b)=δ_ab; hence τ(X†X) is the sum of real
coefficient squares. The four-component module and its normalization are
still representation data.

These are standard real Clifford identities. A solution of the first-order
law satisfies the wave equation, but the second-order equation alone allows
extra Cauchy data and is not an equivalent replacement.

## NP-2. Constitutive coefficients give a Lorentzian characteristic form

Supply smooth real βᵢ(t,x) and an invertible real 3-by-3 V(t,x), and set

\[
M^i=\beta^i I+\sum_{a=1}^3 V_{ai}A_a.
\]

For a cotangent vector (τ,k), write s=τ−β·k and v=Vk. The principal matrix
of ∂t−Mⁱ∂i is sI−v·A. Since (v·A)²=|v|²I and tr(v·A)=0,

\[
\det(sI-v\cdot A)=(s^2-|v|^2)^2.
\]

Thus the characteristic quadratic form can be chosen as

\[
\boxed{q^{-1}(\xi,\xi)=-(\tau-\beta\cdot k)^2+|Vk|^2.}
\]

The invertible change (τ,k)↦(τ−β·k,Vk) proves signature (−,+,+,+).
Real characteristic roots are τ=β·k±|Vk|, each of multiplicity two away from
k=0. Singular V is explicitly excluded: it gives a degenerate quadratic form.
Lower-order terms do not change this principal symbol.

This law supplies a conditional single Lorentzian cone. It does not make
q a pseudo-Kähler metric or identify it with the Hessian lift. V, β, the
choice of four coordinates, and any dimensional calibration remain input.
V=c₀I gives an arbitrary positive speed c₀; the equations have not predicted c.
A characteristic cone alone fixes only a conformal class. The displayed
q⁻¹ uses the supplied normalization of the time derivative.

For comparison, requiring a linear symmetric symbol H(k)=kᵢBᵢ to satisfy
H(k)²=hⁱʲkᵢkⱼI forces {Bᵢ,Bⱼ}=2hⁱʲI by comparing coefficients.
Consequently the Clifford closure is exactly the scalar-square assumption,
not a consequence of thermodynamic stability alone.

## NP-3. The native volume gives the unique constant skew multiplier

Define

\[
\boxed{\mathsf J=A_1A_2A_3
=\Gamma_0\Gamma_1\Gamma_2\Gamma_3=K\otimes R.}
\]

The Clifford relations imply J²=−I, J†=−J and [J,Aᵢ]=0. This J is a
field-component operator, distinct from the tangent complex structure in
GS-1 and from the single return Γ₀. In particular Γ₀ anticommutes with Aᵢ
and cannot substitute for J in the action below.

More precisely, every constant real matrix S commuting with all three Aᵢ
has the form aI+bJ. To prove it explicitly, use

\[
A_3=\begin{pmatrix}I_2&0\\0&-I_2\end{pmatrix},\quad
A_1=\begin{pmatrix}0&K\\K&0\end{pmatrix},\quad
A_2=\begin{pmatrix}0&L\\L&0\end{pmatrix}.
\]

Commuting with A₃ makes S=diag(P,Q); commuting with A₁ makes Q=KPK.
Commuting with A₂ then makes P commute with LK=−R, so P=aI+bR and
Q=aI−bR. Therefore the skew matrices in this commutant are exactly bJ.
Every nonzero one is invertible.

This also proves the variational uniqueness statement. A constant-multiplier
first-order operator S(∂t−Aᵢ∂i) can be formally self-adjoint only if
S†=−S and (SAᵢ)†=−SAᵢ. Since Aᵢ†=Aᵢ, these are precisely
S†=−S and [S,Aᵢ]=0. Thus S=bJ. This statement is restricted to the given
four real commuting fields, constant invertible multipliers, and the same
first-order equation; it does not classify every possible action formulation.

## NP-4. Conservation fixes the symmetric drift

Let ρ(t,x)>0 be a prescribed density. Supply C(t,x) with C†=−C and define

\[
Q=\tfrac12\left[(\partial_t\log\rho)I
-\partial_i M^i-M^i\partial_i\log\rho\right],\qquad
B=\partial_t-M^i\partial_i+Q+C.
\]

Here and below repeated spatial indices are summed. Q†=Q. The formal
adjoint in the spacetime pairing ∫ρ u†v is computed using

\[
(\partial_t)^*=-\partial_t-\partial_t\log\rho,\qquad
(-M^i\partial_i)^*=M^i\partial_i+\partial_iM^i+M^i\partial_i\log\rho.
\]

It follows that B*=−B. Conversely, for the given ρ and Mⁱ, imposing
B*=−B fixes the symmetric part of its zero-order coefficient to Q;
the skew part C remains free. This is selection within a declared conservation
contract, not selection of the density or the principal law.

For a solution Bψ=0 the exact local balance is

\[
\boxed{\partial_t\left(\tfrac12\rho\,\psi^\dagger\psi\right)
-\partial_i\left(\tfrac12\rho\,\psi^\dagger M^i\psi\right)=0.}
\]

To check it, expand the two derivatives, substitute
∂tψ=Mⁱ∂iψ−Qψ−Cψ, and use ψ†Cψ=0. The remaining terms cancel by 2Q's
definition. On a periodic spatial domain or with zero boundary flux,
E(t)=∫ρ ψ†ψ/2 is conserved and positive on nonzero fields.
For decaying fields the same conclusion requires the boundary term to vanish.

E is an energy estimate/norm for this system. No assertion identifies it with
thermodynamic internal energy, observed matter energy or the Hamiltonian of
the action in NP-5. Smooth positive ρ and smooth bounded coefficients on a
regular local patch supply the usual symmetric-hyperbolic Cauchy setting;
global existence, singular backgrounds and arbitrary boundary conditions are
not claimed.

## NP-5. An action, derived by the formal-adjoint condition

Impose the further lower-order compatibility gate

\[
\boxed{[C,\mathsf J]=0.}
\]

Since J is constant and commutes with all Mⁱ and Q, it then commutes with B.
Together with B*=−B and J†=−J this proves

\[
(\mathsf J B)^*=B^*(-\mathsf J)=B\mathsf J=\mathsf J B.
\]

For any fixed nonzero real normalization Z, define

\[
\boxed{S_{\rm prop}[\psi]
=\frac Z2\int \rho\,\psi^\dagger\mathsf J B\psi\,dt\,d^3x.}
\]

For compactly supported field variations, keeping the background data fixed,

\[
\delta S_{\rm prop}
=Z\int\rho\,\delta\psi^\dagger\mathsf J B\psi,
\qquad \delta S_{\rm prop}=0\iff B\psi=0.
\]

The equivalence uses invertibility of J. The term ψ†JQψ vanishes pointwise
because JQ is skew. Q nevertheless appears in the Euler–Lagrange equation
through integration by parts of the variable derivative coefficients and
density. One may therefore omit Q from the written action integrand without
omitting it from the field equation.

This is a concrete inverse variational result: for the stated law its
first-order kinetic multiplier is forced up to overall scale, and its
weighted action is reconstructed. It is not yet the requested derivation of
a unique physical L_λ from uncut UGD. In particular ψ has not been identified
with λ, no rule V[λ] has been supplied, and no metric field equation follows
when V, β and ρ are held fixed.

With constant backgrounds, write B=∂t−G where G=Aᵢ∂i−C. The action's
Hamiltonian is H=Z∫ψ†JGψ/2. This quadratic functional need not be positive;
it is distinct from E. For example on a periodic one-dimensional wave, with
C=0 and G=A₁∂x, use ψ=v cos(kx)±JA₁v sin(kx). These two fields have equal
positive E, while H has opposite signs for k≠0. Positivity of the conserved
norm therefore does not prove a bounded physical Hamiltonian or a viable
quantized matter theory.

## NP-6. A sharp counterexample: conservation does not force this action

At constant coefficients set C=mΓ₀ for real m≠0. It is skew, so the
NP-4 conserved norm still holds. But J anticommutes with Γ₀, hence

\[
[C,J]\ne0,\qquad \tfrac12(JC+(JC)^\dagger)
=\tfrac12(JC+CJ)=0.
\]

The quadratic action with multiplier J loses this entire mass term. Its
Euler–Lagrange equation cannot be Bψ=0. NP-3 shows no other constant
invertible skew multiplier repairs it on the same four real commuting
components. Doubling fields, admitting other variations, or using Grassmann
fields would be different problems and are not ruled out.

This counterexample even retains a gapped wave equation. Let
G_m=Aᵢ∂i−mΓ₀. Since {Γ₀,Aᵢ}=0, G_m²=Δ−m²; solutions obey
(∂t²−Δ+m²)ψ=0. Thus a Lorentzian principal symbol, a positive conserved
norm and a mass gap still do not guarantee the particular variational
completion in NP-5. This exposes the actual additional compatibility gate.

## Relationship to the Einstein target and prior results

The original request has three different questions. GS-1 rules out a direct
Lorentzian pseudo-Kähler metric; NP-2 constructs a Lorentzian propagation
symbol in a separately specified continuum. GS-4 solves one scalar
Monge–Ampère problem; NP-5 reconstructs an action for the selected
four-component propagation law. Neither result supplies a curvature action
or a physical identification between those sectors.

To derive Einstein dynamics, a next theorem must select the coframe/metric
dynamics and permitted variations, then show that the resulting Euler–Lagrange
equation reduces consistently to Einstein's equation. The existing
[thermo classical-field section](../thermodynamic-response-corrections/source/sections/06b_classical_field_equations.tex)
does prove that reduction under its declared action class. That conditional
result is preserved. No value of G, c, Λ, α or ℏ follows from this package.

Clifford linearization, symplectic variational principles, and
symmetric-hyperbolic energy estimates are established mathematical methods.
For a related primary treatment of hyperbolic field systems and their energy
estimates see J. M. M. Senovilla,
[Symmetric hyperbolic systems for a large class of fields in arbitrary dimension](https://arxiv.org/abs/gr-qc/0607120).
The present proof is self-contained and does not identify its coefficient
module with the tensor-field constructions in that paper. Novelty or
experimental discovery is not inferred from a new theorem label.
