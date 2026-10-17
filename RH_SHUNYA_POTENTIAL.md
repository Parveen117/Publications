# SPT-1 — Shunya potential and nonlinear cut transport

17 October 2026. Written continuation at Publications commit
ab727070cde8545d954a33750a9e4a009c43ab1e, following
[CEM-1](RH_CUT_EQUILIBRIUM_MEMORY.md).

The proposed term **shunya potential** is assigned a precise mathematical
object here. In the native retained/complementary coordinates (r,h), put
\[
\chi(y)=y(1-y),\qquad
\boxed{\phi_y(r,h)=\tfrac12\bigl(r^2+\chi(y)h^2\bigr).}
\]
It is conserved by the CEM transport. More strongly, every pointwise
polynomial invariant of that transport is a polynomial in y and phi_y.
Thus this potential generates the complete polynomial invariant algebra
of the specified two-channel flow.

Functions of this potential generate an exactly solvable nonlinear
extension. In particular the exponential potential
F_eta(s)=(exp(eta s)-1)/eta, with F_0(s)=s, changes the transport rate
while preserving phi_y. Repeated cuts still freeze the retained channel,
and their rescaled second response is derived exactly. Positivity of that
response for every polynomial profile remains equivalent to RH.

The TVSP connection is supplied separately and explicitly. Its four
corner labels and Legendre differential structure persist on every
regular source chart and through every available derivative order.
The invariant object across corner representations is the common
contact one-form; U,H,F,G are generally different numerical potentials.
The Maxwell relations require an exact constitutive potential. An
idempotent cut alone does not imply them.

**Status: RH remains unproved.** The new claims are written analytic
derivations, not machine-checked source certificates. No scientific
computation, numerical replay, independent agent or workflow dispatch
is performed. The naming does not claim a new general Hamiltonian,
Legendre or contact-geometry theory.

## Source and notation contract

Retain CEM's actual finite-anchor functional L_a, with a>=3249/4,
the polynomial carrier V=R[y], and its bilinear forms
\[
B_a(f,h)=\mathcal L_a[fh],\qquad
M_a(f,h)=\mathcal L_a[\chi fh].
\tag{1}
\]
The actual source points satisfy 0<|u_mu|<1 and u_mu!=1,
their nonzero conjugate weights have sum |w_mu|<=2,
and their only accumulation point is 0. The weights are not
assumed positive. OQD and CEM give
\[
\operatorname{rad}B_a=\operatorname{rad}M_a=\{0\},\qquad
M_a(f,f)\ge0\ \hbox{for all real polynomial }f
\Longleftrightarrow\mathrm{RH}.
\tag{2}
\]
All products at complex source points remain bilinear.

For the original native state x, use the exact decomposition
\[
x=Jr+k_h,\qquad Jf=(f,f),\qquad k_h=((1-y)h,-yh).
\]
Its invertible polynomial coordinate change puts the CEM cut,
pairing and generator into the form
\[
\Pi=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
D_y=\begin{pmatrix}1&0\\0&\chi\end{pmatrix},\qquad
A_0=\begin{pmatrix}0&2\chi\\-2&0\end{pmatrix}.
\tag{3}
\]
Here A_0 is the same (W-W^(-1))/2 written in these coordinates;
it is not a new source matrix. CEM proves
A_0^2=-4chi I, Pi A_0 Pi=0 and A_0^T D_y+D_y A_0=0.

The symbols below have distinct roles.

| Symbol | Role |
|---|---|
| y | Native source variable; its actual values may be complex |
| a | Fixed finite arithmetic source anchor |
| r,h | Retained and complementary state coordinates |
| t | Declared transport parameter |
| n | Repeated-cut count |
| eta | New scalar deformation of a potential |
| N | Finite derivative or Taylor depth, where stated |
| Lambda, lambda_X | The retained TVSP response section and its components |

Neither eta nor n is silently identified with an earlier lambda
coordinate. The TVSP covariant response tower retains its tensor
types and selected connections.

The term shunya refers to the zero compressed motion Pi A_0 Pi=0,
not to a demand that phi_y have value zero. The zero vector remains
a fixed vector of linear transport; a nonzero cut profile need
not be an ambient stationary state.

## SPT-T1. The quadratic potential generates all polynomial invariants

On R[y,r,h] define the native transport derivation
\[
\mathcal D=2\chi h\,\partial_r-2r\,\partial_h,\qquad
\mathcal D y=0.
\]
Then
\[
\boxed{
\ker\mathcal D=\mathbb R[y,\phi_y],\qquad
\phi_y=\tfrac12(r^2+\chi h^2).
}
\tag{4}
\]
In particular phi_y is conserved along exp(tA_0).

**Proof.** Direct differentiation gives D phi_y=0, proving one
inclusion. For the converse pass temporarily to the coefficient
field K=R(y) and its algebraic extension containing s with
s^2=-chi. Since chi is a nonzero rational function, s!=0.
The invertible linear coordinates
\[
z_+=r+sh,\qquad z_-=r-sh
\]
satisfy
\[
\mathcal D z_+=-2s z_+,\qquad
\mathcal D z_-=2s z_-.
\]
Thus a monomial z_+^j z_-^k has eigenvalue 2s(k-j).
In characteristic zero a polynomial killed by D contains only
j=k terms. It is therefore a polynomial in
z_+z_-=r^2+chi h^2.

Now return to an invariant I in R[y,r,h]. Its expression is
sum_m b_m(y)(r^2+chi h^2)^m over the extension field.
Set h=0. Since I(y,r,0) is an ordinary polynomial in y,r,
uniqueness of its coefficients of r^(2m) forces each b_m(y)
to belong to R[y]. This proves (4). Conversely an invariant
of exp(tA_0) is killed by D on differentiating at t=0; a
polynomial killed by D is constant along every such orbit.
QED.

This is a classification of pointwise polynomial observables
of the specified two-channel flow. It does not classify all
nonlocal functionals coupling different source values.

For a quadratic invariant the statement is especially rigid.
A symmetric polynomial matrix Q obeys
A_0^TQ+QA_0=0 exactly when
\[
\boxed{Q=b(y)D_y.}
\tag{5}
\]
Indeed the off-diagonal entry of Q must vanish and its lower
diagonal entry must equal chi times its upper one. Thus the
native pairing selects phi_y up to a coefficient factor and
the conventional normalization 1/2. Invariance alone does
not select an absolute physical unit.

At the exceptional scalar fibers y=0 or 1 the generator
degenerates and its fiberwise invariant algebra can enlarge.
Equation (4) concerns the polynomial family before specializing.
Neither exceptional fiber is an actual source point in (2).

## SPT-T2. The source potential contains the retained reading and the transverse memory

Define the source shunya potential by
\[
\boxed{
\mathscr S_a[r,h]
:=\mathcal L_a[\phi_y(r,h)]
=\tfrac12 B_a(r,r)+\tfrac12 M_a(h,h).
}
\tag{6}
\]
It is conserved by the full native transport. Its transverse
second variation is precisely CEM's memory:
\[
\boxed{
\left.\partial_\epsilon\partial_\delta
\mathscr S_a[r,\epsilon f+\delta h]\right|_{\epsilon=\delta=0}
=M_a(f,h).
}
\tag{7}
\]
A fixed projection instead satisfies
\[
\boxed{
\mathscr S_a[r,h]-\mathscr S_a[r,0]
=\tfrac12M_a(h,h).
}
\tag{8}
\]
Consequently conservation under full transport does not imply
invariance under deletion of the complementary channel.

**Proof.** Equation (6) is the exact CEM decomposition.
Equation (4), or A_0^T D_y+D_y A_0=0, proves pointwise
conservation, which remains true after applying L_a.
Expand the quadratic form to prove (7) and (8). QED.

The series defining L_a apply to transported profiles: on the
closed unit disk the coefficient matrices and polynomial initial
profiles are bounded, and the exponential series is locally
uniform in t. The absolute weight budget permits source
evaluation and differentiation without positivity.

The condition S_a[r,h]>=0 for all polynomial pairs is equivalent
to RH. One implication follows by setting r=0 and using (2).
Under RH all actual y are in (0,1) with positive weights, so
both terms in (6) are nonnegative and the potential is strictly
positive on every nonzero polynomial pair.

Before that condition is proved, S_a is an invariant quadratic
potential, not an assumed positive energy. In particular (8)
has no unconditional sign.

## SPT-T3. Functions of the invariant generate exact nonlinear transport

Let F(s) be a real polynomial or a real entire function and
write b(s)=F'(s). With the declared canonical matrix
J_0=[[0,1],[-1,0]], define
\[
X_F(q)=2J_0\nabla_q F(\phi_y(q))
=b(\phi_y(q))A_0q,\qquad q=(r,h)^T .
\tag{9}
\]
Then phi_y, F(phi_y), and the source potential S_a are
conserved. The exact flow is
\[
\boxed{
q(t)=\exp\!\left[t\,b(\phi_y(q_0))A_0\right]q_0 .
}
\tag{10}
\]

**Proof.** The chain rule gives (9). The derivative of phi_y
along A_0q vanishes by T1, and multiplication by the scalar
b(phi_y) preserves that identity. Hence b(phi_y) is constant
on each trajectory. Equation (10) solves the resulting linear
equation with that trajectory's constant coefficient, and also
solves (9). Its initial value is q_0. The vector field is locally
Lipschitz, so this solution is unique. For polynomial initial
profiles, the coefficients in (10) are bounded on compact
source and t sets, which justifies applying L_a. QED.

The nonlinear dependence is in the initial-state-dependent rate.
Where b does not vanish, this family keeps the same trajectories
as the base flow, with changed parametrization. It is not a
claim to generate every possible nonlinear physical law.

There is also a converse within the polynomial Hamiltonian
class. A polynomial Hamiltonian Psi(y,r,h) preserves phi_y
exactly when its canonical bracket with phi_y vanishes,
equivalently D Psi=0. By T1, it must have the form
Psi=F(y,phi_y). Thus the specified class of potential-preserving
polynomial Hamiltonian dynamics is exhaustive, not just a list
of examples. Here differentiation of F in its second argument
replaces b in (9).

A useful entire family is
\[
\boxed{
F_\eta(s)=
\begin{cases}(e^{\eta s}-1)/\eta,&\eta\ne0,\\s,&\eta=0,\end{cases}
\qquad F_\eta'(s)=e^{\eta s}.
}
\tag{11}
\]
For it,
\[
q_\eta(t)=e^{\,t e^{\eta\phi_y(q_0)}A_0}q_0,\qquad
\phi_y(q_\eta(t))=\phi_y(q_0).
\tag{12}
\]
The scalar derivative never vanishes, even on the complex source.
Thus this deformation does not create new stationary profiles
by a vanishing rate factor. Its Pi-compressed velocity still
vanishes on the cut h=0.

The invariant shared by different eta flows is phi_y. The
numerical values F_eta(phi_y) generally differ with eta. For
each fixed eta, its own value is conserved along that flow.
This distinction prevents a changing potential function from
being called a numerically unchanged potential.

## SPT-T4. Nonlinear repeated cuts have an exact rescaled response at every depth

Retain CEM's entire scalar
\[
c_t(\chi)=\sum_{m=0}^\infty\frac{(-4\chi)^m t^{2m}}{(2m)!}.
\]
Starting from the cut state (r,0), one full F-flow step followed
by the same cut has retained coordinate
\[
\Psi_\delta(r)=r\,c_{\delta b(r^2/2)}(\chi).
\]
Let r_(n,F)(t) be n compositions of Psi_(t/n), and define
\[
R_{a,n,F}(t;f,r)=\mathcal L_a[f\,r_{n,F}(t)]
\quad(f,r\in\mathbb R[y]).
\]
Then for every integer n>=1,
\[
\boxed{
-\frac n4\,\partial_t^2R_{a,n,F}(0;f,r)
=\mathcal L_a[\chi f r\,b(r^2/2)^2].
}
\tag{13}
\]
Moreover r_(n,F)(t) tends to r locally uniformly in t and
on the closed source disk as n tends to infinity.

**Proof.** Applying (10) to (r,0) and projecting yields Psi.
Its Taylor expansion is
\[
\Psi_\delta(r)
=r-2\chi r\,b(r^2/2)^2\delta^2+O(\delta^4).
\tag{14}
\]
It has zero first delta derivative at zero, is the identity
there, and has input derivative one there. The second
derivative of n compositions with delta=t/n is therefore the
sum of n identical second derivatives divided by n^2:
\[
\left.\partial_t^2 r_{n,F}(t)\right|_{t=0}
=-\frac4n\chi r\,b(r^2/2)^2.
\]
Source evaluation gives (13).

For convergence, initial r is uniformly bounded on |y|<=1.
Choose a larger fixed bounded set of possible retained values.
Analyticity of b and c gives a uniform step bound
|Psi_delta(z)-z|<=C|delta|^2 on that set for sufficiently
small delta. For |t|<=T, the sum of n bounds is CT^2/n.
For large enough n this is smaller than the margin to the
chosen set's boundary, so induction keeps every iterate inside
it. The same bound then proves convergence to r. QED.

For F(s)=s, formula (13) is exactly CEM's bilinear memory M_a(f,r).
For a nonlinear F it is generally nonlinear in r; it must not
be relabeled as that same bilinear form.

For the exponential family (11) the response simplifies to
\[
\boxed{
-\frac n4\,\partial_t^2R_{a,n,F_\eta}(0;f,r)
=\mathcal L_a[\chi f r\,e^{\eta r^2}].
}
\tag{15}
\]
This is exact at every finite cut count. The limiting unscaled
reading is constant in t and conceals the same response.

## SPT-T5. Nonlinear potential deformation retains the full RH sign condition

For a fixed real entire F with F'(0)!=0, define its self-response
\[
\mathscr N_{a,F}[r]
=\mathcal L_a[\chi r^2 F'(r^2/2)^2].
\]
Then
\[
\boxed{
\mathrm{RH}
\Longleftrightarrow
\mathscr N_{a,F}[r]\ge0
\quad\hbox{for every real polynomial }r.
}
\tag{16}
\]
This includes every fixed real eta in (11).

**Proof.** Under RH the source points are real in (0,1) and
the weights are positive. Since F has real coefficients, its
derivative is real on r^2/2. Every summand in N is nonnegative.

Conversely substitute r=epsilon p, with p any fixed real
polynomial and epsilon real tending to zero. Uniformly on
the closed source disk,
\[
F'(\epsilon^2p^2/2)=F'(0)+O(\epsilon^2).
\]
The absolute source budget therefore gives
\[
\boxed{
\lim_{\epsilon\to0}
\frac{\mathscr N_{a,F}[\epsilon p]}{\epsilon^2}
=F'(0)^2 M_a(p,p).
}
\tag{17}
\]
Positivity for every polynomial amplitude implies M_a(p,p)>=0
for every p. Since F'(0)^2>0, equation (2) proves RH. QED.

Thus nonlinear rates and infinitely frequent cuts do not remove
the source positivity obligation. For the exponential family the
small-amplitude coefficient is exactly M_a, independent of eta.
Equation (16) is an equivalence and a diagnostic, not a proof of
its right-hand side. No actual negative source profile is supplied.

## SPT-T6. TVSP corners preserve a common differential structure at every regular layer

Use the original diagram's clockwise labels from the T--V
corner: 1 H, 2 G, 3 F, 4 U, with axes T(+y), V(+x),
S(-y), P(-x). A corner's diagram position is not its natural
independent-variable pair.

For a fixed-amount smooth constitutive source U(S,V), put
T=U_S and P=-U_V. On a regular patch where the corresponding
Legendre coordinate changes are invertible, define
\[
H=U+PV,\qquad F=U-TS,\qquad G=U+PV-TS.
\tag{18}
\]
Their common one-form is
\[
\boxed{
\begin{aligned}
\Theta
&=dU-T\,dS+P\,dV\\
&=dH-T\,dS-V\,dP\\
&=dF+S\,dT+P\,dV\\
&=dG+S\,dT-V\,dP .
\end{aligned}
}
\tag{19}
\]
The pullback of Theta to the constitutive surface is zero.
The four Maxwell identities are:

| Diagram corner | Potential and natural variables | Differential on the source | Maxwell relation |
|---|---|---|---|
| 1 | H(S,P) | dH=T dS+V dP | (partial T/partial P)_S=(partial V/partial S)_P |
| 2 | G(T,P) | dG=-S dT+V dP | (partial S/partial P)_T=-(partial V/partial T)_P |
| 3 | F(T,V) | dF=-S dT-P dV | (partial S/partial V)_T=(partial P/partial T)_V |
| 4 | U(S,V) | dU=T dS-P dV | (partial T/partial V)_S=-(partial P/partial S)_V |

**Proof.** Differentiate (18) with the product rule. The added
terms cancel exactly as shown in (19). Substitution of U_S=T,
U_V=-P makes the pullback zero. Equality of the two mixed
partial derivatives of each C^2 potential gives the table.
Equivalently, exterior differentiation on the source gives
\[
dT\wedge dS=dP\wedge dV.
\tag{20}
\]
Repeated ordinary differentiation preserves these identities
through every order for which the source is differentiable.
On an invertible chart the jet chain rule and local inverse
carry the identities to the other corner charts. QED.

For an infinite smooth source the ordinary identities hold at
every finite jet depth. Analyticity or a controlled completion
is needed to infer nonzero-endpoint values from an infinite
Taylor record. Finite depth does not supply missing derivatives.

This is the exact sense in which the corner structure survives
the response tower. It does not assert U=H=F=G, nor that their
higher derivatives are constant. Equation (18) displays their
differences explicitly.

The native covariant response tower still obeys TLM/NT:
L_1=Lambda, L_(k+1)=nabla L_k, with one extra covector slot
per derivative. Its ordered covariant slots can have the
curvature commutator from NT-8. Maxwell integrability does not
make that selected connection flat. NT-4 already supplies a
stable source with Maxwell identities and nonzero curvature.

The phrase "the cut gives Maxwell" is valid only with its
constitutive exactness hypothesis made explicit. The algebraic
condition Pi^2=Pi alone does not enforce (20). For example a
constant Pi can accompany responses T(S,V)=V and P(S,V)=0;
then T_V+P_S=1, so no potential U generates both responses.
This illustrates the missing hypothesis without altering any
retained source.

The zero one-form on a thermodynamic constitutive surface and
the zero compressed generator Pi A_0 Pi are differently typed
statements. They are not identified merely because both vanish.

## SPT-T7. The shunya potential gives a regular nonlinear TVSP test source

There is an explicit local mathematical adapter. Fix a real
scalar fiber y in (0,1), so chi>0, and normalized affine
coordinates
r=S-S_0, h=V-V_0. For any smooth F define
\[
u_F(S,V)=F\!\left(\tfrac12(r^2+\chi h^2)\right).
\tag{21}
\]
Let q=(r,h)^T and D=diag(1,chi). Its Hessian is
\[
\boxed{
H_F=F'(\phi)D+F''(\phi)(Dq)(Dq)^T.
}
\tag{22}
\]
It is strictly positive exactly when
\[
\boxed{
F'(\phi)>0,\qquad F'(\phi)+2\phi F''(\phi)>0.
}
\tag{23}
\]
For the exponential family F_eta with eta>=0 these inequalities
hold everywhere on this real two-dimensional chart.

**Proof.** Differentiate phi=q^TDq/2 twice to obtain (22).
Conjugate by D^(-1/2). The resulting matrix is
F' I+F'' zz^T, with z=D^(1/2)q and |z|^2=2phi.
Its eigenvalues are F' and F'+2phi F''. At q=0 they coincide.
This proves (23). For (11), they are e^(eta phi) and
e^(eta phi)(1+2eta phi), both positive when eta>=0 and phi>=0.
QED.

The finite Taylor potentials
\[
F_{\eta,N}(s)=
\sum_{m=0}^{N}\frac{\eta^m s^{m+1}}{(m+1)!}
\tag{24}
\]
also obey (23) for eta>=0 and s>=0: their first derivatives
have constant term one and nonnegative coefficients, and their
second derivatives are nonnegative. They converge with all
derivatives on compact sets to (11).

Therefore the four corner constructions and Maxwell identities
apply at every one of these finite potential truncations and
at their entire limit. This is a supplied source family, not
a claim that the existing covariant Lambda tower is a list of
powers phi^n. The actual response tower of each u_F is obtained
by the declared derivatives and connection, as in NT/TLM.

Temperature and pressure for this test source are
T=F' r and P=-F' chi h in the normalized chart. If positive
T,P,S,V are wanted, choose a suitable local patch with r>0,
h<0 and positive translated S,V. Nothing here identifies the
whole transport orbit with a physical thermodynamic process.

Crucially, this adapter fixes a real y in (0,1). It does not
prove that all actual arithmetic source values u_mu are real.
The positive Hessian of this test source cannot substitute for
the still-unproved source positivity in (16).

## SPT-T8. Compatible jets preserve the invariant through the whole retained tower

For the exact nonlinear flow in T3, write the convergent time
jets r(t)=sum r_j t^j and h(t)=sum h_j t^j. Pointwise conservation
is equivalent to the coefficient identities
\[
\boxed{
\sum_{i+j=N}\bigl(r_i r_j+\chi h_i h_j\bigr)
=\begin{cases}r_0^2+\chi h_0^2,&N=0,\\0,&N\ge1.\end{cases}
}
\tag{25}
\]
These identities persist under every compatible truncation.

For the existing TVSP source transport, retain TET's
H(t)=D^2u(gamma(t)) and
\[
w'=-\tfrac12H^{-1}H'w.
\]
Then the corresponding quadratic response potential
I(t)=w(t)^TH(t)w(t)/2 obeys
\[
\boxed{
I(t)=I(0),\qquad
\sum_{i+j+k=N}w_i^T H_j w_k
=\begin{cases}w_0^TH_0w_0,&N=0,\\0,&N\ge1.\end{cases}
}
\tag{26}
\]
For smooth nonanalytic sources the latter statement concerns
derivatives and formal jets, not an asserted convergent endpoint
series.

**Proof.** Expand the conserved quadratic expressions and
collect the indicated coefficient. For the TVSP expression,
\[
\frac{d}{dt}(w^THw)
=w^T(H'-A^TH-HA)w=0,\qquad A=\tfrac12H^{-1}H'.
\]
Symmetry of H,H' gives A^TH+HA=H'. Every coefficient of
total degree N uses only layers through N, so truncation to
that depth preserves it. This is precisely the ordered jet
product rule retained by TLM. QED.

The original layers may change and cancel in these sums.
Conservation of the full scalar does not force each layer to
vanish or each corner value to remain fixed. For general TVSP
source paths the transport is ordered as in TET; replacing it
by an unordered exponential is not justified by this invariant.

Equation (26) reuses TET's compatible response transport.
It does not identify that tangent transport with the distinct
arithmetic CEM flow. The shared quadratic conservation rule,
the explicit test adapter T7, and their different source
hypotheses are all retained.

## Outcome and remaining obligation

The name shunya potential now has a specified mathematical
meaning and a nontrivial uniqueness property:

- phi_y generates every pointwise polynomial invariant of the
  native two-channel flow.
- S_a=L_a[phi_y] combines the retained pairing and the transverse
  memory. Its conservation survives the specified nonlinear
  potential deformations.
- Nonlinear repeated cuts have the exact response (13). Their
  all-profile sign condition still contains the complete RH gate.
- The TVSP corners retain the common differential form (19),
  and their Maxwell relations survive every regular derivative
  depth. Their numerical potential values need not be equal.
- Compatible jets retain the invariant as coefficient
  cancellation, rather than as frozen individual layers.

The next arithmetic obligation is positivity of the actual
M_a or of the equivalent nonlinear self-response for all
profiles. The polynomial invariant classification and the
Maxwell identities do not supply that sign. TPW's two-prime
window and all earlier scope limits remain unchanged.

No universal propagation speed, physical zero-energy vacuum,
black-hole information law, new mass gap, larger positive
support window or all-order RH proof is inferred.

## Sources, attribution and verification

Retained repository sources were read at the parent commit:

- [CEM-1](RH_CUT_EQUILIBRIUM_MEMORY.md): exact cut coordinates,
  source memory, finite detector and arithmetic positivity gate.
- [OQD-1](RH_OBSERVER_QUOTIENT_DESCENT.md): faithful source pairing.
- [NT-1--NT-8](papers/native-thermodynamic-curvature/THEOREM.md):
  supplied potential, six responses, Maxwell differentials,
  selected connection and covariant response tower.
- [TLM-1](papers/native-tvsp-lambda-matrices/THEOREM.md):
  source jet algebra, tensor types and truncation rules.
- [TET-1](papers/native-tvsp-exponential-transport/THEOREM.md):
  ordered metric-compatible transport and source/path distinctions.
- [NJC-1](papers/native-jacobian-cut-certification/THEOREM.md):
  corrected tower hypotheses and limits of cut-only claims.

The supplied original TVSP diagram was read and its labels
retained. The earlier exponential manuscript was also inspected;
its unestablished interpretations are not substituted for
NJC's corrected hypotheses.

Classical attribution:

- [MIT 3.20, Recitation 5: Legendre transforms and Maxwell relations](https://ocw.mit.edu/courses/3-20-materials-at-equilibrium-sma-5111-fall-2003/68bfe3a448119ace697230a1aade3ccc_r5.pdf).
  The four potential differentials and their mixed-derivative
  argument are standard thermodynamics, derived again in T6.
- C. S. Lopez-Monsalvo, F. Nettel, V. Pineda-Reyes and
  L. F. Escamilla-Herrera,
  [Contact polarizations and associated metrics in geometric thermodynamics](https://arxiv.org/abs/2006.12754).
  Legendre changes are changes of contact representation.
  No assertion that every Hessian metric is Legendre invariant
  is made here.

The Hamiltonian first-integral mechanism and polynomial
weight-space proof are standard methods. The development is
their explicit native source application and nonlinear
cut-response formula, with the precise RH and TVSP boundaries.

Verification covers written proofs, retained-source definitions,
text structure, links, exact immutable-head readback and draft
status. No numerical scientific calculation, source evaluation,
independent agent, replay or CI run is added. Reading supplied
files and maintaining their file metadata is document handling.
