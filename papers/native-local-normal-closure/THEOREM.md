# LNC-1 — A local canonical completion of the normal bracket

## What is selected and what is proved

[LCF-1](../native-local-coframe-covariance/THEOREM.md) supplies a coframe density and cancels the mixed local obstruction, but its matter-only normal generator omits geometry-field normal transport. Here a selected local canonical extension supplies that transport and closes all three **classical 1+1-dimensional** secondary-constraint brackets, including the massive LCF matter coupling.

The canonical geometry variables, an additional real scalar `phi`, the action ansatz and its potential are declared inputs. The coefficient relations, brackets, source terms and constraint preservation are derived. Closure does not select the potential: two models have the same algebra and the same initial constraint data but different coframe acceleration. This is a precise check against mistaking an invariant of a chosen model for a uniquely recovered law of nature.

This is not a new discovery of two-dimensional dilaton-gravity constraint structures. The construction is of that established type; for primary context see R. Meyer, *Constraints in two-dimensional Dilaton Gravity with Fermions*, [hep-th/0512267](https://arxiv.org/abs/hep-th/0512267). That work uses fermionic constraints and a local Lorentz sector; we do not import its conclusions for our different canonical matter contract. The calculations below rederive the stated model's brackets explicitly. Neither extra scalar nor canonical Poisson structure is claimed to follow uniquely from the native cut.

Reuse: EMK's associative direction algebra; LCF's coframe and source; RCG's distinction between record invariance and model selection; QGC's distinction between classical feedback and autonomous quantum evolution. All predecessor packets remain frozen.

## 1. Phase space and action

Work locally in one spatial coordinate on smooth fields, with compactly supported smearings and variations, so integrations by parts have no surface terms. Let `e(x)>0` be the coframe density, `pi(x)` its conjugate scalar, and `(phi(x),p_phi(x))` a scalar/density canonical pair:

\[
 \{e(x),\pi(y)\}=\delta(x-y),\qquad
 \{\phi(x),p_\phi(y)\}=\delta(x-y). \tag{1}
\]

For classical commuting complex two-component matter fields, choose
`{psi_a(x),psi_b^*(y)}=delta_ab delta(x-y)/(i eta)` with `eta>0`. This is the classical symplectic contract behind the linear LCF matter action. It is not a fermionic quantization or a physical quantum state-space construction.

Set `C=diag(1,-1)`, `R=[[0,-1],[1,0]]`, and reuse

\[
 D_f=f\partial_x+f'/2,\qquad A_N(e)=C D_N+gRNe,
 \qquad J=i\eta g\psi^*R\psi. \tag{2}
\]

For any real smooth potential `V(phi)`, define the geometry densities

\[
 \boxed{h_g=e\pi p_\phi+\frac{e'}e\phi'-\phi''+e^2V(\phi),}
 \qquad d_g=-e\pi'+p_\phi\phi'. \tag{3}
\]

Define smeared constraints

\[
 \begin{aligned}
 H[N]&=\int N h_g\,dx+H_m[N],&
 H_m[N]&=i\eta\int\psi^*A_N(e)\psi\,dx,\\
 D[f]&=\int f d_g\,dx+D_m[f],&
 D_m[f]&=i\eta\int\psi^*D_f\psi\,dx.
 \end{aligned} \tag{4}
\]

`H_m[N]` is linear in `N`; integration by parts gives its density
`(i eta/2)(psi^* C psi'-psi^*' C psi)+e J`.
The action is

\[
 S=\int ds\left[\int(\pi\dot e+p_\phi\dot\phi+i\eta\psi^*\dot\psi)dx
                 -H[N]-D[\beta]\right]. \tag{5}
\]

`N,beta` are multiplier fields. Their variation imposes `H(x)=D(x)=0`; their primary conjugate momenta, if added, are the usual vanishing multiplier momenta and are not the new result. The operation parameter `s` is not calibrated physical time. Positivity of `e`, boundary conditions and smoothness are local-domain hypotheses, not a proof of global existence.

## 2. Why these derivative coefficients appear

The temporary chart `q=log e`, `p=e pi` is canonical because `pi delta e=p delta q`; it is used only to shorten the analytic calculation. Native arithmetic does not evaluate logarithms. Equations (3) become

\[
 h_g=p p_\phi+q'\phi'-\phi''+U,
 \quad d_g=pq'+p_\phi\phi'-p',\quad U=e^{2q}V(\phi). \tag{6}
\]

The spatial generator gives

\[
 \delta_f q=fq'+f',\quad \delta_f p=(fp)',\quad
 \delta_f\phi=f\phi',\quad\delta_f p_\phi=(fp_\phi)'. \tag{7}
\]

Within the declared ansatz `h=a p p_phi+b q'phi'+c phi''+U`, with constant `a,b,c`, its departure from weight-two covariance is

\[
 \delta_f h-fh'-2f'h=(b+c)f''\phi'. \tag{8}
\]

Direct product differentiation proves (8). Thus this ansatz requires `c=-b`. The normal bracket is

\[
 \{H_g[N],H_g[M]\}
 =ab D_g[\omega]+a(b+c)\int\omega p' dx,
 \qquad \omega=NM'-MN'. \tag{9}
\]

Its derivation is the same functional-derivative calculation below with these coefficients retained. Matching the inherited matter derivative bracket requires `ab=1`. We select the representative `a=b=1,c=-1`. These are conditions within this ansatz, not a uniqueness theorem among all gravitational actions. In particular no relation fixes the function `V`.

## 3. The full classical constraint algebra

For `b(f,h)=fh'-hf'`,

\[
 \boxed{\begin{aligned}
 \{D[f],D[h]\}&=D[b(f,h)],\\
 \{D[f],H[N]\}&=H[b(f,N)],\\
 \{H[N],H[M]\}&=D[NM'-MN'].
 \end{aligned}} \tag{10}
\]

These are exact smeared identities on the declared smooth phase space, before imposing the constraints. The normal smearing `N` has weight minus one: it is the densitized lapse of LCF, with scalar lapse `alpha=Ne`. Equation (10) is not a claim that the scalar-lapse hypersurface algebra has constant structure coefficients in higher dimensions.

### Geometry brackets

The transformations in (7), or `delta e=(fe)'`, `delta pi=f pi'`, directly give the geometry spatial bracket. They also make `h_g` a weight-two density: the `f'' phi'` terms from `q'phi'` and `-phi''` cancel. Therefore
`{D_g[f],H_g[N]}=H_g[fN'-Nf']` after integration by parts.

For the normal bracket the four functional derivatives are

\[
 \begin{aligned}
 (H_g[N])_q&=-(N\phi')'+2NU,& (H_g[N])_p&=Np_\phi,\\
 (H_g[N])_\phi&=-(Nq')'-N''+NU_\phi,&
 (H_g[N])_{p_\phi}&=Np.
 \end{aligned} \tag{11}
\]

Substitution into the canonical bracket cancels every potential term and gives the pointwise bracket integrand

\[
 \omega(p_\phi\phi'+pq')+\omega'p
 =\omega d_g+(\omega p)'. \tag{12}
\]

The last term integrates to zero for the stated support. This proves the third geometry bracket and (9) by keeping the coefficients. Omitting `-phi''` would lose the derivative-of-momentum part of `D_g` and would also break (8).

### Matter and mixed terms

For prescribed geometry, quadratic classical matter functionals satisfy

\[
 \{\psi^*A\psi,\psi^*B\psi\}_m=(i\eta)^{-1}\psi^*[A,B]\psi
\]

with spatial integration understood. Thus LCF's operator identities imply

\[
 \{H_m[N],H_m[M]\}_m=D_m[\omega],
\]
\[
 \{D_m[f],H_m[N]\}_m
 =H_m[b(f,N)]+\int NJ(fe)'dx. \tag{13}
\]

Geometry cancels the last term exactly:

\[
 \{D_g[f],H_m[N]\}=-\int NJ(fe)'dx. \tag{14}
\]

Indeed `(D_g[f])_pi=(fe)'` and `(H_m[N])_e=NJ`. Finally,

\[
 \{H_g[N],H_m[M]\}=-\int NM\,e p_\phi J\,dx, \tag{15}
\]

which is symmetric in `N,M`, so it cancels its antisymmetrized partner. All other cross terms are zero. Adding (12)–(15) proves all of (10), for nonzero `g` and arbitrary smooth `V`. No discarded matter-source term is needed to obtain closure.

## 4. Normal geometry evolution and the retained source

Hamilton variation of (5) gives

\[
 \begin{aligned}
 \dot e&=Ne p_\phi+(\beta e)',\\
 \dot\phi&=Ne\pi+\beta\phi',\\
 \dot\pi&=-N\pi p_\phi+\frac{(N\phi')'}e-2NeV-NJ+\beta\pi',\\
 \dot p_\phi&=(N e'/e)'+N''-Ne^2V'+(\beta p_\phi)',\\
 \dot\psi&=(A_N(e)+D_\beta)\psi.
 \end{aligned} \tag{16}
\]

The geometry now has a normal component `dot e|_N=Ne p_phi`. Applying two normal Hamilton vector fields to `e` yields `-(omega e)'`, with the sign following the convention `X_H F={F,H}`. This is the field-direction bracket absent in LCF's matter-only generator. Equation (10), rather than a literal replacement inside LCF's fixed-background differential operators, is the completion on the enlarged canonical phase space.

Along shift characteristics the first equation is multiplicative in `e`, so a smooth solution with bounded coefficients stays positive over its interval of existence. This is conditional preservation, not a global PDE well-posedness theorem.

The matter term `-NJ` in `dot pi` is exactly LCF's independent-coframe variation. Locally any positive `e` can be realized by a monotone circle-cut lift

\[
 \theta(x)=\theta_0+2\int_{x_0}^x e(y)dy,
 \quad e=\theta'/2.
\]

For compactly supported variations of that lift, the same matter functional obeys
`delta H_m/delta theta=-(NJ)'/2`, recovering LCF's cut-source formula. This is a local realization and a variational identity, not an asserted pointwise canonical equality between `pi` and `p_theta`. Global winding, boundary data and the integration constant remain extra record data; on a periodic cut field the allowed integral of `e` is constrained by winding and is not automatically preserved by arbitrary normal multipliers in this local model.

For fixed test smearings, (10) gives

\[
 \dot D[f]=H[b(f,N)]+D[b(f,\beta)],\qquad
 \dot H[M]=D[b(M,N)]+H[b(M,\beta)]. \tag{17}
\]

Consequently any smooth solution starting on the constraint surface remains there, for admissible multipliers and no boundary flux. Establishing existence of such solutions is separate from this algebraic preservation implication.

## 5. Closure cannot select the physical dynamics

Choose zero matter and spatially homogeneous initial data
`e=1, pi=0, phi=0, p_phi=0`, with `N=1`, `beta=0`. On a local spatial patch both models below satisfy `H=D=0` initially:

| Selected potential | Initial `dot p_phi` | Initial `ddot e` |
|---|---:|---:|
| `V(phi)=0` | `0` | `0` |
| `V(phi)=phi` | `-1` | `-1` |

Both obey exactly (10). The first admits a static solution. The second admits the homogeneous local solution
`phi=pi=0`, `e=sech(s)`, `p_phi=-tanh(s)`; direct differentiation verifies (16) and both constraints. The certificate checks the initial acceleration and canonical constraints using exact rational data, not transcendental evaluation. On a noncompact line these homogeneous examples are local solutions, not claims of finite total integrated energy.

Thus equal closure and equal initial constraint values do not determine equal dynamics. If two closed theory presentations are compared, potential, boundary, canonical and measurement contracts must remain in the record. An invariant under common comparison transport need not be an invariant of every proposed physical evolution.

## 6. Boundary of this result

The LCF normal-bracket obligation is solved **within the selected classical local model**. This is a stronger result than another finite unitary example: the bracket proof holds for arbitrary smooth compact smearings and any smooth potential, including nonzero matter rest coupling.

It is still a classical 1+1 model with an additional scalar, not completed quantum gravity. In particular:

- The extra canonical pair and action ansatz are not derived uniquely from the original cut information.
- The kinematic cross-momentum term is not bounded below; physical-sector positivity needs reduction and analysis.
- Commuting classical matter fields are not physical quantized fermions. Operator ordering, constraint anomalies, domains and a physical inner product remain open.
- A regular spatial interpretation and local cut realization have been selected. Global sheet/winding conditions and the physical 3+1 interacting extension are not proved.
- Neither LCF's new spatial coframe nor this local theory inherits AGD/ZCR's global-mode quantum completion or recovery estimates without a separate argument.

The next concrete target is a physical reduction and quantization of this specified local model, with potential and boundary conditions fixed, testing whether the classical algebra survives and whether a positive physical sector exists. Classical closure alone cannot settle those questions or uniquely fix a potential.

## Certificate scope

The native rational engine checks polynomial functional-derivative identities, total derivative boundary cancellation, matter differential-operator brackets, source cancellation, coefficient failures and potential underdetermination. The native rewrite engine emits serialized algebraic witnesses and rejects a corrupted witness. Python 3.12 handles I/O only. General smooth-field and preservation statements have written proofs; finite certificate counts do not turn them into proof-assistant-checked theorems.
