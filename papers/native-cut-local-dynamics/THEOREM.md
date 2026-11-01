# CLD-1 — Dynamics after the observation cut

## Meaning, lineage and scope

Here **local means relative to the retained information sector after an observation cut**. It is an operational definition, not an assertion about separation of spatial points. A parameter `s` orders operations; identifying it with laboratory time needs a separate calibration. First we give a law with no differentiable parameter at all.

This packet extends, rather than renames, established results: [GE3](../generalized-euler-evolution/GE3_RETURNING_MEMORY.md) already eliminates a fixed hidden block; [NJC-T3](../native-jacobian-cut-certification/THEOREM.md) already constructs the cut-compatible connection; [ZCR-1](../native-zero-reference-recovery/THEOREM.md) already proves fixed-cut leakage and retained error estimates. The new combination tracks changing observation cuts, gives the corresponding closure defect and conditioning budget, and evaluates the defect for moving reflection cuts of the actual [AGD-1](../native-action-geometry-dynamics/THEOREM.md) operator. [IRT-1](../native-invariant-ratio-transport/THEOREM.md)'s stationary reflection invariant is recovered unchanged.

Sections 1–4 concern finite carriers. Section 5 separately proves a common-core identity and a closure criterion in the infinite AGD model. Finite norm estimates are **not** asserted for its unbounded moving-cut connection. Written analytic arguments and native algebraic certificates have different scopes; source hashes and finite examples do not machine-prove the analytic arguments.

## 1. A cut-relative law without a clock

Let `T_j : H_j -> H_{j+1}` be unitary arrows and let `P_j` be orthogonal observation cuts, `Q_j=I-P_j`. Define the typed blocks

\[
 a_j=P_{j+1}T_jP_j,\quad b_j=P_{j+1}T_jQ_j,\quad
 c_j=Q_{j+1}T_jP_j,\quad d_j=Q_{j+1}T_jQ_j.
\]

For `x_j=P_j psi_j`, `y_j=Q_j psi_j`, exact evolution is

\[
 x_{j+1}=a_jx_j+b_jy_j,\qquad y_{j+1}=c_jx_j+d_jy_j.
\]

Successive substitution, with empty products the identity on the appropriate hidden fiber, gives

\[
 x_{n+1}=a_nx_n+b_nd_{n-1}\cdots d_0y_0
 +\sum_{j=0}^{n-1}b_nd_{n-1}\cdots d_{j+1}c_jx_j. \tag{1}
\]

This is GE3's memory identity with the domains now changing with the cuts. Its proof is induction in the hidden-block recurrence; no physical time or exponential hypothesis is needed. Preparation in the first retained sector sets `y_0=0`; it does not remove later return paths.

In particular,

\[
 P_2T_1T_0P_0-(P_2T_1P_1)(P_1T_0P_0)
 =P_2T_1Q_1T_0P_0. \tag{2}
\]

For equal finite retained ranks, `c_j=0` is equivalent to
`T_j P_j T_j^*=P_{j+1}` and then `b_j=0`. Indeed the isometry maps the first retained range into a range of the same dimension, hence onto; unitarity maps the orthogonal complements onto one another. This is exact autonomous cut-relative transport. Equal rank is essential to this implication as stated.

**Exact example.** Set `T_0=T_1=I`, `P_0=P_2=diag(1,0)`, and `P_1=u u^*`, `u=(3/5,4/5)`. Equation (2) separates endpoint amplitude `1` into retained-path amplitude `9/25` and return amplitude `16/25`. With no intermediate physical measurement the endpoint survival is `1`. Actually selecting the retained outcome at the intermediate and final measurements has joint success probability `81/625`. The normalized final retained state is the same, but the outcome ledger is different. Merely writing a block decomposition is not a measurement intervention.

## 2. The moving-cut closure theorem

Let `P(s)` be a continuously differentiable finite orthogonal projection on a compact parameter interval, `H(s)=H(s)^*` continuous, and `eta>0` constant. Write

\[
 K=[P',P],\qquad W'=KW,\quad W(0)=I,\qquad P_0=P(0).
\]

Differentiating `P^2=P` gives `P'P+PP'=P'`, hence `PP'P=QP'Q=0`. It follows that `K^*=-K` and `[K,P]=P'`. Uniqueness of the linear matrix equation then gives `P(s)=W(s)P_0W(s)^*` and `W^*W=I`.

This is the existing NJC connection in a forward-transport convention: for `J=2P-I`, its minimal covariant coefficient is `JJ'/2=[P,P']=-K`. It is not a second independently derived connection.

For `i eta psi'=H psi` and `z=W^*psi`,

\[
 i\eta z'=\widetilde H z,\qquad
 \widetilde H=W^*(H-i\eta K)W. \tag{3}
\]

The connection term is Hermitian and its sign follows directly by differentiating `W^*psi`. Define

\[
 \boxed{\mathfrak D=P'+\frac{i}{\eta}[H,P]},\qquad
 b=Q_0\widetilde H P_0
   =-i\eta W^*Q\mathfrak D P W. \tag{4}
\]

`D` is Hermitian and offdiagonal. To verify the second identity use
`QKP=QP'P`, `Q[H,P]P=QHP`. Consequently the following are equivalent:

1. `D=0` throughout the interval.
2. `P'=-(i/eta)[H,P]`, so `P(s)=U(s,0)P_0U(s,0)^*` for the full propagator.
3. `b=0`, and the moving retained and hidden sectors evolve independently.
4. `Tr(P(s)rho(s))` is constant for every initial density matrix.

For the last equivalence, the density equation gives

\[
 \frac{d}{ds}\operatorname{Tr}(P\rho)=\operatorname{Tr}(\mathfrak D\rho). \tag{5}
\]

Unitary evolution reaches every density matrix at each parameter value. A Hermitian matrix whose expectation vanishes on all pure states is zero. Conservation for one specially prepared state does not establish this all-state conclusion.

The gauge change `W -> W Z`, with unitary differentiable `Z` commuting with `P_0`, gives
`Htilde -> Z^* Htilde Z - i eta Z^* Z'`. Its additional term is block diagonal, and `b -> Z_Q^* b Z_P`. Thus the vanishing and singular values of the cut-motion defect are independent of the frame chosen inside each retained/hidden fiber.

Choosing the cut to follow the full propagator guarantees closure but uses that propagator as input. This is a transport identity, not a new source equation or a prediction of gravity.

## 3. Return memory and a conditioning budget

Set `a=P_0 Htilde P_0`, `d=Q_0 Htilde Q_0`. For `x=P_0z`, `y=Q_0z`,

\[
 i\eta x'=ax+b^*y,\qquad i\eta y'=bx+dy.
\]

If `U_d(s,t)` is the unitary propagator of `d`, variation of constants gives

\[
 x'=-\frac i\eta ax-\frac i\eta b(s)^*U_d(s,0)y_0
 -\frac1{\eta^2}\int_0^s b(s)^*U_d(s,t)b(t)x(t)\,dt. \tag{6}
\]

The last term is returned information. It survives even after a selected initial preparation `y_0=0`. This is the changing-cut extension of GE3/ZCR, not an assumption of Markov evolution in the retained sector.

For a normalized initial retained pure state let `x_iso=U_a(s,0)x_0` and

\[
 \Lambda(s)=\int_0^s\frac{\|b(t)\|}{\eta}\,dt
 =\int_0^s\|Q(t)\mathfrak D(t)P(t)\|\,dt.
\]

Unitarity and the two variation-of-constants formulas give

\[
 \|y(s)\|\le\Lambda(s),\quad
 w(s)=\|x(s)\|^2\ge1-\Lambda(s)^2,\quad
 \|x(s)-x_{iso}(s)\|\le\tfrac12\Lambda(s)^2. \tag{7}
\]

Indeed `||x(t)||<=1` bounds the hidden integral by `Lambda`; the retained difference is bounded by `integral Lambda(t) dLambda(t)=Lambda(s)^2/2`.
For any bounded retained observable `O`,

\[
 |\langle x,Ox\rangle-\langle x_{iso},Ox_{iso}\rangle|
 \le\|O\|\Lambda^2. \tag{8}
\]

If `Lambda<1`, the normalized conditional reading satisfies

\[
 \left|\frac{\langle x,Ox\rangle}{w}
 -\langle x_{iso},Ox_{iso}\rangle\right|
 \le\frac{2\|O\|\Lambda^2}{1-\Lambda^2}. \tag{9}
\]

Proof: divide the numerator error plus `(1-w)` times the isolated expectation by `w`, and use (7)–(8). The same survival and observable bounds hold for mixed initial states supported in `P_0`, by convex decomposition before conditioning. The amplitude bound in (7) is a pure-state statement.

Example: `Lambda<=1/10` guarantees survival at least `99/100`, retained amplitude error at most `1/200`, unconditioned unit-norm observable error at most `1/100`, and conditional error at most `2/99`. Normalization must not silently erase the success weight. This conditional state update is not an autonomous nonlinear evolution of an unknown unmeasured state (QGC-1).

## 4. A moving detector is not necessarily a changing source

In two dimensions let `R=[[0,-1],[1,0]]`, `W=exp(omega s R)` and `P=W P_0W^*`. For `H=0`, the source vector is constant, but (3) gives `Htilde=-i eta omega R`. Initial vector `(1,0)` has moving-frame coordinates `(cos(omega s),-sin(omega s))`; retained weight is `cos^2(omega s)`. At a frame with cosine `3/5` and sine `4/5`, this is `9/25`. Omitting the connection predicts `1` and is rejected.

If instead `H=i eta omega R`, source and cut co-move, `Htilde=0`, and weight remains `1`. Thus a changing reading can originate in the observation cut; it is not by itself evidence of source loss or curvature.

## 5. Application to the full AGD geometry

Use the AGD conventions `L=-i partial_theta`, `p=eta L`,

\[
 H_g=p^2+8-2\cos2\theta,\quad r=(1-\cos\theta)/2,\quad
 v=1+a f(x)r,\quad a=3/5,\quad f=x^2/(1+x^2),
\]
\[
 H=H_g+i\eta C(v\partial_x+v_x/2)+i\eta gRv,
 \qquad C=\operatorname{diag}(1,-1).
\]

Here `R` acts on the matter direction register. IRT proves that reflection `S phi(theta)=phi(-theta)` commutes with this full operator, and `Pi_+=(I+S)/2` reduces it. Consider a continuously differentiable rigidly rotated detector

\[
 W_\alpha=e^{-i\alpha(s)L},\quad P_\alpha=W_\alpha\Pi_+W_\alpha^*,
 \qquad \alpha(0)=0.
\]

On the smooth periodic, spatially compact core used in AGD,
`W_alpha^* F(theta) W_alpha=F(theta+alpha)` and

\[
 \widetilde H_\alpha=W_\alpha^*HW_\alpha-\alpha' p. \tag{10}
\]

The forward connection `W'W^*=-i alpha' L` is odd under reflection, hence is the minimal cut connection on this core. It is unbounded. `P_alpha` is not assumed norm-C1; Sections 2–3 are not invoked as unbounded operator theorems.

Write

\[
 H_f=i\eta\{C(f\partial_x+f'/2)+gRf\}.
\]

Parity separates even and odd terms of (10), giving the exact common-core identity

\[
 \boxed{
 b_\alpha=\Pi_-\left[
 2\sin(2\alpha)\sin(2\theta)
 +\frac a2\sin\alpha\sin\theta\,H_f
 -\alpha' p\right]\Pi_+ .} \tag{11}
\]

Proof: the odd part of `-2cos(2theta+2alpha)` is `2sin(2alpha)sin(2theta)`; the odd part of `r(theta+alpha)` is `sin(alpha)sin(theta)/2`. The kinetic term and all remaining terms are even. Momentum is odd. Sandwiching between opposite parity sectors removes every even term. This proves (11) directly for the infinite operator; finite Fourier fixtures below only check the implementation.

**Closure within this detector family.** At a parameter value, (11) vanishes on every even core state if and only if

\[
 \sin\alpha=0,\qquad \alpha'=0. \tag{12}
\]

To prove necessity choose an interval in `x` on which `f` is nonzero, a nonzero smooth compactly supported profile `chi`, a nonzero matter spinor, and the constant even angular mode. Replace `chi(x)` by `exp(i k x)chi(x)`. In (11) the only term growing with `k` is the spatial derivative coefficient
`(a/2)sin(alpha)sin(theta) i eta C f partial_x`.
Its norm after division by `|k|` has a nonzero limit unless `sin(alpha)=0`. No angular multiplication or momentum term cancels this coefficient. With `sin(alpha)=0`, both trigonometric terms vanish; testing any nonconstant cosine mode forces `alpha'=0`. Sufficiency is immediate. Over a connected interval the closed rigid detector is a constant integer-pi rotation, equivalent to the original reflection cut. This is not a classification of all possible moving cuts.

**Initial leakage coefficient.** Let `alpha'(0)=omega` and initially prepare normalized `cos(n theta)` (`n>=1`) times a normalized smooth matter profile. Strong differentiation of `W_alpha(s)^* exp(-iHs/eta)` on this initial core vector gives
`y'(0)=i omega L psi_0`, since `Pi_- H Pi_+=0`. Thus

\[
 1-w(s)=\omega^2 n^2s^2+o(s^2). \tag{13}
\]

Only differentiability at zero is required: `y(s)=s y'(0)+o(s)` in norm. For `n=1`, `omega=3/5`, the coefficient is `9/25`. It is a short-parameter asymptotic coefficient, not a finite-step probability and not a fitted constant. It describes detector-induced sector change even though IRT's fixed reflection weights remain conserved.

## 6. Certificates and remaining obligation

The unchanged native rational matrix engine checks moving-frame signs, block covariance, coherent return versus actual intermediate selection, exact budgets, Fourier reflection identities and (11) in symmetric finite sections. The unchanged native rewrite engine proves projection-derivative and reflection identities with serialized witness replay. Corrupted witnesses and stronger false conclusions are rejected. `SOURCE_PINS.json` binds the manuscript, verifier, runner and reused sources; `EXPECTED.sha256` binds the generated certificate.

The infinite proof in Section 5 uses principal derivative separation and a smooth-core differentiation argument; a finite truncation is not its proof. No certified numerical experiment establishes the physical interpretation of `eta`, the lapse profile or the chosen AGD action.

**What this closes:** the meaning and exact evolution law of cut-local information, the extra term required when the cut moves, the all-state closure criterion in finite dimension, and the moving-reflection criterion in AGD.

**What remains for quantum gravity:** deriving and preserving a spatial local coframe/source constraint algebra (including AGD's residual `g R N f'`), a controlled physical interacting limit and physical calibration. Operational cut-locality alone does not settle spatial locality, refoliation or Einstein dynamics. The next construction must keep the cut-motion defect explicit while addressing those independent obligations.
