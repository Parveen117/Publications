# GE1: domains, stable exponential change and the native heat dock

3 October 2026. Written proofs plus exact finite controls; Python 3.12 only.

## Premise order and scope

Use the completed cut-complex field earned by RKF F00-E/F00-G, its central
quarter turn `iota`, positive modulus from the cut square, and its factorial
exponential. Finite positive pairings, their null quotients and Cauchy
completion are used in the same order as RH T01 and YM50. No pre-existing
infinite Hilbert space, Haar integral, spectral theorem or abstract semigroup
generation theorem supplies a proof below. A Hilbert-space coordinate
description can be made afterwards; its familiar name is not a new primitive.

This is a scoped implementation of familiar evolution/domain mathematics,
not a priority claim for Fourier modes, Cayley transforms, implicit Euler,
closed generators or near resonance. For the established approximation
literature see A. Gomilko and Y. Tomilov, *On rates in Euler's formula for
C0-semigroups*, [arXiv:1301.4406](https://arxiv.org/abs/1301.4406). That paper's
functional calculus is not used as an input to the elementary proofs here.

There are two different carriers below. The Generalized Euler character
module is declared explicitly. The YM application uses the already derived
quaternion coefficient carrier and its counted reference; it does not identify
that noncommuting carrier with a torus. Seam, ratio, winding and branch-memory
coordinates of the larger Generalized Euler lift are not erased by this
restricted phase-sector construction, nor are their dynamics certified here.
The RH source's existing multiplicative Euler-scale flow is a third carrier.
Its weighted energy is not generally preserved by dilation; see the audited
growth bound in SOURCE_AUDIT.md. An isometry result on one pairing cannot be
transferred to a different weighted energy merely by changing notation.

## GE1-T1. A character core earns its completion and real generator

Let `J` be finite or countable and `K = Z^(J)` the integer sequences with
finite support. Begin with finite formal Laurent sums

\[
 f=\sum_{k\in K}f_k z^k,\quad z^kz^l=z^{k+l},\quad
 (z^k)^\dagger=z^{-k},\quad
 \phi(f)=f_0.
\]

Thus `phi(f^dagger f)=sum_k |f_k|^2` by finite coefficient extraction.
The finite cut-square inequality bounds the pairing. Complete this energy
norm to obtain `H_GE`: explicitly, cut-complex sequences with summable
squared coefficients. Completeness follows by coordinatewise limits and a
finite-tail estimate, as in RH T01. The Laurent algebra acts on this module:
each `z^k` shifts coefficients isometrically, so multiplication by a finite
`p` has bound `sum |p_k|`. Products of two arbitrary completed vectors are
**not** thereby declared to belong to the completion.

Choose real frequencies `omega_j`; their physical selection is not derived.
For a finite-support `k`, put `alpha_k=sum k_j omega_j`. Define

\[
 (U(t)f)_k=\operatorname{Exp}_\Sigma(\iota t\alpha_k)f_k,
 \qquad t\in\widehat{\mathbb R}_\Sigma.
 \tag{1}
\]

**Theorem.** `U` is a strongly continuous group of matching isometries.
Its strong derivative generator is exactly

\[
 (Df)_k=\iota\alpha_k f_k,\qquad
 \mathcal D(D)=\{f:\sum_k\alpha_k^2|f_k|^2<\infty\}.
 \tag{2}
\]

It is closed and skew-adjoint for the earned pairing. Finite Laurent sums
are a graph core: truncation approximates both `f` and `Df`.

**Proof.** The native exponential addition and dagger laws give composition
and modulus one. For strong continuity, truncate to finitely many coefficients;
their continuity is scalar, while the remaining norm is at most twice the
tail norm. For `f` in (2), the scalar bound
`|Exp(iota t alpha)-1| <= |t alpha|` and finite-tail splitting justify
the difference quotient and its limit. Conversely a strong derivative has
coordinate `iota alpha_k f_k`; its norm is finite, so (2) is necessary.
Coordinate limits prove closedness. Testing an adjoint candidate on every
`z^k` gives the same square-summability condition and adjoint `-D`.
Truncation proves the graph-core statement. QED.

The finite-dimensional torus chart is `z^k(theta)=Exp(iota k.theta)`.
There `D=sum omega_j partial_theta_j` and (1) is translation. The algebraic
construction works also for countably many phase labels without asserting
that every informal smooth function on an unspecified infinite torus has
already been included.

## GE1-T2. Complex directions, analytic vectors and the failed forward step

Allow `omega=a+iota b` for this paragraph only. Write `beta_k=k.b` and
`alpha_k=k.a`. The diagonal exponential has maximal domain

\[
 \mathcal D(U_\mathbb C(t))=
 \{f:\sum_k e^{-2t\beta_k}|f_k|^2<\infty\}.
 \tag{3}
\]

On the full bilateral character module it is bounded for a nonzero real
`t` **if and only if `b=0`**. Indeed a nonzero component of `b` gives
arbitrarily growing multipliers along one of `k=+m e_j` or `-m e_j`.
Formula (3) still defines a closed densely defined operator, by coordinate
limits; it is not an everywhere-defined diffusion semigroup.

For the change operator `C=I+D_omega`, on its diagonal maximal domain,
`ker C` consists of square-summable coefficients supported on
`k.omega=iota`. On a smooth finite-torus chart the analogous coefficients
are rapidly decreasing, not necessarily finite sums. In particular
`omega=iota` gives the mode **`z^1`**, not `z^-1`. On this kernel the
admitted exponential is multiplication by `e^-t`. This is one selected
eigenspace, not decay of all states.

The formal Taylor series `sum t^m D^m f/m!` is valid in norm when
`sum |t|^m ||D^m f||/m!` is finite. Membership in all graph-power domains
alone is insufficient. For `J={1}`, `omega=1`, let

\[
 f_k=2^{-\lceil\sqrt{k}\rceil}\ (k\geq1),\qquad f_k=0\ (k\leq0).
 \tag{4}
\]

Every power-weighted coefficient sum converges, so this is smooth in the
torus chart and belongs to every `D^m` domain. At `k=m^2`, the norm of
the `m`th Taylor term is at least

\[
 |t|^m m^{2m}2^{-m}/m!\ \geq (|t|m/2)^m.
 \tag{5}
\]

For every nonzero `t` these terms fail to tend to zero, although (1)
is well defined and norm preserving.

Forward Euler has a separate failure. At cutoff `k<=n`, time `t=1`,
take `n=m^2`. The `k=n` coefficient of `(I+D/n)^n P_n f` has squared
modulus

\[
 2^{m^2-2m}\longrightarrow\infty.
 \tag{6}
\]

Here `f` is the **same** vector (4) for all cutoffs. Thus convergence for
each fixed finite matrix does not justify a joint cutoff/step limit, even
on smooth data. A cutoff projection in this statement keeps a finite
interval of characters including `1,...,n`.

## GE1-T3. A stable reversible exponential approximation

For real frequencies and `h` define, using (2),

\[
 V_h=(I-hD/2)^{-1}(I+hD/2),\qquad
 (V_h f)_k=\frac{1+\iota h\alpha_k/2}{1-\iota h\alpha_k/2}f_k.
 \tag{7}
\]

The right-hand formula defines the bounded extension on every vector;
the unbounded product on the left is initially interpreted on `D(D)`.
`V_h` is an isometry, has inverse `V_-h`, and commutes with character
cutoffs. For `f` in `D(D^3)`,

\[
 \boxed{\|V_{t/n}^{\,n}f-U(t)f\|
 \leq\frac{|t|^3}{12n^2}\|D^3f\|.}
 \tag{8}
\]

If `P_N` is any increasing finite-character exhaustion, then for every
`f in H_GE`, `V_(t/n)^n P_N f -> U(t)f` jointly as `N,n -> infinity`,
without a rate relation between them. This is **strong** convergence,
not operator-norm convergence over arbitrary high frequencies.

**Proof.** Modulus one in (7) proves stability. On a finite core,
`V_s'=D(I-s^2D^2/4)^-1 V_s`; the inverse has matching gain at most one.
Consequently
`||(V_s'-D V_s)f|| <= s^2 ||D^3 f||/4`.
Differentiate `U(h-s)V_s f` and integrate the bound to get
`||(V_h-U(h))f|| <= |h|^3 ||D^3 f||/12`.
Telescoping `n` isometries proves (8); graph truncation extends it to
`D(D^3)`. For arbitrary `f`, fix a finite core vector `p` whose norm
distance is at most `epsilon`. For any cutoff containing `p`, the error
is bounded by `2 epsilon+|t|^3||D^3p||/(12n^2)`. First choose `p`, then
`n`. This is a joint-limit proof, not substitution of a cutoff-dependent
large derivative bound into (8). QED.

## GE1-T4. A finite native core determines a closed heat generator

This theorem applies directly to the YM carrier and is independent of the
commuting phase model. Let `V_0 subset V_1 subset ...` be finite positive
pairing spaces with consistent inclusions, and let `H` be the completion
of their union. Suppose `L` preserves every `V_d`, is symmetric there
and obeys `<p,Lp> >= 0`. Null relations, if present in a polynomial
presentation, are first quotiented out. Set `E_d(t)=Exp_Sigma(-t L_d)`.

**Theorem.** The following constructions agree on the dense core and
extend consistently:

1. `E(t)` is a strongly continuous contraction semigroup on `H`.
2. `L` has a closed graph extension `A=closure(L)`, and this extension
   is self-adjoint and nonnegative for the earned pairing.
3. For every `h>0`, the finite inverses `(I+hL_d)^-1` extend to the
   contraction `B_h=(I+hA)^-1` on all of `H`.
4. The strong generator of `E(t)` is exactly `-A`, with domain `D(A)`.
   Thus its domain is specified, not inferred from a formal power series.

**Proof.** Finite factorial tails construct `E_d`, and differentiating
its squared norm gives `-2<E_d p,L_d E_d p><=0`. This proves contraction
without diagonalizing. Finite uniqueness and invariance make the actions
consistent. Core approximation gives their bounded extension, composition
and strong continuity.

For `y` in a finite core,
`||y||^2 <= <y,(I+hL)y> <= ||y|| ||(I+hL)y||`.
Finite elimination supplies the inverse and its contraction bound.
Its inverses on nested spaces agree. Symmetry makes `L` closable:
if `p_j->0` and `Lp_j->q`, testing against any core vector forces `q=0`.
For `x_j` in the core converging to `x`, set `y_j=B_h x_j`.
Then `y_j->B_h x` and `Ly_j=(x_j-y_j)/h` converges. Hence `B_h x`
belongs to the graph closure and `(I+hA)B_h x=x`. Passing the reverse
identity along graph approximations proves `B_h(I+hA)y=y` on `D(A)`.

Symmetry and nonnegativity persist along graph limits. If `x` belongs
to the adjoint domain with `A^*x=z`, set `y=B_h(x+hz)`. Testing against
`(I+hL)p` shows `<(I+hL)p,x-y>=0` for every core `p`. The range of
`I+hL` on the core is the entire core. Thus `x=y`, proving
`D(A^*)=D(A)` and self-adjointness without an external representation
or generation theorem.

On the core, `E(t)p-p=-integral_0^t E(s)Lp ds`. The integral is the
native count-sum limit of a continuous curve. Graph approximation extends
the identity to `D(A)`, so its strong derivative is `-A`. Conversely,
if `f` has a strong derivative `g` at zero, test its difference quotients
against a core `p`. Symmetry of `E(t)` yields `<p,g>=<-Lp,f>`.
The resulting adjoint-domain test and self-adjointness force
`f in D(A)` and `g=-Af`. QED.

## GE1-T5. Positive implicit exponential with a certified error

Under GE1-T4, for `t>=0`, `n>=1`, and `f in D(A^2)`,

\[
 \boxed{\|B_{t/n}^{\,n}f-E(t)f\|
 \leq\frac{3t^2}{2n}\|A^2f\|.}
 \tag{9}
\]

At `t=0` both sides are interpreted with `B_0=I`. For general `f`,
and any finite core `p`, the right-hand side can be replaced by
`2||f-p||+3t^2||L^2p||/(2n)`. Thus convergence is strong on the
whole completion. Orthogonal projections `P_N` onto the nested `V_N`
satisfy `B_(t/n)^n P_N f -> E(t)f` jointly, with no cutoff/step relation.
The numerical constant is a sufficient bound, not claimed optimal.

**Proof.** On the core the exact inverse identity gives
`B_h f=f-hLf+h^2 B_h L^2f`. Twice integrating the finite factorial
derivative gives
`E(h)f=f-hLf+integral_0^h(h-s)E(s)L^2f ds`.
Their difference has norm at most `3h^2||L^2f||/2`. Telescoping
contractions and using commutation on the invariant core gives (9).

To extend the bound to `D(A^2)`, note that the orthogonal complement
of each finite `V_N` is invariant under the symmetric resolvents and
`E(t)`. Therefore `P_N A^j f=A^j P_N f` for `j<=2` on the indicated
domains, by testing against the finite core. The three tails tend to
zero. Finally use the same fixed-core `2 epsilon` argument as in T3.

There is also the resolvent identity

\[
 \boxed{B_h f=\int_0^\infty e^{-s}E(hs)f\,ds.}
 \tag{10}
\]

Construct it by finite native sums on `[0,S]`; the omitted norm is at
most `e^-S ||f||`. On a finite core, integration by parts in
`e^-s E(hs)f` gives `(I+hL)` times (10) equal to `f`.
Contraction and density extend the identity. If the existing heat action
preserves a specified function order, constants and a bounded reference
functional, (10) preserves the same structures. On YM's uniform function
completion these statements follow directly from positive sums and their
uniform limit. Positivity is not inferred from a coefficient-space norm.

Equations (9)--(10) do **not** select `L`, its protocol, its energy scale
or a clock. They evaluate the evolution of a law that has already been
specified. In particular one cannot replace `D` by a Laplacian and then
call that substitution a derivation of heat from first-order drift.

The reversible midpoint formula from T3 must not be mistaken for a
positive heat step. For the three-record generator with diagonal entries
2 and off-diagonal entries -1, `(I+2L)^-1(I-2L)` contracts the matching
norm but has diagonal entry `-1/7`. The positive resolvent `(I+4L)^-1`
has nonnegative entries and row sums one. Exact inversion verifies both;
norm stability alone is not a positivity certificate.

## GE1-T6. Time averages earn a projection, but not a Born law or gap

For the real character flow of T1 define the native count integral
`M_T f = T^-1 integral_0^T U(t)f dt`. For `T>0`, its multiplier is

\[
 m_T(k)=\begin{cases}1,&\alpha_k=0,\\
 (e^{\iota T\alpha_k}-1)/(\iota T\alpha_k),&\alpha_k\ne0.
 \end{cases}
\]

It contracts the matching norm and converges strongly to the coefficient
projection `P_0` onto `alpha_k=0`. For any finite core `p`,

\[
 \|M_T f-P_0 f\|\leq2\|f-p\|+
 \frac2T\left(\sum_{\alpha_k\ne0}\frac{|p_k|^2}{\alpha_k^2}\right)^{1/2}.
 \tag{11}
\]

**Proof.** The integral average has norm at most one. Integrate each
finite scalar exponential and bound its numerator by two. Finite-tail
approximation proves (11) and strong convergence. QED.

An actual normalized vector and an independently specified orthogonal
projection yield nonnegative weights `||P f||^2/||f||^2`; that algebraic
fact alone does not derive a detector, a measurement law, or the draft's
identification of winding density with probability. Constants belong to
`ker D` even when the flow vector is nonzero. The draft's equation
`ker D = {v=0}` confuses functions with vector fields.

For a decisive mass-gap control, choose two phase frequencies `(1,sqrt(2))`
in the earned radial field. On a character labelled `(p,-q)`, the
separately selected heat generator `L=-D^2` has eigenvalue
`lambda=(p-q sqrt(2))^2`. Its kernel consists only of constants, since
integer squares cannot solve `p^2=2q^2` nontrivially (parity descent).

Define positive integer pairs by

\[
 (p_1,q_1)=(1,1),\quad
 (p_{m+1},q_{m+1})=(p_m+2q_m,p_m+q_m).
\]

Then `p_m^2-2q_m^2=(-1)^m`, so

\[
 \boxed{0<\lambda_m=(p_m-q_m\sqrt2)^2
 =\frac1{(p_m+q_m\sqrt2)^2}\longrightarrow0.}
 \tag{12}
\]

Every finite nonconstant character cutoff has a positive minimum rate,
but the full centered completion has **no positive uniform decay gap**.
The chosen characters have norm one and mean zero; their heat decay
is `e^(-t lambda_m)`, contradicting any proposed fixed positive rate.
This is a counterexample to the inference from finite positivity, not a
claim that every commuting multi-direction heat law lacks a gap.

## Actual YM application and limits

Take `H` to be YM50's counted-reference recognition completion and `V_d`
the restrictions of coefficient polynomials of degree at most `d`, modulo
their null relations. For any fixed protocol tensor `C>=0`, YM51 derives

\[
 L_C=-\sum_{a,b}C_{ab}D_aD_b,\qquad
 \Phi_Q(\bar p L_Cp)=\sum_r w_r\Phi_Q(|D_{v_r}p|^2)\geq0.
\]

It preserves polynomial degree and is symmetric. Therefore T4 and T5
apply to this **actual** generator, including singular `C` and
noncommuting `D_a`. They identify its graph closure, prove the positive
resolvent bridge and supply stable rational approximation. The heat action
agrees with YM51 on polynomials and hence on the completion by density.

For `C=diag(0,1,1)`, the existing YM52 compact-sector gap is `1/2`.
It is not derived anew here. Rational degree-two controls also use a
rotated rank-two tensor, so the verification does not replace the YM
engine by diagonal character multipliers. On finite raw homogeneous
coefficient spaces the already derived tensor energy provides an
independent exact contraction check before the reference null quotient.

T1--T3 treat constant commuting phase directions. Their brackets vanish.
The nonzero native YM brackets come from the quaternion carrier, as
already derived in YM50--52. Connection curvature additionally requires
its connection and frame-subtraction law (EMK-C1); none of these objects
is silently identified with another.

This closes neither the actual interacting infinite-chain row-closure
defect nor spatial lattice removal, four-dimensional Yang--Mills,
asymptotic freedom, Clay's mass gap, RH, or quantum gravity. It does not
replace the transfer steps or alter the constants of YM53--55. A future
use of these approximations inside interacting normalized histories must
separately propagate its error through those histories.
