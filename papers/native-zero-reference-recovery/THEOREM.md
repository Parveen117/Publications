# ZCR-1 — Fixed reference, closed observation and controlled metric recovery

15 October 2026. Continuation of AGD-1 using the exact sector criterion already proved in GCP-1.

**Result.** A reference called zero does not select classical dynamics. The useful invariant conditions are vanishing observation-sector leakage and a scalar geometry reading on that sector. Their exact version is inherited from GCP-T5. This packet supplies a finite-carrier error budget and, more substantially, a normalizable finite-width family in the actual variable-profile AGD model whose matter evolution approaches the flat-metric law with an explicit error. The geometry is allowed to evolve throughout the comparison. The limit changes the dimensionless action scale and preparation together; it is not a relabelling of zero or a unit choice.

## 1. Reuse and scope

| Existing result | Reuse here |
|---|---|
| GCP-T4/T5 | Geometry variance is a compression defect; an invariant scalar-lapse sector gives exact metric reduction |
| GE3-T3 | Eliminating a retained complement produces a return-memory term |
| GE4-T3 | A generator residual controls the difference between full and reduced evolution |
| AGD-T2/T3/T4 | Infinite circle geometry, self-adjoint joint evolution, and spatial locality |
| AGD section 6 | Nonconstant circle lapse has no normalizable sharp-geometry eigenstate |

The new main result is ZCR-T3 below: a quantitative finite-time recovery bound on the full spatially varying AGD operator. The algebraic lemmas and Duhamel argument are adaptations of existing machinery, not claims of new general projection or perturbation theory. All quantities remain dimensionless. The physical calibration of the action scale, clock and lapse dictionary is not derived here.

## 2. What it means for an observation sector to close

First take a finite positive carrier, a Hermitian generator `H`, an orthogonal retained cut `P`, complement `Q=I-P`, and `eta>0`. Set

\[
A=PHP|_P,\quad D=QHQ|_Q,\quad B=QHP,
\quad U(s)=\exp(-\iota sH/\eta).
\]

The exact identity is

\[
\boxed{PH^2P-(PHP)^2=B^\dagger B.}                         \tag{1}
\]

It is the inherited retained-record defect with `H` as the tested operation. Consequently `B=0`, `[H,P]=0`, and
`PU(s)P|_P=exp(-iota sA/eta)` for all real `s` are equivalent. For the reverse implication, differentiate twice at zero and use (1). Hermiticity turns `QHP=0` into a reducing, rather than just one-sided invariant, sector.

For an unbounded self-adjoint operator, the appropriate exact hypothesis is that `P` reduces its unitary group (equivalently its spectral resolution). A formal commutator on an unspecified core is not enough. We do not extend the finite block norm estimates below to arbitrary unbounded blocks.

To obtain a definite metric as well, require the geometry observable to act scalarly on the geometry sector `S`:

\[
V(x)|_S=v_{\rm cl}(x)I_S,\qquad H_g S\subseteq S
\quad\text{with a reducing-domain interpretation}.        \tag{2}
\]

GCP-T5 then supplies factorization into scalar-metric matter propagation and internal geometry evolution. An arbitrary closed sector alone does not imply a scalar metric. Conversely an initially sharp lapse alone does not ensure closure under the generator. The recovered scalar metric can be curved; zero geometry dispersion is not the same statement as zero curvature.

### A change of zero leaves the obstruction unchanged

For real scalars `b,d`,

\[
Q(H+bI)P=QHP,\qquad
(V-dI-(v_{\rm cl}-d)I)P=(V-v_{\rm cl}I)P.                 \tag{3}
\]

Energy-origin shifts only supply a global phase. The geometry equation in (3) changes both the reading and its reference. Changing `V` in the physical propagation law without changing that dictionary would change the model instead of its labels.

A two-state counterexample is `H=Z`, `rho=P_0=(I+Z)/2`. The density is stationary, but `Var_rho(X)=1` and `[X,Z]!=0`. A fixed quantum state need not make its observable algebra classical. With `H=bX`, `V=diag(1,2)` and the same initial `P_0`, the lapse is initially sharp but the sector leaks: `B†B=b²P_0` and the loss probability is `sin²(bs/eta)`.

## 3. Finite-carrier leakage and readout error

**ZCR-T1 (closed-sector error budget).** In the finite setting of section 2, let `ell=||B||/eta`. For normalized initial data in `P`, and real `s`,

\[
\|QU(s)P\|\le \min(1,\ell|s|),
\]
\[
\boxed{\|PU(s)P|_P-e^{-\iota sA/\eta}\|
\le\tfrac12\ell^2s^2.}                                   \tag{4}
\]

The leakage probability is at most `min(1,ell²s²)`. For a bounded retained observable `O=POP`, its **unconditioned** full reading and compressed reading differ by at most

\[
\|O\|\ell^2s^2.                                          \tag{5}
\]

These are upper bounds and may exceed trivial bounds at large times. Conditioning on survival would require an additional probability denominator.

**Proof.** Write `x=PU(s)psi`, `y=QU(s)psi`. The block equations are
`iota eta x'=Ax+B†y`, `iota eta y'=Bx+Dy`, with `y(0)=0`. Variation of constants and `||x||<=1` give `||y(s)||<=ell|s|`. Integrating the retained equation relative to `exp(-iota sA/eta)` then gives `ell integral_0^|s| ell t dt`, proving (4). Since both the exact retained vector and its compressed comparator have norm at most one, the difference of their quadratic readings is at most `2||O||` times their vector difference, which proves (5). QED.

The eliminated equation also gives

\[
x'=-\frac{\iota}{\eta}Ax-
\frac1{\eta^2}\int_0^s B^\dagger e^{-\iota(s-t)D/\eta}B x(t)\,dt.
                                                               \tag{6}
\]

Thus the error term is the retained return information, not an adjustment of the reference value. Equation (6) is the unitary version of the already retained GE3 memory mechanism.

## 4. A finite-width circle preparation with exact moments

Return to the infinite circle carrier of AGD. Write

\[
p=-\iota\eta\partial_\theta,\quad
r(\theta)=\frac{1-\cos\theta}{2},\quad
H_g=p^2+6+4\sin^2\theta.                                  \tag{7}
\]

For each integer `m>=1` use the normalized, smooth, periodic vector

\[
\phi_m(\theta)=
\frac{e^{-\iota m\theta}(1+e^{\iota\theta})^{2m}}
{\sqrt{\binom{4m}{2m}}}
=\frac{2^{2m}\cos^{2m}(\theta/2)}{\sqrt{\binom{4m}{2m}}}.
                                                               \tag{8}
\]

Its Fourier support is exactly `-m,...,m`. The coefficients are
`binom(2m,m+n)/sqrt(binom(4m,2m))`; all quadratic readings are rational. This is a normalizable preparation near the reference cut `theta=0`, not a delta state.

**ZCR-T2 (packet moments).** In (8),

\[
q_m:=\langle r\rangle=\frac1{4m+2},\qquad
\langle r^2\rangle=\frac3{(4m+2)(4m+4)},
\]
\[
\langle L^2\rangle=\frac{m^2}{4m-1},\qquad L=-\iota\partial_\theta,
\]
\[
\boxed{\epsilon_m:=\langle H_g-6\rangle
=\eta^2\frac{m^2}{4m-1}
+\frac{2(4m+1)}{(m+1)(2m+1)}.}                            \tag{9}
\]

In particular `Var(r)>0` for every finite `m`.

**Proof.** Vandermonde's finite coefficient identity gives the normalization in (8). It can also be obtained by taking the constant Fourier coefficient of `(1+e^{iota theta})^{2m}(1+e^{-iota theta})^{2m}`. For the moments, put `y=sin²(theta/2)`. Integrals over the circle become integrals with positive weight proportional to `y^(-1/2)(1-y)^(2m-1/2)` on `[0,1]`. Integration by parts gives the recurrence

\[
\frac{\int_0^1 y^{j+1-1/2}(1-y)^{2m-1/2}dy}
{\int_0^1 y^{j-1/2}(1-y)^{2m-1/2}dy}
=\frac{j+1/2}{j+2m+1}.                                   \tag{10}
\]

The first two moments follow at `j=0,1`. Since
`phi_m'=-m tan(theta/2) phi_m`, integration of its squared derivative gives
`<L²>=m²<y/(1-y)>=m²/(4m-1)`. The latter ratio follows by integrating the derivative of `y^(1/2)(1-y)^(2m-1/2)`; endpoint terms vanish for `m>=1`. Finally `sin²theta=4y(1-y)` gives (9). The variance is the difference of the first two rational moments, strictly positive. QED.

The native certificate checks (9) by independent finite Fourier convolution and number-operator multiplication, rather than numerical quadrature or fitted coefficients. The general `m` proof is (10), not extrapolation from those fixtures.

## 5. Localization persists for a controlled operation interval

Let `phi(s)=exp(-iota sH_g/eta)phi(0)` be a normalized finite-energy geometry state and put
`epsilon=<H_g-6>`, `q(s)=<r>_s`. Then

\[
\boxed{\sqrt{q(s)}\le\sqrt{q(0)}+\sqrt\epsilon\,|s|.}       \tag{11}
\]

**Proof.** Energy conservation and positivity of both terms in (7) give `||p phi(s)||²<=epsilon`. On smooth states the commutator identity yields
`q'=<p r'+r'p>`, where `r'=sin(theta)/2` and `(r')²=r(1-r)<=r`. Thus `|q'|<=2 sqrt(epsilon) sqrt(q)`. Apply the differential inequality to `sqrt(q+delta)` and let `delta` decrease to zero. Form-domain approximation extends it to finite-energy data. QED.

This estimate is intentionally referenced to one well. Small energy alone would not distinguish `theta=0` from the other potential minimum `theta=pi`; the initial `q(0)` does. In particular, centering a packet at the other well has the same energy but approaches lapse `1+af(x)`, not `1`.

The bound is valid for all finite `s` but only useful while its right side is small. It does not claim that the wavefunction is stationary or that tunneling is forbidden forever.

## 6. Full spatial AGD dynamics approaches the fixed flat metric

Keep the **actual variable profile** of AGD,

\[
a=3/5,\quad f(x)=\frac{x^2}{1+x^2},\quad
v(\theta,x)=1+af(x)r(\theta),\quad C=Z,\quad R=XZ.
\]

The full self-adjoint joint operator and reference operator are

\[
H=H_g+\iota\eta C(v\partial_x+v_x/2)+\iota\eta gRv,
\]
\[
H_{\rm ref}=H_g+H_{\rm flat},\qquad
H_{\rm flat}=\iota\eta(C\partial_x+gR).                    \tag{12}
\]

They act on the same joint carrier. The reference geometry evolves under `H_g`; it is not frozen or reset. Let `chi` be any normalized direction-valued matter wavefunction in `H¹(R)`, let `d=||partial_x chi||`, and define

\[
K=d+\tfrac12+|g|.
\]

**ZCR-T3 (finite-width flat-metric recovery).** For initial product data
`Psi_0=phi_m tensor chi`, let `Psi(s)=exp(-iota sH/eta)Psi_0` and
`Psi_ref(s)=exp(-iota sH_ref/eta)Psi_0`. For `s>=0`,

\[
\boxed{\|\Psi(s)-\Psi_{\rm ref}(s)\|
\le E_m(s):=aK\left(s\sqrt{q_m}
+\tfrac12s^2\sqrt{\epsilon_m}\right).}                    \tag{13}
\]

For every bounded matter observable `O`, including any bounded detector effect, the reduced full reading and flat-metric reading obey

\[
\boxed{|\langle I\otimes O\rangle_{\Psi(s)}
-\langle O\rangle_{\chi_{\rm flat}(s)}|
\le2\|O\|E_m(s).}                                       \tag{14}
\]

Here the notation orders geometry before matter; the differential operators in (12) have the same content as AGD's opposite tensor-order convention.

**Proof.** The flat evolution commutes with spatial translation, hence
`||partial_x chi_flat(t)||=d`, and preserves the matter norm. Also `||f||_infinity<=1` and `||f'||_infinity<=1`: use `2|x|<=1+x²` in `|f'|=2|x|/(1+x²)²`. On the reference product field,

\[
(H-H_{\rm ref})/\eta
=\iota a r\,[C(f\partial_x+f'/2)+gRf].                   \tag{15}
\]

Both `C` and `R` have operator norm one. Since `0<=r<=1`, the norm of (15) applied to that product field is at most
`aK ||r phi(t)|| <= aK sqrt(q(t))`.

For smooth product data, differentiate
`exp(-iota(s-t)H/eta)Psi_ref(t)` and integrate. Unitarity of the full group bounds its norm difference by the integral of the preceding residual. Equation (11) then gives (13). The product of a periodic `H²` geometry state and an `H¹` matter state belongs to the joint closure domain by approximation with AGD's smooth core; all terms in (12) converge under these approximations. The reference group preserves these separate regularities. This justifies the differentiation and then extends (13) to the stated `H¹` matter data. Finally the difference of two normalized vector quadratic readings is bounded by twice the observable norm times their vector difference; tracing out geometry leaves the reference matter reading unchanged. This proves (14). QED.

There is no bounded-operator assumption about the full perturbation in (15): its spatial derivative is controlled on the evolving reference states. Conservation of `epsilon_m` is used only along the reference geometry flow, not assumed for geometry alone under the interacting full `H`. The full unitary in the Duhamel comparison retains that interaction and its backreaction. There is also no uniform estimate for arbitrary high-spatial-frequency matter; the explicit `d` is essential.

### An explicit joint limit and a rational budget

Choose integer `m=n²` and dimensionless action scale `eta=1/m`, with `n>=1`. The exact moments imply

\[
q_m\le\frac1{4n^2},\qquad
\epsilon_m\le\frac{13}{3n^2}<\frac{441}{100n^2}.
\]

Indeed the kinetic term in (9) becomes `1/(4m-1)<=1/(3m)`, and the potential term is at most `4/m`. Consequently

\[
\boxed{E_m(s)\le\frac{aK}{n}
\left(\frac{s}{2}+\frac{21s^2}{20}\right).}                \tag{16}
\]

This proves convergence to the scalar flat-metric matter law, uniformly on every fixed finite operation interval, for fixed `d` and `g`. The full geometry wavefunctions need not converge in norm to a delta vector. The comparison error tends to zero for this family of generators on the same joint carrier, and the matter readout has a fixed limiting law.

For an explicit budget take
`n=1000`, `m=1000000`, `eta=1/1000000`, `g=0`, `s<=1`, and `d<=1`. Then `K<=3/2` and

\[
\boxed{E_m(s)\le\frac{279}{200000},\qquad
\text{readout error for }\|O\|\le1\ \le\frac{279}{100000}.} \tag{17}
\]

Normalized matter profiles with `d=1` exist; for example `chi(x)=exp(-|x|)e` for a fixed unit direction vector `e` is in `H¹`, has norm one and derivative norm one. The theorem is independent of a spatial grid. The large `m` in (17) is an analytic preparation budget, not a claim that a million-mode matrix was simulated. Its initial Fourier carrier has `2m+1` modes, and the exact subsequent evolution is infinite.

Equation (17) is a worst-case mathematical bound for selected dimensionless inputs. It is not experimental accuracy, a derived value of a physical constant, or an optimal preparation cost.

## 7. Why fixing zero or narrowing alone is insufficient

At fixed `eta>0`, shrinking (8) indefinitely raises its geometry kinetic energy:
`<p²>=eta²m²/(4m-1)`, which diverges with `m`. Thus the controlled limit in (16) changes preparation and dimensionless action scale together. The estimate cannot justify arbitrary-precision classical recovery by narrowing at fixed `eta`. No claim that (16) is a necessary condition for every possible classical limit is made.

A scalar change of units changes the numerical representatives and their dictionary consistently. It does not change the zero/nonzero status of (1), or a dimensionless accumulated leakage budget `ell |s|` when the generator, parameter and action unit are transformed consistently. Nor does it implement `eta=1/m` inside one already calibrated physical model. If `eta` is eventually interpreted as an effective action ratio, the physical regime making that ratio small has to be identified separately.

The flat-metric matter equation remains a quantum spinor equation. A classical particle trajectory, unique measurement outcome/Born rule, Einstein dynamics, and local gravitational constraint closure do not follow from (13). Those are distinct remaining targets. The explicit variable-profile comparison is the progress here: it turns the proposed fixed reference into a finite-time, finite-precision operational statement with stated preparation costs.

## 8. Coverage and next target

The native certificate covers exact projection defects, origin-shift invariance, counterexamples to insufficient criteria, finite Fourier packet moments, translation to the other well, circle word identities, and rational inequalities in the explicit budget. The analytic proofs of (11)–(16) are written arguments over the already admitted AGD completion; they are not proof-assistant checked by the finite certificate.

The next dynamical target remains the local geometry/coframe action and its constraint algebra. ZCR-1 now supplies a concrete recovery estimate that any proposed extension should preserve or improve. Recovery of a general curved evolving classical field, physical calibration, and a complete gravitational classical limit remain open.
