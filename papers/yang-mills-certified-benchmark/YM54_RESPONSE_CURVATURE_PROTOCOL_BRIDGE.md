# YM-54: response curvature to a counted native heat protocol

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

This chapter connects the existing native thermodynamic response algebra
to YM51's counted-turn heat construction and YM53's interacting chain.
It does not derive another operator algebra. For a declared two-direction
response protocol, the exact curvature marker f and trace-free response
budget tau give

\[
\boxed{\beta_{\rm protocol}\geq\frac{4f^2}{\tau}.}
\]

If every fixed site response has |f|>=f_0>0 and tau<=T, the resulting
rank-two protocol tensors satisfy the uniform second-eigenvalue floor
beta_0=4 f_0^2/T. Therefore YM53 applies in the sufficient window

\[
\boxed{|\theta|<\frac{f_0^2}{600T}.}
\]

This is an explicit response-to-protocol construction and a conditional
uniform interacting bound. The source response, two directions, equal
record counts, freezing at a reference state and heat clock are specified
choices. Neither the algebra nor curvature alone selects them physically.
The response Gram matrix is not identified with the protocol tensor by
notation: their different dimensions and factor of four are proved below.

## Inputs and source order

Use NT-1--NT-4 in [Native thermodynamic curvature](../native-thermodynamic-curvature/THEOREM.md)
and the EMK-C1 derivation contract. The scalar field, dagger, faithful
two-mode sector and positive square roots are already available:

\[
K^2=L^2=1,\ R^2=-1,\ L=RK,\ [K,L]=-2R,
\quad K^\dagger=K,\ L^\dagger=L,\ R^\dagger=-R.
\]

The central scalar iota_Sigma squares to -1 and commutes with these
words; it is distinct from the internal element R. Write Tr=2 sc,
where sc is NT-1's cyclic native scalar part. In the faithful two-mode
representation this is the ordinary matrix trace.

The response element H is NT-1's positive invertible element
(a+c)1/2+(a-c)K/2+bL, with a>0 and ac-b^2>0. Its existing native
factor B obeys B^dagger B=H. The derivations commute with dagger and
obey the declared frame bracket. NT-3 selects A_i=X_i/2, where
X_i=H^{-1}delta_i H, and proves

\[
F_{12}=-\tfrac14[X_1,X_2].                                  \tag{1}
\]

The frame-subtracted curvature is understood when directions do not
commute. A_i=X_i is the different, flat pure-gauge connection. The
selection A=X/2 is made in the response trivialization; it is not
reimposed after every variable gauge change.

Only the value of H and its two first response derivatives at a chosen
reference state enter the construction at one chain site. These data
are frozen when that site's heat acts on the separate compact quaternion
variable. A state-dependent-coefficient diffusion is not asserted.

**Lineage.** Positive-matrix trace metrics and commutator geometry are
established mathematics; see Pennec, Fillard and Ayache, *A Riemannian
Framework for Tensor Computing* (2006), DOI 10.1007/s11263-005-3222-z,
listed on the [author's research page](https://www-sop.inria.fr/members/Xavier.Pennec/ManifoldValuedTensorProcessing.html).
The contribution here is the normalized native word/record bridge into
this programme. All displayed identities are derived from the pinned
native sources; no theorem about that external metric is a proof premise.

## YM54-T1: separate scale and positive shape response

Define

\[
\sigma_i=\tfrac12\operatorname{Tr}X_i,\qquad
\widetilde X_i=X_i-\sigma_i1,\qquad
Y_i=B\widetilde X_iB^{-1}.
\]

In the differentiable determinant chart, Tr X_i=delta_i log det H.
It is not generally zero. For example H=e^x 1 gives X_x=1 and
Tr X_x=2. The native product rule and positivity instead give

\[
X_i^\dagger H=HX_i,\qquad Y_i^\dagger=Y_i,\quad
\operatorname{Tr}Y_i=0.
\]

**Proof.** Since delta_i H is self-dagger,
X_i^dagger H=delta_i H=HX_i. Subtracting the real scalar sigma_i
preserves this identity. Insert H=B^dagger B to prove the conjugated
adjoint identity; cyclicity proves its zero trace. The native basis
therefore gives unique real p_i,q_i with

\[
Y_i=p_iK+q_iL.                                               \tag{2}
\]

Set

\[
G_{ij}=\operatorname{Tr}(\widetilde X_i\widetilde X_j)
      =\operatorname{Tr}(Y_iY_j)=2(p_ip_j+q_iq_j),\qquad
\tau=\operatorname{Tr}G.
\]

Thus G is a positive semidefinite **two-by-two shape Gram matrix**.
It is positive definite exactly when the two shape vectors (p_i,q_i)
are independent. The full response Gram splits as
Tr(X_i X_j)=G_ij+2 sigma_i sigma_j. A nonzero trace/scale response can
therefore increase its size without creating commutator curvature.
One must not replace Tr(X_i X_j) by Tr(X_i^dagger X_j): for
H=diag(2,1), delta H=L, the self-pairings are respectively 1 and 5/4.

This is native positive factorization, not the false assertion that a
real operator is automatically self-adjoint. If H is indefinite the
conclusion fails: H=K, delta H=L gives Tr((H^{-1}delta H)^2)=-2.

## YM54-T2: curvature measures shape area, with a quantitative floor

Let d=p_1q_2-q_1p_2 and Z=BF_{12}B^{-1}. From (1)--(2),

\[
[Y_1,Y_2]=-2dR,\quad Z=\tfrac d2R,\quad
\mathcal E_H(F):=\operatorname{Tr}(Z^\dagger Z)=d^2/2.
\]

The two-dimensional Gram determinant is consequently

\[
\boxed{\det G=4d^2=8\mathcal E_H(F).}                        \tag{3}
\]

For NT's positive oriented factor B, B J_H B^{-1}=R, where
J_H=RH/sqrt(det H). The native oriented marker is
f=-sc(J_H F_{12}), so Z=fR, d=2f and

\[
\boxed{\det G=16f^2.}                                      \tag{4}
\]

Indeed B R B^dagger=(det B)R and det B=sqrt(det H), proving the
quarter-turn conjugation. These are identities in this two-mode response
sector, not statements about every native curvature target. Reversing the
ordered response directions changes f's sign, while its squared magnitude
is unchanged. A reflected factor changes the coordinate d's sign; using
the transformed J_H keeps the intrinsic marker f unchanged.

If f!=0, let 0<g_-<=g_+ be G's eigenvalues in its proved finite
representation. Their sum is tau and product is 16f^2, hence

\[
g_-\geq\frac{16f^2}{\tau}>0.                               \tag{5}
\]

The elementary inequality uses g_+<=tau. Independence at each instance
does not by itself give a uniform positive floor across a family.
Uniform |f|>=f_0>0 and tau<=T are sufficient data for such a floor.
The class is nonempty only if T>=8f_0, since tau>=2 sqrt(det G)=8|f|.
The same argument accepts certified interval bounds on |f| and tau;
a sign-uncertain interval containing zero supplies no positive bound.

## YM54-T3: compactify the shape words and count four actual turns

Inside the existing complexified EMK algebra put

\[
e_1=\iota_\Sigma K,\quad e_2=\iota_\Sigma L,\quad e_3=R.
\]

Direct multiplication gives e_a^2=-1,
e_1e_2=e_3, e_2e_3=e_1, e_3e_1=e_2 and reversed products with the
opposite sign; e_a^dagger=-e_a. This realizes YM50's native quaternion
sector. In particular iota_Sigma Y_i=p_i e_1+q_i e_2 is a valid
imaginary turn word. Y_i itself is not substituted for a skew turn.

Define v_i=(p_i,q_i,0). Choose four equally counted labels (i,s),
i=1,2, s=+1,-1, and the actual readout

\[
\mathcal S_\epsilon\phi(x)=\tfrac14\sum_{i=1}^2\sum_{s=\pm1}
\phi\!\left(\operatorname{Exp}_\Sigma(-s\epsilon\,v_i\cdot e/2)x\right).
                                                                    \tag{6}
\]

The half-angle is YM51's fixed normalization. Thus (6) is precisely
YM51's protocol with two weights 1/2, not an assumed differential heat
law. Its signed first moment is zero. Expanding the paired turns gives
I-epsilon^2 L_C/2 at second order, with

\[
\boxed{C=\tfrac12(v_1v_1^T+v_2v_2^T),\qquad
L_C=-\tfrac12(D_{v_1}^2+D_{v_2}^2).}                        \tag{7}
\]

Constants, positivity, Phi_Q and dagger symmetry are preserved by each
readout. Let V have columns v_1,v_2. Then G=2V^TV and C=VV^T/2.
The nonzero spectra of VV^T and V^TV agree, as follows by applying V
or V^T to a nonzero eigenvector. Therefore, when f!=0,

\[
\boxed{\operatorname{spec}C=\{0,g_-/4,g_+/4\},\quad
\operatorname{tr}C=\tau/4,\quad
\beta_{\rm protocol}=g_-/4\geq4f^2/\tau.}                  \tag{8}
\]

The exact smaller positive eigenvalue is
(tau-sqrt(tau^2-64f^2))/8. Its discriminant is nonnegative by the
two-by-two Gram identity. Rank two, rather than isotropy, is the result.
Its missing direct direction is generated by the nonzero turn bracket.
No claim that connection curvature equals a Yang--Mills gauge field
strength is required or made by this bridge.

YM52 also gives the exact **free** centered rate of this protocol:

\[
\gamma_{\rm free}=\min\{\tau/16,g_-/4\},\qquad
\gamma_{\rm even}=g_-/4.                                   \tag{9}
\]

It is not Tr G/2 or the interacting chain rate. A response matrix and
the evolution energy generator are different objects.

## YM54-T4: refinement and the choices hidden by a tensor alone

For (6), m_2=tau/4 and m_4=(G_11^2+G_22^2)/8. YM51-T2 proves on
degree d

\[
\|(\mathcal S_{\sqrt{2t/n}})^n-E_t^C\|_d
\leq\frac{d^4t^2}{96n}
 \left(\frac{G_{11}^2+G_{22}^2}{8}+\frac{3\tau^2}{16}\right)
\leq\frac{5d^4t^2T^2}{1536n},\quad\tau\leq T.              \tag{10}
\]

The last inequality uses G_11,G_22>=0 and their sum tau. YM51's
coefficient approximation extends this limit to the same uniform and
recognition completions; (10) is finite-degree, not a degree-uniform
operator-norm error at zero time.

Changing the positive factor to B'=UB, U^dagger U=1, conjugates the
Y_i. For real two-mode factors, U is orthogonal and its action on
(iota K,iota L,R) is a rotation O of the compact three-frame. Thus
C'=OCO^T: its eigenvalue floor is independent of the factor. A constant
orthogonal recombination of the two response directions also leaves
VV^T unchanged. A general rescaling or nonorthogonal recombination need
not do so and changes the declared sampling protocol.

Equal C does not erase the labels or determine finite events. The
direction pairs (K,2L) and ((3K+8L)/5,(-4K+6L)/5) have the same
VV^T but different m_4. Their readouts already differ at fourth order
on the scalar quaternion coordinate. The operator limit coincides;
the finite records and higher-moment error do not.

More directly, replace (6) by idle weight 1-rho and active weight rho,
0<rho<=1. The source H, response curvature f and G stay fixed, but
C becomes rho C and every free rate scales by rho. This supplies an
explicit obstruction to deriving a unique clock or gap from response
curvature alone. Equal active counts in (6) are the declared process.

## YM54-T5: pass the derived floor to the interacting chain

At each site of an open finite chain, freeze a response jet satisfying
the same source contract, |f_i|>=f_0>0 and tau_i<=T. Apply (6) at each
site, permitting different positive factors and fixed compact orientations.
Then

\[
\beta_0=4f_0^2/T>0,\qquad b_i(C_i)\geq\beta_0,
\qquad\max_i\operatorname{tr}C_i\leq T/4.                 \tag{11}
\]

These are precisely YM53's rank-two/coarse-mixing and trace hypotheses.
For its declared bounded overlap interaction B_int=theta sum v_bond,
choose integer J>=1 and R>1 with

\[
R\left(\frac{5}{2J}+240\frac{|\theta|}{\beta_0}J\right)<1,
\qquad\gamma=\frac{\beta_0\log R}{10J}.                   \tag{12}
\]

For every finite width and 0<h<=1/beta_0, the actual square-sourced
transfer satisfies YM53's full vacuum-complement bound
||Q_h(S_h/lambda_h)^nQ_h||<=exp(-gamma*n*h). It survives the
operator-norm time limit at each fixed width. With b_*=|theta|(m-1),
the norm refinement budget is at most

\[
b_*t e^{b_*t}\sqrt{3T\,\mathrm{mesh}/4}.                    \tag{13}
\]

The limiting interacting ground-source energy also inherits YM53-T6's
weighted derivative-energy identity and lower bound gamma. This does
not replace that ground weight by the response metric H or free Phi_Q.

J=5 and the same strict-gate choice of R as YM53 prove the displayed
window |theta|<beta_0/2400=f_0^2/(600T). Thus this chapter supplies
one concrete native source of the tensor bounds used there. The
anisotropic infinite-volume/joint-cutoff state and actual row closure
remain separate obligations, unchanged by this substitution.

## YM54-T6: the existing curved thermo fixture gives an exact protocol

Reuse NT-4's dimensionless potential
u(s,v)=10s-10v+s^2+sv+(5/2)v^2+(1/2)s^2v at the origin. Its jet is

\[
H=\begin{pmatrix}2&1\\1&5\end{pmatrix},\quad
H_s=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
H_v=\begin{pmatrix}1&0\\0&0\end{pmatrix}.
\]

For exact rational conjugation use Bhat=[[2,1],[0,3]], which obeys
Bhat^T Bhat=2H. The common factor sqrt(2) cancels from conjugation,
so it is equivalent to the normalized B of T1. Direct native coefficient
evaluation gives

\[
v_1=(1/9,1/3,0),\quad v_2=(2/9,-1/6,0),\quad
G=\begin{pmatrix}20/81&-5/81\\-5/81&25/162\end{pmatrix},
\]
\[
f=-5/108,\quad\tau=65/162,\quad
C=\operatorname{diag}(5/162,5/72,0).                        \tag{14}
\]

The bound (8) gives beta_0=5/234, while this fixture's exact second
eigenvalue is 5/162. The sufficient curvature-budget window is
|theta|<1/112320. For theta=1/262144, J=6, R=5/4, (12) passes and

\[
\gamma=\log(5/4)/2808>0.000079467076.
\]

Copying these frozen jets along a chain gives a concrete family to
which T5 applies; it does not identify a physical material with a
four-dimensional gauge vacuum. Formula (9) gives free full rate
65/2592 and even rate 5/162, distinct from both Tr G/2=65/324 and
the conservative interacting bound. Deleting the cubic potential term
leaves H unchanged at the origin but sets both response derivatives,
f and this constructed protocol tensor to zero. Point eigenvalues of
H alone therefore cannot supply this bridge.

## YM54-T7: scale, degeneracy and operator-source boundaries

1. Multiplying H by a positive scalar function adds only
   (delta_i log s)1 to X_i. The shape words, G, F and constructed C
   are unchanged. This protocol intentionally reads shape; it is not a
   universal dynamics of all scalar response changes.
2. Multiplying both derivatives by k scales G,C,beta by k^2 and
   f by k^2. The squared curvature energy scales by k^4. Simultaneously
   changing t,h by 1/k^2 and theta by k^2 preserves the transfer and
   comparison margins. A physical time unit has not been selected.
3. For shape words K and epsilon L, epsilon>0, every Gram is positive,
   but C=diag(1/2,epsilon^2/2,0) and its smaller positive eigenvalue
   tends to zero. Tr G/2=1+epsilon^2 stays above one. Positivity and
   trace alone cannot certify a uniform gap.
4. Even a fixed nonzero curvature is insufficient without a budget:
   shape words tK and L/t have f=1/2 for every t>0, but the smaller
   positive protocol eigenvalue tends to zero as t grows. Tau grows.
5. Independent full responses 1 and K have positive full Gram but
   commute; their shape Gram has rank one. Zero f refuses this
   certificate. Such a refusal is not a claim that every other possible
   process on the same algebra must be gapless.

The existing operator sources have specific roles, rather than being
interchangeable names for a Hamiltonian. RKF T51's Aghora law carries
odd/even exchange; T56 records why a nonzero odd-signed generator must
not simply be declared positive in the derived square-source sense.
Here the imaginary lift in T3 makes lawful skew turns, and the counted
sum of derivative squares makes their positive heat energy. This is an
explicit construction, not a new assertion about every Aghora operator.
RKF T75's RV 7.59.12 selective release retains information under its
declared tail/Gram conditions. Its infinite-tail premises are not
silently transferred to this new source family. Those operator papers
are audited lineage/context; T1--T6 rely on NT and YM50--53.

## Evidence and claim boundary

The exact certificate checks the two-mode word embedding into the
existing quaternion engine, positive-factor response jets, scale/shape
separation, curvature/Gram minors, actual four-label moments, independent
coefficient generator actions, factor and direction covariance, finite
record distinctions, outward refinement bounds and the NT fixture's
interacting gate. Degenerate/indefinite inputs, false trace/gap claims,
missing compactification and unsupported uniformity claims have explicit
negative controls. No canonical engine or prior certificate is changed.

General arguments are written proofs supported by exact finite controls,
not formal proof-assistant verification or external expert certification.
The model selects a two-direction equal-count response protocol at frozen
reference jets. Physical selection of that protocol, response law,
observer directions, clock, coupling and state remains open, as do a
state-dependent diffusion extension, anisotropic joint/volume limit,
actual row closure, NCG quantum measure, four-dimensional QFT, AF/Clay
and quantum gravity. An algebraic completion alone proves none of these.
