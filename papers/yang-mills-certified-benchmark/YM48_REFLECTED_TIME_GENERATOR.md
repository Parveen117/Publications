# YM-48: reflected time action, its generator, and a nonzero local excitation

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

**Result on the existing heat/positive-functional adapter:** YM-47's joint
history limit constructs a completed positive-form carrier on which literal
forward time translation is a strongly continuous positive self-adjoint
contraction semigroup. Its nonnegative self-adjoint generator has an explicit
resolvent-range domain and graph core. The same gamma from YM-45/46 gives
an all-source gap above the unique vacuum on this carrier.

Local finite-coefficient row observables also have a controlled time-zero
embedding. A normalized fundamental character at one site supplies a
nonzero vacuum-orthogonal vector, so this gap is not a statement about a
one-dimensional vacuum space alone.

The carrier is the **reflected history completion**. Its identification
with a closed instantaneous-row Markov theory is not assumed. A complete
observable representation, native Phi_Sigma/NCG dictionary, physical clock,
relativistic spacetime and four-dimensional Yang--Mills remain open.

## Inputs, reused mathematics and carrier distinction

Keep abs(theta)<1/1680, kappa(a)=theta a, the countable spatial chain,
the full SU(2) coefficient heat family and the admitted finite reference
functional of YM-44/45/46/47. Fix one lawful YM46-T1 cell (or its full-window
construction) and its rate gamma>0. Write Omega for YM-47's positive
functional on finite sums of bounded local history products.

At a finite interval I, use the common positive unit vacuum h_I and
\[
V_I(t)=U_I(t)/\|U_I(t)\|,\qquad
\|Q_I V_I(t)Q_I\|\leq e^{-\gamma t}.
\]
These are already constructed. All readouts below are the vacuum histories
whose remote time endpoints have already been removed.

RKF N10/N11 supply the positive-pairing, null-space and bounded-extension
discipline; N08 supplies the existing resolvent algebra; T28 requires an
honest common carrier and explicit tails. Here that common algebraic carrier
is the space of finite future histories. We pass inequalities of their
positive forms to the limit, rather than claiming operator-norm convergence
of finite Hilbert spaces with different dimensions. WC2's bounded series
calculus is not applied to an unbounded generator. The canonical engine
is unchanged.

Reflection-positive reconstruction has the established
[Osterwalder--Schrader 1973](https://doi.org/10.1007/BF01645738) /
[1975](https://doi.org/10.1007/BF01608978) lineage. The application proof is
given here; no general reconstruction theorem or full relativistic
OS/Wightman axiom package is imported as a premise. No priority claim is made.

## YM48-T1: finite folding identifies the positive history form

Let A_+ be the linear span of finite products of bounded local row
observables at **strictly positive** heat times, together with the constant
one. Let r reflect t to -t, without conjugating coefficients. Define
\[
(F,G)_+=\Omega(\overline{rF}\,G).
\tag{1}
\]
YM-47 gives positivity. For a product
F=F_1(t_1)...F_k(t_k), 0<t_1<...<t_k, define the finite folded source
\[
\Psi_I(F)=V_I(t_1)M_{F_1}V_I(t_2-t_1)M_{F_2}\cdots
V_I(t_k-t_{k-1})M_{F_k}h_I,\qquad \Psi_I(1)=h_I.
\tag{2}
\]
Extend linearly. It is also h_I times the finite conditional expectation
of the future history given the row at zero. This follows by inserting
the finite Doob factors; successive h_I and normalization factors cancel.
It proves independence of a chosen product expansion and
||Psi_I(F)||<=||F||_infinity.

Finite integration, time reversal and the common vacuum give
\[
\Gamma_I(\overline{rF}G)=\langle\Psi_I(F),\Psi_I(G)\rangle,\quad
\Psi_I(\tau_t F)=V_I(t)\Psi_I(F),\quad t\geq0,
\tag{3}
\]
where tau_t moves every time label forward by t. If some times coincide
within a product, combine their functions before applying (2).
These are finite identities, not assertions about an infinite path measure.

YM-47's bounds apply term by term to (3), so its left side converges to
(1). The finite folded sources themselves are not declared to converge in
one finite-volume norm.

## YM48-T2: quotient, completion and actual composition

Positive-form Cauchy--Schwarz follows by expanding
(F+zG,F+zG)_+ for cut-complex z, including the zero-diagonal case.
Its null space N is the radical of the pairing. Thus A_+/N has the
well-defined positive pairing (1); complete its Cauchy sequences in that
norm to obtain H_ref. Put e=[1], so ||e||=1.

Finite contraction in (3), followed by the YM-47 limit, gives
\[
\|\,[\tau_tF]\,\|_{\rm ref}\leq\|\,[F]\,\|_{\rm ref}.
\tag{4}
\]
In particular null histories stay null. Define T(t)[F]=[tau_t F].
It extends uniquely to a contraction on H_ref. The literal identity
tau_s tau_t=tau_(s+t) on the common history core, followed by (4), proves
\[
T(s)T(t)=T(s+t),\quad T(0)=I,\quad T(t)e=e.
\tag{5}
\]
Finite self-adjointness in (3) passes to
(F,T(t)G)_+=(T(t)F,G)_+. Hence the bounded T(t) is self-adjoint.
Also T(t)=T(t/2)^dagger T(t/2), so it is positive in this Hilbert form.

The proof of composition is this norm-controlled translation construction.
Composition is not inferred merely from convergence of selected correlations.
The null space N need not be an ideal for multiplication by an arbitrary
future observable. No such multiplication representation has been smuggled
into the quotient.

## YM48-T3: strong continuity and the all-source gap

Represent F as a finite sum of products, let M_F be the sum of their
coefficient magnitudes times the products of the sup norms, and choose
epsilon>0 below all their nonconstant observation times. Constant terms
are fixed by translation. YM47-T3 gives
\[
\|(V_I(t)-I)V_I(\epsilon)\|\leq\min(1,t/\epsilon).
\]
Every term in (2) starts with this smoothing factor. Therefore, uniformly
in I, and then in its limiting history norm,
\[
\boxed{\|(T(t)-I)[F]\|_{\rm ref}
\leq M_F\min(1,t/\epsilon).}\tag{6}
\]
This proves strong continuity at zero on a dense core. Contraction and a
two-approximation estimate extend it to every x in H_ref, and (5) gives
strong continuity at every t>=0. No full operator-norm continuity at zero
is asserted.

Let Q=I-|e><e|. In the finite form,
<h_I,Psi_I(F)>=Gamma_I(F). Apply the finite all-source gap to
Psi_I(F)-Gamma_I(F)h_I, and pass both squared norms through (3). Obtain
\[
\boxed{\|T(t)Qx\|\leq e^{-\gamma t}\|Qx\|,
\qquad x\in H_{\rm ref}.}\tag{7}
\]
Density extends the bound to every completed source. Invariance of e and
self-adjointness give QT(t)=T(t)Q. This is an operator inequality on the
newly constructed carrier, not just a two-point bound.

## YM48-T4: a constructive resolvent with cut and quadrature tails

For lambda>0 define
\[
R_\lambda x=\int_0^\infty e^{-\lambda t}T(t)x\,dt.
\tag{8}
\]
The integral means strong Riemann cut sums. On compact intervals the
integrand is norm-continuous on each source; the omitted tail has bound
\[
\left\|\int_L^\infty e^{-\lambda t}T(t)x\,dt\right\|
\leq e^{-\lambda L}\|x\|/\lambda.
\tag{9}
\]
Completeness supplies the limit. Thus ||R_lambda||<=1/lambda. Finite
dagger identities and positive forms pass through these sums, so R_lambda
is positive and self-adjoint.

For a core history with the data of (6), left endpoint quadrature on
[0,L] with mesh at most delta has total error at most
\[
\boxed{M_F\left\{\tfrac12 L\delta(\lambda+1/\epsilon)
+e^{-\lambda L}/\lambda\right\}.}\tag{10}
\]
Indeed the integrand has Lipschitz constant at most
M_F(lambda+1/epsilon); integrate the linear error on each cell.
For each finite sum, YM-47 separately controls its history matrix elements.
This specifies both the time-integration error and the underlying
volume/refinement error, without claiming a uniform finite-volume lift.

Absolute norm tails permit regrouping the double integral for
R_lambda R_mu by u=s+t. The elementary scalar inner integral yields
\[
R_\lambda-R_\mu=(\mu-\lambda)R_\lambda R_\mu.
\tag{11}
\]
This realizes the existing N08 resolvent identity on the present source
completion; it does not extend the old mass algebra by notation.

R_lambda is injective: if R_lambda x=0, then
0=<x,R_lambda x>=integral exp(-lambda t)||T(t/2)x||^2 dt.
Strong continuity at zero makes this impossible unless x=0.
Moreover lambda R_lambda x->x as lambda grows. On the core (6) gives
the explicit error M_F/(epsilon lambda), and boundedness extends convergence
to the completion.

The range of each R_lambda is dense. One direct proof is that its
orthogonal complement is ker(R_lambda^dagger)=ker(R_lambda)={0}.
Equation (11) also shows that these ranges agree for all lambda>0.
On QH_ref, (7) improves (9) to
exp(-(lambda+gamma)L)||Qx||/(lambda+gamma) and gives
||R_lambda Q||<=1/(lambda+gamma).

## YM48-T5: generator, domain, self-adjointness and graph core

Define, for any lambda>0,
\[
\mathcal D(H)=\operatorname{Ran}R_\lambda,\qquad
H R_\lambda x=x-\lambda R_\lambda x.
\tag{12}
\]
Injectivity makes this unambiguous. The common-range identity (11), or the
derivative calculation below, makes it independent of lambda.

Translation of (8) gives the exact formula
\[
T(h)R_\lambda x
=e^{\lambda h}\left[R_\lambda x-
\int_0^h e^{-\lambda s}T(s)x\,ds\right].
\]
Strong continuity therefore proves
\[
\lim_{h\downarrow0}\frac{T(h)u-u}{h}=-Hu,\qquad u\in\mathcal D(H).
\tag{13}
\]
Conversely, if this limit exists and equals -v, the semigroup property
gives d(T(s)u)/ds=-T(s)v. Integrate its damped derivative with the tail
(9): u=R_lambda(v+lambda u). Thus (12) is precisely the strong-generator
domain, rather than only a convenient subdomain.

The operator is densely defined and closed. For the latter, if
u_n->u and Hu_n->v, use u_n=R_lambda(Hu_n+lambda u_n) and boundedness
of R_lambda to obtain u=R_lambda(v+lambda u) and Hu=v.
Symmetry follows either from (13) and self-adjoint T, or directly from
self-adjoint R_lambda. Positivity follows from
\[
\langle u,Hu\rangle
=\lim_{t\downarrow0}\{\|u\|^2-\langle u,T(t)u\rangle\}/t\geq0.
\]
It is **self-adjoint**, not merely symmetric: H+lambda is onto by (12).
If v is in the adjoint domain, set w=R_lambda(H^*v+lambda v).
Then <(H+lambda)z,v-w>=0 for every z in D(H). Surjectivity gives v=w
in D(H), and H^*v=Hv.

For any fixed lambda, R_lambda(A_+/N) is a graph core: approximate x
by core histories x_n, then R_lambda x_n->R_lambda x and
H R_lambda x_n=x_n-lambda R_lambda x_n->x-lambda R_lambda x.
The domain is invariant under T(t), because T(t) commutes with R_lambda.

Finally, (7) and (13) give the quadratic gap
\[
\boxed{He=0,\qquad
\langle u,Hu\rangle\geq\gamma\|Qu\|^2,\quad u\in\mathcal D(H).}\tag{14}
\]
Thus ker H is exactly the span of e. We have constructed the normalized
heat-time generator; we have not identified it with an unbounded formal
sum of infinite-volume local terms or with a physical Hamiltonian in
material units.

## YM48-T6: local coefficient energy controls the time-zero boundary

This additional calculation uses the **existing SU(2) coefficient adapter**,
not a new primitive metric. On its finite coefficient algebra, let D_(i,alpha)
be the left-action derivations generated at site i by
-iota sigma_alpha/2. Tensor-product differentiation gives Leibniz;
finite reference invariance gives integration by parts. The representation
Casimir identity, in exactly the heat convention already used by YM-44,
is
\[
L_I=-\sum_{i,\alpha}D_{i,\alpha}^2,\qquad
\langle f,L_I f\rangle=\sum_{i,\alpha}\|D_{i,\alpha}f\|^2.
\tag{15}
\]
This follows on every matrix coefficient by the spin-j identity
sum J_alpha^2=j(j+1)I, then by finite products and linearity. It fixes
the factor 1/4 in the fundamental-coordinate computation below.

For completeness, L_I has domain sum_C C^2|f_C|^2<infinity in the
orthogonal coefficient completion. Coordinate truncation proves
closedness and self-adjointness directly. Since B_I is bounded, a
Neumann inverse for L_I+lambda followed by the bounded perturbation
constructs the inverse of L_I-B_I+E_I+lambda for sufficiently large
lambda. Symmetry and surjectivity, as in T5, show that
\[
H_I=L_I-B_I+E_I,\quad\mathcal D(H_I)=\mathcal D(L_I),
\qquad \|U_I(t)\|=e^{E_I t}.
\]
The scalar exponential follows from the common vacuum, the semigroup
identity and continuity. YM-44's Volterra identity gives the generator
on D(L_I), so the normalized semigroup is precisely V_I.
The same identity applied to h_I gives L_I h_I=(B_I-E_I)h_I.
In particular h_I belongs to D(L_I), and H_I h_I=0.

Multiplication by a finite local coefficient polynomial F preserves
D(L_I). To check rather than assume the domain, approximate h in the
L_I graph norm by finite coefficients, use (15) to control each D h, and use
\[
L_I(Fh)=F L_Ih+(L_IF)h-2\sum_{i,\alpha}(D_{i,\alpha}F)(D_{i,\alpha}h).
\]
All displayed multipliers are bounded. Closure then proves the assertion.
The same argument applies to |F|^2.

Expand the derivative squares of F h_I in (15) and subtract the vacuum
equation tested against |F|^2 h_I. The interaction and vacuum scalar
cancel, leaving the exact identity
\[
\boxed{\langle Fh_I,H_I Fh_I\rangle
=\pi_I\left(\sum_{i,\alpha}|D_{i,\alpha}F|^2\right)
\leq C_F,\quad
C_F=\left\|\sum_{i,\alpha}|D_{i,\alpha}F|^2\right\|_\infty.}\tag{16}
\]
Here pi_I=Phi_I(h_I^2\,dot), and h_I is real and positive. This
normalization-dependent identity is not obtained by discarding B_I.
Its constant depends only on the local polynomial F, not on I.

For u in D(H_I), positive contraction and finite geometric summation give
\[
0\leq\langle u,(I-V_I(t))u\rangle\leq t\langle u,H_Iu\rangle.
\]
For the upper bound put W=V_I(t/n):
sum_(j<n)<u,W^j(I-W)u><=n<u,(I-W)u>, because successive terms
differ by the positive form of W^j(I-W)^2. Let n grow and use the
generator derivative. Consequently
\[
\boxed{\|(V_I(t)-I)Fh_I\|^2\leq2t C_F.}\tag{17}
\]
This supplies the missing uniform boundary continuity for this local
coefficient core, without asserting norm continuity of the full heat family.

## YM48-T7: an isometric row embedding and its actual readouts

For a local coefficient polynomial F, consider [F(epsilon)] in H_ref.
Equations (3), (17), and the finite-volume limit imply
\[
\|[F(\epsilon)]-[F(\delta)]\|^2
\leq2|\epsilon-\delta|C_F.
\]
Thus
\[
JF=\lim_{\epsilon\downarrow0}[F(\epsilon)]
\]
exists. The finite estimate against Fh_I, followed by the volume limit,
gives
\[
\boxed{\|JF\|^2=\pi(|F|^2),\qquad
\langle JF,JG\rangle=\pi(\overline F G),}\tag{18}
\]
where pi is the row restriction of Omega. The construction is linear
and factors through the row null space. It extends isometrically to
H_row^coeff, the positive-form completion of local coefficient polynomials.

The translation bound away from zero and (17) also give
\[
T(t)JF=[F(t)]\quad(t>0),\qquad
\|(T(t)-I)JF\|\leq\sqrt{2tC_F},
\]
\[
\boxed{\langle JF,T(t)JG\rangle
=\Omega(\overline{F(0)}G(t)),\quad t\geq0.}\tag{19}
\]
For t>0, pass the finite inner products using (17) to control the
time-zero endpoints uniformly in I; t=0 follows from (18).
Thus the new action reproduces the actual local row correlations.

This proves an **embedding**, not that J is onto H_ref or that its
closed range is invariant under T(t). Those are substantive
instantaneous-row closure/Markov obligations. It also does not put every
bounded measurable row observable or every JF in the operator domain
D(H); the stated core, continuity and resolvent-domain claims are distinct.

## YM48-T8: a nonzero excitation in the interacting carrier

Use F(x)=chi_(1/2)(x_i)/2=x_0 at one site, with
x=x_0 I+iota sum_(alpha=1)^3 x_alpha sigma_alpha and sum x_alpha^2=1,
including alpha=0. Finite V_I and B_I are invariant under a common left
SU(2) multiplication of all sites: each bridge character is conjugated.
The unique positive normalized vacuum is therefore invariant. Averaging
this symmetry with the already admitted finite reference functional shows
that every one-site marginal of pi_I equals that reference functional.
The argument passes to pi.

In particular pi_I(F)=0, pi_I(F^2)=1/4. The mean vanishes under x->-x;
coordinate symmetry and sum x_alpha^2=1 give the second moment.
Direct left multiplication by exp(-iota t sigma_alpha/2) gives
D_alpha x_0=x_alpha/2, hence
\[
\sum_\alpha|D_\alpha F|^2=(1-F^2)/4,\qquad
\langle Fh_I,H_I Fh_I\rangle=3/16.
\tag{20}
\]
Equation (17)'s preceding quadratic estimate now yields
\[
\|[F(\epsilon)]\|_{\rm ref}^2
=\lim_I\langle Fh_I,V_I(2\epsilon)Fh_I\rangle
\geq\frac14-\frac{3\epsilon}{8}.
\]
At epsilon=1/4,
\[
\boxed{\langle e,[F(1/4)]\rangle=0,\qquad
\|[F(1/4)]\|_{\rm ref}^2\geq5/32>0.}\tag{21}
\]
This is a nonzero vector of the interacting chain's own carrier for
every theta in the stated window, without appending a free spectator.
It is an allowed coefficient source here, not a separately established
gauge-invariant physical observable.
It is not a proof of scattering, of a non-Gaussian continuum, or of the
Clay nontriviality requirements.

## YM48-T9: refusal witnesses and what remains

1. **Reflection null is not a multiplication ideal.** In an independent
   reversible two-state fixture let f=+/-1 and P_s f=q f, 0<q<1.
   The future history F=f(2s)-q f(s) folds to zero, but
   F f(2s) folds to the nonzero constant 1-q^2.
   This rejects an automatic bounded representation of all future
   multiplication operators on H_ref.
2. **A projected time map need not compose.** C P^2 C and (CPC)^2
   differ when the eliminated component returns. T2 avoids this by
   translating full histories and proving contraction on their quotient.
3. **Domain and topology matter.** Unbounded free Casimirs preserve
   ||K_t-I||=1 for t>0. Finite cut generators do not supply a bounded
   all-content generator or justify a norm exponential series in H.
4. **The resolvent needs both errors.** A finite time integral omits
   exp(-lambda L)/lambda even on the vacuum. A coarse quadrature also
   has a nonzero mesh error. Both are present in (10).
5. **The strict window is retained.** No cell with beta>=1, rho<=1 or
   gamma<=0 is promoted through reconstruction. A gap inequality alone
   would allow a vacuum-only space; T8 supplies a separate nonzero vector.
6. **The energy normalization is fixed.** The Casimir convention and
   the vacuum subtraction E_I are required in (16)--(20). The coordinate
   factor 1/4 is checked algebraically. The resolvent sign is H+lambda;
   it is not inferred from a classical-force analogy.

The next identification is whether the instantaneous row completion is
all of H_ref and supports the required bounded local observable action,
with the native measure/NCG intertwiner separately retained. Arbitrary
bounded row time-collision claims, a countably additive path measure,
physical real time, relativistic covariance/locality, the spatial
continuum, AF, Clay and QG are not established here.

## Evidence and reproduction

The written infinite proof is not mechanically formalized. Exact finite
controls check reflected folding against full path enumeration, null
translation and the non-ideal witness, mixed-source gap bounds, two-sided
resolvent identities and generator reconstruction, outward quadrature tails,
the fundamental Casimir/Leibniz identities and weighted polynomial energy
identities, plus the strict parameter and refusal gates. Their test
fixtures do not replace the infinite SU(2) argument.
The path fixture is a four-state product chain with rational transition
weights at ticks of duration log(2). The weighted polynomial fixtures use
positive sources h=1+/-x_0/3 with test potential (Lh)/h to check the
ground-source energy identity; they do not claim to compute the interacting
chain vacuum. Their exact sphere moments follow the coordinate-rotation
recurrence and sum x_alpha^2=1.

All prior certificate bytes remain frozen. The source record pins the
canonical RKF contracts by commit/blob and the inherited local evidence.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym48_reflected_time.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym48_reflected_time.py -v
~~~

Default/--check is read-only; --write explicitly regenerates this chapter's
evidence. CI uses one Python 3.12 job.
