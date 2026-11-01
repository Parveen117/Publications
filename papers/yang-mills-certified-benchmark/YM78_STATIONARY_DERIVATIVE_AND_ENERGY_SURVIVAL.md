# YM-78: stationary derivative fields and local energy survival

Monty Dabas. 5 October 2026. **Python 3.12 only.**

YM77 controlled physical-time increments but did not prove a nonzero
reflected local state in the continuum. This stage obtains an exact,
positive physical Dirichlet energy for an actual local observable,
constructs stationary tempered derivative-field subsequences, and gives
an explicit positive reflected-norm lower bound at each regulator.
The latter bound decays as the spacing goes to zero. A rational
frequency contrast isolates precisely the additional estimate needed
for survival at a fixed physical energy scale.

**This is not a proof of four-dimensional Yang–Mills existence or its
mass gap.** The new field limits can vanish, or have nonzero contact
covariance with zero reflected linear norms. Their identification with
renormalized Yang–Mills fields, full Euclidean covariance, and a positive
gap remain open. The finite certificates check the stated algebra,
geometry and examples; they do not replace the analytic proofs.

## 1. Fixed model and a local source with exact positive energy

Use the unchanged finite periodic SU(2) model of YM76–YM77:
\[
 H_a=L_E-\theta_a\sum_bW_b-E_0,\quad
 \mathsf H_a=\kappa_aH_a,\quad
 \mu_a=h_a^2\,d\Phi,\quad
 \kappa_a\sqrt{\theta_a}=\frac ca,\qquad c=\sqrt{c_Ec_B}>0 .
 \tag{1}
\]
The actual vacuum \(h_a\) is normalized, positive and smooth.
Write \(r_b=\arccos W_b\). The radial identities established in YM76 are
\[
 \Gamma(r_b)=1\ \text{a.e.},\qquad
 \langle Fh_a,\mathsf H_aFh_a\rangle
       =\mu_a\Gamma_{\rm phys}(F),\quad
 \Gamma_{\rm phys}=\kappa_a\Gamma .
 \tag{2}
\]
All observables below are centered when viewed as Hilbert-space sources.
Let \(a\to0\), \(aL\to\infty\), and take even \(L\ge4\).
The coupling remains the original physical calibration, not a parameter
defined from an observed spectral gap.

Fix one orientation, say \(xy\). Select its plaquettes based at
\((2i,2j,k)\), with coordinates interpreted periodically. Their link sets
are pairwise disjoint: distinct selected squares in one plane have no
common edge, and different \(z\)-planes have different links.
Associate to each selected plaquette the spatial block
\[
 B_b=[2ia,(2i+2)a)\times[2ja,(2j+2)a)\times[ka,(k+1)a).
 \tag{3}
\]
The blocks partition a fundamental physical cube, each with volume \(4a^3\).
For a real \(f\in C_c^\infty(\mathbb R^3)\), with support inside the cube, set
\[
 \bar f_b=\frac1{4a^3}\int_{B_b}f(x)\,dx,\quad
 P_af=\sum_b\bar f_b\,1_{B_b},\qquad
 \Psi_a(f)=2\theta_a^{1/4}a^2
               \sum_b\bar f_b(r_b-\mu_a r_b).
 \tag{4}
\]
The factor 2 compensates for the selected density \(1/4\).
Disjoint link support is used only for the gradient identity; no
independence of plaquette random variables is assumed.

### Theorem 78.1 — exact source energy and a quantitative floor

For \(v_a(f)=\Psi_a(f)h_a\), the following are exact:
\[
 \boxed{\Gamma_{\rm phys}(\Psi_a(f))
       =q_a(f):=\langle v_a,\mathsf H_av_a\rangle
       =4ca^3\sum_b\bar f_b^2=c\|P_af\|_2^2.}
 \tag{5}
\]
In particular,
\[
 c\left(\|f\|_2^2-\frac{4a^2}{\pi^2}\|\nabla f\|_2^2\right)_+
       \le q_a(f)\le c\|f\|_2^2,\qquad
 q_a(f)\longrightarrow c\|f\|_2^2 .
 \tag{6}
\]
These statements hold at every positive coupling.

**Proof.** Distinct selected angles have zero carré-du-champ cross
term. Equations (1)–(2) give
\(\Gamma_{\rm phys}\Psi_a=4\kappa_a\sqrt{\theta_a}a^4
\sum_b\bar f_b^2\), proving (5).
Block averaging is an orthogonal \(L^2\) projection. On each rectangular
block, the first nonzero Neumann eigenvalue is \(\pi^2/(4a^2)\).
Thus
\[
 \|f-P_af\|_2^2\le\frac{4a^2}{\pi^2}\|\nabla f\|_2^2,\qquad
 \|P_af\|_2^2=\|f\|_2^2-\|f-P_af\|_2^2 .
\]
This proves (6). \(\square\)

There is now a strictly positive, volume- and coupling-independent
source-energy floor for each nonzero fixed \(f\) at sufficiently small
\(a\). This is a form-energy statement, not a reflected-norm or gap
lower bound.

## 2. Stationary space-time derivative fields

Let \(X_a(t)\), \(t\in\mathbb R\), be the two-sided stationary reversible
diffusion with invariant law \(\mu_a\) and physical generator
\[
 \mathcal G_aF=-\kappa_a L_EF
              +2\kappa_a\Gamma(\log h_a,F).
 \tag{7}
\]
For a real \(J\in C_c^\infty(\mathbb R^4)\), define
\[
 Z_a(J)=-\int_{\mathbb R}
                 \Psi_a(\partial_tJ(t,\cdot))(X_a(t))\,dt.
 \tag{8}
\]
It is the distributional time derivative of the angle observable, or
equivalently of its anchored increment process. The anchor disappears
because the time integral of \(\partial_tJ\) is zero.

For global test functions use only blocks of the fundamental cube and
extend the resulting spatial field by zero outside that cube. In
particular, do not periodically replicate the field on all of
\(\mathbb R^3\): such replication would invalidate a global \(L^2\)
test-function bound. Every fixed compact support lies inside the
growing cube eventually.

### Theorem 78.2 — uniform exponential bounds and stationary limits

At every regulator,
\[
 \boxed{\mathbb E e^{\lambda Z_a(J)}
       \le e^{c\lambda^2\|J\|_{L^2(\mathbb R^4)}^2}.}
 \tag{9}
\]
Consequently,
\[
 \mathbb EZ_a(J)=0,\quad
 \mathbb E|Z_a(J)|^2\le2c\|J\|_2^2,\quad
 \mathbb E|Z_a(J)|^{2p}\le2p!(4c\|J\|_2^2)^p,\quad p\ge1 .
 \tag{10}
\]
The variance bound in (10) is sharper than substituting \(p=1\) in
the general tail-moment estimate.

There are common subsequential laws on
\(\mathcal S'(\mathbb R^4)\), stationary under ordinary space and time
translations, with these same bounds. They are supported locally in
\(H^{-s}_{\rm loc}(\mathbb R^4)\) for every \(s>2\).
They preserve reflection positivity with the **odd** time reflection
\[
 (\Theta_Z Z)(J)=-Z(J^\theta),\qquad
 J^\theta(t,x)=J(-t,x).
 \tag{11}
\]
They need not be nonzero.

**Proof of the exponential estimate.** First replace \(r\) by smooth
radial approximations \(r_\epsilon\) converging uniformly to \(r\), with
\(|r_\epsilon'|\le1\), and center them in the same actual vacuum.
Equation (5) becomes an upper bound with the same constant.
Choose time endpoints outside the support of \(J\), and put
\(F_t=\Psi_a(J_t)\),
\(S=\int\partial_tF_t(X_t)dt\), \(D=\int\mathcal G_aF_t(X_t)dt\).
Forward and reversed Itô formulas yield
\[
 0=S+D+M_+,\qquad 0=-S+D+M_-,
 \qquad Z_a(J)=\tfrac12(M_+-M_-).
 \tag{12}
\]
Each martingale has quadratic variation at most
\(2c\int\|J_t\|_2^2dt\).
The two martingales need not be independent: Cauchy–Schwarz and their
separate exponential martingale estimates give
\[
 \mathbb Ee^{\lambda(M_+-M_-)/2}
 \le(\mathbb Ee^{\lambda M_+})^{1/2}
     (\mathbb Ee^{-\lambda M_-})^{1/2}
 \le e^{c\lambda^2\|J\|_2^2}.
\]
At fixed regulator the smooth approximations in (8) converge uniformly
over diffusion paths. Fatou proves (9) for the actual angles.
Differentiating the moment generating function at zero gives the
variance bound. The tail bound
\(\mathbb P(|Z_a(J)|>u)\le2e^{-u^2/(4c\|J\|_2^2)}\)
and tail integration give (10).

**Tempered tightness and regularity.** The variance bound extends
\(J\mapsto Z_a(J)\) continuously from compact smooth tests to real
\(L^2(\mathbb R^4)\), as a map into \(L^2\) of probability.
Let \(h_{\mathbf n}\) be the four-dimensional Hermite basis and
\(\lambda_{\mathbf n}=1+|\mathbf n|\). For \(r>4\),
\[
 \mathbb E\sum_{\mathbf n}\lambda_{\mathbf n}^{-r}
                     |Z_a(h_{\mathbf n})|^2
       \le2c\sum_{\mathbf n}\lambda_{\mathbf n}^{-r}<\infty.
 \tag{13}
\]
This realizes each field in a negative oscillator Hilbert space.
A bounded ball for exponent \(r\) is compact in the corresponding
space for any exponent \(r'>r\). Markov's inequality therefore gives
uniform tightness in that weaker separable Hilbert space, and hence
a tempered subsequential law. The oscillator exponent convention in
(13) is stated explicitly; it is not an assertion of global
unweighted Sobolev integrability of a stationary field.

For compactly supported \(\chi\), the Fourier coefficient
\(Z_a(\chi e^{i\xi\cdot x})\) has second moment bounded by
\(2c\|\chi\|_2^2\), with the complexification interpreted sesquilinearly.
Integration against \((1+|\xi|^2)^{-s}\) is finite in four dimensions
when \(s>2\). A countable exhaustion and slightly weaker Sobolev
indices give the asserted simultaneous local regularity.

**Translations and reflection.** Time stationarity follows from
stationarity of \(X_a\). Spatial translations by \((2ai,2aj,ak)\)
preserve the selected lattice family and the vacuum. For any fixed
translation choose such lattice translations converging to it.
The \(L^2\) covariance bound controls the test-function error, and
the boundary is absent for fixed compact supports in large boxes.
Spatial translation invariance follows in the limit.

For a bounded continuous positive-time cylinder function \(F(Z_a)\),
(8) makes it a function of the positive-time diffusion path.
Time reversal of \(X_a\) induces exactly (11). Reversibility and
conditional independence of past and future given \(X_a(0)\) imply
\[
 \mathbb E\big[\overline{F(\Theta_ZZ_a)}F(Z_a)\big]
 =\mathbb E\big|\mathbb E[F(Z_a)\mid X_a(0)]\big|^2\ge0 .
 \tag{14}
\]
Weak convergence passes this inequality. Uniform moments extend it
to polynomial cylinders. \(\square\)

For \(n\) tests and \(p=\lceil n/2\rceil\), Hölder and (10) give
\[
 \left|\mathbb E\prod_{i=1}^nZ_a(J_i)\right|
 \le[2p!]^{n/(2p)}(4c)^{n/2}\prod_i\|J_i\|_2 .
 \tag{15}
\]
These estimates also give uniform integrability, so all such moments
pass along the same field-law subsequence. They are sub-Gaussian
estimates, not a Gaussianity or independence claim.

The same proof applies to YM77's full three-orientation observable
\(\Phi_a\), with \(c\) in (9)–(15) replaced by \(4c\).
It preserves that observable's spatial cubic symmetries. The selected
single-orientation field \(\Psi_a\) is not claimed to be rotationally
invariant, or identical to the full field's limit. Neither construction
establishes rotations mixing time with space.

## 3. Actual positive reflected norms at each regulator

Put
\[
 D_a(T,f)=\|(1-e^{-T\mathsf H_a})v_a(f)\|^2,\qquad T>0 .
 \tag{16}
\]
For the anchored \(\Psi_a\)-increment process this is its reflected
two-time covariance, as in YM77.

### Theorem 78.3 — an explicit lower certificate, with its cutoff loss

For any nonnegative self-adjoint operator \(H\), \(v\in D(H)\), and
\(q=\langle v,Hv\rangle>0\), \(m_0=\|v\|^2\), \(m_2=\|Hv\|^2\),
\[
 \boxed{\|(1-e^{-TH})v\|^2
       \ge\frac{T^2q^2}{m_0+Tq+T^2m_2}.}
 \tag{17}
\]
For the actual selected source and \(\theta_a\ge25\), YM76 implies
\[
 m_0\le\frac{9\pi^2}{32a^2}\|f\|_1^2=:K_fa^{-2},\qquad
 m_2\le\frac{10c^2}{a^4}\|f\|_1^2=:R_fa^{-4}.
 \tag{18}
\]
Thus the exact \(q_a\) of (5) gives
\[
 \boxed{D_a(T,f)\ge
 \frac{T^2q_a(f)^2}
      {K_fa^{-2}+Tq_a(f)+T^2R_fa^{-4}}>0}
 \tag{19}
\]
whenever \(P_af\ne0\). Formula (6) makes this strict for every fixed
nonzero \(f\) and all sufficiently small \(a\).
For fixed \(T,f\), the right side of (19) is asymptotic to
\(q_0(f)^2a^4/R_f\), where \(q_0(f)=c\|f\|_2^2>0\).
It therefore proves no positive continuum lower bound.

**Proof.** For \(x\ge0\),
\[
 (1-e^{-x})^2\ge\frac{x^2}{1+x+x^2}.
 \tag{20}
\]
Indeed \(e^x\ge1+x+x^2/2\), and the remaining rational comparison is
\[
 (1+x/2)^2(1+x+x^2)-(1+x+x^2/2)^2
          =x^2(1+x)/4\ge0.
\]
Apply spectral calculus, then weighted Cauchy–Schwarz:
\[
 q^2\le
 \left(\int\frac{E^2}{1+TE+T^2E^2}\,d\nu_v(E)\right)
 \left(\int(1+TE+T^2E^2)\,d\nu_v(E)\right).
\]
This proves (17). For (18), the actual single-angle bounds are
\[
 \|(r_b-\mu_ar_b)h_a\|^2\le\frac{9\pi^2}{8\sqrt{\theta_a}},
 \quad
 \|H_a(r_b-\mu_ar_b)h_a\|^2\le40\sqrt{\theta_a}.
\]
Use the triangle inequality and
\(\sum_b|\bar f_b|\le\|f\|_1/(4a^3)\) in (4), followed by (1).
The radial domain statement in YM76 ensures \(v_a\in D(\mathsf H_a)\).
Substitution proves (19). \(\square\)

The exact rational replay uses \(9<\pi^2<10\). It can replace
\(K_f\) by \(45\|f\|_1^2/16\), and the lower energy in (6) by
\(c(\|f\|_2^2-4a^2\|\nabla f\|_2^2/9)_+\).
When using an energy interval, substitute its lower endpoint in the
numerator and upper endpoint in the denominator.

## 4. Where the positive source energy goes

Fix \(f\ne0\) and write the physical spectral measures
\[
 \nu_a(dE)=\|P_{\mathsf H_a}(dE)v_a(f)\|^2,\qquad
 \eta_a(dE)=E\,\nu_a(dE),\quad \eta_a([0,\infty))=q_a(f).
 \tag{21}
\]
At each regulator \(\eta_a\) has no atom at zero. Its total mass now
converges to the known positive value \(c\|f\|_2^2\).
This fact alone does not locate that mass in a fixed physical band.

### Theorem 78.4 — contact covariance and a rational survival criterion

The stationary derivative covariance in time is the distribution
\[
 C_a(t)=2q_a\delta(t)-\int_{(0,\infty)}Ee^{-E|t|}\,d\eta_a(E).
 \tag{22}
\]
With \(\widehat j(\omega)=\int e^{-i\omega t}j(t)dt\), its multiplier is
\[
 K_a(\omega)=2\int\frac{\omega^2}{E^2+\omega^2}\,d\eta_a(E),
 \quad
 \mathbb E|Z_a(j\otimes f)|^2
       =\frac1{2\pi}\int K_a(\omega)|\widehat j(\omega)|^2d\omega .
 \tag{23}
\]
For a real positive-time test \(j\), the **odd-reflected** linear norm is
\[
 R_a(j,f)=\int E\left|\int_0^\infty e^{-Et}j(t)dt\right|^2d\eta_a(E).
 \tag{24}
\]
For a fixed physical frequency \(s>0\), define
\[
 \Delta_a(s,f)=K_a(2s)-K_a(s)
 =\int\frac{6(E/s)^2}{[(E/s)^2+1][(E/s)^2+4]}\,d\eta_a(E).
 \tag{25}
\]
Then
\[
 \boxed{D_a(1/s,f)\ge\frac{\Delta_a(s,f)}{3s}.}
 \tag{26}
\]
At every regulator there is also an exact all-scale sum rule,
\[
 \boxed{\sum_{k\in\mathbb Z}\Delta_a(2^ks,f)=2q_a(f).}
 \tag{27}
\]
A positive lower bound on \(\Delta_a(s,f)\) at one fixed physical
frequency would therefore prove survival. Equation (27) does not
provide that bound: its mass can move to \(k\to+\infty\) or
\(k\to-\infty\).

**Proof.** The undifferentiated covariance is
\(\int e^{-E|t|}d\nu_a(E)\). Differentiation in both time arguments
gives its negative second derivative, proving (22)–(23).
For disjoint positive and negative times the contact term disappears.
The minus sign in (11) changes the remaining negative kernel into
\(\int Ee^{-E(t+u)}d\eta_a(E)\), which factorizes to (24).

For \(x\ge0\),
\[
 \frac{6x^2}{(x^2+1)(x^2+4)}
 \le\frac{3x}{(1+x)^2}
 \le\frac{3(1-e^{-x})^2}{x}.
 \tag{28}
\]
The value at \(x=0\) is interpreted continuously.
The first inequality reduces to
\(x^2(x-1)^2+4-2x\ge0\): use \(4-2x\ge0\) for \(x\le2\),
and \(x^2(x-1)^2\ge x^2\ge2x\) for \(x\ge2\).
The second follows from \(e^x\ge1+x\).
Integrate (28) against \(\eta_a\), with \(x=E/s\), to obtain (26).
Finally, (25) telescopes over dyadic frequencies; dominated convergence
gives \(K_a(0)=0\) and \(K_a(\infty)=2q_a\), proving (27).
\(\square\)

**The exact endpoint distinction.** Pass along a subsequence for which
\(\eta_a\) converges weakly on the compactification \([0,\infty]\):
\[
 \bar\eta=m_0^{\rm end}\delta_0+\eta_*+m_\infty^{\rm end}\delta_\infty .
 \tag{29}
\]
Here endpoint masses are not the \(m_0\) source norm in (17).
For \(\omega\ne0\), the limiting multiplier is
\[
 K_*(\omega)=2m_0^{\rm end}
       +2\int_{(0,\infty)}\frac{\omega^2}{E^2+\omega^2}\,d\eta_*(E).
 \tag{30}
\]
Mass at infinity is invisible. Mass at zero contributes temporal
white-noise covariance. For \(j\ge0\), nonzero and supported in a
compact subset of \((0,\infty)\), the kernel in (24) is strictly
positive on \((0,\infty)\) and tends to zero at both endpoints.
The same is true of the kernels for (16) and (25).
Consequently, along this subsequence, the following are equivalent:

- the interior measure \(\eta_*\) is nonzero;
- the limiting \(D_a(T,f)\) is positive, for any fixed \(T>0\);
- the limiting contrast \(\Delta_a(s,f)\) is positive, for any fixed \(s>0\);
- the limiting reflected linear norm in (24) is positive.

These statements concern this selected linear observable sector.
Vanishing of all its linear norms alone does not classify the entire
nonlinear reflected algebra.

For sharp controls, consider one-dimensional spectral fixtures with
unit energy measure: \(E_n=n,\ \|v_n\|^2=1/n\), or
\(E_n=1/n,\ \|v_n\|^2=n\). Both have \(q_n=1\) and both have
\(D_n(T)\to0\) and \(\Delta_n(s)\to0\).
The first has \(K_n(\omega)\to0\). The second has
\(K_n(\omega)\to2\) for \(\omega\ne0\).
A stationary Gaussian Ornstein–Uhlenbeck realization of the second
example converges after differentiation to temporal white noise,
whose positive and negative half-lines are independent and whose
reflected Hilbert space is vacuum only.
These are counterexamples to proposed inferences, not substitutes
for actual Yang–Mills vacua.

### A resolvent witness for the next computation

There is a direct connection to YM75's source-resolved resolvent method.
For any nonnegative self-adjoint \(H\) and \(v\in D(H^{1/2})\), define
\[
 B_T=TH(I+TH)^{-1},\quad
 \mathcal R_T(v)=\|B_Tv\|^2 .
\]
Scalar inequalities give the safe comparison
\[
 \mathcal R_T(v)\le\|(1-e^{-TH})v\|^2\le4\mathcal R_T(v).
 \tag{31}
\]
For the upper comparison use
\(1-e^{-x}\le\min(x,1)\le2x/(1+x)\).
Moreover,
\[
 \mathcal R_T(v)
 =\sup_{z\in D(H)}
 \left\{2T\operatorname{Re}\langle H^{1/2}v,H^{1/2}z\rangle
                   -\|(I+TH)z\|^2\right\}.
 \tag{32}
\]
Indeed \((I+TH)D(H)\) is the whole Hilbert space, and
the expression is \(2\operatorname{Re}\langle B_Tv,b\rangle-\|b\|^2\)
with \(b=(I+TH)z\). Completing the square proves (32).
A certified positive value for one witness \(z\), uniformly at a
fixed physical \(T\), would discharge survival without controlling
the divergent unfiltered norm of \(v\).
No such actual uniform witness has yet been supplied.

## 5. Local sources exclude the charged winding obstruction

### Theorem 78.5 — exact neutral-sector reduction

Let \(Z_g\), \(g\in\mathbb Z_2^3\), be the three center-seam
transformations, and let
\[
 P_\chi=\frac18\sum_g\chi(g)Z_g
 \tag{33}
\]
be the character projections. They commute with the actual Hamiltonian.
Uniqueness and positivity put the vacuum in the neutral sector.
For every gauge-invariant observable \(O\) whose link support lies in
a contractible region,
\[
 P_\chi Oh_a=0\quad(\chi\ne0).
 \tag{34}
\]
Thus its full Hamiltonian cyclic subspace is neutral.

**Proof.** A flat center cocycle restricts on a contractible region
to a vertex-sign coboundary. On that region its action is an ordinary
vertex gauge transformation, under which \(O\) is invariant.
This proves (34); commutation with spectral calculus proves the
cyclic-space claim. \(\square\)

For a finite-dimensional retained space generated by such local sources,
the Schur coupling to every charged hidden sector is zero.
Define the hidden inverse directly on the neutral complement.
At a charged eigenvalue the full hidden inverse can fail to exist:
one must not write a full inverse there and cancel an undefined
quantity against a zero coupling.

In the actual side-four electric cut \(R=3\), the full retained
dimension \(241\) becomes neutral dimension \(193=1+192\);
48 winding directions are charged. The neutral hidden electric
floor remains \(9/2\), attained by a contractible six-edge rectangle.
The uniform two-dimensional block and 191 transverse directions
of YM75 remain unchanged. At \(\theta=1/48\), the same rational
inertia calculation gives
\[
 0.49069036\le\Delta_{\rm neutral,bare}\le3.00712209.
 \tag{35}
\]
These are bare Hamiltonian units; physical energies multiply by
\(\kappa_a\). This removes irrelevant charged directions and clarifies scope; it
does not improve the numerical enclosure or extend it to weak coupling.

The neutral sector still contains nonlocal adjoint and double-winding
observables. It is not identified with the local continuum state space.
Nor can YM77's seven charged norm bounds simply be paired to solve
local survival. If \(U=AB\), \(U'=A'B\) share a winding tail, the exact
SU(2) trace identity is
\[
 2W(U)W(U')-W(UU')=W(UU'^{-1})=W(AA'^{-1}).
 \tag{36}
\]
Tail cancellation returns an ordinary local loop when \(A,A'\)
are local alternatives. For one plaquette its centered norm obeys
\(\operatorname{Var}_{\mu_a}(W_b)\le2\mu_a(1-W_b)
\le9/(2\sqrt{\theta_a})\). The charged lower norm does not survive
this localization argument.

## 6. Scope, reproducibility and the remaining estimate

The advance over YM77 is an actual positive local source energy,
a stationary tempered derivative-field construction with a full
test-function moment hierarchy, and a positive finite-regulator
reflected-norm certificate. The latter explicitly loses four powers
of the spacing. The frequency contrast and resolvent witness say
what must improve to retain a nonzero reflected state.

For a chosen physical trajectory, the immediate unresolved target is
one fixed \(f\ne0,s>0\) and \(d>0\) such that
\[
 \inf_{a\text{ sufficiently small}}\Delta_a(s,f)\ge d .
 \tag{37}
\]
Even (37) would establish only a nonzero local linear sector.
It would not supply a uniform positive lower edge for the full
local excitation spectrum, rotations mixing time with space,
Yang–Mills short-distance identification, or all compact simple
gauge groups. Those are separate obligations.

The original problem statement requires a nontrivial local quantum
Yang–Mills theory on four-dimensional space-time with a positive
mass gap, together with the specified axiomatic and short-distance
properties [1]. The temperedness and reflection statements proved
here do not establish that full set of properties.

The deterministic certificate records the selected link geometry,
exact polynomial block averages and projection defects, rational
energy and norm budgets, exponential-comparison remainders,
frequency kernels and telescoping sums, endpoint-collapse controls,
and the neutral electric-cut replay. It pins 281 upstream inputs
and four local proof/code/test/workflow inputs. Upstream source
files and the shared engine are unchanged.

Run from the repository root:

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym78_stationary_derivative_energy_survival.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym78_stationary_derivative_energy_survival.py -v
~~~

The dedicated workflow also replays YM77, YM76 and YM75. A PASS
certifies these finite checks and source identities, not the open
continuum assertions.

## References

1. A. Jaffe and E. Witten, *Quantum Yang–Mills Theory*, official
   Clay Mathematics Institute problem statement, especially §4.
   [Primary text](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf).
2. K. Osterwalder and R. Schrader, *Axioms for Euclidean Green's
   Functions II*, Communications in Mathematical Physics **42**
   (1975), 281–305.
   [Publisher record](https://doi.org/10.1007/BF01608978).
   Cited for the reconstruction framework, not as a claim that
   this stage verifies its full hypotheses.
3. L. Wu, *Forward-backward martingale decomposition and compactness
   results for additive functionals of stationary ergodic Markov
   processes*, Annales de l'Institut Henri Poincaré, Probabilités
   et Statistiques **35** (1999), 121–141.
   [Primary archive](https://www.numdam.org/item/AIHPB_1999__35_2_121_0/).
   The particular time-dependent estimate used here is derived
   directly in (12).
4. [YM75](YM75_RH_SOURCE_RESOLVED_SPECTRAL_TRANSFER.md),
   [YM76](YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS.md),
   [YM77](YM77_PHYSICAL_SMEARING_AND_CONTINUUM_SURVIVAL.md).
   Their exact source identities and analytic hypotheses are retained.
