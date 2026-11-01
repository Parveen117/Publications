# YM-80: bounded phases, energy budgets and ultraviolet escape

Monty Dabas. 5 October 2026. **Python 3.12 only.**

This stage proves an order-\(a\) lower certificate for the reflected norm
of a new bounded, smooth, local phase observable in the actual
untruncated SU(2) model. Its physical energy stays bounded as
\(a\to0\). The proof combines a stronger retained-energy estimate,
an exact phase chain rule and the local third-moment estimate from
YM79. Adjusting the regularizer improves the leading certified
coefficient by more than a factor of five.

The new observable is different from the raw linear source whose
certificate improved from \(a^4\) to \(a^2\) in YM79.
**The order-\(a\) certificate is still not a uniform positive lower
bound.** An exact comparison diffusion shows that the available
moment data cannot yield one. A further actual-model theorem shows
that the explicit amplitude scaling \(\lambda_a^2q_a=Bc/a\)
places a fixed fraction of the phase norm at energies of order \(1/a\).

The remaining task is made precise: control a vacuum-dependent
quadratic form, or establish positive spectral weight in a fixed
physical energy band. Neither a nontrivial four-dimensional
Yang–Mills continuum nor its positive mass gap is proved here.

## 1. Model, support and quantities held fixed

Keep the Hamiltonian, native SU(2) metric and actual positive vacuum
used in YM76–YM79:
\[
 A_\theta=L_E-\theta V,\quad L_E=-\Delta,\quad
 H=A_\theta-E_0,\quad \mathsf H_a=\kappa_a H,\quad
 \mu_a=h_a^2d\Phi,\quad \kappa_a\sqrt{\theta_a}=c/a.
 \tag{1}
\]
Here \(c=\sqrt{c_Ec_B}>0\). The spatial lattice is periodic, with
even \(L\ge4\); the continuum trajectory has \(a\to0\), \(aL\to\infty\)
and \(\theta_a\ge25\). No relation between \(L\) and \(\theta\) beyond
these conditions enters the estimates below.

For real \(f\in C_c^\infty(\mathbb R^3)\), use the selected \(xy\)
plaquettes rooted at \((2i,2j,k)\). Their four-link supports are
disjoint. The associated averaging cells have sides \(2a,2a,a\).
Let \(P_af\) be the block-constant orthogonal projection and set
\[
 d_b=2\theta_a^{1/4}a^2\bar f_b,\qquad
 q_a=\kappa_a\sum_bd_b^2=c\|P_af\|_2^2
       \longrightarrow q_0=c\|f\|_2^2.
 \tag{2}
\]
Only cells meeting the fixed support of \(f\) enter. For \(f\ne0\),
\(q_a>0\) eventually, and \(q_a\le c\|f\|_2^2\).

For a regularizer parameter \(0<\eta\le1\), define
\[
 \epsilon_a=\eta\theta_a^{-1/4},\quad
 \rho_\epsilon(r)=
   \arccos\!\left(\frac{\cos r}{\sqrt{1+\epsilon^2}}\right),\quad
 G_a=\sum_bd_b\{\rho_{\epsilon_a}(r_b)
                         -\mu_a\rho_{\epsilon_a}(r_b)\},
 \quad r_b=\arccos W_b.
 \tag{3}
\]
The function \(\rho_\epsilon(\arccos W)\) is smooth in every link,
including \(W=\pm1\). Thus \(G_a h_a\) and
\(\exp(i\lambda G_a)h_a\) lie in every finite Hamiltonian domain
at each regulator. These observables are gauge invariant and
center neutral, and retain the fixed physical spatial support.

Write \(\Gamma_a(F)=\kappa_a|\nabla F|^2\). For real \(G_a\), put
\[
 q_G=\mu_a\Gamma_a(G_a),\quad u_a=q_G/q_a,\quad
 m_j(G)=\langle G_ah_a,\mathsf H_a^jG_ah_a\rangle.
 \tag{4}
\]
For any vector \(w\), let
\[
 D_T(w)=\|(1-e^{-T\mathsf H_a})w\|^2,\qquad T>0.
 \tag{5}
\]
All times and energies below are physical unless explicitly stated.

## 2. Retained energy and a better regularizer

### Theorem 80.1 — actual smooth-source moment bounds

Set
\[
 \beta_\eta=\frac{1}{1+(183/5)\eta^2},\qquad
 C_\eta=8568\eta^{-2}+419280\eta^{-4}+8+24\eta^{-1}.
 \tag{6}
\]
For every regulator in Section 1,
\[
 0\le\Gamma_a(G_a)\le q_a,\qquad
 \beta_\eta\le u_a\le1,\qquad
 m_3(G)\le C_\eta\,\frac{c^2}{a^2}\,q_a.
 \tag{7}
\]
If \(B^2\ge C_\eta\) and \(r_a=Bc/a\), spectral Cauchy–Schwarz gives
the correlated bound
\[
 \boxed{m_2(G)\le r_a q_a\sqrt{u_a}.}
 \tag{8}
\]

**Proof.** Disjoint link supports and
\[
 \rho_\epsilon'(r)=
 \frac{\sin r}{\sqrt{\sin^2r+\epsilon^2}}
\]
give the pointwise bound and the diagonal formula
\(q_G=\kappa_a\sum_b d_b^2\mu_a(\rho_\epsilon'(r_b))^2\).
YM79, equation (12), proved in the actual vacuum that
\[
 \mu_a\csc^2r_b\le(183/5)\sqrt{\theta_a}.
\]
Put \(X=\epsilon_a^2\csc^2r_b\). The function \(x\mapsto(1+x)^{-1}\)
is convex on \([0,\infty)\), so
\[
 \mu_a(\rho_\epsilon')^2
 =\mu_a(1+X)^{-1}
 \ge(1+\mu_aX)^{-1}\ge\beta_\eta.
 \tag{9}
\]
The poles have zero measure and the displayed Hardy bound gives the
required integrability.

YM79's local third-moment identity and covariant derivative bounds
apply for \(\epsilon\le1\). This holds here because
\(\epsilon^2\le\eta^2/5\le1/5\).
Before substituting \(\eta=1/16\), that proof bounds the bare third
moment by
\[
 [9K\eta^{-2}+Z\eta^{-4}+8+24\eta^{-1}]
       \theta_a\sum_bd_b^2,\quad K=952,\quad Z=419280.
\]
Multiplication by \(\kappa_a^3\), followed by (1) and (2), proves
(7). Finally \(m_2(G)^2\le m_1(G)m_3(G)\) and
\(m_1(G)=u_aq_a\) prove (8). No independence, clustering or
uniform gap enters this argument. \(\square\)

The explicit choice used in this stage is
\[
 \boxed{\eta=\frac7{30},\quad
 \beta=\frac{1500}{4489},\quad B=12000.}
 \tag{10}
\]
Indeed,
\[
 C_{7/30}=\frac{339994914968}{2401}<12000^2,\qquad
 \sqrt\beta\le\frac{29}{50}.
 \tag{11}
\]
The choice is convenient and certified; no claim of a global
optimum over regularizers is made.

For comparison, applying the same Jensen argument to YM79's
\(\eta=1/16\) gives \(\beta_0=1280/1463\) and its existing
\(B_0=2^{18}\). The new leading coefficient obeys
\[
 \left[
 \frac{\beta^{3/2}/12000}{\beta_0^{3/2}/2^{18}}
 \right]^2
 =\frac{2404884362496}{90458382169}>25.
 \tag{12}
\]
Thus it improves the old-regularizer phase certificate by more
than fivefold. This compares two phase certificates; it does not
compare norms of different states.

## 3. A bounded phase with bounded physical energy

### Theorem 80.2 — actual order-\(a\) survival certificate

Define
\[
 U_{a,\lambda}=e^{i\lambda G_a},\quad
 w_{a,\lambda}=(U_{a,\lambda}-\mu_aU_{a,\lambda})h_a,\quad
 e_a=\lambda^2q_a.
 \tag{13}
\]
For \(\lambda\ne0\), \(q_a>0\), and \(T>0\),
\[
 \boxed{
 D_T(w_{a,\lambda})\ge
 \delta_a(T,e_a):=
 \frac{\beta^2T^2e_a^2}
 {1+\beta Te_a+\beta T^2e_a^2+
              \sqrt\beta\,r_aT^2e_a}.}
 \tag{14}
\]
The equal-time centered norm also satisfies
\[
 \boxed{\|w_{a,\lambda}\|^2\ge
 \frac{\beta^2e_a}{\beta e_a+\sqrt\beta\,r_a}.}
 \tag{15}
\]
The right sides may be replaced by the smaller, entirely rational
expressions obtained by putting \(29/50\) in place of \(\sqrt\beta\).

**Proof.** Since \(U\) has modulus one,
\[
 m_0(w)=1-|\mu_aU|^2\le1,\qquad m_1(w)=e_au_a.
 \tag{16}
\]
The ground-state chain rule is
\[
 \mathsf H_a(Uh_a)=
 U\{\lambda^2\Gamma_a(G_a)h_a+
                      i\lambda\,\mathsf H_a(G_ah_a)\}.
 \tag{17}
\]
The two terms inside braces are real and purely imaginary,
respectively. Consequently their norm cross term vanishes:
\[
 m_2(w)=\lambda^4\mu_a\Gamma_a(G_a)^2+\lambda^2m_2(G)
 \le e_a^2u_a+e_ar_a\sqrt{u_a}.
 \tag{18}
\]
Centering does not change positive spectral moments.
The inequality used in the first term is
\(\Gamma_a(G_a)^2\le q_a\Gamma_a(G_a)\); using only
\(\Gamma_a(G_a)^2\le q_a^2\) would discard useful information.

YM78's scalar moment certificate gives
\[
 D_T(w)\ge\frac{T^2m_1(w)^2}
 {m_0(w)+Tm_1(w)+T^2m_2(w)}
 \ge\frac{T^2e_a^2u_a^2}
 {1+(Te_a+T^2e_a^2)u_a+T^2e_ar_a\sqrt{u_a}}.
 \tag{19}
\]
For \(A,C\ge0\),
\[
 \frac{d}{du}\frac{u^2}{1+Au+C\sqrt u}
 =\frac{u(2+Au+\tfrac32C\sqrt u)}
             {(1+Au+C\sqrt u)^2}>0\quad(u>0).
 \tag{20}
\]
We may therefore substitute \(u_a\ge\beta\) in the whole quotient;
substituting separately in a denominator would not be justified.
This proves (14). The spectral inequality
\(m_1^2\le m_0m_2\) gives (15), using the same monotonicity without
the constant term. \(\square\)

The certificate is increasing in \(e>0\): it has the form
\(ke^2/(1+be+de^2)\), whose derivative is
\(ke(2+be)/(1+be+de^2)^2>0\).
Under the guaranteed energy budget \(e_a\le E_*\), the best value
of this certificate is attained at \(e_a=E_*\). This is optimization
under the sufficient budget \(\lambda^2q_a\le E_*\), not under a
hypothetically known exact value of \(u_a\).

There is no need to adjust \(\lambda\) at each cutoff to obtain
the leading coefficient. The fixed amplitude
\[
 \lambda^2=\frac{E_*}{c\|f\|_2^2}
 \tag{21}
\]
satisfies \(e_a\le E_*\), \(e_a\to E_*\), and the actual phase energy
\(m_1(w)\le E_*\). At every fixed \(T>0\),
\[
 \delta_a(T,e_a)\sim\frac{\beta^{3/2}E_*}{Bc}\,a.
 \tag{22}
\]
This is the asymptotic of the lower certificate, not a claim that
the actual norm is asymptotic to a constant times \(a\).

The rational coefficient improves by more than fivefold too:
\[
 \frac{\beta^2/[12000(29/50)]}
      {\beta_0^2/[2^{18}(15/16)]}
 =\frac{12039575625}{2337530036}>5.
 \tag{23}
\]
Changing the regularizer has not established equality of the
raw-angle and smooth-phase continuum limits. For example, the
direct estimate
\[
 \|F_a-G_a\|_\infty\le
       2\epsilon_a\sum_b|d_b|
       \le\frac{\eta\|f\|_1}{a}
 \tag{24}
\]
does not tend to zero. The factor two allows for their distinct
centerings. No such identification is used above.

## 4. Why these moment bounds cannot close the continuum problem

### Proposition 80.3 — sharpness as a statement about the input data

The exponent in (22), and its coefficient with exact square roots,
are sharp for the relaxed moment data used in Theorem 80.2:
\(0\le\Gamma(G)\le q\), \(\beta\le m_1(G)/q\le1\),
\(m_3(G)\le r_a^2q\), (8), and the phase chain rule.
Here \(r_a=Bc/a\); the replacement of \(C_\eta\) by \(B^2\)
is part of this relaxation. This is not a statement about the
Yang–Mills vacuum or sharpness under the smaller exact constant
\(C_\eta\) in (7).

To see this, on the real line take the positive Ornstein–Uhlenbeck
generator on \(L^2(\mu_a)\)
\[
 K_a=-\beta q\,\partial_x^2+\omega_ax\partial_x,\quad
 \mu_a=N(0,\beta q/\omega_a),\quad
 G(x)=x,\quad \omega_a=r_a/\sqrt\beta.
 \tag{25}
\]
Its vacuum is the constant function and
\[
 \Gamma(G)=\beta q,\quad
 m_1(G)=\beta q,\quad m_2(G)=r_aq\sqrt\beta,\quad
 m_3(G)=r_a^2q.
 \tag{26}
\]
For \(e=\lambda^2q\), the phase has Poisson spectral weights at
energies \(n\omega_a\), with mean
\[
 s_a=\frac{\beta^{3/2}e}{r_a}.
\]
Centering removes the zero-energy atom. Its centered norm and
filtered norm are
\[
 m_0(w)=1-e^{-s_a},\qquad
 D_T(w)=1-2e^{-s_a(1-e^{-T\omega_a})}
                +e^{-s_a(1-e^{-2T\omega_a})}.
 \tag{27}
\]
For fixed \(T,e>0\), \(D_T(w)\sim s_a\). Uniformly over
\(0\le e\le E_*\), it is bounded above by
\(\beta^{3/2}E_*/r_a\to0\).
Thus these inputs alone cannot prove uniform positive survival.
Even if the sharper third-moment constant in (7) is retained,
the same construction with \(r_a=\sqrt{C_\eta}\,c/a\) proves the
same order-\(a\) obstruction; it does not saturate the particular
coefficient obtained using \(B\).

The finite certificate uses directed rational enclosures for the rate
\(r_a/\sqrt\beta\) at the selected value \(\beta=1500/4489\),
and also checks the unit-retention comparison. Independent Gaussian,
Hermite and Poisson controls are compared with decimal evaluations.
The tests additionally include the rational parameter \(u=1/4\),
\(\omega=2r\), as a formula check at a different retained-energy
parameter, not as a model satisfying the selected lower bound on \(u\).
They are explicitly labelled comparison diffusions, never actual
interacting Yang–Mills vacua.

## 5. Increasing the amplitude forces ultraviolet norm in the actual model

### Theorem 80.4 — a fixed norm packet at diverging physical energy

Choose the cutoff-dependent amplitude
\[
 e_a=\lambda_a^2q_a=r_a=Bc/a.
 \tag{28}
\]
Let \(\nu_a\) be the spectral measure of the centered actual phase
state \(w_{a,\lambda_a}\). Then
\[
 \boxed{
 \nu_a\!\left([\beta r_a/2,\;8r_a/\beta]\right)
     \ge d_{\rm UV}:=\frac{3\beta^2}{32}>0.}
 \tag{29}
\]
In particular, for every fixed \(t>0\),
\[
 D_t(w_{a,\lambda_a})\ge
 d_{\rm UV}(1-e^{-t\beta r_a/2})^2,\qquad
 \liminf_{a\to0}D_t(w_{a,\lambda_a})\ge d_{\rm UV}.
 \tag{30}
\]

**Proof.** Equations (16)–(18) give
\[
 m_0\le1,\qquad m_1\ge\beta r_a,\qquad m_2\le2r_a^2.
\]
For \(\ell=\beta r_a/2\), spectral Cauchy–Schwarz implies
\[
 m_1\le\ell\,m_0+\sqrt{m_2\,\nu_a([\ell,\infty))},
\]
hence \(\nu_a([\ell,\infty))\ge\beta^2/8\).
For \(M=8r_a/\beta\), Markov's inequality gives
\(\nu_a((M,\infty))\le m_2/M^2\le\beta^2/32\).
Subtracting proves (29), and the increasing filter proves (30).
\(\square\)

This rules out retaining the whole time-zero state's norm in a
strongly continuous physical-time Hilbert-space limit.
Specifically, if the norms and all heat correlations of these
vectors converged to those of a vector \(w\) for a nonnegative
self-adjoint continuum Hamiltonian, then
\(D_t(w_{a,\lambda_a})\to\|(1-e^{-tH})w\|^2\).
Equation (30) would keep the latter at least \(d_{\rm UV}\)
for every \(t>0\), contradicting strong continuity at \(t=0\).

A continuum construction can discard this ultraviolet norm.
The theorem neither excludes a separate surviving low-energy
component nor proves that all actual phase norm escapes.

A comparison with complete escape is explicit: take (25) with
unit retention, \(\omega_a=r_a\), and \(e_a=r_a\).
Its phase spectral law is \(r_a\) times a Poisson variable of
mean one. The centered norm and \(D_T\) tend to \(1-e^{-1}>0\),
but
\[
 \|e^{-sK_a}(1-e^{-TK_a})w_a\|^2\longrightarrow0
 \quad(s,T>0).
 \tag{31}
\]
The time-zero norm is therefore insufficient without control
at positive time separation.

## 6. A positive-time criterion and the precise missing vacuum estimate

### Proposition 80.5 — separated-time survival

For any centered state with \(\nu_a([0,\infty))\le1\), define
\[
 S_{s,T}(w_a)=
 \|e^{-s\mathsf H_a}(1-e^{-T\mathsf H_a})w_a\|^2
 =\int e^{-2sE}(1-e^{-TE})^2\,\nu_a(dE),
 \quad s,T>0.
 \tag{32}
\]
The integrand vanishes at both \(E=0\) and \(E=\infty\).
Along weak convergence of the bounded measures on the compactified
half-line, its limit is positive exactly when the limiting measure
has positive mass in \((0,\infty)\).

More quantitatively, fixed-band mass
\(\nu_a([\ell,M])\ge b>0\), with \(0<\ell<M<\infty\),
implies
\[
 S_{s,T}\ge b\,e^{-2sM}(1-e^{-T\ell})^2.
 \tag{33}
\]
Conversely, \(S_{s,T}\ge d>0\) implies at least \(d/2\) mass in
a fixed band: take
\[
 M=\frac{\log(4/d)}{2s},\qquad
 \ell=\min\{\sqrt d/(2T),M/2\}.
\]
Below \(\ell\) the integrand is at most \(T^2\ell^2\le d/4\);
above \(M\) it is at most \(e^{-2sM}=d/4\).
The surviving band contribution is at least \(d/2\), and the
integrand is at most one.

If \(D_T(w_a)\ge d>0\) and \(m_1(w_a)\le E_*\), a useful direct
bound is
\[
 \boxed{S_{s,T}(w_a)\ge d\,\exp(-2sE_*/d).}
 \tag{34}
\]
Indeed, the filtered measure
\((1-e^{-TE})^2\nu_a(dE)\) has mass \(D_T\) and first moment
at most \(E_*\). Jensen's inequality for \(e^{-2sE}\) gives
\(S_{s,T}\ge D_T\exp(-2sE_*/D_T)\).
This expression increases with \(D_T>0\), proving (34).

For our present certificate \(d=\delta_a\), (34) is positive
at each cutoff but its lower bound still tends to zero.
A uniform positive \(d\), together with a bounded physical
energy, would suffice. Neither the current order-\(a\) bound
nor the high-amplitude ultraviolet packet supplies it.

### Proposition 80.6 — the vacuum quadratic form that remains to be controlled

Write \(\chi=\log h_a\) and let
\[
 \mathcal L_\mu=-\Delta-2\nabla\chi\cdot\nabla
\]
be the positive ground-state generator in bare units.
For smooth real \(F\), weighted integration by parts and the
Riemannian Bochner identity give the exact formula
\[
 \|\mathcal L_\mu F\|_{L^2(\mu_a)}^2
 =\mu_a\!\left[
 \|\operatorname{Hess}F\|_{\rm HS}^2+
 (\operatorname{Ric}-2\operatorname{Hess}\chi)
                              (\nabla F,\nabla F)\right].
 \tag{35}
\]
For completeness, apply the usual identity
\[
 \tfrac12(\Delta+2\nabla\chi\cdot\nabla)|\nabla F|^2
 =\langle\nabla F,\nabla(\Delta+2\nabla\chi\cdot\nabla)F\rangle
  +\|\operatorname{Hess}F\|^2+
  (\operatorname{Ric}-2\operatorname{Hess}\chi)(\nabla F,\nabla F),
\]
integrate against the invariant measure, and integrate the first
term once more. This proves (35), including its signs.

Let \(\phi_b=\rho_{\epsilon_a}(r_b)\) and
\[
 Q_{bc}=
 \langle\mathcal L_\mu\phi_b,\mathcal L_\mu\phi_c\rangle_{L^2(\mu_a)}.
\]
The matrix \(Q\) is positive semidefinite. For distinct selected
blocks, product geometry makes the ordinary Hessian and Ricci
cross terms vanish, so polarization of (35) gives
\[
 \boxed{Q_{bc}=
 -2\mu_a\operatorname{Hess}(\log h_a)
                 (\nabla\phi_b,\nabla\phi_c)\quad(b\ne c).}
 \tag{36}
\]
The exact physical second moment is
\[
 m_2(G)=\kappa_a^2\sum_{b,c}d_bd_cQ_{bc}.
 \tag{37}
\]
A finite-range plaquette potential does not by itself make the
vacuum Hessian in (36) finite range. Disjoint source supports
therefore do not justify discarding these cross terms.

One concrete sufficient next estimate, for a fixed nonzero \(f\),
would be
\[
 \kappa_a^2d^{\mathsf T}Qd\le M_f<\infty
 \quad\hbox{uniformly along the continuum trajectory.}
 \tag{38}
\]
If \(0<q_*\le q_a\le Q_*\) and (38) holds, then every fixed
\(\lambda\ne0\) satisfies the uniform lower bound
\[
 D_T(w_{a,\lambda})\ge
 \frac{\beta^2T^2\lambda^4q_*^2}
 {1+T\lambda^2Q_*+T^2\lambda^4Q_*^2+T^2\lambda^2M_f}>0.
 \tag{39}
\]
Equation (34) would then give positive separated-time survival.
This is a sufficient criterion, not a necessary one, and (38)
has not been proved. Even proving it would establish survival
for this source family, not a gap above every nonvacuum state
or identification of the required Yang–Mills fields.

## 7. Reproduction and scope

The exact certificate replays the retained-energy constants,
regularizer comparisons, phase moment bounds, energy-budget
optimization, ultraviolet packet arithmetic, separated-time
bounds and labelled Gaussian/Poisson controls. Directed rational
bounds are used for exponentials and square roots; independent
decimal tests check selected enclosures.

The written proofs establish the analytic statements. Finite
arithmetic controls do not stand in for universal analytic proofs.
The source manifest pins the 295 upstream inputs through YM79
and the four local manuscript/code/test/workflow inputs.
The shared engine and all earlier stage files remain unchanged.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym80_bounded_phase_energy_certificates.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym80_bounded_phase_energy_certificates.py -v
~~~

The dedicated workflow also replays YM75–YM79 certificates and tests.

| Statement | Status in this stage |
|---|---|
| Actual smooth local phase, bounded physical energy | Proved |
| Order-\(a\) reflected-norm lower certificate | Proved |
| More than fivefold coefficient improvement | Proved for the stated certificates |
| Fixed ultraviolet norm packet at amplitude squared proportional to \(1/a\) | Proved in the actual model |
| Positive survival at fixed physical time separation, uniformly in \(a\) | Open |
| Uniform vacuum quadratic-form estimate (38) | Open; sufficient criterion |
| Four-dimensional local Yang–Mills field construction | Open |
| Full Euclidean covariance and positive continuum mass gap | Open |

Sources: the proofs use the unchanged local geometry and vacuum
estimates in [YM76](YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS.md),
the physical source construction in
[YM77](YM77_PHYSICAL_SMEARING_AND_CONTINUUM_SURVIVAL.md) and
[YM78](YM78_STATIONARY_DERIVATIVE_AND_ENERGY_SURVIVAL.md), and the
smooth third-moment argument in
[YM79](YM79_LOCAL_WITNESSES_AND_SURVIVAL_CERTIFICATES.md).
The target remains the
[Jaffe–Witten problem description](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf).
Clay's [current problem page](https://www.claymath.org/millennium/yang-mills-the-maths-gap/)
lists the mass-gap problem as unsolved (checked 5 October 2026).
