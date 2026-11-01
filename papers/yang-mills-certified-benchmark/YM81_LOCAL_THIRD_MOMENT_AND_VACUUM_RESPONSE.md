# YM-81: local third moments and vacuum response

Monty Dabas. 5 October 2026. **Python 3.12 only.**

This stage resolves the spatial structure of the moment estimate
left open in YM80. For the actual smooth plaquette sources, the
third-moment matrix is exactly nearest neighbor, although the
second-moment matrix contains nonlocal vacuum correlations.
Its spatial variation has a uniform physical bound. The unresolved
part of this sufficient continuum-survival route is reduced to
one zero-momentum vacuum coefficient.

An independent exact identity controls the squared size of a
vacuum-Hessian row from above and below. Its physical size grows
as \(1/a^2\), so a uniform estimate cannot be obtained merely by
bounding every vacuum derivative in absolute value.

A two-resolvent filter then isolates positive physical energies.
It gives an actual order-\(a^3\) lower certificate for a positive-time
phase state and a concrete spectral-band certificate for the
actual side-four Wilson source.

**Uniform continuum survival and the four-dimensional Yang–Mills
mass gap remain unproved.** The new spatial estimate is unconditional;
boundedness of its zero-momentum coefficient is a separate,
unproved hypothesis. The finite side-four band moves with the
physical clock and is not a continuum mass-gap interval.

## 1. Fixed model and notation

Use the unchanged SU(2) Hamiltonian, its positive normalized vacuum,
and the native product metric:
\[
 A_\theta=L_E-\theta V,\quad L_E=-\Delta,\quad
 H=A_\theta-E_0,\quad Hh=0,\quad
 \mu=h^2d\Phi,\quad \chi=\log h.
 \tag{1}
\]
The physical Hamiltonian is \(\mathsf H_a=\kappa_aH\), with
\[
 \kappa_a\sqrt{\theta_a}=c/a,\qquad c=\sqrt{c_Ec_B}>0.
 \tag{2}
\]
For the continuum estimates, take \(\theta_a\ge25\), \(a\to0\),
and \(aL\to\infty\). Each finite-volume vacuum is smooth and
strictly positive.

On link \(e\), let \(D_{e,\alpha}\), \(\alpha=1,2,3\), be the
orthonormal invariant fields with
\[
 [D_{e,\alpha},D_{e,\beta}]
   =\epsilon_{\alpha\beta\gamma}D_{e,\gamma},\quad
 \nabla_{D_{e,\alpha}}D_{e,\beta}
   =\tfrac12\epsilon_{\alpha\beta\gamma}D_{e,\gamma}.
 \tag{3}
\]
Different links commute and have zero cross connection.
This normalization has fundamental Casimir \(3/4\) and
\(\operatorname{Ric}=\tfrac12 I\) on each SU(2) link.

Translation and cubic symmetry make the following quantities
independent of the chosen link or plaquette:
\[
 \tau=\sum_\alpha\|D_{e,\alpha}h\|_2^2,\qquad
 m=\mu W_b.
 \tag{4}
\]
For \(\theta\ge25\), YM76 gives
\[
 m\ge1-\frac9{4\sqrt\theta},\qquad
 \frac{\sqrt\theta}{32}\le\tau\le\frac94\sqrt\theta.
 \tag{5}
\]
All Haar integrals and inner products below use the actual
untruncated model unless explicitly labelled comparisons.

## 2. Theorem 81.1 — an exact vacuum-Hessian row identity

For a fixed link \(e\), define
\[
 \mathcal R_e=
 \sum_f\mu\|\operatorname{Hess}_{f,e}\chi\|_{\rm HS}^2,
 \tag{6}
\]
where the sum is over all links \(f\), including \(f=e\), and the
Hilbert–Schmidt norm sums the three tangent directions on each link.
Then
\[
 \boxed{\mathcal R_e+\frac{\tau}{2}=\frac{3\theta}{2}m.}
 \tag{7}
\]
Consequently, for \(\theta\ge25\),
\[
 \boxed{\frac35\theta\le\mathcal R_e\le\frac32\theta,\qquad
 \frac{3c^2}{5a^2}\le\kappa_a^2\mathcal R_e
                           \le\frac{3c^2}{2a^2}.}
 \tag{8}
\]

**Proof.** The link Casimir commutes with every invariant link
field. Differentiating the ground-state equation gives
\[
 H(D_i h)=\theta(D_iV)h,\qquad i=(e,\alpha).
 \tag{9}
\]
Write \(\psi_i=D_i\chi\), so \(D_i h=h\psi_i\).
The ground-state transform and Haar integration by parts give
\[
 \begin{aligned}
 \mu\,\nabla\psi_i\cdot\nabla\psi_j
 &=\langle D_i h,H D_j h\rangle\\
 &=\theta\mu(\psi_iD_jV)
 =-\frac{\theta}{2}\mu D_iD_jV .
 \end{aligned}
 \tag{10}
\]
The last equality uses
\(\mu(\psi_i g)=-\tfrac12\mu(D_i g)\).
In particular this Gram entry is zero when distinct links
\(i,j\) share no plaquette. This locality concerns the Gram
matrix of derivatives of the score fields, not pointwise
locality of \(\operatorname{Hess}\chi\).

Summing \(i=j=(e,\alpha)\), and using the four plaquettes
incident to a link,
\[
 \sum_\alpha D_{e,\alpha}^2V
   =-\frac34\sum_{b\ni e}W_b,
\]
so
\[
 \sum_{f,\beta,\alpha}
      \mu(D_{f,\beta}D_{e,\alpha}\chi)^2
 =\frac{3\theta}{2}m.
 \tag{11}
\]
For \(f=e\), the ordered derivative differs from the covariant
Hessian by
\(\tfrac12\epsilon_{\beta\alpha\gamma}D_{e,\gamma}\chi\).
The Hessian is symmetric and this connection contribution is
antisymmetric. Their summed cross term is zero, and
\[
 \sum_{\alpha,\beta}
 \left|\tfrac12\epsilon_{\beta\alpha\gamma}
                       D_{e,\gamma}\chi\right|^2
 =\tfrac12|\nabla_e\chi|^2.
\]
For \(f\ne e\), there is no connection contribution.
Since \(\mu|\nabla_e\chi|^2=\tau\), equation (11) proves (7).

Using (5), \(m\le1\), and \(\tau\ge0\),
\[
 \frac32\theta-\frac92\sqrt\theta
 \le\mathcal R_e\le\frac32\theta.
\]
The lower endpoint is at least \(3\theta/5\) for \(\theta\ge25\).
Equation (2) proves the physical bounds. \(\square\)

The vectors \(D_i h\) used in the proof belong to the full link
Hilbert space; individually they need not be gauge invariant.
The row identity concerns the actual gauge-invariant vacuum.
It does not produce physical charged particles or a mass gap.
The squared row norm also supplies neither a sign nor spatial
decay for individual Hessian entries.

As a normalization check, fixed-volume perturbation at \(\theta=0\)
gives
\[
 h=1+\frac{\theta}{3}V+O(\theta^2),\quad
 m=\frac{\theta}{6}+O(\theta^2),\quad
 \tau=\frac{\theta^2}{12}+O(\theta^3),\quad
 \mathcal R_e=\frac{5\theta^2}{24}+O(\theta^3).
 \tag{12}
\]
These coefficients follow from native Haar orthogonality,
\(L_EV=3V\), and the four incident plaquettes. This check is
not the proof of the all-coupling identity.

## 3. Theorem 81.2 — the actual third-moment matrix is nearest neighbor

Take even \(L\ge6\). Select the \(xy\) plaquettes rooted at
\((2i,2j,k)\), and use their disjoint four-link supports.
The coarse translation steps in lattice units are
\[
 s_x=(2,0,0),\qquad s_y=(0,2,0),\qquad s_z=(0,0,1).
\]
Use YM80's smooth regularizer
\[
 \phi_b=\rho_\epsilon(r_b),\quad
 \rho_\epsilon(r)=
 \arccos\!\left(\frac{\cos r}{\sqrt{1+\epsilon^2}}\right),\quad
 \epsilon=\frac7{30}\theta^{-1/4},\quad r_b=\arccos W_b.
 \tag{13}
\]
Define the bare moment matrices
\[
 M^{(j)}_{bc}=
 \langle(\phi_b-\mu\phi_b)h,\,
      H^j(\phi_c-\mu\phi_c)h\rangle,\quad
 Q=M^{(2)},\quad \mathcal T=M^{(3)}.
 \tag{14}
\]
They are real, translation invariant and positive semidefinite.
In fact,
\[
 M^{(1)}=u I,\qquad
 u=\mu(\rho_\epsilon')^2\in[\beta,1],\quad
 \beta=\frac{1500}{4489}.
 \tag{15}
\]

The third-moment matrix has the exact form
\[
 \mathcal T_{bc}=0
 \quad\text{unless }c=b,\ b\pm s_x,\ b\pm s_y,\ b\pm s_z.
 \tag{16}
\]
Write its diagonal as \(A\), and its couplings in the three
positive directions as \(K_x,K_y,K_z\). Then
\[
 \boxed{|K_x|,|K_y|\le\theta/8,\qquad |K_z|\le\theta/2.}
 \tag{17}
\]
Its Fourier symbol and zero-momentum value are
\[
 \widehat{\mathcal T}(k)=
 A+2\sum_{i=x,y,z}K_i\cos k_i,\qquad
 t_0=A+2(K_x+K_y+K_z).
 \tag{18}
\]
Here \(t_0>0\) at every finite regulator. No sign is assumed for
the individual \(K_i\).

**Proof.** Set \(X_b=\nabla\phi_b\) and
\(J_b=X_b\cdot\nabla+\tfrac12\Delta\phi_b\).
These operators are skew on Haar space, and
\(H(\phi_bh)=-2J_bh\).
Polarizing YM79's double-commutator identity gives, for \(b\ne c\),
\[
 \mathcal T_{bc}
 =-2\langle h,[J_b,[H,J_c]]h\rangle
 =-2\theta\,\mu(X_bX_cV).
 \tag{19}
\]
To verify the last equality, the kinetic commutator acts only
on block \(c\) and commutes with \(J_b\); disjoint blocks also
give \([J_b,J_c]=0\). The remaining commutator is with the
multiplication operator \(\theta X_cV\).

A plaquette contributing to \(X_bX_cV\) must contain a link
from each block. Native lattice incidence gives exactly one
such bridging plaquette for each horizontal coarse neighbor,
and four for each vertical neighbor. There are none for other
pairs. Each bridge meets each block in exactly one link.

For a shared link,
\[
 |\nabla_e\phi_b|=\rho_\epsilon'/2\le1/2.
\]
For two distinct links of a plaquette, the mixed Hessian block
of \(W\) has operator norm at most \(1/4\): a native tangent of
length \(|v|\) inserts the SU(2) generator \(i\,v\cdot\sigma/2\),
and the normalized trace of the product is bounded by the
product of these operator norms. Inverse links only change
signs and unitary conjugations. There is no cross-link
connection term. Thus one bridge contributes
\[
 |X_bX_cW|\le(1/4)(1/2)(1/2)=1/16.
\]
Equation (19) proves (16)–(17).

Translation invariance and symmetry give (18).
Moreover, \(\sum_b\phi_b\) is nonconstant: varying a link in
one selected block changes that summand and no other.
The kernel of \(H\) is precisely the vacuum line. Therefore
\[
 t_0=\frac1N\left\|H^{3/2}
        \left(\sum_b(\phi_b-\mu\phi_b)\right)h\right\|^2>0,
 \tag{20}
\]
where \(N=L^3/4\). The same argument shows positivity for every
nonzero linear combination of the selected sources. \(\square\)

At \(L=4\), opposite horizontal coarse neighbors coincide;
two bridging faces contribute to that matrix entry. The same
formulas hold with oriented parallel-bond multiplicities.
The analytic continuum statement uses \(L\ge6\); native
geometry controls separately verify \(L=4,6,8\).

The second-moment matrix \(Q\) is generally not nearest neighbor.
What follows from operator Cauchy–Schwarz is
\[
 \boxed{Q^2\preceq u\mathcal T,\qquad
 \widehat Q(k)^2\le u\,\widehat{\mathcal T}(k).}
 \tag{21}
\]
Indeed the Gram matrix of the vectors
\(H^{1/2}(\phi_b-\mu\phi_b)h\) and
\(H^{3/2}(\phi_b-\mu\phi_b)h\) is
\[
 \begin{pmatrix}uI&Q\\Q&\mathcal T\end{pmatrix}\succeq0.
\]
Its Schur complement proves (21).
This preserves the full actual correlation matrix; it does
not replace it by independent plaquettes.

## 4. Theorem 81.3 — uniform spatial variation and one remaining coefficient

Let \(f\in C_c^\infty(\mathbb R^3)\) be real, with its support
inside the physical torus for all sufficiently large regulators.
Use the cells of sides \(2a,2a,a\), their block projection \(P_a\),
and
\[
 d_b=2\theta_a^{1/4}a^2\bar f_b,\qquad
 G_a=\sum_bd_b(\phi_b-\mu_a\phi_b),\qquad
 q_a=c\|P_af\|_2^2.
 \tag{22}
\]
Then the actual physical third moment satisfies
\[
 \boxed{
 \langle G_ah_a,\mathsf H_a^3G_ah_a\rangle
 \le \kappa_a^2t_{0,a}\,q_a+
                  \frac{c^3}{2}\|\nabla f\|_2^2.}
 \tag{23}
\]
Equivalently, writing \(\sigma_a=t_{0,a}/\theta_a\),
the first term is \(c^2\sigma_aq_a/a^2\).
The existing bound is only
\[
 0<\sigma_a\le C_{7/30}
       =339994914968/2401;
 \tag{24}
\]
no bound \(\sigma_a=O(a^2)\) is established here.

There is also a uniform momentum-variation bound. For allowed
physical torus momenta \(p\), with
\(k=(2ap_x,2ap_y,ap_z)\),
\[
 \boxed{\left|\kappa_a^2
 [\widehat{\mathcal T}_a(k)-t_{0,a}]\right|
                   \le\frac{c^2}{2}|p|^2.}
 \tag{25}
\]
This is a bound on an energy-moment symbol. It is not a particle
dispersion relation or a proof of Lorentz covariance.

**Proof.** Discrete summation by parts in (18) gives the exact identity
\[
 d^{\mathsf T}\mathcal T d
 =t_0\sum_bd_b^2-
        \sum_iK_i\sum_b(d_{b+s_i}-d_b)^2.
 \tag{26}
\]
Cell averaging commutes with coarse translations and is an
\(L^2\) contraction. With physical steps
\(\ell_x=\ell_y=2a,\ell_z=a\),
\[
 \sum_bd_b^2=\sqrt\theta\,a\|P_af\|_2^2,\qquad
 \sum_b(d_{b+s_i}-d_b)^2
 \le\sqrt\theta\,a\,\ell_i^2\|\partial_i f\|_2^2.
 \tag{27}
\]
The last inequality follows by integrating the derivative along
the translation segment. Multiplying (26) by \(\kappa_a^3\)
and using (17) and (2), each direction contributes at most
\(c^3\|\partial_i f\|_2^2/2\). This proves (23).

Finally \(1-\cos x\le x^2/2\), (17), and (18) give
\[
 |\widehat{\mathcal T}(k)-t_0|
 \le\sum_i|K_i|k_i^2\le(\theta a^2/2)|p|^2.
\]
Equation (2) proves (25).
The upper bound in (24) is YM80's third-moment estimate applied
to the constant coefficient vector. \(\square\)

### A sufficient continuum gate

For this moment route, the remaining scalar estimate is
\[
 \boxed{\sup_a\kappa_a^2t_{0,a}\le C_0<\infty.}
 \tag{28}
\]
If it were established, then for fixed nonzero \(f\),
\(q_a\to c\|f\|_2^2>0\), and (23) would bound the third moment
uniformly. For example, choose \(0<q_*\le q_a\le Q_*\) eventually,
put
\[
 S=C_0Q_*+\tfrac12c^3\|\nabla f\|_2^2,\qquad
 M=\sqrt{Q_*S}.
\]
Spectral Cauchy–Schwarz gives \(m_2(G_a)\le M\).
For fixed \(\lambda\ne0\), the centered phase
\(w_a=(e^{i\lambda G_a}-\mu_ae^{i\lambda G_a})h_a\) then has
\[
 m_0(w_a)\le1,\quad m_1(w_a)\ge\beta\lambda^2q_*,
 \quad m_1(w_a)\le E_*:=\lambda^2Q_*,
 \quad m_2(w_a)\le\lambda^4Q_*^2+\lambda^2M.
\]
The scaled-vector resolvent witness consequently gives the
uniform positive floor
\[
 R_T(w_a):=
 \left\|\frac{T\mathsf H_a}{1+T\mathsf H_a}w_a\right\|^2
 \ge
 \frac{T^2\beta^2\lambda^4q_*^2}
 {1+2TE_*+T^2(\lambda^4Q_*^2+\lambda^2M)}>0.
 \tag{29}
\]
Section 5 converts it to positive-time, finite-energy survival.

Equation (28) is sufficient, not necessary. A third moment can
diverge while a fixed positive-energy component survives:
for instance a comparison spectral measure with mass one at
\(E=1\) and mass \(n^{-2}\) at \(E=n\) has bounded first moment
and diverging third moment. The third-moment coefficient is
therefore not the physical mass gap.

Neither translation symmetry nor zero-mean smearing removes
the first term in (26). Zero mean excludes the single mode
\(k=0\), while \(t_0\) is an additive term at every momentum.
The strict finite-regulator positivity in (20) also rules out
an exact zero row-sum identity for this source. It does not
rule out the scaling required in (28).

## 5. Theorem 81.4 — a rational positive-time energy window

For any nonnegative self-adjoint \(K\), define
\[
 F_T(K)=2[(I+TK/3)^{-1}-(I+3TK)^{-1}],\qquad
 f(x)=\frac{16x}{3x^2+10x+3}.
 \tag{30}
\]
The scalar function satisfies \(0\le f\le1\),
\(f(1/x)=f(x)\), and vanishes at both \(0\) and \(\infty\).
For a vector \(v\) of finite first moment, put
\[
R_T=\left\|\frac{TK}{1+TK}v\right\|^2,\quad
 q=\|K^{1/2}v\|^2,\quad Q_T=\|F_T(K)v\|^2.
\]
If \(R_T>0\), then
\[
 \boxed{Q_T\ge\frac{256R_T^3}{(3Tq+R_T)^2}.}
 \tag{31}
\]
For \(\|v\|\le1\), the window state also obeys
\[
 \langle F_Tv,KF_Tv\rangle\le2/T,\qquad
 \langle F_Tv,K^{-1}F_Tv\rangle\le2T,
 \tag{32}
\]
where the inverse is on the positive spectral subspace.

**Proof.** Let \(b(x)=x/(1+x)\). A direct polynomial calculation gives
\[
 \frac{b(x)^3}{f(x)}
 \le\frac{3x+b(x)^2}{16},\qquad x>0.
 \tag{33}
\]
After multiplication by \(16(1+x)^3\), the nonnegative
remainder is \(7x^2+3x\).
Hölder's inequality for the spectral measure gives
\[
 R_T=\int(f^2)^{1/3}(b^3/f)^{2/3}\,d\nu
 \le Q_T^{1/3}\left(\int b^3/f\,d\nu\right)^{2/3}.
\]
At zero the ratios are defined by continuous extension;
integrability at infinity follows from the finite first moment
and (33). Applying (33) proves (31). The right side increases with \(R_T\), so
a certified lower value for \(R_T\) and an upper value for \(q\)
can be inserted.

For (32), \(xf(x)^2\le2\). The needed polynomial is
\[
 \begin{aligned}
 &(3x^2+10x+3)^2-128x^3\\
 &=9x^2(x-4)^2+4x(x-13/4)^2+(71/4)x+9\ge0.
 \end{aligned}
\]
Reciprocal symmetry gives \(f(x)^2/x\le2\). Spectral integration
proves both bounds. \(\square\)

The filter has the semigroup representation
\[
 F_T(K)=\int_0^\infty j_T(t)e^{-tK}\,dt,\quad
 j_T(t)=2[(3/T)e^{-3t/T}-(1/(3T))e^{-t/(3T)}].
 \tag{34}
\]
Here \(\int j_T=0\) and \(\|j_T\|_1\le4\). In the actual vacuum
diffusion, with \(K=\mathsf H_a\) and bounded local observable \(O\),
\(\int j_T(t)O(X_t)\,dt\) is a bounded positive-time functional,
and its OS vector is \(F_T(K)Oh\). It can be approximated by
integrals over \([\varepsilon,R]\subset(0,\infty)\): the discarded
part is bounded uniformly over paths by
\(\|O\|_\infty\int_{(0,\varepsilon)\cup(R,\infty)}|j_T(t)|dt\),
which tends to zero. On each compact interval, the finite-regulator
semigroup's strong continuity gives the corresponding cylinder
approximation in \(L^2\). Thus neither the endpoint at zero nor
the infinite time support requires an extra domain assumption.
The zero integral removes the vacuum component.
The time smearing, rather than a new time-zero local operator,
is what this construction supplies.

Because \(f^2\) vanishes at both spectral endpoints, along
compactified weak convergence of measures of mass at most one,
a positive limiting \(Q_T\) is equivalent to positive limiting
mass in \((0,\infty)\). It is a survival test; positivity at each
cutoff alone is insufficient.

### An actual polynomial lower certificate

For YM80's smooth phase write
\[
 e=\lambda^2q_a,\quad r_a=12000c/a,\quad
 u=q_G/q_a\in[\beta,1].
\]
The phase identities give
\(m_0\le1\), \(m_1=eu\), and
\(m_2\le e^2u+er_a\sqrt u\).
The resolvent witness has the correlated denominator
\[
 D_u=1+(2Te+T^2e^2)u+T^2er_a\sqrt u,\qquad
 R_T\ge T^2e^2u^2/D_u.
 \tag{35}
\]
The factor \(2T\) here is required; YM80's heat-filter denominator
had \(T\). Substitution into (31) gives
\[
 \boxed{Q_T\ge
 \frac{256T^4e^4\beta^4}
 {D_\beta(3D_\beta+Te\beta)^2}.}
 \tag{36}
\]
To justify minimizing over \(u\), each denominator factor of
the form \(c_0+Au+C\sqrt u\), with \(c_0>0\), has logarithmic derivative
\(uD'/D<1\). The numerator has power four and the denominator
has total power three, so the whole quotient is increasing.

For fixed \(T,e>0\), the exact-square-root certificate in (36)
is asymptotic to
\[
 \frac{256\beta^{5/2}e}{9T^2r_a^3}.
 \tag{37}
\]
Replacing \(\sqrt\beta\) by \(29/50\) gives a rational lower
certificate of the same order \(a^3\). This is unconditional
for the actual phase family, but still vanishes.

The relaxed OU comparison from YM80 attains the leading
coefficient in (37): its rate is \(r_a/\sqrt\beta\), its phase
Poisson mean is \(\beta^{3/2}e/r_a\), and its leading positive
spectral atom has precisely the high-energy attenuation of
\(f^2\). Thus this moment information alone does not imply
uniform window survival. This comparison is not an interacting
Yang–Mills vacuum.

### An actual finite-regulator spectral band

For the actual side-four model at \(\theta=1/48\), YM79 proved
for \(v=(W_b-\mu W_b)h\), in bare time,
\[
 R_1(v)\ge32607/250000,\quad
 \|v\|^2\le1,\quad
 \langle v,Hv\rangle=\mu(1-W_b^2)\le1.
\]
Equation (31) therefore gives
\[
 \boxed{Q_1(v)\ge
 \frac{61632531115632}{1063322424390625}
 >\frac{36}{625}.}
 \tag{38}
\]
Outside the band \([1/32,32]\), the window satisfies
\(f(x)^2\le(512/3395)^2\).
Since the total spectral mass is at most one, a single
outside-mass budget, not two separate unit budgets, suffices:
\[
 \begin{aligned}
 \nu_v([1/32,32])
 &\ge Q_1(v)-(512/3395)^2\\
 &>\frac{36}{625}-(512/3395)^2>\frac1{32}.
 \end{aligned}
 \tag{39}
\]
Thus at least \(1/32\) of unnormalized local Wilson norm lies
in this band of the actual untruncated Hamiltonian.
The corresponding physical band is
\([\kappa/32,32\kappa]\), with physical filter time \(1/\kappa\).
This establishes a finite-regulator spectral packet, not a
band fixed along a continuum trajectory or a global gap.

## 6. Reproduction and the remaining task

The packet checks native bridge incidence, SU(2) connection
normalization, Ward constants, finite Fourier identities, matrix
moment inequalities, conditional survival budgets, rational
window inequalities, and the actual side-four spectral anchor.
Comparison matrices and OU models are labelled separately.
The analytic proofs above establish the universal statements;
finite controls replay their constants and test failure cases.

The source manifest preserves 302 upstream inputs through YM80.
Only this stage's manuscript, code, tests, workflow, evidence
and README entries are added or updated.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym81_local_third_moment_vacuum_response.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym81_local_third_moment_vacuum_response.py -v
~~~

| Question | Result |
|---|---|
| How large is a vacuum-Hessian row? | Exact identity and two-sided physical \(a^{-2}\) squared-size bound |
| Are all source moment matrices local? | The smooth-source third moment is nearest neighbor; the second need not be |
| Can spatial variation be controlled uniformly? | Yes: \(c^3\|\nabla f\|_2^2/2\) in (23) |
| What remains in this sufficient moment route? | Uniform boundedness of \(\kappa_a^2t_{0,a}\) |
| Is that remaining bound proved? | No |
| Is there an actual positive-time window certificate? | Yes, order \(a^3\) at fixed energy |
| Is a local finite-regulator spectral packet certified? | Yes, actual Wilson mass \(>1/32\) in the stated band |
| Are continuum Yang–Mills existence and its positive mass gap proved? | No |

Sources within the unchanged programme:
[YM76](YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS.md) supplies
the actual weak-coupling bounds;
[YM79](YM79_LOCAL_WITNESSES_AND_SURVIVAL_CERTIFICATES.md) supplies
the smooth double commutator and actual local Wilson witness;
[YM80](YM80_BOUNDED_PHASE_ENERGY_CERTIFICATES.md) supplies the
phase chain rule, retained energy and vacuum-Hessian formulation.
The target remains the
[Jaffe–Witten problem description](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf).
