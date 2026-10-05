# YM-77: physical smearing, anchored space-time limits and continuum survival

Monty Dabas. 5 October 2026. **Python 3.12 only.**

This stage moves from YM76's microscopic one-plaquette state to observables
smeared over fixed physical regions. It proves uniform physical energy and
time-increment bounds, constructs subsequential anchored random-distribution
limits on three space coordinates and one physical-time coordinate, and
preserves their reflection positivity. It also constructs an explicit global
center-charged state and proves an all-coupling upper bound for the full
finite-volume physical gap.

**The nontrivial four-dimensional Yang–Mills theory and its positive mass gap
are not proved.** The anchored limits may be zero, and their normalization
has not been identified with renormalized Yang–Mills curvature fields.
An explicit reflected-increment norm now isolates a stronger nontriviality
test. The winding-state estimate concerns global center sectors, not a
center-neutral local-observable gap.

## 1. Actual model and the physical limit

Use the unchanged SU(2) Hamiltonian on the periodic spatial cube of side
\(L\ge4\), with lattice spacing \(a\), \(p=3L^3\) links and plaquettes:
\[
 A_\theta=L_E-\theta\sum_bW_b,\quad
 H=A_\theta-E_0,\quad \mathsf H_a=\kappa_a H,\quad
 \mu_a=h_a^2\,d\Phi.
 \tag{1}
\]
The actual ground state \(h_a\) is unique, positive, smooth and normalized.
The physical calibration is fixed throughout:
\[
 \theta_a=\frac{c_B}{c_Eg_a^4},\qquad
 \kappa_a=\frac{c_Eg_a^2}{a},\qquad
 c=\sqrt{c_Ec_B},\qquad \kappa_a\sqrt{\theta_a}=\frac ca.
 \tag{2}
\]
The coupling is not retuned by defining \(a\) from a measured gap.

For the spatial \(\mathbb R^3\) limit take a sequence
\[
 a\longrightarrow0,\qquad aL\longrightarrow\infty.
 \tag{3}
\]
A fixed value of \(aL\) instead gives a spatial torus. Weak bare coupling
\(\theta_a\to\infty\) may be imposed, but Theorems 77.1–77.3 hold for every
positive coupling on any sequence satisfying (2)–(3).

Write \(r_b=\arccos W_b\). YM76 established
\[
 \Gamma(r_b)=1,\qquad
 \Gamma\!\left(\sum_b d_b r_b\right)\le4\sum_b d_b^2,
 \tag{4}
\]
and the ground-state form identity
\[
 \langle Fh_a,\mathsf H_aFh_a\rangle
       =\mu_a\Gamma_{\rm phys}(F),\qquad
 \Gamma_{\rm phys}=\kappa_a\Gamma
 \tag{5}
\]
for real centered \(F\) in the form domain.
The radial gradient identities and bounds hold almost everywhere;
the pole approximation used below preserves their uniform bounds.

## 2. Theorem 77.1 — a local smear with bounded physical energy

Label the three plaquette orientations by \(\omega\).
Let \(f=(f_\omega)\) be a real smooth compactly supported test function on
\(\mathbb R^3\), with
\(\|f\|_2^2=\sum_\omega\int|f_\omega(x)|^2dx\).
For a plaquette based at a lattice site \(x_b\), use its associated spatial
cell \(Q_{a,b}\) and the exact average
\[
 \bar f_{a,b}=a^{-3}\int_{Q_{a,b}}f_\omega(x)\,dx.
 \tag{6}
\]
There is one copy of the cell partition for each orientation.
For sufficiently large boxes the support is away from the periodic boundary.

Define the centered local multiplication observable
\[
 \boxed{
 \Phi_a(f)=\theta_a^{1/4}a^2
       \sum_b\bar f_{a,b}\,[r_b-\mu_a r_b].
 }
 \tag{7}
\]
Then the following bound holds pointwise almost everywhere on the configuration space:
\[
 \boxed{
 \Gamma_{\rm phys}(\Phi_a(f))
 \le4c\,a^3\sum_b\bar f_{a,b}^{\,2}
 \le4c\|f\|_2^2=:C_f.
 }
 \tag{8}
\]
In particular
\[
 \langle\Phi_a(f)h_a,\mathsf H_a\Phi_a(f)h_a\rangle\le C_f.
 \tag{9}
\]

**Proof.** Apply (4) with coefficients
\(\theta_a^{1/4}a^2\bar f_{a,b}\), multiply by \(\kappa_a\), and
use (2). Cellwise Cauchy–Schwarz gives
\(a^3\bar f_{a,b}^{\,2}\le\int_{Q_{a,b}}f_\omega^2\).
Summation gives (8); (5) proves (9). \(\square\)

No lower bound on \(\|\Phi_a(f)h_a\|\) was used. This is an
**unnormalized** physical energy bound, not a normalized Rayleigh bound.
The factor in (7) is selected by (8); it is not yet a justified
renormalization of a curvature polynomial.

For \(\theta_a\ge25\), YM76 supplies only the one-time estimate
\[
 \operatorname{Var}_{\mu_a}\Phi_a(f)
 \le\frac{9\pi^2}{8}a^4
       \left(\sum_b|\bar f_{a,b}|\right)^2
 \le\frac{9\pi^2}{8a^2}\|f\|_1^2.
 \tag{10}
\]
Here \(\|f\|_1=\sum_\omega\int|f_\omega(x)|\,dx\).
Thus the variance can neither be assumed bounded nor bounded away from zero.

There is a useful comparison normalization. The ordinary Riemann smear
\(\widetilde\Phi_a=a\Phi_a\) has uniformly bounded variance by (10) and
energy at most \(a^2C_f\). Its physical-time increments tend to zero in
\(L^2\). Any subsequential law for these smears is time-independent.
If a uniform positive physical gap \(m\) were independently known, its
centered variance would be at most \(a^2C_f/m\) and would vanish.

## 3. Theorem 77.2 — actual anchored space-time compactness

At each finite regulator the ground-state transform gives the stationary
reversible diffusion \(X_a(t)\), \(t\in\mathbb R\), with invariant measure
\(\mu_a\). Its generator is
\[
 \mathcal G_aF
 =-\kappa_aL_EF+2\kappa_a\Gamma(\log h_a,F).
 \tag{11}
\]
Vacuum heat correlations of multiplication observables equal correlations
of this diffusion. Define the anchored process
\[
 Y_a(t,f)=\Phi_a(f)(X_a(t))-\Phi_a(f)(X_a(0)).
 \tag{12}
\]
There is no spectral projection in this definition.

For all real \(\lambda\), \(s,t\), and positive integers \(q\),
\[
 \boxed{
 \mathbb E\exp\{\lambda[Y_a(t,f)-Y_a(s,f)]\}
 \le\exp\{C_f\lambda^2|t-s|\},
 }
 \tag{13}
\]
\[
 \boxed{
 \mathbb E|Y_a(t,f)-Y_a(s,f)|^{2q}
 \le2q!\,[4C_f|t-s|]^q.
 }
 \tag{14}
\]
Here \(2q!\) means \(2\,(q!)\), not \((2q)!\).
Neither inequality assumes a mass gap or a one-time variance bound.

**Proof of the increment bound.** Write the forward martingale for a
finite-regulator observable \(F=\Phi_a(f)\) as
\[
 M_t-M_s=F(X_t)-F(X_s)-\int_s^t\mathcal G_aF(X_u)\,du.
 \tag{15}
\]
Its quadratic variation on \([s,t]\) is at most
\(2C_f(t-s)\), by (8). On a finite time window ending at \(T\),
reversibility gives the corresponding backward martingale
\[
 \widehat M_{T-s}-\widehat M_{T-t}
 =F(X_s)-F(X_t)-\int_s^t\mathcal G_aF(X_u)\,du.
 \tag{16}
\]
Subtracting these expressions cancels the drift exactly:
\[
 F(X_t)-F(X_s)
 =\tfrac12[(M_t-M_s)
          -(\widehat M_{T-s}-\widehat M_{T-t})].
 \tag{17}
\]
Each martingale increment has exponential moment at most
\(\exp(C_f\lambda^2(t-s))\), by the exponential-martingale inequality.
Cauchy–Schwarz applied to (17) proves (13); no independence of the
two martingales is assumed. Optimizing the exponential Markov bound gives
\[
 \mathbb P(|F(X_t)-F(X_s)|\ge u)
 \le2\exp[-u^2/(4C_f|t-s|)].
 \tag{18}
\]
Integrating this tail gives (14).

The nonsmooth radial poles cause no unproved Itô step. Smooth \(r\) to
\(r_\epsilon\), flattened near \(0,\pi\), with
\(|r_\epsilon'|\le1\). The holonomy poles have codimension three,
and \(L_Er=-2\cot r\) is locally \(L^2\).
These approximants converge in the generator graph domain at every
fixed regulator; the squared \(L^2\) error of the singular second
derivative tends to zero with the cutoff radius.
The finite-regulator ground-state drift is smooth and bounded.
Martingale identities pass in \(L^2\); the same pointwise bound (8)
holds throughout the approximation. No uniform bound on that drift
across regulators is needed. \(\square\)

**Spatial distribution statement.** Realize \(\Phi_a\) as the
orientation-valued piecewise constant random distribution with value
\(\theta_a^{1/4}a^{-1}(r_b-\mu_a r_b)\) on its cells.
For every \(s>3/2\), the laws of \(Y_a\) are tight in
\[
 C_{\rm loc}\bigl(\mathbb R;H^{-s}_{\rm loc}(\mathbb R^3;\mathbb R^3)\bigr).
 \tag{19}
\]
Consequently there is a common subsequence with a random-distribution
limit \(Y\), on three space coordinates and physical time, with \(Y(0)=0\).
It has locally Hölder continuous time paths of every exponent
\(\alpha<1/2\), in the indicated weaker local Sobolev spaces.

Here are the compactness details. Choose \(3/2<s_0<s\), localize by
a smooth compactly supported \(\chi\), and use a real orthonormal
trigonometric basis \(e_k\) on a fixed torus containing its support.
The \(L^2\) norms of \(\chi e_k\) are uniformly bounded. Minkowski's
inequality and (14) give, with a finite constant \(B_{\chi,s_0}\),
\[
 \mathbb E\|\chi(Y_a(t)-Y_a(s))\|_{H^{-s_0}}^{2q}
 \le2q!\,[B_{\chi,s_0}|t-s|]^q,\qquad
 \sum_k(1+|k|^2)^{-s_0}<\infty.
 \tag{20}
\]
Kolmogorov's bounds control any time Hölder exponent below
\(1/2-1/(2q)\). Compact embedding from localized \(H^{-s_0}\)
into \(H^{-s}\), combined with time equicontinuity, gives tightness
in continuous-path topology. A countable exhaustion of space and
time gives one common subsequence. Taking arbitrarily large \(q\)
gives the stated Hölder regularity.
A further diagonal extraction over rational \(s_j\downarrow3/2\)
gives a common subsequence for every \(s>3/2\).
Using only a bound in \(H^{-s}\), without the stronger compactly
embedded space, would not suffice for this argument.

The limit has stationary time increments, inherited from stationarity
of \(X_a\). It is not asserted to be a stationary field in time.
Spatial translation invariance follows by approximating each translation
by a lattice translation and using the \(L^2\) continuity in (8), (14)
of translated test functions. The lattice cubic symmetries also pass.
Full four-dimensional Euclidean covariance does not follow.

## 4. Theorem 77.3 — reflection positivity and a nontriviality certificate

For every bounded continuous positive-time cylinder functional \(F\)
of \(Y_a\), with reflection \((\Theta Y)(t)=Y(-t)\),
\[
 \mathbb E[\overline{F(\Theta Y_a)}F(Y_a)]
 =\mathbb E\left|\mathbb E[F(Y_a)\mid X_a(0)]\right|^2\ge0.
 \tag{21}
\]
Indeed, conditional on \(X_a(0)\), the past and future are independent;
reversibility identifies their conditional laws. The anchor is fixed by
reflection and may be included in both cylinder functionals.
Weak path convergence passes (21) to the limit. The uniform higher
moments in (14) also pass polynomial cylinder reflection inequalities
by uniform integrability.

There is also a genuine temporal Hilbert-space construction for this limit.
Stationary increments make the reanchoring action
\[
 (\tau_sY)(t)=Y(t+s)-Y(s)
 \tag{21a}
\]
measure preserving. It obeys \(\Theta\tau_s\Theta=\tau_{-s}\).
For \(s\ge0\), \(U_sF=F\circ\tau_s\) preserves the positive-time
cylinder algebra. On that algebra define
\(B(F,G)=\mathbb E[\overline{F(\Theta Y)}G(Y)]\) and \(q(F)=B(F,F)\).
Reflection positivity gives the quotient Hilbert space after completion.
The identities
\[
 q(U_sF)=B(F,U_{2s}F),\qquad
 q(U_sF)^2\le q(F)\,q(U_{2s}F),\qquad
 q(U_uF)\le\|F\|_{L^2}^2
 \tag{21b}
\]
imply \(q(U_sF)\le q(F)\), by iterating and taking \(2^{-k}\)-th
powers. Null vectors therefore remain null. The induced semigroup is
symmetric and contractive; path continuity and bounded continuous
cylinders give strong continuity on a dense set and hence on the
completion. It has a self-adjoint nonnegative generator \(\mathsf H_\infty\)
and the constant vacuum vector.

This construction uses reanchored shifts, not ordinary time translations
of \(Y(t)\). It may produce only the vacuum; vacuum uniqueness is not
asserted. It does not automatically represent arbitrary cylinder
multiplication by physical operators, or construct covariant local
Yang–Mills fields. The shifted increment
\(Y(T+s,f)-Y(s,f)\) does represent the heat evolution of the state
in (22), with its finite-regulator correlations passing to the same limit.

Reflection positivity and a nonzero stochastic increment alone do not
ensure a nonvacuum quantum state. For example, two-sided Brownian motion
anchored at zero has independent halves and a rank-one reflection form.
Its ordinary increment variance is nonzero, but its reflected
increment pairing is zero.

For \(T>0\) define an actual physical Hilbert-space state
\[
 v_a(f)=\Phi_a(f)h_a,\qquad
 u_a(T,f)=(I-e^{-T\mathsf H_a})v_a(f).
 \tag{22}
\]
Its squared norm is exactly the reflected increment:
\[
 \boxed{
 D_a(T,f):=\mathbb E[Y_a(-T,f)Y_a(T,f)]
 =\|u_a(T,f)\|^2.
 }
 \tag{23}
\]
Spectral calculus and (9) give unconditional estimates
\[
 \boxed{
 \|u_a(T,f)\|^2\le TC_f,\qquad
 \langle u_a(T,f),\mathsf H_a u_a(T,f)\rangle\le C_f.
 }
 \tag{24}
\]
These states have both a uniformly bounded norm and a uniformly bounded
physical energy. They are not claimed to have norm bounded below.
The filter is the state represented by the original Euclidean increment
cylinder; it is not asserted to be a local commuting time-zero operator.

Let \(\sigma_a\) be the unnormalized spectral measure of \(u_a(T,f)\)
for \(\mathsf H_a\). For \(\epsilon,M>0\),
\[
 \boxed{
 \sigma_a([0,\epsilon])\le T^2\epsilon C_f,\qquad
 \sigma_a([M,\infty))\le C_f/M.
 }
 \tag{25}
\]
To verify this, use
\((1-e^{-TE})^2\le TE\), the sharper bound
\((1-e^{-TE})^2\le T^2E^2\) for the low-energy tail,
and \((1-e^{-TE})^2\le1\) for the energy numerator.
Uniform high-energy tails give spectral tightness; the low-energy
estimate implies that every subsequential measure has no atom at zero.
Its support may nevertheless accumulate at zero, and its total mass
may vanish.

There is a useful inverse-energy bound as well:
\[
 \boxed{\int_{(0,\infty)}E^{-1}\,d\sigma_a(E)\le T^2C_f.}
 \tag{25a}
\]
It follows from \((1-e^{-TE})^2/E\le T^2E\) in the source
spectral measure. The inverse is taken as zero on the vacuum, which
the filter annihilates. Cauchy–Schwarz consequently gives
\(\langle u_a,\mathsf H_a u_a\rangle\ge
\|u_a\|^4/(T^2C_f)\) when \(C_f>0\).
For the zero test function all states and bounds are zero; the
survival criterion below always assumes \(C_f>0\).

Now a single quantitative lower certificate would suffice for a
nonzero finite-energy observed state. If along the chosen trajectory
\[
 \boxed{D_a(T,f)\ge d>0,}
 \tag{26}
\]
then the interval
\[
 \boxed{
 \left[\frac{d}{4T^2C_f},\frac{4C_f}{d}\right]
 }
 \tag{27}
\]
contains at least \(d/2\) unnormalized spectral mass and at least half
the mass of the normalized state. The physical Rayleigh quotient is
between \(d/(T^2C_f)\) and \(C_f/d\). Formula (25) proves the band claim by subtracting at
most \(d/4\) at each end. Since \(d\le TC_f\), the band endpoints
are consistently ordered.

The common path subsequence also passes (23) to the limit by fourth-moment
uniform integrability. Thus (26) supplies a nonzero reflected state of
the limit, stronger than nonzero ordinary stochastic variance.
It has **not** been established for the actual family (7).
Nor would it alone bound the spectrum on the entire local field algebra.

Two exact spectral controls explain why (26) is independent. With
energy \(E=n\) and source mass \(1/n\), the source energy equals one,
but its filtered norm
\[
 \frac{(1-e^{-Tn})^2}{n}\longrightarrow0.
 \tag{28}
\]
With energy \(E=1/n\) and source mass \(n\), the source energy is again
one and its ordinary time-increment variance tends to \(2T\), while
\[
 n(1-e^{-T/n})^2\longrightarrow0.
 \tag{29}
\]
These are logical spectral controls, not Yang–Mills vacua.
The first loses norm at high energy; the second loses reflected norm
at zero energy. Energy boundedness alone cannot rule out either loss.

## 5. Theorem 77.4 — seven center channels give an actual global state

At a fixed lattice vertex take the three rooted axial winding holonomies
\(U_1,U_2,U_3\). Their cycles have \(L\) links each and disjoint link
sets. For every nonempty \(S\subset\{1,2,3\}\), let
\[
 W_S=\tfrac12\operatorname{Tr}\!\left(\prod_{i\in S}^{i\ {\rm increasing}}U_i\right).
 \tag{30}
\]
The product is an ordered SU(2) product, not a product of traces.
Each \(W_S\) is smooth and gauge invariant.

We first need a pointwise quaternion statement:
\[
 \boxed{\sum_{\varnothing\ne S\subset\{1,2,3\}} W_S^2\ge\frac15}
 \qquad\text{for every }(U_1,U_2,U_3)\in{\rm SU}(2)^3.
 \tag{31}
\]

**Proof of (31).** Write \(U_i=a_i+v_i\), where \(a_i\in\mathbb R\)
and \(v_i\in\mathbb R^3\). Put
\[
 x_{ij}=\operatorname{Re}(U_iU_j),\quad
 z=\operatorname{Re}(U_1U_2U_3),\quad
 A=\sum_i a_i^2,\quad B=\sum_{i<j}x_{ij}^2,\quad
 C=z^2,\quad S=A+B+C.
 \tag{32}
\]
Let \(G=(v_i\cdot v_j)\), and let \(X\) have off-diagonal entries
\(x_{ij}\) and zero diagonal. Then
\[
 G=I+aa^T-2\operatorname{diag}(a_i^2)-X,\qquad
 \|X\|\le\sqrt{2B},\qquad
 \|aa^T-2\operatorname{diag}(a_i^2)\|\le A.
 \tag{33}
\]
For the last estimate, its upper bound follows from \(aa^T\le AI\).
Its lower bound follows from the exact identity, for real \(y\),
\[
 A\|y\|^2+y^T[aa^T-2\operatorname{diag}(a_i^2)]y
 =\sum_{i<j}(a_i y_j+a_j y_i)^2.
 \tag{34}
\]

Suppose \(0<S\le1/5\). Then
\(A+\sqrt{2B}\le\sqrt{2S}<2/3\), since
\(A+\sqrt{2(S-A)}\) is decreasing in \(A\) on this range.
Thus \(G\ge I/3\), while
\(\operatorname{tr}G=3-A\ge14/5\). Writing its three eigenvalues as
\(1/3+y_i\), \(y_i\ge0\), gives
\[
 \det G\ge\frac{\operatorname{tr}G-2/3}{9}
          \ge\frac{32}{135}>\frac15.
 \tag{35}
\]
But for \(T=v_1\cdot(v_2\times v_3)\), quaternion multiplication gives
\[
 T=-z+a_1x_{23}+a_2x_{13}+a_3x_{12}-2a_1a_2a_3,
 \quad |T|\le\sqrt C+\sqrt{AB}+A\sqrt A.
 \tag{36}
\]
Let \(x=A+B\). For \(S\le1/4\),
\[
 \sqrt{S-x}+\frac x2+x\sqrt S
 \le\sqrt S+
 x\left(\frac12+\sqrt S-\frac1{2\sqrt S}\right)
 \le\sqrt S.
 \tag{37}
\]
Here the first step uses the tangent bound for the square root.
Therefore \(T^2\le S\le1/5\), contradicting \(T^2=\det G\).
If \(S=0\), (33) gives \(G=I\), while (36) gives \(T=0\),
which is already a contradiction. This proves (31). \(\square\)

For each spatial direction, multiply by the central element \(-I\)
all links crossing its periodic seam. Every plaquette intersects such
a seam in zero or two links, so \(V\) is invariant; the electric
operator and Haar measure are invariant as well.
Uniqueness of the positive vacuum makes \(h_a\) invariant under all
three center transformations. The seven \(W_S\) carry the seven
distinct nontrivial characters of this \(\mathbb Z_2^3\) symmetry.
Consequently
\[
 \mu_aW_S=0,\qquad
 \mu_a(W_SW_{S'})=0,\quad
 \langle W_Sh_a,H W_{S'}h_a\rangle=0\quad(S\ne S').
 \tag{38}
\]
Every link in the concatenated word for \(W_S\) occurs once; hence
\[
 \Gamma(W_S)=\frac{|S|L}{4}(1-W_S^2).
 \tag{39}
\]
Combining all seven characters gives an explicit actual centered state
\[
 F=\sum_{S\ne\varnothing}W_S,\qquad
 \boxed{\frac15\le\|Fh_a\|^2\le7,\qquad
         \langle Fh_a,H Fh_a\rangle\le3L.}
 \tag{40}
\]
The total length of the seven words is \(12L\), which explains the
constant \(3L\). There is no vacuum-variance assumption.

The variational principle therefore proves, at every finite coupling,
\[
 \boxed{\Delta_{\rm full,bare}\le15L,\qquad
        \Delta_{\rm full,phys}\le15\kappa_a L.}
 \tag{41}
\]
“Full” means the complete gauge-invariant finite-volume Hilbert space,
including nontrivial global center/electric-flux sectors.
This upper bound does not apply by the same argument to the
center-neutral sector or to a local-observable continuum gap.

## 6. Theorem 77.5 — an actual trajectory with closing global gap

Set \(c_E=c_B=1\) and, for integers \(n\ge2\), choose
\[
 a_n=\frac1n,\qquad L_n=n^2,\qquad g_n^2=n^{-4}.
 \tag{42}
\]
Then
\[
 \theta_n=n^8,\quad \kappa_n=n^{-3},\quad
 a_nL_n=n\longrightarrow\infty,\qquad
 \boxed{0<\Delta_{\rm full,phys,n}\le\frac{15}{n}\to0.}
 \tag{43}
\]
Each fixed regulator still has a positive gap. The upper bound concerns
their limit, without computing a large lattice or replacing it by a toy model.

For the normalized state
\(f_n=Fh_n/\|Fh_n\|\), (40) gives the actual estimate
\[
 0\le1-\langle f_n,e^{-t\mathsf H_n}f_n\rangle
 \le\frac{15t}{n}.
 \tag{44}
\]
Its physical spectral probability measures converge to a zero-energy
point mass. On the same models, YM76's normalized one-plaquette state
still carries at least \(1/3600\) of its weight in the physical band
\[
 [n/45,80n].
 \tag{45}
\]
Thus globally slow and locally ultraviolet probes coexist in the
actual Hamiltonian. Weak bare coupling and plaquette concentration
alone cannot imply a uniform full physical lattice gap.

The trajectory (42) has not been identified with a renormalized
Yang–Mills trajectory satisfying its short-distance requirements.
Moreover these winding probes span the entire periodic box, whose
physical size tends to infinity. They can escape a limit based only
on bounded local regions. Equation (43) therefore does **not** refute
a positive gap in a correctly constructed local Yang–Mills theory.

## 7. What has advanced and what the complete continuum proof still needs

The following distinction is now concrete:

| Object | Result supplied here | Remaining condition |
|---|---|---|
| Physical local smears | Uniform pointwise energy bound (8) | Identify the correct renormalized local fields |
| Anchored space-time distributions | Common subsequential limits, time regularity, spatial translations | Prove a suitable limit is nonzero and restore the required full covariance |
| Reflection form | Positivity survives (21) | A nonzero local reflected norm such as (26) |
| Temporal Hilbert space | Strongly continuous positive-energy evolution under reanchoring | Nontriviality, physical local operators and full spacetime covariance |
| Filtered physical states | Uniform norm, energy and two-sided spectral-tail bounds | Positive survival, then control of the complete relevant state space |
| Global center channels | Explicit nonzero state and gap upper bound (41) | They do not replace the local center-neutral continuum sector |
| Positive continuum mass gap | Not established | Uniform positive lower spectral edge on a complete local field core |

The Jaffe–Witten problem requires a nontrivial local quantum Yang–Mills
theory with the specified axiomatic and short-distance properties and a
finite positive mass gap. Its Section 6.5, footnote 2, explicitly warns
that compactness alone is insufficient. The present anchored construction
establishes several properties of a limit, but not those full requirements.
In particular, reflection positivity alone does not supply Euclidean
covariance, clustering, hierarchy-wide reconstruction bounds, or the
identification of Yang–Mills fields.

The next local calculation is sharply specified: certify a positive
lower bound for (23), for a physical test function and time, on a
justified trajectory; then establish a positive lower spectral edge
on the full local field core. The global center-sector upper bound
cannot substitute for either step.

## 8. Reproducible evidence

The analytic arguments above are written proofs, not finite sampling
arguments. The accompanying exact controls check quaternion identities,
native winding words and center characters, cell-average energy budgets,
increment moment constants, spectral filters and the scaling in (42).
Finite spectral counterexamples are labeled as logical controls; no
interacting vacuum is replaced by them. Large-\(n\) trajectory bounds
evaluate proved formulas and do not diagonalize those large lattices.

The source manifest pins 274 unchanged upstream files and four local
inputs. The result and SHA-256 pin use **Python 3.12 only**.
The targeted suite contains 29 new tests and 61 unchanged YM76/YM75 tests.
The dedicated workflow also replays the unchanged YM76 and YM75 evidence.
No proof-assistant formalization or completed continuum mass-gap proof
is claimed.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym77_physical_smearing_continuum_survival.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym77_physical_smearing_continuum_survival.py -v
~~~

## 9. Sources

- [YM76: local defects and weak-coupling dynamics](YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS.md):
  actual angle domains, the pointwise collective gradient bound and the
  microscopic spectral band.
- [YM75: source-resolved transfer](YM75_RH_SOURCE_RESOLVED_SPECTRAL_TRANSFER.md)
  and [YM74: all-coupling enclosures](YM74_ALL_COUPLING_ELECTRIC_ENCLOSURES.md):
  native conventions, complete electric cuts and exact spectral comparisons.
- Jaffe and Witten, [Quantum Yang–Mills Theory](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf),
  Sections 3–4 and Section 6.5, footnote 2: the continuum target and
  the distinction between compactness and a full construction.
- The forward/backward martingale argument is the established reversible
  diffusion method, often called the Lyons–Zheng decomposition. Its needed
  finite-regulator identity and constants are derived in (15)–(18).

“New” describes this programme's development. No external priority claim
is made for variational principles, reflection positivity, martingale
compactness, or general spectral filtering.
