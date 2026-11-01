# MOP-1 — A common probability model for multiple actual observers

17 October 2026. Written continuation of CCI-1 at
`a395208e6ee9e919e4bab56e38b0f9ded2d7d34f`.

The proposal is to compare observations at several positions through
their spectral coordinates and probability weights. The retained
SE48/SE49 theory already constructs one positive spectral response
when two complete source pencils are positive. LAT already gives
the exact change of probability weights under an observer shift.

This note supplies a new unconditional range for those premises.
Both actual dyadic pencils are positive on an explicit growing
bank. For sufficiently large anchors the prefix grows from RBP's
retained \(N/8\) range to \(N/6\). One resulting finite spectrum
then matches the actual response and its first derivative at every
retained position. Its probability weights obey the same exact
transport law at all those positions, and a uniform bound controls
its response on the whole interval between and beyond the probes.

The finite eigenvalues are those of the acquired interpolating
model. Identifying every one with an actual zeta zero would require
additional source information. The all-order RH condition remains
open. No scientific execution or fresh source evaluation is used.

## 1. Source, spectral coordinate, and observer position

Retain
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q(t)=\mathcal X'(t)/\mathcal X(t).
\]
An actual positive-height zero \(s=\beta+i\gamma\) determines
\[
\delta=\beta-\tfrac12,\qquad
\mu=(\gamma-i\delta)^2,\qquad
(s-\tfrac12)^2=-\mu.
\tag{1}
\]
Thus RH concerns reality and positivity of this actual spectral
coordinate \(\mu\), equivalently \(\Re s=1/2\). It does not assert
that the original zeros lie on the real \(s\)-axis.

At a finite anchor \(a\ge a_0=3249/4\), keep the actual functional
\[
y_a(\mu)=\frac a{a+\mu},\qquad
\mathcal L_a[f]
=\frac1{q(a)}\sum_\mu\frac{\nu_\mu}{a+\mu}f(y_a(\mu)).
\tag{2}
\]
The retained source has \(q(a)>0,\ |y_a(\mu)|<1\), and total
absolute normalized weight at most two. No all-order positive
measure is presumed for (2).

| Object | Role |
|---|---|
| \(\zeta,\xi,\mathcal X\) | Functions defining the retained arithmetic source |
| Actual \(\mu\) in (1) | Source spectral coordinate; it may be nonreal unless RH is proved |
| \(a\), and \(t_i=a(1+x_i)\) below | Positions at which the source is observed |
| Finite \(\mu_j^{(m)}\) | Eigenvalues of a positive model acquired from a proved finite bank |
| \(y_j(t)=t/(t+\mu_j^{(m)})\) | That model's coordinate at observation position \(t\) |
| \(\pi_j(t)\) | Its normalized positive weight at that position |

A model eigenvalue is invariant under a faithful change of its
observer coordinates. Its coordinate \(y_j(t)\) and probability
weight generally change when the observation position changes.
Changing the model depth can also change the model eigenvalues.

## 2. Retained finite positivity and two complementary kernels

Write
\[
b=\sqrt a,\qquad \ell=\log(b/(2\pi)),\qquad
N_*(a)=\left\lfloor\frac{e^{-4}\sqrt a}{4\pi}\right\rfloor.
\tag{3}
\]
RBP proves, for \(b\ge64\) and \(1\le N\le N_*(a)\), simultaneous
reference comparisons
\[
\frac1{2\ell}A_N^*\preceq A_N\preceq2A_N^*,\qquad
\frac1{2\ell}B_N^*\preceq B_N\preceq2B_N^*.
\tag{4}
\]
Here \(A_N=(c_{i+j+1}),\ B_N=(c_{i+j}-c_{i+j+1})\),
\(c_j=\mathcal L_a[y^j]\), and the star uses
\(d\nu_*=dy/(\pi\sqrt{y(1-y)})\) on \((0,1)\).

Its acquired \(N\)-node probability measure
\(\nu_N=\sum_{\ell'=1}^N\omega_{\ell'}\delta_{\tau_{\ell'}}\)
has distinct \(\tau_{\ell'}\in(0,1)\), positive weights summing to
one, and matches \(c_0,\ldots,c_{2N-1}\) exactly. This measure is
used only to estimate the actual source below.

Set
\[
F_a(x)=q(a(1+x))/q(a),\qquad
x_i=2^{-i-1},\quad i=0,\ldots,m-1.
\]
Define the two actual kernels
\[
K_y(x,z)=\mathcal L_a\!\left[
\frac{y}{(1+xy)(1+zy)}\right]
=\frac{F_a(x)-F_a(z)}{z-x},
\]
\[
K_c(x,z)=\mathcal L_a\!\left[
\frac{1-y}{(1+xy)(1+zy)}\right]
=\frac{(1+x)F_a(x)-(1+z)F_a(z)}{x-z}.
\tag{5}
\]
The diagonal entries use the derivative limits.
The complementary numerator \(1-y\) is essential for the
second spectral pencil.

Use the same Newton observers as RLO/RBP:
\[
b_i(y)=\frac{y^i}{\prod_{k=0}^i(1+x_k y)},\qquad
G_m^y=(\mathcal L_a[y b_i b_j])_{i,j<m},\qquad
G_m^c=(\mathcal L_a[(1-y)b_i b_j])_{i,j<m}.
\tag{6}
\]
The real invertible Newton matrix \(D_m\), defined by
\(b_i=\sum_j(D_m)_{ij}(1+x_jy)^{-1}\), satisfies
\[
G_m^y=D_mK_y[\mathbf x]D_m^T,\qquad
G_m^c=D_mK_c[\mathbf x]D_m^T.
\tag{7}
\]
These are identities for the actual bilinear source, before any
positivity conclusion.

## MOP-T1. Both Newton channels have explicit positive margins

Let \(m\le N\le N_*(a)\), with \(b\ge64\), and put \(R=3/2\).
For the two model forms defined by replacing \(\mathcal L_a\)
with integration against \(\nu_N\),
\[
\boxed{
G_m^{y,\mathrm{mod}},G_m^{c,\mathrm{mod}}
\succeq\frac1{12\ell\,49^m}I.
}
\tag{8}
\]
Their actual-source errors satisfy
\[
\boxed{
\begin{aligned}
\|G_m^y-G_m^{y,\mathrm{mod}}\|
 &\le2304mR^{-(2N-2m+1)},\\
\|G_m^c-G_m^{c,\mathrm{mod}}\|
 &\le3840mR^{-(2N-2m+1)}.
\end{aligned}}
\tag{9}
\]
Consequently the single sufficient condition
\[
\boxed{92160\,\ell\,m\,49^mR^{-(2N-2m+1)}\le1}
\tag{10}
\]
gives
\[
\boxed{G_m^y,G_m^c\succeq\frac1{24\ell\,49^m}I.}
\tag{11}
\]

**Reference-form bound.** RBP already proves
\(\lambda_{\min}(A_m^*)\ge24\,49^{-m}\).
The same bound holds for \(B_m^*\). To see this directly, put
\(y=(1+\cos\theta)/2\). The polynomials
\[
\psi_j(y)=\sqrt2\,
\frac{\sin((j+1/2)\theta)}{\sin(\theta/2)}
\]
are orthonormal for \((1-y)d\nu_*\). Their unnormalized first
two polynomials are \(1,4y-1\), and their recurrence is
\(v_{j+1}=(4y-2)v_j-v_{j-1}\).
Their coefficient \(l^1\) norms are at most \(7^j\).
Thus the coefficient matrix \(Q\) obeys
\(Q^TB_m^*Q=I\) and
\(\|Q\|^2\le2\sum_{j<m}49^j<49^m/24\).
Inverting this congruence proves
\[
\lambda_{\min}(B_m^*)\ge24\,49^{-m}.
\tag{12}
\]

**Model comparison.** For a real vector \(v\), use RBP's exact
common denominator:
\[
\sum_{i<m}v_i b_i(y)=p_v(y)/D(y),\qquad
D(y)=\prod_{k<m}(1+x_k y).
\]
At the positive model nodes, \(1\le D\le e<3\).
If \(T\) maps \(v\) to the monomial coefficients of \(p_v\),
RBP's dyadic estimate gives \(\|Tv\|\ge\|v\|/4\).
For either weight \(y\) or \(1-y\), moment matching through
degree \(2m-1\) therefore gives
\[
v^TG_m^{\bullet,\mathrm{mod}}v
\ge\frac1{9}(Tv)^T(A_m\ \text{or}\ B_m)(Tv)
\ge\frac1{12\ell\,49^m}\|v\|^2.
\]
This proves (8). All order comparisons in this paragraph occur
on the acquired real positive model or on the already proved
polynomial forms, not at complex actual source points.

**Actual error.** RBP proves, on \(|y|=R=3/2\),
\[
\prod_{k=0}^i|1+x_ky|\ge1/16.
\]
Hence \(y b_i b_j\) has modulus at most
\(256R^{i+j+1}\), while \((1-y)b_i b_j\) has modulus at most
\(640R^{i+j}\). Cauchy's coefficient bound and matching of the
first \(2N\) moments leave only a Taylor tail. Each remaining
moment difference has modulus at most three. For the complementary
entry its error is at most
\[
3\cdot640\sum_{k\ge2N}R^{i+j-k}
=5760R^{i+j-2N}
\le3840R^{-(2N-2m+1)}.
\]
The corresponding \(y\)-entry bound is RBP's
\(2304R^{-(2N-2m+1)}\). The real symmetric row-sum bound proves
(9). Condition (10) makes both errors at most half the model
margin, proving (11). QED.

## MOP-T2. The actual multi-observer range grows to \(N/6\)

Take \(N=N_*(a)\). Condition (10), and thus both conclusions in
(11), hold throughout either schedule
\[
\boxed{
\begin{array}{ll}
N\ge128,&1\le m\le\lfloor N/8\rfloor,\\
N\ge2048,&1\le m\le\lfloor N/6\rfloor.
\end{array}}
\tag{13}
\]
The second schedule strengthens the retained dyadic prefix
coefficient by a factor \(4/3\), at its stated larger threshold.

**Proof of the first schedule.** Retain
\[
\ell<4+\log(2(N+1))<N/8\qquad(N\ge128).
\]
Since \(49<R^{10}\), the left side of (10) is less than
\[
92160\ell mR^{-(2N-12m+1)}
\le11520\ell N R^{-N/2}
<1440N^2R^{-N/2}.
\]
The last expression decreases for \(N\ge128\).
At \(N=128\), \(R^{64}>2^{32}\) and
\(1440(128)^2<2^{25}\), so it is below \(2^{-7}<1\).

**Proof of the second schedule.** A slightly sharper elementary
comparison is
\[
49<R^{39/4}.
\tag{14}
\]
Indeed \(R^8>25\), so
\(R^{39}>(2/3)25^5=(25/24)50^4>49^4\).
For \(m\le N/6\), the left side of (10) is consequently less than
\[
92160\ell mR^{-2N+(47/4)m-1}
\le15360\ell N R^{-N/24}
<1920N^2R^{-N/24}.
\tag{15}
\]
This decreases for \(N\ge2048\), since \(\log R>1/3\) and
\(2/N-\log R/24<0\). At \(N=2048\),
\(R^{N/24}>R^{84}>2^{42}\), while \(1920N^2<2^{33}\).
Thus (15) is below \(2^{-9}<1\). This proves both schedules. QED.

No arithmetic source value or matrix spectrum is evaluated in
these threshold arguments. The general sufficient gate (10)
remains available outside the two simple schedules.

## MOP-T3. Both actual source pencils acquire one common finite spectrum

At the covered positions
\[
t_i=a(1+x_i),\qquad r_i=q(t_i),\qquad
T=\operatorname{diag}(t_i),
\]
form the actual SE48 pencils
\[
L_{ij}=-q[t_i,t_j],\qquad
J_{ij}=(tq(t))[t_i,t_j],
\]
\[
L_{ii}=-q'(t_i),\qquad
J_{ii}=q(t_i)+t_iq'(t_i).
\tag{16}
\]
Here \(f[t_i,t_j]=(f(t_i)-f(t_j))/(t_i-t_j)\).
Equation (5) gives exactly
\[
L_m=\frac{q(a)}a K_y[\mathbf x],\qquad
J_m=q(a)K_c[\mathbf x].
\]
Therefore (11) implies
\[
\boxed{
D_mL_mD_m^T\succeq
\frac{q(a)}{24a\ell\,49^m}I,\qquad
D_mJ_mD_m^T\succeq
\frac{q(a)}{24\ell\,49^m}I.
}
\tag{17}
\]
Both actual forms are strictly positive.

Apply the existing SE48 finite-pencil construction:
\[
\mathfrak A_m=L_m^{-1}J_m,\qquad
q_m(z)=r^T(J_m+zL_m)^{-1}r.
\tag{18}
\]
The operator is self-adjoint for its acquired positive \(L_m\)
pairing and has \(m\) distinct positive eigenvalues. Its response is
\[
\boxed{
q_m(z)=\sum_{j=1}^m\frac{R_j}{z+\mu_j},
\qquad \mu_j>0,\quad R_j>0,
}
\tag{19}
\]
and
\[
\boxed{q_m(t_i)=q(t_i),\qquad q_m'(t_i)=q'(t_i)
\quad(0\le i<m).}
\tag{20}
\]
Here and below the finite-model label \(\mu_j\) abbreviates
\(\mu_j^{(m)}\). It does not identify the model pole with an
independently certified actual source pole.

For completeness, the realization uses the exact identities
\(J+LT=r1^T\) and \(J+TL=1r^T\).
In the symmetric chart \(L^{-1/2}JL^{-1/2}\), the eigenvalues are
positive. If \(Jv=\mu Lv\), then
\((T+\mu I)Lv=1(r^Tv)\).
The number \(r^Tv\) cannot vanish for a nonzero \(v\), and an
eigenspace of dimension greater than one would contain such a
vector. Hence the eigenvalues are simple and all residues in
(19) are positive. At a retained position,
\((J+t_iL)e_i=r\). Substitution into (18), followed by one
derivative, proves (20).

The spectral-realization theorem is retained from SE48. The new
premise supplied here is unconditional positivity of both actual
pencils on the growing dyadic schedules (13), with explicit
Newton-coordinate margins.

Under an invertible observer-coordinate change, transport \(L,J,r\)
together. The operator changes by similarity, and \(q_m\), its
\(\mu_j\), and its residues \(R_j\) remain unchanged.

## MOP-T4. One probability family matches all retained observer positions

For every real \(t>0\), define from the same model (19)
\[
y_j(t)=\frac{t}{t+\mu_j},\qquad
\pi_j(t)=\frac{R_j}{(t+\mu_j)q_m(t)}.
\tag{21}
\]
These satisfy \(0<y_j(t)<1,\ \pi_j(t)>0,\ \sum_j\pi_j(t)=1\).
The common spectral and residue invariants are
\[
\boxed{
\mu_j=t\,\frac{1-y_j(t)}{y_j(t)},\qquad
R_j=\frac{tq_m(t)\pi_j(t)}{y_j(t)}.
}
\tag{22}
\]

For any \(s,t>0\),
\[
\boxed{
\frac{q_m(s)}{q_m(t)}
=\sum_j\pi_j(t)\frac{t+\mu_j}{s+\mu_j},
}
\tag{23}
\]
\[
\boxed{
\pi_j(s)=
\frac{\pi_j(t)(t+\mu_j)/(s+\mu_j)}
{\sum_k\pi_k(t)(t+\mu_k)/(s+\mu_k)}.
}
\tag{24}
\]
Every multi-step update telescopes to the direct update. With
\(\rho=s/t\),
\[
y_j(s)=\frac{\rho y_j(t)}{1+(\rho-1)y_j(t)}.
\tag{25}
\]
These are LAT's exact positive-model transport rules applied to
the newly acquired multi-observer model.

**Proof.** Positivity and normalization follow directly from
(19). Solving the first expression in (21) and using the second
gives (22). Substitute
\(R_j=(t+\mu_j)q_m(t)\pi_j(t)\) into (19) to get (23);
normalization at \(s\) gives (24). Products of likelihood ratios
\((t+\mu_j)/(s+\mu_j)\) telescope along any chain of positions,
as do the normalizing response ratios. Direct substitution
proves (25). QED.

At the retained positions, the normalization is the actual one,
\(q_m(t_i)=q(t_i)\). Therefore the same positive probabilities give
all cross-position ratios exactly:
\[
\frac{q(t_k)}{q(t_i)}
=\sum_j\pi_j(t_i)\frac{t_i+\mu_j}{t_k+\mu_j},
\]
and also the actual first normalized derivative:
\[
\boxed{
-\frac{t_i q'(t_i)}{q(t_i)}
=\sum_j\pi_j(t_i)y_j(t_i).
}
\tag{26}
\]
The response \(q\) itself is not a probability. Equation (21)
specifies the normalization that produces probabilities in the
acquired finite model.

Thus the different observers do share one finite spectral model.
Their probabilities are related by the measured response ratios;
they are not independently assigned probabilities whose average
is declared to be a new source. At positions outside the retained
set the same family remains an exact model family; its difference
from actual \(q\) is bounded next.

## MOP-T5. The common model has a uniform actual-source response bound

Let
\[
w_i=\sqrt{t_i},\qquad w_*=\sqrt{3a/2}=w_0,\qquad
M_a=2q(3a/2)+\frac{41}{8w_*^2}+\frac1{8w_*}.
\tag{27}
\]
For every covered bank of size \(m\), throughout \(\Re w>1\),
\[
\boxed{
|q(w^2)-q_m(w^2)|
\le M_a\frac{|w+w_*|^2}{|w|}
\prod_{i=0}^{m-1}
\left|\frac{w-w_i}{w+w_i-2}\right|^2.
}
\tag{28}
\]
In particular, for every real \(t\in[a,3a/2]\),
\[
\boxed{
|q(t)-q_m(t)|\le6\sqrt a\,M_a\,36^{-m}.
}
\tag{29}
\]
This bounds the full actual source on a continuous interval,
including points not used to acquire the model.

**Actual-source analytic bound.** The retained SE49 proof gives
\[
F(w)=wq(w^2),\qquad
|F(w)|<41/8+|w|/4\quad(\Re w>1).
\tag{30}
\]
Its Euler, Gamma and completed-source terms are all included.
No RH premise enters this half-plane bound.

**Finite-model analytic bound.** Put \(F_m(w)=wq_m(w^2)\).
Every summand in (19) has positive real part on \(\Re w>0\),
because
\[
\frac{w}{w^2+\mu}
=\frac12\left(\frac1{w-i\sqrt\mu}+
             \frac1{w+i\sqrt\mu}\right).
\]
At \(w_*\), (20) gives \(F_m(w_*)=c=w_*q(3a/2)>0\).
The half-plane Schwarz inequality applied to
\((F_m-c)/(F_m+c)\) yields
\[
|F_m(w)|\le c\,\frac{1+\rho(w)}{1-\rho(w)},\qquad
\rho(w)=\left|\frac{w-w_*}{w+w_*}\right|.
\]
Since
\[
\frac{1+\rho}{1-\rho}
=\frac{(|w+w_*|+|w-w_*|)^2}{4w_*\Re w},
\]
one obtains
\[
|F_m(w)|\le q(3a/2)\frac{(|w|+w_*)^2}{\Re w}.
\]
On \(\Re w>1\), use
\(|w+w_*|^2\ge |w|^2+w_*^2\) to get
\[
\left|\frac{F_m(w)}{(w+w_*)^2}\right|\le2q(3a/2).
\]
Similarly, (30) implies
\[
\left|\frac{F(w)}{(w+w_*)^2}\right|
\le\frac{41}{8w_*^2}+\frac1{8w_*}.
\]
Therefore \((F-F_m)/(M_a(w+w_*)^2)\) is holomorphic and bounded
by one on \(\Re w>1\).

**Interpolation factors.** Equation (20) gives a double zero at
each \(w_i\). Divide the bounded analytic function successively
by the half-plane Schwarz factors
\((w-w_i)/(w+w_i-2)\), twice for each probe. This is the same
finite interpolation argument proved in SE49, now with the
declared moving normalization. Restoring the prefactor proves
(28).

For real \(w\in[b,w_*]\), every \(w_i\) is in that interval.
Since \(w_*<5b/4\) and \(b\ge64\),
\[
|w-w_i|<b/4,\qquad w+w_i-2\ge2b-2\ge3b/2.
\]
Thus each unsquared factor has modulus at most \(1/6\).
Also \((w+w_*)^2/w\le4w_*^2/b=6b\).
Equation (29) follows. QED.

For two covered prefixes \(m,n\) at the same anchor, both include
the first position \(3a/2\) and use the same \(M_a\). Their
response difference therefore obeys
\[
|q_m(t)-q_n(t)|
\le6\sqrt a\,M_a(36^{-m}+36^{-n})
\quad(t\in[a,3a/2]).
\tag{31}
\]
This is controlled agreement of their responses. It does not
identify individual finite eigenvalues across different depths.

## 3. What has advanced, and what still needs proof

| Item | Result |
|---|---|
| Actual \(y\)-kernel dyadic positivity | Retained RBP range and a wider large-anchor schedule |
| Actual complementary kernel | Positive with an explicit source error and margin |
| Both actual spectral pencils | Positive throughout (13) |
| Number of simultaneous actual response/derivative probes | Up to \(\lfloor N_*(a)/6\rfloor\) when \(N_*(a)\ge2048\) |
| Common positive finite spectrum and residues | Acquired from those actual pencils |
| Probabilities across observation positions | One exact normalized transport family |
| Actual response between probes | Uniform bound (29) |
| Identification of all model poles with actual xi zeros | Not established |
| All-order fixed-anchor positivity / RH | Open |

The eventual prefix size now has the explicit asymptotic scale
\[
m_+(a)\sim\frac{e^{-4}}{24\pi}\sqrt a.
\tag{32}
\]
RBP's retained \(N/8\) dyadic schedule had coefficient
\(e^{-4}/(32\pi)\). The \(4/3\) comparison concerns these proved
sufficient ranges, not an optimal constant or an actual count
of positive eigenvalues.

The model poles lie at \(-\mu_j<0\) in the \(t\) coordinate,
and their inverse images are \(s=1/2\pm i\sqrt{\mu_j}\).
These are finite-model spectral locations; their residues need not
be integers. The earlier SE48
pipeline uses separate actual-source zero enclosures when
certifying a particular arithmetic zero; that step is not
replaced by interpolation.

The number of probes grows while their interval \([a,3a/2]\)
moves outward. This does not establish all bank sizes at one
finite anchor. There is a precise limitation in the error factors.
Take \(m=\lfloor N_*(a)/6\rfloor\) and put
\(c=e^{-4}/(24\pi)\). For any fixed \(w\) with \(\Re w>1\),
\[
\boxed{
\lim_{a\to\infty}
\prod_{i<m}\left|\frac{w-w_i}{w+w_i-2}\right|^2
=e^{-4c(\Re w-1)}>0.
}
\tag{33}
\]
Indeed \(m/b\to c\) and summability of \(x_i=2^{-i-1}\) gives
\(\sum_{i<m}(w_i-1)^{-1}=m/(b-1)+O(b^{-1})\to c\).
For fixed \(w\), the logarithm of each squared factor is
\(-4(\Re w-1)/(w_i-1)+O_w((w_i-1)^{-3})\), uniformly as
\(b\to\infty\). The summed remainder is \(O_w(m/b^3)\to0\).
This proves (33): the factors that are small within the moving
interval need not vanish at a fixed old observation position.

Likewise, separately positive banks at different anchors do not
establish positivity of their combined bank: the mixed source
readings and both Schur innovations remain required. The exact
weight transport within a fixed finite model does not supply
those additional actual-source signs.

Consequently CCI's full negative-index count has not been shown
to vanish. BRR's all-width return requirement, PRW's unrestricted
support window, and FCP's restricted support window retain their
previous status. This continuation advances an actual finite
multi-observer range and its shared response, without replacing
the original source by a presumed positive infinite measure.

## References and verification scope

- [RBP-1](RH_RELATIONAL_BANK_POSITIVITY.md): actual polynomial positivity, finite probability carrier, dyadic Newton bounds, and the preceding prefix schedule.
- [RLO-1](RH_RELATIONAL_OBSERVER_CONFLUENCE.md): actual divided-difference kernels and exact Newton congruence.
- [SE48 retained source](papers/native-prime-observation-closure/sources/SE48.md), theorem 1: positive spectral pencil and response/derivative interpolation.
- [SE49 retained source](papers/native-moving-zero-covariance/sources/SE49.md), theorem 2: complete-source half-plane bound and Schwarz-factor interpolation.
- [POC-1](papers/native-prime-observation-closure/THEOREM.md): the same source pencil for general observer positions and the all-order compatibility requirement.
- [LAT-1](papers/native-layer-anchor-transport/THEOREM.md), theorem 2: exact probability transport, spectral poles and residue invariants.
- [CSS-1](RH_CURVATURE_SQUARE_SYMMETRY.md): actual spectral-coordinate reality and the complete curvature/first-derivative record.
- [CCI-1](RH_CONJUGATE_CANCELLATION_INDEX.md): the unresolved actual-source negative-index count.
- [BRR-1](RH_BOUNDARY_OBSERVER_RETURN.md): universal-width source return obligation.

The comparisons, threshold inequalities, spectral realization and
analytic error bounds are written proofs. No source value, finite
spectrum, probability vector or held-out error has been numerically
evaluated in this continuation. No literature-wide novelty,
optimality or machine-checked analytic proof is claimed.

Only this note and its README entry change. No frozen packet,
engine or workflow is altered. There is no scientific computation,
independent agent, replay or workflow dispatch. The documentation
commit uses `[skip ci]`, and the PR remains an unmerged draft.
