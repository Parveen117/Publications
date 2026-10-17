# RCD-1 — Residual completion and dual observer transport

17 October 2026. Written continuation of UCB-1 at
838d6d2b2aa83123ca1d92790e2a83236bdb5068.

The question is whether incomplete observer readings can be completed
at a preferred dual position, with mock-theta completion as a possible
guide. This note first determines the exact omitted source.

For a finite positive actual observer bank, its residual factors into
known interpolation zeros and a filtered actual-source response.
The filter is invariant under the native change of anchor. Every
hypothetical off-critical pole survives the filter, and the original
residual retains that pole's exact integer residue.

On UCB's declared finite range, the filtered response is strictly
positive throughout the observation interval, including the fitted
positions. There the first missing derivative has a strictly positive
enclosure: the model matches q and q' but underestimates q'' at each
fitted position. This gives an explicit excess in the next source
moment above the model probability reading.

Integrating the residual gives a source correction with a quantitative
amplitude bound and an exact additive change of potential under
renormalization. Its reflection law is proved separately from its sign.
The critical line is the fixed set of the dual reflection; reflection
does not force the residual to vanish.

No non-holomorphic mock-modular completion or global RH closure is
established. The quantitative bank and support-window ranges are
unchanged. No scientific computation is used.

## 1. Source, model, and coordinate conventions

Retain the actual source of CWR/SE39:
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q(z)=\frac{\mathcal X'(z)}{\mathcal X(z)}
=\sum_\mu\frac{\nu_\mu}{z+\mu}.
\tag{1}
\]
The series is absolutely and locally uniformly convergent off its
poles. Distinct source values arise from positive-height xi zeros,
\[
\mu=(\gamma-i(\beta-1/2))^2,
\]
with positive integer multiplicity \(\nu_\mu\). The source values
may be nonreal. They form a conjugation-symmetric discrete set,
with the retained inverse-square-height summability.

Let \(A=\{t_1,\ldots,t_m\}\), \(m\ge1\), be a set of distinct
positive source positions with both actual pencils strictly positive:
\[
L_{ij}=-q[t_i,t_j],\qquad
J_{ij}=(tq(t))[t_i,t_j],\qquad u_i=q(t_i).
\]
Their diagonal entries use the derivative limits.
SE48 supplies the unchanged model
\[
q_A(z)=u^T(J+zL)^{-1}u
      =\sum_{j=1}^m\frac{r_j}{z+\eta_j},
\qquad r_j>0,\quad \eta_j>0,
\tag{2}
\]
with distinct \(\eta_j\) and exact matching of q and q' at every
\(t_i\). Here \(\eta_j\) are model eigenvalues, distinct in meaning
from the actual \(\mu\). Set
\[
R_A(z)=q(z)-q_A(z).
\tag{3}
\]

For the quantitative statements, use UCB's common anchor a:
\[
b=\sqrt a,\qquad \ell=\log\frac b{2\pi},\qquad
N_*(a)=\left\lfloor\frac{e^{-4}b}{4\pi}\right\rfloor,
\]
\[
I_a=[a/2,3a/2],\qquad
M(a)=\lfloor N_*(a)/12\rfloor,\qquad
\delta_k(a)=\frac6{\ell(441/4)^k}.
\tag{4}
\]
Assume \(N_*(a)\ge128\), \(A\subset I_a\), and \(m+1\le M(a)\).
This is UCB's extra-position condition; it is not needed for the
purely algebraic residual identities.

Write \(x_i=t_i/a-1\) and \(S_k(x)=\sum_{j=0}^{k-1}x^{2j}\).
UCB-T6 gives, at \(t=a(1+x)\in I_a\) outside A,
\[
R_A(t)\ge q(a)(1+t/a)\delta_{m+1}(a)
\frac{\prod_{i=1}^m(x-x_i)^2}{S_{m+1}(x)}>0.
\tag{5}
\]
At fitted positions \(R_A=R_A'=0\). The actual product in the
source sums is always bilinear; no absolute-square substitution
is made at nonreal \(\mu\).

## RCD-T1. Exact factorization of the complete source residual

Define the monic polynomials and rational filter
\[
\omega_A(z)=\prod_{i=1}^m(z-t_i),\qquad
d_A(z)=\prod_{j=1}^m(z+\eta_j),
\]
\[
p_A(\mu)=\prod_{j=1}^m(\mu-\eta_j),\qquad
\Pi_A(\mu)=\prod_{i=1}^m(\mu+t_i),\qquad
h_A(\mu)=\frac{p_A(\mu)}{\Pi_A(\mu)}.
\tag{6}
\]
Then the following is an identity of meromorphic functions:
\[
\boxed{
R_A(z)=\left(\frac{\omega_A(z)}{d_A(z)}\right)^2
\mathcal T_A(z),\qquad
\mathcal T_A(z)=
\sum_\mu\frac{\nu_\mu h_A(\mu)^2}{z+\mu}.
}
\tag{7}
\]
The filtered response \(\mathcal T_A\) is holomorphic near every
positive fitted position. Formula (7) retains the entire actual
source, including any nonreal spectral values.

**Bilinear residual identity.** Put
\(f(\mu)=((t_i+\mu)^{-1})_{i=1}^m\),
\(P(z)=J+zL\), and \(v(z)=P(z)^{-1}u\).
Direct divided differences of (1) give
\[
L=\sum_\mu\nu_\mu f(\mu)f(\mu)^T,\qquad
J=\sum_\mu\nu_\mu\mu f(\mu)f(\mu)^T,\qquad
u=\sum_\mu\nu_\mu f(\mu).
\]
Consequently expansion of the square yields
\[
\boxed{
R_A(z)=\sum_\mu\frac{\nu_\mu}{z+\mu}
\left(1-(z+\mu)f(\mu)^Tv(z)\right)^2.
}
\tag{8}
\]
Indeed its three terms are q, \(-2u^Tv\), and \(v^TP(z)v\).
Since \(P(z)v=u\), they sum to (3).

**Factorization of the residual observer.** Its numerator after
multiplication by \(\det P(z)\Pi_A(\mu)\) is a polynomial in z
of degree at most m. The exact identity \(P(t_i)e_i=u\) makes it
vanish at every \(z=t_i\), for every regular \(\mu\).
It is therefore \(\omega_A(z)K(\mu)\) for a polynomial K.

Evaluate the numerator at \(\mu=-z\). The term containing
\(z+\mu\) vanishes, while
\(\Pi_A(-z)=(-1)^m\omega_A(z)\). Thus
\[
K(-z)=(-1)^m\det P(z).
\]
Since \(\det P(z)=\det L\,d_A(z)\), this identifies
\(K(\mu)=\det L\,p_A(\mu)\). Therefore
\[
\boxed{
1-(z+\mu)f(\mu)^Tv(z)
=\frac{\omega_A(z)}{d_A(z)}\,h_A(\mu).
}
\tag{9}
\]
Substitute into (8) to prove (7).

The filter tends to one at large \(|\mu|\); its denominator has
zeros only at the negative observation positions, where the
retained actual source has no support. The source summability
therefore makes all these series locally uniformly absolutely
convergent away from their stated poles. The identities extend
meromorphically from regular z. QED.

This is a multipoint rational-interpolation remainder identity,
derived here from the actual source pencils. It is not a new
general spectral theorem. Although a square occurs in (7)--(9),
\(h_A(\mu)^2\) may be complex. Global positivity is not supplied
by the notation.

## RCD-T2. The residual filter survives anchor transport and retains off-line poles

At a positive anchor a define native coordinates
\[
y=\frac a{a+\mu},\qquad
y_j(a)=\frac a{a+\eta_j},\qquad x_i(a)=t_i/a-1.
\]
The same filter has the exact native expression
\[
\boxed{
h_{A,a}(y)=
\frac{\prod_{j=1}^m(1-y/y_j(a))}
     {\prod_{i=1}^m(1+x_i(a)y)},\qquad
h_{A,a}(y_a(\mu))=h_A(\mu).
}
\tag{10}
\]
For \(a'=\rho a>0\), with
\(\phi_\rho(y)=\rho y/(1+(\rho-1)y)\),
\[
\boxed{
h_{A,a'}(\phi_\rho(y))=h_{A,a}(y).
}
\tag{11}
\]
The identity is rational and is read at regular arguments. It
keeps the physical observation set A and its model fixed.

**Proof.** Substitute \(\mu=a(1-y)/y\) in (6); the m factors
\(a/y\) cancel between numerator and denominator. Also
\[
1-\frac{\phi_\rho(y)}{\phi_\rho(y_j)}
=\frac{1-y/y_j}{1+(\rho-1)y},
\]
\[
1+x_i(a')\phi_\rho(y)
=\frac{1+x_i(a)y}{1+(\rho-1)y}.
\]
The common factors cancel between the two products, proving
(11). Composition follows from the same native Mobius law
already proved in LAT. QED.

For a regular anchor with \(q(a)\ne0\), retain the actual
source functional
\[
\mathcal L_a[g]=\frac1{q(a)}
\sum_\mu\frac{\nu_\mu}{a+\mu}
g\!\left(\frac a{a+\mu}\right)
\]
whenever the displayed sum converges. The filtered source
then satisfies
\[
\boxed{
\mathcal T_A(z)=q(a)\,
\mathcal L_a\!\left[
\frac{h_{A,a}(y)^2}{1+\lambda y}\right],
\qquad \lambda=z/a-1.
}
\tag{12}
\]
This follows term by term from
\((a+\mu)/(z+\mu)=(1+\lambda y_a(\mu))^{-1}\).
Thus changing the coordinate and the source weights together
preserves the entire filtered response.

Now suppose a hypothetical actual source value \(\mu_0\) is
nonreal. Because all zeros of \(p_A\) are positive real and all
zeros of \(\Pi_A\) are negative real,
\[
h_A(\mu_0)\ne0.
\]
The filtered response has the nonremovable pole
\[
\boxed{
\operatorname*{Res}_{z=-\mu_0}\mathcal T_A(z)
=\nu_{\mu_0}h_A(\mu_0)^2\ne0.
}
\tag{13}
\]
In the original residual its residue is exactly unchanged:
\[
\boxed{
\operatorname*{Res}_{z=-\mu_0}R_A(z)=\nu_{\mu_0}.
}
\tag{14}
\]
Indeed at \(z=-\mu_0\),
\(\omega_A(z)/d_A(z)=\Pi_A(\mu_0)/p_A(\mu_0)=h_A(\mu_0)^{-1}\),
so the factors in (7) cancel exactly in the residue.

These are statements about any hypothetical nonreal source
value; none is asserted to exist. A finite model or anchor change
cannot erase one if it exists. The filter can annihilate a
positive real source value coinciding with a model eigenvalue;
the noncancellation assertion (13) specifically concerns
nonreal values.

## RCD-T3. The fitted readings conceal a strictly positive second-derivative defect

Under the extra-position range (4)--(5), the filtered response
obeys on all of \(I_a\), including the fitted positions,
\[
\boxed{
\mathcal T_A(t)\ge
q(a)(1+t/a)\delta_{m+1}(a)
\frac{d_A(t)^2}{a^{2m}S_{m+1}(t/a-1)}>0.
}
\tag{15}
\]
At every \(t_i\in A\), the residual has a zero of exactly order
two, and
\[
\boxed{
R_A''(t_i)=
2\left(\frac{\omega_A'(t_i)}{d_A(t_i)}\right)^2
\mathcal T_A(t_i)>0.
}
\tag{16}
\]
An explicit lower bound is
\[
\boxed{
R_A''(t_i)\ge
\frac{2q(a)(1+t_i/a)\delta_{m+1}(a)}
     {a^2 S_{m+1}(x_i)}
\prod_{\substack{k=1\\k\ne i}}^m(x_i-x_k)^2.
}
\tag{17}
\]

**Proof.** Away from the fitted positions, divide (5) by
\(\omega_A(t)^2/d_A(t)^2\), using
\(\omega_A(t)=a^m\prod_i(x-x_i)\), to obtain (15).
Both the filtered series and its lower bound are continuous
at every \(t_i\), so the inequality extends there.
The model polynomial \(d_A(t_i)\) is positive, the observation
root of \(\omega_A\) is simple, and (15) is strictly positive.
Taylor expansion in (7) proves the exact order and (16).
Use \(\omega_A'(t_i)=a^{m-1}\prod_{k\ne i}(x_i-x_k)\)
in (15)--(16) to obtain (17). QED.

There is a matching finite upper enclosure. Choose any
\(t_*\in A\), put \(w_i=\sqrt{t_i}\), \(w_*=\sqrt{t_*}\), and
\[
K_A=q(t_*)+\frac{41}{8w_*^2}+\frac1{8w_*}.
\tag{18}
\]
Then UCB's complete-source Schwarz bound gives
\[
\boxed{
0<R_A''(t_i)\le
\frac{K_A(w_i+w_*)^2}{8w_i^3(w_i-1)^2}
\prod_{k\ne i}
\left|\frac{w_i-w_k}{w_i+w_k-2}\right|^2.
}
\tag{19}
\]
To verify the constant, let \(w\to w_i\) in UCB (28), divide
by \(|w-w_i|^2\), and use
\[
R_A(w^2)=2w_i^2R_A''(t_i)(w-w_i)^2
         +O((w-w_i)^3).
\]
The i-th squared Schwarz denominator tends to
\(4(w_i-1)^2\); the response prefactor tends to
\(K_A(w_i+w_*)^2/w_i\). This is exactly (19).
No second derivative is numerically evaluated.

For the actual source layers and their model versions, define
\[
c_n(t)=\frac{(-t)^n q^{(n)}(t)}{n!\,q(t)},\qquad
c_n^A(t)=\frac{(-t)^n q_A^{(n)}(t)}{n!\,q_A(t)}.
\]
At every fitted point,
\[
\boxed{
c_0=c_0^A=1,\quad c_1=c_1^A,\quad
c_2-c_2^A=\frac{t_i^2}{2q(t_i)}R_A''(t_i)>0.
}
\tag{20}
\]
In the positive model, let
\[
\pi_j(t_i)=\frac{r_j}{(t_i+\eta_j)q(t_i)},\qquad
y_j(t_i)=\frac{t_i}{t_i+\eta_j}.
\]
Their sum of weights is one. Consequently the exact reading is
\[
\boxed{
c_2(t_i)-c_1(t_i)^2
=\operatorname{Var}_{\pi(t_i)}(y_j(t_i))
+\frac{t_i^2}{2q(t_i)}R_A''(t_i).
}
\tag{21}
\]
Only the first term is the variance of the stated finite
probability model. The equality does not presume a positive
all-order probability measure for the actual xi source.

This is an actual next-layer result on the stated range:
matched values and slopes retain a strictly missing second
derivative. The word curvature here, if used, means the curvature
of this scalar response. It is not an identification with CSS's
commutator curvature or with physical spacetime curvature.

## RCD-T4. Integrating the residual gives a controlled source correction

Let
\(\Omega=\mathbb C\setminus(-\infty,0]\), choose \(t_*>0\) with
\(\mathcal X(t_*)>0\), and use the principal logarithms to define
\[
\boxed{
\mathcal X_A^{t_*}(z)=\mathcal X(t_*)\,
\exp\!\left[
\sum_{j=1}^m r_j
\log\frac{z+\eta_j}{t_*+\eta_j}
\right],\qquad z\in\Omega.
}
\tag{22}
\]
This is holomorphic and nonzero on \(\Omega\), with logarithmic
derivative \(q_A\). The residues \(r_j\) are retained; they are
not replaced by one in a characteristic determinant.

Define
\[
\boxed{
\mathcal C_A^{t_*}(z)=
\frac{\mathcal X(z)}{\mathcal X_A^{t_*}(z)}.
}
\tag{23}
\]
It is holomorphic on \(\Omega\), equals one at \(t_*\), and
\[
\boxed{
\frac{(\mathcal C_A^{t_*})'}{\mathcal C_A^{t_*}}=R_A
}
\tag{24}
\]
as a meromorphic identity. Its zeros in \(\Omega\), with their
multiplicities, are exactly those of the actual \(\mathcal X\)
there.

**Proof.** For \(z\in\Omega\), none of \(z+\eta_j\) meets the
principal logarithm's cut. The exponential in (22) is therefore
single-valued, holomorphic, and nonzero. Differentiate it to
obtain (2), and take the quotient to prove (23)--(24).
A nonvanishing holomorphic denominator leaves the numerator's
zeros and multiplicities unchanged. QED.

For the quantitative range (4), choose \(t_*\in A\).
The source and model amplitudes are positive on \(I_a\).
The real reconstruction potential
\[
\mathcal V_A^{t_*}(t)=\log\mathcal C_A^{t_*}(t)
=\int_{t_*}^t R_A(v)\,dv
\tag{25}
\]
is strictly increasing on \(I_a\): its derivative is positive
outside the finite fitted set. At each fitted point it has
\[
(\mathcal V_A^{t_*})'=0,\qquad
(\mathcal V_A^{t_*})''=0,\qquad
(\mathcal V_A^{t_*})'''=R_A''(t_i)>0.
\]
UCB's upper bound gives the explicit amplitude control
\[
\boxed{
|\mathcal V_A^{t_*}(t)|
\le 9\sqrt a\,K_A\,4^{-m}|t-t_*|,
\qquad t\in I_a.
}
\tag{26}
\]
In particular
\[
e^{-\epsilon(t)}\le
\frac{\mathcal X(t)}{\mathcal X_A^{t_*}(t)}
\le e^{\epsilon(t)},\qquad
\epsilon(t)=9\sqrt a\,K_A\,4^{-m}|t-t_*|.
\tag{27}
\]
The potential is negative to the left of \(t_*\), zero at
\(t_*\), and positive to its right.

Changing only the normalization position gives the exact law
\[
\boxed{
\mathcal C_A^{b}(z)=
\frac{\mathcal C_A^{a}(z)}{\mathcal C_A^{a}(b)},\qquad
\mathcal V_A^{b}(t)=
\mathcal V_A^{a}(t)-\mathcal V_A^{a}(b).
}
\tag{28}
\]
The first identity holds for positive normalization positions a,b
with \(\mathcal X(a),\mathcal X(b)>0\); the real potential identity
holds on an interval where the source is positive. Both follow by comparing the constants in (22).
Successive normalizations compose exactly, and the derivative
\(R_A\) is unchanged.

For separately normalized source and model responses there is
also an anchor correction. If
\(F_a(z)=q(z)/q(a)\) and
\(\widetilde F_{A,a}(z)=q_A(z)/q_A(a)\), then
\[
\boxed{
\frac{F_a(z)}{\widetilde F_{A,a}(z)}
=\frac{q(z)/q_A(z)}{q(a)/q_A(a)}.
}
\tag{29}
\]
Thus normalization at a cannot silently discard its residual.
This identity holds where the displayed denominators are nonzero.

The correction (23) uses the actual source; it is an exact
description of what the model omits, not an independent
construction of the full source from finite data. Its zero-free
property on all of \(\Omega\) would itself exclude off-critical
xi zeros and is not proved here. No global logarithm of
\(\mathcal C_A\) is assumed in the presence of such possible
zeros. The potential in (25) is defined on the declared positive
interval, where that issue does not arise.

## RCD-T5. Dual reflection fixes the critical line but does not erase the residual

The residual has real coefficients in the meromorphic sense:
\(R_A(\bar z)=\overline{R_A(z)}\).
Lift it to the original source coordinate by
\[
\mathcal D_A(s)=(s-\tfrac12)
R_A\bigl((s-\tfrac12)^2\bigr).
\tag{30}
\]
Then
\[
\boxed{
\mathcal D_A(1-s)=-\mathcal D_A(s),\qquad
\mathcal D_A(1-\bar s)=-\overline{\mathcal D_A(s)}.
}
\tag{31}
\]
At a regular point on \(\Re s=1/2\), \(\mathcal D_A(s)\) is
purely imaginary, while \(R_A(-\gamma^2)\) is real for
\(s=1/2+i\gamma\). Neither statement forces that value to vanish.

**Proof.** The first reflection changes \(s-1/2\) to its
negative and leaves its square unchanged. The second changes
it to its negative conjugate. Conjugation symmetry of q and
the real finite pencils proves (31). The fixed set of
\(s\mapsto1-\bar s\) is exactly \(\Re s=1/2\), giving the
last assertion wherever the functions are regular. QED.

An actual zero is a logarithmic-derivative pole, so the regular
value statement is not applied at that pole. Equations (13)--(14)
retain the off-line poles under the same reflection laws.
The symmetry therefore identifies a distinguished line without
proving that the actual pole support is confined to it.

The chart \(z=(s-1/2)^2\) sends the critical line to the negative
real axis. A statement that every real-axis reading is intrinsically
incomplete would not be invariant under this change of chart.
The specific finite incompleteness proved here is (16)--(21).
Exact analytic values on an interval, together with the admitted
analytic class, have the separate uniqueness supplied by the
identity theorem; finite reconstruction is a different requirement.

## 2. Theta balance and the type of completion

The ordinary theta representation already gives a precise dual
balance for the actual source. With
\(\psi(x)=\sum_{n\ge1}e^{-\pi n^2x}\), the standard theta integral
recorded in NIST DLMF 25.5.13--14 yields
\[
\xi(s)=\frac12+\frac{s(s-1)}2
\int_1^\infty
\bigl(x^{s/2}+x^{(1-s)/2}\bigr)\psi(x)\,\frac{dx}{x}.
\tag{32}
\]
For \(s=\sigma+i\gamma\), direct algebra gives
\[
\log\frac{|x^{s/2}|}{|x^{(1-s)/2}|}
=(\sigma-\tfrac12)\log x.
\tag{33}
\]
For every \(x>1\), the amplitudes balance precisely at
\(\sigma=1/2\); there their sum is
\[
x^{s/2}+x^{(1-s)/2}
=2x^{1/4}\cos\!\left(\frac\gamma2\log x\right).
\tag{34}
\]
This balance is a property of the paired transform. It does not
supply a zero-exclusion argument away from the balanced line.

Zwegers' mock-theta completion adds a suitable non-holomorphic
real-analytic correction to obtain the required modular law.
That is a different type of completion from (23). Away from
its poles,
\[
\partial_{\bar z}R_A(z)=0,\qquad
\partial_{\bar z}\mathcal C_A^{t_*}(z)=0
\quad(z\in\Omega).
\tag{35}
\]
Thus neither is directly a non-holomorphic correction with
nonzero shadow in the same complex coordinate, or after a
holomorphic change of that coordinate.

A mock-modular identification would require an explicit
upper-half-plane variable, modular group and transformation
law, and a verified non-holomorphic correction/shadow. None is
constructed or assumed here. Modular weight \(1/2\) and the
source location \(\Re s=1/2\) remain distinct quantities.

Sources for this context:
[NIST DLMF 25.5](https://dlmf.nist.gov/25.5#E13);
[Zwegers, Mock Theta Functions, 2002 thesis](https://qseries.org/fgarvan/zwegers-thesis.pdf),
introduction and chapter 1. The theta identity and completion
framework are attributed background, not new claims of this note.

## 3. Established results and remaining obligation

| Object | Result |
|---|---|
| Complete source residual | Exact bilinear square and factorization (7)--(9) |
| Native anchor change | The entire rational residual filter is invariant, (10)--(12) |
| Hypothetical off-critical source pole | Survives filtering; its residue in \(R_A\) stays exactly \(\nu_\mu\) |
| Fitted value and slope | Match exactly, as retained from SE48 |
| First missing derivative | Strict positive lower and upper bounds for \(q''-q_A''\), (17), (19), on UCB's extra-position range |
| Finite probability reading | Exact additional second-moment term, (20)--(21) |
| Completed-source amplitude | Explicit finite product with retained residues and controlled correction, (22)--(27) |
| Change of potential normalization | Adds a constant; its derivative remains unchanged, (28) |
| Critical-line reflection | Exact parity and reality properties, not residual vanishing |
| Mock-modular completion or RH | Not established |

The strongest new finite conclusion is that the apparently exact
value/slope observation retains a quantitatively nonzero next-layer
defect. Removing the known double interpolation factors exposes
a positive actual response throughout the covered interval,
while keeping every possible off-critical pole in its meromorphic
continuation.

The filtered response still contains the original unknown source
support. A proof that it is holomorphic on all of
\(\Omega=\mathbb C\setminus(-\infty,0]\) would, by (13), exclude
every off-critical zero. The local real-interval sign (15) does
not establish that global property. Likewise the zero-free
property of (23) on all of \(\Omega\) remains unproved.

UCB's finite bank-size gate, MOP's dyadic range, BRR's all-width
obligation, and the earlier full/restricted support windows are
unchanged. No conclusion about a physical mass, speed, quantum
gravity, or a new numerical constant is inferred from this
scalar response correction.

## References and verification scope

- [UCB-1](RH_OBSERVER_BANK_UNION.md), T6: actual lower envelope, explicit lower residual, and complete-source upper bound.
- [SE48](papers/native-prime-observation-closure/sources/SE48.md), theorem 1: finite pencils, characteristic spectrum, and value/derivative interpolation.
- [CWR-1](RH_CUT_WAVE_RETURN.md): actual genus-zero source sum, native source functional, and retained complex source budget.
- [LAT-1](papers/native-layer-anchor-transport/THEOREM.md): exact native anchor transport and normalized model invariants.
- [MOP-1](RH_MULTI_OBSERVER_PROBABILITY.md): common positive model probabilities and their observer dependence.
- [CSS-1](RH_CURVATURE_SQUARE_SYMMETRY.md): the separate commutator-curvature identity and the need to keep source and geometric quantities typed.
- [CCI-1](RH_CONJUGATE_CANCELLATION_INDEX.md): the still unresolved actual-source negative-index count.
- [BRR-1](RH_BOUNDARY_OBSERVER_RETURN.md): the remaining all-width source return requirement.

The remainder identity, transport calculations, derivative enclosures,
and amplitude bounds are written proofs. No source value, derivative,
spectrum, probability, or zero is numerically evaluated. No
literature-wide novelty or machine-checked analytic proof is claimed.

Only this note and its README entry change. No scientific computation,
independent agent, source replay, or workflow dispatch is performed.
Frozen sources, certificates, engines, and workflows remain unchanged.
The documentation commit uses [skip ci]. The PR is stacked on
draft #123 and remains unmerged.
