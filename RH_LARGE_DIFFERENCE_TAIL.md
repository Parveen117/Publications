# RH continuation: positivity for large difference order

17 October 2026. Continues [the positive layer sectors](RH_POSITIVE_LAYER_SECTORS.md)
at Publications commit `e93ace29b37edb53605e25125b93b2ab9588a6bb`.

This note proves a second unbounded region: for each fixed derivative
index \(r\), every sufficiently large difference order \(k\) is positive,
with an explicit threshold. Combining it with the retained finite rectangle
settles all \(k\) in the first ten million binomial rows. The remaining
two-index condition is still unbounded, so RH is not concluded.

## Retained source and inputs

Use the unchanged actual source
\[
a=3249/4,\quad \mathcal X(w^2)=\xi(1/2+w),\quad
q=\mathcal X'/\mathcal X,\quad
c_n=\frac{(-a)^nq^{(n)}(a)}{n!\,q(a)},\quad
D_{r,k}=\sum_{j=0}^{k}(-1)^j\binom{k}{j}c_{r+j}.
\tag{1}
\]
We retain \(800<a<900\), \(0<q(a)<1/2\), and the finite-height theorem
with \(T=3\cdot10^{12}\) cited and scoped in the preceding note.

For each nontrivial zero \(\rho=\beta+i\gamma\) with \(\gamma>0\),
including multiplicity, put
\[
\delta=\beta-\tfrac12,\quad
\mu=(\gamma-i\delta)^2,\quad y=\frac{a}{a+\mu}.
\]
The preceding S1 proves the absolutely convergent actual-source identity
\[
q(a)D_{r,k}
=\frac1a\sum_{\gamma>0}y^{r+1}(1-y)^k
=\sum_{\gamma>0}\frac{y^r(1-y)^k}{a+\mu}.
\tag{2}
\]
Every zero, including any hypothetical off-line zero above T, remains
in this sum. The same source argument gives
\[
\sum_{\gamma>T}\frac1{|a+\mu|}
\le2q(a)<1,\qquad |y|<1\quad(\gamma>T).
\tag{3}
\]
At \(\gamma\le T\), all terms in (2) are positive real numbers.

A new classical input is used only to ensure at least one zero in a
distant height interval. We use the coarse consequence
\[
\left|N(x)-\frac{x}{2\pi}\log\frac{x}{2\pi e}\right|
\le\log x+10\qquad(x\ge e).
\tag{4}
\]
Here \(N(x)\) counts all zeros with \(0<\gamma\le x\), with multiplicity,
without a restriction on their real parts. Rosser's bound implies (4):
in the normalization without the \(7/8\) constant, its error is at most
\(0.137\log x+0.443\log\log x+2.463\), and
\(\log\log x\le\log x\) for \(x\ge e\).

For the accessible source and the normalization of \(N\), see
Hasanalizade, Shen and Wong, *Counting zeros of the Riemann zeta function*,
**Journal of Number Theory** 235 (2022), 219–241, equation (1.1) and
the Rosser row of Table 1:
[author-hosted published paper](https://www-math.nsysu.edu.tw/~pjwong/stuff/CountingRiemannZeros.pdf),
[DOI 10.1016/j.jnt.2021.06.032](https://doi.org/10.1016/j.jnt.2021.06.032).
The historical source is Rosser, *Explicit bounds for some functions of
prime numbers*, **American Journal of Mathematics** 63 (1941), 211–232,
[DOI 10.2307/2371291](https://doi.org/10.2307/2371291).
We use the classical Rosser estimate as reproduced in that table,
not a sharper subconvexity-based constant.

## Theorem V1. A zero exists in each required distant interval

For every real \(X\ge1000\), at least one nontrivial zero has
\(X<\gamma\le2X\). No critical-line location is asserted.

**Proof.** Let \(F(x)=\frac{x}{2\pi}\log\frac{x}{2\pi e}\).
For \(x\ge1000\), \(F'(x)=\frac1{2\pi}\log\frac{x}{2\pi}>1/8\).
Equation (4) gives
\[
N(2X)-N(X)>
\frac X8-2\log(2X)-20.
\]
At \(X=1000\), this is positive since \(\log2000<8\).
Its derivative \(1/8-2/X\) is positive for \(X\ge1000\). Hence
\(N(2X)>N(X)\) throughout that range. QED.

## Theorem V2. A moving cutoff separates positive terms from a small remainder

Put \(m=r+1\ge1\). For any integer \(k\ge1\), define
\[
H=(16ak)^{1/3},\qquad X=\sqrt{ak}.
\tag{5}
\]
Assume
\[
H\ge4m,\quad H^2\ge a+1/4,\quad X\ge H,\quad X\ge1000.
\tag{6}
\]
Then
\[
\boxed{
q(a)D_{r,k}>
\frac1{6000(6k)^m}-\exp(-k^{1/3}/8).
}
\tag{7}
\]

**Proof, positive terms above H.** The phase bounds from S3 of the
preceding note apply to every zero with \(\gamma>1\):
\[
|\arg y|<1/\gamma,\qquad
|\arg(1-y)|<2a/\gamma^3.
\tag{8}
\]
These follow from \(|\delta|<1/2\); they do not assume \(\delta=0\).
One can use continuous phases in the right half-planes of
\(\mu\) and \(a+\mu\), so the phase of \(y^m(1-y)^k\), before
reduction modulo \(2\pi\), has absolute value at most
\[
\frac m\gamma+\frac{2ak}{\gamma^3}
\le\frac14+\frac18=\frac38
\qquad(\gamma\ge H).
\tag{9}
\]
Its real part is therefore greater than half its magnitude.
Every term with \(\gamma\ge H\) is positive in real part.
Terms with \(\gamma\le T\) are also positive by the retained
finite-height theorem. Thus only
\[
\mathcal B=\{\rho:T<\gamma<H\}
\]
can make a negative contribution. If \(H\le T\), this set is empty.

**Proof, attenuation on the potentially negative interval.** For
\(\rho\in\mathcal B\), write \(\mu=\gamma^2-\delta^2-2i\gamma\delta\).
Then
\[
1-|1-y|^2
=\frac{a(a+2\gamma^2-2\delta^2)}{|a+\mu|^2}
\ge\frac{a}{2(a+\gamma^2+1/4)}.
\tag{10}
\]
Indeed \(|a+\mu|\le a+\gamma^2+1/4\), and
\[
a+2\gamma^2-2\delta^2
\ge\tfrac12(a+\gamma^2+1/4)
\]
because \(a>800\), \(\gamma>T>1\), and \(\delta^2<1/4\).
Using \(\sqrt{1-u}\le e^{-u/2}\), \(\gamma<H\), and (6), gives
\[
|1-y|^k
\le\exp\left(-\frac{ak}{4(a+H^2+1/4)}\right)
\le\exp\left(-\frac{ak}{8H^2}\right)
<\exp(-k^{1/3}/8).
\tag{11}
\]
For the last inequality,
\[
\frac{ak}{8H^2}
=\frac18\left(\frac a{256}\right)^{1/3}k^{1/3}
>\frac18 k^{1/3},
\]
since \(a>256\). Equations (3) and (11), and \(|y|^r\le1\), bound
the absolute contribution of all terms in \(\mathcal B\) by
\[
\sum_{\mathcal B}
\left|\frac{y^r(1-y)^k}{a+\mu}\right|
<\exp(-k^{1/3}/8).
\tag{12}
\]
This estimates the entire potentially negative interval, including
all multiplicities, with no zero-count truncation.

**Proof, an explicit positive contribution.** V1 supplies a zero with
\(X<\gamma\le2X\). By (6) it lies in the positive phase range.
For this zero,
\[
|y|=\frac a{|a+\mu|}
\ge\frac a{a+\gamma^2+\delta^2}
>\frac1{6k},
\tag{13}
\]
because \(\gamma^2\le4ak\), \(\delta^2<1/4<a\), and \(k\ge1\).
Also \(|\mu|=\gamma^2+\delta^2>ak\), so
\[
|1-y|=\frac{|\mu|}{|a+\mu|}
\ge\frac{|\mu|}{a+|\mu|}
>\frac{k}{k+1},\qquad
|1-y|^k>e^{-1}>1/3.
\tag{14}
\]
Combining (9), (13), and (14), its real contribution to (2) exceeds
\[
\frac1{6a}(6k)^{-m}>\frac1{6000}(6k)^{-m}.
\tag{15}
\]
All other terms outside \(\mathcal B\) have positive real part.
Subtract (12) from (15) to obtain (7). Absolute convergence justifies
the splitting and summing real parts. QED.

The height interval furnishing (15) may contain off-line zeros.
Its contribution is positive because of the kernel's controlled phase,
not because a zero has been presumed to lie on the critical line.

## Theorem V3. Explicit positivity for all sufficiently large k in each row

For \(m\ge1\), set
\[
L_m=\log(2m),\qquad
K(m)=10^9m^3[\log(2m)]^3.
\tag{16}
\]
For the actual source, every pair of nonnegative integers \(r,k\)
satisfying \(k\ge K(r+1)\) obeys
\[
\boxed{
D_{r,k}>\frac1{6000(6k)^{r+1}}>0.
}
\tag{17}
\]

**Proof of the cutoff conditions.** Since \(L_m>1/2\),
\[
k\ge K(m)>\frac{10^9}{8}m^3>10^8m^3.
\]
It follows that \(H^3=16ak>10^{12}m^3\), so \(H>10^4m\).
In particular \(H\ge4m\) and \(H^2>a+1/4\).
Also \(X^2=ak>8\cdot10^{10}m^3\), hence \(X>1000\).
The inequality \(X\ge H\) is equivalent to \(ak\ge256\), which
holds here. Thus V2 applies.

**Proof that the remainder is smaller than half the positive bound.**
For real \(k>0\), put
\[
G_m(k)=\frac{k^{1/3}}8-\log12000-m\log(6k).
\]
At \(k_0=K(m)\), the inequalities
\(\log12000<10\), \(\log6<2\), \(\log10<3\),
\(\log m\le L_m\), and \(\log L_m\le L_m\) give
\[
G_m(k_0)>
125mL_m-10-m(29+6L_m)
=119mL_m-29m-10>0.
\tag{18}
\]
The final inequality uses \(L_m>1/2\) and \(m\ge1\).
Furthermore,
\[
G'_m(k)=\frac1k\left(\frac{k^{1/3}}{24}-m\right)>0
\qquad(k\ge K(m)),
\tag{19}
\]
because \(k^{1/3}\ge1000mL_m>500m\).
Consequently
\[
\exp(-k^{1/3}/8)
<\frac1{12000(6k)^m}
\qquad(k\ge K(m)).
\]
Equation (7) now yields
\(q(a)D_{r,k}>[12000(6k)^m]^{-1}\).
Finally \(q(a)<1/2\) proves (17). QED.

Thus every fixed row has a proved positive tail. The constants are
conservative and no optimal threshold is claimed. The lower bound is
not uniform as either index increases.

## Theorem V4. Ten million complete binomial rows

For the actual source,
\[
\boxed{
D_{r,k}>0
\qquad 0\le r\le9{,}999{,}999,\quad k\ge0.
}
\tag{20}
\]

**Proof.** If \(m=r+1\le10^7\), then \(\log(2m)<20\), so
\[
K(m)<10^9(10^7)^3\,20^3
=8\cdot10^{33}<10^{34}.
\tag{21}
\]
The elementary logarithmic bound follows, for example, from
\(\log2<1\) and \(\log10<5/2\).
The preceding S3 already proves positivity for every
\(0\le r<3\cdot10^{12}\) and \(0\le k\le10^{34}\).
V3 proves positivity for \(k\ge K(m)\).
These ranges overlap and together cover every \(k\ge0\)
for the stated rows. QED.

These are rows of the scalar binomial array \(D_{r,k}\).
Their count is not a count of positive Hankel matrix dimensions
or newly verified zeta zeros.

## Theorem V5. The precise remaining sign condition

Retain the unbounded sector \(r\ge1,\ 0\le k\le60r\) and the rectangle
\(0\le r<T,\ 0\le k\le10^{34}\) from the preceding note.
After adding V3, the only unsettled pairs form
\[
\boxed{
\mathcal U_{\rm new}=
\left\{(r,k)\in\mathbb Z_{\ge0}^2:
\begin{array}{l}
k>60r,\quad k<K(r+1),\\
r\ge3\cdot10^{12}\ \text{or}\ k>10^{34}
\end{array}
\right\}.
}
\tag{22}
\]
There are no pairs of this set with \(r<10^7\).
For each fixed \(r\), its set of possible \(k\) is finite,
but \(r\) itself is unbounded. The retained ILC criterion now gives
\[
\boxed{\mathrm{RH}\ \Longleftrightarrow\
D_{r,k}\ge0\ \hbox{for every }(r,k)\in\mathcal U_{\rm new}.}
\tag{23}
\]

**Proof.** All pairs outside (22) belong to at least one of the three
proved positive regions. Conversely (22) is their exact complement
(the pair \((0,0)\) lies in the rectangle). Apply ILC-T2 and V4. QED.

An actual-source estimate throughout the remaining middle range is
still required. Neither rowwise eventual positivity nor the first
ten million complete rows supplies positivity at every pair.
No full RH proof, all-order Hankel positivity or global eta improvement
is asserted.

## Evidence and execution

| Result | Evidence |
|---|---|
| A zero in every distant interval \((X,2X]\) | Rosser's count bound and V1; no critical-line assumption |
| All terms above a moving cutoff have positive real part | Written phase bound (8)-(9) |
| Potentially negative terms are collectively small | Full source budget (3), attenuation (10)-(12) |
| Every fixed row has an explicit positive tail | V3 |
| First ten million binomial rows are positive for all k | V4, combining the new tail with the preceding rectangle |
| Every pair in the remaining set is nonnegative | Open |
| New numerical evaluation, zero count or machine certificate | None |

This written continuation reuses the same source and existing certificate
lineage in Publications. The additional zero-count theorem is explicitly
cited. No local scientific computation, independent agent or scientific
CI stage is introduced. The documentation commit uses skip-CI to avoid
replaying unchanged checks. Engines, predecessor proofs, numerical
certificates and workflows are unchanged. The new arguments are written
proofs, not newly machine-checked certificates.

Classical zero counting and the canonical product keep their attribution.
No literature-wide novelty claim is made for the explicit consequences
derived here.
