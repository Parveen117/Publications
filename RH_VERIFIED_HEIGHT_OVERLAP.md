# RH continuation: verified-height overlap and a sharper phase tail

17 October 2026. Continues [the quadratic phase-window argument](RH_QUADRATIC_PHASE_WINDOW.md)
at Publications commit `c44bf2f30f601c0514eeb1796fc3598b92cda256`.

This continuation strengthens both sides of the unsettled middle range.
A zero within the already verified height gives a much wider positive
sector. Rechecking the phase-window proof with sharper constants reduces
its sufficient tail threshold to \(100(r+1)^2\log(2r+2)\).
Their overlap proves positivity for every difference order in each
scalar binomial row \(0\le r<10^{17}\).

The result uses the same actual source and the same cited finite-height
theorem. No new zero verification or numerical array calculation is
performed. The remaining row index is unbounded, so RH is not concluded.

## Retained source and evidence

Use
\[
a=3249/4,\quad \mathcal X(w^2)=\xi(1/2+w),\quad
q=\mathcal X'/\mathcal X,\quad
c_j=\frac{(-a)^jq^{(j)}(a)}{j!\,q(a)},
\]
\[
D_{r,k}=\sum_{j=0}^k(-1)^j\binom{k}{j}c_{r+j}.
\tag{1}
\]
For every nontrivial zero \(\rho=\beta+i\gamma\), \(\gamma>0\),
with multiplicity retained, let
\[
\delta=\beta-\tfrac12,\qquad
\mu=(\gamma-i\delta)^2,\qquad y=\frac a{a+\mu}.
\]
The unchanged actual-source identities and bounds are
\[
q(a)D_{r,k}
=\frac1a\sum_{\gamma>0}y^{r+1}(1-y)^k
=\sum_{\gamma>0}\frac{y^r(1-y)^k}{a+\mu},
\tag{2}
\]
\[
800<a<900,\qquad 0<q(a)<1/2,\qquad
\sum_{\gamma>T}|a+\mu|^{-1}<1,\quad T=3\cdot10^{12}.
\tag{3}
\]
These sums converge absolutely. All terms with \(\gamma\le T\) are
positive real numbers. With
\[
u_T=\frac a{a+T^2-1/4},
\]
S1 of [the positive-sector note](RH_POSITIVE_LAYER_SECTORS.md) gives
\[
\left|q(a)D_{r,k}
-\sum_{0<\gamma\le T}
\frac1{a+\gamma^2}
\left(\frac a{a+\gamma^2}\right)^r
\left(\frac{\gamma^2}{a+\gamma^2}\right)^k
\right|<u_T^r.
\tag{4}
\]

The finite-height input is the theorem of Platt and Trudgian,
*The Riemann hypothesis is true up to \(3\cdot10^{12}\)*,
Bulletin of the London Mathematical Society 53 (2021), 792–797:
[Theorem 1](https://arxiv.org/pdf/2004.09765),
[DOI 10.1112/blms.12460](https://doi.org/10.1112/blms.12460).
We use the conservative height \(T\) above. No assertion about the
critical-line location of higher zeros is imported.

C1 of [the preceding phase-window note](RH_QUADRATIC_PHASE_WINDOW.md)
also supplies the retained consequence
\[
X\ge1000\quad\Longrightarrow\quad
\text{some zero has }X<\gamma\le X+100.
\tag{5}
\]
It follows from the coarse Rosser estimate
\[
\left|N(X)-\frac X{2\pi}\log\frac X{2\pi e}\right|
\le\log X+10\qquad(X\ge e).
\]
The Rosser constants and counting normalization are reproduced in
equation (1.1) and Table 1 of Hasanalizade, Shen and Wong,
*Counting zeros of the Riemann zeta function*, Journal of Number Theory
235 (2022), 219–241:
[published author-hosted paper](https://www-math.nsysu.edu.tw/~pjwong/stuff/CountingRiemannZeros.pdf),
[DOI 10.1016/j.jnt.2021.06.032](https://doi.org/10.1016/j.jnt.2021.06.032).
The historical source is Rosser, *Explicit bounds for some functions
of prime numbers*, American Journal of Mathematics 63 (1941), 211–232,
[DOI 10.2307/2371291](https://doi.org/10.2307/2371291).
Statement (5) counts zeros without restricting their real parts.

## Theorem O1. A wide positive sector from the verified height

Set
\[
B_T=\frac{T^2}{10a}.
\tag{6}
\]
Then
\[
10^{21}<B_T<2\cdot10^{21},
\qquad
\boxed{D_{r,k}>0\quad\text{whenever }0\le k\le B_T r.}
\tag{7}
\]
For \(r\ge100\), the proof gives the explicit lower bound
\[
D_{r,k}>
2u_T^r\left(\frac{2^r}{T^2}-1\right)>0.
\tag{8}
\]

**Proof.** Apply (5) with \(X=T/2-100\). It supplies a zero
\[
T/2-100<\gamma_\star\le T/2.
\]
The finite-height theorem places this zero on the critical line.
Since \(T>2000\), it satisfies
\[
T^2/5<\gamma_\star^2\le T^2/4.
\tag{9}
\]
Indeed \(T/2-100>9T/20\), and \((9/20)^2>1/5\).

Let \(v=a/(a+\gamma_\star^2)\). The coarse source bounds imply
\[
\frac1{a+\gamma_\star^2}>\frac1{T^2},\qquad
\frac v{u_T}\ge
\frac{a+T^2-1/4}{a+T^2/4}>\frac72.
\tag{10}
\]
For the last inequality it suffices that \(T^2>20a+2\).
Furthermore,
\[
(1-v)^k
=\exp\{-k\log(1+a/\gamma_\star^2)\}
\ge\exp(-5ak/T^2)
\ge e^{-r/2}
\quad(k\le B_Tr).
\tag{11}
\]
The selected low-zero term in (4) consequently exceeds
\[
\frac{u_T^r}{T^2}
\left(\frac72 e^{-1/2}\right)^r
>\frac{2^r}{T^2}u_T^r.
\tag{12}
\]
Here \(e<3<49/16\) implies \(e^{1/2}<7/4\).
All other low-zero terms are positive, while (4) bounds the full
high-zero remainder. Thus
\[
q(a)D_{r,k}>
u_T^r\left(\frac{2^r}{T^2}-1\right).
\]
For \(r\ge100\), \(2^r\ge2^{100}>10^{30}>T^2\);
use \(2^{10}>10^3\). Division by \(q(a)<1/2\) proves (8).

For \(0\le r<100\), the range \(k\le B_Tr\) has
\(k<2\cdot10^{23}<10^{34}\) and \(r<T\).
It lies in the retained S3 positive rectangle. This includes
\(r=k=0\). The inequalities on \(B_T\) follow from \(800<a<900\)
and \(T^2=9\cdot10^{24}\). QED.

This proof uses one critical-line zero whose existence follows from
the cited count and whose location follows from the cited finite-height
theorem. It does not require its ordinate to be computed or isolated.

## Theorem O2. Recalibrating the phase-window proof

Let
\[
m=r+1\ge10^7,\quad L=\log(2m),\quad
Q_*(m)=100m^2L,\quad K(m)=10^9m^3L^3.
\]
For integers \(r,k\) satisfying
\[
Q_*(m)\le k\le K(m),
\tag{13}
\]
put
\[
p=\frac r{r+k},\quad M=p^r(1-p)^k,\quad
G=\sqrt{\frac{ak}{r}},\quad w=10\sqrt{\frac Lr}.
\]
Then
\[
\boxed{
q(a)D_{r,k}>\frac{pM}{256a},\qquad
D_{r,k}>\frac{pM}{128a}>0.
}
\tag{14}
\]

**Proof.** The algebraic phase and modulus identities used in C2–C4
of the preceding note remain exact. The following checks establish
the required estimates on the enlarged range (13); they do not assume
the preceding, stronger lower bound on \(k\).

First,
\[
L>2,\quad p<1/8,\quad w<1/4,\quad
G^2=\frac{ak}{r}\ge100arL>160000r.
\tag{15}
\]
Thus \(G>400\sqrt r>1000\).
The estimate for \(w\) uses the same decreasing function
\(\log(2m)/(m-1)\) and its bound at \(m=10^7\) as in C2.

For any zero, define
\[
x=\gamma^2+\delta^2,\qquad
t=\frac a{a+x},\qquad s=t/p.
\]

**Positive phase near the maximum.** At a high zero with
\(|s-1|\le w\), the exact logarithmic derivative
\[
\frac{d}{dz}\log
\left(\frac{a^{r+1}z^{2k}}{(a+z^2)^{r+1+k}}\right)
=\frac{2(ka-(r+1)z^2)}{z(a+z^2)}
\]
integrated along \(z=\gamma-iv\), \(0\le |v|\le|\delta|<1/2\),
gives the phase bound
\[
|\arg_{\rm cont}\{y^{r+1}(1-y)^k\}|
\le\frac{2r|s-1|}{\gamma}+\frac{2r}{\gamma^2}+\frac2\gamma.
\tag{16}
\]
This is the derivation in C2, retaining the cancellation in the
numerator; its conditions are \(\gamma>T\), \(p<1/8\), and \(w<1/4\).
Also \(\gamma^2\ge ak/(4r)\) in this window.
The three terms in (16) are bounded by
\[
\frac{40r\sqrt L}{\sqrt{ak}}
\le\frac4{\sqrt a}<\frac16,\qquad
\frac{8r^2}{ak}\le\frac8{100aL}<\frac1{1000},\qquad
\frac2T<\frac1{1000}.
\]
Their sum is less than \(1/4\). Consequently every high-zero term
inside this window has real part greater than half its magnitude.

**Full contribution outside the window.** Write \(n=r+k\) and
\(J(s)=s-1-\log s\). As in C3, the exact modulus identity is
\[
h:=|y|^r|1-y|^k
=t^r(1-t)^k
\left(1-\frac{4\delta^2t^2}{a}\right)^{-n/2}.
\tag{17}
\]
For high zeros, \(t<10^{-22}\), so
\[
\log(h/M)\le-rJ(s)+\frac{nt^2}{a}.
\tag{18}
\]
This follows by taking the tangent of \(\log(1-t)\) at \(p\)
and using \(-\log(1-u)\le2u\).
For \(s\le2\), the correction now satisfies
\[
\frac{nt^2}{a}\le\frac{4r^2}{ak}
\le\frac4{100aL}<\frac1{100}.
\]
For \(s\ge2\), \(J(s)\ge s/10\) and
\[
\frac{nt^2}{a}=\frac{rt s}{a}<\tfrac12 rJ(s).
\]
Therefore the same envelope holds throughout the high-zero range:
\[
h<2M e^{-rJ(s)/2}.
\]
Since \(|s-1|>w\) implies \(J(s)\ge w^2/3\), the complete
source budget (3) yields
\[
\sum_{\substack{\gamma>T\\|s-1|>w}}
\left|\frac{y^r(1-y)^k}{a+\mu}\right|
<2M e^{-16L}.
\tag{19}
\]

**One positive contribution.** By (5), a zero exists with
\(G<\gamma\le G+100\). Let
\[
u=\frac{\gamma^2+\delta^2-G^2}{a+G^2}.
\]
Equation (15) implies
\[
0<u\le\frac{200}{G}+\frac{10001}{G^2}
<\frac1{2\sqrt r}+\frac1{8r}<\frac1{\sqrt r}.
\]
As \(p=a/(a+G^2)\),
\[
s=\frac1{1+u},\quad |s-1|<1/\sqrt r<w,\quad t>p/2.
\tag{20}
\]
On \(1/2\le s\le3/2\), the function
\(g(s)=r\log s+k\log(1-ps)\) obeys \(g'(1)=0\) and
\(|g''(s)|\le5r\). Taylor's theorem gives the slightly weaker
positive-profile bound
\[
t^r(1-t)^k>M e^{-5/2}>M/32.
\tag{21}
\]
For the last inequality use \(e^{5/2}<e^3<27<32\).
The modulus correction in (17) is at least one, and \(|y|\ge t\).
This selected zero has phase cosine greater than \(1/2\) if its
height exceeds \(T\), by (16); otherwise its phase is zero by the
finite-height theorem. Its real contribution to (2) is therefore
greater than
\[
\frac1{2a}|y|\,|y|^r|1-y|^k>\frac{pM}{128a}.
\tag{22}
\]

**Dominance over the entire possibly negative part.** The upper
bound \(k\le K(m)\), together with \(r\ge m/2\), gives
\[
p\ge\frac1{4\cdot10^9m^2L^3}.
\]
Since \(a<900\), \(L<m\), and \(m\ge10^7\),
\[
\frac{p}{256a}>
\frac1{10^{15}m^2L^3}>
2(2m)^{-16}=2e^{-16L}.
\tag{23}
\]
The second comparison follows, for example, from
\(m^{11}>2\cdot10^{15}\) and \(L<m\).
The absolute remainder (19) is thus smaller than \(pM/(256a)\).

All low-zero terms and all high-zero terms inside the window are
positive in real part. The selected term (22) lies inside the window.
Subtract (19) from (22), using (23), to obtain the first bound in (14).
Absolute convergence justifies this splitting, and \(q(a)<1/2\)
gives the second bound. QED.

The logarithmic factor has decreased from a cube to the first power.
The explicit lower bound (14) is claimed on (13), not automatically
on every pair beyond its upper endpoint.

## Theorem O3. A sharper positive tail in every row

For all nonnegative integers \(r,k\),
\[
\boxed{
k\ge100(r+1)^2\log(2r+2)
\quad\Longrightarrow\quad D_{r,k}>0.
}
\tag{24}
\]

**Proof.** For \(r+1\le10^7\), V4 of
[the large-difference note](RH_LARGE_DIFFERENCE_TAIL.md) already
proves positivity for all \(k\). For larger rows and
\(Q_*(r+1)\le k\le K(r+1)\), apply O2.
For \(k>K(r+1)\), apply V3 of that same preceding note.
These cases cover (24). QED.

## Theorem O4. All difference orders below row index \(10^{17}\)

For the actual source,
\[
\boxed{
D_{r,k}>0\qquad 0\le r<10^{17},\quad k\ge0.
}
\tag{25}
\]

**Proof.** For \(1\le r<10^{17}\), let \(m=r+1\le10^{17}\).
Then \(m^2/(m-1)\le2m\) and
\(\log(2m)<44\), using \(\log2<1\) and \(\log10<5/2\).
Consequently
\[
\frac{Q_*(r+1)}r
=\frac{100m^2\log(2m)}{m-1}
\le200m\log(2m)
<8.8\cdot10^{20}<10^{21}<B_T.
\tag{26}
\]
The O1 sector \(k\le B_Tr\) and the O3 tail
\(k\ge Q_*(r+1)\) overlap and together cover all \(k\ge0\).
The remaining row \(r=0\) was already completely covered by V4
of the large-difference note. QED.

The number \(10^{17}\) counts scalar binomial rows, with every
difference order covered analytically in each row. It does not count
newly verified zeros or positive Hankel matrix dimensions. Nothing
like a \(10^{17}\)-entry array is generated.

## Theorem O5. The remaining RH condition is one explicit middle band

Define
\[
\mathcal U_{\rm overlap}=
\left\{(r,k)\in\mathbb Z_{\ge0}^2:
\begin{array}{l}
r\ge10^{17},\\
B_Tr<k<100(r+1)^2\log(2r+2)
\end{array}
\right\}.
\tag{27}
\]
Then the retained ILC-T2 criterion gives
\[
\boxed{
\mathrm{RH}\ \Longleftrightarrow\
D_{r,k}\ge0\quad\text{for every }(r,k)\in\mathcal U_{\rm overlap}.
}
\tag{28}
\]

**Proof.** O4 covers all lower rows, O1 covers the lower side of
every remaining row, and O3 covers its upper side. Their exact
complement is (27). Every previously stated positive sector and
rectangle is contained in this union. ILC-T2 identifies RH with
nonnegativity at every index pair, giving (28). QED.

The interval in each fixed row is finite; the row index is unbounded.
The present estimates do not control every pair in (27).
Thus the work has narrowed the sign problem but has not proved RH,
all-order Hankel positivity, or a global eta improvement.

## Evidence and execution

| Result | Evidence |
|---|---|
| Positive sector \(k\le B_Tr\) | One zero inside the verified height; full remainder bound O1 |
| Sharper tail \(k\ge100(r+1)^2\log(2r+2)\) | Explicit rechecking of phase, modulus, and dominance bounds O2–O3 |
| Every \(k\) in rows \(r<10^{17}\) | Overlap inequality (26) |
| Nonnegativity on the remaining unbounded band | Open |
| New numerical computation or machine certificate | None |

Only this note and the root README entry are changed. The arguments
are written proofs, not newly machine-checked certificates.
No local scientific computation, independent agents, numerical replay,
or new scientific CI stage is used. Source engines, prior certificates,
and workflows are unchanged; the documentation commit uses [skip ci].
The cited classical results retain their attribution. No optimality or
literature-wide novelty claim is made.
