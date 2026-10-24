# RH continuation: a quadratic threshold from phase cancellation

17 October 2026. Continues [the large-difference tail](RH_LARGE_DIFFERENCE_TAIL.md)
at Publications commit `ff440839a76d827f51c01d1a73f74124b418f213`.

The preceding proof treated the two phase contributions separately.
Here their cancellation near the maximum of the actual kernel gives
a stronger result:
\[
D_{r,k}>0\quad\hbox{whenever}\quad
k\ge Q(r+1),\qquad
Q(m)=10^8m^2[\log(2m)]^3.
\]
The power of the row index has fallen from three to two. Combined with
the retained rectangle, this covers every difference order in the first
50 billion scalar binomial rows. The remaining row index is still
unbounded; an RH proof is not concluded.

## Unchanged source and cited inputs

Retain the actual xi source and notation
\[
a=3249/4,\qquad
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,
\]
\[
c_j=\frac{(-a)^j q^{(j)}(a)}{j!\,q(a)},\qquad
D_{r,k}=\sum_{j=0}^k(-1)^j\binom{k}{j}c_{r+j}.
\tag{1}
\]
For each nontrivial zero \(\rho=\beta+i\gamma\), \(\gamma>0\),
including multiplicity, put
\[
\delta=\beta-\tfrac12,\quad
\mu=(\gamma-i\delta)^2,\quad y=\frac a{a+\mu}.
\]
The retained source identities are
\[
q(a)D_{r,k}
=\frac1a\sum_{\gamma>0}y^{r+1}(1-y)^k
=\sum_{\gamma>0}\frac{y^r(1-y)^k}{a+\mu},
\tag{2}
\]
\[
800<a<900,\quad 0<q(a)<1/2,\quad
\sum_{\gamma>T}|a+\mu|^{-1}<1,\qquad T=3\cdot10^{12}.
\tag{3}
\]
The sums converge absolutely. All zeros with \(0<\gamma\le T\)
are on the critical line, so their terms in (2) are positive.

The source and absolute remainder budget are proved in S1 of
[the positive-sector note](RH_POSITIVE_LAYER_SECTORS.md).
The finite-height input is Platt and Trudgian,
*The Riemann hypothesis is true up to \(3\cdot10^{12}\)*,
Bulletin of the London Mathematical Society 53 (2021), 792–797,
[Theorem 1](https://arxiv.org/pdf/2004.09765),
[DOI 10.1112/blms.12460](https://doi.org/10.1112/blms.12460).
Its computation is not rerun here.

The only zero-counting input is the same coarse classical estimate
used in the preceding continuation:
\[
|N(x)-F(x)|\le\log x+10,\qquad
F(x)=\frac{x}{2\pi}\log\frac{x}{2\pi e},\qquad x\ge e.
\tag{4}
\]
It follows from Rosser's constants \(0.137,0.443,2.463\), with threshold 2,
in the normalization without a \(7/8\) constant. These are reproduced
in equation (1.1) and Table 1 of Hasanalizade, Shen and Wong,
*Counting zeros of the Riemann zeta function*, Journal of Number Theory
235 (2022), 219–241:
[published author-hosted paper](https://www-math.nsysu.edu.tw/~pjwong/stuff/CountingRiemannZeros.pdf),
[DOI 10.1016/j.jnt.2021.06.032](https://doi.org/10.1016/j.jnt.2021.06.032).
The original source is Rosser, *Explicit bounds for some functions of
prime numbers*, American Journal of Mathematics 63 (1941), 211–232,
[DOI 10.2307/2371291](https://doi.org/10.2307/2371291).
Here \(N\) counts every nontrivial zero of positive height, with
multiplicity, without assuming its real part is \(1/2\).

## Theorem C1. A zero within a fixed height increment

Every real \(X\ge1000\) has a nontrivial zero with
\[
X<\gamma\le X+100.
\tag{5}
\]
No critical-line assertion is made.

**Proof.** For \(x\ge X\ge1000\), \(2\pi<8\) and
\(X>(2\pi)^2\) imply
\[
F'(x)=\frac{\log(x/(2\pi))}{2\pi}>
\frac{\log X}{16}.
\]
Since \(X+100\le2X\), equation (4) gives
\[
N(X+100)-N(X)>
\frac{100}{16}\log X-2\log(2X)-20.
\]
The error \(2\log(2X)+20\) is less than \(6\log X\):
use \(\log X>6\) and \(2\log2+20<22\).
The right side is therefore greater than \((1/4)\log X>0.
This proves (5). QED.

## Working range and the real comparison kernel

For C2–C4 set
\[
m=r+1\ge10^7,\quad L=\log(2m),\quad
Q(m)\le k\le K(m),
\tag{6}
\]
where
\[
Q(m)=10^8m^2L^3,\qquad K(m)=10^9m^3L^3.
\]
All \(r,k\) under consideration are integers. Define
\[
n=r+k,\quad p=\frac r{r+k},\quad
M=p^r(1-p)^k,\quad
G=\sqrt{\frac{ak}{r}},\quad
w=10\sqrt{\frac Lr}.
\tag{7}
\]
Then
\[
0<p<1/8,\quad 0<w<1/4,\quad
G\ge10^4\sqrt r>1000.
\tag{8}
\]
For the bound on \(w\), the function
\(\log(2m)/(m-1)\) decreases for \(m\ge2\).
At \(m=10^7\), use \(\log(2m)<20\) to obtain
\(L/r<1/1600\). The remaining bounds follow directly from
\(k\ge10^8m^2L^3\), \(L>1\), and \(a>800\).

For each zero define the real comparison variable
\[
x=\gamma^2+\delta^2=|\mu|,\qquad
t=\frac a{a+x},\qquad s=\frac tp.
\tag{9}
\]
The real profile \(t^r(1-t)^k\) has its unique maximum \(M\) at
\(t=p\). At \(\delta=0\), this corresponds exactly to height \(G\).

## Theorem C2. Positive phase in a window about the maximum

For every zero with \(\gamma>T\) and \(|s-1|\le w\),
\[
\Re\{y^{r+1}(1-y)^k\}>
\frac12|y^{r+1}(1-y)^k|.
\tag{10}
\]
Thus each such zero contributes positively to (2), even if it is off
the critical line.

**Proof.** Fix this zero and let \(z=\gamma-iv\), with \(v\) between
0 and \(\delta\). Use continuous logarithms along this segment for
\[
f(z)=\frac{a^{r+1}z^{2k}}{(a+z^2)^{r+1+k}}.
\]
At \(v=0\), \(f(\gamma)>0\), while \(f(\gamma-i\delta)\)
is the expression in braces in (10). Logarithmic differentiation gives
\[
\frac{f'(z)}{f(z)}
=\frac{2(ka-(r+1)z^2)}{z(a+z^2)}.
\tag{11}
\]
Integration along the segment and \(2|\delta|<1\) give an upper bound
on the absolute continuous phase by the supremum of
\[
\frac{|ka-(r+1)z^2|}{|z|\,|a+z^2|}.
\tag{12}
\]
This step retains the cancellation in the numerator.

Since \(x=\gamma^2+\delta^2\),
\[
|ka-(r+1)z^2|
\le |ka-rx|+r|x-z^2|+|z^2|,
\]
\[
|ka-rx|=an\,\frac{|s-1|}{s},\quad
|x-z^2|\le\gamma+\tfrac12,\quad
|z^2|\le\gamma^2+\tfrac14.
\]
Also
\[
|z|\ge\gamma,\quad
|a+z^2|\ge a+\gamma^2-\tfrac14
\ge\tfrac12(a+x),\quad
|a+z^2|>\gamma^2.
\]
Apply the first denominator bound to \(|ka-rx|\) and the second to
the other two terms. The phase is at most
\[
\frac{2r|s-1|}{\gamma}+\frac{2r}{\gamma^2}+\frac2\gamma.
\tag{13}
\]
Because \(t=ps<2p\), \(p<1/8\), and \(a>800\),
\[
\gamma^2=\frac a{ps}-a-\delta^2
\ge\frac a{4p}\ge\frac{ak}{4r}.
\]
Therefore the first two terms in (13) are bounded respectively by
\[
\frac{4r^{3/2}w}{\sqrt{ak}}
=\frac{40r\sqrt L}{\sqrt{ak}}<\frac1{100},
\qquad
\frac{8r^2}{ak}<\frac1{100}.
\]
Here use \(k\ge10^8m^2L^3\), \(m\ge r\), and \(L>1\).
The last term is \(2/\gamma<2/T<1/100\).
The absolute continuous phase is less than \(3/100<1/4\).
Its cosine exceeds \(1/2\), proving (10). QED.

The cancellation in (11) is controlled in a finite window.
It has not been asserted for every height or every index pair.

## Theorem C3. A small combined contribution outside the positive window

For the high zeros outside this window,
\[
\sum_{\substack{\gamma>T\\ |s-1|>w}}
\left|\frac{y^r(1-y)^k}{a+\mu}\right|
<2M e^{-16L}.
\tag{14}
\]

**Proof.** The exact modulus identity
\[
|a+\mu|^2=(a+x)^2-4a\delta^2
\]
gives
\[
h:=|y|^r|1-y|^k
=t^r(1-t)^k
\left(1-\frac{4\delta^2t^2}{a}\right)^{-n/2}.
\tag{15}
\]
For high zeros, \(0<t<a/T^2<10^{-22}\). As
\(-\log(1-u)\le2u\) for \(0\le u\le1/2\), it follows that
\[
\log h\le r\log t+k\log(1-t)+\frac{nt^2}{a}.
\tag{16}
\]
Write
\[
J(s)=s-1-\log s\ge0.
\]
Concavity of \(\log(1-t)\), using its tangent at \(p\), gives
\[
r\log t+k\log(1-t)-\log M\le-rJ(s).
\tag{17}
\]
For \(s\le2\), the correction in (16) is
\[
\frac{nt^2}{a}=\frac{rp s^2}{a}
\le\frac{4r^2}{ak}<\frac1{100}.
\]
For \(s\ge2\), \(J(s)\ge s/10\). To see this, evaluate
\(J(s)-s/10\) at 2, using \(\log2<3/4\), and note that its derivative
\(9/10-1/s\) is positive for \(s\ge2\).
Thus in the second case
\[
\frac{nt^2}{a}=\frac{rt s}{a}
\le\frac{10t}{a}\,rJ(s)<\tfrac12 rJ(s).
\]
Both cases yield the uniform envelope
\[
h<2M\exp\{-rJ(s)/2\}.
\tag{18}
\]
For \(0<w<1/4\), monotonicity on each side of \(s=1\) and
integration of \(J'(s)=(s-1)/s\) imply
\[
|s-1|>w\quad\Longrightarrow\quad J(s)\ge w^2/3.
\]
Consequently
\[
h<2M e^{-rw^2/6}
=2M e^{-(100/6)L}<2M e^{-16L}
\]
outside the window. Multiply by \(|a+\mu|^{-1}\), sum, and apply
the complete source budget (3). This proves (14). QED.

This is an absolute bound on every high-zero term outside the window.
It includes any potentially negative term; no high zero is omitted.

## Theorem C4. A positive contribution near the maximum dominates

In the working range (6),
\[
\boxed{
q(a)D_{r,k}>\frac{pM}{16a},\qquad
D_{r,k}>\frac{pM}{8a}>0.
}
\tag{19}
\]

**Proof of a positive contribution.** C1 supplies a zero with
\(G<\gamma\le G+100\). Its \(|\delta|<1/2\) is retained; if its
height exceeds \(T\), no critical-line location is assumed.
For this zero let
\[
u=\frac{\gamma^2+\delta^2-G^2}{a+G^2}.
\]
Because \(G\ge10^4\sqrt r\),
\[
0<u\le\frac{200}{G}+\frac{10001}{G^2}
<\frac1{10\sqrt r}.
\]
Since \(p=a/(a+G^2)\), it follows that
\[
s=\frac1{1+u},\qquad
|s-1|<\frac1{10\sqrt r}<w,\qquad t>p/2.
\tag{20}
\]
Put \(g(s)=r\log s+k\log(1-ps)\), so \(g'(1)=0\).
On \(1/2\le s\le3/2\), \(p<1/8<1/4\) gives
\[
|g''(s)|=\frac r{s^2}+\frac{kp^2}{(1-ps)^2}\le5r.
\]
Taylor's theorem and (20) imply
\[
t^r(1-t)^k\ge M e^{-1/40}>M/2.
\tag{21}
\]
The factor correcting this profile in (15) is at least one, and
\(|y|\ge t\). If this zero has \(\gamma>T\), C2 bounds the
cosine of the full phase below by \(1/2\). If \(\gamma\le T\),
its phase is zero by the finite-height theorem. In either case,
its real contribution to (2) exceeds
\[
\frac1{2a}|y|\,|y|^r|1-y|^k
>\frac{pM}{8a}.
\tag{22}
\]

**Proof that the remaining possibly negative part is smaller.**
Since \(r=m-1\ge m/2\), \(k\le K(m)\), and \(r<K(m)\),
\[
p=\frac r{r+k}\ge\frac1{4\cdot10^9m^2L^3}.
\tag{23}
\]
Using \(a<900\),
\[
\frac{p}{16a}>\frac1{10^{14}m^2L^3}.
\]
On the other hand, \(L<m\) and \(m\ge10^7\) imply
\[
2e^{-16L}=\frac2{(2m)^{16}}
<\frac1{10^{14}m^2L^3}
<\frac p{16a}.
\tag{24}
\]
For the middle comparison it suffices that
\(m^{11}>2\cdot10^{14}\), which holds throughout this range.
Equations (14) and (24) bound the absolute combined contribution
outside the positive window by \(pM/(16a)\).

All low-zero terms and all high-zero terms inside the window have
positive real parts. The selected term (22) is inside that window.
Subtracting the outside bound from (22) proves the first inequality
in (19). Absolute convergence justifies the splitting. Dividing by
\(0<q(a)<1/2\) proves the second. QED.

## Theorem C5. A quadratic large-difference threshold for every row

For all nonnegative integers \(r,k\),
\[
\boxed{
k\ge10^8(r+1)^2[\log(2r+2)]^3
\quad\Longrightarrow\quad D_{r,k}>0.
}
\tag{25}
\]

**Proof.** If \(m=r+1\le10^7\), every \(k\) was already settled
by V4 in the preceding note. For \(m>10^7\) and
\(Q(m)\le k\le K(m)\), apply C4. For \(k>K(m)\), apply V3
of the preceding note. These ranges cover every pair in (25). QED.

The explicit lower bound (19) was proved on the working range (6);
only positivity is asserted uniformly over the full range (25).

## Theorem C6. Fifty billion complete scalar binomial rows

For the unchanged actual source,
\[
\boxed{
D_{r,k}>0\qquad
0\le r\le49{,}999{,}999{,}999,\quad k\ge0.
}
\tag{26}
\]

**Proof.** For \(m=r+1\le5\cdot10^{10}\),
\(\log(2m)\le\log(10^{11})<28\), using \(\log10<5/2\).
Hence
\[
Q(m)<10^8(5\cdot10^{10})^2\,28^3<10^{34};
\tag{27}
\]
for the last inequality, \(28^3<40000\) suffices.
The retained S3 rectangle covers
\(r<3\cdot10^{12}\) and \(k\le10^{34}\).
C5 covers \(k\ge Q(r+1)\). These intervals overlap and cover
every \(k\ge0\) for the stated rows. QED.

This is a uniform analytic consequence. No fifty-billion-entry array
was computed. The number counts rows of \(D_{r,k}\), not positive
Hankel matrix dimensions, newly checked zeta zeros, or independent
certificates.

## Theorem C7. Updated remaining sign condition

The previously proved sector \(r\ge1,\ k\le60r\) and rectangle
\(r<T,\ k\le10^{34}\) remain valid. With C5, the exact complement
of these three stated positive regions is
\[
\mathcal U_{\rm quad}=
\left\{(r,k)\in\mathbb Z_{\ge0}^2:
\begin{array}{l}
k>60r,\quad k<Q(r+1),\\
r\ge3\cdot10^{12}\ \text{or}\ k>10^{34}
\end{array}
\right\}.
\tag{28}
\]
There are no pairs here with \(r<5\cdot10^{10}\).
The retained ILC-T2 criterion gives
\[
\boxed{
\mathrm{RH}\ \Longleftrightarrow\
D_{r,k}\ge0\quad\hbox{for every }(r,k)\in\mathcal U_{\rm quad}.
}
\tag{29}
\]

**Proof.** Every pair outside (28) lies in one of the three proved
regions; the exceptional sector endpoint \((0,0)\) is already in
the rectangle. Apply the all-pairs equivalence ILC-T2 and C6. QED.

Each fixed row has a finite remaining middle interval, but the row
index is unbounded. The present phase and modulus estimates do not
settle that interval for every row. No full RH proof, all-order
Hankel positivity, or global eta improvement is asserted.

## Provenance and verification scope

| Item | Evidence |
|---|---|
| A zero within 100 units of a prescribed height \(X\ge1000\) | Rosser's cited coarse count and C1 |
| Phase cancellation about the kernel maximum | Exact derivative (11), proved window estimate C2 |
| Full contribution outside the window | Exact modulus (15), envelope and source budget C3 |
| Quadratic threshold | Written proofs C4–C5 |
| Fifty billion complete scalar rows | C6, combining C5 with the retained rectangle |
| Positivity on every remaining pair | Open |
| New numerical evaluation or machine certificate | None |

The arguments are written mathematical proofs, not newly machine-checked
certificates. No local scientific computation, independent agent,
zero-count replay, or new scientific CI stage is used. Source engines,
earlier certificates, and workflows are unchanged. The documentation
commit uses [skip ci]. Only this note and its README entry are added.

All classical source results keep their attribution. The constants are
conservative; no optimality or literature-wide novelty claim is made.
