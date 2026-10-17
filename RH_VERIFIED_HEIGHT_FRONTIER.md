# RH continuation: an exact high-zero envelope and a verified-height frontier

17 October 2026. Continues [the verified-height overlap](RH_VERIFIED_HEIGHT_OVERLAP.md)
at Publications commit `d626f7bf3c2d1487d06f3ede0c93397627869efe`.

The new step is an envelope for the entire unknown high-zero contribution,
using its dependence on the ratio \(k/r\). Below the verified height,
this envelope decreases with zero height. A nearby verified zero can
therefore be compared directly against the complete high-zero remainder.

This yields a row-dependent positive region approaching the verified
height, a uniform sector with kernel peak at most \(0.99T\), and all
difference orders in scalar rows \(0\le r<10^{18}\). It does not
increase the verified zero height \(T\), and does not conclude RH.

## Source, notation, and retained inputs

Retain the actual source
\[
a=3249/4,\quad \mathcal X(w^2)=\xi(1/2+w),\quad
q=\mathcal X'/\mathcal X,
\]
\[
c_j=\frac{(-a)^jq^{(j)}(a)}{j!\,q(a)},\qquad
D_{r,k}=\sum_{j=0}^k(-1)^j\binom{k}{j}c_{r+j}.
\tag{1}
\]
For nontrivial zeros \(\rho=\beta+i\gamma\), \(\gamma>0\),
including multiplicity, put
\[
\delta=\beta-\tfrac12,\quad
\mu=(\gamma-i\delta)^2,\quad
y=\frac a{a+\mu}.
\]
The unchanged actual-source identities and bounds are
\[
q(a)D_{r,k}
=\sum_{\gamma>0}\frac{y^r(1-y)^k}{a+\mu},
\tag{2}
\]
\[
800<a<900,\quad 0<q(a)<1/2,\quad
\sum_{\gamma>T}|a+\mu|^{-1}<1,\quad T=3\cdot10^{12}.
\tag{3}
\]
All sums converge absolutely. The low-height terms \(\gamma\le T\)
are positive real numbers.

The critical-line input below \(T\) is Platt and Trudgian,
*The Riemann hypothesis is true up to \(3\cdot10^{12}\)*,
Bulletin of the London Mathematical Society 53 (2021), 792–797:
[Theorem 1](https://arxiv.org/pdf/2004.09765),
[DOI 10.1112/blms.12460](https://doi.org/10.1112/blms.12460).
The high-zero budget and product identity are from S1 of
[the positive-sector note](RH_POSITIVE_LAYER_SECTORS.md).

C1 of [the quadratic phase-window note](RH_QUADRATIC_PHASE_WINDOW.md)
supplies the following consequence of the cited coarse Rosser estimate:
\[
X\ge1000\quad\Longrightarrow\quad
\text{some zero has }X<\gamma\le X+100.
\tag{4}
\]
The underlying counting estimate counts all zeros, not just
critical-line zeros. Its normalization and Rosser constants appear
in equation (1.1) and Table 1 of Hasanalizade, Shen and Wong,
*Counting zeros of the Riemann zeta function*, Journal of Number Theory
235 (2022), 219–241:
[published author-hosted paper](https://www-math.nsysu.edu.tw/~pjwong/stuff/CountingRiemannZeros.pdf),
[DOI 10.1016/j.jnt.2021.06.032](https://doi.org/10.1016/j.jnt.2021.06.032).
The historical source is Rosser, *Explicit bounds for some functions
of prime numbers*, American Journal of Mathematics 63 (1941), 211–232,
[DOI 10.2307/2371291](https://doi.org/10.2307/2371291).

Finally, O3–O4 of the preceding note retain
\[
D_{r,k}>0\quad\hbox{if}\quad
k\ge Q_*(r+1),\qquad Q_*(m)=100m^2\log(2m),
\tag{5}
\]
and positivity for every \(k\) in all rows \(r<10^{17}\).
These are written predecessor results, not new numerical computations.

## Theorem F1. A decreasing envelope for the full high-zero remainder

Let \(r\ge1\), \(k\ge0\), and
\[
b=k/r,\qquad G=\sqrt{ab}=\sqrt{\frac{ak}{r}}.
\]
Assume \(G\le T\). Set
\[
d_T=\sqrt{(a+T^2)^2-a},\qquad
E_{r,k}(T)=
\left(\frac a{d_T}\right)^r
\left(\frac{T^2}{d_T}\right)^k.
\tag{6}
\]
Then the full high-zero remainder obeys
\[
\boxed{
\left|
\sum_{\gamma>T}\frac{y^r(1-y)^k}{a+\mu}
\right|<E_{r,k}(T).
}
\tag{7}
\]

**Proof.** For a high zero, put \(x=|\mu|=\gamma^2+\delta^2>T^2\).
The exact identity
\[
|a+\mu|^2=(a+x)^2-4a\delta^2
\]
and \(|\delta|<1/2\) imply
\[
|y|^r|1-y|^k
\le
\left[
\frac{a x^b}{((a+x)^2-a)^{(1+b)/2}}
\right]^r.
\tag{8}
\]
Denote the bracket by \(A_b(x)\). Logarithmic differentiation gives
\[
\frac{d}{dx}\log A_b(x)
=\frac{b}{x}
-\frac{(1+b)(a+x)}{(a+x)^2-a}
=\frac{(a+x)(ab-x)-ab}
{x((a+x)^2-a)}.
\tag{9}
\]
Since \(ab=G^2\le T^2<x\), the last numerator is negative.
Thus \(A_b\) decreases on \(x\ge T^2\), including the case \(b=0\),
and
\[
|y|^r|1-y|^k<E_{r,k}(T).
\]
Multiply by \(|a+\mu|^{-1}\), sum over every high zero,
and use the complete source budget (3). This proves (7). QED.

This estimate retains the \(k\)-dependence of the high-zero tail.
It does not assume that individual high-zero phases are positive.

## Theorem F2. Exact comparison with a nearby verified zero

Assume
\[
r\ge1,\qquad 1000\le G=\sqrt{ak/r},\qquad h=G+100\le T.
\tag{10}
\]
Define
\[
P_{r,k}(h)=
\frac1{a+h^2}
\left(\frac a{a+h^2}\right)^r
\left(\frac{h^2}{a+h^2}\right)^k
\]
and
\[
\Delta_T(G)=
\log\frac{d_T}{a+h^2}
+\frac{G^2}{a}
\log\frac{h^2d_T}{T^2(a+h^2)}.
\tag{11}
\]
Then
\[
q(a)D_{r,k}>P_{r,k}(h)-E_{r,k}(T).
\tag{12}
\]
In particular,
\[
\boxed{
r\Delta_T(G)\ge\log(a+h^2)
\quad\Longrightarrow\quad D_{r,k}>0.
}
\tag{13}
\]

**Proof.** Equation (4) supplies a zero with
\(G<\gamma_\star\le G+100=h\le T\).
The finite-height theorem places it on the critical line.
For \(b=k/r\), let
\[
f_b(x)=\log\frac a{a+x}+b\log\frac x{a+x}.
\]
Then
\[
f_b'(x)=\frac{ab-x}{x(a+x)}.
\tag{14}
\]
The real kernel decreases for \(x\ge ab=G^2\).
Hence the selected zero contributes at least \(P_{r,k}(h)\)
to the low-height sum in (2). All other low terms are positive,
and F1 bounds the entire high remainder, proving (12).
A direct ratio gives
\[
\frac{P_{r,k}(h)}{E_{r,k}(T)}
=\frac{\exp(r\Delta_T(G))}{a+h^2}.
\]
Condition (13) makes this ratio at least one.
The strict remainder bound in (12) then proves positivity. QED.

No ordinate within the interval \((G,G+100]\) is numerically evaluated.

## Theorem F3. An explicit row-dependent verified-height frontier

Let
\[
\ell_T=\log(a+T^2),\qquad
\mathcal G_T(r)=
T-100-2\sqrt{1+\frac{T^2\ell_T}{r}}\quad(r\ge1).
\tag{15}
\]
For all integer pairs satisfying
\[
r\ge1,\qquad 1000\le G=\sqrt{ak/r}\le\mathcal G_T(r),
\tag{16}
\]
one has \(D_{r,k}>0\).

**Proof.** In the notation of F2, compare the real profile at \(h^2\)
and \(T^2\). Since \(h=G+100\ge G\), equation (14) gives
\[
f_b(h^2)-f_b(T^2)
=\int_{h^2}^{T^2}\frac{x-G^2}{x(a+x)}\,dx
\ge\frac{(T^2-h^2)^2}{2T^2(a+T^2)}
\ge\frac{(T-h)^2}{4T^2}.
\tag{17}
\]
Here \(a<T^2\), \(x-G^2\ge x-h^2\), and \(T+h\ge T\).
The correction between the complex envelope and the real profile is
\[
e_T(b)=
-\frac{1+b}{2}
\log\left(1-\frac a{(a+T^2)^2}\right).
\]
Using \(-\log(1-u)\le2u\), \(0<u<1/2\), and
\(b=G^2/a\le T^2/a\), one obtains
\[
0<e_T(b)\le\frac{(1+b)a}{(a+T^2)^2}
\le\frac1{a+T^2}<\frac1{T^2}.
\tag{18}
\]
Therefore
\[
\boxed{
\Delta_T(G)>
\frac{(T-G-100)^2-4}{4T^2}.
}
\tag{19}
\]
Condition (16) makes the right side at least \(\ell_T/r\).
Thus \(r\Delta_T(G)>\ell_T\ge\log(a+h^2)\), and F2 applies.
The same condition ensures \(h<T\), as required. QED.

The number 100 comes from the retained zero-count interval.
The additional margin in (15) comes from the envelope estimate.
These margins are sufficient bounds, not intrinsic spectral boundaries.

## Theorem F4. A positive sector reaching \(0.99T\)

For the actual source,
\[
\boxed{
D_{r,k}>0\quad\hbox{whenever}\quad
0\le k\le B_{99}r,\qquad
B_{99}=\frac{9801}{10000}\frac{T^2}{a}.
}
\tag{20}
\]
In particular \(B_{99}>9.801\cdot10^{21}\).

More strongly, for every \(r\ge10^7\),
\[
\boxed{
0\le k\le\frac r a\,\mathcal G_T(r)^2
\quad\Longrightarrow\quad D_{r,k}>0.
}
\tag{21}
\]

**Proof.** First, \(a+T^2<10^{25}\) and \(\log10<5/2\) give
\(\ell_T<63\). For \(r\ge10^7\),
\[
1+\frac{T^2\ell_T}{r}<9\cdot10^{-6}T^2,
\]
since \(T=3\cdot10^{12}\). Also \(100<10^{-3}T\).
Hence
\[
\mathcal G_T(r)>T-10^{-3}T-6\cdot10^{-3}T
=0.993T>0.99T.
\tag{22}
\]
For \(1000\le G\le\mathcal G_T(r)\), apply F3.
For \(0\le G<1000\), one has
\[
k/r=G^2/a<1250<T^2/(10a);
\]
the preceding O1 sector applies. This proves (21) for \(r\ge10^7\),
and therefore (20) on those rows by (22).

Rows \(r<10^7\) are already positive at every \(k\) by the retained
complete-row result. This includes the sector endpoint \(r=k=0\).
Finally \(a<900\) implies
\[
B_{99}>
\frac{9801}{10000}\frac{9\cdot10^{24}}{900}
=9.801\cdot10^{21}.
\]
QED.

No new height verification is being claimed: \(G\) is the peak height
of a comparison kernel, while the cited verified height remains \(T\).

## Theorem F5. Every difference order in rows below \(10^{18}\)

For the unchanged source,
\[
\boxed{
D_{r,k}>0\qquad 0\le r<10^{18},\quad k\ge0.
}
\tag{23}
\]

**Proof.** For \(1\le r<10^{18}\), let \(m=r+1\le10^{18}\).
Then
\[
\frac{Q_*(r+1)}r
=\frac{100m^2\log(2m)}{m-1}
\le200m\log(2m)<9.2\cdot10^{21}<B_{99}.
\tag{24}
\]
For the strict bound use
\(\log(2m)\le\log2+18\log10<46\).
Thus the F4 sector and the retained positive tail (5) overlap,
covering every \(k\) in these rows. Row \(r=0\) was already
completely positive. QED.

This counts scalar binomial rows covered by an analytic argument.
It is not a count of newly checked zeros, positive Hankel dimensions,
or numerical entries evaluated.

## Theorem F6. The remaining sign condition in height coordinates

Define
\[
\mathcal U_{\rm height}=
\left\{(r,k)\in\mathbb Z_{\ge0}^2:
\begin{array}{l}
r\ge10^{18},\\
\displaystyle\frac r a\,\mathcal G_T(r)^2
<k<100(r+1)^2\log(2r+2)
\end{array}
\right\}.
\tag{25}
\]
Then the retained all-pairs criterion ILC-T2 gives
\[
\boxed{
\mathrm{RH}\ \Longleftrightarrow\
D_{r,k}\ge0\quad\hbox{for every }(r,k)\in\mathcal U_{\rm height}.
}
\tag{26}
\]

**Proof.** F5 covers every lower row. In each remaining row,
F4 covers the lower side and (5) covers the upper side.
Their exact complement is (25). Apply ILC-T2. QED.

In terms of peak height \(G=\sqrt{ak/r}\), the remaining band is
\[
\mathcal G_T(r)<G<
\sqrt{\frac{100a(r+1)^2\log(2r+2)}r},
\qquad r\ge10^{18}.
\tag{27}
\]
The lower sufficient frontier has the limit
\[
\lim_{r\to\infty}\mathcal G_T(r)=T-102,
\]
where 102 reflects the coarse margins in this proof. The upper edge
in (27) grows without bound. Consequently the fixed verified-height
input cannot, by this estimate alone, settle the unbounded height band.
One must establish additional actual-source sign control there, or
supply stronger justified input. Merely substituting a larger unverified
\(T\) would not be a valid step.

The remaining band is not asserted to contain a negative entry.
Its nonnegativity at every pair is still open. No complete RH proof,
all-order Hankel positivity, or global eta improvement is concluded.

## Evidence and execution

| Item | Evidence |
|---|---|
| Entire high-zero tail with its \(k\)-dependence | Exact modulus and monotonicity F1 |
| Positive contribution near the comparison peak | Cited count interval plus verified-height theorem, F2 |
| Row-dependent height frontier | Integral profile comparison and correction bound, F3 |
| Uniform sector through peak height \(0.99T\) | F4 |
| All \(k\) in rows \(r<10^{18}\) | Overlap F5 |
| Nonnegativity throughout the unbounded remaining height band | Open |
| New numerical computation or machine certificate | None |

Only this note and the root README entry are changed.
These are written mathematical proofs, not newly machine-checked
certificates. No local scientific computation, independent agent,
zero verification replay, or new scientific CI stage is used.
Engines, numerical certificates, and workflows are unchanged.
The documentation commit uses [skip ci]. Classical results retain
their attribution; no optimality or literature-wide novelty is claimed.
