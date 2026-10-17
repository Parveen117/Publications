# RH continuation: finite observer detection and one Hankel tower

17 October 2026. Continues [the verified-height frontier](RH_VERIFIED_HEIGHT_FRONTIER.md)
at Publications commit `76d09efd30ec4c0727ebcbb74ad2261feba5440c`.

This step addresses the sign information missing from zero counts.
For the actual discrete xi source, any hypothetical off-line zero pair
produces a negative determinant in a finite two-dimensional compression
of the source moment form. The compression uses polynomials of potentially
large degree; it is not a fixed two-layer RH test.

The construction gives an explicit sufficient depth once the pair and
its separation data are specified. It also proves that one unshifted
Hankel tower suffices as an all-order RH criterion. A separate comparison
construction shows that the retained counting, symmetry, finite-height,
and anchor-budget inputs alone permit both positive and indefinite
sources. Additional information about the actual xi source is therefore
necessary to close the remaining sign condition.

No off-line zeta zero, negative actual determinant, or negative actual
binomial reading is exhibited here. No new positive index region is claimed.

## Actual source and the polynomial pairing

Retain
\[
a=3249/4,\quad \mathcal X(w^2)=\xi(1/2+w),\quad
q=\mathcal X'/\mathcal X,\quad
c_n=\frac{(-a)^nq^{(n)}(a)}{n!\,q(a)}.
\]
Let \(S\) be the set of distinct values
\[
u=\frac a{a+(\gamma-i\delta)^2},
\qquad \delta=\beta-\tfrac12,
\]
from nontrivial zeros \(\beta+i\gamma\) with \(\gamma>0\).
Write \(\nu_u\) for the multiplicity and
\[
w_u=\frac{\nu_u u}{a q(a)}.
\]
The unconditional product identity of S1 in
[the positive-sector note](RH_POSITIVE_LAYER_SECTORS.md) gives
\[
c_n=\sum_{u\in S}w_u u^n.
\tag{1}
\]
This is an absolutely convergent complex-weight representation; positivity
of its weights is not assumed.

The set \(S\) is closed under conjugation, satisfies \(0<|u|<1\),
and has no accumulation point except 0. Every set
\(\{u\in S:|u|\ge s\}\), \(s>0\), is finite.
Indeed the zero heights tend to infinity and
\(u=a/\gamma^2+O(\gamma^{-3})\). The finite-height theorem places
all low zeros on the critical line, while every high zero has
\(\Re(a+\mu)>0\) and \(|a+\mu|>a\).

The retained source estimates also imply
\[
\sum_{u\in S}|w_u|\le2,\qquad c_0=1.
\tag{2}
\]
For completeness, low-zero terms in \(q(a)=\sum(a+\mu)^{-1}\)
are positive. At high zeros, S1 proves
\(|(a+\mu)^{-1}|<2\Re(a+\mu)^{-1}\).
Sum these inequalities and divide by \(q(a)>0\) to obtain (2).

For a real polynomial \(p\), define
\[
\mathscr L[p]=\sum_{u\in S}w_u p(u),\qquad
\mathscr L[t^n]=c_n.
\tag{3}
\]
Conjugation symmetry makes this functional real-valued.
It does not yet make \(\mathscr L[p^2]\) nonnegative.

Use the same size convention as ILC-T1:
\[
H_N=(c_{i+j})_{i,j=0}^{N-1},\qquad N\ge1.
\tag{4}
\]
For a real coefficient vector \(v\) representing \(p\),
\(v^\mathsf{T}H_Nv=\mathscr L[p^2]\).

The canonical-product background remains attributed to
[NIST DLMF 25.2(iv)](https://dlmf.nist.gov/25.2#iv).
The existing layer and binomial criteria are in
[ILC-T1–T2](papers/native-infinite-layer-cut/THEOREM.md).

## Theorem D1. A nonreal pair has an indefinite two-dimensional block

For \(z\notin\mathbb R\) and \(A\ne0\), define the real symmetric matrix
\[
P(A,z)=
\begin{pmatrix}
2\Re A&2\Re(Az)\\
2\Re(Az)&2\Re(Az^2)
\end{pmatrix}.
\]
Then
\[
\boxed{\det P(A,z)=-4|A|^2(\Im z)^2<0.}
\tag{5}
\]

**Proof.** Factor the bilinear moment matrix as
\[
P(A,z)=
\begin{pmatrix}1&1\\z&\bar z\end{pmatrix}
\begin{pmatrix}A&0\\0&\bar A\end{pmatrix}
\begin{pmatrix}1&z\\1&\bar z\end{pmatrix}.
\]
Its determinant is
\(|A|^2(\bar z-z)^2=-4|A|^2(\Im z)^2\). QED.

Thus conjugate pairing gives a real matrix but leaves one positive
and one negative eigenvalue. Under any real invertible observer-basis
change \(R\),
\[
\det(R^\mathsf{T}PR)=(\det R)^2\det P<0.
\tag{6}
\]
The negative determinant cannot be removed by such a change of basis.

## Theorem D2. Finite isolation of a hypothetical off-line pair

Suppose the actual source has \(z\in S\setminus\mathbb R\).
Let
\[
s=|z|,\qquad \eta=|\Im z|,\qquad \nu=\nu_z=\nu_{\bar z},
\]
and define the finite real polynomial
\[
Q(t)=
\prod_{\substack{u\in S,\ |u|\ge s\\u\ne z,\bar z}}(t-u).
\tag{7}
\]
Each distinct value occurs once; a single vanishing polynomial factor
annihilates its contribution regardless of its multiplicity.
The empty product is 1. Let \(d=\deg Q\).

After these outer values and the selected pair are removed, the remaining
support is nonempty and has a largest modulus \(s_0\) with \(0<s_0<s\).
This follows from the sole accumulation point being 0.
Write \(Q(t)=\sum_j q_jt^j\), and set
\[
C=2\left(\sum_j|q_j|\right)^2,\qquad
\kappa=\frac{\nu s|Q(z)|^2}{a q(a)}>0.
\tag{8}
\]
For \(n\ge0\), consider the two real observer polynomials
\[
p_{n,0}(t)=t^nQ(t),\qquad p_{n,1}(t)=t^{n+1}Q(t),
\]
and their full-source matrix
\[
\mathcal H_n=
\left(\mathscr L[p_{n,i}p_{n,j}]\right)_{i,j=0}^{1}.
\tag{9}
\]
If
\[
\boxed{
\frac C\kappa\left(\frac{s_0}{s}\right)^{2n}
\le\frac{\eta^2}{4},
}
\tag{10}
\]
then \(\det\mathcal H_n<0\). In particular it suffices to take
\[
\boxed{
n\ge
\max\left(
0,\left\lceil
\frac{\log(4C/(\kappa\eta^2))}
{2\log(s/s_0)}
\right\rceil
\right).
}
\tag{11}
\]
Thus a finite moment matrix \(H_{n+d+2}\) is not positive semidefinite.

**Proof.** The polynomial \(Q\) has real coefficients because its
root set is closed under conjugation. It is nonzero at \(z,\bar z\)
and kills every other support value of modulus at least \(s\).
The selected pair contributes \(P(A_n,z)\), where
\[
A_n=\frac{\nu z^{2n+1}Q(z)^2}{a q(a)},\qquad
|A_n|=\kappa s^{2n}.
\]
All other surviving terms have modulus at most \(s_0\).
Equation (2), \(|u|<1\), and
\(|Q(u)|\le\sum_j|q_j|\) give
\[
\mathcal H_n=P(A_n,z)+R_n,\qquad
|(R_n)_{ij}|\le C s_0^{2n}\quad(0\le i,j\le1).
\tag{12}
\]
Also \(|P(A_n,z)_{ij}|\le2\kappa s^{2n}\).
Expanding the determinant and using (5) yields
\[
\det\mathcal H_n
\le-4\kappa^2s^{4n}\eta^2
+8\kappa C s^{2n}s_0^{2n}
+2C^2s_0^{4n}.
\tag{13}
\]
Let \(v_n=(C/\kappa)(s_0/s)^{2n}\).
After dividing by \(\kappa^2s^{4n}>0\), the right side is
\[
-4\eta^2+8v_n+2v_n^2.
\]
Under (10) it is at most
\(-2\eta^2+\eta^4/8<0\), since \(0<\eta<1\).
This proves the determinant claim; taking logarithms gives (11).

The two polynomials in (9) have degree at most \(n+d+1\).
Their matrix is a real congruence compression of \(H_{n+d+2}\).
A positive semidefinite matrix cannot have a compression with negative
determinant, so this finite source moment matrix fails positivity. QED.

The depth bound depends on the hypothetical pair, the annihilator,
and the separation ratio \(s_0/s\). No actual values of these data
are supplied. The argument neither locates an off-line zero nor
provides a universal fixed matrix size that decides RH.

## Theorem D3. One unshifted Hankel tower suffices for the actual source

The following conditions are equivalent:
\[
\boxed{
\begin{array}{ll}
\text{(i)}&\mathrm{RH},\\
\text{(ii)}&H_N\succ0\quad\text{for every }N\ge1,\\
\text{(iii)}&H_N\succeq0\quad\text{for every }N\ge1,\\
\text{(iv)}&\det H_N>0\quad\text{for every }N\ge1.
\end{array}
}
\tag{14}
\]

**Proof.** Under RH, every support value \(u\) in (1) lies in
\((0,1)\) and every weight \(w_u\) is positive.
There are infinitely many distinct support values.
A nonzero polynomial has only finitely many roots, so
\[
\mathscr L[p^2]=\sum_{u\in S}w_up(u)^2>0
\]
for every nonzero real polynomial. This proves (i) implies (ii);
(ii) implies (iii) directly.

If RH fails, a zero with \(\delta\ne0\) gives
\[
\Im\frac a{a+(\gamma-i\delta)^2}
=\frac{2a\gamma\delta}{|a+(\gamma-i\delta)^2|^2}\ne0.
\]
Theorem D2 then supplies a finite \(H_N\) that is not positive
semidefinite. This proves (iii) implies (i).

Finally, the leading principal minors of \(H_N\) are exactly
\(\det H_1,\ldots,\det H_N\). Sylvester's criterion, equivalently
successive completion of squares, makes (ii) and (iv) equivalent.
QED.

This supplements ILC-T2: for this actual discrete source, one
unshifted Hankel tower is enough for the logical all-order criterion.
The two shifted forms of ILC still provide their finite positive
contraction construction. This theorem does not replace those
finite hypotheses by positivity of a single small \(H_N\).

## Theorem D4. A finite negative moment form gives a finite binomial witness

Suppose a real polynomial \(f\) satisfies \(\mathscr L[f^2]<0\).
Let
\[
\epsilon=-\tfrac12\mathscr L[f^2]>0,\qquad
g(t)=f(t)^2+\epsilon=\sum_{\ell=0}^{D}a_\ell t^\ell,
\]
and
\[
A_g=\sum_{\ell=0}^{D}|a_\ell|\ell(\ell-1).
\tag{15}
\]
For any integer
\[
N\ge\max(D,1),\qquad N>A_g/\epsilon,
\tag{16}
\]
at least one index \(0\le j\le N\) satisfies
\[
\boxed{D_{j,N-j}<0.}
\tag{17}
\]

**Proof.** Since \(\mathscr L[1]=1\),
\(\mathscr L[g]=-\epsilon\), while \(g(t)\ge\epsilon\)
for \(0\le t\le1\).
Expand \(g\) in the Bernstein basis of degree \(N\):
\[
g(t)=\sum_{j=0}^{N}b_{N,j}\binom Nj t^j(1-t)^{N-j},
\]
\[
b_{N,j}=
\sum_{\ell=0}^{D}a_\ell\frac{(j)_\ell}{(N)_\ell},
\tag{18}
\]
with the falling-factorial ratio equal to 1 for \(\ell=0\).
The identity follows from the factorial moments of the binomial
polynomial; it is the degree-elevation identity used in ILC-T2.

For \(0\le j\le N\), compare sampling \(\ell\) labels from \(N\)
with and without replacement. The probabilities that all sampled
labels lie in a fixed \(j\)-element set are respectively
\((j/N)^\ell\) and \((j)_\ell/(N)_\ell\).
Couple the samples until a repeated label occurs. The probability
of a repeat is at most \(\ell(\ell-1)/(2N)\), by the union bound.
Hence
\[
\left|\frac{(j)_\ell}{(N)_\ell}-(j/N)^\ell\right|
\le\frac{\ell(\ell-1)}{2N}.
\]
Equations (15)–(16) imply
\[
|b_{N,j}-g(j/N)|\le A_g/(2N)<\epsilon/2.
\]
Thus every \(b_{N,j}>\epsilon/2>0\).
Apply \(\mathscr L\) to (18):
\[
-\epsilon
=\sum_{j=0}^{N}b_{N,j}\binom Nj D_{j,N-j}.
\]
At least one of these binomial readings must be negative. QED.

A negative eigenvector of the real matrix in D2 supplies the
polynomial \(f=t^nQ(t)(\alpha+\beta t)\) required here.
Thus D2 and D4 give a finite witness construction from a hypothetical
off-line pair, with explicit sufficient bounds at each step.

For the actual source, every such witness would have to lie in the
unsettled height band of F6 in the preceding note, since all other
pairs have already been proved positive. No such actual witness has
been found or evaluated in this continuation.

## Theorem D5. Height counts and the retained structural bounds do not force positivity

There exist two comparison sources with the same positive-height
zero-count function, the same finite-height critical-line data below
\(T=3\cdot10^{12}\), conjugation and reflection symmetry, genus-zero
products, and \(0<q(a)<1/2\), such that:

- the first source has every finite \(H_N\) positive definite;
- the second source has a finite indefinite moment matrix and a
  negative binomial reading.

Both satisfy the same coarse zero-count estimate used in the preceding
positive-region proofs. These comparison sources are not asserted to
be the actual completed zeta function.

**Construction and proof.** List all actual positive zero heights
\(\gamma\), with multiplicity, and choose one height \(\gamma_0>T\).
Project each listed zero onto the critical line and add one further
critical-line zero at height \(\gamma_0\). In the \(t\)-variable define
\[
\mathcal X_+(t)=
\left[\prod_{\gamma>0}(1+t/\gamma^2)\right](1+t/\gamma_0^2).
\tag{19}
\]
At least two copies of the factor at \(\gamma_0\) are present.
Fix \(\delta_0=1/4\), set
\(\mu_0=(\gamma_0-i\delta_0)^2\), and replace exactly two copies:
\[
\mathcal X_-(t)=\mathcal X_+(t)
\frac{(1+t/\mu_0)(1+t/\bar\mu_0)}
{(1+t/\gamma_0^2)^2}.
\tag{20}
\]
The denominator cancels existing factors. Both are entire genus-zero
products. The functions
\(\mathcal X_\pm((s-1/2)^2)\) have conjugation symmetry and are
invariant under \(s\mapsto1-s\).
The replacement in (20) turns two upper zeros at \(1/2+i\gamma_0\)
into the pair \(1/2\pm1/4+i\gamma_0\), and makes the corresponding
lower-half-plane replacement.

The count functions are exactly equal to each other:
\[
N_+(x)=N_-(x)=N_{\xi}(x)+\mathbf1_{\{x\ge\gamma_0\}}.
\tag{21}
\]
Below \(T\), all zero data remain the actual critical-line data.
The products converge because \(\sum\gamma^{-2}<\infty\).
The same zero-count growth gives order \(1/2\) in \(t\):
the bound on \(\sum\log(1+|t|/\gamma^2)\) gives the upper order,
and the zero count gives the matching lower order.

For \(x\ge e\), Rosser's bound for the actual count, with the extra
unit in (21), gives
\[
|N_\pm(x)-F(x)|
\le0.137\log x+0.443\log\log x+3.463
<\log x+10,
\]
where \(F(x)=x\log(x/(2\pi e))/(2\pi)\).
Thus both retain the coarse counting estimate used in C1.

It remains to check the anchor budget rather than assume it.
For an actual high zero put \(S_0=a+\gamma^2\), \(d=\delta^2\).
A direct comparison gives
\[
\frac{1/(a+\gamma^2)}{\Re(1/(a+\mu))}
=
1+\frac{d(3\gamma^2-a+d)}
{S_0(S_0-d)}
\le1+\frac1{\gamma^2}.
\tag{22}
\]
Here \(0\le d<1/4\), \(a>1/4\), and \(\gamma>T\);
the denominator exceeds \(\gamma^4\) and the numerator is
less than \(\gamma^2\) when \(d>0\). For \(d=0\) the ratio is one.
Low actual zeros are already on the line.
Using the retained bound \(q_\xi(a)<49/114\) from SE49, (22) yields
\[
0<q_+(a)\le(1+T^{-2})q_\xi(a)+T^{-2}
<(1+T^{-2})\frac{49}{114}+T^{-2}<\frac12.
\tag{23}
\]
At the replaced height, \(\Re(1/(a+\mu_0))<
1/(a+\gamma_0^2)\), since
\(\delta_0^2(3\gamma_0^2-a+\delta_0^2)>0\).
Every real part remains positive, so
\(0<q_-(a)<q_+(a)<1/2\).
The same elementary source-budget argument as (2) applies to both.

All support values for \(\mathcal X_+\) are in \((0,1)\), with
positive weights and infinite distinct support. Its moment matrices
are positive definite. The source \(\mathcal X_-\) contains the
nonreal pair from (20), so D2 and D4 give its finite negative witnesses.
This completes the construction. QED.

The original Rosser constants and counting normalization are reproduced
in equation (1.1) and Table 1 of Hasanalizade, Shen and Wong,
*Counting zeros of the Riemann zeta function*, Journal of Number Theory
235 (2022), 219–241:
[published author-hosted paper](https://www-math.nsysu.edu.tw/~pjwong/stuff/CountingRiemannZeros.pdf),
[DOI 10.1016/j.jnt.2021.06.032](https://doi.org/10.1016/j.jnt.2021.06.032).
The historical source is Rosser, American Journal of Mathematics
63 (1941), 211–232,
[DOI 10.2307/2371291](https://doi.org/10.2307/2371291).
The anchor estimate used in (23) is retained in
[SE49](papers/native-moving-zero-covariance/sources/SE49.md)
and [the positive-sector note](RH_POSITIVE_LAYER_SECTORS.md).

The comparison models do not retain the actual Euler/Gamma identity
or the actual source's full derivative sequence. They prove a limitation
of the listed structural inputs, not a failure of RH for zeta.
They are mathematical constructions in this note; no repository source
engine or actual source data are replaced.

## Consequence for the next actual-source step

The finite-height frontier and the positive tails remain valid.
This continuation adds an exact detector and a sharper logical target:
\[
\boxed{
\mathscr L[p^2]>0
\quad\text{for every nonzero finite real polynomial }p,
}
\tag{24}
\]
equivalently the single determinant tower in (14).
For the actual source, establishing (24) from its complete arithmetic
formula would prove RH. Equations (19)–(23) explain why height counts,
reflection symmetry, and the anchor budget alone cannot establish it.

A new proof must therefore use additional actual-source structure,
such as a justified positive representation of the combined arithmetic
source, without assuming the positive zero measure that RH would supply.
No such all-order positivity proof is supplied in this note.

## Evidence and execution

| Result | Evidence |
|---|---|
| Isolated nonreal pair gives a negative determinant | Exact identity D1 |
| Finite full-source detection of any hypothetical pair | Annihilator, complete tail bound, and depth estimate D2 |
| One unshifted Hankel tower is equivalent to RH | D3 |
| Finite negative binomial witness follows | Explicit Bernstein degree elevation D4 |
| Counts and the listed structural bounds are insufficient alone | Two comparison products with the same count, D5 |
| Actual off-line zero or actual negative witness | None found or asserted |
| Actual all-order Hankel positivity | Open |

Only this note and its README entry are added. The proofs are written
arguments, not newly machine-checked certificates. No local scientific
computation, independent agent, numerical replay, or new scientific CI
stage is used. Existing engines, certificates, and workflows remain
unchanged; the documentation commit uses [skip ci].
Polynomial moment methods, congruence invariance, Sylvester's criterion,
and Bernstein degree elevation retain their classical status. No
literature-wide novelty or optimal-depth claim is made.
