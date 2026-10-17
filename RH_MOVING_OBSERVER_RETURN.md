# MOR-1 — Moving observers return finite-anchor signs

17 October 2026. Continues NIF-1 and the finite observer detector at
Publications commit 4ed004e7dcfe991eeec1715d15ba0537e2c70087.

NIF-1 proved that every fixed finite source bank approaches a positive
arcsine bank as its anchor increases. That limit discards some source
information. This continuation retains it by transporting the observer
as well as the anchor.

For the actual xi source, a square-root weighted analytic observer
preserves every finite-anchor polynomial pairing exactly. A finite
polynomial approximation has a uniform, source-weighted error bound.
For a fixed target degree and error, its required degree grows at most
linearly with the anchor ratio. A hypothetical negative polynomial
witness therefore persists at all larger anchors in linearly growing
banks. Conversely, positivity along a sequence of sufficiently growing
banks is equivalent to RH.

These are written source-transport and error theorems. They do not
establish positivity of those growing actual banks or exhibit any
negative actual witness.

## Actual source and coordinates

Retain
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q=\mathcal X'/\mathcal X,\qquad a_0=3249/4.
\]
A nontrivial zero \(\beta+i\gamma\), \(\gamma>0\), gives
\[
\delta=\beta-\tfrac12,\qquad
\mu=(\gamma-i\delta)^2,\qquad u_\mu=\Re\mu.
\]
Sum over distinct values of \(\mu\), with positive integer multiplicity
\(\nu_\mu\). The set is closed under conjugation. The retained
unconditional genus-zero product gives
\[
q(a)=\sum_\mu\frac{\nu_\mu}{a+\mu}.
\tag{1}
\]
For a real anchor \(a\ge a_0\), put
\[
\rho=\frac a{a_0}\ge1,\qquad r=1-\rho^{-1},\qquad
f_\rho=\frac{q(a)}{q(a_0)},\qquad
y_a(\mu)=\frac a{a+\mu}.
\tag{2}
\]
The source moments and their unshifted size-\(N\) bank are
\[
c_n(a)=\frac{(-a)^nq^{(n)}(a)}{n!\,q(a)},\qquad
H_N(a)=(c_{i+j}(a))_{0\le i,j<N}.
\tag{3}
\]
The polynomial pairing is
\[
\mathcal L_a[p]=
\frac1{q(a)}\sum_\mu\frac{\nu_\mu}{a+\mu}\,
p(y_a(\mu)),\qquad \mathcal L_a[y^n]=c_n(a).
\tag{4}
\]
For functions holomorphic on a neighborhood of the closed unit disk,
the same absolutely convergent sum defines \(\mathcal L_a\).
Real Taylor coefficients give real readings by conjugation symmetry.

This is a bilinear source functional. No positive Hilbert realization,
positive actual zero measure, or all-order moment positivity is assumed.

The retained product and finite detector are in
[the positive-sector note](RH_POSITIVE_LAYER_SECTORS.md) and
[Theorems D1-D4](RH_FINITE_OBSERVER_DETECTOR.md).
The anchor convention is [LAT-1](papers/native-layer-anchor-transport/THEOREM.md).
The fixed-depth limit is [NIF-1](RH_INFINITY_CUT_FIXED_POINT.md).

## MOR-T1. A uniform absolute source budget survives every outward anchor move

For every real \(a\ge a_0\),
\[
q(a)>0,\qquad 0<|y_a(\mu)|<1,
\]
and
\[
\boxed{
\mathcal Q(a):=\sum_\mu\frac{\nu_\mu}{a+u_\mu}\le2q(a),\qquad
\frac1{q(a)}\sum_\mu\frac{\nu_\mu}{|a+\mu|}\le2.
}
\tag{5}
\]
In particular \(f_\rho>0\). The auxiliary quantity \(\mathcal Q\)
is only an absolute-error majorant; its real support does not replace
the actual complex source in (1)-(4).

**Proof.** Retain the published Platt--Trudgian theorem below
\(T=3\cdot10^{12}\), already used in the source notes.
A zero with \(0<\gamma\le T\) has \(\delta=0\), so
\(u_\mu=\gamma^2>0\) and \(\mu\) is real.
Above \(T\), the critical-strip bound \(|\delta|<1/2\) gives
\[
u_\mu=\gamma^2-\delta^2>0,\qquad
|\Im\mu|=2\gamma|\delta|<\gamma.
\]
Since \(a\ge a_0>1/4\),
\(a+u_\mu>\gamma^2\) for those high zeros.
Thus \(|\Im\mu|/(a+u_\mu)<1/\gamma<1\).

Writing \(v_\mu=\Im\mu\), every real part
\[
\Re\frac1{a+\mu}
=\frac{a+u_\mu}{(a+u_\mu)^2+v_\mu^2}
\]
is positive. For low real terms the following comparison is equality;
for the other terms it gives
\[
\frac1{a+u_\mu}\le2\Re\frac1{a+\mu}.
\]
The sums converge absolutely because \(u_\mu\) grows quadratically
with \(\gamma\) and the retained zero count gives
\(\sum_\mu\nu_\mu/\gamma^2<\infty\).
Conjugation makes (1) real, and summing proves its positivity and the
first bound in (5). Also
\[
\frac1{|a+\mu|}\le\frac1{a+u_\mu},\qquad
|y_a(\mu)|\le\frac a{a+u_\mu}<1.
\]
These prove the remaining claims. QED.

No new zero-height computation is used. The finite-height input is
Platt and Trudgian,
[The Riemann hypothesis is true up to 3 times 10^12](https://arxiv.org/abs/2004.09765).
The canonical-product background remains attributed to
[DLMF 25.2(iv)](https://dlmf.nist.gov/25.2#iv).

## MOR-T2. An explicit moving observer preserves the old pairing exactly

For a real polynomial \(p\), define
\[
\boxed{
(\mathcal O_\rho p)(y)=
(1-r y)^{-1/2}
p\!\left(\frac{y}{\rho(1-r y)}\right).
}
\tag{6}
\]
For \(\rho>1\) use the analytic square root equal to one at \(y=0\).
Its singularity is outside the closed unit disk, at
\(y=1/r=\rho/(\rho-1)>1\). For \(\rho=1\), set
\(\mathcal O_1p=p\).

For any two real polynomials \(p,h\),
\[
\boxed{
f_\rho\,\mathcal L_a[
(\mathcal O_\rho p)(\mathcal O_\rho h)]
=\mathcal L_{a_0}[ph].
}
\tag{7}
\]
Consequently, an observer basis carried by
\(\sqrt{f_\rho}\mathcal O_\rho\) has exactly its old real Gram matrix,
including all positive and negative directions.

The old moments also have the exact return series
\[
\boxed{
c_n(a_0)=
f_\rho\rho^{-n}\sum_{\ell=0}^\infty
\binom{n+\ell}{n}r^\ell c_{n+\ell}(a),
\qquad n\ge0.
}
\tag{8}
\]

**Proof.** At each actual source point, direct algebra gives
\[
1-r y_a=\frac{a_0+\mu}{a+\mu},\qquad
\frac{y_a}{\rho(1-r y_a)}=\frac{a_0}{a_0+\mu}=y_{a_0}.
\tag{9}
\]
The source weights obey
\[
\frac{f_\rho}{q(a)(a+\mu)(1-r y_a)}
=\frac1{q(a_0)(a_0+\mu)}.
\tag{10}
\]
Multiplying (6) for \(p,h\), applying (9)-(10) term by term, and
summing proves (7). The root has a single analytic branch on a
neighborhood of the closed unit disk; its square is exactly
\((1-r y)^{-1}\), so no conjugation or hidden absolute square has
been inserted into the bilinear source pairing.

For (8), apply the same weight relation to \(y_{a_0}^n\):
\[
c_n(a_0)=f_\rho\rho^{-n}
\mathcal L_a[y^n(1-r y)^{-n-1}].
\]
Expand by the binomial series. At each fixed \(\rho\), \(r<1\);
the expansion converges uniformly and absolutely on the closed
unit disk. Equation (5) therefore justifies interchange with the
source sum, yielding (8). QED.

Under an all-order positive realization, (7) becomes an isometry of
the corresponding polynomial pairings. Before that premise is proved,
it is the exact real bilinear identity just stated. In particular it
does not turn an indefinite source into a positive one.

The observer's singularity approaches \(y=1\) at distance
\(1/(\rho-1)\) as the anchor increases. Holding its polynomial degree
fixed would discard precisely the changing resolution needed in (7).

## MOR-T3. Finite moment return has an error controlled by depth divided by anchor ratio

For \(n,K\ge0\), define the finite return
\[
R_{n,K}(a)=f_\rho\rho^{-n}
\sum_{\ell=0}^K\binom{n+\ell}{n}r^\ell c_{n+\ell}(a).
\]
Then
\[
\boxed{
|c_n(a_0)-R_{n,K}(a)|
\le 2^{n+2}\exp\!\left(-\frac{K+1}{2\rho}\right).
}
\tag{11}
\]
For \(\rho=1\) the return is exact. For any fixed \(n\) and
\(0<\epsilon<1\), the sufficient condition
\[
K+1\ge
2\rho\big((n+2)\log2+\log(1/\epsilon)\big)
\tag{12}
\]
makes its error at most \(\epsilon\). Thus fixed-order, fixed-accuracy
return requires at most linear growth in \(\rho\) of the number of
additional retained layers. The constants are not claimed optimal.

**Proof.** We first derive the elementary tail inequality used here
and in MOR-T4. For \(\alpha>0\), \(0\le t<1\), and \(K\ge0\), let
\[
E_{\alpha,K}(t)=\sum_{\ell>K}\frac{(\alpha)_\ell}{\ell!}t^\ell.
\]
When \(t>0\), positivity of the coefficients gives
\[
E_{\alpha,K}(t)
\le t^{(K+1)/2}\sum_{\ell=0}^\infty
\frac{(\alpha)_\ell}{\ell!}(\sqrt t)^\ell
=t^{(K+1)/2}(1-\sqrt t)^{-\alpha}.
\]
Multiplying by \((1-t)^\alpha\) yields
\[
\boxed{
(1-t)^\alpha E_{\alpha,K}(t)
\le (1+\sqrt t)^\alpha t^{(K+1)/2}
\le 2^\alpha t^{(K+1)/2}.
}
\tag{13}
\]
At \(t=0\), both tails vanish and the same inequality holds.

In the absolute source remainder of (8), replace
\(|a+\mu|\) by the smaller positive denominator \(a+u_\mu\).
Set \(t_\mu=(a-a_0)/(a+u_\mu)\le r\).
With \(\alpha=n+1\), the bound (13) leaves the factor
\[
\frac{a_0^n}{q(a_0)(a+u_\mu)^{n+1}}
(1-t_\mu)^{-n-1}
=\frac{a_0^n}{q(a_0)(a_0+u_\mu)^{n+1}}.
\]
Since \(a_0/(a_0+u_\mu)\le1\), summing this factor with multiplicity
is at most \(\mathcal Q(a_0)/q(a_0)\le2\).
For \(\rho>1\), \(t_\mu\le r=1-1/\rho\), and hence
\(t_\mu^{(K+1)/2}\le\exp(-(K+1)/(2\rho))\).
Combining these estimates proves (11); taking logarithms proves (12).
The case \(\rho=1\) has \(r=0\) in (8) and is exact. QED.

The binomial expansion itself is classical
[DLMF 4.6.7](https://dlmf.nist.gov/4.6#E7).
The source-weighted bound is derived above; no numerical replay or
assumed positive actual spectral measure enters it.

## MOR-T4. Finite polynomial observers return a whole Gram matrix with explicit error

Write a real polynomial of degree at most \(d\) as
\(p(y)=\sum_{j=0}^d p_jy^j\), and put
\[
B(p)=\sum_{j=0}^d|p_j|.
\]
For \(K\ge0\), define the finite transported polynomial
\[
\boxed{
P_{\rho,K}[p](y)=
\sum_{j=0}^d p_j\rho^{-j}y^j
\sum_{\ell=0}^K
\frac{(j+1/2)_\ell}{\ell!}(r y)^\ell.
}
\tag{14}
\]
It has degree at most \(d+K\). Each binomial factor is truncated
after order \(K\); this is the specified polynomial, not an implicit
claim about a different total-degree truncation.

For real polynomials \(p,h\) of degree at most \(d\),
\[
\boxed{
\left|\mathcal L_{a_0}[ph]-
f_\rho\mathcal L_a[P_{\rho,K}[p]P_{\rho,K}[h]]\right|
\le
2^{d+3}B(p)B(h)
\exp\!\left(-\frac{K+1}{2\rho}\right).
}
\tag{15}
\]

In particular, for the size-\(N\) monomial observer define the
\((K+N)\)-by-\(N\) real coefficient matrix \(T_{\rho,K,N}\) by
\[
(T_{\rho,K,N})_{j+\ell,j}
=\rho^{-j}\frac{(j+1/2)_\ell}{\ell!}r^\ell,
\quad 0\le j<N,\quad 0\le\ell\le K,
\tag{16}
\]
and set all other entries to zero. Define the returned matrix
\[
G_{\rho,K,N}
=f_\rho\,T_{\rho,K,N}^{\mathsf T}
H_{K+N}(a)T_{\rho,K,N}.
\tag{17}
\]
Then, in the fixed ordinary coefficient norm at \(a_0\),
\[
\boxed{
\|H_N(a_0)-G_{\rho,K,N}\|_2
\le \varepsilon_{\rho,K,N}
:=2^{N+2}N\exp\!\left(-\frac{K+1}{2\rho}\right).
}
\tag{18}
\]
Equivalently,
\[
G_{\rho,K,N}-\varepsilon_{\rho,K,N}I
\preceq H_N(a_0)\preceq
G_{\rho,K,N}+\varepsilon_{\rho,K,N}I.
\tag{19}
\]

**Proof.** For \(0\le t\le1\), introduce the positive majorant
\[
\overline O_{\rho,p}(t)=
(1-rt)^{-1/2}
\sum_{j=0}^d |p_j|
\left(\frac{t}{\rho(1-rt)}\right)^j.
\]
Set
\[
\eta_{\rho,K,d}
=2^{d+1/2}\exp\!\left(-\frac{K+1}{2\rho}\right).
\]
For \(|y|\le1\), apply (13) with \(\alpha=j+1/2\) and
\(t=r|y|\le r\) to each term in (6). It gives
\[
|\mathcal O_\rho p-P_{\rho,K}[p]|
\le \eta_{\rho,K,d}\,\overline O_{\rho,p}(|y|).
\]
Positivity of the binomial coefficients also gives
\[
|\mathcal O_\rho p|,\ |P_{\rho,K}[p]|
\le\overline O_{\rho,p}(|y|).
\]
Consequently the absolute difference of the two products for \(p,h\)
is at most
\(2\eta_{\rho,K,d}\overline O_{\rho,p}(|y|)
\overline O_{\rho,h}(|y|)\).

At an actual source point, let \(z_\mu=a/(a+u_\mu)\).
Then \(|y_a|\le z_\mu<1\), the majorants are increasing, and
\[
1-rz_\mu=\frac{a_0+u_\mu}{a+u_\mu},\qquad
\frac{z_\mu}{\rho(1-rz_\mu)}
=\frac{a_0}{a_0+u_\mu}\le1.
\]
Therefore
\[
\overline O_{\rho,p}(|y_a|)
\overline O_{\rho,h}(|y_a|)
\le \frac{a+u_\mu}{a_0+u_\mu}B(p)B(h).
\]
The full source-weighted majorant is bounded by
\[
\frac{f_\rho}{q(a)}
\sum_\mu\frac{\nu_\mu}{|a+\mu|}
\overline O_{\rho,p}(|y_a|)
\overline O_{\rho,h}(|y_a|)
\le\frac{\mathcal Q(a_0)}{q(a_0)}B(p)B(h)
\le2B(p)B(h).
\]
Using the exact identity (7) now bounds the product error by
\(4\eta_{\rho,K,d}B(p)B(h)\).
Since \(4\cdot2^{d+1/2}<2^{d+3}\), equation (15) follows.
At \(\rho=1\), the approximation is exact and the same upper bound
is valid.

The columns in (16) are exactly the coefficients of (14) for
\(p=y^j\), proving (17). Apply (15) to each pair of these basis
polynomials with \(d=N-1\) and \(B=1\).
Every matrix entry error is at most
\(2^{N+2}\exp(-(K+1)/(2\rho))\).
The spectral norm of a real symmetric \(N\)-by-\(N\) matrix is at
most \(N\) times its largest absolute entry, proving (18)-(19).
QED.

This is a finite source-return certificate template. If an actual
returned matrix has \(G_{\rho,K,N}\succeq
(\varepsilon_{\rho,K,N}+\zeta)I\) for \(\zeta>0\), then
\(H_N(a_0)\succeq\zeta I\).
If only \(H_{K+N}(a)\succeq0\) is known, then \(G_{\rho,K,N}\succeq0\),
so \(H_N(a_0)\succeq-\varepsilon_{\rho,K,N}I\).

The Gram return uses actual moments through degree \(2(K+N-1)\).
It cannot be executed from fewer source layers by declaring the
missing ones zero. At rational \(\rho\), all entries of \(T\) are
rational; the source normalization \(f_\rho\) must still be retained
with its actual enclosure in any future numerical application.
No such new bank or threshold is evaluated in this note.

## MOR-T5. A hypothetical negative witness survives in linearly growing banks

Suppose a real polynomial \(p\) of degree \(d\) satisfies
\[
\mathcal L_{a_0}[p^2]=-\delta<0.
\]
Let \(B=B(p)>0\), and set
\[
\boxed{
\mathcal C(p)=2\log\!\left(\frac{2^{d+4}B^2}{\delta}\right),\qquad
K_\rho=\lceil \mathcal C(p)\rho\rceil.
}
\tag{20}
\]
Then \(\mathcal C(p)>0\), and for every \(\rho\ge1\),
\[
\boxed{
f_\rho\,
\mathcal L_a[P_{\rho,K_\rho}[p]^2]<-\delta/2.
}
\tag{21}
\]
Hence the actual bank \(H_{K_\rho+d+1}(a)\) is indefinite at
every anchor \(a=\rho a_0\).

**Proof.** Equation (5) implies
\(|\mathcal L_{a_0}[p^2]|\le2B^2\), so the argument of the logarithm
in (20) is at least \(2^{d+3}>1\).
By (15) and the choice of \(K_\rho\), the error between the left
side of (21) and \(-\delta\) is strictly less than
\[
2^{d+3}B^2 e^{-\mathcal C(p)/2}=\delta/2.
\]
This proves (21). The polynomial has degree at most \(K_\rho+d\),
and \(f_\rho>0\), so it is a negative direction of the stated Hankel
bank. That bank also has the positive constant direction \(c_0=1\),
making it indefinite. QED.

The normalized margin \(\delta/B^2\), and therefore \(\mathcal C(p)\), is
unchanged by multiplying \(p\) by a nonzero real scalar.
Its value and the degree \(d\) depend on the hypothetical witness.
This is not a universal numerical RH test size.

Theorem D2 of the retained detector supplies such a finite polynomial
if an off-line pair exists. MOR-T5 therefore proves that such a pair
could not disappear from all sufficiently resolved moving banks,
even though NIF-T3 makes each fixed-size bank positive eventually.
No actual negative polynomial, off-line pair, or value of
\(\delta\) is assumed to have been found.

Combining this with NIF-T3 gives a further precise consequence.
If RH fails, define
\[
n_-(a)=\min\{N\ge1:H_N(a)\not\succeq0\}.
\]
Then \(n_-(a)\) is finite at every \(a\ge a_0\), and
\[
n_-(a)\longrightarrow\infty\quad(a\to\infty),\qquad
n_-(\rho a_0)\le\mathcal C(p)\rho+d+2.
\]
The first statement follows because NIF-T3 makes every fixed finite
bank positive for all sufficiently large anchors. The upper bound is
(21) and the ceiling estimate. Thus any hypothetical obstruction
would move to higher observer depth, but this theorem supplies a
linear upper bound on the depth needed to retain one fixed witness.
No numerical lower growth rate is asserted.

## MOR-T6. A growing-depth criterion returns RH to the finite reference anchor

Let \(\rho_j\to\infty\) and let \(M_j\) be positive integers with
\[
\frac{M_j}{\rho_j}\longrightarrow\infty.
\tag{22}
\]
If the actual banks satisfy
\[
H_{M_j}(\rho_j a_0)\succeq0
\quad\text{for every }j,
\tag{23}
\]
then RH holds. The same conclusion applies if these conditions hold
along an unbounded subsequence. Conversely, RH makes every such
bank positive definite.

In particular, set for integers \(m\ge2\)
\[
M(m)=\lceil m\log(2m)\rceil.
\]
Then
\[
\boxed{
\mathrm{RH}
\quad\Longleftrightarrow\quad
H_{M(m)}(m a_0)\succeq0
\ \text{for arbitrarily large integers }m.
}
\tag{24}
\]
The logarithmic schedule is a convenient sufficient schedule;
any choice satisfying (22) gives the corresponding criterion.
No optimal growth rate is asserted.

**Proof.** Fix an arbitrary finite \(N\).
For sufficiently large \(j\), choose \(K_j=M_j-N\ge0\).
Under (23), the returned matrix (17) is positive semidefinite.
Equation (19) therefore gives
\[
H_N(a_0)\succeq
-2^{N+2}N
\exp\!\left(-\frac{M_j-N+1}{2\rho_j}\right)I.
\]
The right-side error tends to zero by (22).
Testing on each real coefficient vector and taking the limit proves
\(H_N(a_0)\succeq0\). Since \(N\) was arbitrary, every finite
unshifted reference bank is positive semidefinite.
The actual-source equivalence in detector Theorem D3 now implies RH.

Conversely, under RH all actual support values lie in \((0,1)\),
all weights in (4) are positive, and there are infinitely many
distinct support points. Every nonzero finite polynomial therefore
has strictly positive squared reading at every anchor. This proves
positive definiteness of all the banks in question.
Finally \(M(m)/m\to\infty\), so applying the first part along the
arbitrarily large integers in (24) proves that equivalence. QED.

NIF-1 proves positivity after increasing the anchor for each fixed
finite size. It does not prove (23) or (24), where the matrix size
increases with the anchor. MOR-T5 gives the corresponding explanation:
a hypothetical finite-anchor negative direction would require
increasing observer resolution to remain visible.

The new actual-source obligation is now explicit. It is sufficient
to prove positivity for the banks in (24), or for any unbounded
sequence satisfying (22)-(23), using the complete arithmetic source.
The fixed-point asymptotic at fixed derivative degree is not such
a proof.

## Evidence, attribution and continuation

| Statement | Evidence or status |
|---|---|
| Uniform outward-anchor source bounds | MOR-T1; retained product and published finite-height theorem |
| Exact moving-observer preservation of old pairings | MOR-T2 |
| Fixed-order moment return with linear-depth accuracy | MOR-T3 |
| Finite whole-bank return with explicit operator error | MOR-T4 |
| Persistence of a hypothetical negative witness | MOR-T5 |
| Growing-depth positivity criterion equivalent to RH | MOR-T6 and retained detector D3 |
| Actual positivity along the required growing banks | Open |
| New negative actual witness or off-line zero | None found or asserted |
| New numerical source bank, improved fixed-anchor positive region, or RH closure | Not claimed |

The generalized binomial expansion, analytic change of variables,
bilinear congruence, and finite-dimensional norm estimates are
classical tools. Their needed identities and bounds are derived here.
The programme contribution is their actual-source observer return,
the explicit source-weighted error, and the resulting quantified
connection between moving anchors and the fixed-anchor RH gate.
No literature-wide novelty or optimality claim is made.

This note follows the user's observer/invariant route while retaining
all source terms. It does not identify a thermodynamic capacity ratio
with a zeta observable without an additional adapter.

Only this note and its README entry are added. No local scientific
computation, independent agent, numerical replay, or new scientific
CI stage is used. These are written proofs, not new machine-checked
certificates. Existing engines, source data, certificates and workflows
are unchanged. The commit includes [skip ci]. The change remains a
draft stacked on PR #97 and is not merged.
