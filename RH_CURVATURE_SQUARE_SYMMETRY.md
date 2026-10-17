# CSS-1 — Curvature-square symmetry and the critical line

17 October 2026. Written analytic continuation of RBP-1.

The question is whether the programme's curvature and squared-generator
identities can supply a symmetry condition that forces the actual Riemann
zeros onto the critical line. This note recovers the precise native
antecedents and derives a condition on the unchanged arithmetic source.

The result is an exact bridge, not a proof of its remaining premise.
A curvature square loses a two-branch record. Reality of that square at
one source anchor can therefore conceal an off-line mode. Reality at two
distinct anchors, or reality of both the square and its first anchor
derivative at one anchor, removes this ambiguity **mode by mode** and is
equivalent to RH. A source moment version and a complete first-derivative
reconstruction are proved below. Reality or symmetry of aggregate
readings alone is insufficient; it already holds unconditionally.

No scientific execution, fresh source evaluation, independent agent or
workflow dispatch is used. The predecessor quantitative bounds and all
frozen packets are unchanged.

## 1. Exact native antecedents and their types

The sources in this table are read at predecessor commit
`5e12d79c34039d6002adc494890ed7f79b017163`.

| Source | Retained identity | Assumptions and scope |
|---|---|---|
| [NT-3](papers/native-thermodynamic-curvature/THEOREM.md) | \(A_i^H=X_i/2,\ X_i=H^{-1}\delta_iH;\ F^H_{ij}=-[X_i,X_j]/4\) | Metric compatibility and the selected strain representative in the declared response trivialization |
| [QG-R1, R3–R5](papers/native-relative-curvature/THEOREM.md) | \(F_{\rm rel}=D_BC+C\wedge C\) | Source and reference connections, derivative terms, and comparison directions are retained |
| [MGH-1, section 4.1](papers/native-matter-geometry-holonomy/THEOREM.md) | \(\Omega_{\rm MGH}^{\dagger}\Omega_{\rm MGH}=4M_o^2\otimes I\) | The stated finite Hermitian matter operator and selected configuration operations |
| [CWR-1](RH_CUT_WAVE_RETURN.md), [SPT-1](RH_SHUNYA_POTENTIAL.md) | \([P,Z]^2=-4y(1-y)I\) | Actual arithmetic source evaluation; no all-order positive pairing assumed |

In MGH, \(G=M_o\otimes X_g,\ J=I\otimes Z_g\), and
\(\Omega_{\rm MGH}=[-\iota G,-\iota J]\).
The connected diagnostic is
\(\mathcal C=4(I_2-I_1^2)\). Its relation to the chosen geometric hinge
is a theorem of that finite model. It is not an identification of an
arithmetic source mode with physical matter or spacetime curvature.

NT-3 does contain an exact selection of one half. For real constant
\(t\), consider its family \(A_i=tX_i\). Its metric-compatibility
defect is
\[
\delta_iH-A_i^\dagger H-HA_i=(1-2t)\delta_iH.
\tag{1}
\]
Thus, wherever some \(\delta_iH\ne0\), compatibility in this family
forces \(t=1/2\). Its curvature is
\(F^{tX}_{ij}=t(t-1)[X_i,X_j]\).
The metric-compatible strain representative is consequently curved
when this commutator is nonzero; symmetry does not mean flatness.
The general compatible connection still has NT-3's additional
metric-skew freedom.

This \(t\) is a transport coefficient. It has not been proved equal to
the real part \(\beta\) of an arithmetic zero. Equating the two labels
would insert the desired conclusion as an extra premise. The following
construction instead uses the actual zero/source map already retained
in the RH programme.

## 2. Keep the actual source and its normalization

For every finite anchor \(a\ge a_0=3249/4\), retain
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,
\]
\[
\delta=\beta-\tfrac12,\qquad
\mu=(\gamma-i\delta)^2,\qquad
u_\mu(a)=\frac a{a+\mu},\qquad
w_\mu(a)=\frac{\nu_\mu}{q(a)(a+\mu)}
        =\frac{\nu_\mu u_\mu}{a q(a)}.
\tag{2}
\]
The sum is over distinct source values arising from positive-height
zeros, with their multiplicities. CWR/OQD retain
\[
q(a)>0,\quad 0<|u_\mu|<1,\quad u_\mu\ne1,\quad
\sum_\mu |w_\mu|\le2,\quad \sum_\mu w_\mu=1.
\tag{3}
\]
The support is infinite, conjugation symmetric, and has only zero as
an accumulation point. Define
\[
\mathcal L_a[f]=\sum_\mu w_\mu f(u_\mu),\qquad
c_j(a)=\mathcal L_a[y^j],\qquad
\mathscr D=a\frac{d}{da}.
\tag{4}
\]
The derivative holds the arithmetic source \(\mu,\nu_\mu\) fixed.
It varies the observation anchor, not a zero or a physical clock.

Products in \(\mathcal L_a\) remain bilinear. In particular they are
not replaced by absolute squares. All polynomial readings in (4)
and their real-anchor derivatives are real by conjugation symmetry,
whether or not RH holds.

For later noncancellation, every individual weight has positive real
part:
\[
\Re u_\mu
=\frac{a(a+\gamma^2-\delta^2)}{|a+\mu|^2}>0,\qquad
\Re w_\mu>0.
\tag{5}
\]
Here \(|\delta|<1/2\) and \(a\ge a_0>1/4\) suffice. This is not a
claim that the complex weights are positive real numbers.

## CSS-T1. The source curvature square and its first derivative retain different records

Use the inherited native cut and reflection
\[
P(y)=\begin{pmatrix}y&1-y\\y&1-y\end{pmatrix},\qquad
Z=\operatorname{diag}(1,-1).
\]
Write
\[
\Omega(y):=[P(y),Z]
 =\begin{pmatrix}0&-2(1-y)\\2y&0\end{pmatrix},\qquad
k(y)=4y(1-y),\qquad t(y)=2y-1.
\]
Then
\[
\boxed{\Omega(y)^2=-k(y)I,\qquad t(y)^2=1-k(y).}
\tag{6}
\]
This is the existing native commutator, now used as a squared-curvature
diagnostic. It is a frozen two-operation commutator; it is not the full
curvature of an unspecified variable connection.

At an actual source mode put \(k_a(\mu)=k(u_\mu(a))\). Its first
anchor derivative satisfies
\[
\boxed{\mathscr D u_\mu=u_\mu(1-u_\mu),\qquad
\mathscr D k_a=k_a(1-2u_\mu).}
\tag{7}
\]
Since \(k_a\ne0\), the pair recovers the source coordinate:
\[
\boxed{
u_\mu=\frac12\left(1-\frac{\mathscr D k_a}{k_a}\right),\qquad
\mu=a\,\frac{k_a+\mathscr D k_a}{k_a-\mathscr D k_a}.
}
\tag{8}
\]
The denominator in the second expression is \(2u_\mu k_a\ne0\).

**Proof.** Multiply the two displayed matrices to obtain (6).
Differentiate \(a/(a+\mu)\) with \(\mu\) fixed:
\(\mathscr D u=a\mu/(a+\mu)^2=u(1-u)\).
The chain rule gives (7), and solving it gives (8). QED.

The square alone identifies \(u\) and \(1-u\). In source coordinates
this identification is
\[
\mu\longleftrightarrow a^2/\mu.
\tag{9}
\]
It is different from the functional-equation/conjugation pairing
\(\mu\leftrightarrow\overline\mu\).
For a real angular carrier \(y=\cos^2\theta\), (6) reads
\(k=\sin^2(2\theta)\); it loses the complementary-angle branch.

The native algebra also has formal metric-skewness
\(\Omega^Tg+g\Omega=0\), \(g=\operatorname{diag}(y,1-y)\).
It does not make \(g(u_\mu)\) a positive Hermitian metric at a complex
source point. Consequently \(-\Omega^2\) cannot be declared positive
by silently replacing it with an already positive
\(\Omega^\dagger\Omega\). MGH's adjoint identity has its own positive
finite-state premises.

## CSS-T2. Exact confinement from two anchors or one first derivative

Let
\[
D_a=|a+\mu|^2,\qquad r_\mu=\gamma^2+\delta^2=|\mu|.
\]
Then
\[
\boxed{
\Im k_a(\mu)=
\frac{8a\gamma\delta(r_\mu^2-a^2)}{D_a^2}.
}
\tag{10}
\]
For each actual mode,
\[
\boxed{k_a(\mu)\in\mathbb R
\quad\Longleftrightarrow\quad
\delta=0\ \text{or}\ a=\gamma^2+\delta^2.}
\tag{11}
\]

If \(\delta\ne0\), the exceptional anchor \(a=r_\mu\), when admitted,
has the explicit readings
\[
u_\mu=\tfrac12+\frac{i\delta}{2\gamma},\qquad
k_a=1+\frac{\delta^2}{\gamma^2}>1,\qquad
\mathscr D k_a=-\,i\frac{\delta}{\gamma}k_a\notin\mathbb R.
\tag{12}
\]
Thus its square is real and positive even though the zero is off-line.

For any two distinct finite anchors \(a,b\ge a_0\), the following are
equivalent:
\[
\boxed{
\begin{aligned}
&\mathrm{RH};\\
&k_a(\mu),k_b(\mu)\in\mathbb R
       &&\text{for every actual source mode }\mu;\\
&k_a(\mu),\mathscr D k_a(\mu)\in\mathbb R
       &&\text{for every actual source mode }\mu;\\
&k_a(\mu)\in(0,1]
       &&\text{for every actual source mode }\mu.
\end{aligned}}
\tag{13}
\]
The third and fourth lines each use just the one fixed anchor \(a\).

**Proof.** Write \(u=v+ih\). From (2),
\[
h=\frac{2a\gamma\delta}{D_a},\qquad
1-2v=\frac{|\mu|^2-a^2}{D_a}.
\]
Since \(\Im(4u(1-u))=4h(1-2v)\), this proves (10).
The anchors and \(\gamma\) are positive, giving (11).
Substitution at \(a=\gamma^2+\delta^2\) gives (12).

A nonzero \(\delta\) cannot have both \(a=r_\mu\) and \(b=r_\mu\).
This proves the two-anchor implication. If the two quantities in
the third line of (13) are real, (8) makes \(u_\mu\) real.
Then \(\Im u_\mu=2a\gamma\delta/D_a\) forces \(\delta=0\).
The last line excludes the only off-line real-square alternative
in (12). Conversely, under RH, \(\mu=\gamma^2>0\),
\(0<u_\mu<1\), and (7) and \(k=1-(2u-1)^2\) prove all three
conditions. QED.

The first-derivative condition is source-matched symmetry under complex
conjugation, not a demand that the derivative vanish. It forces
\(\beta=1/2\); it does not force \(u_\mu=1/2\).
The condition concerns every mode individually. Real values of the
summed curvature moments and their derivatives already follow from
pairing and do not establish (13).

There is also a precise generator control. In SPT's cut coordinates,
\[
\Omega=\begin{pmatrix}0&k/2\\-2&0\end{pmatrix},\qquad
D=\operatorname{diag}(1,k/4),\qquad W=tI+\Omega.
\]
At the off-line exceptional anchor (12), \(D\) is positive and
\(\Omega^\dagger D+D\Omega=0\). Nevertheless \(t=i\delta/\gamma\)
is nonreal, and the original reflection is not self-adjoint for
this positive pairing. Stability of the commutator flow alone
therefore does not establish the native wave's positive unitary
realization. Retaining the cut and complementary branch matters.

## CSS-T3. A source-defined curvature response cannot hide its nonreal poles

Define the actual curvature moments and response
\[
d_n(a)=\mathcal L_a[k(y)^n],\qquad
R_a^{\rm curv}(z)=\sum_{n\ge0}d_n(a)z^n
 =\sum_\mu\frac{w_\mu(a)}{1-zk_a(\mu)}
 \quad (|z|<1/8).
\tag{14}
\]
The indicated disk follows from \(|u_\mu|<1\), hence
\(|k_a(\mu)|<8\), and the total absolute weight budget.

For a distinct nonzero curvature value \(v\), let
\[
W_v(a)=\sum_{\mu:\,k_a(\mu)=v}w_\mu(a).
\]
There are at most two terms, since \(k(y)=4y(1-y)\) is quadratic
and the source map \(\mu\mapsto u_\mu\) is injective. Equation (5)
gives
\[
\Re W_v>0,\qquad W_v\ne0.
\tag{15}
\]
The distinct \(v\)'s have only zero as an accumulation point.
The response therefore extends meromorphically to the finite
\(z\)-plane, with the exact nonzero residues
\[
\boxed{\operatorname*{Res}_{z=1/v}R_a^{\rm curv}(z)=-W_v/v\ne0.}
\tag{16}
\]
Repeated or complementary modes cannot cancel these poles.

**Proof.** On any compact set avoiding the listed poles, the tail
has \(|zk_a(\mu)|\le1/2\) and is normally convergent by (3).
The finitely many remaining terms are rational. Grouping equal
curvature values and taking their residues proves (16). QED.

Thus (13) has an exact scalar-source version: reality of all poles
of \(R_a^{\rm curv}\) and \(R_b^{\rm curv}\), for two distinct
anchors, is equivalent to RH. At one anchor, the stronger condition
that all poles belong to \([1,\infty)\) is equivalent to RH.
These are global analytic conditions on the actual response,
not conclusions from finitely many coefficients.

There is also a real moment-bank formulation. Put
\[
C_N(a)=(d_{i+j}(a))_{0\le i,j<N}.
\tag{17}
\]
At each fixed anchor,
\[
\boxed{
C_N(a)\succeq0\ \text{for every }N
\quad\Longleftrightarrow\quad
k_a(\mu)\in\mathbb R\ \text{for every }\mu.
}
\tag{18}
\]
When these hold, the matrices are in fact positive definite.

To prove the reverse direction in (18), real curvature support is
stable under conjugation. Its grouped weights are consequently real,
and (15) makes them strictly positive. They sum to one. The support
is infinite, so no nonzero polynomial has zero squared reading.

For the forward direction, the finite-pair detector used in
[OQD](RH_OBSERVER_QUOTIENT_DESCENT.md) applies to the curvature
support without assuming positive weights. Here is its decisive
step. If \(v\) is nonreal, choose a real polynomial \(Q\) vanishing
on all other support points of modulus at least \(|v|\), except
\(v,\bar v\). There are only finitely many such points. On the two
observers \(x^nQ(x),x^{n+1}Q(x)\), the leading pair matrix has entries
\(2\Re(A_n v^{i+j})\), where \(A_n=W_vQ(v)^2v^{2n}\ne0\).
Its determinant is
\[
-4|A_n|^2(\Im v)^2<0.
\]
The remaining support has strictly smaller maximum modulus, so its
matrix error is \(o(|A_n|)\). The determinant remains negative for
all sufficiently large \(n\), contradicting all-bank positivity.
The absolute weight bound and bounded polynomial values justify the
tail estimate. This proves (18).

Combining (18) with T2 gives the exact two-anchor criterion
\[
\boxed{\mathrm{RH}\quad\Longleftrightarrow\quad
C_N(a)\succeq0\ \text{and}\ C_N(b)\succeq0
\quad\text{for every }N,\qquad a\ne b.}
\tag{19}
\]
A single all-order curvature bank can miss precisely the exceptional
circle \(\gamma^2+\delta^2=a\): at that anchor the conjugate off-line
pair has one real value \(k>1\) and a positive grouped weight.
This is a characterization of a possible missed mode, not an
assertion that such an arithmetic zero exists.

## CSS-T4. The first derivative restores the complete source pairing

A complete source-level alternative uses the curvature tower and
its anchor derivative, with no zero-by-zero input.

Differentiation of the actual normalized weights gives
\[
\frac{\mathscr D q(a)}{q(a)}=-c_1(a),\qquad
\mathscr D w_\mu=(c_1-u_\mu)w_\mu.
\tag{20}
\]
Consequently for every polynomial \(f\),
\[
\mathscr D\,\mathcal L_a[f]
=\mathcal L_a[(c_1-y)f+y(1-y)f'].
\tag{21}
\]
In particular, putting
\[
e_n=\mathcal L_a[y k^n],\qquad
h_n=\mathcal L_a[t k^n],
\]
one obtains
\[
\boxed{
\mathscr D d_n=(n+c_1)d_n-(2n+1)e_n,\qquad
h_n=\frac{(2c_1-1)d_n-2\mathscr D d_n}{2n+1}.
}
\tag{22}
\]
The derivative of a source average includes the weight derivative.
It is not merely the average of the mode derivative in (7).

**Proof.** Differentiate the convergent source sum for \(q\) to
get \(\mathscr Dq=-q c_1\). Differentiating (2) then gives (20).
Use (7) in the derivative of each weighted evaluation to get (21).
For \(f=k^n\),
\(y(1-y)f'=n k^n(1-2y)\).
This proves the first formula in (22); \(h_n=2e_n-d_n\)
proves the second. All differentiations are locally uniform at
admitted anchors: the source tails are bounded by convergent inverse
squared-height tails, and the extra rational factors and their
derivatives are uniformly bounded on each compact anchor interval.
QED.

Every real polynomial has a unique decomposition
\[
\boxed{f(y)=A(k(y))+t(y)B(k(y)),\qquad A,B\in\mathbb R[k].}
\tag{23}
\]
Indeed, write \(f((1+t)/2)\) as its even part in \(t\) plus
\(t\) times an even polynomial and use \(t^2=1-k\).
Thus the exact data
\[
c_1(a),\quad d_n(a),\quad\mathscr Dd_n(a)\quad(n\ge0)
\tag{24}
\]
determine every original polynomial reading \(\mathcal L_a[f]\)
and every bilinear reading \(\mathcal L_a[fg]\).
The normalization \(d_0=1\) is included. The value \(c_1\) is retained
source data; it is not inferred from the identity
\(\mathscr Dd_0=0\).

There is an explicit positive-symmetry test containing the missing
cross-record. Define
\[
C_N^+=(d_{i+j+1})_{i,j<N},\qquad
E_N=(h_{i+j})_{i,j<N},
\]
\[
\boxed{
\mathcal J_N(a)=
\begin{pmatrix}
C_N&E_N\\
E_N&C_N-C_N^+
\end{pmatrix}.
}
\tag{25}
\]
For coefficient vectors \(\alpha,\beta\) of
\(A(k)=\sum_{j<N}\alpha_j k^j\) and
\(B(k)=\sum_{j<N}\beta_j k^j\), its quadratic reading is
\[
\begin{pmatrix}\alpha\\\beta\end{pmatrix}^{T}
\mathcal J_N
\begin{pmatrix}\alpha\\\beta\end{pmatrix}
=\mathcal L_a[(A(k)+tB(k))^2].
\tag{26}
\]
The \(2N\) polynomials
\[
1,k,\ldots,k^{N-1},\quad t,tk,\ldots,tk^{N-1}
\]
form a basis of the original polynomials of degree below \(2N\).
Their leading degrees are the distinct even and odd degrees
\(0,\ldots,2N-1\). If \(S_N\) is their monomial coefficient matrix,
\[
\boxed{\mathcal J_N=S_N^T H_{2N}(a)S_N,\qquad S_N\ \text{invertible}.}
\tag{27}
\]
OQD's actual-source criterion now gives
\[
\boxed{\mathrm{RH}\quad\Longleftrightarrow\quad
\mathcal J_N(a)\succeq0\quad\text{for every }N}
\tag{28}
\]
at any one finite admitted anchor.

Equations (22), (25), and (27) are a faithful implementation of the
curvature proposal: curvature-square readings plus their correctly
normalized first derivatives restore the complementary branch and
the complete original pairing.

The blocks in (25) are real symmetric unconditionally. The unresolved
requirement is their positivity for every size, including the cross
block \(E_N\). Setting \(E_N=0\) would impose an additional
\(y\leftrightarrow1-y\) symmetry not supplied by the arithmetic
functional equation. For example \(h_0=2c_1-1\) is actual source
data, not a quantity that may be discarded by selecting a frame.

## CSS-T5. The existing finite positive range transfers exactly

This change of record retains RBP's proved finite progress. Let
\[
b=\sqrt a,\quad \ell=\log(b/(2\pi)),\qquad
n_*(a)=\left\lfloor\frac{e^{-4}\sqrt a}{4\pi}\right\rfloor.
\tag{29}
\]
For \(b\ge64\) and integers \(N\ge1\) with \(2N\le n_*(a)\),
RBP's reference-form comparison and (27) give
\[
\boxed{
\frac1{2\ell}\mathcal J_N^*
\preceq\mathcal J_N(a)
\preceq2\mathcal J_N^*,\qquad
\mathcal J_N(a)\succ0.
}
\tag{30}
\]
Here the star uses the retained arcsine source
\(d\nu_*=dy/(\pi\sqrt{y(1-y)})\) on \((0,1)\).

The reference blocks can also be read explicitly. The map
\(y\mapsto k=4y(1-y)\) preserves this arcsine distribution, while
\(t\) changes sign under \(y\mapsto1-y\). Therefore
\[
d_n^*=\frac{\binom{2n}{n}}{4^n},\qquad h_n^*=0,\qquad
\mathcal J_N^*=
\operatorname{diag}(H_N^*,\,H_N^*-A_N^*)\succ0.
\tag{31}
\]
For the distribution statement take \(y=(1+\cos\theta)/2\) with
uniform \(\theta\in(0,\pi)\); then \(k=\sin^2\theta\) has the same
arcsine law. Both reference integrals in (31) are strictly positive
on nonzero polynomials.

This is a congruence of the retained RBP theorem, not a new
numerically certified range or an improved square-root coefficient.
The first block of (30) also supplies finite curvature-bank
positivity at each anchor in that range. Two finite banks do not
satisfy the all-size requirement in (19).

## 3. The remaining obligation, precisely stated

There are now three exact routes to the same actual-source target:

1. Prove modewise reality of the curvature square at two distinct
   anchors, or of its square and first anchor derivative at one anchor.
2. Prove positivity of all actual curvature moment banks at two
   distinct anchors.
3. Prove positivity of every completed block (25) at one anchor,
   retaining the derivative-determined cross-record.

None of these universal premises has been proved here.
The native QG adjoint identity, NT-3's transport coefficient \(1/2\),
functional-equation reflection, conservation of the shunya potential,
and finite RBP positivity do not separately establish them.

The source-specific new content is the exact exceptional-circle
diagnostic, its removal by an anchor derivative or a second anchor,
noncancellation for the folded curvature response, and reconstruction
of the full pairing from the curvature tower and its first variation.
These are written analytic results in the declared programme.
No literature-wide novelty, machine-checked proof, physical curvature
identification, or resolution of RH is claimed.

The quantitative relational range of RBP, PRW's unrestricted support
window, and FCP's separate restricted support window remain unchanged.
All source weights, normalization, complementary records, and
all-order quantifiers are retained.

## References and continuity

- [NT-1–NT-7: native thermodynamic curvature](papers/native-thermodynamic-curvature/THEOREM.md), especially NT-3 and the specific Hessian tangent adapter NT-4.
- [QG-R1: native relative curvature](papers/native-relative-curvature/THEOREM.md), especially the full derivative and reference terms in R3–R5.
- [MGH-1: matter–geometry holonomy](papers/native-matter-geometry-holonomy/THEOREM.md), section 4.1 for the exact \(\Omega^\dagger\Omega\) rule and its model scope.
- [CWR-1: cut-wave return](RH_CUT_WAVE_RETURN.md), for the actual source, native matrices, wave gain, and source meromorphy.
- [SPT-1: shunya potential](RH_SHUNYA_POTENTIAL.md), for the same commutator in retained/complementary coordinates and its conserved bilinear potential.
- [OQD-1: observer quotient descent](RH_OBSERVER_QUOTIENT_DESCENT.md), for source faithfulness, the finite-pair detector, and all-polynomial positivity.
- [RBP-1: relational bank positivity](RH_RELATIONAL_BANK_POSITIVITY.md), for the unchanged explicit finite positive range.
