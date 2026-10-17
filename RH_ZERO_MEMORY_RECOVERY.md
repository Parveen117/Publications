# ZMR-1 — Zero-chart memory, reference motion and complete source recovery

17 October 2026. Written continuation of GRR-1 at
ce872fb61ccb2d1b1254cc325604820fd218a443.

This note tests the proposal that a cut's infinitely close readings,
together with the memory of its moving reference, should recover the
source and settle its sign. The recovery part has an exact realization
in the existing native cut algebra. The sign implication does not
follow from that recovery alone.

There are two concrete results.

First, the complete local germ of one native cut-return, together
with its first anchor derivative as a germ in the same variable,
recovers every actual source moment and the actual meromorphic
response. A fixed-reference constant return by itself can identify
distinct sources: it sees y(1-y), which loses the orientation between
y and 1-y. The reference derivative supplies the missing channel.

Second, the transport preserving the cut-memory form needs an
extra inverse-square-root Jacobian relative to the transport
preserving the ordinary source pairing. Its exact weight has
power 3/2 where the ordinary observer weight has power 1/2.
Both its finite transport and its source-derived connection are
given below.

The reconstructed memory matrices retain the RH-equivalent sign
condition. An explicit finite comparison source has a bounded,
positive-definite constant cut-return and a negative memory
direction. It is not the xi source or an RH counterexample.
It shows why complete reconstruction and carried invariance
cannot replace the actual-source sign proof.

## Retained source and the meaning of the zero chart

Use the actual source of CEM, MOR, LAT and FWT:
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q(z)=\frac{\mathcal X'(z)}{\mathcal X(z)}
=\sum_\mu\frac{\nu_\mu}{z+\mu}.
\]
At a finite anchor \(a\ge a_0=3249/4\), put
\[
y_a(\mu)=\frac a{a+\mu},\qquad
w_\mu(a)=\frac{\nu_\mu}{q(a)(a+\mu)},\qquad
\mathcal L_a[p]=\sum_\mu w_\mu(a)p(y_a(\mu)).
\tag{1}
\]
The actual source values are
\(\mu=(\gamma-i\delta)^2\), with
\(\delta=\beta-1/2\), positive zero height gamma and multiplicity
\(\nu_\mu\). Retain the proved source properties
\[
q(a)>0,\quad \Re\mu>0,\quad
0<|y_a(\mu)|<1,\quad y_a(\mu)\ne1,\quad
\sum_\mu|w_\mu(a)|\le2,\quad\mathcal L_a[1]=1 .
\tag{2}
\]
The support and weights are conjugation symmetric, so readings
of real polynomials are real. This is a bilinear source
functional; its point evaluations are not replaced by absolute
squares.

Write
\[
v(y)=y(1-y),\quad
B_a(f,h)=\mathcal L_a[fh],\quad
\mathcal M_a(f,h)=\mathcal L_a[vfh],\quad
c_n(a)=\mathcal L_a[y^n]
=\frac{(-a)^nq^{(n)}(a)}{n!\,q(a)}.
\tag{3}
\]
CEM's specified native exponential generator is
\[
A(y)=
\begin{pmatrix}0&-2(1-y)\\2y&0\end{pmatrix},
\qquad A^2=-4vI,\qquad PAP=0,
\]
with
\[
P(y)=\begin{pmatrix}y&1-y\\y&1-y\end{pmatrix},
\quad g(y)=\operatorname{diag}(y,1-y),\quad Jf=(f,f).
\]
Its even scalar return is the entire function
\[
C(t,v)=\sum_{k=0}^\infty\frac{(-4v)^kt^{2k}}{(2k)!},
\qquad P e^{tA}P=C(t,v)P.
\tag{4}
\]
The displayed series requires no square-root choice.

The zero chart here is the regular response coordinate \(t=0\),
where the compressed generator is zero. It is not a zero
vector carrying nonzero linear evolution. Three parameters
retain distinct roles:

| Parameter | Role |
|---|---|
| t | CEM's specified exponential-return parameter |
| \(\tau=\log(a/a_{\rm ref})\) | Logarithmic anchor motion |
| \(\lambda=(z-a)/a\) | Local argument of the recovered source response |

The formulas below connect their readings; they do not identify
these parameters with each other or with physical time.

CEM already retains the rescaled memory of repeated cuts:
\[
-\frac n4
\left.\partial_t^2\mathcal L_a[C(t/n,v)^nfh]\right|_{t=0}
=\mathcal M_a(f,h).
\tag{5}
\]
Thus that memory was not absent from the prior construction.
The question here is what one scalar reading determines, and
how this memory changes under anchor motion.

## ZMR-T1. One constant return can lose orientation even with its source scale

Define the normalized and unnormalized returns
\[
C_{a,j}(t)=\mathcal L_a[y^j C(t,v)],\qquad
R_j(a,t)=q(a)C_{a,j}(t),\qquad j=0,1 .
\tag{6}
\]
Both are entire even functions of t, real on the real axis.
Their complete derivatives at zero give
\[
\boxed{
C_{a,j}^{(2k)}(0)=(-4)^k\mathcal L_a[y^jv^k],
\qquad C_{a,j}^{(2k+1)}(0)=0.
}
\tag{7}
\]
The constant return \(C_{a,0}\) sees only polynomials in v.
It need not determine the moments of y.

**Proof of analyticity.** On the source disk \(|y|\le1\),
\(|v|\le2\). The series (4) and all its fixed t-derivatives
converge uniformly on bounded complex t-sets and this disk.
The absolute source budget (2) permits termwise pairing and
differentiation, proving (7). This concerns the bounded native
generator (4); it does not assert differentiability to all
orders of FWT's different, unsmoothed continuous wave.

**Exact comparison example.** Fix \(A\ge a_0\), and consider
the two finite rational sources
\[
q_1(z)=\frac3{z+3A},\qquad
q_2(z)=\frac1{z+A/3}.
\tag{8}
\]
They are comparison sources, not xi. At the common anchor A,
\[
q_1(A)=q_2(A)=\frac3{4A}.
\]
Their sole normalized source coordinates are respectively
\(y_1=1/4\) and \(y_2=3/4\). Both have \(v=3/16\), so
\[
C_{A,0}^{[1]}(t)=C_{A,0}^{[2]}(t)
=\cos(\sqrt3\,t/2)
\quad\hbox{for every complex }t.
\]
Their unnormalized returns also agree. Nevertheless
\(c_1^{[1]}=1/4\) and \(c_1^{[2]}=3/4\), and the sources in
(8) are different. Thus no number of t-derivatives of this
one fixed-reference return separates them. QED.

This is a failure of injectivity for the specified scalar
observation on a larger source class. It does not claim that
two distinct xi functions satisfy the same arithmetic definition.

## ZMR-T2. The first reference derivative recovers the missing channel

Let \(\mathscr D=a\partial_a\), differentiating the anchor with
the source values and multiplicities held fixed. Then
\[
\boxed{
(1+t\partial_t)R_1(a,t)
=\frac t2\partial_tR_0(a,t)-\mathscr D R_0(a,t).
}
\tag{9}
\]
In particular, for \(t\ne0\),
\[
\boxed{
R_1(a,t)=\frac1t\int_0^t
\left[\frac s2\partial_sR_0(a,s)-\mathscr D R_0(a,s)\right]ds,
}
\tag{10}
\]
with the removable value
\(R_1(a,0)=-\mathscr D q(a)\).
The local germs of \(R_0(a,\cdot)\) and
\(\mathscr D R_0(a,\cdot)\) at one anchor determine the whole
germ of \(R_1(a,\cdot)\).

**Proof.** For each fixed source value,
\[
\mathscr D y=v,\qquad
\mathscr D v=v(1-2y),\qquad
\mathscr D\left(\frac{\nu_\mu}{a+\mu}\right)
=-y\,\frac{\nu_\mu}{a+\mu}.
\tag{11}
\]
Also the entire series (4) gives
\(2v\partial_v C=t\partial_t C\). Differentiating the
unnormalized source sum for \(R_0\) therefore yields
\[
\begin{aligned}
\mathscr D R_0
&=-R_1+
 \sum_\mu\frac{\nu_\mu}{a+\mu}
 v(1-2y)\partial_v C(t,v)\\
&=-R_1+\frac t2\partial_tR_0-t\partial_tR_1.
\end{aligned}
\]
The source sums and these derivatives converge locally
uniformly in a regular positive-anchor neighborhood and on
bounded complex t-sets: the source denominators are locally
comparable, y and its indicated derivatives are bounded,
and the entire series is uniformly bounded there.
This justifies the differentiation and proves (9).

The left side of (9) is \(\partial_t(tR_1)\).
Integrate from zero. Analyticity excludes a nonzero integration
constant divided by t, proving (10) and uniqueness. QED.

At the level of derivatives, (9) is the explicit reconstruction
\[
\boxed{
R_1^{(2k)}(a,0)=
\frac{k R_0^{(2k)}(a,0)
-\mathscr D R_0^{(2k)}(a,0)}{2k+1}.
}
\tag{12}
\]
Here superscripts denote t-derivatives only. For k=0 this
recovers the actual first layer
\(c_1=-a q'(a)/q(a)\).
The reference response is derived from the fixed actual source;
no freely chosen connection is used to make a desired reading
constant.

The data in (10) comprise two complete analytic germs in t.
They are not two numbers or finitely many near-zero samples.
No noise-stability claim is made for recovering all derivatives.

## ZMR-T3. The two recovered germs determine every source layer

Every real polynomial in y has a unique representation
\[
p(y)=F(v(y))+yG(v(y)),\qquad F,G\in\mathbb R[v].
\tag{13}
\]
Consequently (7) and (12) recover \(\mathcal L_a[p]\) for every
polynomial p, and therefore every \(c_n(a)\).

**Proof of the algebra.** Since \(y^2=y-v\), reduction gives
existence. More explicitly, put \(U_0(v)=0\), \(U_1(v)=1\) and
\[
U_{n+1}(v)=U_n(v)-vU_{n-1}(v).
\]
Then for \(n\ge1\),
\[
\boxed{y^n=yU_n(v)-vU_{n-1}(v).}
\tag{14}
\]
The n=1 and n=2 cases are immediate, and the quadratic relation
proves the recurrence. For uniqueness, replace y by \(1-y\)
in a vanishing expression (13) and subtract. This gives
\((2y-1)G(v(y))=0\) as a polynomial identity. Thus \(G=0\),
then \(F=0\). QED.

For example, writing \(C_j=C_{a,j}\),
\[
\begin{aligned}
c_0&=C_0(0)=1,\qquad c_1=C_1(0),\\
c_2&=c_1+\frac14C_0''(0),\\
c_3&=c_1+\frac14\bigl(C_0''(0)+C_1''(0)\bigr),\\
c_4&=c_1+\frac14\bigl(C_0''(0)+2C_1''(0)\bigr)
              +\frac1{16}C_0^{(4)}(0).
\end{aligned}
\tag{15}
\]
These are exact source identities, with the second channel
obtainable from the mixed reference readings (12).

The normalized actual response is then recovered by the
convergent germ
\[
\boxed{
\frac{q(a(1+\lambda))}{q(a)}
=\mathcal L_a\!\left[\frac1{1+\lambda y}\right]
=\sum_{n=0}^\infty(-\lambda)^n c_n(a),
\qquad |\lambda|<1.
}
\tag{16}
\]
Indeed the geometric expansion converges uniformly on the
source disk for \(|\lambda|<1\), and (2) justifies pairing it.
Also \(q(a)=R_0(a,0)\), so the scale is retained.
The identity theorem then uniquely determines the actual
meromorphic response on its continuation domain. Its poles
and residues are not removed by this reconstruction.

Integrating q determines \(\mathcal X(z)/\mathcal X(a)\) on
each simply connected zero-free domain reached from a.
An independent value \(\mathcal X(a)\) would be needed to fix
the absolute multiplicative normalization of \(\mathcal X\).
The recovered q already determines its zero locations and
multiplicities. This is exact determination from complete
analytic data, not a finite algorithm for all zeros or a sign proof.

## ZMR-T4. Cut memory has its own source-derived connection

For a differentiable family of admissible analytic profiles
\(f_\tau(y)\), the normalized source satisfies
\[
\boxed{
\frac d{d\tau}\mathcal L_a[f_\tau]
=\mathcal L_a\!\left[
(\partial_\tau+v\partial_y)f_\tau+(c_1-y)f_\tau
\right],\qquad a=a_{\rm ref}e^\tau .
}
\tag{17}
\]
For profiles and derivatives uniformly bounded on the source
support locally in tau, this follows by differentiating (1):
the weight derivative is \((c_1-y)w_\mu\), because
\(\mathscr D q=-q c_1\), and (11) gives the coordinate derivative.

The connections compatible with the ordinary pairing and the
memory pairing, respectively, are
\[
\begin{aligned}
\nabla_\tau^B&=\partial_\tau+v\partial_y+\frac{c_1-y}{2},\\
\boxed{\nabla_\tau^{\mathcal M}
=\partial_\tau+v\partial_y+\frac{1+c_1-3y}{2}.}
\end{aligned}
\tag{18}
\]
They obey
\[
\begin{aligned}
\frac d{d\tau}B_a(f,h)
 &=B_a(\nabla_\tau^B f,h)+B_a(f,\nabla_\tau^B h),\\
\frac d{d\tau}\mathcal M_a(f,h)
 &=\mathcal M_a(\nabla_\tau^{\mathcal M} f,h)
   +\mathcal M_a(f,\nabla_\tau^{\mathcal M} h).
\end{aligned}
\tag{19}
\]

**Proof.** Apply (17) to fh for the first identity.
For the second, apply it to vfh and use
\(v\partial_y v=v(1-2y)\). The extra coefficient on fh is
\(c_1-y+1-2y=1+c_1-3y\). Split it equally between the two
slots, proving (18)--(19). QED.

Thus carrying a profile parallel for B is not automatically
parallel transport for \(\mathcal M\). The difference of the
two connections is multiplication by \((1-2y)/2\).
The distinction comes from the measured factor v in the
specified cut memory, not from adding a new positivity axiom.

Metric compatibility here means the bilinear identities (19).
It holds before positivity is known. TET's separate positive
thermodynamic Hessian cannot be substituted for this actual
arithmetic pairing without proving a matching source identity.

## ZMR-T5. The finite memory transport has an extra Jacobian factor

Fix two anchors \(b=\rho a\), with \(\rho\ge1\), and put
\[
d_\rho(y)=1+(\rho-1)y,\qquad
\phi_\rho(y)=\frac{\rho y}{d_\rho(y)},\qquad
f_\rho=\frac{q(b)}{q(a)},\qquad r=1-\rho^{-1}.
\]
The exact coordinate and weight laws are
\[
y_b=\phi_\rho(y_a),\quad
w_\mu(b)=\frac{w_\mu(a)}{f_\rho d_\rho(y_a)},\quad
v(\phi_\rho(y))=\frac{\rho v(y)}{d_\rho(y)^2},
\quad \phi_\rho'(y)=\frac{\rho}{d_\rho(y)^2}.
\tag{20}
\]
Define the carried memory observer, in the new coordinate Y, by
\[
\boxed{
(\mathsf T^\mathcal M_{b\leftarrow a}p)(Y)
=\sqrt{\frac{f_\rho}{\rho}}\,
 (1-rY)^{-3/2}
 p\!\left(\frac{Y}{\rho(1-rY)}\right).
}
\tag{21}
\]
The branch is the analytic one equal to one at \(Y=0\).
For \(\rho>1\), its singularity is outside the closed unit disk,
at \(Y=1/r>1\); for \(\rho=1\) it is the identity.

For all real polynomial profiles p,h,
\[
\boxed{
\mathcal M_b(\mathsf T^\mathcal M_{b\leftarrow a}p,
             \mathsf T^\mathcal M_{b\leftarrow a}h)
=\mathcal M_a(p,h).
}
\tag{22}
\]
These transports compose exactly on their carried source class.
Their reverse is the inverse on that class. Their infinitesimal
parallel equation is \(\nabla_\tau^\mathcal M p_\tau=0\).

**Proof.** At corresponding source points,
\(1-r y_b=1/d_\rho(y_a)\). Thus the new profile equals
\[
\sqrt{\frac{f_\rho}{\rho}}\,d_\rho(y_a)^{3/2}p(y_a).
\]
The new weight times the new v times the product of two such
profiles is, by (20),
\[
\frac{w_\mu(a)}{f_\rho d_\rho}\,
\frac{\rho v(y_a)}{d_\rho^2}\,
\frac{f_\rho}{\rho}d_\rho^3\,p(y_a)h(y_a)
=w_\mu(a)v(y_a)p(y_a)h(y_a).
\]
Summing proves (22). The transported functions are analytic
near the closed source disk at every fixed outward ratio,
so the absolute source budget justifies the operation.

The coordinate maps compose by multiplication of anchor ratios.
The weight factors telescope. Equivalently, at a source point,
\[
d_\rho(y_a)=\frac{b+\mu}{a+\mu}.
\]
Use the continuous logarithms of \(a+\mu\) along the positive
anchor axis; \(\Re\mu>0\) prevents zeros or branch crossings.
Their differences telescope, proving composition and the
inverse statement on the carried class.

Finally the logarithmic derivative of the displayed profile
factor along a source coordinate is
\[
\frac12(-c_1-1)+\frac32y
=-\frac{1+c_1-3y}{2}.
\]
Together with \(\dot y=v\), this is precisely the parallel
equation from (18). QED.

MOR's ordinary source observer has instead the factor
\[
\mathsf T^B_{b\leftarrow a}p
=\sqrt{f_\rho}(1-rY)^{-1/2}
 p\!\left(\frac{Y}{\rho(1-rY)}\right).
\]
At corresponding points, the additional memory factor is
\[
\boxed{
\frac{\text{memory observer factor}}{\text{ordinary observer factor}}
=\frac{d_\rho(y)}{\sqrt\rho}
=\bigl(\phi_\rho'(y)\bigr)^{-1/2}.
}
\tag{23}
\]
This is the exact Jacobian compensation for the chosen memory.
No ratio of profile values at their zeros is intended in (23).

Every finite bank of carried memory observers has exactly
its old Gram matrix by (22). Its positive and negative directions
are therefore preserved. The transported observers are generally
analytic functions, not polynomials of the old finite degree.
No new finite truncation-error bound or finite-bank range is
asserted by this exact transport identity.

This also clarifies the two different motions: CEM-T4 carries
memory under its exponential at a fixed anchor. Equation (22)
carries the same memory when the anchor itself changes.

## ZMR-T6. The recovered jets give the full memory sign matrix

Define the recovered sequences
\[
m_k=\mathcal L_a[v^{k+1}],\qquad
n_k=\mathcal L_a[yv^{k+1}]\qquad(k\ge0).
\]
They come directly from the response derivatives:
\[
m_k=(-4)^{-k-1}C_{a,0}^{(2k+2)}(0),\qquad
n_k=(-4)^{-k-1}C_{a,1}^{(2k+2)}(0).
\tag{24}
\]
For \(N\ge1\), form
\[
M_N=(m_{i+j})_{0\le i,j<N},\quad
N_N=(n_{i+j})_{0\le i,j<N},\quad
M_N^+=(m_{i+j+1})_{0\le i,j<N}.
\]
The complete memory Gram matrix in the ordered bank
\[
1,v,\ldots,v^{N-1},\quad y,yv,\ldots,yv^{N-1}
\]
is
\[
\boxed{
\mathsf K_N(a)=
\begin{pmatrix}
M_N&N_N\\
N_N&N_N-M_N^+
\end{pmatrix}.
}
\tag{25}
\]
It is congruent to CEM's memory matrix \(D_{2N}(a)\), where
\[
D_{2N}(a)=
(c_{i+j+1}-c_{i+j+2})_{0\le i,j<2N}.
\]
In fact their determinants agree. Hence the precise inherited
condition is
\[
\boxed{
\mathrm{RH}\quad\Longleftrightarrow\quad
\mathsf K_N(a)\succeq0\ \hbox{for every }N,
}
\tag{26}
\]
at any one retained finite anchor.

**Proof.** The first block in (25) pairs \(v^i\) with \(v^j\)
under the v-weighted source; the cross block pairs \(v^i\)
with \(yv^j\). For the last block use \(y^2=y-v\), giving
\(n_{i+j}-m_{i+j+1}\).
The interleaved basis \(v^i,yv^i\), \(0\le i<N\), has successive
degrees \(0,1,\ldots,2N-1\), with leading coefficients
\((-1)^i,(-1)^i\). Its coefficient change from monomials is
invertible with determinant one; grouping the bank only changes
the sign of that determinant. Thus (25) and \(D_{2N}\) are
congruent and have equal determinants.

All even-size conditions imply the odd-size ones by restriction.
CEM-T5 proves that positivity of every \(D_N\) is equivalent
to RH by its retained actual-source finite detector.
This proves (26). QED.

The reconstruction in ZMR-T3 now specifies every entry of the
actual sign test from the zero-chart and reference-memory data.
It does not prove their sign. In particular, positivity of
\(M_N\) alone ignores the orientation and the last block.

If \(M_N\succ0\) for a specified N, its exact remaining
finite Schur condition is
\[
\boxed{
\mathsf K_N\succeq0
\quad\Longleftrightarrow\quad
N_N-M_N^+-N_NM_N^{-1}N_N\succeq0.
}
\tag{27}
\]
This is ordinary block completion with all cross-memory retained.
Neither its hypothesis nor its conclusion is asserted at every N.

The desired factorization into a positive metric would establish
(26) only if it reproduces the actual bilinear source. Writing
\(\mathcal M_a(f,h)=\mathcal G_a(AJf,AJh)/4\), as CEM does,
is exact, but \(\mathcal G_a\) is not yet a positive metric.
Replacing it by an ordinary Euclidean norm would change the form.

## ZMR-T7. Stable scalar returns and exact transport do not force the memory sign

There is a finite comparison source respecting conjugation,
positive integer residues and the same algebraic transport laws,
whose constant native return is bounded and positive definite
on the whole real line while its memory form is indefinite.

Fix \(A\ge a_0\), set
\[
\delta=\frac14,\qquad
\gamma=\sqrt{A-\delta^2},\qquad
\mu=(\gamma-i\delta)^2,\qquad
q_{\rm toy}(z)=\frac1{z+\mu}+\frac1{z+\overline\mu}.
\tag{28}
\]
This is not the Riemann source and does not retain its exact
Euler/Gamma identity. In particular it supplies no actual
off-line zero of xi.

At the anchor A, let
\[
\eta=\frac{\delta}{2\gamma}>0,\qquad
y=\frac12+i\eta,\qquad
v=y(1-y)=\frac14+\eta^2 .
\]
Then
\[
q_{\rm toy}(A)=\frac1A,\qquad
y_A(\mu)=y,\quad y_A(\overline\mu)=\overline y,\qquad
\mathcal L_{\rm toy}[p]=y\,p(y)+\overline y\,p(\overline y).
\tag{29}
\]
Consequently
\[
\boxed{C_{{\rm toy},0}(t)=\cos(2\sqrt v\,t).}
\tag{30}
\]
For real t this is bounded by one and is positive definite:
it is the equal-weight sum of the two real-frequency
characters with frequencies \(2\sqrt v\) and \(-2\sqrt v\).
Its scalar memory is \(\mathcal M_{\rm toy}(1,1)=v>0\).

Nevertheless, with \(c_1=1/2-2\eta^2\), the two-observer
memory bank is
\[
\boxed{
\mathsf K_1=
v\begin{pmatrix}1&c_1\\c_1&c_1-v\end{pmatrix},
\qquad
\det\mathsf K_1=-4\eta^2v^3<0.
}
\tag{31}
\]
Equivalently the explicit polynomial \(p(y)=y-c_1\) has
\[
\boxed{\mathcal M_{\rm toy}(p,p)=-4\eta^2v^2<0.}
\tag{32}
\]

**Proof.** Since \(A=\gamma^2+\delta^2\),
\(A+\mu=2\gamma(\gamma-i\delta)\). Its reciprocal and its
conjugate sum to \(1/A\), and multiplication by A gives
\(y=1/2+i\delta/(2\gamma)\). This proves (29).
The two v-values coincide, and the two weights sum to one,
giving (30).

Also
\[
\mathcal L_{\rm toy}[y]=y^2+\overline y^2=1-2v=c_1,\qquad
\mathcal L_{\rm toy}[y^2]=y^3+\overline y^3=1-3v=c_1-v.
\]
Insert these in the memory bank to obtain (31):
\[
(c_1-v)-c_1^2=v(1-4v)=-4\eta^2v .
\]
Completing the scalar square gives (32). QED.

All derivative, recovery and transport identities above hold
for this finite source by the same exact algebra. Its negative
memory direction remains negative when transported by (21).
Thus the identities do not, on their own, enforce positivity.

The wave in (30) is CEM's constant native cut-return.
It is not FWT's faithful continuous wave \(\mathscr W_a(u)\);
their parameters and generators are different. There is no
contradiction with FWT's actual-wave growth criterion.
For the toy, the folded frequency exceeds one, and the full
source still contains the nonreal pair in (28).

## Outcome of the proposed route

| Question | Established result |
|---|---|
| Was cut memory already retained? | Yes: CEM's exact rescaled second response (5) |
| Does the full germ of one fixed-reference scalar always recover the source? | No: (8) has equal returns and equal source scales but different source moments |
| What additional observation restores the orientation? | The first anchor derivative as a complete germ, through (9)--(12) |
| Can all other source layers then be recovered? | Yes, exactly through the polynomial decomposition and convergent response germ (13)--(16) |
| How is memory carried when the reference moves? | The source-derived connection (18) and the Jacobian-compensated transport (21)--(23) |
| Does faithful recovery establish positivity? | No: the explicit comparison source (28)--(32) retains a negative memory direction |
| What remains for the actual xi source? | All-size positivity in (26), or an independently derived equivalent actual-source factorization |

The useful correction is not to erase reference motion. Its
first mixed response supplies an orientation channel missing
from a constant fixed-reference cut-return. The memory pairing
also requires its own Jacobian compensation when carried
between anchors. Both statements are derived from the retained
actual source and do not assume RH.

The exact recovery is an information statement. Its input is
infinite analytic data, and no uniform stable recovery from
finitely many close measurements is proved. Its output is the
actual signed source, including any obstruction that source
might contain.

The next sign task can now be expressed in the recovered
blocks (25), or in (27) where its hypothesis is available.
No proof that these blocks are positive at every depth has
been obtained. In particular this note does not certify the
GRR first-prime residual condition.

GRR's full-family quarter window and improved mean-zero
first-prime margin are unchanged. PRW's stronger smaller-window
bound, FCP's separate restricted interval, and native finite-bank
ranges are unchanged. The all-width source positivity and
FWT growth obligations, hence RH, remain open.

## References and verification scope

- [CEM-1](RH_CUT_EQUILIBRIUM_MEMORY.md): specified native exponential, zero compressed generator, rescaled memory, its RH-equivalent sign condition, and polynomial observability.
- [MOR-1](RH_MOVING_OBSERVER_RETURN.md): actual source budget and square-root weighted ordinary observer transport.
- [LAT-1](papers/native-layer-anchor-transport/THEOREM.md): actual source hierarchy and coordinate/weight motion.
- [TET-1](papers/native-tvsp-exponential-transport/THEOREM.md): source-derived connections, moving-frame terms, and the distinct positive thermodynamic chart.
- [RFT-1](RH_REFLECTION_TRANSFER.md): regular zero chart, faithful source transfer and the distinction between bilinear adjoint products and positive factorizations.
- [FWT-1](RH_FAITHFUL_WAVE_TRANSPORT.md): the separate continuous wave, actual-source recovery and remaining global growth criterion.
- [GRR-1](RH_GAMMA_REMAINDER_RETURN.md): unchanged local source bounds and the unproved first-prime sufficient certificate.

Analytic jets, polynomial reduction, bilinear metric transport,
Jacobian weights and Schur completion are standard tools.
The new assembly here is the explicit mixed-reference recovery
identity, the memory-specific anchor connection and factor,
and the recovered orientation-sensitive sign matrix.
No literature-wide novelty or optimality is asserted.

Only this note and its README entry change. Verification consists
of written derivation review, document checks, immutable-head
content and draft PR metadata. No scientific computation,
numerical xi evaluation, independent agent, source replay,
finite certificate evaluation, or workflow dispatch is performed.
The displayed comparison sources are checked by written exact
algebra. Frozen sources, engines, certificates and workflows
remain unchanged. The commit includes [skip ci].
The PR is stacked on draft #128 and remains unmerged.
