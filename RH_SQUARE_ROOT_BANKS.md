# SQR-1 — Square-root source banks and positive layer innovations

17 October 2026. Written continuation of the actual source and native
cut-wave programme, based on Publications commit
aa08780e457bae7f6e89cc7ff94d026c8fb822e4.

The preceding reflection-transfer note identifies the actual unshifted
Hankel banks as the folded transfer's reflection forms. This note proves
their positivity through an explicit size proportional to sqrt(a), where
a is the source anchor. QRI previously supplied a logarithmic positive
range. The improvement comes from estimating the full wave generating
function and its matrix coefficients together.

Every newly adjoined layer within this range has an explicitly positive
Schur innovation. The proof keeps the rational, Gamma and complete
prime-power terms of the actual xi source. It assumes neither RH nor
a positive source measure.

**Status: RH remains open in this programme.** The proven positive
range has order sqrt(a), whereas the retained sufficient return
schedule has order sqrt(a) log(a). These are written analytic
derivations, not additional machine-checked source certificates.

## Retained source and the folded reflection bank

For every real a>=a_0=3249/4, retain
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q=\mathcal X'/\mathcal X,\qquad
c_n(a)=\frac{(-a)^nq^{(n)}(a)}{n!\,q(a)},\qquad
H_N(a)=(c_{i+j}(a))_{0\le i,j<N}.
\tag{1}
\]
MOR and CWR give q(a)>0 and the actual, possibly complex, source
functional L_a with L_a[y^n]=c_n and L_a[1]=1. Conjugation makes
these moments real. All source products remain bilinear.

Write
\[
b=\sqrt a,\qquad \ell=\log\frac b{2\pi},\qquad
F(w)=wq(w^2)=\frac12\frac{\xi'}{\xi}(1/2+w).
\tag{2}
\]
Here ell is a scalar logarithm, distinct from L_a. Let
\[
r_k(a)=\mathcal L_a[T_k(2y-1)],\quad k\ge0,\qquad
u(z)=\frac{1-z}{1+z}.
\]
CWR-T2 gives the exact wave identity, initially near z=0,
\[
\boxed{
C_a(z):=1+2\sum_{k\ge1}r_k(a)z^k
=\frac{F(bu(z))}{F(b)}.
}
\tag{3}
\]
Every disk used below lies in an analytic region of this same
right-hand side, so its Taylor coefficients remain the actual r_k.

For reference use the arcsine law and its orthonormal polynomials
\[
d\nu_*(y)=\frac{dy}{\pi\sqrt{y(1-y)}},\qquad
\phi_0=1,\qquad \phi_j=\sqrt2\,T_j(2y-1)\quad(j\ge1).
\]
If Q_N contains their monomial coefficient columns, put
\[
H_N^*=\left(\int_0^1 y^{i+j}\,d\nu_*\right)_{i,j<N},\qquad
\widehat H_N=Q_N^TH_N(a)Q_N,\qquad Q_N^TH_N^*Q_N=I.
\tag{4}
\]
This reference law is an estimating device, not a replacement for
the actual source.

RFT-T3 to T4 retain the native P, Z, g and W=(2P-I)Z, with
\[
\mathsf T=\frac{2I+W+W^{-1}}4=yI,\qquad
\mathcal R=2P-I,\qquad d=(1,1)^T.
\]
For G_a(f,h)=L_a[f^T g h] they give
\[
\mathcal G_a(\mathsf T^i d,\mathcal R\mathsf T^j d)=c_{i+j}(a).
\tag{5}
\]
Thus the positivity below is positivity of the original folded
reflection bank. It concerns unshifted H_N; it does not assert the
separate finite shifted-bank conditions A_N>0 and B_N>0.

| Parameter | Meaning |
|---|---|
| a, b=sqrt(a) | Source anchor and its square root |
| N | Number of polynomial directions in the unshifted bank |
| n | Degree of a newly adjoined monic layer |
| z, varrho | Auxiliary complex wave variable and circle radius |
| rho=a/a_0 | Anchor ratio, used only when comparing with CWR return |

No parameter above is identified with local lambda displacement,
physical velocity, thermodynamic capacity or semigroup time.

## SQR-T1. A uniform bound for the full arithmetic source

For Re(w)>=2, with principal logarithms in the right half-plane,
\[
\boxed{
\left|F(w)-\frac14\log\frac w{2\pi}\right|
\le\frac2{|w|}+3\,2^{-\Re w}.
}
\tag{6}
\]

**Proof.** Logarithmic differentiation of the completed xi formula
in the Euler half-plane gives the same arithmetic identity retained
in RH_ARITHMETIC_GATE, now at complex w:
\[
F(w)=\frac{w}{w^2-1/4}
+\frac{\psi(1/4+w/2)-\log\pi}{4}
-\frac12\sum_{m\ge2}\Lambda(m)m^{-1/2-w}.
\tag{7}
\]
It is analytic for Re(w)>1/2. In that region the prime series is
absolutely convergent on compact subsets.

For completeness the Gamma error has an elementary uniform bound.
Let B_1(t)={t}-1/2, with arbitrary values at integers. Applying
first-order Euler summation to (z+t)^(-1), and using the digamma
series, gives for Re(z)>0
\[
\psi(z)=\log z-\frac1{2z}
+\int_0^\infty\frac{B_1(t)}{(z+t)^2}\,dt.
\tag{8}
\]
Indeed, the finite sum from 0 to M-1 is
\[
\log(z+M)-\log z
+\frac12\left(\frac1z-\frac1{z+M}\right)
-\int_0^M\frac{B_1(t)}{(z+t)^2}\,dt,
\]
and subtracting it from log(M) and taking the limit proves (8).
Absolute integrability justifies that limit.

Write z=R exp(i theta), with |theta|<pi/2. Direct integration gives
\[
\int_0^\infty\frac{dt}{|z+t|^2}
=\frac{|\theta|}{R|\sin\theta|}
\le\frac{\pi}{2R},
\]
where the value at theta=0 is 1/R. Hence
\[
|\psi(z)-\log z|
\le\frac{1/2+\pi/4}{|z|}.
\tag{9}
\]
For z=(w+1/2)/2, Re(w)>=2 implies |z|>=|w|/2. Also
\[
\left|\log z-\log(w/2)\right|
=\left|\log\left(1+\frac1{2w}\right)\right|
\le\frac2{3|w|}.
\]
The Gamma contribution to the error in (6) is consequently at most
(5/12+pi/8)/|w|. The rational term is at most
16/(15|w|), because |w|>=2.

For the complete prime-power term, Lambda(m)<=log(m) and the
decreasing function log(x)/x^2 on [2,infinity) give
\[
S:=\sum_{m\ge2}\frac{\log m}{m^2}
\le\frac{\log2}{4}+\int_2^\infty\frac{\log x}{x^2}\,dx
=\frac{\log2}{4}+\frac{\log2+1}{2}<2.
\]
Writing x=Re(w)>=2, we obtain
\[
\frac12\sum_{m\ge2}\Lambda(m)m^{-1/2-x}
\le\frac S2\,2^{3/2-x}<3\,2^{-x}.
\]
Finally 16/15+5/12+pi/8<2. Combining these estimates proves
(6), including the entire prime-power sum. QED.

The standard inputs to (8) are the
[NIST DLMF digamma series](https://dlmf.nist.gov/5.7#E6) and
[Euler summation](https://dlmf.nist.gov/2.10#i).
The needed bound is derived here rather than invoked as an
asymptotic formula with unspecified constants.

## SQR-T2. The wave remainder on a circle adapted to bank size

Suppose
\[
b\ge64,\qquad 2\le N\le b/8,\qquad
\varrho=1-\frac1N,\qquad M=\frac{1+\varrho}{1-\varrho}=2N-1.
\tag{10}
\]
Define the exact analytic remainder
\[
E_a(z)=C_a(z)-1-\frac1\ell\log\frac{1-z}{1+z}
=2\sum_{k\ge1}e_k(a)z^k,\qquad e_0=0.
\tag{11}
\]
Then
\[
\boxed{
\sup_{|z|\le\varrho}|E_a(z)|
\le\frac{128N}{b\ell}.
}
\tag{12}
\]
The leading wave coefficients separated in (11) are
\[
r_k^{(0)}=
\begin{cases}
-1/(k\ell),&k\text{ odd},\\
0,&k\ge2\text{ even},
\end{cases}
\qquad e_k=r_k-r_k^{(0)}\quad(k\ge1).
\tag{13}
\]

**Proof.** On the closed disk in (10), the Cayley map obeys
\[
\Re u\ge\frac1M,\qquad
\frac1M\le|u|\le M,\qquad
|\log u|\le\log M+\frac\pi2.
\tag{14}
\]
Since M<=b/4, Re(bu)>=4 throughout the disk. Formula (7)
is analytic on a neighborhood of its image. Thus (3) and (11)
are analytic on a neighborhood of this closed disk.

Put
\[
d(w)=F(w)-\tfrac14\log(w/(2\pi)),\qquad d_0=d(b).
\]
At b>=64, (6) gives |d_0|<=3/b, using 3*2^(-b)<=1/b.
Also ell>2 and b ell>=24, so
\[
F(b)=\ell/4+d_0\ge\ell/8>0.
\tag{15}
\]
The logarithms in (14) satisfy
\[
\frac{|\log u|}{\ell}
\le\frac{\log(b/4)+\pi/2}{\ell}<3.
\]
Here log(b/4)=ell+log(pi/2), log(pi/2)<1 and pi/2<2.

For x>=4, 2^(-x)<=1/x. Applying (6) at bu gives
\[
|d(bu)|\le\frac{2M}{b}+\frac{3M}{b}
\le\frac{10N}{b}.
\]
The exact quotient identity is
\[
E_a(z)=
\frac{d(bu)-d_0(1+\log(u)/\ell)}{F(b)}.
\tag{16}
\]
Its numerator has modulus at most 10N/b+12/b<=16N/b.
Equation (15) proves (12). Finally
log((1-z)/(1+z))=-2 sum_(k odd) z^k/k proves (13).
No term of the actual source has been discarded. QED.

## SQR-T3. A matrix remainder estimate without exponential coefficient loss

For the same b and N,
\[
\boxed{
\widehat H_N=I+\frac{\mathsf D_N}{2\ell}+\mathsf R_N,\qquad
\|\mathsf R_N\|_2\le\frac{896N}{b\ell}.
}
\tag{17}
\]
The real symmetric correction is precisely the QRI correction:
\[
D_{ij}=0\quad(i,j\text{ of the same parity}),\qquad
D_{0j}=D_{j0}=-\frac{2\sqrt2}{j}\quad(j\text{ odd}),
\]
\[
D_{ij}=-2\left(\frac1{|i-j|}+\frac1{i+j}\right)
\quad(i,j\ge1\text{ of opposite parity}).
\tag{18}
\]

**Proof.** The product identity 2T_iT_j=T_(i+j)+T_|i-j|
and (13) give (17)-(18), with exact remainder entries
\[
R_{00}=0,\qquad R_{0j}=R_{j0}=\sqrt2\,e_j,\qquad
R_{ij}=e_{|i-j|}+e_{i+j}\quad(i,j\ge1).
\tag{19}
\]
It remains to estimate this whole matrix from the analytic
function E_a instead of bounding each monomial coefficient.

Here is the required elementary matrix estimate. Let
E(z)=2 sum_(k>=1) e_k z^k have real coefficients, e_0=0,
be analytic near |z|<=varrho, and satisfy |E|<=M_E there.
Let D_varrho=diag(1,varrho,...,varrho^(N-1)).
The strictly lower triangular matrix
\[
(A_\varrho)_{ij}=
\begin{cases}
e_{i-j}\varrho^{i-j},&i>j,\\
0,&i\le j
\end{cases}
\]
is the compression to the first N coefficients of multiplication
by E(varrho z)/2 on the unit circle. Parseval's identity gives
||A_varrho||_2<=M_E/2. Thus A=D_varrho^(-1) A_varrho D_varrho
satisfies
\[
\|A+A^T\|_2\le\varrho^{-(N-1)}M_E.
\tag{20}
\]
Its entries are e_|i-j|.

Likewise the matrix B_varrho=(e_(i+j) varrho^(i+j)) has norm
at most M_E/2. To see this directly, put
P_x(z)=sum_(i<N) x_i z^i. For complex vectors x and y,
\[
x^*B_\varrho y
=\frac1{2\pi}\int_0^{2\pi}
\frac{E(\varrho e^{it})}{2}
\,\overline{P_x(e^{it})}\,P_y(e^{-it})\,dt.
\]
Cauchy--Schwarz and Parseval bound this by
(M_E/2)||x||_2||y||_2. Therefore B=D_varrho^(-1) B_varrho
D_varrho^(-1) obeys
\[
\|B\|_2\le\frac12\varrho^{-2(N-1)}M_E.
\tag{21}
\]
Its entries are e_(i+j).

For J=diag(1/sqrt(2),1,...,1), the matrix in (19) is
J(A+A^T+B)J. Since ||J||_2=1,
\[
\|\mathsf R_N\|_2
\le\left(\varrho^{-(N-1)}
+\tfrac12\varrho^{-2(N-1)}\right)M_E<7M_E.
\tag{22}
\]
For the last inequality use varrho=1-1/N,
varrho^(-(N-1))<=e and e+e^2/2<7. Equations (12) and
(22) prove (17). QED.

The auxiliary circle norms in this proof control coefficients.
They do not assume that L_a[p^2] is positive. In particular no
absolute-square substitution is made in the actual complex
source functional.

## SQR-T4. An explicit positivity gate for each admissible bank

For every N>=2 the leading correction satisfies
\[
\boxed{\|\mathsf D_N\|_2\le2\log N+8+2\pi.}
\tag{23}
\]
Consequently, whenever (10) holds, define
\[
d_N(a)=\frac{\log N+4+\pi+896N/b}{\ell}.
\]
Then
\[
\boxed{
\|\widehat H_N-I\|_2\le d_N(a),\qquad
(1-d_N)H_N^*\preceq H_N(a)\preceq(1+d_N)H_N^*.
}
\tag{24}
\]
In particular d_N(a)<1 is an unconditional sufficient positivity
condition for the actual bank.

**Proof.** On indices 1,...,N-1, separate (18) into its difference
and sum parts. The difference part has absolute row sum at most
\[
4\sum_{\substack{1\le k<N\\k\ {\rm odd}}}\frac1k
\le4+2\log N.
\]
This bounds its symmetric operator norm.

For the sum part, the matrix with entries 1/(i+j), i,j>=1,
has norm at most pi, by the weighted Schur estimate with
weights i^(-1/2):
\[
\sum_{j\ge1}\frac1{(i+j)\sqrt j}
\le\int_0^\infty\frac{dx}{(i+x)\sqrt x}
=\frac\pi{\sqrt i}.
\]
The integrand is decreasing; this proves the first inequality.
Dropping the parity restriction provides an entrywise majorant,
so the sum part of (18) has norm at most 2pi.

Finally the row coupling index 0 to the other indices has
squared Euclidean norm
\[
8\sum_{\substack{1\le j<N\\j\ {\rm odd}}}\frac1{j^2}
\le8\left(1+\int_0^\infty\frac{dx}{(2x+1)^2}\right)
=12<16.
\]
The symmetric matrix supported on this row and column has
operator norm less than 4. Adding these three estimates proves
(23). Equation (17) then gives the norm bound in (24).
Conjugating by Q_N^(-1) gives the displayed comparison with
H_N^*. QED.

## SQR-T5. Unconditional square-root positive banks

Set the explicit conservative constant and bank schedule
\[
\boxed{
\kappa=\frac{e^{-9}}{2\pi},\qquad
N_+(a)=\max\{1,\lfloor\kappa\sqrt a\rfloor\}.
}
\tag{25}
\]
For every a>=a_0 and every integer 1<=N<=N_+(a),
\[
\boxed{H_N(a)\succ0.}
\tag{26}
\]
For every N>=2 in this range there is the stronger comparison
\[
\boxed{
\frac1{2\ell}H_N^*\preceq H_N(a)\preceq2H_N^*.
}
\tag{27}
\]

**Proof.** H_1=(c_0)=(1), so only N>=2 needs work.
Then N<=kappa b forces b>=2/kappa>64 and N<=b/8.
Furthermore
\[
\log N\le\log(\kappa b)=\ell-9,
\]
so (24) gives
\[
1-d_N(a)\ge\frac{5-\pi-896\kappa}{\ell}.
\]
Using e>2 and 2pi>6,
kappa<1/3072 and 896kappa<7/24. Since pi<4,
\[
5-\pi-896\kappa>\frac{17}{24}>\frac12.
\]
This proves the lower bound in (27), and d_N<1 gives its
upper bound. H_N^* is positive definite, proving (26). QED.

The constant is not claimed optimal. At the reference anchor this
conservative schedule only asserts H_1>0; the stronger retained
finite certificates there remain available. The advance is an
explicit growing analytic range for large anchors.

For comparison, fix 0<eta<1/2 and put
N_eta(a)=floor(a^(1/2-eta)). The same bound, when a is sufficiently
large, gives
\[
\limsup_{a\to\infty}d_{N_\eta(a)}(a)\le1-2\eta,\qquad
\liminf_{a\to\infty}
\lambda_{\min}(\widehat H_{N_\eta(a)})\ge2\eta.
\tag{28}
\]
Indeed log(N_eta)/ell tends to 1-2eta and N_eta/b tends to zero.
Thus the normalized positivity margin can stay bounded away
from zero at every fixed sub-square-root power.

QRI's earlier remainder bound amplified a scalar coefficient
error by 100^N. Equation (17) instead controls the whole wave
matrix by 896N/(b ell). Combined with the coefficient 2 of
log N in (23), this is what makes (25) possible.

## SQR-T6. Every covered next-layer innovation is strictly positive

For n>=1 let
\[
v_n=(c_n,c_{n+1},\ldots,c_{2n-1})^T,\qquad
s_n(a)=c_{2n}-v_n^T H_n(a)^{-1}v_n.
\tag{29}
\]
Whenever n+1<=N_+(a), this Schur complement exists and satisfies
\[
\boxed{
\frac{2^{-4n}}{\ell}\le s_n(a)\le2^{2-4n}.
}
\tag{30}
\]
In particular det(H_(n+1))/det(H_n)=s_n>0 throughout the
covered chain.

**Proof.** Positive definiteness from T5 permits completing the
square in the coefficients of every monic polynomial of degree n:
\[
s_n(a)=
\min_{p(y)=y^n+\text{terms of degree}<n}\mathcal L_a[p^2].
\tag{31}
\]
For the arcsine reference the minimizing monic polynomial is
2^(1-2n) T_n(2y-1). Its squared norm is
\[
s_n^*=2^{1-4n},
\tag{32}
\]
because the squared norm of T_n is 1/2 for n>=1.
Applying (27) to every such polynomial, and taking minima, gives
s_n^*/(2ell)<=s_n(a)<=2s_n^*. This is exactly (30).
The block determinant identity gives the determinant ratio. QED.

More generally, under the gate (10) and d_(n+1)<1,
the same argument gives
(1-d_(n+1))s_n^*<=s_n(a)<=(1+d_(n+1))s_n^*.
These are actual layer signs within a specified range. They do
not establish acceptance of every subsequent Schur step at one
fixed anchor, nor finite shifted-bank positivity beyond the
separate retained certificates.

## SQR-T7. A sharper location of any hypothetical RH obstruction

If RH fails, the retained detector makes
\[
n_-(a)=\min\{N:H_N(a)\not\succeq0\}
\]
finite for every a>=a_0. The new unconditional positive range gives
\[
\boxed{
n_-(a)>N_+(a),\qquad
\liminf_{a\to\infty}\frac{n_-(a)}{\sqrt a}\ge\kappa
\quad\text{if RH fails}.
}
\tag{33}
\]
CWR-T7 retains the conditional upper bound from any one
hypothetical negative reference polynomial p:
\[
n_-(\rho a_0)\le
2\sqrt\rho\log\rho+C_p\sqrt\rho+1,\qquad \rho\ge1.
\tag{34}
\]
Here C_p=(4/3) log(288B_4(p)^2/delta), with
L_(a_0)[p^2]=-delta<0. No such actual witness is asserted.

**Proof.** Equation (26) excludes every size up to N_+(a).
Divide this exclusion by sqrt(a) and take the lower limit to
obtain (33). Equation (34) is exactly the retained CWR-T7 return
of the hypothetical witness, with its constants unchanged. QED.

Equivalently the conditional depth interval is now
\[
\boxed{
\liminf_{\rho\to\infty}
\frac{n_-(\rho a_0)}{\sqrt\rho}\ge\kappa\sqrt{a_0},\qquad
\limsup_{\rho\to\infty}
\frac{n_-(\rho a_0)}{\sqrt\rho\log\rho}\le2.
}
\tag{35}
\]
This raises the previous logarithmic lower bound to square-root
growth, while preserving the known square-root-log upper bound.

## The remaining source assertion

CWR-T6 proves the exact equivalence
\[
\boxed{
\mathrm{RH}\quad\Longleftrightarrow\quad
H_{\lceil4\sqrt m\log(2m)\rceil}(m a_0)\succeq0
\quad\text{for arbitrarily large integers }m.
}
\tag{36}
\]
At these anchors T5 supplies sizes asymptotic to
kappa sqrt(a_0) sqrt(m), whereas (36) asks for order
sqrt(m) log(m). Their ratio grows like a constant times log(m).

The same distinction is visible directly in CWR's sufficient
return condition: M/sqrt(rho)-2 log(rho) must tend to +infinity.
Our schedule makes the first term tend to kappa sqrt(a_0);
it does not satisfy that condition.

Therefore the written derivation advances the actual positivity
range but does not close the return argument. Proving positivity
at the schedule (36), or proving all unshifted signs at one
fixed anchor by another valid source argument, remains the RH
obligation. A sufficient return upper bound is not a proof that
a logarithmic factor is necessary for every possible method.

## Evidence and continuation scope

| Statement | Evidence or status |
|---|---|
| Explicit full-source half-plane bound | SQR-T1; complete arithmetic source and derived Gamma estimate |
| Bank-dependent analytic wave remainder | SQR-T2; exact CWR generating identity |
| Linear-in-N remainder and logarithmic leading norm | SQR-T3 to T4; coefficient operators and weighted Schur estimate |
| Actual unshifted positivity through explicit kappa sqrt(a) | SQR-T5; written uniform inequalities |
| Strictly positive next-layer innovations in that range | SQR-T6; finite Schur minimization |
| Square-root lower obstruction depth if RH fails | SQR-T7; no negative actual witness asserted |
| Positivity at the CWR square-root-log return schedule | Open |
| New numerical source certificate or physical mass identification | Not supplied |

Source lineage:
[arithmetic source](RH_ARITHMETIC_GATE.md),
[finite detector](RH_FINITE_OBSERVER_DETECTOR.md),
[NIF](RH_INFINITY_CUT_FIXED_POINT.md),
[MOR](RH_MOVING_OBSERVER_RETURN.md),
[QRI](RH_REFLECTION_RATIO_GATE.md),
[CWR](RH_CUT_WAVE_RETURN.md),
[OQD](RH_OBSERVER_QUOTIENT_DESCENT.md), and
[RFT](RH_REFLECTION_TRANSFER.md).

The classical ingredients are Euler summation, the digamma series,
Chebyshev products, Parseval's identity, weighted matrix estimates
and Schur complements. Their needed applications are given above.
No literature-wide novelty or optimality claim is made.

Only this note and its README entry are added. No local scientific
computation, independent agent, numerical source evaluation or
scientific replay is performed. Existing data, engines, certificates
and workflows remain unchanged. These are written proofs, not
new machine-checked certificates. The documentation commit includes
[skip ci]; the change is a draft stacked on PR #102 without a merge.
