# NIF-1 — The actual source at infinite anchor and its limiting cut

17 October 2026. Written continuation of LAT-1, ILC-1, LOT-1 and the
finite observer detector, based on Publications commit
3defc08c2f59c9d30e23664461dd16439b44064d.

The actual completed-zeta response has a scale-normalized limiting
shape. Its layers tend to the arcsine moments, and its fixed finite
cuts tend to explicit Chebyshev cuts. The constant seed readings
1/2 and 1/pi follow from this limit without assuming RH.

This note also states precisely which infinity is taken, how the
limiting law remains unchanged under the retained moving-anchor rule,
and when a cut reading can be independent of the prepared state.
The result is an actual-source asymptotic theorem, not a derivation
of all physical laws or of a dimensional speed of light.

## Coordinates and retained definitions

Retain
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\qquad
a_0=3249/4.
\]
The original source displacement is \(\lambda_0=(t-a_0)/a_0\).
Move the anchor along
\[
a=a_0(1+\Lambda)=a_0e^\tau,\qquad \Lambda>-1.
\tag{1}
\]
Here \(\Lambda\) is the original displacement evaluated at the new
anchor; \(\tau=\log(a/a_0)\) is LAT's logarithmic anchor parameter.
Near this new anchor use a separate local displacement \(z\):
\[
F_a(z)=\frac{q(a(1+z))}{q(a)},\qquad
c_n(a)=\frac{(-a)^nq^{(n)}(a)}{n!\,q(a)},\qquad
F_a(z)=\sum_{n\ge0}(-1)^nc_n(a)z^n.
\tag{2}
\]
The relation to the original displacement is
\(1+\lambda_0=(1+\Lambda)(1+z)\).
Derivative depth \(n\), matrix size \(N\), anchor displacement
\(\Lambda\), and a transport parameter are distinct quantities.
No identification with a thermodynamic response component or a
physical clock is assumed.

For sufficiently large positive \(a\), Theorem NIF-T1 below proves
\(q(a)>0\). At such anchors define the real symmetric forms
\[
H_N(a)=(c_{i+j}(a)),\quad A_N(a)=(c_{i+j+1}(a)),\quad
B_N(a)=(c_{i+j}(a)-c_{i+j+1}(a)),\quad 0\le i,j<N.
\tag{3}
\]
The previously certified finite anchors and their data are unchanged.

## NIF-T1. The actual arithmetic source has a normalized infinity limit

Put
\[
L_a=\log\frac{\sqrt a}{2\pi}.
\]
For each fixed \(0<r<1\), uniformly on \(|z|\le r\) as real
\(a\to+\infty\),
\[
q(a(1+z))=
\frac{L_a+\frac12\log(1+z)}{4\sqrt a\,\sqrt{1+z}}
+O_r(a^{-1}).
\tag{4}
\]
All roots and logarithms in this disk are the branches equal to one
and zero, respectively, at \(z=0\). In particular \(q(a)>0\)
eventually, and
\[
\boxed{
F_a(z)=\frac1{\sqrt{1+z}}
+\frac{\log(1+z)}{2L_a\sqrt{1+z}}
+O_r\!\left(\frac1{\sqrt a\,L_a}\right).
}
\tag{5}
\]
Thus the normalized source germs converge locally uniformly on the
unit disk to
\[
\boxed{F_*(z)=(1+z)^{-1/2}.}
\tag{6}
\]

For every fixed \(n\ge0\), define
\[
m_n=\frac{(1/2)_n}{n!}=\frac{\binom{2n}{n}}{4^n},\qquad
S_n=\sum_{j=0}^{n-1}\frac1{j+1/2},\quad S_0=0.
\]
Then
\[
\boxed{
c_n(a)=m_n\left(1-\frac{S_n}{2L_a}\right)
+O_n\!\left(\frac1{\sqrt a\,L_a}\right).
}
\tag{7}
\]
The error constants here may depend on the fixed order. In particular,
\[
-\frac{a q'(a)}{q(a)}
=c_1(a)=\frac12-\frac1{2L_a}
+O\!\left(\frac1{\sqrt a\,L_a}\right)\longrightarrow\frac12.
\tag{8}
\]
This is a dimensionless logarithmic-response limit.

**Proof.** The retained actual arithmetic identity is
\[
q(t)=\frac1{t-1/4}
+\frac{\psi(1/4+\sqrt t/2)-\log\pi}{4\sqrt t}
-\frac1{2\sqrt t}\sum_{m\ge2}\Lambda_{\rm vM}(m)m^{-1/2-\sqrt t}.
\tag{9}
\]
Here \(\Lambda_{\rm vM}\) is the von Mangoldt function, not the anchor
coordinate. The identity holds in the Euler half-plane
\(\operatorname{Re}(1/2+\sqrt t)>1\), by logarithmic differentiation
of the completed xi formula. See the retained
[arithmetic derivation](RH_ARITHMETIC_GATE.md), equations (2)-(6).

For \(t=a(1+z)\), \(|z|\le r<1\), the principal root satisfies
\(\operatorname{Re}\sqrt t\ge\sqrt{a(1-r)}\). Consequently the entire
disk enters that zero-free half-plane as \(a\) increases. If
\(\sigma=\operatorname{Re}(1/2+\sqrt t)\ge2\), then
\[
\left|\sum_{m\ge2}\Lambda_{\rm vM}(m)m^{-1/2-\sqrt t}\right|
\le 2^{-(\sigma-2)}\sum_{m\ge2}\frac{\log m}{m^2}.
\]
The last series is finite. Thus the complete prime-power contribution
is exponentially small uniformly on the disk; it is bounded, not
silently omitted from the source.

The sectorial digamma asymptotic gives
\[
\psi(1/4+\sqrt t/2)=\log(\sqrt t/2)+O_r(a^{-1/2}).
\]
The rational term in (9) is \(O_r(a^{-1})\). Substitution proves (4).
At \(z=0\), it gives \(q(a)=L_a/(4\sqrt a)+O(a^{-1})\), which is
positive for large \(a\). Division proves (5); the denominator's
relative error is \(O(1/(\sqrt a L_a))\).

The errors in (5) are holomorphic on any smaller fixed disk.
Cauchy's coefficient formula therefore permits each fixed Taylor
coefficient to be taken with the same type of bound. Since
\[
(1+z)^{-\alpha}
=\sum_{n\ge0}(-1)^n\frac{(\alpha)_n}{n!}z^n,\qquad
\log(1+z)(1+z)^{-1/2}
=-\left.\partial_\alpha(1+z)^{-\alpha}\right|_{\alpha=1/2},
\]
and \(\partial_\alpha(\alpha)_n=(\alpha)_n
\sum_{j=0}^{n-1}(\alpha+j)^{-1}\), equation (7) follows.
Its \(n=1\) instance proves (8). QED.

For comparison, keeping the original anchor fixed and simply sending
its argument to infinity gives a different displayed reading:
\[
\boxed{
\lim_{\Lambda\to+\infty}
\frac{\sqrt\Lambda}{\log\Lambda}
\frac{q(a_0(1+\Lambda))}{q(a_0)}
=\frac1{8\sqrt{a_0}\,q(a_0)}.
}
\tag{10}
\]
The unscaled response tends to zero. Equation (10) follows directly
from (4) at \(z=0\) with \(a=a_0(1+\Lambda)\).
The constant depends on the declared normalization; multiplying by
\(8\sqrt{a_0}q(a_0)\) makes the limit one by that normalization.
Neither this operation nor (8) introduces physical units.

The classical inputs are the completed xi definition
[DLMF 25.4.4](https://dlmf.nist.gov/25.4#E4), the full Mangoldt series
[DLMF 27.4.12](https://dlmf.nist.gov/27.4#E12), and the sectorial
digamma expansion and remainder bounds
[DLMF 5.11.2 and Section 5.11(ii)](https://dlmf.nist.gov/5.11).
The source-normalized conclusions above are derived from those inputs.

## NIF-T2. Moving coordinates can preserve the whole limiting law

For any probability measure \(\nu\) on \([0,1]\), write
\[
F(z)=\int_0^1\frac{d\nu(y)}{1+zy}.
\]
The retained LAT anchor transformation is
\[
(\mathcal T_\rho F)(z)=\frac{F(\rho-1+\rho z)}{F(\rho-1)},\quad \rho>0,
\tag{11}
\]
or, in spectral coordinates,
\[
\phi_\rho(y)=\frac{\rho y}{1+(\rho-1)y},\qquad
d\nu_\rho(\phi_\rho(y))
=\frac{d\nu(y)}{[1+(\rho-1)y]F(\rho-1)}.
\tag{12}
\]
These formulas include a weight change as well as coordinate motion.

The probability responses fixed by every \(\mathcal T_\rho\) are
exactly
\[
\boxed{F_\alpha(z)=(1+z)^{-\alpha},\qquad 0\le\alpha\le1.}
\tag{13}
\]
For \(0<\alpha<1\), their unique compact moment measures are
\[
d\nu_\alpha(y)=
\frac{y^{\alpha-1}(1-y)^{-\alpha}}{B(\alpha,1-\alpha)}\,dy.
\tag{14}
\]
The endpoints are \(\nu_0=\delta_0\) and \(\nu_1=\delta_1\).
Their layers are \(m_n(\alpha)=(\alpha)_n/n!\).

The actual source limit in NIF-T1 selects \(\alpha=1/2\), hence
\[
\boxed{
d\nu_*(y)=\frac{dy}{\pi\sqrt{y(1-y)}},\qquad 0<y<1.
}
\tag{15}
\]
In particular,
\[
\boxed{
\lim_{a\to\infty}D_{r,k}(a)
=\int_0^1y^r(1-y)^k\,d\nu_*(y)
=\frac{(2r)!(2k)!}{4^{r+k}r!k!(r+k)!}>0
}
\tag{16}
\]
for every fixed pair of nonnegative integers \(r,k\).

**Proof.** Let \(g(x)=F(x-1)>0\) for \(x>0\).
Being fixed by (11) for every \(\rho\) is equivalent to
\(g(\rho x)=g(\rho)g(x)\). Continuity and \(g(1)=1\) imply
\(g(x)=x^b\): the function \(\log g(e^t)\) is continuous and additive.
Moreover \(b=g'(1)=F'(0)=-\int y\,d\nu(y)\in[-1,0]\).
Writing \(\alpha=-b\) proves the necessity of (13).
Direct substitution proves sufficiency.

For (14), put \(u=y/(1-y)\). The integral of its unnormalized
density divided by \(1+zy\) becomes
\(\int_0^\infty u^{\alpha-1}/[1+(1+z)u]\,du\).
For real \(z>-1\), rescaling \(u\) extracts
\((1+z)^{-\alpha}\), with normalization \(B(\alpha,1-\alpha)\).
Analytic continuation gives the response on its slit domain.
The beta integral gives its moments
\((\alpha)_n/n!\). Compact polynomial density gives uniqueness.

At \(\alpha=1/2\), substitute \(y=(1+\cos\theta)/2\);
the measure is uniform \(d\theta/\pi\) on \(0<\theta<\pi\).
This proves (15) and its moments \(m_n\).
Taking the finite binomial combination in (7) gives (16);
the beta integral and half-integer Gamma values give its factorial
form. QED.

This is an exact version of a distribution remaining invariant while
its coordinates move. Under \(\rho=e^\tau\), individual coordinates
satisfy \(\dot y=y(1-y)\), yet the law (15) is fixed by the combined
coordinate-and-weight transformation (12). Without the weight
renormalization, the pushforward of this measure tends to
\(\delta_1\) as \(\rho\to\infty\); it is not stationary.

The same conclusion is visible in LAT's source hierarchy:
\[
\dot c_n=(n+c_1)c_n-(n+1)c_{n+1}.
\]
The sequences \(m_n(\alpha)\) solve its stationary equations because
\((n+1)m_{n+1}(\alpha)=(n+\alpha)m_n(\alpha)\).
Conversely, a stationary compact probability moment sequence obeys
this recurrence with \(\alpha=c_1\), and therefore belongs to (13).

The beta family, beta integral and classical arcsine law are standard.
See [DLMF 5.12.1](https://dlmf.nist.gov/5.12#E1).
The fixed-point classification here concerns the specified LAT
renormalization, not every possible transport or every physical law.

## NIF-T3. Finite actual cuts approach an explicit limiting cut

For every fixed \(N\), the three matrices in (3) converge entrywise
to \(H_N^*,A_N^*,B_N^*\) obtained by replacing \(c_n(a)\) by \(m_n\).
All three limiting matrices are strictly positive. Their determinants
are
\[
\boxed{
\det H_N^*=2^{-(N-1)(2N-1)},\qquad
\det A_N^*=\det B_N^*=2^{-N(2N-1)}.
}
\tag{17}
\]
Consequently, for each fixed \(N\), both actual forms
\(A_N(a),B_N(a)\) are positive for every sufficiently large \(a\).
This existence statement does not compute a new numerical threshold.

Use Gram--Schmidt with positive leading coefficients to express
\(Y_N(a)=H_N(a)^{-1}A_N(a)\) in its source-orthonormal polynomial
basis, with first vector \(e_0=1\). Denote that real symmetric matrix
by \(J_N(a)\). For each fixed \(N\),
\[
J_N(a)\longrightarrow J_N^*,
\]
where \(J_1^*=(1/2)\), and for \(N\ge2\), \(J_N^*\) has diagonal
\(1/2\), first off-diagonal entry \(1/(2\sqrt2)\), and all subsequent
off-diagonal entries \(1/4\). In particular the five-direction limit is
\[
\boxed{
J_5^*=\frac14
\begin{pmatrix}
2&\sqrt2&0&0&0\\
\sqrt2&2&1&0&0\\
0&1&2&1&0\\
0&0&1&2&1\\
0&0&0&1&2
\end{pmatrix}.
}
\tag{18}
\]
It is an exact limiting matrix, not a new evaluation at the archived
finite anchor.

The limiting finite eigenvalues and seed weights are
\[
y_{j,N}=\frac{1+\cos\theta_{j,N}}2,\quad
\theta_{j,N}=\frac{(2j-1)\pi}{2N},\quad
v_{j,N}=\frac1N,\qquad 1\le j\le N.
\tag{19}
\]
Apply the retained ILC cut construction in this basis:
\[
\Pi_N^*=
\begin{pmatrix}
J_N^*&\sqrt{J_N^*(I-J_N^*)}\\
\sqrt{J_N^*(I-J_N^*)}&I-J_N^*
\end{pmatrix}.
\tag{20}
\]
The actual finite cuts, expressed in the same prescribed
orthonormalization convention, converge in norm to (20) for fixed \(N\).

Let \(E_1=(e_0,0)\) and \(E_2=(0,e_0)\).
Their two-by-two matrix of limiting cut readings is
\[
\boxed{
\big(\langle E_i,\Pi_N^*E_j\rangle\big)_{i,j=1}^2
=
\begin{pmatrix}1/2&s_N\\s_N&1/2\end{pmatrix},\qquad
s_N=\frac1{2N\sin(\pi/(2N))}.
}
\tag{21}
\]
Thus \(s_5=(1+\sqrt5)/10\), and as \(N\to\infty\),
\[
\boxed{
M_*=\begin{pmatrix}1/2&1/\pi\\1/\pi&1/2\end{pmatrix}.
}
\tag{22}
\]
This last two-by-two matrix is a compression of the cut onto the two
seed channels; it is not itself the full idempotent cut.

**Proof.** Equation (7) proves entry convergence at each fixed size.
The limiting forms are the polynomial Gram forms with the strictly
positive densities \(1,y,1-y\) relative to (15). A nonzero polynomial
cannot vanish on an interval, so all three forms are strictly positive.
Continuity of eigenvalues gives eventual positivity of the actual
finite forms. For example, if
\(\epsilon_N(a)=\max_{0\le j\le2N-1}|c_j(a)-m_j|\) and
\(\delta_N=\min\{\lambda_{\min}(A_N^*),\lambda_{\min}(B_N^*)\}>0\),
then \(2N\epsilon_N(a)<\delta_N\) suffices for both signs.

The orthonormal polynomials of (15) are
\[
p_0(y)=1,\qquad p_n(y)=\sqrt2\,T_n(2y-1),\quad n\ge1,
\]
because \(T_n(\cos\theta)=\cos(n\theta)\) and the cosine functions
have the corresponding uniform orthogonality. The cosine recurrence
gives the printed tridiagonal matrix for multiplication by \(y\).
Gram--Schmidt coefficients vary continuously with a positive finite
moment matrix, proving convergence of the source-orthonormal
representations.

The degree-\(n\) monic polynomial has squared norm \(2^{1-4n}\)
for \(n\ge1\), and the constant has squared norm one. The product of
these norms for \(0\le n<N\) gives \(\det H_N^*\) in (17).
The characteristic polynomial of \(J_N^*\) is the degree-\(N\)
monic orthogonal polynomial, whose roots are (19).
Its leading normalization and \(T_N(-1)=(-1)^N\) give
\(\det J_N^*=2^{1-2N}\). Symmetry under \(y\mapsto1-y\) gives
the same determinant for \(I-J_N^*\).
Since \(\det A_N^*=\det H_N^*\det J_N^*\), and similarly for \(B\),
the remaining determinants follow.

The equally weighted nodes (19) integrate every polynomial of degree
at most \(2N-1\): after the cosine substitution, the required identities
are \(N^{-1}\sum_j\cos(k\theta_{j,N})=0\) for \(1\le k<2N\).
They follow by pairing for odd \(k\), and a finite geometric sum for
even \(k\). Thus they give the seed spectral measure of \(J_N^*\).
The ILC projection identity follows from commuting functional calculus.
Continuity of the positive square root gives convergence of the
actual cuts to (20).

The diagonal entries of (21) are the mean of \(y\) and \(1-y\), both
\(1/2\). Its off-diagonal entry is
\[
\frac1N\sum_{j=1}^N\sqrt{y_{j,N}(1-y_{j,N})}
=\frac1{2N}\sum_{j=1}^N\sin\theta_{j,N}
=\frac1{2N\sin(\pi/(2N))}.
\]
This proves (21)-(22), including the exact value at \(N=5\). QED.

The infinite limiting model is independently well defined:
on \(L^2(\nu_*)\), let \(Y_*f(y)=yf(y)\), \(e(y)=1\), and
\(S_*=\sqrt{Y_*(I-Y_*)}\). Then
\[
\Pi_*=\begin{pmatrix}Y_*&S_*\\S_*&I-Y_*\end{pmatrix}
\tag{23}
\]
is the same ILC projection. Its two seed readings are directly
\(\int y\,d\nu_*=1/2\) and
\(\int\sqrt{y(1-y)}\,d\nu_*=1/\pi\).
Polynomial compressions of \(Y_*\), padded as in ILC-T3, converge
strongly to \(Y_*\), and their cuts converge strongly to (23).

The order used here is explicit:
first \(a\to\infty\) at each fixed finite \(N\), then the completion
of the limiting moment law as \(N\to\infty\). Constructing the
actual infinite positive cut at one fixed finite anchor remains the
all-order source problem.

The trigonometric representation used in this proof is classical
[DLMF 18.5.1](https://dlmf.nist.gov/18.5#E1).
No priority claim is made for Chebyshev quadrature or these standard
arcsine moment determinants.

## NIF-T4. A cut can be an exponential limit; universal readings have an exact test

For any orthogonal projection \(\Pi\) on a positive Hilbert carrier,
including (20) or (23), define the specified comparison transport
\[
K=2\Pi-I,\qquad \mathcal B_\sigma=e^{\sigma K}.
\]
Its exact formula is
\[
\boxed{\mathcal B_\sigma=e^\sigma\Pi+e^{-\sigma}(I-\Pi).}
\tag{24}
\]
Consequently
\[
\boxed{
e^{-|\sigma|}\mathcal B_\sigma\longrightarrow
\begin{cases}
\Pi,&\sigma\to+\infty,\\
I-\Pi,&\sigma\to-\infty,
\end{cases}
\quad\text{in operator norm},
}
\tag{25}
\]
with norm error at most \(e^{-2|\sigma|}\).

This construction starts with the cut and defines a new comparison
flow from it. It provides an exact exponential realization of the
cut; it does not derive the cut independently or identify this
nonunitary flow with LAT anchor motion, TET source transport, or
the retained unitary wave \(W=(2\Pi-I)Z\).
The parameter \(\sigma\) is correspondingly separate from those
parameters.

For a bounded self-adjoint observable \(O\) and state \(x\) with
\(\Pi x\ne0\), the normalized reading along its positive branch is
\[
\lim_{\sigma\to\infty}
\frac{\langle\mathcal B_\sigma x,O\mathcal B_\sigma x\rangle}
{\|\mathcal B_\sigma x\|^2}
=
\frac{\langle\Pi x,O\Pi x\rangle}{\|\Pi x\|^2}.
\tag{26}
\]
There is one and the same constant \(c\) for every such input state
if and only if
\[
\boxed{\Pi O\Pi=c\Pi.}
\tag{27}
\]

**Proof.** The orthogonal decomposition into the range and kernel
of \(\Pi\) gives eigenvalues \(+1,-1\) for \(K\), proving (24).
For \(\sigma\ge0\), subtracting \(\Pi\) from the left side of (25)
gives exactly \(e^{-2\sigma}(I-\Pi)\). For \(\sigma\le0\),
subtracting \(I-\Pi\) gives \(e^{2\sigma}\Pi\).
This proves the error bound and both limits.

After dividing \(\mathcal B_\sigma x\) by \(e^\sigma\), its components
are \(\Pi x+e^{-2\sigma}(I-\Pi)x\). Normalizing proves (26).
If (27) holds, the quotient is \(c\). Conversely, constancy of
the quotient implies that the quadratic form of
\(\Pi O\Pi-c\Pi\) vanishes on the range of \(\Pi\).
Polarization makes its associated bounded self-adjoint operator zero,
proving (27). QED.

These readings are unchanged by a consistently carried isometric
frame: replacing \((\Pi,O,x)\) by
\((U\Pi U^*,UOU^*,Ux)\) preserves the identities and the quotient.
Changing the physical preparation or changing the observer
independently is a different operation.

The constants in (22) refer to the specified source seed channels.
They are not common readings of every observable and state.
For an explicit example on (23), let
\[
Vf=(\sqrt{Y_*}f,\sqrt{I-Y_*}f),\qquad
O=\operatorname{diag}(Y_*,Y_*).
\]
Both \(V1\) and \(V\sqrt{2y}\) are unit vectors in the cut range.
Their \(O\) readings are, respectively,
\[
\int y\,d\nu_*=\frac12,\qquad
2\int y^2\,d\nu_*=\frac34.
\]
Thus the cut keeps genuine state information. A nontrivial
all-state constant requires (27), not just an infinite parameter.
Physical identification of any such reading as a speed requires
specified spatial and temporal observables and their units.

## NIF-T5. The infinite limit does not retain every finite-source distinction

There are two explicit reasons the positive limiting cut does not
close the finite-anchor RH condition.

First, finite derivative depth and a large response argument need
not commute. For the limiting finite matrix in (19), define
\[
F_N^*(z)=\frac1N\sum_{j=1}^N\frac1{1+zy_{j,N}},\qquad z\ge0.
\]
Then
\[
\boxed{
\lim_{N\to\infty}\lim_{z\to\infty}\sqrt z\,F_N^*(z)=0,\qquad
\lim_{z\to\infty}\lim_{N\to\infty}\sqrt z\,F_N^*(z)=1.
}
\tag{28}
\]
Thus a finite matrix's infinite-argument response is not a substitute
for the infinite-depth response.

Second, finite nonreal zero insertions preserve every limiting
moment and cut established above. Choose \(u>0\) and \(v>0\), and set
\[
P(t)=(t+u)^2+v^2,\qquad
\widetilde{\mathcal X}(t)=\mathcal X(t)P(t).
\]
Then
\[
\widetilde q(t)=q(t)+\frac{2(t+u)}{(t+u)^2+v^2}.
\tag{29}
\]
The added term is \(O_r(a^{-1})\) uniformly for \(t=a(1+z)\),
\(|z|\le r<1\), so the entire conclusions (5)-(8) and their fixed
finite matrix limits remain the same for this modified source.
Nevertheless the modified entire function has the extra nonreal
zeros \(-u+iv\) and \(-u-iv\).

**Proof.** For each fixed \(N\), every \(y_{j,N}\) is strictly positive,
and
\[
0\le\sqrt z\,F_N^*(z)
\le\frac1{N\sqrt z}\sum_{j=1}^N y_{j,N}^{-1}\longrightarrow0.
\]
For fixed \(z\ge0\), the finite sum is a midpoint Riemann sum in
\(\theta\), so it tends to
\(\int(1+zy)^{-1}d\nu_*(y)=(1+z)^{-1/2}\).
Multiplying by \(\sqrt z\) and taking \(z\to\infty\) proves (28).

Logarithmic differentiation proves (29). On the specified disk,
\(|t|\ge a(1-r)\), so the added rational term is uniformly \(O_r(a^{-1})\).
It can therefore be absorbed into the remainder in (4).
The coefficient argument in NIF-T1 and the finite-dimensional
arguments in NIF-T3 apply unchanged.
The roots of \(P\) are the stated nonreal pair. QED.

The comparison source in (29) is not the Riemann xi function and
does not preserve its exact Euler/Gamma identity. It is not an
RH counterexample. It establishes that this asymptotic limit by
itself does not recover even all finite zero modifications.
The original actual source is retained throughout NIF-T1 to T3.

The exact quantifiers supplied by NIF-T3 are
\[
\boxed{
\text{for every finite }N\text{ there exists }a_\star(N)
\text{ such that }A_N(a),B_N(a)>0
\text{ for all }a\ge a_\star(N).
}
\tag{30}
\]
They do not supply a single finite anchor with positivity at every
\(N\). The latter is the retained RH-equivalent condition.
The positive limit and the two seed constants are therefore an
unconditional development, while the all-order actual-source
positivity problem remains open.

## Proof status, source lineage and next obligation

| Result | Status |
|---|---|
| Actual normalized source limit and fixed-order correction, NIF-T1 | Written proof from the complete arithmetic source |
| Classification of stationary LAT probability responses, NIF-T2 | Written proof; the actual source selects the arcsine member asymptotically |
| Eventual positivity at each fixed finite size, NIF-T3 | Written existence proof; no new numerical anchor threshold evaluated |
| Limiting five-direction matrix, cut and seed readings | Exact written derivation |
| Cut-generated exponential limit and constant-reading criterion, NIF-T4 | Written operator identities for the expressly defined flow |
| Noncommuting limits and asymptotic information loss, NIF-T5 | Written exact comparison arguments |
| Infinite positive actual source cut at a fixed finite anchor | Open |
| RH, physical speed of light, Maxwell/Yang--Mills or quantum-gravity closure | Not established by this note |

The retained source lineage is:
[LAT anchor flow](papers/native-layer-anchor-transport/THEOREM.md),
[ILC moment and cut construction](papers/native-infinite-layer-cut/THEOREM.md),
[LOT carried observers](papers/native-layer-observer-transport/THEOREM.md),
[TET metric-compatible source transport](papers/native-tvsp-exponential-transport/THEOREM.md),
and [the finite observer detector](RH_FINITE_OBSERVER_DETECTOR.md).

The next actual-source obligation is to retain and control the
deviation from this fixed point at unbounded derivative depth.
Convergence for every fixed degree is already proved here; the
uniform finite-anchor sign control needed for RH is not.
For a physical constant, the corresponding task is to identify the
observable and show the scalar-on-cut condition (27), with units.

No local scientific computation, independent agent, numerical replay,
or new scientific CI stage was used. This is a written derivation,
not a proof-assistant-checked certificate. Existing source evaluators,
finite certificates and workflow definitions are unchanged.
The documentation commit uses [skip ci], and the change remains a
draft stacked on PR #96.
