# QRI-1 — Quadrant signs, reflection ratios and growing RH banks

17 October 2026. Written continuation of NT, TLM, VSR, NIF and MOR,
based on Publications commit 6e82f58299f5bf30e72cc15169d0affbb717a461.

The four thermodynamic Maxwell charts give one consistent sign ledger.
Their capacity ratio is a normalized Schur complement, so reversing a
cross-channel orientation leaves that ratio unchanged. The same
algebra applies to an actual RH observer bank, with its own source
pairing and explicitly selected observer directions.

The new actual-source result is stronger than fixed-size convergence:
without assuming RH, the unshifted Hankel banks are positive at every
sufficiently large anchor through a size proportional to log a. The
proof retains the first logarithmic source correction. In an arcsine
Chebyshev observer basis this correction couples only opposite
reflection parities. Its sign disappears from the squared-coupling
ratio, while its magnitude still controls the positivity condition.

No thermodynamic observable is identified with a zeta observable.
The two settings share a proved finite-form identity. The actual
all-order RH gate remains the source-specific condition.

## QRI-T1. All four Maxwell charts and their positive source metrics

Use dimensional source coordinates (S,V), with fixed material amount,
T=U_S, P=-U_V, and
\[
\mathsf H=D^2U=\begin{pmatrix}a&b\\b&c\end{pmatrix},
\qquad a>0,\quad \Delta=ac-b^2>0.
\tag{1}
\]
Thus c>0. Affine calibrations in the native source do not change
these second derivatives. Use E=U+PV for enthalpy, F=U-TS and
G=U+PV-TS. The complete chart ledger is:

| Potential and controls | Exact differential | Maxwell relation |
|---|---|---|
| U(S,V) | dU=T dS-P dV | T_V at S = -P_S at V |
| E(S,P) | dE=T dS+V dP | T_P at S = V_S at P |
| F(T,V) | dF=-S dT-P dV | S_V at T = P_T at V |
| G(T,P) | dG=-S dT+V dP | S_P at T = -V_T at P |

In those respective control orders their Hessians are
\[
D^2U=\begin{pmatrix}a&b\\b&c\end{pmatrix},\qquad
D^2E=\frac1c\begin{pmatrix}\Delta&-b\\-b&-1\end{pmatrix},
\]
\[
D^2F=\frac1a\begin{pmatrix}-1&b\\b&\Delta\end{pmatrix},\qquad
D^2G=-\frac1\Delta\begin{pmatrix}c&b\\b&a\end{pmatrix}.
\tag{2}
\]
Their signatures differ. The source metric transported to the
same four charts, however, is respectively
\[
\mathsf H,\qquad
g_E=\operatorname{diag}(\Delta/c,1/c),\qquad
g_F=\operatorname{diag}(1/a,\Delta/a),\qquad
g_G=\frac1\Delta\begin{pmatrix}c&b\\b&a\end{pmatrix}.
\tag{3}
\]
All four metrics in (3) are positive.

**Proof.** The one common derivative system is
\[
dT=a\,dS+b\,dV,\qquad dP=-b\,dS-c\,dV.
\tag{4}
\]
Solving it in each chosen pair of controls gives (2). Equality of
each potential's mixed derivatives gives the corresponding Maxwell
relation, including its displayed sign.

For completeness, the coordinate Jacobians back to (S,V) are
\[
Q_E=\begin{pmatrix}1&0\\-b/c&-1/c\end{pmatrix},\quad
Q_F=\begin{pmatrix}1/a&-b/a\\0&1\end{pmatrix},\quad
Q_G=\frac1\Delta\begin{pmatrix}c&b\\-b&-a\end{pmatrix}.
\]
Multiplication of Q^T H Q proves (3). In particular a raw Legendre
Hessian is not, in every chart, the positive transported metric.
QED.

The established VTPS order is (V,T,P,S). Reordering the Helmholtz
controls to (V,T) and using the output (P,S) gives exactly VSR-T2:
\[
\mathsf M_{\rm VTPS}
=\frac{\partial(P,S)}{\partial(V,T)}
=\frac1a\begin{pmatrix}-\Delta&-b\\-b&1\end{pmatrix},
\qquad
g_{\rm VTPS}=\operatorname{diag}(\Delta/a,1/a).
\tag{5}
\]
Thus its negative determinant is compatible with a stable source;
it is the response matrix, while g is the acquired metric.
All four charts have now been included without assigning an
unprovided geometric quadrant numbering to the compass.

These are the classical Maxwell and Legendre identities, rederived
from the retained native source. The four relations and the capacity/
compressibility identity are also stated in the Oxford
[Basic Thermodynamics, Handout 6](https://www-thphys.physics.ox.ac.uk/people/AlexanderSchekochihin/A1/2012/handout6_2012.pdf).
Their use here does not establish electromagnetic field equations.

## QRI-T2. What cancels in the capacity ratio, including the five-reading source

At a state with positive temperature and the retained nonzero
entropy/volume calibrations, NT-1 gives
\[
C_V=T/a,\quad C_P=Tc/\Delta,\quad
K_S=Vc,\quad K_T=V\Delta/a.
\]
Here K denotes a bulk modulus, the reciprocal of compressibility.
Consequently
\[
\boxed{
\Gamma_{\rm th}
=\frac{C_P}{C_V}
=\frac{K_S}{K_T}
=\frac{\lambda_V}{\lambda_P}
=\frac{\lambda_S}{z_T}
=\frac{ac}{\Delta}
=\frac1{1-\chi},\qquad
\chi=\frac{b^2}{ac}\in[0,1).
}
\tag{6}
\]
In particular (C_P/C_V)(K_T/K_S)=1 exactly.

Let D=diag(epsilon_S,epsilon_V), with each epsilon equal to +1 or
-1. Changing the two source-axis orientations by congruence sends
H to D H D: a,c stay fixed and b becomes epsilon_S epsilon_V b.
Both chi and Gamma_th are unchanged. Exchanging the two selected
directions likewise preserves this scalar algebraic ratio.

More intrinsically, for two linearly independent selected directions
e,f in a positive form h,
\[
\Gamma(e,f)=
\frac{h(e,e)h(f,f)}
     {h(e,e)h(f,f)-h(e,f)^2}.
\tag{7}
\]
It is invariant under any passive frame change that carries both
directions and the form, and under independent nonzero rescalings
of e,f. Changing the selected physical directions can change it.

There is also an exact cut interpretation. In the original source
coordinates put
\[
P_P=\begin{pmatrix}1&0\\-b/c&0\end{pmatrix}.
\]
Then
\[
P_P^2=P_P,\quad H P_P=P_P^T H,\quad
\boxed{
\frac{\|P_P e_S\|_H^2}{\|e_S\|_H^2}
=\frac{\Delta}{ac}=\Gamma_{\rm th}^{-1}.
}
\tag{8}
\]
Its range satisfies dP=0; its kernel is the volume direction.

**Proof.** Equation (4) gives the four response coefficients, hence
(6). The two orientation factors in b disappear in b^2, proving
the stated sign cancellation. The contractions in (7) are unchanged
under carried frames and have the same total scaling in numerator
and denominator. Direct multiplication gives
H P_P=diag(Delta/c,0); evaluation on e_S proves (8). QED.

On the unchanged TLM five-reading source,
\[
a=2+\mu v,\quad b=1+\mu s,\quad c=5,\quad
\Gamma_\mu(s,v)=
\frac{5(2+\mu v)}{5(2+\mu v)-(1+\mu s)^2}.
\tag{9}
\]
At the calibrated cut,
\[
H_0=\begin{pmatrix}2&1\\1&5\end{pmatrix},\quad
\mathsf M_{\rm VTPS,0}
=\frac12\begin{pmatrix}-9&-1\\-1&1\end{pmatrix},\quad
\Gamma_\mu(0)=10/9.
\]
The retained cut curvature is 5 mu^2/324, so the common ratio
10/9 does not by itself determine the source curvature.
The source-specific reflection from VSR gives
\[
\Gamma_{-\mu}(s,v)=\Gamma_\mu(-s,-v),\qquad
\partial_s\Gamma_\mu(0)=20\mu/81,\quad
\partial_v\Gamma_\mu(0)=-5\mu/81.
\tag{10}
\]
Thus the value is orientation invariant, while it can change with
the material state. The closure product one does not force the
individual ratio to be a constant function.

At every finite jet depth the same identity holds with native
coefficient convolution:
\[
j^N\Gamma=(j^Na\star j^Nc)\star
\bigl(j^Na\star j^Nc-j^Nb\star j^Nb\bigr)^{-1}.
\tag{11}
\]
The inverse exists because its constant coefficient is Delta>0.
For the paired source reflection, degree-alpha coefficients change
by (-1)^|alpha|, and the regular jet matrices are conjugated by
the retained graded sign matrix. Odd derivatives need not vanish.
This is an all-finite-depth algebra identity, not a separate
positive self-adjoint interpretation of the regular jet matrix.

An observer covector and a state transform as O'=OD and x'=Dx,
so O'x'=Ox. Exchanging a scalar pairing's two roles does not
automatically add a minus sign. Antisymmetric oriented forms do
acquire a minus sign on exchange. The Maxwell signs above come
from the specified potential differentials; an observed positive
reading by itself does not establish an underlying negative one.

Maxwell consistency alone also does not imply stability. For
example the integrable constant Hessian [[1,2],[2,1]] obeys all four
mixed-derivative identities wherever those charts are defined,
but Delta=-3 and Gamma=-1/3. This algebraic comparison is not a
replacement for the stable native material source.

## QRI-T3. The reflection-invariant gate for an RH observer bank

For any real symmetric finite form split into two selected banks,
\[
\mathsf G=\begin{pmatrix}\mathsf A&\mathsf C\\
                         \mathsf C^T&\mathsf B\end{pmatrix},
\qquad \mathsf A>0,\quad\mathsf B>0,
\]
define its normalized cross coupling
\[
\mathsf K=\mathsf A^{-1/2}\mathsf C\mathsf B^{-1/2}.
\]
Then
\[
\boxed{
\mathsf G>0\ \Longleftrightarrow\ \|\mathsf K\|_2<1,\qquad
\frac{\det\mathsf G}{\det\mathsf A\,\det\mathsf B}
=\det(I-\mathsf K^T\mathsf K).
}
\tag{12}
\]
When the strict inequality holds, the normalized inverse Schur form is
\[
\mathscr G=(I-\mathsf K^T\mathsf K)^{-1}.
\tag{13}
\]
Its eigenvalues are 1/(1-sigma_j^2), where sigma_j are the singular
values of K, including any zero values on the B side.

Reflecting either bank replaces C by -C and K by -K. Equations
(12)-(13) are unchanged. Independent invertible basis changes within
the selected banks preserve the singular values, since their squares
are the generalized eigenvalues of
C^T A^-1 C versus B. Common positive scaling of the whole form also
cancels. An arbitrary change of the selected subspaces is a different
observation and need not preserve these numbers.

**Proof.** Congruence by diag(A^-1/2,B^-1/2) gives [[I,K],[K^T,I]].
Its Schur complement is I-K^T K. Completing the square proves
positivity equivalence; block elimination proves the determinant
formula. The reflection and generalized-eigenvalue statements follow
by substitution. QED.

For the actual RH source retain
\[
\mathcal X(w^2)=\xi(1/2+w),\quad q=\mathcal X'/\mathcal X,\quad
a_0=3249/4,\quad
c_j(a)=\frac{(-a)^j q^{(j)}(a)}{j!\,q(a)},\quad
H_N(a)=(c_{i+j}(a))_{0\le i,j<N}.
\tag{14}
\]
Write L_a[y^j]=c_j(a) for its polynomial functional; this symbol
does not assume that L_a is positive. At a step where H_n(a)>0
and c_(2n)(a)>0, put
\[
h_n=(c_n,\ldots,c_{2n-1})^T,\qquad
s_n=c_{2n}-h_n^T H_n^{-1}h_n,
\]
\[
\boxed{
\chi_n=\frac{h_n^T H_n^{-1}h_n}{c_{2n}},\qquad
\Gamma_n=\frac{c_{2n}}{s_n}=\frac1{1-\chi_n}.
}
\tag{15}
\]
The next bank is positive exactly when chi_n<1; then Gamma_n>=1.
If chi_n=1 the Schur innovation is zero and the ratio is undefined.
If chi_n>1 the enlarged bank has a negative direction.

Indeed the new residual vector is (-H_n^-1 h_n,1), and
\[
s_n=\min_{\deg p<n}L_a[(y^n-p(y))^2].
\tag{16}
\]
The minimum is valid because the old quadratic part H_n is
positive, even before the new innovation's sign is known.
Once the enlarged bank is positive, its orthogonal projection
onto this residual along the old bank has squared retention
s_n/c_(2n)=Gamma_n^-1 on the new raw monomial direction. This is
precisely the finite-form identity in (8), with different source
and observer data.

At a_0 the preceding positive-layer sector result gives c_j>0
at every depth. Starting at H_1=(1), acceptance of (15) at every
successive n is therefore equivalent to positivity of every H_N,
and hence to RH by the retained detector D3. The universal acceptance
statement is still to be proved.

MOR-T2 gives a further exact invariance: its transported observers
satisfy
\[
f_\rho L_a[(O_\rho p)(O_\rho h)]=L_{a_0}[ph],
\qquad a=\rho a_0,\quad f_\rho>0.
\tag{17}
\]
Hence all the ratios in (12)-(15) are preserved for a co-transported
selected bank, when they are defined. The finite polynomial
replacement has the explicit MOR-T4 error. A fresh monomial bank at
the new anchor is a different selected bank; it is not covered by
an equality that forgets the observer transport.

## QRI-T4. The actual infinity correction is reflection odd

Retain the complete-source expansion in NIF-T1 and set
\[
L=\log\frac{\sqrt a}{2\pi},\qquad
F_a(z)=\frac{q(a(1+z))}{q(a)}.
\]
Here L is a scalar logarithm, distinct from the functional L_a.
For large a let
\[
\mathcal R_a(z)=F_a(z)-(1+z)^{-1/2}
-\frac{\log(1+z)}{2L\sqrt{1+z}},\qquad
E(a)=\sup_{|z|\le1/2}|\mathcal R_a(z)|.
\]
NIF-T1 supplies
\[
\boxed{E(a)=O((\sqrt a\,L)^{-1}).}
\tag{18}
\]
Its proof retains the full rational, Gamma and prime-power source.

For the arcsine reference measure and its orthonormal observers set
\[
d\nu_*(y)=\frac{dy}{\pi\sqrt{y(1-y)}},\qquad
\phi_0(y)=1,\quad
\phi_j(y)=\sqrt2\,T_j(2y-1)\quad(j\ge1).
\tag{19}
\]
Let Q_N contain their monomial coefficient columns. Define
\[
\widehat H_N=Q_N^T H_N(a)Q_N,\qquad
g(y)=\log\frac{1-y}{y},\qquad
(\mathsf D_N)_{ij}=\int\phi_i\phi_j g\,d\nu_*.
\]
Then
\[
\boxed{
\widehat H_N=I+\frac{\mathsf D_N}{2L}+\mathsf R_N,\qquad
\|\mathsf R_N\|_2\le
\varepsilon_N(a):=\frac{2}{99}\,100^N E(a).
}
\tag{20}
\]
This bound holds at every finite N for which (18)'s analytic disk
is in the regular large-anchor regime; its explicit N dependence
will allow N to grow.

The correction has the exact entries
\[
D_{ij}=0\quad\text{if }i,j\text{ have the same parity};
\qquad D_{0j}=D_{j0}=-2\sqrt2/j\quad(j\text{ odd}),
\]
\[
D_{ij}=-2\left(\frac1{|i-j|}+\frac1{i+j}\right)
\quad(i,j\ge1\text{ of opposite parity}).
\tag{21}
\]
In particular, for J_N=diag((-1)^j),
\[
\boxed{
J_N\mathsf D_NJ_N=-\mathsf D_N,\qquad
\|\mathsf D_N\|_2\le9(1+\log N).
}
\tag{22}
\]
Here J_N denotes this parity reflection, not the Jacobi matrix
used for the limiting multiplication operator in NIF.

**Proof.** For 0<alpha<1 the beta family in NIF-T2 has moments
(alpha)_n/n! and density
y^(alpha-1)(1-y)^(-alpha)/B(alpha,1-alpha). Differentiating
its Stieltjes response (1+z)^(-alpha) at alpha=1/2 shows
\[
\int\frac{g(y)}{1+zy}\,d\nu_*(y)
=\frac{\log(1+z)}{\sqrt{1+z}}.
\tag{23}
\]
Differentiation is justified by integrable endpoint logarithms
for alpha in a compact subinterval of (0,1). Therefore Cauchy's
coefficient formula on |z|=1/2 bounds every residual moment by
\[
\left|c_n-\int y^n\,d\nu_*
-\frac1{2L}\int y^n g\,d\nu_*\right|\le E(a)\,2^n.
\tag{24}
\]

For a polynomial p=sum p_k y^k put B_2(p)=sum |p_k|2^k.
All zeros of phi_j, j>=1, lie in (0,1); its coefficients alternate.
Consequently
\[
B_2(\phi_j)=|\phi_j(-2)|=\sqrt2\,T_j(5)
\le\sqrt2\,10^j,\qquad B_2(\phi_0)=1.
\]
The last inequality follows from
T_j(5)=[(5+sqrt(24))^j+(5-sqrt(24))^j]/2 and 5+sqrt(24)<10.
It follows from (24) that
\[
|(\mathsf R_N)_{ij}|\le E(a)B_2(\phi_i)B_2(\phi_j).
\]
For a real vector t, apply Cauchy--Schwarz to
(sum |t_i| B_2(phi_i))^2. Since
\[
\sum_{j<N}B_2(\phi_j)^2
\le1+2\sum_{j=1}^{N-1}100^j
<\frac{2}{99}100^N,
\]
this proves (20).

To compute D, write y=(1+cos theta)/2. The reference measure is
dtheta/pi and g=2 log tan(theta/2). Subtracting the real parts of
the logarithmic power series for 1-r exp(i theta) and
1+r exp(i theta), then taking their integrable Abel limit, gives
\[
g=-4\sum_{\substack{k\ge1\\k\ {\rm odd}}}
       \frac{\cos(k\theta)}k
\]
in the sense of its integrals against any fixed trigonometric
polynomial. Equivalently the regularized series converges to g in
L^1(0,pi); the logarithmic endpoint singularities are integrable.
Cosine orthogonality and
2 cos(i theta)cos(j theta)=cos((i-j)theta)+cos((i+j)theta)
give (21).

For the norm estimate put h_N=sum_(k=1)^N 1/k. Row zero has
absolute sum at most 3h_N. Any other row has the zero-column
contribution at most 3, the |i-j| contribution at most 4h_N,
and the i+j contribution at most 2h_N. Thus every row sum is
at most 9h_N<=9(1+log N). Symmetry bounds the spectral norm
by the maximum absolute row sum, proving (22). QED.

The classical beta and Chebyshev identities used in this proof are
[DLMF 5.12](https://dlmf.nist.gov/5.12) and
[DLMF 18.5.1](https://dlmf.nist.gov/18.5#E1); the logarithmic power
series is [DLMF 4.6.1](https://dlmf.nist.gov/4.6#E1).
The residual estimate and the finite reflection matrix are derived
above for the retained actual source.

The basis reflection satisfies phi_j(1-y)=(-1)^j phi_j(y).
It is a passive reflection of the selected observer basis.
It does not assert L_a[p(y)]=L_a[p(1-y)] for the finite-anchor
actual source. Only the limiting arcsine law has that symmetry;
the displayed first correction is antisymmetric under it.
Nor is this y reflection identified with VSR's material-parameter
reflection without an additional adapter.

## QRI-T5. Positive banks of logarithmically increasing size

Define the explicit bounds
\[
b_N(a)=\frac{9(1+\log N)}{2L},\qquad
d_N(a)=b_N(a)+\varepsilon_N(a).
\tag{25}
\]
If d_N(a)<1, then
\[
\boxed{
(1-d_N)H_N^*\preceq H_N(a)\preceq(1+d_N)H_N^*,
\qquad H_N(a)>0,
}
\tag{26}
\]
where H_N^*=(int y^(i+j) dnu_*) is the arcsine monomial bank.

For each fixed 0<eta<1, put
\[
\boxed{
N_\eta(a)=
\left\lfloor\frac{(1-\eta)\log a}{4\log10}\right\rfloor.
}
\tag{27}
\]
For all sufficiently large real a, every bank with
1<=N<=N_eta(a) is positive definite, without assuming RH.
Uniformly over this range,
\[
d_N(a)=
O_\eta\!\left(\frac{\log\log a}{\log a}
+\frac{a^{-\eta/2}}{\log a}\right)\longrightarrow0.
\tag{28}
\]

**Proof.** Equations (20) and (22) give
||Hhat_N-I||<=d_N. This is equivalent under the invertible Q_N
to (26), since Q_N^T H_N^* Q_N=I. For (27),
\[
100^{N_\eta(a)}\le a^{(1-\eta)/2}.
\]
Combine this with (18), and use L asymptotic to (log a)/2.
The b_N term is O_eta(log log a/log a), while the residual
term is O_eta(a^(-eta/2)/log a). Both bounds increase with N,
so they hold uniformly for all smaller banks. Their sum is
eventually less than one. QED.

This advances NIF's fixed-N statement to an explicit unbounded
growth schedule. The large-anchor threshold has not been
numerically evaluated, and the constant in (27) is not claimed
optimal. Equation (26) is a comparison of forms; it does not
assert an N-independent Euclidean eigenvalue lower bound.

There is a sharper sign-cancellation statement for the reflection
split itself. Reorder the Chebyshev bank into its even and odd
indices. In this order (20) reads
\[
\widehat H_N=
\begin{pmatrix}I+R_{ee}&D_{eo}/(2L)+R_{eo}\\
D_{eo}^T/(2L)+R_{eo}^T&I+R_{oo}\end{pmatrix}.
\tag{29}
\]
For N>=2, if epsilon_N<1, the two diagonal forms are positive.
Their normalized cross coupling in QRI-T3 satisfies
\[
\|\mathsf K_N\|_2\le
k_N:=\frac{b_N+\varepsilon_N}{1-\varepsilon_N}.
\]
If also k_N<1, its ratio form obeys
\[
\boxed{
I\preceq\mathscr G_N\preceq\frac1{1-k_N^2}I.
}
\tag{30}
\]
Along (27) these conditions hold eventually, and
\[
\boxed{
\|\mathscr G_N-I\|_2=
O_\eta\!\left((\log\log a/\log a)^2\right)
}
\tag{31}
\]
uniformly for 2<=N<=N_eta(a).

Indeed ||R_ee||,||R_oo||,||R_eo||<=epsilon_N, and the inverse
square roots of the two diagonal forms have norm at most
(1-epsilon_N)^(-1/2). This proves the cross bound; (12)-(13)
give (30). Equation (28) then proves (31). This is where the
first-order orientation sign cancels and the coupling enters
quadratically. The remaining inequality is its magnitude bound.

For example the first two Chebyshev directions give
\[
\widehat H_2(a)=
\begin{pmatrix}1&-\sqrt2/L\\-\sqrt2/L&1\end{pmatrix}
+O((\sqrt a\,L)^{-1}).
\tag{32}
\]
Reflecting either one of these directions flips the off-diagonal
sign. Their scalar inverse Schur ratio is
\[
\Gamma_{eo}(a)
=1+\frac2{L^2}
+O\!\left(L^{-4}+(\sqrt a\,L)^{-1}\right).
\tag{33}
\]
This is an actual asymptotic RH-source ratio; it is not a physical
heat capacity. It follows from (21), (20), and division by the two
positive diagonal entries, followed by 1/(1-chi).

## QRI-T6. A ratio depends on the selected layer, even at the fixed point

For the monomial extension in (15), the limiting arcsine values are
\[
\boxed{
s_n^*=2^{1-4n},\qquad
\Gamma_n^*=\frac12\binom{4n}{2n}\quad(n\ge1).
}
\tag{34}
\]
Thus Gamma_1^*=3, Gamma_2^*=35, and in general these layer ratios
increase. They are invariant under reflection of the same selected
directions, but they are not one common number for all depths.

More quantitatively, if (26) holds at size N with d_N<1, then for
each 1<=n<N,
\[
\boxed{
\frac{1-d_N}{1+d_N}\Gamma_n^*
\le\Gamma_n(a)\le
\frac{1+d_N}{1-d_N}\Gamma_n^*.
}
\tag{35}
\]
In particular the relative comparison tends to one uniformly through
n<N_eta(a).

**Proof.** For n>=1 the normalized polynomial
phi_n=sqrt2 T_n(2y-1) has leading coefficient 2^(2n-1/2).
Dividing by that coefficient produces the monic orthogonal
residual of y^n. Its squared norm is 2^(1-4n), proving the
first formula. The raw monomial has squared norm
m_(2n)=binom(4n,2n)/2^(4n), proving the second.

For (35), apply (26) first to y^n and then to every polynomial
y^n-p, deg p<n. Taking minima and using (16) gives
(1-d_N)s_n^*<=s_n(a)<=(1+d_N)s_n^* and the analogous
bounds for c_(2n). Division proves (35). QED.

The even/odd bank ratio in (31), the monomial extension ratio in
(34), and the material ratio 10/9 select different directions or
subspaces. A passive frame change carrying those selections preserves
each one's value. Replacing the selection does not. This distinction
specifies what an invariant ratio means in the programme.

## QRI-T7. A quantified depth interval for any hypothetical RH obstruction

If RH fails, retain MOR's definition
\[
n_-(a)=\min\{N:H_N(a)\not\succeq0\}.
\]
The detector and MOR guarantee that this number is finite at
every a>=a_0. QRI-T5 now strengthens its previous merely divergent
lower behavior to
\[
\boxed{
\liminf_{a\to\infty}\frac{n_-(a)}{\log a}
\ge\frac1{4\log10}.
}
\tag{36}
\]
For any one fixed negative polynomial witness supplied by the
conditional detector, MOR-T5 still gives
\[
\boxed{
n_-(\rho a_0)\le\mathcal C(p)\rho+\deg p+2,\qquad \rho\ge1.
}
\tag{37}
\]
Both statements are conditional descriptions of where a hypothetical
failure would remain observable; no such witness has been found.

**Proof.** For each fixed eta in (0,1), QRI-T5 gives
n_-(a)>N_eta(a) eventually. Divide by log a, take the lower
limit, and then let eta decrease to zero to obtain (36).
Equation (37) is the unchanged MOR-T5 witness-return bound. QED.

The logarithmic positive-bank schedule in (27) has
N_eta(a)/(a/a_0)->0. It therefore does not meet MOR-T6's sufficient
schedule M/(a/a_0)->infinity. In particular it does not prove the
existing equivalent target
\[
H_{\lceil m\log(2m)\rceil}(m a_0)\succeq0
\quad\hbox{for arbitrarily large integers }m.
\tag{38}
\]
That growing-bank positivity, or acceptance of every actual
fixed-anchor Schur step, is still the RH obligation.
The relevant obstruction is a cross-coupling magnitude or a
nonpositive innovation. Reversing observer signs preserves that
verdict by real congruence.

## Evidence, attribution and continuation

| Statement | Evidence or status |
|---|---|
| Four Maxwell charts and transported metrics | QRI-T1; classical identities from the retained material source |
| Capacity ratio, cut retention and finite-jet reflection | QRI-T2; explicit source algebra |
| Reflection-invariant RH bank gate and carried-ratio identity | QRI-T3; Schur elimination and MOR-T2 |
| Reflection-odd first actual-source correction with uniform matrix error | QRI-T4; NIF full-source remainder and derived Chebyshev bounds |
| Unconditional actual H_N positivity through logarithmic growth in a | QRI-T5; written uniform proof, threshold not numerically evaluated |
| Quadratic even/odd ratio correction and layer-dependent monomial ratios | QRI-T5 to T6 |
| Logarithmic lower and linear upper depth bounds if RH fails | QRI-T7 and retained conditional detector |
| Fixed-anchor all-order positivity or the MOR sufficient growing schedule | Open |
| Physical identification of the RH ratio with C_P/C_V, Maxwell waves, or gravity | Not established |

Source lineage:
[NT](papers/native-thermodynamic-curvature/THEOREM.md),
[TLM](papers/native-tvsp-lambda-matrices/THEOREM.md),
[VSR](papers/native-vtps-source-reflection/THEOREM.md),
[LCW](papers/native-layer-continuation-window/THEOREM.md),
[actual arithmetic source](RH_ARITHMETIC_GATE.md),
[finite detector](RH_FINITE_OBSERVER_DETECTOR.md),
[NIF](RH_INFINITY_CUT_FIXED_POINT.md), and
[MOR](RH_MOVING_OBSERVER_RETURN.md).

Schur complements, congruence invariance, beta differentiation,
Chebyshev polynomials and Cauchy coefficient bounds are classical
mathematics. Their necessary steps are given here; no general
priority or optimality claim is made. The new programme development
is the actual-source parity correction, its explicit growing-bank
error bound, and the resulting quantified positivity and ratio laws.

Only this note and its README entry are added. No local scientific
computation, agent audit, new source evaluation or scientific replay
is performed. These are written proofs, not new machine-checked
certificates. Existing engines, source data, certificates and
workflows remain unchanged. The documentation commit includes
[skip ci]. The change is a draft stacked on PR #98 and is not merged.
