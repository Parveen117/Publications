# PTE-1 — Prime transport energy and monotone defect flow

17 October 2026. Written continuation at Publications commit
7d04b9faee1ee78df0d00a68bd6956cedc2708cd, following
[CWF-1](RH_CRITICAL_WAVE_FORMS.md).

The square-root prime normalization has an exact unitary transport
interpretation. A geometric resolvent of each primitive prime combines
all its powers into one norm reading. Substitution into the actual
Guinand-Weil form does not give an unconditional norm-square:
it gives a positive difference-energy minus an explicit scalar term.

This note derives that identity, proves a boundary coercivity estimate,
and enlarges the explicitly covered boundary window to L<=1/32 with
margin at least 1/16. The closed energy has a compact inverse on every
finite window. Its exact positivity test is a compact-operator norm
bound. A monotone coupling then counts all negative directions without
cancellation, while the ordinary Fredholm index stays zero.

**Status: RH remains unproved.** The compact norm bound is established
in the stated small window, not for every L. The new window still lies
before the first prime translation. These are written analytic proofs,
not machine-checked source certificates. No scientific computation,
source evaluation, independent agent or replay is performed.

## Source and conventions

Use CWF's real even test space
\[
\mathcal V_L=
\{g\in C_c^\infty((-L,L),\mathbb R):g(-x)=g(x)\},
\qquad
\widehat g(t)=\int_{\mathbb R}g(x)e^{itx}\,dx.
\]
All functions are extended by zero outside their support. Put
\[
T_u g(x)=g(x-u),\qquad
C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,
\qquad
\ell(g)=\int_{\mathbb R}g(x)\cosh(x/2)\,dx.
\tag{1}
\]
For real even g, C_g(u)=(g*g)(u). It is zero for |u|>=2L.

Retain the actual boundary source form from CWF-T1--T3:
\[
\begin{aligned}
\mathscr Q_\xi(g,g)
={}&\ell(g)^2+
\frac1{4\pi}\int_{\mathbb R}|\widehat g(t)|^2\Omega(t)\,dt
-\sum_{2\le n<e^{2L}}a_n C_g(\log n),\\
\Omega(t)&=\Re\psi(1/4+it/2)-\log\pi,\qquad
a_n=\frac{\Lambda(n)}{\sqrt n}.
\end{aligned}
\tag{2}
\]
Here Lambda(p^k)=log p and Lambda(n)=0 otherwise. There is no RH
assumption in (2). Its positivity for every L and every g in V_L
is equivalent to RH by the retained CWF finite-detector argument.

Define
\[
\begin{gathered}
\Omega_0=\Omega(0)
=-\gamma_{\!E}-\frac{\pi}{2}-3\log2-\log\pi<0,\\
k(u)=\frac{e^{-u/2}}{1-e^{-2u}}\quad(u>0),\qquad
S_L=\sum_{2\le n<e^{2L}}a_n,\qquad
c_L=S_L-\frac{\Omega_0}{2}>0 .
\end{gathered}
\tag{3}
\]
Euler's constant gamma_E is unrelated to a zero ordinate.

| Parameter | Meaning here |
|---|---|
| p, n=p^k | Primitive prime and prime power |
| sqrt(n) | Arithmetic radial relabelling |
| a and sqrt(a) | Retained source anchor and bank-depth scale |
| L | Half-width of a compact wave profile |
| lambda | New diagnostic coupling defined in T5 |
| t | Fourier variable in (2); not physical time |

In particular, sqrt(n) and sqrt(a) are not identified by notation.
The coupling in T5 is not silently identified with the native anchor
displacement or the heat parameter in an index formula.

## PTE-T1. Root-normalized prime transport gives an exact prime-power observer

On L2((0,infinity),dx), for s>0 define
\[
(U_s f)(x)=s^{-1/2}f(x/s),\qquad
(Jf)(u)=e^{u/2}f(e^u).
\]
Then
\[
\boxed{
\|U_s f\|_2=\|f\|_2,\quad
U_sU_r=U_{sr},\quad
U_s^*=U_{1/s},\quad
JU_sJ^{-1}=T_{\log s}.
}
\tag{4}
\]
The map J is unitary from multiplicative-coordinate L2(dx) onto
additive-coordinate L2(du).

For a primitive prime p put r_p=p^(-1/2) and define the bounded operator
\[
\boxed{
V_p=(1-r_p^2)^{1/2}(I-r_p T_{\log p})^{-1}.
}
\tag{5}
\]
For real even compact g,
\[
\boxed{
\|V_p g\|_2^2
=\|g\|_2^2+
2\sum_{k\ge1}p^{-k/2}C_g(k\log p).
}
\tag{6}
\]
Only p^k<e^(2L) can contribute to the correlation sum. Consequently
the entire prime part of (2) is exactly
\[
\boxed{
-\sum_{2\le n<e^{2L}}a_n C_g(\log n)
=\frac12\sum_{p<e^{2L}}\log p
 \bigl(\|g\|_2^2-\|V_p g\|_2^2\bigr).
}
\tag{7}
\]
For every p>=e^(2L), the bracket on the right is zero.

**Proof.** Substituting x=sy proves the norm identity in (4);
composition and inversion prove the next two identities.
Substituting x=e^u proves the norm identity for J, and direct
evaluation of JU_s gives the translation formula.

Since r_p<1 and T_log(p) is unitary, the resolvent in (5) equals
the operator-norm convergent series sum_(j>=0) r_p^j T_(j log p).
In its squared norm, diagonal pairs contribute
||g||^2/(1-r_p^2). Grouping the pairs whose indices differ by k
gives 2 r_p^k C_g(k log p)/(1-r_p^2).
The double series is absolutely convergent because each correlation
has absolute value at most ||g||^2.
Multiplication by 1-r_p^2 proves (6).

The correlation vanishes once the two supports are disjoint.
Thus (6) contains only the stated finite set of nonzero correlations.
Grouping the n sum in (2) by primitive primes proves (7).
The last statement follows by the same support argument. QED.

The operator V_p is a resolvent observation, not a unitary map.
Its norm is bounded by sqrt((1+r_p)/(1-r_p)). Exact transport
normalization therefore does not mean that each bracket in (7)
is nonnegative.

The root coordinate r=sqrt(n) rewrites a_n as Lambda(n)/r and
log n as 2 log r. This is an exact change of coordinates.
The first identity (4) explains the square-root normalization with
respect to the specified measure dx; a different measure changes
that normalization rule.

## PTE-T2. The complete source is positive difference-energy minus a scalar

For g in V_L define
\[
\begin{aligned}
\mathcal E_\Gamma[g]
&=\frac12\int_0^\infty k(u)\,
       \|g-T_u g\|_2^2\,du,\\
\mathcal E_L[g]
&=|\ell(g)|^2+\mathcal E_\Gamma[g]
 +\frac12\sum_{2\le n<e^{2L}}a_n
       \|g-T_{\log n}g\|_2^2.
\end{aligned}
\tag{8}
\]
These quantities are finite and nonnegative. The exact identity is
\[
\boxed{\mathscr Q_\xi(g,g)=\mathcal E_L[g]-c_L\|g\|_2^2.}
\tag{9}
\]
In particular, the map
\[
\mathcal A_L g=
\left(
\ell(g),\
\sqrt{k(u)/2}\,(g-T_u g),\
\left\{\sqrt{a_n/2}\,(g-T_{\log n}g)\right\}_{n<e^{2L}}
\right)
\tag{10}
\]
takes values in the corresponding Hilbert direct sum, satisfies
||A_L g||^2=E_L[g], and gives
\[
\boxed{
\mathscr Q_\xi(g,g)
=\|\mathcal A_Lg\|^2-c_L\|g\|_2^2 .
}
\tag{11}
\]
The desired norm-square factorization of Q itself still requires
a lower bound on A_L, not merely the norm preservation of U_p.

**Proof.** The digamma integral, applied at 1/4+it/2 and 1/4,
gives
\[
\Omega(t)-\Omega_0
=2\int_0^\infty k(u)(1-\cos(tu))\,du.
\tag{12}
\]
The integral converges because k(u) is of order 1/u at zero
and decays exponentially at infinity. Its integrand is nonnegative.

Plancherel, followed by Tonelli on that nonnegative integrand, yields
\[
\frac1{4\pi}\int |\widehat g(t)|^2\Omega(t)\,dt
=\frac{\Omega_0}{2}\|g\|_2^2+\mathcal E_\Gamma[g].
\tag{13}
\]
The factor is fixed by
||g-T_u g||^2=2(||g||^2-C_g(u)).
The same identity gives, term by term,
\[
-a_n C_g(\log n)
=\frac{a_n}{2}\|g-T_{\log n}g\|_2^2-a_n\|g\|_2^2.
\]
Substitution into (2) proves (9).
Near u=0 the difference norm is at most u||g'||_2, so (8)
is finite; its tail is integrable because the difference norm is
at most 2||g||_2. The remaining terms form a finite sum.
Equation (10) now proves (11). QED.

An omitted prime power with log n>=2L would add exactly
a_n||g||^2 to E_L and exactly a_n to c_L, because its two
translated supports are disjoint. Thus extending a finite prime
cutoff beyond the support changes both sides of the subtraction
by the same amount. It does not create positivity.

The analytic input for (12) is
[NIST DLMF 5.9.16](https://dlmf.nist.gov/5.9.E16).
The value of psi(1/4) in (3) follows from
[DLMF 5.4.19](https://dlmf.nist.gov/5.4.E19).
No unpaired infinite difference of norms is used.

## PTE-T3. Boundary coercivity enlarges the proven window to 1/32

Define
\[
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+
 \arctan(e^{-v/2}),\qquad v>0.
\tag{14}
\]
Then for every g in V_L,
\[
\boxed{
\mathcal E_\Gamma[g]\ge F(L)\|g\|_2^2,\qquad
\mathscr Q_\xi(g,g)\ge
 \left(F(L)+\frac{\Omega_0}{2}-S_L\right)\|g\|_2^2
 +|\ell(g)|^2 .
}
\tag{15}
\]
In particular,
\[
\boxed{
0<L\le\frac1{32}
\quad\Longrightarrow\quad
\mathscr Q_\xi(g,g)\ge\frac1{16}\|g\|_2^2.
}
\tag{16}
\]

**Proof.** Substituting z=e^(-u/2) in the integral for F gives
2 integral_0^(e^(-v/2)) dz/(1-z^4), which is the expression in (14).

Write I=(-L,L). Symmetrizing the difference integral gives the
exact boundary decomposition
\[
\begin{aligned}
\mathcal E_\Gamma[g]
={}&\frac14\int_I\int_I
 k(|x-y|)|g(x)-g(y)|^2\,dy\,dx\\
&+\frac12\int_I|g(x)|^2
 \bigl(F(L-x)+F(L+x)\bigr)\,dx .
\end{aligned}
\tag{17}
\]
All expressions are understood through nonnegative integrals.
The diagonal singularity is harmless for smooth g.
The two outside-inside portions of the full double integral
produce the second term in (17).

Since k is positive and strictly decreasing, F''=-k'>0.
Convexity gives
(F(L-x)+F(L+x))/2>=F(L).
Dropping the first nonnegative term proves the first bound in (15).
Using (9) and retaining the nonnegative pole square proves the second.

Here is an elementary bound with no numerical source evaluation.
The inequalities
\[
0<\gamma_E<1,\quad 3<\pi<22/7<16/5,\quad
2/3<\log2<7/10,\quad \log\pi<6/5
\]
give Omega_0/2>-3.
For example, log2>2/3 follows from the positive series
log2=2 sum_(j>=0) (1/3)^(2j+1)/(2j+1).
The upper logarithm bounds follow by positive exponential series:
sum_(j=0)^4 (7/10)^j/j!>2 and
sum_(j=0)^3 (6/5)^j/j!>22/7.
The bound gamma_E<1 follows from its harmonic-sum definition.

For L<=1/32, e^(-L/2)>=1-L/2>=63/64. Hence
\[
\begin{aligned}
F(L)
&\ge\operatorname{artanh}(63/64)+\arctan(63/64)\\
&\ge\frac12\log127+\frac{\pi}{4}-\frac1{64}\\
&>\frac73-\frac1{254}+\frac34-\frac1{64}.
\end{aligned}
\tag{18}
\]
The arctangent bound uses that its derivative is at most one.
For the last line use
log127=7log2-log(128/127)>14/3-1/127.
The displayed rational expression minus 3 is
1/12-1/254-1/64, which is greater than 1/16:
after subtracting 1/16 the remainder is 1/192-1/254>0.

Finally 2L<=1/16<log2, so S_L=0. Equations (15)--(18)
prove (16). QED.

This strengthens the explicitly covered boundary range from
CWF-T4's L<=2kappa, kappa=e^(-9)/(2pi), to L<=1/32.
The earlier margin 1/4 remains available in its smaller window;
the new margin is 1/16 over the larger one.

The interval (16) is still before the first prime threshold
L=(log2)/2. It therefore does not prove an interacting-prime
window. It is an analytic improvement for the actual boundary
form, not an evaluation of prime data.

Also, CWF's convergence is for fixed profiles. Equation (16)
does not by itself enlarge SQR's uniform positive finite-bank
sizes. No new lower constant for n_-(a)/sqrt(a) is claimed.

## PTE-T4. Every finite window has a compact positive operator and an exact norm gate

Let H_L be the Hilbert space of even complex L2 functions supported
in [-L,L]. Complexify the real form by
Q(g_1+ig_2)=Q(g_1)+Q(g_2), with real even g_1,g_2.
The energies below use squared absolute values.

Let D_L be the closure of the even smooth compact tests in the
norm (||g||_2^2+E_Gamma[g])^(1/2).
The form E_L extends to a densely defined closed positive form on
H_L. Let G_L be its associated positive self-adjoint operator:
\[
\mathcal E_L[g]=\|G_L^{1/2}g\|_2^2,\qquad g\in D_L.
\]
Then
\[
\boxed{
G_L\ge F(L)I>0,\qquad
G_L^{-1}\text{ is compact},\qquad
K_L:=c_LG_L^{-1}\text{ is positive and compact}.
}
\tag{19}
\]
The source form extends as the closed form of G_L-c_L I, and
\[
\boxed{
\mathscr Q_\xi\ge0\text{ on }\mathcal V_L
\ \Longleftrightarrow\
\|K_L\|\le1.
}
\tag{20}
\]
The number of negative directions is finite and equals
\[
\boxed{
n_-(G_L-c_L I)
=\#\{\text{eigenvalues of }K_L\text{ greater than }1\},
}
\tag{21}
\]
counting multiplicities. If (20) holds, an actual norm-square
factorization is
\[
\boxed{
\mathscr Q_\xi(g,g)
=\|(I-K_L)^{1/2}G_L^{1/2}g\|_2^2 .
}
\tag{22}
\]
Conversely, a Hilbert-space norm-square factorization for every
test implies (20). The factorization (22) is therefore conditional
on precisely the still-required sign bound.

For the explicitly covered range there is a strict contraction:
\[
\boxed{
0<L\le1/32\quad\Longrightarrow\quad
\|K_L\|\le48/49<1.
}
\tag{23}
\]

**Proof.** In Fourier variables, E_Gamma has the nonnegative multiplier
\[
m(t)=\frac{\Omega(t)-\Omega_0}{2},\qquad
\mathcal E_\Gamma[g]=\frac1{2\pi}\int m(t)|\widehat g(t)|^2\,dt.
\]
Multiplication by sqrt(m) is a closed operator on Fourier L2.
Consequently the smooth-test form is closable; its closure is
exactly the domain D_L just defined. Evenness and support are
closed constraints in L2, and the test space is dense in H_L.

The pole functional is bounded on H_L by Cauchy--Schwarz.
Each prime difference has norm at most 2||g||, and only finitely
many enter. These are bounded form additions to E_Gamma.
They preserve closedness and its domain. The bound in (15)
extends by closure, giving G_L>=F(L)I.

To prove compactness, use the digamma asymptotic
Omega(t)=log(|t|/(2pi))+o(1) as |t| tends to infinity.
It follows either from the retained digamma estimates or
[NIST DLMF 5.11.2](https://dlmf.nist.gov/5.11.E2).
Thus m(t) tends to infinity. A bounded set in D_L has uniformly
small Fourier L2 tails:
\[
\frac1{2\pi}\int_{|t|>R}|\widehat g(t)|^2\,dt
\le
\frac{\mathcal E_\Gamma[g]}{\inf_{|t|>R}m(t)}
\longrightarrow0
\]
uniformly on that bounded set. For fixed R, the map from
L2([-R,R]) to L2([-L,L]) with Fourier kernel e^(-itx)
is Hilbert--Schmidt. The restricted Fourier truncations are
therefore compact and converge uniformly in L2 on the form-norm
unit ball. The embedding D_L into H_L is compact.

The positive closed-form representation gives G_L, and its inverse
takes bounded L2 sets into bounded form-norm sets. The preceding
compact embedding proves G_L^(-1) is compact.

By the spectral theorem, the eigenvalues e_j of G_L are positive,
have finite multiplicity and tend to infinity.
Those of K_L are c_L/e_j. The eigenvalues of G_L-c_L I
are e_j-c_L. This proves (20) and (21), including possible zeros.
Density of smooth tests in D_L ensures that a negative form
direction in the closure is approximated by a negative smooth
test. Thus no extra positivity premise is introduced by completion.

Under (20), I-K_L has a bounded nonnegative square root.
For g in D_L,
\[
\|(I-K_L)^{1/2}G_L^{1/2}g\|^2
=\|G_L^{1/2}g\|^2-c_L\|g\|^2,
\]
proving (22). The converse follows because squared Hilbert norms
are nonnegative.

For L<=1/32, T3 gives G_L>= (c_L+1/16)I.
Here S_L=0 and 0<c_L=-Omega_0/2<3, so
||K_L||<=c_L/(c_L+1/16)<48/49. This proves (23). QED.

CWF's window profile is now the attained spectral bottom
\[
\boxed{
\lambda_\xi(L)
=\min\sigma(G_L)-c_L
=c_L\bigl(\|K_L\|^{-1}-1\bigr).
}
\tag{24}
\]
This uses the smooth tests as a form core; the minimizing vector
need not itself be smooth or compactly supported away from the
endpoints.

The Hilbert completion in this construction follows from the
explicit positive energy (8), with its domain and source matching
proved. Positivity of the actual source form is not used to
construct it. Nothing identifies G_L with a physical Hamiltonian
or its eigenvalues with particle masses.

## PTE-T5. A monotone zero-to-infinity flow counts defects that index equality misses

For a fixed finite L define a diagnostic coupling
\[
\boxed{
\tau(\lambda)=\frac{\lambda}{1+\lambda},\qquad
B_L(\lambda)=I-\tau(\lambda)K_L,\qquad 0\le\lambda<\infty,
}
\tag{25}
\]
and set B_L(infinity)=I-K_L.
It is norm-continuous on the compactified interval and
\[
B_L(0)=I,\qquad
\frac{d}{d\lambda}B_L(\lambda)
=-\frac{K_L}{(1+\lambda)^2}\preceq0.
\tag{26}
\]
Every B_L(lambda), including its limiting endpoint, is Fredholm
and has ordinary Fredholm index zero:
\[
\boxed{
\operatorname{ind}B_L(0)
=\operatorname{ind}B_L(\lambda)
=\operatorname{ind}B_L(\infty)=0.
}
\tag{27}
\]
Nevertheless define the unsigned negative count
\[
d_L(\lambda)=n_-(B_L(\lambda)),\qquad d_L(0)=0.
\]
For lambda>0,
\[
\boxed{
d_L(\lambda)
=\#\{\kappa_j(K_L)>1+1/\lambda\},\qquad
d_L(\infty)=\#\{\kappa_j(K_L)>1\}.
}
\tag{28}
\]
It is finite, nondecreasing, and equals its limiting value for
all sufficiently large lambda. Each eigenvalue kappa_j>1 crosses
zero once, at lambda=1/(kappa_j-1), and only in the downward
direction. An eigenvalue kappa_j=1 reaches zero only at infinity
and never produces a negative direction.

Consequently
\[
\boxed{
\begin{aligned}
\mathrm{RH}
&\Longleftrightarrow
d_L(\infty)=d_L(0)=0\quad\text{for every }L>0\\
&\Longleftrightarrow
d_L(\lambda)=0\quad\text{for every finite }L>0,\ \lambda>0 .
\end{aligned}
}
\tag{29}
\]
Equality of the ordinary indices in (27) imposes no corresponding
restriction: that equality holds even if (29) fails.

For this diagnostic reading, the endpoint identification 0~infinity
is compatible exactly when the two counts agree. Monotonicity and
d_L(0)=0 would then force the count to vanish everywhere. Thus the
faithful endpoint-matching condition contains the remaining arithmetic
positivity problem; imposing the matching is not a proof of it.

**Proof.** Compactness of K_L gives norm continuity and the Fredholm
property of identity minus tau K_L. One can see Fredholmness
directly: only finitely many compact eigenvalues can be close
to 1/tau; on their orthogonal complement I-tau K_L is boundedly
invertible. Self-adjointness identifies kernel and cokernel, so
the Fredholm index is zero even when a kernel is present.

The spectral theorem gives the eigenvalue paths
1-tau(lambda)kappa_j, with remaining essential spectral point 1.
Solving their strict negativity condition proves (28).
Only finitely many kappa_j exceed 1, so the count stabilizes.
Solving for a zero gives the stated crossing parameter and
(26) supplies its direction. The case kappa_j=1 is immediate.

By T4, d_L(infinity)=0 exactly when Q is nonnegative on V_L.
CWF-T3 then proves (29). If the final count is positive, at least
one kappa_j>1 gives a crossing at finite lambda. This also proves
the last equivalence. QED.

If I-K_L is invertible, the usual spectral flow along tau in [0,1],
with upward crossings counted positive, is -d_L(infinity).
If the endpoint has a kernel, (28) specifies the strict-negative
count without assigning an artificial negative sign to endpoint
zero modes.

Thus an equality of endpoint readings can be made faithful to the
RH gate by choosing the unsigned defect count, not an index that
already cancels kernel against cokernel. Its required vanishing
has not been proved for every L. The diagnostic path (25) is a
constructed monotone coupling of the actual window operator;
no identification with the original native lambda evolution
is asserted.

For the standard Fredholm distinction, see Dan Freed,
[Lectures 9 and 10, Corollary 9.17](https://abel.math.harvard.edu/~dafr/M392C-2015/Notes/lecture9.pdf).
The ordinary index is locally constant in a continuous Fredholm
family. Eigenvalue crossings do not generally change that index.
The monotone crossing count in this note is a separate observable.
No elliptic-manifold realization or Atiyah--Singer application to
the zeta source is claimed.

## Source connection, closure boundary and a normalization clarification

For any smooth compact even g, CWF supplies the actual source limit
\[
b^2q(b^2)\mathcal L_{b^2}[P_{b,g}^2]
\longrightarrow \mathscr Q_\xi(g,g),
\qquad
\deg P_{b,g}\le\lfloor Lb/2\rfloor.
\tag{30}
\]
Hence a kappa_j(K_L)>1 would, through smooth form-core approximation,
give a negative actual boundary profile and then a negative finite
native bank at every sufficiently large anchor. This is conditional;
no such actual kappa_j or profile is exhibited here.

Conversely, the larger positive boundary window does not give
operator-norm convergence for all changing finite-bank vectors.
The SQR finite-bank range remains its established range.

For clarity in CWF-T1's second-difference estimate, the complete
scaled inequality is
\[
h_b\sum_k
 \bigl|g(kh_b)-2g((k-1)h_b)+g((k-2)h_b)\bigr|
 |\omega^k|
\le C_g h_b^2 .
\tag{31}
\]
The sum without the prefactor h_b is O_g(h_b); the prefactor
belongs to P_(b,g). There are O_g(1/h_b) terms of size O_g(h_b^2),
and the wave weights are uniformly bounded on their support.
After division by |1-omega|^2, (31) gives exactly CWF's bound (7).
This spells out the normalization in that proof; its final
estimate and tail argument do not change.

The strongest closure obtained here is:
\[
\boxed{
\begin{array}{c}
\text{exact root-normalized prime observers}\\
\Downarrow\\
\mathscr Q_\xi=\mathcal E_L-c_L\|\cdot\|^2\\
\Downarrow\\
\mathscr Q_\xi\ge0
\ \Longleftrightarrow\ \|K_L\|\le1
\ \Longleftrightarrow\ d_L(\infty)=0 .
\end{array}
}
\tag{32}
\]
The norm bound is proved for 0<L<=1/32. Its extension to every
finite L remains the RH problem in this representation.

| Established | Not established by it |
|---|---|
| Exact prime resolvent norm including all powers | Nonnegativity of each prime contribution |
| Positive difference-energy with explicit scalar subtraction | Positivity after the subtraction at every width |
| Boundary margin 1/16 for L<=1/32 | Positivity at or beyond the first prime translation |
| Compact positive K_L and finite negative count per window | A bound ||K_L||<=1 for every L |
| Monotone crossing count without cancellation | Vanishing of that count for all actual windows |
| Exact ordinary-index equality at both endpoints | A sign certificate for the source |

## Attribution and execution

The Guinand-Weil formula and Weil criterion remain attributed in
[CWF](RH_CRITICAL_WAVE_FORMS.md). The digamma identities and
asymptotic used here are cited above. Unitary dilation, the
geometric resolvent identity, the closed-form construction and
the compact spectral test are standard mathematical tools;
no literature-wide novelty claim is made for them.

The programme-specific statements here are the exact substitution
into the retained source form, the explicit boundary margin,
and the resulting source-faithful defect count. All proof
assumptions and remaining sign conditions are stated.

Retained programme references:
[CWF](RH_CRITICAL_WAVE_FORMS.md),
[BRR](RH_BOUNDARY_OBSERVER_RETURN.md),
[SQR](RH_SQUARE_ROOT_BANKS.md),
[finite observer detector](RH_FINITE_OBSERVER_DETECTOR.md).

Only this note and its README entry are added. No local scientific
computation, numerical source evaluation, independent agent,
scientific replay or workflow change is made. These written
proofs are not newly machine-checked certificates. The
documentation commit includes [skip ci], is stacked on PR #105,
and is not merged.
