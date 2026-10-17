# LAR-1 — Entire source return with retained coefficient error

17 October 2026. Continue LWP-1 from commit
`3a843da9377f1a03e54869e5f3ea22f82ac42499`.
The actual completed source is unchanged:

\[
\mathcal X(w^2)=\xi(1/2+w),\quad q=\mathcal X'/\mathcal X,
\quad c_n(a)=\frac{(-a)^n q^{(n)}(a)}{n!q(a)},\quad a=w^2.
\]

At the retained real anchors, q(a)>0 and c_0=1. Define
\(D_{r,k}=\sum_{\ell=0}^k(-1)^\ell\binom{k}{\ell}c_{r+\ell}\).
The parameter tau remains dimensionless. No physical clock is selected.

## LAR-T1. Unconditional entire analytic return

For w>4, LAA-T1 supplies an analytic source disk and its Cauchy bound:

\[
|c_n|\le M R^{-n},\qquad R=15/16,\qquad
M=Q_w/q(a),\qquad Q_w=41/(2w)+1/4.
\tag{1}
\]

These hypotheses do not require RH or all-order moment positivity.
Define the actual-source scalar function

\[
\boxed{\mathcal B_a(\tau)=\sum_{j=0}^{\infty}
\frac{(-4\tau^2)^j}{(2j)!}D_{j+1,j}.}
\tag{2}
\]

It is an even entire function, real on the real axis. Write
\(C=M/R\), \(A=4(R+1)/R^2=1984/225\). Then

\[
|\mathcal B_a(\tau)|\le C\cosh(\sqrt A\,|\tau|),\quad
\mathcal B_a(0)=c_1,\quad \mathcal B_a''(0)=-4(c_2-c_3).
\tag{3}
\]

**Proof.** The binomial theorem and (1) give

\[
|D_{j+1,j}|\le M R^{-j-1}(1+R^{-1})^j
=C\bigl((R+1)/R^2\bigr)^j.
\]

Thus the absolute j-th term of (2) is at most
\(C(A|\tau|^2)^j/(2j)!\). The sum is the displayed hyperbolic
cosine and is uniformly bounded on every compact tau disk. The
Weierstrass convergence theorem gives an entire sum; evenness, reality
and the derivatives follow from the coefficients. No positive measure
or infinite Hilbert space was used. QED.

This proves existence for the actual infinite source sequence. It does
not construct a positive realization or assert that all finite source
extensions pass their positivity conditions.

## LAR-T2. Coarse exact rational remainder

Let \(S_J=\sum_{j=0}^{J-1}(-4\tau^2)^jD_{j+1,j}/(2j)!\),
\(x=A|\tau|^2\), and choose J>=1 with
\(r_J=x/((2J+1)(2J+2))<1\). Then

\[
\boxed{|\mathcal B_a-S_J|\le E_J(\tau):=
\frac{C x^J}{(2J)!}\frac1{1-r_J}.}
\tag{4}
\]

**Proof.** The ratio of successive absolute majorant terms, beginning
at index J, is at most r_J and decreases thereafter. Sum the geometric
majorant. For any fixed tau, a sufficiently large J meets this condition;
the ratio restriction belongs to this finite bound, not to (2). QED.

For c_0,...,c_19, J=10 gives exactly the LWP prefix S_10, now with an
unconditional analytic remainder. This generally exceeds LWP's radius
\(c_1|\tau|^{20}/20!\), which uses a positive realization. The two
radii have different hypotheses and must not be substituted for each other.

## LAR-T3. Every retained coefficient can contribute

Given c_0,...,c_K with K>=1, define

\[
P_K(\tau)=\sum_{j=0}^{K-1}\frac{(-4\tau^2)^j}{(2j)!}
\sum_{\ell=0}^{\min(j,K-j-1)}
(-1)^\ell\binom j\ell c_{j+1+\ell}.
\tag{5}
\]

Equivalently \(P_K=\sum_{n=1}^K c_n\kappa_n\), where

\[
\kappa_n(\tau)=\sum_{j=\lceil(n-1)/2\rceil}^{n-1}
\frac{(-4\tau^2)^j}{(2j)!}
(-1)^{n-j-1}\binom j{n-j-1}.
\tag{6}
\]

This polynomial has degree at most 2(K-1); it retains every occurrence
of every known coefficient in the infinite series. In particular, K=19
gives degree at most 36 without obtaining derivatives beyond order 19.
P_K is not a Taylor truncation of degree 2K-2: its higher coefficients
retain only the known contributions.

For each j>=0, put

\[
U_{K,j}=M R^{-j-1}
\sum_{\ell=\max(0,K-j)}^j\binom j\ell R^{-\ell},
\tag{7}
\]

where a sum with lower limit above its upper limit is zero. For any
J>=K satisfying the ratio condition in LAR-T2,

\[
\boxed{|\mathcal B_a-P_K|\le
\sum_{j=0}^{J-1}\frac{(4|\tau|^2)^j}{(2j)!}U_{K,j}
+E_J(\tau).}
\tag{8}
\]

**Proof.** A term uses c_n with n=j+1+ell, 0<=ell<=j. The condition
n<=K gives (5); solving for j at fixed n gives (6). The omitted terms
satisfy ell>=K-j. Apply (1) to these terms to obtain (7). At j>=K no
known coefficient occurs, so the entire remaining tail from J is bounded
by (4). The absolute convergence established in T1 justifies the split.
All unknown contributions vanish for 2j+1<=K; at K=19 the first missing
index is j=10. QED.

No assertion is made that refined and coarse interval enclosures are
always nested. They use different centers and interval arithmetic may
have different dependency loss. The certificate compares their radii and
checks consistency at the declared actual readings.

For interval source data, use the lower endpoint of q(a) in M, evaluate
P_K with outward intervals, and add the symmetric radius (8). This
retains both uncertainty in the known data and the entire unknown tail.
The finite consistency check |c_n|<=MR^(-n) cannot establish the bound for
unknown n; that bound comes from the source theorem LAA-T1.

## LAR-T4. Conditional identification with the positive observer return

If the actual sequence admits its inherited probability measure nu on
[0,1], with c_n=integral y^n dnu, then

\[
\mathcal B_a(\tau)=\int_0^1 y
\cos\bigl(2\tau\sqrt{y(1-y)}\bigr)\,d\nu(y).
\tag{9}
\]

This is the LWP/LOT weighted-channel exponential return of e_+=(1,0).
In particular, for real tau, |B_a(tau)|<=c_1, and the sharper positive
Taylor radius of LWP applies to this same analytic function.

**Proof.** Finite binomial expansion gives
D_(j+1,j)=integral y^(j+1)(1-y)^j dnu. On [0,1], 4y(1-y)<=1,
so the cosine series is uniformly absolutely convergent for tau in
compact sets. Interchange summation and integration to obtain (9).
The inherited weighted generator squares to -4y(1-y) on each channel,
which identifies the seed return. The real-axis bound follows from
|cos|<=1. QED.

For the actual completed source, the all-order criterion remains the
one documented in ILC and its RH connection. T1 does not prove this
additional hypothesis. An explicit counterexample shows why: take
c_0=1, c_1=1/2, c_2=-1/4 and c_n=0 for n>=3. The first A_1,B_1
are both positive and (2) is the entire polynomial
B(tau)=1/2+tau^2/2. But det H_2=-1/2 and B(1)=1>c_1. It cannot be
a unitary seed return with squared seed norm c_1. This is a synthetic
counterexample to an inference, not an assertion about xi.

## Remaining obligation

The actual source now has an unconditional entire scalar return with
finite certified error budgets. A global positive observer realization
still needs the actual all-order positivity theorem. A physical mass,
field equation or spacetime interpretation requires further structure.
No theorem here supplies those additional claims.
