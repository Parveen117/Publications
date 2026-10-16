# Complete native higher-layer evaluation

The source is exactly PIW/SE48's q(t), with t=(s-1/2)^2. Only
the existing native rational F type, the inherited 110-place outward
interval type and exact integer recurrences perform scientific arithmetic.
Python 3.12 handles file I/O and process orchestration only.

## 1. Higher source coefficients at s=29

Expand Q(s+z)=xi'(s+z)/xi(s+z)=sum_k Q_k z^k. From the
completed xi definition,

\[
Q_k=(-1)^k(s^{-k-1}+(s-1)^{-k-1})+G_k+E_k.
\]

The constant G_0=(psi(s/2)-log pi)/2 uses PIW's convergent log,
Machin pi and signed digamma-remainder implementation unchanged.
For k>=1, differentiating the absolutely convergent polygamma series gives

\[
G_k=\frac{(-1)^{k+1}}{2^{k+1}}\zeta(k+1,s/2).
\tag{1}
\]

The named zeta in (1) is the positive-real Hurwitz series; it is
computed below, not imported from a floating-point special-function library.

## 2. Hurwitz sum and a full Euler--Maclaurin remainder

For integer r>=2, v>0, shift M=512, y=v+M and K=32, compute

\[
Z=\sum_{j=0}^{M-1}(v+j)^{-r}+
\frac{y^{1-r}}{r-1}+\frac12 y^{-r}
+\sum_{j=1}^{K}\frac{B_{2j}}{(2j)!}(r)_{2j-1}y^{-r-2j+1}.
\]

Here (r)_n is the rising factorial. The exact enclosure is

\[
\boxed{|\zeta(r,v)-Z|\le
\frac{|B_{2K}|}{(2K)!}(r)_{2K-1}y^{-r-2K+1}.}
\tag{2}
\]

Proof: Euler--Maclaurin with the last j=K term included leaves a
periodic-Bernoulli integral with derivative of order 2K. The Fourier
series gives |Btilde_(2K)(x)|<=|B_(2K)|. Integrating
(r)_(2K)(y+x)^(-r-2K) gives (2). Both signs of this remainder
are retained. The last included term is not itself claimed to be
the first omitted term. Exact Bernoulli numbers come from PIW's
finite binomial recurrence. All rational summands are rounded outward.

## 3. Euler coefficients, combinatorics and omitted primes

Absolute convergence for real s>1 gives

\[
E_k=\frac{(-1)^{k+1}}{k!}
\sum_p(\log p)^{k+1}\sum_{m\ge1}m^k p^{-ms}.
\tag{3}
\]

Define integer polynomials P_0(u)=1 and

\[
P_{k+1}(u)=(1+ku)P_k(u)+u(1-u)P_k'(u).
\]

Induction by u d/du proves sum_(m>=1)m^k u^m=
uP_k(u)/(1-u)^(k+1). Their coefficients are nonnegative, and
P_k(1)=k! for k>=1, with P_0(1)=1. Thus all powers of every
retained prime are summed exactly. If P=p^s and P_k(u)=sum_j a_j u^j,
the rational weight used before interval multiplication is

\[
\frac{\sum_j a_jP^{k-j}}{k!(P-1)^{k+1}}.
\]

This avoids losing precision by separately rounding a huge numerator
and a sub-grid reciprocal, the issue already diagnosed by PIW.

For prime cutoff C=4096 put d=1-(C+1)^(-s). Every omitted positive
coefficient sum in (3), after division by k!, is bounded by

\[
T_k=d^{-k-1}\int_C^\infty (\log x)^{k+1}x^{-s}\,dx,
\]
\[
\int_C^\infty(\log x)^h x^{-s}dx=
C^{1-s}\sum_{j=0}^h\frac{h!}{(h-j)!}
\frac{(\log C)^{h-j}}{(s-1)^{j+1}}.
\tag{4}
\]

Indeed P_k(u)<=k! on [0,1] (also valid for k=0), each denominator
is bounded by d^(k+1), and the omitted-prime sum is at most the
omitted-integer sum. The integrand decreases for s log C>h, which
the implementation checks. The integral test and repeated integration
by parts prove (4). The tail enters with sign (-1)^(k+1); no prime
or prime power is discarded. The execution retains 564 primes,
the last 4093, and all 16 signed tails.

## 4. Exact change to dimensionless lambda layers

Set w=s-1/2, a=w^2. Then

\[
z(\lambda)=w(\sqrt{1+\lambda}-1),\qquad
q(a(1+\lambda))=
\frac{Q(s+z(\lambda))}{2w\sqrt{1+\lambda}}.
\tag{5}
\]

The rational binomial coefficients of (1+lambda)^(1/2) and
(1+lambda)^(-1/2), exact truncated composition, and outward
convolution determine every coefficient through degree 15. The
inner series has zero constant, so no unknown higher source layer
can influence those coefficients. This is exact jet algebra; a
degree-15 Taylor polynomial is not claimed to equal q at s=31.
Instead the finite resolvent in ILC-T4 is evaluated there, then
compared with PIW's independent complete source box.

The normalization c_0=1 is exact. Other coefficients are divided
outwardly by the certified positive q(a). Correlations may widen
the intervals but cannot invalidate enclosure. Native LDL decides
the actual forms' signs; midpoint values never certify positivity.
