# PIW-1 complete source evaluation and outward-error proof

The arithmetic below realizes the actual completed xi response at odd
integer original s. It uses the unchanged native F engine and POC's
110-place directed rational intervals. It does not fit a model to zeta
zeros, assume RH, or remove a source term because its contribution is
small. The old source boxes are imported unchanged; the new jets are
recomputed by the procedure here.

## Elementary constants

For 0<x<1, integrating the finite geometric identity for 1/(1+t^2)
gives the alternating arctangent series, with absolute error at most
x^(2m+1)/(2m+1) after m terms. We take m=96 and use
pi=16 atan(1/5)-4 atan(1/239). The tangent addition law gives exactly
tan(4 atan(1/5)-atan(1/239))=1; elementary angle bounds place that
angle in (0,pi/2), fixing pi/4 rather than another tangent branch.

For any positive rational y, first scale y=2^k z with 1<=z<2. If
v=(z-1)/(z+1), then 0<=v<1/3 and

\[
\log z=2\sum_{j=0}^{m-1}\frac{v^{2j+1}}{2j+1}+E_m,
\qquad
0\le E_m\le\frac{2v^{2m+1}}{(2m+1)(1-v^2)}.
\]

Use m=128. The same formula at z=2 supplies log 2. The remainder
follows by bounding each omitted denominator below by 2m+1 and
summing a positive geometric series. Log pi is enclosed by applying
this monotone rational logarithm to the lower and upper pi endpoints.
No approximate decimal constant is treated as exact.

## Gamma and pole contributions

The completed definition in [DLMF 25.4](https://dlmf.nist.gov/25.4)
gives Q and Q' in PIW equation (12). For v=s/2 set y=v+256. The
digamma expansion and its positive-real signed remainder bound are
[DLMF 5.11.2 and 5.11(ii)](https://dlmf.nist.gov/5.11):

\[
\psi(y)=\log y-\frac1{2y}
 -\sum_{k=1}^{K}\frac{B_{2k}}{2k y^{2k}}+R_K,
\qquad
|R_K|\le\frac{|B_{2K+2}|}{(2K+2)y^{2K+2}},\quad K=24.
\]

The remainder has the sign of -B_(2K+2). For completeness this bound
also follows from the positive-real Binet representation of psi:
expand 1/(1+(t/y)^2) by its exact finite geometric identity inside
the positive Bose integral. The residual denominator is at least one,
so its magnitude is at most the first omitted moment and its sign is
fixed. This justifies a finite enclosure, not merely an asymptotic
approximation.

All Bernoulli numbers are generated as exact fractions from B_0=1 and
sum_(k=0)^n binom(n+1,k)B_k=0 for n>=1. The
[digamma recurrence](https://dlmf.nist.gov/5.5#E2) gives
psi(v)=psi(y)-sum_(k=0)^255 1/(v+k). For s=2n+1 use the exact
[half-integer trigamma identity](https://dlmf.nist.gov/5.15#E4):

\[
\psi'(s/2)=\pi^2/2-4\sum_{k=1}^n(2k-1)^{-2}.
\]

Thus Q retains 1/s+1/(s-1), half the digamma reading, -log pi/2,
and the Euler logarithmic derivative. Q' retains the derivatives of
all four contributions. Coordinate conversion gives the two formulas
in PIW (12), including the w^-3 term in q'.

## Complete Euler contribution and all omitted primes

Absolute convergence of the Euler product for real s>1 allows
termwise logarithmic differentiation. Summing all powers of each
retained prime exactly gives

\[
\frac{\zeta'}\zeta(s)
 =-\sum_p\frac{\log p}{p^s-1},\qquad
\left(\frac{\zeta'}\zeta\right)'(s)
 =\sum_p\frac{(\log p)^2p^s}{(p^s-1)^2}.
\]

The source is the [Euler product](https://dlmf.nist.gov/25.2#E11),
not a finite replacement for it. A finite integer sieve retains every
prime p<=M=4096, namely 564 primes with last prime 4093. For omitted
p>M, 1-p^-s >= c=1-(M+1)^-s. Consequently the two positive tails
are bounded by c^-1 sum_(n>M) log n/n^s and
c^-2 sum_(n>M) (log n)^2/n^s. For s>=3 and M>=2 both summands
decrease on [M,infinity), so the sums are at most their integrals from
M to infinity. Integration by parts gives exactly PIW (13).

Every retained logarithm and both integral bounds are computed outwards.
The negative Euler term is enclosed by [-retained-T_0,-retained],
and its positive derivative by [retained,retained+T_1]. Nothing beyond
the cutoff is assigned zero.

Form p^s/(p^s-1)^2 as an exact native rational before multiplying its
logarithm interval. On a fixed absolute grid, splitting it into a huge
factor and a reciprocal below the grid is still outward-valid but can
produce uselessly wide bounds. A regression control records that loss
of precision. It is not evidence of a negative actual source sign.

## Certification boundary

An interval operation encloses every exact operation on its input
intervals. Integer floor and ceiling implement directed rounding,
including for negative values. Native outward LDL therefore certifies
strict positivity whenever all pivot lower bounds are positive.
The response solve and both Schur evaluations retain all mixed entries.
Their dependency overestimates may widen a box; they cannot justify a
sign when its lower bound does not establish it.

The written formulas and tail arguments justify that the generated
boxes contain the actual source. Agreement with the earlier s=23
archive is an independent numerical cross-check, not a substitute for
that derivation. The native certificate proves the stated finite signs
from these enclosures. It does not turn the infinite RH premise into a
consequence of a finite calculation.
