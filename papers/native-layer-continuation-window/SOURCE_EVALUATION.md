# Actual LCW-1 execution

The original source is unchanged: Xcal(w^2)=xi(1/2+w),
q=Xcal'/Xcal. The anchor is original s=29, w=57/2, a=3249/4.
`sourceLayers(29,17)` independently computes the complete normalized
derivative tower c_0,...,c_17. It reuses the frozen ILC evaluator,
including 564 primes through cutoff 4096, the exact Eulerian-polynomial
prime weights, signed omitted-prime tails, Gamma/pole terms and the
shift-512, 32-term Hurwitz expansion with its explicit remainder.
The last Bernoulli term is included, as documented in ILC; it is not
mislabelled a first omitted term. All 18 moment widths are below 1e-66.

The old 16 moment enclosures overlap the retained ILC certificate.
The next two derivatives are evaluated from this complete source;
they are not obtained by picking canonical positions or a zero measure.
No evaluator or predecessor certificate was edited.

The following decimals are approximations for orientation. The full
outward rational bounds are in CERTIFICATE.json.

| Reading | Approximate value |
|---|---:|
| c_16 | 0.00210082766879004928 |
| Lower admissible c_16 | 0.00210082766670829749 |
| Upper admissible c_16 | 0.00210082767894994720 |
| c_17 | 0.00167651502072235002 |
| Lower admissible c_17 | 0.00167651502047142674 |
| Upper admissible c_17 | 0.00167651502219916656 |
| Canonical position p_16 | 0.17005484065229024 |
| Canonical position p_17 | 0.14523209771876147 |

Certified coarse bounds are 0.170<p_16<0.171 and
0.145<p_17<0.146. The two new Schur pivots exceed 2.509e-13 and
1.476e-12 respectively. Independent 9x9 LDL elimination agrees with
these innovations and certifies both complete enlarged forms.
All 171 binomial readings r+k<=17 are strictly positive.

The even width W_16 lies in (1.224e-11,1.225e-11); the odd width
W_17 lies in (1.727e-12,1.728e-12). The next, unexecuted even
width W_18 lies in (2.144e-13,2.145e-13), and agrees with the
canonical-position width recurrence. Its approximate endpoints are
0.00133946287110567269 and 0.00133946287132015386.
This forecasts the admissible interval for c_18. It does not compute
c_18 or certify that the actual next derivative occupies this interval.

At original s=31, t=3721/4 and lambda=472/3249, the previously
independently evaluated source value is held out. The new response
lies strictly above the old eight-direction response and below that
actual source value. Its error is positive and less than 1.412e-31
(approximately 1.41157535576e-31). Independently, ILC-T4 yields

\[
\frac{(w+w_0)^2}{w}
\left|\frac{w-w_0}{w+w_0-2}\right|^{18}<7.417\times10^{-25},
\quad w=61/2,\quad w_0=57/2.
\]

This analytic bound needs the certified finite source signs and the
retained half-plane majorant, not RH. Both bounds are checked in
native exact arithmetic. The comparison does not identify the fitted
poles with actual xi zeros.

Validation: 67 exact checks, 4 symbolic certificates and 15 negative
controls. Controls include each wrong side of both windows, both even
endpoints, one-positive-innovation failure, a positive finite binomial
triangle with a negative moment determinant, and an admissible chosen
sequence distinct from the source under discussion. Exact uniform
moment fixtures test the arbitrary-order identities at several orders;
the general proof is LCW-T1–T3, not those finite examples.
