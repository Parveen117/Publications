# LAT-1 actual execution

The reference remains original s=29, w_0=57/2, a_0=3249/4.
The unchanged `sourceLayers(29,19)` computes twenty source layers.
It retains all 564 Euler primes through cutoff 4096, full prime-power
weights, omitted-prime tails, pole terms and Gamma/Hurwitz remainder
enclosures. All twenty moment widths are below 1e-62. All eighteen
overlapping old layers agree with the LCW certificate.

The values below are approximate displays, not replacement bounds.
CERTIFICATE.json contains the outward rational enclosures.

| New source reading | Approximate value |
|---|---:|
| c_18 | 0.00133946287113518028 |
| c_19 | 0.00107115930733865748 |
| LCW even canonical position p_18 | 0.13757658102269333 |
| LCW odd canonical position p_19 | 0.13674815771485775 |
| New A pivot | 3.47997157561e-15 |
| New B pivot | 2.19680610251e-14 |

The actual pair passes both sharp windows. The previous c_18 window
is reproduced. Direct LDL independently certifies both complete 10x10
forms; all 210 binomial readings through total depth 19 are positive.

LAT-T4 is executed at K=19 using the exact-source finite model guaranteed
by those forms. No interval midpoint is substituted for that model.
The following tau values are approximate; positivity is checked using
their full rational values and outward matrix elimination.

| Target order m | Entire rho interval | tau_A | tau_B |
|---|---|---:|---:|
| 3 | [17/20,23/20] | 6.66898999778e-5 | 5.89614398673e-5 |
| 4 | [11/12,13/12] | 5.03965293197e-6 | 4.79555857945e-6 |

Both adjusted forms are strictly positive in both rows. The bounds
apply simultaneously to every real rho in the stated interval.
The source coordinate is s=1/2+(57/2)sqrt(rho). The first row contains
rho_27=2809/3249 and rho_31=3721/3249, and the full real interval
between them. The second row gives the stronger four-direction
certificate on a smaller interval. Ten-direction positivity is
certified at the reference point only.

The wider attempted m=3 range rho in [5/6,7/6] is not certified by
this sufficient estimate. This is explicitly an inconclusive estimate,
not evidence that an actual form becomes negative there.

The finite Taylor forecast at s=31 encloses both the independently
evaluated source value and its first derivative from PIW. The separate
ten-direction response has positive residual below 4.160e-35
(approximately 4.15937701688e-35), improving the prior response.
ILC's independent Schwarz bound with exponent 20 is below 9.131e-28.
Both comparisons use exact outward arithmetic. The source and finite
rational model remain distinct.

Validation: 61 exact checks, 3 symbolic certificates and 14 negative
controls. Exact finite measures test anchor composition/inversion,
required weight changes, invariant poles/residues/response, polynomial
transport, moment and binomial evolution. Negative controls reject
missing source layers, a nonpositive surrogate, divergent tail ratios,
uncarried weights and unsupported continuation beyond the certified
interval. The written uniform proof supplies the continuum result;
finite examples alone do not.
