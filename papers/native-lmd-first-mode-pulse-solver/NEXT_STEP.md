# The remaining case needs spatial relations, not only more amplitudes

Author: Monty. 31 October 2026. Post-execution mathematical diagnosis.

The remaining case is h=3/8, t=1/2, k=2 with requested probability
halfwidth 1/200. Its first-mode solver is valid but returns halfwidth about
0.01301658. Even ideal unbinned first-mode conditioning leaves a term between
0.012662092 and 0.012662095 in the chosen conditional-variance certificate.
Finer q/radial cells cannot remove that term.

## 1. Use existing mode moments before constructing another large partition

The earlier characteristic packet already certified I_N(h) for selected
Fourier cutoffs. Applying the same probability transfer factor

    exp(-h k^2/4) a_t(k)^2 / 4

to I(h)-I_N(h) gives the following *ideal mode-truncation budget terms*.
They contain no finite source quantization error and are not lower bounds
on the actual oscillatory approximation error. Display endpoints are outward.

| Retained nonzero Fourier pairs N | Ideal probability-budget interval | Below 0.005 before source bins? |
| --- | --- | --- |
| 0 | [0.021952198, 0.021952202] | No |
| 1 | [0.012662092, 0.012662095] | No |
| 2 | [0.009395529, 0.009395532] | No |
| 4 | [0.006575711, 0.006575714] | No |
| 8 | [0.004438428, 0.004438431] | Yes |
| 16 | [0.002937136, 0.002937140] | Yes |

The native [diagnostic](frontier.cjs) also subtracts the inherited phase
interpolation budget. N=8 is the first **previously tested cutoff** leaving
positive room for finite source compression. It is not proved to be the
smallest integer cutoff, and a working 8-mode expectation solver has not yet
been supplied. N counts spatial Fourier pairs, not qubits or physical
spacetime dimensions.

## 2. Why amplitude-only conditioning has a persistent budget

The first-pair reduction may suggest discarding every Fourier orientation
and keeping only all radii. For this source that is a legitimate conditional
expectation, but the following bound shows it is insufficient for the
requested quadratic certificate in the remaining case.

Write the nth unnormalized sine/cosine pair as independent Gaussians with
variance h/n. Let its radius be R_n. The radii are mutually independent and
independent of q~N(0,h). Conditional angular averaging, using the positive
Wick martingale and the Bessel identity from THEOREM.md, gives

    E[J | q,R_1,...,R_N]
      = exp(q-h/2) product_(n=1)^N exp(-h/(2n)) B_0(R_n).

Each factor has mean one and second moment B_0(h/n). Therefore

    E[E[J | q,R_1,...,R_N]^2]
       = exp(h) product_(n=1)^N B_0(h/n).                (1)

For x>=0, the standard positive Bessel series gives

    B_0(x)=sum_(k>=0) (x^2/4)^k/(k!)^2
           <= sum_(k>=0) (x^2/4)^k/k! = exp(x^2/4).

The increasing conditional-expectation martingale is L2 bounded by E J^2.
Its limit is E[J|q,R_1,R_2,...]. Using sum_(n>=1) 1/n^2=pi^2/6 in (1),

    E[E[J | all amplitudes]^2] <= exp(h+h^2 pi^2/24),
    ||J-E[J | all amplitudes]||_2^2
        >= I(h)-exp(h+h^2 pi^2/24).                      (2)

All finite amplitude-bin schemes capture no more second moment than this
ideal conditioning. For the remaining case, the native outward interval for
the probability transfer factor times the right-hand side of (2) is

    [0.006357675, 0.006357679] > 1/200.                  (3)

Thus adding arbitrarily many radii while discarding their angular relations
cannot make this particular variance-remainder certificate meet the target.
This is not a no-go for amplitude methods using a sharper error theorem,
and (3) does not bound the true probability approximation error from below.
It identifies a limitation of the current certificate and source compression.

## 3. An explicit relative-phase term

With two spatial harmonics, set

    s_2(theta)=q+r_1 cos(theta-alpha_1)
                     +r_2 cos(2 theta-alpha_2),
    delta=alpha_2-2 alpha_1.

A common spatial-origin shift removes alpha_1 but preserves delta. The
spatial integral in the mass is

    Z(r_1,r_2,delta)
      = B_0(r_1) B_0(r_2)
        +2 sum_(m>=1) B_(2m)(r_1) B_m(r_2) cos(m delta),  (4)

where B_m is the standard modified Bessel function I_m. To prove (4), expand
each exponential in its absolutely convergent Laurent series and retain
the constant coefficient. The required generating identity is
[DLMF 10.35.1](https://dlmf.nist.gov/10.35.E1); its use here is not a new
special-function theorem. The Wick prefactor is exp(q-5h/4), because the
diagonal source variance through N=2 is h(1+1+1/2).

For r_1,r_2>0, positivity of the Bessel coefficients gives

    Z(r_1,r_2,0)-Z(r_1,r_2,pi)
        =4 sum_(m positive odd) B_(2m)(r_1) B_m(r_2) > 0. (5)

The two fields have identical radii and different total mass. Already the
leading relation-dependent contribution in (4) is
r_1^2 r_2 cos(delta)/8. Independently deleting both orientations would erase
this contribution; a common observer-coordinate rotation does not erase it.

The next concrete numerical target is a certified source compression for
N=8 that retains enough of these relative-phase relations, with all phase,
tail and partition errors within the remaining budget. A coarser conditional
compression may be used only with its actual residual variance. No supplied
finite encoding can be silently identified with the whole field state.

This diagnosis was performed after the four-case refinement, uses frozen
mode moments and new written inequalities, and did not tune or rerun those
cases. It supplies a better next algorithmic target. It does not close
autonomous gravitational dynamics, local apparatus or the physical 3+1 gate.
