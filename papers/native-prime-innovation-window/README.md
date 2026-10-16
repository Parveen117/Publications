# PIW-1 — Sharp innovation windows and a new complete prime reading

17 October 2026 (Asia/Kolkata). Continues POC-1 in Publications.

The next-source acceptance rule now has a sharp scalar form. If
delta=q(x)-q_N(x), the two Schur innovations satisfy

\[
x\alpha+\beta=\delta(1-\kappa\delta),\qquad \kappa>0.
\]

An admissible derivative exists exactly when 0<delta<1/kappa.
The actual derivative must also lie in the derived strict slope window.
Two dimensionless coordinates record the value position and the slope
position within these windows. Both are invariant under an invertible
old-source basis change when the observation columns move with it.

A determinant normalization removes the double zeros already imposed
by the retained observations and gives an exact Riccati-type inequality
for the remaining actual response. This exposes the source condition
that still needs all-order control; the factorization does not assert
that the actual source automatically satisfies it.

The finite execution advances from eight to nine prime directions.
A new complete xi jet at original s=29 includes Gamma, pole and Euler
terms with explicit outward tails. Both new innovations and both 9x9
forms are positive. The value and slope positions lie respectively in
(0.140,0.141) and (0.142,0.143). At independently evaluated, unfitted
s=31, the response error is positive and less than 3e-17.

- [Four theorems and proofs](THEOREM.md)
- [Complete source evaluation and remainder bounds](SOURCE_EVALUATION.md)
- [Source citations and scope](SOURCE_MAP.md)
- [Usage and replay contract](USAGE.md)
- [Source pins](SOURCE_PINS.json)
- [Exact certificate](CERTIFICATE.json)

```bash
python3.12 papers/native-prime-innovation-window/run.py --check
```

The scientific computation uses the existing native rational matrix and
rewrite engines, with POC's outward interval type. Python 3.12 handles
I/O and replay orchestration only. The old source packets remain frozen.
This certifies a new ninth prime direction; it does not certify the full
old 24-direction bank plus that direction. All-order positivity, global
eta zero, RH and physical identification remain separate open targets.
