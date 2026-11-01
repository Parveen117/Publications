# NAC-1 — Native unit-angle calibration and finite loop curvature

15 October 2026, India. Continuation of native phase observation, QG-R1,
NJC-1 and N4H-1. Original native arithmetic and rewrite engines unchanged.

The proposed relation between exponential transport and a turn reference
has a precise form. A retained generator and angle lift give

\[
 U(\theta)=\exp\!\left(\frac{\theta}{\theta_0}L_0\right),
 \qquad L_0=\theta_0 A,\quad A^2=-I,\quad A^\dagger=-A.
\]

The angle ratio is invariant under a change of angular unit. Native `pi`
specifies the half-turn of this circular sector. The scalar numbers `e`
and `pi` are fixed; the local quantities are the generator, angle and
transport. A bare numerical `e/pi` does not specify a curvature.

For fixed unit circular directions `A,B`, with
`AB+BA=-2 mu I`, define
`H=V^dagger U^dagger V U`, `U=exp(theta A)`, `V=exp(phi B)`.
The exact return defect is

\[
 \boxed{\chi(\theta,\phi)
 =\frac{1-\operatorname{sc}H}{2}
 =(1-\mu^2)\sin^2\theta\sin^2\phi.}
\]

Thus one nonzero-angle calibration determines the coefficient
`K=chi0/(sin^2(theta0) sin^2(phi0))=1-mu^2`, and predicts the return
defect at any other angles **within the same fixed-direction protocol**.
The finite-angle correction is explicit; it is not an extrapolation from
an infinitesimal commutator.

| Calibration data | Exact result |
|---|---:|
| `sin(theta)=4/5`, `sin(phi)=3/5`, perpendicular directions | `chi=144/625` |
| Same angles, direction overlap `mu=3/5` | `chi=2304/15625` |
| Same angles, coincident directions | `chi=0` |
| One-axis defect `kappa=2/5`, twice the retained angle | `kappa2=32/25` |

For fractional rescaling `p/q` in lowest terms, the endpoint rotation
requires winding modulo `q`; arbitrary real rescaling is supplied by the
full angle lift. Angle-dependent directions require their ordered transport
data. These are the records that make the calibration usable.

- [THEOREM.md](THEOREM.md): definitions, proofs, reconstruction domain and lineage.
- [CERTIFICATE.json](CERTIFICATE.json): native algebraic witnesses and exact controls.

Replay from the repository root:

```bash
python3.12 papers/native-angular-calibration/run.py --check
```

Python 3.12 performs I/O only. Scientific arithmetic uses the original native
JS engines. The certificate verifies algebraic identities and exact fixtures;
the completed exponential, angle lift and limit arguments have written proofs.
No new matter coupling, spacetime dynamics or gravitational constant is
postulated by this packet.
