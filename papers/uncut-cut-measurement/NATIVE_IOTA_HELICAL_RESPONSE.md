# v1.8: native iota, helical sheet transport and a derived response

Monty Dabas — 25 September 2026.

**Result.** Starting with the existing cut-complex scalar and the existing
EMK-G3 helical return map, a declared local quadratic response cost produces
one positive operator for both a stationary source response and a native
norm-preserving continuation. On an infinite integer-sheet orbit, an exact
source kernel is `G_s=(2/3)(1/2)^|s|`. Its sheet-local reading is `2/3`; the
three-site sheet-summed reading is `6/7`. Their difference `4/21` is an exact
sum of hidden-sheet contributions. The quotient is correct for its summed
target; it must not be reported as a sheet-local measurement.

This is a constructive mathematical extension under stated constitutive
inputs. It is not a derivation of the physical gravitational field, a
physical clock, detector frequencies or the numerical action scale. The
bounded infinite-operator statements below have written proofs; the exact
computational certificate checks specified algebraic obligations and controls.

## 1. Consume the foundation already present

The source-pinned [F00-E theorem](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/86198d29cbf30390059079f38675c952e506ea9c/theorems/foundation/F00E_NATIVE_EULER_FROM_IOTA_COMPLEX.md)
starts after the cut-complex construction and derives the native exponential
by factorial polynomials and completion. Its existing order is

`oriented cut -> iota_Sigma -> completed C_Sigma -> Exp_Sigma -> Euler flow`.

In particular `iota_Sigma^2=-1`, `iota_Sigma^dagger=-iota_Sigma`, and
`N(a+iota_Sigma b)=a^2+b^2`. Exponential addition, inverse, dagger covariance,
circular and hyperbolic factors are available dependencies. They are not new
claims of this note. F00-E assumes its previously established complete ordered
radial field; the finite F00/F00-E certificate audits the algebra and bounds,
not all universal foundation statements by sampling.

The exact witness uses an **unchanged copy** of the upstream
`NativeCutScalar` implementation. Its source hash, F00-E's source hash and
the consumed EMK-G2/G3 code hashes are recorded in
[NATIVE_HELICAL_SOURCE_PINS.json](NATIVE_HELICAL_SOURCE_PINS.json).
No Python complex number or floating-point proof margin is used.

There is also an explicit scalar-to-pair intertwiner:

`M(a+iota_Sigma b) = [[a,-b],[b,a]]`.

The multiplication rule gives `M(zw)=M(z)M(w)`,
`M(z^dagger)=M(z)^T`, and injectivity follows from the first column.
Consequently `M(iota_Sigma)=[[0,-1],[1,0]]` **in this specified module**.
This closes the unbound represented-quarter-turn step for this construction;
it does not identify every earlier two-component carrier with the native
scalar without an intertwiner.

## 2. Keep EMK geometry and its helical lift typed

EMK-G1 already supplies `g=A(v)^2 du^2+dv^2`, `A>0` even,
and `K=-A''/A`. EMK-G2 supplies the periodic base seam. EMK-G3 supplies
the mapping-torus gluing

`(u,v,f) ~ (u+L,v,rho^(-1)f)`,

where the exact certificate carries phase modulo `KMOD=6` and

`rho(phi,sigma,k)=(phi+alpha mod 6, sigma+beta, k+q)`.

The field amplitude in `C_Sigma` and the fibre's phase coordinate `phi` are
different objects. No physical or scalar-phase identification is assumed.

Choose a base cycle with `N>=3` nodes and follow **one full lifted orbit**.
For `s=Nm+j`, `0<=j<N`, label its state by

`ell_s=(j,rho^m f_0)`.

If `q!=0`, this label map is injective: equality of the integer sheet forces
the same `m`, then the base label forces the same `j`. Increasing `s` through
the seam applies `rho`, and decreasing it applies `rho^-1`. The resulting
orbit is an infinite chain, not a finite ring with the sheet silently wrapped.
Other lifted orbits are outside this single-orbit sector.

The witness uses the **existing G3** values `alpha=3`, `beta=2/5`, `q=2`
and `N=3`. At `s=6`, two base circuits give fibre `(0,4/5,4)` from seed
`(0,0,0)`. Phase and base return; scale and integer sheet do not.

### How the metric can enter the cost

For a supplied scalar Dirichlet cost on the EMK surface,

`E_grad = (1/2) integral [A^-1 N(partial_u psi) + A N(partial_v psi)] du dv`.

This follows from `g^-1=diag(A^-2,1)` and `dV_g=A du dv`. Its local scalar
operator is `Delta_g=A^-2 partial_u^2+partial_v^2+(A'/A)partial_v`.
The global field must in addition respect the mapping-torus gluing.

On the seam the induced metric is `A(0)^2 du^2`. A one-dimensional
Dirichlet cost discretized by piecewise-linear fields at spacing `h=L/N`
has uniform edge weight proportional to `1/[A(0)h]`. Thus a declared
normalization can give the `w=1` seam model below. Selecting this cost,
normalization and single-seam restriction is extra constitutive information.
This note does not solve the full transverse curved surface. For the
quadratic EMK-G1 family, `A(0)=1` for every kappa, so this seam-only
experiment cannot identify transverse metric curvature.

## 3. NH-1: the source operator follows from the local cost

Complete the finitely supported `C_Sigma`-valued fields on the orbit in the
norm square `sum_s N(psi_s)`. Denote this sequence space by `H_orb` and
let `(S psi)_s=psi_(s-1)`. The shift is unitary: reindexing the sum gives
`S^dagger=S^-1`. No continuous Haar measure is required for this countable
completion.

Fix radial constants `w>0`, `eta>0` and a square-summable source `b`. Declare

`E(psi;b) = (w/2) sum_s N(psi_s-psi_(s-1))`
`           + (eta/2) sum_s N(psi_s) - Re sum_s b_s^dagger psi_s`.

Here `eta` is a **grounding stiffness**, not a derived particle mass. It
grounds every site to the fixed background and changes U36's ungrounded
bulk-shift symmetry. A gapped, screened one-dimensional response is the
intended mathematical sector, not an inverse-square gravity model.

With `D=I-S`, variation of this explicitly declared quadratic law gives

`H=w D^dagger D+eta I=(eta+2w)I-w(S+S^dagger)`.

**Proof and consequences.** Expansion of each norm square gives this
operator, and

`<psi,H psi> = w ||D psi||^2+eta ||psi||^2 >= eta ||psi||^2`.

Thus `H` is self-adjoint, bounded by `eta+4w`, and strictly positive.
Writing `H=(eta+2w)[I-w(S+S^dagger)/(eta+2w)]` gives a norm-convergent
Neumann inverse because `2w/(eta+2w)<1`. It follows that

`psi_*=H^-1 b`, `||H^-1||<=1/eta`,

and completing the square gives

`E(psi;b)-E(psi_*;b)=(1/2)<psi-psi_*,H(psi-psi_*)>`.

Therefore the stationary response is unique. With one operator fixed,
source superposition and source reversal follow. Locality, quadratic degree,
weights, grounding and source assignment were declared; iota alone does not
select them. This is the helical-orbit extension of U36's cost argument.

## 4. NH-2: exact integer-sheet Green response

Set `w=1`, `eta=1/2`, and place a unit native-radial source at `s=0`.
Then

`(H psi)_s=(5/2)psi_s-psi_(s-1)-psi_(s+1)`.

The exact solution on the **infinite** lifted orbit is

`G_s=(2/3)(1/2)^|s|`.

**Proof.** For `r=1/2`, `(5/2)r-1-r^2=0`, so the recurrence vanishes
at every nonzero index, on both sides of the source. At zero it is
`(2/3)[5/2-2(1/2)]=1`. The geometric series gives
`sum_s G_s=2` and `sum_s G_s^2=20/27`, so `G` belongs to both the
absolute-sum and square-sum spaces. NH-1's inverse uniqueness now proves
the formula for every integer `s`. No finite truncation is used in this proof.

For general positive `w,eta`, the same calculation gives `G_s=c r^|s|`,
where `0<r<1` solves `(eta+2w)r=w(1+r^2)` and
`c=[eta+2w-2wr]^-1`. This describes the chosen grounded law's decay.
It is not a universal gravitational distance dependence.

## 5. NH-3: native coherent continuation of the same response operator

Use the native quarter-turn and the same `H` to choose the centered generator
`Q=-iota_Sigma H`. This is a continuation policy inherited in form from the
quadratic source programme; multiplying by iota is lawful but not the only
possible policy. Define

`U_tau = sum_(n>=0) (-iota_Sigma tau H)^n/n!`.

**Proof.** The factorial majorant with `||H||<=eta+4w` proves convergence
in operator norm for every finite radial `tau`. The same Cauchy-product
argument as F00-E gives `U_tau U_sigma=U_(tau+sigma)`. Since
`Q^dagger=-Q`, dagger covariance gives `U_tau^dagger=U_(-tau)` and hence
`U_tau^dagger U_tau=I`. Differentiation of the norm-convergent series gives
`dU_tau/dtau=-iota_Sigma H U_tau`. Every series term commutes with `H`.

Thus the affine evolution

`psi(tau)=psi_*+U_tau[psi(0)-psi_*]`

preserves both centered norm and excess source cost. This connects the
native scalar, full lifted orbit, stationary Green response and coherent
continuation in one explicit system. It also proves that the conservative
continuation does not itself relax to equilibrium.

Only after a chosen energy scale and clock adapter does this become an
energy-time equation with laboratory units. No numerical hbar is claimed.
The finite certificate checks coefficients of `U^dagger U` through degree
12 in the exact Laurent shift algebra; the all-orders claim follows from
the written factorial convergence proof, not from that finite check.

## 6. NH-4: exactly what the visible cycle reads

For an absolutely summable field define

`(P psi)_j=sum_(m in Z) psi_(j+Nm)`, `j=0,...,N-1`.

This sums all helical sheet images. It is bounded on the absolute-sum
space; it is **not** a bounded observation on the entire square-sum space.
Reindexing absolutely convergent sums proves `P H=H_cycle P`, where

`H_cycle=(eta+2w)I-w(S_cycle+S_cycle^dagger)`.

So the projected stationary equation is lawful for the **summed** target.
It does not preserve the sheet-local target: `P delta_0=P delta_N`,
although those fields have different value at `s=0`.

For the same source `b=delta_0` and `N=3`, the visible source-site response is

`(P G)_0 = (2/3)[1+2 sum_(m>=1)(1/2)^(3m)] = 6/7`.

Equivalently it is the source diagonal of the inverse of

`H_cycle=[[5/2,-1,-1],[-1,5/2,-1],[-1,-1,5/2]]`.

The sheet-local response is `G_0=2/3`, so

`(P G)_0-G_0=4/21`.

The positive difference is the exact contribution of the other sheets.
These are two different observables of one solution, not two contradictory
predictions for the same detector. Lifting the ring solution periodically
back to the line would also repeat the source; it is not the original
localized-source boundary problem. A physical experiment must specify
which observable its detector reads.

## 7. Evidence and the next physics obligation

The [exact witness](certificates/native_helical_response.py) imports the
unchanged native scalar and G3 transition, checks source hashes, the faithful
scalar-to-pair map, integer-sheet labels, operator factorization and source
recurrence, geometric sums, cycle inverses, quotient intertwining, and unitary
series coefficients. Negative controls change an edge sign, decay factor,
sheet readout, generator quarter-turn and dagger norm. The master source pin
binds this manuscript, the witness, copied scalar and upstream geometry code.

The construction advances a particular mathematical sector. It retains an
already certified native structure instead of replacing it by a bare visible
circle. It uses familiar positive graph operators and factorial-series
methods; no global novelty claim is made for those methods.

The next physical obligation is to identify a real source, sheet-sensitive
observation and energy/clock mapping for this operator, or derive a different
constitutive operator from a stronger native law. The law must be specified
before comparing a held-out measurement. Transverse EMK curvature, a universal
source-mass coupling, Lorentzian spacetime dynamics and nonlinear gravity are
not established by this one-seam calculation.

Reproduce from this directory:

```bash
python certificates/native_helical_response.py
python certificate.py --check
```
