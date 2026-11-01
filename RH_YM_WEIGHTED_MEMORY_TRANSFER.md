# RYM-1 — What the recent RH estimates transfer to Yang–Mills

Author: Monty. 18 October 2026. Written derivation only.

The recent RH comparison hierarchy has a concrete application to the
actual Yang–Mills electric cut: retain a known higher electric energy
in the inverse-memory bound instead of charging every source component
at the lowest hidden energy. On the same side-four lattice, coupling
theta=1/48 and complete 241-dimensional retained space used by YM75,
this gives the improved working enclosure
\[
\boxed{0.49072844\le\Delta_{4^3}(1/48)\le3.00712209.}
\tag{1}
\]
YM75's preceding enclosure was [0.49069036,3.00712209].
This is a modest improvement of a finite-volume bare-energy bound.
It is not a physical particle mass or a continuum mass-gap result.
No simulation, scientific computation or certificate replay is used
to obtain (1); the rational sign proof is written below.

## 1. Source audit and the meaning of prime weights

The RH and latest YM work are on different branches of the same
Publications repository. The RH branch contains the older YM55 files;
its index is not the current YM checkpoint. This audit reads the
following immutable sources, without modifying the YM branch:

- [BRT-T52–T55](https://github.com/Parveen117/Publications/blob/50464900436e768e237808020ab3e63b86e566a2/RH_BOUNDARY_REGULAR_TRIAL.md):
  coherent source comparison and the complete weighted inverse residual.
- [PTE](https://github.com/Parveen117/Publications/blob/50464900436e768e237808020ab3e63b86e566a2/RH_PRIME_TRANSPORT_ENERGY.md):
  prime translations, root normalization and the actual signed source.
- [HCM prime-clock map](https://github.com/Parveen117/Publications/blob/50464900436e768e237808020ab3e63b86e566a2/papers/native-heat-observer-memory/PRIME_CLOCK_MAP.md):
  prime valuations as an optional exact audit of rational gains.
- [YM73](https://github.com/Parveen117/Publications/blob/5290ac92f362941023ff550eb3a8aa09da770302/papers/yang-mills-certified-benchmark/YM73_PHYSICAL_MASS_MEMORY_METRIC.md),
  [YM74](https://github.com/Parveen117/Publications/blob/5290ac92f362941023ff550eb3a8aa09da770302/papers/yang-mills-certified-benchmark/YM74_ALL_COUPLING_ELECTRIC_ENCLOSURES.md),
  [YM75](https://github.com/Parveen117/Publications/blob/5290ac92f362941023ff550eb3a8aa09da770302/papers/yang-mills-certified-benchmark/YM75_RH_SOURCE_RESOLVED_SPECTRAL_TRANSFER.md):
  physical normalization, the complete electric cut, and the previous
  quantitative RH-to-YM transfer.
- [YM91](https://github.com/Parveen117/Publications/blob/5290ac92f362941023ff550eb3a8aa09da770302/papers/yang-mills-certified-benchmark/YM91_CONDITIONAL_MEMORY_FLUX.md),
  [YM92](https://github.com/Parveen117/Publications/blob/5290ac92f362941023ff550eb3a8aa09da770302/papers/yang-mills-certified-benchmark/YM92_WILSON_FLUX_CORRECTOR.md),
  [YM93](https://github.com/Parveen117/Publications/blob/5290ac92f362941023ff550eb3a8aa09da770302/papers/yang-mills-certified-benchmark/YM93_BOCHNER_CONTINUATION.md):
  actual conditional memory, flux correction and the latest
  volume-uniform small-theta window.

These are inspected working proofs and source records, not results
independently validated merely by their presence on GitHub.

In the RH source a prime power n=p^k carries
\[
u_n=\log n,\qquad a_n=\frac{\log p}{p^{k/2}},
\qquad -a_n C_g(u_n).
\tag{2}
\]
They are a translation distance, a source weight and its signed
contribution. PTE constructs a positive difference energy but retains
the scalar subtraction:
\[
q_L[g]=\mathcal E_L[g]-c_L\|g\|^2.
\]
It explicitly does not identify its operator with a physical Hamiltonian
or its eigenvalues with particle masses. Nor does the HCM valuation
record supply such an identification.

The actual YM physical gap uses a different source:
\[
A_\theta=L_E-\theta\sum_pW_p,\qquad H=A_\theta-E_0,\qquad
\Delta_{\rm phys}=\kappa\,\inf\sigma
 (H|_{\Omega^\perp_{\rm phys}}).
\tag{3}
\]
Here p in the sum labels plaquettes, not prime numbers. The physical
clock kappa and vacuum subtraction are essential. For a finite
vacuum-orthogonal observation cut, YM73 also retains
\[
M_0=BD^{-1}B^*,\quad Z=I+BD^{-2}B^*,\quad
K=\kappa Z^{-1/2}(C-M_0)Z^{-1/2}.
\tag{4}
\]
Thus a static returned scalar cannot simply be renamed a mass.
No source intertwiner, physical clock identification or vacuum-sector
identification between (2) and (3) has been supplied here.

The transferable tool is the inverse comparison, with its source and
domain checked afresh. No RH numerical constant is transplanted.

## 2. RYM-T1 — A weighted full-residual bound for the actual electric cut

Use YM74's finite periodic lattice and its complete physical electric
projection P=P_R. Let Q=I-P, let N_plaq be the number of plaquettes,
and put b=|theta|N_plaq. Since |W_p|<=1,
\[
A_\theta=L_E-\theta V,\quad \|V\|\le N_{\rm plaq},\qquad
A_\theta=\begin{pmatrix}C&B\\B^*&D_A\end{pmatrix}.
\]
Let lambda_* be the lowest physical electric energy in Q and d=lambda_*-b.
For a real z<d define
\[
S=D_A-z,\qquad
\mathcal D_z=Q(L_E-b-z)Q,\qquad F=B^* .
\]
Both positive operators have lower bound d-z, and
\[
S\succeq\mathcal D_z,\qquad
S^{-1}\preceq\mathcal D_z^{-1}.
\tag{5}
\]
This is form order. The electric projection reduces L_E and preserves
its domain; the potential is bounded. Inverse order follows from the
variational supremum for the inverse quadratic form. It does not
require D_A to commute with L_E.

For any finite-column trial X with range in the operator domain of S,
write
\[
R=F-SX,\qquad
T_X=F^*X+X^*F-X^*SX.
\]
YM75's exact completion now gives the sharper enclosure
\[
\boxed{
T_X\preceq F^*S^{-1}F
 =T_X+R^*S^{-1}R
 \preceq T_X+
 \sum_{\lambda>R_{\rm cut}}
       \frac{R^*\Pi_\lambda R}{\lambda-b-z}.
}
\tag{6}
\]
Here R_cut denotes the cutoff parameter in P_R, to distinguish it
from the residual R; Pi_lambda are the full physical electric
eigenspace projections in Q. Multiplicity is retained.

The sum converges in finite-matrix norm: the columns of R lie in the
Hilbert space and every inverse weight is at most (d-z)^(-1).
The exact hidden inverse remains infinite dimensional. Compared with
YM75's uniform error, the lower Schur endpoint improves by
\[
\boxed{
\sum_{\lambda>R_{\rm cut}}
\left(\frac1{d-z}-\frac1{\lambda-b-z}\right)
R^*\Pi_\lambda R\succeq0.
}
\tag{7}
\]
Strict improvement in a direction requires a nonzero residual reading
above lambda_*. It is not inferred just from a nonzero residual.

Equations (5)–(7) apply the standard inverse-order and square-completion
principles also used in BRT-T53–T54. This is not a new general positivity
theorem. All mixed source terms remain within each Gram. For finite
electric polynomial trials the residual has finite electric support,
by YM74's multiplication-shell bound; this does not truncate S^(-1).

## 3. RYM-T2 — An actual energy-eight source component

Specialize to YM74–YM75's periodic side-four lattice and R_cut=3.
There are 192 elementary plaquettes and 48 winding four-cycles.
The complete retained orthonormal basis is
\[
1,\qquad \phi_c=2W_c\quad(c\text{ one of those }240\text{ cycles}).
\]
The hidden electric floor is lambda_*=9/2. For each elementary
plaquette p define its spin-one character
\[
a_p=\chi_1(U_p)=4W_p^2-1.
\tag{8}
\]
These are smooth physical vectors, orthonormal in the product Haar
Hilbert space. Each has spin one on four links, so
\[
L_Ea_p=8a_p.
\]
The plaquette holonomy has Haar distribution. Using YM74's moments
Phi(W_p^2)=1/4 and Phi(W_p^4)=1/8 gives directly
Phi((4W_p^2-1)^2)=16/8-8/4+1=1. Distinct plaquettes
have an edge belonging to only one of them; integration over that
edge makes their inner product zero. Let Pi_adj project onto their
192-dimensional span. It is a subprojection of the full energy-eight
eigenspace and lies in Q.

Since the electric operator does not couple P to Q,
F=-theta QVP. Character multiplication gives
\[
W_p\phi_p=2W_p^2=\frac{1+a_p}{2}.
\]
Consequently the exact retained source reading is
\[
\boxed{
\Pi_{\rm adj}F1=0,\quad
\Pi_{\rm adj}F\phi_p=-\frac{\theta}{2}a_p,\quad
\Pi_{\rm adj}F\phi_w=0
\quad(w\text{ a winding cycle}),
}
\]
\[
\boxed{F^*\Pi_{\rm adj}F=\frac{\theta^2}{4}P_{\rm plaq}.}
\tag{9}
\]
P_plaq is the retained projection onto the elementary plaquette traces.

To verify that (9) has no omitted cross term, consider
< a_p,W_q phi_c >. If q and c are distinct four-cycles, an edge in
their symmetric difference has odd fundamental parity. Changing that
link by the central element -I reverses W_q phi_c and leaves a_p
unchanged, so the integral vanishes. If q=c is an elementary plaquette
different from p, an edge of p absent from q instead gives zero by
Haar integration of the nontrivial spin-one matrix coefficient.
Only q=c=p survives. The constant-column assertion follows by the
same parity rule. This proves (9) on the complete retained basis.

Let s=9/2-b-z>0 and t=8-b-z=s+7/2. Retaining just this known
part of the electric hierarchy already gives
\[
S\succeq sI+\frac72\Pi_{\rm adj},\qquad
S^{-1}\preceq s^{-1}I-(s^{-1}-t^{-1})\Pi_{\rm adj}.
\]
At X=0, the whole exact memory therefore satisfies
\[
\boxed{
B(D_A-z)^{-1}B^*
\preceq\frac{BB^*}{s}
-\frac{\theta^2}{4}\left(\frac1s-\frac1t\right)P_{\rm plaq}.
}
\tag{10}
\]
Every other hidden component is still charged at the full uniform
allowance. The subtracted correction is justified by a source
projection, not by deleting a hidden mode. The comparison is strictly
sharper on the plaquette subspace for nonzero theta. It is not
positive definite as an improvement on all 241 retained directions.

## 4. RYM-T3 — A written rational sign test improves the same gap

Set theta=1/48, so b=4 and d=1/2. YM75 supplies
\[
BB^*=\theta^2\operatorname{diag}
 \left(0,\frac{191I+J}{4},48I\right).
\tag{11}
\]
Split the plaquettes into the uniform vector and its 191 transverse
directions. The retained constant/uniform coupling has square
k^2=48theta^2=1/48. By (10), the lower Schur comparison in that
two-dimensional sector is
\[
\begin{pmatrix}
-z&-k\\
-k&3-z-\dfrac{382}{9216s}-\dfrac{1}{9216t}
\end{pmatrix},
\quad s=\frac12-z,\quad t=4-z.
\tag{12}
\]
The old numerator 383 has become 382 at the low denominator, with
the remaining one charged at the correct higher denominator. This
is a bound on the same interacting operator, not a new Hamiltonian.

Choose the rational point
\[
z_1=\frac{2419}{5000}=0.4838,\quad
s=\frac{81}{5000},\quad t=\frac{17581}{5000}.
\]
Its first pivot -z_1 is negative. The second pivot is
\[
P_1=3-z_1-\frac{382}{9216s}
              -\frac1{9216t}+\frac1{48z_1}.
\]
The following written rational inequalities suffice:
\[
\frac{382}{9216s}
=\frac{1910000}{746496}<\frac{2559}{1000},
\quad
\frac1{48z_1}=\frac{5000}{116112}>\frac{43}{1000},
\quad
\frac1{9216t}<\frac1{25000}.
\tag{13}
\]
For the first, the positive cross-product difference is
1910283264-1910000000=283264. For the second it is
5000000-4992816=7184. For the third, t>3 and
3*9216=27648>25000.

It follows that
\[
\boxed{
P_1>
\frac{12581}{5000}-\frac{2559}{1000}
-\frac1{25000}+\frac{43}{1000}
=\frac1{6250}>0 .
}
\tag{14}
\]
All 48 winding pivots are positive, since
\[
\frac{1}{48s}=\frac{5000}{3888}<\frac32,
\qquad 3-z_1>\frac52.
\]
The old transverse penalty is smaller than that winding penalty
because 191/9216<1/48; (10) only reduces it. All 191 transverse
pivots are therefore positive too.

The lower Schur comparison has exactly one negative direction at
z_1. The actual Schur matrix dominates it, and the full hidden block
D_A-z_1 is positive. Schur inertia on the full physical Hilbert
space gives
\[
E_1(A_\theta)\ge z_1.
\tag{15}
\]
Use the unchanged YM75 rational enclosures
\[
c_0\le-0.00692844,\qquad E_0(A_\theta)\ge-0.00712209,
\qquad E_1(A_\theta)\le c_1=3.
\]
The exact expression for the retained ground Ritz value remains
c_0=(3-sqrt(109/12))/2. Thus
\[
\Delta=E_1-E_0\ge z_1-c_0\ge0.49072844,
\qquad
\Delta\le3.00712209,
\]
which proves (1), conditional only on the inspected predecessor
working estimates just cited. No independent new physical assumption
or computed vacuum input is introduced.

## 5. What helped, what has not transferred

| Development | Actual benefit | Scope |
|---|---|---|
| Earlier RH source-matrix work, YM75 | [0.4904,3.07874] became [0.49069036,3.00712209] | Same finite cube, coupling and complete cut |
| Recent weighted inverse hierarchy, RYM-T1–T3 | Lower endpoint now 0.49072844 | Same cube and source; no larger retained space |
| Prime-power grouping in RH | Shows why source correlations must be retained together | Prime labels and RH constants are not YM masses |
| Nested returned memory | Preserves the full source when eliminating in stages | Already a Schur identity; no extra global sign follows |
| YM91–YM93 conditional flux work | Actual nonzero memory is bounded with its density drift and residual | Separate vacuum-weighted space; Haar weights cannot be imported unchanged |

In particular, the energy-eight decomposition above is in the
unshifted Haar electric-cut representation of YM74–YM75.
YM93's conditional inverse acts in L2(nu_Q), with the actual
vacuum-weighted density and centered source. Its leading hidden
plaquette frequencies are 3k/4 with k=2 or 3, but its exact generator
also has the density drift. Treating those free frequencies as
eigenvalues of the interacting conditional generator would omit
precisely the residual retained in YM92–YM93. No improved YM93
memory constant is claimed here.

The inspected YM93 working result is volume-uniform:
\[
0\le\theta\le1/100,\quad
d_{93}(\theta)=\frac12-32\theta-520\theta^2\ge16/125,
\]
with relative dynamic memory at most
418609/2073600<1/4 at the endpoint. Its all-coupling
gradient estimate |grad_e log h|<=4theta does not extend
that gap interval by itself.

In the declared physical calibration,
theta=c_B/(c_E g^4) grows without bound when g tends to zero.
The finite-cube improvement (1), at theta=1/48, does not supply
a volume-uniform estimate there, and the YM93 small-theta result
does not cover the continuum weak-coupling trajectory.

Native curvature has not been identified automatically with
Riemann curvature. YM93 explicitly chooses the round SU(2)
metric and distinguishes its background Ricci tensor from the
vacuum-weighted Ricci tensor. YM91 also supplies an engineered
comparison family where its memory estimate succeeds despite
failure of a worst-case weighted-curvature sufficient test.
That family is not the actual Wilson vacuum.

The useful next target is the complete actual conditional
Poisson/flux tensor, or a source-resolved retained spectral
pencil, with controlled costs along the physical trajectory.
An improved finite-volume lower bound is progress on an
estimate; it does not settle uniformity, finite physical
energy, survival of local states, or the four-dimensional
continuum construction.

## 6. Execution and gate accounting

Only this written note and its README entry are added to the
existing draft PR #138, with [skip ci]. The YM branches and
their evidence packets remain unchanged. No subagents,
numerical simulations, independent scientific computation,
certificate replay or workflow dispatch is performed.
Repository-content readback and actual Actions-run counts
are recorded separately in the PR; those metadata checks
are not mathematical certification.

**Closed in the working derivation:** the source-matched
weighted inverse transfer, the actual energy-eight projection,
and the improved finite-cube enclosure (1).
**Unresolved:** a sharper actual YM93 conditional-memory
estimate, volume-uniform control on the physical weak-coupling
trajectory, surviving finite-energy local states and the
four-dimensional continuum mass gap. RH all-width positivity
also remains open; its established unrestricted interval is
unchanged.
**New independent assumption:** none. A prime-to-physical-mass
identification has not been asserted.
