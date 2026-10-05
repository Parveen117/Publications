# YM85 — Exterior reference frames and the cost of physical response

Monty Dabas. 6 October 2026. Python 3.12 evidence packet.

## What this changes after YM84

YM84 correctly removes *bare* open-path markers by Gauss averaging and
gives a physical-sector projection criterion. Its isolated plaquette
example has a much better physical resampling gap than the unrestricted
one. A full lattice contains additional reference information, however:
exterior links can dress local markers into gauge-invariant observables.
This stage determines exactly what that does to the local criterion.

1. For any finite SU(2) lattice with smooth strictly positive invariant
   density, a connected exterior graph with at least two independent
   cycles generically supplies a complete gauge reference frame. The
   conditional response operator has the **same norm** on the physical
   Hilbert space as on the full space. The optimal integrated physical
   block Poincare floor also equals the unrestricted conditional floor.
2. This applies to every pair of single-link blocks in the actual
   homogeneous cubic vacuum, at every finite coupling and volume.
   Consequently YM84's single-link physical row test has exactly the
   same coefficients as YM83's unrestricted test. Gauge projection
   alone cannot improve this particular test.
3. In that same actual vacuum, an exterior frame produces three physical
   plaquette markers with an exact local energy matrix. Their partial
   block Rayleigh quotients are at least
   \((2/9)\sqrt\theta-1\) for \(\theta\ge25\). This is a bound for
   a stated three-dimensional marker space, not for all physical states.
4. Two of the normalized frame markers are in physical L2 but outside
   the full electric form domain. Smooth cutoffs have an explicit
   logarithmically divergent electric cost. The first marker does have
   finite full electric energy. Thus retaining a response channel and
   producing a finite-energy physical state are different operations.

These results use the actual vacuum measure, without solving for its
density. They identify a missing energy obligation in the gauge-resolved
route. They do not establish its remaining weak-coupling block bounds,
a nonzero continuum limit, or a positive four-dimensional mass gap.
Tree gauge and body-frame coordinates are established constructions;
the contribution here is their explicit use in the YM83--YM84 conditional
norms and energy-domain tests. No priority claim for gauge fixing is made.

## 1. Model, norms and partial forms

Let \(\mathcal V,\mathcal E\) be a finite connected oriented graph,
with one SU(2) variable per edge and normalized product Haar measure
\(\Phi\). Use the gauge action

\[
 U_e\longmapsto g_{s(e)}U_eg_{t(e)}^{-1}.
 \tag{1}
\]

Let \(d\mu=w\,d\Phi\), where w is smooth, strictly positive and
gauge invariant. On this finite compact configuration manifold,

\[
 0<m_\mu:=\min w\le w\le M_\mu:=\max w<\infty.
 \tag{2}
\]

The intended application is the unchanged homogeneous periodic cubic
Hamiltonian, with side L>=4,

\[
 A_\theta=L_E-\theta\sum_p W_p,\qquad
 H=A_\theta-E_0,\qquad w=h^2,\qquad \mu=h^2\Phi.
 \tag{3}
\]

Its unique positive normalized ground state h has all the properties
above. The native link Casimir is 3/4; equivalently the link sphere has
radius 2. The ground-state transform identifies the electric Dirichlet
form \(\mathcal E(f)=\mu\sum_e\Gamma_e(f)\) with the form of H
on vectors hf. Gauss averaging is \(\Pi_G\), with range \(\mathcal H_G\).

For a block B write \(P_B=\mathbb E_\mu[\cdot\mid U_{\bar B}]\)
and \(Q_B=I-P_B\). For disjoint B,C define

\[
 T_{B,C}=P_BP_C-P_{B\cup C},\qquad
 c_{B,C}=\|T_{B,C}\|,\qquad
 c^G_{B,C}=\|T_{B,C}|_{\mathcal H_G}\|.
 \tag{4}
\]

These are norms of the full exterior-conditioned operator, not of an
unconditioned marginal. Define the best physical block floor by

\[
 \rho_B^G=\inf_{\substack{f\in C^\infty\cap\mathcal H_G\\Q_Bf\ne0}}
     \frac{\mu\Gamma_B(f)}{\|Q_Bf\|_\mu^2},
 \qquad \Gamma_B=\sum_{e\in B}\Gamma_e.
 \tag{5}
\]

At fixed exterior z, let \(\lambda_B(z)\) be the optimal Poincare
constant for the positive conditional density \(\mu_B(\cdot\mid z)\),
using the native product-link gradient. Put
\(\rho_B=\operatorname*{ess\,inf}_z\lambda_B(z)\).
This is a conditional-measure constant; it is not an isolated block
Hamiltonian's gap with a different vacuum.

The **partial form domain** for B completes smooth functions using
\(\|f\|_\mu^2+\mu\Gamma_B(f)\). It does not charge derivatives
outside B. It is larger than the full electric form domain. Smooth
invariant exterior cutoffs are admissible in (5) even if their exterior
gradient becomes large along an approximating sequence.

## 2. Theorem 85.1 — A complete frame made only from the exterior

Fix a set S of internal edges. Assume its exterior graph
\((\mathcal V,\mathcal E\setminus S)\) is connected and has cycle rank
at least two. Choose an exterior spanning tree with root o. Let T_v be
its transport from o to v. For every chord or internal edge put

\[
 Z_e=T_{s(e)}U_eT_{t(e)}^{-1}.
 \tag{6}
\]

Every Z_e transforms by the same conjugation at o. Choose two distinct
exterior chords, with holonomies
\(A=(a_0,\mathbf a)\), \(B=(b_0,\mathbf b)\) in unit-quaternion
coordinates. On the set

\[
 r=|\mathbf a|>0,\qquad
 s=|\mathbf b-(\mathbf b\cdot\mathbf n)\mathbf n|>0,
 \qquad \mathbf n=\mathbf a/r,
 \tag{7}
\]

define the oriented body frame

\[
 \mathbf e_1=\mathbf n,\qquad
 \mathbf e_2=\frac{\mathbf b-(\mathbf b\cdot\mathbf n)\mathbf n}{s},
 \qquad \mathbf e_3=\mathbf e_1\times\mathbf e_2,
 \qquad R=(\mathbf e_1,\mathbf e_2,\mathbf e_3)\in SO(3).
 \tag{8}
\]

For \(Z_e=(z_{e,0},\mathbf z_e)\), the canonical coordinates

\[
 \widehat Z_e=(z_{e,0},R^T\mathbf z_e)
 \tag{9}
\]

are gauge invariant. At fixed exterior, the map from the internal
variables to \(\widehat Z_S\) is a product of Haar-preserving native
isometries. Its dependence on the exterior uses no internal edge.

The regular set (7) has full \(\mu\)-measure. Its exterior gauge
stabilizer is just the common central sign, which acts trivially on
every internal edge. Equations (6)--(9) separate gauge orbits there.

**Proof.** Tree transport transforms as
\(T_v\mapsto g_oT_vg_v^{-1}\), proving the common-conjugation rule.
SU(2) conjugation rotates imaginary quaternions in SO(3), so R transforms
by that rotation and (9) is invariant. Two noncollinear imaginary
vectors determine an oriented frame. Their common SU(2) centralizer is
\(\{\pm I\}\); the tree then determines every other vertex gauge
transformation from the root value. A common central sign does nothing
to internal links.

Replacing non-tree edges by (6), while keeping tree variables, preserves
product Haar measure. The two selected exterior chord variables are
independent Haar variables in these coordinates. The event r=0 is a
null set, and, for any r>0, the event s=0 is a circle inside the B
three-sphere and is also null. Absolute continuity gives the assertion
for \(\mu\). No uniform lower bound on r or s is asserted.

Finally the tree sets its own links to identity; the frame fixes the
remaining common conjugation uniquely up to its ineffective sign.
The canonical exterior and the freely varying canonical internal
variables therefore specify the orbit. At fixed exterior, (6) is a
left/right translation and (9) is conjugation, proving the isometry
assertion. QED.

### The exterior condition on the cubic lattice

Deleting any set S of at most four links from the L>=4 periodic cubic
lattice leaves it connected. Each deleted edge has four three-edge
plaquette detours, from the positive and negative directions of the two
transverse axes. Those four detours have disjoint edge sets. At most
three other deleted links can obstruct them, so at least one detour
survives. Replace each deleted edge in an original path by a surviving
detour to obtain an exterior path.

Its cycle rank is

\[
 |\mathcal E\setminus S|-|\mathcal V|+1
   =2L^3-|S|+1\ge125.
 \tag{10}
\]

Thus (7)--(9) apply both to S equal to one or two arbitrary links and
to S equal to all four links of one plaquette. The theorem is not
asserted for arbitrary large blocks whose exterior disconnects or
does not contain the required independent cycles.

## 3. Theorem 85.2 — The physical restriction preserves local worst cases

If the exterior of \(B\cup C\) satisfies Theorem 85.1, then

\[
 \boxed{c^G_{B,C}=c_{B,C}
   =\operatorname*{ess\,sup}_z c_{B,C}(z),}
 \tag{11}
\]

where the rightmost coefficient is the conditional maximal correlation
of B and C at fixed exterior z. If the exterior of B satisfies that
theorem, then

\[
 \boxed{\rho_B^G=\rho_B
   =\operatorname*{ess\,inf}_z\lambda_B(z).}
 \tag{12}
\]

Both identities hold for the actual finite-volume vacuum (3), at every
finite theta. In particular, for every pair of single-link blocks,
the YM84 physical row sum equals the YM83 unrestricted row sum.

**Proof of (11).** Disintegrate over all exterior variables. Conditional
expectation over either block is a fiber operator. On each fiber,
\(P_BP_C-P_{B\cup C}\) has norm equal to conditional maximal
correlation. The norm of this decomposable operator is the essential
supremum of its fiber norms. This also proves the upper bound for its
physical restriction.

For the reverse inequality, use the exterior frame. An arbitrary
function of the canonical internal variables and canonical exterior
is a physical function. Thus the physical Hilbert space contains the
entire conditional fiber L2 space over every regular exterior orbit,
with no remaining stabilizer restriction. Gauge-related exteriors
have unitarily conjugate fiber operators. Disintegration over exterior
orbits therefore retains exactly the same essential supremum.

One can make the norm argument without choosing extremizing vectors
measurably over all exteriors. Pick a regular exterior where the norm
is within epsilon of that supremum and choose a smooth approximate
fiber witness. The positive smooth conditional densities make the
corresponding norm ratio continuous locally. Transport the witness
with (9) and multiply by a smooth invariant exterior cutoff supported
inside that regular neighborhood in orbit space. The ratio remains
within an arbitrarily small error. Such functions are globally smooth
after extending by zero. Let epsilon tend to zero.

**Proof of (12).** Integrating each conditional Poincare inequality
gives \(\rho_B^G\ge\rho_B\). To reverse it, choose a smooth fiber
test whose conditional quotient is arbitrarily close to \(\rho_B\)
at a regular exterior. Such an exterior may be chosen because the
regular set has full measure, and smooth fiber tests are dense in the
conditional form domain. Apply the same frame transport and invariant
exterior localization. Subtract its conditional mean.

Native B derivatives see only the internal function: both the frame
and the cutoff are exterior functions. The isometries (6)--(9) preserve
the conditional gradient. The integrated quotient consequently
approaches the selected conditional quotient. Exterior derivatives
can be large, but (5) does not charge them. This proves the reverse
inequality. QED.

### What does and does not follow

Bare open-path markers in YM84 still have zero Gauss average. The new
markers use *additional exterior information*. The isolated cycle in
YM84 has no such exterior frame, so its physical resampling gap 2 and
its quantum gap bounds remain valid.

Equations (11)--(12) do not say that the physical quantum gap equals
the unrestricted quantum gap. They concern two specified conditional
operator norms and a specified partial energy inequality. Nor do they
show that either row sum is large at weak coupling: YM84's
unconditioned open-path calculation cannot replace the conditioned
coefficient in (11).

The identities do show that simply restricting those local constants
to physical functions gives no improvement whenever the exterior
supplies the frame. An improvement must use information absent from
that worst-case test, such as full electric energy or a sharper
estimate of the physical resampling generator itself.

## 4. Theorem 85.3 — Actual dressed plaquette stiffness

Remove all four edges of an elementary plaquette p and build the
frame from its exterior, as above. Let its holonomy, transported to
the tree root, be \(X=(W_p,\mathbf x)\). Put

\[
 F_k=\mathbf e_k\cdot\mathbf x,\qquad
 \mathbf F=(F_1,F_2,F_3)^T.
 \tag{13}
\]

These bounded functions are gauge invariant and in the partial form
domain for p. Almost everywhere their native local energy is exactly

\[
 \boxed{\Gamma_p(F_i,F_j)=\delta_{ij}-F_iF_j,\qquad
        |\mathbf F|^2=1-W_p^2.}
 \tag{14}
\]

Define \(M=\mu(\mathbf F\mathbf F^T)\),
\(C=\operatorname{Cov}_\mu(\mathbf F)\), and
\(G_{ij}=\mu\Gamma_p(F_i,F_j)\). With
\(q=\mu(1-W_p^2)\) and \(\delta=\mu(1-W_p)\),

\[
 G=I_3-M,\qquad
 0<C\le M\le qI_3\le2\delta I_3,\qquad
 G\ge(1-2\delta)I_3.
 \tag{15}
\]

For nonzero constant \(v\in\mathbb R^3\), if \(2\delta<1\),

\[
 \frac{\mu\Gamma_p(v\cdot\mathbf F)}
      {\operatorname{Var}_\mu(v\cdot\mathbf F)}
 \ge\frac{1-2\delta}{2\delta}.
 \tag{16}
\]

YM76's actual homogeneous-vacuum bound
\(\delta\le9/(4\sqrt\theta)\), \(\theta\ge25\), gives

\[
 \boxed{\frac{\mu\Gamma_p(v\cdot\mathbf F)}
              {\operatorname{Var}_\mu(v\cdot\mathbf F)}
       \ge\frac29\sqrt\theta-1.}
 \tag{17}
\]

For theta>=81 this is at least \(\sqrt\theta/9\).

**Proof.** At fixed exterior the frame is constant. Differentiating
the plaquette holonomy in any one of its four distinct links gives
the native SU(2) gradient of X, up to left/right isometries. A linear
unit-quaternion coordinate has native gradient pairing
\((\delta_{ij}-F_iF_j)/4\). Sum over the four links to get (14).
All matrix inequalities in (15) follow from positivity and
\(\operatorname{tr}M=q\le2\delta\). The conditional density of
the plaquette links is strictly positive; its holonomy can vary over
all SU(2) while holding three links fixed. Therefore every nonzero
v has positive variance. This also proves C is positive definite.
Equations (16)--(17) follow by the Rayleigh quotient. QED.

Exterior cutoffs vanishing near r=0 and s=0 approximate (13) in L2
and in the p-partial form norm. They commute with every p derivative,
so (14) passes to that partial form closure. This does not yet put
all three markers in the full electric form domain.

The inequality is a high-energy lower quotient on a specified marker
space. It is neither a lower bound on the full quantum spectral gap
nor a proof that these vectors have no low-energy spectral component.
A large first spectral moment alone does not locate all spectral mass.

## 5. Theorem 85.4 — Two frame coordinates have infinite electric energy

Let e_B be the exterior chord whose holonomy defines B in (7). The
tree, A and the unrotated target plaquette X do not depend on this
chord variable. Native differentiation of the frame gives exactly

\[
 \boxed{\Gamma_{e_B}(F_2)=\frac{F_3^2}{4s^2},\qquad
        \Gamma_{e_B}(F_3)=\frac{F_2^2}{4s^2}.}
 \tag{18}
\]

Consequently, for **every** finite-volume vacuum (3),

\[
 F_2,F_3\in\mathcal H_G\cap L^\infty(\mu),\qquad
 F_2,F_3\notin\operatorname{Dom}\mathcal E.
 \tag{19}
\]

The first marker F_1 *does* belong to the full form domain. Thus its
centered physical vector h(F_1-\mu F_1) is an actual finite-energy
state and satisfies the lower quotient (17), with full energy in
place of partial energy if desired. This does not prove that the
state has bounded energy in physical continuum units.

**Proof.** Write \(\mathbf b=b_\parallel\mathbf n+s\mathbf e_2\).
At fixed n and x, the Euclidean b-gradient of F_2 is
\(F_3\mathbf e_3/s\). It is tangent to the unit quaternion sphere
because \(\mathbf b\cdot\mathbf e_3=0\). The native metric supplies
the factor 1/4 in (18). Rotating \(\mathbf e_3=\mathbf n\times
\mathbf e_2\) gives the second identity with the same squared norm.

In the tree-Haar coordinates, A and B are independent Haar variables.
At fixed exterior and three target links, X is Haar as the fourth
target link varies. Thus \(\Phi(F_2^2\mid A,B)=
\Phi(F_3^2\mid A,B)=1/4\). For a fixed n, the squared perpendicular
radius \(u=s^2\) under Haar B is uniform on [0,1]. This follows by
splitting the unit three-sphere into its two orthogonal coordinate
planes, or from the Beta(1,1) density. Hence for 0<epsilon<1,

\[
 \int_{s\ge\epsilon}\Gamma_{e_B}(F_2)\,d\Phi
   =\frac1{16}\int_{\epsilon^2}^1\frac{du}{u}
   =\frac18\log\frac1\epsilon,
 \tag{20}
\]

and likewise for F_3. The density floor (2) proves divergence in
the actual measure. Any Sobolev weak derivative would agree with
the displayed classical derivative on every regular open set, so
this divergence rules out full form membership; assigning arbitrary
values on the null singular set cannot repair it.

For F_1 only the normalization \(\mathbf a/|\mathbf a|\) is needed.
Its derivatives are bounded by a finite graph-dependent multiple of
\(1+1/r\). Haar A satisfies

\[
 \int r^{-2}dA=\frac2\pi\int_0^\pi d\alpha=2.
 \tag{21}
\]

The bounds (2) and compactness give finite total gradient energy.
Smooth radial cutoffs also verify weak differentiability: the cutoff
energy near r=0 is O(epsilon), since Haar volume there is
O(epsilon^3). Thus F_1 lies in the full form closure. QED.

### Quantitative cost of regularization

Take smooth invariant cutoffs F_(k,epsilon) that agree with F_k on
\(\{r\ge1/2,\ s\ge\epsilon\}\), and vanish near both singular
sets. They may be chosen to converge to F_k in L2 and in the
p-partial form norm as epsilon tends to zero. For k=2,3,

\[
 \boxed{\mathcal E(F_{k,\epsilon})
       \ge\frac{m_\mu}{16}\log\frac1\epsilon.}
 \tag{22}
\]

Indeed (21) and Markov give \(\Phi(r\ge1/2)\ge1/2\); A is
independent of B and the target Haar holonomy. Apply (20) on that
subset, where the cutoff equals the original marker. The coefficient
\(m_\mu\) is positive at each finite regulator but is not claimed
uniform in coupling or volume. Equation (22) is a fixed-regulator
energy-domain obstruction, not a continuum mass-gap bound.

The same distinction explains the localization argument in (12).
The partial block form does not see the exterior gradient required
to select the frame. The full physical Hamiltonian does see it.

## 6. Consequence for the continuum route

The YM84 physical projection theorem remains valid. For local blocks
with a complete exterior frame, however, its two optimally defined
inputs coincide with the unrestricted conditional inputs. It cannot
gain a stronger single-link estimate solely from the word "physical."

The exact identities (14) and (18) identify where additional information
can enter: the block geometry and the electric cost of the exterior
reference. A next criterion must control those together, or estimate
the full physical resampling generator beyond its pairwise row bound.
No such new uniform weak-coupling criterion is assumed here.

The physical clock is still
\(\kappa\sqrt\theta=c/a\). The finite-energy marker F_1 in this
stage has physical Rayleigh quotient at least
\(\kappa((2/9)\sqrt\theta-1)\), and for theta>=81 at least
\(c/(9a)\). It therefore does not supply the bounded-energy
continuum-state sequence requested by YM84. Collective smearing,
actual-vacuum correlation control and state survival remain required.

## 7. Evidence and reproduction

The exact Python controls verify native quaternion products and gauge
transport, complete exterior frames, disjoint lattice detours, a finite
group analogue of conditional norm preservation, the plaquette energy
matrix, tangent derivatives in (18), Haar singularity constants and
the stated weak-coupling arithmetic. The finite group is a separate
comparison fixture; it is not the SU(2) vacuum. Exact finite checks
do not prove the infinite-dimensional norm identities or the continuum
claim; those identities are established by the written arguments above.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym85_exterior_frames_energy.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym85_exterior_frames_energy.py -v
```

The manifest preserves previous research inputs and evidence. No
independent external review or proof-assistant verification is claimed.

## References and lineage

- [YM84](YM84_GAUGE_RESOLVED_CONTINUUM_SCALING.md): physical projection
  criterion, bare open-path tower and the isolated-cycle comparison.
- [YM83](YM83_ACTUAL_VACUUM_CONDITIONAL_RESPONSE.md): the actual
  conditional measure, block Poincare floors and response row criterion.
- [YM76](YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS.md): native link
  normalization and the volume-uniform actual plaquette defect estimate.
- F. Fuerstenberg, G. Rudolph and M. Schmidt, *Defining relations for
  the orbit type strata of SU(2)-lattice gauge models*, J. Geom. Phys.
  119 (2017), 66--81, [arXiv:1702.01047](https://arxiv.org/abs/1702.01047).
  See its tree reduction and residual common-conjugation action.
- D. M. Grabowska, C. F. Kane and C. W. Bauer, *Fully gauge-fixed SU(2)
  Hamiltonian for quantum simulations*, Phys. Rev. D 111 (2025), 114516,
  [arXiv:2409.10610](https://arxiv.org/abs/2409.10610),
  [DOI](https://doi.org/10.1103/PhysRevD.111.114516). Cited for the
  established maximal-tree and body-frame construction, not for the
  conditional norm or partial/full form conclusions proved here.
