# Native loop returns, Cartan curvature and the cosmological coefficient

Author: Monty Dabas. Development R4, 30 September 2026. Results LR-1–LR-8.

**Result.** A declared quadratic, oriented readout of small native matrix loop
returns has a curvature-square continuum limit. Combine the admitted Lorentz
connection and coframe in a Cartan connection. The real native coefficient
algebra then fixes the relative Einstein and cosmological terms of this
quadratic class:

\[
\boxed{\Lambda=-12\sigma u^2=-\frac{3\sigma}{\ell^2},\qquad
\kappa=-\frac1{2\beta\sigma u^2},\qquad
\kappa\Lambda=\frac6\beta,\qquad \ell=\frac1{2u}>0.}
\tag{1}
\]

Here u is an admitted inverse-length scale, β a nonzero readout coefficient,
and σ the square sign of the translation-like coefficient generators. The
original real four-component realization permits σ=+1, giving Λ<0. One
additional two-state native factor with R²=−I permits σ=−1, giving Λ>0, with
all matrices still real. Eight real components are sufficient, and two copies
are minimal **among copies of the fixed four-component Lorentz module**.

Lorentz invariance below means the connected proper Lorentz group, preserving
the chosen native orientation. There are three invariant symmetric quadratic curvature pairings in
this class. They give a Holst coefficient γ=(α−ζ)/β. The unprojected ordinary
trace has α=ζ and its entire contribution is bulk-inactive; setting ζ=0 by an
even-curvature readout restores an independent Holst term. Requiring the larger
Cartan internal symmetry instead leaves only the bulk-inactive ordinary trace.
Thus the native J readout/Lorentz reduction is substantive input.

This is a constructive extension of MG, not primitive-only selection of the
readout, continuum, scale, module or physical constants. The underlying
curvature-square gravity mechanism is the established MacDowell–Mansouri/Cartan
construction [1–3]. The development here makes its coefficient realization,
loop readout, sign restriction, native-factor extension and thermo variation
interface explicit and reproducible. It does not claim a new discovery of de
Sitter space, the MacDowell–Mansouri action or the Cartan interpretation of
curvature. The proper one-sided seam still has no faithful finite realization
in these matrices.

## LR-1. A precise loop-return process with a curvature-square limit

Admit a smooth real matrix connection 𝓐 on an oriented four-dimensional patch.
Use the convention in which based, positively oriented coordinate squares of
side h have holonomy

\[
U_{\mu\nu}(h)=I+h^2\mathcal F_{\mu\nu}+O(h^3),\qquad
\mathcal F=d\mathcal A+\mathcal A\wedge\mathcal A.
\tag{2}
\]

This is the smooth shrinking-loop interface of SH-2. It requires an identified
connection, orientation and loop scale. A single native sign or a
noncontractible return does not supply those inputs. Every loop used in the
following pairing is based at the **same point**; returns at different points
must first be transported to a common comparison frame.

Let D_μν=U_μν−I. For a fixed internal readout Q define the local oriented
quadratic record

\[
\begin{split}
\mathscr R_h(Q)=\tau\{Q(&D_{01}D_{23}+D_{23}D_{01}
-D_{02}D_{13}-D_{13}D_{02}\\
&+D_{03}D_{12}+D_{12}D_{03})\}.
\end{split}
\tag{3}
\]

τ is the normalized real matrix trace. Substituting (2) proves

\[
h^{-4}\mathscr R_h(Q)
=[\tau(Q\mathcal F\wedge\mathcal F)]_{0123}+O(h).
\tag{4}
\]

The two product orders in (3) are retained because matrix factors need not
commute. The six terms are exactly the coefficient of the wedge square; this
is an oriented pairing of complementary loop planes, not a sum of positive
return-defect norms. Under a uniform smooth expansion on a compact patch,
summing (3) once per coordinate cell is a Riemann sum for
∫τ(Q𝓕∧𝓕), with error O(h). Smooth fields with uniformly bounded derivatives
on a neighborhood of the patch provide such an expansion. This derives the
leading continuum action from the **declared loop functional**. It does not
assert exact finite-lattice diffeomorphism symmetry or convergence of every
finite-lattice stationary point.

More generally apply a constant symmetric internal bilinear pairing to the
two loop defects in each term. The same expansion yields its wedge pairing
on 𝓕. LR-4 classifies the Lorentz-invariant choices on the admitted ten-dimensional
Cartan algebra. Selecting this quadratic, local, oriented process law remains
an assumption; higher loop moments and nonlocal readouts are not excluded by
the native arithmetic alone.

## LR-2. The Cartan connection from the existing native matrices

Retain the MG convention

\[
\{\Gamma_a,\Gamma_b\}=2\eta_{ab}I,\quad
\eta=\operatorname{diag}(-1,1,1,1),\quad
J=\Gamma_0\Gamma_1\Gamma_2\Gamma_3,\quad J^2=-I.
\]

The six bivectors commute with J and the four vectors anticommute with J.
There are three explicit real realizations used below:

| Realization | Vector matrices Y_a | Bivectors B_ab | Readout Ĵ | σ |
|---|---|---|---|---|
| Original four-component module | Γ_a | Γ_aΓ_b | J | +1 |
| Extra native factor K²=I | K⊗Γ_a | I₂⊗Γ_aΓ_b | I₂⊗J | +1 |
| Extra native factor R²=−I | R⊗Γ_a | I₂⊗Γ_aΓ_b | I₂⊗J | −1 |

The latter two act on eight real components. These R and K are exactly the
real two-by-two native coefficient matrices already admitted in NP, not newly
postulated complex scalars. B_ab is defined for a<b in the table and extended
by B_ba=−B_ab, B_aa=0. In all three cases

\[
[Y_a,Y_b]=2\sigma B_{ab},\qquad
[B_{ab},Y_c]=2(\eta_{bc}Y_a-\eta_{ac}Y_b).
\tag{5}
\]

The ten independent generators close. With M_ab=B_ab/2 and P_a=Y_a/2,
[P_a,P_b]=σM_ab, which is the orthogonal Lie algebra with extended internal
metric diag(η,−σ). Thus σ=+1 is the anti-de Sitter algebra so(2,3), and σ=−1
the de Sitter algebra so(1,4). This is an algebraic group identification in a
declared representation, not a selection of the physical spacetime dimension.

Let e^a be an invertible oriented coframe, ω^{ab}=−ω^{ba} an independent
Lorentz connection and u>0 constant. Define

\[
\widehat\Omega=\frac14\omega^{ab}B_{ab},\qquad
\mathcal A=\widehat\Omega+u e^aY_a,
\tag{6}
\]
\[
\mathcal C=\frac14 R^{ab}B_{ab},\quad
\mathcal B=e^a\wedge e^b B_{ab},\quad
\mathcal T=T^aY_a,\quad \rho=\sigma u^2,
\]

where R^{ab}=dω^{ab}+ω^a{}_c∧ω^{cb}, and T^a=de^a+ω^a{}_b∧e^b.
Repeated antisymmetric indices include both orders, as in MG. In the
four-component realization 𝓑=E∧E. Expanding d𝓐+𝓐² using (5) gives

\[
\boxed{\mathcal F=\mathcal C+u\mathcal T+\rho\mathcal B.}
\tag{7}
\]

The normalized trace is tr/4 or tr/8 as appropriate; τ₈(I₂⊗X)=τ₄(X).
Thus adding the native factor does not introduce a hidden factor of two into
the action normalization. The MG thermo adapter can supply e and all its
variations; the connection ω and constant u are still independent admitted
data.

## LR-3. Why the original real four-component sector has only one sign

Suppose four real matrices Y_a on the **same** original module transform as a
Lorentz vector under its six fixed bivectors. The complete solution is

\[
Y_a=r\Gamma_a+sJ\Gamma_a,\qquad r,s\in\mathbb R.
\tag{8}
\]

To see completeness, use the sixteen-element Clifford basis
I,J,Γ_a,JΓ_a,Γ_aΓ_b (a<b). Rotation covariance forces Y_0 to lie in the span
of I,J,Γ_0,JΓ_0; the bivectors have no rotation-invariant component.
Write Y_0=A I+B J+rΓ_0+sJΓ_0. Then [Γ_0Γ_i,Y_0]=2Y_i fixes
Y_i=rΓ_i+sJΓ_i. The relation [Γ_0Γ_i,Y_i]=2Y_0 forces A=B=0.
Hence only the two constants in (8) remain. The executable check also solves
all 64 real unknown entries directly:
the intertwiner constraint matrix has rank 62, and these two independent
solutions exhaust its kernel.

Since JΓ_aJ=Γ_a, cross terms cancel and

\[
Y_aY_b=(r^2+s^2)\Gamma_a\Gamma_b.
\tag{9}
\]

Every nonzero real embedding therefore has a **positive** square factor.
Replacing Γ_a by JΓ_a cannot produce the opposite sign. The obstruction concerns
this fixed Lorentz module and equivariant linear embedding, not every possible
UGD representation or higher-dimensional theory.

In two copies of this module, the independent native R factor commutes with
the original gamma matrices and squares to −I. Then
(R⊗Γ_a)(R⊗Γ_b)=−I₂⊗Γ_aΓ_b. This realizes the missing sign, as (5) shows.
One copy is excluded by (9), and two copies suffice. This is the stated
minimality result among integral copies of the original module, with no claim
that eight is minimal over every representation of every relevant group.

## LR-4. Classifying the quadratic readout and reducing it to MG

For X in the ten-dimensional Cartan algebra, define the even and odd native
readouts

\[
X_+=\frac12(X-\widehat JX\widehat J),\qquad
X_-=\frac12(X+\widehat JX\widehat J).
\tag{10}
\]

They project onto the Lorentz bivectors and translation-like vectors. In
particular 𝓕_+=𝓒+ρ𝓑 and 𝓕_−=u𝓣. Every real, constant, Lorentz-invariant
symmetric bilinear pairing on this algebra is a linear combination of
τ(X_+Y_+), τ(ĴX_+Y_+) and τ(X_−Y_−).

**Proof of completeness.** Under spatial rotations, write the generators as
three boosts, three rotations, three spatial translations and one temporal
translation. Rotation invariance permits seven scalar pairing coefficients:
boost/boost, rotation/rotation, their cross pairing, spatial/spatial,
temporal/temporal, boost/spatial and rotation/spatial. Boost invariance sets
rotation/rotation to minus boost/boost, temporal/temporal to minus
spatial/spatial, and both mixed translation pairings to zero. Three coefficients
remain; the three displayed trace pairings are independent and invariant.
This also follows by solving the 55-unknown symmetric bilinear problem: its
Lorentz-invariance constraints have rank 52 for either σ. □

Consequently the general action in this constant quadratic wedge class is

\[
S=\int\left[\alpha\tau(\mathcal F_+\wedge\mathcal F_+)
+\beta\tau(\widehat J\mathcal F_+\wedge\mathcal F_+)
+\zeta\tau(\mathcal F_-\wedge\mathcal F_-)\right].
\tag{11}
\]

No spacetime Hodge star, inverse coframe, variable coefficients or additional
fields are included. Such additions would define a broader class.

There is also a single matrix stationary equation. Define the linear curvature
response

\[
\mathcal K(\mathcal F)=\alpha\mathcal F_+
+\beta\widehat J\mathcal F_++\zeta\mathcal F_-.
\]

This map is self-adjoint for the trace pairing and takes the Cartan algebra
into itself. With δ𝓕=D_𝓐δ𝓐, compact-support integration by parts gives
δS=−2∫τ(D_𝓐𝓚(𝓕)∧δ𝓐). The trace pairing on the ten generators is
nondegenerate, and u≠0 makes independent δe, δω span every δ𝓐. Consequently

\[
\boxed{D_{\mathcal A}\mathcal K(\mathcal F)=0.}
\tag{11a}
\]

For the unprojected family ζ=α, the Bianchi identity D_𝓐𝓕=0 reduces (11a),
when β≠0, to the native readout law

\[
\boxed{D_{\mathcal A}\{\widehat J,\mathcal F\}=0.}
\tag{11b}
\]

These are conditional field equations of the declared loop process, rather
than a separate insertion of an Einstein tensor. Their even component is
ρ[(α−ζ)I+βĴ]D_Ω𝓑=0. For real coefficients and β≠0, its multiplier is
invertible; MG's torsion argument gives T=0. The odd component then supplies
the coframe equation. The following density reduction identifies its exact
Einstein normalization and cosmological coefficient.

Put p=τ(Ĵ𝓑∧𝓒), h=τ(𝓑∧𝓒) and v=τ(Ĵ𝓑∧𝓑). The MG trace identities give

\[
p=-\frac14\epsilon_{abcd}e^a\wedge e^b\wedge R^{cd},\quad
h=-\frac12e^a\wedge e^b\wedge R_{ab},\quad v=-24\,\mathrm{vol}_e.
\]

Odd Clifford traces vanish, τ(𝓑²)=0 and τ(Ĵ𝓣²)=0. Moreover
τ(𝓣²)=σ T_a∧T^a. Using the directly differentiated Nieh–Yan identity
d(e_a∧T^a)=T_a∧T^a+2h yields the **full density identity**

\[
\begin{split}
\mathcal L={}&\alpha\tau(\mathcal C\wedge\mathcal C)
+\beta\tau(\widehat J\mathcal C\wedge\mathcal C)
+\zeta\rho\,d(e_a\wedge T^a)\\
&+2\rho[\beta p+(\alpha-\zeta)h]+\beta\rho^2v.
\end{split}
\tag{12}
\]

The first two terms have only boundary first variations for a Lorentz
connection; the third is exact. They are bulk-inactive under compactly
supported variations, without being assumed zero pointwise or globally.
For β≠0 and u≠0 the remaining terms are exactly the MG action, with

\[
\kappa=-\frac1{2\beta\rho},\qquad
\gamma=\frac{\alpha-\zeta}{\beta},\qquad
\Lambda=-12\rho.
\tag{13}
\]

All numerical factors in (1) now follow from the chosen matrix normalization.
MG-3–MG-5 imply zero torsion and the **full** vacuum equations
G_μν+Λq_μν=0, including when e is supplied by its sixteen-channel thermo
adapter with complete variations. Adding the independent declared Maxwell
sector gives the MG Einstein–Maxwell equations with the same Λ. This step
does not supply the still-missing curved-spin-matter completion.

For the unprojected two-trace readout
S=∫[ατ(𝓕²)+βτ(Ĵ𝓕²)], one has ζ=α. Therefore γ=0 and

\[
\tau(\mathcal F^2)=\tau(\mathcal C^2)+\rho\,d(e_a\wedge T^a).
\tag{14}
\]

Its ordinary trace is entirely bulk-inactive. It is incorrect to retain its
apparent Holst term while dropping the compensating torsion-square term.
Conversely a readout that first retains only 𝓕_+ sets ζ=0 and permits
γ=α/β. Thus zero Holst coupling follows from the unprojected readout choice,
not from every native quadratic readout. The Λ relation in (13) survives this
three-weight extension.

## LR-5. What symmetry and positivity do, and do not, select

The actions in (11) are Lorentz invariant: Ĵ commutes with the Lorentz
generators, and the even/odd decomposition is preserved. They are differential
four-form actions, with the usual coordinate covariance. They need not be
invariant under the larger internal Cartan group when Ĵ is held fixed.

Requiring invariance under all ten internal generators restricts (11) to
β=0 and α=ζ. One way to see the second restriction is to apply a vector
generator to the pairing of a vector and a bivector: the commutators (5)
relate their trace weights and force equality. A boost/rotation cross pairing
with a vector generator forces β=0. Equivalently the full symmetric-bilinear
constraint matrix has rank 54 out of 55 unknowns for either sign, leaving only
τ(XY). By (14), its action has no gravitational bulk equations.

An explicit control makes the J issue visible. In the original real module set

\[
\mathcal F=\Gamma_0\Gamma_1\,dx^0\wedge dx^1
+\Gamma_2\,dx^2\wedge dx^3,
\qquad \delta\mathcal F=[\Gamma_3,\mathcal F].
\]

Then [δτ(J𝓕∧𝓕)]_0123=4, whereas the ordinary-trace variation is zero.
Every Lorentz-generator variation preserves the J readout. Thus obtaining
gravity here requires retaining the native J readout/Lorentz reduction; full
unbroken Cartan internal invariance of a constant quadratic pairing would
remove it. Promoting J to a transforming field is a possible larger theory,
but requires new constraints and dynamics which are not supplied in this note.

Nor does thermodynamic stability by itself select this readout. Take
𝓕_01=Γ_0Γ_1, 𝓕_23=±Γ_2Γ_3 and other components zero. The two oriented
J densities are ∓2, while the normalized Frobenius defect norms are both 2.
The action is indefinite as a local oriented functional; it is not a positive
entropy-production or return-norm minimization functional. The positive
thermodynamic Hessian used by MG supplies a regular adapter, not a proof of
this gravitational process law.

## LR-6. Native return square and the two cosmological sectors

Equation (13) fixes the cosmological coefficient relative to the soldering
scale. For positive u, set ℓ=1/(2u). Then

| Native realization | σ | Λ | β required if κ>0 is imposed |
|---|---|---|---|
| Original real four-component module | +1 | −3/ℓ² | β<0 |
| Extra factor K²=I | +1 | −3/ℓ² | β<0 |
| Extra factor R²=−I | −1 | +3/ℓ² | β>0 |

The square sign, not the sign of β, determines the sign of Λ in this class:
changing β rescales both the Palatini and volume coefficients together.
Admitting the independent native R factor therefore supplies a real positive-Λ
sector that the original module cannot provide by a Lorentz-equivariant linear
change of its four vector generators.

If one additionally requires both κ>0 and β>0, this class selects σ=−1.
Those are additional sign requirements; positivity of an equilibrium Hessian
does not establish them. The enlargement and readout still need a native
process selector. The numerical magnitude of u, β and any physical unit
identification have not been derived.

## LR-7. Flat combined transport can encode curved Einstein geometry

There are explicit local solutions for both signs. Let n be 3 for σ=+1 and
0 for σ=−1, so η_nn=σ. On a patch x^n>0 define

\[
e^a=\frac\ell{x^n}dx^a,\qquad
q=\frac{\ell^2}{(x^n)^2}\eta_{ab}dx^a dx^b,\qquad
\omega^{ab}=\frac\sigma\ell
(\delta^a_n e^b-\delta^b_n e^a).
\tag{15}
\]

Direct exterior differentiation gives

\[
T^a=0,\qquad R^{ab}=-\frac\sigma{\ell^2}e^a\wedge e^b,
\qquad \operatorname{Ric}(q)=-\frac{3\sigma}{\ell^2}q,
\qquad R(q)=-\frac{12\sigma}{\ell^2}.
\tag{16}
\]

For u=1/(2ℓ), the two even contributions in (7) cancel, so **𝓕=0**.
All contractible combined-connection returns are therefore trivial locally,
while the Levi-Civita metric curvature is nonzero. The n=3 metric is an
anti-de Sitter patch; n=0 gives a de Sitter patch. These are actual smooth
local solutions, checked both by coframe/connection differentiation and an
independent coordinate-metric Ricci calculation. No global topology or
geodesic-completeness conclusion is intended.

At 𝓕=0 the first variation of any quadratic action (11) vanishes. This supplies
an explicit solution-existence result for the declared theory. It does not
prove uniqueness, stability or general Cauchy well-posedness. Nor do the
Einstein equations require 𝓕=0: for a torsion-free Einstein solution, 𝓕 retains
its Weyl-curvature part. Exact nonzero-Weyl algebraic controls show this
distinction for both signs; those latter controls are pointwise tensors, not
additional claimed global solutions.

The interpretation is precise: the combined transport compares geometry to
its admitted de Sitter/anti-de Sitter model, rather than measuring only the
Levi-Civita connection. This standard Cartan distinction [2] provides a
concrete meaning of native return closure **within the chosen adapter**.
It does not turn every flat native return into physical gravitational
curvature or erase SH-1's representation and global-memory ambiguities.

## LR-8. Scale ambiguity, the zero-scale boundary, and reproducibility

For any constant d>0, transform

\[
e\longmapsto d e,\qquad u\longmapsto u/d,\qquad
\omega\longmapsto\omega,
\quad (\alpha,\beta,\zeta)\longmapsto(\alpha,\beta,\zeta).
\tag{17}
\]

The **entire** combined connection (6), its curvature, all its loop returns
and the quadratic readout are unchanged. However

\[
q\longmapsto d^2q,\qquad
\kappa\longmapsto d^2\kappa,\qquad
\Lambda\longmapsto\Lambda/d^2.
\tag{18}
\]

Thus even complete knowledge of these return records cannot separately fix
the physical length scale and dimensionful couplings without an independent
calibration or a law breaking this ambiguity. The relation κΛ=6/β is invariant
under (17), but β itself is an admitted action/readout normalization. Physical
action and length units remain to be supplied; the rational fixtures use
dimensionless coordinates and coefficients. No value of G, c, Λ, α_em or ℏ
is predicted by these control numbers.

For example the internal-unit choice σ=−1, u=1/4, β=8 gives κ=1 and Λ=3/4.
The dilation d=2 yields κ=4 and Λ=3/16 with **identical returns** and the same
product 3/4. This is an explicit nonselection witness, not an empirical fit.

The zero-scale boundary is also substantive. At u=0 with finite fixed
coefficients, 𝓐 contains no coframe, and (11) is only a combination of
curvature characteristic forms. It does not give zero-Λ Einstein gravity.
To keep κ fixed as ρ→0, β must scale as −1/(2κρ); after an explicit
bulk-equivalent subtraction of the divergent Euler term, the local bulk
action has the ordinary zero-Λ limit. This statement concerns local bulk
equations, not convergence of the unmodified global action or boundary
charges. Variable u or variable readout coefficients would require a new
derivation and are not covered by the constant-coefficient theorem.

`loop_curvature.py` adds exact standard-library/Fraction controls:

- the exhaustive 64-unknown Lorentz-vector intertwiner problem, with rank 62;
  all three explicit real representations and their bracket signs;
- the exhaustive 55-unknown symmetric-pairing problem for both signs, with
  Lorentz rank 52 and full Cartan rank 54;
- twelve directly differentiated connection jets and thirty-six three-weight
  density reductions, with all characteristic/exact terms retained;
- twelve based rectangles for affine connections, with exact noncommutative
  Taylor coefficients through order four and inverse-path controls; four
  independent oriented return-density coefficients;
- ninety-six coframe variations of the Cartan-square seed against the
  unchanged, Einstein-tensor-checked MG action;
- twelve independent even/odd splits of the compact matrix field equation,
  including four directly differentiated Bianchi identities;
- eight independent metric/coframe checks of (15), for both cosmological
  signs, plus two nonzero-Weyl algebraic controls;
- the full-group symmetry failure, the positive-norm/oriented-sign control,
  three identical-return scale ambiguities and the zero-scale boundary.

The rank calculations exhaust the unknown coefficients in their fixed finite
linear problems; they are not random sampling. The variable geometry fixtures
support the written analytic proofs and do not replace them. All earlier
GS/NP/CP/MG proof files and result groups remain unchanged and source-pinned.
Verification is read-only by default:

```bash
python papers/ugd-kahler-propagation/verify.py --check
```

The remaining selection problem is now specific: determine the module and its
square sign, the native readout and its symmetry reduction, the scale u and
the normalization β from an independently stated process principle. This
development supplies their algebraic consequences and an explicit continuum
route, without presupposing that the present primitive laws already select
those inputs.

## Primary references

1. S. W. MacDowell and F. Mansouri, “Unified Geometric Theory of Gravity and
   Supergravity,” *Physical Review Letters* **38**, 739 (1977),
   [doi:10.1103/PhysRevLett.38.739](https://doi.org/10.1103/PhysRevLett.38.739).
2. D. K. Wise, “MacDowell–Mansouri gravity and Cartan geometry,”
   [arXiv:gr-qc/0611154](https://arxiv.org/abs/gr-qc/0611154), submitted 2006,
   revised 2009. Primary exposition of the combined connection, model
   curvature and gravity construction.
3. D. K. Wise, “Symmetric Space Cartan Connections and Gravity in Three and
   Four Dimensions,” [arXiv:0904.1738](https://arxiv.org/abs/0904.1738) (2009).
   Primary treatment of invariant pairings and the Immirzi/Cartan interface.

Local dependencies: the unchanged MG proof, thermo adapter and exterior-algebra controls,
NP native matrices, CP thermo chart, the pre-entropy return record and the
SH holonomy/geometry distinction. Their exact file identities are recorded
in `SOURCE_PINS.json`.
