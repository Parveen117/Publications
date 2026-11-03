# Fine-structure selection: native return sensitivity and the dressed source

Monty Dabas · AS-1–AS-5 · 30 September 2026

**Numerical alpha: NOT DERIVED.** This packet tests the proposed route from
native cut selection to an electromagnetic coupling. It finds a precise
remaining freedom even in an exactly closed native aperture, and derives the
source-dependent coefficient that a physical electromagnetic adapter would
need to select. No experimental constant is an input to the verifier.

## Sources and contract

Consume canonical RKF commit `3cc5a33b05c16d59c90994ddda69dedc0d392424`:
its actual R2 paired-depth return solver, Pāṇinian engine, cut arithmetic,
and WC5 source-preserving elimination theorem. The source API is imported,
not copied. The [existing alpha reduction](../thermo-compass-foundations/downstream/FINE_STRUCTURE_R1.md)
supplies the conditional target q_e²/(4πZγ), with its physical normalization
and low-momentum matching requirements. The [CT trajectory](../native-cut-energy-transport/THEOREM.md)
and thermo phase connection are used only within their stated scope.

Two operator layers must remain distinct. The R2 excursion return is a
native inverse corner, not automatically a positive energy Hessian or a
photon propagator. The later positive sourced quadratic model is an explicitly
admitted adapter. No theorem equates its parameters with R2 cell products.
R2 uses L=KR; CF/CT uses RK=-L. We keep R2's convention in AS-1/AS-2.

## AS-1. Exact cut closure leaves a continuum of response sensitivities

In the R2 paired-depth family, let the two positive cell products be a,b.
The completed boundary response is F=1+xL, with

\[
bx^2+x-a=0,\qquad x>0.
\]

R2 proves aperture convergence for every positive periodic profile. Its
normalized response Q=F/2 is an exact dagger-symmetric idempotent precisely
when x=1, or a=b+1. Thus every b>0 yields the same native cut Q=(1+L)/2.

Probe the supplied cells by a common positive pair multiplier z: (a,b)→
(az,bz). This z is a declared dimensionless perturbation of edge products,
not photon momentum, time, charge or an already fixed laboratory input.
The unique positive completed response satisfies

\[
bz x(z)^2+x(z)-az=0. \tag{1}
\]

For z>0 its derivative with respect to x is 1+2bzx>0, so the response is
smooth. At the exact-cut background a=b+1, z=1, differentiation gives

\[
\boxed{x(1)=1,\quad x'(1)=\frac1{1+2b},\quad
x''(1)=-\frac{2b(3+4b)}{(1+2b)^3}.} \tag{2}
\]

Indeed (1+2bzx)x'=a-bx², and
(1+2bzx)x''=-4bxx'-2bz(x')². These prove (2). The slope ranges over all
of (0,1) as the positive b varies, even though Q stays exactly unchanged.
Two fully native exact examples are:

| Cell period | Cut at z=1 | x'(1) | x''(1) |
|---|---|---:|---:|
| (2,1) | (1+KR)/2 | 1/3 | -14/27 |
| (3,2) | (1+KR)/2 | 1/5 | -44/125 |

Hence a rule that selects only exact aperture idempotence has not selected
return stiffness. The second boundary phase y=b/(1+b) distinguishes these
operators; x'(1)=(1-y)/(1+y). This is missing target information, not a
contradiction in native cut closure or a claim that the complete operators
are observationally identical.

Conversely a supplied slope s∈(0,1) reconstructs b=(1-s)/(2s), a=b+1.
Using an observed electromagnetic coupling to choose s would therefore be
an inverse design or fit unless a separate theorem independently identifies
and selects this probe and its slope. Renaming s as alpha is not that theorem.

## AS-2. Full response jets and finite-aperture error are retained

Set t=z-1 and x(1+t)=sum_(n≥0)c_n t^n, c0=1. For n≥1, native scalar
coefficient comparison in (1) gives

\[
(1+2b)c_n=(b+1)\mathbf1_{n=1}
-b\left[\sum_{i=1}^{n-1}c_i c_{n-i}
       +\sum_{i=0}^{n-1}c_i c_{n-1-i}\right]. \tag{3}
\]

This fixes every jet after b is supplied. It does not select b. The code
checks exact residual cancellation through order six, and c1 and 2c2 against
(2), using the canonical native cut-rational arithmetic.

Finite differences of a cutoff response do not automatically equal the
completed derivative. We avoid that interchange. From (1),

\[
x'=\frac{x}{z(1+2bzx)}>0,\qquad x''<0\quad(z>0).
\]

For 0<δ<1, concavity gives

\[
\frac{x(1+\delta)-1}{\delta}\le x'(1)
\le\frac{1-x(1-\delta)}{\delta}. \tag{4}
\]

Call R2's unchanged finite-aperture enclosure at both positive arguments,
retaining its nonnegative-tail contract and refinement error. If L+ and L−
are the respective certified lower endpoints, then
(L+−1)/δ ≤ x'(1) ≤ (1−L−)/δ. Thus the infinite response slope is bracketed
by actual aperture computations without differentiating a truncated series
or assuming an absent tail. For any fixed δ the source aperture error
converges to zero; the secant error vanishes subsequently as δ→0.

The code replays the native paired Schur identities and uses the original
solver for these enclosures. A finite execution supplements the general
convergence and concavity proofs; it is not the proof of every infinite case.

## AS-3. The coupling target must dress its source as well as its stiffness

After an explicit positive self-dagger quadratic representation and removal
of relevant null directions, declare a visible channel v, hidden channels h,
a block form H=[[A,B],[B†,C]], C>0, and a source column q=(qv,qh).
For a scalar source amplitude j, set

\[
E=\tfrac12(v,h)^\dagger H(v,h)-\operatorname{Re}\{j^\dagger q^\dagger(v,h)\}.
\]

WC5's native ordered elimination, applied to both operator and source, gives

\[
Z=A-BC^{-1}B^\dagger,\quad q_{\rm eff}=q_v-BC^{-1}q_h,
\quad h=C^{-1}(q_h j-B^\dagger v). \tag{5}
\]

Completing the square yields

\[
E_{\rm eff}(v;j)=\tfrac12 v^\dagger Zv
-\operatorname{Re}(j^\dagger q_{\rm eff}^\dagger v)
-\tfrac12|j|^2q_h^\dagger C^{-1}q_h.
\]

When Z is invertible, the full source response is

\[
\boxed{q^\dagger H^{-1}q
=q_{\rm eff}^\dagger Z^{-1}q_{\rm eff}
+q_h^\dagger C^{-1}q_h.} \tag{6}
\]

The first term carries the retained response; the second is the hidden-source
contact contribution. Omitting qh while retaining the Schur correction to Z
can change even the pole residue. Both terms are basis-covariant when the
source and quadratic form are transported together. In particular h=T h'
gives B'=BT, C'=T†CT, qh'=T†qh; a visible rescaling v'=t v transports
Z'=Z/t² and qeff'=qeff/t for real nonzero t, preserving qeff²/Z.

Consider the exact scalar illustration, for k>0,

\[
H(k)=\begin{pmatrix}\zeta k^2+B^2/C&B\\B&C\end{pmatrix},
\quad C>0,\ \zeta>0. \tag{7}
\]

Its positive Schur complement is Z(k)=ζk² and its source response is

\[
q^T H(k)^{-1}q=\frac{(q_v-Bq_h/C)^2}{\zeta k^2}+\frac{q_h^2}{C}. \tag{8}
\]

For B=1,C=2,qv=2,qh=1, the residue is 9/(4ζ) and the contact term is 1/2.
Ignoring the hidden source gives 4/ζ instead. The verifier checks (8) against
full native matrix inversion at nonzero k, not an inverse of the singular
k=0 matrix. The parameter k in this illustration is a supplied spectral
coordinate; identifying it with physical momentum is a further adapter.

If an electromagnetic matching theorem supplies one transverse massless
photon channel, the normalized electron source, a dimensionless action and
its low-momentum limit, the candidate alpha is that channel's pole residue
divided by 4π in the convention of the earlier note. Equation (8) alone
supplies neither a photon nor an electron, and its regular contact term is
not to be included as the massless pole residue. Charged quantum effects and
scale matching are not performed by this finite quadratic calculation.

## AS-4. The current geometrical and process laws do not remove this freedom

There are two separate non-selection statements, not a universal impossibility
theorem for RKF.

1. R2 exact cut selection leaves b free by AS-1. Even a full chosen boundary
   response does not reconstruct individual opening and closing amplitudes.
2. In the admitted quadratic field adapter, multiplying the photon quadratic
   cost by a positive factor λ while keeping the normalized matter/source
   charge fixed leaves the homogeneous field equations unchanged, but changes
   qeff²/Zγ by 1/λ. This is a different constitutive model, not a field-coordinate
   change: a coordinate change transports charge too and preserves the ratio.

The newer CT clock does not close this gap. Under W→λW, H→λH and
p→λp, its state law -H^-1p, response attenuation ratios and internal clock
are unchanged. Likewise J_H=RH/sqrt(detH) and the thermo connection
H^-1 dH/2 are unchanged for constant λ>0. Hessian density and numerical
energies do change; those are not claimed to be identical data. A calibrated
energy datum can constrain this scale, but none has been identified with the
physical photon action in the present construction.

In the phase-field model, scaling electric and magnetic quadratic response
matrices by the same λ leaves its source-free propagation equation unchanged.
A measured or derived common propagation speed alone therefore does not fix
the source coupling. An independent normalization law must remove these
freedoms relative to the matter/action standard; choosing λ=1 is a convention
or model choice until that relative identification is proved.

## AS-5. A concrete acceptance test for a future alpha prediction

A candidate derivation on this route must provide, before matching the observed
constant: a selected native interaction profile and vacuum; the charged source
identified with the electron; the photon channel and dimensionless action
normalization; source-preserving elimination with controlled tails; and the
low-momentum matching map. AS-1's b must be selected if that return family is
used, rather than solved backwards from a target number.

The immediate mathematical progress is the exact sensitivity map (2)–(4) and
the source-complete residue (5)–(8). The immediate unresolved law is one that
fixes interaction strength relative to the charged source and action standard.
Neither aperture idempotence nor the new internal clock currently does so.

**Result of this attempt: no numerical alpha prediction.** No proximity search
among numerical expressions, no measured alpha input and no adjustment to its
value is performed. This scoped non-selection result identifies a missing law;
it does not say that no richer native law can supply it. Computational PASS,
written proofs, physical matching and independent validation remain distinct.
