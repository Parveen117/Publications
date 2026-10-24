# SRT-1 — Graded axis-response transport through a singular chart

## 1. Native contract and what is new

The question is whether a cut can carry the symmetric response structure behind
the TVSP compass, including axis information that a diagonal reading omits.
We give a positive reconstruction theorem through a specified singular chart,
and state exactly which additional record it consumes.

Reuse NT-1's native response element

    h = (a+c)1/2 + (a-c)K/2 + b L,
    K^2=L^2=1, R^2=-1, L=RK, KL=-LK,

represented faithfully by the symmetric matrix H=[[a,b],[b,c]]. All coefficient
arithmetic in the certificate uses the existing native cut/rational engine.
The real symmetric sector is sufficient here. It is not an identification of
every physical constitutive tensor with a two-component response.

CRC-1 already proves record sufficiency by kernel preservation, NT-6 distinguishes
descent from compression, and NJC-1/NT-8 state the covariant response tower.
We reuse these; their general intertwining statements are not new results here.
The additional result is an explicit analytic divisibility criterion, the exact
jet depth required for singular response reconstruction, and a signed thermo/EM
adapter that retains the potential differential.

The parameter t below is a dimensionless chart parameter, not assumed to be
physical time, a quantum state index or the unique lambda of another tower.
Weights and a comparison frame are declared input. No spacetime curvature
singularity or unique gravitational transport is inferred from chart rank loss.

## 2. Axis information is not determined by the diagonal

Use J=L=[[0,1],[1,0]], P=(1+J)/2 and Q=1-P. J flips two normalized readings;
P extracts their symmetric diagonal component. J is invertible; P is not.
For any response H, the full record is

    H = PHP + PHQ + QHP + QHQ.                         (SRT.1)

Retaining PHP alone discards the other blocks. In particular

    H_+ = [[5/2,1],[1,3/2]],
    H_- = [[3/2,1],[1,5/2]]

are positive definite, have the same determinant 11/4, the same ratio
1-b^2/(ac)=11/15, and the same diagonal compression PHP=3P. Their labelled
axis responses differ. J transposes their axes: H_-=J^T H_+ J.

The two unitary transports I and J both satisfy PU=UP and UP=P. They agree
on every vector in the retained diagonal but act oppositely on Q. Therefore
even diagonal intertwining plus unitarity does not determine the full response
transport. Equation (SRT.1) identifies the omitted record explicitly.

For a surjective cut A:V->W, a symmetric form H descends as H=A^T C A iff
ker A is contained in ker H. Proof: necessity follows by inserting a null
vector of A. Conversely define C(Ax,Ay)=H(x,y); annihilation of ker A makes
this independent of both representatives and surjectivity gives uniqueness.
This is the bilinear specialization of the existing observer-kernel criterion.
For A=[1,0], descent requires H12=H22=0. Symmetry H12=H21 alone is insufficient.
In particular a positive definite H cannot descend losslessly through a cut
with a nonzero kernel. Keeping an augmented record changes that conclusion.

## 3. Singular chart theorem and the required tower depth

Let H(t) be a real analytic symmetric d by d matrix near t=0. Fix integers
w_i>=0 and the native diagonal scaling

    A(t)=diag(t^w_1,...,t^w_d),
    C(t)=A(t)^(-T) H(t) A(t)^(-1),  t != 0.           (SRT.2)

This is the response transformation under y=A(t)x at fixed t. The chart loses
rank at zero if any w_i>0. The theorem concerns a removable singularity of C,
not inversion of A(0).

**Theorem.** C has a real analytic extension through t=0 iff for every pair
i,j the Taylor coefficients of H_ij of orders below w_i+w_j vanish. Equivalently
H_ij is divisible by t^(w_i+w_j) in the analytic germ ring. The extension is
unique and symmetric, and its normalized jets satisfy

    C_ij^(k)(0)/k!
       = H_ij^(k+w_i+w_j)(0)/(k+w_i+w_j)!,  k>=0.    (SRT.3)

**Proof.** For t nonzero each entry is exactly H_ij(t)/t^(w_i+w_j).
An analytic extension exists precisely when this Taylor series contains no
negative powers. Removing the vanishing prefix gives a convergent series on
the same sufficiently small neighborhood. Coefficient comparison proves
(SRT.3). Any continuous extension must take the same limit, so it is unique.
Symmetry holds entrywise. This proves necessity as well as sufficiency.

To reconstruct only C(0), the required record consists of the coefficients
at orders w_i+w_j, together with the verified vanishing conditions and weights.
To reconstruct derivatives through order k, the corresponding entries need
coefficients through w_i+w_j+k. A uniform sufficient depth is 2 max(w_i)+k.
For unrestricted symmetric families this depth is necessary in the worst case:
changing an entry first at its required order changes C_ij(0) while leaving all
lower coefficients unchanged. Symmetry alone imposes no relation fixing it.

An entire analytic germ C requires the entire corresponding germ H, not just
these finite boundary jets. Finite measurements do not certify exact analyticity
or exact vanishing orders. No noise-stable recovery claim is made.

If C(0) is positive definite, continuity gives a positive definite C on a small
two-sided neighborhood. Then H=A^T C A is positive definite for t!=0, while H(0)
may be degenerate. Positivity of H away from zero alone does not force positive
definiteness at the reconstructed boundary: C=diag(1,t^2) is a counterexample.

## 4. Exact passage example: first and second response grades

Choose w=(0,1), so A=P_axis+t Q_axis with P_axis=(1+K)/2, Q_axis=(1-K)/2.
For arbitrary constants a,b,c let

    H(t) = [[a, b t], [b t, c t^2]].                 (SRT.4)

Then C(t)=[[a,b],[b,c]] for all nonzero t and has exactly that extension at zero.
The retained record is

    a=H11(0),  b=H12'(0),  c=H22''(0)/2.           (SRT.5)

The value H(0) alone sees only a. Its apparent loss of cross and second-axis
response does not prevent recovery when the specified first/second jets were
retained before discarding the family. No information is created at zero.

For a=c=2, b=1 the reconstructed determinant is 3 and C is positive definite.
Taking b=-1 gives the same collapsed H(0) and the same capacity-style ratio
Delta/(ac)=3/4, but a different signed cross-response. The first jet distinguishes
them. Taking c=3 instead of c=2 leaves every zero/first jet unchanged but changes
the reconstructed second-axis response; the second jet is necessary.

This is a verified sense in which symmetric information can pass a singular
chart: it passes in a graded record. If H(t)=I instead, then the same A gives
C22=t^-2 and there is no finite extension. The chosen scaling weights matter:
the same raw family with w=(0,0) has limit diag(a,0), not the limit in (SRT.5).

The incoming chart parameter and its orientation also matter. Replacing t by
-t reverses the raw first jet; the chart A must be transformed consistently.
Transport of a signed cross-response cannot be reconstructed from squared ratios.

## 5. Changing frames: retain the derivative term

On t!=0 put Gamma=A_dot A^-1. Differentiating (SRT.2) gives

    C_dot = A^-T H_dot A^-1 - Gamma^T C - C Gamma.    (SRT.6)

Proof: differentiate both inverse factors in A^-T H A^-1 and collect terms.
For (SRT.4), Gamma=diag(0,1/t), so all apparent 1/t terms cancel and C_dot=0.
Omitting the frame terms would infer a spurious changing response. For variable
analytic C the same calculation recovers its actual derivative. Equation
(SRT.6) is a transformation identity, not a selected evolution equation.

If U_t(x)=x^T H(t)x/2, then at fixed t it equals y^T C(t)y/2 with y=A(t)x.
The recovered potential is finite at zero whenever C extends. When t varies,
the total differential also contains the partial-t term. A frozen-t reciprocity
relation alone never supplies that evolution or the missing derivative term.

## 6. What symbol replacement can faithfully transport

For a C2 potential U(q1,q2), put p_i=partial_i U. The response is symmetric:
partial_2 p1=partial_1 p2. At fixed t, a constant invertible linear chart
q_new=Aq and p_new=A^-T p preserves p.dq and transforms the Hessian as in
(SRT.2). This is a genuine potential/reciprocity-preserving adapter. Arbitrary
nonlinear coordinate changes require their additional Hessian terms.

The four corner functions follow from the same differential by subtracting
selected conjugate products. For S subset {1,2},

    U_S=U-sum_(i in S) p_i q_i,
    dU_S=sum_(i not in S) p_i dq_i - sum_(i in S) q_i dp_i.

These are always identities on the original state space. Treating each U_S as
a potential in the exchanged independent variables additionally requires its
Legendre chart to be regular. A label at a corner is not proof of invertibility.

| Adapter | Coordinates q | Conjugates p | Exact differential |
|---|---|---|---|
| Thermodynamic energy | (S,V) | (T,-P) | dU=T dS-P dV |
| Reversible EM energy, scalar component | (D,B) | (E,H_em) | dw=E dD+H_em dB |
| Same signed TVSP compass convention | (D,H_em) | (E,-B) | dPsi=E dD-B dH_em |

The last row retains T->E, V->H_em, S->D, P->B, but requires the partial Legendre
transform Psi=w-B H_em, solving H_em=w_B for B. It is not mere unsigned renaming.
For the quadratic EM energy w=(a D^2+2b DB+c B^2)/2, c!=0, its response is

    Hess(Psi) = [[a-b^2/c, b/c], [b/c, -1/c]].        (SRT.7)

Both cross derivatives agree. For a=c=2,b=1 this gives
[[3/2,1/2],[1/2,-1/2]]. The negative entry is required by this mixed potential;
its Hessian need not be positive even when the energy Hessian is positive.
Do not confuse a partial Legendre change of potential with the positive-form
congruence of section 3. Units, conjugate pairing and signs are part of the adapter.

The general reciprocal response identities and singular extension theorem are
label-independent. Thermal capacities and mechanical responses are differently
typed constrained derivatives in NT-1; this packet does not replace them with
an unsigned list called Cp,Cv,Cs,Ct. It preserves their source definitions.

## 7. Electromagnetic field equations require a differential transport too

The thermodynamic Maxwell relations are cross-response identities. Physical
electromagnetism also has flux/source equations. In a chosen oriented chart write
the field strength F and excitation G as 2-forms, current j as a 3-form:

    dF=0,  dG=j,  G=chi(F).

For a linear constitutive map chi, a candidate map T2 on 2-forms and T3 on
3-forms transports these equations if

    d' T2 = T3 d,  chi' T2 = T2 chi,  j'=T3 j.      (SRT.8)

Proof: substitute F'=T2 F and G'=T2 G in each target equation. These conditions
give a sufficient structure-preserving transport; they do not guarantee faithful
recovery if a T has a kernel. Current conservation also needs the next-degree
chain relation and d^2=0. The finite complex fixture checks precisely this
algebraic part; it is not a discretization/continuum proof of physical Maxwell PDEs.

Symmetry of a two-by-two energy Hessian does not supply T2,T3, the spacetime
exterior derivative, a six-component constitutive tensor or boundary sources.
A symbol substitution alone cannot establish (SRT.8). In particular this paper
does not yet connect the two-dimensional singular response adapter to the full
electromagnetic 2-form bundle.

This distinction is established in the primary electromagnetic literature:
Y. N. Obukhov and F. W. Hehl, *On the boundary-value problems and the validity of
the Post constraint in modern electromagnetism*,
[arXiv:0707.1112](https://arxiv.org/abs/0707.1112), separates flux equations,
source equations and constitutive response. We use that physical dictionary;
the finite native algebraic implication in (SRT.8) is proved above.

## 8. Certification and the next obligation

The analytic proof in section 3 covers all dimensions and all stated analytic
germs. The native certificate checks exact polynomial instances, required jet
orders, signed cross-response recovery, axis/diagonal ambiguity, frame terms,
Legendre signs and differential-chain controls. Serialized native word witnesses
replay the selected identities; a corrupted witness must fail. No analytic
continuation theorem is replaced by a finite sampling certificate.

The result supplies a concrete record-completion rule: retain each response
entry at the sum of its two axis grades. These grades describe vanishing order,
not quantum energy levels. They need not equal the generations of the nonlinear
TVSP lift or the covariant NT-8 tower without a separate intertwining theorem.

The next obligation is to derive admissible weights and a faithful transport
from the actual local cut/matter dynamics, then test that transport against
(SRT.8) or the gravitational constraints. This packet does not choose those
dynamics, quantize the readings, or solve quantum gravity. It establishes when
the proposed information passage is mathematically possible and what it costs
in retained response data.
