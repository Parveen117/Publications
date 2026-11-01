# YM-50: native counted-record reference functional and the heat bridge

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

**Result:** on the already derived unit-quaternion EMK sector, a specified
finite counted-record readout of the canonical stationary-path functional
has a positive, normalized, invariant limit. An explicit finite polynomial
equation supplies its convergence tail. Its values equal the formerly
admitted SU(2) reference integral on coefficient polynomials and their
uniform completion. A symmetric native turn protocol also constructs the
same free heat action, with an explicit finite-content refinement error.
The recognition pairing, bounded observables and declared interacting
transfers are intertwined on that completion.

This closes a **specified compact reference/heat representation bridge**.
It does not identify the entire native path algebra with functions on a
compact group. It does not close RH's general UGD completion, select a
physical probability state, choose a clock/coupling trajectory, or identify
the chain with four-dimensional Yang--Mills. The protocol and sector are
declared; their invariance and representation properties are proved.

## Source audit and premise order

The source audit re-read RH-Framework's current T01 ledger, T01-A/B/C,
E1/E3, E4, E5 and contamination audit at fd104f46. T01-E5C/E6 remain open
there. YM-27 correctly recorded that earlier boundary; this chapter does
not rewrite it or claim to solve those more general UGD tasks.

Use the cut field and radial completion already earned by F00-E, the
associative EMK sector and native coefficient pairing of NCG-1/3, and
YM-F1's quaternion/fabric identification. Write

\[
e_ae_b=-\delta_{ab}+\epsilon_{abc}e_c,\qquad
e_a^\dagger=-e_a,\qquad
\mathcal Q_\Sigma=\{x=x_0+\sum_{a=1}^3x_ae_a:x^\dagger x=1\}.
\tag{1}
\]

The four-component coefficient energy sum x_a^2 is a derived scalar-part
identity. The unit sphere is its coordinate description, not a primitive
measure space. A_coeff denotes cut-complex polynomial readouts restricted
to (1); A_unif is their uniform completion. Begin with rational coefficients
and extend them in the already constructed cut field.

No Haar integral, infinite Hilbert space, spectral theorem, or classical
heat equation generates the constructions T1--T6. Finite tensor energies,
finite elimination, polynomial differentiation and cut-Cauchy limits are
used explicitly. The later comparison to existing SU(2) objects is marked
as a representation theorem. This is a native-source implementation of
invariant averaging and diffusion, not a priority claim for those general
mathematical ideas. Haar's 1933 paper, DOI 10.2307/1968346, is historical
lineage, not an imported proof premise.

## YM50-T1: an actual finite Phi_Sigma readout

In the quaternion sector choose the rational unit turns
q_a=(3+4e_a)/5, a=1,2,3. Use a twelve-letter protocol alphabet:
six distinct idle records and the six records q_a,q_a^dagger.
Distinct records remain labelled even when they have the same endpoint.

For a polynomial p define the finite turn average

\[
(Pp)(x)=\frac12p(x)+\frac1{12}
  \sum_{a=1}^3\{p(q_ax)+p(q_a^\dagger x)\}.
\tag{2}
\]

Thus P^k p(1) is exactly a sum over 12^k labelled word records divided
by their count. Product order is the order in (2); word order is never
silently commuted. Define

\[
\varphi_N(p)=\frac1N\sum_{k=0}^{N-1}(P^kp)(1),\qquad N\geq1.
\tag{3}
\]

This has an explicit realization in the native finite path algebra.
Take the finite atlas labelled by k=0,...,N-1 and an (N-1)-letter record.
At label (k,w), read only its first k letters and retain the remaining
labels. Let D_(p,N) be the diagonal stationary-path kernel with value
p(endpoint_k(w)). There are M_N=N12^(N-1) labels. RH T01-B's
stationary zero-residue extraction and T01-E3's normalized finite counts
give the exact identities

\[
\boxed{\varphi_N(p)=M_N^{-1}\Phi_\Sigma(D_{p,N}),\quad
\varphi_N(\bar p\,r)=M_N^{-1}
 \Phi_\Sigma(D_{p,N}^\dagger\star D_{r,N}).}
\tag{4}
\]

The padding replicates the kth endpoint exactly 12^(N-1-k) times, proving
(4). Diagonal composition also gives D_p D_r=D_(pr), dagger preservation,
phi_N(1)=1 and phi_N(|p|^2)>=0. If p is pointwise nonnegative, every
finite readout is nonnegative. Finite cut-square Cauchy--Schwarz gives
|phi_N(p)|<=||p||_infinity.

The normalization and record protocol are specified data, not measured
detector probabilities. Equation (4) is the bridge to the existing
Phi_Sigma; no new functional is silently renamed as that raw extraction.

## YM50-T2: finite coefficient energy and the stationary sector

For homogeneous degree d use the finite symmetric tensor coefficient
space. If p(x)=sum_(|m|=d) p_m x^m, its coefficient energy is

\[
\|p\|_d^2=\sum_{|m|=d}\frac{m_0!m_1!m_2!m_3!}{d!}|p_m|^2.
\tag{5}
\]

This is NCG's sum-of-squares tensor energy: the symmetric tensor entry
p_m/binomial(d;m) appears binomial(d;m) times. Quaternion multiplication
by a unit element is orthogonal for the four coefficient squares.
Its d-fold tensor action, and its symmetric restriction, preserve (5).
No invariant integral is used to define this finite energy.

Let U_a p(x)=p(q_a x), and let dagger here mean the adjoint for (5).
Then U_a^dagger=U_(q_a^dagger), P is self-adjoint and contractive, and

\[
\langle p,(I-P)p\rangle_d
 =\frac1{12}\sum_a\|(I-U_a)p\|_d^2.
\tag{6}
\]

Hence the fixed space consists exactly of polynomials fixed by all three
turns. Its determination is algebraic:

* zeta=(3+4 iota)/5 has infinite order. The integer polynomials
  T_0=2,T_1=X,T_(n+1)=XT_n-T_(n-1) satisfy
  T_n(zeta+zeta^(-1))=zeta^n+zeta^(-n). If zeta^n=1, their monic
  leading term at X=6/5 would give, after multiplication by 5^n,
  6^n=0 modulo 5, an impossibility.
* Left multiplication by e_a has square -I. The finite projectors
  (I +/- iota E_a)/2 decompose its cut-complexified coefficient
  space. On the d-fold tensor, U_a has factors zeta^(d-2j).
  Infinite order shows that its fixed vectors have zero infinitesimal
  charge. Thus a common fixed polynomial is annihilated by all three
  left quaternion derivations.
* At a unit x, the three e_a x are an orthonormal basis of its tangent
  directions. This follows directly from (1), or from
  sum_a(e_a x)(e_a x)^T=|x|^2 I-xx^T.
  The rational unit chart
  ((1-|v|^2),2v)/(1+|v|^2), with v replaced by t v, connects 1 to
  every unit point except -1. The derivative of a fixed polynomial
  along this chart is zero. A rational function over a characteristic-zero
  cut field with zero derivative is constant. The remaining point follows
  by applying any q_a to -1. Therefore the polynomial is constant on (1).

Homogeneity now identifies the fixed space: it is zero for odd d, and
span{v_d} for even d, where v_d(x)=(sum x_a^2)^(d/2).
For the latter statement scale a nonzero x by its positive native norm;
for odd d compare x and -x. Polynomial equality follows over the infinite
cut field. Let Pi_d be the orthogonal projection onto this explicitly
known space using (5); Pi_d=0 for odd d. Its entries are rational.

## YM50-T3: a constructive Cauchy tail without a spectral theorem

Finite elimination solves, uniquely on the orthogonal complement,

\[
u_d=(I-P+\Pi_d)^{-1}(I-\Pi_d)p,\qquad
(I-P)u_d=p-\Pi_dp .
\tag{7}
\]

To check invertibility, pair a kernel vector with itself in (6).
The terms from I-P and Pi_d are nonnegative; their simultaneous
vanishing forces the vector to be both fixed and orthogonal to the
fixed space. It is zero. Finite Gaussian elimination supplies the inverse.

The finite telescoping identity gives

\[
\frac1N\sum_{k=0}^{N-1}P^kp-\Pi_dp
 =\frac{u_d-P^Nu_d}{N}.
\]

Every term of P is evaluation after a unit turn. Therefore P contracts
the uniform cut bound, and

\[
\boxed{\left\|\frac1N\sum_{k=0}^{N-1}P^kp-\Pi_dp\right\|_\infty
\leq\frac{2B_p}{N},\qquad
B_p=\sum_m D_\Sigma((u_d)_m).}
\tag{8}
\]

For a sum of homogeneous polynomials, add their B_p. The inverse in (7)
and this bound are finite rational computations on every rational input.
There is no asserted degree-independent bound on B_p.

Since Pi_dp is constant on (1), (3) has a cut-Cauchy limit, denoted Phi_Q(p).
The constants obtained from different polynomial presentations agree,
because all finite evaluations of a zero function vanish. Positivity,
normalization, the pairing identity and ||Phi_Q||<=1 pass from (4).
For f in A_unif and a polynomial p,

\[
|\varphi_N(f)-\Phi_Q(f)|
\leq2\|f-p\|_\infty+2B_p/N.
\tag{9}
\]

This constructs the uniform extension with an explicit two-part tail.

## YM50-T4: invariance, uniqueness and exact moment values

Left and right multiplication by any unit quaternion preserve (5) and
the radial fixed vector v_d. Thus their actions preserve the coefficient
of Pi_d p, proving left and right invariance of Phi_Q. The same argument
works for every orthogonal change of the four native coefficients.

Uniqueness is finite: if a normalized linear functional nu on polynomial
readouts is invariant under all q_a and q_a^dagger, (7) gives
nu(p)=nu(Pi_dp)=Phi_Q(p). No measure uniqueness theorem is assumed.
Consequently another finite protocol with the same fixed spaces yields
the same limiting polynomial functional.

The moment values are

\[
\boxed{\Phi_Q(x^m)=0\ \text{if any }m_a\text{ is odd};\qquad
\Phi_Q\!\left(\prod_{a=0}^3x_a^{2n_a}\right)
=\frac{\prod_a(2n_a-1)!!}{4\cdot6\cdots(2n+2)},\
n=\sum_a n_a.}
\tag{10}
\]

The empty product is one. Odd moments vanish directly because Pi_d's
radial vector has only even exponents. For even moments, differentiating
invariance under coefficient-plane rotations relates the four raised
moments: M(n+e_a)/(2n_a+1) has the same value for every a.
The polynomial identity sum x_a^2=1 says their sum is M(n).
Hence M(n+e_a)=(2n_a+1)M(n)/(2n+4), proving (10) from M(0)=1.
These are consequences of counted native turns, not imported sphere
integrals. Independent finite coefficient projections check the formula.

## YM50-T5: the reference and recognition representation bridge

Only now introduce the existing SU(2) matrix chart of NCG-3/YM-F1.
It is a product/dagger-preserving identification of (1), with
half-fundamental character x_0. The previously admitted reference
functional Phi_ref is normalized and invariant under the six turns.
The uniqueness argument above therefore proves

\[
\boxed{\Phi_Q(p)=\Phi_{\rm ref}(p),\qquad p\in A_{\rm coeff},
\quad\text{and hence on }A_{\rm unif}.}
\tag{11}
\]

The reference measure is a comparison object in (11); it was not used
to generate Phi_Q. Independent finite products of the counted protocol
give Phi_(Q,I) and satisfy the same identity at every finite set of sites.
Independence of those reference copies is the declared product protocol,
not a conclusion that physical interacting states factorize.

The native recognition energy is Phi_Q(bar p p). Its finite cut-square
bound follows from (4); remove its null radical and take recognition-Cauchy
classes as in T01-B/N11. This defines the completed native coefficient
module. Equation (11) gives an onto isometry W_I to the old coefficient
completion, because both have the same dense coefficient core. Moreover

\[
W_I M_p^{\rm native}=M_p^{\rm ref}W_I,\qquad
\|M_p^{\rm native}\|\leq\|p\|_\infty.
\tag{12}
\]

The bound is first Phi_Q(|pf|^2)<=||p||_infinity^2 Phi_Q(|f|^2);
positivity makes it descend through the null space and extend.
The resulting Hilbert-space description is thus a proved representation
of this constructed positive completion.

This is not an isometry of the raw full record/path algebra: distinct
records with the same quaternion endpoint can have positive native path
energy but identical readouts. Their path memory has not been erased by
an algebraic identity. Also, convergence of phi_N(p) does not assert that
the endpoint functions on successive word atlases are recognition-Cauchy
in the pre-existing refinement module; see T8.

## YM50-T6: the free heat action from native symmetric turns

Define D_a p by left substitution exp_Sigma(-u e_a/2)x at u=0.
The exponential here is F00-E's convergent factorial construction,
applied to a native square-minus-one element. On degree d, each D_a
is skew for (5), with coefficient-energy bound ||D_a||<=d/2.
This follows by summing its d tensor-factor actions, each of bound 1/2.
Set L=-sum_a D_a^2; it is positive on this finite energy and
||L||<=3d^2/4.

Use the six-turn average

\[
S_\epsilon p=\frac16\sum_a
 \{p(\operatorname{Exp}_\Sigma(-\epsilon e_a/2)x)
       +p(\operatorname{Exp}_\Sigma(\epsilon e_a/2)x)\}.
\]

Each S_epsilon is a positive, unital uniform contraction and preserves
Phi_Q. Finite factorial series define E_d(t)=Exp_Sigma(-tL).
Differentiating its coefficient energy proves contraction:
d||E_d(t)v||_d^2/dt=-2 sum_a||D_a E_d(t)v||_d^2<=0.
The native series also proves composition and self-adjointness.

The symmetric fourth-order Taylor remainder and the second-order
remainder for E_d(t/n), with epsilon=sqrt(6t/n), give

\[
\|S_\epsilon-I+(t/n)L\|_d\leq3d^4t^2/(32n^2),\qquad
\|E_d(t/n)-I+(t/n)L\|_d\leq9d^4t^2/(32n^2).
\]

These estimates use the norm-one turn action and contractive E_d inside
their integral remainders; there is no omitted growing exponential.
Telescoping n contractions proves

\[
\boxed{\|S_{\sqrt{6t/n}}^{\,n}-E_d(t)\|_d
\leq 3d^4t^2/(8n).}
\tag{13}
\]

Evaluation at a unit x has bound one in (5), by its symmetric tensor
formula. Thus (13) also controls the uniform error on each homogeneous
polynomial. Sum the degree bounds for a general polynomial.
The limit is consistent under the unit-norm relation because every
finite turn respects it. Uniform approximation and contraction extend
the limit, positivity, invariance and composition to A_unif.

Direct quaternion multiplication gives the polynomial identity

\[
L=\tfrac14\{\mathsf E(\mathsf E+2)-|x|^2\Delta_{\rm coeff}\},
\quad \mathsf E=\sum x_a\partial_a.
\tag{14}
\]

The coefficient differential expressions on the right are derived
computational descriptions of the left-turn sum. A homogeneous harmonic
polynomial of degree ell has L-value ell(ell+2)/4=j(j+1), j=ell/2.
The decomposition into |x|^(2k) times harmonics follows by finite
elimination from
Delta_coeff(|x|^(2k)h_ell)=4k(ell+k+1)|x|^(2k-2)h_ell.
These are exactly the coefficient factors of the old free heat family.
The signed axis relabelling between NCG-3/YM-F1 and its Pauli chart
does not change the sum of three derivative squares.

Hence (12) intertwines this native turn limit with the existing K_t.
On the recognition completion the action extends contractively by
Phi_Q invariance and finite cut-square averaging; core approximation
proves strong continuity. The factor 6 and the generators e_a/2 fix the
heat-parameter convention. This constructs one specified isotropic
turn law; it does not derive the physical clock or prove that every
native process must have a second-order isotropic generator.

## YM50-T7: declared interactions and limiting readouts are preserved

YM-F1 already writes the bridge multiplier as the native identity
coefficient of a relative quaternion, and its weight as the factorial
exponential of the face-residue sum. Such bounded exponentials are
uniform limits of polynomials. On each finite chain, (12) and T6 imply

\[
W_I K_t^{\rm native}=K_t^{\rm ref}W_I,\qquad
W_I\operatorname{Exp}_\Sigma(aB)^{\rm native}
 =\operatorname{Exp}_\Sigma(aB)^{\rm ref}W_I,
\]
\[
\boxed{W_I S_a^{\rm native}=S_a^{\rm ref}W_I,\qquad
S_a=K_{a/2}\operatorname{Exp}_\Sigma(aB)K_{a/2}.}
\tag{15}
\]

Thus operator norms, positive vacuum selection, normalized transfers
and local A_unif readouts agree on the same finite chain. Norm limits
and vacuum vectors already constructed by YM-44/45 are transported
through W_I. There is no assertion that the different W_I themselves
have a common infinite-volume operator limit.

Instead, equality of finite local histories and the explicit YM-47
tails transport their scalar limits. Their reflected pairings agree;
the null quotient and completion therefore give an isometry of the
closed A_unif-generated history sectors. It intertwines literal time
translation and the time-zero observables of YM-48/49. No identification
with the larger arbitrary-bounded-history completion is added.

The chain, face law, independent reference copies and trajectory
kappa(a)=theta a remain declared choices; the transferred gap still
uses abs(theta)<1/1680. T1--T6 replace the reference/heat starting
objects by the specified native construction; they do not derive the
physical interacting state from unconstrained native primitives.
The NCG classical action's quantum measure, general UGD state-locked
atlas, actual row-closure defect and four-dimensional continuum remain open.

## YM50-T8: refusal witnesses and evidence scope

1. A finite protocol must explore the full compact sector. Retaining
   only the e_1 turn leaves x_2^2+x_3^2 invariant; starting at 1 its
   value stays zero, whereas (10) gives 1/2. A single-axis native
   average cannot certify the SU(2) reference functional.
2. Finite depth is not exact invariance: P x_0=(4/5)x_0. Its finite
   Cesaro average at 1 is nonzero, while Phi_Q(x_0)=0. The tail in
   (8) is required.
3. Scalar readout convergence is not raw-record vector convergence.
   On the nested full-word atlas put f_n(w)=x_0(endpoint_n(w)).
   Direct degree-two recursion gives
   P x_0=(4/5)x_0 and
   P x_0^2=1/4+(43/75)(x_0^2-1/4).
   Consequently the squared difference of f_(n+1) and the refined f_n is
   1/10-(1/50)(43/75)^n, bounded below by 2/25, with limit 1/10.
   Their recognition presentations are
   not Cauchy. The new state pairing is constructed from limits of
   readouts, rather than falsely invoking E4 on these endpoint vectors.
4. Different full paths can have the same endpoint. The empty word and
   q_a q_a^dagger share endpoint 1 but are distinct retained records.
   The difference of their unit path coefficients has native energy 2.
   An endpoint observer is not a faithful full-memory quotient.
5. Wrong normalization makes the constant readout grow like the record
   count; dropping padding weights changes (3). Averaging only final
   eigenvalues does not prove the product/dagger readout identities.
6. The heat scaling is fixed: a turn of size sqrt(t/n) instead of
   sqrt(6t/n) gives generator L/6. A single-axis protocol also gives
   a different generator. Symmetric-turn convergence is finite-content;
   (13) has d^4 and is not a content-uniform operator-norm estimate.
7. This compact bridge is not the entire RH E5C/E6 claim, nor proof of
   physical state/clock/trajectory selection, countable path-measure
   extension, NCG quantum-field identification, AF, Clay or QG.

The written general arguments are not mechanically formalized.
Exact controls compare raw finite word counts with polynomial transfer,
check coefficient energies and fixed spaces, solve (7), verify (8),
derive moments independently by finite projection, test left/right
invariance and the native recognition pairing, and check the turn
generator and heat normalization. The origin ledger separates native
inputs, constructed objects, comparison identities and open selections.
No general RKF/RH engine or old evidence is modified.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym50_native_reference.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym50_native_reference.py -v
~~~

Default/--check is read-only. Only --write regenerates the new evidence.
CI uses one Python 3.12 job.
