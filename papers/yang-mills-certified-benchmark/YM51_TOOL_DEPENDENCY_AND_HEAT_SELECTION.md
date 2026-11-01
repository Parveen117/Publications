# YM-51: tool dependency audit and the native heat-selection boundary

Monty Dabas. 3 October 2026. Verification: **Python 3.12 only**.

## Audit result

The existing native curvature calculus does not require replacing every
mathematical tool. It requires preserving the carrier, composition,
observer, pairing and completion assumptions of each tool. This chapter
audits the tools on the current YM-43--50 route, not every theorem in the
framework or every branch of physics.

The decisive remaining choice is the symmetric second moment of the
native turn protocol. We construct different positive heat laws on the
same quaternion algebra, positive reference and derivative frame. They
have the same action on every linear coefficient, but different quadratic
decay. A full-support family has decay tending to zero although the
native order brackets remain fixed. Curvature and the reference pairing
therefore do not select the heat law or its gap.

An explicit additional invariance condition forces isotropy up to one
positive scale. Ordinary covariance, which transforms the protocol data
along with the frame, does not force isotropy. This is a selection
classification and a non-selection result, not a new physical clock or
an extension of the existing interacting gap window.

## A. Existing tools and their exact status

The machine-readable companion is
[YM51_DEPENDENCY_LEDGER.json](certificates/YM51_DEPENDENCY_LEDGER.json).
Its dependency checks are provenance/scope checks, not a proof assistant.

| Tool | Existing derivation or representation | Boundary relevant to this route |
| --- | --- | --- |
| Scalars and lawful products | F00-E; RKF N01; EMK associative sector and path dagger | Not every partial UGD numeral product is an associative scalar field |
| Tensor, dual and contraction | Extra R7; central-field typed slots, dual evaluation and induced Lie action | Cross-carrier contraction needs a bridge; a metric is not an endomorphism; no arbitrary typed/noncommutative tensor product is supplied |
| Derivatives and connection curvature | EMK-C1: inner derivations, Leibniz, frame subtraction, gauge covariance and Bianchi | Connection coefficients and the physical direction frame remain choices |
| Observation and memory | RKF N03/N06/N09; Extra R8--R13 | Compression is not generally multiplicative; descent and target faithfulness are separate; moving cuts need derivative terms |
| Spectral and curvature readouts | Extra R11; marked recovery on its specified family | A spectrum, trace or zero mean is not a full-curvature or future-completeness certificate |
| Metric and adjoint | NT positive response pairing; Extra R13/R36 protocol-derived information/response metrics | These are particular pairings; identifying one with a physical spacetime metric needs a separate map |
| Riemann specialization | Extra R10; NT-4's explicit Hessian connection intertwiner | Tangent identification, metric compatibility and torsion-free connection laws are required |
| Reference functional and Hilbert representation | YM50 counted Phi_Sigma readout, Phi_Q, null quotient and completion | Specified compact sector; not the full raw-path/refinement module or general RH E5C/E6 |
| Free heat action | YM50 symmetric native turns and finite-content convergence | Protocol and heat parameter were specified; this chapter classifies their second-moment freedom |
| Interacting time action and observables | YM43--49; fixed law, trajectory and abs(theta)<1/1680 | The old proofs apply to that family, not automatically to the alternative generators below |
| Clocks and physical interpretation | Morphic clock calculus; Extra R39--R41 constructed relational clocks | A clock must resolve the process; a native tick does not select physical seconds or the YM clock coupling |

The historical sources already contain corrections, rather than an
instruction to discard their mathematics:

* EMK-C1 distinguishes raw [delta_i,delta_j] from connection curvature
  [nabla_i,nabla_j]-nabla_[delta_i,delta_j]. Its noncommuting inner frame
  with A=0 has zero connection curvature.
* R8 proves [PAP,PBP]=P[A,B]P-PAQBP+PBQAP. R9 retains the odd/odd
  commutator in a balanced even readout. Neither supports universal
  observation-induced flatness.
* R7 proves that the continuous two-mode R and K flows have no common
  nonzero fixed symmetric metric. Finite involutions, flows, alternating
  forms and moving pairings are distinct structures.
* NT-3 gives F^(qX)=q(q-1)[X_i,X_j] for X=H^-1 delta H. The strain
  choice q=1/2 and the flat pure-gauge choice q=1 are different
  connections on the same response data.
* R12's supplied Gaussian experiment is not the last metric result.
  R13 constructs binary record/noise/event ingredients for specified
  protocols; R36 derives a phase response metric from its actual source.
  Neither is a universal spacetime-metric selector.
* Morphic Calculus's certification ledger retains open/overbroad claims;
  its chosen-clock derivative is not used as an unconditional theorem.
  R41 constructs relational evolution for a declared process V; it
  does not select the YM heat process.

These source identities are reused unchanged. Familiar tensor, connection,
quadratic-form and diffusion ideas retain their established mathematical
lineage. The contribution here is the audited dependency map, an explicit
native protocol comparison, its quantitative refinement bound and the
resulting selection/observation tests.

## B. Contract for the comparison

Use exactly YM50's unit-quaternion EMK carrier Q_Sigma, coefficient
algebra A_coeff, uniform completion A_unif, positive reference Phi_Q
and finite homogeneous coefficient energy ||.||_d. Put

\[
D_a p(x)=\left.\frac{d}{du}
 p(\operatorname{Exp}_\Sigma(-u e_a/2)x)\right|_{u=0},
\qquad D_v=\sum_{a=1}^3v_aD_a.
\tag{1}
\]

Quaternion multiplication derives [D_1,D_2]=D_3 and its cyclic variants.
Each D_v is skew on the finite coefficient energy and on the Phi_Q
polynomial pairing. The former follows from the native coefficient
norm; the latter follows by differentiating Phi_Q's turn invariance.
The coefficient bound is ||D_v||_d<=d|v|/2.

Here |v|^2=sum v_a^2 is the scalar-part norm of the imaginary quaternion,
not a physical spatial metric. The directions a label an internal
quaternion frame. Holding this entire frame and any independently
specified connection fixed holds its curvature fixed in all comparisons.
No free-turn commutator is relabelled as NCG gauge or Riemann curvature.

## YM51-T1: the declared counted protocol produces a symmetric tensor

Choose finitely many native imaginary coefficient vectors v_r and
nonnegative rational weights w_r with sum w_r=1. Rational weights are
realized by repeated labelled records, with equal counts for the two
opposite turns. Zero vectors are idle records. Define

\[
S_\epsilon^{\mathcal P}p(x)=
\sum_r\frac{w_r}{2}
\{p(\operatorname{Exp}_\Sigma(-\epsilon v_r/2)x)
 +p(\operatorname{Exp}_\Sigma(\epsilon v_r/2)x)\}.
\tag{2}
\]

This is an actual finite normalized diagonal-record readout of the
type in YM50-T1. It preserves constants, order positivity and Phi_Q;
it contracts both the uniform norm and the Phi_Q norm. Dagger and
opposite-turn symmetry make its action self-adjoint for the latter.

The protocol selects

\[
C_{\mathcal P}=\sum_rw_r v_rv_r^T\succeq0,\quad
m_2=\sum_rw_r|v_r|^2=\operatorname{tr}C,\quad
m_4=\sum_rw_r|v_r|^4,
\]
\[
\boxed{L_C=-\sum_rw_rD_{v_r}^2
       =-\sum_{a,b}C_{ab}D_aD_b.}
\tag{3}
\]

For polynomial f,
Phi_Q(bar f L_C f)=sum_r w_r Phi_Q(|D_(v_r)f|^2)>=0.
This derives the energy form from the counted protocol. It does not
deduce its weights from curvature. We compare finite protocols of this
form, not every possible stochastic, nonlocal or memory-bearing law.

## YM51-T2: a native refinement bound for the whole protocol class

On degree d define E_C(t)=Exp_Sigma(-tL_C) by the existing finite
factorial series. Energy differentiation proves contraction, and the
series proves composition. For t>=0 and n>=1,

\[
\boxed{\left\|(S_{\sqrt{2t/n}}^{\mathcal P})^n-E_C(t)\right\|_d
\leq\frac{d^4t^2}{96n}(m_4+3m_2^2).}
\tag{4}
\]

**Proof.** Symmetric Taylor expansion eliminates odd terms. Its
fourth-order remainder is bounded by
epsilon^4 sum_r w_r ||D_(v_r)||_d^4/24, because the finite turn action
inside the remainder is norm one. Substitution of epsilon^2=2t/n
gives d^4 m_4 t^2/(96n^2).
Also ||L_C||_d<=d^2 m_2/4. The second-order remainder of the
contractive E_C(t/n) is at most d^4 m_2^2 t^2/(32n^2).
Subtract the common I-(t/n)L_C term and telescope n contractions.
The sum is (4). These one-parameter remainders can be obtained by
integrating the convergent finite-matrix factorial series; no compact
group integral or spectral theorem supplies the construction.

Unit-point evaluation has norm one in YM50's tensor energy. Therefore
(4), multiplied by ||p||_d and summed across degrees when necessary,
controls polynomial uniform error. Every finite turn respects sum x_i^2=1;
the limiting action is consistent on A_coeff. Uniform approximation and
contraction extend it to A_unif. The positive reference pairing then
extends it contractively to its recognition completion, where core
approximation gives strong continuity. Different polynomial protocols
with the same C have exactly the same limiting action on these carriers.

YM50 is recovered using v_a=sqrt(3)e_a, w_a=1/3. Then C=I, m_2=3,
m_4=9, and (4) is precisely 3d^4t^2/(8n), with the original turn size
sqrt(6t/n). The square root is in the already completed radial field.
The controls also use entirely rational finite protocols.

Equal C does not identify finite steps or full retained histories.
For example one unit-axis turn and a protocol with weight 1/4 on twice
that axis plus 3/4 idle have the same C. On x_0 their one-step readouts
are cos_Sigma(epsilon/2) and 3/4+cos_Sigma(epsilon)/4.
Their difference, second minus first, starts at epsilon^4/128.
The continuum equality has not erased the finite record distinction.

## YM51-T3: curvature, reference and linear probes do not select heat

For 0<eta<=1 put

\[
C_\eta=\operatorname{diag}(3-2\eta,\eta,\eta).
\tag{5}
\]

For rational eta this is realized by v_a=3e_a with weights
(3-2eta)/9, eta/9, eta/9 and an idle weight 2/3. Every axis occurs
with positive weight. All protocols have the same quaternion algebra,
derivative brackets, Phi_Q and trace m_2=3.

On every linear coefficient x_mu,

\[
\boxed{L_Cx_\mu=\tfrac14\operatorname{tr}(C)x_\mu,\qquad
E_{C_\eta}(t)x_\mu=e^{-3t/4}x_\mu.}
\tag{6}
\]

**Proof.** On linear coordinates the anticommutator of the native
generators is D_aD_b+D_bD_a=-delta_ab I/2. Contract with the
symmetric C. Its off-diagonal terms cancel and its diagonal terms
give (6).

Nevertheless let

\[
f=x_2^2+x_3^2-\tfrac12\sum_{\mu=0}^3x_\mu^2.
\]

Direct native differentiation and YM50's moments give

\[
D_1f=0,\quad D_2^2f=D_3^2f=-f,\quad
\Phi_Q(f)=0,\quad \Phi_Q(f^2)=1/12,
\]
\[
\boxed{L_{C_\eta}f=2\eta f,\qquad
E_{C_\eta}(t)f=e^{-2\eta t}f.}
\tag{7}
\]

Thus even exact observation of all linear-coefficient decay cannot
distinguish this family. Its centered all-source decay rate cannot
exceed 2eta. Taking eta down to zero disproves a uniform positive
rate across these protocols despite full axis support at every eta>0.
At eta=0 the nonzero centered f is stationary.
This does not contradict YM45's fixed isotropic small-bridge gap.
It changes the free protocol; that theorem is not asserted uniformly
over (5). No interacting-chain gap for these alternatives is claimed.

## YM51-T4: the exact extra symmetry assumption

For a native unit quaternion g let R_g describe conjugation of its
three imaginary coefficients, and U_g p(x)=p(gx). Native multiplication
gives orthogonality of R_g and

\[
U_g^{-1}D_vU_g=D_{R_gv},\qquad
\boxed{U_g^{-1}L_CU_g=L_{R_gCR_g^T}.}
\tag{8}
\]

Equation (8) is covariance. It holds for anisotropic C as well as
isotropic C, with the protocol tensor transformed.

For comparison, impose the additional active invariance condition
U_g^-1 L_C U_g=L_C with the same C. For this constant symmetric
second-order class it is equivalent to R_g C R_g^T=C. Indeed,

\[
\boxed{C_{ab}=-2\,(L_C(x_ax_b))(1),\qquad a,b=1,2,3,}
\tag{9}
\]

so C->L_C is injective. At x=1 the imaginary coordinates vanish and
D_i x_a=-delta_ia/2; the two product-rule terms prove (9).

Only four native conjugations are needed to enforce full isotropy:
g=e_1,e_2,e_3 and g=(1+e_1+e_2+e_3)/2. The first three flip two
imaginary axes and force every off-diagonal C entry to vanish. The
last cyclically permutes the axes and forces the three diagonal
entries to agree. Consequently

\[
\boxed{\text{fixed-law invariance under these conjugations}
\iff C=cI,\quad c\geq0.}
\tag{10}
\]

Conversely cI is invariant under every native unit conjugation by
orthogonality. The invariance requirement is a stated physical/model
selection hypothesis, not a consequence of covariance or a redefinition
of curvature. Positivity alone also allows c=0.
These are constant internal conjugations; neither local gauge Ward laws
nor physical spatial isotropy are derived by this argument.

## YM51-T5: a lawful observer of the missing tensor

Equation (9) gives six fixed real scalar channels, one for each
a<=b, that recover the entire symmetric C. Six are minimal among
fixed linear readouts on the unrestricted six-dimensional symmetric
coefficient space: fewer rows have a nonzero kernel. The same lower
bound holds on positive definite C, since a sufficiently small
perturbation of I in any such kernel remains positive definite.
This specializes the existing RKF/R11 target-recovery argument.
It is not a universal minimum for protocols with extra prior constraints.

For diagonal C_eta, the response in (7) already distinguishes eta.
The one-parameter witness and the unrestricted six-channel bank have
different reconstruction scopes. They recover the limiting second
moment, not every microscopic word probability, history or physical
state. Equation (4) supplies the refinement error before treating
finite observations as continuum readouts.

## YM51-T6: isotropy still leaves the clock and interaction choices

If the invariance hypothesis gives C=cI, then

\[
L_C=cL_I,\qquad E_C(t)=E_I(ct).
\tag{11}
\]

All positive c preserve the same reference and native frame brackets.
Fixing tr(C)=3 or linear-coefficient decay 3/4 sets c=1 as a heat
parameter convention. It does not derive physical seconds, energy,
h-bar or a coupling trajectory.

The relative interacting strength also changes when the clock changes.
For the already declared finite-chain expression
S_a=K^C_(a/2) Exp_Sigma(a theta V) K^C_(a/2), write b=ca.
Then it is exactly
K^I_(b/2) Exp_Sigma(b (theta/c) V) K^I_(b/2).
Thus theta and a cannot both be held fixed while calling c a harmless
relabeling of the same interacting law. This is an algebraic parameter
identity; all step/window hypotheses of any imported gap estimate must
still be checked. No new gap window is inferred here.

## C. Consequence for the programme

The audit leaves the proven algebra, tensor rules, curvature/observer
corrections and scoped completions in place. It rules out three shortcuts:
using curvature alone to choose a Laplacian; treating covariance as
fixed-law isotropy; and using agreement of linear decay to certify the
full evolution or its gap.

For the current YM family the mathematical next operator problem remains
YM49's actual infinite-chain row-closure defect. For a physical selection
claim, the next separate gate is to justify fixed-law internal isotropy
or identify and retain a non-isotropic C, then connect the remaining
clock and interaction ratio to a declared native process/calibration.
Existing relational-clock and response constructions are candidate
inputs; they have not been proved to select this process.

General typed-module extensions, moving observers, RH E5C/E6, NCG quantum
measure, physical spacetime/4D Yang--Mills, AF, Clay and QG remain open at
their previously stated scopes. The word curvature has not changed
these target requirements. No novelty claim for general invariant
quadratic forms or diffusion is made.

The written general proof is not mechanically formalized or externally
expert-certified. Exact controls check finite coefficient actions,
native conjugation covariance, the invariant-tensor nullspace,
quadratic recovery, the full-support slow-mode family, independent
turn counts and outward refinement enclosures. Ledger mutations reject
missing parents, unmarked choices and promotion of physical selections.
They check the recorded contract, not all mathematical truth.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym51_heat_selection.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym51_heat_selection.py -v
~~~

Only --write regenerates this new chapter's evidence. Upstream research,
canonical engines and the Extra Ideas R1--R46 register remain unchanged.
