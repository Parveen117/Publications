# YM-56: a Hessian lambda family, preserved symmetry and an exact compact gap

Monty Dabas. 4 October 2026. Runtime: **Python 3.12 only**.

This chapter develops the symmetry/curvature question using the existing
NT and YM50--55 engines. A cubic potential realizes the response directions
K and lambda L after removal of scalar strain. Its selected native curvature
is lambda R/2 at the reference state. The prescribed equal-count compact
protocol has exact centered free rate

\[
 \boxed{\gamma(\lambda)=\min\{(1+\lambda^2)/8,
                                  \min(1,\lambda^2)/2\}.}
\]

An exact involution symmetry survives for every lambda. At lambda=0,
additional stationary modes appear; for nonzero lambda the constants are
the only stationary sector. Thus this example isolates loss of control of
modes, rather than universal symmetry breaking, as the gap-closing event.

Lambda here is a declared dimensionless coefficient of a constitutive
family. It is not identified with the unit-typed thermodynamic lambda
responses, YM55's interaction ratio, or a physical particle mass.
All response jets are frozen at the reference state before constructing
the heat protocol. No state-dependent diffusion is asserted.

## Inputs and attribution

Use NT-1--4, YM52-T4, YM54-T1--3 and YM55 with their existing contracts.
In the NT cut chart K=diag(1,-1), R=[[0,-1],[1,0]], L=RK, so
K,L are self-dagger, R is skew, and [K,L]=-2R. This chart differs by
basis choice from the earlier conversation's swap-matrix K. The central
iota used in the compact lift commutes with these internal words.

The exact compact spectral formula is inherited from YM52, which credits
Lauret, *The smallest Laplace eigenvalue of homogeneous 3-spheres*,
Theorems 1.1--1.2, https://arxiv.org/abs/1801.04259v3. The present additions
are the integrable constitutive family, its symmetry/kernel interpretation,
its explicit curvature-budget specialization, and a direct fixed-profile YM55
application. Neither the general spectral formula nor the even/odd
decomposition is claimed as new mathematics.

## YM56-T1: an integrable positive Hessian realizes the shape family

In centered dimensionless coordinates choose
\[
 u_\lambda(x,y)=10x-10y+\tfrac12(x^2+y^2)
              +\tfrac{\lambda+2}{6}x^3+\tfrac\lambda2xy^2.
\]
Its Hessian is
\[
 H_\lambda=\begin{pmatrix}1+(\lambda+2)x&\lambda y\\
                         \lambda y&1+\lambda x\end{pmatrix}.
\]
At the origin H=I. It is positive in a neighborhood for every fixed real
lambda. More quantitatively, for |lambda|<=B and
|x|,|y|<=r=1/[4(B+1)], the maximum absolute row sum of H-I is at most
(2B+2)r=1/2. Symmetry gives H>=I/2 there. Temperature u_x and pressure
-u_y remain positive after possibly shrinking this chart. Positive entropy
and volume origins can be restored by translations.

At the origin the two full response derivatives are
\[
 X_x=H_x=\operatorname{diag}(\lambda+2,\lambda),\qquad
 X_y=H_y=\lambda L.
\]
Their scalar parts are 1+lambda and zero. Hence their trace-free shapes
are Y_x=K and Y_y=lambda L, with positive factor B_factor=I. The
mixed third derivatives agree because these jets came from a potential.
The scalar part is essential for this integrability; simply taking
H_x=K,H_y=lambda L would generally violate it.

For NT's selected strain connection A_i=X_i/2,
\[
 \boxed{F_{xy}=-\tfrac14[X_x,X_y]=\tfrac\lambda2R,
             \quad f=\lambda/2,\quad \tau=2(1+\lambda^2).}       \tag{1}
\]
All statements in (1) are at the reference state. The Hessian metric is
locally a legitimate NT tangent adapter; at that state its Gaussian
curvature is -lambda/2 in NT's convention. Ordinary exact-differential
Maxwell identities continue to hold at nonzero lambda.

## YM56-T2: order survives as bracket control while symmetry is preserved

YM54's four equally counted signed labels, with compact turn directions
v_1=(1,0,0), v_2=(0,lambda,0), give
\[
 C_\lambda=\tfrac12\operatorname{diag}(1,\lambda^2,0),\qquad
 \mathcal L_\lambda=-\tfrac12(D_1^2+\lambda^2D_2^2).             \tag{2}
\]
These coefficients follow from actual label second moments, not an
identification of the response H with the energy operator.

On the compact quaternion carrier let Q act by pullback of conjugation
by e_1. Then Q^2=I and QD_1Q=D_1, QD_2Q=-D_2, QD_3Q=-D_3.
Thus Q commutes with both (2) and the symmetric finite-turn readouts
for every lambda; the second pair of labels is merely exchanged.
Antipodal parity is a separate symmetry, also present for every lambda.
No physical gauge-symmetry identification is made for Q.

Meanwhile [D_1,lambda D_2]=lambda D_3. The first two directions and
their bracket span the internal tangent directions precisely when
lambda is nonzero. YM52's energy is
\[
 \mathcal E_\lambda(\psi)=\tfrac12\|D_1\psi\|^2
                           +\tfrac{\lambda^2}{2}\|D_2\psi\|^2.
\]
The source all-content result proves constants are its only stationary
sector for lambda!=0. At zero, nonconstant centered functions invariant
under the first turn remain stationary. No claim that symmetry implies
zero curvature or zero gap is valid even in this family.

The family is even in lambda at the heat-generator level but odd in its
oriented curvature marker: C_{-lambda}=C_lambda, f_{-lambda}=-f_lambda.
Squaring the turn directions loses that orientation information.

## YM56-T3: exact opening rate and sharp witnesses

The ordered eigenvalues of C are 0, min(1,lambda^2)/2,
max(1,lambda^2)/2. Substitution in YM52-T4 proves the opening formula.
The odd antipodal sector has rate (1+lambda^2)/8 and the centered even
sector has rate min(1,lambda^2)/2. These are observer sectors, not the
even/odd splitting under Q.

For an explicit witness, let q_0,...,q_3 be quaternion coordinates and
\[
 w_j=\sum_{k\in\{1,2,3\}\setminus\{j\}}q_k^2
                              -\tfrac12\sum_{k=0}^3q_k^2.
\]
Then Phi(w_j)=0, Phi(w_j^2)=1/12, D_jw_j=0 and
D_k^2w_j=-w_j for k!=j. Thus w_1 has rate lambda^2/2 and w_2
has rate 1/2. Each linear coordinate has rate (1+lambda^2)/8.
Together with YM52's block lower bounds, these witnesses prove sharpness
on the completion, not just a numerical truncation.

In particular gamma=lambda^2/2 for |lambda|<=1/sqrt(3). At lambda=0,
w_1 witnesses zero gap above constants. A rank-one generator still has
positive eigenvalues above its larger kernel; these statements differ.
At lambda=1 the full rate is 1/4, while the even rate is 1/2.

## YM56-T4: recover GE2's curvature-budget law and specialize its plateau

[GE2-T6](../generalized-euler-evolution/GE2_CLOCK_FREE_RECOGNITION_GENERATOR.md)
already gives the normalized curvature law. The following calculation
recovers it with chi^2 equal to GE2's kappa and relates it to this
chapter's explicit Hessian family; it is not a newly discovered general law.

For any source admitted by YM54 with tau>0, put
\[
 \chi=8|f|/\tau,\qquad 0\leq\chi\leq1.
\]
The bound follows from det G=16f^2 and tr G=tau. The nonzero eigenvalues
of C are g_-/4,g_+/4, where
g_\pm=(tau +/- sqrt(tau^2-64f^2))/2. Hence
\[
 \boxed{\frac{16\gamma_{\rm full}}\tau
 =\min\{1,\ 2(1-\sqrt{1-\chi^2})\}.}                       \tag{3}
\]
Zero f includes the rank-one case; tau=0 is handled separately as zero
protocol, with no division by tau. At fixed tau the full rate reaches
its maximum tau/16 as soon as chi^2>=3/4. Thus a whole anisotropic
range attains the same full rate; it does not select a unique symmetry.
The even rate reaches its unique shape-balanced maximum tau/8 at chi=1.

For the cubic family, chi=2|lambda|/(1+lambda^2), so the full-rate
plateau after normalization by tau is 1/3<=lambda^2<=3.
This is a substitution into the existing spectral theorem, with an
explicit response-curvature interpretation. It is neither a unique
physical clock nor a bound determined by curvature without a budget.
GE2's intrinsic record clock divides the old rate by m=tau/4; thus
gamma_rec=4 gamma_full/tau and (3) is exactly 4 gamma_rec.
Indeed lambda and 1/lambda have the same chi; their unnormalized clocks
and rates differ unless their trace budgets are matched.

## YM56-T5: fixed-profile interacting inheritance with a sharper floor

At each chain site choose and freeze a member of T1 with
\[
 0<\ell\leq|\lambda_i|\leq B<\infty.
\]
Allow fixed compact orientations and use the same profile in every finite
interval. Formula (2) directly gives
\[
 \boxed{\beta=\tfrac12\min(1,\ell^2),\qquad
 M=\tfrac12(1+B^2),\qquad |\theta|<\min(1,\ell^2)/4800.}     \tag{4}
\]
These are YM55's second-eigenvalue, trace and interaction hypotheses.
Thus its positive/reflection-positive local-history joint limit and
centered correlation bound apply, with gamma_chain=beta log R/(10J)
for a passing YM55 gate. This is an application of that written result,
not a new proof of the joint-limit theorem.

For comparison, the curvature-budget-only route uses f0=ell/2 and
T0=2(1+B^2), giving beta_curv=ell^2/[2(1+B^2)]. Direct knowledge of
this family permits the stronger beta in (4); no upstream estimate is
changed. For ell=1/2, B=2, beta=1/8, M=5/2 and theta=1/40000
satisfy (4). The generic curvature-only floor is 1/40 and its sufficient
window would refuse this theta. A refusal is not proof of gaplessness.

The free formula in T3 is not the interacting chain's exact gap. Nor is
the free degeneracy at lambda=0 automatically an interacting phase
transition. The sufficient interacting window collapses as ell tends
to zero; this comparison theorem alone decides no phase at its boundary.

## YM56-T6: the remaining coercivity question is quantitative

Each member with lambda!=0 is gapped above constants, but lambda_n=1/n
gives gamma_n=1/(2n^2) for n>=2. At fixed nonzero curvature, taking
shape directions tK,L/t also makes the small protocol eigenvalue tend
to zero while the response budget diverges, as YM54 already exhibits.
Even rank two at every site is insufficient for a uniform free-chain
gap if |lambda_i| tends to zero: a centered w_1 localized at that site
has Rayleigh quotient lambda_i^2/2. It also rules out a uniform centered
time-correlation decay rate for that free profile.

A proposed physical mass-gap argument must therefore supply a uniform
energy inequality on the physical vacuum complement, with appropriate
gauge constraints, through the relevant spatial cutoff limit. Symmetry,
nonzero local curvature and finite-cutoff positivity do not alone supply
that inequality. For a selected physical Hilbert space the target has
the form <psi,H_phys psi> >= Delta ||psi||^2 for psi orthogonal to the
vacuum, with Delta>0 uniform in the required approximation.

The present chapter fixes no physical seconds, hbar, gauge action or mass.
YM55 retains a spatial chain lattice; this work does not turn it into a
four-dimensional Yang--Mills continuum, solve row closure, or settle the
Jaffe--Witten problem. The new result is a concrete integrable source,
an exact symmetry-preserving opening law, and its scoped chain inheritance.

## Reproduction and evidence

Run with Python 3.12 (other versions and optimized assertion removal are
refused):

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym56_lambda_curvature_gap.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym56_lambda_curvature_gap.py -v
```

The certificate differentiates the cubic independently, reuses native
response and quaternion engines, checks exact polynomial sharp modes,
finite-spin energy inequalities, symmetry, orientation loss, the
normalized-rate identity and YM55's stricter profile gate. Rational
controls include zero, signs, rank loss, trace scaling and vanishing
uniform floors. Source hashes bind the inherited files and new inputs.
General statements rely on the displayed proofs and source theorems;
finite controls are not proof-assistant verification or external review.
