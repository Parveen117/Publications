# YM-61: independent sources, higher response, and collective spectral probes

Monty Dabas. 4 October 2026. **Python 3.12 only.**

The response diagram becomes a calculation engine when its directions
are independent source variables, its mixed derivatives have a specified
potential, and its inverse steps preserve the native algebra. This
chapter develops that mechanism for the actual YM59--60 scalar model.
It computes spatial correlations beyond a single bond and uses them
together in collective spectral probes.

For any spatially averaged history limit constructed in YM60, write
\(c_r(\theta)=\langle q_0\cdot q_r\rangle\), with c_0=1. At unit
neighbour rate the new exact responses are
\[
 \begin{array}{ll}
 c_1=\theta/2+(1123/11844)\theta^3+O(|\theta|^5),&
 c_2=(11/24)\theta^2+O(\theta^4),\\
 c_3=(49/108)\theta^3+O(|\theta|^5),&
 c_4=(4853/10368)\theta^4+O(\theta^6).
 \end{array}                                                \tag{1}
\]
Here c_1 is inherited from YM60. The other three are new. The error
constants are finite and uniform over the constructed history limits
at each fixed separation; they are not numerically evaluated. No
convergence radius uniform in separation is asserted.

There is an exact cancellation in the leading collective response.
Let m_L be the best variational upper constraint from linear coordinate
readouts on L consecutive sites, with L fixed. For every k>0,
\[
 \boxed{m_L(\theta)=\frac{1+2k}{8}
       -\frac12\cos\!\left(\frac{\pi}{L+1}\right)|\theta|
       +O_L(\theta^2),\qquad L\ge2.}                        \tag{2}
\]
The coefficient of abs(theta) is independent of k, although the free
kinetic operator and the individual source susceptibility depend on k.
This is a model identity in the retained heat clock, not a dimensional
physical constant. m_L is an upper spectral constraint, not a claimed
mass eigenvalue or a positive lower gap.

## Carrier, lineage and scope

Keep the native compact quaternion carrier Q, reference Phi, derivative
algebra and extensive-rate clock from YM50--60. On the open N-site chain,
\[
 H_{N,k}=-\tfrac12\sum_i A_i^2-\tfrac k2\sum_i(Z_i-Z_{i+1})^2,
 \quad A_i=D_{1,i},\ Z_i=D_{2,i},\quad
 A_N(\theta)=H_{N,k}-\theta\sum_i q_i\cdot q_{i+1}.           \tag{3}
\]
Use N>=2 and k>0. The actual positive vacuum is h_N, not a chosen
polynomial or Gibbs ansatz. Auxiliary sources below are set to zero
after taking the response; the final correlations refer to (3).

Mixed-potential derivatives, resolvent perturbation, linked clusters
and variational structure-factor bounds are established methods. For
background on the latter, see [Feynman, Phys. Rev. 94 (1954), 262](https://doi.org/10.1103/PhysRev.94.262);
for an example of high-order linked-cluster work, see
[Coester et al., Phys. Rev. B 94 (2016), 125109](https://arxiv.org/abs/1605.05848).
Neither paper's model, gap or convergence theorem is imported here.
The model-specific content is the native source recursion, its parity
selection, exact site elimination and resulting collective coefficients.
Written limiting proofs and exact finite checks are distinct. This is
not proof-assistant certification or a physical 4D Yang--Mills construction.

## YM61-T1: source reciprocity and a multivariable Poisson tower

For fixed r>=1 introduce
\[
 W_{N,r}=\sum_{i=1}^{N-r}q_i\cdot q_{i+r},\qquad
 A_N(\theta,\eta)=H_{N,k}-\theta V_N-\eta W_{N,r}.           \tag{4}
\]
An empty sum is zero. More generally finitely many independent scalar
sources may be included. Finite-volume vacuum simplicity from YM60
persists for every real choice of these bounded potentials. If
A=H-sum_a g_a O_a, E is its simple ground energy, K=A-E, and
Q_h=I-|h><h|, then
\[
 \begin{split}
 -\partial_a E&=\langle O_a\rangle_h,\\
 -\partial_a\partial_b E
 &=2\operatorname{Re}\langle Q_h O_a h,
                   (K|_{h^\perp})^{-1}Q_h O_b h\rangle.
 \end{split}                                                \tag{5}
\]
Differentiate the eigenvalue equation with <h,partial_a h>=0 to
obtain (5). Thus the response matrix is symmetric and positive
semidefinite. The Maxwell-type mixed-derivative reciprocity concerns
these scalar sources; it is not a derivation of electromagnetic Maxwell
equations. Noncommuting native derivatives are fully compatible with
this integrability of a scalar energy potential.

Splitting an open chain costs at most abs(theta)+r abs(eta) in crossing
scalar bonds, while the crossing kinetic square is positive. The
superadditivity argument of YM59 therefore gives
\[
 e_r(\theta,\eta)=\lim_{N\to\infty}E_N(\theta,\eta)/N       \tag{6}
\]
for all real theta,eta. It is jointly concave. Its restriction eta=0
is the existing e(theta).

Use the same native carré du champ Gamma and local inverse R=H^-1
on centered sources. With U=theta V+eta W, the homogeneous logarithmic
recursion is
\[
 H F_1=U,\qquad
 G_n=\sum_{a+b=n}\Gamma(F_a,F_b),\quad
 E_n=-\Phi(G_n),\quad H F_n=G_n-\Phi(G_n).                  \tag{7}
\]
Each F_n is homogeneous of source degree n and has native coordinate
degree at most 2n. At fixed n,r, each local term has interval span at
most n(r+1); this coarse bound allows both source spans and intervening
kinetic edges. There are finitely many relative and physical-boundary
patterns. Centering and R preserve their supports, by YM59-T1. The
coefficient-norm argument used in YM59 then gives, for every fixed p,
\[
 \left|E_N(\theta,\eta)/N-P_{p,N}(\theta,\eta)\right|
 \le C_{p,r,k}(|\theta|+|\eta|)^{p+1},\quad
                                      |\theta|+|\eta|\le1.  \tag{8}
\]
The coefficients of P_(p,N) have a bulk limit with O(1/N) boundary
corrections. Thus (8) passes to e_r. The constants may grow with p and r;
this proof does not sum the infinite tower.

For orientation, the bulk quadratic source-response coefficients at zero are
\[
 \chi_r=\frac1{2\omega_r},\qquad
 \omega_1=\omega=\frac{1+3k}{4},\quad
 \omega_r=\nu=\frac{1+2k}{4}\ (r\ge2).                    \tag{9}
\]
They are the limits of the normalized finite-volume Hessian entries;
no differentiability of the bulk energy on a neighbourhood is assumed.
Different separations have zero mixed quadratic coefficient. Indeed,
H(q_i dot q_j)=omega_(j-i)(q_i dot q_j) in the bulk, its reference
variance is 1/4, and different endpoint pairs are orthogonal by site
parity. These are susceptibilities to the independent sources, not
claims that c_r(theta) starts linearly in the original nearest-bond theta.

## YM61-T2: parity determines the first possible spatial response

At each site the antipodal map q_i -> -q_i commutes with H, its centered
inverse and site averaging. Give every scalar bond its own source
label temporarily. A coefficient of the energy can be nonzero only
if every site has even total incidence in its list of inserted bonds.
This follows by changing the sign of just that site's coordinates.
Kinetic coupling does not remove the rule: it preserves each site's parity.

Consider the coefficient linear in eta for a source bond from 0 to r.
Its odd endpoint incidences must be completed by nearest-neighbour
source insertions. Modulo two, the only finite chain solution has each
of the r bonds between 0 and r inserted an odd number of times and
all other nearest bonds an even number of times. Consequently
\[
 [\theta^j\eta]e_r=0\quad\text{if }j<r
                       \text{ or } j-r\text{ is odd}.       \tag{10}
\]
At j=r the r path bonds occur exactly once. Together with the eta bond
they form a labelled cycle; r=1 is the two-parallel-edge case.
There is one such cycle per translate. This identifies the leading
coefficient a_r without a factorial ambiguity:
\[
                   [\theta^r\eta]e_r=-a_r.                \tag{11}
\]
The rule gives an earliest allowed order, not a general proof that its
coefficient must be nonzero. Nonzero values are checked here for r=1,...,4.

For the original eta=0 histories, Hellmann--Feynman and the spatial
average give
\[
 -\partial_{\eta,-}e_r(\theta,0)
       \le c_r(\theta)\le
 -\partial_{\eta,+}e_r(\theta,0).                          \tag{12}
\]
There are at most r exterior-crossing readouts in the finite spatial
average; their contribution vanishes as N grows. Concave secants
justify (12) even at a nondifferentiability point.

To turn (11) into a uniform asymptotic for all selected history limits,
use (8) with p=2r+3 and eta step h=abs(theta)^(r+2). Its secant remainder
is O(h)+O(abs(theta)^(2r+4)/h)=O(abs(theta)^(r+2)). The polynomial
coefficient rule (10) then gives
\[
                 c_r(\theta)=a_r\theta^r
                                   +O(|\theta|^{r+2}).      \tag{13}
\]
For c_1's cubic term use the higher-order YM60 secant argument. These
are fixed-order bulk response jets. Smoothness of every derivative on
a whole coupling interval, or convergence of their infinite series,
does not follow merely from the existence of all these asymptotics.

## YM61-T3: an exact cycle recursion with justified site elimination

Let C be the r path edges plus the independently labelled closing edge.
For S a subset of these labels, let d_i(S) be its incidence at site i,
and let P_S average precisely those sites with d_i(S)=2. Set f_empty=1.
For a nonempty proper subset S define
\[
 s_S=P_S\sum_{b\in S}v_b f_{S\setminus\{b\}},\qquad
 f_S=R s_S,\qquad
 a_r=\Phi\!\left(\sum_{b\in C}v_b f_{C\setminus\{b\}}\right).
                                                               \tag{14}
\]
All operators here are the actual bulk marginal of H on these finite
sites, including the external endpoint noises. Every proper subset of
a cycle has an odd-incidence endpoint. Therefore s_S is centered,
and every energy coefficient on a proper subset is zero. Expanding the
finite-volume eigenvalue equation with intermediate normalization
Phi(psi)=1 gives (14), initially without P_S; the usual normalization
subtractions vanish for this particular squarefree cycle coefficient.
This does not authorize dropping them for repeated-edge coefficients
such as YM59's quartic energy.

Here is why P_S is exact. A site with two inserted incident cycle edges
will never be touched by any remaining scalar insertion. Its averaging
projection commutes with H and R. On removing a last edge b, every site
already completed in S without b is disjoint from b; thus its projection
also commutes with multiplication by v_b. One may move all completed-site
projections back through the remaining calculation without changing
the final mean. Averaging a site with only one inserted edge would erase
an active endpoint and is wrong.

After this elimination, each active site has coordinate degree one;
inactive sites have degree zero. On the active set T,
\[
 H|_T=\frac{1+2k}{8}|T|\,I
               +k\sum_{\{i,i+1\}\subset T}Z_iZ_{i+1}.      \tag{15}
\]
The mixed terms must be retained. They commute with each other on this
multilinear space, but eliminating them changes the native inverse.
The implementation uses the already certified polynomial derivative,
reference pairing and exact Krylov solver, rather than a new operator.

The exact results at k=1 are

| Separation r | a_r | Proper-subset inversions | Largest Krylov dimension |
|---|---:|---:|---:|
| 1 | 1/2 | 2 | 1 |
| 2 | 11/24 | 6 | 1 |
| 3 | 49/108 | 14 | 2 |
| 4 | 4853/10368 | 30 | 2 |

There is also an independent closed calculation at r=2:
\[
 a_2(k)=\frac1{4\omega\nu}+\frac1{8\omega^2}.               \tag{16}
\]
Use Phi[(q_0 dot q_2)(q_0 dot q_1)(q_1 dot q_2)]=1/16 and sum
the six labelled insertion orders with their correct endpoint rates.
At k=1, (16) is 11/24. Incorrectly assigning the adjacent-bond rate to
the separated endpoints gives 3/8 instead. Factorizing c_2 as c_1^2
would give 1/4. Neither substitution reproduces the actual response.

## YM61-T4: collect the responses into a spectral matrix problem

In any YM60 bulk history limit, let
\(F_b=\sum_{j=1}^L b_j q_j^a\), for real b and a fixed quaternion
coordinate. The exact symmetry and form identities of YM60 give
\[
 \begin{split}
 G_L{}_{ij}&=c_{|i-j|},\\
 M_L&=\frac18\big[(1+2k)I-kc_1 J_1\big],\\
 \langle F_b^2\rangle&=\tfrac14 b^T G_L b,\\
 \int\lambda\,d\rho_{F_b}(\lambda)&\le\tfrac14 b^T M_L b.
 \end{split}                                                \tag{17}
\]
J_r is the L by L matrix with ones at |i-j|=r and zeros elsewhere.
The last inequality allows loss of first moment in the history limit;
it is not silently replaced by equality. G_L is positive semidefinite
because it is an actual covariance matrix. If b^T G_L b>0, Jensen yields
\[
 \frac{\langle F_b(0)F_b(t)\rangle}{\langle F_b^2\rangle}
 \ge \exp\!\left[-t\frac{b^T M_L b}{b^T G_L b}\right].       \tag{18}
\]
Near theta=0 at each fixed L, G_L=I+O(theta) is positive definite, and
\[
 m_L=\lambda_{\min}(B_L),\qquad
                  B_L=G_L^{-1/2}M_L G_L^{-1/2}.             \tag{19}
\]
This is an optimized variational constraint on spectral thresholds,
not the restriction of the full Hamiltonian to an invariant coordinate
subspace. Multiplication by the interaction mixes higher degrees.
Padding b by a zero shows the exact nesting m_(L+1)<=m_L wherever
these minima are defined. The minimizing readout can depend on theta.

For equal-amplitude phase readouts one obtains the finite-window form
\[
 S_L(p)=1+2\sum_{r=1}^{L-1}(1-r/L)c_r\cos(rp),\qquad
 R_L(p)=\frac{1+2k-2k(1-1/L)c_1\cos p}{8S_L(p)}.            \tag{20}
\]
Its use requires S_L(p)>0. Complex phase readouts use the corresponding
Hermitian correlation; alternatively split into real sine and cosine
readouts. No infinite structure-factor sum or summability is assumed.
L=2 with p=0 or pi recovers YM60's aligned two-site result.

## YM61-T5: a universal first response and noncommuting higher ones

Put rho=(1+2k)/8 and a_1=2/(1+3k). Since c_2=O(theta^2),
\[
 G_L=I+a_1\theta J_1+O_L(\theta^2),\quad
 M_L=\rho I-\tfrac{k a_1}{8}\theta J_1+O_L(\theta^3).
\]
Consequently
\[
 \boxed{B_L=\rho I-\tfrac14\theta J_1+O_L(\theta^2).}       \tag{21}
\]
The cancellation is a_1(rho+k/8)=1/4. Thus the first response is the
same adjacency matrix for every positive k, in the stated normalization.
The sine vector u_j=sin(pi j/(L+1)) has the largest J_1 eigenvalue
2 cos(pi/(L+1)); the other sine modes give the rest of the spectrum.
For theta<0 the lowest mode is its staggered version. This proves (2).

At k=1 the higher matrix response through cubic order is
\[
 \begin{split}
 B_L={}&\tfrac38 I-\tfrac14\theta J_1
 +\theta^2\left[\tfrac18J_1^2-\tfrac{11}{64}J_2\right]\\
 &+\theta^3\left[-\tfrac{1123}{23688}J_1-\tfrac{49}{288}J_3
       +\tfrac{77}{768}(J_1J_2+J_2J_1)-\tfrac1{16}J_1^3\right]
 +O_L(\theta^4).                                         \tag{22}
 \end{split}
\]
For finite windows J_1 and J_2 do not generally commute. The two
orders in their anticommutator are essential to a symmetric B_L.
This supplies a concrete higher-response role for noncommutativity,
with both its location and coefficient fixed by the calculation.

The optimal branch through quadratic order follows without diagonalizing
the full interacting Hamiltonian. Write z=pi/(L+1). The lowest
first-response mode is simple on each side of theta=0. Since
J_1^2=J_2+D with D=diag(1,2,...,2,1), ordinary finite-matrix perturbation
applied to (22) gives, for L>=2,
\[
 m_L=\tfrac38-\tfrac12\cos z\,|\theta|
 +\left[\tfrac5{32}+
        \frac{3-11/(L+1)}{16}\sin^2 z\right]\theta^2
 +O_L(|\theta|^3).                                        \tag{23}
\]
For L=2 this is 3/8-|theta|/4+theta^2/8. For L=3 it is
3/8-(sqrt(2)/4)|theta|+(21/128)theta^2. These expansions do not
assert differentiability of the optimized even branch at zero.

An entirely rational four-site readout uses b=(1,2,2,1), staggered at
negative theta. Directly from (1) and (17), its upper constraint is
\[
 R_b=\tfrac38-\tfrac25|\theta|+\tfrac{73}{400}\theta^2
                  -\tfrac{1293529}{11844000}|\theta|^3
                  +O(\theta^4).                           \tag{24}
\]
It improves YM60's two-site bound for all sufficiently small nonzero
coupling, since the leading improvement is -(3/20)abs(theta).
The size of that interval is not numerically certified here.

## YM61-T6: the missing uniform step remains explicit

Equations (13), (22) and (23) hold for every fixed r or L. Their proof
does not give constants uniform as r,L grow. One therefore cannot sum
(13) into an infinite structure factor, interchange an arbitrary tower
limit with a volume limit, or infer exponential spatial clustering.
Likewise the limit of the coefficient in (2) is -1/2, but that fact
alone does not give a uniform O(theta^2) error at L=infinity.

There is a valid weaker conclusion. Define m_*(theta) as the infimum
of all nonzero finite-coordinate Rayleigh quotients for a selected
history state, allowing null vectors to be omitted. At theta=0 it equals
rho. From (2), taking the small-coupling limsup for each fixed L first,
\[
 \limsup_{t\downarrow0}\frac{m_*(t)-\rho}{t}\le-\tfrac12.    \tag{25}
\]
The estimate holds for any selection of the constructed states, because
the fixed-L errors above are uniform over those selections. It is an
upper constraint, not a matching derivative or an actual mass formula.

If a genuine mixing gap Delta>0 were proved for the full reconstructed
vacuum, (17) would necessarily imply the matrix inequalities
\[
                         \Delta G_L\le M_L                \tag{26}
\]
for every L, using the spectral first-moment inequality. Conversely,
checking (26) only on coordinate probes would control only that family;
it would not prove the weighted Poincare inequality for all readouts.
The missing lower-gap step is still control of the actual interacting
vacuum uniformly in width, or uniform decay for the full centered
history algebra. A convergent local source construction or another
suitable uniform argument would be needed to complete it.

The diagram organizes a hierarchy of exact operations here: scalar
potential, reciprocal response, parity-selected mixed insertions, native
inverse, collective covariance and spectral constraints. This is a
specified working instance of the proposed organizing idea. It does
not establish that placement alone determines every law of physics.
Physical gauge constraints and the four-dimensional continuum remain
separate obligations.

## Exact evidence

The certificate records all 52 proper-subset inverse witnesses for
r=1,...,4, checks the independent triangular-cycle formula and source
susceptibilities, and verifies the matrix congruence through cubic
order without floating-point square roots. It also tests elimination,
parity, noncommuting boundary matrices, source normalization, rational
collective readouts and invalid covariance refusal. The written proof,
not numerical sampling, supplies the uniform fixed-order and history
statements.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym61_source_response_collective.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym61_source_response_collective.py -v
```

Default and --check are read-only; only --write updates this chapter's
result and hash. The manifest preserves 160 upstream files. The dedicated
Python 3.12 workflow also replays YM60 and YM59.
