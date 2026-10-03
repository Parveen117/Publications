# GE4: recover the smallest observer-memory model from native response

3 October 2026. Python 3.12 only.

GE3 showed that an observer's first derivative does not determine its
retained evolution. GE4 gives a constructive recovery test. In the native
quadratic example, the known generator ceiling and **three** initial
response derivatives certify the entire one-memory law. More generally,
a finite moment-rank test gives the smallest number of linear memory
coordinates. An explicit error bound covers approximate two-state closure.

## Admitted source and provenance

Work on a fixed finite native coefficient core with its already derived
positive pairing. Let L be self-adjoint, 0<=L<=M I, M>0, and normalize
a declared probe u to norm one. No primitive Hilbert space is assumed.
The native heat and clock are those of GE1–GE3. Write

\[
c(\tau)=\langle u,e^{-\tau L}u\rangle,\qquad
m_k=(-1)^kc^{(k)}(0)=\langle u,L^ku\rangle,\quad m_0=1.       \tag{1}
\]

A constant nonzero multiplicative readout gain cancels by normalizing at
zero; additive backgrounds must first be removed.
All derivatives and M must use the same native clock units. This packet
does not infer physical seconds or certify differentiation of noisy data.
Its examples are exact arithmetic fixtures, not experimental measurements.

The input is assumed to come from **one fixed admitted source**. Exact
finite moment equalities then imply the vector relations proved below.
Finite matrices alone do not certify that an unknown physical process is
linear, stationary, self-adjoint or in this source class.

## GE4-T1: moment rank is the minimal linear memory count

Let V_n=span(u,Lu,...,L^(n-1)u) and form

\[
H_n=(m_{i+j})_{0\le i,j<n},\quad
N_n=(m_{i+j+1})_{0\le i,j<n}.                               \tag{2}
\]

H_n is the Gram matrix of these native word probes, so its rank equals
dim V_n. If H_r is positive definite and rank H_(r+1)=r, then
V_(r+1)=V_r and L V_r is contained in V_r. Self-adjointness makes V_r
reducing. In the word basis u,Lu,...,L^(r-1)u its pairing and generator are

\[
G=H_r,\qquad \widehat L=H_r^{-1}N_r.                        \tag{3}
\]

Thus G Lhat=N_r is symmetric positive and M G-N_r is positive. Starting
from e_0 and observing e_0^dagger G reconstructs c(tau) exactly for all
tau. There is no need to introduce square roots to implement this basis.
The full linear state dimension is r; retaining u as the visible direction
leaves **r-1 orthogonal memory coordinates**.

**Proof of minimality and uniqueness.** Any other admitted realization
with these moments has the same Gram matrix on its word vectors. It
needs at least r dimensions. The map between corresponding word vectors
preserves their pairing; the flat-rank relation makes it intertwine L.
It therefore identifies the minimal cyclic realizations by a pairing
isometry. This proves uniqueness only up to that change of hidden basis.
Extra modes orthogonal to the cyclic probe sector remain unidentifiable.
The minimum concerns autonomous linear positive-pairing realizations,
not arbitrary nonlinear encodings or externally prescribed time laws.

All moments through m_(2r) suffice to check the stated exact flatness
condition. A merely small determinant is not an exact rank certificate.
This is classical Gram/Krylov realization algebra instantiated in the
native carrier, not a claim of new general realization theory.

## GE4-T2: the native ceiling makes three derivatives sufficient

Put a=m_1 and q=m_2-a^2. If q=0, then ||(L-a)u||^2=0, so the probe
already closes: c(tau)=exp(-a tau) and no memory coordinate is needed.
Suppose q>0. The source inequalities imply 0<a<M and
a^2<m_2<=Ma. Define entirely from m_1,m_2,m_3 and the known ceiling

\[
r=\frac{Ma-m_2}{M-a},\quad
\Xi=-m_3+(M+2r)m_2-(2Mr+r^2)a+Mr^2.                         \tag{4}
\]

Then 0<=r<a<M and

\[
\Xi=\langle(L-r)u,(M-L)(L-r)u\rangle\ge0.                   \tag{5}
\]

If Xi=0, positivity of M-L implies (M-L)(L-r)u=0. Hence
L^2u=(M+r)Lu-Mr u and the cyclic source has exactly two dimensions.
The full visible response and its GE3 memory kernel are

\[
c(\tau)=w e^{-r\tau}+(1-w)e^{-M\tau},\quad
w=\frac{M-a}{M-r},\qquad
K(s)=q e^{-\delta s},\quad \delta=M+r-a.                     \tag{6}
\]

**Proof.** Positivity of L and M-L gives m_2<=Ma; the variance gives
m_2>a^2. Equation (5) is its polynomial expansion. A zero quadratic
form for a positive finite operator implies that it annihilates the
vector, by a native positive-factor decomposition. This proves the
degree-two word relation. The projectors (M-L)/(M-r) and (L-r)/(M-r)
on this cyclic space give (6), using their initial weights. They are
positive there since the relation has roots r and M. In the orthogonal
word basis u,v=(L-a)u the pairing is diag(1,q) and the generator is

\[
\widehat L=\begin{pmatrix}a&q\\1&\delta\end{pmatrix}.         \tag{7}
\]

Eliminating the second coordinate by GE3 gives K(s)=q exp(-delta s).
The sign of a normalized hidden coordinate is a basis choice; the
observed moments do not identify that sign.

Xi=0 is a sufficient ceiling-based certificate, not a necessary condition
for every two-state source: two rates strictly below M can have Xi>0 and
still close exactly in two dimensions. The flat-rank test in T1 handles
that case. A positive Xi alone does not prove that extra memory is needed.

## GE4-T3: approximate closure with a proved error

For Xi>=0, retain q>0 and the same r. Whether or not closure is exact,
compress to V_2=span(u,(L-a)u). Its native pairing is G=diag(1,q), and
its generator still has the form (7), now with

\[
\delta=\frac{m_3-2am_2+a^3}{q}.                             \tag{8}
\]

This is the actual positive compression, so 0<=Lhat<=M in pairing G.
Let chat(tau)=e_0^dagger G exp(-tau Lhat)e_0. Then

\[
\boxed{|c(\tau)-\widehat c(\tau)|
 \le\tau\sqrt{M\Xi/q}.}                                   \tag{9}
\]

**Proof.** Set E=(M-L)(L-r)u. Since 0<=M-L<=M,
||E||^2<=M Xi. For the orthonormal hidden direction v/sqrt(q), the
only leakage from V_2 is (I-P_2)Lv/sqrt(q); its norm is bounded by
sqrt(M Xi/q), because subtracting a linear combination of u and Lu
does not change the perpendicular component of L^2u. The visible
direction has zero leakage. If J is the isometric inclusion from the
two-coordinate pairing, then ||LJ-JLhat|| has this bound. Integrating
the derivative of exp(-(tau-s)L)J exp(-s Lhat), and using both native
contractions, gives (9). No general unbounded-generator block is used.

For a degree-d GE2 rational history with fixed shape, budget tau and
maximum turn size epsilon, its scalar readout therefore differs from
chat(tau) by at most

\[
K_d\epsilon^2\tau+\tau\sqrt{M\Xi/q},\qquad
K_d=d^2(2d^2+1)/96.                                        \tag{10}
\]

The first term is GE2's coefficient-norm bound. Equation (10) uses a
probe normalized in that coefficient pairing; in another pairing use
its appropriate fixed-core norm conversion for the first term. The
native examples below have a common scalar pairing on their invariant
diagonal sectors, so the same unit-probe bound applies there.

The denominator q makes conditioning explicit. Small q requires its
own positive lower bound before using (9) with uncertain inputs. Small
Xi for exact moments is an approximation certificate; it is not by
itself a noisy-data derivative or exact-rank certificate.

## GE4-T4: native recovery, an extra-memory example and blindness

**Recover GE3.** Its degree-two ceiling is M=d^2/4=1. The first three
normalized moments are

\[
(m_1,m_2,m_3)=(3/4,5/8,9/16).
\]

Equations (4) give q=1/16, r=1/2 and Xi=0. Thus precisely one memory
coordinate is sufficient and necessary, delta=3/4, and the whole GE3
law and K(s)=exp(-3s/4)/16 are recovered from these data. The ordinary
flat-rank method also gives rank H_2=rank H_3=2, using m_4=17/32.
The ceiling test saves that fourth derivative in this case.

Even with M=1, the first two derivatives do not suffice: the positive
matrix [[3/4,1/4],[1/4,1/2]] obeys that ceiling and has the same m_1,m_2,
but a different m_3 and response. This is a positive-core comparison,
not another claim about the fixed GE3 polynomial sector.

**A genuinely larger native probe sector.** For
C=diag(1/6,1/3,1/2), let
f_j=x_0^2+x_j^2-rho/2, j=1,2,3. The unchanged native operator satisfies
L_C f_j=(1-C_jj)f_j. These three centered words are mutually orthogonal
with equal Phi norm. Their equal sum, normalized by its own norm, has
moments m_k=((5/6)^k+(2/3)^k+(1/2)^k)/3. Here

\[
a=2/3,\quad q=1/54,\quad r=11/18,\quad
\Xi=5/972,\quad M\Xi/q=5/18.
\]

H_3 is positive definite and rank H_4=3: two memory coordinates are
necessary and sufficient. Equations (3) reconstruct the exact three-state
law. The two-state compression remains a controlled approximation, not
an exact reconstruction. Native polynomial substitution verifies this
source; it is not a fitted three-exponential data set.

**Why the ceiling is load-bearing.** The positive comparison matrix

\[
L_\eta=\begin{pmatrix}3/4&1/4&0\\1/4&3/4&\eta\\0&\eta&3/4\end{pmatrix},
\qquad \eta=1/8,
\]

with probe (1,0,0), has the same first three moments as GE3, but a
different fourth moment and a three-dimensional cyclic source. It is
positive, but violates L<=I. It is an admissible positive-core
comparison, not another realization of GE3 under its native ceiling.
Dropping that ceiling would invalidate the three-derivative conclusion.

**What remains invisible.** Appending an orthogonal positive mode with
arbitrarily small rate epsilon changes the ambient gap while leaving
every moment and readout of the original probe unchanged. The augmented
matrix can still obey the same ceiling. Minimal observable memory is
therefore not the total hidden-state count or a global mass-gap proof.
This is a positive-core nonidentifiability example; no additional mode
of the fixed YM carrier is asserted to exist.

## Certification and lineage

The proof uses frozen GE1–GE3 and the native YM coefficient engine at
Publications commit `e470af355c187b5d9fb051a07025180afd44c7ba`.
The earlier [NI return-identification result](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-return-identification/THEOREM.md)
was reviewed at that historical commit: it identifies a paired-cell
parameter from two symmetric probes and is a different inverse problem.
It was not found at current main and is not represented as a current
local predecessor. R11's static marker-rank recovery motivates the
comparison; GE4 uses temporal word moments and a declared generator bound.

Gram/Krylov realization and minimal-state reconstruction have established
lineage, including B. L. Ho and R. E. Kalman,
[Effective construction of linear state-variable models from input/output functions](https://doi.org/10.1524/auto.1966.14.112.545)
(1966). Publisher bibliographic metadata was checked for attribution;
no external realization theorem is a proof premise. The native increment
is the ceiling-based three-derivative recovery, its explicit approximation
bound, and its docking to the already certified rational-arrow clock.
No general priority claim is made for the elementary moment identities.

Written fixed-core proofs and exact finite controls are separate evidence;
there are 14 new adversarial/exact tests, with 62 tests across GE1–GE4.
This is not formal proof-assistant or independent expert certification.
The actual interacting row, moving/unbounded observers, experimental
derivative inference, physical source/clock selection and 4D/Clay remain
open. The reconstructed visible sector alone cannot close those gates.
