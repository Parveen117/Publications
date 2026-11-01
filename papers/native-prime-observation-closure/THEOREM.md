# POC-1 — Equal-reading centres and prime-sampled source recovery

16 October 2026. Continues POS-1 in Publications and reuses the pinned
SE26, SE48 and SE49 RH source results. The geometry and the analytic
sampling claims have different carriers; no unproved identification
between them is made.

## Source contract

POS-1 classifies elementary additive observers of its declared integral
seam/shear algebra. It does not embed a prime as a unique point of a
metric plane. To make the three-point proposal testable, in T1--T2 we
choose the existing TET tangent frame with lattice Z^2, retained base
metric H_0=[[2,1],[1,5]], and squared-distance readings. In native NM1
coordinates the same statements are transported by S=I+R, as in POS-T4.
Calling a point prime-content means gcd of its adapted coordinates is
prime, using POS-T2's orbit invariant. This is a stated point model,
not the definition of a prime observer module.

T3--T5 instead use the actual completed RH response in SE48/SE49,
including its Gamma, prime and pole contributions. Sampling at a prime
integer s does not mean truncating the Euler product to that prime.
The common terminology does not supply a map from a circumcentre to a
zeta zero or from either object to mass.

## POC-T1. Exact equal-reading centre and failure of automatic prime closure

Let H be a real symmetric positive definite 2x2 matrix and let a_0,a_1,a_2
be noncollinear real columns. Define R_i(x)=(x-a_i)^T H(x-a_i). Set

\[
A=2\begin{pmatrix}(a_1-a_0)^T H\\(a_2-a_0)^T H\end{pmatrix},
\qquad b=\begin{pmatrix}a_1^THa_1-a_0^THa_0\\
a_2^THa_2-a_0^THa_0\end{pmatrix}.
\]

Then there is exactly one equal-reading centre, and it is

\[
\boxed{x_*=A^{-1}b,\qquad R_0(x_*)=R_1(x_*)=R_2(x_*).}
\tag{1}
\]

**Proof.** Expanding R_i-R_0 cancels x^T Hx and leaves
R_i-R_0=b_i-A_i x. The two displacement rows are independent;
positive H is invertible, so A is invertible. This proves existence
and uniqueness. QED.

The comparison cut C(r_0,r_1,r_2)=(r_1-r_0,r_2-r_0) has kernel the
whole common-reading line, not only zero. Retaining the mean alongside
the contrasts recovers all three readings. This is the three-channel
specialization of SE26-T1's common-mode/contrast distinction. Equality
of readings by itself supplies neither a zero radius nor an elementary
additive observer.

For the declared prime-content inputs

\[
a_0=(2,0),\qquad a_1=(0,3),\qquad a_2=(5,0),\qquad H=H_0,
\]

the exact result is

\[
A=\begin{pmatrix}-2&26\\12&6\end{pmatrix},\quad
b=\binom{37}{42},\quad \det A=-324,\quad
\boxed{x_*=\left(\frac{145}{54},\frac{44}{27}\right).}
\tag{2}
\]

The centre is not a lattice point, so it has no integer coordinate
content to classify as prime. Thus prime-content input does not imply
prime-content centre in this source-metric realization. This is a
counterexample to that specific automatic closure rule, not a claim
that every possible native reading must be squared distance.

More generally a common invertible frame map B transports points to
Ba_i, the metric to B^-T H B^-1, and the centre to Bx_*. All readings
are unchanged. In particular a TET path co-transports this construction
lawfully; holding the metric fixed under a nonisometric shear would be
a different observation contract. The POS holonomy boundary remains.

## POC-T2. Local prime centres have a determinant gate

Suppose H and all a_i are integral in the adapted lattice, so A,b in
(1) are integral, and det A != 0. For a prime p not dividing det A,
every depth r>=1 has exactly one solution to

\[
A x=b\pmod {p^r}.
\tag{3}
\]

The solutions are compatible under reduction, define one point of
Z_p^2, and agree with the p-adic reading of the rational centre A^-1 b.
For coprime good moduli the solutions combine by POS's CRT map.

**Proof.** The adjugate formula gives the inverse modulo p^r because
det A is a unit. Uniqueness makes all reductions compatible. Since
the rational inverse has denominators dividing det A before cancellation,
it is p-adically integral and represents those same solutions. The CRT
maps commute with integral A, proving the last assertion. QED.

At a bad prime, rank compatibility must be tested; a solution at the
first depth need not extend. For (2), modulo 2 there is no solution.
Modulo 3 the solutions are (0,2),(1,0),(2,1), but modulo 9 there are
none: six times the first equation plus the second gives
162 x_2=264, impossible modulo 9. Modulo 5 and 7 the unique solutions
are (0,2) and (1,5); their CRT combination modulo 35 is (15,12).

For every good prime the ambient observer V_p is elementary by POS-T1.
That is a statement about the chosen ambient module, not about its
centre point selecting p. The same rational centre has a consistent
reading at every good prime, with infinitely many possible choices.
This proves an exact local compatibility rule without asserting the
missing prime-centre selection law.

## POC-T3. Prime-indexed actual source pencils

Use exactly SE48's meromorphic response

\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q(t)=\mathcal X'(t)/\mathcal X(t),\qquad
F(w)=wq(w^2)=\tfrac12\xi'(1/2+w)/\xi(1/2+w).
\]

Let p_1=3,p_2=5,... be the increasing odd primes and set

\[
w_j=p_j-\tfrac12,\qquad t_j=w_j^2,\qquad u_j=q(t_j).
\tag{4}
\]

Thus the actual original s-readings occur at s=p_j. Define

\[
L_{ij}=-q[t_i,t_j],\quad J_{ij}=(tq(t))[t_i,t_j],\quad
L_{ii}=-q'(t_i),\quad J_{ii}=u_i+t_iq'(t_i).
\tag{5}
\]

For any finite N for which these actual matrices are positive definite,
the inherited SE48-T1 response

\[
q_N(z)=u^T(J_N+zL_N)^{-1}u
       =\sum_{k=1}^N\frac{a_k}{z+\mu_k},\qquad a_k,\mu_k>0
\tag{6}
\]

matches both actual q and q' at every retained t_j. Its poles are
finite-model poles, not independently certified zeros of xi.

For clarity the algebra is independent of integer spacing. With
T=diag(t_i), direct divided-difference arithmetic gives
J+LT=u1^T and J+TL=1u^T. The positive symmetric chart
L^-1/2 J L^-1/2 gives (6). If Jv=mu Lv, then
(T+mu I)Lv=1(u^Tv). A nonzero v cannot have u^Tv=0, so the
eigenvalues are simple and all residues positive. At each node
(J+t_iL)e_i=u, giving both interpolation identities by differentiating
the inverse. This is the existing SE48 theorem specialized to (4),
not a new spectral theorem.

Each next prime observation is accepted precisely by SE49-T1's two
Schur innovations. For x=t_{N+1}, r=q(x), d=q'(x), let

\[
\ell_i=(u_i-r)/(x-t_i),\quad j_i=(xr-t_i u_i)/(x-t_i).
\]

Then positivity extends if and only if

\[
\boxed{\alpha=-d-\ell^TL_N^{-1}\ell>0,\qquad
\beta=r+xd-j^TJ_N^{-1}j>0.}
\tag{7}
\]

This follows by two exact triangular congruences. Primality of the
new node does not make either inequality automatic.

## POC-T4. Prime sampling has sufficient analytic density

For a=0 or 1 and Re w>a define

\[
B_{N,a}(w)=\prod_{j=1}^N
\left|\frac{w-w_j}{w+w_j-2a}\right|^2.
\tag{8}
\]

For every compact K in Re w>a, with alpha_K=min_K Re w-a>0,
there is a finite constant C_K,a independent of N such that

\[
\boxed{B_{N,a}(w)\le C_{K,a}
       [\log(p_N+1)]^{-4\alpha_K}\longrightarrow0
       \quad\text{uniformly on }K.}
\tag{9}
\]

This is analytic uniqueness density, not topological density of the
prime nodes in the plane.

**Proof of the arithmetic input.** For X>=2 the finite Euler product
over p<=X contains every term 1/n for n<=floor X, so
product_{p<=X}(1-1/p)^-1 >= sum_{n<=X}1/n >= log(floor X+1).
The logarithm's terms beyond 1/p have total at most
sum_{n>=2}1/[n(n-1)]=1. Hence for integral X>=2,

\[
\sum_{p\le X}\frac1p\ge\log\log(X+1)-1.
\tag{10}
\]

Removing prime 2 costs 1/2, and 1/(p-1/2-a)>=1/p for the odd
primes when a=0,1. Thus the cumulative reciprocal sampling depth is
at least log log(p_N+1)-3/2. This elementary proof of prime-reciprocal
divergence is credited to classical Euler-product reasoning; no prime
number theorem, RH or prime-gap assumption is used.

**Proof of (9).** Write z=w-a and b_j=w_j-a. On a fixed compact,
for sufficiently large b_j,

\[
\log\left|\frac{b_j-z}{b_j+z}\right|^2
=-4\operatorname{Re}(z)/b_j+O_K(b_j^{-3}).
\]

The error sum converges since the b_j are a separated increasing
subsequence of positive half-integers. Finitely many initial factors
are at most one, because
|w+w_j-2a|^2-|w-w_j|^2=4(Re w-a)(w_j-a)>0.
At a zero of a factor the bound is immediate. Sum the expansion and
use (10), absorbing the finite initial sum and the bounded cubic tail
into C_K,a. This proves (9). QED.

Increasing observation count alone is insufficient. For geometrically
spaced w_j=2^j and w=1/8, the product in (8) at a=0 remains above
1/2: each factor is at least 1-2^(-j-1), and every finite product
is at least 1 minus the sum of those deficits, which is greater than
1/2. Repeated copies of old readings likewise add no new independent
conditions. On Re w=0 every factor in B_N,0 has modulus one, so
(9) is not a bound uniform up to the critical-line boundary.

## POC-T5. RH criterion on the prime bank, with a finite replay

For the actual source (4)--(5), the following are equivalent:

1. RH.
2. Both L_N and J_N are positive definite for every finite prime-bank N.
3. The initial two signs and both actual innovations (7) are positive
   at every successive prime observation.

Under these equivalent hypotheses q_N converges locally uniformly to
actual q on C minus the nonpositive real axis. With
c=(5/2)q(25/4)>0 and rho(w)=|(w-5/2)/(w+5/2)|,

\[
|q(w^2)-q_N(w^2)|
\le\frac{4c}{|w|(1-\rho(w))^2}B_{N,0}(w),\qquad \Re w>0.
\tag{11}
\]

This is an equivalent criterion and a conditional convergence theorem.
The all-order signs are not established here.

**Proof.** (2) iff (3) follows from (7). Under RH, the inherited
SE39/SE48 source product gives q(t)=sum_gamma m_gamma/(t+gamma^2).
The matrices in (5) are the corresponding Cauchy Grams with positive
weights m_gamma and m_gamma gamma^2. A nonzero finite rational
combination of the distinct node functions cannot vanish at infinitely
many distinct gamma^2; hence both Grams are strictly positive. This
proves (1) implies (2), with the same actual source as SE48/SE49.

For the reverse direction, put F_N(w)=wq_N(w^2). Equation (6) gives
positive real part on Re w>0, and F_N(5/2)=c. Thus
h_N=(F_N-c)/(F_N+c) maps that half-plane to the unit disk and
|h_N(w)|<=rho(w). For M>N, (h_M-h_N)/2 is bounded by one
and has double zeros at all w_1,...,w_N. Repeated half-plane
Schwarz division proves its bound B_N,0. Inverting the Cayley
transform gives exactly the finite-pair estimate

\[
|q_M(w^2)-q_N(w^2)|
\le\frac{4c}{|w|(1-\rho(w))^2}B_{N,0}(w).
\tag{12}
\]

By T4 these responses have a holomorphic locally uniform limit in
Re w>0. To identify it with the actual source without assuming RH,
reuse SE49-T2's independently proved Euler-half-plane majorant:

\[
|q(w^2)-q_N(w^2)|
\le \frac{|w+5/2|^2}{|w|} B_{N,1}(w),\qquad\Re w>1.
\tag{13}
\]

That proof uses only finite pencil positivity, the common anchor,
actual holomorphy in the Euler-product region and interpolation double
zeros. None of its steps requires consecutive integer probes. Its
bounded analytic function (F-F_N)/(w+5/2)^2 has norm at most one;
division by the selected prime-node Schwarz factors gives (13).
Thus T4 identifies the limit with actual q on Re w>1. The meromorphic
identity theorem extends this equality to Re w>0. Any xi zero with
Re s>1/2 would give a nonremovable logarithmic-derivative pole there,
a contradiction. Reflection of the actual xi source excludes left-side
zeros. Therefore RH follows. Letting M tend to infinity in (12)
proves (11). QED.

The usable unconditioned finite bound is (13) whenever that finite
prime pencil has been certified. Bound (11) for the actual full source
requires the all-order premise. Equal readings, a circumcentre, a
formal zero contrast or an arbitrarily assigned positive model cannot
substitute for the actual innovations (7).

### Executed finite source scope

The pinned SE49 archive contains complete outward enclosures for the
24 consecutive original s-readings 3,...,26. This packet extracts

\[
s\in\{3,5,7,11,13,17,19,23\}
\tag{14}
\]

and replays both 8x8 principal pencils with exact native F arithmetic
and directed rational interval rounding. Every LDL pivot is strictly
positive, so all eight nested prime-bank prefixes are certified from
the inherited source enclosures. This is a new selection/replay of
existing source data; it does not enlarge the old V_24 or recompute
the original special-function tails.

The smallest certified pivot lower bounds exceed 10^-19 for L and
10^-15 for J. These are coordinate-dependent residual pivots, not
invariant spectral gaps or physical energy thresholds.

The response at the held-out original s=27, w=53/2, is enclosed using
the same eight-prime pencil and compared with SE49's independent
held-out actual response box. The certificate stores the full rational
error enclosure and the separate finite analytic bound (13). Neither
number is a global eta estimate or a zero-height localization bound.

At that one held-out point the direct absolute error is less than
10^-15. The independent general analytic upper bound (13) is exactly
103938208563200/38778237372401894277. Both numbers refer to the same
eight-prime bank; the direct source-box comparison is sharper there.

The verifier also retains SE48's finite-to-global masking lesson in
an explicit rational control: add epsilon times
[(t+4-iota)^-1+(t+4+iota)^-1], epsilon=10^-8, to the toy response
1/(t+2)+2/(t+5)+3/(t+11). Both three-node prime pencils remain
positive, while the added poles -4+/-iota persist. This is a toy
countermodel to a finite-sign inference, not a claim about actual xi.

MP's observer-marker completeness theorem and recursive-residue theorem
are cited for the relevant boundaries: repeated matching measurements
need not be complete, and finite/selected-state closure does not imply
uniform reconstruction. The source-specific positive-real and Schwarz
argument above supplies the appropriate analytic convergence mechanism
under its explicit premise; it does not assume the stronger operator-
norm reconstruction hypotheses of that separate MP theorem.

## Result and next exact obligation

The automatic prime-centre rule fails for the declared squared-distance
model, while its correct prime-local compatibility is determined by
det A and lifting. On the actual RH response, the complete prime bank
does have sufficient analytic sampling density. RH is equivalently
reduced to both actual Schur signs at every prime observation. The
eight-prime finite replay is positive; the universal signs remain open.

No improved global eta, newly certified zeta zero, physical mass,
prime-selected spacetime lattice or full quantum-gravity closure is
claimed. A next arithmetic advance must control the next and ultimately
all actual prime innovations with their full source and tails retained.
