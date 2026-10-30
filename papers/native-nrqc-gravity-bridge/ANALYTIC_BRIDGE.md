# NRQC bridge: a current–density covariance benchmark from LMD
Author: Monty. 30 October 2026.

Status: written derivation and capability assessment, conditional on the frozen
LMD results. No new calculator execution, benchmark, certificate replay,
hardware run or continuum/gravity closure is reported.

## Sources and the question

- Publications parent: `35d63e4ac0d8cf161c549654ce6b430d68afbde4`.
- [LMD-1](../native-logarithmic-material-density/THEOREM.md), T3–T6:
  renormalized density, current Stone generators, mass subspaces and moments.
- [MPQ-1](../native-material-profile-quantization/THEOREM.md), section 2:
  the original separating chart and inverse.
- NRQC implementation contract:
  [EMK_GEOMETRY.md](https://github.com/Parveen117/quantum-classical-provisional/blob/fdd39048a6ddaa2dd78aa3197c274d7766f31a0e/native_relational_qc/EMK_GEOMETRY.md).
  Its geometry query takes signed Pauli words and returns the finite pure-state
  covariance tensor. The complex rank and real tangent rank are distinct.
- Standard background, not a novelty claim: Kolodrubetz et al.,
  [Geometry and non-adiabatic response in quantum and classical systems](https://arxiv.org/abs/1602.01062),
  Physics Reports 697 (2017).

Question: can NRQC advance the existing gravity problem rather than produce
an unrelated large-qubit example? First derive a target from the actual LMD
operators, then identify exactly what the present calculator can represent.

## 1. A specific target on the existing physical space

Keep LMD's positive orientation sector, dimensionless 0<h<1, covariance
C(z,y)=K(z-y)=1-log(2|sin(pi(z-y))|), and normalized vacuum psi=1.
Do not reinterpret h as a numerical SI value of Planck's constant. Set

    a=1+h/2,
    v(z)=sin(2 pi z)/(2 pi),   w(z)=cos(2 pi z),
    R=R_a[v],                J=J_h(w)=M(w).

These are the existing local material current and a signed window reading of
the positive Wick measure. A signed w does not assert a positive J.
LMD-T3 and T5 put the vacuum in both operator domains. Define the quadratic
form tensor by vectors, with inner products conjugate-linear in the first slot:

    Q_ab = <(H_a-<H_a>)psi, (H_b-<H_b>)psi>, H=(R,J).

No assertion that psi lies in the domains of all products H_a H_b is needed.
This is the same covariance/Gram construction used by the finite-cut adapter,
now for declared continuum operators. Positivity is the ordinary Gram
identity; it is not a new positivity theorem. The tensor describes these
response directions. It is not a spacetime metric or Riemann tensor.

## 2. Derivation of all four entries

Both means vanish. E J=integral w=0. The real half-density representation
makes R1 purely imaginary, while self-adjointness makes <1,R1> real; it is
therefore zero.

For this v, LMD equation (17) gives

    ||b_C||^2=1,             ||S_v||_HS^2=1/8.

Indeed |v_1|^2=1/(16 pi^2); substitution into the two stated Fourier sums
gives the factors above. LMD equation (16) therefore gives

    A(h)=Q_RR=h(1+h/2)^2/4+h^2/16
             =h/4+5h^2/16+h^3/16.                     (1)

LMD's Wick second-moment formula gives

    B(h)=Q_JJ
        =integral_0^1 integral_0^1
          cos(2 pi z) cos(2 pi y) exp[h K(z-y)] dz dy.  (2)

The subtraction of the squared mean is zero. This is an analytic integral,
not a computed decimal or an exactly admitted calculator coefficient.
It is finite because the singularity has exponent h<1.

For the off-diagonal entry use the corrected current identity, not an
unjustified product-domain commutator:

    (i/h)[R,J]=J_h(v w').

Let P_L be the spectral projection of total mass M(S1) onto [0,L].
LMD proves that it reduces R, while J and J_h(vw') are bounded there.
Hence the displayed commutator holds on Dom R in this subspace. Put
psi_L=P_L 1. These vectors converge to 1 in the graph norms of R and J:
P_L commutes with R and J, and the two vacuum graph norms are finite.
The right-side expectation converges as well, by domination by total mass.

Thus the limit of the bounded-mass weak commutator gives

    2 i Im Q_RJ = -i h E J_h(v w')
                = -i h integral v w'.

Here integral v w'=-1/2. Since R1 is purely imaginary and J1 is real,
Re Q_RJ=0. Consequently

    Q(h) = [ A(h),    i h/4 ]
           [ -i h/4, B(h)  ].                         (3)

The sign follows the stated generator order (R,J) and LMD's convention.
Reversing the order reverses the antisymmetric component.

## 3. A strict finite-h distinction from the one-qubit limit

The first Wiener-chaos component of M(w) has squared norm

    h <w,Cw> = h/4,

because w is the n=1 cosine and lambda_1=1/2. Orthogonality of Wiener
chaoses at finite cutoff and L2 convergence in LMD-T3 imply

    B(h) >= h/4.

Together with (1), this gives the explicit strict bound

    det Q=A(h)B(h)-h^2/16
         >= h^3(5+h)/64 > 0,    0<h<1.               (4)

This is a consequence of inherited Gaussian moments and the selected
observables, not a universal lower bound for arbitrary gravity models.
Both real tangent rank and complex Gram rank are 2 at finite h.

Expansion under the integral in (2), justified by the integrable logarithmic
kernel bound for 0<h<=h0<1, gives B(h)=h/4+O(h^2). Therefore

    Q(h)/h -> (1/4) [ 1,  i ]
                     [ -i, 1 ]     as h -> 0.        (5)

The limiting complex rank is 1, while its real metric still has rank 2.
Rank need not be continuous in this limit.

For a pure state on a single qubit, all centered vectors
(H_a-<H_a>)psi belong to the one-dimensional complex space psi-perp.
Their Gram rank is at most 1, for any Hermitian controls. Thus a single
pure qubit cannot realize (3) exactly for any 0<h<1, even though it
realizes (5). Mixed-state or enlarged-system encodings are outside that
specific one-pure-qubit obstruction.

## 4. Exact finite moment representation, and the calculator boundary

The existing NRQC zero-state geometry query with H=(X,Y) gives

    Q_0 = [1, i; -i, 1].

Multiplying both controls by sqrt(h)/2 matches the leading tensor h Q_0/4,
not the finite-h tensor. The rescaling is an analytic comparison; the
current geometry JSON accepts signed Pauli words, not arbitrary scalars.
The request below is a source-derived fixture, NOT a newly executed run:

    {
      "schema": "nrqc-calculator-v1",
      "n_qubits": 1,
      "operations": [],
      "backend": "stabilizer",
      "queries": [{
        "id": "semiclassical_current_density_shape",
        "kind": "geometry",
        "generators": [
          {"paulis": {"0": "X"}},
          {"paulis": {"0": "Y"}}
        ]
      }]
    }

There is also an explicit exact covariance representation in dimension 3.
Let c=h/4, alpha=sqrt(A), beta=c/sqrt(A),
gamma=sqrt(B-c^2/A)>0, and choose orthonormal |0>,|1>,|2>. Put

    H_R=alpha (|1><0|+|0><1|),
    H_J=i beta |1><0|-i beta |0><1|
          +gamma (|2><0|+|0><2|).

In state |0>, both means vanish and the two centered vectors are
alpha|1> and i beta|1>+gamma|2>. Their Gram matrix is exactly (3).
Embed these three levels in two qubits, leaving |11> unused.
This is elementary Gram realization of two moments. It is NOT an
intertwiner of LMD's full operator algebra, spectra, flows or dynamics.

With basis |q1 q0>, |0>=|00>, |1>=|01>, |2>=|10>, the controls are

    H_R=(alpha/2)(I_1+Z_1) X_0,
    H_J=(beta/2)(I_1+Z_1) Y_0
           +(gamma/2) X_1 (I_0+Z_0).                 (6)

NRQC can evaluate the constituent Pauli correlations. But its present
geometry API does not accept these weighted sums directly. The coefficients
contain B(h); their exact membership in the admitted coefficient field
is not established for general h. A future implementation must distinguish
symbolic coefficients from certified interval approximations. It must not
silently round them or label the continuum integral as an exact gate result.

Feeding A and B into (6) merely reproduces supplied moments. It does not
calculate the missing continuum answer or demonstrate computational advantage.

## 5. Which research obstruction remains?

The original MPQ inverse uses

    j=epsilon exp(s), y=d/j, theta=Y^{-1}(y),
    m=j/a(theta), P=integral m, r=theta+log|m/P|.

LMD replaces j dz by a random measure and s by a distribution. The pointwise
division d/j and the subsequent nonlinear inverse have not thereby been
defined. The response tensor (3) does not supply those operations.
The frozen 3D scalar shear obstruction also remains unchanged.

This calculation provides a source-specific diagnostic for a future
finite representation: matching only the leading one-qubit tensor omits
strict finite-h fluctuations. It does not close the original
geometric-observable reconstruction, independent unreduced ordering,
3+1 constraints, a continuum limit, or physical gravity validation.

Next implementation gate: select a bounded or smeared geometric observable
from the actual LMD inverse problem, construct its finite representation,
and bound the state, observable and cutoff errors separately. Only after
that map exists should NRQC and the existing native calculation be timed
on the SAME observable, state, cutoff, tolerance and preprocessing costs.
For this two-moment analytic target the existing formulas are already more
direct; no NRQC speed or memory advantage is established.

## 6. Review and negative controls

Written checks performed for this derivation:
- n=1 Fourier normalization in both LMD sums;
- off-diagonal sign from the corrected mass-subspace commutator;
- graph-norm passage for the unbounded operators;
- Wick first-chaos lower bound and strict determinant bound;
- h -> 0 asymptotics, real/complex rank distinction;
- direct action of (6) on |00> and preservation of the imaginary cross term.

These are analytic checks, not a fresh test-suite pass or independent review.
Negative comparisons to retain in an eventual implementation:
uncorrected a=1 is a different density-flow law; dropping the imaginary
entry loses the signed response; replacing B by h/4 omits higher chaoses;
a one-pure-qubit exact finite-h claim contradicts (4); matching (6) alone
does not simulate LMD dynamics. No numeric execution or CI dispatch was
performed for this note. Frozen predecessor packets remain untouched.
