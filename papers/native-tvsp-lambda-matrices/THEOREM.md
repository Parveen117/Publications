# TLM-1 — TVSP-generated layers in every matrix entry

16 October 2026. Parent: `35d63e4ac0d8cf161c549654ce6b430d68afbde4`.

## Construction and source contract

An observer's scalar reading is enriched by the local TVSP response that
generated it. Layers are derivatives of that response, with the tensor slots
and connection retained. They are not independently assigned interior matrices.
This packet supersedes the unpublished generic directed-dash recursion; that
recursion is not a dependency or part of the present result.

The source is NT-1/3/4/8 in `native-thermodynamic-curvature/THEOREM.md`, with
the tower hypotheses corrected by NJC-T7. All older files remain frozen.
Jet algebras and their regular representations are established mathematics.
The development here is their explicit TVSP source contract, native adapter,
curvature readout and exact equal-value/different-geometry example.

Fix a normalized, affine, oriented two-coordinate chart `(s,v)`, a point p,
positive S and V, and a supplied smooth energy u. Put

    T=u_s, P=-u_v, H=D²u=[[a,b],[b,c]], Delta=ac-b²,
    a>0, Delta>0,
    Lambda=(lambdaP,lambdaV,lambdaS,zT,zS,lambdaT)
          =(-S Delta/c, -S a, -V c, -V Delta/a, P/c, P a/Delta).       (1)

S and V are calibrated affine coordinates, not confused with their translated
increments s,v. These are NT's six normalized responses; they do not silently
replace the older uniform entropy-scaled four-response convention.

The context C includes this source germ, units, affine chart, observer point,
orientation, and the selected strain connection

    A_i = (1/2) H^(-1) partial_i H.                                  (2)

It also specifies a flat connection on the six normalized scalar response
components E and the Hessian tangent connection (2). The flat E connection
is an allowed choice, not something forced by six component labels.

For a scalar field f, a contextual entry is

    mathfrak f_N = [ C_N | j_p^N f ],
    C_N = (declared context, j_p^N H, j_p^N Lambda),
    O_p(mathfrak f_N)=f(p).                                         (3)

The source germ is a mathematical specification; its finite displayed record
contains the indicated jets. A numerical implementation can generate the record
from that specification, but cannot infer an unknown germ from one reading.
Context can be shared once by a matrix instead of duplicated in every cell.
Distinct local contexts must be explicitly carried to a common comparison
chart before the operations below; there is no implicit multiplication of
geometries belonging to unrelated sources.

We distinguish the response section Lambda, an individual response lambdaX,
the derivative order N, and the scalar family parameter mu used below. N=5,
f(p)=5 and mu=5 are three different assertions. No equality between them is
built into the notation.

## TLM-T1. Source-generated entry algebra and ordinary recovery

In one fixed context, let f range over scalar fields obtained from the source
responses, explicitly admitted coordinates, and regular scalar operations.
The enclosing algebra of N-jets is

    J_N = C[epsilon_s,epsilon_v]/(monomials of total degree N+1),
    j_p^N f = sum_(|alpha|<=N) (partial^alpha f(p)/alpha!) epsilon^alpha.

The epsilons are central formal increments. The adapter uses the native exact
cut field Q+iota Q as the coefficient subfield. This is a finite local record,
not a convergence assertion for a Taylor series.

For compatible jets x,y define

    (x star y)_alpha = sum_(beta+gamma=alpha) x_beta y_gamma,          (4)
    (x+y)_alpha=x_alpha+y_alpha,    (x*)_alpha=conjugate(x_alpha).

Then J_N is an associative unital algebra with involution. Taylor lifting
preserves products, sums and conjugation. The degree-zero observer pi_N is
a unital algebra homomorphism, and truncation to lower depth is likewise a
homomorphism. A scalar jet is invertible iff its constant coefficient is
nonzero; its inverse is the truncated jet of the reciprocal germ.

Proof. The multiindex Leibniz rule, divided by alpha!, is exactly (4).
Either association of three factors sums the ordered product
x_beta y_gamma z_eta over beta+gamma+eta=alpha. Terms above N cannot feed
back into smaller degree. The constant coefficient gives the observer and
identity statements; conjugation reverses products after matrix amplification.
For x=x_0(1+h), h has positive degree and h^(N+1)=0, so
x^(-1)=x_0^(-1) sum_(k=0)^N (-h)^k. The converse follows by taking degree zero.

Source admissibility is separate from closure of this enclosing algebra.
Multiplying Hessian matrices need not produce a symmetric positive Hessian of
another energy. No theorem declares every operator product a new thermodynamic
state, nor every arbitrary list of layer coefficients a source-derived record.

Differentiation has the typed rule

    partial_i : J_(N+1) -> J_N,
    (partial_i x)_alpha=(alpha_i+1)x_(alpha+e_i).                     (5)

It is well defined because derivatives of a function vanishing through order
N+1 vanish through order N. It is not an endomorphism derivation on J_N:
epsilon_s^(N+1)=0 there, whereas its formal derivative
(N+1)epsilon_s^N is nonzero. A displayed derivative consumes one retained layer.

## TLM-T2. One finite matrix retains all layers

For an m-by-n matrix M of contextual entries in one aligned context,

    (M star N)_ij=sum_k M_ik star N_kj,
    (M*)_ij=(M_ji)*,
    Pi(M)= [pi_N(M_ij)].                                            (6)

Composition is associative, dagger reverses products, and

    Pi(M star N)=Pi(M)Pi(N),    Pi(M*)=Pi(M)*.                        (7)

There is also a faithful single-matrix encoding. Order monomials by total
degree, with s before v at degree one. Their number is d=(N+1)(N+2)/2.
For M in M_(m,n)(J_N), define the block matrix

    Reg_N(M)_(alpha,beta) = M_(alpha-beta) if alpha>=beta componentwise,
                            0 otherwise.                           (8)

The blocks are m-by-n. It represents multiplication on J_N-valued columns,
so Reg(M star N)=Reg(M)Reg(N). The constant-input block column consists of
all M_alpha and is an explicit decoder. Thus the encoding is faithful.
For a square matrix its diagonal blocks all equal M_0, hence

    det(t I-Reg_N(M)) = det(t I-M_0)^d.                              (9)

Proof. Matrix multiplication in (8) sums M_(alpha-gamma)N_(gamma-beta)
over beta<=gamma<=alpha, exactly (4). The zero input monomial proves
faithfulness. Increasing degree makes (8) block lower triangular, proving (9).

This gives both visibility and a warning about a specific observer: matrix
entries/action recover the layers, but this characteristic polynomial cannot
see them. Equal spectra are not equal layered records or equal geometry.

Matched source unitary frames obey, at every retained depth,

    M_(out<-in)=G_out* star M star G_in,
    N_(next<-out) star M_(out<-in)=G_next* star (N star M) star G_in.  (10)

This follows by inserting G_out star G_out*=I. Actual smooth unitary germs
supply such jets; independent changes to individual entries need not define
one composable frame change. Nonlinear chart changes require the actual jet
chain rule and tensor/connection transport. A bare Hessian recomputed after
a nonlinear coordinate change is not automatically the same metric.

The coefficientwise involution in J_N must not be confused with ordinary
conjugate transpose of Reg_N. For N>=1 a formal coordinate is a nonzero
self-adjoint nilpotent. In a positive bounded-operator *-representation its
image would be a self-adjoint nilpotent, hence zero. Both coordinates, and
thus the whole positive-degree ideal, are killed. Consequently there is no
faithful positive bounded-operator *-representation of the full jet algebra.
Reg_N is an algebra encoding, not a bank of new self-adjoint quantum operators.
Positivity of the actual H(s,v) on its regular source patch remains valid.

## TLM-T3. Relation to the existing lambda response tower

The existing covariant tower is retained with its actual tensor type:

    L_1=Lambda,  L_(k+1)=nabla L_k,
    L_(k+1) in E tensor (T*)^(tensor k).                             (11)

In the chosen flat E frame and Hessian tangent connection, its first layers are

    (L_2)^A_i = partial_i Lambda^A,
    (L_3)^A_ij = partial_i partial_j Lambda^A
                  - Gamma^r_ij partial_r Lambda^A.                 (12)

More generally differentiate a tensor component and subtract one Christoffel
term on every lower tangent slot. Thus j^N Lambda and j^N H reconstruct all
L_1 through L_(N+1). For N=0 only L_1 is claimed. For k>=2, order k uses
derivatives of H through at most k-1 and derivatives of Lambda through k.

Proof. H^(-1) and its derivatives are determined recursively on Delta>0;
Gamma=(1/2)H^(-1)dH has order one in H. Induction in the stated tensor product
rule proves the order bound and reconstruction. NT-8 then gives, on any layer T,

    ([nabla_i,nabla_j]T)(v_1,...,v_k)
       = -sum_r T(v_1,...,F^H_ij v_r,...,v_k),                     (13)

since F^E=0 in this chosen fibre frame. Other source choices retain the
additional F^E T term. Higher ordered covariant slots need not be symmetric;
their curvature defect is recovered using connection data, not discarded.

This is how the finite scalar jet storage interfaces with the actual TVSP
tensor tower. Layer k is not Lambda^k, a new scalar angle, or a physical
spatial cell. A cut grading may be inherited only when its specified base
and fibre actions preserve the source/connection hypotheses of NJC-T7.
No arbitrary potential is declared cut-equivariant just by drawing a compass.

## TLM-T4. Minimal context layer for the selected curvature

Put X_i=H^(-1)partial_i H. NT-3 proves for (2)

    F^H_sv = -[X_s,X_v]/4,
    kappa = (H F^H_sv)_12 / Delta.                                 (14)

Indices in the numerator are one-based. In the Hessian tangent adapter this
is Gaussian curvature, with NT's orientation convention. Formula (14)
depends only on j^1H at the observer point. Thus a first context layer
suffices; order zero is insufficient over the admitted source class, by T5.
This is minimal universal derivative order relative to an H-jet readout,
with the connection selection fixed and no independent disclosure of the
source parameter. It is not a claim that each component of the first jet
is independently necessary.

For completeness, the cancellation behind (14) is

    partial_s X_v-partial_v X_s=-[X_s,X_v],
    d(X/2)+(X/2) wedge (X/2)=-[X_s,X_v]/4.

It follows from differentiating H^(-1)H=I. H=D²u additionally makes the
selected compatible connection torsion free, identifying it with the
Hessian metric connection. This is a reuse of NT, not a new curvature law.

The connection selection is indispensable. For example the same constant
H=I admits the metric-compatible connections A=0 and A_s=0,A_v=sR,
R=[[0,-1],[1,0]], with curvature 0 and R. The latter is not the selected
torsion-free Hessian connection. Even every H layer cannot choose between
independently admitted transport laws. Context must name the law it uses.

## TLM-T5. An actual visible 5 with source-generated curved layers

Use S=1+s,V=1+v and the parameterized NT-4 energy

    u_mu=10s-10v+s²+sv+(5/2)v²+(mu/2)s²v,
    H_mu=[[2+mu v,1+mu s],[1+mu s,5]],
    a=2+mu v, b=1+mu s, Delta=5a-b².                               (15)

The common chart |s|,|v|<1/20, |mu|<=5 is stable: a>=7/4 and
Delta>=35/4-25/16=115/16>0. The source's units and coordinates have already
been normalized; mu is an independent dimensionless deformation control.
It is not equal by definition to any component of Lambda.

All six response values at p=(0,0) are, for every mu,

    Lambda(p)=(-9/5,-2,-5,-9/2,2,20/9),
    H_mu(p)=[[2,1],[1,5]].                                         (16)

The positive mechanical reading f=-lambdaS=Vc=5(1+v) has f(p)=5.
Its entire scalar germ, including all its derivatives, is identical for
every mu. The source context is nevertheless different:

    partial_s H=[[0,mu],[mu,0]],
    partial_v H=[[mu,0],[0,0]].                                     (17)

Explicit multiplication in (14) gives throughout the chart

    F_sv = mu²/(4 Delta²) [[5b,25],[-5a,-5b]],
    kappa = 5mu²/(4 Delta²),
    kappa(p)=5mu²/324.                                              (18)

In particular:

| Deformation mu | Observer f(p) | Source curvature kappa(p) |
|---|---|---|
| 0 | 5 | 0 |
| 1 | 5 | 5/324 |
| 5 | 5 | 125/324 |
| -5 | 5 | 125/324 |

Proof. Directly differentiate (15), substitute the six definitions (1),
invert H with determinant Delta and use (14). The inverse is
Delta^(-1)[[5,-b],[-b,a]]. The commutator yields (18), including its sign.
In this family partial_v lambdaV(p)=-mu; hence the resolved source reading
recovers mu=-partial_v lambdaV(p), and

    K([C_1|j^1 f]) = (5/324)(partial_v lambdaV(p))².                 (19)

Equation (19) is calibrated for this family at this observer point. General
source records use (14). Mu and -mu have equal curvature but different first
layers: the layer retains orientation information that this scalar curvature
forgets. Flatness at mu=0 follows from constant H in (15). It does not mean
Lambda=0, f=0, or that all higher response layers vanish.

The first-layer response matrix can literally be stored as one matrix:

    Reg_1(H_mu) =
      [ H0   0    0 ]
      [ Hs   H0   0 ] ,
      [ Hv   0    H0 ]
    H0=[[2,1],[1,5]], Hs=[[0,mu],[mu,0]], Hv=[[mu,0],[0,0]].          (20)

It is a six-by-six matrix with the same characteristic polynomial
(t²-7t+9)^3 for every mu. Its lower blocks distinguish every mu and recover
(18). Matrix dimension six is storage size, not physical spacetime dimension.

There is no universal function of f(p), or even of the entire germ of this f,
which gives the selected source curvature: mu=0 and mu=1 have identical f
but different kappa. This proves the need to retain additional source response
information. The added context is not extracted for free from a numeral.

## Evidence and continuation

Written proofs establish the arbitrary finite-depth statements. The native
certificate checks source differentiation, rational response layers, ordered
matrix products, exact ordinary recovery, one-matrix decoding, differentiated
connection curvature, and failures of stronger interpretations. Python 3.12
does I/O only; unchanged native F/Cut and Paninian engines do the arithmetic
and symbolic reductions. Source pins distinguish this adapter from a second
shared operator engine.

No physical space quantization, intrinsic geometry of every bare scalar,
continuum quantum representation, new gravitational constraints, or 3+1
quantum-gravity completion is inferred. The present advance is a usable
source-generated matrix record and an exact curvature-sensitive reading.
Its next physical use requires an actual family of geometry observables,
admissible preparation/measurement and domains realizing the retained layers.
