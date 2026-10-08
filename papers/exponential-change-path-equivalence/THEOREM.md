# ECP1 — Exponential change, ordered paths and the flat sector

8 October 2026. Written proofs and scoped exact computational certificate.

## Premises and relation to the programme

Work in a unital associative algebra of lawful endomorphisms over a characteristic-zero field. Exponentials first mean formal factorial series in central parameters. Finite matrices give convergent analytic models; arbitrary unbounded generators require independently established domains. Smooth vector fields below act on smooth functions by Lie differentiation and their local flows, not by an unrestricted Taylor series on every smooth function.

Sources: Publications `papers/emk-ugd-algebra/CONNECTION_CALCULUS.md` (EMK-C1), `papers/native-thermodynamic-curvature/THEOREM.md` (NT-2/3), `papers/generalized-euler-evolution/THEOREM.md` and its source audit (domains); Extra Ideas at `772eac6ca90fe508209c88518c32025f2300d726`, `uncut/up4/UP4_COST_OF_COMMUTING_CUTS.md` and `02-relational-response/NATIVE_CURVATURE_FLOW_R45.md`. These supply the native premise order and prior results. Lie brackets, commuting-flow criteria and connection curvature are established mathematics; no priority claim is made for them. ECP1 supplies a linked repair of the uploaded path-equivalence argument, formal word replay and negative controls.

## T1. The exponential retains ordering

Put U(t)=Exp(tA), V(s)=Exp(sB). In formal series,

    U(t)V(s)-V(s)U(t) = ts[A,B] + terms of total degree >=3.
    L(t,s)=U(t)V(s)U(-t)V(-s)=I+ts[A,B]+terms of total degree >=3.

Every remainder monomial contains both t and s. Multiplying the four factorial series proves this, without changing the order of a word. Moreover U(t)V(s)=V(s)U(t) as formal series iff [A,B]=0: necessity reads the ts coefficient; sufficiency reorders A^m B^n. This is not an identification of a finite loop with its leading term.

For finite matrices in any submultiplicative norm, let r=2(|t| ||A||+|s| ||B||). The total-degree tail after degree two in L is bounded by exp(r) r^3/6, since the product's degree-n coefficient norm sum is bounded by r^n/n!. This coarse bound is sufficient, not optimal. Thus the formal leading defect is also an analytic small-loop statement here.

## T2. Lie derivatives and response symmetry

Let X,Y be smooth vector fields on a finite-dimensional smooth manifold and D_X,D_Y their Lie derivatives on scalar functions. Leibniz and differentiation give

    [D_X,D_Y] = D_[X,Y].

For smooth f on a compact neighborhood admitting the local flows, two-variable Taylor expansion gives

    (U_X(t)U_Y(s)-U_Y(s)U_X(t))f
        = ts[D_X,D_Y]f + O(|t|^2|s|+|t||s|^2).

The remainder is locally uniform; derivatives on that neighborhood bound it. Only this finite-order expansion is used. The local flows commute on all functions iff [X,Y]=0. Necessity follows from the mixed derivative. For sufficiency, the flow of X preserves Y when [X,Y]=0 (differentiate its pushforward); it therefore maps integral curves of Y to integral curves with the same initial point. Local ODE uniqueness proves commutation.

For a smooth potential E, define T=D_X E, P=-D_Y E, B1=D_Y T and B2=-D_X P. Then

    B2-B1 = [D_X,D_Y]E.

If [D_X,D_Y]=alpha D_X+beta D_Y, the right side is alpha T-beta P. Reciprocal response is a zero-defect sector for that potential. It is weaker than operator commutation: one potential may be blind to a nonzero bracket.

Exact model: X=partial_x, Y=exp(a x)partial_y. Then [X,Y]=aY and

    U_X(t)U_Y(s)f=f(x+t,y+s exp(a(x+t))),
    U_Y(s)U_X(t)f=f(x+t,y+s exp(a x)).

The family has commuting flows at a=0 and generally order-dependent flows at a!=0. For E=y the response defect is a exp(a x); for E=x it vanishes even at a!=0.

## T3. Native connection curvature requires frame subtraction

Let derivations delta_i satisfy [delta_i,delta_j]=c_ij^k delta_k with central constant coefficients. On the regular module define nabla_i=delta_i+A_i (left multiplication). Leibniz expansion proves

    ([nabla_i,nabla_j]-c_ij^k nabla_k)f=F_ij f,
    F_ij=delta_i A_j-delta_j A_i+[A_i,A_j]-c_ij^k A_k.

Thus a noncommuting frame with A=0 has F=0. The flow bracket alone is not this coefficient curvature. In commuting coordinate directions the leading exponential loop of the covariant generators is F; in a noncommuting frame one must also compare/subtract the base-frame displacement. Orientation and transport signs must be declared.

For an invertible response H put X_i=H^-1 delta_i H. Differentiating H^-1 H=I and antisymmetrizing gives

    delta_i X_j-delta_j X_i-c_ij^k X_k=-[X_i,X_j].

For central constant q, substitution into the definition gives

    F^(qX)_ij=q(q-1)[X_i,X_j].

Hence q=1 is flat even for noncommuting X_i. The NT-3 self-dagger response and strain selection give q=1/2 and F=-[X_i,X_j]/4. These selection hypotheses are inputs; ECP1 does not derive H from the uploaded scalar y0. Symmetric point responses, closed finite walks and this connection flatness remain separate predicates.

## T4. Exact repair of the draft's affine equation

The draft's T(y)=y'+y0 is affine, not a linear derivation. Its two displayed calculations have difference y0'-y0''. They establish only

    path equality iff y0''=y0', hence y0=C1 exp(x)+C2

on a connected interval. Solving the scalar ODE does not construct a multidirectional connection or imply curvature. In particular the exponential solution is a zero-defect solution of this equality.

To represent affine operations faithfully, use pairs (y,c) and homogeneous operators

    T_hat(y,c)=(Dy+c y0,0), D_hat(y,c)=(Dy,0),
    S_hat(y,c)=(y+c y0,c).

Then T_hat^2(y,1)=(D^2y+D y0,0), whereas
T_hat D_hat S_hat(y,1)=(D^2y+D^2y0,0). Adding the same y0 readout reproduces the draft's two expressions exactly. These are words of different lengths, not an inverse four-step flow loop. This preserves the actual draft without pretending its defect is automatically NT-3 curvature.

The closed affine translation step y -> Dy+y0 itself can instead be homogenized as (Dy+c y0,c); its products have their own explicitly computable residues. Neither homogenization licenses a generator that changes the dimension of its carrier.

## Evidence boundary

General formal and smooth statements have the written proofs above. The executable verifier replays exact free-word coefficients through degree four, rational finite-matrix loops, the Lie-response witness, affine draft expressions, pure-gauge cancellation and rational-curved counterexample. Deliberate false replacements must be rejected. Finite controls do not certify all unbounded analytic domains, independently review the paper, prove universal physics, or supply a proof-assistant formalization. Flatness is not rationality: rational matrices already carry nonzero curvature. Exponential evolution alone is not a curvature law.

Run with Python 3.12 and sympy: `python3.12 verify.py`. The certificate hashes this theorem and the verifier; it is scoped computational evidence, not external validation.
