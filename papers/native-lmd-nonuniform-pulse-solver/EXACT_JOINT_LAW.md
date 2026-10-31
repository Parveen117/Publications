# Post-run analytic extension: total mass and its first spatial harmonic

Author: Monty. 31 October 2026.

Status: written derivation under a newly credited external theorem. This
extension is not used by the frozen first-mode numerical campaign. Its
positive-cell solver remains to be implemented and qualified separately.

## 1. External input and exact convention match

Reda Chhaibi and Joseph Najnudel, *On the circle, Gaussian Multiplicative
Chaos and Beta Ensembles match exactly*, [arXiv:1904.00578v3](https://arxiv.org/pdf/1904.00578v3),
20 February 2025, Theorem 2.1 (pages 6-7), identifies circular GMC as a
random probability measure times a dependent mass factor. Its independent
Verblunsky coefficients a_j are rotationally invariant, with
|a_j|^2~Beta(1,B(j+1)/2); for B>2 the mass factor is

    W=product_(j>=0) (1-2/[B(j+1)])/(1-|a_j|^2).

The covariance and Wick convention are in equations (1.2.1)-(1.2.2),
pages 4-5. This measure-valued theorem is the new analytic input, not a
consequence of matching a finite list of moments. It is credited to those
authors; its proof is not supplied by the present numerical tests.

For this LMD application put alpha=h/2, gamma=sqrt(alpha), B=2/alpha=4/h.
Their real field G has covariance -2 log|chord| and coefficient gamma, so
gamma G has precisely the LMD nonzero-mode covariance -h log|chord|.
Their measure is normalized by dtheta/(2 pi) and Wick-subtracted by half
the exponent variance. LMD-T3 identifies this regularization with its own
Fourier limit. The independent LMD zero mode still contributes
Q=exp(q-h/2), q~N(0,h). No additional smooth covariance or radius factor
is inserted. The working range 0<h<1/2 is strictly subcritical here.

## 2. Isolate one spectral coefficient without losing its mass dependence

Let U_j=1-|a_j|^2. Then U_j~Beta((j+1)/alpha,1), independently. The degree-one
orthogonal polynomial is z-conj(a_0), so its orthogonality to constants
gives integral z mu(dz)=conj(a_0). Only its phase convention changes under
conjugation. Consequently, for the *same* source measure,

    J=integral M_h(dz)=QW,
    L=integral exp(2 pi i z) M_h(dz)=QW conj(a_0),
    W=(1-alpha)U_0^(-1) T,
    T=product_(j>=1) (1-alpha/(j+1)) U_j^(-1).          (1)

The tail T is independent of U_0 and the uniform phase of a_0. Total mass
W itself is not independent of U_0. Keeping that dependence is essential.

The tail has an explicit one-variable law. For a finite product T_m over
n=j+1=2,...,m, elementary Beta integration gives

    E T_m^p=product_(n=2)^m (1-alpha/n)^p/(1-alpha p/n).

For p=i y this is the characteristic function of log T_m. Gamma products
and the Gamma-ratio asymptotic give its limit

    E exp(i y log T)=Gamma(2-i alpha y)/Gamma(2-alpha)^(i y).

Here log T_m converges almost surely: the centered independent log terms
have summable variances, and their mean corrections are O(n^(-2)). Since
|T_m^(i y)|=1, dominated convergence justifies the characteristic limit.
Uniqueness of characteristic functions, not a finite-moment inference,
therefore identifies

    T =_law G^(-alpha)/Gamma(2-alpha),   G~Gamma(2,1). (2)

The Gamma-ratio step uses the standard leading ratio asymptotic, e.g.
[DLMF 5.11.13](https://dlmf.nist.gov/5.11.E13). Explicitly, the residual
ratio is Gamma(m+1-alpha)^p Gamma(m+1)^(1-p)/Gamma(m+1-alpha p), tending to one.
Independence in (1) allows G to be taken independent of U_0 and the phase.

## 3. Exact finite-variable joint source representation

Write U_0=V^alpha, where V is uniform on (0,1), and let Phi be independent
uniform angle. With q,G,V,Phi independent, equations (1)-(2) give

    J =_law Q (GV)^(-alpha)/Gamma(1-alpha),
    L =_joint-law J sqrt(1-V^alpha) exp(i Phi).         (3)

The joint-law notation means the same J and V occur on both right-hand
sides. This is a representation of (J,L), not a replacement for every
field observable or a physical preparation of those four variables.

For the first-harmonic source family, (3) gives the exact weighted law

    J_(e,beta) =_law J[1+e sqrt(1-V^alpha)cos(Phi-beta)]. (4)

It includes all spatial modes through the imported joint measure law.
It is not the first Fourier cutoff of the Gaussian field. The angle may
then be integrated using the ordinary Bessel function C=J_0:

    chi_e(a)=E_(q,G,V) exp(i a J)
                  C(a e J sqrt(1-V^alpha)).             (5)

As a uniform check, the transformation X=GV,Y=G(1-V) has joint density
exp(-X-Y) on X,Y>0: the Jacobian cancels the shape-two Gamma factor G.
Thus X~Exp(1), and setting e=0 recovers the preceding full mass law exactly.

## 4. Moment and dependence checks, proved analytically

Direct integration in G and V gives

    E J=1,
    E J^2=I(h),
    E |L|^2=I(h) alpha/(1-alpha)=I(h)h/(2-h),
    E J_e^2=I(h)[1+(e^2/2)alpha/(1-alpha)].             (6)

The uniform phase makes E L=E(JL)=0 and gives the same variance to its real
and imaginary parts. More generally the source pair has

    E[J_(e,beta) J_(f,delta)]
      =I(h)[1+(ef/2)alpha/(1-alpha)cos(beta-delta)].    (7)

This is a joint *source-mass* identity. It is not a simultaneous sharp
reading of the conjugate detector or a joint geometric POVM.

A useful analytic negative control is to draw J independently of the
normalized harmonic radius R=sqrt(1-V^alpha). Its unweighted E R^2 is
alpha/(1+alpha), while the actual mass-squared weighted ratio in (6) is
alpha/(1-alpha). The independent shortcut therefore gives the wrong
weighted second moment. The earlier uniform Exp(1) sample cannot simply
be paired with an independent Beta-distributed shape.

## 5. Concrete next numerical contract

Factor J=Q A(G) B_0(V) and |L|=Q A(G) B_1(V), where

    A(g)=g^(-alpha)/Gamma(2-alpha),
    B_0(v)=(1-alpha)v^(-alpha),
    B_1(v)=B_0(v)sqrt(1-v^alpha).

Both A and B_0 have mean one. Their complete second moments are

    E A^2=Gamma(2-2alpha)/Gamma(2-alpha)^2,
    E B_0^2=(1-alpha)^2/(1-2alpha),
    E B_1^2=alpha(1-alpha)/(1-2alpha).

Gamma(2,1) probability cells and the A first masses are Gamma CDF differences
of shapes 2 and 2-alpha. Uniform V cells have probability v-u and B_0 first
mass v^(1-alpha)-u^(1-alpha). The remaining positive first mass is

    b_1(u,v)=(1-alpha) integral_u^v x^(-alpha)sqrt(1-x^alpha) dx.

Its certified evaluation, including the singular zero endpoint, is the new
special-function obligation. Conditional averaging on q,G,V cells while
retaining Phi then captures the second moment

    S_q S_A [S_B0+(e^2/2)S_B1].                        (8)

Subtract (8) from the full moment (6) and use the same conditional Taylor
bound. There is no omitted Fourier-mode term under this new representation.
L2 conditional convergence follows from the finite moments above and
exhausting partitions; fixed undivided tails are not an arbitrary-accuracy
claim. The numerical first masses, finite error budget, runtime and target
probabilities have not been executed for this representation in this packet.

Next execution should freeze those cells and precision, recover e=0,
verify (6), reject the independent-shape shortcut, and overlap the retained
first-mode certificates. This is the next implementation target. The
existing campaign and its unresolved tolerance result remain untouched.
