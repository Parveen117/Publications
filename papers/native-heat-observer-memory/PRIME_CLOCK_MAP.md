# HCM-3 — Coprime gain recovery, prime audits and rational reference correction

15 October 2026. Extends the actual HCM/YM54 source; no new heat law is chosen.

## 1. The source gives a primitive exponent pair

The two spectral transport gains, not the raw readings q and h, are

    A(t)=exp(-5t/72),  B(t)=exp(-5t/162).

Set u(t)=exp(-5t/648). Then A=u^9 and B=u^4, with 0<u<=1 for finite heat time
t>=0. These powers are coprime: gcd(9,4)=1 and 9-2*4=1. The important property
is coprimality, not that 9 or 4 is prime (neither is). Prime numbers themselves
have exactly two positive divisors, 1 and themselves.

The exact map and its inverse are

    u -> (u^9,u^4),       u=A/B^2,
    A^4=B^9.

This uses standard Bezout/monomial algebra; no novelty for that general method
is claimed. Its application here identifies a root-free common parameter of
the already published heat source.

## 2. Root-free compatibility theorem

For nonzero complex A,B, the condition A^4=B^9 is necessary and sufficient for
a unique nonzero u satisfying A=u^9 and B=u^4. The inverse is u=A/B^2.

Proof. Necessity is immediate. If A^4=B^9, then (A/B^2)^4=B, and
(A/B^2)^9=A*A^8/B^18=A. Uniqueness follows from u^9/u^8=u.
For positive rational A,B the recovered u is positive rational. Restricting
u to (0,1] is an additional heat-time admissibility check; an arbitrary
compatible complex pair is not an actual dissipative heat observation.

More generally, for gcd(m,n)=1 and integers r,s with rm+sn=1, the map
u -> (u^m,u^n) on nonzero complex values has inverse u=A^r B^s on its
compatibility locus A^n=B^m. Negative powers require nonzero inputs.
If gcd(m,n)=d>1, every d-th root of unity gives an ambiguity u -> omega u.
The primitive exponent condition is exactly what removes that ambiguity.

Zero is excluded from the inverse. It occurs only in the infinite-time limit
of this heat law; declaring its lost phase or scale recoverable would be false.

## 3. What matching primes contributes

A positive rational has a unique finite prime-valuation record

    u=product_p p^(k_p),    k_p in Z.

Therefore the source gain pair has valuations

    v_p(A)=9k_p,   v_p(B)=4k_p.

Conversely, rational gains are compatible exactly when
4v_p(A)=9v_p(B) for every prime. Coprimality then gives
k_p=v_p(A)-2v_p(B). Matching only the set of primes is insufficient; their
integer multiplicities must match. For example 1/2 and 1/4 have the same prime
support but fail A^4=B^9.

The valuations are an exact audit representation, not the required recovery
algorithm: A/B^2 already reconstructs u without factorization. The certificate
factors small declared rational examples and checks both routes. This does not
claim efficient factorization of arbitrary large integers or a general exact
representation of irrational/transcendental observations by finite prime lists.

## 4. The corresponding phase statement and its limit

On the unit circle, A=u^9 and B=u^4 give

    arg(A)=9 arg(u),  arg(B)=4 arg(u) modulo 2pi,
    arg(u)=arg(A)-2 arg(B) modulo 2pi.

Thus the paired phase record determines the common phase without a root-branch
ambiguity. Either channel alone may be ambiguous. This recovers an endpoint
phase, not an integer path-winding ledger; distinct lifted paths can have the
same endpoints.

In the actual heat source u>0, so this common parameter's phase is zero.
The changing phase of z=q+i h from HCM-2 is a different quantity, arising from
the relative channel mixture. The complex monomial theorem is an algebraic
extension, not a derived phase law for that measured z.

One may assign an angle theta_p to each prime and define a multiplicative
phase readout exp(i sum k_p theta_p). Such a choice is extra structure, not a
canonical prime-to-phase law. It can discard information: assigning both 2
and 3 the quarter-turn i identifies them. Any root-of-unity assignment also
annihilates some nontrivial prime powers. A finite phase label cannot replace
the exact valuation record for unrestricted positive rationals.

## 5. Rational observation schedules and premeasurement correction

In HCM's retained frame let Pslow and Pfast be the fixed complementary EMK
spectral projectors. Its exact transfer at parameter u is

    T(u)=u^4 Pslow+u^9 Pfast,
    T(u)^(-1)=u^(-4)Pslow+u^(-9)Pfast.

Both matrices are rational whenever 0<u<=1 is rational. For u=1/2, the
mathematical clock is t*=(648/5)log 2, and

    T(1/2)=[[521,372],[372,304]]/12800,
    T(1/2)^(-1)=[[4864,-5952],[-5952,8336]]/25.

If x(t*)=(q(t*),h(t*)), the corrected functional is

    q(0)=(4864/25)q(t*)-(5952/25)h(t*).

This may be implemented as a calibrated mixed readout before recording, if the
instrument can access that linear functional, or after acquiring both channels.
It cannot reconstruct arbitrary hidden-channel information from q(t*) alone.
The amplitude/gain calibration and observer accessibility are additional inputs;
the desired answer cannot be used to fit them.

Choosing a rational u is a mathematical sampling contract, not an assertion
that an irrational physical time can be set with infinite precision. These
are declared heat ticks, not calibrated seconds. Exactness of the rational
certificate and uncertainty of a physical implementation are separate.

## 6. Correction and clock-error budgets

For the retained Euclidean pairing (equivalently its common reference Gram
scale), the inverse operator norm is u^(-9). If the correct calibrated transfer
is used, additive measurement error e produces reconstruction error at most
u^(-9)||e||. At u=1/2 the worst-case amplification is 512, attained on Pfast.
Thus algebraic recoverability does not imply numerically stable recovery late
in the heat process.

Suppose the actual and assumed parameters are u and uhat in (0,1], with
|log(uhat/u)|<=epsilon. The spectral decomposition gives

    ||T(uhat)^(-1)T(u)-I|| <= exp(9epsilon)-1.

For observed y=T(u)x0+e and reconstructed xhat=T(uhat)^(-1)y,

    ||xhat-x0|| <= [exp(9epsilon)-1]||x0|| + uhat^(-9)||e||.

Proof: the two relative eigenvalues are exp(-4delta), exp(-9delta), where
delta=log(uhat/u); bound their deviations from 1 and then the noise term.
If calibrated gains Ahat,Bhat are used, the root-free estimator satisfies
log(uhat/u)=log(Ahat/A)-2log(Bhat/B). Hence absolute log-gain error bounds
epsilon_A,epsilon_B give epsilon_A+2epsilon_B. The root-free formula is an
exact algebraic identity, not a claim of optimal statistical conditioning.

## 7. Result and open boundaries

The source now has one recoverable common gain parameter, an optional prime
audit, a precise paired-phase theorem, rational transfer certificates and a
quantitative correction budget. Standard unique prime factorization and Bezout
identities are used explicitly. No relation to RH zeta zeros is proved by this
prime appearance. No full UGD numeral homomorphism, source-independent phase
law, recovery through total information loss or physical instrument design is
claimed. All source-specific calculations remain in the HCM free coefficient
block and its declared observer.
