"""Independent floating-point cross-check of the operator entries used by
verify.py, straight from the source formula (ASM (2)):

    q_L[g] = l(g)^2 + (Omega_0/2) N_g + int_0^inf k(u)(N_g - C_g(u)) du
             - sum_{n in N_L} a_n C_g(log n),

evaluated by adaptive quadrature on explicit polynomial profiles, with no
use of the Taylor kernel, the Legendre energy identity, the boundary
decomposition or any piecewise formula of the verifier.  It is a sanity
check of the analytic pipeline, not part of the certificate: it uses
binary floats and SciPy and is accurate to roughly 1e-9.

Usage: python3 crosscheck_quadrature.py L_numer L_denom N [matrices.json]
Prints the symmetric N x N matrix <psi_i, T psi_j> (unnormalized
Legendre P_(2j)(x/L), N <= 6) and the maximum deviation from the
reference block stored in crosscheck_reference.json for that width
(entries of the verifier's own --matrices exploration run), or from
the 'A' entries of a --matrices file given as a fourth argument.
"""
import sys, json, math
import numpy as np
from numpy.polynomial import legendre as Lg
from scipy.integrate import quad

Lnum, Lden, N = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
L = Lnum/Lden
gamma_E = 0.57721566490153286
Omega0 = -gamma_E - math.pi/2 - 3*math.log(2) - math.log(math.pi)

def prime_powers():
    out = []
    for n in range(2, 64):
        p = min(d for d in range(2, n+1) if n % d == 0)
        m = n
        while m % p == 0: m //= p
        if m == 1 and math.log(n) < 2*L: out.append((n, math.log(p)/math.sqrt(n)))
    return out
PRIMES = prime_powers()

def psi(j):
    c = np.zeros(2*j+1); c[-1] = 1.0
    return lambda x: Lg.legval(x/L, c) if abs(x) < L else 0.0

def k(u): return math.exp(-u/2)/(1-math.exp(-2*u))

def corr(f, g, u):
    # C_{f,g}(u) = int f(x) g(x-u) dx over the overlap of the supports
    lo, hi = max(-L, -L+u), min(L, L+u)
    if hi <= lo: return 0.0
    return quad(lambda x: f(x)*g(x-u), lo, hi, limit=200, epsabs=1e-13, epsrel=1e-13)[0]

def pair(f, g):
    # polarization of q_L: <f, T g>
    lf = quad(lambda x: f(x)*math.cosh(x/2), -L, L, limit=200, epsabs=1e-13)[0]
    lg = quad(lambda x: g(x)*math.cosh(x/2), -L, L, limit=200, epsabs=1e-13)[0]
    Nfg = quad(lambda x: f(x)*g(x), -L, L, limit=200, epsabs=1e-13)[0]
    def integrand(u):
        c = 0.5*(corr(f, g, u) + corr(g, f, u))
        return k(u)*(Nfg - c)
    # the kernel is integrable at 0 after the difference; split at 2L (no overlap beyond)
    gam = quad(integrand, 0, 2*L, limit=400, epsabs=1e-12)[0] + Nfg*quad(k, 2*L, np.inf, limit=200)[0]
    prime = sum(a*0.5*(corr(f, g, u) + corr(g, f, u)) for n, a in PRIMES for u in [math.log(n)])
    return lf*lg + Omega0/2*Nfg + gam - prime

A = np.zeros((N, N))
for i in range(N):
    for j in range(i, N):
        A[i, j] = A[j, i] = pair(psi(i), psi(j))
np.set_printoptions(precision=12, linewidth=160)
print('L =', Lnum, '/', Lden, ' active prime powers:', [n for n, _ in PRIMES])
print(A)
import os
if len(sys.argv) > 4:
    d = json.load(open(sys.argv[4]))
    ref = np.array([[float(x[0]) for x in row[:N]] for row in d['A'][:N]])
else:
    key = f'{Lnum}/{Lden}'
    refs = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'crosscheck_reference.json')))
    if key not in refs:
        print('no reference block for this width; add a --matrices file as fourth argument'); sys.exit(0)
    assert N <= 6
    ref = np.array([[float(x) for x in row[:N]] for row in refs[key]['A_block_6'][:N]])
print('max |A - A_ref| =', np.abs(A-ref).max(), ' max |A_ref| =', np.abs(ref).max())
