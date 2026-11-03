"""Exact chart controls; consume the pinned original thermo implementation."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import importlib.util
import json
import sys


def run():
    if sys.version_info[:2] not in ((3, 11), (3, 12)) or not __debug__:
        raise RuntimeError('Python 3.11/3.12 without -O required')
    source = Path(__file__).resolve().parents[1]/'thermo-compass-foundations/model.py'
    spec = importlib.util.spec_from_file_location('ct_thermo', source)
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    counts = dict(trajectory_fixtures=0, moving_principal_seams=0,
                  exact_isotropic_steps=0, clock_tail_ledgers=0, rejected_mutations=0)
    rows = []
    directions = [(Q(1), Q(0)), (Q(0), Q(1)), (Q(3, 5), Q(4, 5)),
                  (Q(4, 5), Q(-3, 5)), (Q(-1), Q(0)), (Q(-3, 5), Q(-4, 5))]
    kr = lambda i, j: Q(i == j)
    for nu, root, eps, e in product([Q(4, 3), Q(3, 2), Q(5, 3)],
                                  [Q(1, 4), Q(1, 2), Q(1), Q(2)],
                                  [Q(0), Q(1, 4), Q(1)], directions):
        r = root**nu.denominator
        h = root**int((nu-2)*nu.denominator)  # kappa=1/nu
        delta = nu-2
        x = tuple(r*u for u in e)
        H = m.matrix([[h*(kr(i,j)+delta*e[i]*e[j])+eps*kr(i,0)*kr(j,0)
                       for j in range(2)] for i in range(2)])
        m.require_stable(H)
        inv = m.inverse(H)
        p = ((h+eps)*x[0], h*x[1])
        X = tuple(-sum(inv[i][j]*p[j] for j in range(2)) for i in range(2))
        assert tuple(sum(H[i][j]*X[j] for j in range(2)) for i in range(2)) == tuple(-v for v in p)
        det = m.det(H)
        rho = sum(x[i]*X[i] for i in range(2))/r**2
        omega = (x[0]*X[1]-x[1]*X[0])/r**2
        assert rho == -h*(h+eps)/det < 0
        assert omega == h*eps*delta*e[0]*e[1]/det
        energy_rate = sum(p[i]*X[i] for i in range(2))
        assert energy_rate == -sum(X[i]*H[i][j]*X[j] for i,j in product(range(2), repeat=2)) < 0
        # The local variational characterization, independent rational perturbation.
        eta = (Q(2, 3), Q(-5, 7))
        J = lambda u: sum(p[i]*u[i] for i in range(2))+sum(u[i]*H[i][j]*u[j] for i,j in product(range(2), repeat=2))/2
        assert J(tuple(X[i]+eta[i] for i in range(2)))-J(X) == sum(eta[i]*H[i][j]*eta[j] for i,j in product(range(2),repeat=2))/2 > 0
        # Differentiate the exact Hessian along X, not along a fixed radial ray.
        def third(i,j,k):
            return h*delta/r*(kr(i,j)*e[k]+kr(i,k)*e[j]+kr(j,k)*e[i]+(delta-2)*e[i]*e[j]*e[k])
        dH = [[sum(third(k,i,j)*X[k] for k in range(2)) for j in range(2)] for i in range(2)]
        b = (H[0][0]-H[1][1], 2*H[0][1])
        db = (dH[0][0]-dH[1][1], 2*dH[0][1])
        c,d = e
        direct_db = (delta*h*(delta*rho*(c*c-d*d)-4*omega*c*d),
                     2*delta*h*(delta*rho*c*d+omega*(c*c-d*d)))
        assert db == direct_db
        if b[0]**2+b[1]**2:
            seam_rate = (b[0]*db[1]-b[1]*db[0])/(2*(b[0]**2+b[1]**2))
            assert seam_rate == (b[0]*direct_db[1]-b[1]*direct_db[0])/(2*(b[0]**2+b[1]**2))
            counts['moving_principal_seams'] += 1
        rows.append(dict(s=str(x[0]),v=str(x[1]),ds=str(X[0]),dv=str(X[1]),rho=str(rho),omega=str(omega)))
        counts['trajectory_fixtures'] += 1
    # nu=3/2, kappa=2/3: r=t², p=t e. Two attenuation steps compose exactly.
    for t,a,b in product([Q(1,2),Q(1),Q(2)], [Q(1,2),Q(2,3)], [Q(1,3),Q(3,4)]):
        energy = lambda z: Q(2,3)*z**3
        radius = lambda z: z*z
        determinant = lambda z: Q(1,2)/z**2
        assert radius(a*t)/radius(t) == a*a
        assert energy(a*t)/energy(t) == a**3
        assert determinant(a*t)/determinant(t) == a**-2
        assert (b*(a*t))**2 == ((a*b)*t)**2
        released = (energy(t)-energy(a*t))+(energy(a*t)-energy(a*b*t))
        assert released == energy(t)-energy(a*b*t) > 0
        counts['exact_isotropic_steps'] += 1
    # A genuinely rotating, anisotropic finite arrow; no numerical inverse fit.
    eps = Q(7,5)
    x0,x1 = (Q(3,5),Q(4,5)),(Q(1,20),Q(3,80))
    r0,r1,h0,h1,a = Q(1),Q(1,16),Q(1),Q(4),Q(3,16)
    assert sum(z*z for z in x0)==r0*r0 and sum(z*z for z in x1)==r1*r1
    p0,p1 = ((h0+eps)*x0[0],h0*x0[1]),((h1+eps)*x1[0],h1*x1[1])
    assert p1 == tuple(a*z for z in p0)
    assert x0[0]*x1[1]-x0[1]*x1[0] < 0
    W0,W1 = Q(2,3)+eps*x0[0]**2/2,Q(2,3)*Q(1,4)**3+eps*x1[0]**2/2
    assert W0==Q(689,750) and W1==Q(73,6000) and W1<W0
    counts['exact_anisotropic_steps']=1
    # Native positive log's rational Cayley body and certified positive tail.
    def body(N): return 2*sum((Q(1,3)**(2*j+1)/Q(2*j+1) for j in range(N+1)),Q(0))
    def bound(N): return 2*Q(1,3)**(2*N+3)/(Q(2*N+3)*(1-Q(1,9)))
    for n in range(9):
        correction = body(n+1)-body(n)
        assert correction > 0
        assert body(n+12)-body(n) < bound(n)
        assert body(n+1)+bound(n+1) < body(n)+bound(n)
        assert body(n)+correction == body(n+1)
        # A finite future body provides an exact retained-tail ledger witness.
        tail = body(n+12)-body(n)
        assert body(n)+tail == body(n+1)+(tail-correction)
        assert body(n)+tail != body(n+1)+(tail+correction)
        assert body(n) != body(n+1)  # a zero-tail mutation is refuted
        counts['clock_tail_ledgers'] += 1
    # Clock rescaling and omitted clock-curvature term, exact isotropic fixture.
    # x'= -2x, x''=4x in sigma, at sigma'=2 and sigma''=3.
    x=Q(3,5)
    assert Q(2)**2*(4*x)+Q(3)*(-2*x) == 10*x
    assert Q(2)**2*(4*x) != 10*x
    counts['rejected_mutations'] = 3  # tail sign, zero tail, omitted clock term
    return dict(counts=counts,native_generator_rows=rows,
                benchmark=dict(a='1/2',r_before='1',r_after=str(Q(1,2)**2),
                               W_before='2/3',W_after=str(Q(2,3)*Q(1,2)**3),
                               det_before='1/2',det_after=str(Q(1,2)/Q(1,2)**2)),
                anisotropic_benchmark=dict(epsilon=str(eps),a=str(a),state_before=[str(z) for z in x0],state_after=[str(z) for z in x1],W_before=str(W0),W_after=str(W1)),
                log2_enclosure=dict(lower=str(body(8)),upper=str(body(8)+bound(8)),tail_bound=str(bound(8))))

if __name__ == '__main__':
    print(json.dumps(run(),sort_keys=True))
