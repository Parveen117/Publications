"""Exact controls consuming original thermo and CID modules without edits."""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
from math import factorial
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[2]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def run():
    if sys.version_info[:2] not in ((3, 11), (3, 12)) or not __debug__:
        raise RuntimeError('Python 3.11/3.12 without -O required')
    thermo = load('cf_thermo', 'papers/thermo-compass-foundations/model.py')
    cid = load('cf_cid', 'papers/curvature-information-duality/certificates/cid1_curvature_information_duality.py')
    directions = [(Q(1), Q(0)), (Q(0), Q(1)), (Q(3, 5), Q(4, 5)),
                  (Q(4, 5), Q(-3, 5)), (Q(-1), Q(0)), (Q(-3, 5), Q(-4, 5))]
    counts = dict(energy_jets=0, full_christoffel_comparisons=0,
                  thermo_closures=0, radial_rotation=0, cid_ledgers=0)
    signs = set()
    native_seams = []
    kr = lambda i, j: Q(i == j)
    for nu, root, epsilon, e in product(
            [Q(4, 3), Q(3, 2), Q(5, 3)],
            [Q(1, 4), Q(1, 2), Q(1), Q(2)],
            [Q(0), Q(1, 4), Q(1)], directions):
        # r=root^denominator makes every power in the rational benchmark exact.
        r = root ** nu.denominator
        h = root ** int((nu-2)*nu.denominator)  # kappa=1/nu
        delta = nu-2
        q = tuple(r*x for x in e)
        hessian = thermo.matrix([[h*(kr(i, j)+delta*e[i]*e[j])
                                  + epsilon*kr(i, 0)*kr(j, 0)
                                  for j in range(2)] for i in range(2)])

        def third(i, j, k):
            return h*delta/r*(kr(i, j)*e[k]+kr(i, k)*e[j]+kr(j, k)*e[i]
                               +(delta-2)*e[i]*e[j]*e[k])

        def fourth(i, j, k, l):
            return h*delta/r**2*(
                kr(i, j)*kr(k, l)+kr(i, k)*kr(j, l)+kr(i, l)*kr(j, k)
                +(delta-2)*(kr(i, j)*e[k]*e[l]+kr(i, k)*e[j]*e[l]
                            +kr(i, l)*e[j]*e[k]+kr(j, k)*e[i]*e[l]
                            +kr(j, l)*e[i]*e[k]+kr(k, l)*e[i]*e[j])
                +(delta-2)*(delta-4)*e[i]*e[j]*e[k]*e[l])

        dh = tuple(thermo.matrix([[third(i, j, k) for k in range(2)]
                                 for j in range(2)]) for i in range(2))
        ddh = tuple(tuple(thermo.matrix([[fourth(i, j, k, l) for l in range(2)]
                                        for k in range(2)]) for j in range(2))
                    for i in range(2))
        # Independent Taylor expansion of the energy validates every derivative jet.
        z = {(1, 0): 2*q[0]/r**2, (0, 1): 2*q[1]/r**2,
             (2, 0): 1/r**2, (0, 2): 1/r**2}
        power, jet, choose = {(0, 0): Q(1)}, {}, Q(1)
        radial_energy = root**int(nu*nu.denominator)/nu
        for degree in range(5):
            if degree:
                nxt = {}
                for (i, j), a in power.items():
                    for (k, l), b in z.items():
                        if i+j+k+l <= 4:
                            key = (i+k, j+l)
                            nxt[key] = nxt.get(key, Q(0))+a*b
                power = nxt
                choose *= (nu/2-degree+1)/degree
            for key, value in power.items():
                jet[key] = jet.get(key, Q(0))+radial_energy*choose*value
        jet[(2, 0)] = jet.get((2, 0), Q(0))+epsilon/2
        for order in (2, 3, 4):
            for indices in product(range(2), repeat=order):
                i, j = indices.count(0), indices.count(1)
                actual = jet.get((i, j), Q(0))*factorial(i)*factorial(j)
                expected = (hessian[indices[0]][indices[1]] if order == 2
                            else third(*indices) if order == 3 else fourth(*indices))
                assert actual == expected
        thermo.require_stable(hessian)
        determinant = thermo.det(hessian)
        expected_det = h*h*(nu-1)+epsilon*h*(1+delta*e[1]**2)
        assert determinant == expected_det > 0
        # Independent directional identities from homogeneity and radial differentiation.
        for j, k in product(range(2), repeat=2):
            lhs = sum(e[i]*third(i, j, k) for i in range(2))
            assert lhs == delta/r*(hessian[j][k]-epsilon*kr(j, 0)*kr(k, 0))
            for l in range(2):
                assert sum(e[i]*fourth(i, j, k, l) for i in range(2)) == (delta-1)/r*third(j, k, l)
        f = thermo.curvature(hessian, *dh)
        direct = thermo.curvature_direct(hessian, dh, ddh)
        assert f == direct
        curvature = thermo.gaussian_curvature(hessian, *dh)
        closed = epsilon*h*h*delta*delta*(e[0]**2-(nu-1)*e[1]**2)/(4*r*r*determinant**2)
        assert curvature == closed
        signs.add((curvature > 0)-(curvature < 0))
        if epsilon == 0:
            assert f == thermo.ZERO
        counts['full_christoffel_comparisons'] += 1

        temperature = 100+h*q[0]+epsilon*q[0]
        pressure = 100-h*q[1]
        volume = 100+q[1]
        assert temperature > 0 and pressure > 0 and volume > 0 and 100+q[0] > 0
        responses = thermo.responses(hessian, temperature, volume)
        assert responses['C_P']/responses['C_V'] == responses['K_S']/responses['K_T']
        assert responses['C_P'] >= responses['C_V']
        counts['thermo_closures'] += 1
        formation = root**int(nu*nu.denominator)/nu+epsilon*q[0]**2/2
        assert formation > 0
        gradient = (h*q[0]+epsilon*q[0], h*q[1])
        radial_work = sum(gradient[i]*e[i] for i in range(2))
        assert radial_work == h*r+epsilon*r*e[0]**2 > 0
        assert sum(hessian[i][j]*e[i]*e[j] for i, j in product(range(2), repeat=2)) == h*(nu-1)+epsilon*e[0]**2
        counts['energy_jets'] += 1

        a, b = hessian[0][0]-hessian[1][1], 2*hessian[0][1]
        norm2 = a*a+b*b
        if norm2:
            # Angle derivative independently evaluated from Hessian third jets.
            ar = sum(e[k]*(dh[k][0][0]-dh[k][1][1]) for k in range(2))
            br = 2*sum(e[k]*dh[k][0][1] for k in range(2))
            angle_derivative = (a*br-b*ar)/(2*norm2)
            closed_angle = epsilon*delta*delta*h*(2*e[0]*e[1])/(2*r*norm2)
            assert angle_derivative == closed_angle
            if epsilon and e[0]*e[1]:
                assert angle_derivative != 0  # dropped radial rotation must fail
            counts['radial_rotation'] += 1
        if len(native_seams) < 12 or (nu == Q(3, 2) and root == 1 and epsilon == 1):
            native_seams.append(dict(a=str(a), b=str(b), squared_norm=str(norm2)))

    assert signs == {-1, 0, 1}
    ledger_rows = []
    probabilities = (Q(1, 4), Q(1, 4), Q(1, 8), Q(1, 8), Q(1, 4))
    for t, e in product([Q(1, 2), Q(1), Q(2)], directions):
        f = (-e[1], e[0])
        vectors = [tuple(x/t for x in e), tuple(-x/t for x in e),
                   tuple(2*x/t for x in f), tuple(-2*x/t for x in f), (Q(0), Q(0))]
        stats = tuple(tuple(v[i] for v in vectors) for i in range(2))
        covariance = cid.covariance(probabilities, stats)
        expected = [[(e[i]*e[j]/2+f[i]*f[j])/t**2 for j in range(2)] for i in range(2)]
        assert covariance == expected and cid.is_pd(covariance)
        for blocks in [[(0, 1), (2, 3), (4,)], [(0,), (1,), (2,), (3,), (4,)],
                       [(0,), (1, 2, 3, 4)]]:
            _, _, recognized, discarded = cid.coarse_grain(probabilities, stats, blocks)
            assert cid.madd(recognized, discarded) == covariance
            assert cid.is_psd(recognized) and cid.is_psd(discarded)
            if len(blocks) == 3:
                assert recognized == [[0, 0], [0, 0]] and discarded == covariance
                assert recognized != covariance  # erased-memory control
            if len(blocks) == 5:
                assert discarded == [[0, 0], [0, 0]]
            counts['cid_ledgers'] += 1
        ledger_rows.append(dict(t=str(t), determinant=str(cid.det(covariance))))

    # Full formation energy decreases inward while determinant increases.
    benchmark = [dict(r=str(t**4), formation=str(Q(2, 3)*t**6),
                      determinant=str(Q(1, 2)/t**4), radial_response=str(t**2))
                 for t in (Q(1), Q(1, 2))]
    assert Q(benchmark[1]['formation']) < Q(benchmark[0]['formation'])
    assert Q(benchmark[1]['determinant']) > Q(benchmark[0]['determinant'])
    assert Q(benchmark[1]['radial_response']) < Q(benchmark[0]['radial_response'])
    return dict(counts=counts, curvature_signs=sorted(signs), native_seams=native_seams,
                cid_examples=ledger_rows, fixed_outcome_probabilities=[str(p) for p in probabilities],
                benchmark=benchmark, scope=dict(physical_gravity=False, SI_c=False,
                global_fixed_statistic_Fisher_model=False, infinite_total_information=False))


if __name__ == '__main__':
    print(json.dumps(run(), sort_keys=True))
