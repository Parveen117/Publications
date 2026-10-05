"""Direct transfer of the six thermo lambda/z definitions. Python 3.12.

The primary check differentiates the four potential functions and evaluates
Jacobian ratios. Closed Hessian formulas provide a separate comparison.
Only point values and first derivatives of response jets are consumed.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SPEC = importlib.util.spec_from_file_location('nsr2_previous', HERE/'nsr1_response.py')
previous = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(previous)
Jet = previous.Jet

NAMES = ('lambda_q', 'lambda_y', 'lambda_x', 'z_p', 'z_x', 'lambda_p')


def diff(jet, axis):
    """Differentiate the stored finite Taylor polynomial once."""
    if axis not in (0, 1):
        raise ValueError('Only two independent chart variables')
    out = {}
    for key, value in jet.c.items():
        if key[axis]:
            target = list(key)
            target[axis] -= 1
            out[tuple(target)] = key[axis]*value
    return Jet(out)


def jacobian(f, g):
    return f.derivative(0)*g.derivative(1)-f.derivative(1)*g.derivative(0)


def constrained(f, variable, held):
    denominator = jacobian(variable, held)
    if not denominator:
        raise ValueError('The requested constrained chart is singular')
    return jacobian(f, held)/denominator


def chart(energy, x0, y0):
    x, y = previous.variables(x0, y0)
    p, q = diff(energy, 0), diff(energy, 1)
    return dict(U=energy, x=x, y=y, p=p, q=q,
                Fx=energy-p*x, Fy=energy-q*y, G=energy-p*x-q*y)


def derivative_readouts(energy, x0, y0):
    h, _, _ = previous.metric_jets(energy)
    previous.check_positive(h)
    d = chart(energy, x0, y0)
    triples = (('G', 'x', 'q'), ('Fx', 'x', 'y'), ('Fy', 'y', 'x'),
               ('G', 'y', 'p'), ('U', 'q', 'x'), ('Fx', 'q', 'p'))
    return {name: constrained(d[f], d[v], d[held])
            for name, (f, v, held) in zip(NAMES, triples)}


def formula_jets(energy, x0, y0):
    """Closed-form values and first derivatives from the four-order energy jet."""
    h, _, _ = previous.metric_jets(energy)
    previous.check_positive(h)
    x, y = previous.variables(x0, y0)
    p, q = diff(energy, 0), diff(energy, 1)
    a, b, c = diff(p, 0), diff(p, 1), diff(q, 1)
    delta = a*c-b*b
    return dict(zip(NAMES, (-x*delta/c, -x*a, -y*c, -y*delta/a, q/c, q*a/delta)))


def formula_readouts(energy, x0, y0):
    return {k: v.derivative() for k, v in formula_jets(energy, x0, y0).items()}


def thermal_readouts(energy, entropy, volume):
    d = chart(energy, entropy, volume)
    pressure = -d['q']
    return dict(lambda_P=constrained(d['G'], d['x'], pressure),
                lambda_V=constrained(d['Fx'], d['x'], d['y']),
                lambda_S=constrained(d['Fy'], d['y'], d['x']),
                z_T=constrained(d['G'], d['y'], d['p']),
                z_S=constrained(d['U'], pressure, d['x']),
                lambda_T=constrained(d['Fx'], pressure, d['p']))


def polynomial_fixture(x0, y0, seed):
    x, y = previous.variables(x0, y0)
    dx, dy = x-x0, y-y0
    a, b, c = F(2+seed % 4), F(seed % 3-1), F(3+seed % 5)
    p, q = F(4+seed % 2), F((seed % 3-1)*3)
    third = [F((seed*(k+2)+k) % 5-2) for k in range(4)]
    fourth = [F((seed*(k+3)+k) % 7-3) for k in range(5)]
    u = p*dx+q*dy+(a*dx*dx+2*b*dx*dy+c*dy*dy)/2
    u += sum(third[k]*dx**(3-k)*dy**k for k in range(4))
    u += sum(fourth[k]*dx**(4-k)*dy**k for k in range(5))
    return u


def signed_pair(sign):
    if sign not in (-1, 1):
        raise ValueError('Expected a signed-coupling comparison')
    x, y = previous.variables(1, 1)
    dx, dy = x-1, y-1
    return 2*dx+3*dy+dx*dx+sign*dx*dy+F(5, 2)*dy*dy


def check_identities(energy, x0, y0):
    values = derivative_readouts(energy, x0, y0)
    assert values == formula_readouts(energy, x0, y0)
    h, _, _ = previous.metric_jets(energy)
    a, b, c, delta = h[0][0], h[0][1], h[1][1], previous.nt.det(h)
    q = energy.derivative(1)
    lq, ly, lx, zp, zx, lp = (values[k] for k in NAMES)
    assert lq*c/delta == ly/a == -x0
    assert lx*zx == lp*zp == -y0*q
    assert lq*lx == ly*zp
    assert lq-ly == x0*b*b/c
    assert zp-lx == y0*b*b/a
    assert lp-zx == q*b*b/(c*delta)
    if x0 and y0:
        assert (lq/ly)*(lx/zp) == 1
    return dict(values=values, gamma_1_extension=delta/(a*c),
                gamma_2_extension=a*c/delta,
                literal_lambda_ratios_defined=bool(x0 and y0))


def source_digest():
    path = HERE/'NSR2_SOURCE_PINS.json'
    raw = path.read_bytes()
    for relative, expected in json.loads(raw)['sha256'].items():
        actual = hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError('Source mismatch: '+relative)
    return hashlib.sha256(raw).hexdigest()


def evidence():
    previous.require_runtime()
    digest = source_digest()
    generic = []
    for seed in range(36):
        x, y = (F(0), F(1), F(3, 2))[seed % 3], (F(0), F(2))[seed % 2]
        data = check_identities(polynomial_fixture(x, y, seed), x, y)
        generic.append(dict(seed=seed, x=x, y=y, **data))

    # Reproduce the original definitions with P=-U_V, not q=U_V.
    s, v = previous.variables(1, 1)
    ds, dv = s-1, v-1
    thermo_u = 10*ds-10*dv+ds*ds+ds*dv+F(5, 2)*dv*dv+ds*ds*dv/2
    thermo = thermal_readouts(thermo_u, 1, 1)
    expected_thermo = dict(zip(('lambda_P','lambda_V','lambda_S','z_T','z_S','lambda_T'),
                              (F(-9,5),F(-2),F(-5),F(-9,2),F(2),F(20,9))))
    assert thermo == expected_thermo
    generic_thermo = derivative_readouts(thermo_u, 1, 1)
    assert generic_thermo['z_x'] == -thermo['z_S']
    assert generic_thermo['lambda_p'] == -thermo['lambda_T']
    assert thermo['lambda_S']*thermo['z_S'] == thermo['lambda_T']*thermo['z_T'] == -10

    bi = check_identities(previous.born_infeld(F(1,2),F(1,2)), F(1,2), F(1,2))
    expected_bi = dict(zip(NAMES, (F(-3,8),F(-2,5),F(-2,5),F(-3,8),F(5,8),F(2,3))))
    assert bi['values'] == expected_bi
    assert (bi['gamma_1_extension'],bi['gamma_2_extension']) == (F(15,16),F(16,15))

    # Fully explicit conservative circuit model, in declared normalized units.
    q0, phi0, capacity, inductance = F(2), F(3), F(5), F(7)
    charge, flux = previous.variables(q0, phi0)
    circuit_u = charge*charge/(2*capacity)+flux*flux/(2*inductance)
    circuit = check_identities(circuit_u, q0, phi0)
    assert tuple(circuit['values'][k] for k in NAMES) == (
        -q0/capacity,-q0/capacity,-phi0/inductance,-phi0/inductance,phi0,phi0)

    # Same six point properties, opposite mixed response, distinct first tower.
    pair = []
    for sign in (-1, 1):
        energy = signed_pair(sign)
        properties = check_identities(energy, F(1), F(1))
        slope = formula_jets(energy, 1, 1)['z_x'].derivative(0)
        assert slope == F(sign, 5)
        pair.append(dict(sign=sign, mixed=energy.derivative(0,1),
                         properties=properties, next_z_x_derivative=slope))
    assert pair[0]['properties']['values'] == pair[1]['properties']['values']

    # A zero local coordinate is not a license to delete the physical origin.
    origin_energy = signed_pair(1)
    physical = derivative_readouts(origin_energy, 1, 1)
    changed_origin = derivative_readouts(origin_energy, 0, 0)
    assert all(changed_origin[k] == 0 for k in NAMES[:4])
    assert all(physical[k] != 0 for k in NAMES[:4])

    # Re-run the exact predecessor evidence without writing any old file.
    old = previous.evidence()
    saved_old = json.loads((HERE/'NSR1_RESULT.json').read_text())
    assert old == saved_old

    return previous.encode(dict(
        status='PASS_NSR2_DIRECT_CONJUGATE_PROPERTY_TRANSFER', runtime='Python 3.12',
        source_manifest_sha256=digest, generic_cases=generic,
        original_thermodynamic_properties=thermo,
        electromagnetic_Born_Infeld_properties=bi, linear_circuit_properties=circuit,
        signed_coupling_pair=pair,
        origin_control=dict(physical=physical,redefined_local=changed_origin),
        previous_NSR1_evidence_unchanged=True,
        scope='Exact four-potential derivative transfer; empirical validation and priority are separate'
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = evidence()
    raw = (json.dumps(data, indent=2, sort_keys=True)+'\n').encode()
    digest = hashlib.sha256(raw).hexdigest()
    output, expected = HERE/'NSR2_RESULT.json', HERE/'EXPECTED_NSR2.sha256'
    if args.write:
        output.write_bytes(raw)
        expected.write_text(digest+'\n')
    elif output.read_bytes() != raw or expected.read_text().strip() != digest:
        raise RuntimeError('NSR2 evidence differs from the committed result')
    print(json.dumps(dict(status=data['status'],sha256=digest,
                          generic_cases=len(data['generic_cases']),NSR1_preserved=True)))


if __name__ == '__main__':
    main()
