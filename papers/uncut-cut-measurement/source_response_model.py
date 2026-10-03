"""Exact static source-cost candidate; physical interpretation is conditional.

No distance, mass, clock, spatial dimension or gravity datum is inferred here.
Weights, graph, quadratic mismatch cost and source coupling are supplied.
"""
from fractions import Fraction as Q

import gravity_model as gm
import interaction_model as im


def column(values):
    return gm.matrix([[Q(x)] for x in values])


def vector(a):
    a = gm.matrix(a)
    if len(a[0]) != 1:
        raise ValueError('Expected a column')
    return tuple(row[0] for row in a)


def network(n, edges, boundary):
    if type(n) is not int or n < 2:
        raise ValueError('At least two vertices required')
    boundary = tuple(boundary)
    if (not boundary or len(set(boundary)) != len(boundary)
            or any(type(v) is not int or not 0 <= v < n for v in boundary)):
        raise ValueError('Distinct grounded boundary vertices required')
    edges = tuple((u, v, Q(c)) for u, v, c in edges)
    im.links_for(n, [(u, v) for u, v, _ in edges])
    if any(c <= 0 for _, _, c in edges):
        raise ValueError('All link stiffnesses must be positive')
    if any(not set(comp).intersection(boundary)
           for comp in im.components(n, [(u, v) for u, v, _ in edges])):
        raise ValueError('Every component must meet the grounded boundary')
    interior = tuple(v for v in range(n) if v not in boundary)
    if not interior:
        raise ValueError('At least one free vertex required')
    lap = [[Q(0) for _ in range(n)] for _ in range(n)]
    for u, v, c in edges:
        lap[u][u] += c
        lap[v][v] += c
        lap[u][v] -= c
        lap[v][u] -= c
    h = gm.matrix([[lap[u][v] for v in interior] for u in interior])
    return {'n': n, 'edges': edges, 'boundary': boundary,
            'interior': interior, 'h': h}


def solve(net, sources):
    sources = tuple(Q(x) for x in sources)
    if len(sources) != len(net['interior']):
        raise ValueError('One source value per free vertex required')
    values = vector(gm.mul(gm.inverse(net['h']), column(sources)))
    field = [Q(0)] * net['n']
    for v, x in zip(net['interior'], values):
        field[v] = x
    return tuple(field)


def currents(net, field):
    field = tuple(Q(x) for x in field)
    if len(field) != net['n']:
        raise ValueError('One field value per vertex required')
    return tuple(c*(field[u]-field[v]) for u, v, c in net['edges'])


def cut_flux(net, field, selected):
    selected = set(selected)
    if not selected <= set(net['interior']):
        raise ValueError('Cut region must contain only free vertices')
    return sum(((int(u in selected)-int(v in selected))*j
                for (u, v, _), j in zip(net['edges'], currents(net, field))), Q(0))


def energy(net, field, sources):
    field, sources = tuple(map(Q, field)), tuple(map(Q, sources))
    if len(field) != net['n'] or len(sources) != len(net['interior']):
        raise ValueError('Wrong field/source length')
    if any(field[v] for v in net['boundary']):
        raise ValueError('Grounded field must vanish on boundary')
    return (sum((c*(field[u]-field[v])**2/2 for u, v, c in net['edges']), Q(0))
            - sum((q*field[v] for v, q in zip(net['interior'], sources)), Q(0)))


def stationary_cost(h, b):
    h, b = gm.matrix(h), column(b)
    return -gm.mul(gm.transpose(b), gm.mul(gm.inverse(h), b))[0][0]/2


def eliminate(h, b, retained):
    """Stationary elimination, not an exact quotient of the whole dynamics."""
    h, b = gm.matrix(h), tuple(map(Q, b))
    n, retained = len(h), tuple(retained)
    if (len(h[0]) != n or len(b) != n or not retained
            or len(set(retained)) != len(retained)
            or any(type(v) is not int or not 0 <= v < n for v in retained)):
        raise ValueError('Invalid retained coordinates')
    hidden = tuple(i for i in range(n) if i not in retained)
    if not hidden:
        return {'h': gm.matrix([[h[i][j] for j in retained] for i in retained]),
                'b': tuple(b[i] for i in retained), 'constant': Q(0)}
    def block(rows, cols):
        return gm.matrix([[h[i][j] for j in cols] for i in rows])
    a, c, d = block(retained, retained), block(retained, hidden), block(hidden, hidden)
    inv = gm.inverse(d)
    bh = column([b[i] for i in hidden])
    return {'h': gm.sub(a, gm.mul(gm.mul(c, inv), gm.transpose(c))),
            'b': vector(gm.sub(column([b[i] for i in retained]), gm.mul(gm.mul(c, inv), bh))),
            'constant': -gm.mul(gm.transpose(bh), gm.mul(inv, bh))[0][0]/2}


def native_generator(h):
    h = gm.matrix(h)
    if h != gm.transpose(h):
        raise ValueError('Symmetric source-cost Hessian required')
    return im.kron(h, ((0, -1), (1, 0)))


def affine_step(h, source, state, step=Q(1, 2)):
    """(I-hG)x'=(I+hG)x-2h Zb, with G=H tensor R."""
    h, step = gm.matrix(h), Q(step)
    if not step:
        raise ValueError('Nonzero protocol step required')
    g = native_generator(h)
    z = im.kron(gm.identity(len(h)), ((0, -1), (1, 0)))
    x, b = column(state), column(source)
    if len(x) != len(g) or len(b) != len(g):
        raise ValueError('Two coordinates per free vertex required')
    unit = gm.identity(len(g))
    rhs = gm.sub(gm.mul(gm.add(unit, gm.scale(g, step)), x),
                 gm.scale(gm.mul(z, b), 2*step))
    return vector(gm.mul(gm.inverse(gm.sub(unit, gm.scale(g, step))), rhs))


def paired_cost(h, source, state):
    k = im.kron(h, gm.identity(2))
    x, b = column(state), column(source)
    return (gm.mul(gm.transpose(x), gm.mul(k, x))[0][0]/2
            - gm.mul(gm.transpose(b), x)[0][0])


def radial(capacities, source=1):
    """Exact path/lumped-shell stiffness model; shell symmetry is a hypothesis."""
    capacities = tuple(map(Q, capacities))
    if not capacities or any(c <= 0 for c in capacities):
        raise ValueError('Positive shell stiffnesses required')
    q = Q(source)
    drops = tuple(q/c for c in capacities)
    field = [Q(0)] * (len(drops)+1)
    for i in reversed(range(len(drops))):
        field[i] = field[i+1]+drops[i]
    return {'field': tuple(field), 'drops': drops,
            'flux': tuple(c*d for c, d in zip(capacities, drops))}


def source_probe_force(capacities, source=1, probe=1, spacing=1):
    """Inward signed finite-difference cross-energy response, conditional units."""
    spacing = Q(spacing)
    if spacing <= 0:
        raise ValueError('Positive independent displacement calibration required')
    return tuple(-Q(probe)*d/spacing for d in radial(capacities, source)['drops'])


def source_only_cross_energy(h, source, probe):
    source, probe = tuple(map(Q, source)), tuple(map(Q, probe))
    both = tuple(a+b for a, b in zip(source, probe))
    if len(source) != len(probe):
        raise ValueError('Equal source dimensions required')
    return stationary_cost(h, both)-stationary_cost(h, source)-stationary_cost(h, probe)


def nonnegative_interval(value, error):
    value, error = Q(value), Q(error)
    if value <= 0 or error < 0 or value-error <= 0:
        raise ValueError('Positive magnitude interval separated from zero required')
    return value-error, value+error


def scale_test(first, first_error, second, second_error, predicted_ratio, prediction_error=0):
    """Bounded-error rejection only; overlap does not validate a physical model."""
    lo1, hi1 = nonnegative_interval(first, first_error)
    lo2, hi2 = nonnegative_interval(second, second_error)
    ratio = Q(predicted_ratio)
    prediction_interval = nonnegative_interval(ratio, prediction_error)
    interval = (lo2/hi1, hi2/lo1)
    separated = (interval[1] < prediction_interval[0]
                 or prediction_interval[1] < interval[0])
    return {'ratio_interval': interval, 'prediction': ratio,
            'prediction_interval': prediction_interval,
            'verdict': 'REJECTED' if separated else 'UNRESOLVED'}
