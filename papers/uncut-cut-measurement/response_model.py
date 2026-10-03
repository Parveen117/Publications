"""Native represented response tensors and two explicitly different cut ledgers.

No probability/readout law, physical curvature, or spacetime is postulated.
Complex tensors are stored as exact (real, imaginary) rational matrices.
"""
from fractions import Fraction as Q

import gravity_model as gm
import grading_model as gr
import interaction_model as im


def zero(rows, columns=None):
    return gm.matrix([[0] * (rows if columns is None else columns) for _ in range(rows)])


def quarter_turn(m):
    im.links_for(m, ())
    return im.kron(gm.identity(m), ((0, -1), (1, 0)))


def tensor(response, z, metric=None):
    """Y^T H Y - i Y^T H Z Y; input columns are calibrated responses."""
    response, z = gm.matrix(response), gm.matrix(z)
    n = len(response)
    h = gm.identity(n) if metric is None else gm.matrix(metric)
    if (len(z), len(z[0]), len(h), len(h[0])) != (n, n, n, n):
        raise ValueError('Response, quarter-turn and metric must share a carrier')
    if gm.mul(z, z) != gm.scale(gm.identity(n), -1):
        raise ValueError('The represented quarter-turn must square to minus identity')
    if h != gm.transpose(h) or gm.mul(gm.transpose(z), h) != gm.scale(gm.mul(h, z), -1):
        raise ValueError('Require a symmetric metric and a skew-adjoint quarter-turn')
    # Positivity of H is a hypothesis; default H=I and transported H are positive.
    return (gm.congruence(response, h),
            gm.scale(gm.congruence(response, gm.mul(h, z)), -1))


def tensor_add(a, b):
    return tuple(gm.add(x, y) for x, y in zip(a, b))


def tensor_sub(a, b):
    return tuple(gm.sub(x, y) for x, y in zip(a, b))


def projector(cut):
    """Euclidean row-space projector, allowing redundant or zero readout rows."""
    cut = gm.matrix(cut)
    reduced, pivots = gm.rref(cut)
    if not pivots:
        return zero(len(cut[0]))
    basis = reduced[:len(pivots)]
    return gm.mul(gm.mul(gm.transpose(basis),
                         gm.inverse(gm.mul(basis, gm.transpose(basis)))), basis)


def split_tensor(response, cut, z):
    """Real-coordinate cut; keeps the imaginary cross-seam term explicit."""
    response, cut, z = gm.matrix(response), gm.matrix(cut), gm.matrix(z)
    if len(cut[0]) != len(response):
        raise ValueError('Cut acts on response outputs, not on parameter labels')
    p = projector(cut)
    q = gm.sub(gm.identity(len(p)), p)
    full = tensor(response, z)
    visible = tensor(gm.mul(p, response), z)
    hidden = tensor(gm.mul(q, response), z)
    cross_operator = gm.add(gm.mul(gm.mul(p, z), q), gm.mul(gm.mul(q, z), p))
    cross = (zero(len(response[0])), gm.scale(gm.congruence(response, cross_operator), -1))
    return {'full': full, 'visible': visible, 'hidden': hidden, 'cross': cross,
            'projector': p, 'global_ledger_closes': gm.commutator(p, z) == zero(len(p))}


def paired_repair(cut, z):
    """Smallest row-space extension invariant under Z, using Z^2=-I."""
    cut, z = gm.matrix(cut), gm.matrix(z)
    n = len(cut[0])
    if (len(z), len(z[0])) != (n, n) or gm.mul(z, z) != gm.scale(gm.identity(n), -1):
        raise ValueError('Quarter-turn must act on the readout carrier')
    reduced, pivots = gm.rref(cut + gm.mul(cut, z))
    observer = reduced[:len(pivots)] if pivots else (tuple(Q(0) for _ in range(n)),)
    return {'observer': observer, 'rank': len(pivots),
            'extra_scalar_channels': len(pivots) - gr.rank(cut)}


def order_response(u, v):
    """Signed difference UV - VU: chronological V,U minus U,V."""
    return gm.commutator(u, v)


def mean_covariance(points, weights):
    points = gm.matrix(points)  # one vector per ensemble entry
    weights = tuple(Q(p) for p in weights)
    if len(points) != len(weights) or any(p <= 0 for p in weights) or sum(weights) != 1:
        raise ValueError('A finite ensemble needs positive weights summing to one')
    d = len(points[0])
    mean = tuple(sum(p*x[j] for p, x in zip(weights, points)) for j in range(d))
    cov = gm.matrix([[sum(p*(x[i]-mean[i])*(x[j]-mean[j])
                          for p, x in zip(weights, points))
                      for j in range(d)] for i in range(d)])
    return mean, cov


def conditional_ledger(points, weights, labels):
    """CID-1 total covariance applied to signed native response vectors.

    This is conditioning on ensemble labels, not a real-coordinate projection.
    """
    points = gm.matrix(points)
    weights, labels = tuple(Q(p) for p in weights), tuple(labels)
    mean, total = mean_covariance(points, weights)
    if len(labels) != len(points):
        raise ValueError('One cut label is required for each ensemble entry')
    keys = tuple(dict.fromkeys(labels))
    masses, means, within = [], [], zero(len(mean))
    for key in keys:
        indices = [i for i, a in enumerate(labels) if a == key]
        mass = sum(weights[i] for i in indices)
        mu, cov = mean_covariance([points[i] for i in indices],
                                  [weights[i]/mass for i in indices])
        masses.append(mass)
        means.append(mu)
        within = gm.add(within, gm.scale(cov, mass))
    _, between = mean_covariance(means, masses)
    return {'mean': mean, 'total': total, 'recognized': between, 'discarded': within,
            'block_labels': keys, 'block_weights': tuple(masses), 'block_means': tuple(means)}


def target_repair(cut, target):
    """Inherited linear target-faithfulness rank formula (Spectral II/T32)."""
    cut, target = gm.matrix(cut), gm.matrix(target)
    if len(cut[0]) != len(target[0]):
        raise ValueError('Input cut and response target must share a source')
    reduced, pivots = gm.rref(cut + target)
    return {'rank': len(pivots), 'extra_scalar_channels': len(pivots)-gr.rank(cut),
            'observer': reduced[:len(pivots)] if pivots else (reduced[0],)}
