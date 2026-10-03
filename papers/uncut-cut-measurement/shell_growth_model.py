"""Exact graph-distance shell capacity, with no physical distance identification."""
from fractions import Fraction as Q
from itertools import product

import source_response_model as sr


def _edges(n, pairs):
    return tuple((u, v, Q(1)) for u, v in sorted({tuple(sorted((u, v)))
            for u, v in pairs if u != v}))


def torus2(side):
    """Four-regular finite torus, used as a candidate graph, not spacetime."""
    if type(side) is not int or side < 7:
        raise ValueError('Side must be at least seven')
    return side*side, _edges(side*side, ((side*x+y, side*((x+1)%side)+y)
                           for x in range(side) for y in range(side))) + _edges(
        side*side, ((side*x+y, side*x+(y+1)%side)
                    for x in range(side) for y in range(side)))


def circulant(size, jumps):
    if (type(size) is not int or size < 7 or not jumps or
            any(type(j) is not int or j <= 0 or j*2 >= size for j in jumps)
            or len(set(jumps)) != len(jumps)):
        raise ValueError('Positive distinct short jumps on a finite cycle required')
    return size, _edges(size, ((u, (u+j)%size)
                               for u in range(size) for j in jumps))


def lattice_ball(dimension, radius):
    if (type(dimension) is not int or not 1 <= dimension <= 4
            or type(radius) is not int or radius < 1):
        raise ValueError('Supported lattice dimensions 1..4 and positive radius')
    vertices = [x for x in product(range(-radius, radius+1), repeat=dimension)
                if sum(abs(t) for t in x) <= radius]
    vertices.sort(key=lambda x: (sum(abs(t) for t in x), x))
    labels = {x: i for i, x in enumerate(vertices)}
    pairs = ((labels[x], labels[x[:k]+(x[k]+1,)+x[k+1:]])
             for x in vertices for k in range(dimension)
             if x[:k]+(x[k]+1,)+x[k+1:] in labels)
    return len(vertices), _edges(len(vertices), pairs)


def tree_ball(degree, radius):
    if (type(degree) is not int or degree < 3
            or type(radius) is not int or radius < 1):
        raise ValueError('Degree >=3 and positive radius required')
    depths, edges = [0], []
    node = 0
    while node < len(depths):
        if depths[node] >= radius:
            node += 1
            continue
        for _ in range(degree if node == 0 else degree-1):
            child = len(depths)
            depths.append(depths[node]+1)
            edges.append((node, child, Q(1)))
        node += 1
    return len(depths), tuple(edges)


def distances(n, edges, root=0):
    if type(n) is not int or n < 1 or not 0 <= root < n:
        raise ValueError('Finite rooted graph required')
    neighbors = [[] for _ in range(n)]
    for u, v, c in edges:
        if not 0 <= u < n or not 0 <= v < n or u == v or Q(c) <= 0:
            raise ValueError('Simple positive graph required')
        neighbors[u].append(v)
        neighbors[v].append(u)
    d, queue = [None]*n, [root]
    d[root] = 0
    for u in queue:
        for v in neighbors[u]:
            if d[v] is None:
                d[v] = d[u]+1
                queue.append(v)
    if any(k is None for k in d):
        raise ValueError('Connected graph required')
    return tuple(d)


def shell_profile(n, edges, root=0, radius=3):
    d = distances(n, edges, root)
    if type(radius) is not int or radius < 1 or max(d) < radius:
        raise ValueError('Declared radius must fit the graph')
    counts = tuple(sum(int(k == j) for k in d) for j in range(radius+1))
    capacities = tuple(sum(Q(c) for u, v, c in edges
                           if min(d[u], d[v]) <= j < max(d[u], d[v]))
                       for j in range(radius))
    return {'distances': d, 'shell_sizes': counts, 'capacities': capacities,
            'ball_sizes': tuple(sum(counts[:j+1]) for j in range(radius+1))}


def grounded_ball(n, edges, root=0, radius=3):
    d = distances(n, edges, root)
    if radius <= 0 or max(d) < radius:
        raise ValueError('A reachable grounded shell is required')
    boundary = [i for i, k in enumerate(d) if k >= radius]
    return sr.network(n, edges, boundary)


def radial_equitable(n, edges, root=0, radius=3):
    """Check whether every vertex at one shell sees equal inward/outward weights."""
    d = distances(n, edges, root)
    if radius <= 0 or max(d) < radius:
        raise ValueError('A reachable grounded shell is required')
    weights = [[Q(0), Q(0)] for _ in range(n)]
    for u, v, c in edges:
        c = Q(c)
        if d[u] < d[v]:
            weights[u][1] += c
            weights[v][0] += c
        elif d[v] < d[u]:
            weights[v][1] += c
            weights[u][0] += c
    return all(len({tuple(weights[i]) for i, k in enumerate(d) if k == j}) == 1
               for j in range(1, radius))


def stationary_shell_means(n, edges, root=0, radius=3):
    net = grounded_ball(n, edges, root, radius)
    source = [Q(0)]*len(net['interior'])
    source[net['interior'].index(root)] = Q(1)
    field = sr.solve(net, source)
    d = distances(n, edges, root)
    means = tuple(sum((field[i] for i, k in enumerate(d) if k == j), Q(0))
                  /sum(int(k == j) for k in d) for j in range(radius+1))
    return {'field': field, 'means': means,
            'mean_drops': tuple(means[j]-means[j+1] for j in range(radius))}


def lattice_shell_formula(dimension, shell):
    """Exact outward channel count of the unweighted L1 lattice ball."""
    from math import comb
    if dimension not in (1, 2, 3, 4) or type(shell) is not int or shell < 0:
        raise ValueError('Dimension 1..4 and nonnegative shell required')
    if shell == 0:
        return 2*dimension
    def sphere(dim):
        return sum(2**k*comb(dim, k)*comb(shell-1, k-1)
                   for k in range(1, min(dim, shell)+1))
    return dimension*(sphere(dimension)+sphere(dimension-1))
