"""Finite exact shell-capacity evidence, not a physical area measurement."""
from collections import Counter
from fractions import Fraction as Q

import shell_growth_model as sh
import source_response_model as sr


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def build_checks():
    counts = {'same_size_valence_weight_countermodels': 0,
              'exact_lattice_profiles': 0,
              'lattice_link_counts': 0,
              'tree_profiles': 0,
              'tree_link_counts': 0,
              'stationary_source_solutions': 0,
              'independent_ball_flux_checks': 0,
              'radial_tree_checks': 0,
              'nonradial_checks': 0}
    torus, cycle = sh.torus2(7), sh.circulant(49, (1, 2))
    profiles = []
    for n, links in (torus, cycle):
        degrees = Counter(v for u, w, _ in links for v in (u, w))
        require(n == 49 and len(links) == 98 and set(degrees.values()) == {4}
                and all(c == 1 for _, _, c in links), 'same-carrier graph control')
        profiles.append(sh.shell_profile(n, links, radius=3))
        counts['same_size_valence_weight_countermodels'] += 1
    require(profiles[0]['capacities'] == (4, 12, 20), 'torus capacity')
    require(profiles[1]['capacities'] == (4, 6, 6), 'circulant capacity')
    for dimension in range(1, 5):
        for radius in range(1, 5):
            n, links = sh.lattice_ball(dimension, radius)
            profile = sh.shell_profile(n, links, radius=radius)
            expected = tuple(sh.lattice_shell_formula(dimension, j) for j in range(radius))
            require(profile['capacities'] == expected, 'lattice edge count vs binomial identity')
            require(sum(profile['shell_sizes']) == n, 'lattice shell partition')
            counts['exact_lattice_profiles'] += 1
            counts['lattice_link_counts'] += radius
    for degree in range(3, 6):
        for radius in range(1, 5):
            n, links = sh.tree_ball(degree, radius)
            profile = sh.shell_profile(n, links, radius=radius)
            require(profile['capacities'] == tuple(degree*(degree-1)**j
                                                   for j in range(radius)), 'tree capacity')
            require(sum(profile['shell_sizes']) == n, 'tree shell partition')
            counts['tree_profiles'] += 1
            counts['tree_link_counts'] += radius
    named = [('torus49', *torus), ('cycle49', *cycle),
             ('lattice3', *sh.lattice_ball(3, 3)),
             ('tree4', *sh.tree_ball(4, 3))]
    responses = {}
    for name, n, links in named:
        profile = sh.shell_profile(n, links, radius=3)
        net = sh.grounded_ball(n, links, radius=3)
        source = [Q(0)]*len(net['interior'])
        source[net['interior'].index(0)] = Q(1)
        data = sh.stationary_shell_means(n, links, radius=3)
        distances = profile['distances']
        for j in range(3):
            ball = (i for i, k in enumerate(distances) if k <= j)
            require(sr.cut_flux(net, data['field'], ball) == 1, 'independent ball current')
            counts['independent_ball_flux_checks'] += 1
        require(sr.energy(net, data['field'], source) == sr.stationary_cost(net['h'], source),
                'stationary cost')
        counts['stationary_source_solutions'] += 1
        responses[name] = data
        if name == 'tree4':
            require(sh.radial_equitable(n, links, radius=3), 'tree equitability')
            require(data['mean_drops'] == tuple(1/c for c in profile['capacities']), 'radial drops')
            counts['radial_tree_checks'] += 1
        else:
            require(not sh.radial_equitable(n, links, radius=3), 'false radial equitability')
            require(data['mean_drops'][1] != 1/profile['capacities'][1], 'false local inverse')
            counts['nonradial_checks'] += 1
    cubic = responses['lattice3']
    require(cubic['mean_drops'] == (Q(1, 6), Q(13, 378), Q(5, 378)), 'cubic mean')
    lines, links = sh.lattice_ball(1, 4)
    require(sh.shell_profile(lines, links, radius=4)['capacities'] == (2, 2, 2, 2), 'upper-bound control')
    return {
        'counts': counts,
        'U39_equal_local_contract_different_capacity': {
            'carrier_count_each': 49, 'degree_each': 4,
            'unit_links_each': 98,
            'periodic_square_shell_capacity': list(map(int, profiles[0]['capacities'])),
            'cycle_two_jumps_shell_capacity': list(map(int, profiles[1]['capacities']))},
        'U40_declared_3_direction_count': {
            'first_shell_capacities': [6, 30, 78],
            'identity': 'C_3(n)=12*n^2+12*n+6 for integer n>=0',
            'physical_dimension_selected': False,
            'exact_pure_n_squared_at_finite_radii': False},
        'U41_cubic_local_response_obstruction': {
            'mean_drops': list(map(str, cubic['mean_drops'])),
            'inverse_aggregate_capacities': ['1/6', '1/30', '1/78'],
            'shell_field_constant': False},
        'upper_bound_saturation_control': '1D constant capacity 2 obeys a quadratic upper bound',
        'source_type': 'Exact rational finite graph constructions; no physical distance or empirical data',
        'native_axioms_force_quadratic_capacity': False,
        'physical_inverse_square_derived': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(build_checks(), indent=2))
