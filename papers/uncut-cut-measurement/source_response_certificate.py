"""Finite exact evidence for a conditional source-response candidate, U36-U38."""
from fractions import Fraction as Q
from itertools import combinations

import gravity_model as gm
import interaction_model as im
import source_response_model as sr


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def subsets(items):
    items = tuple(items)
    return (s for k in range(1, len(items)+1) for s in combinations(items, k))


def build_checks():
    counts = dict(weighted_graphs=0, source_solutions=0, cut_flux_ledgers=0,
                  stationary_eliminations=0, energy_completions=0,
                  native_affine_laws=0, radial_profiles=0,
                  pair_energy_force_checks=0, nonlinear_scaling_controls=0,
                  bounded_error_ratio_checks=0)
    for n in range(2, 5):
        possible = tuple(combinations(range(n), 2))
        for mask in range(1, 1 << len(possible)):
            links = tuple(e for i, e in enumerate(possible) if mask & (1 << i))
            if len(im.components(n, links)) != 1:
                continue
            for weighted in [False, True]:
                edges = [(u, v, Q(i+1, i+2) if weighted else Q(1))
                         for i, (u, v) in enumerate(links)]
                net = sr.network(n, edges, [n-1])
                h, m = net['h'], n-1
                counts['weighted_graphs'] += 1
                source_list = [tuple(Q(i == j) for i in range(m)) for j in range(m)]
                source_list += [tuple(Q((-1)**i*(i+1)) for i in range(m))]
                for b in source_list:
                    field = sr.solve(net, b)
                    counts['source_solutions'] += 1
                    # Edge-current summation is independent of the matrix solve.
                    for selected in subsets(range(m)):
                        require(sr.cut_flux(net, field, selected) == sum(b[i] for i in selected), 'cut flux')
                        counts['cut_flux_ledgers'] += 1
                    delta = tuple(Q(i+1) for i in range(m))+(Q(0),)
                    shifted = tuple(x+y for x, y in zip(field, delta))
                    require(sr.energy(net, shifted, b)-sr.energy(net, field, b)
                            == sr.energy(net, delta, [0]*m), 'completion of cost')
                    counts['energy_completions'] += 1
                    for retained in subsets(range(m)):
                        reduced = sr.eliminate(h, b, retained)
                        answer = sr.vector(gm.mul(gm.inverse(reduced['h']), sr.column(reduced['b'])))
                        require(answer == tuple(field[i] for i in retained), 'stationary cut response')
                        require(sr.stationary_cost(reduced['h'], reduced['b'])+reduced['constant']
                                == sr.energy(net, field, b), 'hidden source energy')
                        counts['stationary_eliminations'] += 1
                source = tuple(Q(i == 0) for i in range(2*m))
                state = tuple(Q(i-1, 3) for i in range(2*m))
                g, j = sr.native_generator(h), im.native_grading(m)
                require(gm.transpose(g) == gm.scale(g, -1), 'anti-self-dagger')
                require(gm.mul(gm.mul(j, g), j) == gm.scale(g, -1), 'cut-odd')
                after = sr.affine_step(h, source, state)
                require(sr.paired_cost(h, source, state) == sr.paired_cost(h, source, after), 'affine cost')
                centre = sr.vector(gm.mul(gm.inverse(im.kron(h, gm.identity(2))), sr.column(source)))
                require(sr.affine_step(h, source, centre) == centre, 'source fixed point')
                counts['native_affine_laws'] += 1
    for length in range(2, 10):
        for alpha in range(4):
            capacities = tuple(Q(n**alpha) for n in range(1, length+1))
            net = sr.network(length+1, [(i, i+1, c) for i, c in enumerate(capacities)], [length])
            for source in [Q(1, 3), Q(1), Q(4)]:
                formula = sr.radial(capacities, source)
                field = sr.solve(net, [source]+[0]*(length-1))
                require(field == formula['field'], 'radial formula vs matrix inverse')
                require(set(formula['flux']) == {source}, 'radial source flux')
                counts['radial_profiles'] += 1
            source, probe = [Q(2)]+[Q(0)]*(length-1), Q(3)
            energies = []
            for i in range(length):
                p = [Q(0)]*length
                p[i] = probe
                energies.append(sr.source_only_cross_energy(net['h'], source, p))
            energies.append(Q(0))  # Grounded endpoint has zero cross term.
            force = sr.source_probe_force(capacities, 2, 3)
            for i in range(length):
                require(-(energies[i+1]-energies[i]) == force[i], 'source-only cross-energy force')
                counts['pair_energy_force_checks'] += 1
    for p in [2, 3, 4]:
        for n in range(1, 13):
            for amplitude in [1, 2, 3]:
                source, c, drop = Q(amplitude**(p-1)), Q(n**(2*(p-1))), Q(amplitude, n*n)
                require(c*drop**(p-1) == source, 'nonlinear radial source balance')
                counts['nonlinear_scaling_controls'] += 1
    for ratio in [Q(1, 4), Q(4)]:
        for magnitude in [Q(1, 2), Q(1), Q(3)]:
            radius = Q(1, 1000)
            for s1 in [-1, 0, 1]:
                for s2 in [-1, 0, 1]:
                    gate = sr.scale_test(magnitude+s1*radius, radius,
                                         ratio*magnitude+s2*radius, radius, ratio)
                    require(gate['verdict'] == 'UNRESOLVED', 'compatible ratio falsely rejected')
                    counts['bounded_error_ratio_checks'] += 1
    wrong_distance = sr.scale_test(1, Q(1, 100), Q(1, 2), Q(1, 100), Q(1, 4))
    wrong_source = sr.scale_test(1, Q(1, 100), 2, Q(1, 100), 4)
    require(wrong_distance['verdict'] == wrong_source['verdict'] == 'REJECTED', 'rejection controls')
    net = sr.network(4, [(0, 1, 1), (1, 2, 4), (2, 3, 9)], [3])
    screened_h = gm.add(net['h'], gm.identity(3))
    screened_field = sr.vector(gm.mul(gm.inverse(screened_h), sr.column([1, 0, 0])))+(Q(0),)
    screened_flux = sr.currents(net, screened_field)
    require(1 > screened_flux[0] > screened_flux[1] > screened_flux[2] > 0, 'onsite flux leakage')
    return {
        'counts': counts,
        'W45_inverse_square_source_one': list(map(str, sr.source_probe_force([1, 4, 9, 16]))),
        'W46_same_native_rules_other_exponents': [0, 1, 2, 3],
        'W47_cubic_cost_quartic_capacity': 'Same 1/n^2 drop, source times 4 gives response times 2',
        'W48_onsite_flux': list(map(str, screened_flux)),
        'W49_naive_cut_stiffness': '1 versus correct effective 36/49',
        'W50_posthoc_radius': 'A 1/n profile becomes 1/r^2 under fitted r=sqrt(n); forbidden calibration',
        'W51_no_relaxation': 'Affine Cayley preserves excess cost; a nonstationary state does not settle',
        'distance_control': wrong_distance['verdict'],
        'source_control': wrong_source['verdict'],
        'data_type': 'Synthetic exact rational fixtures, not laboratory observations',
        'native_primitive_selects_quadratic_cost_or_quadratic_shell_growth': False,
        'physical_inverse_square_law_derived_without_additional_hypotheses': False,
        'physical_gravity_test_passed': False,
        'quantum_gravity_established': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(build_checks(), indent=2))
