"""Exhaustive v0.2 evidence, kept separate from the v0.1 domains."""
from itertools import product

import finite_model as fm
import family_model as fam


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def implies(fine, coarse):
    return all(fine[x] != fine[y] or coarse[x] == coarse[y]
               for x in range(len(fine)) for y in range(len(fine)))


def brute_closed(cut, family):
    return all(cut[x] != cut[y] or all(cut[t[x]] == cut[t[y]] for t in family)
               for x in range(len(cut)) for y in range(len(cut)))


def enumerated_words(n, family):
    """Direct execution oracle, independent of refinement and pair-graph search."""
    result = []
    for depth in range(n + 1):
        for word in product(range(len(family)), repeat=depth):
            images = []
            for x in range(n):
                current = x
                for a in word:
                    current = family[a][current]
                images.append(current)
            result.append((word, tuple(images)))
    return tuple(result)


def check_case(cut, family, candidates, words, counts):
    n = len(cut)
    result = fam.family_closure(cut, family)
    h = result['cut']
    require(result['depth'] <= n - len(set(cut)), 'U8 depth bound')
    require(implies(h, cut) and brute_closed(h, family), 'U8 invalid closure')
    signatures = tuple(tuple(cut[images[x]] for _, images in words) for x in range(n))
    oracle = fm.canonical(signatures)
    require(h == oracle, 'U8 direct word oracle disagrees')
    for d in candidates:
        if implies(d, cut) and brute_closed(d, family):
            require(implies(d, h), 'U8 closure is not coarsest')
            counts['closed_refinement_comparisons'] += 1
    require(fam.family_closure(h, family)['cut'] == h, 'U11 non-idempotent closure')

    k = fm.repair_size(cut, h)
    memory = fm.repair_channel(cut, h)
    joint = fam.joint_cut(cut, memory)
    require(len(set(memory)) == k and brute_closed(joint, family),
            'U9 attaining repair failed')
    require(implies(joint, h) and implies(h, joint), 'U9 repair changed partition')
    if k > 1:
        for smaller in product(range(k - 1), repeat=n):
            require(not brute_closed(tuple(zip(cut, smaller)), family),
                    'U9 smaller memory closes the family')
            counts['smaller_memory_maps_rejected'] += 1

    for subset in ((),) + tuple((t,) for t in family):
        weaker = fam.family_closure(cut, subset)['cut']
        require(implies(h, weaker) and fm.repair_size(cut, weaker) <= k,
                'U9 added-arrow monotonicity failed')
        if len(subset) == 1:
            require(weaker == fm.stable_history(cut, subset[0])['cut'],
                    'U8 single-arrow specialization disagrees with U6')
        counts['subfamily_comparisons'] += 1

    redundant = tuple(range(n))
    if family:
        # Apply first arrow, then last arrow, by direct array lookup.
        redundant = tuple(family[-1][family[0][x]] for x in range(n))
    enlarged = family + (tuple(range(n)), redundant)
    require(fam.family_closure(cut, enlarged)['cut'] == h,
            'U10 identity/composite changed closure')
    require(fam.family_closure(cut, tuple(reversed(family)))['cut'] == h,
            'U10 generator-list order changed closure')
    counts['generator_presentations_checked'] += 2

    for word, images in words:
        require(fam.word_map(n, family, word) == images, 'Word execution mismatch')
        for x in range(n):
            label = h[x]
            for a in word:
                label = result['transitions'][a][label]
            require(label == h[images[x]], 'U10 observed word composition failed')
        counts['word_descent_checks'] += 1

    for x, y in product(range(n), repeat=2):
        expected = next((w for w, images in words if cut[images[x]] != cut[images[y]]), None)
        actual = fam.distinguishing_word(cut, family, x, y)
        require(actual == expected, 'U12 shortest-word oracle disagrees')
        require((actual is None) == (h[x] == h[y]), 'U12 false indistinguishability')
        if actual is not None:
            require(len(actual) <= n - len(set(cut)), 'U12 length bound failed')
        counts['diagnostic_pairs_checked'] += 1
    counts['cut_family_pairs_checked'] += 1
    return result


def check_joint_cuts(cuts, closures, counts):
    for c, d in product(cuts, repeat=2):
        joined = fam.joint_cut(c, d)
        hc, hd = closures[c], closures[d]
        require(closures[joined]['cut'] == fam.joint_cut(hc['cut'], hd['cut']),
                'U11 closure of joint cuts failed')
        counts['joint_cut_comparisons'] += 1
        if implies(d, c):
            require(implies(hd['cut'], hc['cut']), 'U11 refinement monotonicity failed')
            r = fm.factor_map(hd['cut'], hc['cut'])
            for td, tc in zip(hd['transitions'], hc['transitions']):
                require(all(r[td[b]] == tc[r[b]] for b in set(hd['cut'])),
                        'U11 refinement transition diagram failed')
            counts['refinement_intertwiners_checked'] += 1


def build_checks():
    counts = {
        'exhaustive_cut_family_pairs': 0,
        'cut_family_pairs_checked': 0,
        'four_state_fixtures': 0,
        'closed_refinement_comparisons': 0,
        'smaller_memory_maps_rejected': 0,
        'subfamily_comparisons': 0,
        'generator_presentations_checked': 0,
        'word_descent_checks': 0,
        'diagnostic_pairs_checked': 0,
        'joint_cut_comparisons': 0,
        'refinement_intertwiners_checked': 0,
    }
    domains = []
    for n in range(1, 4):
        cuts = tuple(fm.partitions(n))
        arrows = tuple(fm.all_endomaps(n))
        domains.append({'carrier_size': n, 'cut_partitions': len(cuts),
                        'ordered_two_arrow_families': len(arrows) ** 2,
                        'word_oracle_max_depth': n})
        for family in product(arrows, repeat=2):
            words = enumerated_words(n, family)
            closures = {c: check_case(c, family, cuts, words, counts) for c in cuts}
            check_joint_cuts(cuts, closures, counts)
            counts['exhaustive_cut_family_pairs'] += len(cuts)

    # Explicit four-state checks supplement, but do not enlarge, the exhaustive domain.
    cut, A, B = (0, 0, 0, 1), (0, 2, 2, 3), (0, 0, 3, 3)
    cuts4 = tuple(fm.partitions(4))
    for family in ((A, B), ((1, 2, 3, 3),)):
        check_case(cut, family, cuts4, enumerated_words(4, family), counts)
        counts['four_state_fixtures'] += 1
    identities = (tuple(range(4)), tuple(range(4)))
    words4 = enumerated_words(4, identities)
    closures4 = {c: check_case(c, identities, cuts4, words4, counts) for c in cuts4}
    counts['four_state_fixtures'] += len(cuts4)
    check_joint_cuts(cuts4, closures4, counts)

    ha = fam.family_closure(cut, (A,))['cut']
    hb = fam.family_closure(cut, (B,))['cut']
    both = fam.family_closure(cut, (A, B))
    joined = fam.joint_cut(ha, hb)
    costs = [fm.repair_size(cut, h) for h in (ha, hb, both['cut'])]
    require(costs == [1, 2, 3] and costs[2] > costs[0] * costs[1],
            'W5 separate-memory product control failed')
    require(not brute_closed(joined, (A, B)), 'W5 naive join falsely closes')
    word = fam.distinguishing_word(cut, (A, B), 0, 1)
    require(word == (0, 1), 'W5 shortest diagnostic word changed')
    require(all(fam.distinguishing_word(cut, (t,), 0, 1) is None for t in (A, B)),
            'W5 an unmixed continuation unexpectedly separates')
    ab = tuple(B[A[x]] for x in range(4))
    expanded = fam.family_closure(cut, (A, B, ab))
    require(expanded['cut'] == both['cut'] and expanded['depth'] == 1 and
            both['depth'] == 2, 'W5 generator-dependent depth control failed')

    coarse, fine, transition = (0, 0, 0), (0, 0, 1), (0, 2, 2)
    hc = fam.family_closure(coarse, (transition,))['cut']
    hd = fam.family_closure(fine, (transition,))['cut']
    finer_costs = [fm.repair_size(coarse, hc), fm.repair_size(fine, hd)]
    require(finer_costs == [1, 2], 'W6 initial-cut cost control failed')

    return {
        'status': 'PASS_FINITE_CHECKS',
        'exhaustive_domains': domains,
        'additional_four_state_scope': 'W5; the U6 sharp-depth chain; identity pair on all 15 cuts',
        'counts': counts,
        'controls': {
            'W5_mixed_words': {
                'initial_cut': cut, 'arrows': {'A': A, 'B': B},
                'single_closures': [ha, hb], 'naive_join': joined,
                'family_closure': both['cut'],
                'class_counts': [len(set(p)) for p in both['levels']],
                'memory_sizes_A_B_family': costs,
                'family_binary_bits': fm.fixed_binary_bits(costs[2]),
                'diagnostic_candidates': [0, 1], 'shortest_word': ['A', 'B'],
                'unmixed_words_cannot_distinguish': True,
                'depth_before_after_naming_composite': [both['depth'], expanded['depth']],
            },
            'W6_initial_detail_can_increase_prediction_memory': {
                'coarse_cut': coarse, 'fine_cut': fine, 'transition': transition,
                'closures': [hc, hd], 'memory_sizes': finer_costs,
            },
        },
        'scope': {
            'general_proofs': 'FAMILY_CONTINUATION.md U8-U12 under declared finite hypotheses',
            'all_finite_words_admitted': True,
            'arrow_labels_known_to_observer': True,
            'physical_memory_acquisition_derived': False,
            'minimum_physical_action_derived': False,
            'quantum_or_spacetime_identification': False,
        },
    }
