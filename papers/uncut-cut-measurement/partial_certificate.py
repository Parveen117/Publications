"""Independent small-carrier oracles for partial continuation statements."""
from itertools import product

import finite_model as fm
import family_model as fam
import partial_model as pm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def refines(fine, coarse):
    return all(fine[x] != fine[y] or coarse[x] == coarse[y]
               for x in range(len(fine)) for y in range(len(fine)))


def closed(cut, arrows):
    # Direct pairwise gate, independent of decoder/signature construction.
    for x, y in product(range(len(cut)), repeat=2):
        if cut[x] != cut[y]:
            continue
        for t in arrows:
            u, v = t[x], t[y]
            if (u is None) != (v is None):
                return False
            if u is not None and cut[u] != cut[v]:
                return False
    return True


def distances(cut, arrows, availability):
    """Backward distance relaxation over all pairs, no partition or forward BFS."""
    n = len(cut)
    inf = n * n + 2
    d = [[0 if cut[x] != cut[y] else inf for y in range(n)] for x in range(n)]
    if availability:
        for x, y in product(range(n), repeat=2):
            if any((t[x] is None) != (t[y] is None) for t in arrows):
                d[x][y] = min(d[x][y], 1)
    for _ in range(n * n):
        following = [row[:] for row in d]
        for x, y in product(range(n), repeat=2):
            for t in arrows:
                u, v = t[x], t[y]
                if u is not None and v is not None:
                    following[x][y] = min(following[x][y], 1 + d[u][v])
        if following == d:
            return d, inf
        d = following
    raise RuntimeError("Independent pair distances did not stabilize")


def execute(arrows, word, x):
    for a in word:
        if x is None:
            return None
        x = arrows[a][x]
    return x


def verify_witness(cut, arrows, x, y, witness, distance, inf, availability):
    require((witness is None) == (distance == inf), "U15 missing/spurious witness")
    if witness is None:
        return
    require(len(witness['word']) == distance, "U15 witness is not shortest")
    u, v = x, y
    for a in witness['prefix']:
        u, v = arrows[a][u], arrows[a][v]
        require(u is not None and v is not None, "U15 prefix is not commonly admissible")
    if witness['kind'] == 'output':
        require(witness['prefix'] == witness['word'] and cut[u] != cut[v],
                "U15 output witness fails")
    else:
        require(availability and witness['kind'] == 'availability',
                "U15 safe-only search fabricated an availability observation")
        a = witness['arrow']
        require(witness['word'] == witness['prefix'] + (a,), "U15 malformed query")
        require((arrows[a][u] is None) != (arrows[a][v] is None),
                "U15 availability witness fails")


def build_checks():
    counts = {key: 0 for key in (
        'cut_two_partial_arrow_cases', 'closed_refinement_comparisons',
        'diagnostic_pair_mode_checks', 'smaller_memory_maps_rejected',
        'word_descent_checks', 'totalization_comparisons')}
    domains = []
    for n in range(1, 4):
        cuts = tuple(fm.partitions(n))
        maps = tuple(product((None,) + tuple(range(n)), repeat=n))
        domains.append({'carrier_size': n, 'cut_partitions': len(cuts),
                        'partial_endomaps': len(maps), 'ordered_arrow_pairs': len(maps) ** 2})
        for arrows in product(maps, repeat=2):
            valid = {p: closed(p, arrows) for p in cuts}
            for cut in cuts:
                result = pm.partial_closure(cut, arrows)
                h = result['cut']
                require(refines(h, cut) and valid[h], 'U14 invalid closure')
                require(result['depth'] <= n - len(set(cut)), 'U14 depth bound')
                for a, t in enumerate(arrows):
                    require((pm.descended_partial(cut, t) is not None) == closed(cut, (t,)),
                            'U13 exact-domain quotient gate mismatch')
                    require(result['transitions'][a] is not None, 'U13 final quotient missing')
                for p in cuts:
                    if refines(p, cut) and valid[p]:
                        require(refines(p, h), 'U14 closure is not coarsest')
                        counts['closed_refinement_comparisons'] += 1

                # Direct tagged word signatures, including every possible distinguishing depth.
                words = tuple(w for length in range(n - len(set(cut)) + 1)
                              for w in product(range(2), repeat=length))
                signatures = []
                for x in range(n):
                    signature = []
                    for w in words:
                        y = execute(arrows, w, x)
                        signature.append(() if y is None else (cut[y],))
                    signatures.append(tuple(signature))
                require(refines(h, signatures) and refines(signatures, h), 'U14 word oracle')

                k = fm.repair_size(cut, h)
                memory = fm.repair_channel(cut, h)
                require(closed(tuple(zip(cut, memory)), arrows), 'U14 repair does not close')
                if k > 1:
                    for m in product(range(k - 1), repeat=n):
                        require(not closed(tuple(zip(cut, m)), arrows), 'U14 smaller repair works')
                        counts['smaller_memory_maps_rejected'] += 1

                tc, ta = pm.totalize(cut, arrows)
                total = fam.family_closure(tc, ta)['cut']
                require(fm.canonical(total[:-1]) == h and total[-1] not in total[:-1],
                        'U14 tagged sink equivalence failed')
                counts['totalization_comparisons'] += 1

                for availability in (False, True):
                    d, inf = distances(cut, arrows, availability)
                    for x, y in product(range(n), repeat=2):
                        witness = pm.diagnostic(cut, arrows, x, y,
                                                observe_availability=availability)
                        verify_witness(cut, arrows, x, y, witness, d[x][y], inf, availability)
                        if availability:
                            require((witness is None) == (h[x] == h[y]), 'U15 quotient mismatch')
                            if witness is not None:
                                require(len(witness['word']) <= n - len(set(cut)), 'U15 depth')
                        elif witness is not None:
                            require(len(witness['word']) <= n * n - 1, 'U15 safe pair bound')
                        counts['diagnostic_pair_mode_checks'] += 1

                # Domain and endpoint agree after each short generated composition.
                for length in range(3):
                    for w in product(range(2), repeat=length):
                        actual = pm.word_map(n, arrows, w)
                        for x in range(n):
                            label = h[x]
                            for a in w:
                                value = result['transitions'][a][label]
                                if not value:
                                    label = None
                                    break
                                label = value[0]
                            require(label == (None if actual[x] is None else h[actual[x]]),
                                    'U13 word domain/endpoint descent failed')
                            counts['word_descent_checks'] += 1
                counts['cut_two_partial_arrow_cases'] += 1

    c7, a7 = (0, 0, 0, 1), ((0, None, 3, 3),)
    require(pm.diagnostic(c7, a7, 0, 1) is None and
            pm.diagnostic(c7, a7, 1, 2) is None and
            pm.diagnostic(c7, a7, 0, 2)['word'] == (0,), 'W7 nontransitivity')
    c8, a8 = (0, 0), ((0, None),)
    require(pm.diagnostic(c8, a8, 0, 1) is None, 'W8 invented executable result')
    require(pm.diagnostic(c8, a8, 0, 1, observe_availability=True)['kind'] == 'availability',
            'W8 missing availability distinction')
    require(fm.repair_size(c8, pm.partial_closure(c8, ((0, 1),))['cut']) == 1 and
            fm.repair_size(c8, pm.partial_closure(c8, a8)['cut']) == 2,
            'W8 domain restriction memory increase failed')
    c9, a9 = (0, 0, 0), ((1, 2, None),)
    r9 = pm.partial_closure(c9, a9)
    require(r9['depth'] == 2 and fm.repair_size(c9, r9['cut']) == 3, 'W9 delayed permission')
    require(pm.diagnostic(c9, a9, 0, 1) is None and
            pm.diagnostic(c9, a9, 0, 1, observe_availability=True)['prefix'] == (0,),
            'W9 query confused with executable output')
    return {
        'domains': domains, 'checks': counts,
        'controls': {
            'W7_safe_nontransitive': {'cut': c7, 'arrows': a7,
                                      'safe_indistinguishable_pairs': [[0, 1], [1, 2]],
                                      'distinguished_pair': [0, 2], 'word': [0]},
            'W8_availability_not_execution': {'cut': c8, 'arrows': a8,
                                               'safe_witness': None, 'availability_query': 0,
                                               'memory_total_identity': 1, 'memory_restricted': 2},
            'W9_delayed_availability': {'cut': c9, 'arrows': a9,
                                        'class_counts': [len(set(p)) for p in r9['levels']],
                                        'depth': r9['depth'], 'memory_labels': 3}},
        'scope': 'All cuts and ordered pairs of partial endomaps through three states; W7 is a named four-state control. Availability queries require a separately justified interface; no physical instrument is certified.'}
