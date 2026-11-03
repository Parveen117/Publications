"""Partial continuation: domain-aware quotients and common-domain experiments.

None is an undefined *transition*, never a readout value or a physical result.
Descended arrows use () for disabled and (next_label,) for enabled.
"""
from collections import deque

import finite_model as fm


def _family(cut, transitions):
    fm._same_nonempty(cut)
    family = tuple(tuple(t) for t in transitions)
    for t in family:
        fm._same_nonempty(cut, t)
        if any(y is not None and (type(y) is not int or not 0 <= y < len(cut))
               for y in t):
            raise ValueError("An arrow value must be a carrier index or None")
    return family


def _candidate(cut, x):
    if type(x) is not int or not 0 <= x < len(cut):
        raise ValueError("Candidate must belong to the presented carrier")


def _word(family, word):
    result = tuple(word)
    if any(type(a) is not int or not 0 <= a < len(family) for a in result):
        raise ValueError("Word contains an undeclared arrow")
    return result


def descended_partial(cut, transition):
    """Exact domain and successor descent; None means no quotient exists."""
    (t,) = _family(cut, (transition,))
    target = tuple(() if y is None else (cut[y],) for y in t)
    return fm.factor_map(cut, target)


def partial_closure(cut, transitions):
    """Coarsest refinement retaining both domains and all future readouts."""
    family = _family(cut, transitions)
    current = fm.canonical(cut)
    levels = [current]
    for _ in range(len(cut)):
        signatures = tuple(
            (current[x],) + tuple(() if t[x] is None else (current[t[x]],)
                                  for t in family)
            for x in range(len(cut)))
        following = fm.canonical(signatures)
        if following == current:
            return {"cut": current, "depth": len(levels) - 1,
                    "levels": tuple(levels),
                    "transitions": tuple(descended_partial(current, t) for t in family)}
        current = following
        levels.append(current)
    raise RuntimeError("Finite partial refinement did not stabilize")


def word_map(n, transitions, word):
    """Execute in written order; propagate inapplicability without execution."""
    if type(n) is not int or n < 1:
        raise ValueError("Carrier size must be a positive integer")
    family = _family(tuple(range(n)), transitions)
    word = _word(family, word)
    positions = tuple(range(n))
    for a in word:
        positions = tuple(None if x is None else family[a][x] for x in positions)
    return positions


def totalize(cut, transitions):
    """Mathematical comparison only: fresh tagged output for an absorbing sink."""
    family = _family(cut, transitions)
    n = len(cut)
    outputs = tuple((1, value) for value in cut) + ((0,),)
    arrows = tuple(tuple(n if y is None else y for y in t) + (n,)
                   for t in family)
    return outputs, arrows


def diagnostic(cut, transitions, x, y, *, observe_availability=False):
    """Shortest output experiment, or explicitly tagged availability witness.

    Default mode only traverses arrows enabled for both candidates. Availability
    mode also allows a final domain query after a common admissible prefix. It
    does NOT authorize execution of the final arrow at a disabled candidate.
    Length counts arrow positions, including a terminal availability query.
    """
    family = _family(cut, transitions)
    _candidate(cut, x)
    _candidate(cut, y)
    queue = deque([(x, y, ())])
    seen = {(x, y)}
    while queue:
        u, v, word = queue.popleft()
        if u is None or v is None:
            # Only asymmetric undefined pairs are queued, at their true BFS depth.
            return {"kind": "availability", "word": word,
                    "prefix": word[:-1], "arrow": word[-1]}
        if cut[u] != cut[v]:
            return {"kind": "output", "word": word, "prefix": word}
        for a, t in enumerate(family):
            pair = (t[u], t[v])
            if None in pair:
                if not observe_availability or pair == (None, None):
                    continue
            if pair not in seen:
                seen.add(pair)
                queue.append((*pair, word + (a,)))
    return None
