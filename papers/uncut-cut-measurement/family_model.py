"""Finite continuation families and distinguishing experiments (v0.2)."""
from collections import deque

import finite_model as fm


def _family(cut, transitions):
    fm._same_nonempty(cut)
    result = tuple(tuple(t) for t in transitions)
    for transition in result:
        fm._transition(cut, transition)
    return result


def joint_cut(*cuts):
    fm._same_nonempty(*cuts)
    return fm.canonical(tuple(zip(*cuts)))


def closed_for_family(cut, transitions):
    family = _family(cut, transitions)
    return all(fm.descended_transition(cut, t) is not None for t in family)


def family_closure(cut, transitions):
    """Coarsest refining observation on which every admitted arrow descends.

    Empty families are allowed: only the initial observation is required.
    Arrow order in the family names the generators, not a sequence of execution.
    """
    family = _family(cut, transitions)
    current = fm.canonical(cut)
    levels = [current]
    for _ in range(len(cut)):
        signatures = tuple((current[x],) + tuple(current[t[x]] for t in family)
                           for x in range(len(cut)))
        following = fm.canonical(signatures)
        if following == current:
            return {"cut": current, "depth": len(levels) - 1,
                    "levels": tuple(levels),
                    "transitions": tuple(fm.descended_transition(current, t)
                                         for t in family)}
        current = following
        levels.append(current)
    raise RuntimeError("Finite family refinement did not stabilize")


def word_map(n, transitions, word):
    """Compose in execution order: word (a,b) means T_b composed with T_a."""
    if type(n) is not int or n < 1:
        raise ValueError("Carrier size must be a positive integer")
    family = _family(tuple(range(n)), transitions)
    positions = tuple(range(n))
    for a in word:
        if type(a) is not int or a < 0 or a >= len(family):
            raise ValueError("Word contains an undeclared arrow")
        positions = tuple(family[a][x] for x in positions)
    return positions


def distinguishing_word(cut, transitions, x, y):
    """Shortest word that separates the initial candidates, or None.

    Breadth-first search uses at most |X|^2 ordered pairs. This is a model
    prediction of an experiment, not access to unmeasured future observations.
    """
    family = _family(cut, transitions)
    if any(type(i) is not int or i < 0 or i >= len(cut) for i in (x, y)):
        raise ValueError("Candidates must belong to the presented carrier")
    queue = deque([(x, y, ())])
    seen = {(x, y)}
    while queue:
        u, v, word = queue.popleft()
        if cut[u] != cut[v]:
            return word
        for a, transition in enumerate(family):
            pair = (transition[u], transition[v])
            if pair not in seen:
                seen.add(pair)
                queue.append((*pair, word + (a,)))
    return None
