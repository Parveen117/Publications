"""Exact finite controls for the conditional RKF coupling reduction.

This file has no experimental alpha input and predicts no physical constant.
The witness matrix is a declared arithmetic fixture, not a selected vacuum.
"""
from fractions import Fraction as Q
import json


def tr(a):
    return [list(r) for r in zip(*a)]


def mul(a, b):
    return [[sum((x*y for x, y in zip(row, col)), Q(0))
             for col in zip(*b)] for row in a]


def inv(a):
    n = len(a)
    m = [[Q(x) for x in row] + [Q(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for c in range(n):
        pivot = next(i for i in range(c, n) if m[i][c])
        m[c], m[pivot] = m[pivot], m[c]
        d = m[c][c]
        m[c] = [x/d for x in m[c]]
        for i in range(n):
            if i != c:
                d = m[i][c]
                m[i] = [x-d*y for x, y in zip(m[i], m[c])]
    return [r[n:] for r in m]


def schur(a):
    b = [a[0][1:]]
    c = [row[1:] for row in a[1:]]
    return Q(a[0][0]) - mul(mul(b, inv(c)), tr(b))[0][0]


H = [[Q(5), Q(1), Q(1)],
     [Q(1), Q(2), Q(0)],
     [Q(1), Q(0), Q(3)]]
q = Q(2)
Z = schur(H)
gamma = q*q/Z
checks = {}

assert Z == Q(25, 6) and gamma == Q(24, 25)
checks['schur_coefficient'] = str(Z)

# Independent solution of the full sourced system.
J = [[q], [Q(0)], [Q(0)]]
state = mul(inv(H), J)
assert state == [[Q(12, 25)], [Q(-6, 25)], [Q(-4, 25)]]
assert mul(tr(J), state)[0][0] == gamma
checks['full_source_solution_agrees'] = True

# Completion of squares over a finite exact audit lattice.
n = 0
for x in map(Q, range(-2, 3)):
    for h1 in map(Q, range(-2, 3)):
        for h2 in map(Q, range(-2, 3)):
            v = [[x], [h1], [h2]]
            full = mul(mul(tr(v), H), v)[0][0]
            reduced = Z*x*x + 2*(h1+x/2)**2 + 3*(h2+x/3)**2
            assert full == reduced
            n += 1
checks['completion_of_squares_cases'] = n

# Arbitrary invertible hidden-coordinate change h = T h_new.
T = [[Q(1), Q(1)], [Q(0), Q(2)]]
P = [[Q(1), Q(0), Q(0)],
     [Q(0), T[0][0], T[0][1]],
     [Q(0), T[1][0], T[1][1]]]
assert schur(mul(mul(tr(P), H), P)) == Z
checks['hidden_basis_invariance'] = True

# Visible coordinate change x_new = t*x transports both source and cost.
for t in (Q(-3), Q(1, 2), Q(2), Q(5, 3)):
    assert (q/t)**2 / (Z/t**2) == gamma
checks['visible_normalization_cases'] = 4

# Eliminate h1 then h2; compare with simultaneous elimination.
assert Q(5) - Q(1, 2) - Q(1, 3) == Z
checks['sequential_elimination_agrees'] = True

# Sensitivity controls: deleting memory and changing physical stiffness
# must change the coupling, unlike a mere coordinate change.
assert q*q/Q(5) != gamma
assert q*q/schur([[2*x for x in row] for row in H]) == gamma/2
checks['memory_deletion_detected'] = True
checks['physical_stiffness_change_detected'] = True

print(json.dumps({
    'status': 'PASS_EXACT_FINITE_COUPLING_CONTROLS',
    'Z_eff': str(Z),
    'gamma_q_squared_over_Z_eff': str(gamma),
    'checks': checks,
    'alpha_physical_status': 'NOT_DERIVED',
    'matrix_status': 'DECLARED_FIXTURE_NOT_NATIVE_VACUUM_SELECTION'
}, indent=2))
