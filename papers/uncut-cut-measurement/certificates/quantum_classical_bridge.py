"""Exact Fraction oracle for QB-1/QB-3 on one grounded native graph."""

from fractions import Fraction as F


def c(re=0, im=0):
    return F(re), F(im)


ZERO, ONE, I = c(), c(1), c(0, 1)


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def neg(a):
    return -a[0], -a[1]


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conj(a):
    return a[0], -a[1]


def norm2(a):
    return mul(a, conj(a))[0]


def mm(a, b):
    return [[sum_c(mul(a[i][k], b[k][j]) for k in range(2))
             for j in range(2)] for i in range(2)]


def sum_c(values):
    out = ZERO
    for v in values:
        out = add(out, v)
    return out


def dagger(a):
    return [[conj(a[j][i]) for j in range(2)] for i in range(2)]


def mv(a, v):
    return [sum_c(mul(a[i][j], v[j]) for j in range(2)) for i in range(2)]


def density(v):
    n = sum(norm2(a) for a in v)
    assert n > 0
    return [[tuple(q / n for q in mul(v[i], conj(v[j])))
             for j in range(2)] for i in range(2)]


def diagonal(a):
    return [[a[i][j] if i == j else ZERO for j in range(2)]
            for i in range(2)]


def evolve(u, rho):
    return mm(mm(u, rho), dagger(u))


def run():
    # The same two-site graph has one unit link and one unit grounding per site.
    h = [[F(2), F(-1)], [F(-1), F(2)]]
    assert [sum(row) for row in h] == [F(1), F(1)]
    assert h[0][0] * h[1][1] - h[0][1] ** 2 == F(3)
    source = [F(1), F(2)]
    star = [F(4, 3), F(5, 3)]
    assert all(sum(h[i][j] * star[j] for j in range(2)) == source[i]
               for i in range(2))
    s, v = F(2), F(1)
    t = h[0][0] * s + h[0][1] * v - source[0]
    p = source[1] - h[1][0] * s - h[1][1] * v
    cv = t / h[0][0]
    cp = t / (h[0][0] - h[0][1] ** 2 / h[1][1])
    ks = v * h[1][1]
    kt = v * (h[1][1] - h[0][1] ** 2 / h[0][0])
    assert (t, p, cv, cp, ks, kt) == (F(2), F(2), F(1), F(4, 3), F(2), F(3, 2))
    assert cp / cv == ks / kt == F(4, 3)

    a, b = (F(-2, 5), F(4, 5)), (F(2, 5), F(1, 5))
    u = [[a, b], [b, a]]
    left = [[add(ONE if i == j else ZERO, neg(mul(I, c(h[i][j]))))
             for j in range(2)] for i in range(2)]
    right = [[add(ONE if i == j else ZERO, mul(I, c(h[i][j])))
              for j in range(2)] for i in range(2)]
    identity = [[ONE, ZERO], [ZERO, ONE]]
    assert mm(left, u) == right  # Native Cayley defining equation at h=1.
    assert mm(dagger(u), u) == identity
    prob = [[norm2(u[i][j]) for j in range(2)] for i in range(2)]
    assert prob == [[F(4, 5), F(1, 5)], [F(1, 5), F(4, 5)]]
    assert [sum(row) for row in prob] == [F(1), F(1)]
    assert [sum(prob[i][j] for i in range(2)) for j in range(2)] == [F(1), F(1)]

    psi0 = [ONE, ZERO]
    rho0 = density(psi0)
    psi1, psi2 = mv(u, psi0), mv(u, mv(u, psi0))
    assert sum(norm2(x) for x in psi1) == sum(norm2(x) for x in psi2) == F(1)
    rho1 = evolve(u, rho0)
    assert rho1 == density(psi1)
    memory = [[add(rho1[i][j], neg(diagonal(rho1)[i][j]))
               for j in range(2)] for i in range(2)]
    coherent = diagonal(evolve(u, rho1))[1][1][0]
    erased = diagonal(evolve(u, diagonal(rho1)))[1][1][0]
    classical_markov = sum(prob[1][j] * prob[j][0] for j in range(2))
    residue = diagonal(evolve(u, memory))[1][1][0]
    assert coherent == norm2(psi2[1]) == F(16, 25)
    assert erased == classical_markov == F(8, 25)
    assert coherent - erased == residue == F(8, 25)

    # Same first-step populations, opposite relative phase: no closed
    # autonomous update on site probabilities can recover the next readout.
    flipped = [psi1[0], neg(psi1[1])]
    assert [norm2(x) for x in flipped] == [norm2(x) for x in psi1]
    assert norm2(mv(u, flipped)[1]) == F(0) != coherent

    return {
        "grounded_graph_vertices": 2,
        "source_response": [str(x) for x in star],
        "thermo_ratio": str(cp / cv),
        "one_step_transition": [[str(x) for x in row] for row in prob],
        "coherent_site_2_after_two_steps": str(coherent),
        "phase_erased_site_2_after_two_steps": str(erased),
        "seam_residue": str(residue),
        "same_diagonal_different_future_control": True,
        "physical_detector_frequency_law_derived": False,
        "physical_clock_and_hbar_derived": False,
    }


if __name__ == "__main__":
    print("PASS_EXACT_QUANTUM_CLASSICAL_BRIDGE", run())
