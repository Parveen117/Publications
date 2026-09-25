"""GP-1/GP-2 exact coherent inverse and phase-erased no-go controls."""

from fractions import Fraction as F

if __package__:
    from . import quantum_classical_bridge as qc
else:
    import quantum_classical_bridge as qc


def sub(a, b):
    return qc.add(a, qc.neg(b))


def div(a, b):
    n = qc.norm2(b)
    assert n > 0
    return tuple(x / n for x in qc.mul(a, qc.conj(b)))


def inv2(a):
    det = sub(qc.mul(a[0][0], a[1][1]), qc.mul(a[0][1], a[1][0]))
    return [[div(a[1][1], det), div(qc.neg(a[0][1]), det)],
            [div(qc.neg(a[1][0]), det), div(a[0][0], det)]]


def cayley_checks(h, u, green, step=F(1)):
    lhs = [[sub(qc.ONE if i == j else qc.ZERO, qc.mul(qc.I, qc.c(step * h[i][j])))
            for j in range(2)] for i in range(2)]
    rhs = [[qc.add(qc.ONE if i == j else qc.ZERO, qc.mul(qc.I, qc.c(step * h[i][j])))
            for j in range(2)] for i in range(2)]
    assert qc.mm(lhs, u) == rhs
    identity = [[qc.ONE, qc.ZERO], [qc.ZERO, qc.ONE]]
    assert qc.mm(qc.dagger(u), u) == identity
    assert [[sum(h[i][k] * green[k][j] for k in range(2))
             for j in range(2)] for i in range(2)] == [[F(1), F(0)], [F(0), F(1)]]

    plus = [[qc.add(u[i][j], qc.ONE if i == j else qc.ZERO)
             for j in range(2)] for i in range(2)]
    minus = [[sub(u[i][j], qc.ONE if i == j else qc.ZERO)
              for j in range(2)] for i in range(2)]
    recovered = [[qc.mul(qc.c(0, step), q) for q in row]
                 for row in qc.mm(plus, inv2(minus))]
    assert recovered == [[qc.c(green[i][j]) for j in range(2)] for i in range(2)]
    return [[qc.norm2(u[i][j]) for j in range(2)] for i in range(2)]


def run():
    graph_a = [[F(2), F(-1)], [F(-1), F(2)]]
    graph_b = [[F(11, 12), F(-5, 12)], [F(-5, 12), F(11, 12)]]
    assert [sum(r) for r in graph_a] == [F(1), F(1)]
    assert [sum(r) for r in graph_b] == [F(1, 2), F(1, 2)]
    assert -graph_a[0][1] == F(1) and -graph_b[0][1] == F(5, 12)
    green_a = [[F(2, 3), F(1, 3)], [F(1, 3), F(2, 3)]]
    green_b = [[F(11, 8), F(5, 8)], [F(5, 8), F(11, 8)]]
    ua = [[(F(-2, 5), F(4, 5)), (F(2, 5), F(1, 5))],
          [(F(2, 5), F(1, 5)), (F(-2, 5), F(4, 5))]]
    ub = [[(F(4, 25), F(22, 25)), (F(11, 25), F(-2, 25))],
          [(F(11, 25), F(-2, 25)), (F(4, 25), F(22, 25))]]
    pa = cayley_checks(graph_a, ua, green_a)
    pb = cayley_checks(graph_b, ub, green_b)
    assert ua != ub
    phase = (F(4, 5), F(-3, 5))
    assert qc.norm2(phase) == F(1)
    assert ub == [[qc.mul(phase, ua[i][j]) for j in range(2)] for i in range(2)]
    for psi in ([qc.ONE, qc.ZERO], [qc.ONE, qc.ONE], [qc.ONE, qc.I]):
        rho = qc.density(psi)
        assert qc.evolve(ua, rho) == qc.evolve(ub, rho)
    assert pa == pb == [[F(4, 5), F(1, 5)], [F(1, 5), F(4, 5)]]
    assert green_a[0][1] == F(1, 3) != green_b[0][1] == F(5, 8)

    source = [F(1), F(2)]
    field_a = [sum(green_a[i][j] * source[j] for j in range(2)) for i in range(2)]
    field_b = [sum(green_b[i][j] * source[j] for j in range(2)) for i in range(2)]
    assert field_a == [F(4, 3), F(5, 3)]
    assert field_b == [F(21, 8), F(27, 8)]
    assert -green_a[0][1] == -F(1, 3) and -green_b[0][1] == -F(5, 8)

    # A coherent map fixes h*H, not H when native h is unknown.
    scale = F(2)
    assert [[scale * x for x in row] for row in graph_a] != graph_a
    scaled_graph = [[scale * x for x in row] for row in graph_a]
    scaled_green = [[x / scale for x in row] for row in green_a]
    assert cayley_checks(scaled_graph, ua, scaled_green, F(1, 2)) == pa
    assert scaled_green != green_a

    return {
        "same_dephased_transition": [[str(x) for x in row] for row in pa],
        "distinct_coherent_maps": True,
        "same_quantum_density_channel_up_to_global_phase": True,
        "global_phase_relation": [str(x) for x in phase],
        "cross_response_a": str(-green_a[0][1]),
        "cross_response_b": str(-green_b[0][1]),
        "source_field_a": [str(x) for x in field_a],
        "source_field_b": [str(x) for x in field_b],
        "full_coherent_cayley_inversion_passed": True,
        "unknown_step_scale_ambiguous": True,
        "universal_gravitational_coupling_derived": False,
    }


if __name__ == "__main__":
    print("PASS_EXACT_GRAVITY_PHASE_TRANSLATOR", run())
