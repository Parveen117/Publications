"""TF-1/TF-2 exact rational witness; the general proof is in the manuscript."""

from fractions import Fraction as F


def schur_one_hidden(k, b):
    zz = k[2][2]
    assert zz > 0
    d = [[k[i][j] - k[i][2] * k[2][j] / zz for j in range(2)]
         for i in range(2)]
    source = [b[i] - k[i][2] * b[2] / zz for i in range(2)]
    return d, source


def run():
    k = [[F(3), F(-1), F(-1)],
         [F(-1), F(4), F(-1)],
         [F(-1), F(-1), F(2)]]
    b = [F(0), F(10), F(0)]
    # Triangle unit links plus independent positive grounding on nodes 0,1.
    for i in range(3):
        assert all(k[i][j] == -1 for j in range(3) if i != j)
    assert [sum(row) for row in k] == [F(1), F(2), F(0)]
    d, source = schur_one_hidden(k, b)
    assert d == [[F(5, 2), F(-3, 2)], [F(-3, 2), F(7, 2)]]
    assert source == [F(0), F(10)]
    determinant = d[0][0] * d[1][1] - d[0][1] ** 2
    assert k[0][0] > 0 and k[0][0] * k[1][1] - k[0][1] ** 2 > 0
    assert determinant == F(13, 2) and k[2][2] * determinant == F(13)

    # The hidden state is exactly the stationary minimizer at each fixed y.
    for s, v in [(F(1), F(1)), (F(2), F(1)), (F(1, 2), F(3, 2))]:
        z = -(k[2][0] * s + k[2][1] * v - b[2]) / k[2][2]
        y = [s, v]
        full = [s, v, z]
        full_e = sum(full[i] * k[i][j] * full[j] / 2
                     for i in range(3) for j in range(3)) - sum(
                         b[i] * full[i] for i in range(3))
        projected_e = sum(y[i] * d[i][j] * y[j] / 2
                          for i in range(2) for j in range(2)) - sum(
                              source[i] * y[i] for i in range(2))
        assert full_e == projected_e

    # Basis change z_new=2*z cannot alter visible equilibrium energy.
    scale = [F(1), F(1), F(1, 2)]
    changed_k = [[scale[i] * k[i][j] * scale[j] for j in range(3)]
                 for i in range(3)]
    changed_b = [scale[i] * b[i] for i in range(3)]
    assert schur_one_hidden(changed_k, changed_b) == (d, source)

    s = v = F(1)
    t = d[0][0] * s + d[0][1] * v - source[0]
    p = source[1] - d[1][0] * s - d[1][1] * v
    assert (t, p) == (F(1), F(8))
    assert d[0][1] == d[1][0]  # Maxwell T_V=-P_S, as P_S=-D_VS.
    cv = t / d[0][0]
    cp = t / (d[0][0] - d[0][1] ** 2 / d[1][1])
    ks = v * d[1][1]
    kt = v * (d[1][1] - d[0][1] ** 2 / d[0][0])
    assert (cv, cp, ks, kt) == (F(2, 5), F(7, 13), F(7, 2), F(13, 5))
    assert cp / cv == ks / kt == F(35, 26)

    # A measured pressure contaminated by +epsilon*S breaks exact closure.
    epsilon = F(1, 7)
    assert d[0][1] != -( -d[1][0] + epsilon)
    return {
        "graph_vertices": 3,
        "hidden_coordinates": 1,
        "effective_hessian_determinant": "13/2",
        "temperature_at_unit_chart": str(t),
        "pressure_at_unit_chart": str(p),
        "response_ratio": str(cp / cv),
        "maxwell_perturbation_rejected": True,
        "hidden_basis_invariant": True,
        "physical_entropy_volume_identification_derived": False,
    }


if __name__ == "__main__":
    print("PASS_EXACT_WITNESS TF-1 TF-2; Maxwell negative control rejected", run())
