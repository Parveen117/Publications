"""Exact finite witness for the represented-iota -> EMK local-geometry bridge.

This certificate proves a typed finite adapter on one rational response plane.
It does not identify the adapter with physical gravity, spacetime curvature,
Planck-scale dynamics, or a physical double helix.
"""
from fractions import Fraction as F


def matrix(rows):
    rows = tuple(tuple(F(x) for x in row) for row in rows)
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("nonempty rectangular matrix required")
    return rows


def transpose(a):
    return tuple(zip(*matrix(a)))


def mul(a, b):
    a, b = matrix(a), matrix(b)
    if len(a[0]) != len(b):
        raise ValueError("incompatible matrix product")
    bt = transpose(b)
    return tuple(tuple(sum(x*y for x, y in zip(row, col)) for col in bt)
                 for row in a)


def scale(a, q):
    q = F(q)
    return tuple(tuple(q*x for x in row) for row in matrix(a))


def identity(n):
    return matrix([[int(i == j) for j in range(n)] for i in range(n)])


def congruence(y, h):
    return mul(mul(transpose(y), h), y)


def det2(a):
    a = matrix(a)
    if (len(a), len(a[0])) != (2, 2):
        raise ValueError("2x2 determinant only")
    return a[0][0]*a[1][1] - a[0][1]*a[1][0]


def rank(a):
    a = [list(row) for row in matrix(a)]
    row = 0
    for col in range(len(a[0])):
        pivot = next((k for k in range(row, len(a)) if a[k][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        d = a[row][col]
        a[row] = [x/d for x in a[row]]
        for k in range(len(a)):
            if k != row and a[k][col]:
                q = a[k][col]
                a[k] = [x-q*y for x, y in zip(a[k], a[row])]
        row += 1
        if row == len(a):
            break
    return row


def response_tensor(y, z):
    """Q = G + i A in exact real/imaginary storage."""
    y, z = matrix(y), matrix(z)
    h = identity(len(y))
    assert mul(z, z) == scale(h, -1)
    assert transpose(z) == scale(z, -1)
    g = congruence(y, h)
    a = scale(congruence(y, z), -1)
    return g, a


def plane_invariant(y, z):
    g, a = response_tensor(y, z)
    assert transpose(g) == g
    assert transpose(a) == scale(a, -1)
    det_g, det_a = det2(g), det2(a)
    if det_g <= 0:
        raise ValueError("nondegenerate positive response plane required")
    gram_margin = det_g - det_a
    assert det_a >= 0 and gram_margin >= 0
    return {
        "G": g,
        "A": a,
        "det_G": det_g,
        "det_A": det_a,
        "chi": det_a/det_g,
        "gram_margin": gram_margin,
    }


def emk_local_two_jet(chi, length_squared):
    """Normalized reflection-symmetric seam two-jet from a supplied scale."""
    chi, length_squared = F(chi), F(length_squared)
    if not (F(0) <= chi <= F(1)) or length_squared <= 0:
        raise ValueError("admissible chi and positive length^2 required")
    kappa = chi/length_squared
    w0, w1, w2 = F(1), F(0), 2*kappa
    curvature = w1*w1/(4*w0*w0) - w2/(2*w0)
    assert curvature == -kappa
    return {
        "length_squared": length_squared,
        "kappa": kappa,
        "W_0": w0,
        "W_1": w1,
        "W_2": w2,
        "K_0": curvature,
    }


def strings(a):
    if isinstance(a, tuple):
        return [strings(x) for x in a]
    if isinstance(a, F):
        return str(a)
    return a


def run():
    z = matrix([
        [0, -1, 0, 0],
        [1,  0, 0, 0],
        [0,  0, 0, -1],
        [0,  0, 1,  0],
    ])
    i4 = identity(4)
    assert mul(z, z) == scale(i4, -1)
    assert transpose(z) == scale(z, -1)

    y = matrix([
        [1, 0],
        [0, 1],
        [0, 1],
        [0, 0],
    ])
    base = plane_invariant(y, z)
    assert base["G"] == matrix([[1, 0], [0, 2]])
    assert base["A"] == matrix([[0, 1], [-1, 0]])
    assert base["chi"] == F(1, 2)
    assert base["gram_margin"] == 1

    b = matrix([[2, 1], [0, 1]])
    relabelled = plane_invariant(mul(y, b), z)
    assert det2(b) == 2
    assert relabelled["det_G"] == 4*base["det_G"]
    assert relabelled["det_A"] == 4*base["det_A"]
    assert relabelled["chi"] == base["chi"]

    reverse = matrix([[1, 0], [0, -1]])
    reversed_plane = plane_invariant(mul(y, reverse), z)
    assert reversed_plane["G"] == base["G"]
    assert reversed_plane["A"] == scale(base["A"], -1)
    assert reversed_plane["chi"] == base["chi"]

    local = emk_local_two_jet(base["chi"], 4)
    assert local["kappa"] == F(1, 8)
    assert local["W_2"] == F(1, 4)
    assert local["K_0"] == F(-1, 8)

    p_complex = matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ])
    assert mul(p_complex, z) == mul(z, p_complex)
    visible = plane_invariant(mul(p_complex, y), z)
    assert visible["chi"] == 1 != base["chi"]

    c = matrix([[1, 0, 0, 0], [0, 0, 1, 0]])
    cz = mul(c, z)
    assert rank(c) == 2 and rank(c + cz) == 4
    paired_extra = rank(c + cz) - rank(c)
    assert paired_extra == 2

    local_no_sheet = (base["chi"], local["K_0"], 0)
    local_with_sheet = (base["chi"], local["K_0"], 1)
    assert local_no_sheet[:2] == local_with_sheet[:2]
    assert local_no_sheet[2] != local_with_sheet[2]

    z_bad = identity(4)
    assert mul(z_bad, z_bad) != scale(i4, -1)

    return {
        "protocol": "QUANTUM_EMK_LOCAL_GEOMETRY_V1",
        "native_iota_representation_fixture": {
            "Z_squared_minus_identity": True,
            "Z_skew_adjoint": True,
            "wrong_identity_representation_rejected": True,
        },
        "response_plane": {
            "G": strings(base["G"]),
            "A": strings(base["A"]),
            "det_G": str(base["det_G"]),
            "det_A": str(base["det_A"]),
            "chi": str(base["chi"]),
            "gram_margin": str(base["gram_margin"]),
            "GL2_label_invariance": True,
            "orientation_reversal_flips_A_not_chi": True,
        },
        "emk_local_two_jet": {
            "supplied_length_squared": str(local["length_squared"]),
            "kappa": str(local["kappa"]),
            "W_0": str(local["W_0"]),
            "W_1": str(local["W_1"]),
            "W_2": str(local["W_2"]),
            "K_0": str(local["K_0"]),
            "global_warp_selected": False,
        },
        "cut_controls": {
            "Z_invariant_lossy_cut_chi": str(visible["chi"]),
            "uncut_chi": str(base["chi"]),
            "Z_invariance_alone_preserves_chi": False,
            "paired_repair_extra_scalar_channels_on_fixture": paired_extra,
        },
        "helical_guard": {
            "same_local_two_jet_different_sheet_memory": True,
            "sheet_memory_is_local_metric_curvature": False,
            "physical_double_helix_identified": False,
        },
        "scope": {
            "represented_iota_to_response_complex_structure": True,
            "response_plane_invariant_to_emk_local_two_jet_adapter": True,
            "physical_length_scale_derived": False,
            "physical_gravity_identified": False,
            "spacetime_metric_derived": False,
            "quantum_gravity_established": False,
        },
    }


if __name__ == "__main__":
    print("PASS_QUANTUM_EMK_LOCAL_GEOMETRY", run())
