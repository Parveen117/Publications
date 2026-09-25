"""Exact Spin(2)/SO(2) quarter-turn and flat spin-structure controls."""
from fractions import Fraction
import json


def cmul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def cpow(a, n):
    result = (1, 0)
    for _ in range(n):
        result = cmul(result, a)
    return result


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def run():
    i = (0, 1)
    # At deficit pi: vector angle pi and spin half-angle pi/2.
    vector_pi = cpow(i, 2)
    spin_pi = i
    require(vector_pi == (-1, 0) and cmul(spin_pi, spin_pi) == vector_pi,
            "cone half-angle covering identity")
    # Spin -1 corresponds to a full vector turn (and to the boundary
    # delta=2*pi, excluded from the nondegenerate alpha>0 cone family).
    spin_two_pi = cpow(i, 2)
    vector_two_pi = cpow(i, 4)
    require(spin_two_pi == (-1, 0) and vector_two_pi == (1, 0),
            "double-cover kernel")
    # One flat spin deck sign: every plaquette returns +1, but the
    # horizontal fundamental cycle returns -1 on spin and +1 on vectors.
    edge = lambda a: -1 if a == 2 else 1
    squares = [edge(a)*edge(a) for a in range(3) for _ in range(3)]
    cycles = [edge(0)*edge(1)*edge(2) for _ in range(3)]
    require(squares == [1]*9 and cycles == [-1]*3,
            "flat spin monodromy")
    require(all(x*x == 1 for x in cycles), "vector return is identity")
    # For f=1+u², f'(v)-f'(u)=2(v-u); around a unit-width
    # rectangle the integrated connection equals 2 times area.
    u, v, width = Fraction(1, 3), Fraction(2, 3), Fraction(3, 5)
    area = (v-u)*width
    integral = (2*v-2*u)*width
    require(integral == 2*area and integral / area == 2,
            "local curvature coefficient")
    return {"cone_deficit_pi_vector_minus_identity": True,
            "cone_deficit_pi_spin_phase": "i (up to orientation and spin-structure sign)",
            "spin_minus_one_vector_identity": True,
            "flat_spin_plaquettes": len(squares),
            "flat_spin_negative_cycles": len(cycles),
            "smooth_adapter_rectangle_area": "1/5",
            "smooth_adapter_curvature_coefficient": 2,
            "spin_minus_one_selects_cone_deficit": False,
            "dimensionful_hbar_selected_from_sign": False}


if __name__ == "__main__":
    print(json.dumps(run(), sort_keys=True))
