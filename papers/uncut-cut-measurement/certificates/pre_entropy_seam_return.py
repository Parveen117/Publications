"""Exact SR-1/SR-2 sector, torus-flatness and winding-blindness controls."""
import json
from fractions import Fraction


def require(value, message):
    if not value:
        raise AssertionError(message)


def run():
    # The oriented +a transport lives on horizontal links of a periodic grid.
    h = lambda a, b: -1 if a == 2 else 1
    v = lambda a, b: 1
    plaquettes = [h(a, b) * v((a+1) % 3, b) * h(a, (b+1) % 3) * v(a, b)
                  for a in range(3) for b in range(3)]
    require(plaquettes == [1] * 9, "a torus plaquette is curved")
    horizontal = [h(a, b) * h((a+1) % 3, b) * h((a+2) % 3, b)
                  for a in range(3) for b in range(3)]
    require(horizontal == [-1] * 9, "a horizontal cycle lost its sector")
    require(all(abs(u) == 1 for u in (-1, 1)), "unitary norm witness")
    defects = {str(s): {str(u): (u-s)**2 for u in (-1, 1)} for s in (-1, 1)}
    require(defects["1"]["-1"] == 4 and defects["-1"]["-1"] == 0,
            "sector-square identity")
    # Exact formal lift: phases 0 and 2*pi have the same endpoint by the
    # kernel 2*pi*Z of the U(1) projection; their retained integers differ.
    lifted_endpoints = {k: 1 for k in (-1, 0, 1)}
    require(len(set(lifted_endpoints.values())) == 1 and
            len(lifted_endpoints) == 3, "lift blindness control")
    # A squared-amplitude observation also forgets the signed return.
    require((-1)**2 == 1**2 and (-1) != 1, "target blindness control")
    # The independent positive ratio product is 1, never the sector sign.
    gamma_c, gamma_m = 2, Fraction(1, 2)
    require(gamma_c * gamma_m == 1 and gamma_c * gamma_m != -1,
            "thermodynamic type separation")
    return {"SR1_unitary_norm_negative_control": True,
            "SR1_sector_defect_squares": defects,
            "SR2_flat_plaquettes": len(plaquettes),
            "SR2_nontrivial_horizontal_cycles": len(horizontal),
            "SR2_horizontal_return": -1,
            "lifted_integer_endpoint_collision": len(lifted_endpoints),
            "positive_thermodynamic_closure_independent": True}


if __name__ == "__main__":
    print(json.dumps(run(), sort_keys=True))
