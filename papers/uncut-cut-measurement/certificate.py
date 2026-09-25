"""Exact finite evidence. Default/check mode never rewrites approved evidence."""
import argparse
import hashlib
import json
from itertools import combinations, permutations, product
from pathlib import Path

import finite_model as fm
import family_certificate
import partial_certificate
import gravity_certificate
import selection_certificate
import grading_certificate
import interaction_certificate
import response_certificate
import identification_certificate
import source_response_certificate

HERE = Path(__file__).resolve().parent
BOUND_SOURCES = (
    "README.md", "MANUSCRIPT.md", "CLAIMS.md", "RESEARCH_PROGRAMME.md",
    "SOURCE_PINS.json", "finite_model.py", "certificate.py",
    "tests/test_finite_model.py",
    "FAMILY_CONTINUATION.md", "family_model.py", "family_certificate.py",
    "tests/test_family_model.py",
    "NOVELTY_AND_LINEAGE.md", "ADMISSIBLE_CONTINUATION.md",
    "partial_model.py", "partial_certificate.py", "tests/test_partial_model.py",
    "GRAVITY_BEFORE_CURVATURE.md", "GRAVITY_SOURCE_PINS.json",
    "gravity_model.py", "gravity_certificate.py", "tests/test_gravity_model.py",
    "SECTOR_AND_COFRAME_SELECTION.md", "SELECTION_SOURCE_PINS.json",
    "selection_model.py", "selection_certificate.py", "tests/test_selection_model.py",
    "NATIVE_GRADING_AND_SEAM_MEMORY.md", "GRADING_SOURCE_PINS.json",
    "grading_model.py", "grading_certificate.py", "tests/test_grading_model.py",
    "NATIVE_INTERACTION_AND_LAWFUL_CUTS.md", "INTERACTION_SOURCE_PINS.json",
    "interaction_model.py", "interaction_certificate.py", "tests/test_interaction_model.py",
    "NATIVE_RESPONSE_TENSOR_AND_CUT_LEDGER.md", "RESPONSE_SOURCE_AUDIT.md",
    "RESPONSE_SOURCE_PINS.json", "response_model.py", "response_certificate.py",
    "tests/test_response_model.py",
    "../curvature-information-duality/README.md",
    "../curvature-information-duality/LINEAGE.md",
    "../curvature-information-duality/certificates/qth1_quantum_recognition_information.py",
    "../curvature-information-duality/certificates/QTH1_RESULT.json",
    "../curvature-information-duality/certificates/EXPECTED_QTH1.sha256",
    "../curvature-information-duality/tests/test_qth1.py",
    "PROBE_CALIBRATION_AND_NATIVE_IDENTIFICATION.md", "IDENTIFICATION_SOURCE_PINS.json",
    "identification_model.py", "identification_certificate.py", "tests/test_identification_model.py",
    "SOURCE_RESPONSE_AND_GRAVITY_TESTS.md", "GRAVITY_EXPERIMENT_CONTRACT.md",
    "SOURCE_RESPONSE_PINS.json", "source_response_model.py",
    "source_response_certificate.py", "tests/test_source_response_model.py",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def equal_implies(fine, coarse):
    """Independent pairwise equality oracle for the enumerated carriers."""
    return all(fine[x] != fine[y] or coarse[x] == coarse[y]
               for x in range(len(fine)) for y in range(len(fine)))


def brute_decoder_exists(cut, target):
    cut_labels, target_labels = tuple(dict.fromkeys(cut)), tuple(dict.fromkeys(target))
    for outputs in product(target_labels, repeat=len(cut_labels)):
        decoder = dict(zip(cut_labels, outputs))
        if all(decoder[a] == v for a, v in zip(cut, target)):
            return True
    return False


def finite_checks():
    counts = {
        "cut_target_pairs": 0,
        "smaller_memory_maps_rejected": 0,
        "refinement_target_triples": 0,
        "common_target_triples": 0,
        "cut_transition_pairs": 0,
        "closed_refinement_comparisons": 0,
        "natural_transition_comparisons": 0,
        "represented_block_pairs": 0,
        "symmetry_permutations": 0,
    }
    domain = []
    for n in range(1, 5):
        cuts = tuple(fm.partitions(n))
        domain.append({"carrier_size": n, "cut_partitions": len(cuts),
                       "target_partitions": len(cuts), "endomaps": n ** n})
        for cut, target in product(cuts, repeat=2):
            actual = fm.factor_map(cut, target)
            require((actual is not None) == brute_decoder_exists(cut, target),
                    f"U1 decoder mismatch: {cut}, {target}")
            memory = fm.repair_channel(cut, target)
            k = fm.repair_size(cut, target)
            require(len(set(memory)) == k, "U2 constructed alphabet has wrong size")
            require(equal_implies(tuple(zip(cut, memory)), target),
                    "U2 constructed channel fails to recover target")
            # If a smaller nonempty alphabet worked it could be embedded in k-1.
            # For k=1, a map from nonempty X into an empty alphabet is impossible.
            if k > 1:
                for smaller in product(range(k - 1), repeat=n):
                    require(not equal_implies(tuple(zip(cut, smaller)), target),
                            "U2 false minimum: a smaller channel works")
                    counts["smaller_memory_maps_rejected"] += 1
            bits = fm.fixed_binary_bits(k)
            require(2 ** bits >= k and (bits == 0 or 2 ** (bits - 1) < k),
                    "U2 binary length is not minimal")
            counts["cut_target_pairs"] += 1

        for coarse, fine in product(cuts, repeat=2):
            if equal_implies(fine, coarse):
                for target in cuts:
                    coarse_sets = fm.target_sets(coarse, target)
                    fine_sets = fm.target_sets(fine, target)
                    for x in range(n):
                        require(fine_sets[fine[x]] <= coarse_sets[coarse[x]],
                                "U3 set enclosure failed")
                    kgc = fm.repair_size(coarse, target)
                    kgd = fm.repair_size(fine, target)
                    kdc = fm.repair_size(coarse, fine)
                    require(kgd <= kgc <= kdc * kgd, "U3 repair bound failed")
                    counts["refinement_target_triples"] += 1
            components = fm.common_components((coarse, fine))
            for target in cuts:
                individual = (equal_implies(coarse, target) and
                              equal_implies(fine, target))
                require(individual == equal_implies(components, target),
                        "U4 common-target characterization failed")
                counts["common_target_triples"] += 1

        for transition in fm.all_endomaps(n):
            closed = {cut: equal_implies(cut, tuple(cut[y] for y in transition))
                      for cut in cuts}
            for cut in cuts:
                observed = fm.descended_transition(cut, transition)
                require((observed is not None) == closed[cut], "U5 descent failed")
                if observed is not None:
                    require(all(observed[cut[x]] == cut[transition[x]]
                                for x in range(n)), "U5 decoder diagram failed")
                result = fm.stable_history(cut, transition)
                stable = result["cut"]
                require(result["depth"] <= n - len(set(cut)), "U6 depth bound failed")
                require(equal_implies(stable, cut), "U6 final cut is not a refinement")
                require(closed[stable], "U6 final cut is not closed")
                # Direct length-(n+1) signatures, independently of the stopping rule.
                signatures = []
                for x in range(n):
                    values, current = [], x
                    for _ in range(n + 1):
                        values.append(cut[current])
                        current = transition[current]
                    signatures.append(tuple(values))
                require(equal_implies(stable, signatures) and
                        equal_implies(signatures, stable), "U6 history oracle mismatch")
                for fine in cuts:
                    if not (equal_implies(fine, cut) and closed[fine]):
                        continue
                    require(equal_implies(fine, stable), "U6 refinement is not coarsest")
                    counts["closed_refinement_comparisons"] += 1
                    if closed[cut]:
                        r = fm.factor_map(fine, cut)
                        tf = fm.descended_transition(fine, transition)
                        require(all(r[tf[b]] == observed[r[b]] for b in set(fine)),
                                "U5 lawful refinement does not intertwine")
                        counts["natural_transition_comparisons"] += 1
                counts["cut_transition_pairs"] += 1

    # Direct full matrix multiplication compared with the retained/memory corners.
    matrices = tuple(product((-1, 0, 1), repeat=4))
    for u, v in product(matrices, repeat=2):
        U = (u[:2], u[2:])
        V = (v[:2], v[2:])
        VU = tuple(tuple(sum(V[i][k] * U[k][j] for k in range(2))
                         for j in range(2)) for i in range(2))
        require(VU[0][0] - V[0][0] * U[0][0] == V[0][1] * U[1][0],
                "U7 composition residue failed")
        counts["represented_block_pairs"] += 1

    transition = (0, 0, 1)
    symmetries = []
    for s in permutations(range(3)):
        if all(s[transition[x]] == transition[s[x]] for x in range(3)):
            symmetries.append(s)
        counts["symmetry_permutations"] += 1
    require(symmetries == [(0, 1, 2)], "S0 negative control failed")
    return domain, counts


def controls():
    c, g = (0, 0, 1, 1), (0, 1, 0, 1)
    require(fm.factor_map(c, g) is None, "W1 hidden bit became recoverable")
    require(fm.repair_size(c, g) == 2 and fm.repair_size(c, c) == 1,
            "W1 target-relative cost failed")
    coarse = (0, 0, 0, 0)
    require(fm.repair_size(coarse, g) == 2 and
            fm.repair_size(coarse, c) * fm.repair_size(c, g) == 4,
            "U3 strict inequality control failed")
    states = ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0))
    cuts = tuple(tuple(x[j] for x in states) for j in range(3))
    pair_witnesses = {f"{i},{j}": fm.joint_fibre((cuts[i], cuts[j]), (1, 1))
                      for i, j in combinations(range(3), 2)}
    require(all(pair_witnesses.values()), "W2 pairwise consistency failed")
    require(not fm.joint_fibre(cuts, (1, 1, 1)), "W2 false global state")
    require(fm.common_components((c, g)) == (0, 0, 0, 0) and
            len(set(zip(c, g))) == 4, "W3 common/joint distinction failed")
    history = fm.stable_history((0, 0, 0, 1), (1, 2, 3, 3))
    require(history["depth"] == 2 and len(set(history["cut"])) == 4,
            "W4 delayed split failed")
    require(fm.repair_size((0, 0, 0, 1), history["cut"]) == 3,
            "W4 memory minimum failed")
    return {
        "W1_hidden_bit": {"readout_only_recovery": False,
                           "hidden_target_labels": 2, "visible_target_labels": 1},
        "U3_strict_product": {"actual": 2, "product_bound": 4},
        "W2_pairwise_not_global": {"states": states, "pair_witness_indices": pair_witnesses,
                                   "global_fibre": []},
        "W3_common_not_joint": {"common_components": 1, "joint_distinct_readouts": 4},
        "W4_delayed_history": {"class_counts": [len(set(x)) for x in history["levels"]],
                               "stable_depth": 2, "memory_labels": 3, "binary_bits": 2},
        "U7_memory_return": {"U": [[0, 1], [1, 0]], "V": [[0, 1], [1, 0]],
                             "visible_step_product": 0, "visible_composite": 1,
                             "memory_return": 1},
        "S0_transitive_not_nontrivial_symmetry": {"transition": [0, 0, 1],
                                                  "commuting_permutations": [[0, 1, 2]]},
    }


def build():
    domain, counts = finite_checks()
    body = {
        "protocol": "UNCUT_CUT_MEASUREMENT_V1_0",
        "status": "PASS_FINITE_CHECKS",
        "arithmetic": "exact integer and rational arithmetic with declared finite coverage; no floating point",
        "domains": domain,
        "represented_blocks": {"shape": [2, 2], "entry_alphabet": [-1, 0, 1]},
        "checks": counts,
        "controls": controls(),
        "family_continuation": family_certificate.build_checks(),
        "admissible_continuation": partial_certificate.build_checks(),
        "conditional_gravity_bridge": gravity_certificate.build_checks(),
        "sector_and_coframe_constraints": selection_certificate.build_checks(),
        "native_grading_and_seam_memory": grading_certificate.build_checks(),
        "native_interaction_and_lawful_cuts": interaction_certificate.build_checks(),
        "native_response_tensor_and_cut_ledger": response_certificate.build_checks(),
        "calibrated_probe_identification": identification_certificate.build_checks(),
        "conditional_source_response_and_gravity_tests": source_response_certificate.build_checks(),
        "source_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest()
                          for p in BOUND_SOURCES},
        "scope": {
            "general_written_proofs": "MANUSCRIPT.md U1-U7, FAMILY_CONTINUATION.md U8-U12, ADMISSIBLE_CONTINUATION.md U13-U16, GRAVITY_BEFORE_CURVATURE.md U17-U20, SECTOR_AND_COFRAME_SELECTION.md U21-U23, NATIVE_GRADING_AND_SEAM_MEMORY.md U24-U26, NATIVE_INTERACTION_AND_LAWFUL_CUTS.md U27-U29, NATIVE_RESPONSE_TENSOR_AND_CUT_LEDGER.md U30-U32, PROBE_CALIBRATION_AND_NATIVE_IDENTIFICATION.md U33-U35 and SOURCE_RESPONSE_AND_GRAVITY_TESTS.md U36-U38 under stated hypotheses; no priority claim for standard methods",
            "availability_interface_physically_established": False,
            "native_gravity_law_physically_identified": False,
            "spacetime_dimension_derived": False,
            "transport_sector_selected_by_native_law": False,
            "coframe_profile_constrained_within_supplied_homogeneous_class": True,
            "homogeneity_derived_from_native_law": False,
            "native_grading_repair_on_declared_edge_assembly": True,
            "native_compatible_inter_edge_candidate_constructed": True,
            "physical_inter_edge_law_selected_from_primitive": False,
            "genuine_lossy_quotient_preserves_declared_interaction": True,
            "model_response_transcript_physically_realized": False,
            "native_response_tensor_and_exact_cut_seam_constructed": True,
            "native_tensor_identified_with_quantum_information_metric": False,
            "antisymmetric_response_identified_with_connection_curvature": False,
            "declared_native_pair_identified_with_augmented_readout": True,
            "bounded_error_probe_rejection_constructed": True,
            "physical_probe_universality_established": False,
            "native_coupling_values_selected_from_primitive": False,
            "explicit_affine_source_cost_candidate_constructed": True,
            "conditional_inverse_square_shell_sector_constructed": True,
            "quadratic_cost_and_shell_growth_selected_from_primitive": False,
            "native_conservative_flow_proves_relaxation": False,
            "physical_inverse_square_derived_without_extra_hypotheses": False,
            "physical_gravity_experiment_passed": False,
            "post_cut_smooth_chart_is_supplied": True,
            "formal_proof_assistant_verification": False,
            "uncut_identified_with_finite_carrier": False,
            "physical_quantum_classical_derivation": False,
            "physical_spacetime_derivation": False,
            "minimum_nature_action_derived": False,
            "empirical_validation": False,
            "independent_review": False,
        },
    }
    return (json.dumps(body, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Read-only verification (default)")
    mode.add_argument("--write", action="store_true", help="Explicitly update evidence and pin")
    args = parser.parse_args()
    data = build()
    digest = hashlib.sha256(data).hexdigest()
    certificate, pin = HERE / "CERTIFICATE.json", HERE / "EXPECTED.sha256"
    if args.write:
        certificate.write_bytes(data)
        pin.write_text(digest + "\n", encoding="ascii")
        print("WROTE_FINITE_EVIDENCE", digest)
    else:
        require(certificate.read_bytes() == data, "Certificate or source drift: review the change")
        require(pin.read_text(encoding="ascii").strip() == digest, "Certificate pin mismatch")
        print("PASS_FINITE_CHECKS", digest)


if __name__ == "__main__":
    main()
