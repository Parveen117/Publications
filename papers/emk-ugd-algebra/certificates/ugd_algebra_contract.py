#!/usr/bin/env python3
"""Audit the existing UGD carry/equality contract without rewriting its source."""
import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def state(value):
    return {'K': value.K, 'digits': [[n, a, s] for n, (a, s) in sorted(value.digits.items())],
            'ledger': value.ledger, 'total_seam': value.total_seam()}


def build():
    if not __debug__:
        raise RuntimeError('Do not disable source/certificate assertions')
    assert sys.version_info[:2] in [(3, 11), (3, 12)]
    pins_bytes = (HERE / 'UGDA0_SOURCE_PINS.json').read_bytes()
    pins = json.loads(pins_bytes)
    for row in pins['local_inputs']:
        source = (ROOT / row['path']).resolve()
        assert source.is_relative_to(ROOT), 'Source outside publication root'
        data = source.read_bytes()
        assert sha(data) == row['sha256'] and blob(data) == row['git_blob_sha'], row['path']
    source = HERE / 'ugd1_numerals.py'
    spec = importlib.util.spec_from_file_location('ugd_a0_pinned_source', source)
    ugd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ugd)
    digits = {s: ugd.digit(0, 0, s, K=4) for s in [-1, 0, 1]}
    x, y, z = [digits[s] for s in [-1, 1, 1]]
    left = ugd.add_digits(ugd.add_digits(x, y), z)
    right = ugd.add_digits(x, ugd.add_digits(y, z))
    assert left != right
    assert left.digits == {0: (0, 1)} and left.ledger == 0
    assert right.digits == {0: (0, 0)} and right.ledger == 1
    assert left.total_seam() == right.total_seam() == 1
    failures = []
    for a, b, c in itertools.product([-1, 0, 1], repeat=3):
        l = ugd.add_digits(ugd.add_digits(digits[a], digits[b]), digits[c])
        r = ugd.add_digits(digits[a], ugd.add_digits(digits[b], digits[c]))
        assert l.total_seam() == r.total_seam() == a + b + c
        if l != r:
            failures.append({'inputs': [a, b, c], 'left': state(l), 'right': state(r)})
    samples = [ugd.Numeral({0: (0, s)}, ledger=l, K=4)
               for s, l in itertools.product([-1, 0, 1], range(-2, 3))]
    for a, b in itertools.product(samples, repeat=2):
        out = ugd.add_digits(a, b)
        assert out.total_seam() == a.total_seam() + b.total_seam()
    bad = ugd.add_digits_broken_carry(digits[1], digits[1])
    assert bad.total_seam() != 2
    return {
        'schema': 'publications.ugd-a0.v1',
        'status': 'PASS_UGD_A0_ASSOCIATIVE_CARRIER_GATE_AUDIT',
        'source_pins_sha256': sha(pins_bytes),
        'strict_state_addition_associative': False,
        'charge_target_addition_associative': True,
        'full_UGD_algebra_certified': False,
        'written_results': ['A0.1 strict-state nonassociativity', 'A0.2 total-charge conservation',
                            'A0.3 associative-algebra import requires a new contract'],
        'exact_witness': {'inputs': [-1, 1, 1], 'phase_exponents': [0, 0, 0],
                          'left': state(left), 'right': state(right)},
        'finite_checks': {'balanced_triples': 27, 'strict_associativity_failures': len(failures),
                          'seam_ledger_conservation_pairs': len(samples) ** 2},
        'all_strict_failure_witnesses': failures,
        'negative_controls': {
            'strict_associativity_promotion_rejected': left != right,
            'charge_quotient_as_full_state_equality_rejected': left != right and
                 left.total_seam() == right.total_seam(),
            'ledger_erasure_rejected': bad.total_seam() != 2},
        'boundary': {'historical_source_modified': False, 'finite_no_go_not_all_UGD_presentations': True,
                     'linked_seam_rule_derived': False, 'infinite_numerals_certified': False,
                     'formal_proof_assistant_verified': False, 'physical_validation': False}}


def main():
    if not __debug__:
        raise RuntimeError('Do not disable source/certificate assertions')
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--check', action='store_true')
    modes.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = build()
    data = (json.dumps(result, sort_keys=True, indent=2, ensure_ascii=False) + '\n').encode()
    result_path, expected = HERE / 'UGDA0_RESULT.json', HERE / 'EXPECTED_UGDA0.sha256'
    if args.write:
        result_path.write_bytes(data)
        expected.write_text(sha(data) + '\n')
    else:
        assert result_path.read_bytes() == data and expected.read_text().strip() == sha(data)
    print(json.dumps({'status': result['status'], **result['finite_checks'],
                      'certificate_sha256': sha(data)}))


if __name__ == '__main__':
    main()
