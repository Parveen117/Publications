"""GitHub-only incremental execution; Python performs provenance and file I/O."""
from pathlib import Path
import ast
import hashlib
import json
import os
import subprocess
import sys
import urllib.request

if sys.version_info[:2] != (3, 12) or os.environ.get('GITHUB_ACTIONS') != 'true':
    raise SystemExit('Use GitHub Actions with Python 3.12')
packet = Path(__file__).resolve().parent
root = packet.parents[1]
baseline = '0c98adf69682f6fee48be85dcbcb7667d77dea2e'
tested = '724c8a6a58a4a6955a53dbc4999f1630f0a64e27'
tested_tree = '1d24f4d40640c951823efcf287487c01f362109b'
repo = 'Parveen117/Publications'
stage = packet.name
def git(*args):
    return subprocess.check_output(['git', *args], cwd=root, text=True).strip()
def api(endpoint):
    request = urllib.request.Request('https://api.github.com/repos/' + repo + endpoint,
        headers={'Authorization': 'Bearer ' + os.environ['GITHUB_TOKEN'],
                 'Accept': 'application/vnd.github+json', 'X-GitHub-Api-Version': '2022-11-28'})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)
subprocess.run(['git', 'fetch', '--no-tags', '--depth=2', 'origin', baseline], cwd=root, check=True)
assert git('rev-parse', baseline + '^{tree}') == '72d69b98bb974ab09a2fd45714d92101cad372ac'
assert git('rev-parse', baseline + '^') == tested
assert git('rev-parse', tested + '^{tree}') == tested_tree
old_packet = root / 'papers/native-layer-analytic-return'
old_outputs = {'papers/native-layer-analytic-return/' + name
               for name in ('CERTIFICATE.json', 'EXPECTED.sha256', 'CI_PROVENANCE.json')}
assert set(git('diff', '--name-only', tested, baseline).splitlines()) == old_outputs
prior = json.loads((old_packet / 'CI_PROVENANCE.json').read_text())
assert prior['source_head_sha'] == tested and prior['source_tree'] == tested_tree
assert prior['run_id'] == '38997057551' and len(prior['reports']) == 59
assert prior['outputs']['CERTIFICATE.json'] == '8a2550206c26f2dc78d9cdc5b742fda5e39e9baf42f6a63465f5ec5273d4bedf'
for name, expected in prior['outputs'].items():
    assert hashlib.sha256((old_packet / name).read_bytes()).hexdigest() == expected
assert hashlib.sha256((old_packet / 'SOURCE_PINS.json').read_bytes()).hexdigest() == prior['source_pins_sha256']
run = api('/actions/runs/38997057551')
assert run['conclusion'] == 'success' and run['head_sha'] == tested
job = api('/actions/jobs/117291825302')
assert job['conclusion'] == 'success'
for name in ('Replay inherited native cut action', 'Replay source-bound quantum geometry chain'):
    assert next(s for s in job['steps'] if s['name'] == name)['conclusion'] == 'success'
allowed = {'README.md', 'MATHEMATICAL_FOUNDATION.md', 'QUANTUM_GRAVITY_RESEARCH.md',
           'RESEARCH_CLOSURES.md', 'scripts/replay_quantum_geometry.py',
           '.github/workflows/quantum-geometry-propagation.yml'}
changed = git('diff', '--name-only', baseline, 'HEAD').splitlines()
assert all(p in allowed or p.startswith('papers/' + stage + '/') for p in changed), 'Predecessor source changed'
tree_ast = ast.parse((root / 'scripts/replay_quantum_geometry.py').read_text())
stages = next(ast.literal_eval(n.value) for n in tree_ast.body if isinstance(n, ast.Assign)
              and any(isinstance(t, ast.Name) and t.id == 'stages' for t in n.targets))
assert stages[-1] == stage and stages[:-1] == [r['stage'] for r in prior['reports']]
outputs = ('SOURCE_PINS.json', 'CERTIFICATE.json', 'EXPECTED.sha256', 'CI_PROVENANCE.json')
present = [(packet / name).exists() for name in outputs]
assert all(present) or not any(present), 'Partial generated-output set'
generating = not any(present)
if generating:
    inputs = json.loads((packet / 'PIN_INPUTS.json').read_text())
    assert inputs['parent'] == baseline
    pins = {'protocol': 'VSR1_SOURCE_PINS_V1', 'repository_parent': baseline,
            'files': {name: hashlib.sha256((packet / name).read_bytes()).hexdigest()
                      for name in inputs['files']}}
    (packet / 'SOURCE_PINS.json').write_text(json.dumps(pins, indent=2) + '\n')
proc = subprocess.run([sys.executable, str(packet / 'run.py'), '--write' if generating else '--check'],
                      cwd=root, text=True, capture_output=True)
if proc.returncode:
    print(proc.stdout)
    print(proc.stderr, file=sys.stderr)
    raise SystemExit(proc.returncode)
report = json.loads(proc.stdout)
if generating:
    provenance = {
        'protocol': 'VSR1_GITHUB_INCREMENTAL_V1', 'repository': repo,
        'run_id': os.environ['GITHUB_RUN_ID'], 'run_attempt': os.environ['GITHUB_RUN_ATTEMPT'],
        'source_head_sha': os.environ['SOURCE_HEAD_SHA'], 'checkout_commit': git('rev-parse', 'HEAD'),
        'source_tree': git('rev-parse', 'HEAD^{tree}'),
        'predecessor_commit': baseline, 'predecessor_tested_commit': tested,
        'predecessor_run_id': 38997057551, 'reused_stages': 59, 'executed_stages': 1,
        'predecessor_sources_unchanged': True,
        'outputs': {name: hashlib.sha256((packet / name).read_bytes()).hexdigest()
                    for name in outputs if name != 'CI_PROVENANCE.json'},
        'report': report,
        'scope': '59 immutable passing stages reused; new stage executed on GitHub; import commit not freshly replayed'}
    (packet / 'CI_PROVENANCE.json').write_text(json.dumps(provenance, indent=2) + '\n')
print(json.dumps({'passed': True, 'reused_stages': 59, 'executed_stages': 1,
                  'total_stages': 60, 'new_stage': stage, **report}))
