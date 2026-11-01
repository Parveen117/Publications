#!/usr/bin/env python3
"""Build the checked-in native edition; verification never runs this writer."""
import os
from pathlib import Path
import subprocess

home = Path(__file__).resolve().parent
source = home / 'source'
pdf = home / 'pdf'
scratch = home.parents[1] / 'tmp' / 'pdfs' / 'native-thermo'
for directory in (source, pdf, scratch):
    directory.mkdir(parents=True, exist_ok=True)
markdown = (home / 'THEOREM.md').read_text()
markdown = markdown.replace('\\[\n', '\n\n$$\n').replace('\\]\n', '$$\n\n')
subprocess.run([
    'pandoc',
    '--from=markdown+tex_math_single_backslash+tex_math_dollars', '--standalone',
    '--shift-heading-level-by=-1', '--variable=fontsize:11pt',
    '--variable=geometry:margin=0.88in', '--variable=colorlinks:true',
    '--variable=linkcolor:blue', '--variable=urlcolor:blue',
    '--include-in-header=' + str(source / 'preamble.tex'),
    '--output=' + str(source / 'main.tex'),
], input=markdown, text=True, check=True)
env = dict(os.environ, SOURCE_DATE_EPOCH='1790985600', FORCE_SOURCE_DATE='1')
for _ in range(2):
    with (scratch / 'build-output.txt').open('w') as output:
        result = subprocess.run([
            'pdflatex', '-interaction=nonstopmode', '-halt-on-error',
            '-output-directory=' + str(scratch), str(source / 'main.tex'),
        ], env=env, stdout=output, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError((scratch / 'build-output.txt').read_text()[-2500:])
log = (scratch / 'main.log').read_text()
bad = [line for line in log.splitlines() if 'Overfull' in line or 'undefined references' in line]
if bad:
    raise RuntimeError('Review manuscript layout/reference warnings: ' + '\n'.join(bad))
(pdf / 'native-thermodynamic-curvature.pdf').write_bytes((scratch / 'main.pdf').read_bytes())
print(pdf / 'native-thermodynamic-curvature.pdf')
