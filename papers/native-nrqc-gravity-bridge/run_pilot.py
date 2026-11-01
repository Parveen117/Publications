"""I/O orchestration only: source arithmetic and postprocessing stay in native JS.

The user-authorized calculator execution calls the existing, pinned NRQC engine.
No source engine is modified, no hardware/API is contacted, and no CI is dispatched.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
PIN='fdd39048a6ddaa2dd78aa3197c274d7766f31a0e'

def run(cmd,**kwargs):
    return subprocess.run(cmd,check=True,text=True,capture_output=True,**kwargs)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--nrqc',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    nrqc=args.nrqc.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    assert run(['git','rev-parse','HEAD'],cwd=nrqc).stdout.strip()==PIN
    assert not run(['git','status','--porcelain','--untracked-files=no'],cwd=nrqc).stdout.strip()
    source=run(['node',str(HERE/'source_interval.cjs')]).stdout
    (out/'source_intervals.json').write_text(source)
    generators=[{'paulis':w} for w in [
        {'0':'X'},{'1':'Z','0':'X'},{'0':'Y'},
        {'1':'Z','0':'Y'},{'1':'X'},{'1':'X','0':'Z'}]]
    states={'vacuum':[], 'flipped':[['x',[0]]],
            'bell':[['h',[0]],['cx',[0,1]]],
            'complex_product':[['h',[0]],['h',[1]],['s',[0]]]}
    execution=[]
    for state,operations in states.items():
        for backend in ['stabilizer','observer']:
            request={'schema':'nrqc-calculator-v1','n_qubits':2,
                     'backend':backend,'operations':operations,
                     'queries':[{'id':'pauli_basis','kind':'geometry','generators':generators}]}
            stem=state+'_'+backend
            req=out/(stem+'.request.json');res=out/(stem+'.result.json')
            req.write_text(json.dumps(request,indent=2)+'\n')
            before=time.perf_counter()
            result=run([sys.executable,'-m','native_relational_qc','calculate',str(req),
                        '--timeout','20','--out',str(res)],cwd=nrqc)
            execution.append({'state':state,'backend':backend,
                              'command_wall_seconds':time.perf_counter()-before})
    # Small local timing sample. Both cold process commands include startup and I/O.
    timings=[]
    for repeat in range(5):
        for engine in (['native','nrqc'] if repeat%2==0 else ['nrqc','native']):
            if engine=='native':
                cmd=['node',str(HERE/'verify_pilot.cjs'),'--baseline',str(out/'source_intervals.json')]
                cwd=HERE
            else:
                cmd=[sys.executable,'-m','native_relational_qc','calculate',
                     str(out/'vacuum_stabilizer.request.json'),'--timeout','20']
                cwd=nrqc
            before=time.perf_counter();run(cmd,cwd=cwd)
            timings.append({'repeat':repeat,'engine':engine,'wall_seconds':time.perf_counter()-before})
    source_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                   for name in ['source_interval.cjs','verify_pilot.cjs','run_pilot.py']}
    native_path=HERE.parent/'native-common-cone-limit'/'native_engine'/'native_operator.cjs'
    source_hashes['../native-common-cone-limit/native_engine/native_operator.cjs']=hashlib.sha256(native_path.read_bytes()).hexdigest()
    meta={'nrqc_commit':PIN,'python':platform.python_version(),
          'node':run(['node','--version']).stdout.strip(),
          'source_sha256':source_hashes,
          'execution':execution,'timings':timings,
          'timing_boundary':'cold local subprocess wall time; shared continuum preparation excluded from both; NRQC weighted assembly measured separately; no scaling inference',
          'memory_comparison':'not measured; no memory advantage claim',
          'hardware_called':False}
    (out/'execution.json').write_text(json.dumps(meta,indent=2)+'\n')
    verified=run(['node',str(HERE/'verify_pilot.cjs'),str(out)])
    (out/'verification.json').write_text(verified.stdout)
    files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.glob('*.json')) if p.name!='manifest.json'}
    (out/'manifest.json').write_text(json.dumps(files,indent=2)+'\n')
    print(verified.stdout)

if __name__=='__main__':main()
