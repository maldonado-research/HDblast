#!/usr/bin/env python3
"""Isolated normal/optimized replay and four bounded negative controls."""
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

def need(ok, marker):
    if not ok: raise ValueError(marker)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    need(not args.output.exists(), 'FRESH_OUTPUT_REQUIRED')
    folder = Path(__file__).resolve().parent
    producer = folder/'exact_angular_structure.py'
    controls = {}
    expected = {'cross_block':'PRIMITIVE_BLOCK_INVARIANCE',
                'tt_boundary':'TT_NEUMANN_LINEAR_PART',
                'tt_bulk_parity':'TT_BULK_SOURCE_PARITY_metric_forcing',
                'reflection':'ANGULAR_REFLECTION_EXCHANGES_TT'}
    with tempfile.TemporaryDirectory(prefix='angular-structure-') as tmp:
        outputs = []
        for optimized in (False, True):
            output = Path(tmp)/('optimized.json' if optimized else 'normal.json')
            command = [sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])
            command += [str(producer),'--repo-root',str(args.repo_root),'--output',str(output)]
            run = subprocess.run(command, capture_output=True, text=True)
            need(run.returncode == 0, 'REPLAY_FAILED_'+run.stderr)
            outputs.append(output.read_bytes())
        need(outputs[0] == outputs[1], 'NORMAL_OPTIMIZED_DIFFER')
        for name in ('ANGULAR_STRUCTURE_NORMAL.json','ANGULAR_STRUCTURE_OPTIMIZED.json'):
            need((folder/name).read_bytes() == outputs[0], 'RETAINED_OUTPUT_DIFFERS_'+name)
        for name, marker in expected.items():
            output = Path(tmp)/(name+'.json')
            command = [sys.executable,'-I','-S','-B','-O',str(producer),
                       '--repo-root',str(args.repo_root),'--output',str(output),'--mutation',name]
            run = subprocess.run(command,capture_output=True,text=True)
            need(run.returncode != 0 and marker in run.stderr and not output.exists(), 'CONTROL_FAILED_'+name)
            controls[name] = {'exit_code':run.returncode,'expected_marker':marker,
                              'observed_marker':run.stderr.strip().splitlines()[-1],
                              'no_artifact_written':True}
        result = {'status':'PASS_ALL_CONTROLS','normal_optimized_byte_identical':True,
                  'isolated_without_site_packages':True,
                  'retained_outputs_equal_to_replay':True,
                  'output_sha256':hashlib.sha256(outputs[0]).hexdigest(),
                  'producer_sha256':hashlib.sha256(producer.read_bytes()).hexdigest(),
                  'controls':controls}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'output':str(args.output)}))

if __name__ == '__main__': main()
