#!/usr/bin/env python3
"""Authenticate this payload and reproduce both exact-state verifier routes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

def require(ok,message):
    if not ok: raise RuntimeError(message)

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def authenticate(root,expected_manifest_sha256):
    manifest=root/'MANIFEST.json'
    require(digest(manifest)==expected_manifest_sha256,'Manifest differs from externally pinned SHA256')
    data=json.loads(manifest.read_text())
    expected={entry['path']:entry for entry in data['files']}
    actual={str(p.relative_to(root)):p for p in root.rglob('*') if p.is_file() and p!=manifest}
    require(set(actual)==set(expected),'Payload inventory differs from manifest')
    for name,p in actual.items():
        require(not p.is_symlink(),'Payload symlink forbidden: '+name)
        entry=expected[name]
        require(p.stat().st_size==entry['bytes'],'Payload size differs: '+name)
        require(digest(p)==entry['sha256'],'Payload SHA256 differs: '+name)
    return {'manifest_sha256':digest(manifest),
            'payload_file_count':len(actual),'payload_bytes':sum(p.stat().st_size for p in actual.values())}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-directory',type=Path,required=True)
    parser.add_argument('--expected-manifest-sha256',required=True,
                        help='Manifest SHA256 from separately authenticated public freeze/review receipt')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    output=args.output_directory.resolve()
    require(not output.exists(),'Fresh replay directory required')
    require(root!=output and root not in output.parents,'Replay outputs must be outside authenticated payload')
    expected_manifest=args.expected_manifest_sha256.lower()
    require(len(expected_manifest)==64 and all(ch in '0123456789abcdef' for ch in expected_manifest),
            'Expected manifest SHA256 must have 64 hexadecimal characters')
    before=authenticate(root,expected_manifest)
    require(sys.version.split()[0]=='3.12.14','Use pinned Python 3.12.14')
    import sympy
    require(sympy.__version__=='1.14.0','Use pinned SymPy 1.14.0')
    output.mkdir(parents=True)
    routes=[('primary','verify_incoming_state.py'),('independent','verify_exact_state.py')]
    env=dict(os.environ)
    env['PYTHONDONTWRITEBYTECODE']='1'
    replay_rows=[]
    for route,script in routes:
        mode_science=[]
        for label,optimized in [('NORMAL',False),('OPTIMIZED',True)]:
            dest=output/f'{route}_{label}.json'
            cmd=[sys.executable,'-B']+(['-O'] if optimized else [])
            cmd+=[str(root/route/script),'--output',str(dest)]
            run=subprocess.run(cmd,check=True,capture_output=True,text=True,timeout=120,env=env)
            (output/f'{route}_{label}.stdout.log').write_text(run.stdout)
            (output/f'{route}_{label}.stderr.log').write_text(run.stderr)
            actual=json.loads(dest.read_text())
            expected=json.loads((root/route/f'EXACT_STATE_{label}.json').read_text())
            require(actual==expected,f'{route} {label} receipt differs from archived receipt')
            require(dest.read_bytes()==(root/route/f'EXACT_STATE_{label}.json').read_bytes(),
                    f'{route} {label} receipt bytes differ from archived receipt')
            science=dict(actual)
            science.pop('python_optimization')
            mode_science.append(science)
            replay_rows.append({'route':route,'mode':label,'receipt_sha256':digest(dest),
                                'archived_receipt_exact_match':True,
                                'check_count':actual['check_count'],
                                'mutation_count':actual['mutation_count']})
        require(mode_science[0]==mode_science[1],route+' scientific fields depend on optimization')
    primary=json.loads((output/'primary_NORMAL.json').read_text())
    independent=json.loads((output/'independent_NORMAL.json').read_text())
    names={'canonical_density':'R_c','normalized_phase_density':'R_A',
           'canonical_pressure':'P_c','normalized_phase_pressure':'P_A',
           'canonical_continuity_work':'work_c','normalized_phase_continuity_work':'work_A',
           'canonical_time_integrated_pressure':'int_pressure_c',
           'normalized_phase_time_integrated_pressure':'int_pressure_A'}
    from fractions import Fraction
    for row in primary['sensitivity_rows']:
        secondary=independent['scientific']['rational_sensitivity'][str(row['K'])]
        for p_name,i_name in names.items():
            require(Fraction(row['uniform_norm_coefficients_upper'][p_name])==Fraction(secondary[i_name]),
                    f'Independent sensitivity disagreement: K={row["K"]}, {p_name}')
    after=authenticate(root,expected_manifest)
    require(before==after,'Payload changed during replay')
    receipt={'status':'PASS_FRESH_EXACT_INCOMING_STATE_REPLAY',
             **after,'routes':replay_rows,'independent_rational_coefficients_match':True,
             'independent_rational_coefficient_comparison_count':24,
             'scientific_comparison_excluded_fields':['python_optimization'],
             'archived_receipt_comparison':'all fields and exact bytes, pinned runtime',
             'physical_source_evaluations':0,'saved_array_decodes':0,
             'actual_incoming_state_error':'NOT_ENCLOSED',
             'full_twelve_case_pressure_contact_certificate':'UNRESOLVED',
             'metric_calibration':'FAIL','higher_dimensional_Big_Bang':'NOT_ESTABLISHED'}
    (output/'REPLAY_RECEIPT.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__=='__main__':main()
