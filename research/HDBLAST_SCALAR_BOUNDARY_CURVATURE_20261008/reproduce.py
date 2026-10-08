"""Reproduce the scalar-boundary checkpoint into a fresh local directory."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--numerical',action='store_true')
    args=ap.parse_args()
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((ROOT/'MANIFEST.json').read_text())
    for row in manifest['files']:
        data=(ROOT/row['path']).read_bytes()
        if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:
            raise ValueError('Delivery integrity failure: '+row['path'])
    calls=[]
    def run(name,script,*extra,optimized=False):
        argv=[sys.executable,'-I','-B']+(['-O'] if optimized else [])+[str(script),*map(str,extra)]
        with (out/(name+'.log')).open('w') as log:
            result=subprocess.run(argv,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=False)
        calls.append(dict(name=name,exit_code=result.returncode))
        if result.returncode:
            raise RuntimeError('Failed '+name+'; inspect '+str(out/(name+'.log')))
    run('finite_exact',ROOT/'independent-math/verify_independent_finite_detuning.py')
    run('primary_translations',ROOT/'literature/verify_prior_art_translations.py')
    constants=ROOT/'numerical-audit/uniform-constants/evaluate_uniform_constants.py'
    theorem=ROOT/'branch-theorem/UNIFORM_BRANCH_THEOREM.md'
    proof_pin=hashlib.sha256(theorem.read_bytes()).hexdigest()
    for suffix,opt in [('normal',False),('optimized',True)]:
        run('uniform_'+suffix,constants,'--theorem',theorem,'--theorem-sha256',proof_pin,'--out',out/('uniform_'+suffix+'.json'),optimized=opt)
        run('independent_uniform_'+suffix,ROOT/'independent-math/check_uniform_constants_independently.py',optimized=opt)
    if (out/'uniform_normal.json').read_bytes()!=(out/'uniform_optimized.json').read_bytes():
        raise ValueError('Interpreter mode changed exact constants')
    if (out/'independent_uniform_normal.log').read_bytes()!=(out/'independent_uniform_optimized.log').read_bytes():
        raise ValueError('Interpreter mode changed independent exact checks')
    reference=ROOT/'reference/HDBLAST-APS-research-review-20261008/general-scalar-response'
    replay=out/'replay'
    shutil.copytree(ROOT/'numerical-audit/replay',replay)
    run('original_algebra',replay/'verify_response_algebra.py')
    if args.numerical:
        # The preserved historical producer prepends this path. Refuse an
        # unrelated import directory rather than silently using its contents.
        if Path('/private/tmp/aps-research-deps').exists():
            raise RuntimeError('Use a clean environment without /private/tmp/aps-research-deps')
        run('family_36',replay/'family_benchmark.py')
        bootstrap=out/'matched_bootstrap.py'
        bootstrap.write_text("from pathlib import Path\nimport runpy,sys\np=Path(__file__).resolve().parent/'replay'\nsys.path.insert(0,str(p))\nrunpy.run_path(str(p/'matched_pair_benchmark.py'),run_name='__main__')\n")
        run('matched_12',bootstrap)
        independent=out/'finite-detuning-review'
        independent.mkdir()
        shutil.copy2(ROOT/'finite-detuning-review/verify_finite_detuning.py',independent)
        run('independent_72',independent/'verify_finite_detuning.py')
        diagnostics=json.loads((independent/'FINITE_DETUNING_DIAGNOSTICS.json').read_text())
        summary=diagnostics['summary']
        runs=[r for m in diagnostics['models'] for r in m['runs']]
        if not (len(runs)==72 and all(r['root_success'] for r in runs)
                and summary['all_C2_negative']
                and summary['all_stable_centered_errors_decrease']
                and summary['max_junction_residual']<1e-10
                and summary['max_H_squared_settings_difference']<1e-10):
            raise ValueError('Independent numerical replay diagnostics failed; not an ODE error certificate')
    run('saved_arithmetic_audit',ROOT/'numerical-audit/audit_results.py','--reference',reference,'--replay',replay,'--out',out/'REPLAY_AUDIT.json')
    run('conditional_remainder',ROOT/'numerical-audit/conditional_remainder.py','--replay-dir',replay,'--out',out/'CONDITIONAL_REMAINDER.json')
    receipt=dict(status='PASS_REPRODUCTION',numerical_integrations_rerun=120 if args.numerical else 0,manifest_files_verified=len(manifest['files']),calls=calls,scope='Mathematical proofs are reviewed in the pinned text. Arithmetic and floating-point replay do not replace those proofs, establish external novelty, or demonstrate a Big Bang origin.')
    (out/'REPRODUCTION.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    main()
