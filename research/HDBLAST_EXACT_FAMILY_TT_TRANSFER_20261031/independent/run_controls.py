"""Run genuine copied-code and copied-source corruptions in both Python modes.

No mutated source is accepted as physical input.  Authentication rejections and
mathematical formula rejections are reported as separate categories.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def need(value,label):
    if not value:
        raise RuntimeError(label)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(),'OUTPUT_EXISTS_OR_SYMLINK')
    need(args.output.parent.is_dir(),'OUTPUT_PARENT_MISSING')
    here=Path(__file__).resolve().parent
    code=(here/'check_independent_transfer_algebra.py').read_text()
    pin_bytes=(here/'SOURCE_PINS.json').read_bytes()
    pins=json.loads(pin_bytes)
    root=args.repo_root.resolve(strict=True)
    mutants=[
        ('canonical_factor_two','actual = (2*b/density)','actual = (b/density)','CANONICAL_FULL_SECOND_MAP'),
        ('canonical_density','actual = (2*b/density)','actual = (2*b)','CANONICAL_FULL_SECOND_MAP'),
        ('pullback_first_velocity','t*t*(a+c*y)','t*(a+c*y)','PULLBACK_ANGULAR_VELOCITY'),
        ('projection_angular_weight',"sphere_mean(integrand,xs) == 0, 'SPHERE_HESSIAN_PROJECTION'",
         "sphere_mean(integrand*(1+xs[0]**2),xs) == 0, 'SPHERE_HESSIAN_PROJECTION'",'SPHERE_HESSIAN_PROJECTION'),
        ('radial_conversion_sign',"s.Symbol('Hy'):-Uz/rho","s.Symbol('Hy'):Uz/rho",'UNIVERSAL_TT_COORDINATE_CONVERSION'),
        ('screen_clock_power',"h = 3*s.Symbol('Ux')/(2*ell*rho**2*c)",
         "h = 3*s.Symbol('Ux')/(2*ell*rho*c)",'SCREEN_CLOCK_FACTOR'),
        ('direct_Ricci_tensor_curvature','-8*h/(r**2*c**2)','-6*h/(r**2*c**2)','DIRECT_WARPED_RICCI_TT_22'),
    ]
    before={row['path']:hashlib.sha256((root/row['path']).read_bytes()).hexdigest() for row in pins['sources']}
    cases=[]
    with tempfile.TemporaryDirectory(prefix='hdblast-independent-exact-TT-controls-') as temp_name:
        temp=Path(temp_name)
        for optimized in (False,True):
            mode='optimized' if optimized else 'normal'
            for name,old,new,error in mutants:
                need(code.count(old)==1,'MUTATION_ANCHOR:'+name)
                case=temp/(mode+'-'+name)
                case.mkdir()
                script=case/'check_independent_transfer_algebra.py'
                script.write_text(code.replace(old,new))
                (case/'SOURCE_PINS.json').write_bytes(pin_bytes)
                cases.append((mode,name,'mathematical_formula',script,root,case/'success.json',error))

        # Build one sparse public-source checkout per mode. It contains every
        # selected source plus the complete retained 51-member proof capsule.
        relbase=Path(pins['retained_capsule']['path'])
        manifest=json.loads((root/relbase/'MANIFEST.json').read_text())
        relative={Path(row['path']) for row in pins['sources']}
        relative.add(relbase/'MANIFEST.json')
        relative.update(relbase/row['path'] for row in manifest['files'])
        source_targets=[
            ('selected_source_bytes',Path(pins['sources'][0]['path']),'SELECTED_SOURCE_HASH:'),
            ('retained_transitive_proof_bytes',relbase/'sol/linear_audit/LINEAR_VARIABLE_WAVE_AUDIT.md','RECOVERED_SOURCE_HASH:'),
        ]
        for optimized in (False,True):
            mode='optimized' if optimized else 'normal'
            for name,target,error in source_targets:
                case=temp/(mode+'-'+name)
                case.mkdir()
                sparse=case/'repo'
                for rel in relative:
                    destination=sparse/rel
                    destination.parent.mkdir(parents=True,exist_ok=True)
                    shutil.copyfile(root/rel,destination)
                mutated=sparse/target
                data=mutated.read_bytes()
                need(bool(data),'EMPTY_MUTATION_TARGET')
                mutated.write_bytes(bytes([data[0]^1])+data[1:])
                script=case/'check_independent_transfer_algebra.py'
                script.write_text(code)
                (case/'SOURCE_PINS.json').write_bytes(pin_bytes)
                cases.append((mode,name,'source_authentication',script,sparse,case/'success.json',error))

        def run(case):
            mode,name,category,script,source_root,output,error=case
            command=[sys.executable,'-I','-B']
            if mode=='optimized': command.append('-O')
            command += [str(script),'--repo-root',str(source_root),'--output',str(output)]
            result=subprocess.run(command,text=True,capture_output=True,timeout=240)
            need(result.returncode!=0,'CORRUPTED_EXECUTION_ACCEPTED:'+mode+':'+name)
            need(error in result.stderr,'WRONG_REJECTION:'+mode+':'+name+':'+result.stderr[-1000:])
            need(not output.exists(),'CORRUPTED_SUCCESS_ARTIFACT:'+mode+':'+name)
            return {'mode':mode,'name':name,'category':category,'expected_failure':error,
                    'actual_returncode':result.returncode,
                    'actual_stderr_witness':[line.strip() for line in result.stderr.splitlines() if error in line][-1],
                    'executed_code_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),
                    'returncode_nonzero':True,'success_artifact_absent':True}

        with ThreadPoolExecutor(max_workers=4) as pool:
            results=list(pool.map(run,cases))
    for rel,digest in before.items():
        need(hashlib.sha256((root/rel).read_bytes()).hexdigest()==digest,'ORIGINAL_SOURCE_CHANGED:'+rel)
    report={'schema':'hdblast-independent-exact-TT-controls-v1',
            'original_code_sha256':hashlib.sha256(code.encode()).hexdigest(),
            'mathematical_formula_subprocess_rejections':sum(r['category']=='mathematical_formula' for r in results),
            'source_authentication_subprocess_rejections':sum(r['category']=='source_authentication' for r in results),
            'all_success_artifacts_absent':True,'original_selected_sources_unchanged':True,
            'controls':sorted(results,key=lambda r:(r['category'],r['name'],r['mode'])),
            'scope':'actual copied-code/source mutations; formula tests do not certify analytic existence'}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='controls'},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
