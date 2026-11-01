#!/usr/bin/env python3
"""Strict subprocess controls: both modes, relocated sources/code, no writes there.

Two fresh relocated positive campaigns make 146 new actual oracle calls total.
Twenty strict rejection subprocesses use pinned raw rows and make no new calls.
"""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
PRODUCER_SHA='7832767495d93e779459019221bedf773266598549046190424a2ed02da79265'
REPORT_SHA='42949cc9f4dc1c2377ac5b6497d79712c0c0d8550a524ad04d83cb799f7df475'
def need(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--raw-report',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(),'FRESH_CONTROL_OUTPUT')
    need(sha(HERE/'check_shell_budget.py')==PRODUCER_SHA and sha(args.raw_report)==REPORT_SHA,'CONTROL_BASELINE_PINS')
    pins=json.loads((HERE/'SOURCE_PINS.json').read_bytes());positives=[];rejects=[]
    with tempfile.TemporaryDirectory(prefix='shell-budget-controls-') as temp:
        tmp=Path(temp);code=tmp/'relocated-component';shutil.copytree(HERE,code)
        repo=tmp/'relocated-source';(repo/'research').mkdir(parents=True)
        for name in pins['capsules']:shutil.copytree(args.repo_root/'research'/name,repo/'research'/name)
        geo=repo/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
        member=geo/'numeric/actual_first_jets.py';original=member.read_bytes()
        freeze=code/'FREEZE.json';freeze_original=freeze.read_bytes()
        pilot=code/'protocol/TRACE_COARSE_PILOT.json';pilot_original=pilot.read_bytes()
        def command(mode,out,raw=True,mutation=None):
            cmd=[sys.executable,'-I','-S','-B']
            if mode=='optimized':cmd.append('-O')
            cmd += [str(code/'check_shell_budget.py'),'--repo-root',str(repo),'--output-dir',str(out)]
            if raw:cmd += ['--raw-report',str(args.raw_report.resolve())]
            if mutation:cmd += ['--inject-error',mutation]
            return cmd
        def reject(mode,case,tags,mutation=None):
            out=tmp/(mode+'-'+case)
            proc=subprocess.run(command(mode,out,True,mutation),cwd=tmp,capture_output=True,text=True)
            message=proc.stderr.strip()
            need(proc.returncode==1 and not out.exists() and not proc.stdout and any(message=='FAIL: '+tag for tag in tags),
                 'STRICT_REJECTION:'+mode+':'+case+':'+message)
            rejects.append({'mode':mode,'case':case,'returncode':1,'output_directory_absent':True,
                            'failure':message,'new_actual_oracle_calls':0})
        for mode in ['normal','optimized']:
            out=tmp/(mode+'-fresh-positive')
            proc=subprocess.run(command(mode,out,False),cwd=tmp,capture_output=True,text=True)
            report=out/'ACTUAL_SHELL_BUDGET.json'
            need(proc.returncode==0 and not proc.stderr and report.is_file() and sha(report)==REPORT_SHA,'RELOCATED_FRESH_POSITIVE:'+mode+':'+proc.stderr)
            info=json.loads(proc.stdout);need(info['actual_new_oracle_calls']==73,'POSITIVE_ORACLE_ACCOUNT')
            positives.append({'mode':mode,'code_relocated':True,'source_capsules_relocated':True,
                              'fresh_oracle_campaign':True,'new_actual_oracle_calls':73,'report_sha256':REPORT_SHA})
            member.write_bytes(original+b"\nraise RuntimeError('SCIENTIFIC_IMPORT_BEFORE_AUTHENTICATION')\n")
            reject(mode,'source_member',['MEMBER:HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py'])
            member.write_bytes(original)
            extra=geo/'UNLISTED_CONTROL_MEMBER';extra.write_text('unlisted\n')
            reject(mode,'unlisted_source_member',['CLOSED_INVENTORY']);extra.unlink()
            member.unlink();member.symlink_to(args.repo_root.resolve()/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py')
            reject(mode,'source_symlink',['MEMBER:HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py']);member.unlink();member.write_bytes(original)
            freeze.write_bytes(freeze_original+b'\n');reject(mode,'freeze',['FREEZE']);freeze.write_bytes(freeze_original)
            pilot.write_bytes(pilot_original+b'\n');reject(mode,'pilot',['PROTOCOL_MEMBER:protocol/TRACE_COARSE_PILOT.json']);pilot.write_bytes(pilot_original)
            for mutation in ['error_gain','kernel_order','clock_scale','P_coefficient']:
                reject(mode,mutation,['FROZEN_SCIENTIFIC_WITNESS','FROZEN_MAIN_SIGNS'],mutation)
            reject(mode,'gauss_coefficient',['TRACE_FORCE_IDENTITIES_COMPATIBLE'],'gauss_coefficient')
    result={'schema':'HDBLAST_SHELL_BUDGET_STRICT_SUBPROCESS_CONTROLS_V1','status':'PASS_STRICT_RELOCATED_FRESH_POSITIVES_AND_NEGATIVE_CONTROLS',
            'producer_sha256':PRODUCER_SHA,'primary_report_sha256':REPORT_SHA,
            'strict_rejection_subprocesses':len(rejects),'positive_subprocesses':len(positives),
            'new_actual_oracle_calls_in_positive_subprocesses':sum(p['new_actual_oracle_calls'] for p in positives),
            'new_actual_oracle_calls_in_rejection_subprocesses':0,'inprocess_mathematical_mutations_counted_as_subprocesses':False,
            'positives':positives,'rejections':rejects,
            'witness_role':'The pinned scientific digest locks the already frozen, independently reconstructable mathematical report. Mutations must fail its intended witness; source and coefficient corruptions fail their earlier named authentication or identity gates.'}
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'strict_rejections':len(rejects),'positive_subprocesses':len(positives),'extra_actual_oracle_calls':146,'sha256':sha(args.output)}))
if __name__=='__main__':
    try:main()
    except (ValueError,OSError) as exc:
        print('FAIL: '+str(exc),file=sys.stderr);raise SystemExit(1)
