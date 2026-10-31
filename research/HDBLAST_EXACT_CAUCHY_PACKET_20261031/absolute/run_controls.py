#!/usr/bin/env python3
"""Execute actual corrupted reports/sources against the production checker."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent

def need(ok,why):
    if not ok:raise RuntimeError(why)

def norm(row,weights):
    ss=sum(w*x*x for w,x in zip(weights,row));den=1<<128;k=math.isqrt(ss.numerator*den*den//ss.denominator)
    if Q(k*k,den*den)<ss:k+=1
    return Q(k,den)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--report',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args();need(not args.output_dir.exists(),'fresh controls directory required');args.output_dir.mkdir(parents=True)
    original=json.loads(args.report.read_bytes());cases=[]
    def altered(name,expected,mutate,witness):
        value=copy.deepcopy(original);mutate(value);path=args.output_dir/(name+'.json');path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
        cases.append({'name':name,'expected':expected,'report':path,'root':args.repo_root,'witness':witness})
    def node(r,identifier):return next(n for n in r['nodes'] if n['id']==identifier)
    altered('coverage_gap','exact parent-child subset identity',lambda r:node(r,'s5L')['r'].__setitem__(1,str(Q(node(r,'s5L')['r'][1])-Q(1,1024))),'end moved left by1/1024')
    altered('coverage_overlap','exact parent-child subset identity',lambda r:node(r,'s5L')['r'].__setitem__(1,str(Q(node(r,'s5L')['r'][1])+Q(1,1024))),'end moved right by1/1024')
    altered('missing_component','all23 complete3-jet table',lambda r:node(r,'s5L')['raw_component_caps'][2].pop(),'rr row reduced from23 to22 components')
    altered('wrong_gram','exact all23 Gram weights',lambda r:r['Gram'].__setitem__(10,str(Q(r['Gram'][10])*2)),'positive Gram component10 doubled')
    altered('wrong_label','exact all23 labels',lambda r:r['labels'].__setitem__(10,'changed:w:2'),'component10 registration label changed')
    altered('invalid_parent_min','valid parent subset component minimum',lambda r:node(r,'s5L')['effective_component_caps'][2].__setitem__(9,str(Q(node(r,'s5L')['effective_component_caps'][2][9])+1)),'effective rr component9 increased by1 without justified parent/raw minimum')
    altered('understated_leaf_H','all23 node Gram aggregation',lambda r:node(r,'s5L')['absolute_H_caps'].__setitem__(2,str(Q(node(r,'s5L')['absolute_H_caps'][2])/2)),'one completed leaf rr norm halved')
    altered('understated_global_H','complete all-leaf absolute aggregation',lambda r:r['candidate_initial_absolute_H_caps'].__setitem__('u_rr',str(Q(r['candidate_initial_absolute_H_caps']['u_rr'])/2)),'global rr upper halved without changing leaves')
    altered('changed_initial_error','frozen complete initial errors',lambda r:r['initial_error_H_caps'].__setitem__('e0_rr',str(Q(r['initial_error_H_caps']['e0_rr'])/2)),'same-candidate initial rr error charge halved')
    altered('changed_velocity','unchanged complete velocity',lambda r:r['velocity']['absolute_H_caps'].__setitem__('v',str(Q(r['velocity']['absolute_H_caps']['v'])/2)),'complete velocity absolute charge halved')
    altered('omitted_tail','complete ordered leaf roster',lambda r:r['leaves'].remove('s1'),'analytic tail leaf[1,65/64] removed from coverage roster')
    def raw_understatement(r):
        n=node(r,'s5L');par=node(r,'s5');k=2;j=9
        need(Q(n['raw_component_caps'][k][j])>0,'nonzero raw perturbation witness')
        n['raw_component_caps'][k][j]=str(Q(n['raw_component_caps'][k][j])/2)
        n['effective_component_caps']=[[str(min(Q(x),Q(y))) for x,y in zip(rr,old)] for rr,old in zip(n['raw_component_caps'],par['effective_component_caps'])]
        n['absolute_H_caps']=[str(norm(list(map(Q,row)),list(map(Q,r['Gram'])))) for row in n['effective_component_caps']]
    altered('raw_formula_understatement','completed exact numerical raw cap s5L',raw_understatement,'raw completed w:0 rr cap halved and all dependent local minima/norms recomputed; numerical evaluator must reject')
    caprel=Path('research/HDBLAST_EXACT_PACKET_INITIAL_20261031')
    for name,extra in [('source_formula_tamper',False),('source_unlisted_member',True)]:
        root=args.output_dir/name;shutil.copytree(args.repo_root/caprel,root/caprel)
        if extra:
            (root/caprel/'unlisted.py').write_text('x=1\n');expected='closed source capsule';witness='unlisted source member added to closed capsule'
        else:
            path=root/caprel/'checks/producer.py';data=path.read_text();needle='scales = [Q(1),self.ell,2*self.ell**2]'
            need(data.count(needle)==1,'source perturbation anchor');path.write_text(data.replace(needle,'scales = [Q(1),self.ell,3*self.ell**2]'))
            expected='member digest';witness='rr jet factorial conversion changed2 to3 on a copied producer'
        cases.append({'name':name,'expected':expected,'report':args.report,'root':root,'witness':witness})
    def execute(case,optimized):
        mode='optimized' if optimized else 'normal';log=args.output_dir/(case['name']+'-'+mode);output=log.with_suffix('.success.json')
        cmd=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[str(HERE/'check_refinement.py'),'--repo-root',str(case['root']),'--report',str(case['report']),'--output',str(output)]
        proc=subprocess.run(cmd,capture_output=True,timeout=180)
        log.with_suffix('.stdout').write_bytes(proc.stdout);log.with_suffix('.stderr').write_bytes(proc.stderr)
        stderr=proc.stderr.decode(errors='replace')
        need(proc.returncode!=0 and not output.exists() and case['expected'] in stderr,'corruption escaped expected failure: '+case['name']+' '+mode+' '+stderr[-1000:])
        return {'case':case['name'],'mode':mode,'returncode':proc.returncode,'successful_output_absent':not output.exists(),'failure':case['expected'],'nonzero_mutation_witness':case['witness'],'stderr_sha256':hashlib.sha256(proc.stderr).hexdigest()}
    with ThreadPoolExecutor(max_workers=4) as pool:
        jobs=[pool.submit(execute,c,optimized) for c in cases for optimized in (False,True)]
        receipts=[job.result() for job in jobs]
    result={'status':'PASS_REAL_CORRUPTION_SUBPROCESSES','case_count':len(cases),'execution_count':len(receipts),'receipts':receipts,'checks_use_runtime_need_not_assertions':True}
    (args.output_dir/'CONTROLS.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':result['status'],'cases':len(cases),'executions':len(receipts)}))
if __name__=='__main__':main()
