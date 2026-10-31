#!/usr/bin/env python3
"""Portable CLI driver for real semantic-corruption subprocesses."""
import argparse,copy,hashlib,json,shutil,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path

def need(v,why):
 if not v:raise RuntimeError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-root',required=True,type=Path);ap.add_argument('--support-root',required=True,type=Path)
 ap.add_argument('--report',required=True,type=Path);ap.add_argument('--report-sha256',required=True);ap.add_argument('--output-dir',required=True,type=Path);a=ap.parse_args()
 need(not a.output_dir.exists(),'FRESH_CONTROL_OUTPUT_DIRECTORY');a.output_dir.mkdir(parents=True)
 checker=Path(__file__).with_name('check_report.py');checker_sha=sha(checker)
 def run(report,out,opt=False,root=None,support=None,expected=None):
  args=[sys.executable]+(['-O'] if opt else [])+[str(checker),'--source-root',str(root or a.source_root),'--support-root',str(support or a.support_root),
   '--report',str(report),'--report-sha256',expected or sha(report),'--output',str(out)]
  return subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 normals=[]
 for opt in (False,True):
  out=a.output_dir/('SEMANTICS_OPTIMIZED.json' if opt else 'SEMANTICS.json');p=run(a.report,out,opt,expected=a.report_sha256)
  need(p.returncode==0,p.stderr);normals.append(out)
 need(normals[0].read_bytes()==normals[1].read_bytes(),'NORMAL_OPTIMIZED_BYTE_IDENTITY')
 original=json.loads(a.report.read_text());c=original['campaigns'][0];seeds=[tuple(map(Q,p)) for p in c['seeds']]
 gap=next(i for i,n in enumerate(c['leaves']) if n['depth']>0 and Q(n['r'][1])<seeds[n['seed']][1])
 overlap=next(i for i,n in enumerate(c['leaves']) if n['depth']>0 and Q(n['r'][0])>seeds[n['seed']][0])
 def shift(d,i,sign):
  n=d['campaigns'][0]['leaves'][i];lo,hi=map(Q,n['r']);w=hi-lo;n['r']=[str(lo+sign*w),str(hi+sign*w)]
 def plus_component(d):
  n=d['campaigns'][0]['leaves'][gap];n['component_caps'][2][9]=str(Q(n['component_caps'][2][9])+1)
 cases=[('gap',lambda d:shift(d,gap,1),'COMPLETE_LEAF_COVER_NO_GAP_OR_OVERLAP'),
 ('overlap',lambda d:shift(d,overlap,-1),'COMPLETE_LEAF_COVER_NO_GAP_OR_OVERLAP'),
 ('dropped_component',lambda d:d['campaigns'][0]['leaves'][gap]['component_caps'][0].pop(),'FULL_3_BY_23'),
 ('wrong_Gram',lambda d:d['campaigns'][0]['Gram'].__setitem__(0,'2'),'SAME_LABEL_ORDER_AND_POSITIVE_GRAM'),
 ('inherited_cap',plus_component,'EXACT_PARENT_COMPONENT_INTERSECTION'),
 ('global_cap',lambda d:d['campaigns'][0]['complete_absolute_H_caps'].__setitem__('u_rr','0'),'EXACT_MAX_LEAF_GLOBAL_CAPS'),
 ('omitted_tail',lambda d:d['campaigns'][0].pop('tail_certificate'),'tail_certificate'),
 ('endpoint_trace',lambda d:d['campaigns'][0]['endpoint_rows'][6].__setitem__('quadratic_exterior_minus_interior_trace',[1,0,0]),'EXACT_QUADRATIC_MATCHING_TRACES'),
 ('dropped_endpoint',lambda d:d['campaigns'][0]['endpoint_rows'].pop(),'COMPLETE_ENDPOINT_COMPONENTS'),
 ('tail_recurrence',lambda d:d['campaigns'][0]['tail_certificate']['polynomials_ascending'][5].__setitem__(4,'361'),'TAIL_EXACT_POLYNOMIAL_RECURRENCE'),
 ('leaf_Gram_norm',lambda d:d['campaigns'][0]['leaves'][gap]['absolute_H_caps'].__setitem__(2,'0'),'EXACT_LEAF_POSITIVE_GRAM_SQRT'),
 ('denominator',lambda d:d['campaigns'][0]['nominal_denominators'].__setitem__('rho',['-1','1']),'NOMINAL_DENOMINATOR_INCLUSION_AND_NONZERO')]
 records=[]
 for name,mutate,marker in cases:
  d=copy.deepcopy(original);mutate(d);report=a.output_dir/(name+'_input.json');write(report,d)
  for opt in (False,True):
   out=a.output_dir/(name+('_O' if opt else '')+'_unexpected.json');p=run(report,out,opt)
   need(p.returncode!=0 and marker in p.stderr and not out.exists(),'EFFECTIVE_SEMANTIC_REJECTION:'+name+':'+p.stderr)
   records.append(dict(case=name,optimized=opt,expected_rejection=marker,rehashed_altered_report=True,nonzero_exit=True,no_success_artifact=True))
 # Alter a copied selected helper while the accepted report digest remains fixed.
 support=a.output_dir/'support_disk_control';support.mkdir()
 for rel in [*original['own_sha256'],'FREEZE.json','source-pins.json']:
  p=support/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.support_root/rel,p)
 helper=support/'scripts/pinned_independent_exact.py';helper.write_bytes(helper.read_bytes()+b'\n# deliberate source mutation\n')
 for opt in (False,True):
  out=a.output_dir/('source_disk'+('_O' if opt else '')+'_unexpected.json');p=run(a.report,out,opt,support=support,expected=a.report_sha256)
  marker='SOURCE_HASH:scripts/pinned_independent_exact.py';need(p.returncode!=0 and marker in p.stderr and not out.exists(),'EFFECTIVE_SOURCE_DISK_REJECTION')
  records.append(dict(case='source_disk',optimized=opt,expected_rejection=marker,nonzero_exit=True,no_success_artifact=True))
 need(sha(checker)==checker_sha,'CHECKER_UNCHANGED')
 report=dict(status='PASS_INDEPENDENT_ABSOLUTE_REPORT_SEMANTIC_CONTROLS',checker_sha256=checker_sha,accepted_report_sha256=a.report_sha256,
   canonical_semantic_report_sha256=sha(normals[0]),normal_optimized_byte_identical=True,controls=records,
   rehashed_semantic_corruption_subprocesses=24,disk_corruption_subprocesses=2,
   all_rejected_nonzero_without_success_artifact=True)
 write(a.output_dir/'CONTROLS.json',report);print(json.dumps({'status':report['status'],'rejected_subprocesses':len(records),'canonical_semantic_report_sha256':report['canonical_semantic_report_sha256']}))

if __name__=='__main__':main()
