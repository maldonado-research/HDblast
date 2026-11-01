#!/usr/bin/env python3
"""Independent portable exact-rational aggregation; no producer imports."""
import argparse,hashlib,json,sys
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
PINS_SHA='9de076a120a0eb5b84afc45e626f145f129a0c1bddb87a1a0921ae9c3441155e'
FREEZE_SHA='67f33c24285cf252b7dde4f40da5b71b1fafde4c41bc1836effba7b6541ed3bc'
REFERENCE_SHA='71d200264d7ba2bc7cfed87e8c1abd6978ab8cd92571c10e1e34a55670cef3f5'
PRODUCER_SHA='7f96dd6689b382642bc1a5b2ff9dfceb76dc4d5f018bfa1dbc4557bf55e1221e'
REPORT_SHA='3a4597ca94e42397a35f46c92dc99a30788f0b00afae3c316c7a58c1868fc4c5'
class Failure(Exception):pass
def need(b,witness):
    if not b:raise Failure(witness)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def read(p):return json.loads(p.read_bytes())
def regular(p):return p.is_file() and not p.is_symlink() and all(not x.is_symlink() for x in p.absolute().parents)
def pair(x):
    need(isinstance(x,list) and len(x)==2,'INTERVAL_PAIR');a,b=map(F,x);need(a<=b,'INTERVAL_ORDER');return a,b
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)
def inv(a):need(a[0]>0 or a[1]<0,'NONZERO_DENOMINATOR');return 1/a[1],1/a[0]
def div(a,b):return mul(a,inv(b))
def sq(a):return F(0) if a[0]<=0<=a[1] else min(x*x for x in a),max(x*x for x in a)
def scale(c,a):return mul((F(c),F(c)),a)
def enc(a):return [str(x) for x in a]
def contains(stored,exact):
    a=pair(stored);return a[0]<=exact[0]<=exact[1]<=a[1]
def authenticate(repo,pins):
    result={}
    for label,h in sorted(pins['capsules'].items()):
        root=repo/'research'/label;manifest=root/'MANIFEST.json'
        need(root.is_dir() and not root.is_symlink() and regular(manifest) and sha(manifest)==h,'SOURCE_MANIFEST_BYTES')
        seen={'MANIFEST.json'}
        for row in read(manifest)['files']:
            rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'SOURCE_SAFE_MEMBER')
            q=root/rel;need(regular(q) and q.stat().st_size==row['bytes'] and sha(q)==row['sha256'],'SOURCE_MEMBER_BYTES')
            need(pins['all_manifest_member_pins'][label+'/'+str(rel)]=={'bytes':row['bytes'],'sha256':row['sha256']},'SOURCE_FLAT_PIN');seen.add(str(rel))
        actual=set()
        for q in root.rglob('*'):
            need(not q.is_symlink(),'SOURCE_SYMLINK')
            if q.is_file():actual.add(str(q.relative_to(root)))
            else:need(q.is_dir(),'SOURCE_SPECIAL_FILE')
        need(actual==seen,'SOURCE_CLOSED_INVENTORY');result[label]={'manifest_sha256':h,'members_including_manifest':len(seen)}
    return result
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--report',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
    need(sys.version_info[:3]==(3,12,14),'PYTHON_VERSION')
    repo=a.repo_root.absolute();report=a.report.absolute();out=a.output_dir.absolute()
    need(all(not p.is_symlink() for p in [repo,*repo.parents,report,*report.parents,out,*out.parents]),'PATH_SYMLINK')
    repo=repo.resolve(strict=True);out=out.resolve()
    need(not out.exists(),'FRESH_OUTPUT_DIRECTORY')
    need(all(base!=out and base not in out.parents for base in [repo,HERE]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
    own=sha(Path(__file__));local={}
    for name,h in [('SOURCE_PINS.json',PINS_SHA),('FREEZE.json',FREEZE_SHA),('FROZEN_REFERENCE.json',REFERENCE_SHA),('producer.py',PRODUCER_SHA)]:
        p=HERE/name;need(regular(p) and sha(p)==h,'LOCAL_'+name);local[name]=h
    pins,freeze,reference=map(read,[HERE/'SOURCE_PINS.json',HERE/'FREEZE.json',HERE/'FROZEN_REFERENCE.json'])
    pilot=HERE/'PILOT_ONLY.json';need(regular(pilot) and sha(pilot)==freeze['exploratory_pilot_sha256'],'SEPARATED_PILOT')
    before=authenticate(repo,pins)
    need(regular(report),'REPORT_REGULAR');r=read(report);need(report.read_bytes()==canonical(r),'REPORT_CANONICAL_BYTES')
    need(r.get('producer_sha256')==PRODUCER_SHA,'PRODUCER_IDENTITY')
    need(r.get('source_pins_sha256')==PINS_SHA and r.get('freeze_sha256')==FREEZE_SHA and r.get('scientific_reference_sha256')==REFERENCE_SHA
         and r.get('source_authentication')==before and r.get('frozen_campaign')==freeze
         and r.get('first_derivative_zero') is True and r.get('scope')==reference['scope']
         and r.get('derivative_convention')==reference['derivative_convention'],'REPORT_METADATA')
    need(r.get('actual_oracle_calls')==9 and len(r.get('rows',[]))==9 and r.get('selected_target')==r['rows'][-1],'EXACT_CAMPAIGN')
    expected_m={k+s for k in ['A','B','C','D','E','L','M'] for s in ['','y','t','yy','yt','tt']}
    expected_s={'s','sy','st','syy','syt','stt'}
    expected_b={'rho','c','theta','a','ay','p','py','Vp','Vpp','Vppp','ell'}
    computed=[];lo,hi=F(31,25),F(63,50)
    for i,row in enumerate(r['rows']):
        xb=[lo+(hi-lo)*i/8,lo+(hi-lo)*(i+1)/8] if i<8 else [F(5,4)]*2
        need(list(map(F,row['scaled_x']))==xb and row['kind']==('continuous_cell' if i<8 else 'selected_target'),'FIXED_TIME_COVER')
        need(row['scaled_r']==['0','0'],'PHYSICAL_SHELL')
        m,s,b=row['actual_metric'],row['actual_scalar'],row['actual_background'];frozen=reference['rows'][i]
        need(set(m)==expected_m and set(s)==expected_s and set(b)==expected_b,'FULL_JET_INVENTORY')
        need(F(b['ell'])==F(1,10**6),'ELL_UNITS')
        need(m==frozen['actual_metric'],'METRIC_JET_REFERENCE')
        need(s==frozen['actual_scalar'],'SCALAR_JET_REFERENCE')
        need(b==frozen['actual_background'],'BACKGROUND_JET_REFERENCE')
        for k,v in {**m,**s,**b}.items():
            if k!='ell':pair(v)
        psi,psit,rho,c=pair(s['s']),pair(s['st']),pair(b['rho']),pair(b['c'])
        K=div(sq(psit),scale(3,sq(rho)))
        S=div(scale(8,sq(psi)),scale(3,mul(sq(rho),sq(c))))
        Q=sub(S,K);ratio=div(S,K);direct=div(add(neg(sq(psit)),div(scale(8,sq(psi)),sq(c))),scale(3,sq(rho)))
        need(contains(row['mean_K2'],K),'SECOND_DERIVATIVE_FACTOR')
        need(contains(row['mean_S2'],S),'ANGULAR_GRADIENT_FACTOR')
        need(contains(row['mean_Q2'],Q) and contains(row['direct_mean_Q2'],direct),'Q2_AGGREGATION')
        need(contains(row['spatial_to_temporal_ratio'],ratio),'RATIO_AGGREGATION')
        need(K[0]>0 and S[0]>0 and Q[1]<0 and ratio[1]<1 and direct[1]<0,'INDEPENDENT_STRICT_SIGN')
        d=pair(row['outward_integer_Q2']);q=pair(row['mean_Q2']);need(d[0]<=q[0]<=q[1]<=d[1]<0,'OUTWARD_INTEGER_DISPLAY')
        computed.append({'mean_K2':K,'mean_S2':S,'mean_Q2':Q,'spatial_to_temporal_ratio':ratio,'direct_mean_Q2':direct})
    continuous={k:(min(x[k][0] for x in computed[:8]),max(x[k][1] for x in computed[:8])) for k in ['mean_K2','mean_S2','mean_Q2','spatial_to_temporal_ratio']}
    need(all(contains(r['continuous_interval_hulls'][k],v) for k,v in continuous.items()),'CONTINUOUS_HULL_AGGREGATION')
    need(all(r[k]==reference[k] for k in reference['scientific_reference_keys']),'UNCHANGED_SCIENTIFIC_VALUES')
    need(sha(report)==REPORT_SHA,'REPORT_INTEGRITY')
    need(authenticate(repo,pins)==before and sha(Path(__file__))==own and all(sha(HERE/name)==h for name,h in local.items())
         and sha(pilot)==freeze['exploratory_pilot_sha256'],'AUDIT_POSTCHECK')
    result={'schema':'HDBLAST_PORTABLE_INDEPENDENT_SHELL_GRADIENT_AGGREGATION_V1','status':'PASS_INDEPENDENT_EXACT_RATIONAL_AGGREGATION',
        'report_sha256':REPORT_SHA,'producer_sha256':PRODUCER_SHA,'auditor_sha256':own,
        'source_pins_sha256':PINS_SHA,'freeze_sha256':FREEZE_SHA,'scientific_reference_sha256':REFERENCE_SHA,
        'scientific_or_producer_imports':False,'authenticated_source_members_before_and_after':True,
        'fixed_cells_verified':8,'selected_target_verified':True,'all_nine_strict_signs':True,
        'independent_exact_aggregate_rows':[{k:enc(v) for k,v in x.items()} for x in computed],
        'continuous_independent_exact_hulls':{k:enc(v) for k,v in continuous.items()},
        'limitation':'Checks source identity, fixed coverage, retained first jets and independent nonlinear aggregation; inherited PDE/oracle error theorems remain premises.'}
    out.mkdir(parents=True,exist_ok=False);target=out/'AUDIT.json';target.write_bytes(canonical(result))
    print(json.dumps({'status':'PASS','witness':'INDEPENDENT_EXACT_RATIONAL_AGGREGATION','audit_sha256':sha(target)},sort_keys=True))
if __name__=='__main__':
    try:main()
    except Failure as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
