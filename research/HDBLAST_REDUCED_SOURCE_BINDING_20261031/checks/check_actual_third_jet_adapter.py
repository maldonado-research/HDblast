#!/usr/bin/env python3
"""Exact coefficient and actual-input audit of the additive third-jet adapter."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
ADAPTER_PIN='4847c79f034a09f423cc020ab3bef86e03d729970084a9bb088e94caa3b9b18a'
MUTATIONS=('third_radial_sign','omit_middle_radial_term','omit_background_radial_term',
           'wrong_background_radial_sign','omit_yyt_rho_term','wrong_ytt_sign',
           'lose_third_time_ell','raw_coefficient_as_derivative','omit_new_radial_slot',
           'replace_actual_fourth_caps_by_zero')
EXPECTED_FAILURE={
 'third_radial_sign':'ordinary_physical_conversion_yyy',
 'omit_middle_radial_term':'ordinary_physical_conversion_yyy',
 'omit_background_radial_term':'ordinary_physical_conversion_yyy',
 'wrong_background_radial_sign':'ordinary_physical_conversion_yyy',
 'omit_yyt_rho_term':'ordinary_physical_conversion_yyt',
 'wrong_ytt_sign':'ordinary_physical_conversion_ytt',
 'lose_third_time_ell':'ordinary_physical_conversion_ttt',
 'raw_coefficient_as_derivative':'raw_normalized_Taylor_factorial_(3, 0)',
 'omit_new_radial_slot':'all_eight_new_pure_radial_slots',
 'replace_actual_fourth_caps_by_zero':'retained_actual_third_slot_Attt',
}


def need(ok,why):
    if not ok:raise ValueError(why)


class P:
    def __init__(self,v=0):
        self.d={k:F(x) for k,x in v.items() if x} if isinstance(v,dict) else ({():F(v)} if v else {})
    @staticmethod
    def cv(v):return v if isinstance(v,P) else P(v)
    @staticmethod
    def var(n):return P({((n,1),):1})
    def __add__(self,v):
        d=dict(self.d)
        for k,x in P.cv(v).d.items():d[k]=d.get(k,F(0))+x
        return P(d)
    __radd__=__add__
    def __neg__(self):return P({k:-x for k,x in self.d.items()})
    def __sub__(self,v):return self+-P.cv(v)
    def __rsub__(self,v):return P.cv(v)+-self
    def __mul__(self,v):
        d={}
        for k,x in self.d.items():
            for j,y in P.cv(v).d.items():
                ex=dict(k)
                for n,p in j:ex[n]=ex.get(n,0)+p
                key=tuple(sorted((n,p) for n,p in ex.items() if p))
                d[key]=d.get(key,F(0))+x*y
        return P(d)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:
            need(len(self.d)==1,'inverse of nonmonomial')
            (k,x),=self.d.items();return P({tuple((v,p*n) for v,p in k):x**n})
        out=P(1)
        for _ in range(n):out=out*self
        return out
    def __truediv__(self,v):return self*P.cv(v)**-1
    def __rtruediv__(self,v):return P.cv(v)*self**-1
    def deriv(self,rules):
        out=P()
        for key,c in self.d.items():
            for name,n in key:
                if name not in rules:continue
                ex=dict(key);ex[name]-=1
                out+=P({tuple(sorted((v,p) for v,p in ex.items() if p)):c*n})*rules[name]
        return out
    def encode(self):return [{'coefficient':str(c),'powers':dict(k)} for k,c in sorted(self.d.items())]


def load_adapter():
    p=HERE/'actual_first_third_jets.py';raw=p.read_bytes()
    need(p.is_file() and not p.is_symlink() and hashlib.sha256(raw).hexdigest()==ADAPTER_PIN,'adapter authentication')
    for path in (p,Path(__file__)):
        need(not any(isinstance(x,ast.Assert) for x in ast.walk(ast.parse(path.read_text()))),'assertion-free checker')
    spec=importlib.util.spec_from_file_location('authenticated_new_third_jet_adapter',p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def calculate(root,mutation):
    new=load_adapter();checks=[]
    def eq(name,left,right):
        remainder=P.cv(left)-right
        need(not remainder.d,'MATH_MISMATCH:'+name+':'+json.dumps(remainder.encode(),sort_keys=True,separators=(',',':')))
        checks.append({'name':name,'zero_collected_coefficients':True})
    ell,rho,a,ay=(P.var(n) for n in ('ell','rho','a','ay'))
    f=lambda i,j:P.var('F_'+str(i)+'_'+str(j))
    def Dr(v):
        rules={'rho':-ell*rho*rho*a,'a':-ell*rho*ay}
        for key in P.cv(v).d:
            for name,_ in key:
                if name.startswith('F_'):
                    _,i,j=name.split('_');rules[name]=f(int(i),int(j)+1)
        return P.cv(v).deriv(rules)
    def Dt(v):
        rules={}
        for key in P.cv(v).d:
            for name,_ in key:
                if name.startswith('F_'):
                    _,i,j=name.split('_');rules[name]=f(int(i)+1,int(j))/ell
        return P.cv(v).deriv(rules)
    Dy=lambda v:-Dr(v)/(ell*rho)
    class Mock:
        def __init__(self,value):self.value=value
        def val(self):return self.value
        def derivative(self,i,j):return f(i,j)
    actual=new.physical_through_three(Mock(f(0,0)),Mock(rho),Mock(a),ay,ell,mutation)
    raw={'':f(0,0)}
    for suffix in ('y','t','yy','yt','tt','yyy','yyt','ytt','ttt'):
        v=f(0,0)
        for axis in suffix:v=Dy(v) if axis=='y' else Dt(v)
        raw[suffix]=v
    for suffix in raw:eq('ordinary_physical_conversion_'+(suffix or 'value'),actual[suffix],raw[suffix])
    for word in ('yt','ty','yyt','yty','tyy','ytt','tyt','tty'):
        v=f(0,0)
        for axis in word:v=Dy(v) if axis=='y' else Dt(v)
        eq('physical_derivative_commutation_'+word,v,raw[''.join(sorted(word,reverse=True))])
    expected_yyy=-f(0,3)/(ell*rho)**3-3*a*f(0,2)/(ell*rho)**2+(ay-a*a)*f(0,1)/(ell*rho)
    eq('independent_pure_radial_third_formula',actual['yyy'],expected_yyy)
    # Wrong omitted terms have exact generic nonzero witnesses.
    witnesses={'middle_radial':(-3*a*f(0,2)/(ell*rho)**2).encode(),
               'background_radial':((ay-a*a)*f(0,1)/(ell*rho)).encode()}
    need(all(witnesses.values()),'nonzero product-rule witnesses')
    adapter=new.Adapter(root);I,T=adapter.I,adapter.T
    if mutation=='raw_coefficient_as_derivative':T.derivative=lambda self,i,j:self.c[i,j]
    fixture=T({(3,0):F(2,7),(1,2):F(3,11),(0,3):F(5,13)})
    for i,j,value in ((3,0,F(12,7)),(1,2,F(6,11)),(0,3,F(30,13))):
        interval=fixture.derivative(i,j)
        need(interval.lo<=value<=interval.hi,'MATH_MISMATCH:raw_normalized_Taylor_factorial_'+str((i,j)))
        checks.append({'name':'raw_normalized_Taylor_factorial_'+str((i,j)),'expected_exact_derivative':str(value),'contained':True})
    if mutation=='replace_actual_fourth_caps_by_zero':
        adapter.master_caps[4]=F(0)
        for i,j in new.KEYS:
            if i+j==4:adapter.gauge_caps[i,j]=[F(0)]*3
        adapter.global_third,adapter.coefficient_boxes=adapter._global_third()
    previous=json.loads((Path(root)/'research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json').read_text())
    need(previous['status']=='PASS_COMPLETE_ORIGINAL_ACTUAL_FULL23_SCALED_BULK_SOURCE_TIME_DERIVATIVE'
         and previous['input_sha256']['master4_report']==new.SELECTED['HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json']
         and previous['input_sha256']['gauge4_report']==new.SELECTED['HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/GAUGE4.json'],
         'same original actual higher-jet reference dependencies')
    for name,row in previous['new_actual_field_third_jet_boxes'].items():
        need(adapter.global_third[name].row()==row,'MATH_MISMATCH:retained_actual_third_slot_'+name)
        checks.append({'name':'retained_actual_third_slot_'+name,'complete_interval_bytes_equal':True})
    # Coefficient order four is unavailable to this lifted calculation and
    # deliberately unused by amplitude outputs through order three.
    sentinel,_=adapter._global_third(unused_fourth_sentinel=F(7,11))
    for name,interval in adapter.global_third.items():
        need(interval.row()==sentinel[name].row(),'MATH_MISMATCH:unused_coefficient_order_four_'+name)
        checks.append({'name':'unused_coefficient_order_four_'+name,'nonzero_sentinel_invariant':True})
    if mutation=='omit_new_radial_slot':adapter.global_third.pop('Ayyy')
    need(len(adapter.global_third)==32 and all(k+'yyy' in adapter.global_third for k in ('A','B','C','D','E','L','M','s')),
         'MATH_MISMATCH:all_eight_new_pure_radial_slots')
    checks.append({'name':'all_eight_new_pure_radial_slots','metric':7,'scalar':1})
    cells=[]
    for label,x,r,clip in [('interior_cell',I(F(1,4),F(1,4)+F(1,4096)),I(F(5,4),F(5,4)+F(1,4096)),False),
                           ('physical_shell',I(F(3,2)),I(0),False),
                           ('clipped_causal_edge',I(F(5,2),F(3)),I(F(5),F(11,2)),True)]:
        m,s,bg=adapter.actual_jet(x,r,clip_causal=clip)
        oldm,olds,oldbg=adapter.old.actual_jet(x,r,clip_causal=clip)
        need(all(m[k].row()==v.row() for k,v in oldm.items()) and all(s[k].row()==v.row() for k,v in olds.items()),
             'MATH_MISMATCH:retained_actual_lower_two_jets_'+label)
        need(set(m)=={p+suffix for p in 'ABCDELM' for suffix in raw} and set(s)=={'s'+suffix for suffix in raw},'complete actual same-field interface')
        cells.append({'label':label,'scaled_x':x.row(),'scaled_r':r.row(),'retained_lower_two_jet_bytes_equal':True,
                      'actual_domain_clipped_to_nonempty_intersection':clip,
                      'actual_metric_jet_slots':len(m),'actual_scalar_jet_slots':len(s),'ordinary_third':new.encode({k:v for k,v in {**m,**s}.items() if k[-3:] in new.THIRD_SUFFIXES})})
        checks.append({'name':'retained_actual_lower_two_jets_'+label,'unchanged':True})
    invalid=[]
    for label,x,r in [('outside_time',I(4),I(0)),('unclipped_crossing',I(F(5,2),3),I(5,F(11,2))),('empty_causal_intersection',I(3),I(6))]:
        try:adapter.actual_jet(x,r)
        except ValueError as err:
            need(any(term in str(err) for term in ('nonempty cell','whole cell inside')),'unexpected domain failure')
            invalid.append(label)
        else:raise ValueError('MATH_MISMATCH:invalid_domain_accepted_'+label)
    checks.append({'name':'actual_domain_guard_rejections','rejected':invalid})
    adapter.readback()
    need(hashlib.sha256((HERE/'actual_first_third_jets.py').read_bytes()).hexdigest()==ADAPTER_PIN,'adapter unchanged after readback')
    return {'status':'PASS_ACTUAL_FIRST_FIELDS_COHERENT_ORDINARY_THIRD_JET_ADAPTER',
            'source_authentication':adapter.before,'new_adapter_sha256':ADAPTER_PIN,
            'arithmetic':'exact rational Laurent identities and outward rational interval enclosures; no point sampling proof',
            'domain':'0<=x<=3,r>=0,x+r<=8; physical-side derivatives and one-sided shell traces',
            'ell':str(adapter.ell),'checks':checks,'independent_physical_conversion_coefficients':{k:v.encode() for k,v in raw.items()},
            'nonzero_chain_rule_witnesses':witnesses,'actual_raw_master_scaled_order_caps':new.encode(adapter.master_caps),
            'actual_raw_gauge_scaled_mixed_caps':{str(i)+','+str(j):new.encode(v) for (i,j),v in sorted(adapter.gauge_caps.items())},
            'actual_coefficient_jet_boxes_through_three':adapter.coefficient_boxes,
            'actual_global_ordinary_third_jet_boxes':new.encode(adapter.global_third),'cell_interface_checks':cells,
            'scope':{'same_actual_first_fields':True,'actual_metric_slots_through_order_three':70,'actual_scalar_slots_through_order_three':10,
                     'new_pure_radial_third_metric_and_scalar_slots':8,'third_jets_are_uniform_coarse_magnitudes':True,
                     'retained_lower_two_jets_unchanged':True,'derivatives_of_independent_interval_selections':False,
                     'new_actual_higher_regularity_theorem':False,'earlier_raw_higher_jet_certificates_replayed_here':False,
                     'actual_background_master_gauge_and_corner_premises_retained':True,'improved_error_bound':False,
                     'second_order_solution_or_nonlinear_family_proved':False,'physical_force_or_cosmological_evidence':False}}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--mutation',choices=('none',)+MUTATIONS,default='none');ap.add_argument('--run-controls',action='store_true');args=ap.parse_args()
    need(not args.output.exists(),'fresh output required')
    report=calculate(args.repo_root.resolve(),args.mutation)
    if args.run_controls:
        controls=[]
        with tempfile.TemporaryDirectory(prefix='actual-third-jet-controls-') as tmp:
            for mutation in MUTATIONS:
                path=Path(tmp)/(mutation+'.json');cmd=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])
                cmd += [str(Path(__file__).resolve()),'--repo-root',str(args.repo_root.resolve()),'--output',str(path),'--mutation',mutation]
                done=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
                marker='ValueError: MATH_MISMATCH:'+EXPECTED_FAILURE[mutation]
                need(done.returncode!=0 and not path.exists() and marker in done.stderr,'corruption not mathematically rejected '+mutation+':'+done.stderr[-1000:])
                controls.append({'mutation':mutation,'returncode':done.returncode,'successful_output_absent':True,'expected_identity':EXPECTED_FAILURE[mutation],
                                 'mathematical_failure':done.stderr.split('ValueError: MATH_MISMATCH:',1)[1].splitlines()[0]})
        report['rejected_corrupted_executions']=controls
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'check_groups':len(report['checks']),'corruptions':len(report.get('rejected_corrupted_executions',[]))}))


if __name__=='__main__':main()
