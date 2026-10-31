#!/usr/bin/env python3
"""Independent exact arithmetic audit of postcontact *report data*.

This program imports no scientific producer or numerical model. It does not
authenticate the reported source values as evaluations of the physical model.
"""
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

GRID=1<<112
ROOT=Path(__file__).resolve().parent

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(v):
    if isinstance(v,Box):return list(map(str,v.ends))
    if isinstance(v,F):return str(v)
    if isinstance(v,dict):return {k:enc(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return list(map(enc,v))
    return v

class Box:
    """Rational interval with separately implemented directed lattice rounding."""
    def __init__(self,a=0,b=None):
        a=F(a); b=a if b is None else F(b)
        if a>b:raise ValueError('ORDER')
        self.ends=(F((a*GRID).__floor__(),GRID),F((b*GRID).__ceil__(),GRID))
    @staticmethod
    def val(x):return x if isinstance(x,Box) else Box(x)
    def __add__(self,x):
        x=self.val(x);return Box(self.ends[0]+x.ends[0],self.ends[1]+x.ends[1])
    __radd__=__add__
    def __neg__(self):return Box(-self.ends[1],-self.ends[0])
    def __sub__(self,x):return self+-self.val(x)
    def __rsub__(self,x):return self.val(x)+-self
    def __mul__(self,x):
        x=self.val(x);z=[a*b for a in self.ends for b in x.ends];return Box(min(z),max(z))
    __rmul__=__mul__
    def reciprocal(self):
        a,b=self.ends
        if a<=0<=b:raise ValueError('ZERO_DENOMINATOR')
        return Box(1/b,1/a)
    def __truediv__(self,x):return self*self.val(x).reciprocal()
    def __rtruediv__(self,x):return self.val(x)*self.reciprocal()
    def __pow__(self,n):
        if n<0:return self.reciprocal()**(-n)
        if n==0:return Box(1)
        if n%2:
            return self*(self**(n-1))
        a,b=(self**(n//2)).ends
        return Box(0 if a<=0<=b else min(a*a,b*b),max(a*a,b*b))
    def mag(self):return max(map(abs,self.ends))

def box(v):return Box(*v)
def rat(v):return F(v)

def audit(d):
    checks={}
    def check(name,value):
        checks[name]=bool(value)
        if not value:raise ValueError(name)
    def same(name,value,key):check(name,enc(value)==d[key])
    T,R=F(5,4),F(1,8)
    check('fixed_actual_postcontact_point',list(map(F,d['point']))==[T,R])
    bulk=d['bulk_cap_cover'];shell=d['shell_history_cover'];curve=d['characteristic_segments']
    check('bulk_cover_64_cells',len(bulk)==64)
    for k,row in enumerate(bulk):
        i,j=divmod(k,8)
        check('bulk_cell_%02d_geometry'%k,list(map(F,row['x']))==[T*i/8,T*(i+1)/8] and list(map(F,row['r']))==[(T+R)*j/8,(T+R)*(j+1)/8])
    f=max(box(v['F_T']).mag() for v in bulk)
    check('bulk_cap_exact_max',f==F(d['bulk_source_cap']))
    check('shell_cover_8_cells',len(shell)==8)
    for i,row in enumerate(shell):
        check('shell_cell_%02d_geometry'%i,list(map(F,row['x']))==[(T-R)*i/8,(T-R)*(i+1)/8])
    g=max(box(v['G_T']).mag() for v in shell)
    check('shell_cap_exact_max',g==F(d['shell_source_cap']))
    knots=sorted({T*i/32 for i in range(33)}|{T-R})
    check('characteristic_33_segments',len(curve)==33)
    integral=Box(0);length=F(0)
    for k,(lo,hi) in enumerate(zip(knots,knots[1:])):
        row=curve[k]
        check('segment_%02d_geometry'%k,list(map(F,row['x']))==[lo,hi] and list(map(F,row['positive_r']))==[R+T-hi,R+T-lo] and list(map(F,row['reflected_r']))==sorted([abs(R-T+lo),abs(R-T+hi)]))
        check('segment_%02d_length'%k,F(row['length'])==hi-lo)
        check('segment_%02d_nonnegative_physical_side'%k,min(map(F,row['reflected_r']))>=0)
        integral=integral+(hi-lo)*(box(row['positive_F'])+box(row['reflected_F']))/2
        length+=hi-lo
    check('characteristic_complete_length',length==T)
    check('reflection_join_is_a_knot',T-R in knots)
    same('integral_dyadic_aggregation',integral,'bulk_source_derivative_integral')
    seed=integral-box(d['shell_endpoint_G'])
    same('correct_minus_neumann_endpoint_seed',seed,'zero_bulk_feedback_seed_x')
    op=d['operator'];caps=op['caps']
    for k in ('A','B','D'):
        check('operator_%s_cap_exact_max'%k,box(op['coefficient_intervals'][k]).mag()==F(caps[k]))
    sigma=F(caps['A'])+F(caps['B'])+T*F(caps['D'])/2
    q=T*sigma
    check('feedback_sigma_exact',sigma==F(op['sigma']))
    check('feedback_q_exact',q==F(op['q']))
    check('feedback_strict_contraction',0<=q<1)
    norm=(2*g+T*f)/(1-q);tail=q*norm
    check('feedback_norm_exact',norm==F(d['full_target_X_upper']))
    check('complete_feedback_tail_exact',tail==F(d['complete_feedback_Ux_error']))
    ux=seed+Box(-tail,tail)
    same('actual_full_Ux_reconstructed',ux,'actual_full_U_Tx_interval')
    b={k:(box(v) if isinstance(v,list) else F(v)) for k,v in d['apex_background'].items()}
    m={k:box(v) for k,v in d['apex_actual_first_m'].items()}
    rho,c=b['rho'],b['c'];radius=rho*c
    scale=F(3,2)/(b['ell']*rho**2*c)
    check('screen_scale_positive',scale.ends[0]>0)
    same('screen_scale_reconstructed',scale,'screen_scale')
    fixed= ((m['D']+4*m['M']-m['L'])*m['Mt']-2*m['M']*m['Lt'])/(4*rho) -m['E']*(m['D']+m['L'])/(2*radius) +(m['C']/radius-m['My']/2)*(radius*(m['Ct']-b['theta']*m['C'])/(2*rho)-radius*m['Ey']/2+m['B']/2)
    same('all_fixed_quadratic_terms_reconstructed',fixed,'full_fixed_quadratic_F')
    quadratic=4*fixed/radius
    linear=F(3,2)*ux/(b['ell']*rho**2*c)
    same('quadratic_screen_conversion_reconstructed',quadratic,'screen_quadratic')
    same('linear_units_and_second_derivative_normalization',linear,'screen_linear')
    h=linear+quadratic
    same('total_screen_reconstructed',h,'actual_full_screen_h')
    same('opposite_orientation_reconstructed',-h,'actual_full_screen_minus')
    sign='positive' if h.ends[0]>0 else ('negative' if h.ends[1]<0 else 'inconclusive')
    check('reported_sign_exact',sign==d['sign'])
    check('sign_status_consistent',(sign!='inconclusive')==(d['status']=='CERTIFIED_CONDITIONAL_POSTCONTACT_SCREEN_SIGN'))
    charges=d['screen_error_charges']
    check('charged_feedback_scale_error',scale.ends[1]*tail==F(charges['feedback_scale_error_upper']))
    check('charged_seed_width',seed.ends[1]-seed.ends[0]==F(charges['seed_Ux_interval_width']))
    check('charged_quadratic_width',quadratic.ends[1]-quadratic.ends[0]==F(charges['quadratic_interval_width']))
    return {'checks':checks,'checks_count':len(checks),'h':h,'h_display_only':list(map(float,h.ends)),'q':q,'tail':tail,'bulk_cap':f,'shell_cap':g,'status':'PASS_REPORT_DATA_ALGEBRA_ONLY'}

def controls(d):
    variants={}
    a=copy.deepcopy(d);a['shell_endpoint_G']=['1','1'];variants['wrong_neumann_sign_or_endpoint_data']=a
    a=copy.deepcopy(d);a['characteristic_segments'][0]['reflected_r']=['-9/8','-139/128'];variants['unreflected_negative_side']=a
    a=copy.deepcopy(d);a['screen_quadratic']=['0','0'];variants['dropped_fixed_quadratic']=a
    a=copy.deepcopy(d);a['operator']['q']='0';variants['feedback_disabled']=a
    a=copy.deepcopy(d);a['screen_linear']=[str(F(v)/2) for v in a['screen_linear']];variants['extra_half_normalization']=a
    result={}
    for name,a in variants.items():
        try:audit(a)
        except ValueError as exc:result[name]={'rejected':True,'first_failed_witness':str(exc)}
        else:result[name]={'rejected':False}
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--producer-dir',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists() or args.output.is_symlink():raise ValueError('FRESH_OUTPUT')
    normal=args.producer_dir/'POSTCONTACT_NORMAL.json';optimized=args.producer_dir/'POSTCONTACT_OPTIMIZED.json'
    pins={p.name:digest(p) for p in (normal,optimized)}
    d=json.loads(normal.read_text());result=audit(d)
    result['normal_optimized_canonical_equal']=d==json.loads(optimized.read_text())
    if not result['normal_optimized_canonical_equal']:raise ValueError('NORMAL_OPTIMIZED_DIFFER')
    result['independence_scope']='Report algebra only: no model-source numerical function is imported, and source-cell values are accepted solely as report data.'
    result['inputs_before']=pins;result['negative_report_data_controls']=controls(d)
    result['inputs_after']={p.name:digest(p) for p in (normal,optimized)}
    if result['inputs_after']!=pins:raise ValueError('REPORT_INPUTS_CHANGED')
    if not all(v['rejected'] for v in result['negative_report_data_controls'].values()):raise ValueError('NEGATIVE_REPORT_CONTROL_ACCEPTED')
    result['source_sha256']=digest(__file__)
    args.output.write_text(json.dumps(enc(result),indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'checks_count':result['checks_count'],'h_display_only':result['h_display_only'],'canonical_equal':result['normal_optimized_canonical_equal']}))
if __name__=='__main__':main()
