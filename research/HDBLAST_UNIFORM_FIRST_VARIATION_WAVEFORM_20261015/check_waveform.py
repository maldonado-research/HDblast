"""Meaningful finite interval checks; distinct from the analytic PDE argument."""
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import tempfile

ROOT=Path(__file__).resolve().parent


def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/(name+'.py'))
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def independent_exp_negative(x):
    # Positive Taylor bounds for exp(x), then reciprocate. No production import.
    x=F(x);n=512+2*(x.numerator//x.denominator)
    term=total=F(1)
    for k in range(1,n+1):
        term*=x/k;total+=term
    tail=term*x/(n+1)/(1-x/(n+2))
    return 1/(total+tail),1/total


def run():
    iv=module('waveform_intervals')
    checks=[]
    def check(name,condition):
        if not condition:
            raise RuntimeError(name)
        checks.append(name)
    def contains(bounds,lo,hi=None):
        return F(bounds[0],iv.SCALE)<=lo<=F(bounds[1],iv.SCALE) and (
            hi is None or F(bounds[0],iv.SCALE)<=hi<=F(bounds[1],iv.SCALE))
    check('negative_ceiling',iv.ceil_div(-5,4)==-1)
    for x in [F(-1,3),F(0),F(1,7),F(9,5)]:
        check('rational_outward_'+str(x),contains(iv.interval(x),x))
    a,b=iv.interval(F(-3,7),F(2,5)),iv.interval(F(-5,9),F(1,3))
    check('signed_interval_product',all(contains(iv.mul(a,b),x*y) for x in [F(-3,7),F(2,5)] for y in [F(-5,9),F(1,3)]))
    check('square_through_zero',contains(iv.square(a),0,F(9,49)))
    for x in [F(0),F(1,8),F(1),F(4),F(16),F(64),F(128),F(256)]:
        expected=independent_exp_negative(x)
        check('exp_bracket_'+str(x),contains(iv.exp_negative(x),*expected))
    check('exp_zero_exact',iv.exp_negative(0)==(iv.SCALE,iv.SCALE))
    check('exp_tail_nonzero_upper',iv.exp_negative(256)[0]==0<iv.exp_negative(256)[1])
    q,d=iv.profile_cell(F(3,2),F(3,2))
    check('bump_peak_value',q==(iv.SCALE,iv.SCALE) and d==(0,0))
    q1,d1=iv.profile_cell(F(5,4),F(5,4))
    q2,d2=iv.profile_cell(F(7,4),F(7,4))
    check('bump_symmetry',q1==q2 and d1==iv.neg(d2))
    check('rising_derivative_positive',d1[0]>0)
    check('falling_derivative_negative',d2[1]<0)
    check('flat_support_endpoints',iv.profile_cell(1,1)==((0,0),(0,0)) and iv.profile_cell(2,2)==((0,0),(0,0)))
    early=iv.profile_cell(F(21,20),F(21,20))
    check('nonzero_early_interior',early[0][0]>0 and early[1][0]>0)
    edge=iv.profile_cell(1,F(4097,4096))
    check('flat_edge_envelope',edge[0][0]==0 and 0<edge[0][1]<iv.SCALE//10**12 and edge[1][0]==0<edge[1][1])
    # Verify containment at rational subpoints of cells, including support splits.
    for index,(left,right) in enumerate([(F(1),F(9,8)),(F(5,4),F(4,3)),(F(7,5),F(8,5)),(F(15,8),F(2)),(F(9,10),F(11,10))]):
        qrange,drange=iv.profile_cell(left,right)
        valid=True
        for k in range(9):
            point=left+(right-left)*k/8
            qp,dp=iv.profile_cell(point,point)
            valid=valid and qrange[0]<=qp[0]<=qp[1]<=qrange[1] and drange[0]<=dp[0]<=dp[1]<=drange[1]
        check('continuous_cell_subpoints_'+str(index),valid)
    c=iv.sharp_waveform_error();ell,T=c['ell'],c['T'];r,s,det=c['r'],c['s'],c['determinant'];A,V=c['A_upper'],c['V_upper'];m=F(1,1000)
    check('majorant_coupled_rows',(1-r)*A-T*V==1+T and (1-s)*V-r*A/m==(1+T)/m)
    check('uniform_output_sources',c['eta_upper']==ell*(10**6*T*A+(1+250*T)*V))
    check('tail_output_sources',c['eta_tail_upper']==ell*(10**6*T*A+ell/m**2+ell/m+r*A/m+s*V+250*T*V))
    check('published_global_error',c['eta_upper']<c['published_eta'])
    check('published_tail_error',0<c['eta_tail_upper']<c['published_tail'])
    check('support_aware_norm_error',c['eta_upper']**2+c['eta_tail_upper']**2<c['norm_error_upper']**2)
    for bad in [(-1,1),(1,4),(2,1)]:
        try:iv.profile_cell(*bad)
        except ValueError:rejected=True
        else:rejected=False
        check('invalid_cell_'+str(bad),rejected)
    generator=module('compute_waveform')
    with tempfile.TemporaryDirectory(prefix='hdblast-waveform-check-') as temp:
        out=Path(temp)/'fresh'
        result=generator.calculate(out,cells=32)
        check('small_generator_status',result['status']=='PASS')
        check('physical_gain_negative',result['physical_gain_interval'][0]<result['physical_gain_interval'][1]<0)
        check('coefficient_scaling',result['quadrature']['C_interval'][0]==result['physical_norm_prefactor_interval'][0]*result['quadrature']['actual_H_integral_interval'][0]/ell)
        check('continuous_cell_count',result['cell_data']['rows']==34)
        check('cell_export_hash',sha256((out/'WAVEFORM_CELLS.csv').read_bytes()).hexdigest()==result['cell_data']['sha256'])
        check('physical_sign_witnesses',result['sign_witnesses']['rising']['ell_S_interval'][1]<0<result['sign_witnesses']['falling']['ell_S_interval'][0])
    return {'status':'PASS','finite_check_groups':len(checks),'checks':checks,
            'scope':'Finite interval/arithmetic and continuous-export checks; not formal proof of the retained analytic PDE or exact-family hypotheses.'}


if __name__=='__main__':
    try:
        print(json.dumps(run(),indent=2))
    except Exception as error:
        print(json.dumps({'status':'FAIL','failure':str(error),'exception':type(error).__name__},indent=2))
        raise SystemExit(1)
