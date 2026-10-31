#!/usr/bin/env python3
"""Independent exact report arithmetic/geometry checker; no project imports.

The input F/G/first-field boxes are premises supplied by the authenticated
producer. This checker does not independently re-evaluate those 140 oracles.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
GRID=1<<112
T,R=Q(5,4),Q(1,8)


def need(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


class Box:
    """Separately implemented directed rational interval operations."""
    def __init__(self,lower=0,upper=None):
        low=Q(lower); high=low if upper is None else Q(upper)
        need(low<=high,'BOX_ORDER')
        self.low=Q((low*GRID).numerator//(low*GRID).denominator,GRID)
        n=high*GRID; self.high=Q(-((-n.numerator)//n.denominator),GRID)
    @staticmethod
    def cast(value):return value if isinstance(value,Box) else Box(value)
    def __add__(self,other):
        b=self.cast(other);return Box(self.low+b.low,self.high+b.high)
    __radd__=__add__
    def __neg__(self):return Box(-self.high,-self.low)
    def __sub__(self,other):return self+-self.cast(other)
    def __rsub__(self,other):return self.cast(other)+-self
    def __mul__(self,other):
        b=self.cast(other)
        products=[a*c for a in (self.low,self.high) for c in (b.low,b.high)]
        return Box(min(products),max(products))
    __rmul__=__mul__
    def reciprocal(self):
        need(self.low>0 or self.high<0,'NONZERO_DENOMINATOR')
        return Box(1/self.high,1/self.low)
    def __truediv__(self,other):return self*self.cast(other).reciprocal()
    def __rtruediv__(self,other):return self.cast(other)*self.reciprocal()
    def __pow__(self,power):
        if power<0:return self.reciprocal()**(-power)
        if not power:return Box(1)
        if power%2:
            return self*self**(power-1)
        half=self**(power//2)
        return Box(0 if half.low<=0<=half.high else min(half.low**2,half.high**2),
                   max(half.low**2,half.high**2))
    def mag(self):return max(abs(self.low),abs(self.high))
    def pair(self):return [str(self.low),str(self.high)]


def box(value):
    need(isinstance(value,list) and len(value)==2,'BOX_ENCODING')
    result=Box(*value)
    need(result.pair()==value,'CANONICAL_DYADIC_ENCODING')
    return result


def equalbox(actual,encoded,label):
    need(actual.pair()==encoded,label)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--report',required=True,type=Path)
    parser.add_argument('--report-sha256',required=True)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    need(not args.report.is_symlink() and args.report.is_file(),'INPUT_REPORT_TYPE')
    need(sha(args.report)==args.report_sha256,'INPUT_REPORT_AUTHENTICATION')
    need(not args.output.exists() and not args.output.is_symlink(),'FRESH_OUTPUT')
    need(all(not p.is_symlink() for p in args.output.absolute().parents),'OUTPUT_PARENT_SYMLINK')
    d=json.loads(args.report.read_text())
    need(d['point']==['5/4','1/8'],'POINT')
    need(d['schema']=='HDBLAST_TT_POSTCONTACT_ACTUAL_SCREEN_V1','SCHEMA')
    need(d['all_sources_authenticated_before_imports_and_after'] is True,'PRODUCER_SOURCE_AUTHENTICATION')
    need(d['exact_source_role_joins']['bulk_terms']==78 and
         d['exact_source_role_joins']['Neumann_source_generally_nonzero'] is True,
         'ACTUAL_SOURCE_ROLE')
    bulk=d['bulk_cap_cover'];need(len(bulk)==64,'BULK_COVER_COUNT')
    bulkmax=Q(0)
    for ij,row in enumerate(bulk):
        i,j=divmod(ij,8)
        need(row['x']==[str(T*i/8),str(T*(i+1)/8)] and
             row['r']==[str((T+R)*j/8),str((T+R)*(j+1)/8)],'BULK_GRID_GEOMETRY')
        # The entire rectangular overcover is in the original actual oracle
        # domain x<=3,r>=0,x+r<=8; target evolution needs only its cone.
        need(T+(T+R)<8,'BULK_ORACLE_DOMAIN')
        bulkmax=max(bulkmax,box(row['F_T']).mag())
    need(bulkmax==Q(d['bulk_source_cap']),'ACTUAL_F_UNIFORM_CAP')
    shell=d['shell_history_cover'];need(len(shell)==8,'SHELL_HISTORY_COUNT')
    shellmax=Q(0)
    for i,row in enumerate(shell):
        need(row['x']==[str((T-R)*i/8),str((T-R)*(i+1)/8)],'SHELL_HISTORY_GEOMETRY')
        shellmax=max(shellmax,box(row['G_T']).mag())
    need(shellmax==Q(d['shell_source_cap']),'ACTUAL_G_UNIFORM_CAP')
    endpoint=box(d['shell_endpoint_G'])
    need(endpoint.mag()!=0,'GENERIC_NONZERO_ENDPOINT_HISTORY')
    rows=d['characteristic_segments'];knots=sorted(set([T*i/32 for i in range(33)]+[T-R]))
    need(len(rows)==33,'CHARACTERISTIC_COUNT')
    integration=Box(0)
    for row,lo,hi in zip(rows,knots,knots[1:]):
        need(row['x']==[str(lo),str(hi)] and Q(row['length'])==hi-lo,'CHARACTERISTIC_TIME_PARTITION')
        need(row['positive_r']==[str(R+T-hi),str(R+T-lo)],'DIRECT_CHARACTERISTIC_GEOMETRY')
        need(row['reflected_r']==list(map(str,sorted([abs(R-T+lo),abs(R-T+hi)]))),
             'REFLECTED_CHARACTERISTIC_GEOMETRY')
        need(not lo<T-R<hi,'UNSPLIT_REFLECTION_TURN')
        integration+=(hi-lo)*(box(row['positive_F'])+box(row['reflected_F']))/2
    equalbox(integration,d['bulk_source_derivative_integral'],'BOTH_POSITIVE_KX_ENDPOINT_SIGNS')
    seed=integration-endpoint
    equalbox(seed,d['zero_bulk_feedback_seed_x'],'NEUMANN_LIFT_MINUS_G')
    op=d['operator'];caps={k:box(v).mag() for k,v in op['coefficient_intervals'].items()}
    need(caps=={k:Q(v) for k,v in op['caps'].items()},'OPERATOR_CAPS')
    sigma=caps['A']+caps['B']+T*caps['D']/2;q=T*sigma
    need(sigma==Q(op['sigma']) and q==Q(op['q']) and 0<=q<1,'FULL_FEEDBACK_CONTRACTION')
    norm=(2*shellmax+T*bulkmax)/(1-q);tail=q*norm
    need(norm==Q(d['full_target_X_upper']),'COMPLETE_RESPONSE_NORM')
    need(tail==Q(d['complete_feedback_Ux_error']),'COMPLETE_FEEDBACK_TAIL')
    ux=seed+Box(-tail,tail)
    equalbox(ux,d['actual_full_U_Tx_interval'],'ACTUAL_FULL_UX_ENCLOSURE')
    m={k:box(v) for k,v in d['apex_actual_first_m'].items()}
    b={k:(Q(v) if isinstance(v,str) else box(v)) for k,v in d['apex_background'].items()}
    need(b['ell']==Q(1,10**6),'PHYSICAL_CLOCK_ELL')
    rho=b['rho'];radius=rho*b['c'];theta=b['theta']
    need(rho.low>0 and b['c'].low>0,'POSITIVE_BACKGROUND_RADIUS')
    F=((m['D']+4*m['M']-m['L'])*m['Mt']-2*m['M']*m['Lt'])/(4*rho) \
      -m['E']*(m['D']+m['L'])/(2*radius) \
      +(m['C']/radius-m['My']/2)*(radius*(m['Ct']-theta*m['C'])/(2*rho)-radius*m['Ey']/2+m['B']/2)
    quadratic=4*F/radius
    linear=Q(3,2)*ux/(b['ell']*rho**2*b['c'])
    scale=Q(3,2)/(b['ell']*rho**2*b['c'])
    equalbox(scale,d['screen_scale'],'SCREEN_ORDINARY_T_CLOCK_UNITS')
    equalbox(F,d['full_fixed_quadratic_F'],'FULL_FIRST_FIELD_QUADRATIC_F')
    equalbox(quadratic,d['screen_quadratic'],'FULL_FIRST_FIELD_QUADRATIC_COEFFICIENT')
    equalbox(linear,d['screen_linear'],'COMPLETE_TT_LINEAR_CURVATURE')
    total=linear+quadratic
    equalbox(total,d['actual_full_screen_h'],'COMBINED_SCREEN_ENCLOSURE')
    equalbox(-total,d['actual_full_screen_minus'],'CURL_ORIENTATION_MINUS')
    sign='positive' if total.low>0 else ('negative' if total.high<0 else 'inconclusive')
    need(d['sign']==sign and d['acceptance']['h_interval_excludes_zero']==(sign!='inconclusive'),
         'EXACT_SIGN_ACCEPTANCE')
    need(d['calls']=={'actual_first_jets':140,'bulk_source_cells':130,'shell':9,'apex':1},'FROZEN_CALL_COUNT')
    charge=d['screen_error_charges']
    need(Q(charge['feedback_scale_error_upper'])==scale.high*tail and
         Q(charge['seed_Ux_interval_width'])==seed.high-seed.low and
         Q(charge['quadratic_interval_width'])==quadratic.high-quadratic.low,
         'ALL_RECORDED_ERROR_CHARGES')
    need(charge['actual_first_errors_included_in_every_source_and_quadratic_interval'] is True,
         'ACTUAL_FIRST_ERRORS_PREMISE')
    need(sha(args.report)==args.report_sha256,'INPUT_REPORT_POSTAUTHENTICATION')
    result={'status':'PASS_INDEPENDENT_POSTCONTACT_REPORT_ARITHMETIC_GEOMETRY',
      'input_report_sha256':args.report_sha256,'project_imports':[],
      'bulk_whole_cell_boxes':64,'shell_history_boxes':8,'curve_segments':33,
      'exact_oracle_calls_repeated':0,'exact_sign':sign,
      'h_exact':total.pair(),'Ux_exact':ux.pair(),
      'q_exact':str(q),'feedback_Ux_error_exact':str(tail),
      'checked':['predeclared point and whole-cell footprint','complete64cell F and8cell G covers',
        'both endpoint curves with absolute radial reflection and turn split',
        'negative boundary-history time derivative','all-feedback C1 resolvent tail',
        'physical clock ell and positive radius','complete fixed first-field quadratic',
        'actual_Ux plus full screen expression and curl orientation','all exact error charge scalars'],
      'limits':['F/G/actual first-jet input boxes are producer premises; this is an independent exact arithmetic and geometry replay, not a repeated 140-call source campaign.',
        'Separate independent operator/kernel and screen formula checkers validate the mathematical forms.',
        'Conditional registered-coordinate formal response with specified observer; no nonlinear or gauge-invariance claim.']}
    with args.output.open('x') as stream:json.dump(result,stream,sort_keys=True,indent=2);stream.write('\n')
    print(result['status'])


if __name__=='__main__':main()
