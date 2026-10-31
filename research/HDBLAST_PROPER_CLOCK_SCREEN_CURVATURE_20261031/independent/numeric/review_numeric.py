#!/usr/bin/env python3
"""Independent stdlib-only proper-clock numeric aggregation.

The complete actual first-jet oracle remains an explicitly inherited input.
This checker imports no scientific producer and reconstructs every operation
after those authenticated raw boxes using Fraction arithmetic.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
GRID=1<<112
ELL=F(1,10**6)
TSTAR=F(5,4)
RSTAR=F(1,8)
CAPS={
 'HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031':'fb3046fa483cd781060e167324c8002a241898f82627172bcde8187836a3234e',
 'HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031':'e05bf3bb3ca31fe19321568886da29e917b4f1b784239daf69e05c733b71a672',
 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016':'3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
 'HDBLAST_FULL_ANGULAR_RESPONSE_20261026':'a35f84119e84bf2ea8a712ded86ff79977d3ee91ac25f91ccbd202d379c859ac',
 'HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026':'204e647b9f68f2af169b052deabee7c5d7f4bd71195fae31cefbfc2d10265777',
 'HDBLAST_EXACT_PACKET_INITIAL_20261031':'8201452fa5bc745a9ccfd64f67f5f2a15c564a3364af619873825c3a7a6015aa',
}
BACKGROUND_REPORT='abd05a8a93fd2f543a112c344ca7f92cf95ce6d56e379cb80ba1e7d9dfb94882'
FAMILY_JOIN_REPORT='a2aae6d7bd1e522e1cf176349c675fc87b281bf0e8296b6b22e904fc3b49eeee'
TT_REPORT='4a5290e9b3907237a85cb17b72e1262561679d9acb0d5548a2832039416de486'
FORMULA_API='97e8c9be237e1ba47df0b494b29c2b5ebea96e5b11b012a56fdfd9f569d445d1'
FORMULA_REPORT='cb56638ff696f724161150c828b5e8228b22a9796404621dacc2d938ea06f966'
CANONICAL_REPORT='1f2f66bd8e0ba58b1d61db03f07b974e75a73340778ac40a46d0d30fbb7d5216'
PRODUCER_SOURCE='bcab444c1c17f7c62104e510c93aad011cfdae718230dab404f003849aa0e811'
PRODUCER_PINMAP='3bf9ae71f20d406c2620af085327c4564281231877f80e218ec0a2677674cca2'
PRODUCER_FREEZE='66d9f8c7b87ddabab2c052ca6dc331a4827350f1dbe94b9a859565c66a338f5c'
ESSENTIAL_COMPONENT_FILES={
 'background':{
  'BACKGROUND_NORMAL.json':(14234,BACKGROUND_REPORT),
  'DERIVATION.md':(6265,'8cc42cce86d6828fb9b4607b67b82867433d6dddc1bf3b894c7bcdc0058e48ac'),
  'background_certificate.py':(11741,'5fc645c25fe823a9b0a7027daa678b09bffaae34586767f15a4097b6515a4fee'),
  'check_background_certificate.py':(1272,'d200d67247cc5c1f3278a6ad467f2da6e1ca7d68b90d898f79050a0b4eda1b97'),
  'STAGE_FREEZE.json':(2019,'3558841c6b7c433eb5ab55b3b85f5de9a16c9f490917f07b72c75376b7b7e8f9'),
  'SOURCE_PINS.json':(6450,'cf57ec1264106c5758eae6865d07915ad08b577854be31c08073262298ab0ac2'),
  'LOCAL_PINS.json':(527,'03745524e3f7d0875b7b881398437f7c408e4f52e8def2dc433c46a9fe82b26f'),
  'check_family_domain_join.py':(5671,'a9cf7f3131e674f4d617d9c5f995f507cb0cbe1ae8e161061fcc74bbaa010b7b'),
  'FAMILY_DOMAIN_SOURCE_PINS.json':(1226,'a7d2193e97180161107fb07cfd1da6c3d4d0714e8d8aae087342f2cdd0c65a1e'),
  'FAMILY_DOMAIN_JOIN_NORMAL.json':(3782,FAMILY_JOIN_REPORT),
 },
 'tt':{
  'TT_NORMAL.json':(13774,TT_REPORT),
  'check_combined_tt.py':(21752,'83a8e7a132ff180c5f3c144442525f4f835fe9191297be75f34a67328b2450c6'),
 },
}


def need(ok,msg):
    if not ok:raise ValueError(msg)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):return json.loads(Path(path).read_bytes())


def safe_file(path):
    path=Path(path)
    need(path.is_file() and not path.is_symlink(),'FILE_TYPE:'+str(path))
    need(all(not p.is_symlink() for p in path.absolute().parents),'ANCESTOR_SYMLINK')
    return path


def auth_tree(root,pin):
    root=Path(root)
    need(root.is_dir() and not root.is_symlink(),'ROOT_TYPE')
    need(sha(safe_file(root/'MANIFEST.json'))==pin,'MANIFEST:'+root.name)
    rows=read(root/'MANIFEST.json')['files'];seen={'MANIFEST.json'}
    if isinstance(rows,dict):rows=[{'path':k,**v}for k,v in rows.items()]
    for row in rows:
        rel=Path(row['path'])
        need(not rel.is_absolute() and '..'not in rel.parts and str(rel)not in seen,'MANIFEST_MEMBER_PATH')
        p=safe_file(root/rel)
        need(p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'MANIFEST_MEMBER_BYTES:'+str(rel))
        seen.add(str(rel))
    for p in root.rglob('*'):
        need(not p.is_symlink(),'TREE_SYMLINK')
        if p.is_file():need(str(p.relative_to(root))in seen,'UNLISTED_MEMBER:'+str(p))
    return {'manifest_sha256':pin,'members_including_manifest':len(seen)}


class I:
    """Independent reference grid intervals; every primitive rounds outward."""
    def __init__(self,lo=0,hi=None):
        lo,hi=F(lo),F(lo if hi is None else hi)
        need(lo<=hi,'INTERVAL_ORDER')
        self.lo=F((lo*GRID).__floor__(),GRID)
        self.hi=F((hi*GRID).__ceil__(),GRID)
    @staticmethod
    def cv(x):return x if isinstance(x,I) else I(x)
    def __add__(self,x):
        x=self.cv(x);return I(self.lo+x.lo,self.hi+x.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,x):return self+-self.cv(x)
    def __rsub__(self,x):return self.cv(x)+-self
    def __mul__(self,x):
        x=self.cv(x);p=[a*b for a in (self.lo,self.hi)for b in (x.lo,x.hi)]
        return I(min(p),max(p))
    __rmul__=__mul__
    def inv(self):
        need(self.lo>0 or self.hi<0,'ZERO_DENOMINATOR')
        return I(1/self.hi,1/self.lo)
    def __truediv__(self,x):return self*self.cv(x).inv()
    def __rtruediv__(self,x):return self.cv(x)*self.inv()
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        if n==0:return I(1)
        if n%2:return self*self**(n-1)
        a=self**(n//2)
        return I(0 if a.lo<=0<=a.hi else min(a.lo*a.lo,a.hi*a.hi),max(a.lo*a.lo,a.hi*a.hi))
    def mag(self):return max(abs(self.lo),abs(self.hi))
    def row(self):return [str(self.lo),str(self.hi)]


def box(row):
    need(isinstance(row,list) and len(row)==2,'RAW_BOX_SHAPE')
    return I(F(row[0]),F(row[1]))


def boxes(obj):return {k:box(v)if isinstance(v,list)and len(v)==2 else F(v) for k,v in obj.items()}


def same(got,want,name):
    if isinstance(got,I):need(got.row()==want,'INTERVAL_RECONSTRUCTION:'+name)
    else:need(got==F(want),'RATIONAL_RECONSTRUCTION:'+name)


def source_terms(m,b,c,s):
    rho,aa,ay=(b[k] for k in ('rho','a','ay'))
    sp=aa*m['B']-(c*c*m['Dt']-s*s*m['At'])/(2*rho)-c*s*m['Dy']-s*s*m['By']
    sq=c*c*(m['Bt']/rho-aa*(m['A']+m['D'])-m['Dy']/2)+c*s*(m['At']/rho-2*aa*m['B'])+s*s*m['Ay']/2
    M=max((aa*s).mag()+2*(aa*c).mag(),2*(aa*c).mag())
    K=max(2*(ay*c*s).mag(),((2*aa*aa+ay)*c*c).mag())
    return sp,sq,max(sp.mag(),sq.mag()),M,K,(aa*s).mag()


def point_terms(m,b):
    rho,R,theta=b['rho'],b['rho']*b['c'],b['theta']
    eta=R*(m['Ct']-theta*m['C'])/(2*rho)-R*m['Ey']/2-m['B']/2
    Ft=((m['D']+4*m['M']-m['L'])*m['Mt']-2*m['M']*m['Lt'])/(4*rho)-m['E']*(m['D']+m['L'])/(2*R)+(m['C']/R-m['My']/2)*(R*(m['Ct']-theta*m['C'])/(2*rho)-R*m['Ey']/2+m['B']/2)
    Fy=((-m['A']+4*m['M']-m['L'])*m['My']-2*m['M']*m['Ly'])/4+m['C']*(m['A']-m['L'])/(2*R)+(m['E']/R-m['Mt']/(2*rho))*eta
    At=-(m['Mtt']+theta*m['Mt'])/(2*rho**2)+m['Et']/(rho*R)-(m['D']+m['L'])/(2*R**2)
    Ay=-(m['Myy']+b['a']*m['My'])/2+m['Cy']/R-theta*m['Mt']/(2*rho**2)+theta*m['E']/(rho*R)+(m['L']-m['A'])/(2*R**2)
    K=-(m['Myt']+theta*m['My'])/(2*rho)+(m['Ct']+theta*m['C'])/(2*rho*R)+m['Ey']/(2*R)-m['B']/(2*R**2)
    return Ft,Fy,{'At':At,'Ay':Ay,'K':K}


def encode(x):
    if isinstance(x,I):return x.row()
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v)for k,v in x.items()}
    if isinstance(x,list):return [encode(v)for v in x]
    return x


def interval_engine_controls():
    need(I(F(1,3)).lo<=F(1,3)<=I(F(1,3)).hi,'ROUND_ONE_THIRD')
    a=I(-2,3);b=I(-5,7)
    same(a*b,['-15','21'],'FOUR_PRODUCT_ENDPOINTS')
    same(a**2,['0','9'],'SQUARE_ZERO_CROSSING')
    same(-a,['-3','2'],'NEGATIVE_ENDPOINT_ORDER')
    need((I(2,3)/I(5,7)).lo<=F(2,7)and(I(2,3)/I(5,7)).hi>=F(3,5),'RECIPROCAL_ROUNDING')
    try:I(-1,1).inv()
    except ValueError:pass
    else:raise ValueError('ZERO_DIVISION_CONTROL')
    return 6


def independent_domain_checks(bg):
    need(bg['schema']=='proper-clock-background-certificate-v1','BACKGROUND_SCHEMA')
    need(bg['boost']['chi_uniform_all_shooting_roots']is True,'ALL_SHOOTING_ROOTS')
    need(F(bg['ell'])==ELL,'ELL')
    need(F(bg['shoot']['t_dagger'])==ELL*TSTAR and F(bg['shoot']['z_dagger'])==ELL*RSTAR,'COMMON_APEX')
    rho_lo=F(bg['actual_background_local']['rho'][0])
    b_lo,b_hi=map(F,bg['actual_background_local']['rho_a'])
    need(rho_lo>129 and 0<b_lo<=b_hi<20,'WHOLE_TUBE_RHO_AND_ACCELERATION')
    need(list(map(F,bg['actual_background_local']['scaled_z']))==[F(0),F(1,4)],'WHOLE_RADIAL_BACKGROUND_BOX')
    need(list(map(F,bg['compact_tube']['enclosing_scaled_x']))==[F(0),F(8,5)],'AUXILIARY_HORIZON')
    tau_cap=F(1,5500);H=ELL*F(8,5);chi_cap=20*H
    # cosh x=sum x^(2n)/(2n)! <=sum x^(2n)=1/(1-x^2).
    need(0<=chi_cap<1,'RAPIDITY_GEOMETRIC_MAJORANT')
    c_cap=1/(1-chi_cap*chi_cap)
    actual_time_upper=tau_cap*c_cap/129
    family_time_upper=ELL*F(3,2)
    tau_target=F(bg['shoot']['tau_dagger'][1])
    need(0<tau_target<tau_cap<F(bg['compact_tube']['proper_tau'][1]),'TARGET_IN_RESTRICTED_CLOCK_SLAB')
    need(actual_time_upper<family_time_upper,'ACTUAL_FAMILY_COORDINATE_TIME_MARGIN')
    zi_lo,zi_hi=map(F,bg['shoot']['scaled_z_i'])
    need(F(1249,10000)<=zi_lo<=zi_hi<=RSTAR<F(1251,10000),'SOURCE_TUBE_LAUNCH')
    need(20*(ELL*TSTAR)**2/(2*ELL)<F(1,10000),'NARROW_RADIAL_SHIFT')
    need(all(F(bg['gaussian'][k])>0 for k in ('angular_gram_relative_to_kappa_lower','coordinate_determinant_dt_dz_over_dtau_dzi_lower','radial_gram_lower')),'INVERTIBLE_GAUSSIAN_MARGIN')
    return {'restricted_tau_upper':tau_cap,'family_time_upper':family_time_upper,'restricted_coordinate_time_upper':actual_time_upper,'strict_time_margin':family_time_upper-actual_time_upper,'whole_tube_cosh_upper':c_cap,'target_clock_headroom':tau_cap-tau_target,'initial_variation_identity':'zero on initial open unperturbed collar; inherited local completion identity, not inferred from overlapping interval boxes'}


def authenticate_components(repo,bgroot,ttroot,formularoot):
    old={name:auth_tree(repo/'research'/name,pin)for name,pin in CAPS.items()}
    new={}
    for label,root in (('background',bgroot),('tt',ttroot)):
        rows={}
        for rel,(size,pin)in ESSENTIAL_COMPONENT_FILES[label].items():
            p=root/rel
            # Public reports are flat. The private original TT capsule used
            # final/; accepting those same frozen bytes preserves replay there.
            if label=='tt'and rel=='TT_NORMAL.json'and not p.exists():p=root/'final'/rel
            need(sha(safe_file(p))==pin and p.stat().st_size==size,'ESSENTIAL_COMPONENT_FILE:'+label+'/'+rel)
            rows[rel]={'bytes':size,'sha256':pin}
        new[label]={'authentication':'explicit essential public files','files':rows}
    tt_report=ttroot/'TT_NORMAL.json'
    if not tt_report.exists():tt_report=ttroot/'final/TT_NORMAL.json'
    for path,pin,name in ((bgroot/'BACKGROUND_NORMAL.json',BACKGROUND_REPORT,'BACKGROUND'),
                          (bgroot/'FAMILY_DOMAIN_JOIN_NORMAL.json',FAMILY_JOIN_REPORT,'FAMILY_JOIN'),
                          (tt_report,TT_REPORT,'TT'),
                          (formularoot/'proper_clock_coefficient.py',FORMULA_API,'FORMULA_API'),
                          (formularoot/'FINAL_FORMULA_NORMAL.json',FORMULA_REPORT,'FORMULA_REPORT')):
        need(sha(safe_file(path))==pin,'COMPONENT_PIN:'+name)
    new['geometry']={'authentication':'explicit essential public files','files':{
        'proper_clock_coefficient.py':{'bytes':(formularoot/'proper_clock_coefficient.py').stat().st_size,'sha256':FORMULA_API},
        'FINAL_FORMULA_NORMAL.json':{'bytes':(formularoot/'FINAL_FORMULA_NORMAL.json').stat().st_size,'sha256':FORMULA_REPORT}}}
    bg=read(bgroot/'BACKGROUND_NORMAL.json');tt=read(tt_report)
    join=read(bgroot/'FAMILY_DOMAIN_JOIN_NORMAL.json')
    need(join['status']=='PASS_SELECTED_EXACT_FAMILY_PROPER_CLOCK_DOMAIN_JOIN','JOIN_STATUS')
    need(join['frozen_background_sha256']==BACKGROUND_REPORT,'JOIN_BACKGROUND_BINDING')
    need(join['initial_observer_data']['exact_canonical_data_equal_background_near_all_launch_labels']is True and F(join['initial_observer_data']['first_position_system_W0'])==0,'ZERO_INITIAL_OBSERVER_VARIATION')
    for rel,pin in bg['input_sha256'].items():
        # The background report uses keys relative to the repository.
        need(sha(safe_file(repo/rel))==pin,'BACKGROUND_INPUT_PIN:'+rel)
    for rel,pin in join['source_sha256'].items():need(sha(safe_file(repo/rel))==pin,'JOIN_SOURCE_PIN:'+rel)
    need(tt['all_inputs_authenticated_before_and_after']is True and tt['new_source_oracle_calls']==0,'TT_AUTHENTICATED_INHERITED_BOX_REPORT')
    need(tt['enclosure']['tt_term_sign']=='negative','TT_SIGN')
    return old,new,bg,tt,join


def verify_aggregation(r,bg,tt,freeze):
    need(r['schema']=='HDBLAST_FULL_PROPER_CLOCK_COEFFICIENT_V1','REPORT_SCHEMA')
    need(r['status']=='PASS_CONDITIONAL_FULL_PROPER_CLOCK_SCREEN_COEFFICIENT','REPORT_STATUS')
    need(r['point']=={'scaled_x':'5/4','scaled_r':'1/8'}and F(r['ell'])==ELL,'REPORT_EVENT_AND_UNITS')
    need(r['all_inputs_authenticated_before_and_after']is True,'PRODUCER_POSTCHECK')
    need(r['initial_variation_exact_zero_from_open_collar']is True,'ZERO_INITIAL_DATA')
    need(r['frozen_campaign']==freeze,'FROZEN_CAMPAIGN_REPORT')
    need(freeze['time_cells']==64 and list(map(F,freeze['scaled_x']))==[F(0),TSTAR]
         and list(map(F,freeze['scaled_r']))==[F(1249,10000),F(1251,10000)]
         and freeze['no_adaptive_refinement']is True,'FROZEN_CAMPAIGN')
    same(F(r['restricted_tau_upper']),'1/5500','RESTRICTED_CLOCK')
    same(F(r['family_time_upper']),str(ELL*F(3,2)),'FAMILY_TIME')
    same(F(r['tube_rho_lower']),'129','TUBE_RHO_LOWER')
    T=F(bg['shoot']['tau_dagger'][1]);same(T,r['proper_T_upper'],'T')
    c=I(1,F(bg['boost']['cosh_chi_dagger'][1]));s=I(F(bg['boost']['sinh_chi_dagger'][0]),0)
    same(c,r['common_cosh_cap_to_target'],'WHOLE_PATH_COSH')
    same(s,r['common_sinh_cap_to_target'],'WHOLE_PATH_SINH')
    tube_chi=max(abs(F(v))for v in bg['compact_tube']['rapidity'])
    tube_cosh=1/(1-tube_chi**2/2)
    tube_t=F(1,5500)*tube_cosh/129
    for k,v in [('tube_chi_abs_upper',tube_chi),('tube_cosh_upper',tube_cosh),
                ('restricted_background_t_upper',tube_t),('family_time_margin',ELL*F(3,2)-tube_t),
                ('proper_clock_headroom',F(1,5500)-T)]:same(v,r[k],k)
    need(T<F(1,5500)and tube_t<ELL*F(3,2),'FAMILY_DOMAIN')
    rows=r['cells'];need(len(rows)==64,'SOURCE_CELL_COUNT')
    integral=F(0);M=K=acap=F(0)
    for j,row in enumerate(rows):
        xl,xh=TSTAR*j/64,TSTAR*(j+1)/64
        need(list(map(F,row['x']))==[xl,xh],'SOURCE_CELL_TIME_COVERAGE:'+str(j))
        need(row['r']==I(F(1249,10000),F(1251,10000)).row(),'SOURCE_CELL_RADIAL_COVERAGE:'+str(j))
        need(set(row['metric_m1'])=={'A','B','D','At','Ay','Bt','By','Dt','Dy'},'SOURCE_RAW_METRIC_KEYS')
        need(set(row['background'])=={'rho','a','ay'},'SOURCE_RAW_BACKGROUND_KEYS')
        m,b=boxes(row['metric_m1']),boxes(row['background'])
        need(b['rho'].lo>0,'SOURCE_POSITIVE_RHO')
        sp,sq,S,mm,kk,ass=source_terms(m,b,c,s)
        for key,v in [('SP',sp),('SQ',sq),('source_cap',S),('M_cell',mm),('K_cell',kk)]:same(v,row[key],key+':'+str(j))
        weight=ELL*(xh-xl)*b['rho'].hi;same(weight,row['proper_length_upper'],'PROPER_LENGTH:'+str(j))
        integral+=weight*S;M=max(M,mm);K=max(K,kk);acap=max(acap,ass)
    for k,v in [('source_integral',integral),('M',M),('K',K),('as',acap)]:same(v,r[k],k)
    den=1-M*T-K*T*T/2;denx=1-acap*T
    need(den>0 and denx>0,'TRIANGULAR_DENOMINATORS')
    W=integral/den;Y=T*W;X=Y/denx
    for k,v in [('den',den),('den_X',denx),('velocity_cap_W',W),('position_X_cap',X),('position_Y_cap',Y)]:same(v,r[k],k)
    m,b=boxes(r['apex_metric_m2']),boxes(r['apex_background'])
    need(len(m)==42 and set(b)=={'rho','c','theta','a','ay','p','py','Vp','Vpp','Vppp','ell'},'APEX_COMPLETE_TWO_JET')
    need(b['ell']==ELL,'APEX_ELL')
    Ft,Fy,w=point_terms(m,b)
    same(Ft,r['Ft'],'APEX_FT');same(Fy,r['Fy'],'APEX_FY')
    need(set(r['point_weyl_invariants'])==set(w),'WEYL_INVARIANT_KEYS')
    for k,v in w.items():same(v,r['point_weyl_invariants'][k],'WEYL_'+k)
    cd=box(bg['boost']['cosh_chi_dagger']);sd=box(bg['boost']['sinh_chi_dagger'])
    same(cd,r['background_cosh_at_target'],'APEX_COSH');same(sd,r['background_sinh_at_target'],'APEX_SINH')
    R=b['rho']*b['c']
    quadratic=4*(cd*Ft+sd*Fy)/R
    tilt=-4*(I(-X,X)*(cd*w['At']+sd*w['K'])+I(-Y,Y)*(sd*w['Ay']+cd*w['K']))/R
    linear=box(tt['enclosure']['metric_linear_second_TT_coefficient'])
    total=linear+quadratic+tilt
    for k,v in [('boosted_quadratic',quadratic),('full_geodesic_tilt_interval',tilt),('metric_linear_TT',linear),('full_coefficient',total),('opposite_orientation_coefficient',-total)]:same(v,r[k],k)
    need(total.hi<0,'FULL_SIGN_NOT_NEGATIVE')
    scale=10**12;di=[F((total.lo*scale).__floor__(),scale),F((total.hi*scale).__ceil__(),scale)]
    need(list(map(F,r['outward_decimal']))==di,'OUTWARD_DECIMAL_EXACT')
    same(di[1],r['qualitative_sign_upper'],'QUALITATIVE_UPPER')
    same(-di[1]/4,r['ineffective_quadratic_sign_constant'],'INEFFECTIVE_C2_SIGN_CONSTANT')
    need(r['actual_oracle_calls']==65,'ORACLE_COUNT')
    return encode({'source_integral':integral,'M':M,'K':K,'as':acap,'den':den,'den_X':denx,
                   'velocity_cap_W':W,'position_X_cap':X,'position_Y_cap':Y,'point_weyl_invariants':w,
                   'Ft':Ft,'Fy':Fy,'boosted_quadratic':quadratic,'full_geodesic_tilt_interval':tilt,
                   'metric_linear_TT':linear,'full_coefficient':total,'outward_decimal':r['outward_decimal']})


def authenticate_producer(root,report):
    pins={'check_proper_clock.py':PRODUCER_SOURCE,'SOURCE_PINS.json':PRODUCER_PINMAP,'FREEZE.json':PRODUCER_FREEZE}
    for rel,pin in pins.items():need(sha(safe_file(root/rel))==pin,'FROZEN_PRODUCER_FILE:'+rel)
    declared=read(root/'SOURCE_PINS.json')
    expected={'background_report':BACKGROUND_REPORT,'formula_api':FORMULA_API,'freeze':PRODUCER_FREEZE,'tt_report':TT_REPORT}
    need(declared['schema']=='PROPER_CLOCK_NUMERIC_PINS_V1'and declared['old_capsules']==CAPS and declared['new_inputs']==expected,'PRODUCER_INPUT_PINMAP')
    need(report['source_pins_sha256']==PRODUCER_PINMAP and report['producer_source_sha256']==PRODUCER_SOURCE and report['new_input_sha256']==expected,'REPORT_INPUT_BINDING')
    return pins,read(root/'FREEZE.json')


def bind_raw_oracle_inputs(candidate,canonical):
    need(len(candidate['cells'])==len(canonical['cells']),'CANDIDATE_CELL_COUNT')
    for i,(row,base)in enumerate(zip(candidate['cells'],canonical['cells'])):
        need(row['metric_m1']==base['metric_m1']and row['background']==base['background'],'CANONICAL_RAW_ORACLE_CELL:'+str(i))
    need(candidate['apex_metric_m2']==canonical['apex_metric_m2']and candidate['apex_background']==canonical['apex_background'],'CANONICAL_RAW_ORACLE_APEX')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--producer-root',type=Path,required=True)
    ap.add_argument('--report',type=Path,required=True)
    ap.add_argument('--background-root',type=Path,required=True)
    ap.add_argument('--tt-root',type=Path,required=True)
    ap.add_argument('--formula-root',type=Path,required=True)
    ap.add_argument('--candidate-report',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(),'FRESH_OUTPUT')
    need(all(not p.is_symlink()for p in args.output.absolute().parents),'OUTPUT_ANCESTOR_SYMLINK')
    for root in (args.repo_root,args.producer_root,args.background_root,args.tt_root,args.formula_root):
        need(root.resolve()!=args.output.resolve()and root.resolve()not in args.output.resolve().parents,'OUTPUT_OUTSIDE_INPUT_SOURCE')
    own=sha(__file__)
    need(sha(safe_file(args.report))==CANONICAL_REPORT,'CANONICAL_REPORT_PIN')
    canonical=read(args.report)
    producer_pins,freeze=authenticate_producer(args.producer_root,canonical)
    old,new,bg,tt,join=authenticate_components(args.repo_root,args.background_root,args.tt_root,args.formula_root)
    domain=independent_domain_checks(bg)
    controls=interval_engine_controls()
    r=read(safe_file(args.candidate_report))if args.candidate_report else canonical
    bind_raw_oracle_inputs(r,canonical)
    authenticate_producer(args.producer_root,r)
    need(r['source_authentication']==old,'REPORTED_OLD_SOURCE_AUTHENTICATION')
    aggregate=verify_aggregation(r,bg,tt,freeze)
    old2,new2,_,_,_=authenticate_components(args.repo_root,args.background_root,args.tt_root,args.formula_root)
    authenticate_producer(args.producer_root,canonical)
    need(old==old2 and new==new2 and sha(__file__)==own and sha(args.report)==CANONICAL_REPORT,'POST_AUTHENTICATION')
    result=encode({'status':'PASS_INDEPENDENT_FULL_PROPER_CLOCK_NUMERIC_RECONSTRUCTION',
                   'canonical_report_sha256':CANONICAL_REPORT,'checker_sha256':own,
                   'producer_file_pins':producer_pins,'old_source_authentication':old,'new_component_authentication':new,
                   'all_inputs_authenticated_before_and_after':True,'production_scientific_imports':False,
                   'shared_actual_jet_oracle_reproved':False,'source_cells_reconstructed':64,'raw_apex_metric_two_jet_slots':42,
                   'interval_engine_nonzero_controls':controls,'independent_restricted_domain':domain,
                   'reconstructed':aggregate,'candidate_report_sha256':sha(args.candidate_report)if args.candidate_report else None,
                   'scope':'Exact independent aggregation of authenticated inherited first-jet boxes and TT enclosure; retained conditional exact-family/regularity premise and ineffective small-parameter radius.'})
    with args.output.open('x')as f:json.dump(result,f,sort_keys=True,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'outward_decimal':aggregate['outward_decimal'],'output_sha256':sha(args.output)}))


if __name__=='__main__':main()
