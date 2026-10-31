#!/usr/bin/env python3
"""Frozen bounded postcontact actual TT screen coefficient; stdlib only."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
T,R=Q(5,4),Q(1,8)
CAPS={
 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016':'3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
 'HDBLAST_FULL_ANGULAR_RESPONSE_20261026':'a35f84119e84bf2ea8a712ded86ff79977d3ee91ac25f91ccbd202d379c859ac',
 'HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026':'204e647b9f68f2af169b052deabee7c5d7f4bd71195fae31cefbfc2d10265777',
 'HDBLAST_EXACT_PACKET_INITIAL_20261031':'8201452fa5bc745a9ccfd64f67f5f2a15c564a3364af619873825c3a7a6015aa'}
LOCAL={
 'FREEZE.json':'60d2a032c3fe64caa5018c44d8ee693fbb0e57d58254d1fe05b219684c5a4241',
 'STRUCTURE_NORMAL.json':'0e0152ecc7c18d6c895baaab83e42d11b5b6d564609d2398bd08deb4ca3a30a0',
 'check_tt_block.py':'b3ee4150d134b5ab2391a102755a569a19877f0457a42e397c38fc796cb71ac3',
 'SOURCE_PINS.json':'92d315be572092b9c4e39ac2c6ed7ce2cccff81edeecb0a03bb3ff5e4ebee7e2',
 'REFLECTED_TT_PROOF.md':'b2c3f8266e2440ae09bdc70f4f54d2778fb1f21f947885726df5762a299e67dd'}

def need(ok,message):
    if not ok:raise ValueError(message)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_bytes())
def encode(x):
    if hasattr(x,'lo') and hasattr(x,'hi'):return [str(x.lo),str(x.hi)]
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def authenticate(repo):
    result={}
    for label,pin in CAPS.items():
        root=repo/'research'/label
        need(root.is_dir() and not root.is_symlink(),'ROOT_TYPE:'+label)
        need(sha(root/'MANIFEST.json')==pin,'MANIFEST:'+label)
        seen={'MANIFEST.json'}
        for row in read(root/'MANIFEST.json')['files']:
            rel=Path(row['path'])
            need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'SAFE_MEMBER:'+label)
            seen.add(str(rel));path=root/rel
            need(path.is_file() and not path.is_symlink() and path.stat().st_size==row['bytes'] and sha(path)==row['sha256'],'MEMBER_BYTES:'+label+'/'+str(rel))
        for path in root.rglob('*'):
            need(not path.is_symlink(),'TREE_SYMLINK:'+label)
            if path.is_file():need(str(path.relative_to(root)) in seen,'UNLISTED_MEMBER:'+label)
        result[label]={'manifest_sha256':pin,'members_including_manifest':len(seen)}
    for name,pin in LOCAL.items():need(sha(HERE/name)==pin,'LOCAL_PIN:'+name)
    return result

def source_role_algebra(maps):
    rows=maps['maps'];tt=rows['metric_forcing']['ab']
    need(len(tt)==7 and tt[5]==tt[6] and len(tt[5])==78,'IDENTICAL_COMPLETE_78_TERM_TT_ROWS')
    def normalized(terms,shift=None,factor=Q(1)):
        out={}
        for row in terms:
            powers=dict(row['powers'])
            for k,v in (shift or {}).items():powers[k]=powers.get(k,0)+v
            key=tuple(sorted((k,v) for k,v in powers.items() if v))
            out[key]=out.get(key,Q(0))+factor*Q(row['coefficient'])
        return {k:v for k,v in out.items() if v}
    need(normalized(tt[5])==normalized(rows['QEred']['ab'][5],{'ell':2,'rho':2},Q(-2)),'SCALED_SOURCE_MINUS_TWO_ELL2_RHO2')

class Context:
    def __init__(self,repo):
        self.repo=repo;self.before=authenticate(repo);self.ownpin=sha(__file__)
        self.g=repo/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
        self.a=repo/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
        self.s=repo/'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026'
        sys.path.insert(0,str(self.g/'numeric'));sys.path.insert(0,str(self.g/'projection'))
        self.jets=load('postcontact_actual_jets',self.g/'numeric/actual_first_jets.py')
        self.coeff=load('postcontact_actual_source_coefficients',self.g/'projection/source_coefficients.py')
        self.shell=load('postcontact_shell_sources',self.g/'numeric/shell_source_values.py')
        self.screen=load('postcontact_screen_coefficient',self.s/'theory/screen_coefficient.py')
        self.dy=load('postcontact_dyadic',self.g/'numeric/fast_dyadic.py')
        self.op=load('postcontact_linear_operator',self.g/'projection/inputs/linear_operator.py')
        self.alg=load('postcontact_exact_laurent',self.a/'structure/exact_angular_structure.py')
        self.adapter=load('postcontact_source_adapters',self.g/'projection/inputs/source_adapters.py')
        self.maps=read(self.g/'projection/GEOMETRIC_COEFFICIENT_MAPS.json')
        source_role_algebra(self.maps)
        self.terms=self.maps['maps']['metric_forcing']['ab'][5]
        self.I=self.jets.I;self.ell=self.jets.ELL;self.data,self.errors,_=self.jets.inputs()
        need(self.ell==Q(1,10**6) and self.jets.core.GRID==1<<112,'ARITHMETIC_UNITS')
        self.calls=0;self.source_calls=0;self.bulk_audit_done=False
        self.exact_joins()
    def exact_joins(self):
        symbols=self.alg.symbols('ell rho a h c Ht p py b d V V1 V2 sigma1 sigma2')
        matrices=self.op.matrices()
        for name,m in matrices.items():
            need(all(m[i][j]==0 and m[j][i]==0 for i in (21,22) for j in range(23) if j!=i),'TT_PRIMITIVE_CLOSURE:'+name)
            need(m[21][21]==m[22][22],'TT_ORIENTATIONS:'+name)
        m=self.op.assemble(**symbols);z=symbols
        expected={'A':-3*z['ell']*z['Ht'],'B':-3*z['ell']*z['rho']*z['a'],
          'D':z['ell']*z['ell']*(-6/(z['c']*z['c'])+z['rho']*z['rho']*(-2*z['a']*z['a']+2*z['h']*z['h']-2*z['b']-8*z['d']+Q(4,3)*z['V']))}
        for name,v in expected.items():need(m[name][21][21]==m[name][22][22]==v,'TT_EXACT_OPERATOR:'+name)
        need(all(not m['M'][i][j] and not m['M'][j][i] for i in (21,22) for j in range(23)),'TT_ZERO_M')
        need(m['PN'][21][21]==m['PN'][22][22]==1 and not m['PD'][21][21] and not m['PD'][22][22],'NEUMANN_PROJECTORS')
        sh=self.alg.symbols('a p sigmap sigmapp sigmappp h X V L d v')
        boundary=self.op.boundary_forcing(z['ell'],z['rho'],self.adapter.shell_sources(**sh))
        g=-Q(8,3)*z['ell']*z['rho']*(sh['a']*sh['h']+sh['sigmap']*sh['v']/3)*sh['d']
        need(boundary[21]==boundary[22]==g and bool(g),'NONZERO_NEUMANN_PHYSICAL_SOURCE')
        initial=read(HERE/'STRUCTURE_NORMAL.json')
        need(initial['completed_registered_old_new_TT_initial_functions_zero'] is True,'EXACT_REGISTERED_TT_INITIAL_ZERO')
        self.joins={'TT_scalar_operator_exact':True,'TT_rows_columns_close':True,'both_orientations_identical':True,'TT_M_zero':True,'Neumann_source_generally_nonzero':True,'bulk_scaled_source_identity':'metric_forcing_T=-2 ell^2 rho^2 QEred_T','bulk_terms':78,'registered_initial_zero_report_sha256':LOCAL['STRUCTURE_NORMAL.json']}
    def actual(self,x,r):
        self.calls+=1
        need(self.calls<=140,'FROZEN_ORACLE_BUDGET')
        return self.jets.actual_jet(x,r,mode='actual')
    def forcing(self,x,r):
        m,s,b=self.actual(x,r);raw={**m,**s,**b}
        converted={k:self.dy.D(v) if hasattr(v,'lo') else v for k,v in raw.items()}
        value=self.coeff.evaluate_encoded(self.terms,converted);self.source_calls+=1
        if not self.bulk_audit_done:
            reference=self.coeff.evaluate_encoded(self.terms,raw)
            allrows=self.coeff.evaluate_sources(m,s,b,names=('metric_forcing',),path=self.g/'projection/GEOMETRIC_COEFFICIENT_MAPS.json')
            need(value.row()==reference.row()==allrows['metric_forcing']['ab'][5].row()==allrows['metric_forcing']['ab'][6].row(),'DYADIC_FULL_ORIGINAL_SOURCE_IDENTITY')
            self.bulk_audit_done=True
        return self.I(value.lo,value.hi)
    def boundary(self,x):
        m,s,b=self.actual(x,self.I(0));eta0=self.I(*self.data['shell']['boxes']['eta_b']['exact'])
        raw=self.shell.raw_fields(m['A'],m['D'],m['E'],m['L'],m['M'],s['s'],b['a'],b['p'],4*(1+eta0))
        tensor=self.shell.projected(raw)['ab']
        need(tensor[5].row()==tensor[6].row(),'NUMERICAL_SHELL_ORIENTATION')
        return 2*self.ell*b['rho']*tensor[5]
    def operator(self):
        x,r=self.I(0,T),self.I(0,T+R)
        rho,eta,p,a,*_=self.jets.core.background(r,self.data,'none')
        rho,eta,p,a=[v.val() for v in (rho,eta,p,a)]
        c,Ht,_=self.jets.temporal(x);c,Ht=c.val(),Ht.val()
        ay=-1/rho**2-p**2/3;b=-(ay+a**2);d=1/rho**2-a**2
        V=self.jets.core.poly(self.jets.core.U,eta);h=Ht/rho
        A=-3*self.ell*Ht;B=-3*self.ell*rho*a
        D=self.ell**2*(-6/c**2+rho**2*(-2*a**2+2*h**2-2*b-8*d+Q(4,3)*V))
        caps={k:v.mag() for k,v in {'A':A,'B':B,'D':D}.items()}
        sigma=caps['A']+caps['B']+T*caps['D']/2;q=T*sigma
        need(0<=q<1,'FULL_FEEDBACK_CONTRACTION')
        return {'coefficient_intervals':{'A':A,'B':B,'D':D},'caps':caps,'sigma':sigma,'q':q}
    def finish(self):
        need(authenticate(self.repo)==self.before and sha(__file__)==self.ownpin,'POST_AUTHENTICATION')

def campaign(ctx):
    op=ctx.operator();bulk=[];f=Q(0)
    for i in range(8):
        for j in range(8):
            xb=[T*i/8,T*(i+1)/8];rb=[(T+R)*j/8,(T+R)*(j+1)/8]
            value=ctx.forcing(ctx.I(*xb),ctx.I(*rb));f=max(f,value.mag())
            bulk.append({'x':xb,'r':rb,'F_T':value})
        print(json.dumps({'stage':'source_cover','rows_complete':i+1,'rows':8}),flush=True)
    shell=[];g=Q(0)
    for i in range(8):
        xb=[(T-R)*i/8,(T-R)*(i+1)/8];value=ctx.boundary(ctx.I(*xb));g=max(g,value.mag())
        shell.append({'x':xb,'G_T':value})
    endpoint=ctx.boundary(ctx.I(T-R));curve=[];integral=ctx.I(0)
    knots=sorted(set([T*i/32 for i in range(33)]+[T-R]))
    need(len(knots)==34,'DECLARED_CHARACTERISTIC_SPLIT')
    for i,(lo,hi) in enumerate(zip(knots,knots[1:])):
        plus=[R+T-hi,R+T-lo]
        minus=sorted([abs(R-T+lo),abs(R-T+hi)])
        vp=ctx.forcing(ctx.I(lo,hi),ctx.I(*plus));vm=ctx.forcing(ctx.I(lo,hi),ctx.I(*minus))
        integral+=(hi-lo)*(vp+vm)/2
        curve.append({'x':[lo,hi],'positive_r':plus,'reflected_r':minus,'positive_F':vp,'reflected_F':vm,'length':hi-lo})
        if (i+1)%8==0:print(json.dumps({'stage':'characteristic_integrals','segments_complete':i+1,'segments':33}),flush=True)
    need(sum((row['length'] for row in curve),Q(0))==T,'CURVE_COMPLETE_LENGTH')
    seed=integral-endpoint;sigma,q=op['sigma'],op['q']
    norm=(2*g+T*f)/(1-q);tail=q*norm
    ux=seed+ctx.I(-tail,tail)
    m,s,b=ctx.actual(ctx.I(T),ctx.I(R))
    screen=ctx.screen.coefficients(m,b,ux)
    scale=Q(3,2)/(b['ell']*b['rho']**2*b['c'])
    need(scale.lo>0,'POSITIVE_SCREEN_SCALE')
    total=screen['plus'];sign='positive' if total.lo>0 else ('negative' if total.hi<0 else 'inconclusive')
    need(ctx.calls==140 and ctx.source_calls==130,'EXACT_FROZEN_CALL_COUNT')
    return {'status':'CERTIFIED_CONDITIONAL_POSTCONTACT_SCREEN_SIGN' if sign!='inconclusive' else 'VALID_POSTCONTACT_ENCLOSURE_SIGN_INCONCLUSIVE',
       'sign':sign,'point':[T,R],'operator':op,'bulk_cap_cover':bulk,'bulk_source_cap':f,
       'shell_history_cover':shell,'shell_source_cap':g,'shell_endpoint_G':endpoint,
       'characteristic_segments':curve,'zero_bulk_feedback_seed_x':seed,
       'bulk_source_derivative_integral':integral,'full_target_X_upper':norm,
       'complete_feedback_Ux_error':tail,'actual_full_U_Tx_interval':ux,
       'apex_actual_first_m':m,'apex_actual_scalar':s,'apex_background':b,
       'screen_scale':scale,'full_fixed_quadratic_F':screen['F'],
       'screen_quadratic':screen['quadratic'],'screen_linear':screen['linear'],
       'actual_full_screen_h':total,'actual_full_screen_minus':screen['minus'],
       'screen_error_charges':{'feedback_scale_error_upper':scale.hi*tail,'seed_Ux_interval_width':seed.hi-seed.lo,'quadratic_interval_width':screen['quadratic'].hi-screen['quadratic'].lo,'actual_first_errors_included_in_every_source_and_quadratic_interval':True},
       'calls':{'actual_first_jets':ctx.calls,'bulk_source_cells':ctx.source_calls,'shell':9,'apex':1},
       'acceptance':{'all_exact_joins':True,'complete_feedback_q_lt_one':True,'h_interval_excludes_zero':sign!='inconclusive'},
       'limits':['Conditional actual formal target with retained source/regularity/constraint premises.','Specified registered slicing observer and screen; no gauge-invariance claim.','No W1/W2 accuracy comparison or C2 improvement is asserted.','No finite-amplitude nonlinear shell-contact solution or cosmological-origin claim.','Closed upstream sources are authenticated; inherited PDE existence premises are not reproved.']}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    repo=args.repo_root.resolve();output=args.output
    need(not output.exists() and not output.is_symlink(),'FRESH_OUTPUT')
    target=output.resolve();need(repo!=target and repo not in target.parents,'OUTPUT_OUTSIDE_SOURCE_ROOT')
    for parent in output.parents:need(not parent.is_symlink(),'OUTPUT_PARENT_SYMLINK')
    ctx=Context(repo);result=campaign(ctx);ctx.finish()
    result.update({'schema':'HDBLAST_TT_POSTCONTACT_ACTUAL_SCREEN_V1','freeze_sha256':LOCAL['FREEZE.json'],'producer_sha256':ctx.ownpin,'source_authentication':ctx.before,'local_source_pins':LOCAL,'exact_source_role_joins':ctx.joins,'all_sources_authenticated_before_imports_and_after':True})
    with output.open('x') as handle:json.dump(encode(result),handle,sort_keys=True,indent=2);handle.write('\n')
    print(json.dumps({'status':result['status'],'h_display_only':[float(result['actual_full_screen_h'].lo),float(result['actual_full_screen_h'].hi)],'source_calls':ctx.calls,'report_sha256':sha(output)}),flush=True)

if __name__=='__main__':main()
