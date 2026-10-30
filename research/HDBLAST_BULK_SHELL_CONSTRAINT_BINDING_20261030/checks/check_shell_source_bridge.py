#!/usr/bin/env python3
"""Exact physical-shell source, gauge, normalization and Robin binding."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace

sys.dont_write_bytecode = True
PINS = {
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/source_adapters.py':'0233ff25b9873dfc4060ba0d4963f50d0b198dbe71368625017840d225b19f8a',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py':'a008be0d89e2ac44940c273cfd2cbc5ba0d09de8fc172d5d81ef29121ab60189',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/shell_source_values.py':'ed9c6e084eae776509481fee4f490cea349eb5c1665139e286bd0a58be483d7f',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py':'a68fa8025d2170efca6c3e85d443910a796d0c30a9b7dd17ed988cd9e00bbb00',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json':'34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615',
 'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/boundary-review/BOUNDARY_COMPATIBILITY_AUDIT.md':'581bcadecd44c408484d20632d7d3e524fb0bd098370d1fd24da2744bd6c1cf1',
}
PREFIX='research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/'

def need(ok, reason):
    if not ok: raise ValueError(reason)

class P:
    """Collected rational Laurent polynomials with dynamic generators."""
    def __init__(self, value=0):
        if isinstance(value,P): self.t=dict(value.t)
        elif isinstance(value,dict): self.t={k:F(v) for k,v in value.items() if v}
        else: self.t={} if not value else {():F(value)}
    @classmethod
    def var(cls,name): return cls({((name,1),):F(1)})
    def __add__(self,b):
        if isinstance(b,S): return NotImplemented
        out=dict(self.t)
        for k,v in P(b).t.items():out[k]=out.get(k,F(0))+v
        return P(out)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.t.items()})
    def __sub__(self,b):return self+-P(b)
    def __rsub__(self,b):return P(b)+-self
    def __mul__(self,b):
        if isinstance(b,S): return NotImplemented
        out={}
        for k,v in self.t.items():
            for l,w in P(b).t.items():
                d=dict(k)
                for n,e in l:d[n]=d.get(n,0)+e
                key=tuple(sorted((n,e) for n,e in d.items() if e))
                out[key]=out.get(key,F(0))+v*w
        return P(out)
    __rmul__=__mul__
    def inv(self):
        need(len(self.t)==1,'NONMONOMIAL_DENOMINATOR')
        (k,v),=self.t.items();return P({tuple((n,-e) for n,e in k):1/v})
    def __truediv__(self,b):
        if isinstance(b,S):return NotImplemented
        return self*P(b).inv()
    def __rtruediv__(self,b):return P(b)*self.inv()
    def __pow__(self,n):
        need(isinstance(n,int),'INTEGER_POWER')
        if n<0:return self.inv()**(-n)
        out=P(1)
        for _ in range(n):out=out*self
        return out
    def __bool__(self):return bool(self.t)
    def __eq__(self,b):return self.t==P(b).t
    def encoded(self):return [[[[n,e] for n,e in k],str(v)] for k,v in sorted(self.t.items())]

class S:
    def __init__(self,value=0):
        self.c=tuple(P(v) for v in value) if isinstance(value,(list,tuple)) else (P(value),P(),P())
        need(len(self.c)==3,'SERIES_ORDER')
    @classmethod
    def cv(cls,v):return v if isinstance(v,cls) else cls(v)
    def __add__(self,b):return S([x+y for x,y in zip(self.c,self.cv(b).c)])
    __radd__=__add__
    def __neg__(self):return S([-v for v in self.c])
    def __sub__(self,b):return self+-self.cv(b)
    def __rsub__(self,b):return self.cv(b)+-self
    def __mul__(self,b):
        b=self.cv(b);return S([sum((self.c[j]*b.c[k-j] for j in range(k+1)),P()) for k in range(3)])
    __rmul__=__mul__
    def inv(self):
        a,b,c=self.c;return S([1/a,-b/a**2,b*b/a**3-c/a**2])
    def __truediv__(self,b):return self*self.cv(b).inv()
    def __rtruediv__(self,b):return self.cv(b)*self.inv()
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        out=S(1)
        for _ in range(n):out=out*self
        return out

def authenticate(root):
    for rel,pin in PINS.items():
        p=root/rel
        need(p.is_file() and not p.is_symlink(),'SOURCE_TYPE:'+rel)
        need(hashlib.sha256(p.read_bytes()).hexdigest()==pin,'SOURCE_SHA256:'+rel)
    return dict(PINS)

def load(name,path):
    sp=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod);return mod

def numeric_functions(root,mutation):
    source=(root/PREFIX/'numeric/shell_source_values.py').read_text()
    if mutation=='numeric_sigma2':source=source.replace('sigma2=4*(1+eta0)','sigma2=2*(1+eta0)')
    tree=ast.parse(source)
    allowed={'SCALAR','VECTOR','YYK','YHESS'}
    body=[x for x in tree.body if (isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id in allowed for t in x.targets)) or (isinstance(x,ast.FunctionDef) and x.name in ('raw_fields','projected','actual_shell_sources'))]
    env={'F':F};exec(compile(ast.fix_missing_locations(ast.Module(body=body,type_ignores=[])),'authenticated_numeric_shell_functions','exec'),env)
    return env

Y2=(F(1,6),F(1,2),F(1))
GRAD2=(F(4,3),F(2),F(-4))
HESS2=(F(8),F(4),F(24))
YHESS=(F(-4,9),F(-4,3),F(-8,3),F(-1,4),F(1,6),F(2,3),F(2,3))
VECTOR=(F(1,4),F(1,2))

def project_angular_scalar(expression):
    out=[P(),P(),P()]
    for key,value in P(expression).t.items():
        d=dict(key);ang=(d.pop('Y',0),d.pop('G',0),d.pop('H2',0))
        row={(2,0,0):Y2,(0,1,0):GRAD2,(0,0,1):HESS2}.get(ang)
        need(row is not None,'UNEXPECTED_ANGULAR_SCALAR_MONOMIAL:'+str(ang))
        coeff=P({tuple(sorted(d.items())):value})
        out=[a+c*coeff for a,c in zip(out,row)]
    return tuple(out)

def encoded(terms,values):
    out=P()
    for term in terms:
        v=P(F(term['coefficient']))
        for name,power in term['powers'].items():v=v*values[name]**power
        out=out+v
    return out

CONTROLS={
 'omit_unit_normal_second':'AFFINE_METRIC_00',
 'halve_parameter_variation':'ADAPTER_QK00_0',
 'numeric_time_sign':'NUMERIC_RAW_00',
 'omit_tensor_TT':'ADAPTER_QKab_5',
 'wrong_sigma3':'ADAPTER_Qp_0',
 'numeric_sigma2':'ACTUAL_NUMERIC_ONE_COPY_00_0',
 'omit_quadratic_reference_connection':'GENERIC_QC_ACTUAL_0',
 'positive_time_angular_trace':'GENERIC_QC_ACTUAL_0',
 'omit_anisotropic_trace':'GENERIC_QC_ACTUAL_0',
 'coordinate_instead_of_frame_normal_derivative':'QC_MAP_SHELL_0',
 'positive_QR_trace':'GN_LORENTZ_TRACE_0',
 'wrong_inward_sign':'BOUNDARY_K:0',
 'missing_ell':'BOUNDARY_K:0',
 'wrong_two_copy_factor':'ACTUAL_NUMERIC_TWO_COPY_00_0',
 'off_shell_scalar_term_omitted':'OFFSHELL_SCALAR_CODAZZI',
 'wrong_coordinate_Robin':'COORDINATE_COVECTOR_ROBIN',
 'wrong_frame_Robin':'PARALLEL_FRAME_COVECTOR_ROBIN',
 'equate_lowerings_without_gauge':'QC_LOWERING_DIFFERENCE_0',
}

def audit(root,mutation=None):
    before=authenticate(root)
    ad=load('bridge_adapters',root/PREFIX/'projection/inputs/source_adapters.py')
    op=load('bridge_linear_operator',root/PREFIX/'projection/inputs/linear_operator.py')
    num=numeric_functions(root,mutation)
    g={n:P.var(n) for n in ('h','X','V','L','d','v','a','p','eta','rho','ell','c','theta')}
    h,X,V,L,d,v,a,p,eta,rho,ell,c,theta=[g[n] for n in g]
    sigma1=-2*p;sigma2=4*(1+eta);sigma3=P(4)
    beta=a*h+sigma1*v/3
    a0=a*h*h/4-sigma1*h*v/6-sigma2*v*v/6
    p0=p*h*h/4+sigma2*h*v/2+sigma3*v*v/2
    rows=[]
    def eq(x,y,name):
        remainder=P(x)-P(y)
        need(not remainder,'NONZERO:'+name+':'+json.dumps(remainder.encoded(),separators=(',',':')))
        rows.append(name)
    def equal_rows(x,y,prefix):
        need(len(x)==len(y),'ROW_SHAPE:'+prefix)
        for i,(u,w) in enumerate(zip(x,y)):eq(u,w,prefix+'_'+str(i))
    # Independent affine normalization expansion, with X2=2 times eps² coefficient.
    ep=S([0,1,0]);ninv=S([1,-h/2,3*h*h/8])
    if mutation=='omit_unit_normal_second':ninv=S([1,-h/2,0])
    raw_expected={}
    for key,m,sign in (('00',-2*X,-1),('0a',V,0),('ab_trace',2*L,1),('ab_hessian',2*d,0)):
        my=2*a*m+beta*sign
        metric=S([sign,m,0]);metric_y=S([2*a*sign,my,0])
        tension=S([6*a,sigma1*v,sigma2*v*v/2])
        residual=ninv*metric_y/2-tension*metric/6
        eq(residual.c[0],0,'BACKGROUND_PHYSICAL_'+key)
        eq(residual.c[1],0,'FIRST_PHYSICAL_'+key)
        target={'00':-a0+2*beta*X,'0a':-beta*V,'ab_trace':a0-2*beta*L,'ab_hessian':-2*beta*d}[key]
        eq(2*residual.c[2],target,'AFFINE_METRIC_'+key)
        raw_expected[key]=target
    vy=p*h/2-sigma2*v/2
    pres=ninv*S([p,vy,0])+S([sigma1,sigma2*v,sigma3*v*v/2])/2
    eq(pres.c[0],0,'BACKGROUND_SCALAR_PHYSICAL')
    eq(pres.c[1],0,'FIRST_SCALAR_PHYSICAL')
    eq(2*pres.c[2],p0,'AFFINE_SCALAR')
    raw_expected['scalar']=p0
    physical={'00':tuple(z*raw_expected['00'] for z in Y2),
              '0a':tuple(z*raw_expected['0a'] for z in VECTOR),
              'ab':tuple(z*raw_expected['ab_trace']+w*raw_expected['ab_hessian'] for z,w in zip((*Y2,0,0,0,0),YHESS)),
              'scalar':tuple(z*p0 for z in Y2)}
    source=ad.shell_sources(a=a,p=p,sigmap=sigma1,sigmapp=sigma2,
                           sigmappp=3 if mutation=='wrong_sigma3' else sigma3,h=h,X=X,V=V,L=L,d=d,v=v)
    if mutation=='halve_parameter_variation':source={k:tuple(z/2 for z in val) for k,val in source.items()}
    if mutation=='omit_tensor_TT':source['QKab']=source['QKab'][:5]+(P(),P())
    for target,key in (('00','QK00'),('0a','QK0a'),('ab','QKab'),('scalar','Qp')):
        equal_rows(source[key],physical[target],'ADAPTER_'+key)
    raw=num['raw_fields'](h,2*X if mutation=='numeric_time_sign' else -2*X,V,2*L,2*d,v,a,p,sigma2)
    for key in raw_expected:eq(raw[key],raw_expected[key],'NUMERIC_RAW_'+key)
    numeric=num['projected'](raw)
    for key in physical:equal_rows(numeric[key],physical[key],'NUMERIC_PHYSICAL_'+key)
    def identity_interval(lo=0,hi=None):
        if hi is not None:need(lo==hi,'SYMBOLIC_EXACT_SHELL_INPUT')
        return lo
    def actual_jet(x,r):
        need(r==0,'ACTUAL_WRAPPER_SHELL_ONLY')
        return {'A':h,'D':-2*X,'E':V,'L':2*L,'M':2*d},{'s':v},{'a':a,'p':p}
    data={'shell':{'boxes':{'eta_b':{'exact':[eta,eta]}}}}
    num['j']=SimpleNamespace(I=identity_interval,actual_jet=actual_jet,inputs=lambda:(data,None,None))
    actual=num['actual_shell_sources'](P.var('x'))
    if mutation=='wrong_two_copy_factor':actual['two_copy_outward_hessian']={key:tuple(z/2 for z in val) for key,val in actual['two_copy_outward_hessian'].items()}
    for key in physical:
        equal_rows(actual['one_copy_outward_hessian'][key],physical[key],'ACTUAL_NUMERIC_ONE_COPY_'+key)
        equal_rows(actual['two_copy_outward_hessian'][key],tuple(2*z for z in physical[key]),'ACTUAL_NUMERIC_TWO_COPY_'+key)
        equal_rows(actual['two_copy_inward_hessian'][key],tuple(-2*z for z in physical[key]),'ACTUAL_NUMERIC_TWO_COPY_INWARD_'+key)
        equal_rows(actual['coordinate_physical_forcing'][key],tuple((-1 if key=='scalar' else -2)*z for z in physical[key]),'ACTUAL_NUMERIC_PHYSICAL_FORCING_'+key)
    # Independently expand full-metric-lowered normal C_y, including reference connection.
    Y,G,H2=(P.var(n) for n in ('Y','G','H2'))
    H=h*Y;scalar=v*Y;tau=(2*X+6*L-16*d)*Y
    trA2=(4*X*X+12*L*L-64*L*d)*Y*Y-2*V*V*G+4*d*d*H2
    if mutation=='positive_time_angular_trace':trA2=trA2+4*V*V*G
    if mutation=='omit_anisotropic_trace':trA2=trA2-4*d*d*H2
    B=(a*h+sigma1*v/3)*Y
    Hy=-4*a*H+2*a*tau+4*sigma1*scalar/3
    traceMy=2*a*tau+4*B;traceAMy=2*a*trA2+B*tau
    gyy=1+ep*H
    tangent_inv_trace=S([4,-tau,trA2])
    tangent_dy_trace=S([8*a,traceMy-2*a*tau,2*a*trA2-traceAMy])
    reference=a*tangent_inv_trace
    if mutation=='omit_quadratic_reference_connection':reference=S([4*a,-a*tau,0])
    Cup=ep*Hy/(2*gyy**2)-tangent_dy_trace/(2*gyy)+reference
    Clow=gyy*Cup
    eq(Clow.c[0],0,'NORMAL_GAUGE_CONSTANT')
    eq(Clow.c[1],0,'NORMAL_GAUGE_FIRST')
    qc_expected=project_angular_scalar(2*Clow.c[2])
    maps=json.loads((root/PREFIX/'projection/GEOMETRIC_COEFFICIENT_MAPS.json').read_text())
    need(maps['gauge_convention']=='QC=second variation of full metric-lowered gauge covector; QC_background also retained','MAP_GAUGE_CONVENTION')
    values={n:P.var('free_'+n) for n in maps['input_symbols']}
    values.update(A=h,B=P(),C=P(),D=-2*X,E=V,L=2*L,M=2*d,
                  Ay=-4*a*h+2*a*(2*X+6*L-16*d)+4*sigma1*v/3,
                  Dy=-beta,Ly=beta,Ey=P(),My=P(),Bt=P(),Ct=P(),
                  rho=rho,c=c,theta=theta,a=a,p=p)
    # Before imposing the first normal harmonic row, the universal full-
    # metric lowering and background lowering differ by 2 H C1_y.
    Hy0=P.var('Hy');Hy_free=Hy0*Y
    Cup_free=ep*Hy_free/(2*gyy**2)-tangent_dy_trace/(2*gyy)+reference
    Clow_free=gyy*Cup_free
    C1free=Hy_free/2+2*a*H-a*tau-2*sigma1*scalar/3
    eq(Clow_free.c[1],C1free,'GENERIC_FIRST_NORMAL_HARMONIC')
    generic_values=dict(values,Ay=Hy0)
    qc_generic=tuple(encoded(row,generic_values) for row in maps['maps']['QC_actual']['y'])
    qc_background=tuple(encoded(row,generic_values) for row in maps['maps']['QC_background']['y'])
    equal_rows(qc_generic,project_angular_scalar(2*Clow_free.c[2]),'GENERIC_QC_ACTUAL')
    equal_rows(qc_background,project_angular_scalar(2*Cup_free.c[2]),'GENERIC_QC_BACKGROUND')
    difference=project_angular_scalar(2*H*C1free)
    for i,(actual_q,background_q,cross) in enumerate(zip(qc_generic,qc_background,difference)):
        eq(actual_q-background_q,0 if mutation=='equate_lowerings_without_gauge' else cross,'QC_LOWERING_DIFFERENCE_'+str(i))
    if mutation=='coordinate_instead_of_frame_normal_derivative':values['Dy']=-beta-4*a*X
    qc_map=tuple(encoded(row,values) for row in maps['maps']['QC_actual']['y'])
    equal_rows(qc_map,qc_expected,'QC_MAP_SHELL')
    equal_rows(source['QC_y'],qc_expected,'QC_ADAPTER_SHELL')
    traceQR=tuple((-1 if mutation!='positive_QR_trace' else 1)*x+3*y for x,y in zip(source['QK00'],source['QKab'][:3]))
    for i,(gn,tr,qc) in enumerate(zip(source['GN'],traceQR,qc_expected)):eq(gn,-2*tr-2*qc,'GN_LORENTZ_TRACE_'+str(i))
    # Independent all23 physical/source packing, outward-y to inward-r signs.
    forcing=op.boundary_forcing(1 if mutation=='missing_ell' else ell,rho,source)
    if mutation=='wrong_inward_sign':forcing=[-z if i in (0,1,2) else z for i,z in enumerate(forcing)]
    expected={}
    for j,l in enumerate((0,2,4)):
        expected['K:'+str(l)]=-ell*rho*source['GN'][j]
        expected['N:'+str(l)]=2*ell*rho*physical['00'][j]
        expected['w:'+str(l)]=ell*rho*physical['scalar'][j]
        expected['L:'+str(l)]=P()
    for j,l in enumerate((2,4)):
        expected['u:'+str(l)]=2*ell*rho*physical['0a'][j]
        expected['v:'+str(l)]=P()
    for key,z in zip(('J:0','J:2','J:4','H:2','H:4','T:+','T:-'),physical['ab']):expected[key]=2*ell*rho*z
    need(len(expected)==23 and set(expected)==set(op.LABELS),'FULL23_LABELS')
    for key,z in zip(op.LABELS,forcing):eq(z,expected[key],'BOUNDARY_'+key)
    # Exact off-shell Codazzi scalar cancellation BEFORE physical P1=0.
    divK,dtrK,nphi,dphi=(P.var(n) for n in ('divK','dtrK','nphi','dphi'))
    divR_minus_dtrR=divK-dtrK+sigma1*dphi/2
    Pphysical=nphi+sigma1/2
    right=divR_minus_dtrR-(0 if mutation=='off_shell_scalar_term_omitted' else 1)*Pphysical*dphi
    eq(right,divK-dtrK-nphi*dphi,'OFFSHELL_SCALAR_CODAZZI')
    # Induced lower-covector Robin conversion, not additional primal data.
    Ca,Cay,Chat,Chat_y=(P.var(n) for n in ('Ca','Cay','Chat','Chat_y'))
    sym=(Cay-a*Ca-a*Ca)/2
    coordinate=(Cay-(a if mutation=='wrong_coordinate_Robin' else 2*a)*Ca)/2
    eq(sym,coordinate,'COORDINATE_COVECTOR_ROBIN')
    frame_sym=(rho*Chat_y+a*rho*Chat-2*a*rho*Chat)/(2*rho)
    frame_expected=(Chat_y-(2*a if mutation=='wrong_frame_Robin' else a)*Chat)/2
    eq(frame_sym,frame_expected,'PARALLEL_FRAME_COVECTOR_ROBIN')
    need(authenticate(root)==before,'SOURCE_CHANGED_AFTER_READBACK')
    return {'status':'PASS_EXACT_PHYSICAL_SHELL_SOURCE_AND_NORMALIZATION_BINDING',
      'source_SHA256':before,'named_checks':len(rows),'checks':rows,
      'normalization':'second parameter derivatives, one-copy outward QR/QP; ordinary physical frame y,t slots; every PDE source uses stated ell*rho factor',
      'physical_modal_source_coefficients_per_copy':15,'full23_primal_forcing_entries':23,
      'Dirichlet_modal_slots':5,'Robin_modal_slots':18,
      'normal_gauge_map_rows':3,'off_shell_scalar_codazzi_checked_before_P1_zero':True,
      'fullmetric_vs_background_lowering_checked_with_free_Hy':True,
      'source_closure':'two authenticated imported implementation modules use only standard-library imports; numeric wrapper is authenticated AST without actual-jet imports',
      'scope':{'actual_first_trace_interface_mapping':True,'exact_free_shell_profile_identities':True,
       'first_physical_shell_rows_and_normal_harmonic_row_retained_analytic_inputs':True,
       'normal_differentiation_of_shell_only_equations':False,
       'actual_first_wave_PDE_or_global_subsidiary_proved_here':False,
       'natural_QE_Codazzi_five_rows_checked_here':False,
       'finite_amplitude_nonlinear_or_cosmological_realization':False},
      'normal_gauge_reference_profiles':[x.encoded() for x in qc_expected]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--mutation',choices=CONTROLS);ap.add_argument('--run-controls',action='store_true');args=ap.parse_args()
    need(not args.output.exists(),'REFUSE_OVERWRITE')
    report=audit(args.repo_root.resolve(),args.mutation)
    if args.run_controls:
        rejected=[]
        with tempfile.TemporaryDirectory() as tmp:
            for mutation,reason in CONTROLS.items():
                out=Path(tmp)/(mutation+'.json')
                command=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])+[str(Path(__file__).resolve()),'--repo-root',str(args.repo_root.resolve()),'--output',str(out),'--mutation',mutation]
                r=subprocess.run(command,capture_output=True,text=True)
                need(r.returncode!=0 and 'NONZERO:'+reason in r.stderr and not out.exists(),'CONTROL_NOT_REJECTED:'+mutation+':'+r.stderr[-1000:])
                rejected.append({'mutation':mutation,'expected_nonzero_identity':reason,'returncode':r.returncode,'successful_report_absent':True})
        report['rejected_corrupted_executions']=rejected
    report['checker_SHA256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'checks':report['named_checks'],'corruptions':len(report.get('rejected_corrupted_executions',[]))}))

if __name__=='__main__':main()
