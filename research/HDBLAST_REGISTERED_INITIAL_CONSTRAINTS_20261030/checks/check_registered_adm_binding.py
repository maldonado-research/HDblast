#!/usr/bin/env python3
"""Exact all-point ADM-to-registered full23 assembly binding; stdlib only.

The independent producer expands the exact rational lapse/shift prescription.
The implementation side compiles only the authenticated _registered assembly
suffix after its physical profiles have been formed. Those profiles are free
indeterminates here; their independent physical binding is a separate audit.
"""
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

PINS = {
 'research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/theory/CONSTRAINT_FORMULA_NOTE.md':'b53ccfd7c638c9b98df23803b50a6821ae46416adfb95f8a68a9e35b48b93f4f',
 'research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/reference/FORMULAS.json':'f3f037f890655ab05f936491d378a494b56e10f0c8a3d32dda0f1b8b2ea7556f',
 'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md':'bc3ae85c3166575ed15b017468a3e9668615e15dfaca34bb2a13afe27ea28f92',
 'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md':'2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/source_adapters.py':'0233ff25b9873dfc4060ba0d4963f50d0b198dbe71368625017840d225b19f8a',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py':'a008be0d89e2ac44940c273cfd2cbc5ba0d09de8fc172d5d81ef29121ab60189',
 'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py':'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1',
}
CONTROLS = {
 'omit_lapse_value':'REGISTERED_VALUE_N:0',
 'omit_mixed_conversion':'MIXED_VELOCITY',
 'halve_lapse_cross':'LAPSE_SECOND',
 'omit_fy':'RADIAL_SHIFT_SECOND',
 'halve_angular':'ANGULAR_SHIFT_SECOND',
 'omit_lapse_velocity_cross':'TIME_VELOCITY',
 'missing_ell_scaling':'REGISTERED_VELOCITY_K:0',
 'changed_mean_override':'REGISTERED_VALUE_N:0',
 'wrong_vector_projection':'REGISTERED_VELOCITY_u:2',
 'wrong_trace_reversal':'TRACE_REVERSED_NORMAL',
 'use_Enn_alone':'PHYSICAL_HAMILTONIAN_FACTOR',
 'wrong_future_normal':'INITIAL_TIME_NORMAL_FACTOR',
 'use_longitudinal_first_jet':'FIRST_HARMONIC_INITIAL_RADIAL',
 'wrong_initial_QC_radial_sign':'INITIAL_QC_DIRECT_y:0',
}

def need(ok, why):
    if not ok: raise ValueError(why)

class L:
    """Collected exact Laurent polynomial with dynamic named generators."""
    def __init__(self, terms=None):
        if isinstance(terms, L): self.t=dict(terms.t)
        elif isinstance(terms,dict): self.t={k:F(v) for k,v in terms.items() if v}
        else: self.t={} if not terms else {():F(terms)}
    @classmethod
    def cv(cls,x): return x if isinstance(x,cls) else cls(x)
    @classmethod
    def var(cls,name): return cls({((name,1),):F(1)})
    def __add__(self,x):
        if isinstance(x,S): return NotImplemented
        t=dict(self.t)
        for k,v in self.cv(x).t.items():t[k]=t.get(k,F(0))+v
        return L(t)
    __radd__=__add__
    def __neg__(self):return L({k:-v for k,v in self.t.items()})
    def __sub__(self,x):return self+-self.cv(x)
    def __rsub__(self,x):return self.cv(x)+-self
    def __mul__(self,x):
        if isinstance(x,S): return NotImplemented
        t={}
        for k,v in self.t.items():
            for l,w in self.cv(x).t.items():
                key=dict(k)
                for name,power in l:key[name]=key.get(name,0)+power
                key=tuple(sorted((n,p) for n,p in key.items() if p))
                t[key]=t.get(key,F(0))+v*w
        return L(t)
    __rmul__=__mul__
    def inverse(self):
        need(len(self.t)==1,'NONMONOMIAL_DENOMINATOR')
        (k,v),=self.t.items()
        return L({tuple((n,-p) for n,p in k):1/v})
    def __truediv__(self,x):
        if isinstance(x,S):return NotImplemented
        return self*self.cv(x).inverse()
    def __rtruediv__(self,x):return self.cv(x)*self.inverse()
    def __pow__(self,n):
        need(isinstance(n,int),'INTEGER_POWER')
        if n<0:return self.inverse()**(-n)
        out=L(1)
        for _ in range(n):out=out*self
        return out
    def __eq__(self,x):return self.t==self.cv(x).t
    def __bool__(self):return bool(self.t)
    def encode(self):return [[[[n,p] for n,p in k],str(v)] for k,v in sorted(self.t.items())]
    def packet_degree(self):
        packet={'U','Uy','T','Tp','v','f0','f2','f4','fy0','fy2','fy4','q','Pi'}
        return [sum(p for n,p in k if n in packet) for k in self.t]

class S:
    def __init__(self,x=0):
        self.c=tuple(L.cv(y) for y in x) if isinstance(x,(tuple,list)) else (L.cv(x),L(0),L(0))
        need(len(self.c)==3,'SERIES_ORDER')
    @classmethod
    def cv(cls,x):return x if isinstance(x,cls) else cls(x)
    def __add__(self,x):return S([a+b for a,b in zip(self.c,self.cv(x).c)])
    __radd__=__add__
    def __neg__(self):return S([-a for a in self.c])
    def __sub__(self,x):return self+-self.cv(x)
    def __rsub__(self,x):return self.cv(x)+-self
    def __mul__(self,x):
        x=self.cv(x)
        return S([sum((self.c[j]*x.c[k-j] for j in range(k+1)),L(0)) for k in range(3)])
    __rmul__=__mul__
    def inverse(self):
        a,b,c=self.c
        return S([1/a,-b/a**2,b*b/a**3-c/a**2])
    def __truediv__(self,x):return self*self.cv(x).inverse()
    def __rtruediv__(self,x):return self.cv(x)*self.inverse()
    def __pow__(self,n):
        if n<0:return self.inverse()**(-n)
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
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod


def registered_assembly(root, adapter,op,mutation):
    path=root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py'
    tree=ast.parse(path.read_text())
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='Context')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='_registered')
    start=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='zero' for t in n.targets))
    suffix=fn.body[start:]
    # Bind only physical profile inputs, not the prefix calculation. No AST
    # simplification rewrites the implementation expressions under test.
    names=['self','u','uy','ut','vt','fs','fy','K','Pi','rho','p','a','packet']
    new=ast.FunctionDef(name='assembly',args=ast.arguments(posonlyargs=[],args=[ast.arg(arg=n) for n in names],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=ast.parse('ip=1/p\ndy=lambda value:uy').body+suffix,decorator_list=[])
    module=ast.fix_missing_locations(ast.Module(body=[new],type_ignores=[]))
    namespace={'F':F};exec(compile(module,str(path),'exec'),namespace)
    ctx=SimpleNamespace(adapter=adapter,op=op,ell=L.var('ell'),mean_ids=[0,3,6,9,12])
    sym={name:L.var(name) for name in ('U','Uy','T','Tp','v','q','Pi','rho','p','a','packet')}
    fs=[L.var('f'+str(l)) for l in (0,2,4)]
    fy=[L.var('fy'+str(l)) for l in (0,2,4)]
    values,vel=namespace['assembly'](ctx,sym['U'],sym['Uy'],sym['T'],sym['Tp'],fs,fy,sym['q'],sym['Pi'],sym['rho'],sym['p'],sym['a'],sym['packet'])
    if mutation=='omit_lapse_value':values[op.IDX['N:0']]=L(0)
    if mutation=='missing_ell_scaling':vel=[x/L.var('ell') for x in vel]
    if mutation=='changed_mean_override':values[op.IDX['N:0']]*=2
    if mutation=='wrong_vector_projection':vel[op.IDX['u:2']]+=L.var('ell')*sym['U']**2
    return values,vel,hashlib.sha256(ast.dump(ast.Module(body=suffix,type_ignores=[]),include_attributes=False).encode()).hexdigest()


def check(root,mutation=None):
    before=authenticate(root)
    ad=load('binding_source_adapters',root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/source_adapters.py')
    op=load('binding_linear_operator',root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py')
    checks=[]
    def ck(condition,name):need(condition,name);checks.append(name)
    def zero(value,name):ck(not L.cv(value),name)
    def eq(left,right,name):zero(left-right,name)
    s={n:L.var(n) for n in ('U','Uy','T','a','f','fy','Ui','fi','q','rho','Tp','Pi','p','ell')}
    U,Uy,T,a,f,fy,Ui,fi,q,rho,Tp,Pi,p,ell=[s[n] for n in s]
    ep=S([0,1,0]);A=1-4*ep*U;B=1+2*ep*U+ep**2*f/2;P=1+ep*U
    Ay=-4*ep*Uy;By=2*ep*Uy+ep**2*fy/2;Py=ep*Uy
    Ai=-4*ep*Ui;Bi=2*ep*Ui+ep**2*fi/2;Pj=ep*Ui
    Kyy=-2*ep*T/rho+ep**2*q/2;Kab=ep*T/rho+ep**2*q/2
    trK=Kyy/A+3*Kab/B
    dotN=rho**2*P**2*trK
    gamma_y=Ay/(2*A**2)-3*By/(2*A*B)-3*a/A
    reference_y=-S(a)/P**2-3*S(a)/B
    dotBy=rho**2*P**2*(gamma_y-reference_y-S(a)/A-Py/(A*P))
    sphere_bracket=Ai/(2*A)+Bi/(2*B)+Pj/P
    dotBi=-P**2/B*sphere_bracket
    n2=4*rho**2*q-24*rho*U*T
    by1=-6*rho**2*(Uy+4*a*U)
    by2=rho**2*(-74*U*Uy-F(3,2)*fy-194*a*U**2-3*a*f)
    bi2=22*U*Ui-fi/2
    if mutation=='halve_lapse_cross':n2+=12*rho*U*T
    if mutation=='omit_fy':by2+=F(3,2)*rho**2*fy
    if mutation=='halve_angular':bi2-=11*U*Ui
    for key,series,expected in [('LAPSE',dotN,(0,rho*T,n2)),('RADIAL_SHIFT',dotBy,(0,by1,by2)),('ANGULAR_SHIFT',dotBi,(0,0,bi2))]:
        for j,(coeff,want) in enumerate(zip(series.c,expected)):
            eq(coeff*(2 if j==2 else 1),want,key+'_'+('CONSTANT','FIRST','SECOND')[j])
    # Substitute the prescribed derivatives into the independent contracted
    # Christoffel ADM identities. This is a generic all-point coefficient test.
    Ct=trK/(rho*P)-dotN/(rho**3*P**3)
    Cy=-dotBy/(rho**2*P**2)+gamma_y-reference_y-S(a)/A-Py/(A*P)
    Ci=-dotBi/(rho**2*P**2)-sphere_bracket/(rho**2*B)
    for kind,series in [('t',Ct),('y',Cy),('a',Ci)]:
        for j,coeff in enumerate(series.c):zero(coeff,'INITIAL_HARMONIC_'+kind+'_ORDER_'+str(j))
    # Generic polynomial exact C=0 supplies scalar rows l=0,2,4 for t,y,
    # and vector-gradient rows l=2,4 for angular directions.
    for kind,degrees in [('t',(0,2,4)),('y',(0,2,4)),('a',(2,4))]:
        for degree in degrees:zero({'t':Ct,'y':Cy,'a':Ci}[kind].c[2],'INITIAL_HARMONIC_MODAL_'+kind+':'+str(degree))
    # Independent conversion of g=-N²dt²+gamma(dx+beta dt)², beta=0.
    N=rho*P
    gtt=-N**2;gtt2=2*gtt.c[2]/rho**2
    eq(gtt2,-2*U**2,'TIME_LAPSE_VALUE')
    dotgtt=-2*N*dotN
    time2=-8*rho*q+44*U*T
    if mutation=='omit_lapse_velocity_cross':time2+=4*U*T
    eq(2*dotgtt.c[2]/rho**2,time2,'TIME_VELOCITY')
    dotgty=A*dotBy
    mixed2=rho*(-26*U*Uy-F(3,2)*fy-2*a*U**2-3*a*f)
    if mutation=='omit_mixed_conversion':mixed2=by2/rho
    eq(2*dotgty.c[2]/rho,mixed2,'MIXED_VELOCITY')
    eq(2*(B*dotBi).c[2],22*U*Ui-fi/2,'ANGULAR_VELOCITY')
    eq(2*(2*N*Kyy).c[2],2*rho*q-8*U*T,'RADIAL_SPATIAL_VELOCITY')
    eq(2*(2*N*Kab).c[2],2*rho*q+4*U*T,'ANGULAR_SPATIAL_VELOCITY')
    scalarPi=ep*Tp/rho+ep**2*Pi/2
    eq(2*(N*scalarPi).c[2],rho*Pi+2*U*Tp,'SCALAR_VELOCITY')
    # Explicit probability-normalized scalar/vector product identities.
    Y2=(F(1,6),F(1,2),F(1));vector=(F(1,4),F(1,2))
    ck(tuple(ad.Y2)==Y2,'Y2_SCALAR_PRODUCT_TABLE')
    ck(ad.vector_project(1)==vector,'YDY_VECTOR_PRODUCT_TABLE')
    labels=[f'{field}:{l}' for field in ('K','L','N','J','w') for l in (0,2,4)]+[f'{field}:{l}' for field in ('v','u') for l in (2,4)]+['H:2','H:4','T:+','T:-']
    ck(op.LABELS==labels,'ALL23_CANONICAL_ORDER')
    actual_u,actual_v,suffixhash=registered_assembly(root,ad,op,mutation)
    fs=[L.var('f'+str(l)) for l in (0,2,4)];fys=[L.var('fy'+str(l)) for l in (0,2,4)]
    independent_u=[];independent_v=[]
    for label in labels:
        field,degree=label.split(':')
        if degree in ('+','-') or field in ('H','v'):
            uv,vv=L(0),L(0)
        elif field=='u':
            j=(2,4).index(int(degree));uv=L(0);vv=22*U**2*vector[j]-fs[j+1]/2
        else:
            j=(0,2,4).index(int(degree));c=Y2[j]
            if field=='K':uv=L(0);vv=(2*rho*q-8*U*T)*c
            elif field=='L':uv=L(0);vv=rho*((-26*U*Uy-2*a*U**2)*c-F(3,2)*fys[j]-3*a*fs[j])
            elif field=='N':uv=-2*U**2*c;vv=(-8*rho*q+44*U*T)*c
            elif field=='J':uv=fs[j];vv=(2*rho*q+4*U*T)*c
            elif field=='w':uv=-F(3,2)*fys[j]/p;vv=(rho*Pi+2*U*Tp)*c
            else:raise ValueError('UNKNOWN_EXPECTED_LABEL')
        independent_u.append(uv);independent_v.append(ell*vv)
    for i,label in enumerate(labels):
        eq(actual_u[i],independent_u[i],'REGISTERED_VALUE_'+label)
        eq(actual_v[i],independent_v[i],'REGISTERED_VELOCITY_'+label)
    # Adapter is checked independently as a third join, including its ordinary
    # t velocities followed by every-slot ell conversion and manual label map.
    data=ad.initial_frame_data(rho=rho,a=a,p=p,u=U,uy=Uy,ut=T,v=L.var('v'),vt=Tp,Fmom=3*q,Gmom=p*Pi+L.var('Fy'),Fmom_y=L.var('Fy'),f=fs,fy=fys)
    packed_u,packed_v=op.pack_initial(data,ell)
    for i,label in enumerate(labels):
        eq(packed_u[i],independent_u[i],'ADAPTER_VALUE_'+label)
        eq(packed_v[i],independent_v[i],'ADAPTER_VELOCITY_'+label)
    for i in (0,3,6,9,12):
        eq(actual_u[i],packed_u[i],'MEAN_OVERRIDE_VALUE_'+labels[i])
        eq(actual_v[i],packed_v[i],'MEAN_OVERRIDE_VELOCITY_'+labels[i])
    for label in ('H:2','H:4','T:+','T:-'):
        zero(actual_u[labels.index(label)],'ZERO_TENSOR_VALUE_'+label)
        zero(actual_v[labels.index(label)],'ZERO_TENSOR_VELOCITY_'+label)
    # Direct eight-row QC initial source-map calculation, with the corrected
    # first harmonic one-jet. This is separate from the generic ADM identity
    # above, and it does not replace the analytic first gauge-wave initial match.
    m={key:L(0) for key in ('A','B','C','D','E','L','M','Ay','At','By','Bt','Cy','Ct','Dy','Dt','Ey','Et','Ly','Lt','My','Mt')}
    bt=-6*rho*(Uy+4*a*U)
    if mutation=='use_longitudinal_first_jet':bt=L(0)
    m.update(A=-4*U,D=-2*U,L=2*U,Ay=-4*Uy,Dy=-2*Uy,Ly=2*Uy,
             At=-4*T,Dt=-2*T,Lt=2*T,Bt=bt)
    zero((-4*Uy-2*Uy-3*2*Uy)/2-bt/rho+a*(-16*U-2*U-6*U),'FIRST_HARMONIC_INITIAL_RADIAL')
    zero(-(-4*T-2*T+3*2*T)/(2*rho),'FIRST_HARMONIC_INITIAL_TIME')
    zero((4*U-2*U-2*U)/(2*rho),'FIRST_HARMONIC_INITIAL_ANGULAR')
    qc=ad.gauge_source(m,{'rho':rho,'c':L(1),'a':a,'theta':L(0)})
    if mutation=='wrong_initial_QC_radial_sign':qc['y']=tuple(-x for x in qc['y'])
    # Ordinary physical-frame DC formula at t=0 (tanh t=0, sech t=1):
    # Cy=.5(Ky+Ny-3Jy)-Lt/rho+a(4K+N-3J)-lambda*v/rho,
    # C0=Ly-(Kt+Nt+3Jt)/(2rho)+5aL,
    # Ca=vy-ut/rho+5av+(-K+N-J+2beta H)/(2rho).
    # Here first y is outward radial coordinate y; full23 initial angular
    # v,H,L vanish. No reversed r derivative remains in these physical rows.
    direct={'y':[],'0':[],'a':[]}
    for j,degree in enumerate((0,2,4)):
        c=Y2[j]
        cy=(-4*U*Uy*c-3*fys[j])/2-rho*((-26*U*Uy-2*a*U**2)*c-F(3,2)*fys[j]-3*a*fs[j])/rho+a*(-2*U**2*c-3*fs[j])
        ctime=-((2*rho*q-8*U*T)+(-8*rho*q+44*U*T)+3*(2*rho*q+4*U*T))*c/(2*rho)
        direct['y'].append(cy);direct['0'].append(ctime)
        zero(cy+qc['y'][j],'INITIAL_QC_DIRECT_y:'+str(degree))
        zero(ctime+qc['0'][j],'INITIAL_QC_DIRECT_t:'+str(degree))
    for j,degree in enumerate((2,4)):
        ca=-(22*U**2*vector[j]-fs[j+1]/2)/rho+(-2*U**2*Y2[j+1]-fs[j+1])/(2*rho)
        direct['a'].append(ca)
        zero(ca+qc['a'][j],'INITIAL_QC_DIRECT_a:'+str(degree))
    # Universal five-dimensional trace reversal, future timelike n. These
    # are normal projections, not direct differentiated Christoffel evolution.
    X=[L.var('X'+str(i)) for i in range(5)]
    sig=[-1,1,1,1,1]
    E=[[L(0) for _ in range(5)] for _ in range(5)]
    E[0][0]=X[0]
    for i in range(1,5):E[0][i]=E[i][0]=X[i]/2
    trace=sum((sig[i]*E[i][i] for i in range(5)),L(0))
    eq(trace,-X[0],'REDUCED_NORMAL_TRACE')
    Hnn=E[0][0]+trace/2
    if mutation=='wrong_trace_reversal':Hnn=E[0][0]+trace
    eq(Hnn,X[0]/2,'TRACE_REVERSED_NORMAL')
    physicalH=2*Hnn
    if mutation=='use_Enn_alone':physicalH=2*E[0][0]
    eq(physicalH,X[0],'PHYSICAL_HAMILTONIAN_FACTOR')
    for i in range(1,5):eq(E[0][i],X[i]/2,'TRACE_REVERSED_MOMENTUM_'+str(i))
    normal=rho
    if mutation=='wrong_future_normal':normal=-rho
    eq(normal*X[0],rho*X[0],'INITIAL_TIME_NORMAL_FACTOR')
    # Structural collar statement: every coefficient lies in the ideal
    # generated by first packet variables and the physical correction fields.
    # Physical audit supplies f,fy,q,Pi in the packet-jet ideal. Derivation of
    # an ideal element remains an element of its differential ideal; finite
    # initial rows alone are not an all-orders corner/evolution proof.
    for i,label in enumerate(labels):
        ck(all(d>=1 for d in actual_u[i].packet_degree()),'PACKET_IDEAL_VALUE_'+label)
        ck(all(d>=1 for d in actual_v[i].packet_degree()),'PACKET_IDEAL_VELOCITY_'+label)
    ck(authenticate(root)==before,'SOURCE_READBACK_UNCHANGED')
    return {
      'status':'PASS_EXACT_REGISTERED_ADM_ASSEMBLY_BINDING',
      'method':'collected rational Laurent polynomials and truncated parameter series, no point samples',
      'parameter_convention':'g=g0+epsilon*g1+epsilon^2*g2/2',
      'derivative_convention':'ordinary t adapter velocities; every slot multiplied by ell for x=t/ell',
      'labels':labels,'source_sha256':before,'registered_assembly_suffix_AST_sha256':suffixhash,
      'registered_assembly_entries_checked':46,'adapter_entries_checked':46,
      'five_mean_overrides_checked_in_both_arrays':True,
      'zero_initial_TT_and_HessTF_slots':True,
      'initial_harmonic_modal_rows':8,
      'direct_initial_QC_source_map_rows':8,
      'generic_harmonic_modal_rows_are_consequences_of_exact_ADM_identity':True,
      'initial_subsidiary_velocity':'conditional universal trace-reversal implication, not an independent reduced-acceleration or QC_x calculation',
      'corner_scope':'finite registered one-jet belongs to packet differential ideal; smooth flat support and homogeneous recursion remain retained analytic premises',
      'retained_prerequisites':['independent physical-profile binding aa,dd,3K=Fmom,pPi=Gmom-D_yFmom and f/f_y','actual first harmonic wave initial one-jet matches exact ADM prescription','reduced Einstein accelerations for subsidiary normal derivative','physical scalar/reduced Bianchi and shell Codazzi/subsidiary uniqueness beyond this initial audit'],
      'not_proved':['global subsidiary identity from these finite rows alone','exact finite-amplitude nonlinear family','physical relative accuracy','higher-dimensional cosmological origin'],
      'checks':checks,'named_checks':len(checks),
      'direct_initial_QC_profiles':{key:[x.encode() for x in values] for key,values in qc.items()},
      'independent_initial_DC_profiles':{key:[x.encode() for x in values] for key,values in direct.items()},
      'independent_expected_registered_values':[x.encode() for x in independent_u],
      'independent_expected_registered_scaled_velocities':[x.encode() for x in independent_v],
    }


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--mutation',choices=CONTROLS);ap.add_argument('--run-controls',action='store_true');args=ap.parse_args()
    need(not args.output.exists(),'REFUSE_OVERWRITE')
    report=check(args.repo_root.resolve(),args.mutation)
    if args.run_controls:
        results=[]
        with tempfile.TemporaryDirectory() as tmp:
            for name,reason in CONTROLS.items():
                target=Path(tmp)/(name+'.json');cmd=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])+[str(Path(__file__).resolve()),'--repo-root',str(args.repo_root.resolve()),'--output',str(target),'--mutation',name]
                result=subprocess.run(cmd,capture_output=True,text=True)
                need(result.returncode!=0 and reason in result.stderr and not target.exists(),'CONTROL_NOT_REJECTED:'+name)
                results.append({'mutation':name,'expected_reason':reason,'returncode':result.returncode,'report_absent':True})
        report['rejected_corrupted_executions']=results
    report['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'checks':report['named_checks'],'controls':len(report.get('rejected_corrupted_executions',[]))}))

if __name__=='__main__':main()
