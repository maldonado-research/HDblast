#!/usr/bin/env python3
"""Exact registered first-wave pullback/operator binding, standard library only.

The independent side forms coordinate Lie derivatives and contracts the full
parallel-frame connection. No retained harmonic-row or wave-matrix formula is
used to produce the expected covector. Authentic source ASTs supply the
implementation expressions being tested.
"""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

WAVE='research/HDBLAST_FIRST_GAUGE_WAVE_20261016'
MANIFEST_PIN='adb6e3a7e00a8e477af80a49bc6b7a0f4037fafd0dd42ab8175c8b55e3e4465c'
PINS={
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py':'a68fa8025d2170efca6c3e85d443910a796d0c30a9b7dd17ed988cd9e00bbb00',
 WAVE+'/operator/check_gauge_operator.py':'163553f2bdb9469b3ae40bd8f1c9bcf01077b690abac9a2a7b2513619ac99011',
 WAVE+'/coefficients/check_global_gauge_coefficients.py':'5c7e7f51b0257ef6d0d0ee861407f52d8b6a9a480e61367cc84dbc091aa8c42e',
 'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/first-fields/FIRST_HARMONIC_FIELDS_UNPUBLISHED.md':'d15f8c65ffffd0085d195064ccea3c61bb4e719695652d0dd73472c8a06c224d',
}


def need(ok,reason):
    if not ok:raise ValueError(reason)


class P:
    """Collected sparse rational Laurent ring, with named generators."""
    def __init__(self,x=0):
        if isinstance(x,P):self.t=dict(x.t)
        elif isinstance(x,dict):self.t={tuple(k):F(v) for k,v in x.items() if v}
        else:self.t={():F(x)} if x else {}
    @classmethod
    def var(cls,n):return cls({((n,1),):F(1)})
    def __add__(self,b):
        t=dict(self.t)
        for k,v in P(b).t.items():t[k]=t.get(k,F(0))+v
        return P(t)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.t.items()})
    def __sub__(self,b):return self+-P(b)
    def __rsub__(self,b):return P(b)+-self
    def __mul__(self,b):
        t={}
        for k,v in self.t.items():
            for l,w in P(b).t.items():
                m=dict(k)
                for n,e in l:m[n]=m.get(n,0)+e
                m=tuple(sorted((n,e) for n,e in m.items() if e))
                t[m]=t.get(m,F(0))+v*w
        return P(t)
    __rmul__=__mul__
    def __pow__(self,n):
        need(isinstance(n,int),'INTEGER_POWER')
        if n<0:
            need(len(self.t)==1,'NONMONOMIAL_INVERSE')
            (k,v),=self.t.items()
            return P({tuple((s,e*n) for s,e in k):v**n})
        r=P(1)
        for _ in range(n):r=r*self
        return r
    def __truediv__(self,b):return self*P(b)**-1
    def __rtruediv__(self,b):return P(b)*self**-1
    def __bool__(self):return bool(self.t)
    def __eq__(self,b):return self.t==P(b).t
    def partial(self,n):
        out={}
        for k,v in self.t.items():
            d=dict(k);e=d.get(n,0)
            if e:
                d[n]=e-1
                kk=tuple(sorted((s,p) for s,p in d.items() if p))
                out[kk]=v*e
        return P(out)
    def derivation(self,rules):
        names={n for k in self.t for n,e in k}
        return sum((self.partial(n)*rules.get(n,P()) for n in names),P())
    def spacetime_derivative(self,direction):
        table=rules(direction)
        for monomial in self.t:
            for name,power in monomial:
                if name.split('_')[0] in ('f','b','d','W','J') and '_' in name:
                    need(name in table,'JET_ORDER_EXCEEDED:'+name)
        return self.derivation(table)
    def dx(self):return self.spacetime_derivative('x')
    def dr(self):return self.spacetime_derivative('r')
    def derivative(self,ix,ir):
        out=self
        for _ in range(ix):out=out.dx()
        for _ in range(ir):out=out.dr()
        return out
    def val(self):return self
    def sub(self,table):
        out=P()
        for k,v in self.t.items():
            term=P(v)
            for n,e in k:term*=P(table.get(n,P.var(n)))**e
            out+=term
        return out
    def encode(self):return [[[[n,e] for n,e in k],str(v)] for k,v in sorted(self.t.items())]


def field(name,ix=0,ir=0):return P.var(name+'_'+str(ix)+str(ir))


def rules(direction):
    ell,rho,a,ay,c,H,Rr,theta,p,py=map(P.var,('ell','rho','a','ay','c','H','Rr','theta','p','py'))
    out={}
    for n in ('f','b','d','W','J'):
        for i in range(4):
            for j in range(4-i):out[n+'_'+str(i)+str(j)]=field(n,i+(direction=='x'),j+(direction=='r'))
    if direction=='x':out.update(c=ell*H*c,H=ell*(1-H*H),theta=-F(3,2)*ell*H*theta)
    else:out.update(rho=-ell*rho*rho*a,a=-ell*rho*ay,ay=-ell*rho*P.var('ayy'),R=Rr,Rr=P.var('Rrr'),p=-ell*rho*py,py=-ell*rho*P.var('pyy'))
    return out


def normal(poly):
    """Use cosh(t)^(-2)=1-tanh(t)^2 after Laurent collection."""
    out=P();S=1-P.var('H')**2
    for k,v in P(poly).t.items():
        powers=dict(k);e=powers.get('c',0)
        if e<=-2:
            m=(-e)//2;powers['c']=e+2*m
            key=tuple(sorted((n,p) for n,p in powers.items() if p))
            out+=P({key:v})*S**m
        else:out+=P({k:v})
    return out


def authenticate(root):
    raw=(root/WAVE/'MANIFEST.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest()==MANIFEST_PIN,'FIRST_WAVE_MANIFEST')
    seen={'MANIFEST.json'}
    for row in json.loads(raw)['files']:
        rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'WAVE_MANIFEST_PATH')
        seen.add(str(rel));q=root/WAVE/rel
        need(q.is_file() and not q.is_symlink(),'WAVE_MEMBER_TYPE')
        b=q.read_bytes();need(len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],'WAVE_MEMBER:'+str(rel))
    need(all(not q.is_symlink() and (not q.is_file() or q.relative_to(root/WAVE).as_posix() in seen) for q in (root/WAVE).rglob('*')),'WAVE_CLOSED_TREE')
    for rel,pin in PINS.items():
        q=root/rel;need(q.is_file() and not q.is_symlink(),'SOURCE_TYPE:'+rel)
        need(hashlib.sha256(q.read_bytes()).hexdigest()==pin,'SOURCE_PIN:'+rel)
    return {'first_wave_manifest_sha256':MANIFEST_PIN,'first_wave_closed_files':len(seen),'files':dict(PINS)}


def function_ast(root,fn):
    path=root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py'
    tree=ast.parse(path.read_text());node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name==fn)
    env={'ELL':P.var('ell')};exec(compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),str(path),'exec'),env)
    return env[fn]


def implementation_matrices(root):
    rho,a,ay,H,ell,R,theta,Rr=map(P.var,('rho','a','ay','H','ell','R','theta','Rr'))
    path=root/WAVE/'operator/check_gauge_operator.py'
    tree=ast.parse(path.read_text());fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='build')
    targets={'A','B','D','force'}
    body=[x for x in fn.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id in targets for t in x.targets)]
    need(len(body)==4,'ORIGINAL_OPERATOR_ASSIGNMENTS')
    u=R*theta*field('W');uz=u.dr()/ell
    env={'rho':rho,'a':a,'ay':ay,'H':H,'lam':8,'S':1-H*H,'u':u,'uz':uz}
    exec(compile(ast.fix_missing_locations(ast.Module(body=body,type_ignores=[])),str(path),'exec'),env)
    old={k:env[k] for k in targets}
    path=root/WAVE/'coefficients/check_global_gauge_coefficients.py'
    tree=ast.parse(path.read_text());fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='run')
    start=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='a_factor' for t in x.targets))
    stop=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='P0' for t in x.targets))
    env={'ELL':ell,'hh':H,'aa':a,'rr':rho,'ay':ay,'ss':P.var('c')**-2,'zz':P(),'CC':R,'qq':theta,'Cz':Rr/ell,'injection':'none'}
    exec(compile(ast.fix_missing_locations(ast.Module(body=fn.body[start:stop+1],type_ignores=[])),str(path),'exec'),env)
    return old,{k:env[k] for k in ('A','B','D','Ps','P0')}


def angular_data():
    Y=P.var('Y');V=[P.var('Y'+str(i)) for i in range(3)]
    h=[[P.var('H'+str(min(i,j))+str(max(i,j))) for j in range(3)] for i in range(3)]
    h[2][2]=-8*Y-h[0][0]-h[1][1]
    def angular(poly,k,wrong=False):
        rules_={'Y':V[k]}
        for i in range(3):rules_['Y'+str(i)]=h[k][i]
        for i,j in ((0,0),(0,1),(0,2),(1,1),(1,2)):
            rules_['H'+str(i)+str(j)]=-(i==k)*V[j]-(j==k)*V[i]-(1 if wrong else 2)*(i==j)*V[k]
        return P(poly).derivation(rules_)
    return Y,V,h,angular


def audit(root,mutation=None):
    before=authenticate(root);ell,rho,a,ay,c,H,R,B,k,theta,p=map(P.var,('ell','rho','a','ay','c','H','R','Bcoef','k','theta','p'))
    W=field('W');f,b,d=[field(x) for x in ('f','b','d')]
    u=R*theta*W;v=B*theta*(W.dr()/ell+k*W)
    amplitudes=function_ast(root,'frame_amplitudes')
    injection={'wrong_scalar_shift':'wrong_scalar_shift','drop_radial_rho':'drop_radial_rho'}.get(mutation,'none')
    raw_f=rho*f if mutation=='raw_radial_normalization' else f
    amps,scalar=amplitudes(W,[raw_f,b,d],rho,p,a,R,B,k,c,H,theta,injection)
    if mutation=='swap_clock_term':amps['L']=2*(u+a*f+b.dx()/ell)
    Y,V,hess,angular=angular_data();sig=[1,-1,1,1,1]
    if mutation=='wrong_time_contraction':sig[1]=1
    q=[[P() for _ in range(5)] for _ in range(5)]
    q[0][0]=amps['A']*Y;q[0][1]=q[1][0]=amps['B']*Y;q[1][1]=amps['D']*Y
    for i in range(3):
        q[0][i+2]=q[i+2][0]=amps['C']*V[i]
        q[1][i+2]=q[i+2][1]=amps['E']*V[i]
        for j in range(3):q[i+2][j+2]=(amps['L']*Y if i==j else P())+amps['M']*hess[i][j]
    checks=[]
    def eq(name,left,right=0):
        residual=normal(P(left)-right)
        need(not residual,'NONZERO:'+name+':terms='+str(len(residual.t)))
        checks.append(name)
    # Independent coordinate Lie derivative at a sphere-normal point. Metric
    # values and first derivatives suffice; no retained pullback formula used.
    g=[P(1),-rho*rho]+[rho*rho*c*c]*3
    xi=[f*Y,b*Y]+[d*z for z in V]
    def coord(expr,i):return -P(expr).dr()/(ell*rho) if i==0 else P(expr).dx()/ell if i==1 else angular(expr,i-2)
    scale=[P(1),1/rho]+[1/(rho*c)]*3
    for i in range(5):
        for j in range(i,5):
            lie=sum((xi[z]*coord(g[i],z) for z in range(5)),P()) if i==j else P()
            lie+=g[j]*coord(xi[j],i)+g[i]*coord(xi[i],j)
            longitudinal=(-4*u if i==0 else 2*u)*g[i]*Y if i==j else P()
            eq('COORDINATE_LIE_FRAME_'+str(i)+str(j),q[i][j],(lie+longitudinal)*scale[i]*scale[j])
    eq('SCALAR_LIE_PULLBACK',scalar,v+p*f)
    # Exact variable-rho ordinary-derivative conversion, separately compiled
    # from implementation. This differentiates only a formal smooth function.
    physical=function_ast(root,'physical_derivatives')
    J=field('J');derivatives=physical(J,rho,a)
    if mutation=='omit_variable_radial_conversion':derivatives['yy']=J.dr().dr()/(ell*rho)**2
    Dy=lambda z:-P(z).dr()/(ell*rho)
    Dt=lambda z:P(z).dx()/ell
    for suffix,expected in [('',J),('y',Dy(J)),('t',Dt(J)),('yy',Dy(Dy(J))),('yt',Dy(Dt(J))),('tt',Dt(Dt(J)))]:
        eq('PHYSICAL_DERIVATIVE_'+suffix,derivatives[suffix],expected)
    # Full parallel-frame connection; intrinsic sphere connection is zero at
    # this arbitrary normal-frame point. The quadrupole's angular jets remain
    # free subject to their exact covariant trace/third-derivative identities.
    connection={}
    connection[1,0,1]=a;connection[1,1,0]=a
    for i in range(3):
        z=i+2;connection[z,0,z]=a;connection[z,1,z]=H/rho
        connection[z,z,0]=-a;connection[z,z,1]=H/rho
    def frame_derivative(expr,z):return Dy(expr) if z==0 else Dt(expr)/rho if z==1 else angular(expr,z-2,mutation=='wrong_hessian_divergence')/(rho*c)
    dq={}
    for z in range(5):
        for i in range(5):
            for j in range(5):
                dq[z,i,j]=frame_derivative(q[i][j],z)-sum((connection.get((z,i,s),P())*q[s][j]+connection.get((z,j,s),P())*q[i][s] for s in range(5)),P())
    C=[sum((sig[j]*(dq[j,i,j]-dq[i,j,j]/2) for j in range(5)),P()) for i in range(5)]
    old,new=implementation_matrices(root)
    if mutation=='wrong_radial_drift':new['B'][0][0]-=ell*rho*a
    if mutation=='wrong_tangential_drift':new['B'][1][1]+=2*ell*rho*a
    if mutation=='missing_ell_time_factor':new['A']=[[z/ell for z in row] for row in new['A']]
    for n,factor in (('A',ell),('B',ell),('D',ell*ell)):
        for i in range(3):
            for j in range(3):eq('MATRIX_'+n+'_'+str(i)+str(j),new[n][i][j],factor*old[n][i][j])
    force=[new['Ps']*W.dr()+new['P0']*W,P(),P()]
    if mutation=='wrong_forcing_sign':force[0]=-force[0]
    if mutation=='omit_source_ell_squared':force[0]/=ell**2
    if mutation=='omit_source_rho_squared':force[0]/=rho**2
    for i in range(3):eq('SOURCE_RECONSTRUCTION_'+str(i),force[i],ell*ell*old['force'][i])
    fields=[f,b,d]
    residual=[z.dx().dx()-z.dr().dr()-sum((new['A'][i][j]*fields[j].dx()+new['B'][i][j]*fields[j].dr()+new['D'][i][j]*fields[j] for j in range(3)),P()) for i,z in enumerate(fields)]
    eq('GAUGE_OPERATOR_RADIAL',ell*ell*rho*rho*C[0],(force[0]-residual[0])*Y)
    eq('GAUGE_OPERATOR_TIME',ell*ell*rho*C[1],residual[1]*Y)
    afactor=ell*ell*rho*c if mutation=='wrong_angular_frame_factor' else ell*ell*rho/c
    for i in range(3):eq('GAUGE_OPERATOR_ANGULAR_'+str(i),afactor*C[i+2],-residual[2]*V[i])
    # Same operator, zero gauge Cauchy functions. Spatial derivatives and
    # spatial derivatives of their zero velocities vanish as function jets.
    initial={'c':1,'H':0,'theta':1}
    for n in ('f','b','d'):
        for ix,ir in ((0,0),(1,0),(0,1),(0,2),(1,1)):initial[n+'_'+str(ix)+str(ir)]=0
    initial_source=force[0].sub(initial)
    initial['f_20']=P() if mutation=='zero_initial_acceleration' else initial_source
    initial['b_20']=initial['d_20']=P()
    uy=Dy(u).sub(initial);u0=u.sub(initial);ut=Dt(u).sub(initial)
    expected_values={'A':-4*u0,'B':P(),'C':P(),'D':-2*u0,'E':P(),'L':2*u0,'M':P()}
    expected_t={'A':-4*ut,'B':-6*rho*(uy+4*a*u0),'C':P(),'D':-2*ut,'E':P(),'L':2*ut,'M':P()}
    for name in amps:
        eq('INITIAL_VALUE_'+name,amps[name].sub(initial),expected_values[name])
        eq('INITIAL_TIME_JET_'+name,Dt(amps[name]).sub(initial),expected_t[name])
    eq('INITIAL_SCALAR_VALUE',scalar.sub(initial),v.sub(initial))
    eq('INITIAL_SCALAR_TIME',Dt(scalar).sub(initial),Dt(v).sub(initial))
    eq('INITIAL_RADIAL_ACCELERATION',initial_source/ell**2,-6*rho*rho*(uy+4*a*u0))
    # Physical shell traces: tangential b,d values are left free.
    shell={'f_00':0,'f_10':0,'b_01':0,'d_01':0}
    eq('SHELL_MIXED_TIME',amps['B'].sub(shell))
    eq('SHELL_MIXED_ANGULAR',amps['C'].sub(shell))
    eq('SHELL_SCALAR_PULLBACK',scalar.sub(shell),v)
    need(authenticate(root)==before,'SOURCE_READBACK_CHANGED')
    return {'status':'PASS_EXACT_FIRST_WAVE_IMPLEMENTATION_BINDING',
      'arithmetic':'Q Laurent polynomials, actual scaled jet derivations, cosh^(-2)=1-tanh^2 quotient; no point samples',
      'source_authentication':before,'checks':checks,'named_checks':len(checks),
      'coordinate_Lie_metric_components':15,'source_registered_frame_sectors':7,
      'independent_DC_method':'full5Dparallel-frameconnectioncontraction;quadrupoleangularjetsandcovariantthirdderivative',
      'normalization_identities':['ell^2 rho^2 DC_y=(source_f-residual_f)Y','ell^2 rho DC_0=residual_b Y','ell^2 rho/c DC_a=-residual_d D_aY'],
      'initial_specialization':'same audited operator fixes f_tt=-6rho^2(u_y+4au) and b_tt=d_tt=0; compiled reconstruction gives corrected mixed first time jet',
      'source_roles':{'metric_forcing':'notusedhere','QC_evolving':'notcomputedhere','actual_wave_PDE':'analytic premise; this finite audit binds its operator coefficients'},
      'retained_prerequisites':['actual smooth first gauge wave solves the registered PDE with zero Cauchy data','retained actual background and coupled first reconstruction','physical-half-line boundary trace regularity','existing scalar/Bianchi/Codazzi/subsidiary uniqueness chain'],
      'not_proved':['new PDE solve','global nonlinear family','small physical errors','cosmological origin','external specialist approval']}


CONTROLS={
 'wrong_scalar_shift':'SCALAR_LIE_PULLBACK',
 'drop_radial_rho':'COORDINATE_LIE_FRAME_00',
 'raw_radial_normalization':'COORDINATE_LIE_FRAME_00',
 'swap_clock_term':'COORDINATE_LIE_FRAME_22',
 'wrong_time_contraction':'GAUGE_OPERATOR_RADIAL',
 'wrong_hessian_divergence':'GAUGE_OPERATOR_ANGULAR_0',
 'wrong_radial_drift':'MATRIX_B_00',
 'wrong_tangential_drift':'MATRIX_B_11',
 'missing_ell_time_factor':'MATRIX_A_00',
 'wrong_forcing_sign':'SOURCE_RECONSTRUCTION_0',
 'omit_source_ell_squared':'SOURCE_RECONSTRUCTION_0',
 'omit_source_rho_squared':'SOURCE_RECONSTRUCTION_0',
 'wrong_angular_frame_factor':'GAUGE_OPERATOR_ANGULAR_0',
 'omit_variable_radial_conversion':'PHYSICAL_DERIVATIVE_yy',
 'zero_initial_acceleration':'INITIAL_TIME_JET_B',
}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--mutation',choices=CONTROLS);ap.add_argument('--run-controls',action='store_true');a=ap.parse_args()
    need(not a.output.exists(),'REFUSE_OVERWRITE')
    report=audit(a.repo_root.resolve(),a.mutation)
    if a.run_controls:
        rejected=[]
        with tempfile.TemporaryDirectory() as tmp:
            for name,reason in CONTROLS.items():
                target=Path(tmp)/(name+'.json');cmd=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])+[str(Path(__file__).resolve()),'--repo-root',str(a.repo_root.resolve()),'--output',str(target),'--mutation',name]
                r=subprocess.run(cmd,capture_output=True,text=True)
                need(r.returncode!=0 and 'NONZERO:'+reason+':' in r.stderr and not target.exists(),'CONTROL_NOT_REJECTED:'+name+':'+r.stderr[-500:])
                rejected.append({'mutation':name,'expected_identity':reason,'returncode':r.returncode,'report_absent':True})
        report['rejected_corrupted_executions']=rejected
    report['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'checks':report['named_checks'],'controls':len(report.get('rejected_corrupted_executions',[]))}))


if __name__=='__main__':main()
