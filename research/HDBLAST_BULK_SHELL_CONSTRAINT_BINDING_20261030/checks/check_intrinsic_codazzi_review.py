#!/usr/bin/env python3
"""Independent exact coordinate-FRW divergence and natural Codazzi map review.

No scientific producer or arithmetic module is imported. Angular jets are
formal, with only trace(Hess Y)=-8Y and div(Hess Y)=-6dY imposed.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


class P:
    def __init__(self, value=0):
        if isinstance(value, P): self.t = dict(value.t)
        elif isinstance(value, dict): self.t = {k:F(v) for k,v in value.items() if v}
        else: self.t = {():F(value)} if value else {}
    @staticmethod
    def var(name): return P({((name,1),):F(1)})
    def __add__(self,b):
        out=dict(self.t)
        for k,v in P(b).t.items(): out[k]=out.get(k,F(0))+v
        return P(out)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.t.items()})
    def __sub__(self,b):return self+-P(b)
    def __rsub__(self,b):return P(b)+-self
    def __mul__(self,b):
        out={}
        for k,v in self.t.items():
            for l,w in P(b).t.items():
                e=dict(k)
                for n,p in l:e[n]=e.get(n,0)+p
                key=tuple(sorted((n,p) for n,p in e.items() if p))
                out[key]=out.get(key,F(0))+v*w
        return P(out)
    __rmul__=__mul__
    def __pow__(self,n):
        if not isinstance(n,int):raise ValueError('INTEGER_POWER')
        if n<0:
            if len(self.t)!=1:raise ValueError('NONMONOMIAL_INVERSE')
            (k,v),=self.t.items()
            return P({tuple((s,p*n) for s,p in k):v**n})
        r=P(1)
        for _ in range(n):r=r*self
        return r
    def __truediv__(self,b):return self*P(b)**-1
    def __rtruediv__(self,b):return P(b)*self**-1
    def __bool__(self):return bool(self.t)
    def __eq__(self,b):return self.t==P(b).t
    def diff(self,name):
        out={}
        for k,v in self.t.items():
            e=dict(k);n=e.get(name,0)
            if n:
                e[name]=n-1;key=tuple(sorted((s,p) for s,p in e.items() if p))
                out[key]=v*n
        return P(out)
    def names(self):return set(n for k in self.t for n,p in k)
    def substitute(self,table):
        out=P()
        for k,v in self.t.items():
            term=P(v)
            for n,p in k:term=term*table.get(n,P.var(n))**p
            out=out+term
        return out
    def angle_coefficient(self,angle_names,wanted):
        wanted=tuple(sorted((n,p) for n,p in wanted.items() if p))
        out={}
        for k,v in self.t.items():
            ang=tuple((n,p) for n,p in k if n in angle_names)
            if ang==wanted:
                key=tuple((n,p) for n,p in k if n not in angle_names)
                out[key]=out.get(key,F(0))+v
        return P(out)
    def encoded(self):return [[[[n,p] for n,p in k],str(v)] for k,v in sorted(self.t.items())]


PINS={
 'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/boundary-review/BOUNDARY_COMPATIBILITY_AUDIT.md':'581bcadecd44c408484d20632d7d3e524fb0bd098370d1fd24da2744bd6c1cf1',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/geometry/GEOMETRIC_SOURCE_DERIVATION.md':'599c59092ac15c0a375b758e6462fe73aa7969e11965fcef4de6c1b74c60bfe0',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json':'34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/source_adapters.py':'0233ff25b9873dfc4060ba0d4963f50d0b198dbe71368625017840d225b19f8a',
}


def need(ok,message):
    if not ok:raise ValueError(message)


def auth(root):
    for rel,pin in PINS.items():
        path=root/rel
        need(path.is_file() and not path.is_symlink(),'SOURCE_TYPE:'+rel)
        need(hashlib.sha256(path.read_bytes()).hexdigest()==pin,'SOURCE_SHA:'+rel)


def tangent_dt(expr):
    expr=P(expr);out=P()
    for n in expr.names():
        if n in ('h','X','V','L','d','v'):out=out+expr.diff(n)*P.var(n+'t')
        elif n in ('rho','c','theta','a','p','sig2','sig3'):
            if n=='c':out=out+expr.diff(n)*P.var('c')*P.var('theta')
            if n=='theta':out=out+expr.diff(n)*(1-P.var('theta')**2)
        else:raise ValueError('UNEXPECTED_TANGENT_DERIVATIVE:'+n)
    return out


def coordinate_divergence(mutation):
    rho,c,theta,Y,q00,q0,qa,qb,q00t,q0t,qat,qbt=[P.var(n) for n in
        ('rho','c','theta','Y','q00','q0','qa','qb','q00t','q0t','qat','qbt')]
    yi=[P.var('Y'+str(i)) for i in range(1,4)]
    h=[[P.var('H'+str(min(i,j))+str(max(i,j))) for j in range(1,4)] for i in range(1,4)]
    h[2][2]=-8*Y-h[0][0]-h[1][1]
    t=[[[P.var('T'+str(k)+str(min(i,j))+str(max(i,j))) for j in range(1,4)] for i in range(1,4)] for k in range(1,4)]
    # Exact sphere identities: D_i tr(HessY)=-8Y_i and div(HessY)=-6dY.
    for k in range(3):t[k][2][2]=-8*yi[k]-t[k][0][0]-t[k][1][1]
    t[2][2][0]=t[2][0][2]=-6*yi[0]-t[0][0][0]-t[1][1][0]
    t[2][2][1]=t[2][1][2]=-6*yi[1]-t[0][0][1]-t[1][1][1]
    t[1][1][2]=t[1][2][1]=2*yi[2]+t[2][0][0]+t[2][1][1]-t[0][0][2]
    inv=[-1/rho**2]+[1/(rho*c)**2]*3
    metric=[-rho**2]+[(rho*c)**2]*3
    dg=[[[P() for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for i in range(1,4):dg[0][i][i]=2*(rho*c)**2*theta
    Gamma=[[[inv[a]*(dg[b][a][cc]+dg[cc][a][b]-dg[a][b][cc])/2
             for cc in range(4)] for b in range(4)] for a in range(4)]
    qr=[[P() for _ in range(4)] for _ in range(4)]
    qr[0][0]=rho**2*q00*Y**2
    for i in range(3):
        qr[0][i+1]=qr[i+1][0]=rho**2*c*q0*Y*yi[i]
        for j in range(3):qr[i+1][j+1]=(rho*c)**2*(qa*Y**2*int(i==j)+qb*Y*h[i][j])
    dqr=[[[P() for _ in range(4)] for _ in range(4)] for _ in range(4)]
    dqr[0][0][0]=rho**2*q00t*Y**2
    for i in range(3):
        dqr[0][0][i+1]=dqr[0][i+1][0]=rho**2*c*(q0t+theta*q0)*Y*yi[i]
        for j in range(3):dqr[0][i+1][j+1]=(rho*c)**2*((qat+2*theta*qa)*Y**2*int(i==j)+(qbt+2*theta*qb)*Y*h[i][j])
    for k in range(3):
        dqr[k+1][0][0]=2*rho**2*q00*Y*yi[k]
        for i in range(3):
            dqr[k+1][0][i+1]=dqr[k+1][i+1][0]=rho**2*c*q0*(yi[k]*yi[i]+Y*h[k][i])
            for j in range(3):dqr[k+1][i+1][j+1]=(rho*c)**2*(2*qa*Y*yi[k]*int(i==j)+qb*(yi[k]*h[i][j]+Y*t[k][i][j]))
    rows=[]
    for a in range(4):
        div=P()
        for b in range(4):
            cov=dqr[b][b][a]
            if mutation!='omit_shell_connection':
                for dd in range(4):cov=cov-Gamma[dd][b][b]*qr[dd][a]-Gamma[dd][b][a]*qr[b][dd]
            sg=-inv[b] if mutation=='positive_temporal_divergence' and b==0 else inv[b]
            div=div+sg*cov
        trder=sum((inv[b]*dqr[a][b][b] for b in range(4)),P())
        if a==0:trder=trder+sum((-2*theta*inv[b]*qr[b][b] for b in range(1,4)),P())
        rows.append((div-trder)/(rho if a==0 else rho*c))
    angle=set(n for row in rows for n in row.names() if n.startswith(('Y','H','T')))
    return rows,{'Y':Y,'yi':yi,'H':h,'angles':angle,'metric':metric,'inverse':inv,'Gamma':Gamma}


def evaluate_encoded(terms,values):
    # Separate evaluator from the implementation's projection/source_coefficients.py.
    result=P()
    for term in terms:
        item=P(F(term['coefficient']))
        for name,power in sorted(term['powers'].items()):item=item*values[name]**int(power)
        result=result+item
    return result


def calculate(root,mutation=None):
    auth(root)
    checks=[]
    def zero(name,value):
        value=P(value)
        need(not value,'NONZERO:'+name+':'+json.dumps(value.encoded()[:3],separators=(',',':')))
        checks.append(name)
    rows,info=coordinate_divergence(mutation)
    rho,c,theta,q00,q0,qa,qb,qat,qbt,q0t=[P.var(n) for n in
        ('rho','c','theta','q00','q0','qa','qb','qat','qbt','q0t')]
    Y=info['Y'];yi=info['yi'];h=info['H'];GY=sum((x*x for x in yi),P())
    T=3*qa-8*qb
    generic_scalar=(-3*qat+8*qbt-theta*(3*q00+T)-8*q0/c)/rho
    generic_G=q0/(rho*c)
    generic_vector=(2*q00-4*qa+10*qb)/(rho*c)-q0t/rho-4*theta*q0/rho
    generic_H=qb/(rho*c)
    zero('INTRINSIC_TEMPORAL_FROM_CONNECTION',rows[0]-generic_scalar*Y**2-generic_G*GY)
    for i in range(3):zero('INTRINSIC_ANGULAR_FROM_CONNECTION_'+str(i),rows[i+1]-generic_vector*Y*yi[i]-generic_H*sum((h[i][j]*yi[j] for j in range(3)),P()))
    # Extract the primitive coefficients from the independently contracted
    # tensor, rather than using the compact expressions to produce QE targets.
    scalar=rows[0].angle_coefficient(info['angles'],{'Y':2})
    grad=rows[0].angle_coefficient(info['angles'],{'Y1':2})
    vec=rows[1].angle_coefficient(info['angles'],{'Y':1,'Y1':1})
    hes=rows[1].angle_coefficient(info['angles'],{'H12':1,'Y2':1})
    zero('SCALAR_PRIMITIVE_EXTRACTION',rows[0]-scalar*Y**2-grad*GY)
    zero('VECTOR_PRIMITIVE_EXTRACTION',rows[1]-vec*Y*yi[0]-hes*sum((h[0][j]*yi[j] for j in range(3)),P()))
    # Independent original physical shell QR, not imported adapter formulas.
    a,p,hv,X,V,L,d,v,sig2=[P.var(n) for n in ('a','p','h','X','V','L','d','v','sig2')]
    sig1=-2*p
    if mutation=='wrong_scalar_shell_sign':sig1=2*p
    beta=a*hv+sig1*v/3
    A0=a*hv**2/4-sig1*hv*v/6-sig2*v*v/6
    qprofiles={'q00':-A0+2*beta*X,'q0':-beta*V,'qa':A0-2*beta*L,'qb':-2*beta*d}
    if mutation=='omit_anisotropic_shell_source':qprofiles['qb']=P()
    qprofiles.update({n+'t':tangent_dt(qprofiles[n]) for n in ('q00','q0','qa','qb')})
    scalar,grad,vec,hes=[x.substitute(qprofiles) for x in (scalar,grad,vec,hes)]
    scalars=[scalar/6+F(4,3)*grad,scalar/2+2*grad,scalar-4*grad]
    vectors=[vec/4+hes,vec/2-2*hes]
    report=json.loads((root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json').read_text())
    values={n:P.var('free_'+n) for n in report['input_symbols']}
    values.update({n:P.var(n) for n in ('rho','c','theta','a','ay','p','py','Vp','Vpp','Vppp','ell')})
    amplitudes={'A':hv,'D':-2*X,'E':V,'L':2*L,'M':2*d,'s':v}
    for n,x in amplitudes.items():
        values[n]=x;values[n+'t']=tangent_dt(x);values[n+'tt']=P.var('free_'+n+'tt')
    for n in ('B','C'):
        for suffix in ('','t','tt'):values[n+suffix]=P()
    values.update({'Dy':-beta,'Ly':beta,'Ey':P(),'My':P(),
                   'Dyt':-tangent_dt(beta),'Lyt':tangent_dt(beta),'Eyt':P(),'Myt':P(),
                   'sy':p*hv/2-sig2*v/2,'syt':p*P.var('ht')/2-sig2*P.var('vt')/2})
    # First normal gauge Ay and every unused/free normal jet remain FREE.
    if mutation=='coordinate_instead_of_frame_normal':values['Dy']=values['Dy']-2*a*values['D']
    family='QEred' if mutation=='reduced_instead_of_natural' else 'QE'
    produced={field:[evaluate_encoded(row,values) for row in report['maps'][family][field]] for field in ('y0','ya')}
    if mutation=='half_second_variation':produced={k:[x/2 for x in vv] for k,vv in produced.items()}
    for degree,actual,want in zip((0,2,4),produced['y0'],scalars):zero('NATURAL_CODAZZI_TIME_'+str(degree),actual-want)
    for degree,actual,want in zip((2,4),produced['ya'],vectors):zero('NATURAL_CODAZZI_ANGULAR_'+str(degree),actual-want)
    # Independent check of the same five exported rows BEFORE P1=0. This
    # keeps the actual scalar-normal profile free, not an arbitrary interval.
    offvalues=dict(values);offvalues['sy']=P.var('free_sy');offvalues['syt']=P.var('free_syt')
    Pfirst=offvalues['sy']-p*hv/2+sig2*v/2
    offproduced={field:[evaluate_encoded(row,offvalues) for row in report['maps']['QE'][field]] for field in ('y0','ya')}
    for degree,actual,want,coef in zip((0,2,4),offproduced['y0'],scalars,(F(1,6),F(1,2),F(1))):
        zero('OFFSHELL_MAP_CODAZZI_TIME_'+str(degree),actual-want+2*Pfirst*P.var('vt')/rho*coef)
    for degree,actual,want,coef in zip((2,4),offproduced['ya'],vectors,(F(1,4),F(1,2))):
        zero('OFFSHELL_MAP_CODAZZI_ANGULAR_'+str(degree),actual-want+2*Pfirst*v/(rho*c)*coef)
    for field in ('y0','ya'):
        for j,actual in enumerate(produced[field]):
            need(not any(n.startswith('free_') for n in actual.names()),'FREE_NORMAL_JET_SURVIVED:'+field+str(j))
    # Exact off-shell matter cancellation, before imposing P1=0. This uses
    # physical K=R+sigma*gamma/6 and nphi=-sigma'/2+P on a4d shell.
    divR,dtrR,sp,Pval,dphi=[P.var(n) for n in ('divR','dtrR','sp','Pval','dphi')]
    geometric=divR-dtrR-sp*dphi/2
    natural=geometric-(-sp/2+Pval)*dphi
    if mutation=='omit_scalar_codazzi_cancellation':natural=geometric
    zero('OFFSHELL_SCALAR_CODAZZI',natural-(divR-dtrR-Pval*dphi))
    # Second variation of P*dphi: background dphi0=0 and P0=P1=0.
    P0,P1,P2,d0,d1,d2=[P.var(n) for n in ('P0','P1','P2','d0','d1','d2')]
    second=P2*d0+2*P1*d1+P0*d2
    zero('SECOND_SCALAR_CROSS_VANISHES',second.substitute({'P0':P(),'P1':P(),'d0':P()}))
    # Derive the raised unit-normal rescaling and its Hessian cross terms in
    # ordinary epsilon coefficients. The affine boundary has g_ya=0.
    def mulseries(xx,yy):return [sum((xx[j]*yy[k-j] for j in range(k+1)),P()) for k in range(3)]
    invN=[P(1),-hv/2,F(3,8)*hv**2]
    normalized=mulseries(mulseries(invN,invN),[P(1),hv,P()])
    for k,row in enumerate(normalized):zero('RAISED_UNIT_NORMAL_SERIES_'+str(k),row-(1 if k==0 else 0))
    Ezero,Efirst,Ehessian=[P.var(n) for n in ('E0na','E1na','QEna')]
    second=2*mulseries(invN,[Ezero,Efirst,Ehessian/2])[2]
    zero('NORMAL_PULLBACK_SECOND_CROSS_TERMS',second-(Ehessian-hv*Efirst+F(3,4)*hv**2*Ezero))
    zero('NORMAL_VARIATION_LOWER_RESIDUALS_VANISH',second.substitute({'E0na':P(),'E1na':P()})-Ehessian)
    auth(root)
    return {'status':'PASS_INDEPENDENT_COORDINATE_CODAZZI_REVIEW',
      'arithmetic':'own exact rational Laurent ring; no producer or scientific code imports',
      'source_SHA256':PINS,'checks':checks,'check_count':len(checks),
      'coordinate_tensor_dimension':4,'angular_identity_premises':['trace HessY=-8Y','div HessY=-6dY'],
      'natural_Codazzi_modal_rows':5,'off_shell_natural_map_rows_with_P1_free':5,'normal_gauge_Ay_used':False,
      'first_normal_jets_left_free':['Ay','By','Cy','Ayt','Byt','Cyt','all unused yy jets'],
      'bulk_PDE_or_harmonic_identity_used':False,
      'independent_scalar_projected_divergence':[x.encoded() for x in scalars],
      'independent_gradient_projected_divergence':[x.encoded() for x in vectors],
      'scope':'All-point natural source binding on coherent first physical shell traces. Existing geometric Codazzi, smooth trace and actual first-boundary matching proofs are retained. No subsidiary uniqueness, finite-amplitude nonlinear or cosmological realization.'}


CONTROLS={
 'positive_temporal_divergence':'INTRINSIC_TEMPORAL_FROM_CONNECTION',
 'omit_shell_connection':'INTRINSIC_TEMPORAL_FROM_CONNECTION',
 'wrong_scalar_shell_sign':'NATURAL_CODAZZI_TIME_0',
 'omit_anisotropic_shell_source':'NATURAL_CODAZZI_TIME_0',
 'coordinate_instead_of_frame_normal':'NATURAL_CODAZZI_TIME_0',
 'reduced_instead_of_natural':'NATURAL_CODAZZI_TIME_0',
 'half_second_variation':'NATURAL_CODAZZI_TIME_0',
 'omit_scalar_codazzi_cancellation':'OFFSHELL_SCALAR_CODAZZI',
}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--mutation',choices=CONTROLS);ap.add_argument('--run-controls',action='store_true');args=ap.parse_args()
    need(not args.output.exists(),'REFUSE_OVERWRITE')
    result=calculate(args.repo_root.resolve(),args.mutation)
    if args.run_controls:
        controls=[]
        with tempfile.TemporaryDirectory() as folder:
            for mutation,reason in CONTROLS.items():
                output=Path(folder)/(mutation+'.json')
                command=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])+[str(Path(__file__).resolve()),'--repo-root',str(args.repo_root.resolve()),'--output',str(output),'--mutation',mutation]
                run=subprocess.run(command,capture_output=True,text=True)
                need(run.returncode!=0 and 'NONZERO:'+reason+':' in run.stderr and not output.exists(),'CONTROL_NOT_REJECTED:'+mutation+':'+run.stderr[-1500:])
                controls.append({'mutation':mutation,'expected_nonzero_identity':reason,'returncode':run.returncode,'successful_report_absent':True})
        result['rejected_corrupted_executions']=controls
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'check_count':result['check_count'],'controls':len(result.get('rejected_corrupted_executions',[]))}))


if __name__=='__main__':main()
