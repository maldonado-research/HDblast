#!/usr/bin/env python3
"""Exact coherent ordinary-jet witnesses for insufficient pointwise premises.

No parent geometry checker is executed. Covariant derivatives use a separately
constructed upper connection table and apply it to every tensor slot.
"""
import argparse
import ast
from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import product
import json
from math import comb
from pathlib import Path
import sys
sys.dont_write_bytecode = True
CAPSULE = 'research/HDBLAST_NATURAL_FIELD_BINDING_20261030'
MANIFEST_PIN = '5121e4304fad7ab5856977b38c6ebd9fdc205900b543bc1c3d97384d246f378f'
MAP = 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
MAP_PIN = '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
HELPER_PINS = {'codazzi_laurent.py':'178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2','angular_generic.py':'fd14d8cc3ac6f58c206b784c5c57d15430f63f38ca091ac03c7769f9bfc89df9'}

def need(ok, label):
    if not ok: raise ValueError(label)

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def authenticate(root):
    capsule = root/CAPSULE
    need(digest(capsule/'MANIFEST.json') == MANIFEST_PIN, 'NATURAL_MANIFEST_PIN')
    manifest = json.loads((capsule/'MANIFEST.json').read_text())
    files = {}
    for row in manifest['files']:
        path = capsule/row['path']
        need(path.is_file() and not path.is_symlink(), 'CAPSULE_FILE_TYPE:'+row['path'])
        need(path.stat().st_size == row['bytes'] and digest(path) == row['sha256'], 'CAPSULE_FILE_PIN:'+row['path'])
        files[row['path']] = row['sha256']
    need(digest(root/MAP) == MAP_PIN, 'GEOMETRIC_MAP_PIN')
    for name,pin in HELPER_PINS.items():
        need(digest(capsule/'checks'/name) == pin, 'HELPER_PIN:'+name)
    allowed = {'argparse','ast','fractions','hashlib','importlib','itertools','json','math','pathlib','sys','functools'}
    for path in (Path(__file__).resolve(), *(capsule/'checks'/n for n in HELPER_PINS)):
        tree = ast.parse(path.read_text())
        for n in ast.walk(tree):
            need(not isinstance(n,ast.Assert), 'ASSERT_FORBIDDEN:'+path.name)
            if isinstance(n,ast.Import): need(all(x.name.split('.')[0] in allowed for x in n.names), 'IMPORT_FORBIDDEN:'+path.name)
            if isinstance(n,ast.ImportFrom): need(n.module and n.module.split('.')[0] in allowed, 'IMPORT_FORBIDDEN:'+path.name)
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Name): need(n.func.id not in {'eval','exec','compile','__import__'}, 'DYNAMIC_CODE_FORBIDDEN:'+path.name)
    return {'natural_manifest_sha256':MANIFEST_PIN,'natural_files':files,'map_sha256':MAP_PIN,'helpers':HELPER_PINS,'checker_sha256':digest(Path(__file__))}

def load(root,name):
    path = root/CAPSULE/'checks'/(name+'.py')
    spec = importlib.util.spec_from_file_location('witness_'+name,path)
    mod = importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

def run(root, mutation=None):
    before = authenticate(root)
    alg = load(root,'codazzi_laurent'); sphere = load(root,'angular_generic'); P = alg.P
    var = P.var
    rho,a,p,c,theta,Vp,Vpp,Vppp = [var(n) for n in ('rho','a','p','c','theta','Vp','Vpp','Vppp')]
    ay = -rho**-2-p*p/3; py = Vp-4*a*p; V = 6*rho**-2+p*p/2-6*a*a
    def dy(x):
        rules = {'rho':a*rho,'a':ay,'p':py,'Vp':Vpp*p,'Vpp':Vppp*p}
        for lead in ('A','d','f'):
            for k in range(5): rules[lead+'y'*k] = var(lead+'y'*(k+1))
        return sum((x.diff(n)*rules.get(n,P()) for n in x.names()),P())
    def dt(x):
        rules = {'c':c*theta,'theta':1-theta*theta}
        return sum((x.diff(n)*rules.get(n,P()) for n in x.names()),P())
    def quotient(x):
        out = P()
        for mon,coef in x.t.items():
            z = dict(mon);power = z.pop('theta',0)
            need(power>=0,'NEGATIVE_THETA_POWER')
            for k in range(power//2+1):
                w = dict(z);w['theta']=power%2;w['c']=w.get('c',0)-2*k
                out += P({tuple(sorted((n,e) for n,e in w.items() if e)):coef*comb(power//2,k)*(-1)**k})
        return out
    class AP:
        def __init__(self,v=0): self.p = sphere.clean(v) if isinstance(v,dict) else ({sphere.ZERO:P.cv(v)} if v else {})
        @classmethod
        def cv(cls,x):return x if isinstance(x,cls) else cls(x)
        def __add__(self,x):return AP(sphere.add(self.p,self.cv(x).p))
        __radd__=__add__
        def __neg__(self):return AP(sphere.scale(-1,self.p))
        def __sub__(self,x):return self+-self.cv(x)
        def __rsub__(self,x):return self.cv(x)+-self
        def __mul__(self,x):return AP(sphere.mul(self.p,self.cv(x).p))
        __rmul__=__mul__
        def ordinary(self,d):return AP({e:(dy(v) if d==0 else dt(v)) for e,v in self.p.items()})
        def angular(self,d):return AP(sphere.diff(self.p,d))
    I=range(6);J=range(4);sign=(1,-1,1,1,1,1)
    K={ij:AP(v) for ij,v in sphere.P.items()}; X=[AP(v) for v in sphere.X];Y=AP(sphere.Y)
    dY={i:AP(v) for i,v in sphere.grad(sphere.Y).items()}
    gamma={(u,d,v):AP() for u,d,v in product(I,repeat=3)}
    gamma[1,1,0]=AP(a); gamma[0,1,1]=AP(a)
    for i,k in product(J,repeat=2):
        gamma[k+2,i+2,0]=K[i,k]*a;gamma[0,i+2,k+2]=-K[i,k]*a
        gamma[k+2,i+2,1]=K[i,k]*(theta/rho);gamma[1,i+2,k+2]=K[i,k]*(theta/rho)
    def cov(t,idx,d):
        val=t[idx]
        if d==0: out=val.ordinary(0)
        elif d==1: out=val.ordinary(1)*(1/rho)
        else:
            out=sum((K[d-2,k]*val.angular(k) for k in J),AP())
            for s,v in enumerate(idx):
                if v>=2:
                    repl=list(idx);repl[s]=d;out+=X[v-2]*t[tuple(repl)]
            out=out*(1/(rho*c))
        for s,v in enumerate(idx):
            for u in I:
                repl=list(idx);repl[s]=u;out-=gamma[u,d,v]*t[tuple(repl)]
        return out
    point={'rho':Q(1),'c':Q(1),'a':Q(0),'theta':Q(0),'p':Q(1),'Vp':Q(0),'Vpp':Q(0),'Vppp':Q(0),'A':Q(0),'Ay':Q(1),'Ayy':Q(1),'Ayyy':Q(0),'f':Q(0),'fy':Q(1),'fyy':Q(0),'fyyy':Q(0),'fyyyy':Q(0),'fyyyyy':Q(0),'d':Q(0),'dy':Q(0),'dyy':Q(0),'dyyy':Q(1),'dyyyy':Q(0),'dyyyyy':Q(0)}
    def ev(x):
        total=Q(0)
        for mon,coef in quotient(P.cv(x)).t.items():
            v=coef
            for n,e in mon:need(n in point,'UNASSIGNED_JET:'+n);v*=point[n]**e
            total+=v
        return total
    def evAP(x):return sphere.clean({e:ev(v) for e,v in x.p.items()})
    def expectAP(label,left,right=0):
        diff=evAP(left-AP.cv(right));need(not diff,'NONZERO:'+label+':'+json.dumps([[list(k),str(v)] for k,v in diff.items()],separators=(',',':')))
    def encodeAP(x):return [[list(e),str(v)] for e,v in sorted(evAP(x).items())]
    def wave(z):
        dz={(d,b):cov(z,(b,),d) for d,b in product(I,repeat=2)}
        rough={b:sum((sign[k]*cov(dz,(k,b),k) for k in I),AP()) for b in I}
        return {b:rough[b]+z[b,]*(Q(2,3)*V)+(z[b,]*p*p if b==0 else AP()) for b in I}
    q={(i,j):AP() for i,j in product(I,repeat=2)};q[0,0]=Y*var('A')
    dq={(d,i,j):cov(q,(i,j),d) for d,i,j in product(I,repeat=3)}
    connection={(k,i,j):(dq[i,j,k]+dq[j,i,k]-dq[k,i,j])*Q(1,2) for k,i,j in product(I,repeat=3)}
    C={b:sum((sign[k]*connection[b,k,k] for k in I),AP()) for b in I}
    diff={(b,):2*sum((sign[k]*q[b,k]*C[k] for k in I),AP()) for b in I}
    for b in I:expectAP('LOWERING_VALUE:'+str(b),diff[b,])
    lowerwave=wave(diff)
    for b in I:expectAP('LOWERING_WAVE:'+str(b),lowerwave[b],3*Y*Y if b==0 else AP())
    need(bool(evAP(lowerwave[0])),'LOWERING_WAVE_WITNESS_ZERO')
    # The source maps are substituted before differentiating, preserving the
    # actual ordinary-jet relation of their composed coefficients.
    report=json.loads((root/MAP).read_text());vals={n:P() for n in report['input_symbols']}
    vals.update({n:var(n) for n in ('rho','c','a','theta','p','Vp','Vpp','Vppp','A','Ay','Ayy')});vals.update(ay=ay,py=py,ell=P(1))
    maps={name:{field:[alg.encoded(rows,vals) for rows in terms] for field,terms in report['maps'][name].items()} for name in ('QC_actual','QC_background')}
    mapdiff={field:[a-b for a,b in zip(maps['QC_actual'][field],maps['QC_background'][field])] for field in ('y','0','a')}
    mapped={(b,):AP() for b in I}
    for field,b in (('y',0),('0',1)):
        mapped[b,]=sum((AP(base[()])*cf for base,cf in zip(sphere.SCALARS,mapdiff[field])),AP())
    for i in J:mapped[i+2,]=sum((AP(base[i,])*cf for base,cf in zip(sphere.GRADS,mapdiff['a'])),AP())
    for b in I:
        residual=mapped[b,]-diff[b,]
        need(not residual.p,'MAP_LOWERING_SYMBOLIC:'+str(b))
    mapwave=wave(mapped)
    for b in I:expectAP('MAP_LOWERING_WAVE:'+str(b),mapwave[b],lowerwave[b])
    if mutation=='exchange_QC_under_wave':
        for b in I:expectAP('INVALID_POINTWISE_QC_EXCHANGE:'+str(b),mapwave[b])
    # q=L_{f(y)Y partial_y}g0, eta=L_xi phi0, plus d(y)Y in q_00.
    # The added d has zero jets through order two and d_yyy=1 at the point.
    f=var('f');d=var('d')
    q={(i,j):AP() for i,j in product(I,repeat=2)}
    q[0,0]=Y*(2*var('fy'));q[1,1]=Y*(-2*a*f+d)
    for i in J:
        q[0,i+2]=q[i+2,0]=dY[i,]*(f/(rho*c))
        for j in J:q[i+2,j+2]=Y*K[i,j]*(2*a*f)
    g={(i,j):AP() for i,j in product(I,repeat=2)};g[0,0]=AP(1);g[1,1]=AP(-1)
    for i,j in product(J,repeat=2):g[i+2,j+2]=K[i,j]
    dq={(k,i,j):cov(q,(i,j),k) for k,i,j in product(I,repeat=3)}
    connection={(k,i,j):(dq[i,j,k]+dq[j,i,k]-dq[k,i,j])*Q(1,2) for k,i,j in product(I,repeat=3)}
    ric={(i,j):sum((sign[k]*(cov(connection,(k,i,j),k)-cov(connection,(k,i,k),j)) for k in I),AP()) for i,j in product(I,repeat=2)}
    eta=Y*p*f;de={(i,):cov({():eta},(),i) for i in I}
    E={(i,j):ric[i,j]-(de[j,]*p if i==0 else AP())-(de[i,]*p if j==0 else AP())-Q(2,3)*(eta*Vp*g[i,j]+q[i,j]*V) for i,j in product(I,repeat=2)}
    dE={(k,i,j):cov(E,(i,j),k) for k,i,j in product(I,repeat=3)}
    for i,j in product(I,repeat=2):expectAP('POINT_E1_ZERO:'+str((i,j)),E[i,j])
    expectAP('PROLONGATION_YY',dE[0,0,0],Y*Q(1,2))
    expectAP('PROLONGATION_00',dE[0,1,1],Y*Q(-1,2))
    C={b:sum((sign[k]*connection[b,k,k] for k in I),AP()) for b in I}
    delta={b:sum((sign[i]*sign[j]*q[i,j]*(-dE[i,j,b]+Q(1,2)*dE[b,i,j]) for i,j in product(I,repeat=2)),AP())-sum((sign[k]*C[k]*E[k,b] for k in I),AP()) for b in I}
    for b in I:expectAP('POINT_DELTA_B:'+str(b),delta[b],-Y*Y*Q(1,2) if b==0 else AP())
    if mutation=='point_E1_implies_prolongation':
        expectAP('INVALID_POINT_E1_PROLONGATION',2*delta[0])
    need(authenticate(root)==before,'AUTHENTICATED_INPUT_CHANGED')
    return {'status':'PASS_EXACT_COHERENT_POINTWISE_COUNTEREXAMPLES','authentication':before,
        'background_point':{'rho':'1','c':'1','a':'0','theta':'0','p':'1','Vp':'0','Vpp':'0','Vppp':'0','ay':'-4/3','py':'0','V':'13/2'},
        'lowering':{'profiles':'Only A nonzero: A=0, Ay=1, Ayy=1, Ayyy=0; all t jets zero; other profiles identically zero.','value_equal':True,'L0_difference_y_angular':encodeAP(lowerwave[0]),'L0_difference_y_scalar_modes':['1/4','3/2','3'],'identity':'QC_actual-QC_background=(A Ay+8a A^2)Y^2 e^y; L0 difference=3Y^2 e^y at point.','map_substitution_verified':True},
        'first_residual':{'profiles':'q=L_xi g0+d(y)Y dt_frame^2; eta=L_xi phi0; xi=f(y)Y partial_y; f=0, fy=1, fyy=fyyy=fyyyy=0; d=dy=dyy=0,dyyy=1.','E1_all_slots_zero':True,'nabla_y_E1_yy':encodeAP(dE[0,0,0]),'nabla_y_E1_00':encodeAP(dE[0,1,1]),'two_deltaB_y_angular':encodeAP(2*delta[0]),'two_deltaB_y_scalar_modes':['-1/12','-1/2','-1'],'q_yy_value':'2Y','pointwise_prolongation_invalid':True},
        'coherence':'All covariant derivatives formed from commuting ordinary y/t derivatives and the full physical connection; no independent covariant C1 or E1 prolongations introduced.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--mutation',choices=('exchange_QC_under_wave','point_E1_implies_prolongation'));args=ap.parse_args()
    need(not args.output.exists(),'REFUSE_OVERWRITE');result=run(args.repo_root.resolve(),args.mutation)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(result['status'])
if __name__=='__main__':main()
