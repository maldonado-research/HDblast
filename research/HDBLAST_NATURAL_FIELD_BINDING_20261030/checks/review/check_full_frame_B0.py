#!/usr/bin/env python3
"""Independent full frame contraction on arbitrary full23 metric amplitudes.

The parent checker's geometry and Bianchi calculation are never imported.
The shared helpers implement exact Laurent arithmetic and intrinsic sphere
polynomials only. Coordinate Christoffels for the warped physical frame are
contracted as a separately constructed coefficient table.
"""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
PINS = {'codazzi_laurent.py':'178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2',
        'angular_generic.py':'fd14d8cc3ac6f58c206b784c5c57d15430f63f38ca091ac03c7769f9bfc89df9'}

def need(ok, why):
    if not ok: raise ValueError(why)

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def authenticate():
    hashes = {}
    for name, pin in PINS.items():
        path = PARENT / name
        need(path.is_file() and not path.is_symlink(), 'HELPER_TYPE:'+name)
        hashes[name] = digest(path)
        need(hashes[name] == pin, 'HELPER_PIN:'+name)
    for path in (Path(__file__).resolve(), *(PARENT/name for name in PINS)):
        need(not any(isinstance(n, ast.Assert) for n in ast.walk(ast.parse(path.read_text()))), 'ASSERT:'+path.name)
    hashes[Path(__file__).name] = digest(Path(__file__))
    return hashes

def load(name):
    spec = importlib.util.spec_from_file_location('review_'+name, PARENT/(name+'.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def run(mutation=None):
    before = authenticate()
    alg = load('codazzi_laurent'); sphere = load('angular_generic'); P = alg.P
    I = range(6); J = range(4); sign = (1,-1,1,1,1,1)
    rho, c, a, theta = [P.var(n) for n in ('rho','c','a','theta')]
    H = theta/rho
    def add(*xs): return sphere.add(*xs)
    def scale(s,x): return sphere.scale(s,x)
    def mul(x,y): return sphere.mul(x,y)
    def ordinary(x, axis):
        return {e:sum((v.diff(n)*P.var(n+axis) for n in v.names()),P()) for e,v in x.items()}
    coefficients = {prefix+str(l):P.var(prefix+str(l)) for prefix in ('K','N','W','J') for l in (0,2,4)}
    coefficients.update({prefix+str(l):P.var(prefix+str(l)) for prefix in ('v','u','h') for l in (2,4)})
    coefficients.update(Tplus=P.var('Tplus'), Tminus=P.var('Tminus'))
    T = {key:{} for key in product(I,repeat=2)}
    for j,l in enumerate((0,2,4)):
        base = sphere.SCALARS[j][()]
        T[0,0] = add(T[0,0],scale(coefficients['K'+str(l)],base))
        T[0,1] = add(T[0,1],scale(coefficients['W'+str(l)],base))
        T[1,1] = add(T[1,1],scale(coefficients['N'+str(l)],base))
        for i,k in product(J,repeat=2):
            T[i+2,k+2] = add(T[i+2,k+2],scale(coefficients['J'+str(l)],sphere.TENSORS[j][i,k]))
    T[1,0] = T[0,1]
    for j,l in enumerate((2,4)):
        for i in J:
            T[0,i+2] = T[i+2,0] = add(T[0,i+2],scale(coefficients['v'+str(l)],sphere.GRADS[j][i,]))
            T[1,i+2] = T[i+2,1] = add(T[1,i+2],scale(coefficients['u'+str(l)],sphere.GRADS[j][i,]))
        for i,k in product(J,repeat=2):
            T[i+2,k+2] = add(T[i+2,k+2],scale(coefficients['h'+str(l)],sphere.HESS[j][i,k]))
    for j,label in enumerate(('Tplus','Tminus')):
        for i,k in product(J,repeat=2):
            T[i+2,k+2] = add(T[i+2,k+2],scale(coefficients[label],sphere.TT[j][i,k]))
    # Independent connection table Gamma^upper_(direction,slot).
    gamma = {(u,d,v):{} for u,d,v in product(I,repeat=3)}
    gamma[1,1,0] = {sphere.ZERO:a}
    gamma[0,1,1] = {sphere.ZERO:a}
    for i,k in product(J,repeat=2):
        projector = sphere.P[i,k]
        gamma[k+2,i+2,0] = scale(a,projector)
        gamma[0,i+2,k+2] = scale(-a,projector)
        gamma[k+2,i+2,1] = scale(H,projector)
        gamma[1,i+2,k+2] = scale(H,projector)
    if mutation == 'time_connection_sign': gamma[0,1,1] = {sphere.ZERO:-a}
    def frame_derivative(d,i,j):
        if d == 0: out = ordinary(T[i,j],'y')
        elif d == 1: out = scale(1/rho,ordinary(T[i,j],'t'))
        else:
            k = d-2
            out = add(*(mul(sphere.P[k,v],sphere.diff(T[i,j],v)) for v in J))
            if i>=2: out = add(out,mul(sphere.X[i-2],T[d,j]))
            if j>=2: out = add(out,mul(sphere.X[j-2],T[i,d]))
            out = scale(1/(rho*c),out)
        return add(out,*(scale(-1,mul(gamma[u,d,i],T[u,j])) for u in I),
                       *(scale(-1,mul(gamma[u,d,j],T[i,u])) for u in I))
    full = {b:add(*(scale(sign[k],add(frame_derivative(k,k,b),scale(-F(1,2),frame_derivative(b,k,k)))) for k in I)) for b in I}
    rows = {'y':[],'0':[],'a':[]}
    for l in (0,2,4):
        K,N,W,J0 = [coefficients[prefix+str(l)] for prefix in ('K','N','W','J')]
        dy = lambda x:sum((x.diff(n)*P.var(n+'y') for n in x.names()),P())
        dt = lambda x:sum((x.diff(n)*P.var(n+'t') for n in x.names()),P())
        lam=l*(l+2);v=coefficients.get('v'+str(l),P());u=coefficients.get('u'+str(l),P())
        rows['y'].append((dy(K)+dy(N)-3*dy(J0))/2-dt(W)/rho+a*(4*K+N-3*J0)-3*H*W-lam*v/(rho*c))
        rows['0'].append(dy(W)-(dt(K)+dt(N)+3*dt(J0))/(2*rho)+5*a*W-3*H*(N+J0)-lam*u/(rho*c))
        if l:
            h=coefficients['h'+str(l)];beta=2-F(2,3)*lam
            rows['a'].append(dy(v)-dt(u)/rho+5*a*v-4*H*u+(-K/2+N/2-J0/2+beta*h)/(rho*c))
    expected = {0:add(*(scale(cf,base[()]) for cf,base in zip(rows['y'],sphere.SCALARS))),
                1:add(*(scale(cf,base[()]) for cf,base in zip(rows['0'],sphere.SCALARS)))}
    for i in J: expected[i+2]=add(*(scale(cf,base[i,]) for cf,base in zip(rows['a'],sphere.GRADS)))
    residual_counts=[]
    for b in I:
        residual=add(full[b],scale(-1,expected[b]));need(not residual,'FULL_FRAME_B0_NONZERO:'+str(b)+':'+str(list(residual.items())[:1]))
        residual_counts.append({'ambient_covector_slot':b,'angular_terms':len(full[b]),'residual_terms':0})
    # Modal uniqueness is established by a nonsingular diagonal Gram matrix.
    for label,basis in (('scalar',sphere.SCALARS),('gradient',sphere.GRADS)):
        for i,j in product(range(len(basis)),repeat=2):
            inner=sphere.inner(basis[i],basis[j]);need(bool(inner) if i==j else not inner,'GRAM:'+label+':'+str((i,j)))
    # The arbitrary TT amplitudes and all their first derivatives disappear.
    for component in full.values():
        for value in component.values():
            need(not any(n.startswith('Tplus') or n.startswith('Tminus') for n in value.names()),'TT_NOT_IN_KERNEL')
    need(authenticate()==before,'HELPERS_CHANGED')
    return {'status':'PASS_INDEPENDENT_FULL_FRAME_B0_EIGHT_ROWS','authentication':before,
            'arbitrary_metric_amplitudes':20,'arbitrary_first_amplitude_derivatives':40,
            'modal_rows':8,'full_ambient_covector_rows':residual_counts,
            'expected_modal_rows':{field:[x.encode() for x in values] for field,values in rows.items()},
            'independence':'Separate explicit upper connection table; parent Bianchi/geometry checker never imported. Shared exact arithmetic and sphere primitives only.',
            'TT_kernel_verified':True,'background_equations_required':False,
            'scope':'The operator holds for arbitrary physical tensor amplitudes and derivatives, hence for the registered natural QE after substitution. It does not independently rederive QE or the first residual correction.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--mutation',choices=('time_connection_sign',));args=ap.parse_args()
    need(not args.output.exists(),'REFUSE_OVERWRITE')
    report=run(args.mutation);args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(report['status'])
if __name__=='__main__':main()
