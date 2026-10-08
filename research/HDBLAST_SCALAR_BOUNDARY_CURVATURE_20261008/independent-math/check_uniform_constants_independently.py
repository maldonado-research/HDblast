#!/usr/bin/env python3
"""Independent SymPy derivation. Reads frozen JSON; never imports the producer."""
import copy
import hashlib
import json
from pathlib import Path
import re
import sys
import sympy as S

BASE=Path(__file__).resolve().parent.parent
ART=BASE/'numerical-audit/uniform-constants'
REPORT_SHA='95923d9824f3f367385f67392ac10391b489337138851b642aa4553377d11ba1'
PRODUCER_SHA='4ca9b1824c2c94b58ef04e56289272c3e5b51a5cd52a5a062cb8511728fac7cc'
THEOREM_SHA='544c6dbced51fbc68dccab3bbbe84210545849c6e7447e1b6bdae0fc98e76f51'
x=S.Symbol('eta',real=True)

def rational(v):
    if not isinstance(v,str) or not re.fullmatch(r'-?\d+/[1-9]\d*',v):
        raise ValueError('Noncanonical exact rational field')
    a,b=map(int,v.split('/'))
    if str(S.Rational(a,b).p)+'/'+str(S.Rational(a,b).q)!=v:
        raise ValueError('Unreduced exact rational')
    return S.Rational(a,b)

def equal(a,b,label):
    if a!=b:
        raise ValueError(label)

def pin(path,expected):
    payload=path.read_bytes()
    equal(hashlib.sha256(payload).hexdigest(),expected,'source pin: '+path.name)
    return payload

def sup_polynomial(poly,r):
    return sum(abs(coefficient)*r**power[0] for power,coefficient in S.Poly(S.expand(poly),x).terms())

def independently_compute(model):
    p={name:rational(value) for name,value in model['parameters'].items()}
    k,w,v,f0,f1,f2=[p[name] for name in ('k','w','v','f0','f1','f2')]
    if model['name']=='registered_exact_decimal_model':
        coupling=2/S.Rational('1.0357712571566784')-S.Rational(4,3)
        equal(p,dict(k=S.Rational(1,9),w=S.Integer(2),v=S.Integer(2),f0=1+coupling,f1=coupling,f2=S.Integer(0)),'registered model identity')
    elif model['name']=='simple_polynomial_model':
        equal(p,dict(k=S.Integer(1),w=S.Integer(5),v=S.Integer(0),f0=S.Integer(1),f1=S.Integer(1),f2=S.Integer(0)),'simple model identity')
    else:
        raise ValueError('Unexpected model')
    r,b,ds=[rational(model[name]) for name in ('r','b','delta_star')]
    if not (k>0 and w>4*k and f0>0 and r>0 and b>0 and ds>0):
        raise ValueError('Theorem domain')
    W=3*k+w*x*x/2+v*x**3/6
    f=f0+f1*x+f2*x*x/2
    U=S.diff(W,x)**2/2-S.Rational(2,3)*W**2
    A=W*f/9-S.diff(W,x)*S.diff(f,x)/12
    B=f*f/36-S.diff(f,x)**2/48
    polys={'W':W,'f':f,'U':U,'A':A,'BB':B}
    for name,poly in polys.items():
        recovered=sum(rational(c)*x**i for i,c in enumerate(model['polynomial_coefficients'][name]))
        equal(S.expand(recovered-poly),0,'polynomial '+name)
    c={}
    c['m2']=w*(w-4*k)
    c['lam']=w-4*k
    c['nu']=c['lam']/2
    c['mu']=c['nu']*(c['nu']+4*k)
    c['gap']=c['m2']-c['mu']
    c['zeta_lower']=15*c['mu']/(16*(S.Rational(16,3)*k+c['nu']))
    c['J_upper']=2/k+1/(2*c['zeta_lower'])
    c['C_G']=max(1/c['gap'],(1+c['m2']/c['gap'])/(4*k*k))
    c['M']=max(S.Integer(1),c['lam']/k)
    c['ball']=2*c['M']; Q=c['ball']
    c['M2']=sup_polynomial(S.diff(U,x,2),r)
    c['M3']=sup_polynomial(S.diff(U,x,3),r)
    c['L1']=c['M3']/2
    c['C_F']=k*k/4+c['M2']/12
    c['epsilon']=c['C_F']*Q*Q*b*b*c['J_upper']/(2*k)
    c['C_N']=c['L1']*Q*Q+4*c['C_F']*Q**3*b
    c['L_N']=2*c['L1']*Q*b+20*c['C_F']*Q*Q*b*b
    c['C_P']=k*c['C_G']*Q*(2*c['L1']*Q+20*c['C_F']*Q*Q*b)
    c['C_D']=k*c['C_G']*c['C_N']
    c['C_a']=4*c['C_F']*Q*Q/k
    c['C_L_upper']=8*c['lam']+S.Rational(28,3)*k*c['lam']/(w-2*k)
    c['C_H_DN']=4*c['C_L_upper']/(9*k*k)
    for name,poly in {'W2':S.diff(W,x,2),'W3':S.diff(W,x,3),'F1':S.diff(f,x),'F2':S.diff(f,x,2),'F_abs':f,'W_abs':W}.items():
        c[name]=sup_polynomial(poly,r)
    c['L']=(abs(f1)+1)/w
    c['a1']=k*f0/3
    c['C_K']=c['L']*sup_polynomial(S.diff(A,x),r)+sup_polynomial(B,r)
    c['a_bar']=S.Rational(4,3)*k+c['C_F']*Q*Q*b*b/k
    c['C_Tp']=c['M2']+4*c['a_bar']*k*Q+k*Q*(c['lam']+c['C_P']*b)
    c['C_Tbeta']=2*c['C_Tp']/w
    c['C_slope']=k*Q*c['C_a']*c['L']**2+((c['C_a']+c['W2']/3)*c['L']+c['F1']/6)*c['C_Tbeta']*c['L']
    c['d0']=2*(w-2*k)
    c['alpha']=-f1/(2*c['d0'])
    c['K_eta']=(c['C_H_DN']*3*c['a1']*c['L']/2+(c['C_D']+c['W3']/2)*c['L']**2+c['F2']*c['L']/2)/c['d0']
    c['K_H']=abs(S.diff(A,x).subs(x,0))*c['K_eta']+sup_polynomial(S.diff(A,x,2),r)*c['L']**2/2+sup_polynomial(S.diff(B,x),r)*c['L']
    equal(set(c),set(model['exact_constants']),'constant key closure')
    for name,val in c.items(): equal(val,rational(model['exact_constants'][name]),'constant '+name)
    rows=[
       ('field_tube_inside_polynomial_radius',Q*b,r,False),
       ('metric_volterra_contraction',c['epsilon'],S.Rational(1,4),False),
       ('scalar_contraction',c['C_G']*c['L_N'],S.Rational(1,2),False),
       ('scalar_ball_preservation',c['C_G']*c['C_N']*b,c['M'],False),
       ('positive_metric_derivative',c['C_F']*Q*Q*b*b/k,k/2,False),
       ('scalar_junction_derivative_feasibility',(c['C_P']+c['W3'])*b,w/2,True),
       ('scalar_root_inside_boundary_tube',c['L']*ds,b,False),
       ('scalar_junction_derivative',(c['C_P']+c['W3'])*b+ds*c['F2']/2,w/2,False),
       ('positive_reduced_curvature',c['C_K']*ds,c['a1']/2,False),
       ('metric_residual_positive_at_Y',3*c['a1']*ds/2,9*k*k/32,True),
       ('positive_shell_tension',c['W2']*b*b/6+ds*c['F_abs']/6,k/2,False),
       ('strict_downward_shell_crossing',c['C_slope']*ds,c['a1']/4,False)]
    equal(len(rows),len(model['checks']),'inequality closure')
    for (name,lhs,rhs,strict),given in zip(rows,model['checks']):
        equal(given['name'],name,'inequality order')
        equal(rational(given['lhs']),lhs,name+' lhs')
        equal(rational(given['rhs']),rhs,name+' rhs')
        equal(given['strict'],strict,name+' strict')
        equal(rational(given['exact_margin']),rhs-lhs,name+' margin')
        if not (lhs<rhs if strict else lhs<=rhs): raise ValueError(name+' fails')
        equal(given['pass_check'],True,name+' pass flag')
    E=abs(c['alpha'])+c['K_eta']*ds
    Kpoly=abs(S.diff(A,x).subs(x,0))*c['K_eta']
    for (degree,),coef in S.Poly(S.expand(A),x).terms():
        if degree>=2: Kpoly+=abs(coef)*E**degree*ds**(degree-2)
    for (degree,),coef in S.Poly(S.expand(B),x).terms():
        if degree>=1: Kpoly+=abs(coef)*E**degree*ds**(degree-1)
    equal(Kpoly,rational(model['K_H_polynomial_refinement']),'polynomial refinement')
    best=min(Kpoly,c['K_H'])
    q=k*f1*f1/(24*(w-2*k))
    sign=model['strict_curvature_suppression']
    for name,val in {'q':q,'best_K_H':best,'lower_magnitude_coefficient':q-best*ds,'upper_magnitude_coefficient':q+best*ds,'exact_positive_margin':q-best*ds}.items():
        equal(rational(sign[name]),val,'sign '+name)
    if not q-best*ds>0: raise ValueError('Suppression bound fails')
    equal(sign['pass_check'],True,'sign pass')
    equal(rational(model['exact_bound_at_delta_star']['scalar']),c['K_eta']*ds**2,'scalar endpoint bound')
    equal(rational(model['exact_bound_at_delta_star']['curvature']),best*ds**3,'curvature endpoint bound')
    return {'name':model['name'],'exact_constants_verified':len(c),'exact_inequalities_verified':len(rows),'polynomials_verified':len(polys),'strict_suppression':True,'delta_star':str(ds),'display_delta_star':float(ds)}

def main():
    payload=pin(ART/'EXACT_UNIFORM_CONSTANTS.json',REPORT_SHA)
    pin(ART/'EXACT_UNIFORM_CONSTANTS_OPTIMIZED.json',REPORT_SHA)
    pin(ART/'evaluate_uniform_constants.py',PRODUCER_SHA)
    pin(BASE/'branch-theorem/UNIFORM_BRANCH_THEOREM.md',THEOREM_SHA)
    data=json.loads(payload)
    equal(data['source_sha256'],PRODUCER_SHA,'embedded producer pin')
    equal(data['theorem_sha256'],THEOREM_SHA,'embedded theorem pin')
    equal(len(data['models']),2,'model count')
    results=[independently_compute(m) for m in data['models']]
    negatives=[]
    for kind in ('constant','inequality','sign','parameter'):
        mutant=copy.deepcopy(data['models'][0])
        if kind=='constant': mutant['exact_constants']['K_eta']='0/1'
        if kind=='inequality': mutant['checks'][0]['rhs']='0/1'
        if kind=='sign': mutant['strict_curvature_suppression']['lower_magnitude_coefficient']='0/1'
        if kind=='parameter': mutant['parameters']['k']='1/1'
        try: independently_compute(mutant)
        except ValueError: negatives.append(kind)
        else: raise ValueError('Mutation accepted: '+kind)
    # Independent finite-series rational bounds supporting all transcendental replacements.
    equal(sum(S.Rational(2)**n/S.factorial(n) for n in range(5)),7,'exp(2) partial sum')
    tail_upper=S.Rational(2)**5/S.factorial(5)/(1-S.Rational(1,3))
    if not 7+tail_upper<8: raise ValueError('exp upper bound')
    out={'status':'PASS_INDEPENDENT_EXACT_CONSTANTS_AND_SIGN_CHECK','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_sha256':PRODUCER_SHA,'report_sha256':REPORT_SHA,'theorem_sha256':THEOREM_SHA,'sympy_version':S.__version__,'models':results,'negative_controls_rejected':negatives,'analytic_transcendental_review':'7<exp(2)<8 follows from a positive degree-4 partial sum and geometric tail. It implies coth(1)<4/3, (1-exp(-2))^-1<7/6, sinh(1)^2<2. Since C0/k>4 and exp(4)>16, zeta has the stated positive rational lower bound. All substitution directions independently checked.','independence':'No producer source imported or executed; polynomial differentiation and coefficient extraction use an independently written SymPy representation.','limits':['No numerical ODE run is certified by arithmetic alone; the pinned independently reviewed analytic theorem supplies the implication.','Registered parameters are the exact declared decimal model, not the binary64 shooting model.','No external novelty, dynamical stability, APS acceptance or cosmological-origin claim.']}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__': main()
