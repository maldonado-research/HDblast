#!/usr/bin/env python3
"""Pure symbolic preparation of contact and baseline omitted-band inventories."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

def require(c,m):
    if not c:raise RuntimeError(m)

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    require(not args.output.exists(),'Refusing overwrite')
    root=Path(__file__).resolve().parent
    ip=root/'reference_inputs/METRIC_CONTACT_ALGEBRA.json'
    data=json.loads(ip.read_text())
    L,v=s.symbols('L v');hj=s.symbols('h0:8');local={'L':L,'v':v,'pi':s.pi,**{str(h):h for h in hj}}
    D=lambda f:s.expand(L**2*s.diff(f,L)+sum(hj[j+1]*s.diff(f,hj[j]) for j in range(7)))
    terms=lambda f:[{'coefficient':str(co),'L_power':po[0],'h_index':next(j for j,p in enumerate(po[1:]) if p)} for po,co in s.Poly(s.expand(f),L,*hj).terms()]
    inventory={name:{int(n):s.sympify(expr,locals=local) for n,expr in rows.items()} for name,rows in data['coefficient_inventory'].items() if name in ('q','rho','p')}
    for name,rows in inventory.items():
        require(all(n>=5 and n%2==1 for n in rows),'Nonintegrable contact moment')
        for expr in rows.values():require(all(sum(po[1:])==1 for po,co in s.Poly(expr,L,*hj).terms()),'Contact not linear in metric jets')
    out={name:{str(n):terms(expr) for n,expr in rows.items()} for name,rows in inventory.items()}
    qder=[inventory['q']]
    for j in range(2):
        nxt={}
        for n,co in qder[-1].items():
            nxt[n]=nxt.get(n,0)+D(co)
            nxt[n+2]=nxt.get(n+2,0)-2*n*L**3*co
        qder.append({n:s.expand(co) for n,co in nxt.items()})
    out['q_prime']={str(n):terms(expr) for n,expr in qder[1].items()}
    out['q_second']={str(n):terms(expr) for n,expr in qder[2].items()}
    baselines={name:s.factor(s.sympify(data['closed_expressions'][name],locals=local).subs(v,1)-s.sympify(data['closed_expressions'][name],locals=local)) for name in ('Q0','rho0','p0')}
    def ratio(expr):
        num,den=s.fraction(s.factor(expr*s.pi**2))
        m=0
        while s.simplify(num.subs(v,1))==0:
            num=s.cancel(num/(1-v));m+=1
        require(s.Poly(den,v).degree()>=0,'Unsupported rational denominator')
        require(all(sum(po[2:])<=1 for po,co in s.Poly(num,L,v,*hj).terms()),'Nonlinear baseline contact')
        numerator=[]
        for po,co in s.Poly(num,L,v,*hj).terms():
            ids=[j for j,power in enumerate(po[2:]) if power]
            numerator.append({'coefficient':str(co),'L_power':po[0],'v_power':po[1],'h_index':ids[0] if ids else None})
        return {'one_minus_v_power':m,'numerator':numerator,'denominator_coefficients':[str(c) for c in s.Poly(den,v).all_coeffs()],'expression':str(s.factor(expr))}
    baseline_specs={name:ratio(f) for name,f in baselines.items()}
    Dv=lambda f:s.factor(D(f)-L*v*(1-v*v)*s.diff(f,v))
    f=-2*hj[0]*L**2*baselines['Q0']
    reference_q=[]
    for j in range(3):reference_q.append(ratio(f));f=Dv(f)
    # Cross-check differentiated moment inventories against differentiation
    # of a formal fixed-k integrand, with w'=2L^3/w independently applied.
    w=s.symbols('w')
    exact=sum(co/w**n for n,co in inventory['q'].items())
    checks=[]
    for j in range(3):
        target=sum(co/w**n for n,co in qder[j].items())
        require(s.factor(exact-target)==0,'Derivative moment inventory mismatch')
        checks.append({'name':'q_contact_derivative_'+str(j),'status':'PASS'})
        exact=s.expand(D(exact)+2*L**3/w*s.diff(exact,w))
    report={'status':'PASS','scope':'Exact symbolic moment inventories and baseline difference factorization; no physical evaluation.','contact_inventory':out,'baseline_specs':baseline_specs,'q_reference_derivative_specs':reference_q,'checks':checks,'source_contact_sha256':hashlib.sha256(ip.read_bytes()).hexdigest(),'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'moment_inequality':'For n>=5, integral_K^infinity k²/(k²+2L²)^(n/2) dk <= K^(3-n)/(n-3).','coefficient_derivative_rule':'D=L² partial_L+sum h_(j+1) partial_hj; differentiated denominator contributes -2nL³ c_n at power n+2.'}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','contacts':{n:len(x) for n,x in out.items()},'checks':len(checks),'baseline_one_minus_v_powers':{name:spec['one_minus_v_power'] for name,spec in baseline_specs.items()}}))

if __name__=='__main__':main()
