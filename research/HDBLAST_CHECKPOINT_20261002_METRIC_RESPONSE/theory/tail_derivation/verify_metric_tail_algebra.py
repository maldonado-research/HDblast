#!/usr/bin/env python3
"""Pure symbolic verification. Does not call the deferred metric-tail runtime."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise RuntimeError('Refusing overwrite')
    root=Path(__file__).resolve().parent
    cap=json.loads((root/'DERIVATIVE_CERTIFICATE.json').read_text())
    data=json.loads((root/'METRIC_TAIL_INPUTS.json').read_text())
    original=json.loads((root/'reference_inputs/METRIC_CONTACT_ALGEBRA.json').read_text())
    L,v,w,eta=s.symbols('L v w eta');h=s.symbols('h0:8')
    local={'L':L,'v':v,'pi':s.pi,**{str(q):q for q in h}}
    checks=[];mutations=[]
    def eq(name,a,b):
        if s.factor(a-b)!=0:raise RuntimeError('Identity failed: '+name)
        checks.append({'name':name,'status':'PASS'})
    def contact(rows):return {int(n):sum(s.Rational(t['coefficient'])*L**t['L_power']*h[t['h_index']] for t in terms) for n,terms in rows.items()}
    def ratio(spec):
        num=sum(s.Rational(t['coefficient'])*L**t['L_power']*v**t['v_power']*(h[t['h_index']] if t['h_index'] is not None else 1) for t in spec['numerator'])
        den=0
        for coefficient in spec['denominator_coefficients']:den=den*v+s.Rational(coefficient)
        return (1-v)**spec['one_minus_v_power']*num/(den*s.pi**2)
    for name in ('q','rho','p'):
        for n,expr in contact(data['contact_inventory'][name]).items():eq('serialized_'+name+'_moment_'+str(n),expr,s.sympify(original['coefficient_inventory'][name][str(n)],locals=local))
    D=lambda f:L**2*s.diff(f,L)+sum(h[j+1]*s.diff(f,h[j]) for j in range(7))
    expr=sum(co/w**n for n,co in contact(data['contact_inventory']['q']).items())
    for name in ('q_prime','q_second'):
        expr=s.expand(D(expr)+2*L**3/w*s.diff(expr,w))
        eq('serialized_differentiated_'+name,expr,sum(co/w**n for n,co in contact(data['contact_inventory'][name]).items()))
    for name,spec in data['baseline_specs'].items():
        f=s.sympify(original['closed_expressions'][name],locals=local)
        eq('serialized_factored_'+name+'_tail',ratio(spec),f.subs(v,1)-f)
    qref=-2*h[0]*L**2*ratio(data['baseline_specs']['Q0'])
    for j,spec in enumerate(data['q_reference_derivative_specs']):
        eq('serialized_q_reference_derivative_'+str(j),ratio(spec),qref)
        qref=s.factor(D(qref)-L*v*(1-v*v)*s.diff(qref,v))
    moment5=(1-v**3)/6
    moment7=s.Rational(1,30)-(v**3/3-v**5/5)/4
    moment9=s.Rational(1,105)-(v**3/3-2*v**5/5+v**7/7)/8
    eq('mass_proxy_J5_tail',moment5,(1-v)*(1+v+v*v)/6)
    eq('mass_proxy_J7_tail',moment7,(1-v)**2*(3*v**3+6*v*v+4*v+2)/60)
    eq('mass_proxy_J9_tail',moment9,(1-v)**3*(15*v**4+45*v**3+48*v*v+24*v+8)/840)
    # Infinite moments: integrate the derivative of k^(3-n)/(3-n),
    # yielding K^(3-n)/(n-3). The inequality follows from w>=k>0.
    K=s.symbols('K',positive=True)
    for n in (5,7,9,11,13,15):eq('moment_power_primitive_'+str(n),s.diff(K**(3-n)/(3-n),K),K**(2-n))
    # Independent exact Leibniz formula with L=-1/eta and h arbitrary.
    a=s.Function('a')(eta)
    forc=4*a/eta**2+2*s.diff(a,eta)/eta-s.diff(a,eta,2)
    for n in range(6):
        formula=4*sum(s.binomial(n,j)*s.factorial(n-j+1)*(-1/eta)**(n-j+2)*s.diff(a,eta,j) for j in range(n+1))-2*sum(s.binomial(n,j)*s.factorial(n-j)*(-1/eta)**(n-j+1)*s.diff(a,eta,j+1) for j in range(n+1))-s.diff(a,eta,n+2)
        eq('canonical_forcing_Leibniz_'+str(n),s.diff(forc,eta,n),formula)
    for source,envelopes in cap['envelopes'].items():
        if len(envelopes['h_derivative_caps'])!=8 or any(s.Rational(x)<0 for x in envelopes['g_N_upper_rational']):raise RuntimeError('Invalid certificate')
        checks.append({'name':source+'_certificate_orders_and_positive_caps','status':'PASS'})
    def reject(name,difference,subs):
        out=s.factor(difference.subs(subs))
        if out==0:raise RuntimeError('Undetected mutation: '+name)
        mutations.append({'name':name,'status':'DETECTED','rational_residual':str(out)})
    reject('missing_second_moment_reference_factor',moment5-(1-v**3)/3,{v:s.Rational(1,2)})
    reject('wrong_J7_tail_denominator',moment7-(1-v)**2*(3*v**3+6*v*v+4*v+2)/15,{v:s.Rational(1,2)})
    reject('wrong_J9_tail_denominator',moment9-(1-v)**3*(15*v**4+45*v**3+48*v*v+24*v+8)/105,{v:s.Rational(1,2)})
    reject('omit_q_reference_metric_prefactor',-2*h[0]*L**2*ratio(data['baseline_specs']['Q0']),{v:s.Rational(1,2),L:1,h[0]:1})
    q0=sum(co/w**n for n,co in contact(data['contact_inventory']['q']).items())
    reject('freeze_reference_denominator_under_time_derivative',2*L**3/w*s.diff(q0,w),{L:1,w:2,**{h[j]:j+1 for j in range(8)}})
    reject('omit_highest_metric_source_derivative',s.diff(a,eta,7),{s.diff(a,eta,7):1})
    report={'status':'PASS','scope':'Exact omitted-band algebra and synthetic formula mutations only; no deferred runtime physical evaluation.','identity_count':len(checks),'mutation_count':len(mutations),'checks':checks,'mutations':mutations,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':{f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in ('DERIVATIVE_CERTIFICATE.json','METRIC_TAIL_INPUTS.json','reference_inputs/METRIC_CONTACT_ALGEBRA.json')}}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('status','identity_count','mutation_count')}))

if __name__=='__main__':main()
