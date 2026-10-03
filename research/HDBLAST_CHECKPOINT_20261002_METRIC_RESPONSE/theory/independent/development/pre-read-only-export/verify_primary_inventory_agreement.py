#!/usr/bin/env python3
"""Exact independent primary-inventory comparison; no physical evaluation."""
import argparse,contextlib,hashlib,io,json
from pathlib import Path
import sympy as s
p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
here=Path(__file__).resolve().parent
ns={'__file__':str(here/'explore_contacts.py')}
with contextlib.redirect_stdout(io.StringIO()):exec(compile((here/'explore_contacts.py').read_text(),ns['__file__'],'exec'),ns)
L,W,V=s.symbols('L omega v',positive=True);hj=s.symbols('h0:7',real=True)
local={'L':L,'omega':W,'v':V,**{str(t):t for t in hj}}
source=here/'reference_inputs/PRIMARY_FORMULAS.json'
data=json.loads(source.read_text())
expr={key:s.sympify(value,locals=local) for key,value in data['expressions'].items()}
translate={ns['L']:L,ns['w']:W,ns['M']:2*L*L,**dict(zip([ns[n] for n in ('h','h1','h2','h3','h4','h5')],hj))}
convert=lambda x:s.expand(x.subs(translate))
checks=[]
def eq(label,left,right=0):
    residue=s.cancel(left-right)
    if residue!=0:raise RuntimeError(label+': '+str(s.factor(residue)))
    checks.append(label)
poly=lambda rows:sum(s.sympify(coef,locals=local)*W**int(power) for power,coef in rows.items())
for name,reference in [('baseline_scaled_density_subtraction',ns['R']),('baseline_scaled_pressure_subtraction',ns['P']),('variance_subtraction_geo_minus_mass',ns['dS']-ns['dSs']),('density_subtraction_geo_minus_mass',ns['dR']-ns['dRs']),('pressure_subtraction_geo_minus_mass',ns['dP']-ns['dPs'])]:
    eq('primary full Laurent inventory '+name,poly(data['coefficients'][name]),convert(reference))
q0=expr['physical_baseline_Q_times_2pi2'];r0=expr['scaled_baseline_rho_times_2pi2'];p0=expr['scaled_baseline_p_times_2pi2']
def measure_derivative(reference):
    value=convert(reference);total=0
    for power in range(-15,2):
        coef=value.coeff(W,power)
        if coef:total+=coef*(2*L**2)**s.Rational(3+power,2)*V**2/(1-V**2)**s.Rational(5+power,2)
    return total
bare_r_prime=2*L**4*V**3/(1-V**2)**3+3*L**4*V/(2*(1-V**2)**2)
bare_p_prime=2*L**4*V**3/(3*(1-V**2)**3)-L**4*V/(2*(1-V**2)**2)
bare_q_prime=V/(1-V**2)**2
for name,closed,integrand,bare,scale in [('Q',q0,ns['S'],bare_q_prime,L**2),('rho',r0,ns['R'],bare_r_prime,1),('p',p0,ns['P'],bare_p_prime,1)]:
    eq('independent primary baseline primitive derivative '+name,s.diff(closed,V),bare-measure_derivative(integrand)/scale)
    eq('independent primary baseline primitive zero endpoint '+name,closed.subs(V,0))
def integrated_contact(reference):
    value=convert(reference);total=0
    for n in range(5,16,2):
        primitive=(2*L**2)**s.Rational(3-n,2)*sum((-1)**ell*s.binomial((n-5)//2,ell)*V**(2*ell+3)/s.Integer(2*ell+3) for ell in range((n-5)//2+1))
        total+=value.coeff(W,-n)*primitive
    return total
cq=integrated_contact(ns['contacts']['Q']);cr=integrated_contact(ns['contacts']['rho']);cp=integrated_contact(ns['contacts']['p'])
eq('primary finite-K Q bridge convention agreement',expr['Cq_times_2pi2'],-cq)
eq('primary finite-K rho bridge convention agreement',expr['Er_times_2pi2'],cr+(hj[2]+4*L*hj[1])*L**2*q0/2)
eq('primary finite-K p bridge convention agreement',expr['Ep_times_2pi2'],cp-hj[2]*L**2*q0/2)
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
receipt={'status':'PASS','identity_count':len(checks),'identities':checks,'scope':'Exact algebra only; no source/mode/response/quadrature evaluation','primary_inventory_sha256':sha(source),'independent_expanded_chain_rule_sha256':sha(here/'explore_contacts.py'),'source_sha256':sha(Path(__file__)),'sympy':s.__version__,'physical_results_computed':False}
Path(args.output).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('status','identity_count','primary_inventory_sha256','physical_results_computed')}))
