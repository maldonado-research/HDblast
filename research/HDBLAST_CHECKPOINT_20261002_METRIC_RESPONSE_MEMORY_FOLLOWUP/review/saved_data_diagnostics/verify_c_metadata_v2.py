#!/usr/bin/env python3
"""Exact background derivative proof only; no physical profiles or modes."""
import hashlib
import json
from pathlib import Path
import sys
import sympy as S

eta=S.symbols('eta',negative=True);a=S.symbols('a',positive=True)
n=S.symbols('n',integer=True,nonnegative=True)
checks=[]
def exact(name,left,right):
    if S.simplify(left-right)!=0:raise RuntimeError(name)
    checks.append({'name':name,'exact':True})
formula=2*S.factorial(n+1)*a**(n+2)
exact('C0_equals_a_double_prime_over_a',formula.subs(n,0),S.diff(-1/eta,eta,2).subs(eta,-1/a)/a)
exact('all_n_derivative_recurrence',a*a*S.diff(formula,a),formula.subs(n,n+1))
for order in range(8):exact('direct_eta_derivative_'+str(order),S.diff(2/eta**2,eta,order),formula.subs({n:order,a:-1/eta}))
wrong=S.factorial(n+2)*a**(n+2)
witnesses=[]
for order,expected in [(1,2*a**3),(2,12*a**4)]:
    gap=S.simplify(wrong.subs(n,order)-formula.subs(n,order));exact('wrong_formula_witness_'+str(order),gap,expected)
    if gap==0:raise RuntimeError('Wrong formula witness disappeared')
    witnesses.append({'n':order,'old_minus_correct':str(gap),'nonzero_for_a_positive':True})
base=Path(__file__).resolve().parent
old=base/'development/original-strict-wrapper/raw_metric_audit.py';new=base/'raw_metric_audit_v2.py'
if new.read_text()!=old.read_text().replace('math.factorial(n+2)*a**(n+2)','2*math.factorial(n+1)*a**(n+2)'):raise RuntimeError('Unexpected V2 raw auditor change')
checks.append({'name':'only_C_metadata_expression_changed','exact':True})
print(json.dumps({'proof_status':'EXACT_METADATA_IDENTITY_PROOF','python_optimization':sys.flags.optimize,'formula':'C_n=2*(n+1)!*a^(n+2)','identity_count':len(checks),'checks':checks,'wrong_formula_witnesses':witnesses,'original_auditor_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),'v2_auditor_sha256':hashlib.sha256(new.read_bytes()).hexdigest(),'physical_evaluations':0,'scientific_calibration_status':'FAIL_UNCHANGED'},indent=2))
