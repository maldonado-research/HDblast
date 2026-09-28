#!/usr/bin/env python3
"""Independent exact initial-data compatibility identities and saved-result checks.

Requires SymPy. Does not import producers, solve shooting equations, or evolve.
"""
from pathlib import Path
import hashlib,json
import sympy as S

r,v,q,s,u,up,upp,e,b,bp=S.symbols('r v phi s U Up Upp epsilon bump bump_phi',real=True)
sig,sig1,sig2=S.symbols('sigma sigma1 sigma2',real=True)
rho2=(1-v*v)/r-r*s*s/6-r*u/3
phi2=up-4*v*s/r+e*b
D=1-v*v+r*r*(s*s/12-u/6)
def dy(expr):
    return S.diff(expr,r)*v+S.diff(expr,v)*rho2+S.diff(expr,q)*s+S.diff(expr,s)*phi2+S.diff(expr,u)*up*s+S.diff(expr,b)*bp*s
checks=[]
def check(name,expr,zero=True):
    out=S.factor(S.expand(expr));ok=(out==0) if zero else (out!=0)
    checks.append(dict(name=name,pass_check=bool(ok),negative_control=not zero))
    if not ok:raise RuntimeError(name+': '+str(out))

H=-2*r*r*u+6-12*v*v+6*v*v-6*r*rho2-r*r*s*s
M=S.Integer(0)  # A_t=1, B_t=phi_t=0, A_z=B_z=v.
aa=r*rho2-3+3*v*v+S.Rational(2,3)*r*r*u
bb=r*rho2+3-3*v*v+r*r*s*s/2-r*r*u/3
ff=r*r*phi2+4*r*v*s-r*r*up
check('full initial Hamiltonian after proper-distance constraint ODE',H)
check('initial momentum for equal A_z and B_z',M)
check('D transport from differentiated definition',dy(D)+2*v*D/r-e*r*r*s*b/6)
check('integrating-factor ledger',dy(r*r*D)-e*r**4*s*b/6)
check('A acceleration equals minus two D',aa+2*D)
check('B acceleration equals four D',bb-4*D)
check('scalar acceleration is initial bump',ff-e*r*r*b)

# Differentiating the actual junction twice in time, then inserting t=0 data.
# B_t=phi_t=0 eliminates all products of their first time derivatives.
second_A=r*dy(aa)-r*(sig*bb+sig1*ff)/6
second_B=r*dy(bb)-r*(sig*bb+sig1*ff)/6
second_F=r*dy(ff)+r*(sig1*bb+sig2*ff)/2
junction={sig:6*v/r,sig1:-2*s}
check('A second corner including nonzero bump',second_A.subs(junction))
check('B general second corner',second_B.subs(junction)-(-12*v*D+e*r**3*s*b))
check('scalar general second corner',second_F.subs(junction)-(-4*r*s*D+e*r*r*b*(2*v+r*sig2/2)+e*r**3*s*bp))
check('B corner on bump-free neighborhood',second_B.subs(junction).subs({b:0,bp:0})+12*v*D)
check('scalar corner on bump-free neighborhood',second_F.subs(junction).subs({b:0,bp:0})+4*r*s*D)
check('reject wrong sign in D transport',dy(D)-2*v*D/r-e*r*r*s*b/6,False)
check('reject omitting scalar boundary acceleration',second_F.subs(junction).subs({b:0,bp:0}),False)
check('b value zero alone does not remove its spatial derivative',second_F.subs(junction).subs(b,0)+4*r*s*D,False)
check('reject replacing B corner coefficient twelve by eight',second_B.subs(junction).subs({b:0,bp:0})+8*v*D,False)

# Static-neighborhood implication: D=0 eliminates 1-v^2 and turns the
# Hamiltonian ODE into the static Einstein radial equation.
static_v2=1+r*r*(s*s/12-u/6)
check('D zero gives the static radial Einstein ODE',rho2.subs(v*v,static_v2)+r*s*s/4+r*u/6)

base=Path(__file__).resolve().parents[1]/'initial_data'
original=json.loads((base/'CONSTRAINT_SEED_RESULTS.json').read_text())
balanced=json.loads((base/'BALANCED_CONSTRAINT_SEED_RESULTS.json').read_text())
checks_saved=[]
for row in balanced['rows']:
    for key,rr in [('nominal',row),('refined',row['refined'])]:
        ok=rr['bump_at_brane']==0 and rr['phi_b']+1>.85 and max(abs(x) for x in rr['junction_residual'])<1e-10
        checks_saved.append(dict(epsilon=rr['epsilon'],version=key,pass_check=bool(ok),
            open_scalar_gap_to_bump=float(rr['phi_b']+1-.85),direct_D=rr['mass_defect_b'],
            integral_D=rr['mass_defect_integral'],second_corner_B=rr['second_corner_B']))
        if not ok:raise RuntimeError('Saved candidate fails claimed gap/junction bounds')
        if not abs(rr['mass_defect_b'])>100*abs(rr['mass_defect_integral']):
            raise RuntimeError('Unexpected balance ledger precision relationship')

result={'status':'PASS_INDEPENDENT_IDENTITIES_AND_SAVED_DATA_CHECKS',
        'positive_symbolic_checks':sum(not x['negative_control'] for x in checks),
        'negative_controls':sum(x['negative_control'] for x in checks),'symbolic_checks':checks,
        'saved_balanced_candidates_checked':len(checks_saved),'saved_data_checks':checks_saved,
        'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [base/'constraint_seed.py',base/'balanced_constraint_seed.py',base/'verify_constraint_seed_algebra.py',
             base/'CONSTRAINT_SEED_RESULTS.json',base/'BALANCED_CONSTRAINT_SEED_RESULTS.json']},
        'scope':'Exact conditional identities and read-only saved floating-result checks; no new shooting or evolution and no interval existence proof.'}
Path(__file__).with_name('SEED_COMPATIBILITY_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','positive_symbolic_checks','negative_controls','saved_balanced_candidates_checked']}))
