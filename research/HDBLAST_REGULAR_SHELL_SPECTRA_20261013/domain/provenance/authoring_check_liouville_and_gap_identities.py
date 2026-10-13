"""Exact algebra checks supplementing, not replacing, the domain/Hardy proof."""
from pathlib import Path
import hashlib
import json
import sympy as s

y=s.symbols('y', positive=True)
rho=s.Function('rho')(y)
p=s.Function('p')(y)
F=s.Function('F')
x=s.Function('x')(y)
psi=F(x)
alpha=1/(rho**2*p**2)
w=1/(rho**4*p**2)
q=2/(3*rho**2)
# A fixed momentum sign can be chosen because p has no zeros. Replacing p by
# |p| in h changes only its constant sign; all logarithmic derivatives agree.
h=1/(rho**s.Rational(3,2)*p)
u=psi/h
replace_x={s.diff(x,y,2):s.diff(rho,y)/rho**2,s.diff(x,y):-1/rho}
tau=(-s.diff(alpha*s.diff(u,y),y)+q*u)/w
transformed=s.simplify((h*tau).subs(replace_x))
logder=s.simplify(-rho*s.diff(h,y)/h)
V=s.simplify(q/w-rho*s.diff(logder,y)+logder**2)
expected=-s.diff(F(x),x,2)+V*psi
checks={}
checks['unitary_bulk_measure']=s.simplify(h**2-w*rho)==0
checks['liouville_no_first_derivative']=s.simplify(transformed-expected)==0
checks['liouville_log_derivative']=s.simplify(logder-(s.Rational(3,2)*s.diff(rho,y)+rho*s.diff(p,y)/p))==0
flux=s.simplify((alpha*s.diff(u,y)).subs(replace_x))
checks['shell_flux_sign']=s.simplify(flux-h*(-s.diff(F(x),x)+logder*psi))==0

A,L,Ap,Lp,P,R,Upp=s.symbols('a ell a_prime ell_prime p rho U_second', real=True)
V_coeff=R**2*(s.Rational(3,4)*A**2+2*A*L+L**2-s.Rational(3,2)*Ap-Lp+s.Rational(2,3)*P**2)
V_background=s.expand(V_coeff.subs({Lp:Upp-4*Ap-4*A*L-L**2}).subs(Ap,-1/R**2-P**2/3))
checks['background_potential_simplification']=s.simplify(V_background-(-s.Rational(5,2)+R**2*(s.Rational(3,4)*A**2+6*A*L+2*L**2-Upp-P**2/6)))==0
Ac,C=s.symbols('A_c C', nonzero=True, real=True)
cone={rho:y,p:Ac*y}
cone_V=s.simplify(V.subs(cone).doit())
checks['cone_essential_threshold']=s.limit(cone_V,y,0)==s.Rational(25,4)
checks['shifted_threshold']=s.Rational(25,4)-4==s.Rational(9,4)
U1=s.symbols('U_prime', nonzero=True, real=True)
checks['jacobi_normalization_cancellation']=s.simplify((-U1/60)**2/(U1/5)**2)==s.Rational(1,144)

g=s.Function('g')(y)
checks['power_hardy_integration_identity']=s.simplify(s.diff(y**3*g**2,y)-3*y**2*g**2-2*y**3*g*s.diff(g,y))==0
# Check the Riccati/Picone integrand identity independently with generic P0.
P0=s.Function('P0')(y)
f=s.Function('f')(y)
a=s.Function('alpha')(y)
ww=s.Function('w')(y)
qq=4*ww+s.diff(a*s.diff(P0,y),y)/P0
lhs=a*s.diff(P0*f,y)**2+(qq-4*ww)*(P0*f)**2
rhs=a*P0**2*s.diff(f,y)**2+s.diff(a*P0*s.diff(P0,y)*f**2,y)
checks['picone_integrand_identity']=s.simplify(lhs-rhs)==0

out=Path(__file__).with_name('EXACT_ALGEBRA_CHECKS.json')
result={
 'scope':'Exact algebra only; closure, self-adjointness, Hardy coercivity and essential-spectrum arguments are mathematical proofs in the accompanying note.',
 'checks':checks,
 'passed':sum(checks.values()),
 'failed':sum(not v for v in checks.values()),
 'cone_potential_exact_model':str(cone_V),
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'sympy_version':s.__version__,
}
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
if not all(checks.values()):
    raise SystemExit(1)
