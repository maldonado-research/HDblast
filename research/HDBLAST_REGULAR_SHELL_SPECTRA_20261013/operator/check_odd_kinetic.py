"""Second-order moving-boundary EH+GHY kinetic audit in a homogeneous chart.

Direct level-set GHY density, radial ADM divergence flux, moving integration
limit, and shell volume are expanded independently. No imported project code.
Time test functions may have compact support, so time-corner variations vanish.
"""
from pathlib import Path
import argparse
import importlib.util
_support_path = Path(__file__).resolve().with_name('portable_support.py')
_support_spec = importlib.util.spec_from_file_location('hdblast_operator_portable_support', _support_path)
if _support_spec is None or _support_spec.loader is None:
    raise RuntimeError('Unable to load the bundled portable receipt support')
_support_module = importlib.util.module_from_spec(_support_spec)
_support_spec.loader.exec_module(_support_module)
write_result = _support_module.write_result
import hashlib
import json
import sympy as s

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=HERE/'results/odd-kinetic-check.json',
                    help='Deterministic JSON receipt path (default: beside script under results).')
args=parser.parse_args()
y,t,e,z=s.symbols('y t epsilon zeta')
q=s.Function('Q')(t)
rho=s.Function('rho')(y)
nf=s.Function('N')(y)
ps=s.Function('psi')(y)
xf=s.Function('X')(y)
sig=s.Function('sigma')(y)
R,a,N,psi,X,sigma=s.symbols('rho a N psi X sigma')
Q,Qt,Qtt=s.symbols('Q Qt Qtt')
checks={}
record={}

def coeff(expr,k):
    return s.simplify(s.diff(expr,e,k).subs(e,0)/s.factorial(k))
def at_graph_second(expr):
    return s.simplify(coeff(expr,2)+z*q*s.diff(coeff(expr,1),y)
        +z*z*q*q*s.diff(coeff(expr,0),y,2)/2)
def jets(expr):
    return s.expand(expr).subs({rho:R,s.diff(rho,y):a*R,nf:N,ps:psi,
        xf:X,sig:sigma,q:Q,s.diff(q,t):Qt,s.diff(q,t,2):Qtt},simultaneous=True)
def time_ibp(expr):
    """Density divided by sqrt(-gamma)=exp(3t); remove Qtt exactly modulo a
    time derivative. All examples have its coefficient proportional to Q."""
    ex=s.expand(jets(expr))
    f=s.diff(ex,Qtt)
    if f.has(Qtt,Qt):
        raise ValueError('The time-IBP coefficient has unexpected derivatives')
    return s.simplify(ex-f*Qtt-(s.diff(f,Q)*Qt+3*f)*Qt)
def kinetic_time_coefficient(expr):
    return s.simplify(s.diff(time_ibp(expr),Qt,2)/2)
def ok(name,expr):
    v=s.simplify(expr)
    checks[name]=bool(v==0)
    if v!=0: record['FAILED_'+name]=str(v)
    print(('PASS ' if v==0 else 'FAIL ')+name,flush=True)
def detect(name,expr):
    checks[name]=bool(s.simplify(expr)!=0)
    print(('PASS ' if checks[name] else 'FAIL ')+name,flush=True)

# Exact nonlinear completion in longitudinal gauge, used only to order e^2.
n=s.exp(e*nf*q)
A2=rho**2*s.exp(2*e*ps*q)
A4=A2**2
f_t=e*z*s.diff(q,t)
# D=1-n^2 f_t^2/A^2. Its square root and inverse to the required order.
sqrtD=1-e**2*z**2*s.diff(q,t)**2/(2*rho**2)
invD=1+e**2*z**2*s.diff(q,t)**2/(2*rho**2)

# sqrt(-h_graph) K_graph / exp(3t), from K=div(normal), with
# n^y=1/(n sqrtD), n^t=n f_t/(A^2 sqrtD).
ghy=sqrtD/n*(s.diff(A4*invD,y)
    +s.diff(n**2*A2*f_t*invD,t)+3*n**2*A2*f_t*invD)
radial_ghy=s.diff(A4,y)/n
# EH decomposition contains -d_mu(sqrt(-h) h^{mu nu} d_nu n).
# Integrating on y<f(t) gives +sqrt(-h) h^{tt} n_t f_t at its upper edge.
eh_tangential_flux=-A2*s.diff(n,t)*f_t
graph_correction=at_graph_second(ghy-radial_ghy+eh_tangential_flux)
graph_no_flux=at_graph_second(ghy-radial_ghy)
graph_kin=kinetic_time_coefficient(graph_correction)
ok('GHY_graph_plus_EH_flux_one_side',graph_kin+3*a*R**2*z**2)
ok('all_lapse_and_conformal_cross_terms_cancel_in_graph_correction',
   s.diff(graph_kin,N)**2+s.diff(graph_kin,psi)**2)
detect('CONTROL_omit_EH_tangential_divergence_detected',
       kinetic_time_coefficient(graph_no_flux)-graph_kin)

# Intrinsic R[h] from conformal transformation; this density is the radial ADM
# term, before any time integration by parts on the moving integration domain.
Rh=(12+6*e*ps*(s.diff(q,t,2)+3*s.diff(q,t))
       +6*e**2*ps**2*s.diff(q,t)**2)/A2
adm_intrinsic=n*A4*Rh/2
adm_scalar_tangent=n*A2*e**2*xf**2*s.diff(q,t)**2/2
bulk_kin=kinetic_time_coefficient(coeff(adm_intrinsic+adm_scalar_tangent,2))
ok('direct_ADM_bulk_kinetic_before_constraint',
       bulk_kin-R**2*(-3*N*psi-3*psi**2+X**2/2))
ok('direct_ADM_bulk_kinetic_after_constraint',
       bulk_kin.subs(N,-2*psi)-R**2*(3*psi**2+X**2/2))
moving_domain=z*q*coeff(adm_intrinsic,1)
domain_kin=kinetic_time_coefficient(moving_domain)
ok('moving_ADM_domain_first_order_curvature_one_side',domain_kin+3*R**2*z*psi)

# Potential dependence of sigma(phi) cannot affect derivative terms at this
# order: its linear perturbation only multiplies q, while sqrtD-1 is order 2.
shell=-sig*A4*sqrtD
shell_kin=kinetic_time_coefficient(at_graph_second(shell))
ok('shell_volume_kinetic',shell_kin-sigma*R**2*z**2/2)
total_boundary_kin=2*(graph_kin+domain_kin)+shell_kin
ok('full_two_side_boundary_before_odd_continuity',
   total_boundary_kin-R**2*(-6*a*z*z-6*z*psi+sigma*z*z/2))
odd_kin=s.simplify(total_boundary_kin.subs({psi:-a*z,sigma:6*a}))
ok('full_odd_boundary_time_kinetic_coefficient',odd_kin-3*a*R**2*z**2)
P=s.symbols('P')
ok('full_odd_boundary_physical_K_is_18cP_squared',
       (2*odd_kin).subs(z,-P/(a*R**2))-6*P**2/(a*R**2))
detect('CONTROL_omit_moving_domain_detected',
       (2*graph_kin+shell_kin).subs({psi:-a*z,sigma:6*a})-odd_kin)

record.update({
 'graph_correction_density_after_time_IBP':str(time_ibp(graph_correction)),
 'graph_time_kinetic_one_copy':str(graph_kin),
 'moving_domain_time_kinetic_one_copy':str(domain_kin),
 'shell_volume_time_kinetic_one_shell':str(shell_kin),
 'two_side_boundary_time_kinetic_before_continuity':str(total_boundary_kin),
 'odd_boundary_K':str(2*odd_kin),
})
out={'date_utc':'2026-10-13','sympy':s.__version__,
 'n_checks':len(checks),'n_pass':sum(checks.values()),'checks':checks,'expressions':record,
 'scope':'Second-order local kinetic action including moving-domain, EH divergence, graph GHY, and shell. Compact-support time variations avoid corner terms; global symplectic completeness is not proved.'}
raise SystemExit(0 if write_result(out, 11, args.output, Path(__file__)) else 1)
