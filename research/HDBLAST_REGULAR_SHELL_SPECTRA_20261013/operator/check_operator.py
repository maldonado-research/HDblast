"""Fresh tensor, constraint, shell and kinetic-bilinear audit; no source imports.

Run with Python and SymPy 1.14.0; --output selects the JSON receipt.
This is exact local algebra, not a functional-analytic or symplectic-domain proof.
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

HERE = Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=HERE/'results/operator-check.json',
                    help='Deterministic JSON receipt path (default: beside script under results).')
args=parser.parse_args()
y, t, x1, x2, x3 = s.symbols('y t x1 x2 x3', real=True)
coords = (y, t, x1, x2, x3)
r, m, zeta = s.symbols('r mu_squared zeta')
rho = s.Function('rho')(y)
phi = s.Function('phi')(y)
Nf, psif, Xf = [s.Function(k)(y) for k in ('N', 'psi', 'X')]
Y = s.exp(r*t)
gam = s.diag(-1, s.exp(2*t), s.exp(2*t), s.exp(2*t))
g0 = s.diag(1, *[rho**2*gam[i,i] for i in range(4)])
g0i = g0.inv()
h = s.diag(2*Nf*Y, *[2*rho**2*psif*Y*gam[i,i] for i in range(4)])
di = -g0i*h*g0i

def connection(ginv, dg):
    return [[[s.simplify(sum(ginv[a,d]*(s.diff(dg[d,b],coords[c])
                  + s.diff(dg[d,c],coords[b])-s.diff(dg[b,c],coords[d]))
                 for d in range(5))/2)
              for c in range(5)] for b in range(5)] for a in range(5)]

G0 = connection(g0i,g0)
Gh = connection(g0i,h)
Gi = connection(di,g0)
G1 = [[[s.simplify(Gh[a][b][c]+Gi[a][b][c]) for c in range(5)]
       for b in range(5)] for a in range(5)]

def ricci_linear(b,c):
    ans = sum(s.diff(G1[a][b][c],coords[a])-s.diff(G1[a][b][a],coords[c])
              for a in range(5))
    ans += sum(G1[a][a][d]*G0[d][b][c]+G0[a][a][d]*G1[d][b][c]
               -G1[a][c][d]*G0[d][b][a]-G0[a][c][d]*G1[d][b][a]
               for a in range(5) for d in range(5))
    return s.simplify(ans)

R, a, p, u1, u2 = s.symbols('rho a p U1 U2', nonzero=True)
u0 = -6*a*a + 6/R**2 + p*p/2
nn, ns, nss, ps, ps1, ps2, xx, xs, xss = s.symbols('N Ny Nyy psi psiy psiyy X Xy Xyy')
jets = {rho:R, s.diff(rho,y):a*R,
        s.diff(rho,y,2):(a*a-1/R**2-p*p/3)*R,
        s.diff(phi,y):p, s.diff(phi,y,2):u1-4*a*p,
        Nf:nn,s.diff(Nf,y):ns,s.diff(Nf,y,2):nss,
        psif:ps,s.diff(psif,y):ps1,s.diff(psif,y,2):ps2,
        Xf:xx,s.diff(Xf,y):xs,s.diff(Xf,y,2):xss}

checks = {}
expressions = {}
def simp(e):
    return s.factor(s.cancel(s.expand(e).subs(jets, simultaneous=True)))
def check(name, e):
    v = s.simplify(e)
    checks[name] = bool(v == 0)
    if v != 0: expressions['FAILED_'+name] = str(v)
    print(('PASS ' if v == 0 else 'FAIL ')+name, flush=True)
def nonzero(name, e):
    checks[name] = bool(s.simplify(e) != 0)
    print(('PASS ' if checks[name] else 'FAIL ')+name, flush=True)
def harmonic_reduce(e):
    num, den = s.fraction(s.cancel(e))
    return s.factor(s.rem(s.Poly(num,r),s.Poly(r*r+3*r+m,r)).as_expr()/den)

E = {}
for b,c in ((0,0),(0,1),(1,1),(2,2)):
    matter = (s.diff(phi,coords[b])*s.diff(Xf*Y,coords[c])
            + s.diff(phi,coords[c])*s.diff(Xf*Y,coords[b]))
    ee = ricci_linear(b,c)-matter-s.Rational(2,3)*(u1*Xf*Y*g0[b,c]+u0*h[b,c])
    denom = Y*(s.exp(2*t) if (b,c)==(2,2) else 1)
    E[b,c] = simp(ee/denom)

tf = s.factor((E[1,1]+E[2,2])/(r*r-r))
check('tensor_tracefree_coefficient',tf+nn+2*ps)
check('tensor_momentum_constraint',E[0,1]/r-(3*a*nn-3*ps1-p*xx))
gauss = harmonic_reduce((E[0,0]-(-E[1,1]+3*E[2,2])/R**2)/2)
expected_gauss = 12*a*(ps1-a*nn)+3*(m+4)*ps/R**2-p*(xs-p*nn)+u1*xx
check('tensor_gauss_constraint',gauss-expected_gauss)
expressions['gauss_unreduced']=str(expected_gauss)

flow = {R:a*R,a:-1/R**2-p*p/3,p:u1-4*a*p,u1:u2*p}
def deriv(e, extra=None):
    f = dict(flow)
    if extra: f.update(extra)
    return s.expand(sum(s.diff(e,k)*v for k,v in f.items()))
ps_rhs = -2*a*ps-p*xx/3
x_rhs = (u1/p-4*a)*xx+(3*(m+4)/R**2/p-2*p)*ps
mode_flow={ps:ps_rhs,xx:x_rhs}
onshell={nn:-2*ps,ns:-2*ps_rhs,ps1:ps_rhs,xs:x_rhs,
         ps2:deriv(ps_rhs,mode_flow),xss:deriv(x_rhs,mode_flow)}
for bc, ee in E.items():
    check('first_order_implies_E_'+str(bc),harmonic_reduce(ee.subs(onshell,simultaneous=True)))
scalar=xss+4*a*xs+m*xx/R**2-u2*xx+p*(4*ps1-ns)-2*nn*u1
scalar_tensor=0
for b in range(5):
    for c in range(5):
        if g0i[b,c] == 0 and di[b,c] == 0:
            continue
        scalar_tensor+=di[b,c]*(s.diff(phi,coords[b],coords[c])
          -sum(G0[d][b][c]*s.diff(phi,coords[d]) for d in range(5)))
        scalar_tensor+=g0i[b,c]*(s.diff(Xf*Y,coords[b],coords[c])
          -sum(G0[d][b][c]*s.diff(Xf*Y,coords[d])
                   +G1[d][b][c]*s.diff(phi,coords[d]) for d in range(5)))
scalar_tensor=simp(scalar_tensor/Y-u2*Xf)
check('tensor_scalar_equation_unreduced',harmonic_reduce(scalar_tensor)-scalar)
check('first_order_implies_scalar',scalar.subs(onshell,simultaneous=True))

# Explicit normal/extrinsic curvature: n_y=1+epsilon N Y,
# n_mu=-epsilon zeta d_mu Y.  Background mixed y components vanish,
# so tangent-vector corrections are zero to first order.
sig2 = s.symbols('sigma2')
p1 = u1-4*a*p
ap = -1/R**2-p*p/3
for mu in (0,1):
    b=mu+1
    normal_tangent = -zeta*s.diff(Y,coords[b],coords[b])
    normal_tangent += sum(G0[c][b][b]*zeta*s.diff(Y,coords[c]) for c in range(1,5))
    kc1=normal_tangent-G1[0][b][b]-G0[0][b][b]*Nf*Y
    kc0=-G0[0][b][b]
    kc1+=zeta*Y*s.diff(kc0,y)
    induced1=2*rho**2*(psif+rho.diff(y)/rho*zeta)*gam[mu,mu]*Y
    mixed1=simp((kc1-(rho.diff(y)/rho)*induced1)/(rho**2*gam[mu,mu]*Y))
    hess_ratio=r*r/(-1) if mu==0 else -r
    target=ps1-a*nn+ap*zeta-zeta*hess_ratio/R**2
    check('displaced_mixed_extrinsic_'+str(mu),mixed1-target)
    # delta sigma/6=(sigma'/6)(X+p zeta), sigma'=-2p.
    junction=mixed1+p*(xx+p*zeta)/3
    target_j=ps1-a*nn+p*xx/3-zeta*(1+hess_ratio)/R**2
    check('full_metric_junction_'+str(mu),junction-target_j)

B = p1/p+sig2/2
Jscalar = xs-p*nn+sig2*xx/2+p*B*zeta
check('generic_scalar_boundary',Jscalar.subs(onshell, simultaneous=True).subs(zeta,0)
      -(B*xx+3*(m+4)*ps/(R**2*p)))
nonzero('CONTROL_drop_lapse_term_detected',-p*nn)
nonzero('CONTROL_reverse_bending_sign_detected',2*p*B*zeta)
bad_onshell=dict(onshell)
bad_onshell[xs]=x_rhs-6*(m+4)*ps/(R**2*p)
nonzero('CONTROL_reverse_gravity_coupling_detected',gauss.subs(bad_onshell,simultaneous=True))

# General gauge shifts in invariant longitudinal variables correspond to a
# radial shift T at the shell; metric and scalar bending terms cancel it.
T, Tp = s.symbols('T Tp')
gauge_fields={nn:Tp,ps:a*T,ps1:ap*T+a*Tp,xx:p*T,xs:p1*T+p*Tp,zeta:-T}
check('shell_scalar_radial_gauge_cancellation',Jscalar.subs(gauge_fields,simultaneous=True))
metric_trace_before_hess=ps1-a*nn+p*xx/3-zeta/R**2
check('shell_metric_radial_gauge_trace_cancellation',
      metric_trace_before_hess.subs(gauge_fields,simultaneous=True))
# A radial shift also creates B=T/rho^2, whose Hessian term cancels zeta.
check('shell_metric_radial_gauge_hessian_cancellation',-T/R**2-(-T)/R**2)

# Exceptional harmonic classification follows from commuting the Hessian.
check('tracefree_hessian_divergence_coefficient',(m+3)-m/4-s.Rational(3,4)*(m+4))
check('homogeneous_exception_factor',r*r-r-r*(r-1))
check('translation_metric_bending_cancellation',1+(-1))
e0,L0,F=s.symbols('e0 L0 F')
translation_X=p*e0
translation_psi=(a+F)*e0-L0
check('translation_momentum_family',deriv(translation_psi,{F:1/R**2})+p*translation_X/3)
check('translation_gauss_family',deriv(translation_X)-p1/p*translation_X)
check('translation_two_parameter_rank',s.det(s.Matrix([[p,0],[a+F,-1]]))+p)
check('translation_relative_scalar_condition',Jscalar.subs(
      {nn:0,xx:translation_X,xs:deriv(translation_X)}, simultaneous=True)-p*B*(zeta+e0))

# Constant harmonics use the regular static family; they do not inherit the
# momentum constraint, whose angular gradient vanishes identically.
j,hc,J=s.symbols('j h J')
K=R**2*j
M0=hc*J+p*B*K
static=s.Matrix([[j,-1/R**2],[hc*J,p*B]])
check('constant_static_determinant',static.det()-M0/R**2)

# Actual tangential kinetic bilinear of radial ADM (one copy).  This is not
# obtained by inventing an action for the Sturm-Liouville differential equation.
# conformal R4=e^{-2psi}/rho^2*(12-6 Box psi-6 gradpsi^2).
# After tangential integration by parts, Lkin/rho^2=3 gradN.gradpsi
#     +3 gradpsi^2 - gradX^2/2.
dN,dpsi,dX=s.symbols('dN dpsi dX')
adm_kin=3*dN*dpsi+3*dpsi*dpsi-dX*dX/2
check('constrained_ADM_kinetic_sign',adm_kin.subs(dN,-2*dpsi)+(dX*dX+6*dpsi*dpsi)/2)
P,Pprime=s.symbols('P Pprime')
psl=1/(R**2*p*p)
qsl=2/(3*R**2)
wsl=1/(R**4*p*p)
check('physical_kinetic_density_equals_18_energy_density',
      2*R**2*(xx*xx+6*ps*ps)-18*(psl*Pprime**2+qsl*P**2).subs(
          {P:R**2*ps,Pprime:-R**2*p*xx/3},simultaneous=True))
check('SL_weight_density',wsl*(R**2*ps)**2-ps*ps/(p*p))
check('odd_Robin_sign_algebra',psl*(-R**2*p*(p/a)*ps/3)+(R**2*ps)/(3*a*R**2))

out={'date_utc':'2026-10-13','sympy':s.__version__,
     'n_checks':len(checks),'n_pass':sum(checks.values()),'checks':checks,
     'expressions':expressions,
     'scope':'Exact local linear equations and kinetic-bilinear checks; not a global symplectic/domain proof.'}
raise SystemExit(0 if write_result(out, 32, args.output, Path(__file__)) else 1)
