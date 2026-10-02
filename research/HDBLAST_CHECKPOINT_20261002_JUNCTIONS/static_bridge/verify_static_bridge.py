#!/usr/bin/env python3
"""Exact algebra and reanalysis of pinned archived static data; no new ODE/BVP.

Run from any directory. Outputs are confined to this scratch staging directory.
The sources are hash-pinned on every run. No quantum expectation is evaluated.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import numpy as np
import sympy as s

HERE = Path(__file__).resolve().parent
PINNED_COMMIT = 'e17a01b428bb8049e919c42376ab0359e41d613c'
MODEL_REL = Path('hdblast/checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922')
def detect_repo():
    for root in [Path.cwd(), *Path.cwd().parents, HERE, *HERE.parents]:
        for candidate in (root, root/'HDblast'):
            if (candidate/MODEL_REL/'static_branch/PLUS_BRANCH_RESULTS.json').is_file():
                return candidate
    raise SystemExit('Cannot find HDBLAST sources. Pass --repo PATH to a checkout or extracted source tree.')

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo', type=Path, help='HDBLAST checkout or extracted source tree; otherwise detect from current/script ancestors')
parser.add_argument('--output', type=Path, default=HERE, help='output directory; default is the script directory')
args = parser.parse_args()
REPO = args.repo.resolve() if args.repo else detect_repo().resolve()
OUTPUT = args.output.resolve()
OUTPUT.mkdir(parents=True, exist_ok=True)
MODEL = REPO/'hdblast/checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922'
GI = REPO/'research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant'
sources = [MODEL/'matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md',
           MODEL/'static_branch/solve_plus_branch.py',
           MODEL/'static_branch/PLUS_BRANCH_RESULTS.json',
           MODEL/'static_branch/PLUS_BRANCH_PROFILES.npz',
           MODEL/'static_branch/INDEPENDENT_BRANCH_REVIEW.md',
           GI/'gi_core.py', GI/'spectrum.py', GI/'SPECTRUM_RESULTS.json',
           GI/'README.md', GI/'VERIFICATION.md',
           REPO/'research/HDBLAST_CHECKPOINT_20260927/RESULTS_SUMMARY.json',
           REPO/'research/HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW/COMMON_ACTION_SMOOTH_FRW.md']
source_pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
if source_pins['source_commit'] != PINNED_COMMIT:
    raise SystemExit('SOURCE_PINS.json identifies an unexpected source commit')
provenance = {}
for path in sources:
    rel = path.relative_to(REPO).as_posix()
    actual = {'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size}
    if actual != source_pins['sources'].get(rel):
        raise SystemExit('Pinned source mismatch: '+rel)
    provenance[rel] = actual
checks = []
negative = []

def check(name, condition):
    ok = bool(condition)
    checks.append({'name': name, 'passed': ok})
    if not ok:
        raise AssertionError(name)

def zero(expr):
    return s.simplify(expr) == 0

phi, delta, c, kap, dens, pressure, current = s.symbols('phi delta c kap dens pressure current')
W = 1-phi+phi**3/3
U = s.diff(W,phi)**2/2-s.Rational(2,3)*W**2
sig = 2*W+delta*(1+c*phi)
check('model_plus_vacuum_U', U.subs(phi,1) == -s.Rational(2,27))
check('model_plus_vacuum_U_phi', s.diff(U,phi).subs(phi,1) == 0)
check('model_plus_vacuum_U_phiphi', s.diff(U,phi,2).subs(phi,1) == s.Rational(28,9))
check('constant_phi_failure_delta_c_over_two', s.diff(sig,phi).subs(phi,1)/2 == delta*c/2)
ks = (sig+kap*dens)/6
k0 = (sig-kap*(2*dens+3*pressure))/6
check('static_metric_compatibility_p_minus_rho', zero((ks-k0).subs(pressure,-dens)))
q,w,u = s.symbols('q w u')
constraint = q**2-w**2/12+u/6
static_h2 = (sig+kap*dens)**2/36-(s.diff(sig,phi)+kap*current)**2/48+U/6
check('static_constraint_with_both_sources', zero(constraint.subs({q:ks,w:-(s.diff(sig,phi)+kap*current)/2,u:U})-static_h2))

# Exact moving-endpoint column of the vacuum junction Jacobian.
sig0,sig1,sig2,u1,h2 = s.symbols('sig0 sig1 sig2 u1 h2')
q0,w0 = sig0/6,-sig1/2
q_y = -w0**2/4-u/6-q0**2
e1_y = q_y-sig1*w0/6
vac_h2 = sig0**2/36-sig1**2/48+u/6
check('jacobian_endpoint_metric_column_is_minus_H2', zero(e1_y+vac_h2))
e2_y = u1-4*q0*w0+sig2*w0/2
check('jacobian_endpoint_scalar_column_is_B_w', zero(e2_y-(u1/w0-4*q0+sig2/2)*w0))

# Lorentzian local action mapped to its Euclidean density W_E.
H,x,alpha,gamma = s.symbols('H x alpha gamma', nonzero=True)
V,F = s.Function('V')(x),s.Function('F')(x)
Vol = 8*s.pi**2/(3*H**4)
WE = V-12*F*H**2-144*alpha*H**4-24*gamma*H**4
rhoE = s.simplify(-H*s.diff(Vol*WE,H)/(4*Vol))
QE = 2*s.diff(WE,x)
check('local_static_density', zero(rhoE-(V-6*F*H**2)))
check('local_static_mass_source', zero(QE-(2*s.diff(V,x)-24*H**2*s.diff(F,x))))
check('constant_R2_Euler_no_static_sources', not rhoE.has(alpha,gamma) and not QE.has(alpha,gamma))
check('static_common_action_cross_derivative', zero(s.diff(rhoE,x)-QE/2+H*s.diff(QE,H)/8))
generalWE=s.Function('WE')(H,x)
generalrho=generalWE-H*s.diff(generalWE,H)/4
generalQ=2*s.diff(generalWE,x)
check('general_static_common_action_cross_derivative', zero(s.diff(generalrho,x)-generalQ/2+H*s.diff(generalQ,H)/8))
gbar,phistar,m0=s.symbols('gbar phistar m0')
mass=m0**2+gbar**2*(phi-phistar)**2
check('paired_scalar_source_xphi_Q_over_two', zero(s.diff(mass,phi)/2-gbar**2*(phi-phistar)))

# Continuum trace check specialized from the inherited covariant a2(r).
h=s.symbols('h')
a2=(h-2*H**2)**2/2-H**4/15
anomaly=a2/(16*s.pi**2)
check('static_heat_kernel_trace_anomaly', zero(anomaly-(h**2/(32*s.pi**2)-h*H**2/(8*s.pi**2)+29*H**4/(240*s.pi**2))))
ell,lam,Lambda,r=s.symbols('ell lam Lambda r', positive=True)
deg=(ell+1)*(ell+2)*(2*ell+3)/6
check('S4_low_harmonic_degeneracies', deg.subs(ell,0)==1 and deg.subs(ell,1)==5)
weights=[1,-3,3,-1]
for n in range(3):
    check(f'PV_moment_{n}', sum(weights[j]*(x+j*Lambda**2)**n for j in range(4)).expand()==0)
    check(f'PV_reference_moment_{n}', sum(weights[j]*(r+j*Lambda**2)**n for j in range(4)).expand()==0)
check('PV_log_leading_tail_lambda_minus_three', zero(sum(weights[j]*(x+j*Lambda**2)**3/3 for j in range(4))+2*Lambda**6))
check('PV_Q_leading_tail_lambda_minus_four', zero(-sum(weights[j]*(x+j*Lambda**2)**3 for j in range(4))-6*Lambda**6))

# Formal small-detuning/small-source expansion. Frozen rho,J are shell values.
eps,tau,nu,d,a=s.symbols('eps tau nu d a')
eta=a*eps
scalar_lead=(s.Rational(14,9)+2)*a+nu/2
asol=-9*nu/64
check('regular_AdS_scalar_leading_match', zero(scalar_lead.subs(a,asol)))
metric= s.Rational(2,3)+eps*tau+eps*d*eta+2*eta**2+s.Rational(2,3)*eta**3
scalar=eps*nu+4*eta+2*eta**2
Ueta=U.subs(phi,1+eta)
expanded=s.series(metric**2/36-scalar**2/48+Ueta/6,eps,0,3).removeO().expand().subs(a,asol)
expected=eps*tau/27+eps**2*(tau**2/36-d*nu/192+nu**2/384)
check('frozen_source_H2_second_order', zero(expanded-expected))
check('source_free_registered_second_order', zero((expected.subs({nu:d})/eps).subs(eps,1)-(tau/27+tau**2/36-d**2/384)))
t=s.symbols('t')
potential_paired=expected+eps**2*t*asol/27
check('action_paired_linear_potential_second_order', zero(potential_paired.subs(d,nu-t)-(eps*tau/27+eps**2*(tau**2/36-nu**2/384))))
check('frozen_current_linear_cross_term_cancels', zero(s.diff((expected.subs(nu,d+t)).coeff(eps,2),t).subs(t,0)))

def reject(name,residual):
    detected=not zero(residual)
    negative.append({'name':name,'detected':detected,'residual':str(s.simplify(residual))})
    if not detected:
        raise AssertionError(name)
reject('drop_scalar_source_factor_half', s.diff(mass,phi)-gbar**2*(phi-phistar))
reject('drop_Fphi_current', s.diff(rhoE,x)-(2*s.diff(V,x))/2+H*s.diff(2*s.diff(V,x),H)/8)
correct_alpha_rho = -H*s.diff(Vol*(-144*alpha*H**4),H)/(4*Vol)
reject('use_R2_as_static_vacuum_density', 144*alpha*H**4-correct_alpha_rho)
reject('freeze_density_when_absorbing_linear_potential', (expected-potential_paired).coeff(eps,2))
reject('flip_scalar_source_sign', scalar_lead.subs(a,9*nu/64))

archive=json.loads((GI/'SPECTRUM_RESULTS.json').read_text())
row=next(a for a in archive['details'] if a['branch']=='plus' and a['delta']==0.001)
J=np.array(row['static_jacobian']['J'],dtype=float)
forcing=np.diag([1/6,-1/2])
susceptibility=np.linalg.solve(J,forcing)
check('archived_Jacobian_nonsingular_float', np.linalg.det(J)!=0 and np.linalg.cond(J)<20)
check('archived_inverse_algebra_reconstruction', np.max(np.abs(J@susceptibility-forcing))<2e-15)
check('archived_endpoint_metric_column_consistency', abs(J[0,1]+row['shell']['H2_brane'])<1e-10)
check('archived_endpoint_scalar_column_consistency', abs(J[1,1]-row['shell']['B']*row['shell']['phi_y_b'])<1e-10)
try:
    head_result = subprocess.run(['git','-C',str(REPO),'rev-parse','HEAD'],text=True,capture_output=True)
    head = head_result.stdout.strip() if head_result.returncode == 0 else None
except FileNotFoundError:
    head = None
out={'status':'EXACT_ALGEBRA_AND_ARCHIVED_DATA_REANALYSIS_ONLY',
     'source_commit':PINNED_COMMIT, 'observed_git_head':head,
     'all_source_hashes_match_pin':True,'checks':checks,'negative_controls':negative,
     'archived_response':{'coordinates':['log10(abs(eta_h))','y_b'],
                          'forcing_columns':['kappa5_squared*rho','kappa5_squared*J_phi'],
                          'vacuum_delta':row['delta'],'J_raw':J.tolist(),
                          'determinant':float(np.linalg.det(J)),
                          'condition_number':float(np.linalg.cond(J)),
                          'du_d_source':susceptibility.tolist(),
                          'unresolved_components':['dy_b_d_kappa_squared_J_phi: numerical-floor dominated; driven by archived J[0,0]'],
                          'source_row_shell':row['shell'],
                          'scope':'Matrix inversion of archived vacuum data only; no state expectation, sensitivity integration, mode evolution, or BVP.'},
     'sources':provenance}
(OUTPUT/'STATIC_BRIDGE_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'checks_passed':len(checks),
                  'negative_controls_detected':len(negative),'source_commit':PINNED_COMMIT,
                  'all_source_hashes_match_pin':True,
                  'archived_jacobian_condition':out['archived_response']['condition_number']},indent=2))
