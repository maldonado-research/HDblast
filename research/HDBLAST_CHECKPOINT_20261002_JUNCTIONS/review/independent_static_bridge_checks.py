#!/usr/bin/env python3
"""Independent static bridge algebra and archived matrix audit; no integration."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo', type=Path, default=Path('/workspace/HDblast'))
args = parser.parse_args()
REPO = args.repo.resolve()
records, negatives = [], []


def exact(name, expression):
    z = sp.factor(sp.expand(expression))
    if z != 0:
        raise AssertionError(f'{name}: {z}')
    records.append({'name':name, 'passed':True, 'residual':'0'})


def reject(name, expression):
    z = sp.factor(sp.expand(expression))
    if z == 0:
        raise AssertionError(f'Undetected mutation: {name}')
    negatives.append({'name':name, 'rejected':True, 'residual':str(z)})


z = sp.symbols('cosh_y_over_9', real=True)
f = sp.gegenbauer(14,2,z)/sp.gegenbauer(14,2,1)
exact('regular cone polynomial solves linear scalar equation',
      (z*z-1)*sp.diff(f,z,2)+5*z*sp.diff(f,z)-252*f)
exact('regular cone polynomial unit vertex normalization', f.subs(z,1)-1)
growth = sp.Rational(14,9)
exact('regular AdS asymptotic scalar logarithmic derivative', growth**2+sp.Rational(4,9)*growth-sp.Rational(28,9))

epsilon, tau, nu, d, b, a, source_j = sp.symbols('epsilon tau nu d eta2 eta1 source_j')
eta = a*epsilon+b*epsilon**2
phi = 1+eta
W = 1-phi+phi**3/3
U = (phi**2-1)**2/2-sp.Rational(2,3)*W**2
sigma_rho = 2*W+epsilon*tau+epsilon*d*eta
sigma_phi_j = 2*(phi**2-1)+epsilon*nu
H2 = sigma_rho**2/36-sigma_phi_j**2/48+U/6
expanded = sp.series(H2,epsilon,0,3).removeO().expand()
expected_unmatched = epsilon*tau/27+epsilon**2*(tau**2/36+d*a/27-nu*a/6-nu**2/48)
exact('static H squared expansion including unknown second scalar coefficient', expanded-expected_unmatched)
exact('unknown second scalar coefficient cancels through second order', sp.diff(expanded,b))
first_scalar = (growth+2)*a+nu/2
a_solution = -9*nu/64
exact('linear normal scalar matching', first_scalar.subs(a,a_solution))
frozen = expanded.subs(a,a_solution)
frozen_expected = epsilon*tau/27+epsilon**2*(tau**2/36-d*nu/192+nu**2/384)
exact('frozen source expansion', frozen-frozen_expected)
potential = frozen+epsilon**2*source_j*a_solution/27
exact('paired linear potential special result', potential.subs(d,nu-source_j)
      -(epsilon*tau/27+epsilon**2*(tau**2/36-nu**2/384)))
potential_formula = epsilon*tau/27+epsilon**2*(tau**2/36-nu**2/384)
reject('apply potential formula to frozen shell density', frozen-potential_formula)

H,x = sp.symbols('H x', positive=True)
alpha,gamma = sp.symbols('alpha gamma', real=True)
V,F = sp.Function('V')(x),sp.Function('F')(x)
volume = 8*sp.pi**2/(3*H**4)
WE = V-12*F*H**2-144*alpha*H**4-24*gamma*H**4
rho = sp.simplify(-H*sp.diff(volume*WE,H)/(4*volume))
Q = 2*sp.diff(WE,x)
exact('Euclidean metric scaling reproduces Lorentzian local density', rho-(V-6*F*H**2))
exact('Euclidean mass variation reproduces local variance', Q-(2*sp.diff(V,x)-24*H**2*sp.diff(F,x)))
exact('static mixed source identity', sp.diff(rho,x)-Q/2+H*sp.diff(Q,H)/8)
reject('replace local current by fixed-H density derivative', Q/2-sp.diff(rho,x))
G = sp.Function('WE_general')(H,x)
rho_general = G-H*sp.diff(G,H)/4
Q_general = 2*sp.diff(G,x)
exact('static mixed source identity for arbitrary common W', sp.diff(rho_general,x)-Q_general/2+H*sp.diff(Q_general,H)/8)

h = sp.symbols('mass_offset', real=True)
R = 12*H**2
riemann2, ricci2 = 24*H**4,36*H**4
a2 = (h-R/6)**2/2+(riemann2-ricci2)/180
exact('reference de Sitter heat-kernel anomaly', a2/(16*sp.pi**2)
      -(h*h/2-2*h*H**2+sp.Rational(29,15)*H**4)/(16*sp.pi**2))

# Endpoint derivative checks derive from the continuum radial equations.
sig0,sig1,sig2,bulkU,bulkUphi = sp.symbols('sigma sigma_phi sigma_phiphi U U_phi')
q, w = sig0/6,-sig1/2
qy = -w*w/4-bulkU/6-q*q
wy = bulkUphi-4*q*w
e1y = qy-sig1*w/6
e2y = wy+sig2*w/2
H2_boundary = q*q-w*w/12+bulkU/6
exact('vacuum endpoint metric Jacobian column', e1y+H2_boundary)
exact('vacuum endpoint scalar Jacobian column', e2y-(bulkUphi-4*q*w+sig2*w/2))

archive = REPO/'research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant/SPECTRUM_RESULTS.json'
rows=json.loads(archive.read_text())
row=next(z for z in rows['details'] if z['branch']=='plus' and z['delta']==.001)
matrix=np.array(row['static_jacobian']['J'])
forcing=np.diag([1/6,-1/2])
response=np.linalg.solve(matrix,forcing)
residual=np.max(np.abs(matrix@response-forcing))
if residual >= 2e-15:
    raise AssertionError(f'archived matrix inversion: {residual}')
records.append({'name':'archived matrix algebra reconstruction', 'passed':True, 'max_absolute_residual':float(residual)})
endpoint=np.array([-row['shell']['H2_brane'],row['shell']['B']*row['shell']['phi_y_b']])
raw_difference=2e-6*matrix[0,0]
out={'scope':'independent static algebra and archived matrix reanalysis; no BVP, spectral sum or quantum expectation evaluated',
     'checks_passed':len(records),'negative_controls_rejected':len(negatives),
     'checks':records,'negative_controls':negatives,
     'archived_matrix':{'source':str(archive),'source_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
                        'coordinates':['log10(abs(eta_h))','y_b'],
                        'centered_difference_step':1e-6,'J':matrix.tolist(),
                        'forcing_diagonal':[1/6,-1/2],'response_matrix':response.tolist(),
                        'condition_number':float(np.linalg.cond(matrix)),
                        'endpoint_column_absolute_error':(matrix[:,1]-endpoint).tolist(),
                        'J00_underlying_residual_difference':float(raw_difference),
                        'roundoff_warning':'Tiny y_b response to scalar source is driven entirely by tiny J00; no component error certificate or derivative-refinement study is supplied.'},
     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'sympy_version':sp.__version__}
(HERE/'INDEPENDENT_STATIC_BRIDGE_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
print(f'{len(records)} static checks passed; {len(negatives)} negative controls rejected')
