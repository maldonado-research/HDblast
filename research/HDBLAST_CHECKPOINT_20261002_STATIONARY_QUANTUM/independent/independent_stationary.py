#!/usr/bin/env python3
"""Independent acceleration-equation shooting. No producer modules imported.

The registered root and source experiment must be frozen before invoking main.
The callable source formula follows a resolvent finite part and direct metric
variation, independently of determinant-series evaluation. A separate legacy
proper-time source replay may verify selected endpoints after this root solve.
"""
import argparse
import datetime
import hashlib
import json
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.special import digamma
import scipy

DELTA = .001
DETUNING_C = .5975949350280132
INHERITED_ETA_H = -1.3366923651588084e-23
INHERITED_YB = 30.276680395885563
COEFFICIENTS = (-2/27, 0., 14/9, 50/27, -1/6, -4/9, -2/27)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def potential(eta, derivative=0):
    # Explicit differentiated Horner polynomial; avoid phi_h rounding to 1.
    value = 0.
    for j in range(6, derivative-1, -1):
        factor = math.factorial(j)/math.factorial(j-derivative)
        value = value*eta + COEFFICIENTS[j]*factor
    return value


def tension(eta):
    return 2/3 + 2*eta*eta + 2*eta**3/3 + DELTA*(1+DETUNING_C*(1+eta))


def tension_phi(eta):
    return 4*eta + 2*eta*eta + DELTA*DETUNING_C


def cone_state(ell, y):
    eta = -10.**ell
    U, Up, Upp = [potential(eta, k) for k in range(3)]
    a = -U/36
    bb = U*U/4320-Up*Up/750
    cc = Up/10
    dd = Up*(Upp/280+U/630)
    return [y+a*y**3+bb*y**5, 1+3*a*y*y+5*bb*y**4,
            eta+cc*y*y+dd*y**4, 2*cc*y+4*dd*y**3]


def acceleration(y, state):
    R, P, eta, w = state
    return [P, -R*(w*w/4+potential(eta)/6), w,
            potential(eta, 1)-4*P*w/R]


def sources(z, x, reference):
    require(z > 0 and x > 0 and reference > 0, 'Positive source domain required')
    u, v = x/z, reference/z
    nu = np.sqrt(complex(2.25-u))
    psi_sum = digamma(1.5+nu)+digamma(1.5-nu)
    require(abs(psi_sum.imag) < 1e-12, 'Unexpected imaginary resolvent')
    L = psi_sum.real-math.log(v)
    Q = z/(16*math.pi**2)*((u-2)*L-u+v+4/3)
    D = u*u/2-2*u+29/15
    # Direct radius variation of fully normalized W; no trace reconstruction.
    T = u*(u-2)*L-u*u+4*u/3+2*u*v-2*v-v*v/2-D
    rho = z*z/(64*math.pi**2)*T
    return rho, Q


def four_derivatives(z, x, reference, b):
    """The four derivatives of rho(z,phi), J(z,phi), fixed reference."""
    with mp.workdps(50):
        z, x, reference, b = map(mp.mpf, map(str, (z, x, reference, b)))
        u, v = x/z, reference/z
        nu = mp.sqrt(mp.mpf(9)/4-u)
        L = mp.digamma(mp.mpf(3)/2+nu)+mp.digamma(mp.mpf(3)/2-nu)-mp.log(v)
        if abs(nu) < mp.mpf('1e-20'):
            Pu = -mp.polygamma(2, mp.mpf(3)/2)
        else:
            Pu = (mp.polygamma(1, mp.mpf(3)/2-nu)-mp.polygamma(1, mp.mpf(3)/2+nu))/(2*nu)
        D = u*u/2-2*u+mp.mpf(29)/15
        T = u*(u-2)*L-u*u+4*u/3+2*u*v-2*v-v*v/2-D
        Tu = 2*(u-1)*L+u*(u-2)*Pu-3*u+mp.mpf(10)/3+2*v
        Tv = -u*(u-2)/v+2*u-2-v
        C = 1/(16*mp.pi**2)
        Q = C*z*((u-2)*L-u+v+mp.mpf(4)/3)
        Qx = C*(L+(u-2)*Pu-1)
        Qz = C*(-2*L-u*(u-2)*Pu+u-mp.mpf(2)/3)
        rhox = C*z*Tu/4
        rhoz = C*z*(2*T-u*Tu-v*Tv)/4
        xphi = b*mp.sqrt(reference*x)
        xphiphi = b*b*reference/2
        rhophi = xphi*rhox
        Jz = xphi*Qz/2
        Jphi = (xphiphi*Q+xphi*xphi*Qx)/2
        pairing = rhox-Q/2+z*Qz/4
        values = dict(rho_z=rhoz, rho_phi=rhophi, J_z=Jz, J_phi=Jphi,
                      Q_x=Qx, Q_z=Qz, common_action_pairing=pairing)
        return {k: float(mp.re(v)) for k, v in values.items()}


def background(ell, yb, setting):
    y0 = setting['y0']
    require(yb > y0 and math.isfinite(ell), 'Invalid shooting coordinates')
    result = solve_ivp(acceleration, (y0, yb), cone_state(ell, y0),
                       method='DOP853', rtol=setting['rtol'],
                       atol=[1e-14,1e-14,1e-100,1e-100],
                       max_step=setting['max_step'])
    require(result.success, result.message)
    R, P, eta, w = result.y
    require(np.all(R > 0), 'Interior radius ceased to be positive')
    U = np.array([potential(e) for e in eta])
    defect = P*P-1-R*R*(w*w/12-U/6)
    relative = np.abs(defect)/np.maximum(1, P*P)
    end = tuple(float(v[-1]) for v in (R,P,eta,w))
    return dict(state=end, Hamiltonian_max=float(np.max(np.abs(defect))),
                Hamiltonian_relative_max=float(np.max(relative)),
                Hamiltonian_endpoint=float(defect[-1]), steps=int(len(R)),
                nfev=result.nfev)


def endpoint(coordinates, gamma, b, reference, eta_ref, setting):
    back = background(*coordinates, setting)
    R,P,eta,w = back['state']
    z = 1/(R*R)
    mass_factor = 1+b*(eta-eta_ref)/2
    require(mass_factor > 0, 'Registered positive mass-factor domain violated')
    x = reference*mass_factor*mass_factor
    xphi = b*reference*mass_factor
    xphiphi = b*b*reference/2
    rho, Q = sources(z, x, reference)
    J = xphi*Q/2
    E1 = P/R-tension(eta)/6-gamma*rho/6
    E2 = w+tension_phi(eta)/2+gamma*J/2
    projected = (tension(eta)+gamma*rho)**2/36-(tension_phi(eta)+gamma*J)**2/48+potential(eta)/6
    anomaly = ((x-reference)**2/2-2*z*(x-reference)+29*z*z/15)/(16*math.pi**2)
    back.update(ell=float(coordinates[0]), y_b=float(coordinates[1]), eta=eta,
                phi=1+eta, R=R, P=P, w=w, H2=z, x=x, r=reference, gamma=gamma,
                b=b, mass_factor=mass_factor, x_phi=xphi, x_phiphi=xphiphi,
                rho=rho, p=-rho, Q=Q, J=J, E1=E1, E2=E2,
                projected_H2=projected, endpoint_constraint_residual=z-projected,
                trace_residual=-4*rho+x*Q-anomaly)
    return back


def solve_case(start, gamma, b, reference, eta_ref, setting):
    calls = []
    def objective(coordinates):
        data = endpoint(coordinates, gamma, b, reference, eta_ref, setting)
        calls.append({k:data[k] for k in ('ell','y_b','E1','E2','H2','eta')})
        return np.array([data['E1'],data['E2']])/DELTA
    t0 = time.monotonic()
    result = root(objective, start, method='hybr', options={'xtol':5e-11,'maxfev':120,'eps':1e-8})
    end = endpoint(result.x, gamma, b, reference, eta_ref, setting)
    end['four_source_derivatives'] = four_derivatives(end['H2'],end['x'],reference,b)
    controls = dict(reverse_current_defect=-gamma*end['J'],
                    omit_current_defect=-gamma*end['J']/2,
                    wrong_pressure_temporal_defect=gamma*end['rho'])
    # Wrong scalar variation when reference is allowed to track x.
    h=end['x']-reference; z=end['H2']
    Wr=-(h*h/2-2*z*h+29*z*z/15)/(32*math.pi**2*reference)
    controls['tracking_reference_scalar_defect']=gamma*Wr*end['x_phi']/2
    control_checks=[]
    for name,defect in controls.items():
        applicable=bool(gamma != 0 and (name == 'wrong_pressure_temporal_defect' or b != 0))
        control_checks.append(dict(name=name,applicable=applicable,
                                   normalized_defect=abs(defect)/DELTA,
                                   detected=bool(abs(defect)/DELTA > 2e-12) if applicable else None))
    end.update(root_success=bool(result.success), root_message=str(result.message),
               root_nfev=int(result.nfev), root_iterations=calls, seconds=time.monotonic()-t0,
               negative_controls=controls, negative_control_checks=control_checks, setting=setting,
               residual_gate=bool(max(abs(end['E1']),abs(end['E2']))/DELTA <= 2e-12),
               Hamiltonian_gate=bool(end['Hamiltonian_relative_max'] <= 2e-10))
    return end


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--registration-sha256',required=True)
    parser.add_argument('--amplitudes',required=True,help='Frozen comma-separated nonnegative gamma grid')
    parser.add_argument('--slopes',default='-1,0,1')
    parser.add_argument('--reference',type=float,required=True)
    parser.add_argument('--eta-ref',type=float,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(hashlib.sha256(args.registration.read_bytes()).hexdigest()==args.registration_sha256,
            'Registration hash mismatch')
    require(not args.output.exists(), 'Refusing to overwrite prior experiment output')
    amplitudes=[float(s) for s in args.amplitudes.split(',')]
    slopes=[float(s) for s in args.slopes.split(',')]
    require(amplitudes==sorted(set(amplitudes)) and amplitudes[0]>=0,'Need increasing unique amplitudes')
    settings=[dict(rtol=1e-11,y0=1e-3,max_step=.1),dict(rtol=3e-14,y0=1e-4,max_step=.05)]
    out=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
             registration_sha256=args.registration_sha256,
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             versions=dict(numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__),
             model=dict(reference=args.reference,eta_ref=args.eta_ref,amplitudes=amplitudes,slopes=slopes),
             independence='Second-order acceleration bulk, separately coded resolvent/radius-variation source; no producer imports.',
             runs=[],failures=[])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    def persist():
        args.output.write_text(json.dumps(out,indent=2)+'\n')
    persist()
    for setting in settings:
        for b in slopes:
            start=[math.log10(-INHERITED_ETA_H),INHERITED_YB]
            for gamma in amplitudes:
                try:
                    run=solve_case(start,gamma,b,args.reference,args.eta_ref,setting)
                    out['runs'].append(run)
                    print(json.dumps({k:run[k] for k in ('gamma','b','ell','y_b','H2','phi','E1','E2','Hamiltonian_relative_max','root_success')}),flush=True)
                    if run['residual_gate'] and run['Hamiltonian_gate'] and all(c['detected'] for c in run['negative_control_checks'] if c['applicable']):
                        start=[run['ell'],run['y_b']]
                    else:
                        out['failures'].append(dict(gamma=gamma,b=b,setting=setting,error='Residual, Hamiltonian, or negative-control gate failed'))
                except Exception as error:
                    out['failures'].append(dict(gamma=gamma,b=b,setting=setting,error=repr(error)))
                    print('FAILED',gamma,b,repr(error),flush=True)
                persist()
    gates=[]
    length=len(amplitudes)*len(slopes)
    if len(out['runs'])==2*length:
        for first,second in zip(out['runs'][:length],out['runs'][length:]):
            phi_error=abs(first['phi']-second['phi'])
            H2_relative=abs(first['H2']/second['H2']-1)
            coords=max(abs(first['ell']-second['ell']),abs(first['y_b']-second['y_b']))
            gates.append(dict(gamma=first['gamma'],b=first['b'],phi_absolute_difference=phi_error,
                              H2_relative_difference=H2_relative,coordinate_max_absolute_difference=coords,
                              passed=phi_error<=2e-9 and H2_relative<=2e-8 and coords<=2e-7))
    out['refinement_gates']=gates
    out['status']='PASS' if len(gates)==length and not out['failures'] and all(g['passed'] for g in gates) else 'FAIL'
    out['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    persist()
    require(out['status']=='PASS','Independent stationary root gates failed')


if __name__=='__main__':
    main()
