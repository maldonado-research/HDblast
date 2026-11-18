"""Independent momentum transport and proper-time energy ledgers on rolling B1 walls."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import brentq

CUTS = [1.45, 1.486, 1.50, 1.8, 2.15]
COUPLING = 100.
PHI_STAR = .5


def correction_terms(event):
    """Reconstruct the five terms in the source's first-order coefficient."""
    v, acc, jerk, h, dh = (event[k] for k in ('v_fit', 'A', 'B', 'h0_fit', 'h1_fit'))
    pi = math.pi
    terms = {
        'H_squared': (9*pi/4+15/(2*pi))*h*h,
        'H_prime': (3*pi/2-15/(8*pi))*dh,
        'acceleration_squared': (45/(32*pi)-pi/8)*(acc/v)**2,
        'jerk': (pi/8-45/(96*pi))*jerk/v,
        'H_acceleration': (45/(8*pi))*h*acc/v,
    }
    return {k: x/event['q'] for k, x in terms.items()}


def history(path):
    with np.load(path/'rolling_timeseries.npz', allow_pickle=False) as z:
        d = dict(zip(map(str, z['cols']), z['rec'].T, strict=True))
    summary = json.loads((path/'rolling_summary.json').read_text())
    event, = summary['matter']['events']
    if not np.all(np.diff(d['H0tau'])>0) or np.any(np.diff(d['T'])<=0):
        raise ValueError('Nonmonotone history')
    if not (d['H0tau'][0] < CUTS[0] and CUTS[-1] < d['H0tau'][-1]):
        raise ValueError('Proper cuts outside recorded history')
    return d, summary, event


def transport(d, event, nodes):
    """One passing cohort, Gaussian integral over rescaled physical momentum."""
    s = d['H0tau']
    phi = CubicHermiteSpline(s, d['phi_b'], d['v_over_H0'])
    A = CubicHermiteSpline(s, d['A_b'], d['H_over_H0'])
    cross = brentq(lambda t: phi(t)-PHI_STAR, s[0], s[-1])
    vstar = float(phi(cross,1))
    q = COUPLING*abs(vstar)
    Astar = float(A(cross))
    width = 3*abs(vstar)/math.sqrt(q)
    finish = brentq(lambda t: phi(t)-PHI_STAR-width, cross, s[-1])
    slope_min = float(min(phi(np.linspace(cross,finish,2001),1)))
    if slope_min <= .05*vstar:
        raise ValueError('Nonmonotone production progress or active low-speed switch')
    x, w = np.polynomial.legendre.leggauss(nodes)
    # z=sqrt(pi/q) k in [0,8], Gaussian tail beyond cut below 2e-27 by number.
    z, weight = (x+1)*4, w*4
    k = math.sqrt(q/math.pi)*z
    target = weight*(q/math.pi)**1.5*z*z*np.exp(-z*z)/(2*math.pi**2)
    def parameters(t):
        disp = phi(t)-PHI_STAR
        a = np.exp(A(t)-Astar)
        mass = COUPLING*abs(disp)
        omega = np.hypot(k/a,mass)
        return disp,a,mass,omega
    def rhs(t,state):
        disp,a,m,omega = parameters(t)
        # The published ramp stops after the scalar traverses width.
        u = float(np.clip(disp/width,0.,1.))
        inject = 30*u*u*(1-u)**2*phi(t,1)/width
        decay = m*m/(8*math.pi)
        N = state[:-1]
        dN = target*inject - decay*N/omega
        d_rad_comoving = a*decay*N.sum()
        return np.append(dN,d_rad_comoving)
    args = dict(method='DOP853', rtol=2e-10, atol=2e-13, max_step=.0035, dense_output=True)
    before = solve_ivp(rhs,(cross,finish),np.zeros(nodes+1),**args)
    after = solve_ivp(rhs,(finish,CUTS[-1]),before.y[:,-1],**args)
    if not before.success or not after.success:
        raise ValueError('Independent transport failed')
    factor = event['number_factor']
    outputs = {}
    for t in CUTS:
        R = X = 0.
        if t > cross:
            y = (before if t <= finish else after).sol(t)
            _,a,_,omega = parameters(t)
            R, X = float(y[-1]/a**4), float(np.dot(y[:-1],omega)/a**3)
        outputs[str(t)] = dict(R_with_F=factor*R, X_with_F=factor*X,
                               R_without_F=R, X_without_F=X)
    return dict(s_cross=cross, v_cross=vstar, q=q, s_ramp_end=finish, min_v_while_ramping=slope_min,
                gaussian_normalization_ratio=float(target.sum()/(q**1.5/(8*math.pi**3))),
                cuts=outputs)


def ledger(d, summary):
    """Read only sampled physical densities and sources; integrate in s=H0 tau."""
    s = d['H0tau']
    Rhat = d['R']*summary['rho_b']/summary['matter']['b']
    Xhat = d['rho_hat']
    a = np.exp(d['A_b']-d['A_b'][0])
    h = d['H_over_H0']
    # Rebuild decay power from mass and gas number. Q_hat is a separate stored diagnostic.
    Q = (COUPLING*(d['phi_b']-PHI_STAR))**2*d['n_hat']/(8*math.pi)
    work = (d['j_hat']+d['jprod_hat'])*d['v_over_H0']
    pressure_work = 3*h*d['p_hat']
    cols = dict(R=Rhat, X=Xhat, a=a, rad_source4=a**4*Q, decay3=a**3*Q,
                scalar_work3=a**3*work, chi_pressure3=a**3*pressure_work,
                radiation_pressure3=a**3*h*Rhat)
    interpolants = {key:PchipInterpolator(s, vals) for key, vals in cols.items()}
    def integrate(row,key,t):
        times = np.r_[s[s < t],t]
        return float(simpson(np.r_[cols[key][:len(times)-1],row[key]],x=times))
    out = {}
    for t in CUTS:
        row = {key:float(interp(t)) for key,interp in interpolants.items()}
        # Smooth interpolation across a newly opened cohort can leak tiny densities
        # into the pre-crossing interval. The registered empty state is exact there.
        if t <= summary['matter']['events'][0]['s_star']:
            for key in row.keys()-{'a'}:
                row[key] = 0.
        rad = integrate(row,'rad_source4',t)
        decay = integrate(row,'decay3',t)
        work = integrate(row,'scalar_work3',t)
        p_chi = integrate(row,'chi_pressure3',t)
        p_rad = integrate(row,'radiation_pressure3',t)
        R, X, a_t = row['R'], row['X'], row['a']
        out[str(t)] = dict(R=R,X=X,R_from_integrated_decay=rad/a_t**4,
                           U_R=a_t**3*R,U_chi=a_t**3*X,total_scalar_work=work,
                           integrated_decay=decay,chi_pressure_work=p_chi,
                           radiation_pressure_work=p_rad,
                           chi_residual=a_t**3*X-work+p_chi+decay,
                           total_residual=a_t**3*(R+X)-work+p_chi+p_rad)
    return out


def percent(x,base):
    return 100*(x/base-1) if base else None


def analyze(base):
    out = {}
    for grid in ('coarse','fine'):
        branches = {}
        for case in ('quiet','0.001','-0.001'):
            d, summary, event = history(base/grid/case)
            ref = transport(d,event,48)
            check = transport(d,event,32)
            terms = correction_terms(event)
            if abs(sum(terms.values())-event['D_over_q'])>1e-12:
                raise ValueError('First-order coefficient reconstruction failed')
            num = ref['cuts']['2.15']
            lower = check['cuts']['2.15']
            terms_error = max(abs(num[k]/lower[k]-1) for k in num)
            energy = ledger(d,summary)
            last = energy['2.15']
            errors = dict(legendre_32_vs_48_relative_max=terms_error,
                          R_relative_transport_minus_recorded=num['R_with_F']/last['R']-1,
                          X_relative_transport_minus_recorded=num['X_with_F']/last['X']-1,
                          R_relative_ledger_minus_recorded=last['R_from_integrated_decay']/last['R']-1,
                          total_relative_ledger_defect=last['total_residual']/last['total_scalar_work'])
            if max(map(abs, errors.values()))>.002:
                raise ValueError(f'Budget check failed on {grid}/{case}: {errors}')
            branches[case] = dict(event={k:event[k] for k in ('q','s_star','D_over_q',
                                'number_factor','n_hat_at_crossing','screen_pass','G_min_screen')},
                                D_over_q_terms=terms,transport=ref,ledger=energy,checks=errors)
        contrasts = {}
        quiet = branches['quiet']
        for case in ('0.001','-0.001'):
            item = branches[case]
            e, eq = item['event'],quiet['event']
            result = dict(q_pct=percent(e['q'],eq['q']),
                          q_to_3_over_2_pct=percent(e['q']**1.5,eq['q']**1.5),
                          F_pct=percent(e['number_factor'],eq['number_factor']),
                          target_yield_pct=percent(e['n_hat_at_crossing'],eq['n_hat_at_crossing']),
                          delta_D_over_q_components={k:x-quiet['D_over_q_terms'][k]
                                                     for k,x in item['D_over_q_terms'].items()})
            cuts = {}
            for t in map(str,CUTS):
                c,b = item['ledger'][t],quiet['ledger'][t]
                r,u = item['transport']['cuts'][t],quiet['transport']['cuts'][t]
                cuts[t] = dict(dR=c['R']-b['R'], dX=c['X']-b['X'],
                               recorded_R_pct=percent(c['R'],b['R']), recorded_X_pct=percent(c['X'],b['X']),
                               transport_R_pct=percent(r['R_with_F'],u['R_with_F']),
                               exact_fixed_wall_R_no_D_pct=percent(c['R']/e['number_factor'],b['R']/eq['number_factor']),
                               transport_R_no_D_pct=percent(r['R_without_F'],u['R_without_F']),
                               transport_X_no_D_pct=percent(r['X_without_F'],u['X_without_F']),
                               R_per_target_pct=percent(c['R']/e['n_hat_at_crossing'],b['R']/eq['n_hat_at_crossing']))
            result['cuts'] = cuts
            contrasts[case] = result
        out[grid] = dict(branches=branches,contrasts=contrasts)
    return out


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',type=Path)
    print(json.dumps(analyze(parser.parse_args().runs),indent=2,allow_nan=False))
