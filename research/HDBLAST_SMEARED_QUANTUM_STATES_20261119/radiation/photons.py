"""Physical photon thermal spectrum and a conserved-radiation test for a closed bounce."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]


def dimensionless_density(y, tail=42):
    """rho_excess * a_physical^4 in natural units, both transverse polarizations."""
    if y < 0 or not math.isfinite(y):
        raise ValueError('Temperature times radius must be finite and nonnegative')
    if y == 0:
        return 0.
    labels = np.arange(2, max(3, math.ceil(tail*y)+2), dtype=float)
    boltzmann = np.exp(-labels/y)
    return float(np.sum(labels*(labels**2-1)*boltzmann/(-np.expm1(-labels/y)))/math.pi**2)


def y_for_d(d, calibration):
    if d == 0:
        return 0.
    target = 3*d/calibration**2
    upper = max(1.,2*(15*target/math.pi**2)**.25)
    while dimensionless_density(upper) < target:
        upper *= 2
    return brentq(lambda y:dimensionless_density(y)-target, 1e-5, upper,xtol=1e-13)


def evolve(branch, frac, *, direct=False):
    """frac is d/dcrit, where r_gamma=d/a^4. a=H0 a_phys, r=rho/V0."""
    initial_q = branch.fitted.sign*branch.fitted.dq
    ui = float(branch.potential(initial_q)[0])
    d = frac/(4*ui)
    a0 = math.sqrt((1+math.sqrt(1-frac))/(2*ui))
    def unpack(y):
        q,v,third,fourth,fifth = y
        if direct:
            return q,v,third,fourth/third,fifth
        a = a0*math.exp(third)
        return q,v,a,fourth,d/a**4
    def rhs(t,y):
        q,v,a,h,r = unpack(y)
        if not branch.fitted.theta[0] <= branch.top+q <= branch.fitted.theta[-1]:
            raise ArithmeticError("ODE stage outside the calibrated potential")
        u,force = branch.potential(q)
        acceleration = u-v*v/3-r
        return [v,force-3*h*v,y[3] if direct else h,
                a*acceleration if direct else 1/a**2-v*v/2-2*r, -4*h*r if direct else 0.]
    def turnaround(t,y): return y[3]
    turnaround.terminal=True
    turnaround.direction=-1
    def passage(t,y): return branch.fitted.sign*(branch.top+y[0]-branch.fitted.tmin)
    passage.direction=1
    def guard(t,y):
        theta=branch.top+y[0]
        return min(theta-branch.fitted.theta[0],branch.fitted.theta[-1]-theta)-.045
    guard.terminal=True
    guard.direction=-1
    ans=solve_ivp(rhs,(0.,3.),[initial_q,0.,a0 if direct else 0.,0.,d/a0**4],
                  method='RK45' if direct else 'DOP853',rtol=2e-12,atol=2e-14,
                  max_step=.003 if direct else .012,dense_output=True,events=[turnaround,passage,guard])
    if not ans.success or len(ans.t_events[0])!=1 or ans.t[-1]<1e-7:
        raise ArithmeticError('Expanding bounce and regular first maximum required')
    def record(t):
        q,v,a,h,r=unpack(ans.sol(t))
        u=float(branch.potential(q)[0])
        return dict(x=float(t),a=float(a),q=float(q),h=float(h),r_photon=float(r),
                    matter_fraction=float(r/(r+u+v*v/6)),
                    constraint=float(h*h+1/a**2-r-u-v*v/6),
                    dilution_residual=float(r*a**4-d))
    samples=[record(x) for x in np.linspace(0,ans.t[-1],201)]
    start=samples[0]
    residual=max(abs(row['constraint']) for row in samples)
    identity=max(abs(row['matter_fraction']-d/(row['a']**2*(1+row['a']**2*row['h']**2))) for row in samples)
    above=max(row['matter_fraction']-d/row['a']**2 for row in samples)
    if max(residual,identity,above,max(abs(row['dilution_residual']) for row in samples))>2e-9:
        raise ArithmeticError('Radiation balance or closed-universe identity failed')
    return dict(initial=start,first_crossing=record(ans.t_events[1][0]) if len(ans.t_events[1]) else None,
                first_maximum=samples[-1], sampled_maximum_fraction=max(row['matter_fraction'] for row in samples),
                curve=[dict(x=row['x'],a_over_ai=row['a']/a0,photon_fraction=row['matter_fraction'],
                            ceiling=d/row['a']**2,
                            comoving_injection_needed_for_half_over_ai2=(
                                 .5*(row['a']/a0)**2*(1+row['a']**2*row['h']**2)-d/a0**2))
                       for row in samples[::5]],
                max_constraint=residual,max_fraction_identity_error=identity,
                domain_lower_bound_sampled=min(min(branch.top+row['q']-branch.fitted.theta[0],branch.fitted.theta[-1]-branch.top-row['q']) for row in samples))


def calculate(prior):
    spec=importlib.util.spec_from_file_location('prior_round_solver',prior/'homogeneous/exchange.py')
    earlier=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(earlier)
    critical_spec=importlib.util.spec_from_file_location('critical_bounce',Path(__file__).with_name('critical.py'))
    critical=importlib.util.module_from_spec(critical_spec)
    critical_spec.loader.exec_module(critical)
    branches=earlier.branches(prior/'prior/gravity')
    outputs={}
    for key,branch in branches.items():
        initial_q=branch.fitted.sign*branch.fitted.dq
        ui=float(branch.potential(initial_q)[0])
        hpl=float(branch.fitted.Hpl)
        critical_d=1/(4*ui)
        y_critical=y_for_d(critical_d,hpl)
        if abs(dimensionless_density(y_critical)*hpl**2/3/critical_d-1)>2e-13:
            raise ArithmeticError('Finite-volume thermal threshold mismatch')
        a_critical=1/math.sqrt(2*ui)
        thermal=[]
        for y in [0.,1/(2*math.pi),1.,3.,y_critical]:
            d=dimensionless_density(y)*hpl**2/3
            z=1-4*ui*d
            a=math.sqrt((1+math.sqrt(max(z,0)))/(2*ui))
            thermal.append(dict(aT=y,energy_coefficient=dimensionless_density(y),d=d,
                              fraction_critical=d/critical_d,launch_fraction=d/a**2,
                              T_over_H0=y/a,thermal_sum_tail_test=abs(dimensionless_density(y,60)-dimensionless_density(y))))
        runs={}
        for frac in [0.,.01,.25,.75,.95,.99,.9999]:
            runs[f'{frac:g}']=evolve(branch,frac)
        differences=[]
        for frac in [.25,.99]:
            independent=evolve(branch,frac,direct=True)
            for stage in ('first_maximum','first_crossing'):
                saved=runs[f'{frac:g}'][stage]
                check=independent[stage]
                if (saved is None)!=(check is None):raise ArithmeticError('Event sets differ')
                if saved:
                    differences.extend(abs(saved[k]-check[k]) for k in ('x','a','h','q','r_photon','matter_fraction'))
        if max(differences)>1e-8:raise ArithmeticError('Independent solution differs')
        # An interior change in event ordering occurs when first maximum and first passage coincide.
        order=[bool(v['first_crossing']) for v in runs.values()]
        if order[-1] or not order[0]: raise ArithmeticError('Unexpected crossing sequence')
        # Use only the bracket found by adjacent sampled preparations and certify the local change.
        cases=[float(k) for k in runs]
        bracket=next((l,h) for l,h in zip(cases,cases[1:]) if runs[f'{l:g}']['first_crossing'] is not None and runs[f'{h:g}']['first_crossing'] is None)
        # Earlier passages may be followed by returns in loosely radiation-loaded runs,
        # so solve for the speed-qualified coalescence of these two first events locally.
        lo,hi=bracket
        for unused in range(32):
            middle=(lo+hi)/2
            if evolve(branch,middle)['first_crossing'] is None:hi=middle
            else:lo=middle
        boundary=evolve(branch,(lo+hi)/2)
        distance=abs(branch.top+boundary['first_maximum']['q']-branch.fitted.tmin)
        if distance>2e-7:raise ArithmeticError('Ordering boundary misses the minimum')
        outputs[key]=dict(ui=ui,H0_over_Mpl=hpl,critical_check=critical.critical_check(branch,evolve(branch,1.-1.e-8)['first_maximum']['x']),critical_d=critical_d,critical_aT=y_critical,
                          critical_T_over_H0=y_critical/a_critical,thermal_states=thermal,runs=runs,
                          transition_bracket=[lo,hi],transition_field_error=distance,
                          initial_fraction_at_order_transition=boundary['initial']['matter_fraction'],
                          max_independent_difference=max(differences))
    # Analytic high-temperature leading and first curvature corrections, and low-temperature L=2.
    y=40.
    high=dimensionless_density(y)
    expansion=(math.pi**2/15)*y**4-y**2/6+y/(2*math.pi**2)-11/(240*math.pi**2)
    ylow=.04
    low=dimensionless_density(ylow)
    leading=6/math.pi**2*math.exp(-2/ylow)
    if abs(high/expansion-1)>3e-13 or abs(low/leading-1)>1e-9:raise ArithmeticError('Independent spectral limits failed')
    return dict(branches=outputs,spectrum_checks=dict(high_aT=y,relative_high_T_error=high/expansion-1,
                                  low_aT=ylow,relative_low_T_error=low/leading-1))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior',type=Path,default=ROOT/'prior')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior),indent=2,allow_nan=False)+'\n')
