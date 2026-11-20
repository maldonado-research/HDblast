"""Radiation production required in a dissipative, closed Einstein–scalar benchmark."""
import argparse
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
RATES = (0,2,4,6,8,10,12,15,20,24,30,40,50,64,80,96,128,160,192,224,
         256,272,288,304,320,352)


def window(s):
    """Compact smooth dissipative weight in directed displacement s."""
    z = (s-.2)/.6
    return math.exp(4-1/(z*(1-z))) if 0<z<1 else 0.


def calculate_trajectory(branch, gamma, direct=False):
    qi = branch.fitted.sign*branch.fitted.dq
    displacement = branch.fitted.tmin-branch.top-qi
    ui = float(branch.potential(qi)[0])
    ai = 1/math.sqrt(ui)
    def unpack(y):
        q,v,third,fourth,fifth,*_ = y
        if direct:
            a,h = third,fourth/third
            return q,v,a,h,fifth*(ai/a)**4
        return q,v,ai*np.exp(third),fourth,fifth

    def rhs(t, y):
        q,v,a,h,r = unpack(y)
        u,force = branch.potential(q)
        rate = gamma*window((q-qi)/displacement)
        source = rate*v*v/3
        injection = (a/ai)**4*source
        derivatives = ((a*h, a*(u-v*v/3-r), injection) if direct
                       else (h, 1/a**2-v*v/2-2*r, source-4*h*r))
        return (v,force-(3*h+rate)*v,*derivatives,source,4*h*r,h*v*v,injection)

    def turnaround(t,y):
        return y[3]
    turnaround.direction = -1
    turnaround.terminal = True
    def passage(t,y):
        return (y[0]-qi)/displacement-1
    passage.direction = 1
    def outside(t,y):
        theta = branch.top+y[0]
        return min(theta-branch.fitted.theta[0],branch.fitted.theta[-1]-theta)-.045
    outside.terminal = True
    outside.direction = -1
    initial = (qi,0,ai if direct else 0,0,0,0,0,0,0)
    sol = solve_ivp(rhs,(0.,24.),initial,
                    method='RK45' if direct else 'DOP853',dense_output=True,
                    rtol=2e-12 if direct else 2e-10, atol=2e-14 if direct else 3e-12,
                    max_step=.008 if direct else .035, events=(turnaround,passage,outside))
    if not sol.success or len(sol.t_events[0])!=1 or sol.t[-1]<1e-6:
        raise ArithmeticError(f'No regular first turnaround in the fit for gamma={gamma}')
    def variables(t):
        y=sol.sol(t)
        q,v,a,h,r=unpack(y)
        u,_=branch.potential(q)
        parent=u+v*v/6
        total=parent+r
        return y,q,v,a,h,r,parent,total
    def fraction(t):
        *_,r,_,total=variables(t)
        return float(r/total)
    def fraction_rate(t):
        _,q,v,_,h,r,_,total=variables(t)
        source=gamma*window((q-qi)/displacement)*v*v/3
        return float((source-4*h*r+(r/total)*h*(v*v+4*r))/total)

    end=float(sol.t[-1])
    first=float(sol.t_events[1][0]) if len(sol.t_events[1]) else None
    # Locate stationary points and use their continuous dense solutions, with
    # both endpoints included. Each accepted step contributes four search intervals.
    mesh=np.unique(np.r_[(sol.t[:-1,None]+np.diff(sol.t)[:,None]*np.arange(4)/4).ravel(),end])
    y,q,v,a,h,r,parent,total=variables(mesh)
    shape=np.array([window(s) for s in (q-qi)/displacement])
    rates=(gamma*shape*v*v/3-4*h*r+r/total*h*(v*v+4*r))/total
    stationary=[]
    for i in np.flatnonzero(rates[:-1]*rates[1:]<0):
        stationary.append(float(brentq(fraction_rate,mesh[i],mesh[i+1],xtol=3e-13)))
    candidates=sorted({0.,end,*stationary,*([first] if first is not None else [])})
    peak=max(candidates,key=fraction)
    first_peak=max((x for x in candidates if first is None or x<=first),key=fraction)
    half=None
    for lo,hi in zip(candidates,candidates[1:]):
        if fraction(lo)<.5<=fraction(hi):
            half=float(brentq(lambda t:fraction(t)-.5,lo,hi,xtol=3e-13))
            break
    def record(t):
        y,q,v,a,h,r,parent,total=variables(t)
        s=(q-qi)/displacement
        rate=gamma*window(s)
        fraction=float(r/total)
        w_phi=float(v*v/(3*parent)-1)
        dF=float(rate*v*v/(3*total)+h*fraction*(1-fraction)*(3*w_phi-1))
        return dict(x=t,q=float(q),v=float(v),a=float(a),h=float(h),r=float(r),
                    directed_displacement_fraction=float(s),instantaneous_Q=float(rate*v*v/3),
                    parent_equation_of_state=w_phi,fraction_rate=dF,
                    scalar_density=float(parent),fraction=fraction,
                    total_comoving_input_over_ai4=float(y[8]),
                    required_comoving_input_for_half_over_ai4=float(.5*(a/ai)**4*(h*h+1/a**2)),
                    integrated_Q=float(y[5]),direct_to_friction_ratio=(float(
                        rate/(3*h)) if h>1e-7 else None))
    source,dilution,pressure,comoving=y[5:]
    balances=dict(closed_constraint=float(max(abs(h*h+1/a**2-total))),
                  parent=float(max(abs(parent-ui+pressure+source))),
                  radiation=float(max(abs(r-source+dilution))),
                  comoving=float(max(abs((a/ai)**4*r-comoving))),
                  radiation_lower_violation=float(max(0.,-min(r))))
    if max(balances.values())>2e-7:
        raise ArithmeticError(f'Constraint or energy budget failed at gamma={gamma}: {balances}')
    return dict(gamma=float(gamma),first_turnaround=record(end),
                first_minimum=record(first) if first is not None else None,
                maximum_before_turnaround=record(peak),maximum_through_first_minimum=record(first_peak),
                first_half=record(half) if half is not None else None,
                balances=balances)


def calculate(prior):
    spec=importlib.util.spec_from_file_location('archived_round_solver',prior/'homogeneous/exchange.py')
    archived=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(archived)
    old=json.loads((prior/'homogeneous/results.json').read_bytes())['branches']
    branches=archived.branches(prior/'prior/gravity')
    output=dict(rate_profile='Gamma/H0=gamma*b((s-.2)/.6), s=(Theta-Theta_i)/(Theta_min-Theta_i)',
                peak_dimensionless_gamma_scan=list(RATES),initial_radiation=0.,branches={})
    for key,branch in branches.items():
        @lru_cache(None)
        def run(gamma):
            return calculate_trajectory(branch,gamma)
        scan={str(g):run(g) for g in RATES}
        zero=scan['0']
        historical=old[key]['runs']['inverse_f_seed_0']
        mismatch=max(abs(zero[here][field]-historical[there][field])
                     for here,there in [('first_turnaround','at_stop'),('first_minimum','first_minimum_crossing')]
                     for field in ('x','a','q','v','h'))
        if mismatch>2e-8 or any(zero[k]['r']!=0. for k in ('first_turnaround','first_minimum')):
            raise ArithmeticError('Unsourced historical control failed')
        crossings={}
        for label in ('maximum_before_turnaround','maximum_through_first_minimum'):
            brackets=[(left,right) for left,right in zip(RATES,RATES[1:])
                      if (scan[str(left)][label]['fraction']-.5)*(scan[str(right)][label]['fraction']-.5)<0]
            if len(brackets)!=2:
                raise ArithmeticError('The scan has an unexpected number of threshold brackets')
            rows=[]
            for lo,hi in brackets:
                root=float(brentq(lambda g:run(g)[label]['fraction']-.5,lo,hi,xtol=2e-7))
                peak=run(root)[label]
                if label=='maximum_before_turnaround' and (peak['instantaneous_Q']!=0 or
                        abs(peak['parent_equation_of_state']-1/3)>1e-7 or abs(peak['fraction_rate'])>1e-8):
                    raise ArithmeticError('Full threshold peak missed its source-free stationarity condition')
                rows.append(dict(bracket=[lo,hi],gamma=root,peak=peak))
            crossings[label]=rows
        # Resolve the local efficiency maximum around the largest scanned value.
        label='maximum_before_turnaround'
        best=max(range(len(RATES)),key=lambda i:scan[str(RATES[i])][label]['fraction'])
        opt=minimize_scalar(lambda g:-run(float(g))[label]['fraction'],
                            bounds=(RATES[best-1],RATES[best+1]),method='bounded',
                            options={'xatol':5e-6})
        if not opt.success:
            raise ArithmeticError(opt.message)
        optimum=run(float(opt.x))
        independent=[]
        for g in [0.,*[row['gamma'] for row in crossings[label]],float(opt.x),352.]:
            primary=run(g)
            direct=calculate_trajectory(branch,g,direct=True)
            errors={stage:{k:abs(primary[stage][k]-direct[stage][k]) for k in
                          ('x','a','fraction','r','integrated_Q','total_comoving_input_over_ai4')}
                    for stage in ('first_turnaround','first_minimum','maximum_before_turnaround')}
            if max(val for row in errors.values() for val in row.values())>1e-5:
                raise ArithmeticError(f'Independent trajectories disagree: {errors}')
            independent.append(dict(gamma=g,absolute_discrepancies=errors))
        output['branches'][key]=dict(scan=scan,half_thresholds=crossings,
                    local_efficiency_optimum=optimum,independent_formulation=independent,
                    unsourced_historical_maximum_difference=mismatch)
        print(f'{key} positive-transfer sweep and independent checks passed',flush=True)
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=HERE.parent/'prior')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior),indent=2,allow_nan=False)+'\n')
