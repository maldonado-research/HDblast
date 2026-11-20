"""Conservative release into the same spherical Maxwell triplet at a timelike wall."""
import importlib.util
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar


def load_solution():
    path = Path(__file__).with_name('loading.py')
    spec = importlib.util.spec_from_file_location('release_wall_model',path)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    return model.solution


def brief(q):
    keys=('t','a','H','psi','pi','V','mu','K','X','Bmag','sheet','e','rad')
    row={key:float(q[key]) for key in keys}
    row['sheet_fraction']=q['sheet']/q['e'] if q['e']>0 else 0.
    return row


def samples(sol,lo,hi,order=8):
    """Gauss quadrature on every original accepted solver interval in the window."""
    boundaries=np.unique(np.r_[lo,np.concatenate([p.t for p in sol.parts]),hi])
    knots=boundaries[(boundaries>=lo)&(boundaries<=hi)]
    nodes,weights=np.polynomial.legendre.leggauss(order)
    half=np.diff(knots)[:,None]/2
    times=knots[:-1,None]+half*(nodes+1)
    return (half*weights).ravel(),sol.trace_many(times.ravel())


def account(sol,lo,hi):
    weights,z=samples(sol,lo,hi)
    start,end=sol.trace(lo),sol.trace(hi)
    a,h,pi,k,dk,x,b,i,o,e,es,pressure=map(z.get,['a','H','pi','K','Kprime','X','Bmag','I','O','e','sheet','p'])
    kd=dk*pi
    transfer=.5*kd*(b**2-x**2)
    def integral(rate):
        return float(weights@rate)
    body_heat=integral(a**3*transfer)
    emitted=integral(2*a**3*o**2)
    intercepted=integral(2*a**3*i**2)
    gauge_work=integral(a**3*h*es)
    scalar_work=integral(3*a**3*h*(pressure-es/3))
    def delta(observable):
        return float(observable(end)-observable(start))
    scalar_delta=delta(lambda q:q['a']**3*(q['e']-q['sheet']))
    sheet_delta=delta(lambda q:q['a']**3*q['sheet'])
    normalized_delta=delta(lambda q:q['a']**4*q['sheet']/q['K'])
    normalized_rhs=integral(a**4/k*(2*(i*i-o*o)-kd*x*x))
    checks=dict(scalar_budget=abs(scalar_delta+scalar_work+body_heat),
                sheet_budget=abs(sheet_delta+gauge_work-body_heat+emitted-intercepted),
                normalized_gauge_budget=abs(normalized_delta-normalized_rhs))
    if max(checks.values())>5e-7:
        raise ArithmeticError(checks)
    return dict(from_time=lo,to_time=hi,scalar_to_sheet=body_heat,bulk_input=intercepted,
                outgoing_to_both_copies=emitted,sheet_pressure_work=gauge_work,
                scalar_pressure_work=scalar_work,scalar_internal_change=scalar_delta,
                sheet_internal_change=sheet_delta,scaled_gauge_change=normalized_delta,
                scaled_gauge_flux_integral=normalized_rhs,minimum_Kdot=float(min(kd)),errors=checks)


def limit_at_crossing(sol,exit_state,cross):
    weights,z=samples(sol,exit_state['t'],cross['t'])
    monotone=float(min(z['pi']*z['Kprime']))
    if monotone < -1e-10 or max(abs(z['I']))>1e-12:
        return None
    required=cross['a']**4*.9*cross['e']/(exit_state['a']**4*exit_state['sheet'])
    chosen=cross['K']/exit_state['K']
    if required>chosen*(1+1e-8):
        raise ArithmeticError('Postpulse gauge seed bound')
    return dict(measured_seed_density=exit_state['sheet'],minimum_seed_density=exit_state['sheet']*required/chosen,
                required_K_ratio=required,realized_K_ratio=chosen,minimum_Kdot=monotone)


def summarize(sol,late,detailed=False):
    ex=sol.events['exit'][0]
    maximum=next((q for q in sol.events['maximum'] if q['t']>ex['t']),sol.trace(sol.end))
    downs=sol.events['down_end']
    crossing=sorted([q for q in sol.events['sheet90'] if ex['t']<q['t']<maximum['t']],key=lambda r:r['t'])
    cuts=[ex]+crossing+[maximum]
    epochs=[]
    for a,b in zip(cuts,cuts[1:]):
        mid=sol.trace((a['t']+b['t'])/2)
        if mid['H']>0 and mid['sheet']>.9*mid['e']:
            epochs.append(dict(start=brief(a),end=brief(b),efolds=math.log(b['a']/a['a'])))
    # Search all locally bracketing cells, with both expansion endpoints eligible.
    times=np.linspace(ex['t'],maximum['t'],401)
    dat=sol.trace_many(times)
    score=dat['sheet']/dat['e']
    peaks=[times[0],times[-1]]
    def fraction(t):
        q=sol.trace(t)
        return q['sheet']/q['e']
    for idx in np.where((score[1:-1]>=score[:-2])&(score[1:-1]>score[2:]))[0]+1:
        fit=minimize_scalar(lambda t:-fraction(t),bounds=(times[idx-1],times[idx+1]),method='bounded')
        peaks.append(float(fit.x))
    peak=max(peaks,key=fraction)
    native=np.unique(np.concatenate([part.t for part in sol.parts]))
    allpoints=sol.trace_many(native)
    available=allpoints['mu']
    fried=allpoints['H']**2-allpoints['e']*(2+allpoints['e'])+1/allpoints['a']**2-available/allpoints['a']**4
    if min(available)<0 or max(abs(allpoints['rad']))>1e-12 or max(abs(fried))>2e-7:
        raise ArithmeticError(('Gravity or independent fluid check', min(available),max(abs(fried))))
    last=sol.trace(sol.end)
    row=dict(dose=sol.config['dose'],late_K_increment=late,bounce=brief(sol.events['bounce'][0]),
             exit=brief(ex),descent_end=brief(downs[0]) if downs else None,
             expansion_maximum=brief(maximum),peak=brief(sol.trace(peak)),expanding_epochs=epochs,
             terminal=brief(last),minimum_mass=float(min(available)),
             post_exit_min_velocity_during_expansion=float(min(dat['pi'])),
             max_friedmann_error=float(max(abs(fried))),max_fluid_density=float(max(abs(allpoints['rad']))))
    turning=np.flatnonzero((native[:-1]>=ex['t'])&(allpoints['pi'][:-1]>0)&(allpoints['pi'][1:]<=0))
    row['first_scalar_reversal']=None
    if len(turning):
        j=turning[0]
        time=brentq(lambda t:sol.trace(t)['pi'],native[j],native[j+1])
        turn=sol.trace(time)
        row['first_scalar_reversal']=brief(turn)
        if .36<turn['psi']<.38:
            row['first_scalar_reversal']['net_force_per_step_slope']=.111-.5*late*(turn['Bmag']**2-turn['X']**2)
    if detailed:
        targets={'expansion_maximum':maximum,'terminal':last}
        if downs:
            targets['descent_end']=downs[0]
        if epochs:
            targets['radiation_entry']=epochs[0]['start']
            row['crossing_seed_bound']=limit_at_crossing(sol,ex,epochs[0]['start'])
        row['exchange']={name:account(sol,ex['t'],q['t']) for name,q in targets.items()}
    return row


def calculate(solve=None):
    solve=solve or load_solution()
    choices=dict(no_late_coupling=(.123,0.),small=(.123,1500.),below_window=(.123,1850.),
                 lower_success=(.123,1900.),primary=(.123,2000.),larger=(.123,2150.),
                 above_window=(.123,2400.),stalled=(.123,8000.),
                 weaker_pulse=(.12275,2000.),stronger_pulse=(.12325,2000.))
    rows={}
    for name,(dose,late) in choices.items():
        sol=solve(dose=dose,damping=0.,late_gauge=late,tol=2e-12)
        rows[name]=summarize(sol,late,detailed=name in ('primary','no_late_coupling'))
    ref=summarize(solve(damping=0.,late_gauge=2000.,tol=2e-13),2000.,detailed=True)
    base=rows['primary']
    physical=('primary','larger','lower_success')
    if {n for n,x in rows.items() if x['expanding_epochs']}!=set(physical) or rows['stalled']['descent_end']:
        raise ArithmeticError('Changed conservative transfer classification')
    errs=dict(refined_efolds=abs(base['expanding_epochs'][0]['efolds']-ref['expanding_epochs'][0]['efolds']),
              refined_peak=abs(base['peak']['sheet_fraction']-ref['peak']['sheet_fraction']),
              refined_return=abs(base['exchange']['descent_end']['scalar_to_sheet']-
                                 ref['exchange']['descent_end']['scalar_to_sheet']))
    if max(errs.values())>5e-6:
        raise ArithmeticError(errs)
    return dict(cases=rows, refinement_errors=errs)
