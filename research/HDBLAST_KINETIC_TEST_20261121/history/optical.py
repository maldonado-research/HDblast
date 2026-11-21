"""Integrate number-weighted event clocks on native directional-birth histories."""
import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import gammaincc

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import event


def stage(t,rest,clock,when,final_rest,converted_number):
    stop = int(np.searchsorted(t,when))
    if not 0<stop<len(t):
        raise ValueError('Stage outside trace')
    r = np.r_[rest[:stop],np.interp(when,t[stop-1:stop+1],rest[stop-1:stop+1])]
    c = np.r_[clock[:stop],np.interp(when,t[stop-1:stop+1],clock[stop-1:stop+1])]
    mass = np.diff(r)
    born = np.flatnonzero(mass>0)
    if r[0]!=0 or np.any(mass<0) or len(born)==0:
        raise ValueError('Invalid rest source')
    weights = mass[born]/r[-1]
    aged = c[-1]-.5*(c[born]+c[born+1])
    if np.any(aged<0) or abs(weights.sum()-1)>1e-12:
        raise ValueError('Invalid cohort probabilities')
    qrows = {}
    for p in (.05,.5,.95):
        j = min(int(np.searchsorted(r,p*r[-1])),len(r)-1)
        qrows[str(p)] = dict(birth_wall_time=float(t[j]),event_depth_to_stage=float(c[-1]-c[j]))
    return dict(wall_time=when,rest_at_stage=float(r[-1]),number_in_whole_wall=float(r[-1]*converted_number),
                fraction_of_run_final_births=float(r[-1]/final_rest),clock_since_first_birth=float(c[-1]),
                mean_events_per_particle=float(weights@aged),fraction_zero=float(weights@np.exp(-aged)),
                fraction_less_than_two=float(weights@gammaincc(2,aged)),
                fraction_less_than_five=float(weights@gammaincc(5,aged)),cohort_quantiles=qrows)


def evaluate(native,name,meta,lookup):
    raw = np.memmap(native/(name+'.bin'),dtype=np.float64).reshape(-1,len(meta['columns']))
    columns = {k:raw[:,j] for j,k in enumerate(meta['columns'])}
    t,r,D,theta,u = [columns[key] for key in ('t','a','D','theta','u')]
    factor = meta['runs'][name]['channel_fraction']
    mu = meta['daughter_mass_L']
    emitted = 2*math.pi**2*meta['action']/mu
    # K = m² Gamma/n, with n the density of the opposite flavor and sign in its rest frame.
    coeff = lookup.evaluate(theta,u)
    rate = factor*(meta['action']/(6*mu**3))*(D/r**3/np.hypot(1.,u))*coeff
    clock = cumulative_trapezoid(rate,t,initial=0)
    events = {n:meta['runs'][name][n] for n in ['first_hot_time','doubling_time','maximum_tau']}
    for p in (.1,.5,.9):
        events['births_'+str(round(100*p))+'_percent'] = float(np.interp(p*D[-1], D, t))
    work_rate = np.abs(meta['charge_to_mass']*columns['X']*D*u/np.hypot(1.,u))
    work = cumulative_trapezoid(work_rate,t,initial=0)
    work_rel_error = float(abs(work[-1]/meta['runs'][name]['abs_electric']-1))
    if work_rel_error>.003:
        raise ValueError('Absolute electric work differs from recorded solver ledger')
    for p in (.1,.5,.9):
        events['absolute_work_'+str(round(100*p))+'_percent'] = float(np.interp(p*work[-1],work,t))
    rows = {k:stage(t,D,clock,t0,D[-1],emitted) for k,t0 in events.items()}
    # Counting all encounters gives a separate integral identity for the weighted mean
    # even if the individual counters have a non-Poisson distribution.
    encounters = cumulative_trapezoid(D*rate,t,initial=0)
    for row in rows.values():
        when = row['wall_time']
        mean = float(np.interp(when,t,encounters)/row['rest_at_stage'])
        if not math.isclose(mean,row['mean_events_per_particle'],rel_tol=1e-6,abs_tol=1e-8):
            raise ValueError('Integrated event budget disagrees with cohort average')
        row['distribution_free_zero_collision_lower_bound'] = max(0.,1.-mean)
        row['absolute_work_fraction_at_stage'] = float(np.interp(when,t,work)/work[-1])
    pairs = np.divide(encounters,D,out=np.zeros_like(encounters),where=D>0)
    crossed = np.flatnonzero(pairs>=1)
    if len(crossed)==0:
        raise ValueError('Mean event count never reaches one')
    i = crossed[0]
    rows['first_mean_one_event'] = dict(
        wall_time=float(np.interp(1,pairs[i-1:i+1],t[i-1:i+1])),
        birth_fraction=float(np.interp(1,pairs[i-1:i+1],D[i-1:i+1])/D[-1]),
        absolute_work_fraction=float(np.interp(1,pairs[i-1:i+1],work[i-1:i+1])/work[-1]),
        absolute_work_relative_ledger_error=work_rel_error)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--step',type=float,default=.125,help='log10 temperature grid step')
    p.add_argument('--dy',type=float,default=.1,help='scaled rapidity grid step')
    p.add_argument('--order',type=int,default=16)
    args = p.parse_args()
    meta = json.loads((args.native/'metadata.json').read_bytes())
    ts,ys,coeff = event.make_table(args.step,args.dy,args.order,args.jobs)
    lookup = event.Lookup(ts,ys,coeff)
    out = dict(action=meta['action'],daughter_mass_L=meta['daughter_mass_L'],
               log_temperature_step=args.step,scaled_rapidity_step=args.dy,quadrature_order=args.order,
               runs={})
    for name in meta['runs']:
        out['runs'][name] = evaluate(args.native,name,meta,lookup)
        print(name, out['runs'][name]['first_hot_time']['fraction_zero'], flush=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
