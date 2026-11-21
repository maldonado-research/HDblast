"""Test the axial stress at the first post-birth reversal of the reference charge current."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy

import numpy as np
from scipy.integrate import cumulative_trapezoid

HERE=Path(__file__).resolve().parent
STRESS=runpy.run_path(str(HERE/'stress.py'))
EVENT=runpy.run_path(str(HERE.parent/'history/event.py'))
History=STRESS['B']['History']


def first_reversal(h):
    t=h.data['t']; u=h.data['u'];D=h.data['D']
    transitions=np.flatnonzero((u[:-1]<0.)&(u[1:]>=0.)&(D[:-1]>0.))
    if not len(transitions):raise ValueError('No post-birth negative-to-positive sign change')
    i=int(transitions[0])
    when=float(np.interp(0.,u[i:i+2],t[i:i+2]))
    return when,i


def integrated_incidences(h,start,stop,panels,order):
    edges=np.linspace(start,stop,panels+1)
    t=(edges[:-1]+edges[1:])/2
    a,D,theta,u=[h.value(t,k) for k in ('a','D','theta','u')]
    # Direct invariant-mass quadrature at every temporal midpoint avoids a second lookup grid.
    thermal=np.array([EVENT['rate'](T,np.arcsinh(v),order=order) for T,v in zip(theta,u)])
    nproper_over_mass2=h.branch*h.action*D/(6*h.mu**3*a**3*np.hypot(1.,u))
    return float(np.diff(edges)@(D*nproper_over_mass2*thermal))


def full_clock(h,half,turn,panels,order):
    t0=h.data['t'][0]
    lo,hi=sorted((turn,half))
    first=integrated_incidences(h,t0,lo,panels,order)
    second=integrated_incidences(h,lo,hi,max(64,panels//4),order)
    amounts={lo:first,hi:first+second}
    return {'reversal':amounts[turn]/h.value(turn,'D'),'half_absolute_work':amounts[half]/h.value(half,'D'),
            'panels_first_segment':panels,'panels_second_segment':max(64,panels//4),'mass_quadrature_order':order}


def run_case(path,metadata,history,name):
    spec=metadata['runs'][name]
    h=History(path/(name+'.bin'),metadata,spec)
    t,i=first_reversal(h); d=h.data; half=history['runs'][name]['absolute_work_50_percent']
    work=cumulative_trapezoid(abs(d['D']*d['u']/np.hypot(1.,d['u'])*d['d']),d['t'],initial=0.)
    fraction=float(np.interp(t,d['t'],work)/spec['abs_electric'])
    clocks=[full_clock(h,half['wall_time'],t,N,k) for N,k in ((384,16),(1536,24))]
    C=max(c['reversal'] for c in clocks)
    error=abs(clocks[0]['reversal']/clocks[-1]['reversal']-1)
    half_error=abs(clocks[-1]['half_absolute_work']/half['mean_events_per_particle']-1)
    if error>.001 or half_error>.001:
        raise ValueError(f'Thermal participation clock underresolved in {name}')
    marker={'wall_time':t,'mean_events_per_particle':C,'distribution_free_zero_collision_lower_bound':max(0,1-C)}
    nominal=STRESS['stress'](h,marker,0.,0.)
    padded=STRESS['stress'](h,marker,.01,500.)
    result={'time':t,'native_step_bracket':d['t'][i:i+2].tolist(),'u_bracket':d['u'][i:i+2].tolist(),
            'temperature_over_m_at_zero':float(h.value(t,'theta')),
            'fraction_of_total_births':float(h.value(t,'D')/spec['end_rest']),
            'absolute_work_fraction':fraction,'relative_clock_mesh_change':error,'relative_half_work_frozen_clock_difference':half_error,
            'half_absolute_work_reference_incidences':half['mean_events_per_particle'],'two_clock_integrations':clocks,
            'cell_enclosure':nominal,'padded_cell_enclosure':padded}
    if t<half['wall_time']:
        # Incidences accrued to t cannot exceed those of the future half-work stage.
        # Round the legacy endpoint and age-integral budget upward before resizing it.
        half_C=max(half['mean_events_per_particle'],1-half['distribution_free_zero_collision_lower_bound'])
        future=half_C*half['rest_at_stage']/h.value(t,'D')
        marker.update(mean_events_per_particle=future,distribution_free_zero_collision_lower_bound=max(0,1-future))
        result['grant_all_later_half_work_events']=STRESS['stress'](h,marker,0.,0.)
    if nominal['longitudinal_stress_floor_over_fluid_upper']<1.9 or padded['longitudinal_stress_floor_over_fluid_upper']<1.7:
        raise ValueError(f'Original stress incompatibility failed in {name}')
    return name,result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--traces',type=Path,required=True);p.add_argument('--history',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=2);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    meta=json.loads((args.traces/'metadata.json').read_text());clock=json.loads(args.history.read_text())
    if clock['action']!=meta['action'] or clock['daughter_mass_L']!=meta['daughter_mass_L']:
        raise ValueError('Different model parameters in trace and reference event clocks')
    with ThreadPoolExecutor(args.jobs) as threads:
        runs=dict(threads.map(lambda n:run_case(args.traces,meta,clock,n),meta['runs']))
    args.output.write_text(json.dumps({'runs':runs},indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
