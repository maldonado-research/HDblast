"""Extremal electric current and longitudinal stress from compulsory ballistic births."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np
from scipy.integrate import cumulative_trapezoid

HERE=Path(__file__).resolve().parent
B=runpy.run_path(str(HERE/'ballistic.py'))
Q_BIRTH=math.sqrt(1.05**2-1)


def cheapest(cost,weights,floor,available):
    """Minimize cost for surviving mass >=floor and find minimum hits for a cost budget."""
    order=np.argsort(cost)
    values=cost[order]; w=weights[order]
    cumulative=np.cumsum(w); partial=np.cumsum(w*values)
    portion=np.maximum(0.,np.minimum(w,floor-(cumulative-w)))
    minimum=float(np.sum(portion*values))
    if available>=partial[-1]:
        return minimum,0.
    j=int(np.searchsorted(partial,available,side='right'))
    before=cumulative[j-1] if j else 0.
    past=partial[j-1] if j else 0.
    allowed=before+(available-past)/values[j]
    return minimum,max(0.,1.-allowed)


def stress(history,marker,uncertainty,extra_impulse):
    data=history.data;end=marker['wall_time']
    t=np.append(data['t'][data['t']<end],end)
    force=np.interp(t,data['t'],data['d']*data['a']) # native d=Q*X
    increments=np.diff(t)*np.maximum(abs(force[1:]),abs(force[:-1]))
    impulse=cumulative_trapezoid(force,t,initial=0.)
    # Intersect forward and backward reachability intervals within each native time cell.
    low=np.maximum(impulse[:-1],impulse[1:])-increments
    high=np.minimum(impulse[:-1],impulse[1:])+increments
    radius=history.value(end,'a')
    zlo=(impulse[-1]-high)/radius-Q_BIRTH-extra_impulse
    zhi=(impulse[-1]-low)/radius+Q_BIRTH+extra_impulse
    if np.any(zlo>zhi) or np.max(data['R'])>1.05+1e-12:
        raise ValueError('Invalid reach bounds on the birth momenta')
    sign_uncertain=(zlo<=0)&(zhi>=0)
    min_z=np.where(sign_uncertain,0.,np.minimum(abs(zlo),abs(zhi)))
    transverse2=Q_BIRTH**2
    vmax=zhi/np.sqrt(1+zhi*zhi+np.where(zhi<0.,transverse2,0.))
    vmin=zlo/np.sqrt(1+zlo*zlo+np.where(zlo>0.,transverse2,0.))
    pz_min=min_z**2/np.sqrt(1+min_z**2+transverse2)
    energy_less_transverse=(1+min_z**2)/np.sqrt(1+min_z**2+transverse2)
    d=history.value(t,'D'); w=np.diff(d)/d[-1]
    if np.min(w)<0. or abs(w.sum()-1)>1e-12: raise ValueError('Invalid birth weights')
    # Take the most permissive fluid moments across both native endpoints, with
    # Cartesian bounding of theta and |u| to allow their simultaneous variation.
    i=np.searchsorted(data['t'],end)
    i=min(i,len(data['t'])-1)
    us=data['u'][i-1:i+1]; ts=data['theta'][i-1:i+1]
    umax=np.max(np.abs(us)); gmax=np.hypot(1.,umax); theta=np.max(ts)
    h=2.5*theta+np.hypot(1.,1.5*theta)
    eupper=h*gmax-theta/gmax; zzupper=(h*umax**2+theta)/gmax
    transverse_lower=2*np.min(ts)/gmax
    fluid_v=us/np.hypot(1.,us); vlo=float(np.min(fluid_v)); vhi=float(np.max(fluid_v))
    budget=max(marker['mean_events_per_particle'],1.-marker['distribution_free_zero_collision_lower_bound'])
    budget*=1.+uncertainty
    f=max(0.,1.-budget)
    minv,req_vlo=cheapest(1+vmin,w,f,1+vhi)
    maxv,req_vhi=cheapest(1-vmax,w,f,1-vlo)
    zz,req_zz=cheapest(pz_min,w,f,zzupper)
    et,req_tr=cheapest(energy_less_transverse,w,f,eupper-transverse_lower)
    # The pressure and current constraints apply to the SAME unhit set. A positive
    # sum of their costs gives a necessary joint condition. Set its scale to the
    # measured mean lab energy to keep the current term in stress units.
    if abs(vhi)<abs(vlo):
        joint_cost=pz_min+eupper*(1+vmin)
        joint_budget=zzupper+eupper*(1+vhi)
    else:
        joint_cost=pz_min+eupper*(1-vmax)
        joint_budget=zzupper+eupper*(1-vlo)
    joint_minimum,req_joint=cheapest(joint_cost,w,f,joint_budget)
    return dict(joint_witness_floor=joint_minimum,joint_witness_fluid_upper=joint_budget,
                minimum_distinct_hit_fraction_for_joint_witness=req_joint,
                upper_collision_incidences=budget,mandatory_unhit_fraction=f,
                additional_impulse_uncertainty=extra_impulse,extra_relative_event_allowance=uncertainty,
                maximum_native_cell_impulse=float(np.max(increments)/radius),
                allowed_current_velocity_low=minv-1,allowed_current_velocity_high=1-maxv,
                fluid_current_velocity_low=vlo,fluid_current_velocity_high=vhi,
                current_gap=max(0.,minv-1-vhi,vlo-(1-maxv)),
                unhit_longitudinal_stress_floor=zz,fluid_longitudinal_stress_upper=float(zzupper),
                longitudinal_stress_floor_over_fluid_upper=zz/zzupper,
                allowed_transverse_stress_upper_if_total_energy_matches=eupper-et,
                fluid_transverse_stress_lower=float(transverse_lower),
                fluid_total_energy_upper=float(eupper),
                minimum_hit_fraction_for_current=max(req_vlo,req_vhi),
                minimum_hit_fraction_for_longitudinal_stress=req_zz,
                minimum_hit_fraction_for_transverse_stress_at_fixed_total_energy=req_tr,
                minimum_collision_incidences_for_both_current_and_stress=max(req_vlo,req_vhi,req_zz,req_tr,req_joint))


def study(path,clock):
    m=json.loads((path/'metadata.json').read_text());r=json.loads(clock.read_text())
    out={}
    for name,run in m['runs'].items():
        history=B['History'](path/(name+'.bin'),m,run)
        result={}
        for field in ('absolute_work_10_percent','absolute_work_50_percent','absolute_work_90_percent','first_hot_time'):
            marker=r['runs'][name][field]
            result[field]=dict(tau=marker['wall_time'],cell_enclosure=stress(history,marker,0.,0.),
                               sensitivity=stress(history,marker,.01,500.))
        out[name]=result
    for name in ('dedicated_native','dedicated_half'):
        sharp=out[name]['absolute_work_50_percent']['cell_enclosure']
        loose=out[name]['absolute_work_50_percent']['sensitivity']
        if sharp['current_gap']<.065 or sharp['longitudinal_stress_floor_over_fluid_upper']<1.13:
            raise ValueError('Primary moment deficit changed')
        if loose['current_gap']<.01 or loose['longitudinal_stress_floor_over_fluid_upper']<1.05:
            raise ValueError('Padded moment deficit changed')
    return dict(runs=out)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--traces',type=Path,required=True)
    parser.add_argument('--history',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(study(args.traces,args.history),indent=2,sort_keys=True,allow_nan=False)+'\n')
