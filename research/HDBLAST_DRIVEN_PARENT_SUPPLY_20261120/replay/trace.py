"""Reconstruct the sourced homogeneous parent on the measured wall and scalar trajectory."""
import argparse
from bisect import bisect_right
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def step(s):
    if s<=0: return 0.,0.
    if s>=1: return 1.,0.
    z=1/(1-s)-1/s
    v=math.exp(-abs(z))
    return (1 if z>=0 else v)/(1+v), v/(1+v)**2*((1-s)**-2+s**-2)


def cubic(y0,y1,d0,d1,h,x):
    slope=h*d0
    quadratic=3*(y1-y0)-h*(2*d0+d1)
    cubic_coeff=h*(d0+d1)-2*(y1-y0)
    return y0+x*(slope+x*(quadratic+x*cubic_coeff)),(slope+x*(2*quadratic+3*x*cubic_coeff))/h


def propagate(chi,mom,mass,volume,minimum,dt):
    phase=mass*dt
    s=math.sin(phase)
    c=math.cos(phase)
    inertia=mass*volume
    return minimum+c*(chi-minimum)+s*mom/inertia, c*mom-s*inertia*(chi-minimum)


def independent(rows,metrics,parameters,n,spacing):
    data=rows[::spacing]
    if data[-1] is not rows[-1]:
        data.append(rows[-1])
    times=[row['t'] for row in data]
    def fields(t):
        j=min(len(data)-2,max(0,bisect_right(times,t)-1))
        l,r=data[j:j+2]
        h=r['t']-l['t']
        u=(t-l['t'])/h
        a,adot=cubic(l['a'],r['a'],l['a']*l['H'],r['a']*r['H'],h,u)
        psi,pi=cubic(l['psi'],r['psi'],l['pi'],r['pi'],h,u)
        half=parameters['width']/2
        x=(psi-parameters['onset'])/half
        rise,d1=step(x); fall,d2=step(x-1)
        v=parameters['g']*(rise-fall)
        vdot=parameters['g']*(d1-d2)*pi/half
        gain,d3=step((psi-parameters['gate'])/parameters['gate_width']+.5)
        dm=parameters['m_high']-parameters['m_low']
        mass=2*parameters['daughter']*(parameters['m_low']+dm*gain)
        mdot=2*parameters['daughter']*dm*d3*pi/parameters['gate_width']
        return a**3,adot/a,mass,mdot,v,vdot

    events={metrics[s+'_time']:s for s in ('drive_end','before_mass','opening')}
    endpoints=sorted(set(times+list(events)))
    chi=mom=wdrive=wmass=wpress=0.
    checks={}
    target=.0008/(2*parameters['daughter']*parameters['m_low'])
    for begin,end in zip(endpoints,endpoints[1:]):
        dt=(end-begin)/n
        for i in range(n):
            volume,H,mass,mdot,v,vdot=fields(begin+(i+.5)*dt)
            x,p=propagate(chi,mom,mass,volume,v,dt/2)
            delta=x-v
            wdrive-=dt*volume*mass**2*delta*vdot
            wmass+=dt*volume*mass*mdot*delta**2
            wpress+=dt*3*H*volume*.5*((p/volume)**2-(mass*delta)**2)
            chi,mom=propagate(chi,mom,mass,volume,v,dt)
        if end in events or end==endpoints[-1]:
            volume,H,mass,mdot,v,_=fields(end)
            energy=.5*((mom/volume)**2+(mass*(chi-v))**2)*volume
            action=energy/mass
            name=events.get(end,'last_sample')
            checks[name]=dict(t=end,chi=chi,chi_dot=mom/volume,mass=mass,action=action,
                              ratio=action/target, pressure_work=wpress, drive_work=wdrive,
                              mass_work=wmass, comoving_energy=energy,
                              ledger_residual=energy-wdrive-wmass+wpress)
    return checks


def calculate(doc):
    metrics=doc['runs']['primary']
    columns=doc['trace_columns']
    rows=[dict(zip(columns,r,strict=True)) for r in metrics['trace']]
    config=doc['parameters']['primary']
    params=dict(g=config['drive_amplitude'],width=config['drive_width'],onset=config['drive_start'],
                gate=config['decay_threshold'],gate_width=config['decay_width'],daughter=config['daughter_mass'],
                m_low=config['initial_mass_ratio'],m_high=config['final_mass_ratio'])
    if rows[0]['chi'] != 0. or rows[0]['chi_dot'] != 0.:
        raise ValueError('Comparison requires an initially empty homogeneous field.')
    full=independent(rows,metrics,params,32,1)
    half=independent(rows,metrics,params,16,1)
    sparse=independent(rows,metrics,params,32,2)
    errors={s:(c['ratio']/metrics[s+'_target_ratio']-1) for s,c in full.items()
            if s!='last_sample'}
    native = rows[-1]
    final_delta = {key:full['last_sample'][key]-native[key] for key in
                   ('chi','chi_dot','drive_work','mass_work')}
    checks=dict(
        relative_event_ratio_discrepancies=errors, final_minus_native=final_delta,
        max_ratio_change_half_steps=max(abs(c['ratio']/half[k]['ratio']-1) for k,c in full.items()),
        max_ratio_change_half_source=max(abs(c['ratio']/sparse[k]['ratio']-1) for k,c in full.items()),
        max_absolute_ledger_residual=max(abs(c['ledger_residual']) for c in full.values()))
    if max(map(abs,errors.values()))>2e-4 or checks['max_absolute_ledger_residual']>2e-8:
        raise ValueError(checks)
    phase_multiplier=(metrics['before_mass_radius']/metrics['opening_radius'])**3 * (
                      metrics['before_mass_mass']/metrics['opening_mass'])
    initial=full['before_mass']['ratio']
    bound=[initial*phase_multiplier,initial/phase_multiplier]
    if not bound[0]<=full['opening']['ratio']<=bound[1]:
        raise ValueError('Opening action violates the phase-independent finite-time bound.')
    target=.0008/(2*params['daughter']*params['m_low'])
    free_ratios=[row['a']**3*row['rho_chi']/row['M']/target for row in rows
                 if metrics['drive_end_time']<=row['t']<=metrics['before_mass_time']]
    return dict(parameters=params,checks=checks,epochs=full,
                pre_mass_to_opening_phase_independent_bounds=bound,
                postdrive_before_mass_sampled_ratio_range=[min(free_ratios),max(free_ratios)],
                native_end_wkb_diagnostics=dict(H_over_M=metrics['opening_H_over_M'],
                                                abs_Mdot_over_M2=metrics['opening_mass_WKB']))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=HERE.parent/'coupled/PARENT.json')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(json.loads(args.source.read_text())),
                                      indent=2,sort_keys=True,allow_nan=False)+'\n')
