"""How much gauge contrast fits in a positive, slow-response-matched action?"""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=48
MARGIN=.05
EDGE=.2


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def sum_inertia(k):
    return math.fsum(l-1/l for l in range(REFERENCE+1,k+1))


def taper(u):
    """Unit plateau with smooth rising/falling edges, each one fifth of the interval."""
    if u<=0 or u>=1:
        return 0.,0.
    if EDGE<=u<=1-EDGE:
        return 1.,0.
    t=min(u,1-u)/EDGE
    if t<1e-8:
        return 0.,0.
    if 1-t<1e-8:
        return 1.,0.
    y=expit(1/(1-t)-1/t)
    dydt=y*(1-y)*(1/t**2+1/(1-t)**2)
    return float(y),float(dydt/EDGE if u<.5 else -dydt/EDGE)


def geometry(branch,points):
    qstart=points['switch_start']['q']
    width=points['switch_end']['q']-qstart
    y0=[points['switch_start'][k] for k in ('velocity','a','h','x','eta')]
    def rhs(u,y):
        v,a,h,x,eta=y
        time=width/v
        force=branch.potential(qstart+width*u)[1]
        return time*np.array([force-3*h*v,a*h,1/a**2-v*v/2,1.,1/a])
    sol=solve_ivp(rhs,(0.,1.),y0,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01,dense_output=True)
    if not sol.success: raise ArithmeticError(sol.message)
    err=float(max(abs(sol.y[:,-1]-[points['switch_end'][k] for k in ('velocity','a','h','x','eta')])))
    if err>2e-8:
        raise ArithmeticError('Field-coordinate solution missed switch end')
    grid=np.linspace(0.,1.,513)
    a=CubicSpline(grid,sol.sol(grid)[1])
    mid=(grid[1:]+grid[:-1])/2
    spline_error=float(np.max(abs(a(mid)-sol.sol(mid)[1])))
    if spline_error>2e-10:raise ArithmeticError('Radius interpolation differs')
    return a,width,err,spline_error


def probe_reach(branch,archive,matched,gauge,free):
    pts=archive['points']
    a,width,err,interpolation=geometry(branch,pts)
    eps=branch.fitted.Hpl**2
    ai=pts['launch']['a']
    ui=pts['launch']['parent_over_V0']
    qstart=pts['switch_start']['q']
    area=float(a.integrate(0.,1.))
    weighted=quad(lambda u:float(a(u))*taper(u)[0],0.,1.,points=[EDGE,1-EDGE],epsabs=1e-13)[0]
    def flat_slope(q,q0,dq,delta):
        u=(q-qstart)/width
        t,dt=taper(u)
        if t==0. and dt==0.:
            return 0.,0.
        scale=delta/(width*weighted)
        return scale*float(a(u))*t,scale*(float(a(u,1))*t+float(a(u))*dt)/width
    def bump_slope(u):
        return gauge.raw_bump(float(u))/(gauge.NORM*abs(width))
    opt=minimize_scalar(lambda u:-bump_slope(u)/float(a(u)),bounds=(.35,.65),method='bounded',
                        options={'xatol':2e-13})
    maximum_ratio=-opt.fun
    primary_factor=1/(maximum_ratio*abs(width)*area)
    taper_factor=weighted/area
    rows={}
    for k in [49,50,51,52,53,54,56,60,64,72,80,96,128,192,256,384,512,640,768]:
        c=eps*sum_inertia(k)/(2*math.pi**2)
        scale=math.sqrt((1-MARGIN)/c)
        cap=abs(width)*area*scale
        # A necessary kinematic time to span this contrast on any first expansion
        # with nonnegative photon energy and positive scalar inertia. It uses a_end
        # to bound a throughout, so it can only weaken the local margin condition.
        delta=3.
        a_by_contrast=abs(delta)/(abs(width)*scale)
        a_at_peak=abs(delta)*bump_slope(.5)/scale
        def min_time(radius):
            return math.acosh(max(1.,radius/ai))/math.sqrt(ui)
        rows[str(k)]=dict(c=c,upper_contrast_on_reference=cap,original_contrast_limit=scale/maximum_ratio,
                          tapered_contrast_limit=abs(width)*weighted*scale,
                          strong_minimum_end_radius_any_profile=a_by_contrast,
                          strong_earliest_end_any_profile=min_time(a_by_contrast),
                          strong_minimum_midpoint_radius_original=a_at_peak,
                          strong_earliest_midpoint_original=min_time(a_at_peak))
        if abs(cap*primary_factor-rows[str(k)]['original_contrast_limit'])>1e-12:
            raise ArithmeticError('Original profile capacity mismatch')
    integral=quad(lambda u: flat_slope(qstart+width*u,0.,0.,1.)[0]*width,0.,1.,
                  points=[EDGE,1-EDGE],epsabs=1e-13)[0]
    if abs(integral-1)>5e-12 or not 0<primary_factor<taper_factor<1:
        raise ArithmeticError('Taper normalization or upper capacity failed')
    # Fixed field-only profile from the old trajectory, retained unchanged when the
    # matched photons alter the time history, expansion, or direction of the field.
    trials={}
    coasts={}
    for delta,cutoffs in [(.1,[48,64,96,128,192]),(3.,[48,49,52,54,56,60,64,72,80,96])]:
        for k in cutoffs:
            record=matched.run(branch,k,delta,profile=flat_slope)
            trials[f'{k}/{delta:g}']=dict(status=record['status'],minimum_z=record['minimum_bare_inertia'],
                                          end=record['end'],errors=record['errors'])
            if delta==3. and k in (48,56,60,64) and record['status']=='first_minimum':
                expected=dict(record['end'],fraction_of_matter=record['end']['photon_fraction'])
                coasts[str(k)]=free.coast(branch,record['switch_end'],expected)
    checks={}
    for k,delta in [(64,.1),(60,3.),(72,3.)]:
        direct=matched.run(branch,k,delta,profile=flat_slope,direct=True)
        trial=trials[f'{k}/{delta:g}']
        difference=max(abs(trial['end'][f]-direct['end'][f]) for f in ('photon_fraction','a','x','progress','z'))
        if direct['status']!=trial['status'] or difference>2e-7:raise ArithmeticError('Tapered trajectory depends on coordinates')
        checks[f'{k}/{delta:g}']=difference
    return dict(absolute_field_width=abs(width),mean_reference_radius=area,original_capacity_fraction=primary_factor,
                tapered_capacity_fraction=taper_factor,geometric_replay_error=err,maximum_radius_interpolation_error=interpolation,
                cutoffs=rows,tapered_profile_runs=trials,strong_ungated_intervals=coasts,independent_trajectory_error=checks)


def calculate(prior):
    gauge=module('reference_gauge_solver',prior/'gauge/gauge.py')
    matched=module('matched_dynamics',ROOT/'matched/matched.py')
    free=module('unforced_coast',prior/'feedback/coast.py')
    archive=json.loads((prior/'gauge/GAUGE.json').read_bytes())['branches']
    return dict(reference_cutoff=REFERENCE,reserve=MARGIN,taper_end_width=EDGE,
                branches={k:probe_reach(b,archive[k],matched,gauge,free) for k,b in matched.branches(prior).items()})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior',type=Path,default=ROOT/'prior')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior),indent=2,allow_nan=False)+'\n')
