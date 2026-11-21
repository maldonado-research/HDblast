"""Test the compulsory unscattered birth component against a locally thermal transverse distribution."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad

HERE=Path(__file__).resolve().parent
B=runpy.run_path(str(HERE/'ballistic.py'))
MARKERS=('absolute_work_10_percent','absolute_work_50_percent','absolute_work_90_percent','first_hot_time')


def spherical_probability(theta,q):
    # Independent direction average in a spherical, unboosted Maxwell-Jüttner distribution.
    # A longitudinal boost cancels out of this transverse number marginal by parity.
    def density(x):
        return (1+theta*x)*math.sqrt(x*(2+theta*x))*math.exp(-x)
    def cap(x):
        p2=theta*x*(2+theta*x)
        if p2<=q*q:
            return 1.
        y=q*q/p2
        return y/(1.+math.sqrt(1.-y))
    xq=q*q/((1.+math.hypot(1.,q))*theta)
    edges=np.geomspace(xq,80.,32)
    ref=quad(lambda x:density(x)*cap(x),0.,xq,epsabs=1e-18)[0]
    for low,high in zip(edges[:-1],edges[1:]):
        ref+=quad(lambda x:density(x)*cap(x),low,high,epsabs=1e-18,epsrel=3e-11)[0]
    return ref/quad(density,0.,80.,epsabs=1e-15,epsrel=3e-12)[0]


def study(native,history):
    metadata=json.loads((native/'metadata.json').read_text())
    budget=json.loads(history.read_text())
    if (budget['action'],budget['daughter_mass_L'])!=(metadata['action'],metadata['daughter_mass_L']):
        raise ValueError('Event clock and wall history use different mass or density')
    result={}
    for name,spec in metadata['runs'].items():
        trace=B['History'](native/(name+'.bin'),metadata,spec)
        rows={}
        for stage in MARKERS:
            count=budget['runs'][name][stage]
            end=count['wall_time']
            bound=B['bath_independent_floor'](trace,end)
            thermal=bound['maxwell_juttner_tube_fraction']
            numerical=spherical_probability(bound['theta'],bound['transverse_tube_radius'])
            if not math.isclose(numerical,thermal,rel_tol=2e-9):
                raise ArithmeticError('Independent angular integral disagrees with Bessel marginal')
            encounter=count['mean_events_per_particle']
            floor=count['distribution_free_zero_collision_lower_bound']
            if abs(floor-max(0.,1-encounter))>3e-6:
                raise ArithmeticError('Event-budget identity disagrees with integrated birth clock')
            rows[stage]=dict(
                time=end,temperature_over_mass=bound['theta'],tube_radius_over_mass=bound['transverse_tube_radius'],
                maxwell_juttner_tube_fraction=thermal,independent_spherical_integral=numerical,
                thermal_mean_collision_incidences_per_particle=encounter,
                minimum_tube_fraction_given_thermal_event_budget=floor,
                total_variation_floor_given_thermal_event_budget=max(0.,floor-thermal),
                minimum_incidences_for_exact_thermal_tube=1-thermal,
                all_target_rate_upper_bound=bound['global_pair_kernel_upper'],
                minimum_tube_fraction_for_any_target_momenta=bound['mean_survival_floor'],
                any_target_minimum_over_thermal=bound['mandatory_to_thermal_tube_population_ratio'],
                any_target_total_variation_floor=bound['total_variation_floor'],
                any_target_relative_entropy_floor_nats=bound['relative_entropy_floor_nats'],
            )
        d=trace.data
        work_rate=np.abs(metadata['charge_to_mass']*d['D']*d['u']/np.hypot(1.,d['u'])*d['X'])
        work=cumulative_trapezoid(work_rate,d['t'],initial=0.)
        early={}
        # Extend the original clock markers to round early work fractions. No thermal-event input.
        for frac in (.001,.01,.05):
            end=float(np.interp(frac*spec['abs_electric'],work,d['t']))
            bound=B['bath_independent_floor'](trace,end)
            thermal=bound['maxwell_juttner_tube_fraction']
            numerical=spherical_probability(bound['theta'],bound['transverse_tube_radius'])
            if not math.isclose(numerical,thermal,rel_tol=2e-9):
                raise ArithmeticError('Early independent spherical integral disagrees')
            early[str(frac)]=dict(time=end,fraction_of_absolute_electric_work=frac,
                                 fraction_of_all_run_births=float(trace.value(end,'D')/spec['end_rest']),
                                 **bound)
        result[name]=dict(stages=rows,early_field_work=early)
    return {'runs':result,'event_markers':'Fractions of whole-history absolute electrical work, or the reference first hot crossing.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--traces',type=Path,required=True)
    p.add_argument('--history',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(study(args.traces,args.history),indent=2,sort_keys=True,allow_nan=False)+'\n')
