"""Independent checks on the event kernel and refined event-budget landmarks."""
import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import event
import optical


def mom(theta,uniform):
    # Rest-frame d^3p Boltzmann measure, sampled in sqrt((E/m-1)/theta).
    t = np.linspace(0,math.sqrt(55.),16001)
    x = t*t
    pdf = 2*x*(1+theta*x)*np.sqrt(2+theta*x)*np.exp(-x)
    cumulative = cumulative_trapezoid(pdf,t,initial=0.)
    x = np.interp(uniform,cumulative/cumulative[-1],t)**2
    return 1+theta*x,np.sqrt(theta*x*(2+theta*x))


def sample(theta,r,seed):
    seq = qmc.Sobol(d=5,scramble=True,seed=seed).random_base2(17)
    a,p = mom(theta,seq[:,0])
    b,q = mom(theta,seq[:,1])
    i,j = 2*seq[:,2]-1,2*seq[:,3]-1
    transverse = p*q*np.sqrt((1-i*i)*(1-j*j))*np.cos(2*np.pi*seq[:,4])
    # E1*E2-pz1*pz2 after opposite boosts, directly in native rest-frame variables.
    dot = (a*b-i*p*j*q)*math.cosh(2*r) + (i*p*b-j*q*a)*math.sinh(2*r) - transverse
    invariant = np.maximum(4.,2*(1+dot))
    beta = np.sqrt((invariant-4)/invariant)
    h = .30*(invariant+2)
    flux_xs = 6*np.pi*beta*h*h/((invariant-8300.**2)**2+(h*beta)**2)
    return float(np.mean(flux_xs/(a*b))/math.cosh(r))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native',type=Path,required=True)
    parser.add_argument('--history',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=2)
    args=parser.parse_args()
    meta = json.loads((args.native/'metadata.json').read_bytes())
    ref = json.loads(args.history.read_bytes())
    samples = {}
    for theta,r in [(.01,.2),(.001,9.),(10.,4.),(3000.,.7),(3000.,.02)]:
        draws = np.array([sample(theta,r,s) for s in (79,191,223,421,563,691,761,881)])
        direct = event.rate(theta,r,32)
        err = float(np.std(draws,ddof=1)/math.sqrt(len(draws))/direct)
        diff = float(draws.mean()/direct-1)
        if abs(diff)>max(.003,6*err):
            raise ValueError('Independent sample and analytic integral disagree')
        samples[f'theta_{theta}_r_{r}'] = dict(invariant_integral=direct,samples=draws.tolist(),
              relative_sampling_standard_error=err,signed_relative_difference=diff)
    # Refine both the lookup and the invariant-mass quadrature.
    lookup = event.Lookup(*event.make_table(.0625,.05,24,args.jobs))
    full = {}
    for name in meta['runs']:
        rows = optical.evaluate(args.native,name,meta,lookup)
        full[name] = {}
        for stage in ('first_hot_time','absolute_work_50_percent','first_mean_one_event','doubling_time'):
            old = ref['runs'][name][stage]
            vals = {k:v for k,v in rows[stage].items() if k in
                    ('fraction_zero','mean_events_per_particle','absolute_work_fraction')}
            full[name][stage] = {k:dict(value=v,relative_change=v/old[k]-1) for k,v in vals.items()}
    limits = {}
    constant = lambda eta:1.
    for theta,r in [(1e-7,0.),(1e-7,.7),(1e4,0.),(1e4,2.)]:
        exact = ((4*math.sqrt(theta/math.pi) if r==0 else 2*math.sinh(r))
                 if theta<1 else math.cosh(r)*(1+math.tanh(r)**2))
        numeric = event.rate(theta,r,32,cross_section=constant)
        if abs(numeric/exact-1)>.0001:
            raise ValueError('Asymptotic constant-cross-section limit failed')
        limits[f'theta_{theta}_r_{r}'] = dict(quadrature=numeric,asymptotic=exact)
    checks = dict(sobol_samples_per_scramble=2**17,independent_pair_sampling=samples,
        refined_table_logT_step=.0625,refined_table_scaled_rapidity_step=.05,
        refined_quadrature_order=24,refined_clock_comparison=full,constant_sigma_limits=limits)
    args.output.write_text(json.dumps(checks,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
