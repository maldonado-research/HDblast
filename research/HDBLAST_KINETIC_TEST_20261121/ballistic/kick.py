"""Minimum additional transverse impulse for a thermal target on the fixed early wall histories."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import kve


def impulse(theta, tube, fraction):
    z = 1. / theta
    norm = z*z*kve(2, z)
    def tails(y):
        x = math.hypot(z, y)
        return (x*x/norm)*math.exp(z-x)*kve(2, x)
    def density(y):
        x = math.hypot(z, y)
        return y*x*math.exp(z-x)*kve(1, x)/norm
    def integral(func, lo, hi):
        return quad(func, lo, hi, epsabs=1e-13, epsrel=2e-12)[0]
    def quantile(f):
        hi = 2.
        while tails(hi)>1.-f:
            hi *= 2.
        return brentq(lambda y: tails(y)-(1.-f), 0., hi, xtol=1.e-14)
    if not 0.<fraction<1. or theta<=0. or tube<0.:
        raise ValueError("Expected positive temperature, nonnegative tube and interior fraction")
    lo, hi = tube/theta, quantile(fraction)
    norm_test = integral(density, 0., np.inf)
    if not math.isclose(norm_test, 1., rel_tol=2e-12):
        raise ValueError("Transverse density does not integrate to one")
    moment = 0. if lo>=hi else integral(lambda y:(y-lo)**2*density(y), lo, hi)
    check = 0. if lo>=hi else integral(lambda y:2.*(y-lo)*(tails(y)-(1.-fraction)), lo, hi)
    if not math.isclose(check,moment,rel_tol=2.e-11,abs_tol=2e-18):
        raise ValueError("Truncated-moment integrations disagree")
    # Check the analytic tail against a direct density integration at both cutoffs.
    for y,f in [(quantile(fraction*.5),fraction*.5),(hi,fraction)]:
        if not math.isclose(integral(density,0.,y),f,rel_tol=4e-12):
            raise ValueError("Thermal quantile failed")
    return dict(thermal_cutoff_p_perp_over_m=theta*hi,
                rms_extra_impulse_lower_over_m_all_present=theta*math.sqrt(moment),
                supremum_extra_impulse_lower_over_m=max(0.,theta*hi-tube),
                fraction_above_half_quantile_lower=fraction/2.,
                half_quantile_extra_impulse_over_m=max(0.,theta*quantile(fraction/2.)-tube))


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--tube',type=Path,required=True)
    cli.add_argument('--native',type=Path,required=True)
    cli.add_argument('--output',type=Path,required=True)
    args=cli.parse_args()
    data=json.loads(args.tube.read_bytes())['runs']
    meta=json.loads((args.native/'metadata.json').read_bytes())
    col={k:i for i,k in enumerate(meta['columns'])}
    results={}
    for name,run in data.items():
        row=run['early_field_work']['0.01']
        f=row['mean_survival_floor']
        end=row['time']
        a=np.memmap(args.native/(name+'.bin'),dtype='=f8',mode='r').reshape(-1,len(col))
        if len(a)!=meta['runs'][name]['rows']:
            raise ValueError("Incomplete trace")
        first=int(np.searchsorted(a[:,col['D']],0.,side='right'))-1
        after=int(np.searchsorted(a[:,col['t']],end))
        if first<0 or after<=first or after>=len(a):
            raise ValueError("Birth window outside trajectory")
        span=end-float(a[first,col['t']])
        sub=a[first:after+1,:]
        if np.any(np.diff(sub[:,col['a']])<0.):
            raise ValueError("Peak-force lower bound assumes expansion")
        output=impulse(row['theta'],row['transverse_tube_radius'],f)
        upperX=max(float(max(abs(a[first:after,col['X']]))),
                   abs(float(np.interp(end,sub[:,col['t']],sub[:,col['X']]))))
        force=output['supremum_extra_impulse_lower_over_m']/span
        results[name]=dict(marker_time=end,theta=row['theta'],unhit_probability_floor=f,
            birth_tube_over_m=row['transverse_tube_radius'],earliest_birth_left_knot=float(a[first,col['t']]),
            maximum_available_time_over_L=span,peak_extra_force_lower_L_over_m=force,
            equivalent_peak_extra_electric_field_lower=force/meta['charge_to_mass'],
            sampled_peak_parallel_electric_field=upperX,**output)
    args.output.write_text(json.dumps(results,sort_keys=True,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
