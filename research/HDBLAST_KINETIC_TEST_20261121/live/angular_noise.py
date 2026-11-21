"""Single-snapshot angular-pressure calibration for an independent isotropic relativistic gas."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import brentq
from scipy.special import kv, roots_legendre

THETA=2500.
REPLICATES=4096
U,W=roots_legendre(96)
U=(1+U)/2
W=W/2


def massless_tail(n,q,upper,observable="pressure"):
    # Gil-Pelaez integration for sum_i y_i(u_i²-q), y~Gamma(3,1).
    directions=U*U if observable=="pressure" else 2*U-1
    def integrand(t):
        if t==0.:
            return 3*n*((1/3 if observable=="pressure" else 0)-q)
        cf=W@((1-1j*t*(directions-q))**-3)
        return float(np.imag(cf**n)/t)
    val,err=quad(integrand,0.,2.,epsabs=2e-13,epsrel=2e-13,limit=300)
    # For n>=450 the tail beyond t=2 is numerically negligible; extend it as a check.
    rem=quad(integrand,2.,10.,epsabs=2e-15)[0]
    p=.5+(val+rem)/np.pi if upper else .5-(val+rem)/np.pi
    return dict(probability=float(p),integration_error_reported=float(err/np.pi),extended_integral=float(rem/np.pi))


def tail(n,upper,seed):
    q=1/3+(.05 if upper else -.05)
    a=U*U-q
    lam=brentq(lambda v:float(W@(a/(1-v*a)**4)),-1.,1.)
    normal=float(W@((1-lam*a)**-3))
    direction=np.linspace(0,1,16385)
    cdf=cumulative_trapezoid((1-lam*(direction**2-q))**-3,direction,initial=0.)
    cdf/=cdf[-1]
    rng=np.random.default_rng(seed)
    delta=1/THETA
    mass_norm=.5*delta**2*kv(2,delta)
    known=[]; fitted=[]
    for start in range(0,REPLICATES,32):
        shape=(min(32,REPLICATES-start),n)
        u=np.interp(rng.uniform(size=shape),cdf,direction)
        u*=np.where(rng.uniform(size=shape)<.5,-1.,1.)
        alignment=u*u-q
        y=rng.gamma(3.,size=shape)/(1-lam*alignment)
        e=np.hypot(y,delta)
        weight=y*y/e
        zz=np.sum(weight*u*u,axis=1); trace=weight.sum(axis=1)
        flow=(y*u).sum(axis=1); total=e.sum(axis=1)
        pzz=.5*(np.sqrt((total+zz)**2-4*flow**2)-(total-zz))
        lab=zz/trace; frame=pzz/(pzz+trace-zz)
        # Exact likelihood ratio from tilted massless to finite-mass rest MJ.
        logL=(n*(math.log(normal)-math.log(mass_norm))-lam*np.sum(y*alignment,axis=1)
              -np.sum(delta**2/(e+y),axis=1))
        L=np.exp(logL)
        known.extend(L*((lab>q) if upper else (lab<q)))
        fitted.extend(L*((frame>q) if upper else (frame<q)))
    def moments(values):
        return dict(probability=float(np.mean(values)),standard_error=float(np.std(values,ddof=1)/math.sqrt(len(values))))
    return dict(seed=seed,threshold=q,importance_tilt=lam,exact_massless_known_frame=massless_tail(n,q,upper),
                massive_known_frame=moments(known),massive_fitted_landau_frame=moments(fitted))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    rows={}
    # 450 is the observed count of a 512-particle history at its first hot crossing.
    for n in (256,450,512,1024):
        low=tail(n,False,17263+n);high=tail(n,True,23119+n)
        probs={}
        for key in ('massive_known_frame','massive_fitted_landau_frame'):
            probs[key]={'probability':low[key]['probability']+high[key]['probability'],
                        'standard_error':math.hypot(low[key]['standard_error'],high[key]['standard_error'])}
        exact=low['exact_massless_known_frame']['probability']+high['exact_massless_known_frame']['probability']
        diff=abs(probs['massive_known_frame']['probability']-exact)
        if diff>6*probs['massive_known_frame']['standard_error']+2e-8:
            raise ValueError('Importance sampler and massless angular inversion disagree')
        dispersion={}
        for gas_fraction in (1.,.975):
            threshold=1-.9/gas_fraction
            p=2*massless_tail(n,threshold,True,observable='velocity')['probability']
            dispersion[str(gas_fraction)]=dict(fixed_gas_fraction_of_wall=gas_fraction,
                                               max_absolute_axial_momentum_over_energy=threshold,
                                               massless_failure_probability=p)
        rows[str(n)]=dict(dispersion_snapshot=dispersion,lower=low,upper=high,two_tailed=probs,exact_massless_known_frame_tail=exact,
                         exact_massless_known_frame_sd=4/math.sqrt(45*(3*n+1)))
    out=dict(temperature_over_mass=THETA,replicates_per_count_per_tail=REPLICATES,pressure_tolerance=.05,counts=rows)
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
