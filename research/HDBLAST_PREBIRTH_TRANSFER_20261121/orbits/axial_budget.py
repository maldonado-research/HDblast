"""Amplitude-independent energy floor for a linearly polarized wave and finite daughter-induced modes."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import kv

def envelope(history, marker, charge, modes):
    # Interpolation and inequalities are for the recorded piecewise-linear histories.
    time = history["t"]
    end = marker["time"]
    t = np.r_[time[time<end],end]
    z = {key:np.interp(t,time,history[key]) for key in ("a","X","D","R")}
    a,x,d,r=(z[key] for key in ("a","X","D","R"))
    if (np.any(x>=0.) or np.any(np.diff(d)<-1e-18) or np.any(np.diff(a)<0.) or
            d[0]!=0. or not np.isclose(d[-1],marker["D"],rtol=2e-6)):
        raise ValueError("This positive endpoint bound expects a negative electric history and increasing D and a")
    K=2.54
    dt=np.diff(t)
    fraction=d/d[-1]
    # Lower integral of |axial E| weighting cohorts actually present at each time.
    baseline=charge/a[-1]*np.sum(dt*a[:-1]*np.minimum(-x[:-1],-x[1:])*fraction[:-1])
    birth=np.sqrt(max(0.,np.max(r)**2-1.))
    # Each mode sees |source| <= 2 * Q/K * D. Passive damping and propagation are contractions
    # in its Euclidean norm. Bound the whole cell by its right accumulated integral.
    bound=np.cumsum(dt*d[1:]/a[:-1])
    error=(2*np.sqrt(2)*charge**2/(K*a[-1]) *
           np.sum(dt*bound*fraction[1:]/a[:-1]))
    central=charge/a[-1]*trapezoid(a*(-x)*fraction,t)
    dc=cumulative_trapezoid(d/a,t,initial=0.)
    central_error=2*np.sqrt(2)*charge**2/(K*a[-1])*trapezoid(dc*fraction/a,t)
    # For the equal-quantile orbit grids, the actual number-weighted rest deviates from
    # D by <= D(T)/(2Nc) whenever quantile birth times interpolate the same monotone D(t).
    quantile={}
    for cohorts in (18,36,48,72):
        amount=.5/cohorts
        low=np.maximum(0.,fraction[:-1]-amount)
        up=np.minimum(1.,fraction[1:]+amount)
        base_q=charge/a[-1]*np.sum(dt*a[:-1]*np.minimum(-x[:-1],-x[1:])*low)
        u_q=np.cumsum(dt*d[-1]*up/a[:-1])
        e_q=2*np.sqrt(2)*charge**2/(K*a[-1])*np.sum(dt*u_q*up/a[:-1])
        quantile[str(cohorts)]={"axial_push_lower":float(base_q),
                                "induced_per_mode_upper":float(e_q)}
    theta=marker["theta"]
    gamma_th=3*theta+kv(1,1/theta)/kv(2,1/theta)
    def floor(p,mu,n):
        return float(np.hypot(1.,max(0.,p-birth-n*mu)) / gamma_th)
    return dict(thermal_mean_gamma=float(gamma_th),reference_signed_axial_push_lower=float(baseline),
        birth_momentum_upper=float(birth),induced_axial_impulse_per_harmonic_upper=float(error),
        midpoint_check=dict(axial_push=float(central),induced_per_harmonic=float(central_error)),
        energy_over_thermal_lower={str(n):floor(baseline,error,n) for n in modes},
        equal_quantile_energy_over_thermal_lower={str(c):{str(n):floor(q["axial_push_lower"],
            q["induced_per_mode_upper"],n) for n in modes} for c,q in quantile.items()})

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--history",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    args=p.parse_args(); meta=json.loads((args.history/"markers.json").read_text())
    rows={}
    for name, spec in meta["runs"].items():
        with np.load(args.history/(name+".npz")) as h:
            rows[name]={key:envelope(h,m,meta["charge_to_mass"],(0,1,3,7))
                        for key,m in spec["markers"].items()}
    args.output.write_text(json.dumps(rows,indent=2,sort_keys=True,allow_nan=False)+'\n')
if __name__=="__main__":
    main()
