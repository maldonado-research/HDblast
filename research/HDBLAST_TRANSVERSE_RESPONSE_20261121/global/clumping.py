"""Bound all spatial Maxwell modes and binary exposure under a common local density cap."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
TOOLS = runpy.run_path(str(ROOT / "history/budget.py"))
TUBE = TOOLS["TUBE"]
K = 2.54


class ThermalCost:
    """Monotone lower Riemann sum of squared minimum impulse over quantile mass."""
    def __init__(self, theta, bins):
        self.theta = theta
        grid = np.r_[0., np.geomspace(1.e-3, 36.*theta, bins), TUBE]
        self.p = np.unique(grid)
        self.prob = np.maximum.accumulate(TOOLS["thermal_cdf"](self.p, theta))
        self.cost = np.maximum(0., self.p-TUBE)**2
        self.lower = np.r_[0., np.cumsum(np.diff(self.prob)*self.cost[:-1])]

    def __call__(self, f):
        idx = np.clip(np.searchsorted(self.prob, f, side='right')-1,0,len(self.p)-1)
        return np.sqrt(self.lower[idx]+np.maximum(0.,f-self.prob[idx])*self.cost[idx])

    def independent_integral(self, f):
        upper = float(TOOLS["inverse_cdf"](np.array([f]), self.theta)[0])
        if upper <= TUBE:
            return 0.
        def integrand(p):
            return 2*(p-TUBE)*max(0., f-float(TOOLS["thermal_cdf"](p, self.theta)))
        return math.sqrt(quad(integrand,TUBE,upper,epsabs=1.e-13,epsrel=1.e-10)[0])


def suffix(v):
    return np.cumsum(v[::-1])[::-1]


def measure(data, meta, case, marker, bins):
    spec = meta["runs"][case]
    end = spec["markers"][marker]["time"]
    j = int(np.searchsorted(data["t"], end))
    t,a,d = (data[k][:j+1] for k in ("t","a","D"))
    dt = np.diff(t)
    if np.any(dt<=0) or np.any(np.diff(a)<0) or np.any(np.diff(d)<0) or d[0]!=0:
        raise ValueError("Need monotone expanding birth history starting empty")
    q=meta["charge_to_mass"]
    hazard=spec["channel_fraction"]*meta["action"]*TOOLS["pair_kernel_upper"]()/(6*meta["daughter_mass_L"]**3)
    depth=suffix(hazard*dt*d[1:]/a[:-1]**3)
    births=np.diff(d)/d[-1]
    births[-1]=0.  # The marker lies in the last cell. Exclude those births from all floors.
    thermal=ThermalCost(min(data['theta'][j-1:j+1]),bins)
    # Source bound accumulated up to each right endpoint. It is also an upper
    # bound at every interior point of that cell in conformal field units a²F.
    max_field_at_cap_one=q/K*np.cumsum(dt*d[1:]/a[:-1])
    # C bounds each sign/flavor's pointwise number density divided by its reference
    # spatial mean. Evolving all modes adds sqrt(C) from the source and sqrt(C)
    # from sampling the field at particles. Integrate surviving particles after cut.
    up=q/a[-2] * suffix(dt/a[:-1]*np.sqrt(d[1:]/d[-2])*max_field_at_cap_one)
    seed_exposure=q*a[0]**2/a[-2]*suffix(dt/a[:-1]*np.sqrt(d[1:]/d[-2]))
    out={}
    for tangent_factor, sector in [(math.sqrt(2.),"general_fields"),(1.,"perpendicular_neutral_sector")]:
        upper=tangent_factor*up
        for eps in (0.,.01,.05):
            def check(C):
                floor=suffix(births*np.exp(-C*depth))
                mass=np.maximum(0.,floor-eps)
                required=thermal(mass)
                # A positive gap excludes this C. Endpoints near zero can make
                # a relative gap ill-conditioned; decide only by the signed gap.
                diff=required - C*upper
                h=int(np.argmax(diff))
                return float(diff[h]),h, float(floor[h]),float(required[h])
            lo,hi=1.,2.
            if check(lo)[0]<=0:
                excluded=False
            else:
                excluded=True
                while check(hi)[0]>0: hi*=2
                for _ in range(43):
                    mid=.5*(lo+hi)
                    if check(mid)[0]>0: lo=mid
                    else:hi=mid
            h0=check(1.)
            gap,i,born,cost=check(lo)
            precise=thermal.independent_integral(max(0.,born-eps))
            if precise+2e-8<cost or abs(precise-cost)>max(2e-7, .01*precise):
                raise ArithmeticError("Quantile cost quadrature")
            seed_limits={}
            for cap in (1.,1.25,1.5,2.,4.,10.):
                unhit=suffix(births*np.exp(-cap*depth))
                need=thermal(np.maximum(0.,unhit-eps))
                scale=tangent_factor*math.sqrt(cap)*seed_exposure
                deficit=need-cap*upper
                seed_limits[f"{cap:g}"]=float(max(0.,np.max(deficit/scale)))
            out[f"{sector}_eps_{eps:g}"]=dict(
                initial_rms_field_lower_by_density_cap=seed_limits,
                exclusion_at_unit_density=excluded, excluded_cap_lower=lo if excluded else None,
                last_excluded_cohort_time=float(t[i]),last_excluded_cohort_age=float(end-t[i]),
                mandatory_recent_unhit_share=born, required_impulse_recent_rms_lower=cost,
                required_impulse_quadrature=precise, max_field_impulse_recent_rms_at_cap=lo*float(upper[i]),
                max_global_impulse_rms_at_unit_cap=float(upper[0]), unit_cap_maximum_gap=h0[0],
                unit_cap_witness_recent_rms_required=h0[3], unit_cap_witness_recent_rms_available=float(upper[h0[1]]))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--history",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True);args=p.parse_args()
    meta=json.loads((args.history/"markers.json").read_bytes()); rows={}
    for case,spec in sorted(meta["runs"].items()):
        with np.load(args.history/(case+".npz")) as data:
            rows[case]={}
            for marker in spec["markers"]:
                fine=measure(data,meta,case,marker,24576)
                coarse=measure(data,meta,case,marker,12288)
                for k,v in fine.items():
                    v["coarser_thermal_quadrature_cap"]=coarse[k]["excluded_cap_lower"]
                rows[case][marker]=fine
    args.output.write_text(json.dumps(rows,indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=="__main__":
    main()
