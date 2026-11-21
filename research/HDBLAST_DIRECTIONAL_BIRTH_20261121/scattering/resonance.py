"""Evaluate scalar and vector elastic poles using the thermal momentum-drag kernel."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import kv, roots_laguerre, roots_legendre

VECTOR_MASS = 8300.0       # mediator mass / charged daughter mass
VECTOR_WIDTH = 0.30        # on-shell width / mediator mass


def survivor(z):
    return (z**4*kv(2,z) + 3*z**3*kv(3,z))/24


def distribution(cells=4096):
    edges=np.geomspace(.001,45.,cells+1)
    tail=survivor(edges)
    probability=tail[:-1]-tail[1:]
    if np.min(probability)<0:
        raise ArithmeticError("Bessel tail not monotone")
    return edges,probability,float(1-tail[0]+tail[-1])


def response(r, width):
    """Adaptive, independent reference integration; r=M_res/T."""
    def kernel(z):
        weight=(z**3*kv(2,z)+z**4*kv(1,z))/24
        return weight/(1+((r*r/(z*z)-1)/width)**2)
    return sum(quad(kernel,a,b,epsabs=1e-12,epsrel=1e-12)[0] for a,b in [(0,r),(r,max(r+64,90))])


def bounds(rmin,rmax,width,grid):
    """Enclose the thermal integral for every pole position r in one interval."""
    z,w,missing=grid
    low=(rmin/z[1:])**2
    high=(rmax/z[:-1])**2
    far=np.maximum(np.abs(low-1),np.abs(high-1))
    near=np.maximum(0.,np.maximum(low-1,1.-high))
    lower=float(w @ (1/(1+(far/width)**2)))
    upper=float(w @ (1/(1+(near/width)**2)))+missing
    return lower,upper


def global_scalar(grid, width, bins):
    # Cover [0,infinity] including the two infinite-dynamic-range tails.
    edges=np.concatenate(([0.],np.geomspace(.01,100.,bins+1),[np.inf]))
    intervals=[bounds(lo,hi,width,grid)[1] for lo,hi in zip(edges,edges[1:])]
    best=minimize_scalar(lambda logr:-response(math.exp(logr),width),bounds=(-2,3),method="bounded")
    return {"upper_response":max(intervals),"maximizing_ratio_estimate":math.exp(best.x),
            "optimized_response_estimate":-best.fun}


def direct_momentum_ratio(r,width,order=96):
    # Independent massless two-stream CM transport average in x=p1/T, y=p2/T, cos(angle).
    x,w=roots_laguerre(order)
    c,v=roots_legendre(order)
    p=x[:,None];q=x[None,:]
    total=0.
    for angle,angular_weight in zip(c,v):
        k2=2*(1-angle)*p*q
        transfer=p/q+q/p-2*angle
        probability=1/(1+((r*r/k2-1)/width)**2)
        total+=angular_weight/24*float(w @ (p*p*q*q*transfer*probability) @ w)
    return total


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--history",type=Path,required=True)
    p.add_argument("--sweep",type=Path,required=True,help="Prior kinetics/SWEEP.json")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    data=json.loads(args.history.read_text())['runs']
    sweep=json.loads(args.sweep.read_text())['runs']['mass_4e12']
    grid=distribution()
    scalar=global_scalar(grid,.5,2048)
    finer=global_scalar(distribution(8192),.5,4096)
    scalar["finer_upper_response"]=finer["upper_response"]
    strict_cap=sweep['runs']['native']['continuum_pointwise_S0_upper']
    floor=data['mass_4e12']['maximum_action_floor']
    action=5*floor
    scalar["max_rate_over_required_at_selected_action"]=5*.25*finer["upper_response"]
    scalar["max_rate_over_required_at_sharp_sweep_upper"]=strict_cap/floor*.25*finer["upper_response"]
    checks={}
    for name,history in data.items():
        bins=history["bins"]
        lower,upper=[],[]
        sampled=[]
        for r in bins:
            a,b=bounds(VECTOR_MASS/r['max_theta'],VECTOR_MASS/r['min_theta'],VECTOR_WIDTH,grid)
            lower.append(.75*action*a/r['floor'])
            upper.append(.75*action*b/r['min_floor'])
            sampled.append(.75*action*response(VECTOR_MASS/r['theta_at_floor'],VECTOR_WIDTH)/r['floor'])
        checks[name]=dict(minimum_ratio_lower=min(lower),maximum_ratio_upper=max(upper),
                           lowest_sampled_ratio=min(sampled),highest_sampled_ratio=max(sampled),
                           min_theta=min(r['min_theta'] for r in bins),max_theta=max(r['max_theta'] for r in bins),
                           max_abs_drift=max(r['max_abs_drift'] for r in bins),duration=history['duration'],
                           max_thermal_fugacity_upper=max(math.pi*history['desired_drag']*action/
                               (2*r['min_floor']*history['mu']*r['min_theta']) for r in bins))
    moments=[quad(lambda z:z**i*(z**3*kv(2,z)+z**4*kv(1,z))/24,0,100,epsabs=4e-12)[0]
             for i in (0,2)]
    direct={}
    for ratio,width in [(1.50,.30),(3.80,.30),(3.16,.50)]:
        a=direct_momentum_ratio(ratio,width)
        b=response(ratio,width)
        direct[str(ratio)]={"original_three_dimensional":a,"reduced_bessel":b,"absolute_difference":abs(a-b)}
    if max(v["absolute_difference"] for v in direct.values())>2e-5:
        raise ArithmeticError("Direct thermal momentum integral and Bessel kernel disagree")
    if not (abs(moments[0]-1)<1e-10 and abs(moments[1]-20)<1e-9 and
            all(r['minimum_ratio_lower']>1 for r in checks.values()) and
            scalar["max_rate_over_required_at_sharp_sweep_upper"]<1):
        raise ArithmeticError("Thermal controls or interval test failed")
    args.output.write_text(json.dumps(dict(kernel_moments=moments,momentum_controls=direct,scalar=scalar,
        trial_action=action,prior_sharp_shell_ceiling=strict_cap,vector=dict(
        mediator_over_daughter_mass=VECTOR_MASS,width_over_pole_mass=VECTOR_WIDTH,
        effective_alpha=3*VECTOR_WIDTH,runs=checks)),sort_keys=True,indent=2)+'\n')


if __name__=="__main__":
    main()
