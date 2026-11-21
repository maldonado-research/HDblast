"""Two-sided canonical-momentum fold certificate on the accepted-step replay.

Input rows are native-endian float64, with t,dt,k0,k1,d0,d1,drest,f0,f1,
abs_kdot,kmid,P/(12*kmid**2). d is signed db/dtau. The temporal
interpolant is constant in P, linear in k, and linear in signed d.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def triangle(lo, hi):
    """Integral of max(1-abs(x),0) over [lo,hi]."""
    left = np.maximum(-1, lo)
    right = np.minimum(0, hi)
    negative = np.maximum(0, right-left)*(1+(left+right)/2)
    left = np.maximum(0, lo)
    right = np.minimum(1, hi)
    return negative+np.maximum(0, right-left)*(1-(left+right)/2)


def positive_angular_integral(c, lo, hi, T, divisions):
    """Integrate in |d|. The outer 1/|d| part is exact."""
    a = np.maximum(lo,c+T)
    result = .5*T*np.log1p(np.maximum(0,hi-a)/a)
    lower = np.maximum(lo,np.maximum(0,c-T))
    upper = np.minimum(hi,c+T)
    inner = upper>lower
    if np.any(inner):
        first, last, t, cc = lower[inner],upper[inner],T[inner],c[inner]
        step = (last-first)/divisions
        part = np.zeros_like(t)
        for i in range(divisions):
            small = first+i*step
            large = first+(i+1)*step
            val = triangle((cc-small)/t,(cc+small)/t)*t/(2*large)
            val = np.maximum(val,1-(cc+large)/t)
            part += step*np.maximum(0,val)
        result[inner] += part
    return result


def angular(c, d0, d1, T, divisions):
    small = np.minimum(abs(d0),abs(d1))
    large = np.maximum(abs(d0),abs(d1))
    cross = d0*d1<0
    lo = np.where(cross,0,small)
    area = positive_angular_integral(c,lo,large,T,divisions)
    area[cross] += positive_angular_integral(c[cross],lo[cross],small[cross],T[cross],divisions)
    diff = abs(d1-d0)
    result = area/np.where(diff>0,diff,1)
    fixed = diff==0
    if np.any(fixed):
        cc,dd,t = c[fixed],large[fixed],T[fixed]
        result[fixed] = np.where(dd==0,np.maximum(0,1-cc/t),
                                 triangle((cc-dd)/t,(cc+dd)/t)*t/(2*np.where(dd>0,dd,1)))
    return result


def paired_bound(width, J, m, M, variation, margin, active):
    """Excess in a strip in (z,rho^2) for a contiguous time group."""
    # In a smooth region the normalized strip width is M ell^2. Jumps in Z
    # widen it by 4 V ell and reduce the guaranteed covered z length by V.
    good = active & (m>0) & (M>0)
    reduced_width = np.maximum(0,width-variation/np.where(good,m,1))
    effective_J = np.maximum(0,J-6*variation)
    opt = effective_J/(3*np.where(good,M,1))
    root = np.hypot(reduced_width-2*opt,np.sqrt(12.)*opt)
    denom = reduced_width+4*opt+root
    ell = np.minimum(margin,2*opt*reduced_width/np.where(denom>0,denom,1))
    ell = np.maximum(0,ell)
    excess = m*(reduced_width-2*ell)*(2*ell*effective_J-3*M*ell**2)
    return np.where(good,np.maximum(0,excess),0.)


def measure(arr, mu, action_factor=1., divisions=256, max_level=13):
    _,dt,k0,k1,d0,d1,_,f0,f1,_,kmid,T = arr.T
    if (not np.all(np.isfinite(arr)) or np.any(dt<=0) or np.any(T<=0) or
        np.any(np.minimum(k0,k1)<=0) or np.any(f1<=f0)):
        raise ValueError("Expected positive finite shell, source and duration")
    time_error = abs(arr[1:,0]-(arr[:-1,0]+dt[:-1]))
    if np.any(time_error>1e-12) or np.any(f0[1:]!=f1[:-1]):
        raise ValueError("Input steps do not cover a continuous delivery interval")
    if np.any(k0[1:]!=k1[:-1]):
        raise ValueError("Shell radius is discontinuous across input steps")
    c = ((k1-k0)/dt)/mu
    k0, k1, kmid, d0,d1,T = [v/mu for v in (k0,k1,kmid,d0,d1,T)]
    T = T*action_factor
    weights = f1-f0
    P = 12*kmid**2*T
    total = np.sum(P*dt)/np.sum(weights)
    if not np.allclose(P*dt/total,weights,rtol=1e-10,atol=0.):
        raise ValueError("Demand and cumulative birth fractions disagree")
    tlo = P/(12*np.maximum(k0,k1)**2)
    angular_count = P*dt*angular(abs(c),d0,d1,tlo,divisions)
    dmin = np.minimum(abs(d0),abs(d1))
    dmax = np.maximum(abs(d0),abs(d1))
    kmin = np.minimum(k0,k1)
    kmax = np.maximum(k0,k1)
    active = (d0*d1>0) & (kmin>0) & (dmin>abs(c))
    safe_min = np.where(active,dmin,1.)
    safe_max = np.where(active,dmax,1.)
    signed_term = np.sign(d0)*c*(d1-d0)/dt
    term0 = signed_term*kmin/safe_max**3
    term1 = signed_term*kmax/safe_min**3
    zprime_min = 1-(c/safe_min)**2+np.minimum(term0,term1)
    zprime_max = 1-(c/safe_max)**2+np.maximum(term0,term1)
    width = .5*(safe_max+safe_min)*dt
    J = P/(2*kmax*safe_max)
    margin = kmin-kmax*abs(c)/safe_min
    jumps = np.zeros_like(dt)
    jumps[1:] = abs(k0[1:]*c[1:]/np.where(d0[1:]!=0,d0[1:],1) -
                    k1[:-1]*c[:-1]/np.where(d1[:-1]!=0,d1[:-1],1))
    continuous = np.ones_like(active)
    continuous[1:] = (k0[1:]==k1[:-1]) & (np.sign(d0[1:])==np.sign(d1[:-1]))
    # Time is partitioned into dyadic groups. Each node keeps the stronger of
    # its own strip certificate and the sum of its two disjoint children's.
    levels = []
    for level in range(max_level+1):
        group = 2**level
        starts = np.arange(0,len(dt),group)
        def reduce(op,v):
            return op.reduceat(v,starts)
        lo = reduce(np.minimum,zprime_min)
        hi = reduce(np.maximum,zprime_max)
        jump = jumps.copy()
        jump[starts] = 0.
        connected = continuous.copy()
        connected[starts] = True
        pair = paired_bound(reduce(np.add,width),reduce(np.minimum,J),
                            np.maximum(lo,-hi),np.maximum(abs(lo),abs(hi)),
                            reduce(np.add,jump),reduce(np.minimum,margin),
                            reduce(np.logical_and,active & connected))
        if level==0:
            previous = np.maximum(angular_count,pair)
            pair_only = pair
        else:
            previous = np.maximum(np.add.reduceat(previous,np.arange(0,len(previous),2)),pair)
            pair_only = np.maximum(np.add.reduceat(pair_only,np.arange(0,len(pair_only),2)),pair)
        levels.append(float(np.sum(previous)/total))
    return dict(covered_fraction=float(np.sum(weights)), num_steps=len(weights),
                angular_lower_fraction=float(np.sum(angular_count)/total),
                fold_only_lower_fraction=float(np.sum(pair_only)/total),
                combined_lower_fraction=levels[-1], bounds_by_dyadic_depth=levels)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--mu',type=float,required=True)
    parser.add_argument('--action-factor',type=float,default=1.)
    parser.add_argument('--divisions',type=int,default=256)
    args=parser.parse_args()
    a=np.fromfile(args.input,dtype='=f8').reshape(-1,12)
    print(json.dumps(measure(a,args.mu,args.action_factor,args.divisions),indent=2))
