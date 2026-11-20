"""Positive finite photon sectors with two identical low-frequency response coefficients."""
import argparse
from fractions import Fraction
import importlib.util
import json
from math import pi
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
KREF, KMAX = 48, 192
CUTS = (48, 64, 96, 192)


def determinant(columns):
    a, b, c = columns
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - b[0]*(a[1]*c[2]-a[2]*c[1])
            + c[0]*(a[1]*b[2]-a[2]*b[1]))


def mixtures():
    """Signed cofactors yield two exactly normalized convex cutoff mixtures."""
    columns = []
    for k in CUTS:
        k2 = sum(Fraction(j*j-1, j) for j in range(KREF+1, k+1))
        k4 = sum(Fraction(j*j-1, j**3) for j in range(KREF+1, k+1))
        columns.append((1, k2, k4))
    null = [(-1)**i*determinant(columns[:i]+columns[i+1:]) for i in range(4)]
    for row in range(3):
        if sum(cols[row]*n for cols, n in zip(columns, null)) != 0:
            raise ArithmeticError('Rational null-vector construction failed')
    # A places weight at 64/192, B at 48/96.
    if null[1] < 0:
        null = [-v for v in null]
    out = []
    for sign in (1, -1):
        p = [max(sign*v, 0) for v in null]
        out.append([v/sum(p) for v in p])
    for row in range(3):
        a, b = [sum(p[i]*columns[i][row] for i in range(4)) for p in out]
        if a != b:
            raise ArithmeticError('Low-frequency moments differ')
    return np.array(out, dtype=float)


def moment_limits(labels, moment_targets, born_shells, count):
    """Obtain a feasible extremum and numerical dual enclosure for the Born energy."""
    cuts = np.arange(KREF,KMAX+1)
    cut = (labels <= cuts[:,None])
    high = labels > KREF
    weights = labels*(labels**2-1)
    objective = cut @ born_shells
    matrix = np.array([np.ones(len(cuts)),
                       cut @ (high*weights/labels**2),
                       cut @ (high*weights/labels**4)])[:count+1]
    scales = matrix.max(axis=1)
    matrix = matrix/scales[:,None]
    target = np.r_[1.,moment_targets][:count+1]/scales
    bounds = []
    for direction in (1, -1):
        c = direction*objective/max(objective)
        opt = linprog(c, A_eq=matrix, b_eq=target, bounds=(0,None), method='highs')
        if not opt.success:
            raise ArithmeticError(opt.message)
        active = np.nonzero(opt.x>1e-12)[0]
        if len(active)>count+1 or max(abs(matrix@opt.x-target))>1e-10 or min(opt.x)<-1e-10:
            raise ArithmeticError('Moment-constrained primal solution is infeasible')
        dual = opt.eqlin.marginals
        # Because the probabilities sum to one, this shift in the intercept
        # bounds the objective even in the presence of a small reduced-cost error.
        deficit = max(0., float(np.max(dual@matrix-c)))
        certificate = dual@target-deficit
        gap = (c@opt.x-certificate)*max(objective)
        if abs(gap)>2e-7:
            raise ArithmeticError('Primal and dual bounds disagree')
        bounds.append(dict(total_D_Born=float(objective@opt.x),
                           cuts=cuts[active].tolist(),probabilities=opt.x[active].tolist(),
                           signed_numerical_duality_gap_D=float(gap)))
    return dict(minimum=bounds[0],maximum=bounds[1])


def dispersive_bounds(labels, weights, charges2, targets):
    """Bounds for the high-band real susceptibility below its first possible pole."""
    U0,U1=targets
    left,right=1./labels[-1]**2,1./labels[0]**2
    mean=U1/U0
    w_right=(mean-left)/(right-left)
    measure=charges2*(weights/labels**2)
    rows=[]
    for omega in (20., 60., 90.):
        z=omega**2/4.
        def kernel(t):
            return 1/(1-z*t)
        lower=kernel(mean)
        upper=(1-w_right)*kernel(left)+w_right*kernel(right)
        actual=measure@kernel(1./labels**2)/U0
        if np.any(actual<lower-1e-12) or np.any(actual>upper+1e-12):
            raise ArithmeticError('Susceptibility outside positive moment bounds')
        rows.append(dict(conformal_frequency=omega,lower=lower,upper=upper,
                         normalized_A=float(actual[0]),normalized_B=float(actual[1])))
    return rows


def spectrum(roll, gauge, labels, charges, rtol):
    """Independent covariances n, Re(m), Im(m) from a vacuum Gaussian."""
    shape = (3, len(charges), len(labels))
    y0 = np.zeros(shape)
    delta = .1
    def evolve(u, y):
        n, real, imag = y.reshape(shape)
        s = charges*(delta*gauge.raw_bump(u)/gauge.NORM)
        freq = 2*roll.eta(u,1)*labels
        return np.stack((2*s*real, freq*imag+s*(2*n+1), -freq*real)).ravel()
    ans = solve_ivp(evolve,(0.,1.),y0.ravel(),method='DOP853',rtol=rtol,
                    atol=rtol*1e-5,max_step=.012,t_eval=[1.])
    if not ans.success:
        raise ArithmeticError(ans.message)
    n, real, imag = ans.y[:,-1].reshape(shape)
    err = max(np.max(abs(real**2+imag**2-n*(n+1))), -np.min(n))
    if err>2e-10:
        raise ArithmeticError('Gaussian purity or positive occupation failed')
    return n, float(err)


def calculate(prior):
    spec=importlib.util.spec_from_file_location('frozen_gauge',prior/'gauge/gauge.py')
    gauge=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gauge)
    hist=json.loads((prior/'gauge/GAUGE.json').read_bytes())
    branches, old = gauge.load(prior)
    labels=np.arange(2,KMAX+1)
    high=labels>KREF
    weights=labels*(labels**2-1)
    probs=mixtures()
    charges2=probs@(labels<=np.array(CUTS)[:,None])
    targets=np.column_stack([charges2@(high*weights/labels**p) for p in (2,4)])
    if max(abs(targets[0]-targets[1]))>1e-9 or np.any(np.diff(charges2,axis=1)>1e-14):
        raise ArithmeticError('Rounded profiles differ from exact moment construction')
    third=charges2@(high*weights/labels**6)
    if abs(third[0]-third[1])<1e-12:
        raise ArithmeticError('The third moment did not separate the profiles')
    out=dict(common_low_band=KREF,highest_shell=KMAX,cuts=list(CUTS),
             mixture_A=probs[0].tolist(),mixture_B=probs[1].tolist(),
             high_band_matched_coefficients=targets[0].tolist(),high_band_third_moments=third.tolist(),
             calibration='For L>48: sum L(L^2-1) h_L^2/L^(2j), j=1,2. Both profiles have h_L=1 for L<=48.',
             branches={})
    out['high_band_dispersive_bounds']=dispersive_bounds(labels[high],weights[high],charges2[:,high],targets[0])
    for key,b in branches.items():
        roll=gauge.Roll(b,old[key],mesh=4001)
        test=np.concatenate([np.sqrt(charges2),np.ones((1,len(labels)))])
        n,err=spectrum(roll,gauge,labels,test,2e-12)
        refined,err2=spectrum(roll,gauge,labels,test,3e-13)
        old_case=hist['branches'][key]['primary_spectrum']
        born=np.asarray(old_case['occupation_born'])[:len(labels)]*weights/pi**2
        uniform_ref=np.asarray(old_case['occupation'])[:len(labels)]
        bins=refined*weights/pi**2
        d=bins.sum(axis=1)
        d_high=bins[:,high].sum(axis=1)
        unity_error=float(max(abs(refined[2]-uniform_ref)))
        repeat_error=float(np.max(abs((n-refined)@weights))/pi**2)
        if unity_error>2e-10 or repeat_error>3e-7:
            raise ArithmeticError('Independent reference or resolution check failed')
        results={}
        for idx,name in enumerate(('A','B')):
            results[name]=dict(D=float(d[idx]),high_band_D=float(d_high[idx]),
                               D_Born=float(charges2[idx]@born),
                               high_band_D_Born=float(charges2[idx,high]@born[high]),
                               photon_over_original_parent=float(d[idx]*roll.points['first_minimum']['photon_to_parent_per_D']))
        bounds={str(i):moment_limits(labels,targets[0],born,i) for i in (1,2)}
        for limit in bounds.values():
            for z in results.values():
                if not limit['minimum']['total_D_Born']-1e-7 <= z['D_Born'] <= limit['maximum']['total_D_Born']+1e-7:
                    raise ArithmeticError('Explicit example violates the moment bounds')
        out['branches'][key]=dict(profiles=results,ratio_total_A_over_B=float(d[0]/d[1]),
                                  ratio_high_A_over_B=float(d_high[0]/d_high[1]),born_moment_ranges=bounds,
                                  maximum_purity_error=max(err,err2),
                                  unweighted_historical_occupation_error=unity_error,
                                  resolution_absolute_D_error=repeat_error)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--prior',type=Path,default=ROOT/'prior')
    opts=p.parse_args()
    opts.output.write_text(json.dumps(calculate(opts.prior),indent=2,allow_nan=False)+'\n')
