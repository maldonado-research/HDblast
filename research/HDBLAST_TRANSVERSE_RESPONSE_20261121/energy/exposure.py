"""Bound integrated transverse Lorentz-force exposure on an expanding shell."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import kve


TUBE = math.sqrt(1.05**2 - 1)
TOLERANCES = (0., .001, .01, .05)


def root_mean_square(theta, probability):
    """Lowest possible RMS impulse over all born tags for a transverse thermal marginal."""
    if probability <= 0:
        return 0.
    z = 1. / theta
    denominator = kve(2, z)
    def cdf(x):
        # x=p_perp/(m*Theta). The survival probability is w^2 K2(w/Theta)/K2(1/Theta).
        r = math.hypot(z, x)
        log_tail = 2 * math.log(r/z) + math.log(kve(2, r)/denominator) + z - r
        return -math.expm1(min(log_tail, 0.))
    b = TUBE / theta
    if cdf(b) >= probability:
        return 0.
    quantile = brentq(lambda x: cdf(x)-probability, b, 80.)
    squared, _ = quad(lambda x: 2 * (x-b)*(probability-cdf(x)), b, quantile,
                      epsabs=1.e-14, epsrel=2.e-11)
    return theta * math.sqrt(squared)


def birth_operator(dt, weight, birth_fraction):
    """Spectral norm enclosure for the positive matrix M_ij = l_i l_j d_min(i,j)."""
    l = np.sqrt(dt) * weight
    def multiply(x):
        y = l * x
        right = np.cumsum(y[::-1])[::-1] - y
        return l * (np.cumsum(birth_fraction * y) + birth_fraction * right)
    x = l.copy()
    for _ in range(100):
        y = multiply(x)
        ratio = np.divide(y, x, where=(x > 0), out=np.zeros_like(x))
        active = x > 0
        low = float(np.min(ratio[active]))
        high = float(np.max(ratio[active]))
        if high-low <= 2.e-11 * high:
            return low*(1.-2.e-12), high*(1.+2.e-12)
        x = y/np.linalg.norm(y)
    raise ValueError('Positive operator norm failed to converge')


def minimum_exposure(x, required):
    """Minimize squared slope of a cumulative field path constrained above each witness."""
    x = np.array(x, dtype=float)
    required = np.array(required, dtype=float)
    if x[0] != 0. or required[0] != 0. or np.any(np.diff(x)<=0.) or np.any(np.diff(required)<-1.e-13):
        raise ValueError('Nonmonotone cumulative field constraint')
    support=[]
    for i in range(len(x)):
        while len(support)>=2:
            a,b=support[-2:]
            if (required[b]-required[a])*(x[i]-x[b]) > (required[i]-required[b])*(x[b]-x[a]):
                break
            support.pop()
        support.append(i)
    active=np.array(support)
    dx=np.diff(x[active])
    slopes=np.diff(required[active])/dx
    energy=float(np.dot(dx,slopes**2))
    if np.any(slopes< -1.e-13) or np.max(required-np.interp(x,x[active],required[active]))>1.e-12:
        raise ValueError('Constructed envelope fails an original constraint')
    # Dual is a nonnegative sum of the original constraints. For energy rho, each
    # feasible path obeys sum lambda_j g_j <= sqrt(rho * ||sum lambda_j 1_[0,x_j]||²).
    multipliers=slopes-np.append(slopes[1:],0.)
    if np.min(multipliers)<-1.e-12:
        raise ValueError('Negative certificate weight')
    certified=0.
    if energy>0.:
        lhs=np.dot(multipliers,required[active[1:]])
        kernel_norm=np.dot(dx,np.cumsum(multipliers[::-1])[::-1]**2)
        certified=float(lhs**2/kernel_norm)
        if abs(energy-certified)>1.e-10*energy:
            raise ValueError('Primal and dual exposure disagree')
    return certified, len(support)-1


def check_envelope():
    # Compare the constructive solver to independent numerical convex optimization.
    from scipy.optimize import minimize
    x=np.linspace(0,1,19)**1.3
    g=np.exp(x)-1.+.22*np.sin(5*x)**2
    g=np.maximum.accumulate(g)
    v,_=minimum_exposure(x,g)
    dx=np.diff(x)
    objective=lambda w: np.dot(w*w,dx)
    linear=np.tril(np.ones((len(dx),len(dx))))*dx
    reference=minimize(objective,np.full_like(dx,10.),jac=lambda w:2*w*dx,method='SLSQP',
                       constraints=[{'type':'ineq','fun':lambda w:linear@w-g[1:],'jac':lambda w:linear}],
                       options={'ftol':1.e-13,'maxiter':1500})
    if not reference.success or abs(reference.fun/v-1)>2.e-8:
        raise ValueError('Convex-program check failed')
    return float(abs(reference.fun/v-1))


def check():
    # With constant birth flux on [0,1] the continuum kernel min(s,t) has norm 4/pi^2.
    exact = 4 / math.pi**2
    estimates = []
    for count in (128, 512, 2048):
        step = np.full(count, 1./count)
        lo, hi = birth_operator(step, np.ones_like(step), np.arange(1, count+1)/count)
        estimates.append(hi)
        if not hi >= exact or hi - lo > 1.e-9:
            raise ValueError('Constant-injection operator failed')
    if not estimates[0] > estimates[1] > estimates[2] or estimates[-1]/exact-1 > .0006:
        raise ValueError('Uniform-injection endpoint bound failed to converge')
    # Unequal time cells, redshift and injection, compared to the dense symmetric formulation.
    rng = np.random.default_rng(240117)
    dt = rng.uniform(.002,.03,96)
    w = np.linspace(.7, 1., len(dt))
    d = np.cumsum(rng.uniform(size=len(dt))); d /= d[-1]
    y = np.sqrt(dt)*w
    dense = y[:,None]*np.minimum.outer(d,d)*y[None,:]
    large = float(np.linalg.eigvalsh(dense)[-1])
    lo, hi = birth_operator(dt,w,d)
    if not (lo <= large <= hi):
        raise ValueError('Dense norm outside positive-vector bracket')
    return {'constant_injection_exact':exact,'constant_injection_upper':estimates,
            'nonuniform_dense_absolute_error':hi-large,
            'obstacle_optimizer_relative_error':check_envelope()}


def assess(data, meta, name, stage, row, budget):
    time = row['marker_time']
    t = data['t']
    j = int(np.searchsorted(t,time))
    if not t[j-1] <= time <= t[j]:
        raise ValueError('Unbracketed observation')
    d = data['D'][:j+1]
    a = data['a'][:j+1]
    if d[j-1] <= 0. or np.any(np.diff(d) < 0.) or np.any(np.diff(a) < 0.):
        raise ValueError('The causal enclosure needs monotone production and expansion')
    # Include the entire terminal cell, enlarge the fainter earlier field weights and
    # overestimate the cumulative incidence of births everywhere in each time cell.
    dt = np.diff(t[:j+1])
    fractions = np.minimum(1.,d[1:]/d[-2])
    weights = np.minimum(1.,a[1:]/a[-2])
    lo, hi = birth_operator(dt,weights,fractions)
    duration = time-t[0]
    q = meta['charge_to_mass']
    endpoints = np.append(t[:j],time)
    bg = np.interp(endpoints,t,data['X'])**2 + np.interp(endpoints,t,data['B'])**2
    reference = float(np.trapezoid(bg,endpoints))
    spec=meta['runs'][name]
    factor=spec['channel_fraction']*meta['action']*budget['pair_kernel_upper']()/(6.*meta['daughter_mass_L']**3)
    _, _, suffix, cutoff_theta, endpoint_j = budget['enclosed_suffix'](data,time,factor)
    if endpoint_j != j or abs(cutoff_theta-row['min_endpoint_theta'])>1.e-12 or abs(suffix[0]-row['guaranteed_unhit_share'])>1.e-12:
        raise ValueError('Birth-cohort witness differs from independently supplied summary')
    # Lookback variable xi=int (a/a_f)^2 dt. Time beyond the marker can be granted
    # with zero physical field. Right-endpoint weights overestimate available impulse.
    xi=np.append(0.,np.cumsum((weights**2*dt)[::-1]))
    out={}
    for tolerance in TOLERANCES:
        rms = root_mean_square(row['min_endpoint_theta'],max(0.,row['guaranteed_unhit_share']-tolerance))
        # Transverse growth of the momentum magnitude uses the Lorentz force bound |E_perp|+|v_parallel||B_perp|.
        # (|E|+|B|)^2 <= 2 (|E|^2+|B|^2).
        minimum = rms*rms/(2*q*q*hi)
        goal=np.maximum(0.,suffix-tolerance)
        lower_impulse=np.maximum(0.,budget['inverse_cdf'](goal,cutoff_theta)-TUBE)
        field_requirement=np.append(0.,lower_impulse[::-1]/(math.sqrt(2.)*q))
        envelope, _=minimum_exposure(xi,field_requirement)
        one=float(np.max(field_requirement[1:]**2/xi[1:]))
        if envelope<one*(1.-1.e-10):
            raise ValueError('Joint witness weaker than the best individual window')
        joint=max(minimum,envelope)
        out[str(tolerance)] = {'all_windows_field_squared_integral_lower':envelope,
                               'all_windows_rms_field_lower':math.sqrt(envelope/duration),
                               'all_windows_to_recorded_axial_energy_lower':envelope/reference,
                               'best_single_window_field_squared_integral_lower':one,
                               'single_spatial_harmonic_mean_energy_ratio_lower':joint/(2*reference),
                               'strongest_integrated_energy_ratio':joint/reference,'rms_transverse_impulse_lower_over_mc':rms,
                               'time_integral_transverse_field_strength_squared_lower':minimum,
                               'rms_transverse_field_strength_lower':math.sqrt(minimum/duration),
                               'integrated_transverse_to_recorded_axial_energy_lower':minimum/reference}
    return {'time':time,'duration':duration,'endpoint_operator_lower_L':lo,
            'endpoint_operator_upper_L':hi,'all_oldest_particle_time_upper_L':float(dt.sum()),
            'axial_field_squared_integral':reference,'force_exposure':out}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    meta=json.loads((args.history/'markers.json').read_bytes())
    budgets=json.loads((args.history/'BUDGET.json').read_bytes())['native_cell_enclosures']
    result={}
    budget=runpy.run_path(str(Path(__file__).resolve().parents[1]/'history/budget.py'))
    for name, stages in budgets.items():
        with np.load(args.history/(name+'.npz'),allow_pickle=False) as source:
            result[name]={stage:assess(source,meta,name,stage,row,budget) for stage,row in stages.items()}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps({'controls':check(),'runs':result},indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
