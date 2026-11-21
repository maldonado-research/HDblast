"""Optimal coarse momentum counts after unrestricted turns at fixed particle energy."""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import linprog
from scipy.special import kve

# Frozen diagnostic cells from the two-polarization experiment.
ENERGY = np.array([0., 1., 2., 3., 5., np.inf])
COSINE = np.linspace(-1., 1., 5)
# A separate invariant-energy measurement in each scattering-center frame.
SPECTRUM = np.array([0., .125, .25, .375, .5, .75, 1., 1.25, 1.5, 2., 2.5,
                     3., 3.5, 4., 5., 6., 8., 10., 15., np.inf])
FRAMES = (0., .25, .5, .75, 1.)


def rest_density(theta, x):
    gamma = 1. + theta*x
    momentum = math.sqrt(theta*x*(2. + theta*x))
    return math.exp(-x)*gamma*momentum/kve(2, 1./theta)


def integral(function, lo, hi):
    if hi <= lo:
        return 0.
    return quad(function, lo, hi, epsabs=8.e-13, epsrel=8.e-12, limit=120)[0]


@lru_cache(maxsize=None)
def target_cells(theta, reference):
    drift = math.tanh(reference)
    rows = []
    for low, high in zip(ENERGY[:-1], ENERGY[1:]):
        def tilt(x):
            gamma = 1 + theta*x
            velocity = math.sqrt(theta*x*(2 + theta*x))/gamma
            return drift*velocity*rest_density(theta,x)
        mass = integral(lambda x: rest_density(theta,x),low,high)
        current = integral(tilt,low,high)
        rows.append(.5*mass*np.diff(COSINE) + .25*current*np.diff(COSINE**2))
    cells = np.array(rows)
    if abs(cells.sum() - 1.) > 2.e-12 or cells.min() <= 0:
        raise ValueError('Invalid moving-stream normalization')
    return cells


def energy_intervals(theta, difference):
    """Image of each (rest energy, cosine) cell in the elastic frame."""
    boost, speed = math.cosh(difference), math.sinh(difference)
    gamma = 1 + theta*ENERGY
    limits = []
    for low, high in zip(gamma[:-1], gamma[1:]):
        for mu0, mu1 in zip(COSINE[:-1], COSINE[1:]):
            # gamma*S + sqrt(gamma**2-1)*T is convex if T<0.
            stationary = 1/math.sqrt(1 - (mu0*speed/boost)**2) if mu0 < 0 else 1.
            closest = np.clip(stationary, low, high)
            def transform(g, mu):
                return boost*g + speed*math.sqrt((g-1)*(g+1))*mu
            minimum = transform(closest,mu0)
            maximum = max(transform(low,mu1),transform(high,mu1)) if np.isfinite(high) else np.inf
            limits.append((minimum,maximum))
    return np.array(limits)


def energy_spectrum(theta, reference, difference):
    """Integrate the thermal count in a frame moving axially relative to the wall."""
    thermal_boost = math.cosh(difference)
    thermal_drift = math.sinh(difference)
    wall_beta = math.tanh(reference)
    scale = theta*thermal_boost
    def cdf(edge):
        if np.isinf(edge):
            return 1.
        if edge == 0:
            return 0.
        if difference == 0:
            return integral(lambda x: rest_density(theta,x),0.,edge)
        threshold = 1. + scale*edge
        speed_edge = math.acosh(threshold)
        x0 = 2*math.sinh((speed_edge-difference)/2)**2/theta
        x1 = 2*math.sinh((speed_edge+difference)/2)**2/theta
        full = integral(lambda x: rest_density(theta,x),0.,min(x0,100.)) if speed_edge > difference else 0.
        def partial(x):
            g = 1. + theta*x
            p = math.sqrt(theta*x*(2 + theta*x))
            cut = max(-1., min(1., (threshold-thermal_boost*g)/(thermal_drift*p)))
            angular = .5*(1+cut) + .25*wall_beta*p/g*(cut*cut-1.)
            return rest_density(theta,x)*angular
        return full + integral(partial,x0,min(x1,100.))
    probabilities = np.diff([cdf(edge) for edge in SPECTRUM])
    if abs(probabilities.sum()-1.) > 2.e-12 or probabilities.min() < 0:
        raise ValueError('Invalid elastic-frame spectrum')
    return probabilities


def tv(first, second):
    return float(.5*np.abs(first-second).sum())


def minimum_tv(frame_energy, limits, target):
    """Hall's deficit, with an independent fractional transport LP as a primal witness."""
    for edge in np.unique(limits[np.isfinite(limits)]):
        if np.any(frame_energy == edge):
            raise ValueError('A sampled energy lands on a bin reachability boundary')
    reachable = ((frame_energy[:,None] >= limits[:,0]) &
                 (frame_energy[:,None] <= limits[:,1]))
    bits = np.left_shift(1,np.arange(len(target)))
    masks, counts = np.unique(reachable @ bits,return_counts=True)
    fractions = counts/len(frame_energy)
    if masks[0] == 0:
        raise ValueError('A particle could not keep its present direction')
    size = 1 << len(target)
    contained = np.zeros(size)
    contained[masks] = fractions
    proposed = np.zeros(size)
    proposed[bits] = target
    # After these subset sums, contained[S] is the supply with no escape from S.
    for bit in bits:
        grouped = contained.reshape(-1,2*bit)
        grouped[:,bit:] += grouped[:,:bit]
        grouped = proposed.reshape(-1,2*bit)
        grouped[:,bit:] += grouped[:,:bit]
    deficit = contained-proposed
    obstruction = int(np.argmax(deficit))
    allowed = (masks[:,None] & bits) != 0
    groups, bins = np.nonzero(allowed)
    edges = len(groups)
    matrix = np.zeros((len(masks)+len(target),edges))
    matrix[groups,np.arange(edges)] = 1.
    matrix[len(masks)+bins,np.arange(edges)] = 1.
    optimum = linprog(-np.ones(edges),A_ub=matrix,b_ub=np.r_[fractions,target],
                      bounds=(0.,None),method='highs')
    if not optimum.success or abs(1. + optimum.fun - deficit[obstruction]) > 2.e-8:
        raise ValueError('Primal allocation and bottleneck certificate disagree')
    allocation = np.zeros(allowed.shape)
    allocation[groups,bins] = optimum.x
    leftover = fractions-allocation.sum(axis=1)
    allocation[np.arange(len(masks)),allowed.argmax(axis=1)] += leftover
    attainable = allocation.sum(axis=0)
    if abs(tv(attainable,target)-deficit[obstruction]) > 2.e-8:
        raise ValueError('Attainable cell counts differ from the sharp bound')
    # Several equally good allocations or cuts can exist. Publish the stable value.
    return float(max(0.,deficit[obstruction]))


@lru_cache(maxsize=None)
def frame_target(theta, reference, fraction):
    difference = (1.-fraction)*reference
    return energy_intervals(theta,difference),energy_spectrum(theta,reference,difference)


def compare(path):
    with np.load(path,allow_pickle=False) as data:
        p = data['p']
        sign = data['sign']
        theta = float(data['theta'])
        drift = float(data['drift'])
        stored_tv = float(data['reference_tv'])
    reference = math.asinh(abs(drift))
    gamma = np.hypot(1.,np.linalg.norm(p,axis=1))
    aligned = np.sign(drift)*sign*p[:,2]
    particle_count = len(sign)
    if np.any(np.abs(sign)!=1) or np.sum(sign)!=0:
        raise ValueError('Unequal charge weights')
    ref_gamma = math.cosh(reference)*gamma-math.sinh(reference)*aligned
    ref_parallel = math.cosh(reference)*aligned-math.sinh(reference)*gamma
    cosine = ref_parallel/np.hypot(ref_parallel,np.linalg.norm(p[:,:2],axis=1))
    observed = np.histogram2d((ref_gamma-1)/theta,cosine,bins=(ENERGY,COSINE))[0]/particle_count
    target = target_cells(theta,reference)
    original_tv = tv(observed,target)
    if abs(original_tv-stored_tv)>1.e-10:
        raise ValueError('Initial 20-bin result differs from source')
    rows = []
    for fraction in FRAMES:
        center = fraction*reference
        frame_gamma = math.cosh(center)*gamma-math.sinh(center)*aligned
        limits,thermal = frame_target(theta,reference,fraction)
        sharp = minimum_tv(frame_gamma,limits,target.ravel())
        scaled_energy = (frame_gamma-1)/(theta*math.cosh(reference-center))
        sampled = np.histogram(scaled_energy,bins=SPECTRUM)[0]/particle_count
        spec_tv = tv(sampled,thermal)
        if sharp>original_tv+2.e-8 or abs(sampled.sum()-1)>1.e-12:
            raise ValueError('A feasible unturned state was lost')
        if fraction==1. and abs(sharp-tv(observed.sum(axis=1),target.sum(axis=1)))>2.e-12:
            raise ValueError('Rest-frame marginal test')
        rows.append(dict(rapidity_fraction=fraction,minimum_20_cell_tv=sharp,
                         invariant_19_bin_tv=spec_tv,observed_energy_bins=sampled.tolist(),
                         thermal_energy_bins=thermal.tolist()))
    return dict(snapshot=path.name,particles=particle_count,original_20_cell_tv=original_tv,
                frame_tests=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshots',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    cases = [f'{model}_half_{pulse}_{power}.npz'
             for model in ('dedicated','shared') for pulse in ('pressure_fit','small_pulse')
             for power in (12,14)]
    results = [compare(args.snapshots/case) for case in cases]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(
        energy_bin_edges_in_theta_times_relative_lorentz_factor=[*SPECTRUM[:-1],'infinity'],
        runs=results),sort_keys=True,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
