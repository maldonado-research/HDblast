"""Certify a residual invariant-energy discrepancy for every axial elastic frame."""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import kve

# Thresholds on gamma_elastic/(theta * Gamma_relative) are dimensionless.
CUTS = np.linspace(.125, 12., 96)
COSINE_GRID = np.linspace(-1., 1., 1025)
SPEED_GRID = np.linspace(-1., 1., 129)


def gamma_three_cdf(z):
    """CDF of the thermal radial variable q in the zero-mass limit."""
    # For z above 100 the missing probability is below 2e-40.
    z = np.minimum(z,100.)
    return -np.expm1(-z)-np.exp(-z)*(z+.5*z*z)


def thermal_bounds(theta, reference_beta, left, right, cuts):
    """Enclose the massive target CDF throughout one boost-velocity interval."""
    mu_left, mu_right = COSINE_GRID[:-1], COSINE_GRID[1:]
    angular_weights = .5*(mu_right-mu_left)+.25*reference_beta*(mu_right**2-mu_left**2)
    epsilon = 1./theta
    density_tv = .5*math.expm1(epsilon)+epsilon/8
    products = np.stack((left*mu_left,left*mu_right,right*mu_left,right*mu_right))
    low_multiplier = 1 + products.min(axis=0)
    high_multiplier = 1 + products.max(axis=0)
    # True normalized energy exceeds q*(1+b*mu) by 0..epsilon. A zero
    # multiplier contributes unit CDF at all positive thresholds.
    upper_argument = np.divide(cuts[:,None],low_multiplier,
                               out=np.full((len(cuts),len(mu_left)),np.inf),where=low_multiplier>0.)
    lower_argument = np.maximum(cuts[:,None]-epsilon,0.)/high_multiplier
    lower = gamma_three_cdf(lower_argument) @ angular_weights - density_tv
    upper = gamma_three_cdf(upper_argument) @ angular_weights + density_tv
    return np.maximum(0.,lower),np.minimum(1.,upper)


def validate_target_bounds(theta, drift, grid):
    """Check independent massive quadrature at interior points of signed intervals."""
    beta = abs(drift)/np.hypot(1.,drift)
    norm = kve(2,1./theta)
    for b in (-.91,-.33,0.,.36,.89):
        index = min(len(grid)-1,np.searchsorted(SPEED_GRID,b,side='right')-1)
        lower,upper = grid[index]
        boost = 1/math.sqrt(1-b*b)
        rapidity = abs(math.atanh(b))
        for index_cut in (1,7,23,47):
            cut = CUTS[index_cut]
            edge = theta*cut
            lab_rapidity = math.acosh(boost*edge)
            lo = 2*math.sinh((lab_rapidity-rapidity)/2)**2/theta
            hi = 2*math.sinh((lab_rapidity+rapidity)/2)**2/theta
            def density(x):
                gamma = 1.+theta*x
                return math.exp(-x)*gamma*math.sqrt(theta*x*(2.+theta*x))/norm
            def partial(x):
                gamma = 1.+theta*x
                momentum = math.sqrt(theta*x*(2.+theta*x))
                mu = float(np.clip((edge-gamma)/(abs(b)*momentum),-1.,1.))
                tilt = math.copysign(beta,b)*momentum/gamma
                return density(x)*(.5*(1+mu)+.25*tilt*(mu*mu-1.))
            exact = quad(density,0.,lo,epsabs=3.e-12)[0]
            if b:
                exact += quad(partial,lo,hi,epsabs=3.e-12)[0]
            if exact < lower[index_cut]-1.e-10 or exact > upper[index_cut]+1.e-10:
                raise ValueError('Independent massive target falls outside its certified bounds')


@lru_cache(maxsize=None)
def certified_target_grid(theta, drift):
    beta = abs(drift)/np.hypot(1.,drift)
    grid = [thermal_bounds(theta,beta,left,right,CUTS)
            for left,right in zip(SPEED_GRID[:-1],SPEED_GRID[1:])]
    validate_target_bounds(theta,drift,grid)
    return grid


def sign_certificate(rest_energy, axial, target):
    lower_bounds, witness = [], []
    for idx,(left,right) in enumerate(zip(SPEED_GRID[:-1],SPEED_GRID[1:])):
        first = rest_energy + left*axial
        last = rest_energy + right*axial
        # For any b in this interval all energies below the cutoff at both
        # endpoints remain below there. The reverse inclusion bounds the CDF above.
        minimum = np.minimum(first,last)
        maximum = np.maximum(first,last)
        observed_lower = np.searchsorted(np.sort(maximum),CUTS,side='right')/len(maximum)
        observed_upper = np.searchsorted(np.sort(minimum),CUTS,side='right')/len(minimum)
        thermal_lower, thermal_upper = target[idx]
        differences = np.stack((observed_lower-thermal_upper,thermal_lower-observed_upper))
        which, cutoff = np.unravel_index(differences.argmax(),differences.shape)
        # Absolute margin against roundoff of the quadrature bounds.
        certified = max(0.,float(differences[which,cutoff])-1.e-10)
        witness.append(dict(normalized_energy_cut=float(CUTS[cutoff]),population_excess=bool(which==0),
                            guaranteed_difference=certified))
        lower_bounds.append(certified)
    return dict(minimum_over_all_axial_frames=min(lower_bounds),witness_by_speed_cell=witness)


def compare(path):
    with np.load(path,allow_pickle=False) as snapshot:
        p = snapshot['p']
        signs = snapshot['sign']
        drift = float(snapshot['drift'])
        theta = float(snapshot['theta'])
    gamma = np.hypot(1.,np.linalg.norm(p,axis=1))
    along = np.sign(drift)*signs*p[:,2]
    boost = np.hypot(1.,drift)
    rest_energy = (boost*gamma-abs(drift)*along)/theta
    axial = (boost*along-abs(drift)*gamma)/theta
    target = certified_target_grid(theta,drift)
    signs_found = [sign_certificate(rest_energy[signs==sign],axial[signs==sign],target)
                   for sign in (-1.,1.)]
    averaged = sum(r['minimum_over_all_axial_frames'] for r in signs_found)/2.
    return dict(snapshot=path.name,minimum_charge_resolved_tv=averaged,
                negative=signs_found[0],positive=signs_found[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshots',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    results = [compare(args.snapshots/f'{model}_half_{pulse}_{power}.npz')
               for model in ('dedicated','shared') for pulse in ('pressure_fit','small_pulse')
               for power in (12,14)]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(scatter_speed_cell_edges=SPEED_GRID.tolist(),
        normalized_energy_thresholds=CUTS.tolist(),angular_cells=len(COSINE_GRID)-1,
        runs=results),sort_keys=True,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
