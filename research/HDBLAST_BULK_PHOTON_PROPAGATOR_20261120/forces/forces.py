"""Static potential and force of conserved wall test charges on three frozen tuned leaves."""
import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import gamma, roots_legendre

HERE = Path(__file__).resolve().parent
TOP = 16.
EDGES = (0., .015625, .03125, .0625, .125, .25, .5, 1., 1.5, 2., 3., 4., 6., 8., 12., TOP)
DISTANCES = (2., 3., 4., 6., 8., 16., 32., 64., 128., 256.)


def sibling(name):
    spec = importlib.util.spec_from_file_location(name, HERE.parent/name/(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def quadrature(phi, order, spec):
    grid, weights = roots_legendre(order)
    nodes = np.concatenate([(a+b+(b-a)*grid)/2 for a,b in zip(EDGES[:-1],EDGES[1:])])
    coefficients = np.concatenate([(b-a)*weights/2 for a,b in zip(EDGES[:-1],EDGES[1:])])
    # The low end is open, as all Gauss nodes lie inside their panels. The
    # spectral measure on strictly positive mass excludes the discrete zero pole.
    rho = spec.spectral_density(phi, nodes, outer=7.0, rtol=9e-11)
    return nodes, coefficients*rho*spec.norm(phi)


def excess(nodes, coefficients, radius, force):
    x = nodes*radius
    kernel = np.exp(-x)
    if force:
        kernel *= 1+x
    return float(coefficients @ kernel)


def tail(radius, force, moment):
    # The total moment N int rho(m)/m^2 dm is J/N. The maximum of the
    # potential (force) kernel times m^2 is attained below 2/R ((1+sqrt3)/R).
    x = TOP*radius
    assert x >= 1+math.sqrt(3)
    return max(sys.float_info.min, moment*TOP**2*(1+x if force else 1)*math.exp(-x))


def calculate():
    spec = sibling('retarded')
    eucl = sibling('euclidean')
    arch = json.loads((HERE.parent/'euclidean'/'EUCLIDEAN.json').read_text())
    rows = {}
    exponent = 2*spec.S-1
    for phi in (0., .5, -.5):
        n = eucl.Leaf(phi).norm
        moment = arch['leaves'][str(phi)]['correction_over_norm']
        coarse = quadrature(phi, 8, spec)
        fine = quadrature(phi, 16, spec)
        # The spectral routine is most delicate near the zero pole. An independent
        # outgoing-wave radius and integration direction verify the low-mass scale.
        m = np.array([.001,.005,.05,.25,1.,4.,12.])
        reverse = spec.continuum(phi,m,outer=9.,rtol=4e-11,forward=True)
        primary = spec.spectral_density(phi,m,outer=7.,rtol=9e-11)
        independent = float(max(abs(primary/reverse-1)))
        if independent > 3e-5:
            raise ArithmeticError(f'Continuum boundary check diverged at phi={phi}: {independent}')
        large = n*spec.threshold(phi)*gamma(exponent)
        answer = {}
        for force, key in [(False,'potential'),(True,'force')]:
            asymptotic = large*(1+exponent if force else 1)
            values = {}
            for r in DISTANCES:
                eta = excess(*fine,r,force)
                simple = excess(*coarse,r,force)
                remainder = tail(r,force,moment)
                values[f'{r:g}'] = dict(excess_over_zero_mode=eta,upper_bound_on_high_mass_tail=remainder,
                    integration_order_difference=abs(simple-eta),long_distance_ratio=eta*r**exponent/asymptotic)
                if eta<=0 or abs(simple-eta)>5e-6*eta or remainder>1e-6*eta:
                    raise ArithmeticError(f'Finite-grid refinement or tail control failed at {phi}, {r}, {key}')
            def crossing(data, include_tail):
                def residual(r):
                    return excess(*data,r,force)+(tail(r,force,moment) if include_tail else 0)-.01
                return float(brentq(residual,.9,25.,xtol=2e-12))
            rlo = crossing(fine,False)
            rhi = crossing(fine,True)
            coarse_shift = abs(crossing(coarse,False)-rlo)
            if coarse_shift>1e-5 or rhi-rlo>1e-4:
                raise ArithmeticError('One-percent radius did not resolve')
            if abs(values['256']['long_distance_ratio']-1)>.04:
                raise ArithmeticError('Power-law asymptote was not approached')
            answer[key] = dict(samples=values,leading_power=exponent,asymptotic_coefficient=asymptotic,
                    radius_at_one_percent_truncated=rlo,radius_at_one_percent_with_maximum_tail=rhi,
                    radius_quadrature_refinement=coarse_shift)
        rows[f'{phi:g}'] = dict(one_copy_norm=n,finite_spectral_moment=moment,
            check_other_direction_max_relative=independent,
            short_distance_potential_over_coulomb_coefficient=2*n/(math.pi*(1-phi**2)**2),**answer)
    return dict(radius_unit='wall width L0, with brane warp equal to one',
                highest_integrated_mass=TOP,panels=list(EDGES),gauss_orders=[8,16],leaves=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(),indent=2,allow_nan=False)+'\n')
