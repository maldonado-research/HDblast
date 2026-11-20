"""Continuum excitation by a small prescribed Gaussian displacement of a tuned wall."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
from scipy.special import gamma, roots_legendre

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def table(phi, outgoing, box, selected=False, n_segments=240):
    norm = outgoing.norm(phi)
    wall_weight = (1-phi**2)**2
    prefactor = wall_weight**2/(4*norm)
    yb = -math.atanh(phi)
    pb = outgoing.primitive(yb)

    if selected:
        def weight(y):
            return math.exp(-2*(outgoing.primitive(y)-pb)/3)
        norm = quad(weight,yb,45.,epsabs=2e-13,epsrel=2e-13)[0]
        wall_weight = 1.
        prefactor = 1/(4*norm)

    def tail_integrand(y):
        # The integrated sech^4 tail equals (1-tanh y)^2(2+tanh y)/3.
        exp2 = math.exp(-2*y)
        t = (1-exp2)/(1+exp2)
        tail = (2*exp2/(1+exp2))**2*(2+t)/3
        w = (4*exp2/(1+exp2)**2)**2
        warp = -(outgoing.primitive(y)-pb)/3
        return tail**2/(w*math.exp(2*warp))

    inverse_moment = None
    if selected:
        # Compact trial flux j(y)=1-(y-yb) on [yb,yb+1] in the dual
        # variational principle gives an upper bound on the full Q=1 propagator.
        def dual_density(y):
            w = weight(y)
            return 1/w + (1-(y-yb))**2/w**2
        tail_moment = quad(dual_density,yb,yb+1,epsabs=2e-13)[0]
    else:
        inverse_moment = quad(tail_integrand,yb,40.,epsabs=2e-13,epsrel=2e-13)[0]/norm**2
        tail_moment = inverse_moment
    # Interpolate on log mass. The outgoing solver is also evaluated on
    # the nested, finer grid used for every reported central value.
    mass = np.geomspace(.001,20.,2*n_segments+1)
    density = np.empty_like(mass)
    for inds in np.array_split(np.arange(len(mass)),6):
        density[inds] = outgoing.spectral_density(phi,mass[inds],weight="selected" if selected else "wall")
    power = 1. if selected else 2*outgoing.S-2
    coefficient = 1/(4*outgoing.K**3*norm**2) if selected else outgoing.threshold(phi)

    def interpolant(stride):
        curve = PchipInterpolator(np.log(mass[::stride]),np.log(density[::stride]))
        def evaluate(m):
            m = np.asarray(m)
            vals = coefficient*m**power
            in_grid = m >= mass[0]
            if np.any(m > mass[-1]):
                raise ValueError('Spectral table ended before pulse cutoff')
            vals[in_grid] = np.exp(curve(np.log(m[in_grid])))
            return vals
        return evaluate

    coarse, fine = interpolant(2), interpolant(1)

    def integrate(k, width, spectral, points):
        # v^2=T(omega-k). The Gaussian weights the same v interval at each T.
        v_max = math.sqrt(8.)
        nodes, weights = roots_legendre(points)
        v = (nodes+1)*v_max/2
        masses = v*np.sqrt(2*k/width+(v/width)**2)
        jac = 2*(k+v**2/width)/np.sqrt(2*k*width+v**2)
        argument = v*v
        plus = argument+2*k*width
        kernel = 2*math.pi*(argument**2*np.exp(-argument**2)+plus**2*np.exp(-plus**2))
        return prefactor*v_max/2*float(np.dot(weights,jac*spectral(masses)*kernel))

    slow_constant = prefactor*(math.pi*coefficient/2)*gamma((power+5)/4)
    rows = {}
    independent = box.table(phi,k_values=(1.,),selected=True) if selected else box.table(phi)
    for k, widths in [(0.25,(.5,1.,2.)), (1.,(.5,1.,2.,4.,16.,64.,256.,1024.,4096.)),
                      (2.,(.5,1.,2.))]:
        for width in widths:
            result = integrate(k,width,fine,256)
            mesh = integrate(k,width,coarse,256)
            gaussian = integrate(k,width,fine,384)
            cut = math.sqrt((8/width)*(2*k+8/width))
            # m² X(omega-k) decreases beyond the chosen cutoff; the
            # positive-frequency sum kernel there is bounded by 2X(omega-k).
            kernel_tail = 4*math.pi*64*math.exp(-64)
            tail = prefactor*tail_moment*(cut**2+int(selected))*kernel_tail
            key = f'{k:g}/{width:g}'
            row = dict(k=k,width=width,energy_coefficient=result,
                       slow_leading=slow_constant*(2*k/width)**((power+1)/2),
                       spectral_mesh_relative=abs(mesh/result-1),
                       quadrature_relative=abs(gaussian/result-1),
                       omitted_high_mass_bound=tail)
            if k == 1. and width >= 16.:
                # The same positive spectrum determines the static potential.
                # Use x=mr and integrate to 50, well inside the tabulated band.
                r = math.sqrt(width/(2*k))
                nodes, weights = roots_legendre(160)
                x = 25*(nodes+1)
                potential = norm*25/r*float(weights@(fine(x/r)*np.exp(-x)))
                predicted_ratio = (wall_weight/norm)**2*math.pi*gamma((power+5)/4)/(8*gamma(power+1))
                row['static_separation'] = r
                row['static_potential_excess'] = potential
                row['pulse_to_potential'] = result/potential
                row['pulse_to_potential_asymptote'] = predicted_ratio
                # Both (m²+1) exp(-mr) and m² exp(-mr) decrease above x=50.
                row['static_tail_bound'] = norm*tail_moment*((50/r)**2+int(selected))*math.exp(-50)
            if key in independent['estimates']:
                reference = independent['estimates'][key]['energy']
                row['independent_box_relative'] = abs(result/reference-1)
                if row['independent_box_relative'] > 2e-4:
                    raise ArithmeticError((key,row,reference))
            rows[key] = row
    return dict(half_norm=norm,bulk_weight_at_wall=wall_weight,
                inverse_mass_squared_moment=inverse_moment,
                tail_moment_bound=tail_moment,
                spectral_grid_min=float(mass[0]),spectral_grid_max=float(mass[-1]),
                smallest_mass_threshold_ratio=float(density[0]/(coefficient*mass[0]**power)),
                continuum_threshold_power=power,fixed_k_slow_duration_power=-(power+1)/2,
                cases=rows,independent_box=independent)


def calculate():
    outgoing = load('wave_response',ROOT/'retarded/retarded.py')
    box = load('finite_box',ROOT/'pulse/box.py')
    return dict(normalization='continuum energy / incident energy per squared displacement in wall-width units, leading small-amplitude coefficient',
                branches={f'{phi:g}':table(phi,outgoing,box) for phi in (0.,.5,-.5)},
                selected_branches={f'{phi:g}':table(phi,outgoing,box,selected=True) for phi in (0.,.5,-.5)})


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(),indent=2,allow_nan=False)+'\n')
