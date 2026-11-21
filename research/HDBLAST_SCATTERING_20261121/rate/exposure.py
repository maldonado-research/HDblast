"""Collision counts from the first-encounter source and the arbitrary-target rate cap."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import io
import json
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np
from scipy.optimize import brentq
from scipy.stats import qmc


POLE = 8300.
WIDTH = .30
HERE = Path(__file__).resolve().parent


def kernel_upper():
    """Piecewise bound on m² times the total cross section and Møller speed."""
    pole_squared = POLE**2
    cut = pole_squared/2
    low = 24*np.pi*WIDTH**2*(cut+2)**2/(cut*(pole_squared-cut)**2)
    d = WIDTH**2*(1-8/pole_squared)
    high = (24*np.pi/pole_squared * WIDTH**2*(1+4/pole_squared)**2
            * (1+np.sqrt(1+d))/(2*d))
    return max(low, high)


def pair_tags(h, marker, power, scramble, directions):
    # Two independent three-dimensional phase coordinates in one product quadrature.
    coords = qmc.Sobol(6, scramble=True, seed=scramble).random_base2(power)
    columns = []
    for sign, part in [(1,coords[:,:3]), (-1,coords[:,3:])]:
        birth = np.interp(part[:,0]*marker['D'], h['D'], h['t'])
        def at(key):
            return np.interp(birth, h['t'], h[key])
        cosine = directions(at('c'), at('d'), at('target'), part[:,1])
        momentum = sign*at('a')*np.sqrt(np.maximum(0., at('R')**2-1.))
        lateral = momentum*np.sqrt(np.maximum(0., 1-cosine**2))
        phase = 2*np.pi*part[:,2]
        columns.extend((birth, lateral*np.cos(phase), lateral*np.sin(phase), momentum*cosine))
    return np.column_stack(columns)


def cap(h, marker, prefactor):
    end = marker['time']
    inside = h['t'] < end
    t = np.r_[h['t'][inside], end]
    a = np.interp(t, h['t'], h['a'])
    d = np.interp(t, h['t'], h['D'])
    if d[0] != 0 or np.any(np.diff(d)<0) or np.any(np.diff(a)<0):
        raise ValueError('Nonmonotone source or radius')
    dt = np.diff(t)
    increments = kernel_upper()*prefactor*d[1:]*dt/a[:-1]**3
    depths = np.cumsum(increments[::-1])[::-1]
    weights = np.diff(d)/d[-1]
    return dict(first_cohort_depth_cap=float(depths[0]),
                mean_events_cap=float(weights@depths),
                hit_probability_cap=float(weights@(-np.expm1(-depths))))


def budget(t, a, z, kernel, prefactor, final_rest, initial_distance):
    mid = .5*(t[1:]+t[:-1])
    density = prefactor*final_rest/a**3
    dt = np.diff(t)
    seed = dt*density*kernel
    growth = dt*density*z*kernel_upper()
    # Solve Y'=R0+lambda*Y exactly on each midpoint cell and transport each
    # cell's seed forward through the subsequent cells.
    later_growth = np.cumsum(growth[::-1])[::-1]-growth
    def marked_cap(scale):
        x = scale*growth
        exprel = np.divide(np.expm1(x), x, out=np.ones_like(x), where=(x!=0))
        return float(np.sum(scale*seed*exprel*np.exp(scale*later_growth)))
    h_upper = marked_cap(1.)
    seed_total = float(np.sum(seed))
    # A second evaluation advances the cumulative marked count chronologically.
    h_forward = 0.
    for initial, gain in zip(seed, growth):
        exprel = np.expm1(gain)/gain if gain else 1.
        h_forward = np.exp(gain)*h_forward+initial*exprel
    if abs(h_forward-h_upper)>2.e-13:
        raise ValueError('Disagreeing comparison-equation integrations')
    result = dict(free_pair_encounters_per_particle=seed_total,
                  affected_fraction_upper_estimate=h_upper,
                  encounters_per_particle_upper_estimate=2*h_upper-seed_total,
                  centroid_of_free_pair_events_wall_time=float(np.dot(seed,mid)/seed_total))
    if initial_distance is not None:
        # This is sensitivity to a uniform multiplier of m² sigma v_M on a
        # frozen history, not an assertion the given pole permits such a value.
        def multiplier_for(goal):
            if goal<=0: return 0.
            return float(brentq(lambda x:marked_cap(x)-goal,0.,16.))
        result.update(collisionless_20_cell_tv=initial_distance,
                      minimum_20_cell_tv_estimate=max(0.,initial_distance-h_upper),
                      minimum_rate_multiplier_for_exact_20_cell_match=multiplier_for(initial_distance),
                      minimum_rate_multiplier_for_five_point_20_cell_match=multiplier_for(initial_distance-.05))
    return result


def kernel_controls(exe, root):
    """Head-on resonance and zero comoving velocities in a field-free box."""
    p = np.sqrt((POLE/2)**2-1.)
    table = np.zeros((2,8))
    table[:,0] = [0.,1.]
    table[:,1] = 1.
    table.tofile(root/'control-field')
    resonance = 24*np.pi/(POLE**2*np.sqrt(1-4/POLE**2))
    errors = {}
    for name, pair, predicted in [
            ('rest',[0.,0.,0.,0.,0.,0.,0.,0.],0.),
            ('head_on_resonance',[0.,0.,0.,p,0.,0.,0.,-p],resonance),
            ('parallel',[0.,0.,0.,.75,0.,0.,0.,.75],0.)]:
        np.array([pair], dtype='=f8').tofile(root/'control-pair')
        values = np.loadtxt(io.BytesIO(subprocess.check_output(
            [str(exe),str(root/'control-field'),str(root/'control-pair')])))
        error = abs(values[0]-predicted)
        if error>2.e-15*max(1.e-5,predicted) or values[2]!=1:
            raise ValueError('Failed analytic kernel control '+name)
        errors[name]=float(values[0])
    return errors


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent', type=Path, required=True)
    p.add_argument('--fields', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=2)
    opt = p.parse_args()
    study = opt.parent.resolve()
    history = study/'parent/parent/history'
    inherited = runpy.run_path(str(study/'orbits/pairs.py'))
    inputs = inherited['read_inputs'](history, opt.fields)
    meta = json.loads((history/'markers.json').read_text())
    reference = json.loads((study/'feedback/FEEDBACK.json').read_text())
    distances = {}
    for row in reference['runs']:
        if row['steps']==512 and row['pairs']==4096 and not row['feedback']:
            case = {'two_pressure_fit':'two_pressure_fit',
                    'equal_early_and_late_4pct':'early_late_4pct'}.get(row['case'])
            if case:
                distances[(row['history'],case)] = row['coarse_cell_total_variation']
    tasks = []
    basic = {}
    for label,(h,_,marker) in inputs.items():
        prefactor = meta['runs'][label]['channel_fraction']*meta['action']/(6*meta['daughter_mass_L']**3)
        basic[label] = cap(h, marker, prefactor)
        for case in ('two_pressure_fit','early_late_4pct','axial'):
            settings = [(14,256,398141)]
            if label.endswith('half'):
                settings += [(13,192,589031), (15,384,635291)]
            for power,steps,seed in settings:
                tasks.append((label,case,power,steps,seed,prefactor))
    with tempfile.TemporaryDirectory(prefix='pair-clock-') as tmp:
        root = Path(tmp)
        exe = root/'pair_clock'
        subprocess.run(['g++','-std=c++17','-O3',str(HERE/'pair_clock.cpp'),'-o',str(exe)],check=True)
        controls = kernel_controls(exe, root)
        def calculate(task):
            label,case,power,steps,scramble,prefactor = task
            h,f,mark = inputs[label]
            tab = inherited['field_table'](h,f,mark,2,2,'early','late_birth',steps)
            if case=='two_pressure_fit':
                a1,lx,ly = inherited['FIT'][label.split('_')[0]+'_half'][1]
            elif case=='early_late_4pct':
                a1,lx,ly = .04/np.sqrt(2),0.,.04/np.sqrt(2)
            else:
                a1,lx,ly = 0.,0.,0.
            tab[:,4:6] = a1*tab[:,4:6]+lx*tab[:,6:8]
            tab[:,6:8] *= ly
            tags = pair_tags(h, mark, power, scramble, inherited['directions'])
            with tempfile.TemporaryDirectory(dir=root) as wd:
                folder = Path(wd)
                tab.tofile(folder/'field')
                tags.tofile(folder/'tags')
                raw = subprocess.check_output([str(exe),str(folder/'field'),str(folder/'tags')])
            rates = np.loadtxt(io.BytesIO(raw))
            if rates.shape != (steps,3) or not np.all(np.isfinite(rates)):
                raise ValueError('Incomplete kernel trace')
            t = tab[:,0]
            mid = .5*(t[:-1]+t[1:])
            a = np.interp(mid, h['t'], h['a'])
            z = np.interp(mid, h['t'], h['D'])/mark['D']
            discrepancy = float(np.max(np.abs(rates[:,2]-z**2)))
            if discrepancy>.004 or np.max(rates[:,1])>kernel_upper():
                raise ValueError('Live-pair fraction or pointwise invariant kernel check failed')
            results = budget(t,a,z,rates[:,0],prefactor,mark['D'],distances.get((label,case)))
            if not 0<=results['free_pair_encounters_per_particle']<=basic[label]['mean_events_cap']:
                raise ValueError('Free encounter count outside global cap')
            return dict(history=label, case=case, pairs=len(tags), steps=steps, scramble=scramble,
                        live_pair_fraction_max_integration_error=discrepancy,
                        sampled_maximum_kernel=float(np.max(rates[:,1])), **results)
        with ThreadPoolExecutor(opt.jobs) as pool:
            rows = list(pool.map(calculate, tasks))
    output = dict(kernel_upper_bound=kernel_upper(), kernel_controls=controls, arbitrary_target_upper_bounds=basic, runs=rows)
    opt.output.write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
