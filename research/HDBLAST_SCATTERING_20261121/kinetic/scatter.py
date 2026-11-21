"""Relax archived charge streams at fixed volume with relativistic binary collisions."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import brentq
from scipy.special import kve
from scipy.stats import qmc


ENERGIES = np.array([0., 1., 2., 3., 5., np.inf])
COSINES = np.linspace(-1., 1., 5)
# Checkpoints measured in sixteenths of one collision per particle, on average.
CHECKPOINTS = (0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
HERE = Path(__file__).resolve().parent


def gamma(p):
    return np.sqrt(1. + np.einsum('ij,ij->i', p, p))


def enthalpy(theta):
    return kve(1, 1./theta)/kve(2, 1./theta) + 4.*theta


def equilibrium(p):
    mean_energy = gamma(p).mean()
    mean_p = p.mean(axis=0)
    magnitude = np.linalg.norm(mean_p)

    def residual(t):
        h = enthalpy(t)
        boost = np.hypot(1., magnitude/h)
        return h*boost - t/boost - mean_energy

    theta = brentq(residual, 1e-5, mean_energy)
    return theta, mean_p/enthalpy(theta)


def target_cells(theta, u, axis=None):
    """Normalized boosted Juttner counts on a common lab slice, in rest coordinates."""
    projected = np.linalg.norm(u) if axis is None else np.dot(u, axis)
    beta = projected/np.hypot(1., np.linalg.norm(u))
    norm = kve(2, 1./theta)
    cells = []
    for lo, hi in zip(ENERGIES[:-1], ENERGIES[1:]):
        def density(x):
            return np.exp(-x)*(1. + theta*x)*np.sqrt(theta*x*(2. + theta*x))/norm

        def odd(x):
            return np.exp(-x)*theta*x*(2. + theta*x)/norm

        even_int = quad(density, lo, hi, epsabs=1e-12)[0]
        odd_int = quad(odd, lo, hi, epsabs=1e-12)[0]
        cells.append(.5*even_int*np.diff(COSINES)+.25*beta*odd_int*np.diff(COSINES**2))
    arr = np.array(cells)
    if not np.isclose(arr.sum(), 1., atol=1e-10) or arr.min() < 0.:
        raise ValueError('Invalid thermal probability')
    return arr


def histogram(p, theta, u, axis=None):
    magnitude = np.linalg.norm(u)
    direction = u/magnitude if magnitude else np.array([0., 0., 1.])
    boost = np.hypot(1., magnitude)
    g = gamma(p)
    along = p @ direction
    perpendicular = p - along[:, None]*direction
    rest_along = boost*along - magnitude*g
    rest_p = np.hypot(rest_along, np.linalg.norm(perpendicular, axis=1))
    rest_kinetic = rest_p**2 / (1. + np.hypot(1., rest_p))
    if axis is not None:
        rest_along = perpendicular @ axis + rest_along*np.dot(direction,axis)
    cosine = np.divide(rest_along, rest_p, out=np.zeros_like(rest_p), where=rest_p>0)
    obs = np.histogram2d(rest_kinetic/theta, np.clip(cosine,-1.,1.), (ENERGIES,COSINES))[0]/len(p)
    if abs(obs.sum()-1.)>1e-13:
        raise ArithmeticError('Unbinned particle')
    return obs


def sample_juttner(theta, u, power, seed):
    """Lab-count sample via the exact linear tilt of the rest-frame cosine density."""
    q = qmc.Sobol(3, scramble=True, seed=seed).random_base2(power)
    x_grid = np.linspace(0.,np.sqrt(62.),65537)**2
    y = theta*x_grid
    prob = np.exp(-x_grid)*(1.+y)*np.sqrt(y*(2.+y))
    cdf = cumulative_trapezoid(prob,x_grid,initial=0.)
    x = np.interp(q[:,0],cdf/cdf[-1],x_grid)
    rest_g = 1.+theta*x
    rest_p = np.sqrt(theta*x*(2.+theta*x))
    size = np.linalg.norm(u)
    boost = np.hypot(1.,size)
    tilt = (size/boost)*rest_p/rest_g
    mu = (4.*q[:,1]-2.+tilt)/(1.+np.sqrt((1.-tilt)**2+4.*tilt*q[:,1]))
    azimuth = 2.*np.pi*q[:,2]
    axis = u/size if size else np.array([0.,0.,1.])
    cross_axis = np.array([0.,0.,1.]) if abs(axis[2])<.75 else np.array([1.,0.,0.])
    first = np.cross(axis,cross_axis)
    first /= np.linalg.norm(first)
    second = np.cross(axis,first)
    transverse = rest_p*np.sqrt(np.maximum(0.,1.-mu**2))
    return (np.outer(transverse*np.cos(azimuth),first)
            + np.outer(transverse*np.sin(azimuth),second)
            + np.outer(boost*rest_p*mu+size*rest_g,axis))


def tv(observed, ref):
    return float(.5*np.abs(observed-ref).sum())


def split(p):
    return p[::2],p[1::2]


def moments(p):
    return np.r_[gamma(p).mean(),p.mean(axis=0)]


def reference_info(p, theta, drift, kind):
    orig = [(theta,np.array([0.,0.,sign*drift])) for sign in (1,-1)]
    fit = [equilibrium(row) for row in split(p)] if kind=='same' else [equilibrium(p)]*2
    archived = target_cells(*orig[0])
    fitted = [target_cells(t,u) for t,u in fit]
    eq_distances = []
    for power, seed in [(16,81),(17,137)]:
        eq_cells = sum(histogram(sample_juttner(t,u,power,seed+i),*orig[i])
                       for i,(t,u) in enumerate(fit))/2.
        eq_distances.append(tv(eq_cells,archived))
    return orig, fit, archived, fitted, eq_distances


def simulate(exe, trial, temp):
    name, case, power, kind, cone, seed, p, theta, drift = trial
    if not np.isfinite(p).all() or len(p)%4:
        raise ValueError('Invalid input shape')
    orig, fit, archived, fitted, predicted = reference_info(p,theta,drift,kind)
    e0 = moments(p)[0]
    beam_axis = np.array([0., 0., 1.])
    beam_cells = [target_cells(t,u,beam_axis) for t,u in fit]
    baseline = [moments(s) for s in (*split(p),p)]
    schedule = CHECKPOINTS + (1024,2048) if cone==30. else CHECKPOINTS
    rows = []
    with tempfile.TemporaryDirectory(dir=temp) as wd:
        d=Path(wd)
        (d/'schedule').write_text(' '.join(map(str,schedule)))
        p.astype('=f8').tofile(d/'p')
        result = subprocess.check_output([str(exe),str(d/'p'),str(d/'snapshot'),
                                          kind,str(cone),str(seed),str(d/'schedule')],text=True)
        for line in result.splitlines():
            raw, exposure, accepted, proposed, cross_sign, error, majorant = line.split()
            tick=int(raw)
            after = np.fromfile(d/('snapshot'+raw),dtype='=f8').reshape(-1,3)
            if (kind=='opposite' and int(accepted)!=int(cross_sign)) or (kind=='same' and int(cross_sign)!=0):
                raise ArithmeticError('A forbidden collision channel occurred')
            if after.shape!=p.shape or not np.isfinite(after).all():
                raise ArithmeticError('Missing particle record')
            signs=split(after)
            count=sum(histogram(s,*par) for s,par in zip(signs,orig))/2.
            fit_tv=[tv(histogram(s,*par),cells) for s,par,cells in zip(signs,fit,fitted)]
            beam_tv = [tv(histogram(s,*par,beam_axis),cells) for s,par,cells in zip(signs,fit,beam_cells)]
            g = gamma(after)
            pressure = np.einsum('ij,ik,i->jk',after,after,1./g) / g.sum()
            current = [moments(s) for s in (*signs,after)]
            changes=[float(np.max(np.abs(now-then))/e0) for now,then in zip(current,baseline)]
            protected=changes if kind=='same' else changes[-1:]
            if max(protected)>2e-11 or float(error)>2e-11 or float(majorant)>1.+5e-9:
                raise ArithmeticError(('Conservation/majorant error',protected,error,majorant))
            mean_v=[np.mean(s[:,2]/gamma(s)) for s in signs]
            rows.append(dict(collisions_per_particle=tick/16.,exposure=float(exposure),
                             accepted=int(accepted),proposed=int(proposed),opposite_sign=int(cross_sign),
                             archived_tv=tv(count,archived), archived_energy_tv=tv(count.sum(1),archived.sum(1)),
                             archived_direction_tv=tv(count.sum(0),archived.sum(0)), fitted_tv_per_sign=fit_tv,
                             fitted_beam_axis_tv_per_sign=beam_tv,
                             wall_pressure_eigenvalues_over_energy=np.linalg.eigvalsh(pressure).tolist(),
                             mean_axial_v=mean_v, max_abs_four_moment_change_over_initial_energy=changes,
                             max_pair_energy_residual=float(error),max_acceptance=float(majorant)))
    if len(rows)!=len(schedule):
        raise ArithmeticError('Incomplete schedule')
    return dict(history=name,case=case,particles_per_sign=2**power,channels=kind,cone_degrees=cone,
                collision_seed=seed,archived_temperature=theta,archived_signed_drift=float(drift),
                conserved_equilibrium=[dict(theta=t,u=u.tolist()) for t,u in fit],
                predicted_equilibrium_archived_tv_qmc=predicted,initial_four_moments=[m.tolist() for m in baseline],
                samples=rows)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--particles',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--jobs',type=int,default=2)
    args=ap.parse_args()
    work=[]
    for n in ('dedicated_half','shared_half'):
        for case in ('pressure_fit','small_pulse'):
            for power,seed in ((12,7741),(12,12591),(13,31831),(14,84779)):
                with np.load(args.particles/f'{n}_{case}_{power}.npz') as sample:
                    p,theta,drift=sample['p'],float(sample['theta']),float(sample['drift'])
                    if not np.array_equal(sample['sign'],np.tile([1.,-1.],2**power)):
                        raise ValueError('Invalid charge ordering')
                    obs=sum(histogram(s,theta,np.array([0.,0.,sign*drift])) for sign,s in zip((1.,-1.),split(p)))/2.
                    if abs(tv(obs,target_cells(theta,np.array([0.,0.,drift])))-sample['reference_tv'])>1e-12:
                        raise ValueError('Failure to reproduce the archived baseline')
                channels = ['same', 'mixed']
                if case=='small_pulse' and seed in (7741,84779):
                    channels.append('opposite')
                for kind in channels:
                    cones=(180.,30.) if power in (12,14) and seed!=12591 and kind!='opposite' else (180.,)
                    for cone in cones:
                        work.append((n,case,power,kind,cone,seed,p,theta,drift))
    # Both tags of this control start in the same drifting equilibrium.
    up=np.array([0.,0.,3.])
    control=np.stack([sample_juttner(5.,up,13,s) for s in (880,991)],axis=1).reshape(-1,3)
    with tempfile.TemporaryDirectory(prefix='elastic-study-') as name:
        wd=Path(name)
        exe=wd/'elastic'
        subprocess.run(['g++','-O3','-ffp-contract=off','-std=c++17',str(HERE/'elastic.cpp'),'-o',str(exe)],check=True)
        with ThreadPoolExecutor(args.jobs) as pool:
            results=list(pool.map(lambda x:simulate(exe,x,wd),work))
        control_result=simulate(exe,('equilibrium_control','two_comoving_tags',13,'mixed',180.,771,control,5.,3.),wd)
    if max(max(r['fitted_tv_per_sign']) for r in control_result['samples'])>.04:
        raise ArithmeticError('Drifting-equilibrium control moved beyond counting tolerance')
    out=dict(runs=results,equilibrium_control=control_result)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out, allow_nan=False,sort_keys=True,indent=2)+'\n')


if __name__=='__main__':
    main()
