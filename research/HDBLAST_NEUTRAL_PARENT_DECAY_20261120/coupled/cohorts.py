"""Resolve classical collisionless daughter packets on the prescribed zero-drag fluid trajectory."""
import argparse
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
TIMES = [10.4,10.5,10.7,11.,11.2,11.5,12.,12.4]


def projection(s, endpoint, stride, angles, q):
    mu,w=leggauss(angles)
    a=s['a'][endpoint]
    field=s['a']**2*s['B']/2
    left=np.arange(0,endpoint,stride)
    right=np.minimum(endpoint,left+stride)
    gain=s['f'][right]-s['f'][left]
    live=gain>0.
    left,right,gain=left[live],right[live],gain[live]
    if len(gain)==0:
        raise ValueError('There are no daughter particles at the comparison time')
    weights=(gain/gain.sum())[:,None]*w[None,:]/2
    r_birth=(s['R'][left]+s['R'][right])/2
    p_birth=np.sqrt(np.maximum(0,r_birth**2-1))
    a_birth=(s['a'][left]+s['a'][right])/2
    a_field=(field[left]+field[right])/2
    parallel=(a_birth[:,None]*p_birth[:,None]*mu[None,:]+q*(field[endpoint]-a_field[:,None]))/a
    perp2=(a_birth*p_birth/a)[:,None]**2*(1-mu**2)[None,:]
    total=np.sqrt(1.+parallel**2+perp2)
    def average(x):
        return float(np.sum(weights*x))
    velocity=average(parallel/total)
    emean=average(total)
    moment=average(parallel)
    tpar=average(parallel**2/total)
    tperp=average(.5*perp2/total)
    # Positive null projections avoid subtracting nearly equal beam energy
    # and momentum densities when almost every particle co-propagates.
    ep=total+parallel
    em=(1.+perp2)/ep
    # For a beam progressing in the negative direction, prefer the other root.
    em=np.where(parallel<0.,total-parallel,em)
    ep=np.where(parallel<0.,(1.+perp2)/em,ep)
    disc=np.sqrt(average(ep**2/total)*average(em**2/total))
    v_landau=2.*moment/(emean+tpar+disc)
    lorentz=1./np.sqrt((1.-v_landau)*(1.+v_landau))
    energy_rest=.5*(emean-tpar+disc)
    number_rest=lorentz*(1.-v_landau*velocity)
    return dict(mean_current_per_charge_density=velocity,mean_wall_energy_per_mass=emean,
                mean_wall_momentum_per_mass=moment,landau_velocity=float(v_landau),
                mean_landau_kinetic_energy_per_mass=float(energy_rest/number_rest-1.),
                longitudinal_pressure_per_lab_particle=tpar,
                transverse_pressure_per_lab_particle=tperp)


def calculate():
    spec=importlib.util.spec_from_file_location('local_decay_solver',HERE/'decay.py')
    solver=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(solver)
    case=json.loads((HERE/'cases.json').read_bytes())['no_explicit_drag']
    with TemporaryDirectory(prefix='hdblast-vlasov-packets-') as tmp:
        exe=Path(tmp)/'wall'
        solver.compile_solver(exe)
        run=solver.run(exe,case|{'trace':2})
    s=dict(zip(solver.COLUMNS,np.array(run['trace']).T,strict=True))
    if tuple(sorted(s['t'])) != tuple(s['t']):
        raise ValueError('Expected an increasing proper-time trace')
    result=[]
    largest_velocity_change=0.
    largest_angular_change=0.
    for t in TIMES:
        i=int(np.argmin(np.abs(s['t']-t)))
        if abs(s['t'][i]-t)>1.e-10:raise ValueError('Comparison clock changed')
        main=projection(s,i,1,24,case['Q'])
        coarse=projection(s,i,2,24,case['Q'])
        fewer_angles=projection(s,i,1,12,case['Q'])
        neutral=projection(s,i,1,24,0.)
        f=s['f'][i]
        b=case['D']/case['R_low']*(2./s['a'][i])**3
        fluid_v=s['current'][i]/(b*case['Q']*f)
        if abs(neutral['mean_current_per_charge_density'])>1.e-13:
            raise ValueError('Isotropic unforced packets acquired an electric current')
        err=abs(main['mean_current_per_charge_density']-coarse['mean_current_per_charge_density'])
        ang_err=abs(main['mean_current_per_charge_density']-fewer_angles['mean_current_per_charge_density'])
        largest_velocity_change=max(largest_velocity_change,err)
        largest_angular_change=max(largest_angular_change,ang_err)
        result.append(dict(t=t,f=float(f),exact_packets_on_given_field=main,fluid_current_per_charge_density=float(fluid_v),
                      fluid_mean_rest_kinetic_per_mass=float(s['intrinsic_random'][i]/(b*f)),
                      half_time_resolution=coarse,absolute_angular_change=ang_err))
    if largest_velocity_change>.022 or largest_angular_change>1.e-5:
        raise ValueError('Birth quadrature no longer resolves the reported current differences')
    late=next(row for row in result if row['t']==12.)
    if abs(late['fluid_current_per_charge_density']-late['exact_packets_on_given_field']['mean_current_per_charge_density'])<.6:
        raise ValueError('The reported closure disagreement disappeared')
    return dict(source='Fresh zero-explicit-drag fluid run; second-order birth-time interpolation',
                early_cells=case['cells'],parameters=case,birth_time_spacing=.001,
                fine_angular_nodes=24,comparison_angular_nodes=12,
                largest_current_change_halving_time_resolution=largest_velocity_change,
                largest_current_change_halving_angle_resolution=largest_angular_change,
                comparison=result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(),indent=2,sort_keys=True,allow_nan=False)+'\n')
