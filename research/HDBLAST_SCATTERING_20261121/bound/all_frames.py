"""Cover every inertial boost for energy-preserving turns, including oblique boosts."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import quad
from scipy.special import gammainc, kve


def check_primitive(executable):
    """Independent angular integration, including the light-speed and clipped corners."""
    arguments = [(x,s,h,l) for x in (.125,1.25,8.)
                 for s,h in [(0.,0.),(.003,-.73),(.25,.64),(.75,0.),(1.,.91),(1.2,-.62)]
                 for l in (-.2,.41,1.,1.72)]
    queries = ''.join(' '.join(map(str,row))+'\n' for row in arguments)
    values = subprocess.check_output([executable,'--check-cdf'],input=queries,text=True)
    errors = []
    for (x,s,h,l),value in zip(arguments,values.split(),strict=True):
        def density(mu):
            multiplier = l + s*mu
            cdf = gammainc(3.,x/multiplier) if multiplier>0. else 1.
            return .5*(1+h*mu)*cdf
        points = [-l/s] if s and abs(l/s)<1. else []
        exact = quad(density,-1.,1.,points=points,epsabs=1.e-13,epsrel=1.e-13)[0]
        errors.append(abs(float(value)-exact))
    if max(errors)>5.e-10:
        raise ValueError('Primitive disagrees with independent angular integral')
    return len(arguments)


def exact_target(theta, drift_beta, velocity, cutoff):
    """Integral of the massive thermal density, evaluated in the stream rest frame."""
    eps = 1./theta
    speed = float(np.linalg.norm(velocity))
    slope = drift_beta*velocity[2]/speed if speed else 0.
    root = math.sqrt(cutoff**2-(1.-speed**2)*eps**2)
    low = (cutoff**2-eps**2)/(root+speed*cutoff)
    high = (root+speed*cutoff)/(1.-speed**2) if speed<1. else np.inf
    norm = theta**2/kve(2,eps)
    def radial(q):
        return norm*q*q*math.exp(eps-math.hypot(q,eps))
    def partial(q):
        energy = math.hypot(q,eps)
        mu = min(1.,max(-1.,(cutoff-energy)/(speed*q)))
        return radial(q)*(.5*(1.+mu)+.25*slope*(q/energy)*(mu*mu-1.))
    result = quad(radial,0.,low,epsabs=1.e-12)[0]
    if speed:
        result += quad(partial,low,min(high,80.),epsabs=1.e-12)[0]
    return result


def check_witness(certificate, rest, theta, beta):
    center = np.asarray(certificate['weakest_center'])
    width = certificate['weakest_halfwidth']
    cutoff = certificate['weakest_threshold']
    bound = certificate['lower_bound']
    sign = certificate['weakest_sign_of_difference']
    # Any point inside the intersection of this box with the unit velocity ball
    # is a valid witness check. Include all physical corners and the nearest point.
    nearest = np.sign(center)*np.maximum(0.,np.abs(center)-width)
    positions = [nearest,center]
    for bits in range(8):
        shifts = np.array([1. if bits & (1<<dim) else -1. for dim in range(3)])
        positions.append(center+width*shifts)
    checked = 0
    for velocity in positions:
        if np.linalg.norm(velocity)>1.:
            continue
        energies = rest[:,0]+rest[:,1:] @ velocity
        observed = np.count_nonzero(energies<=cutoff)/len(energies)
        thermal = exact_target(theta,beta,velocity,cutoff)
        if sign*(observed-thermal)<bound-1.e-9:
            raise ValueError('Direct massive CDF violates the accepted cube witness')
        checked += 1
    if checked==0:
        raise ValueError('No physical velocity checked in the limiting cube')
    return checked


def run_case(path, executable, folder):
    with np.load(path,allow_pickle=False) as sample:
        momenta = sample['p']
        charge = sample['sign']
        theta,drift = float(sample['theta']),float(sample['drift'])
    boost = math.hypot(1.,drift)
    beta = abs(drift)/boost
    parallel = np.sign(drift)*charge*momenta[:,2]
    gamma = np.hypot(1.,np.linalg.norm(momenta,axis=1))
    # For each sign align its target drift with the third axis. The choice of
    # transverse basis is immaterial since the certificate covers every boost.
    rest = np.column_stack((boost*gamma-abs(drift)*parallel,
                            momenta[:,:2],boost*parallel-abs(drift)*gamma))/theta
    descriptions = []
    for sign in (-1,1):
        data = rest[charge==sign]
        if len(data)*2!=len(charge):
            raise ValueError('Unequal charge weights')
        temp = folder/'four-vectors.bin'
        data.tofile(temp)
        line = subprocess.check_output([executable,temp,str(theta),str(beta),'.05'],text=True)
        bound,visited,accepted,outside,depth,x,y,z,radius,cut,direction = line.split()
        record = dict(charge=sign,lower_bound=float(bound),boxes_visited=int(visited),
                      boxes_certified=int(accepted),boxes_outside=int(outside),maximum_depth=int(depth),
                      weakest_center=list(map(float,(x,y,z))),weakest_halfwidth=float(radius),
                      weakest_threshold=int(cut)*.125,weakest_sign_of_difference=int(direction))
        record['independent_massive_cdf_checks'] = check_witness(record,data,theta,beta)
        descriptions.append(record)
    return dict(snapshot=path.name,per_charge_certificates=descriptions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshots',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='hdblast-all-frames-') as temporary:
        folder = Path(temporary)
        executable = folder/'all_frames'
        subprocess.run(['g++','-std=c++17','-O3',str(Path(__file__).with_suffix('.cpp')),
                        '-o',executable],check=True)
        controls = check_primitive(executable)
        results = [run_case(args.snapshots/f'{model}_half_{pulse}_{power}.npz',executable,folder)
                   for model in ('dedicated','shared') for pulse in ('pressure_fit','small_pulse')
                   for power in (12,14)]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(requested_cdf_gap=.05,primitive_checks=controls,
        runs=results),sort_keys=True,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
