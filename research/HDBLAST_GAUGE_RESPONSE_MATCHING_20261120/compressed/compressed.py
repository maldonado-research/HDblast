"""Positive Gaussian quadrature for a finite photon band with backreaction on a closed sphere."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

ROOT=Path(__file__).resolve().parents[1]
SPLIT=48
TOP=192
ORDERS=(1,2,4,8,16,32,48,64,96)


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def gaussian_band(count):
    """Compress the measure (L-1/L) at t=1/L^2; return frequencies and Hamiltonian weights."""
    labels=np.arange(SPLIT+1,TOP+1,dtype=float)
    measure=labels-1/labels
    nodes=(labels[0]/labels)**2  # t divided by its largest exact value
    mass=math.fsum(measure)
    current=np.sqrt(measure/mass)
    previous=np.zeros_like(current)
    previous_norm=0.
    diagonal=[]
    off_diagonal=[]
    basis=[]
    for index in range(count):
        basis.append(current)
        work=nodes*current-previous_norm*previous
        center=float(current@work)
        diagonal.append(center)
        if index==count-1:
            break
        work-=center*current
        # Prevent spurious repeated quadrature nodes from roundoff at larger orders.
        Q=np.stack(basis)
        for _ in range(2):
            work-=Q.T@(Q@work)
        norm=float(np.linalg.norm(work))
        off_diagonal.append(norm)
        previous,current,previous_norm=current,work/norm,norm
    vals,vecs=eigh_tridiagonal(diagonal,off_diagonal,lapack_driver='stev')
    inertia_weights=mass*vecs[0]**2
    effective=labels[0]/np.sqrt(vals)
    order=np.argsort(effective)
    effective=effective[order]
    inertia_weights=inertia_weights[order]
    if min(inertia_weights)<=0 or effective[0]<labels[0]*(1-1e-13) or effective[-1]>labels[-1]*(1+1e-13):
        raise ArithmeticError('Quadrature lost positivity or the physical band')
    # All 2m moments are checked on a dimensionless abscissa, without vanishing powers of t.
    errors=[]
    for j in range(2*count):
        original=float(measure@nodes**j)
        matched=float(inertia_weights@((labels[0]/effective)**(2*j)))
        errors.append(abs(matched-original)/original)
    if max(errors)>1e-11:
        raise ArithmeticError('Gaussian rule missed a susceptibility moment')
    dispersion={}
    for omega in (16.,32.,64.,80.,96.):
        exact=float(np.sum(measure/(1-(omega/(2*labels))**2)))
        approx=float(np.sum(inertia_weights/(1-(omega/(2*effective))**2)))
        dispersion[str(int(omega))]=dict(exact=exact,approximate=approx,relative_error=approx/exact-1)
    return effective,inertia_weights*effective**2,dict(effective_frequencies=effective.tolist(),
                   inertia_weights=inertia_weights.tolist(),max_relative_moment_error=max(errors),
                   static_frequency_response=dispersion)


def simplify(result):
    if result['status']!='first_minimum' or result['switch_end'] is None:
        raise ArithmeticError('Completed outgoing interval required')
    fields=('x','a','h','v','photon_density','scalar_density','photon_fraction',
            'reference_band_fraction','reference_band_comoving','added_band_comoving','comoving_energy')
    out={stage:{key:float(result[stage][key]) for key in fields} for stage in ('end','switch_end')}
    out['energy_residuals']=result['errors']
    if abs(out['end']['comoving_energy']-out['switch_end']['comoving_energy'])>1e-8:
        raise ArithmeticError('Outgoing comoving energy changed on a constant plateau')
    for row in (out['end'],out['switch_end']):
        if not 0<=row['reference_band_fraction']<=row['photon_fraction']+1e-10:
            raise ArithmeticError('Nonpositive added-band energy')
        if abs(row['reference_band_comoving']+row['added_band_comoving']-row['comoving_energy'])>1e-9:
            raise ArithmeticError('Band energies do not add')
    return out


def calculate(prior):
    matched=module('positive_photon_dynamics',ROOT/'matched/matched.py')
    old=json.loads((prior/'feedback/FEEDBACK.json').read_bytes())['branches']
    rules={m:gaussian_band(m) for m in ORDERS}
    cases={}
    for name,branch in matched.branches(prior).items():
        rows={}
        for delta in (.1,3.):
            ref=simplify(matched.run(branch,TOP,delta,match=False))
            saved=old[name]['runs'][f'{TOP}/{delta:g}']['end']
            checks={k:abs(saved[a]-ref['end'][k]) for a,k in
                    [('x','x'),('a','a'),('h','h'),('v','v'),('photon_density','photon_density'),
                     ('fraction_of_matter','photon_fraction'),('comoving_excess','comoving_energy')]}
            if max(checks.values())>3e-8:
                raise ArithmeticError('Full-band calculation disagrees with archived feedback')
            reduced={}
            for m,(labels,weights,rule) in rules.items():
                results=simplify(matched.run(branch,TOP,delta,match=False,high_band=(labels,weights)))
                a,e=results['end'],ref['end']
                for key in ('x','a','h','v','photon_fraction','reference_band_fraction','added_band_comoving'):
                    results.setdefault('first_minimum_fractional_errors',{})[key]=(a[key]-e[key])/e[key]
                reduced[str(m)]=results
            direct=simplify(matched.run(branch,TOP,delta,match=False,direct=True,high_band=rules[8][:2]))
            differences={k:abs(v-direct['end'][k]) for k,v in reduced['8']['end'].items()}
            if max(differences.values())>2e-7:
                raise ArithmeticError('Quadrature and Bogoliubov descriptions differ')
            rows[f'{delta:g}']=dict(full=ref,archived_differences=checks,rules=reduced,
                                      independent_n8_coordinate_differences=differences)
        cases[name]=rows
    return dict(reference_cutoff=SPLIT,target_cutoff=TOP,orders=list(ORDERS),
                weight='positive linear-response inertia at inverse squared frequency; dynamical energy weight = inertia times squared mode label',
                rules={str(m):r for m,(_,_,r) in rules.items()},branches=cases)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior',type=Path,default=ROOT/'prior')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior),indent=2,allow_nan=False)+'\n')
