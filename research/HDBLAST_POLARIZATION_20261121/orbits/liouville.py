"""Bound the fine-grained support volume swept out by monoenergetic charged births."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import kve


Q=1.e7


def upper_sum(h,marker,electric_norm):
    """Enclose the integrals on linear history cells, including the partial last cell."""
    end=marker['time']
    stop=h['t'].searchsorted(end)
    if not 0<stop<len(h['t']) or h['D'][0]!=0.:raise ValueError('Incomplete birth history')
    def knots(row):return np.r_[row[:stop],np.interp(end,h['t'],row)]
    t=knots(h['t'])
    t[-1]=end
    a,R,field=(knots(x) for x in [h['a'],h['R'],electric_norm])
    if np.min(R)<0. or np.min(a)<=0. or np.min(field)<0.:raise ValueError('Invalid history')
    maximum=lambda seq:np.maximum(seq[:-1],seq[1:])
    mass2=np.maximum(0.,R*R-1.)
    a3=a**3; mass3=mass2**1.5
    rest_variation=np.sum(maximum(a3)*abs(np.diff(mass3))+maximum(mass3)*abs(np.diff(a3)))/3.
    sweep=.5*Q*np.sum(np.diff(t)*maximum(a3)*maximum(mass2)*maximum(field))
    theta=marker['theta']; gamma=np.hypot(1.,np.interp(end,h['t'],h['u']))
    normalization=a[-1]**3*gamma*theta*kve(2,1./theta)
    probability=(rest_variation+sweep)/normalization
    return dict(time=end,theta=theta,drift_gamma=float(gamma),
                birth_radius_volume_term=float(rest_variation),electric_volume_term=float(sweep),
                thermal_volume_scale=float(normalization),thermal_mass_in_support_upper=float(min(1.,probability)),
                total_variation_lower=float(max(0.,1.-probability)))


def bank_norm(waves,source):
    """Worst RMS field of any launch with unit total emitted energy in the specified bank."""
    names=source['profile_order'];gram=np.array(source['emitted_energy_gram'])
    upper=np.zeros_like(waves['t'])
    for k in np.unique(waves['kappa']):
        rows=np.flatnonzero(waves['kappa']==k)
        ids=[names.index(p) for p in waves['profile'][rows]]
        L=np.linalg.cholesky(gram[np.ix_(ids,ids)])
        white=solve_triangular(L,waves['E'][rows],lower=True)
        upper=np.maximum(upper,np.linalg.norm(white,axis=0))
    return upper


def load(path):
    with np.load(path) as z:return dict(z)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('history','fields','scan','source','pairs','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    pairs=json.loads(args.pairs.read_text())['runs']
    inputs={};bank_cases={};example_cases=[]
    markers=json.loads((args.history/'markers.json').read_text())['runs']
    for bank,folder,metadata in [('primary',args.fields,'FIELDS.json'),('storage',args.scan,'STORAGE_SCAN.json')]:
        tables=json.loads((folder.parent/metadata).read_text())
        for label in sorted(tables):
            history='dedicated_native' if label=='dedicated_spatial' else label
            h=load(args.history/(history+'.npz')); marker=markers[history]['markers']['0.01']
            waves=load(folder/(label+'.npz'))
            if not np.array_equal(h['t'],waves['t']):raise ValueError('History grid mismatch')
            g=bank_norm(waves,tables[label])
            bank_cases[bank+'_'+label]={str(eps):upper_sum(h,marker,abs(h['X'])+eps*g)
                                      for eps in (0.,.02,1.)}
            inputs[('scan_' if bank=='storage' else '')+label]=(h,marker,waves)
    for case in pairs:
        label=case['history']
        if case['case']=='evolving_mean_source':
            h=load(args.source/(label+'.npz'))
            marker={'time':float(h['marker_time']),'theta':float(np.interp(h['marker_time'],h['t'],h['theta']))}
            e1,e2=h['mode_E']
        elif case['case'] in ('uniform_local_pressure_match','pooled_pressure_match') or (
              case['case']=='distributed_storage' and case['k2']==515):
            h,marker,waves=inputs[label]
            def coeff(j):
                r=np.flatnonzero((waves['kappa']==case['k'+j])&(waves['profile']==case['profile'+j]))
                if len(r)!=1:raise ValueError('Ambiguous trace')
                return case['amplitude'+j]*waves['E'][r[0]]
            e1,e2=coeff('1'),coeff('2')
        else:continue
        # All selected modes are horizontal. Equal-k pairs in this list use the same
        # normalized spatial polynomial at fixed angular phase; unequal k are orthogonal.
        overlap=np.cos(case['phase2']) if case['k1']==case['k2'] else 0.
        rms=np.sqrt(np.maximum(0.,h['X']**2+e1**2+e2**2+2.*overlap*e1*e2))
        example_cases.append(dict(case=case['case'],history=label,**upper_sum(h,marker,rms)))
    # A force-free gas injected on an expanding momentum sphere fills exactly the ball.
    # The numerator and the starting volume vanish separately for a stationary shell.
    h=dict(t=np.array([0.,1.,2.]),a=np.ones(3),D=np.array([0.,1.,1.]),
           R=np.sqrt(np.array([1.,5.,5.])),u=np.zeros(3))
    test=upper_sum(h,dict(time=1.,theta=1.),np.zeros(3))
    if abs(test['birth_radius_volume_term']-8./3.)>1e-12 or test['electric_volume_term']!=0.:
        raise ValueError('Exact free expansion sphere benchmark failed')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(source_banks=bank_cases,particle_cases=example_cases,
                                          force_free_ball_check=test),allow_nan=False,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':main()
