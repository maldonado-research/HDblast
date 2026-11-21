"""Energy lower bounds for arbitrary launched angular patterns within the measured harmonic band."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import kv, kve


def cumulative_floor(cost, weights, pstar, temperature, uniform):
    """Exhibit single momentum cutoffs where the CDF cannot match the drifting thermal target."""
    cut=temperature*np.geomspace(.002,14.,420)
    cut=cut[cut>pstar]
    w=np.hypot(1.,cut)
    thermal=1.-w*w*kv(2,w/temperature)/kv(2,1./temperature)
    # A cohort above the cutoff must pick up more than cut-pstar in transverse impulse.
    # Cost decreases strictly with later birth, so only the tail of the cohort sum contributes.
    order=np.argsort(cost); c=cost[order]; mass=weights[order]
    p0=np.r_[0.,np.cumsum(mass)]
    p2=np.r_[0.,np.cumsum(mass*c*c)]
    def witness(amplitude):
        threshold=(cut-pstar)/amplitude
        j=np.searchsorted(c,threshold,side='right')
        fraction=p0[j].copy()
        if uniform:
            fraction-=p2[j]/threshold**2
        difference=fraction-thermal
        n=int(np.argmax(difference))
        return dict(guaranteed_excess_cdf=float(max(0.,difference[n])),cutoff_in_mc=float(cut[n]),
                    certified_fraction_below_cutoff=float(fraction[n]),thermal_fraction_below_cutoff=float(thermal[n]))
    lo,hi=1.e-12,1.
    while witness(hi)['guaranteed_excess_cdf']>.01: hi*=2.
    for _ in range(50):
        mid=.5*(lo+hi)
        if witness(mid)['guaranteed_excess_cdf']>.01: lo=mid
        else: hi=mid
    return dict(necessary_launch_rms_for_one_percentage_point=lo,
                limiting_witness=witness(lo), witnesses={str(v):witness(v) for v in (.001,.005,.01,.02)})


def bounds(history, fields, marker):
    modes=fields['kappa']
    if not np.array_equal(modes,np.arange(2,41)) or not np.array_equal(fields['t'],history['t']):
        raise ValueError('Expected consecutive true-collar harmonics and matching history knots')
    t=history['t']
    j=t.searchsorted(marker['time'])
    nodes=np.r_[t[:j],marker['time']]
    dt=np.diff(nodes)
    a=np.interp(nodes,t,history['a']); born=np.interp(nodes,t,history['D']); x=np.interp(nodes,t,history['X'])
    if not (np.all(dt>0.) and np.all(np.diff(a)>=0.) and np.all(np.diff(born)>=0.) and np.all(x<0.) and born[0]==0.):
        raise ValueError('Unexpected interval or reference axial force')
    def endpoint_envelope(coeff):
        amp=np.array([np.interp(nodes,t,v) for v in coeff])
        return np.maximum(abs(amp[:,:-1]),abs(amp[:,1:]))
    electric=endpoint_envelope(fields['E_cos']); magnetic=endpoint_envelope(fields['B_sin'])
    count=2*(modes**2-1)
    # Either helicity has identical electric transfer and equal magnetic magnitude.
    pointwise=np.sqrt((count/3.)@(electric**2+magnetic**2))
    derivative_weight=count*(modes**2-2)/3.
    grad=np.sqrt(derivative_weight@electric**2)+np.sqrt(derivative_weight@magnetic**2)
    # The integrated causal displacement, summed over all daughters then present.
    travel=np.cumsum(dt*born[1:]/a[:-1])
    born_force=born[1:]*(electric.max(axis=0)+magnetic.max(axis=0)) + travel*grad
    universal=born[1:]*pointwise
    loss= np.minimum(universal,born_force)
    charge=1.e7
    normalize=charge/(marker['a']*marker['D'])
    inherited=normalize*np.dot(dt*a[:-1],born[:-1]*np.minimum(-x[:-1],-x[1:]))
    pstar=float(np.sqrt(max(history['R'][:j+1].max()**2-1.,0.)))
    thermal=float(kve(1,1./marker['theta'])/kve(2,1./marker['theta'])+3.*marker['theta'])
    output={}
    for name,force in [('arbitrary_birth_positions',universal),('uniform_birth_positions_finite_travel',loss)]:
        cost=float(normalize*np.dot(dt*a[1:],force))
        deficit=inherited-pstar-thermal
        output[name]=dict(axial_loss_per_unit_outer_rms_upper=cost,
            minimum_launch_rms_to_avoid_excess=float(deficit/cost),
            launch_squared_fraction_required=float((deficit/cost)**2),
            mean_energy_over_thermal_lower={str(alpha):float(max(1.,inherited-pstar-alpha*cost)/thermal)
                                         for alpha in (.005,.01,.02,.05,.1,.16)})
    def tail(v): return np.cumsum(v[::-1])[::-1]
    eta=np.r_[0.,np.cumsum(dt/a[:-1])]
    stage=dt*a[1:]
    local=electric.max(axis=0)+magnetic.max(axis=0)
    everywhere=np.sqrt((count/3.)@electric**2)+np.sqrt((count/3.)@magnetic**2)
    worst=charge/marker['a']*tail(stage*everywhere)
    cohort=charge/marker['a']*(tail(stage*local)+tail(stage*eta[1:]*grad)-eta[:-1]*tail(stage*grad))
    birth_mass=np.diff(born)/marker['D']
    if not np.isclose(birth_mass.sum(),1.): raise ValueError('Birth weights are inconsistent')
    direct={
        'arbitrary_birth_positions':cumulative_floor(worst,birth_mass,pstar,marker['theta'],False),
        'uniform_birth_positions_finite_travel':cumulative_floor(np.minimum(cohort,worst),birth_mass,pstar,marker['theta'],True)}
    return dict(transverse_cdf_test=direct, marker_time=float(nodes[-1]), real_angular_coefficients=int(count.sum()),kappa_min=int(modes.min()),
                kappa_max=int(modes.max()),max_comoving_travel=float(np.sum(dt/a[:-1])),
                inherited_axial_impulse_lower=float(inherited),max_birth_momentum=pstar,thermal_mean_gamma=thermal,
                bounds=output)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history',type=Path,required=True);p.add_argument('--geometry',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    meta=json.loads((args.history/'markers.json').read_bytes())
    results={}
    for label in ('dedicated_native','dedicated_half','shared_native','shared_half','dedicated_spatial'):
        history_label='dedicated_native' if label=='dedicated_spatial' else label
        with np.load(args.history/(history_label+'.npz')) as h, np.load(args.geometry/(label+'.npz')) as f:
            results[label]={key:bounds(h,f,meta['runs'][history_label]['markers'][key]) for key in ('0.001','0.01')}
    for label in ('dedicated','shared'):
        for marker in ('0.001','0.01'):
            native,half=[results[label+'_'+grid][marker] for grid in ('native','half')]
            for rule in half['bounds']:
                ratio=native['bounds'][rule]['minimum_launch_rms_to_avoid_excess']/half['bounds'][rule]['minimum_launch_rms_to_avoid_excess']
                if abs(ratio-1.) > .02: raise ValueError('Unresolved native time step')
    args.output.write_text(json.dumps(results,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
