"""Measure which launch waveforms remain independent in the cohort electric impulse."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import kv


def cohort_operator(h, waves, marker):
    """Mid-cell cohorts, trapezoidal electromagnetic impulse and the reference axial gamma."""
    end=marker['time']; stop=h['t'].searchsorted(end)
    t=np.r_[h['t'][:stop],end]
    def knots(arr):
        return np.r_[arr[:stop],np.interp(end,h['t'],arr)]
    a=knots(h['a']); mass=np.diff(knots(h['D']))/marker['D']
    if not np.isclose(mass.sum(),1.) or np.min(mass)<-1e-13:raise ValueError('Birth measure')
    def operator(f):
        cells=.5*np.diff(t)*(a[:-1]*f[:-1]+a[1:]*f[1:])
        suffix=np.r_[np.cumsum(cells[::-1])[::-1],0.]
        return .5*(suffix[1:]+suffix[:-1])*1e7/marker['a']
    field=np.array([operator(knots(row)) for row in waves['E']]).T
    axial=operator(knots(h['X']))
    weight=mass/np.hypot(1.,axial)
    return field,mass,weight


def optimize(J, metric, weights):
    """Best equal two-dimensional response, per unit total launch energy."""
    L=np.linalg.cholesky(metric)
    F=solve_triangular(L,J.T,lower=True).T*np.sqrt(weights)[:,None]
    _,s,right=np.linalg.svd(F,full_matrices=False)
    s1,s2=map(float,s[:2])
    if not s2>0 or not s1>0: raise ValueError('Launch map has rank below two')
    best=s1*s2/math.hypot(s1,s2)
    mix=solve_triangular(L.T,right[:2].T,lower=False).T*(best/np.array([s1,s2]))[:,None]
    for row in mix:
        if row[np.argmax(abs(row))]<0:row*=-1.
    measured=J@mix.T
    measured=measured.T@(weights[:,None]*measured)
    defect=float(np.max(abs(measured-best**2*np.eye(2)))/best**2)
    if not np.isclose(np.trace(mix@metric@mix.T),1.,rtol=1e-10) or defect>3e-5:
        raise ValueError(('Optimal two-polarization balance',defect,s))
    return s,mix,best


def localized_cdf(h, waves, marker, source, names, bins, centers=1):
    """Causal upper support in each birth-mass and radial-volume cell of S3."""
    end=marker['time']; j=h['t'].searchsorted(end); t=np.r_[h['t'][:j],end]
    def nodes(f):return np.r_[f[:j],np.interp(end,h['t'],f)]
    a=nodes(h['a']); d=nodes(h['D']); dt=np.diff(t); wt=np.diff(d)/marker['D']
    if np.any(np.diff(a)<0.) or np.min(wt)<0.:raise ValueError('Expansion or births')
    quantile=np.arange(1,bins+1,dtype=float)/bins
    # Each rigidly rotated harmonic has P(|z1|²>u)=1-u. A union over
    # at most M centers bounds P(max |z1|²>u) by M(1-u), for any placement.
    radius_squared=1. - (1. - quantile)/centers
    # Distance to the maximum-amplitude Hopf circle is arccos |z1|.
    travel=float(np.sum(dt/a[:-1])); r=np.cos(np.maximum(0.,np.arccos(np.sqrt(radius_squared))-travel))
    dual_e=[]; dual_b=[]; shape=[]; selected=sorted(k for k in np.unique(waves['kappa']) if k>=40)
    for k in selected:
        ix=np.flatnonzero(waves['kappa']==k);ids=[names.index(label) for label in waves['profile'][ix]]
        L=np.linalg.cholesky(source[np.ix_(ids,ids)])
        shape.append((k-1.)*r**(2*(k-2)))
        for array,terms in [('E',dual_e),('B',dual_b)]:
            values=np.array([nodes(row) for row in waves[array][ix]])
            white=solve_triangular(L,values,lower=True)
            p=np.sum(white*white,axis=0);terms.append(np.maximum(p[1:],p[:-1]))
    force=np.sqrt(np.array(shape).T@np.asarray(dual_e))+np.sqrt(np.array(shape).T@np.asarray(dual_b))
    contributions=(1.e7/marker['a'])*(dt*a[1:])[None,:]*force
    caps=np.cumsum(contributions[:,::-1],axis=1)[:,::-1].ravel()
    weights=np.tile(wt/bins,bins)
    order=np.argsort(caps);cum=np.cumsum(weights[order])
    cap=caps[order]; pstar=np.sqrt(np.max(h['R'][:j+1])**2-1.)
    theta=marker['theta']; rows={}
    for eps in (.02,.1,.5,1.):
        cutoff=pstar+eps*cap
        z=np.hypot(1.,cutoff)
        thermal=np.clip(1.-z*z*kv(2.,z/theta)/kv(2.,1./theta),0.,1.)
        witness=int(np.argmax(cum-thermal))
        rows[str(eps)]=dict(cutoff_in_mc=float(cutoff[witness]),
                            guaranteed_mass_below_cutoff=float(cum[witness]),
                            thermal_mass_below_cutoff=float(thermal[witness]),
                            guaranteed_transverse_cdf_excess=float(cum[witness]-thermal[witness]))
    return dict(radius_volume_bins=bins,kappas=list(map(int,selected)),maximum_angular_travel=travel,
                launch_amplitudes=rows)


def case(base,history,label,table,data_dir,include_optimal):
    name='dedicated_native' if label=='dedicated_spatial' else label
    marker=base['runs'][name]['markers']['0.01']
    with np.load(history/(name+'.npz')) as h, np.load(data_dir/'traces'/(label+'.npz')) as waves:
        if not np.array_equal(h['t'],waves['t']):raise ValueError('Field and birth clocks differ')
        J,mass,stress_weight=cohort_operator(h,waves,marker)
        drift=np.interp(marker['time'],h['t'],h['u']); target=marker['theta']/np.hypot(1.,drift)
        source=np.asarray(table['emitted_energy_gram']); names=table['profile_order']
        metrics={}; axes=[]; coeffs=[]; summed_E=[]; summed_B=[]
        for k in np.unique(waves['kappa']):
            ix=np.flatnonzero(waves['kappa']==k)
            if len(ix)<2:continue
            labels=waves['profile'][ix].tolist()
            ids=[names.index(label) for label in labels]
            G=source[np.ix_(ids,ids)]
            s,mixture,best=optimize(J[:,ix],G,stress_weight)
            # The leading and subleading directions of the unweighted momentum second moment.
            mom,_,gain=optimize(J[:,ix],G,mass)
            row=dict(profiles=labels,
                direct_pressure_singular_values_first_two=s[:2].tolist(),
                second_over_first_pressure_amplitude=float(s[1]/s[0]),
                maximum_equal_lateral_pressure_per_unit_launch_energy=best**2,
                equal_pressure_target_launch_amplitude=math.sqrt(target)/best,
                weakest_rms_impulse_max_at_unit_energy=gain)
            if include_optimal:
                row['unit_launch_coefficients_for_two_orthogonal_polarizations']=mixture.tolist()
                for q in range(2):
                    axes.append(k); coeffs.append(q)
                    summed_E.append(mixture[q]@waves['E'][ix]);summed_B.append(mixture[q]@waves['B'][ix])
            if not include_optimal and s[1]/s[0]<1.e-6:
                # Avoid attributing physical precision to nearly parallel templates.
                row=dict(profiles=labels,leading_direct_pressure_singular_value=float(s[0]),
                         second_over_first_below_one_part_per_million=True)
            metrics[str(k)]=row
        result=dict(target_lateral_pressure_per_lab_particle=float(target),modes=metrics)
        if not include_optimal:
            result['localized_high_mode_cdf']=localized_cdf(h,waves,marker,source,names,256)
            result['rotated_center_gaps']={}
            for centers in (2,3,10):
                rotated=localized_cdf(h,waves,marker,source,names,256,centers)
                result['rotated_center_gaps'][str(centers)]={
                    key:max(0.,row['guaranteed_transverse_cdf_excess'])
                    for key,row in rotated['launch_amplitudes'].items()}
            coarse=localized_cdf(h,waves,marker,source,names,128)
            result['coarse_localized_gaps']={key:row['guaranteed_transverse_cdf_excess'] for key,row in coarse['launch_amplitudes'].items()}
        fields=dict(t=waves['t'],kappa=np.array(axes),axis=np.array(coeffs),E=summed_E,B=summed_B) if axes else None
    return result,fields


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history',type=Path,required=True)
    p.add_argument('--fields',type=Path,required=True)
    p.add_argument('--scan',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();out=args.output_dir
    (out/'optimal_traces').mkdir(parents=True,exist_ok=True)
    meta=json.loads((args.history/'markers.json').read_text()); results={}
    for typ,folder,product in [('primary',args.fields,'FIELDS.json'),('preparation_scan',args.scan,'STORAGE_SCAN.json')]:
        inputs=json.loads((folder/product).read_bytes());results[typ]={}
        for label,table in inputs.items():
            row,fields=case(meta,args.history,label,table,folder,typ=='primary')
            if fields: np.savez_compressed(out/'optimal_traces'/(label+'.npz'),**fields)
            results[typ][label]=row
    (out/'IMPULSE.json').write_text(json.dumps(results,indent=2,allow_nan=False,sort_keys=True)+'\n')

if __name__=='__main__': main()
