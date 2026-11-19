"""Consume generated tables: Ward/closure checks and held-out high-band errors."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import cumulative_simpson,simpson

sys.path.insert(0,str(Path(__file__).resolve().parent))
from source import SourceTable,read_native
from quadrature import make_bath,high_moments,LOW,HIGH

HERE = Path(__file__).resolve().parent


def verify(tables,native,output):
    tables,native,output = map(Path,(tables,native,output))
    output.mkdir(parents=True,exist_ok=True)
    report = dict(schema=1,units='hatted physical coefficients per y²',tables={},comparisons={},rules={})
    manifest = dict(schema=1,approximation='fixed-trajectory order-y² positive high-band bath',
        native_inputs={name:hashlib.sha256((native/name).read_bytes()).hexdigest() for name in ('FORK.json','WALLS.zip')},runs={})
    for path in sorted(tables.glob('*.json')):
        if path.name=='SOURCES.json':continue
        row = json.loads(path.read_bytes())
        if all(k in row for k in ('grid','arm','bath','files')):
            row.pop('seconds',None)
            manifest['runs'][f"{row['grid']}_{row['arm']}_{row['bath']}"] = row
    paths = sorted(tables.glob('*.npz'))
    assert paths
    for path in paths:
        grid,arm,bath,scheme = path.stem.split('_')
        q,weights,_,(_,loga,mass) = read_native(native,grid,arm)
        metadata = manifest['runs'][f'{grid}_{arm}_{bath}']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==metadata['files'][scheme]['sha256']
        with np.load(path) as z:
            s,a,h,v = (z[k] for k in ('s','a','h','v'))
            ward = z['drho']+3*h*(z['rho']+z['p'])-z['j']*v
            work = cumulative_simpson(a**3*(z['j']*v-3*h*z['p']),x=s,initial=0)
            work_error = float(np.max(abs(z['energy_comoving']-work)))
            component_error = max(float(np.max(abs(z[key]-sum(z[key+'_'+c] for c in ('parent','daughter','interaction','counterterm'))))) for key in ('rho','p','j','drho','dp','dj'))
            post = s>=metadata['gate_off']+1e-10
            daughter_conformal = z['rad'][post]*a[post]**4
            redshift_error = float(np.max(abs(daughter_conformal-daughter_conformal[0]))/max(abs(daughter_conformal[0]),1e-30))
            closed_interaction = max(float(np.max(abs(z[key+'_'+c][post]))) for key in ('rho','p','j') for c in ('interaction','counterterm'))
            ct_pressure_error = float(np.max(abs(z['p_counterterm']+z['rho_counterterm']/3)))
            daughter_pressure_error = float(np.max(abs(z['p_daughter']-z['rho_daughter']/3)))
            minimum_daughter = float(np.min(z['rad']))
            daughter_scale = max(float(np.max(z['rad'])),1.)
            assert minimum_daughter>=-1e-12*daughter_scale
            assert component_error==0 and closed_interaction==0
            assert redshift_error<2e-12 and ct_pressure_error<1e-11 and daughter_pressure_error<1e-11
            assert np.max(abs(ward))<2e-8 and work_error<2e-6
        table = SourceTable(path)
        mid=(s[:-1]+s[1:])/2
        val,deriv=table.evaluate(mid),table.evaluate(mid,1)
        midpoint_ward=deriv[:,0]+3*loga(mid,1)*(val[:,0]+val[:,1])-val[:,2]*mass(mid,1)/100
        assert np.max(abs(midpoint_ward))<2e-5
        report['tables'][path.stem]=dict(max_abs_nodal_ward=float(np.max(abs(ward))),
            max_abs_interpolated_ward=float(np.max(abs(midpoint_ward))),max_abs_integrated_comoving_work_error=work_error,
            component_sum_error=component_error,after_closure_interaction_and_ct=closed_interaction,
            after_closure_conformal_daughter_energy_relative_drift=redshift_error,
            static_pressure_identity_error=ct_pressure_error,daughter_pressure_identity_error=daughter_pressure_error,
            minimum_daughter_density=minimum_daughter)
    cases=sorted(set('_'.join(path.stem.split('_')[:2]) for path in paths))
    rows=[]
    fields=('rho','p','j','rad')
    for case in cases:
        for scheme in ('bare','static'):
            with np.load(tables/f'{case}_N64_{scheme}.npz') as z:
                s=z['s'];a=z['a'];ref=np.array([z[k] for k in fields]);index=int(np.argmin(abs(s-2.1)))
            with np.load(tables/f'{case}_low_{scheme}.npz') as z:
                np.testing.assert_array_equal(s,z['s']);low=np.array([z[k] for k in fields])
            high=ref-low
            denominator=np.maximum(np.max(abs(high),axis=1),1e-30)
            l1den=np.maximum(simpson(abs(high),x=s,axis=1),1e-30)
            family={}
            for path in sorted(tables.glob(f'{case}_N*_{scheme}.npz')):
                bath=path.stem.split('_')[2]
                with np.load(path) as z:
                    np.testing.assert_array_equal(s,z['s']);test=np.array([z[k] for k in fields])
                error=test-ref
                linf=np.max(abs(error),axis=1)/denominator
                l1=simpson(abs(error),x=s,axis=1)/l1den
                endpoint=(test[:,index]-ref[:,index])/np.maximum(abs(high[:,index]),1e-30)
                family[bath]=dict(high_band_normalized_linf=dict(zip(fields,map(float,linf))),
                    high_band_normalized_l1=dict(zip(fields,map(float,l1))),
                    readout_high_band_relative_error=dict(zip(fields,map(float,endpoint))),
                    reference_high_daughter_comoving=float(high[3,index]*a[index]**3),
                    high_daughter_comoving=float((test[3,index]-low[3,index])*a[index]**3))
                for fi,field in enumerate(fields):
                    rows.append([case,scheme,bath,field,linf[fi],l1[fi],endpoint[fi]])
            report['comparisons'][f'{case}_{scheme}']=family
    q=read_native(native,'fine','quiet')[0]
    for order in (1,2,4,8,16,32,64):
        omega,g2,method=make_bath(q,order)
        highomega,highg2=omega[:,64:],g2[:,64:]
        c,d=np.array([high_moments(k) for k in q]).T
        moment=float(np.max(abs(np.array([np.sum(highg2/highomega**2,axis=1)/c-1,np.sum(highg2/highomega**4,axis=1)/d-1]))))
        assert np.min(highg2)>0 and moment<3e-13
        report['rules'][str(order)]=dict(method=method,min_strength=float(np.min(highg2)),max_relative_contact_and_kinetic_moment_error=moment)
    with (output/'COMPARISON.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(['case','scheme','bath','field','high_normalized_linf','high_normalized_l1','readout_high_relative_error']);writer.writerows(rows)
    (output/'CHECKS.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (output/'SOURCES.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(f"Verified {len(paths)} tables, {len(cases)} wall cases")
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables',type=Path,required=True)
    parser.add_argument('--native',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=HERE)
    args=parser.parse_args()
    verify(args.tables,args.native,args.output)
