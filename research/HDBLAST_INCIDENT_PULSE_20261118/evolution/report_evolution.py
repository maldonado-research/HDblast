"""Recompute the evolution comparisons from recorded wall histories in standard Python."""
import argparse
import csv
import hashlib
import json
import math
import zipfile
from pathlib import Path


def read_series(path):
    with path.open() as f:
        rows=list(csv.DictReader(f))
    return {key:[float(row[key]) for row in rows] for key in rows[0]}


def integrate(x,y):
    return math.fsum((b-a)*(c+d)/2 for a,b,c,d in zip(x,x[1:],y,y[1:]))


def subtract(x,y):
    return [a-b for a,b in zip(x,y,strict=True)]


def maximum(values):
    return max(abs(x) for x in values)


def norm_on(T,f,stop,start=0):
    return maximum(value for t,value in zip(T,f) if start-1e-12<=t<=stop+1e-12)


def nearest(T,f,t):
    j=round(t/(T[1]-T[0]))
    if abs(T[j]-t)>1e-12:
        raise ValueError('Expected common diagnostic node')
    return f[j]


def summarize(root):
    data={}
    histories={}
    for label in ('g1','g2','g4','g8'):
        where=root/label
        meta=json.loads((where/'metadata.json').read_text())
        case={}
        for name,entry in meta['cases'].items():
            path=where/(name+'.csv')
            if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']:
                raise ValueError(path)
            case[name]=read_series(path)
        mass_dir=root.parent/'mass_runs'/label
        mass_meta=json.loads((mass_dir/'metadata.json').read_text())
        for source,name in [('quiet','mass_quiet'),('mass_matched','mass_matched')]:
            path=mass_dir/(source+'.csv')
            if hashlib.sha256(path.read_bytes()).hexdigest()!=mass_meta['cases'][source]['sha256']:
                raise ValueError(path)
            case[name]=read_series(path)
        quiet=case['quiet']
        T=quiet['T']
        histories[label]=case
        stats={}
        for name,series in case.items():
            if series['T']!=T:
                raise ValueError('Unmatched clocks')
            reference=case['mass_quiet'] if name.startswith('mass_') else quiet
            diff=subtract(series['phi_wall'],reference['phi_wall'])
            peak_index=max(range(len(diff)),key=lambda i:abs(diff[i]))
            profile={
                'initial_cut_mismatch_B_Q_u':meta['cases'].get(name,{}).get('initial_cut_mismatch_B_Q_u'),
                'peak_difference_from_quiet':diff[peak_index],
                'T_at_peak_node':T[peak_index],
                'last_difference_from_quiet':diff[-1],
                'maximum_absolute_shell_scalar':maximum(series['phi_wall']),
                'maximum_relative_H_constraint':max(series['H_constraint_rel']),
                'maximum_relative_M_constraint':max(series['M_constraint_rel']),
                'peak_quiet_difference_by_T':{str(t):norm_on(T,diff,t) for t in (.005,.05,.08,.095,.1)},
                'delta_phi_wall_samples':{str(t):nearest(T,diff,t) for t in (.10,.105,.11,.12,.13,.14,.15,.16,.18)},
            }
            stats[name]=profile
        local_pair=subtract(case['delayed']['phi_wall'],case['uncompensated']['phi_wall'])
        negative_pair=subtract(case['delayed_negative']['phi_wall'],case['uncompensated_negative']['phi_wall'])
        positive=subtract(case['delayed']['phi_wall'],quiet['phi_wall'])
        negative=subtract(case['delayed_negative']['phi_wall'],quiet['phi_wall'])
        even=[(a+b)/2 for a,b in zip(positive,negative)]
        odd=[(a-b)/2 for a,b in zip(positive,negative)]
        quiet_drift=subtract(quiet['phi_wall'],[quiet['phi_wall'][0]]*len(T))
        data[label]={
            'dz':meta['spacing'],'dt':meta['dt'],'series_count':len(case),'cases':stats,
            'quiet_max_abs_scalar_drift':maximum(quiet_drift),
            'matched_sign_even_to_odd_suprema_0100_0150':norm_on(T,even,.15,.1)/norm_on(T,odd,.15,.1),
            'negative_delayed_minus_uncompensated_max_by_T':{str(t):norm_on(T,negative_pair,t) for t in (.1,.13,.14,.145,.15,.152,.16,.18)},
            'delayed_minus_uncompensated_max_by_T':{str(t):norm_on(T,local_pair,t) for t in (.1,.13,.14,.145,.15,.152,.16,.18)},
        }
    halves={}
    for label in ('g4','g8'):
        folder=root.parent/'half_runs'/label
        half_meta=json.loads((folder/'metadata.json').read_text())
        D={}
        for name,entry in half_meta['cases'].items():
            path=folder/(name+'.csv')
            if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']:
                raise ValueError(path)
            D[name]=read_series(path)
        T=D['quiet']['T']
        if T!=histories[label]['quiet']['T'] or D['quiet']['phi_wall']!=histories[label]['quiet']['phi_wall']:
            raise ValueError('Half-amplitude control differs')
        half={};curves={}
        for sign,full in [('positive','delayed'),('negative','delayed_negative')]:
            field=D['delayed_'+sign]['phi_wall']
            curves[sign]=subtract(field,D['quiet']['phi_wall'])
            full_curve=subtract(histories[label][full]['phi_wall'],D['quiet']['phi_wall'])
            half[sign]={
               'peak_absolute_delta_phi_before_controls':norm_on(T,curves[sign],.15),
               'max_relative_amplitude_doubling_discrepancy':norm_on(T,[2*h-f for h,f in zip(curves[sign],full_curve)],.15)/norm_on(T,full_curve,.15),
               'matched_minus_uncompensated_supremum_to_015':norm_on(T,subtract(field,D['uncompensated_'+sign]['phi_wall']),.15),
               'initial_cut_mismatch_B_Q_u':half_meta['cases']['delayed_'+sign]['initial_cut_mismatch_B_Q_u']}
        even=[(a+b)/2 for a,b in zip(curves['positive'],curves['negative'])]
        odd=[(a-b)/2 for a,b in zip(curves['positive'],curves['negative'])]
        half['even_to_odd_suprema_0100_0150']=norm_on(T,even,.15,.1)/norm_on(T,odd,.15,.1)
        halves[label]=half
    comparisons={}
    for name in histories['g1']:
        grid_errors=[]
        for low,high in zip(('g1','g2','g4'),('g2','g4','g8')):
            coarse,fine=histories[low], histories[high]
            peak_responses=[]
            for d in (coarse,fine):
                reference=d['mass_quiet'] if name.startswith('mass_') else d['quiet']
                peak_responses.append(subtract(d[name]['phi_wall'],reference['phi_wall']))
            a,b=peak_responses
            if len(a)!=len(b[::2]):
                raise ValueError('Unmatched sampling cadence')
            gap=subtract(a,b[::2])
            grid_errors.append(maximum(gap))
        comparisons[name]={'supremum_successive_grid_delta_phi':grid_errors,
                           'error_reduction_factors':[a/b if b else None for a,b in zip(grid_errors,grid_errors[1:])]}
    return {'time_window':[0,.18],'grids':data,'wall_response_convergence':comparisons,'half_amplitudes':halves}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path,help='wall_histories.zip or extracted runs directory')
    args=parser.parse_args()
    if args.source.is_dir():
        result=summarize(args.source)
    else:
        with zipfile.ZipFile(args.source) as archive:
            result=summarize(zipfile.Path(archive)/'runs')
    print(json.dumps(result,indent=2,allow_nan=False))
