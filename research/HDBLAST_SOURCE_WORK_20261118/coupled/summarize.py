"""Summarize paired native-shell evolutions on a common proper-time cut."""
import argparse
import csv
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np

COLS = ('T','H0tau','A_b','phi_b','H_over_H0','v_over_H0','lapse','rho_hat','p_hat','j_hat','Q_hat','rad_hat',
        'comoving_number','comoving_energy','scalar_work','pressure_work','max_physical_frequency',
        'min_uncertainty','min_occupation','rad','chi','X','R','Hshell','Mshell','Hmax_near','Mmax_near')
BASES = ('quiet','plus','minus','far')
B=6.49e-6
RB=7.838285538073243


def compact(work, target, nodes=24, prep=.70):
    """Store wall observables independently of mutable solver caches."""
    cuts={}
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED) as out:
        for grid in ('coarse','fine'):
            for y in (0,1):
                for feed in (0,1):
                    name=f'{grid}/y{y}_b{feed}'
                    folder=work/grid/f'k{nodes}_y{y}_b{feed}_t{prep:.3f}'
                    prepare=json.loads((folder/'prepare.json').read_text())
                    pre=np.load(str(folder/'fork.npz')+'.history.npz',allow_pickle=False)
                    precols=dict(zip(pre['cols'],pre['rec'].T))
                    def numerical(data, k):
                        key = k if k!='comoving_number' or k in data else 'n_hat'
                        if key not in data and key in ('Hshell','Mshell','Hmax_near','Mmax_near'):
                            return np.full(len(data['T']),np.nan)
                        return np.asarray(data[key],dtype=float)
                    def store_history(path, data):
                        buf=io.StringIO()
                        writer=csv.writer(buf,lineterminator='\n')
                        writer.writerow(COLS)
                        for row in np.array([numerical(data,k) for k in COLS]).T:
                            writer.writerow([format(value,'.17g') for value in row])
                        out.writestr(path,buf.getvalue())
                    store_history(name+'/pre.csv',precols)
                    for branch in BASES:
                        path=folder/branch
                        with np.load(path/'wall_timeseries.npz',allow_pickle=False) as f:
                            d=dict(zip(f['cols'],f['rec'].T))
                        store_history(name+'/'+branch+'.csv',d)
                        if branch!='quiet':
                            cuts[name+'/'+branch]=json.loads((path/'match.json').read_text())
                    cuts[name+'/prepare']=prepare
        out.writestr('initial_data.json',json.dumps(cuts,indent=2)+'\n')


def read(archive, name):
    with archive.open(name+'.csv') as f:
        m=np.genfromtxt(f,delimiter=',',names=True,dtype=float)
    return {k:m[k] for k in m.dtype.names}


def summarize(archive, cut):
    cases={}
    init=json.loads(archive.read('initial_data.json'))
    for grid in ('coarse','fine'):
        for y in (0,1):
            for b in (0,1):
                prefix=f'{grid}/y{y}_b{b}'
                for arm in BASES:
                    name=prefix+'/'+arm
                    pre=read(archive,prefix+'/pre')
                    hist=read(archive,name)
                    d={k:np.r_[pre[k],hist[k]] for k in COLS}
                    s=d['H0tau']
                    if min(np.diff(s))<=0 or not s[0] <= cut <= s[-1]:
                        raise ValueError(name+' bad observer clock')
                    def at(k):return float(np.interp(cut,s,d[k]))
                    vals={k:at(k) for k in ('T','phi_b','H_over_H0','rad_hat','rho_hat','p_hat','j_hat',
                           'Q_hat','comoving_number','comoving_energy','scalar_work','pressure_work')}
                    positive=np.flatnonzero(d['phi_b']>=.5)
                    if not len(positive) or positive[0]==0:
                        raise ValueError(name+' no rising crossing')
                    i=positive[0]
                    vals['crossing_s']=float(np.interp(.5,d['phi_b'][i-1:i+1],s[i-1:i+1]))
                    mask=(s<=cut)
                    vals['max_physical_frequency_through_cut']=float(max(max(d['max_physical_frequency'][mask]),at('max_physical_frequency')))
                    vals['min_uncertainty_through_cut']=float(min(d['min_uncertainty'][mask]))
                    vals['min_occupation_through_cut']=float(min(d['min_occupation'][mask]))
                    denom=np.maximum(1,np.abs(d['scalar_work']))
                    residual=d['comoving_energy']-d['scalar_work']+d['pressure_work']
                    vals['max_scaled_ward_error']=float(np.max(abs(residual[mask])/denom[mask]))
                    vals['ward_error_at_cut']=vals['comoving_energy']-vals['scalar_work']+vals['pressure_work']
                    vals['recorded_constraint_max_postfork_through_cut']={k:float(np.nanmax(d[k][mask & (d['T'] >= 1.11-1e-10)])) for k in ('Hshell','Mshell','Hmax_near','Mmax_near')}
                    h=d['H_over_H0']
                    fraction=d['rad']/h**2
                    vals['rad_fraction']=float(np.interp(cut,s,fraction))
                    vals['max_rad_fraction_through_cut']=float(max(max(fraction[mask]),vals['rad_fraction']))
                    # B1 puts the mixed RX term in its massive-matter contribution.
                    reassign=fraction+RB**2*d['R']*d['X']/(18*h**2)
                    vals['rad_fraction_with_mixed_term']=float(np.interp(cut,s,reassign))
                    vals['restart_phi']=float(d['phi_b'][np.argmin(abs(d['T']-1.11))])
                    vals['restart_s']=float(s[np.argmin(abs(d['T']-1.11))])
                    if vals['max_scaled_ward_error']>2e-5 or vals['min_uncertainty_through_cut']<.248:
                        raise ValueError(name+' failed energy or uncertainty checks')
                    if vals['max_physical_frequency_through_cut']>=B**(-1/3):
                        raise ValueError(name+' violates nominal per-mode screen')
                    cases[name]=vals
    contrasts={}
    for grid in ('coarse','fine'):
        for y in (0,1):
            for arm in BASES:
                if arm=='quiet':continue
                label=f'{grid}/y{y}/{arm}'
                out={}
                for key in ('phi_b','H_over_H0','rho_hat','rad_hat','comoving_number','crossing_s','rad_fraction'):
                    ref=[cases[f'{grid}/y{y}_b{b}/quiet'][key] for b in (0,1)]
                    p=[cases[f'{grid}/y{y}_b{b}/{arm}'][key] for b in (0,1)]
                    diff=[v-q for q,v in zip(ref,p)]
                    out[key]={'passive_delta':diff[0],'coupled_delta':diff[1],'feedback_on_contrast':diff[1]-diff[0],
                              'quiet_feedback':ref[1]-ref[0],
                              'passive_percent':100*diff[0]/ref[0] if ref[0] else None,
                              'coupled_percent':100*diff[1]/ref[1] if ref[1] else None}
                contrasts[label]=out
    return dict(cut_s=cut,B=B,M5_over_H0=B**(-1/3),cases=cases,contrasts=contrasts,preparations={k:v for k,v in init.items() if k.endswith('/prepare')})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build-from',type=Path)
    p.add_argument('--archive',type=Path,default=Path(__file__).with_name('WALLS.zip'))
    p.add_argument('--cut',type=float,default=2.10)
    args=p.parse_args()
    if args.build_from:
        compact(args.build_from,args.archive)
    with zipfile.ZipFile(args.archive) as archive:
        output=summarize(archive,args.cut)
    print(json.dumps(output,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
