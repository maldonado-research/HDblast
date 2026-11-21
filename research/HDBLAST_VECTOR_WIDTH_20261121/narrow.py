"""Rerun the crossed-stress generator on one frozen snapshot at narrower timelike widths."""
import argparse
import hashlib
import importlib.util
import json
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


def load(path):
    spec=importlib.util.spec_from_file_location(path.stem,path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def stats(values):
    values=np.asarray(values)
    return dict(mean=float(values.mean()),batch_se=float(values.std(ddof=1)/np.sqrt(len(values))),batches=values.tolist())


def census(first,second,cut,save):
    ea=np.sqrt(1+(first*first).sum(axis=1));eb=np.sqrt(1+(second*second).sum(axis=1))
    ia=np.flatnonzero(4*ea*eb.max()>=cut);ib=np.flatnonzero(4*eb*ea.max()>=cut)
    count=0; saved=[]
    for start in range(0,len(ia),96):
        ind=ia[start:start+96]
        s=2*(1+ea[ind,None]*eb[ib]-first[ind]@second[ib].T)
        rows,cols=np.nonzero(s>=cut)
        count+=len(rows)
        if save and len(rows):
            saved.append(np.column_stack((ind[rows],ib[cols])))
    return count,np.vstack(saved) if saved else np.zeros((0,2),dtype=int)


def enumerate_tail(tags,pair_indices,generator,kernel,normalization,nodes,azimuths):
    points,weights=np.polynomial.legendre.leggauss(nodes)
    result=np.zeros((2,3))
    for sign,(first,second) in enumerate(tags):
        indices=pair_indices[sign]
        subtotal=np.zeros((2,3))
        for lo in range(0,len(indices),256):
            take=indices[lo:lo+256];n=len(take)
            point=(np.arange(n)+.5)/n
            phase=((lo+np.arange(n)+1)*np.sqrt(2.))%1
            # The center of the index bin selects each explicitly enumerated pair exactly once.
            u=np.column_stack([point,point,phase,np.full(n,.5),np.full(n,.5)])
            values,_=generator.integrate_pairs(first[take[:,0]],second[take[:,1]],u,(points+1)/2,weights/2,azimuths,sign!=0,kernel)
            if np.any(values[0]!=0): raise ArithmeticError('Inconsistent upper-tail classification')
            subtotal+=values[1]
        result+=(1. if sign==0 else .5)*subtotal/(len(first)*len(second))
    return dict(pole=float(result[0,1]*normalization),crossed=float(result[1,1]*normalization))


def run(args):
    generator=load(args.study/'crossed-stress/generator.py')
    partial=load(args.study/'unitarity/partial_wave.py')
    with np.load(args.snapshot) as f:
        p=f['p']; sign=f['sign']
    with np.load(args.history) as f: history=dict(f)
    plus=p[sign==1];minus=p[sign==-1]
    if len(plus)!=16384 or len(minus)!=16384: raise ValueError('Use the finest t=7.05 snapshot')
    if not args.snapshot.name.startswith('t7.05000-p14-n8192-r91283-particle-pole'): raise ValueError('Expected dedicated particle-current pulse')
    tags=((plus,minus),(plus,plus),(minus,minus)); avg=np.sqrt(1.+(p*p).sum(axis=1))
    x,w=np.polynomial.legendre.leggauss(10)
    rows=[]
    for width in (.065,.043,.035):
        mod=load(args.study/'crossed-exploration/check.py'); mod.WIDTH=width;mod.G2=12*np.pi*width
        cut=brentq(lambda s:float(partial.partial_magnitude(s))*width/.30-.5,mod.M**2,1e14)
        kernel=vars(mod);kernel['unitarity_cut']=np.sqrt(cut)/mod.M
        generator.analytic_check(kernel)
        counts=[];saved=[];exact=width<=.043
        for first,second in tags:
            cnt,tail=census(first,second,cut,exact);counts.append(cnt);saved.append(tail)
        measure=generator.measurement(args.snapshot,history,kernel,(x+1)/2,w/2,8,args.power,args.batches)
        norm=measure['normalized_target_density']/(2*measure['mean_gamma'])
        full_tail={}
        if exact:
            for n,az in ((10,8),(14,16),(18,32)):
                full_tail[f'{n}x{az}']=enumerate_tail(tags,saved,generator,kernel,norm,n,az)
        summarized={}
        for label in ('pole','crossed'):
            low=np.array(measure['below_cut'][label+'_derivative']['batches'])
            sampled_high=np.array(measure['above_or_equal_cut'][label+'_derivative']['batches'])
            high=np.full(len(low),full_tail['18x32'][label]) if exact else sampled_high
            summarized[label]=dict(below=stats(low),above=stats(high),total=stats(low+high),
                                   sampled_above=stats(sampled_high))
        pl=np.array(summarized['pole']['below']['batches']);xl=np.array(summarized['crossed']['below']['batches'])
        pa=np.array(summarized['pole']['above']['batches']);xa=np.array(summarized['crossed']['above']['batches'])
        rows.append(dict(width=width,g_squared=mod.G2,unitarity_s=cut,sqrt_s_over_M=kernel['unitarity_cut'],
                         above_ordered_pair_counts=counts,pairs_per_charge_product=len(plus)**2,
                         quadrature_above_exhaustive=full_tail,derivatives=summarized,
                         full_minus_pole=dict(below=stats(xl-pl),above=stats(xa-pa),total=stats(xl+xa-pl-pa)),
                         full_rate_per_wall_time=measure['crossed_born']['encounters_per_particle_per_wall_time'],
                         pole_rate_per_wall_time=measure['isolated_pole']['encounters_per_particle_per_wall_time']))
        print('Finished width',width,'Qdot',summarized['crossed']['total']['mean'],flush=True)
    return dict(snapshot=args.snapshot.stem,momenta_sha256=hashlib.sha256(np.asarray(p,dtype='<f8').tobytes()).hexdigest(),
                initial_Q=measure['initial_wall_shear_over_energy'],gamma_max=float(avg.max()),all_pair_s_ceiling=float(4*avg.max()**2),
                per_channel_draws_per_batch=2**args.power,batches=args.batches,results=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ('study','snapshot','history','output'):
        parser.add_argument('--'+arg,type=Path,required=True)
    parser.add_argument('--power',type=int,default=15)
    parser.add_argument('--batches',type=int,default=8)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='narrow-vector-') as td:
        with np.load(args.snapshot) as source:
            if 'p_final' in source:
                compact=Path(td)/args.snapshot.name
                momenta=source['p_final']
                np.savez(compact,p=momenta,sign=np.tile([1.,-1.],len(momenta)//2))
                args.snapshot=compact
        payload=run(args)
    args.output.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
