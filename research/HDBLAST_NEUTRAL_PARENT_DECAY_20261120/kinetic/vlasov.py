"""Self-consistent coherent Maxwell response of a collisionless, continuously born gas."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent


def background(path,linear=False):
    data=json.loads(path.read_bytes())
    cfg=data['parameters']['primary']
    run=data['runs']['primary']
    src=dict(zip(data['trace_columns'],np.asarray(run['trace']).T))
    start_i=np.flatnonzero((src['psi'][:-1]<cfg['psi_star'])&(src['psi'][1:]>=cfg['psi_star']))[0]
    start=float(src['t'][start_i])
    end=run['maximum_tau']
    assert src['f'][start_i]==0 and src['R'][start_i]<=1 and end>start
    K=src['K'][start_i]
    A=src['a'][start_i]**2*src['B'][start_i]/2
    P=src['a'][start_i]**2*src['X'][start_i]
    ix=(src['t']>=start)&(src['t']<=end+.2)
    t=src['t'][ix]
    if np.ptp(src['K'][ix])>1e-7:raise ValueError('K changes on integration interval')
    psi=CubicHermiteSpline(t,src['psi'][ix],src['pi'][ix])
    opening=brentq(lambda tau:psi(tau)-cfg['psi_star'],t[0],t[1])
    tt=np.insert(t,1,opening)
    frac=PchipInterpolator(tt,np.insert(src['f'][ix],1,0.))
    mass=PchipInterpolator(tt,np.insert(src['R'][ix],1,1.))
    radius=CubicHermiteSpline(t,src['a'][ix],(src['a']*src['H'])[ix])
    potential=PchipInterpolator(t,src['V'][ix])
    def sample(times):
        density=.5*psi(times,1)**2+potential(times)
        if linear:
            return np.c_[times,radius(times),np.interp(times,tt,np.insert(src['f'][ix],1,0.)),
                         np.interp(times,tt,np.insert(src['R'][ix],1,1.)),density]
        return np.c_[times,radius(times),frac(times),mass(times),density]
    return dict(start=start,end=end,A=A,P=P,K=K,D=8*cfg['D']/cfg['R_low'],Q=cfg['Q'],
                opening=opening,sample=sample,
                parent_bulk_integral_change=float(np.interp(end,t,src['comoving_wave_to_EM'][ix])-src['comoving_wave_to_EM'][start_i]),
                parent_scalar_transfer_change=float(np.interp(end,t,src['comoving_scalar_to_EM'][ix])-src['comoving_scalar_to_EM'][start_i]))


def spectrum(shells,z,Q,angles):
    w,Ab,ashell=shells.T
    x=Q*(z['A']-Ab)/z['a']
    s=ashell/z['a']
    u,gu=np.polynomial.legendre.leggauss(angles)
    weight=w[:,None]*gu[None,:]/(2*w.sum())
    par=x[:,None]+s[:,None]*u
    p2=x[:,None]**2+s[:,None]**2+2*(x*s)[:,None]*u
    gamma=np.sqrt(1+np.maximum(p2,0.))
    v=par/gamma
    mean=lambda b:float(np.sum(weight*b))
    idx=np.argsort(gamma.ravel())
    quantile=np.interp([.1,.5,.9,.99],np.cumsum(weight.ravel()[idx]),gamma.ravel()[idx])
    mean_gamma=mean(gamma)
    t01=mean(par)
    t11=mean(par*par/gamma)
    transverse=mean((s[:,None]**2)*(1-u**2)/(2*gamma))
    disc=math.sqrt((mean_gamma+t11)**2-4*t01*t01)
    beta=2*t01/(mean_gamma+t11+disc)
    rest11=(t11-2*beta*t01+beta**2*mean_gamma)/(1-beta*beta)
    return dict(energy_frame_velocity=beta,parallel_share_of_pressure_in_energy_frame=rest11/(rest11+2*transverse),
                mean_gamma=mean_gamma,mean_signed_velocity=mean(v),mean_absolute_parallel_velocity=mean(abs(v)),
                positive_parallel_fraction=mean(par>0),mean_parallel_momentum=mean(par),
                excess_above_same_mean_momentum=(mean_gamma-math.hypot(1.,mean(par)))/mean_gamma,
                gamma_quantiles_10_50_90_99=quantile.tolist(),
                current_agreement=abs(w.sum()*mean(v)-z['current']))


def run(temp,exe,bg,bins,substeps,Q):
    grid=np.linspace(bg['start'],bg['end'],bins*substeps+1)
    path=temp/'input'
    interval=max(1,(bins*substeps)//512)
    with path.open('w') as fh:
        fh.write(' '.join(map(str,(bins,substeps,interval,Q,bg['K'],bg['D'],bg['A'],bg['P'])))+'\n')
        np.savetxt(fh,bg['sample'](grid),fmt='%.17g')
    subprocess.run([str(exe),str(path),str(temp/'trace'),str(temp/'shells')],check=True)
    tr=np.genfromtxt(temp/'trace',names=True)
    shells=np.atleast_2d(np.loadtxt(temp/'shells'))
    total=tr['field']+tr['daughter']+tr['parent']+tr['scalar']
    frac=tr['field']/total
    norm=tr['field'][0]
    m=tr[-1]
    sp=spectrum(shells,m,Q,16)
    sp8=spectrum(shells,m,Q,8)
    def avg(arr,dt):
        first=np.searchsorted(tr['t'],m['t']-dt)
        t=np.r_[m['t']-dt,tr['t'][first:]]
        v=np.r_[np.interp(t[0],tr['t'],arr),arr[first:]]
        return float(np.trapezoid(v,t)/dt)
    res=dict(bins=bins,substeps=substeps,Q=Q,turn_converted=float(m['f']),initial_comoving_field=float(norm),
             turn_comoving_field=float(m['field']),turn_field_share_of_added=float(frac[-1]),
             turn_daughter_share_of_added=float(m['daughter']/total[-1]),
             turn_net_field_loss_fraction=float(1.-m['field']/norm),
             last_quarter_mean_field_share=avg(frac,.25),last_half_mean_field_share=avg(frac,.5),
             max_normalized_total_balance=float(max(abs(tr['net_residual']))/norm),
             max_normalized_field_balance=float(max(abs(tr['field_residual']))/norm),
             turn_current_work_over_initial_field=float(m['fluxwork']/norm),
             positive_current_work_over_initial_field=float(m['forward_work']/norm),
             returned_current_work_over_initial_field=float(m['returned_work']/norm),
             final_quarter_field_minmax=[float(fn(frac[tr['t']>=m['t']-.25])) for fn in (min,max)],
             terminal_potential=float(m['A']),terminal_momentum=float(m['P']),
             born_kinetic_in_proper_comoving_units=float(m['birthkin']),
             terminal_spectrum=sp,angular_change_of_mean_gamma=abs(sp['mean_gamma']/sp8['mean_gamma']-1),
             angular_change_of_signed_current=abs(sp['mean_signed_velocity']-sp8['mean_signed_velocity']))
    res['trace_columns']='t field_share net_field_loss signed_fraction normalized_total_balance'.split()
    res['snapshots']=[[float(tr['t'][i]),float(frac[i]),float(1-tr['field'][i]/norm),float(tr['current'][i]),float(tr['net_residual'][i]/norm)] for i in np.linspace(0,len(tr)-1,41,dtype=int)]
    return res


def calculate(source):
    bg=background(source)
    with tempfile.TemporaryDirectory() as directory:
        temp=Path(directory)
        exe=temp/'kinetic'
        subprocess.run(['g++','-std=c++17','-O3',str(HERE/'vlasov.cpp'),'-o',str(exe)],check=True)
        setup=[('zero',512,8,0.),('coarse',256,32,bg['Q']),('medium',512,32,bg['Q']),
               ('fine',1024,32,bg['Q']),('finer_birth',2048,16,bg['Q']),('finer_time',1024,64,bg['Q']),('reference',4096,32,bg['Q']),
               ('more_births',8192,16,bg['Q']),('smaller_step',4096,64,bg['Q'])]
        runs={name:run(temp,exe,bg,bins,substeps,Q) for name,bins,substeps,Q in setup}
        runs['linear_birth_history']=run(temp,exe,background(source,True),2048,16,bg['Q'])
    ts=np.linspace(bg['start'],bg['end'],32769)
    eta=float(np.trapezoid(1/bg['sample'](ts)[:,1],ts))
    A0,P0=bg['A'],bg['P']
    exactA=A0*math.cos(2*eta)+P0*math.sin(2*eta)/2
    exactP=P0*math.cos(2*eta)-2*A0*math.sin(2*eta)
    control=max(abs(runs['zero']['terminal_potential']-exactA),abs(runs['zero']['terminal_momentum']-exactP))
    if control>5e-8 or abs(runs['more_births']['turn_net_field_loss_fraction']-runs['reference']['turn_net_field_loss_fraction'])>.002:
        raise ValueError('Free-wave or refined kinetic comparison failed.')
    for r in runs.values():
        if r['terminal_spectrum']['current_agreement']>1e-6 or r['max_normalized_total_balance']>.001:
            raise ValueError('Current angular integral or energy budget failed.')
    return dict(background={k:v for k,v in bg.items() if k!='sample'},zero_charge_analytic_field_error=control,runs=runs)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=HERE.parent/'coupled/DECAY.json')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(args.source),indent=2,allow_nan=False)+'\n')
