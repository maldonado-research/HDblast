"""Positive squared-amplitude replay and phase bounds for squeezed parent populations."""
import argparse
import importlib.util
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import simpson


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("moving_mode_transfer", HERE/"transfer.py")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
PHASE = -.5
BRANCHES = ("quiet","plus","minus","far")
OUTPUTS = ("vacuum","squeezed","incoherent_excess","excess","cross_real","cross_imag","profile_squeezed")


def measure(g, eta, q, w, weights, af):
    # Pair total spatial momentum +q. The two terms multiply a_q and a†_-q.
    # This uses the phase convention chi(initial)>0 in each incoming mode.
    phase = np.exp(1j*np.outer(w,eta))
    fplus = (eta[1]-eta[0])*(phase@g)
    fminus = (eta[1]-eta[0])*(phase.conj()@g)
    vacuum = abs(fminus)**2
    squeezed = abs(np.sqrt(2)*fminus.conj()+np.exp(1j*PHASE)*fplus)**2
    profile_angle = -2.2+3*(q/8)**2
    profile = abs(np.sqrt(2)*fminus.conj()+np.exp(1j*profile_angle)*fplus)**2
    incoherent = abs(fplus)**2+vacuum
    cross = 2*np.sqrt(2)*fplus*fminus
    excess = squeezed-vacuum
    if not np.allclose(excess,incoherent+(np.exp(1j*PHASE)*cross).real,atol=2e-15,rtol=1e-11):
        raise ValueError("Gaussian expectation disagrees with direct amplitude")
    if not np.allclose(profile-vacuum,incoherent+(np.exp(1j*profile_angle)*cross).real,atol=2e-15,rtol=1e-11):
        raise ValueError("Quadratic-phase expectation disagrees with direct amplitude")
    densities = [vacuum,squeezed,incoherent,excess,cross.real,cross.imag,profile]
    phase_space = weights[:,None]*(w[:,None]**2-q**2)/(8*np.pi**2)
    return np.stack([np.vstack([np.sum(v*phase_space,axis=0) for v in densities]),
                     np.vstack([np.sum(v*phase_space*(w/af)[:,None],axis=0) for v in densities])])


def integrated(curves,q,limits,frequency_order):
    eta,s,a,v,dv,_=curves
    lo,hi=limits
    active=(s>=lo)&(s<=hi)
    f=np.where(active,np.sin(np.pi*(s-lo)/(hi-lo))**4,0.)
    every=(len(eta)-1)//2048
    keep=np.flatnonzero(active[::every])
    inds=np.arange(max(0,keep[0]-1),keep[-1]+2)*every
    eta=eta[inds]
    g=f[inds,None]*v[inds]
    values=np.zeros((2,len(OUTPUTS),len(q)))
    for i,k in enumerate(q):
        if k<8:
            w,c=base.nodes(k,8,frequency_order)
            values[:,:,i]=measure(g[:,i:i+1],eta,q[i:i+1],w,c,a[-1])[:,:,0]
    results={}
    for dest,low,high in [('screen',8.,60.),('full',60.,768.)]:
        panels=[base.nodes(x,min(x+16.,high),frequency_order) for x in np.arange(low,high,16.)]
        w,c=(np.concatenate(vs) for vs in zip(*panels))
        values+=measure(g,eta,q,w,c,a[-1])
        avg=3*simpson(values*q**2,x=q,axis=-1)/8**3
        results[dest]={k:dict(zip(OUTPUTS,map(float,avg[j]),strict=True)) for j,k in enumerate(('number','energy'))}
        results[dest]["per_q"] = values.copy()
    return results


def run(repo):
    ref=json.loads(base.pinned(repo/base.REF,base.REF_SHA))['matter']['events'][0]['A_star']
    report={'squeeze_occupation':1,'phase_radians':PHASE,'profile_radians':'-2.2 + 3*(q/8)^2','refinements':{}}
    with zipfile.ZipFile(io.BytesIO(base.pinned(repo/base.WALL,base.WALL_SHA))) as archive:
        for steps,intervals,order in ((8192,32,12),(16384,64,16),(32768,128,16)):
            q=np.linspace(0.,8.,intervals+1)
            arms={}
            for arm in BRANCHES:
                curves=base.modes(base.historical(archive,arm,ref),q,steps)
                arms[arm]=dict(a_cut=float(curves[2][-1]),wronskian=curves[-1],
                               windows={key:integrated(curves,q,support,order) for key,support in base.WINDOWS.items()})
            contrasts={}
            for arm in BRANCHES[1:]:
                contrasts[arm]={}
                for key in base.WINDOWS:
                    contrasts[arm][key]={}
                    for band in ("full","screen"):
                        row={}
                        for j,kind in enumerate(("number","energy")):
                            quiet=arms['quiet']['windows'][key][band][kind]
                            pulse=arms[arm]['windows'][key][band][kind]
                            norm=quiet['incoherent_excess']
                            # Report both extensive comoving count/energy per input parent
                            # and locally measured excess per fixed physical volume.
                            ratio=(arms['quiet']['a_cut']/arms[arm]['a_cut'])**3
                            diff=lambda name: (pulse[name]-quiet[name])/norm
                            c=complex(diff('cross_real'),diff('cross_imag'))
                            di=diff('incoherent_excess')
                            qi=arms['quiet']['windows'][key][band]['per_q'][j]
                            pi=arms[arm]['windows'][key][band]['per_q'][j]
                            qr=qi[OUTPUTS.index('cross_real')]+1j*qi[OUTPUTS.index('cross_imag')]
                            pr=pi[OUTPUTS.index('cross_real')]+1j*pi[OUTPUTS.index('cross_imag')]
                            avg_amplitude=lambda cross: 3*simpson(abs(cross)*q**2,x=q)/(8**3*norm)
                            extent=avg_amplitude(pr-qr)
                            local_extent=avg_amplitude(ratio*pr-qr)
                            absolute=di+diff('vacuum')
                            local_di=(ratio*pulse['incoherent_excess']-norm)/norm
                            local_abs=local_di+(ratio*pulse['vacuum']-quiet['vacuum'])/norm
                            bounds=lambda x,r: [100*(x-r),100*(x+r)]
                            ds=diff('excess')
                            if abs(ds-(di+(np.exp(1j*PHASE)*c).real))>1e-12:raise ValueError('Contrast identity')
                            profile_absolute=diff('profile_squeezed')
                            profile_local=(ratio*pulse['profile_squeezed']-quiet['profile_squeezed'])/norm
                            local_vacuum=(ratio*pulse['vacuum']-quiet['vacuum'])/norm
                            row[kind]=dict(quadratic_phase_excess_percent=100*(profile_absolute-diff('vacuum')),
                                quadratic_phase_absolute_percent=100*profile_absolute,
                                quadratic_phase_local_absolute_percent=100*profile_local,
                                quadratic_phase_local_excess_percent=100*(profile_local-local_vacuum),
                                incoherent_percent=100*di,squeezed_percent=100*ds,
                                interval_percent=[100*(di-abs(c)),100*(di+abs(c))],
                                interval_all_phases_percent=bounds(di,extent),
                                absolute_interval_all_phases_percent=bounds(absolute,extent),
                                local_interval_all_phases_percent=bounds(local_di,local_extent),
                                local_absolute_interval_all_phases_percent=bounds(local_abs,local_extent),
                                change_in_absolute_squeezed_percent=100*diff('squeezed'),
                                local_squeezed_percent=100*(ratio*pulse['excess']-quiet['excess'])/norm,
                                local_incoherent_percent=100*(ratio*pulse['incoherent_excess']-norm)/norm)
                        contrasts[arm][key][band]=row
            for arow in arms.values():
                for window in arow['windows'].values():
                    for spectrum in window.values():
                        spectrum.pop('per_q')
            report['refinements'][str(intervals)]=dict(arms=arms, contrasts_percent=contrasts)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository',type=Path)
    print(json.dumps(run(parser.parse_args().repository),indent=2,allow_nan=False))
