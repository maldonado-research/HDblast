"""Gaussian oscillator evolution and a switched detector readout on authenticated wall data."""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np

ZIP_HASH = "9a7069e44c379e68555511e82fe841e0f01a72b23e77fa5ec69caa1373beaf70"
SUMMARY_HASH = "bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f"
CUT = 0.15
WINDOWS = ((0.04, 0.09), (0.10, 0.15))
GAPS = (20., 50., 100., -50.)
WIDTHS = (0., .025)
POINTS = (0., 10., 25., 50., 100.)
INTERVALS = ((0., 5.), (5., 18.), (18., 80.), (80., 320.))
G = 100.


def check(ok, msg):
    if not ok:
        raise ValueError(msg)


def get_series(zf, name, rho):
    data = zf.read(name)
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    def col(key):
        return np.array([float(row[key]) for row in rows])
    t = col("T")
    lapse = np.exp(col("B_wall")) / rho
    clock = np.r_[0., np.cumsum(np.diff(t) * (lapse[1:] + lapse[:-1]) / 2)]
    check(np.all(np.diff(clock)>0) and clock[-1] >= CUT, "Bad wall clock")
    loga = col("A_wall") - float(rows[0]["A_wall"])
    mass = G*(.5-col("phi_wall"))
    check(np.min(mass)>0, "Unexpected crossing")
    # d/ds at the initial slice; at a(0)=1 it also equals d/d dimensionless conformal time.
    dV0 = 2*(mass[0]**2 * float(rows[0]["A_T_wall"])
             - mass[0]*G*float(rows[0]["phi_T_wall"]))/lapse[0]
    return {"s":clock,"loga":loga,"mass":mass,"dV0":dV0,"sha":hashlib.sha256(data).hexdigest()}


def grid(n):
    x,w = np.polynomial.legendre.leggauss(n)
    k = list(POINTS)
    wk = [0.]*len(POINTS)
    for lo,hi in INTERVALS:
        k.extend(((hi+lo)+(hi-lo)*x)/2)
        wk.extend(w*(hi-lo)/2)
    return np.array(k), np.array(wk)


def gaussian_tail(series,K=320.,ell=.025):
    # Enclose total variation of V=a^2 m^2 on segmentwise-linear (log a,m).
    # Includes all saved nodes; the bound through the measurement cut is tighter.
    a = np.exp(series["loga"])
    m = series["mass"]
    tv = 2*max(a)**2*(max(m)**2*sum(abs(np.diff(series["loga"])))+
                       max(m)*sum(abs(np.diff(m))))
    dmax = abs(series["dV0"])/(4*K*K)
    # Frequency is always >= k; allow the larger initial WKB energy as well.
    norm = (math.sqrt(1+(m[0]/K)**2)*(1+dmax*dmax/(2*K*K))
            *math.exp(tv/K**2)/min(a)**2)
    l1_switch = 3*(WINDOWS[1][1]-WINDOWS[1][0])/8
    tail = l1_switch**2 * norm * max(a)**2/(4*math.pi**2*ell**2)*math.exp(-(ell*K/max(a))**2)
    return tail


def propagate(series, k, steps):
    h = CUT/steps
    s = (np.arange(steps)+.5)*h
    a = np.exp(np.interp(s,series['s'],series['loga']))
    m = np.interp(s,series['s'],series['mass'])
    w0 = np.hypot(k,series['mass'][0])
    v0 = (2*w0)**-.5
    v = np.stack([v0,v0]).astype(complex)
    p = -1j*w0*v
    p[1] -= series['dV0']/(4*w0*w0)*v0
    integral = np.zeros((2,len(WINDOWS),len(GAPS),len(WIDTHS),len(k)),complex)
    checkpoints = {}
    # Accumulate with one-step compensation to retain tiny pulse-minus-quiet differences.
    correction = np.zeros_like(integral)
    for j,tm in enumerate(s):
        omega = np.hypot(k,a[j]*m[j])
        angle = omega*h/(2*a[j])
        co,si = np.cos(angle),np.sin(angle)
        mid = co*v+(si/omega)*p
        pmid = co*p-(si*omega)*v
        for wi,(start,end) in enumerate(WINDOWS):
            if not start<tm<end:
                continue
            switch = math.sin(math.pi*(tm-start)/(end-start))**4
            phases = (h*switch/a[j])*np.exp(-1j*np.array(GAPS)*tm)
            smearing = np.exp(-.5*(np.array(WIDTHS)[:,None]*k/a[j])**2)
            entry = (mid[:,None,None,:]*phases[None,:,None,None]*smearing[None,None,:,:]
                     -correction[:,wi])
            tmp = integral[:,wi]+entry
            correction[:,wi] = (tmp-integral[:,wi])-entry
            integral[:,wi] = tmp
        v = co*mid+(si/omega)*pmid
        p = co*pmid-(si*omega)*mid
        if j+1 in (round(.09/h),round(.12/h),steps):
            checkpoints[format((j+1)*h,'.2f')] = (v[:,:len(POINTS)].copy(),p[:,:len(POINTS)].copy())
    return integral,checkpoints,float(np.max(abs(v*p.conjugate()-p*v.conjugate()-1j)))


def responses(raw,k,w):
    # For the unfiltered pointlike operator we retain only explicitly bounded band contributions.
    masks = ((k>5)&(k<18),(k>5)&(k<80),k>0)
    labels = ('band_5_18','band_5_80','gaussian_width_0.025')
    results={}
    for label,mask,wi in zip(labels,masks,(0,0,1)):
        weights = (w*k*k/(2*np.pi**2))*mask
        prob = np.sum(abs(raw[:,:,:,wi,:])**2*weights,axis=-1)
        results[label]=prob
    return results


def difference(p,q):
    # In the shared canonical conformal phase plane |v_q p_p-p_q v_p|^2.
    return {time: (abs(q[time][0]*pp-q[time][1]*vp)**2).tolist()
            for time,(vp,pp) in p.items()}


def run_batch(zf,rho,grid_name,steps,order):
    k,w=grid(order)
    quiet=get_series(zf,'runs/'+grid_name+'/quiet.csv',rho)
    names=[('plus','runs/'+grid_name+'/delayed.csv'),
           ('minus','runs/'+grid_name+'/delayed_negative.csv')]
    if grid_name in ('g4','g8'):
        names += [('half_plus','half_runs/'+grid_name+'/delayed_positive.csv'),
                  ('half_minus','half_runs/'+grid_name+'/delayed_negative.csv')]
    raw,ref,wr=propagate(quiet,k,steps)
    baseline=responses(raw,k,w)
    result={'steps':steps,'quadrature_per_panel':order,'quiet_probability':{n:a.tolist() for n,a in baseline.items()},
            'quiet_wronskian_error':wr,'quiet_gaussian_tail_bound':gaussian_tail(quiet),'runs':{}}
    for name,file in names:
        path=get_series(zf,file,rho)
        vals,pts,wr=propagate(path,k,steps)
        # Evaluate differences before integrating over k; this is resilient to a near-quiet signal.
        delta={}
        for label,mask,idx in (('band_5_18',(k>5)&(k<18),0),('band_5_80',(k>5)&(k<80),0),
                               ('gaussian_width_0.025',k>0,1)):
            term=np.real((vals[:,:,:,idx]-raw[:,:,:,idx])*(vals[:,:,:,idx]+raw[:,:,:,idx]).conjugate())
            delta[label]=np.sum(term*(w*k*k/(2*np.pi**2)*mask),axis=-1).tolist()
        result['runs'][name]={'probability_difference':delta,'n_relative':difference(pts,ref),
                             'wronskian_error':wr,'gaussian_tail_bound':gaussian_tail(path),
                             'csv_sha256':path['sha']}
    return result


def switch_transform(nu,l,r):
    # Exact Fourier transform of the compact sin^4 window including its proper-time offset.
    width=r-l
    def basis(freq):
        return width*np.exp(-1j*freq*(l+width/2))*np.sinc(freq*width/(2*np.pi))
    f=2*np.pi/width
    # Expand cos(j*f*(s-l)) with the appropriate phase at the lower endpoint.
    cos1=(np.exp(-1j*f*l)*basis(nu-f)+np.exp(1j*f*l)*basis(nu+f))/2
    cos2=(np.exp(-2j*f*l)*basis(nu-2*f)+np.exp(2j*f*l)*basis(nu+2*f))/2
    return (3*basis(nu)-4*cos1+cos2)/8


def adiabatic_amplitude(series,k,steps):
    """WKB transport with initial coefficients fixed to the actual common instantaneous state."""
    h=CUT/steps
    s=(np.arange(steps)+.5)*h
    a=np.exp(np.interp(s,series['s'],series['loga']))
    m=np.interp(s,series['s'],series['mass'])
    w0=np.hypot(k,series['mass'][0])
    d0=series['dV0']/(4*w0*w0)
    coeff=1j*d0/(2*w0)
    angle=np.zeros_like(k)
    result=np.zeros_like(k,dtype=complex)
    start,end=WINDOWS[1]
    for j,tm in enumerate(s):
        omega=np.hypot(k,a[j]*m[j])
        inc=omega*h/a[j]
        if start<tm<end:
            f=np.exp(-1j*(angle+inc/2))/np.sqrt(2*omega)
            v=(1+coeff)*f-coeff*f.conjugate()
            switch=math.sin(math.pi*(tm-start)/(end-start))**4
            result+=h*switch*np.exp(-1j*GAPS[0]*tm)*v/a[j]
        angle+=inc
    return result


def source_controls(zf,rho,steps=2400,order=32):
    """Hold either the mass history or the geometry fixed at the same wall proper times."""
    q=get_series(zf,'runs/g8/quiet.csv',rho)
    p=get_series(zf,'runs/g8/delayed.csv',rho)
    k0,w0=np.polynomial.legendre.leggauss(order)
    k=(23+13*k0)/2
    w=w0*13/2*k*k/(2*np.pi**2)
    clocks=np.unique(np.r_[q['s'],p['s']])
    clocks=np.r_[clocks[clocks<CUT],CUT]
    def hybrid(geometry,mass):
        return {'s':clocks,'loga':np.interp(clocks,geometry['s'],geometry['loga']),
                'mass':np.interp(clocks,mass['s'],mass['mass']),'dV0':q['dV0']}
    prob={}
    for name,path in [('quiet',q),('pulse',p),('mass_only',hybrid(q,p)),('geometry_only',hybrid(p,q))]:
        raw,_,_=propagate(path,k,steps)
        z=raw[0,1,0,0]
        wkb=adiabatic_amplitude(path,k,steps)
        prob[name]={'exact':float(np.dot(w,abs(z)**2)), 'WKB_history':float(np.dot(w,abs(wkb)**2))}
    delta={kind: {name:p[kind]-prob['quiet'][kind] for name,p in prob.items()}
           for kind in ('exact','WKB_history')}
    for entry in delta.values():
        entry['interaction']=entry['pulse']-entry['mass_only']-entry['geometry_only']
    maximum_screen_frequency={}
    for nm,path in [('quiet',q),('pulse',p)]:
        # Use a segmentwise upper bound, allowing minima of a and maxima of m at different times.
        a=np.exp(np.interp(clocks,path['s'],path['loga']))
        m=np.interp(clocks,path['s'],path['mass'])
        bound=np.hypot(18/np.minimum(a[:-1],a[1:]),np.maximum(m[:-1],m[1:]))
        maximum_screen_frequency[nm]=float(bound.max())
    return {'probabilities':prob, 'signed_differences':delta,
            'max_physical_frequency_upper_bound_k18':maximum_screen_frequency}


def flat_control():
    k=np.array(POINTS)
    hist={'s':np.array([0.,CUT]),'loga':np.zeros(2),'mass':np.array([50.,50.]),'dV0':0.}
    raw,_,wr=propagate(hist,k,1200)
    omega=np.hypot(k,50.)
    gap=np.array(GAPS)[:,None]
    errs=[]
    for j,(start,end) in enumerate(WINDOWS):
        exact=switch_transform(omega+gap,start,end)/np.sqrt(2*omega)
        errs.append(float(np.max(abs(raw[0,j,:,0]-exact))))
    check(max(errs)<1e-13,'Minkowski switching test failed')
    return {'max_abs_mode_amplitude_error':max(errs),'wronskian_error':wr}


def expanding_control():
    # Exact de Sitter, k=0, a=exp(s) and constant mass yields v_ss+v_s+m^2 v=0.
    m=50.
    k=np.zeros(1)
    lam=math.sqrt(m*m-.25)
    v0=1/math.sqrt(2*m)
    p0=-1j*m*v0
    sine=(p0+v0/2)/lam
    cminus=(v0+1j*sine)/2
    cplus=(v0-1j*sine)/2
    exact=np.array([(cminus*switch_transform(np.array(GAPS)+lam-1.5j,start,end)
                     +cplus*switch_transform(np.array(GAPS)-lam-1.5j,start,end))
                    for start,end in WINDOWS])
    path={'s':np.array([0.,CUT]),'loga':np.array([0.,CUT]),'mass':np.array([m,m]),'dV0':2*m*m}
    errors=[]
    for steps in (1200,2400):
        raw,_,_=propagate(path,k,steps)
        errors.append(float(np.max(abs(raw[0,:,:,0,0]-exact))))
    check(errors[1] < 5e-10 and 3.9<errors[0]/errors[1]<4.1,'Expanding-background control failed')
    return {'time_steps':[1200,2400],'maximum_complex_amplitude_errors':errors}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',type=Path)
    p.add_argument('summary',type=Path)
    args=p.parse_args()
    check(hashlib.sha256(args.archive.read_bytes()).hexdigest()==ZIP_HASH,'Unrecognized pulse histories')
    data=args.summary.read_bytes()
    check(hashlib.sha256(data).hexdigest()==SUMMARY_HASH,'Unrecognized wall summary')
    rho=json.loads(data)['rho_b']
    out={'windows':WINDOWS,'gaps':GAPS,'representative_k':POINTS,'preparations':['instantaneous','WKB_slope'],
         'archive_sha256':ZIP_HASH,'summary_sha256':SUMMARY_HASH,'flat_control':flat_control(),'expanding_control':expanding_control(),'grids':{}}
    with zipfile.ZipFile(args.archive) as zf:
        out['source_controls']=source_controls(zf,rho)
        for name in ('g1','g2','g4','g8'):
            out['grids'][name]=run_batch(zf,rho,name,2400,32)
        for title,steps,order in [('time_refined',4800,32),('momentum_refined',2400,64)]:
            out[title]=run_batch(zf,rho,'g8',steps,order)
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
