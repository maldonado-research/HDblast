"""Reconstruct a bulk annulus from consistent spherical timelike junction traces."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
KAPPA = 6.
DEPTH = .004


def load_wall():
    spec = importlib.util.spec_from_file_location('curved_wall_trace',ROOT/'loading/loading.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sources(fields, plus, minus):
    r,L,A = fields
    scale = np.exp(2*L)
    # J+ J- = (d_x f)^2 - (d_t f)^2.
    nr = plus[0]*minus[0]/scale
    na = plus[2]*minus[2]/scale
    mag2 = (2*A/r**2)**2
    boxr = 4*r+2*(1-nr)/r-KAPPA*r*mag2/3
    return np.array([scale*boxr,
                     scale*(4-1.5*boxr/r-KAPPA*na/(6*r**2)+KAPPA*mag2/3),
                     -(.5/r)*(plus[0]*minus[2]+minus[0]*plus[2])+4*scale*A/r**2])


def reconstruct(solution, start, end, dx, depth=DEPTH):
    layers = round(depth/dx)
    if abs(layers*dx-depth)>1e-12:
        raise ValueError('Require an integral number of inward null cells')
    # Two guard points at each end for the independent fourth-order constraint.
    intervals = round((end-start)/dx)
    idx = np.arange(-layers-2,intervals+layers+3)
    t = start+dx*idx
    row = solution.trace_many(t)
    a = row['a']
    f = np.array([a,np.zeros_like(a),row['A']])
    fx = np.array([-a*(1+row['e']),-1+2*row['e']+3*row['p'],-a*row['Y']])
    ft = np.array([a*row['H'],np.zeros_like(a),a*row['X']])
    jp, jm = fx+ft, fx-ft
    min_strip_mass = np.full(intervals+1,np.inf)
    min_strip_F = np.full(intervals+1,np.inf)
    def sample_geometry(row,p,m,guard):
        area=row[0,guard]
        norm=p[0,guard]*m[0,guard]*np.exp(-2*row[1,guard])
        np.minimum(min_strip_F,norm,out=min_strip_F)
        np.minimum(min_strip_mass,area**2*(1+area**2-norm),out=min_strip_mass)

    for level in range(layers):
        guard=layers-level+2
        sample_geometry(f,jp,jm,slice(guard,-guard))
        src = sources(f,jp,jm)
        # Integrate J+ and J- from opposite feet along exact null lines.
        # Iterate trapezoidal nonlinear right-hand sides at their intersection.
        left = jp[:,2:] + .5*dx*src[:,2:]
        right = jm[:,:-2] + .5*dx*src[:,:-2]
        base = f[:,1:-1]+.25*dx*(jp[:,1:-1]+jm[:,1:-1])
        guess = src[:,1:-1]
        for __ in range(3):
            pnew, mnew = left+.5*dx*guess, right+.5*dx*guess
            new = base + .25*dx*(pnew+mnew)
            guess = sources(new,pnew,mnew)
        f, jp, jm = new, pnew, mnew
    r,L,A = f
    eL = np.exp(L)
    valI, valW = jp[2]/(2*eL*r), -jm[2]/(2*eL*r)
    F = jp[0]*jm[0]/eL**2
    mass = r**2*(1+r**2-F)
    fixed = slice(2,-2)
    def deriv(x):
        return (x[:-4]-8*x[1:-3]+8*x[3:-1]-x[4:])/(12*dx)
    sample_geometry(f,jp,jm,fixed)
    sr = sources(f,jp,jm)[0, fixed]
    cplus = sr+2*deriv(jp[0])-2*jp[0,fixed]*jp[1,fixed]+KAPPA*jp[2,fixed]**2/(3*r[fixed])
    cminus = sr-2*deriv(jm[0])-2*jm[0,fixed]*jm[1,fixed]+KAPPA*jm[2,fixed]**2/(3*r[fixed])
    zp=1+abs(sr)+abs(2*deriv(jp[0]))+abs(2*jp[0,fixed]*jp[1,fixed])+KAPPA*jp[2,fixed]**2/(3*r[fixed])
    zm=1+abs(sr)+abs(2*deriv(jm[0]))+abs(2*jm[0,fixed]*jm[1,fixed])+KAPPA*jm[2,fixed]**2/(3*r[fixed])
    # A third check uses the derivative of the reconstructed geometric mass.
    rdot=(jp[0]-jm[0])/(2*eL)
    normal=-(jp[0]+jm[0])/(2*eL)
    pressure=valI**2+valW**2-.5*(2*A/r**2)**2
    flux=valI**2-valW**2
    mdot=-(2*KAPPA/3)*r**3*(rdot*pressure+normal*flux)
    return dict(t=start+dx*np.arange(intervals+1), r=r[fixed], lapse=eL[fixed],
                I=valI[fixed], W=valW[fixed], M=mass[fixed], F=F[fixed],
                Cplus=cplus,Cminus=cminus,Cscaled=np.maximum(abs(cplus)/zp,abs(cminus)/zm),
                Mres=deriv(mass)/eL[fixed]-mdot[fixed], strip_min_M=min_strip_mass,
                strip_min_F=min_strip_F)


def block(sol, start, stop, dx):
    n=round((stop-start)/dx)
    pieces=[]
    for first in range(0,n,2000):
        last=min(n,first+2000)
        series = reconstruct(sol,start+first*dx,start+last*dx,dx)
        if pieces:
            series={k:v[1:] for k,v in series.items()}
        pieces.append(series)
    return {k:np.concatenate([part[k] for part in pieces]) for k in pieces[0]}


def stats(row):
    t,r,el,I,W=(row[k] for k in ['t','r','lapse','I','W'])
    def integrate(values):
        return float(np.sum((values[:-1]+values[1:])*.5*np.diff(t)))
    return dict(t0=float(t[0]),t1=float(t[-1]),minimum_mass=float(min(row['M'])),
                time_minimum_mass=float(t[np.argmin(row['M'])]),
                minimum_F=float(min(row['F'])),minimum_mass_on_strip=float(min(row['strip_min_M'])),
                minimum_F_on_strip=float(min(row['strip_min_F'])),
                peak_absolute_emitted=float(max(abs(I))),time_peak_emitted=float(t[np.argmax(abs(I))]),
                peak_absolute_return=float(max(abs(W))),
                emitted_integrated=integrate(r**3*el*I**2),
                returned_integrated=integrate(r**3*el*W**2),
                max_null_constraint=max(float(max(abs(row[k]))) for k in ('Cplus','Cminus')),
                normalized_constraint=float(max(row['Cscaled'])),
                max_mass_flux_residual=float(max(abs(row['Mres']))))


def calculate():
    sol=load_wall().solution(dose=.123,calibration_dose=.12)
    start,split=-.008,.024
    # Equal spacings inside each piece, overlapping guard curves discarded at joins.
    cutoff=sol.config['cutoff']
    late_end=min(18.,np.floor((sol.end-DEPTH-.001)/.002)*.002)
    levels={}
    previous=None
    for scale in [1,2,4]:
        early=block(sol,start,split,.00004/scale)
        late=block(sol,split,late_end,.0004/scale)
        quiet=(early['t']>cutoff+DEPTH)&(early['t']<=split)
        signal_after=stats({k:v[quiet] for k,v in early.items()})
        post=(early['t']>=0.)
        data=dict(early=stats(early),late=stats(late),early_after_both_wavefronts=signal_after,
                  after_direct_emission=stats({k:v[post] for k,v in early.items()}))
        if previous is not None:
            prev_early,prev_late=previous
            data['shared_grid_max_I_change']={
                label:float(max(abs(now['I'][::2]-old['I'])))
                for label,now,old in [('early',early,prev_early),('late',late,prev_late)]}
        previous=early,late
        if scale==4:
            ii=np.arange(len(early['t']))
            sample_idx=ii[(ii%50==0)|((ii%25==0)&(early['t']>=-.005)&(early['t']<=.01))]
            data['late_decimated']=[{k:float(late[k][i]) for k in ('t','I','W','M')}
                                    for i in range(0,len(late['t']),1000)]
            data['source_samples']=[dict(t=float(early['t'][i]),I=float(early['I'][i]),W=float(early['W'][i]),
                                         mass=float(early['M'][i])) for i in sample_idx]
            for event,items in sol.events.items():
                if not items:
                    continue
                ti=items[-1]['t']
                if not split <= ti <= late['t'][-1]:
                    continue
                ii=int(np.argmin(abs(late['t']-ti)))
                data.setdefault('late_samples',{})[event]={k:float(late[k][ii]) for k in ('t','r','lapse','I','W','M','F')}
        levels[str(scale)]=data
    finer=block(sol,start,split,.000005)
    launch_fine=stats(finer)
    post=finer['t']>=0.
    post_launch=stats({k:v[post] for k,v in finer.items()})
    delayed_fraction=(post_launch['emitted_integrated']+levels['4']['late']['emitted_integrated'])/(
        launch_fine['emitted_integrated']+levels['4']['late']['emitted_integrated'])
    # Before the launch all fields vanish and the enclosed mass equals four.
    before=block(sol,-.032,-.012,.0001)
    vacuum_error=float(max(abs(before['M']-4)))
    if (launch_fine['minimum_mass_on_strip']<=0. or levels['4']['late']['minimum_mass_on_strip']<=0. or
            launch_fine['minimum_F_on_strip']<=0. or levels['4']['late']['minimum_F_on_strip']<=0. or
            vacuum_error>1e-7 or max(abs(before['I']))+max(abs(before['W']))>1e-14 or
            launch_fine['normalized_constraint']>2e-3):
        raise ArithmeticError('Curvature, mass, constraint, or prelaunch vacuum control failed')
    return dict(depth=DEPTH,launch_fine=launch_fine,post_direct_emission=post_launch,
                delayed_emitted_energy_fraction=float(delayed_fraction),
                vacuum_enclosed_mass_error=vacuum_error,
                source_clock='d proper time = exp(L) dt',
                wall_config=sol.config,wall_end=sol.end,levels=levels)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    path=parser.parse_args().output
    path.write_text(json.dumps(calculate(),sort_keys=True,indent=2)+'\n')
