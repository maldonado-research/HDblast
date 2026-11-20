"""Finite switched-window Yukawa response at order g² in an external FLRW history."""
import math

import numpy as np
from scipy.special import expit
from scipy.integrate import simpson


def smooth(x):
    x=np.asarray(x)
    y=np.zeros_like(x)
    y[x>=1]=1.
    use=(x>0)&(x<1)
    y[use]=expit(1/(1-x[use])-1/x[use])
    return y


def response(psi,radius,turn,mu,*,scenario='rolling',ramp=.45,dt=.003,nv=1800,vmax=5.,phase_cap=1.8):
    start=9.
    # Resolve the largest phase advance as the momentum cutoff is increased.
    maxfreq=2*mu*math.hypot(1.,vmax*float(radius(turn)/radius(start)))
    dt=min(dt,phase_cap/maxfreq)
    t=np.linspace(start,turn,1+math.ceil((turn-start)/dt))
    h=t[1]-t[0]
    mid=(t[1:]+t[:-1])/2
    window=smooth((mid-start)/ramp)*smooth((turn-mid)/ramp)
    a=radius(mid)
    ref=float(radius(turn))
    if scenario=='rolling':
        R=.95+.1*smooth((psi(mid)-.77)/.03+.5)
    else:
        R=np.full_like(mid,1.05 if scenario=='static_open' else .95)
        a=np.full_like(mid,ref)
    M=2*mu*R
    # Uniform comoving momentum k = mu * a(turn) * v. The dimensionless
    # momentum cutoff vmax is increased for independent UV-tail checks.
    v=np.linspace(0,vmax,nv+1)
    k=mu*ref*v
    density=np.zeros_like(v)
    for i in range(0,len(k),64):
        ks=k[i:i+64,None]
        p=ks/a[None,:]
        E=np.hypot(mu,p)
        mismatch=2*E-M
        increments=h*mismatch
        phase=np.cumsum(increments,axis=1)-.5*increments
        # Exact integral of the midpoint frequency on each cell.
        oscillation=np.exp(1j*phase)*np.sinc(increments/(2*math.pi))
        pref=window*p/(a**1.5*E*np.sqrt(M))
        amp=h*np.sum(pref*oscillation,axis=1)
        density[i:i+len(amp)]=ks[:,0]**2*np.abs(amp)**2/(2*math.pi**2)
    probability=float(simpson(density,x=k))
    markov=float(h*np.sum(window**2*M*np.maximum(0.,1.-1./R**2)**1.5)/(8*math.pi))
    return dict(coherent_probability_per_g2=probability,local_first_order_per_g2=markov,
                ratio_to_local=probability/markov if markov else None,
                probability_at_g_0p01=probability*1e-4,
                weight_omitted_by_halving_cutoff=1.-float(simpson(density[:nv//2+1],x=k[:nv//2+1]))/probability)


def calculate(psi,a,turn):
    rows={}
    for mu in (8.,32.,128.):
        values={sc:response(psi,a,turn,mu,scenario=sc) for sc in ('rolling','static_closed','static_open')}
        if mu in (8.,128.):
            ref=values['rolling']['coherent_probability_per_g2']
            finer=response(psi,a,turn,mu,dt=.0015,nv=3600,phase_cap=.9)['coherent_probability_per_g2']
            wider=response(psi,a,turn,mu,dt=.003,nv=3600,vmax=10.)['coherent_probability_per_g2']
            values['relative_finer_time_and_momentum']=abs(finer/ref-1)
            values['relative_extended_cutoff']=abs(wider/ref-1)
        rows[str(mu)]=values
    for mass, v in rows.items():
        for protocol in ('rolling','static_open','static_closed'):
            if not 0 < v[protocol]['probability_at_g_0p01'] < .001:
                raise ArithmeticError('Weak yield left leading-order regime.')
        for key in ('relative_finer_time_and_momentum','relative_extended_cutoff'):
            if key in v and v[key]>.001:
                raise ArithmeticError('Finite-window integration did not converge.')
    if abs(rows['128.0']['rolling']['ratio_to_local']-1)>.02:
        raise ArithmeticError('Expected heavy-mass continuum check failed.')
    closed_short=response(psi,a,turn,32.,scenario='static_closed',ramp=.225)
    main_short=response(psi,a,turn,32.,ramp=.225)
    return dict(proper_time_start=9.,proper_time_end=turn,endpoint_window_width=.45,
                coupling_squared_factored_out=True,spin_and_antiparticle_normalization='one Dirac pair per decay; spin sum included',
                cases=rows,short_window_edge_at_mu32=dict(ramp=.225,
                    static_closed=closed_short,rolling=main_short))
