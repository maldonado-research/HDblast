"""Independent minimal, free scalar modes on the three rolling, coupled B1 histories."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicSpline

def trajectories(root, grid, branch):
    stem=root/grid
    cp=np.load(stem/'checkpoint.npz')
    with np.load(stem/branch/'rolling_timeseries.npz') as data:
        c={k:data['rec'][:,i] for i,k in enumerate(data['cols'])}
    hist=cp['hist']; use=hist[:,0]<c['H0tau'][0]-1e-10
    s=np.r_[hist[use,0],c['H0tau']]
    phi=CubicSpline(s,np.r_[hist[use,1],c['phi_b']])
    h=CubicSpline(s,np.r_[hist[use,2],c['H_over_H0']])
    int_h=h.antiderivative()
    quiet=json.loads((stem/'quiet/rolling_summary.json').read_text())
    split=c['H0tau'][0]; offset=c['A_b'][0]-quiet['matter']['events'][0]['A_star']
    lna=lambda t:offset+int_h(t)-int_h(split)
    drift=max(abs((c['A_b']-quiet['matter']['events'][0]['A_star'])-lna(c['H0tau'])))
    return phi,h,lna,float(drift)

def compute(root,grid,branch,start,step,rtol,order):
    phi,h,lna,drift=trajectories(root,grid,branch)
    k=np.arange(0,60+step/2,step)
    n=len(k)
    def omega2(t):
        return k*k*np.exp(-2*lna(t))+10000*(phi(t)-.5)**2-2.25*h(t)**2-1.5*h(t,1)
    def W(t,which):
        # Differentiation only at adiabatic endpoints. Evolution uses full omega^2.
        f=omega2(t); om=np.sqrt(f); eps=.001
        if which==0: return om
        lo=np.sqrt(omega2(t-eps)); hi=np.sqrt(omega2(t+eps))
        d=(hi-lo)/(2*eps); dd=(hi-2*om+lo)/eps**2
        return om-dd/(4*om**2)+3*d*d/(8*om**3)
    def basis(t,which):
        freq=W(t,which)
        slope=(W(t+.001,which)-W(t-.001,which))/.002
        f=1/np.sqrt(2*freq); fp=(-1j*freq-.5*slope/freq)*f
        return f,fp,float(max(abs(slope)/freq**2))
    def rhs(t,u):
        return np.concatenate((u[n:],-omega2(t)*u[:n]))
    f,fp,in_ad=basis(start,order)
    ends=[1.95,2.02,2.10,2.15]
    sol=solve_ivp(rhs,(start,ends[-1]),np.r_[f,fp],method='DOP853',t_eval=ends,rtol=rtol,atol=rtol/50)
    if not sol.success: raise ValueError(sol.message)
    event=json.loads((root/grid/branch/'rolling_summary.json').read_text())['matter']['events'][0]
    a_star=float(np.exp(lna(event['s_star'])))
    chosen=k<=8
    event_n=float(simpson(k[chosen]**2*event['number_factor']*np.exp(-np.pi*(k[chosen]/a_star)**2/event['q'])/(2*np.pi**2),x=k[chosen]))
    event_leading=event_n/event['number_factor']
    # Subdivide at all polynomial critical points of phi and ln(a) plus a fine grid.
    # On each subinterval, max |m| and min a each occur at an endpoint.
    roots=np.r_[phi.derivative().roots(extrapolate=False),h.roots(extrapolate=False)]
    result={}
    for i,t in enumerate(ends):
        X,P=sol.y[:n,i],sol.y[n:,i]
        f,fp,out_ad=basis(t,order)
        occ=np.abs(f*P-fp*X)**2
        wr=max(abs((1j*(np.conj(X)*P-np.conj(P)*X)).real-1))
        a=np.exp(lna(t)); mass=100*abs(phi(t)-.5)
        dens=k*k*occ/(2*np.pi**2)
        num=lambda lim,skip=1: float(simpson(dens[k<=lim][::skip],x=k[k<=lim][::skip]))
        rho=float(simpson(dens*np.sqrt((k/a)**2+mass**2),x=k)/a**3)
        knots=np.unique(np.r_[np.linspace(start,t,2001),roots[(roots>start)&(roots<t)]])
        ak=np.exp(lna(knots)); mk=np.abs(100*(phi(knots)-.5))
        Emax=float(np.max(np.hypot(8/np.minimum(ak[:-1],ak[1:]),np.maximum(mk[:-1],mk[1:]))))
        band_rho=float(simpson(dens[chosen]*np.sqrt((k[chosen]/a)**2+mass**2),x=k[chosen])/a**3)
        result[str(t)]=dict(comoving_N_8=num(8),N8_half_grid=num(8,2),N8_phys_density=num(8)/a**3,
                            N8_rho_quasiparticle=band_rho,band_E_upper_H0=Emax,
                            comoving_N=num(60),N_45=num(45),N_30=num(30),N_half_grid=num(60,2),
                            physical_rho_quasiparticle=rho,physical_n=num(60)/a**3,
                            a=float(a),mass_H0=float(mass),wronskian_max=float(wr),endpoint_adiabaticity=out_ad)
    return dict(grid=grid,branch=branch,start=start,k_step=step,rtol=rtol,basis_order=order,max_lna_drift=drift,
                event_a_star=a_star,event_comoving_N_8=event_n,event_leading_comoving_N_8=event_leading,
                M5_over_H0=6.49e-6**(-1/3),start_adiabaticity=in_ad,ends=result)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--grid',required=True);parser.add_argument('--branch',required=True)
    parser.add_argument('--start',type=float,default=.80);parser.add_argument('--step',type=float,default=.25)
    parser.add_argument('--rtol',type=float,default=2e-9);parser.add_argument('--order',type=int,choices=(0,2),default=2)
    a=parser.parse_args()
    result=compute(a.runs,a.grid,a.branch,a.start,a.step,a.rtol,a.order)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__': main()
