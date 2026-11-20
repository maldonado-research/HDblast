"""Finite-band Yukawa response of closed-wall Dirac modes on a prescribed parent history."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from scipy.special import expit

HERE=Path(__file__).resolve().parent


def smooth(x):
    x=np.asarray(x); z=np.zeros_like(x,dtype=float)
    z[x>=1]=1
    live=(x>0)&(x<1)
    z[live]=expit(1/(1-x[live])-1/x[live])
    return z


class Parent:
    def __init__(self,doc,case):
        r=doc['cases'][case]; t,a,H,psi,pi,chi,chidot=np.asarray(r['rows']).T
        self.radius=CubicHermiteSpline(t,a,a*H)
        self.psi=CubicHermiteSpline(t,psi,pi)
        self.mu=r['mu']; self.opening=r['opening'];self.turn=r['turn']
        def mass(t):
            x=(self.psi(t)-r['gate'])/r['width']+.5
            return 2*self.mu*(r['initial_ratio']+(r['final_ratio']-r['initial_ratio'])*smooth(x))
        self.mass=mass
        def rhs(t,y):
            return (y[1], -3*self.radius(t,1)/self.radius(t)*y[1]-mass(t)**2*y[0])
        # Both moving minima have already returned to zero at the first archived slice.
        ode=solve_ivp(rhs,(t[0],t[-1]),(chi[0],chidot[0]),rtol=3e-10,atol=1e-14,
                      method='DOP853',dense_output=True,max_step=min(.04,.5/self.mu))
        if not ode.success:raise RuntimeError(ode.message)
        self.osc=ode.sol
        rebuilt=self.osc(t)
        scaled_error=np.hypot(rebuilt[0]-chi,(rebuilt[1]-chidot)/mass(t))
        amp=np.hypot(chi,chidot/mass(t))
        self.max_history_error=float(np.max(scaled_error/amp))

    def cells(self,start,end,ramp,dt):
        count=math.ceil((end-start)/dt)
        h=(end-start)/count
        t=start+(np.arange(count)+.5)*h
        window=smooth((t-start)/ramp)*smooth((end-t)/ramp) if ramp else np.ones_like(t)
        a=self.radius(t); M=self.mass(t); chi,vel=self.osc(t)
        rho=(vel**2+(M*chi)**2)/2
        beta2=np.maximum(0.,1-(2*self.mu/M)**2)
        local=math.pi/4*h*np.sum(a**3*rho*beta2**1.5*window**2)
        friction=math.pi/4*h*np.sum(a**3*vel**2*beta2**1.5*window**2)
        return h,a,window*chi,float(local),float(friction)


def rotate(u,v,p,m,dt):
    en=np.hypot(p,m); c=np.cos(dt*en); s=-1j*np.sin(dt*en)/en
    return c*u+s*(m*u+p*v),c*v+s*(p*u-m*v)


def shells(parent,start,end,ramp,nband,dt,lam=1.,nonlinear=False):
    h,radius,drive,local,friction=parent.cells(start,end,ramp,dt)
    n=np.arange(nband); k=n+1.5; deg=2*(n+1)*(n+2)
    initial=.5*np.arctan2(k/parent.radius(start),parent.mu)
    u,v=-np.sin(initial).astype(complex),np.cos(initial).astype(complex)
    up,vp,um,vm=u,v,u,v
    amp=np.zeros(nband,dtype=complex)
    for a,f in zip(radius,drive):
        p=k/a
        en=np.hypot(p,parent.mu)
        co=np.cos(.5*h*en); si=-1j*np.sin(.5*h*en)/en
        mid_u,mid_v=co*u+si*(parent.mu*u+p*v),co*v+si*(p*u-parent.mu*v)
        du=(parent.mu*mid_u+p*mid_v)/en
        dv=(p*mid_u-parent.mu*mid_v)/en
        sinc=np.sinc(en*h/np.pi)
        # Exact first derivative of the same frozen-cell unitary used for finite lambda.
        amp+=-1j*h*f*((1+sinc)*mid_u*mid_v-(1-sinc)*du*dv)
        u,v=co*mid_u+si*en*du,co*mid_v+si*en*dv
        if nonlinear:
            up,vp=rotate(up,vp,p,parent.mu+lam*f,h)
            um,vm=rotate(um,vm,p,parent.mu-lam*f,h)
    p1=np.abs(amp)**2
    count=float(deg@p1); connected=p1
    summary=dict(pair_response_per_lambda2=count, markov_all_momenta_per_lambda2=local,
                 ratio_to_markov=count/local if local else None,
                 local_velocity_damping_pairs_per_lambda2=friction,
                 ratio_to_velocity_damping=count/friction if friction else None,
                 upper_half_band_share=float((deg[nband//2:]@p1[nband//2:])/count),
                 max_mode_response=float(np.max(p1)),largest_weighted_shell=int(np.argmax(p1*deg)))
    if nonlinear:
        connected=(np.abs(v*up-u*vp)**2+np.abs(v*um-u*vm)**2)/(2*lam**2)
        final=.5*np.arctan2(k/parent.radius(end),parent.mu)
        co,si=np.cos(final),np.sin(final)
        zero=np.abs(co*u+si*v)**2
        plus=np.abs(co*up+si*vp)**2;minus=np.abs(co*um+si*vm)**2
        norm=max(float(np.max(np.abs(abs(a)**2+abs(b)**2-1))) for a,b in ((u,v),(up,vp),(um,vm)))
        summary.update(
            lambda_used=lam,exact_echo_pairs_per_lambda2=float(deg@connected),
            max_exact_common_basis_pair_probability=float(max(plus.max(),minus.max(),zero.max())),
            common_basis_pairs_at_plus=float(deg@plus),common_basis_pairs_at_minus=float(deg@minus),
            common_basis_gravity_only_pairs=float(deg@zero),
            common_basis_sign_even_excess_per_lambda2=float(deg@((plus+minus)/2-zero))/lam**2,
            max_unitarity_error=norm)
        if norm>5e-10:raise ValueError('Lost spinor norm.')
    return summary


def static_check():
    # Exact Fourier transform of w=sin²(pi t/T) on [0,T] for a constant-radius test.
    a,mu,mass,T,A=150.,5.,12.,60.,.003
    n=np.arange(2304);p=(n+1.5)/a;E=np.hypot(mu,p)
    def box(w):
        return T*np.exp(.5j*w*T)*np.sinc(w*T/(2*np.pi))
    def hann(w):
        omega=2*np.pi/T
        return .5*box(w)-.25*(box(w+omega)+box(w-omega))
    amp=(A/2)*(p/E)*(hann(2*E-mass)+hann(2*E+mass))
    finite=float((2*(n+1)*(n+2))@np.abs(amp)**2)
    rho=.5*mass**2*A**2
    local=math.pi/4*a**3*rho*(1-(2*mu/mass)**2)**1.5*(3*T/8)
    class Static:
        def __init__(self):
            self.mu=mu
        def radius(self,t):
            return a
        def cells(self,start,end,ramp,dt):
            count=math.ceil(T/dt);h=T/count
            tm=(np.arange(count)+.5)*h
            return h,np.full(count,a),A*np.cos(mass*tm)*np.sin(np.pi*tm/T)**2,local,local
    discrete=shells(Static(),0.,T,0.,len(n),.005)['pair_response_per_lambda2']
    if abs(discrete/finite-1)>.0005:raise ValueError('Static Fourier and spinor solutions disagree.')
    return dict(radius=a,mu=mu,mass=mass,duration=T,comoving_amplitude=A,shell_count=len(n),
                numerical_mode_response=discrete,solver_to_fourier_ratio=discrete/finite,
                exact_switched_pair_response=finite,local_infinite_volume=local,ratio=finite/local)


def calculate(doc):
    reference=static_check()
    if abs(reference['ratio']-1)>.002:raise ValueError('Golden-rule normalization check failed.')
    cases={}
    for name in ('stable','heavy_stable'):
        bg=Parent(doc,name);o=bg.opening;turn=bg.turn
        settings={'short_crossing':(o-.5,o+.5,.125),
                  'early_closed':(o-1.25,o-.25,.125),
                  'through_turnaround':(o-.5,turn,.45),
                  'gentle_edges':(o-.5,turn,.9)}
        probe=np.linspace(o,turn,math.ceil((turn-o)*bg.mu*10))
        values=bg.osc(probe)[0]
        crossings=np.flatnonzero(values[1:]*values[:-1]<0)
        roots=[brentq(lambda t:bg.osc(t)[0],probe[i],probe[i+1]) for i in (crossings[0],crossings[-1])]
        settings['free_zeros_after_opening']=(roots[0],roots[1],0.)
        results={}
        for label,(lo,hi,ramp) in settings.items():
            dense=label in ('short_crossing','through_turnaround','free_zeros_after_opening')
            entry=shells(bg,lo,hi,ramp,512,.002,lam=1.,nonlinear=dense)
            entry.update(start=lo,end=hi,ramp=ramp,shell_count=512)
            if dense:
                fine=shells(bg,lo,hi,ramp,512,.001)
                wide=shells(bg,lo,hi,ramp,1024,.002)
                entry.update(halved_step_relative_change=fine['pair_response_per_lambda2']/entry['pair_response_per_lambda2']-1,
                             doubled_band_relative_change=wide['pair_response_per_lambda2']/entry['pair_response_per_lambda2']-1)
            if dense and (abs(entry['halved_step_relative_change'])>.003 or
                          abs(entry['doubled_band_relative_change'])>.10 or
                          abs(entry['exact_echo_pairs_per_lambda2']/entry['pair_response_per_lambda2']-1)>.02):
                raise ValueError(f'Convergence check failed for {name}/{label}: {entry}')
            results[label]=entry
        cases[name]=dict(mu=bg.mu,opening=o,turn=turn,relative_history_error=bg.max_history_error,
                         windows=results)
    return dict(static_control=reference,cases=cases)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,default=HERE/'HISTORY.json')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    result=calculate(json.loads(args.data.read_text()))
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
