"""Timelike spherical junction with characteristic Maxwell input and a surface kinetic coupling."""
import argparse
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('previous_scalar_wall',ROOT/'prior/radiation_exit.py')
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
PRIMARY_DOSE = .123
CALIBRATION_DOSE = .12
PRIMARY_WIDTH = .001
PRIMARY_CAP = 40.
PRIMARY_G = math.sqrt(8)


def step(x):
    if x<=0: return 0.
    if x>=1: return 1.
    return float(expit(1/(1-x)-1/x))


def material(psi,width,strength,cap):
    # Smooth C3 finite caps in BOTH directions for the extended trajectory.
    signed=strength*psi
    s=abs(signed)
    log_cap=math.log(cap)
    onset=log_cap-.2
    u=(s-onset)/.4
    if u<=0:
        value,gradient=s,1.
    elif u>=1:
        value,gradient=log_cap,0.
    else:
        value=onset+.4*(u-2.5*u**4+3*u**5-u**6)
        gradient=(1-u)**3*(1+3*u+6*u*u)
    k=width*math.exp(math.copysign(value,signed))
    return k,k*strength*gradient


@lru_cache(maxsize=16)
def waveform(cap,return_width=0.):
    stop=3+.25*math.log(cap-1)+.0625
    def shape(t):
        if not 0<t<stop: return 0.
        x=4*(t-3)
        return math.sqrt(expit(x))*expit(-x)*step(t/.35)*step((stop-t)/.125)
    knots=[.35,3.,stop-.125]
    if not return_width:
        norm=math.sqrt(2*quad(lambda t:shape(t)**2,0.,stop,points=knots,
                              epsabs=2e-14,epsrel=4e-13)[0])
        return stop,lambda t:shape(t)/norm
    area=quad(shape,0.,stop,points=knots,epsabs=2e-14,epsrel=4e-13)[0]
    def lobe(v):
        return math.exp(4.-1/(v*(1-v))) if 0<v<1 else 0.
    lobe_area=quad(lobe,0.,1.,epsabs=2e-14)[0]
    def dual(t):
        return shape(t) if t<=stop else -area*lobe((t-stop)/return_width)/(return_width*lobe_area)
    end=stop+return_width
    norm=math.sqrt(2*quad(lambda t:dual(t)**2,0.,end,points=knots+[stop],
                          epsabs=3e-14,epsrel=4e-13)[0])
    return end,lambda t:dual(t)/norm


class WallSolution:
    def __init__(self,parts,describe,events,config):
        self.parts=parts
        self.describe=describe
        self.events=events
        self.config=config
        self.end=float(parts[-1].t[-1])
    def trace(self,t):
        if t<0:
            a=math.sqrt(4-t*t)
            return dict(t=t,a=a,H=-t/a**2,mu=4.,A=0.,X=0.,Y=0.,I=0.,O=0.,e=0.,p=0.,
                        psi=0.,pi=0.,K=self.config['width'],rad=0.,sheet=0.,
                        Bmag=0.,acc_over_a=-4/a**4,Kprime=0.,V=0.,Tnn=0.,odd_squared=0.,flux=0.)
        if t>self.end+1e-11: raise ValueError('Beyond generated solution')
        for part in self.parts:
            if t<=part.t[-1]+1e-14:
                return self.describe(t,part.sol(t))
        return self.describe(t,self.parts[-1].sol(self.end))
    def trace_many(self,times):
        # Dense interpolation is vectorized; nonlinear material laws are pointwise.
        times=np.asarray(times)
        output=[None]*len(times)
        done=np.zeros(len(times),dtype=bool)
        for j in np.flatnonzero(times<0):
            output[j]=self.trace(times[j]);done[j]=True
        for part in self.parts:
            idx=np.flatnonzero((~done)&(times<=part.t[-1]+1e-12))
            if len(idx):
                states=part.sol(times[idx])
                for j,y in zip(idx,states.T):
                    output[j]=self.describe(times[j],y)
                done[idx]=True
        if not np.all(done): raise ValueError('Beyond generated solution')
        return {key:np.array([r[key] for r in output]) for key in output[0]}


def solution(dose=PRIMARY_DOSE,width=PRIMARY_WIDTH,cap=PRIMARY_CAP,
             g=PRIMARY_G,calibration_dose=CALIBRATION_DOSE, potential=True,damping=4.,
             stop_time=24.,tol=2e-10,scalar_jump=None,coherence=1.,late_gauge=0.,return_width=0.):
    if not 0<=coherence<=1 or dose<0 or width<=0:
        raise ValueError("Invalid input dose, coherence fraction or pulse scale")
    # Calibration stays constant during an amplitude scan with a fixed material.
    strength=g/(width*math.sqrt(calibration_dose))
    stop,signal=waveform(cap,return_width)
    cutoff=stop*width
    def variables(t,y):
        a,h,mu,c,momentum,psi,pi,rad=y[:8]
        k,kprime=material(psi,width,strength,cap)
        if late_gauge:
            activation,derivative=prior.smooth_step((psi-prior.EXIT)/prior.WIDTH)
            k+=late_gauge*activation
            kprime+=late_gauge*derivative/prior.WIDTH
        v,dv,gamma=prior.properties(psi,damping) if potential else (0.,0.,0.)
        beam=math.sqrt(dose/width)*signal(t/width)
        i=math.sqrt(coherence)*beam
        odd_squared=(1-coherence)*beam**2
        x=momentum/k
        o=x-i
        normal=o-i
        magnetic=2*c/a
        sheet=.5*k*(x*x+magnetic**2)
        e=.5*pi**2+v+rad+sheet
        p=.5*pi**2-v+(rad+sheet)/3
        pressure=i*i+o*o+2*odd_squared-.5*magnetic**2
        flux=i*i-o*o
        accel=(1+e)*(1-2*e-3*p)-1-mu/a**4-2*pressure
        return k,kprime,v,dv,gamma,i,x,o,normal,magnetic,sheet,e,p,pressure,flux,accel,odd_squared
    def describe(t,y):
        a,h,mu,c,momentum,psi,pi,rad=y[:8]
        k,dk,v,_,_,i,x,o,n,b,sh,e,p,pressure,flux,acc,odd=variables(t,y)
        return dict(t=float(t),a=float(a),H=float(h),mu=float(mu),A=float(a*c),X=float(x),
                    Y=float(n),I=float(i),O=float(o),e=float(e),p=float(p),psi=float(psi),
                    pi=float(pi),K=k,rad=float(rad),sheet=float(sh),Bmag=float(b),acc_over_a=float(acc),Kprime=float(dk),V=float(v),Tnn=float(pressure),odd_squared=float(odd),flux=float(flux))
    def rhs(t,y):
        a,h,mu,c,momentum,psi,pi,rad=y[:8]
        k,dk,v,dv,gamma,i,x,o,n,b,sh,e,p,pressure,flux,acc,odd=variables(t,y)
        return [a*h,acc-h*h,-4*a**4*(pressure*h+flux*(1+e)),x-h*c,
                -2*h*momentum-4*k*c/a**2-2*n,pi,.5*dk*(x*x-b*b)-dv-(3*h+gamma)*pi,
                gamma*pi*pi-4*h*rad,2*(i*i+odd),2*(o*o+odd),2*a**3*flux,3*a**3*h*p,a**4*gamma*pi*pi]
    def event(fun,direction=0,terminal=False):
        fun.direction=direction;fun.terminal=terminal
        return fun
    def rad90(t,y):
        return y[7]-.9*variables(t,y)[11]
    events={
        'bounce':event(lambda t,y:y[1],1),
        'maximum':event(lambda t,y:y[1],-1),
        'plateau':event(lambda t,y:y[5]-.012,1),
        'exit':event(lambda t,y:y[5]-.36,1),
        'down_end':event(lambda t,y:y[5]-.38,1),
        'rad90':event(rad90),
        'inner':event(lambda t,y:y[0]-1.6,-1,True),
        'outer':event(lambda t,y:y[0]-16.,1,True)}
    if late_gauge or damping==0.:
        def sheet90(t,y):
            v=variables(t,y)
            return v[10]-.9*v[11]
        events['sheet90']=event(sheet90)
    ends=sorted(set([0.,min(stop_time,cutoff),min(stop_time,cutoff+5*width),stop_time]))
    parts=[]
    init=np.zeros(13);init[:3]=[2.,0.,4.]
    if scalar_jump is not None:
        init[1]=scalar_jump['H']
        init[2]=scalar_jump['mu']
        init[6]=math.sqrt(2*scalar_jump['e'])
        init[10]=8*scalar_jump['e']
    recorded={key:[] for key in events}
    for lo,hi in zip(ends,ends[1:]):
        max_step=min(.025, width/12) if lo<cutoff+1e-14 else .025
        if g==0:
            max_step=min(max_step,width)
        if dose>0 and lo>=cutoff and init[6]<0 and not potential:
            max_step=min(max_step,width/cap)
        options=dict(rtol=tol,atol=tol*.001,dense_output=True,max_step=max_step,
                     first_step=min(width/25,max_step,hi-lo),events=list(events.values()))
        # A pulse too weak to climb the potential can roll back across the short
        # kinetic transition. Switch to an implicit method if an explicit trial
        # stage overflows on this stiff branch. The failed stage is discarded.
        try:
            with np.errstate(over='raise',invalid='raise',divide='raise'):
                integrate=solve_ivp(rhs,(lo,hi),init,method='DOP853',**options)
        except FloatingPointError:
            integrate=solve_ivp(rhs,(lo,hi),init,method='Radau',**options)
        if not integrate.success: raise ArithmeticError(integrate.message)
        parts.append(integrate);init=integrate.y[:,-1]
        for name,ts,ys in zip(events,integrate.t_events,integrate.y_events):
            for t,y in zip(ts,ys):
                if t>1e-9 and (name not in ('sheet90','rad90') or variables(t,y)[11]>1e-10):
                    recorded[name].append(describe(t,y))
        if integrate.status==1: break
    return WallSolution(parts,describe,recorded,dict(dose=dose,width=width,g=g,cap=cap,
                         cutoff=cutoff,calibration_dose=calibration_dose,potential=potential,damping=damping,coherence=coherence,late_gauge=late_gauge,return_width=return_width))


def summarize(sol,*,save_trace=False):
    config=sol.config
    all_t=np.concatenate([part.t for part in sol.parts])
    all_y=np.concatenate([part.y for part in sol.parts],axis=1)
    trace=sol.trace_many(all_t)
    a,h,m,e,p=map(trace.get,['a','H','mu','e','p'])
    residual=h*h-e*(2+e)+1/a**2-m/a**4
    angular=np.sqrt(1+a*a-m/a**2+(a*h)**2)/a-(1+e)
    law=a**3*e+all_y[11]-all_y[10]
    e_rad=a**4*trace['rad']-all_y[12]
    inner=1+a*a-m/a**2
    minima=dict(mass=float(min(m)),wall_f=float(min(inner)),surface_energy=float(min(e)),
                surface_null_enthalpy=float(min(e+p)),min_gauge_weight=float(min(trace['K'])))
    errors=dict(friedmann=float(max(abs(residual))),angular=float(max(abs(angular))),
                surface_first_law=float(max(abs(law))),radiation_work=float(max(abs(e_rad))))
    if max(errors.values())>2e-6 or min(m)<-1e-8 or min(inner)<=0 or min(e+p)<-2e-10:
        raise ArithmeticError((minima, errors))
    final=sol.trace(sol.end)
    selected_events={key:items[:1] for key,items in sol.events.items() if items}
    epochs=[]
    boundaries={0.,sol.end}
    for key in ('rad90','bounce','maximum'):
        boundaries.update(row['t'] for row in sol.events[key])
    cuts=sorted(boundaries)
    for lo,hi in zip(cuts,cuts[1:]):
        mid=sol.trace((lo+hi)/2)
        if mid['H']>0 and mid['rad']>.9*mid['e']:
            start,stop=sol.trace(lo),sol.trace(hi)
            epochs.append(dict(start=lo,end=hi,start_radius=start['a'],
                               end_radius=stop['a'],efolds=math.log(stop['a']/start['a'])))
    row=dict(config=config,events=selected_events,final=final,epochs=epochs,minima=minima,errors=errors,
             incident_dose=float(all_y[8,-1]),returned_dose=float(all_y[9,-1]))
    if save_trace:
        row['trajectory']=[sol.trace(t) for t in np.linspace(0.,sol.end,161)]
    return row


def limiting_jump(dose,g,calibration,cap):
    # The chart is a local orthonormal frame with time scaled by K0.
    stop,shape=waveform(cap)
    scaled_g=g*math.sqrt(dose/calibration)
    def drive(t,y):
        p,psi,vel,_,_=y
        k,dk=material(psi,1.,scaled_g,cap)
        x=p/k
        i=shape(t)
        return [4*i-2*x,vel,.5*dk*x*x,2*i*i,2*(x-i)**2]
    inside=solve_ivp(drive,(0.,stop),np.zeros(5),method='DOP853',
                     rtol=2e-12,atol=1e-14,max_step=.025)
    tail=solve_ivp(drive,(stop,stop+30*cap),inside.y[:,-1],method='DOP853',
                   rtol=2e-12,atol=1e-14,max_step=cap/5)
    if not inside.success or not tail.success: raise ArithmeticError('Sheet limit failed')
    p,psi,vel,got,returned=tail.y[:,-1]
    if max(abs(got-1),abs(returned+.5*vel*vel-1))>2e-8:raise ArithmeticError('Sheet energy budget')
    stored=dose*.5*vel**2
    back=dose*returned
    h=-(dose+back)
    mu=4+16*(h*h-stored*(2+stored))
    return dict(incident=dose,returned=back,e=stored,H=h,mu=mu,efficiency=stored/dose)


def calculate_boundary():
    controls={
        'primary':{},'smaller_dose':dict(dose=.12),'larger_dose':dict(dose=.124),
        'weak_pulse':dict(dose=.10), 'no_pulse':dict(dose=0.),
        'zero_coupling':dict(g=0.), 'no_potential':dict(potential=False),
        'no_damping':dict(damping=0.), 'even_99_percent':dict(coherence=.99),
        'bipolar_short':dict(return_width=40.,damping=0.,late_gauge=2000.),
        'bipolar_long':dict(return_width=500.,damping=0.,late_gauge=2000.),
        'even_98_percent':dict(coherence=.98),'even_90_percent':dict(coherence=.90),
        'uncorrelated':dict(coherence=.50),'odd_only':dict(coherence=0.),
        'refined':dict(tol=2e-11)}
    outputs={name:summarize(solution(**cfg),save_trace=name=='primary') for name,cfg in controls.items()}
    jump=limiting_jump(PRIMARY_DOSE,PRIMARY_G,CALIBRATION_DOSE,PRIMARY_CAP)
    thin=summarize(solution(dose=0.,scalar_jump=jump))
    grid={}
    for width in (.002,.001,.0005,.00025,.000125):
        result=outputs['primary'] if width==PRIMARY_WIDTH else summarize(solution(width=width))
        grid[str(width)]={key:dict(t=ev[0]['t'],a=ev[0]['a'],mu=ev[0]['mu']) for key,ev in result['events'].items()}
    if (not outputs['primary']['epochs'] or not outputs['refined']['epochs'] or
            any(outputs[name]['epochs'] for name in ('no_pulse','zero_coupling','no_potential','no_damping'))):
        raise ArithmeticError('Radiation or control classification changed')
    signed_areas={}
    for name in ('bipolar_short','bipolar_long'):
        w=controls[name]['return_width']
        t0,_=waveform(PRIMARY_CAP)
        t1,dual=waveform(PRIMARY_CAP,w)
        signed_areas[name]=quad(dual,0.,t1,points=[.35,3.,t0-.125,t0],
                                epsabs=2e-13,epsrel=1e-12)[0]
        if abs(signed_areas[name])>3e-12:
            raise ArithmeticError('Incident bipolar area failed to cancel')
    return dict(runs=outputs,signed_template_areas=signed_areas,
                thin_limit=dict(jump=jump,instantaneous=thin,finite_width=grid))


def companion(name):
    spec = importlib.util.spec_from_file_location('curved_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def calculate():
    return dict(boundary=calculate_boundary(),return_field=companion('return_field').calculate(),
                forward=companion('forward').calculate())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(),indent=2,allow_nan=False)+'\n')
