"""First-scatter survival of ballistic birth cohorts on fixed fluid-generated fields."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import qmc
from scipy.special import kve

HERE=Path(__file__).resolve().parent
MARKERS=('half_work','nine_tenths_work','first_hot')


def inverse_angle(uniform,c,d,target):
    """Quantile for p(xi) proportional to min(1,abs(c-d*xi)/target)."""
    xi=2*uniform-1
    # A constant weight affects total accepted rest, not the conditional direction.
    active=(target>0)&(d!=0)&(np.abs(c)<np.abs(d)+target)
    if not np.any(active):
        return xi
    f=target[active]; cc=c[active]; dd=d[active]; h=uniform[active]
    cuts=np.column_stack((-np.ones_like(f),np.clip((cc-f)/dd,-1,1),
                          np.clip(cc/dd,-1,1),np.clip((cc+f)/dd,-1,1),np.ones_like(f)))
    cuts.sort(axis=1)
    lo=cuts[:,:-1]; hi=cuts[:,1:]; delta=hi-lo
    g0=np.minimum(1.,np.abs(cc[:,None]-dd[:,None]*lo)/f[:,None])
    g1=np.minimum(1.,np.abs(cc[:,None]-dd[:,None]*hi)/f[:,None])
    sat=np.abs(cc[:,None]-dd[:,None]*(lo+hi)/2)>=f[:,None]
    g0[sat]=1.; g1[sat]=1.
    areas=.5*(g0+g1)*delta
    total=np.cumsum(areas,axis=1)
    x=h*total[:,-1]
    k=np.minimum(3,(total<=x[:,None]).sum(axis=1))
    idx=np.arange(len(f)); scale=delta[idx,k]
    past=np.where(k>0,total[idx,np.maximum(k-1,0)],0.)
    area=np.divide(x-past,scale,out=np.zeros_like(past),where=scale>0)
    start=g0[idx,k]; slope=g1[idx,k]-start
    end=np.sqrt(np.maximum(0.,start**2+2*slope*area))
    u=np.divide(2*area,start+end,out=np.zeros_like(area),where=start+end>0)
    xi[active]=lo[idx,k]+scale*np.clip(u,0.,1.)
    return xi


class History:
    def __init__(self, fname, meta, run):
        self.mu=meta['daughter_mass_L']; self.action=meta['action']
        self.branch=run['channel_fraction']; self.hot=run['first_hot_time']
        raw=np.memmap(fname,mode='r',dtype=np.float64).reshape(-1,len(meta['columns']))
        ci={key:i for i,key in enumerate(meta['columns'])}
        times=raw[:,ci['t']]; rest=raw[:,ci['D']]
        lo=max(0,np.searchsorted(rest,0.,side='right')-1)
        hi=np.searchsorted(times,self.hot,side='right')+1
        table=raw[lo:hi]
        self.data={k:np.array(table[:,i]) for k,i in ci.items()}
        data=self.data
        data['I']=cumulative_trapezoid(data['a']*meta['charge_to_mass']*data['X'],data['t'],initial=0.)
        v=data['u']/np.hypot(1.,data['u'])
        work=cumulative_trapezoid(np.abs(meta['charge_to_mass']*data['D']*v*data['X']),data['t'],initial=0.)
        work_hot=np.interp(self.hot,data['t'],work)
        self.events={'half_work':np.interp(.5*work_hot,work,data['t']),
                     'nine_tenths_work':np.interp(.9*work_hot,work,data['t']),'first_hot':self.hot}
        self.native_work=work_hot
        self.birth_D,idx=np.unique(data['D'],return_index=True)
        self.birth={k:data[k][idx] for k in ('t','a','I','R','c','d','target')}

    def births(self,rest,uniform):
        # Invert the cumulative *accepted* charge rest. A stores both signs/axes.
        shape=np.shape(rest)
        b={k:np.interp(rest.ravel(),self.birth_D,col) for k,col in self.birth.items()}
        xi=inverse_angle(uniform.ravel(),b['c'],b['d'],b['target'])
        ap=b['a']*np.sqrt(np.maximum(0.,(b['R']-1.)*(b['R']+1.)))
        return np.stack((b['t'],ap*xi-b['I'],ap*np.sqrt(np.maximum(0.,1-xi**2))),axis=-1).reshape(*shape,3)

    def value(self,t,col):
        return np.interp(t,self.data['t'],self.data[col])


def max_pair_rate_bound():
    """Global bound on m² sigma v_M for the stated pole, including threshold."""
    pole2=8300.**2
    gamma=.30
    cut=pole2/2
    below=24*math.pi*gamma**2*(cut+2)**2/(cut*(pole2-cut)**2)
    # Ratio maximize x/[(x-1)^2+d*x^2] = 1/[2(sqrt(1+d)-1)].
    d=gamma**2*(1-8/pole2)
    above=(24*math.pi/pole2)*gamma**2*(1+4/pole2)**2*(math.sqrt(1+d)+1)/(2*d)
    return max(below,above)


def bath_independent_floor(history,end):
    d=history.data
    before=d['t']<end
    t=np.append(d['t'][before],end)
    a=np.append(d['a'][before],history.value(end,'a'))
    rest=np.append(d['D'][before],history.value(end,'D'))
    density=history.branch*history.action*rest/(6*history.mu**3*a**3)
    depth=cumulative_trapezoid(density*max_pair_rate_bound(),t,initial=0.)
    # Trapezoidal quadrature over cumulative births, normalized at the endpoint.
    survival=np.exp(depth-depth[-1])
    integral=np.sum(np.diff(rest)*.5*(survival[1:]+survival[:-1]))
    # A monotone cell enclosure avoids treating an endpoint trapezoid as a certified floor.
    # On the expanding interval N is nondecreasing and a is nondecreasing. Use the
    # largest number and smallest radius in each cell, and place its births at its left end.
    if np.any(np.diff(a)<0) or np.any(np.diff(rest)<0):
        raise ValueError('Monotone birth and expansion interval required')
    cell_upper=history.branch*history.action*rest[1:]*max_pair_rate_bound()*np.diff(t)/(6*history.mu**3*a[:-1]**3)
    depth_upper=np.cumsum(cell_upper[::-1])[::-1]
    floor=float(np.diff(rest)@np.exp(-depth_upper)/rest[-1])
    # All unscattered particles still lie in a tube of bounded physical p_perp/m.
    R=max(1.,np.max(d['R'][before]),history.value(end,'R'))
    tube_radius=math.sqrt((R-1.)*(R+1.))
    theta=history.value(end,'theta')
    # Stable evaluation avoids subtracting nearly equal K2 terms when theta is large.
    z,w=np.polynomial.legendre.leggauss(12)
    b=1.+.5*(R-1.)*(z+1.)
    thermal=.5*(R-1.)*float(w@(b*b*np.exp((1.-b)/theta)*kve(1,b/theta)))/(theta*kve(2,1/theta))
    binary_kl=floor*math.log(floor/thermal)+(1.-floor)*(math.log1p(-floor)-math.log1p(-thermal)) if floor>thermal else 0.
    return dict(global_pair_kernel_upper=max_pair_rate_bound(),earliest_birth_depth_upper=float(depth_upper[0]),
                mean_survival_endpoint_estimate=float(integral/rest[-1]),mean_survival_floor=floor,
                theta=float(theta),transverse_tube_radius=tube_radius,maxwell_juttner_tube_fraction=thermal,
                mandatory_to_thermal_tube_population_ratio=floor/thermal,
                total_variation_floor=max(0.,floor-thermal),relative_entropy_floor_nats=binary_kl)

def measure(history,end,exe,scratch,tag_power,target_power,steps,seed):
    tags=qmc.Sobol(d=2,scramble=True,seed=seed).random_base2(tag_power)
    targets=qmc.Sobol(d=3,scramble=True,seed=seed+1).random_base2(target_power)
    amount=history.value(end,'D')
    tags_b=history.births(amount*tags[:,0],tags[:,1])
    edges=np.linspace(history.data['t'][0],end,steps+1)
    dt=np.diff(edges); t=(edges[1:]+edges[:-1])/2
    a=history.value(t,'a'); rest=history.value(t,'D')
    impulse=history.value(t,'I')
    density=history.branch*history.action*rest/(6*history.mu**3*a**3)
    with open(scratch,'wb') as out:
        tags_b.tofile(out)
        np.cos(2*np.pi*targets[:,2]).tofile(out)
        for begin in range(0,steps,64):
            sl=slice(begin,min(steps,begin+64))
            masses=rest[sl,None]*targets[:,0]
            samples=history.births(masses,np.broadcast_to(targets[:,1],masses.shape))
            states=np.stack((t[sl],dt[sl],a[sl],impulse[sl],density[sl]),axis=1)
            np.column_stack((states,samples.reshape(len(states),-1))).tofile(out)
    response=subprocess.check_output([str(exe),str(scratch),str(steps),str(len(tags)),str(len(targets))],text=True)
    tau=np.fromstring(response,sep='\n')
    if len(tau)!=len(tags) or np.min(tau)<0 or not np.all(np.isfinite(tau)):raise ValueError('Invalid hazard output')
    a=history.value(end,'a'); I=history.value(end,'I')
    z=(tags_b[:,1]+I)/a; perp=tags_b[:,2]/a; E=np.sqrt(1+z*z+perp**2)
    survival=np.exp(-tau)
    return dict(tau=float(end),born_rest=float(amount),mean_depth=float(np.mean(tau)),
                mean_survival=float(np.mean(survival)),mean_one_or_more=float(np.mean(-np.expm1(-tau))),
                depth_quantiles=np.quantile(tau,[.1,.5,.9]).tolist(),
                survivor_kinetic_energy_fraction=float(np.mean(survival*(E-1))/np.mean(E-1)),
                ballistic_mean_v=float(np.mean(z/E)),fluid_v=float(history.value(end,'u')/np.hypot(1.,history.value(end,'u'))),
                ballistic_mean_gamma=float(np.mean(E)),transverse_pressure_share=float(np.mean(perp**2/E)/np.mean((perp**2+z*z)/E)))


def synthetic_checks(exe,folder):
    cases={'rest':(0.,0.,0.,0.,.4,1.,0.,0.),'resonant':(4150.,0.,4150.,0.,0.,1.,0.,0.),
           'translated_and_redshifted':(-40000.,.21,-700.,.18,0.,3.,2000.,.432),
           'parallel':(8500.,0.,-8500.,0.,0.,1.,0.,0.)}
    out={}
    for name,(z,q,z2,q2,birth,a,impulse,phi) in cases.items():
        path=folder/(name+'.bin')
        tag=[birth,z,q]; target=[0.,z2,q2]; dt=1.
        with path.open('wb') as f:
            np.asarray([*tag,phi,.5,dt,a,impulse,1.,*target]).tofile(f)
        measured=float(subprocess.check_output([str(exe),str(path),'1','1','1'],text=True))
        z=(z+impulse)/a; z2=(z2+impulse)/a; q=q/a; q2=q2/a
        product=math.sqrt(1+z*z+q*q)*math.sqrt(1+z2*z2+q2*q2)
        s=max(4.,2*(1+product+z*z2+q*q2*phi))
        beta=math.sqrt(1-4/s)
        # Compute the flux and cross section separately, including the finite threshold.
        g=.30*(s+2)
        sigma=12*math.pi*g*g/(s*((s-8300**2)**2+(beta*g)**2))
        expected=(1-birth)*sigma*.5*s*beta/product
        if not math.isclose(measured,expected,rel_tol=1e-10,abs_tol=1e-20):
            raise ValueError(f'Pair kernel check failed: {name}')
        out[name]={'measured':measured,'expected':expected}
    h=(np.arange(64)+.5)/64
    quant=inverse_angle(h,np.zeros_like(h),np.ones_like(h),np.ones_like(h))
    error=float(np.max(np.abs(quant-np.sign(2*h-1)*np.sqrt(np.abs(2*h-1)))))
    if error>1e-12:raise ValueError('Angular birth quantiles disagree')
    out['angular_quantile_max_abs_error']=error
    return out


def angular_controls():
    z,w=np.polynomial.legendre.leggauss(16)
    controls={}
    for b2 in (0.,.5,1.):
        A=3.-2*b2; B=b2*b2
        density=(A+B*z*z)/(2*(A+B/3))
        p2=(3*z*z-1)/2
        controls[str(b2)]=dict(normalization=float(w@density),P1=float(w@(density*z)),
                               P2=float(w@(density*p2)),expected_P2=2*B/(15*(A+B/3)))
    return controls


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--traces',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=2)
    args=p.parse_args()
    meta=json.loads((args.traces/'metadata.json').read_text())
    series={'coarse':(7,7,768,731),'fine':(8,8,1536,931),'repeat':(8,8,1536,1133)}
    series['dense']=(9,10,2048,3091)
    results={}
    with tempfile.TemporaryDirectory(prefix='hdblast-ballistic-') as folder:
        folder=Path(folder); exe=folder/'hazard'
        subprocess.run(['g++','-O3','-std=c++17',str(HERE/'hazard.cpp'),'-o',str(exe)],check=True)
        tests=synthetic_checks(exe,folder)
        histories={name:History(args.traces/(name+'.bin'),meta,run) for name,run in meta['runs'].items()}
        jobs=[(name,marker,res) for name in histories for marker in MARKERS for res in (series if name.endswith('native') else ('fine',))]
        def evaluate(job):
            name,marker,res=job; hist=histories[name]
            r=measure(hist,hist.events[marker],exe,folder/('_'.join(job)+'.dat'),*series[res])
            print(*job,r['mean_survival'],r['mean_depth'],flush=True)
            return job,r
        with ThreadPoolExecutor(args.jobs) as pool:
            for (name,marker,res),r in pool.map(evaluate,jobs):
                results.setdefault(name,{}).setdefault(marker,{})[res]=r
    report={'parameters':{'daughter_mass_L':meta['daughter_mass_L'],'action':meta['action'],
                          'charge_to_mass':meta['charge_to_mass'],'pole_mass_over_m':8300.,'width_over_mass':.3},
            'any_target_bounds':{name:{marker:bath_independent_floor(hist,end) for marker,end in hist.events.items()}
                                 for name,hist in histories.items()},
            'kernel_controls':tests,'angular_controls':angular_controls(),'meshes':{k:dict(tags=2**a,targets=2**b,time_midpoints=c,sobol_seed=d)
                                                          for k,(a,b,c,d) in series.items()},
            'runs':results}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
