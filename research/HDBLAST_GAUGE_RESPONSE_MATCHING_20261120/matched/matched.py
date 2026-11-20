"""Conservative first-descent calculation with the leading finite-band inertia matched."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
NORM = quad(lambda u: math.exp(-1/(u*(1-u))), 0., 1., epsabs=2e-16)[0]
KREF = 48
RESERVE = .05


def branches(prior):
    spec=importlib.util.spec_from_file_location('calibrated_parent', prior/'homogeneous/exchange.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.branches(prior/'prior/gravity')


def switch(q,q0,dq,delta):
    width=.6*dq
    z=((q-q0)/dq-.2)/.6
    if z<=0. or z>=1.:
        return 0.,0.
    f=z*(1.-z)
    g=delta*math.exp(-1/f)/(width*NORM)
    return g,g*(1-2*z)/(width*f*f)


def run(branch,K,delta,*,match=True,direct=False,reserve=RESERVE,profile=switch,high_band=None):
    """Match a finite sharp band to K=48. Stop while bare inertia still positive."""
    q0=branch.fitted.sign*branch.fitted.dq
    dq=branch.fitted.tmin-branch.top-q0
    ui=float(branch.potential(q0)[0])
    ai=1/math.sqrt(ui)
    eps=branch.fitted.Hpl**2
    coef=eps/(3*math.pi**2)
    ls=np.arange(2,K+1,dtype=float)
    w=ls*(ls**2-1.)
    if high_band is not None:
        if match or K<=KREF:
            raise ValueError('A compressed high band requires K>48 and match=False')
        nodes,weights=(np.asarray(part,dtype=float) for part in high_band)
        if (nodes.ndim!=1 or weights.shape!=nodes.shape or not nodes.size
                or not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(weights))
                or np.any(nodes<=KREF) or np.any(nodes>K) or np.any(np.diff(nodes)<=0)
                or np.any(weights<=0)):
            raise ValueError('Invalid quadrature nodes or energy weights')
        low=ls<=KREF
        ls=np.r_[ls[low],nodes]
        w=np.r_[w[low],weights]
    high=ls[ls>KREF]
    c=eps/(2*math.pi**2)*math.fsum(high-1/high) if match else 0.
    count=len(ls)
    retained=ls<=KREF
    def unpack(y):
        q=y[0]
        a=y[2] if direct else ai*math.exp(y[2])
        h=y[3]/a if direct else y[3]
        g,gq=profile(q,q0,dq,delta)
        z=1.-c*g*g/a**2
        # Positive extension for trial Runge--Kutta stages below the terminal surface.
        # Only z>=reserve contributes to the saved solution.
        safe=max(.5*reserve,z)
        v=y[1] if direct else y[1]/(a**3*safe)
        az_a=2*(1.-z)
        zq=-2*c*g*gq/a**2
        first=y[6:6+count]
        second=y[6+count:6+2*count]
        third=y[6+2*count:]
        if direct:
            n=(first+second)*.25
            real=(first-second)*.25
            imag=third*.5
        else:
            n,real,imag=first,second,third
        r=coef*(w@n)/a**4
        source=2*coef*g*(w@real)/a**4
        return q,v,a,h,g,z,safe,az_a,zq,n,real,imag,r,source
    def derivatives(x,y):
        q,v,a,h,g,z,safe,az_a,zq,n,real,imag,r,S=unpack(y)
        if not branch.fitted.theta[0] <= branch.top+q <= branch.fitted.theta[-1]:
            raise ArithmeticError('Integration left the calibrated potential')
        u,force=branch.potential(q)
        s=g*v
        omega=ls/a
        if direct:
            # Physical scalar velocity and scale-factor acceleration, with canonical
            # photon quadratures X,Y shifted by their incoming unit vacuum variances.
            vp=(force-3*S-(3*z+az_a)*h*v-zq*v*v/2)/safe
            app=a*(u-(z/3+az_a/12)*v*v-r)
            dx=2*s*(1+y[6:6+count])+2*omega*y[6+2*count:]
            dy=-2*s*(1+y[6+count:6+2*count])-2*omega*y[6+2*count:]
            dz=omega*(y[6+count:6+2*count]-y[6:6+count])
            geometry=[v,vp,a*h,app]
            waves=(dx,dy,dz)
        else:
            pp=a**3*(force+zq*v*v/2-3*S)
            hp=1/a**2-(z/2+az_a/12)*v*v-2*r
            geometry=[v,pp,h,hp]
            waves=(2*s*real,2*omega*imag+s*(2*n+1),-2*omega*real)
        source=v*S
        work=h*(z+az_a/6)*v*v+source
        return np.concatenate((geometry,[a**4*source,work],*waves))
    def crossing(x,y): return (y[0]-q0)/dq-1.
    crossing.terminal=True
    crossing.direction=1
    def turning(x,y): return y[3]
    turning.terminal=True
    turning.direction=-1
    def lower_inertia(x,y):
        a=y[2] if direct else ai*math.exp(y[2])
        g,_=profile(y[0],q0,dq,delta)
        return 1.-c*g*g/a**2-reserve
    lower_inertia.terminal=True
    lower_inertia.direction=-1
    def switch_end(x,y): return (y[0]-q0)/dq-.8
    switch_end.direction=1
    y0=np.zeros(6+3*count)
    y0[0]=q0
    y0[2]=ai if direct else 0.
    sol=solve_ivp(derivatives,(0.,4.),y0,method='DOP853',rtol=2e-11 if direct else 2e-10,
                  atol=8e-15 if direct else 7e-14,max_step=.0025 if direct else .005,
                  dense_output=True,events=(crossing,turning,lower_inertia,switch_end))
    labels=['first_minimum','first_turnaround','kinetic_reserve']
    stamps=[(name,float(times[0])) for name,times in zip(labels,sol.t_events[:3]) if len(times)]
    if not sol.success or len(stamps)!=1: raise ArithmeticError('Missing terminal event: '+sol.message)
    state,xend=stamps[0]
    def obs(x):
        y=sol.sol(x)
        q,v,a,h,g,z,safe,az_a,zq,n,real,imag,r,S=unpack(y)
        u=float(branch.potential(q)[0])
        rho=u+z*v*v/6
        pressure=(z+az_a/3)*v*v/6-u
        low=coef*(w[retained]@n[retained])/a**4
        return dict(reference_band_fraction=low/(rho+r),reference_band_comoving=low*a**4,
                    added_band_comoving=(r-low)*a**4,x=float(x),q=q,progress=(q-q0)/dq,a=a,h=h,v=v,z=z,photon_density=r,
                    scalar_density=rho,scalar_pressure=pressure,photon_fraction=r/(rho+r),
                    comoving_energy=r*a**4, work_residual=r*a**4-y[4],
                    scalar_residual=rho-ui+y[5],constraint=h*h+1/a**2-rho-r,
                    purity_residual=float(max(abs(real**2+imag**2-n*(1+n))/(1+abs(n*(1+n))))))
    def inertia(x):
        y=sol.sol(x)
        a=y[2] if direct else ai*math.exp(y[2])
        g,_=profile(y[0],q0,dq,delta)
        return 1-c*g*g/a**2
    grid=np.linspace(0.,xend,161)
    snapshots=[obs(x) for x in grid]
    errors={field:max(abs(p[field]) for p in snapshots) for field in
            ['constraint','work_residual','scalar_residual','purity_residual']}
    if max(errors.values())>8e-7:
        raise ArithmeticError('Conservative or covariance check failed: '+str(errors))
    zs=[p['z'] for p in snapshots]
    minima=[float(zs[-1]),1.]
    if c>0.:
        for i in range(1,len(grid)-1):
            if zs[i]<zs[i-1] and zs[i]<=zs[i+1]:
                opt=minimize_scalar(inertia,bounds=(grid[i-1],grid[i+1]),method='bounded',
                                    options={'xatol':2e-13})
                if not opt.success: raise ArithmeticError('Missing inertia extremum')
                minima.append(float(opt.fun))
    if min(minima)<reserve-1e-7: raise ArithmeticError('Crossed the kinetic margin')
    return dict(cutoff=K,delta=delta,added_inertia_coefficient=c,status=state,
                minimum_bare_inertia=min(minima),end=obs(xend),errors=errors,
                switch_end=obs(sol.t_events[3][0]) if len(sol.t_events[3]) else None)


def calculate(prior):
    spec=importlib.util.spec_from_file_location('archived_coast',prior/'feedback/coast.py')
    free=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(free)
    data={}
    archived=json.loads((prior/'feedback/FEEDBACK.json').read_bytes())['branches']
    for key,b in branches(prior).items():
        rows={}
        for delta in (.01,.1,.25,.5,1.,3.):
            limits={48,49,52,64,96,128,192}
            if delta==3.:
                limits.update(range(49,55))
            if delta in (.01,.1):
                limits.update([256,384,512])
            if delta==.1:
                limits.update([640,768])
            for limit in sorted(limits):
                label=f'{limit}/{delta:g}'
                rows[label]=run(b,limit,delta)
        uncalibrated={}
        for limit in (48,64,96,128,192,256,384,512):
            label=f'{limit}/0.1'
            bare=run(b,limit,.1,match=False)
            if bare['status']!='first_minimum': raise ArithmeticError('Missing bare first crossing')
            uncalibrated[label]=bare
            if label in archived[key]['runs']:
                expected=archived[key]['runs'][label]['end']
                if max(abs(bare['end'][field]-expected[field]) for field in ('x','a','h','v','photon_density'))>2e-8:
                    raise ArithmeticError('Unmatched solver disagrees with its archive')
        for delta in (.01,.1,.25,.5,1.,3.):
            label=f'48/{delta:g}'
            prior_end=archived[key]['runs'][label]['end']
            if max(abs(rows[label]['end'][field]-prior_end[field]) for field in ('x','a','h','v','photon_density'))>2e-8:
                raise ArithmeticError('Reference changed')
        # Evaluate both sides of the first terminal cutoff with several positive margins.
        reserve_scan={}
        for limit in ((52,53,54) if b.fitted.sign<0 else (50,51,52)):
            for margin in (.02,.05,.10):
                record=rows[f'{limit}/3'] if margin==RESERVE else run(b,limit,3.,reserve=margin)
                reserve_scan[f'{limit}/{margin:g}']=dict(status=record['status'],minimum_inertia=record['minimum_bare_inertia'],
                                                       progress=record['end']['progress'],x=record['end']['x'])
        control={}
        last_admissible=53 if b.fitted.sign<0 else 51
        for limit,delta in [(48,3.),(49,3.),(last_admissible,3.),(last_admissible+1,3.),
                            (64,.1),(192,.01),(128,.5),(512,.1)]:
            alt=run(b,limit,delta,direct=True)
            base=rows[f'{limit}/{delta:g}']
            if alt['status']!=base['status']: raise ArithmeticError('Independent stop differs')
            fields=('x','v','progress','a','h','photon_density','photon_fraction','z')
            changes={k:abs(alt['end'][k]-base['end'][k]) for k in fields}
            if max(changes.values())>2e-6: raise ArithmeticError('Independent states differ')
            control[f'{limit}/{delta:g}']=changes
        coasts={}
        for limit in range(48,last_admissible+1):
            label=f'{limit}/3'
            trajectory=rows[label]
            reference=trajectory['end']
            expected=dict(reference,fraction_of_matter=reference['photon_fraction'])
            if trajectory['status']!='first_minimum' or trajectory['switch_end'] is None:
                raise ArithmeticError('Requested completed strong run is unavailable')
            coasts[label]=free.coast(b,trajectory['switch_end'],expected)
        data[key]=dict(runs=rows,uncalibrated=uncalibrated,reserve_scan=reserve_scan,
                       free_intervals=coasts,independent_coordinate_differences=control)
    return dict(reference_cutoff=KREF,reserve=RESERVE,branches=data)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=ROOT/'prior')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior),indent=2,allow_nan=False)+'\n')
