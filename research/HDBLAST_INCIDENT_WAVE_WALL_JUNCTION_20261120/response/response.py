"""Couple a transverse wave on a tuned bulk to a passive, pressure-driven boundary."""
import argparse
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('fixed_scattering',ROOT/'scattering/scattering.py')
scattering = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scattering)


@lru_cache(maxsize=12)
def geometry(law,dx,length,map_width):
    n=round(length/dx)+1
    x=np.linspace(0.,length,n)
    mid=(x[:-1]+x[1:])/2
    widths=np.full(n,dx)
    widths[[0,-1]]=dx/2
    knots=np.r_[np.linspace(-.95,0.,200,endpoint=False),np.linspace(0.,length,3*n)]
    def y_rhs(z,y):
        return np.exp(-scattering.original().primitive(y)/3)
    def trace(grid):
        sol=solve_ivp(y_rhs,(grid[0],grid[-1]),[0.],t_eval=grid,method='DOP853',
                      rtol=3e-13,atol=1e-14)
        if not sol.success:raise ArithmeticError(sol.message)
        return sol.y[0]
    y=np.r_[trace(knots[knots<=0][::-1])[::-1][:-1],trace(knots[knots>=0])]
    warp=-scattering.original().primitive(y)/3
    logw=3*warp if law=='selected' else warp-4*(np.logaddexp(y,-y)-math.log(2))
    weight=CubicSpline(knots,logw)
    mass=widths*np.exp(weight(x))
    spring=np.exp(weight(mid))/dx
    shape=np.maximum(1-x/map_width,0.)
    midshape=np.maximum(1-mid/map_width,0.)
    nr=np.count_nonzero(shape)
    er=np.count_nonzero(midshape)
    return (x,mid,widths,mass,spring,weight,shape[:nr]**4,-4*shape[:nr]**3/map_width,
            midshape[:er]**4,-4*midshape[:er]**3/map_width)


def evolve(law, incident, dx=.05, *, length=60., end=36., fixed=False,
           damping=.8, frequency=1., map_width=3., rtol=2e-8, snapshots=False):
    k=1.
    x,mid,widths,mass,spring,weight,h,hprime,hm,hprime_m=geometry(law,dx,length,map_width)
    # Use an integral number of cells throughout the run.
    dx_actual=x[1]-x[0]
    if abs(dx_actual-dx)>1e-10:raise ValueError('Choose compatible length and grid spacing')
    n=len(x)
    nr,er=len(h),len(hm)
    seed=scattering.prepare_packet(law,'balanced')
    initial_q,initial_vel=scattering.initial_data(x,seed)
    support=scattering.shapes(x)[0][0]
    # Keep both continuum constraints exact also on the numerical quadrature.
    for data in (initial_q,initial_vel):
        data-=((mass@data)/(mass@support))*support
    initial_E=.5*(mass@(initial_vel**2+k*k*initial_q**2)+spring@np.diff(initial_q)**2)
    scale=math.sqrt(incident/initial_E)
    q0,v0=initial_q*scale,initial_vel*scale
    N=seed['half_norm']
    def coefficients(b):
        wz=np.exp(weight(x[:nr]+b*h))
        ratio=widths[:nr]*wz*(1+b*hprime)/mass[:nr]
        perturbed=np.exp(weight(mid[:er]+b*hm))/((1+b*hprime_m)*dx)
        return wz,ratio,perturbed
    def rhs(t,state):
        q=state[:n]; v=state[n:2*n]
        b,vel,heat,rc,rs,test_b,test_v,test_heat,tc,ts=state[2*n:]
        wz,mratio,spring_near=coefficients(b)
        pi=v.copy()
        pi[:nr]/=mratio
        dq=np.diff(q)
        flux=spring*dq
        flux[:er]=spring_near*dq[:er]
        mixed=np.zeros(nr+1)
        mixed[:nr]=wz*vel*h*pi[:nr]
        flux[:nr]+=.5*(mixed[:-1]+mixed[1:])
        qt=pi.copy()
        grad=np.empty(nr)
        grad[0]=dq[0]/dx
        grad[1:]=.5*(dq[:nr-1]+dq[1:nr])/dx
        qt[:nr]+=vel*h*grad/(1+b*hprime)
        vt=np.empty(n)
        vt[0]=flux[0];vt[-1]=-flux[-1];vt[1:-1]=np.diff(flux)
        vt/=mass
        vt-=k*k*q
        vt[:nr]-=(mratio-1)*k*k*q[:nr]
        pressure=.5*wz[0]*((1-vel*vel)*pi[0]**2-k*k*q[0]**2)
        acc=0. if fixed else -frequency**2*b-damping*vel-pressure
        # In the fixed branch, test_b supplies the first displacement driven by
        # the fixed-wall pressure. The complex moment uses a rotating phase.
        test_acc=-frequency**2*test_b-damping*test_v-pressure
        phase_real=-math.sin(k*t)*k*wz[0]*q[0]
        phase_imag=-math.cos(k*t)*k*wz[0]*q[0]
        return np.r_[qt,vt,vel,acc,damping*vel*vel,phase_real*vel,phase_imag*vel,
                     test_v,test_acc,damping*test_v**2,phase_real*test_v,phase_imag*test_v]
    sample=np.linspace(0.,end,round(end/.1)+1)
    sol=solve_ivp(rhs,(0.,end),np.r_[q0,v0,np.zeros(10)],method='DOP853',
                  t_eval=sample,max_step=min(.09,2.3*dx),rtol=rtol,atol=rtol*1e-4)
    if not sol.success:raise ArithmeticError(sol.message)
    max_res=0.; moment_error=0.; far=0.; summary=[]
    for index,t in enumerate(sample):
        state=sol.y[:,index]
        q=state[:n];v=state[n:2*n]
        b,vel,heat,rc,rs,test_b,test_v,test_heat,tc,ts=state[2*n:]
        wz,mratio,spring_near=coefficients(b)
        m=mass.copy();m[:nr]*=mratio
        stiffness=spring.copy();stiffness[:er]=spring_near
        pi=v.copy();pi[:nr]/=mratio
        energy=.5*(m@(pi*pi+k*k*q*q)+stiffness@np.diff(q)**2)
        oscillator=.5*(vel*vel+frequency**2*b*b)
        residual=abs(energy+oscillator+heat-incident)/incident
        max_res=max(max_res,residual)
        far=max(far,abs(q[-1]),abs(pi[-1]))
        # Integral tails of q and pi vanish outside the light cone of the seed.
        # Restore the half-line zero-mode norm beyond the finite checking grid.
        norm=N+np.sum(m[:nr]-mass[:nr])
        moment=np.exp(-1j*k*t)*(mass@v+1j*k*(m@q))
        moment_error=max(moment_error,abs(moment-(rc+1j*rs)))
        pole=abs(moment)**2/(2*norm)
        if index%50==0 or index==len(sample)-1:
            summary.append(dict(time=float(t),wall_position=float(b),wall_speed=float(vel),
                         field_energy=float(energy),wall_mechanical_energy=float(oscillator),
                         dissipated_energy=float(heat),projected_energy=float(pole),
                         moment_source_energy=float((rc*rc+rs*rs)/(2*norm))))
    if max(abs(sol.y[2*n]))>=.6 or max(abs(sol.y[2*n+1]))>=.5 or far>1e-6:
        raise ArithmeticError('Boundary speed, coordinate range or domain is inadequate')
    final=summary[-1]
    return dict(incident=incident,spacing=dx,duration=end,frequency=frequency,damping=damping,
                fixed=fixed,peak_absolute_displacement=float(max(abs(sol.y[2*n]))),
                peak_wall_speed=float(max(abs(sol.y[2*n+1]))),final=final,
                peak_relative_balance_error=float(max_res),
                peak_absolute_moment_difference=float(moment_error),far_field_amplitude=float(far),
                **(dict(linear_dissipated_coefficient=float(test_heat/incident**2),
                        linear_projected_coefficient=float((tc*tc+ts*ts)/(2*N*incident**3))) if fixed else {}),
                **({'times':summary} if snapshots else {}))


def compatibility_check():
    path = Path(__file__).with_name('compatibility.py')
    definition = importlib.util.spec_from_file_location('frozen_wall_compatibility', path)
    module = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(module)
    return module.calculate()


def calculate():
    results={}
    for law in ('wall','selected'):
        runs={}
        for dx in (.05,.025):
            for e in (.01,.04,.16):
                runs[f'{dx}/{e}']=evolve(law,e,dx,snapshots=(dx==.025 and e==.04))
            runs[f'{dx}/fixed']=evolve(law,1.,dx,fixed=True)
        controls={'finest':evolve(law,.04,.00625),'fine':evolve(law,.04,.0125), 'coarse':evolve(law,.04,.1),
                  'later':evolve(law,.04,.025,end=44.), 'other_map':evolve(law,.04,.025,map_width=4.5),
                  'undamped':evolve(law,.04,.025,damping=0.), 'weak':evolve(law,.0025,.025)}
        primary=runs['0.025/0.04']
        half=primary['final']
        finer=controls['fine']['final']
        finest=controls['finest']['final']
        diagnostics={}
        for key in ('dissipated_energy','projected_energy','moment_source_energy'):
            new=(4*finest[key]-finer[key])/3
            old=(4*finer[key]-half[key])/3
            diagnostics[key]=dict(extrapolated=new,previous_extrapolated=old,
                                  relative_grid_shift=abs(new/old-1),
                                  relative_map_change=abs(controls['other_map']['final'][key]/half[key]-1),
                                  relative_time_change=abs(controls['later']['final'][key]/half[key]-1))
        base=runs['0.025/fixed']
        weak=controls['weak']
        pred_d=base['linear_dissipated_coefficient']*weak['incident']**2
        pred_p=base['linear_projected_coefficient']*weak['incident']**3
        tests=dict(weak_dissipation_over_leading=weak['final']['dissipated_energy']/pred_d,
                   weak_projection_over_leading=weak['final']['projected_energy']/pred_p,
                   fine_source_over_projection=finest['moment_source_energy']/finest['projected_energy'])
        if any(abs(v-1)>.025 for v in tests.values()):
            raise ArithmeticError(('Weak-amplitude or moment check',law,tests))
        if (max(v['relative_grid_shift'] for v in diagnostics.values())>.006 or
                max(v['relative_map_change'] for v in diagnostics.values())>.012):
            raise ArithmeticError(('Mesh or map control failed',law,diagnostics))
        results[law]=dict(runs=runs,controls=controls,refinement=diagnostics,tests=tests)
    return dict(model='fixed bulk with unit effective wall mass, positive spring and passive damping in conformal coordinates',
                branches=results, compatibility=compatibility_check())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(),indent=2,allow_nan=False)+'\n')
