"""Finite transverse photon oscillators with dynamical parent and closed Einstein gravity."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp

ROOT = Path(__file__).resolve().parents[1]
NORM = quad(lambda z: math.exp(-1/(z*(1-z))),0.,1.,epsabs=2e-16)[0]


def load_prior(prior):
    path=prior/'homogeneous/exchange.py'
    spec=importlib.util.spec_from_file_location('finite_parent_model',path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.branches(prior/'prior/gravity')


def profile(z):
    if z <= 0 or z >= 1:
        return 0.
    return math.exp(-1/(z*(1-z)))/NORM


def run(branch,limit,delta,*,coupled=True,coordinates=False,refine=False):
    q0=branch.fitted.sign*branch.fitted.dq
    qmin=branch.fitted.tmin-branch.top
    distance=qmin-q0
    ui=float(branch.potential(q0)[0])
    ai=1/math.sqrt(ui)
    ls=np.arange(2,limit+1,dtype=float)
    weights=ls*(ls**2-1)
    coeff=branch.fitted.Hpl**2/(3*math.pi**2)
    count=len(ls)
    def fields(y):
        n=y[7:7+count]
        real=y[7+count:7+2*count]
        imag=y[7+2*count:]
        if coordinates:
            # Real canonical covariances X=2<I²Lq²>,Y=2<p²/(I²L)>,Z=<pq+qp>.
            # Use shifted variances to keep a weak excitation numerically visible.
            return .25*(n+real),.25*(n-real),.5*imag
        return n,real,imag
    def terms(y):
        q,v,log_a,h=y[:4]
        a=ai*math.exp(log_a)
        if not branch.fitted.theta[0] <= branch.top+q <= branch.fitted.theta[-1]:
            raise ArithmeticError('Field left the archived potential grid')
        u,force=branch.potential(q)
        n,real,imag=fields(y)
        density_factor=coeff/a**4
        rho=density_factor*(weights@n)
        grad=delta*profile(((q-q0)/distance-.2)/.6)/(.6*distance)
        exchange=2*density_factor*grad*(weights@real)
        return a,float(u),float(force),rho,grad,exchange,n,real,imag
    def rhs(x,y):
        q,v,log_a,h=y[:4]
        a,u,force,r,grad,per_velocity,n,real,imag=terms(y)
        s=grad*v
        omega=ls/a
        if coordinates:
            # X'=2s X+2ω Z, Y'=−2s Y−2ω Z, Z'=ω(Y−X).
            # Z=2 Im <bb>, and each diagonal variable has vacuum value one.
            X,Y,Z=y[7:7+count],y[7+count:7+2*count],y[7+2*count:]
            wave=[2*s*(1+X)+2*omega*Z,-2*s*(1+Y)-2*omega*Z,omega*(Y-X)]
        else:
            wave=[2*s*real,2*omega*imag+s*(2*n+1),-2*omega*real]
        ray=1/a**2-v*v/2-(2*r if coupled else 0.)
        return np.concatenate(([v,force-3*h*v-(3*per_velocity if coupled else 0.),h,ray,
                                a**4*v*per_velocity,h*v*v+(v*per_velocity if coupled else 0.),s],
                               *wave))
    def crossing(x,y): return (y[0]-q0)/distance-1
    crossing.direction=1
    crossing.terminal=True
    def switch_end(x,y): return (y[0]-q0)/distance-.8
    switch_end.direction=1
    def switch_start(x,y): return (y[0]-q0)/distance-.2
    switch_start.direction=1
    def midpoint(x,y): return (y[0]-q0)/distance-.5
    midpoint.direction=1
    def turn(x,y): return y[3]
    turn.direction=-1
    turn.terminal=True
    def direction_change(x,y): return y[1]/distance
    direction_change.direction=0
    def outside(x,y):
        theta=branch.top+y[0]
        return min(theta-branch.fitted.theta[0],branch.fitted.theta[-1]-theta)-.02
    outside.direction=-1
    outside.terminal=True
    y0=np.zeros(7+3*count)
    y0[0]=q0
    ans=solve_ivp(rhs,(0.,5.),y0,method='DOP853',rtol=2e-11 if refine else 2e-10,
                  atol=5e-15 if refine else 1.5e-13,max_step=(.004 if refine else .008),
                  events=[crossing,switch_end,switch_start,turn,outside,midpoint,direction_change],dense_output=True)
    if not ans.success or len(ans.t_events[4]) or not (len(ans.t_events[0]) or len(ans.t_events[3])):
        raise ArithmeticError('Trajectory ended outside the specified first-event test')
    def log_i(q):
        position=((q-q0)/distance-.2)/.6
        if position<=0:
            return 0.
        if position>=1:
            return delta
        return delta*quad(profile,0.,position,epsabs=5e-14,epsrel=5e-14)[0]
    extrema=[x for x in ans.t_events[6] if x>1e-9]
    endpoints=[0.]+[log_i(ans.sol(x)[0]) for x in extrema]
    variations=np.r_[0.,np.cumsum(np.abs(np.diff(endpoints)))]
    def record(x):
        y=ans.sol(x)
        a,u,force,r,grad,per_v,n,real,imag=terms(y)
        q,v,log_a,h,C,W=y[:6]
        parent=u+v*v/6
        value=log_i(q)
        i=np.searchsorted(extrema,x,side='right')
        variation=variations[i]+abs(value-endpoints[i])
        if abs(value-y[6])>2e-8:
            raise ArithmeticError('Coefficient differs from its signed trajectory integral')
        return dict(x=float(x),a=a,h=h,q=q,v=v,progress=(q-q0)/distance,photon_density=r,
                    parent_density=parent,photon_over_parent=r/parent,
                    fraction_of_matter=r/(r+parent),comoving_excess=a**4*r,
                    integrated_source_comoving=C,total_variation_log_I=float(variation), constraint=h*h+1/a**2-parent-(r if coupled else 0.),
                    source_residual=a**4*r-C,parent_residual=parent-ui+W,
                    purity_residual=float(np.max(abs(real**2+imag**2-n*(n+1))/(1+abs(n*(n+1))))),
                    highest_proper_frequency_over_Mpl=limit*branch.fitted.Hpl/a)
    end=record(ans.t[-1])
    stages={key:(record(float(ev[0])) if len(ev) else None) for key,ev in zip(
        ['first_crossing','switch_end','switch_start','first_maximum'],ans.t_events[:4])}
    if len(ans.t_events[5]):
        mid=float(ans.t_events[5][0])
        stages['midpoint']=record(mid)
        y=ans.sol(mid)
        a,u,force,r,grad,per_v,n,real,imag=terms(y)
        high=ls>=.8*limit
        eta_s=a*grad*y[1]
        estimate=eta_s**2*np.sum(weights[high]/ls[high]**2)/4
        stages['midpoint']['high_band_energy_over_adiabatic_leading']=float((weights[high]@n[high])/estimate) if estimate>0 else None
    else:
        stages['midpoint']=None
    times=np.linspace(0.,ans.t[-1],151)
    sampled=[record(x) for x in times]
    errors={k:max(abs(r[k]) for r in sampled+[end]) for k in ('constraint','source_residual','parent_residual','purity_residual')}
    max_var=end['total_variation_log_I']
    bound=coeff*sum(weights)*math.sinh(max_var)**2/(ai**2)
    if coupled and any(row['fraction_of_matter']>bound+2e-9 for row in sampled):
        raise ArithmeticError('Finite band exceeds the integrated squeeze upper bound')
    if max(errors.values())>4e-8:
        raise ArithmeticError('Finite oscillator or background balance failed: '+str(errors))
    return dict(cutoff=limit,delta=delta,coupled=coupled,
                field_reversals_before_stop=[dict(x=float(x),progress=float((ans.sol(x)[0]-q0)/distance)) for x in ans.t_events[6] if x>1e-9],fraction_ceiling_from_total_variation=bound,stages=stages,end=end,checks=errors,
                sampled_peak_fraction=max(r['fraction_of_matter'] for r in sampled),
                sampled_peak_fraction_during_switch=max((r['fraction_of_matter'] for r in sampled if .2<=r['progress']<=.8),default=0.))


def calculate(prior,quick=False):
    spec=importlib.util.spec_from_file_location('finite_budget',Path(__file__).with_name('budget.py'))
    budget=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(budget)
    coast_spec=importlib.util.spec_from_file_location('photon_free_interval',Path(__file__).with_name('coast.py'))
    coast=importlib.util.module_from_spec(coast_spec)
    coast_spec.loader.exec_module(coast)
    outputs={}
    for name,b in load_prior(prior).items():
        rows={}
        for n in ([32,64] if quick else [16,32,48,64,96,128,192]):
            for delta in ([.1,.5] if quick else [.01,.1,.25,.5,1.]):
                rows[f'{n}/{delta:g}']=run(b,n,delta)
        if not quick:
            for n in (16,24,32,48,64,96,128,192):
                for delta in (1.5,2.,3.):
                    rows[f'{n}/{delta:g}']=run(b,n,delta)
            # Selected after scanning for an off-switch fraction near one-half.
            rows['32/1.535']=run(b,32,1.535)
        control={}
        for delta in (0.,.1,.5):
            control[f'passive/{delta:g}']=run(b,64,delta,coupled=False)
        primary=rows['64/0.1']
        control['alternative_covariance']=run(b,64,.1,coordinates=True,refine=True)
        control['inverse_coupling']=run(b,64,-.1)
        for key in ('alternative_covariance','inverse_coupling'):
            alternate=control[key]
            for field in ('x','a','v','photon_density','parent_density','fraction_of_matter'):
                if abs(primary['end'][field]-alternate['end'][field])>2e-8:
                    raise ArithmeticError('Independent photon covariance or duality mismatch')
        if not quick:
            control['strong_alternative_covariance']=run(b,48,3.,coordinates=True,refine=True)
            check=control['strong_alternative_covariance']['end']
            ref=rows['48/3']['end']
            if max(abs(check[field]-ref[field]) for field in ('x','a','v','photon_density','fraction_of_matter'))>2e-8:
                raise ArithmeticError('Strong coupled calculation differs in canonical coordinates')
        ai=1/math.sqrt(float(b.potential(b.fitted.sign*b.fitted.dq)[0]))
        free_intervals={}
        for key in ('48/3','128/3','48/0.1','32/1.5','32/1.535'):
            if key in rows:
                stages=rows[key]['stages']
                free_intervals[key]=coast.coast(b,stages['switch_end'],stages['first_crossing'])
        outputs[name]=dict(runs=rows,controls=control,free_intervals=free_intervals,
                           finite_band_envelope=budget.budgets(ai,b.fitted.Hpl))
    return dict(model='finite conformal harmonics 2 through cutoff, both polarizations; instantaneous oscillator zero point removed at all times',
                branches=outputs)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior',type=Path,default=ROOT/'prior')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--quick',action='store_true')
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior,args.quick),indent=2,allow_nan=False)+'\n')
