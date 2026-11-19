"""Closed, homogeneous Einstein--two-scalar evolution in the archived hilltop units."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent


def branches(gravity):
    spec = importlib.util.spec_from_file_location('constraint', gravity / 'initial/constraint.py')
    constraint = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(constraint)
    return constraint.load_branches(gravity / 'previous')


def simulate(branch, offset, seed, coupled, *, direct=False, field_cap=None, horizon=5.):
    """seed = initial daughter density divided by initial scalar potential.
    q=Theta-Theta_top, v=dq/dx, N=log(a/a_i), h=H/H0, r=rho_chi/V0,
    V0=3 Mpl^2 H0^2. Direct check evolves the daughter velocity and a''.
    """
    pot = branch.fitted
    ftop = pot.spline(branch.top)[1]

    def terms(q):
        theta = branch.top + q
        u, du, _ = pot.at(theta)
        u -= branch.du_top*q
        du -= branch.du_top
        if q == 0.:
            du = 0.  # The fitted hilltop is stationary by construction.
        if coupled:
            f = pot.spline(theta)[1]
            df = pot.spline(theta, 1)[1]
            return u, du, ftop/f, -df/f
        return u, du, 1., 0.

    ui, dui, zi, betai = terms(offset)
    ri = seed*ui
    ai = 1/np.sqrt(ui+ri)
    # P here is the dimensionless conserved charge a^3 Z d(chi/Mpl)/dx.
    charge = ai**3*np.sqrt(6*ri*zi)

    def observables(y):
        q, v, third, fourth, w = y[:5]
        if direct:
            a, h = third, fourth/third
            n = np.log(a/ai)
        else:
            n, h = third, fourth
            a = ai*np.exp(n)
        u, du, z, beta = terms(q)
        r = z*w*w/6 if direct else ri*np.exp(-6*n)*zi/z
        source = -r*beta*v
        return a,h,n,u,du,z,beta,r,source

    def rhs(x, y):
        q,v,third,fourth,w = y[:5]
        a,h,n,u,du,z,beta,r,source = observables(y)
        vp = -3*h*v - 3*du + 3*r*beta
        if direct:
            ap = fourth
            app = a*(u-v*v/3-2*r)
            wp = -(3*h+beta*v)*w
            return v,vp,ap,app,wp,source,6*h*r, h*v*v, np.exp(6*n)*source
        return (v,vp,h,1/a**2-v*v/2-3*r,0.,
                source,6*h*r,h*v*v,np.exp(6*n)*source)

    def maximum(x,y):
        return y[3]
    maximum.direction = -1.
    maximum.terminal = True

    def crossing(x,y):
        return pot.sign*(branch.top+y[0]-pot.tmin)
    crossing.direction = 1.

    def outside(x,y):
        th = branch.top+y[0]
        return min(th-pot.theta[0],pot.theta[-1]-th)-.05
    outside.terminal = True
    outside.direction = -1.

    initial = [offset,0.,ai if direct else 0.,0.,np.sqrt(6*ri/zi) if direct else 0.,0.,0.,0.,0.]
    def cap(x,y):
        return field_cap-abs(y[0]) if field_cap is not None else 1.
    cap.terminal = True
    cap.direction = -1.

    result = solve_ivp(rhs,(0.,horizon),initial,events=(maximum,crossing,outside,cap),
                       method='DOP853' if not direct else 'RK45',
                       max_step=.015 if not direct else .006, dense_output=True,
                       rtol=2e-11 if not direct else 1e-12,
                       atol=1e-14 if not direct else 3e-15)
    if not result.success:
        raise ValueError(result.message)

    def record(x):
        y = result.sol(x)
        q,v,_,_,w,int_source,int_dilution,int_friction,int_weighted = y
        a,h,n,u,du,z,beta,r,source = observables(y)
        constraint = h*h + 1/a**2 - u-v*v/6-r
        return dict(x=x,q=q,theta=branch.top+q,v=v,a=a,h=h,N=n,u=u,
                    r=r,r_fraction=r/(u+v*v/6+r),Z=z,dlogZ_dq=beta,
                    source_to_daughter=source,dilution_rate=6*h*r,
                    daughter_density_rate=source-6*h*r,
                    integrated_source=int_source,integrated_dilution=int_dilution,
                    comoving_stiff_gain_factor=zi/z if ri else None,
                    Hprime=1/a**2-v*v/2-3*r,acceleration_over_a=u-v*v/3-2*r,
                    residuals=dict(friedmann=constraint,
                        daughter=r-ri-int_source+int_dilution,
                        parent=u+v*v/6-ui+int_friction+int_source,
                        weighted_daughter=r*np.exp(6*n)-ri-int_weighted,
                        direct_charge=a**3*z*w-charge if direct else None))

    samples = [record(float(x)) for x in np.linspace(0.,result.t[-1],251)]
    audit = {k:max(abs(s['residuals'][k]) for s in samples) for k,v in samples[0]['residuals'].items() if v is not None}
    ending = next((tag for i,tag in ((0,'first_maximum'),(2,'fit_cut'),(3,'field_cap'))
                   if len(result.t_events[i])), 'time_limit')
    if result.t[-1] < .12:
        raise ValueError('Probe time exceeds the computed trajectory.')
    return dict(offset=offset,seed_r_over_ui=seed,coupled=coupled,charge=charge,
                initial_parent_acceleration=-3*dui+3*ri*betai,initial=record(0.),
                probe=record(.12), first_minimum_crossing=record(float(result.t_events[1][0])) if len(result.t_events[1]) else None,
                stop=ending,at_stop=record(float(result.t[-1])),maximum_sampled_abs_residual=audit)


def calculate(gravity):
    data = {}
    for name,branch in branches(gravity).items():
        pot = branch.fitted
        runs = {}
        for seed in (0.,1e-5,1e-3,1e-2):
            for coupled in ((True,) if seed == 0 else (True, False)):
                key = f'{"inverse_f" if coupled else "constant"}_seed_{seed:g}'
                runs[key] = simulate(branch,pot.sign*pot.dq,seed,coupled)
        direct = simulate(branch,pot.sign*pot.dq,.01,True,direct=True)
        reduced = runs['inverse_f_seed_0.01']
        err = {stage:{k:abs(direct[stage][k]-reduced[stage][k])
                      for k in ('x','q','v','a','h','r','integrated_source')}
               for stage in ('probe','first_minimum_crossing','at_stop') if reduced[stage] and direct[stage]}
        if any(r['stop']!='first_maximum' for r in runs.values()):
            raise ValueError('A run failed to reach the first maximum within the fit.')
        if max(max(run['maximum_sampled_abs_residual'].values()) for run in runs.values())>2e-8:
            raise ArithmeticError('Constraint or component balance failed.')
        if max(max(row.values()) for row in err.values())>1e-8 or max(direct['maximum_sampled_abs_residual'].values())>2e-8:
            raise ArithmeticError('Independent trajectory check failed.')
        crest = {label:simulate(branch,0.,charge_seed,variable,field_cap=.01,horizon=.8)
                 for label,charge_seed,variable in [('charged',.01,True),('constant_metric',.01,False),('uncharged',0.,True)]}
        if (crest['charged']['stop'] != 'field_cap' or crest['charged']['at_stop']['q'] <= 0
            or any(crest[key]['at_stop']['q'] != 0 for key in ('constant_metric','uncharged'))):
            raise ArithmeticError('Exact crest response controls failed.')
        crossing = runs['inverse_f_seed_0.01']['first_minimum_crossing']
        if crossing is None or crossing['integrated_source']*pot.pm<=0:
            raise ArithmeticError('First-passage source sign changed.')
        data[name]=dict(theta_top=branch.top,theta_min=pot.tmin,H0_over_Mpl=pot.Hpl,
                        seed_parent_abs_displacement=pot.dq,archived_detuning=.001,
                        ftop=pot.spline(branch.top)[1],runs=runs,exact_crest_controls=crest,
                        direct_velocity_scale_factor_check=dict(differences=err,
                            residuals=direct['maximum_sampled_abs_residual']))
    return dict(model='closed two-derivative four-dimensional EFT; homogeneous canonical parent and shift-symmetric daughter',
                normalization='x=H0*t; a=H0*a_physical; fields in reduced Planck units; densities in V0=3*Mpl^2*H0^2; Z=ftop/f or 1',
                branches=data)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gravity',type=Path,default=HERE.parent/'prior/gravity')
    parser.add_argument('--output',type=Path,default=HERE/'results.json')
    args = parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.gravity),indent=2,allow_nan=False)+'\n')
