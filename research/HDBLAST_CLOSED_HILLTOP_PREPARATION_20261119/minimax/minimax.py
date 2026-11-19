"""Sharp all-time spatial RMS in unstable closed-sphere harmonics, on the linear hilltop."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar


class Transfer:
    """A starts at (1,0); B starts at (0,1). All selected modes are unstable throughout."""
    def __init__(self, h, mu2, n, stop):
        self.h, self.mu2, self.n = h, mu2, n
        if mu2 <= n * (n + 2) * h**2:
            raise ValueError('Maximum-principle condition fails at the bounce')
        def rhs(x, values):
            log_a, rate, j = values
            ch = math.cosh(h*x)
            pos = mu2 - n*(n+2)*h**2 / ch**2
            return [rate, pos - 3*h*math.tanh(h*x)*rate - rate**2,
                    math.exp(-2*log_a)/ch**3]
        sol = solve_ivp(rhs, (0,stop), [0.,0.,0.], dense_output=True, method='DOP853',
                        max_step=.08, rtol=2e-12, atol=2e-14)
        if not sol.success:
            raise ArithmeticError(sol.message)
        self.solution = sol.sol

    def at(self, time):
        log_amplitude, _, ratio = self.solution(time)
        a = math.exp(log_amplitude)
        return a, a * ratio

    def endpoint_basis_check(self, end):
        # Positive u with u(0)=1,u(T)=0 is integrated backwards, avoiding cancellation.
        def rhs(s,y):
            x = end-s
            force = self.mu2 - self.n*(self.n+2)*self.h**2 / math.cosh(self.h*x)**2
            return [y[1], 3*self.h*math.tanh(self.h*x)*y[1] + force*y[0]]
        backward = solve_ivp(rhs,(0,end),[0.,1.],dense_output=True,method='DOP853',
                             max_step=.08,rtol=3e-12,atol=2e-14)
        if not backward.success:
            raise ArithmeticError(backward.message)
        t = np.linspace(0,end,401)
        u = backward.sol(end-t)[0]/backward.y[0,-1]
        a_end,b_end = self.at(end)
        w = np.array([self.at(x)[1]/b_end for x in t])
        # Reconstruct A(t) = u(t) + w(t) A(T) in a well-conditioned chart.
        expected = np.array([self.at(x)[0] for x in t])
        err = float(np.max(abs((u+w*a_end)/expected-1)))
        if err>1e-7 or min(np.min(u),np.min(w))<-1e-10 or max(u+w)>1+1e-9:
            raise ArithmeticError('Two-endpoint interpolation check failed')
        return err


def minimum_velocity_variance(a, b, kappa, D):
    """Minimum per-component Var(v0) if both endpoint variances must be <=D."""
    if b==0.:
        return kappa**2/(4*D)
    s = kappa*b/2
    if s > D:
        return math.inf
    if (s/D)**2 <= 1-1/a**2:
        return (kappa*a)**2/(4*D)
    # Both endpoints saturate. Rationalize 1-sqrt(1-(s/D)^2).
    denom = 1+math.sqrt(max(0.,1-(s/D)**2))
    return D*((a-1)**2+2*a*(s/D)**2/denom)/b**2


def budget_horizon(mode, kappa, width, fraction, unconstrained):
    D = width**2/9
    def energy(t):
        a,b = mode.at(t)
        return 1.5*minimum_velocity_variance(a,b,kappa,D)
    if energy(0.) > fraction:
        return dict(status='exceeds budget at bounce',time=0.)
    end = np.nextafter(unconstrained,0.)
    # Compute exact endpoint on the RMS frontier to avoid ulp rounding of kappa*B/(2D).
    a,_ = mode.at(unconstrained)
    b = 2*D/kappa
    frontier = 1.5*D*(1+a*a)/b**2
    if frontier <= fraction:
        return dict(status='quantum RMS limit',time=unconstrained)
    def residual(t):
        if t >= end:
            return frontier-fraction
        return energy(t)-fraction
    return dict(status='kinetic budget limit',time=brentq(residual,0,end,xtol=1e-13))


def check_optimization():
    # Directly minimize over covariance matrices at given endpoint variances.
    # Fixing x sets optimal y=min(D,A²x+s²/x) by differentiating w.r.t. y.
    err = 0.
    for a in (1.01,1.5,8.,64.):
        for ratio in (.03,.4,.88,.999,.99999):
            D,b = 1.,.7
            s=ratio*D
            kap=2*s/b
            def objective(x):
                y=min(D,a*a*x+s*s/x)
                return (y+a*a*x-2*a*math.sqrt(max(0.,x*y-s*s)))/b**2
            ans=minimize_scalar(objective,bounds=(s*s/D,D),method='bounded',
                                options={'xatol':1e-14})
            benchmark=min(ans.fun,objective(D),objective(s*s/D))
            exact=minimum_velocity_variance(a,b,kap,D)
            err=max(err,abs(benchmark/exact-1))
    if err>1e-8:
        raise ArithmeticError('Endpoint covariance optimization disagrees')
    return err


def branch_result(branch):
    mu2 = branch['negative_mass_squared_over_Htop_squared']
    growth = branch['limiting_linear_growth_rate']
    h=(mu2-growth**2)/(3*growth)
    modes={n:Transfer(h,mu2,n,20.) for n in range(2,9)}
    outcomes={}
    for key,case in branch['detunings'].items():
        time=case['primary']['x_match']
        kappa=2*case['archive_displacement_scale']**2*h**3
        def minimum_rms(t, band):
            return math.sqrt(kappa/2*sum((n+1)**2*modes[n].at(t)[1] for n in band))
        entry={'time_required':time,'kappa':kappa,'n2_B':modes[2].at(time)[1],
               'band_results':{},'budget_tests':{}}
        for tag,band in [('n2',[2]), ('n2_to_8',list(modes))]:
            limits={}
            for d in (.0001,.001,.01):
                horizon=brentq(lambda t:minimum_rms(t,band)-d,0,time,xtol=2e-15)
                limits[str(d)]=horizon
                if tag=='n2':
                    entry['budget_tests'][str(d)]={str(f):budget_horizon(modes[2],kappa,d,f,horizon)
                                                    for f in (.01,.0001)}
            entry['band_results'][tag]=dict(minimax_rms_at_required_time=minimum_rms(time,band),
                                             rms_horizons=limits)
        entry['interpolation_check_max_relative_error']=max(m.endpoint_basis_check(time) for m in modes.values())
        outcomes[key]=entry
    return outcomes


def run(prior, output):
    previous=json.loads((prior/'hilltop/HILLTOP.json').read_text())['branches']
    result=dict(direct_covariance_minimum_max_relative_error=check_optimization(),
                branches={key:branch_result(data) for key,data in previous.items()})
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.prior,args.output)
