"""Check the norm of a proposed bulk photon and the fixed-sector energy of an added excitation."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp

HERE = Path(__file__).resolve().parent
SHELLS = np.array([2., 3., 4., 6., 8., 12., 16., 32., 64., 128., 256.])


def module_at(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def profile_integral(phi):
    # Integral from the shell to the throat, for G=(1-phi^2)^2, per copy.
    return (1+phi)**2*(2-phi)/3


def sech2(y):
    z = math.exp(-2*abs(y))
    return 4*z/(1+z)**2


def primitive(y):
    logcosh = abs(y)+math.log1p(math.exp(-2*abs(y)))-math.log(2.)
    return y+(2/3)*logcosh+math.tanh(y)**2/6


def static_patch(prior):
    source = prior/'origin'
    background = module_at(source/'hdblast_background.py','registered_static_chart')
    ledger = json.loads((source/'REGISTERED_SHELL_FLOAT_SOLUTION.json').read_bytes())
    end, pole = ledger['v_b'],ledger['phi_h']
    # Integrate from the regular cone jets. The fourth component is the warp-squared integral.
    eps = .001
    def rhs(v,state):
        rho,p,s = state[:3]
        rp = math.sqrt(1+rho*rho*(s*s/12-background.U(p)/6))
        return [rp,s,background.U1(p)-4*rp*s/rho,rho*rho]
    def solve(start,step):
        initial = [*background.series(pole,start),start**3/3-background.U(pole)*start**5/90]
        sol = solve_ivp(rhs,(start,end),initial,method='DOP853',rtol=2e-13,atol=1e-16,max_step=step)
        if not sol.success:
            raise ArithmeticError(sol.message)
        rho,p,s,area=sol.y[:,-1]
        return dict(rho_b=float(rho),phi_b=float(p),gravitational_M4_squared=float(2*area/rho**2))
    values=solve(eps,.012)
    refined=solve(eps/2,.006)
    error=max(abs(values[k]-refined[k]) for k in values)
    # Pole field datum is stored as a float, with a tiny distance from its AdS fixed point.
    if (abs(refined['gravitational_M4_squared']-2.0677328)>3e-6 or
        abs(refined['rho_b']/ledger['rho_b']-1)>1e-6 or error>1e-4):
        raise ArithmeticError('Static slice calibration changed')
    return dict(geometry='single registered detuned de Sitter exterior, in L0 units',
                cone_to_shell_proper_length=end,two_copy_constant_gauge_norm_times_g5_squared=2*end,
                H_shell=1/ledger['rho_b'],vector_continuum_threshold_over_H=.5,
                tensor_continuum_threshold_over_H=1.5,
                static_refined=refined,static_refinement_max_absolute=error)


def kinetic_data(fitted, theta):
    phi,f = fitted.spline(theta)
    dphi,d_f = fitted.spline(theta,1)
    integral=profile_integral(phi)
    return phi,f,integral,.5*dphi*(1-phi*phi)/integral


def location(fitted, state):
    theta=state['theta']
    phi,f,integral,slope=kinetic_data(fitted, theta)
    y=-np.arctanh(phi)
    numerical=quad(lambda u:sech2(y+u)**2,0,np.inf,epsabs=2e-13)[0]
    gravity=2*quad(lambda u:math.exp(-2*(primitive(y+u)-primitive(y))/3),0,np.inf,
                   epsabs=2e-12)[0]
    if max(abs(numerical-integral),abs(gravity-f))>1e-9:
        raise ArithmeticError('Localized or gravitational norm quadrature failed')
    return dict(theta=theta,phi_b=float(phi),M4_squared_rolling_leaf=float(f),
                half_norm_times_g5_squared=float(integral),
                d_log_I_d_theta=float(slope),d_log_I_d_x=float(slope*state['v']),
                norm_integral_error=float(abs(numerical-integral)),
                M4_integral_error=float(abs(gravity-f)))


def excess(branch, archive, step):
    """Energy difference of a single initial excitation, measured with the original Maxwell stress."""
    start=archive['initial']
    crossing=archive['first_minimum_crossing']
    fit=branch.fitted
    initial=np.array([start['theta'],0.,start['a'],0.,0.])
    def derivative(t,state):
        theta,v,a,h,eta=state
        potential,force=branch.potential(theta-branch.top)
        return np.array([v,force-3*h*v,a*h,1/a**2-v*v/2,1/a])
    timed=solve_ivp(derivative,(0.,crossing['x']),initial,method='DOP853',
                    rtol=2e-13,atol=2e-15,max_step=step,dense_output=True)
    if not timed.success:
        raise ArithmeticError(timed.message)
    end=timed.y[:,-1]
    trajectory_error=max(abs(end[:4]-[crossing['theta'],crossing['v'],crossing['a'],crossing['h']]))
    initial_norm=kinetic_data(fit,initial[0])[2]
    def coefficients(x):
        theta,v,a,h,eta=timed.sol(x)
        _,_,norm,slope=kinetic_data(fit,theta)
        return norm/initial_norm,slope*v,a,eta
    def mixing(t,columns):
        alpha,beta=columns.reshape(2,-1)
        _,rate,_,eta=coefficients(t)
        factor=rate*np.exp(2j*SHELLS*eta)
        return np.r_[factor*beta.conj(),factor*alpha.conj()]
    control=dict(method='DOP853',rtol=2e-12,atol=1e-14,max_step=step)
    zero=np.zeros_like(SHELLS)
    unit=np.ones_like(SHELLS)
    sol=solve_ivp(mixing,(0.,crossing['x']),np.r_[unit,zero].astype(complex),**control)
    def canonical(t,columns):
        q,p=columns.reshape(2,-1)
        i2,_,a,_=coefficients(t)
        return np.r_[p/(a*i2),-i2*SHELLS**2*q/a]
    # Arbitrary constant I(0) is absorbed into the charge convention once, at launch.
    canon=solve_ivp(canonical,(0.,crossing['x']),np.r_[1/np.sqrt(2*SHELLS),-1j*np.sqrt(SHELLS/2)],**control)
    if not sol.success or not canon.success:raise ArithmeticError('Gauge ODE failure')
    alpha,beta=sol.y[:,-1].reshape(2,-1)
    q,p=canon.y[:,-1].reshape(2,-1)
    i2,final_rate,a1,_=coefficients(crossing['x'])
    measured=(abs(p)**2/i2+i2*SHELLS**2*abs(q)**2)/SHELLS
    energy=1+2*abs(beta)**2
    variation=abs(math.log(i2)/2)
    # No interior turning point of I on the first monotone leg.
    speeds=timed.sol(np.linspace(0.,crossing['x'],513))[1]
    if min(branch.fitted.sign*speeds)<-1e-10:raise ArithmeticError('Monotone bound used on a turning trajectory')
    bound=math.cosh(2*variation)
    if (max(energy)>bound+1e-10 or trajectory_error>1e-8 or
        max(abs(measured-energy))>2e-9 or max(abs(abs(alpha)**2-abs(beta)**2-1))>2e-10):
        raise ArithmeticError('Canonical energy, phase norm, monotone variation, or background control failed')
    # One integration by parts in conformal time yields an instantaneous endpoint term.
    endpoint_coefficient=(a1*final_rate)**2/2
    measured_coefficient=float(2*SHELLS[-1]**2*abs(beta[-1])**2)
    if abs(measured_coefficient/endpoint_coefficient-1)>.01:
        raise ArithmeticError('Large-L finite-time endpoint coefficient changed')
    return dict(crossing_x=crossing['x'],first_crossing_eta=end[4],log_I_change=math.log(i2)/2,
                added_energy_over_adiabatic_upper_bound=bound,
                moving_endpoint_high_L_coefficient=endpoint_coefficient,
                largest_shell_measured_coefficient=measured_coefficient,
                round_redshift_E_adiabatic_over_E_bounce=initial[2]/a1,
                canonical_energy_comparison_max_error=float(max(abs(measured-energy))),
                Wronskian_max_error=float(max(abs(abs(alpha)**2-abs(beta)**2-1))),
                archived_background_max_error=float(trajectory_error),
                added_energy=[dict(L=int(l),energy_over_adiabatic=float(e),
                                   energy_over_bounce=float(e*initial[2]/a1)) for l,e in zip(SHELLS,energy)])


def calculate(prior):
    loader=module_at(prior/'homogeneous/exchange.py','rolling_parent')
    background=loader.branches(prior/'prior/gravity')
    archived=json.loads((prior/'homogeneous/results.json').read_bytes())['branches']
    data={}
    for key,branch in background.items():
        archive=archived[key]['runs']['inverse_f_seed_0']
        base=excess(branch,archive,.004)
        fine=excess(branch,archive,.002)
        change=max(abs(x['energy_over_adiabatic']-y['energy_over_adiabatic']) for
                   x,y in zip(base['added_energy'],fine['added_energy']))
        if change>2e-9:raise ArithmeticError('Mode refinement changed the energy')
        data[key]=dict(bounce=location(branch.fitted,archive['initial']),
                       first_minimum=location(branch.fitted,archive['first_minimum_crossing']),
                       response=fine,added_energy_refinement_max_shift=change)
    k=5/9
    exponent=(1+4/k)/2
    return dict(static_cone=static_patch(prior),infinite_tuned_leaf=dict(
                throat_curvature_k=k,minimal_photon_normalizable=False,
                localized_G='(1-phi^2)^2',vacuum_zero_mode_constant=True,
                canonical_zero_mode_asymptotic_power=exponent,
                inverse_square_potential_coefficient=exponent*(exponent+1),
                continuum_lower_edge_mass_squared=0.),
                branches=data)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=HERE.parent/'prior')
    parser.add_argument('--output',type=Path,required=True)
    opts=parser.parse_args()
    opts.output.write_text(json.dumps(calculate(opts.prior),indent=2,allow_nan=False)+'\n')
