"""Homogeneous canonical wall scalar, powered by incident outgoing Vaidya radiation."""
import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import expit


DURATION = .01
LOW = .004
HIGH = .012
INNER = 1.6


def bump(x):
    if 0 < x < 1:
        return math.exp(4 - 1 / (x * (1 - x)))
    return 0.


BUMP_AREA = quad(bump, 0., 1., epsabs=1e-14, epsrel=1e-13)[0]


def potential(psi, height):
    if psi <= LOW:
        return 0., 0.
    if psi >= HIGH:
        return height, 0.
    x = (psi - LOW) / (HIGH - LOW)
    slope = -1 / x + 1 / (1 - x)
    rise = expit(slope)
    derivative = expit(-abs(slope)) * expit(abs(slope)) * (1 / x ** 2 + 1 / (1 - x) ** 2)
    return height * rise, height * derivative / (HIGH - LOW)


def run(kick=.485, height=.111, tolerance=2e-11):
    def sources(t, state):
        a, h, mass, psi, pi = state[:5]
        pot, gradient = potential(psi, height)
        e, p = .5*pi*pi + pot, .5*pi*pi - pot
        j = kick * bump(t / DURATION) / (DURATION * BUMP_AREA)
        Q = j * pi
        B = 1 + e + h
        return e, p, gradient, Q, B

    def rhs(t, state):
        a, h, mass, psi, pi = state[:5]
        e, p, gradient, Q, B = sources(t, state)
        j = kick * bump(t / DURATION) / (DURATION * BUMP_AREA)
        accel = (1 + e) * (1 - 2*e - 3*p) - 1 - mass/a**4 - Q
        return [a*h, accel-h*h, -2*a**4*B*Q, pi, -3*h*pi - gradient + j,
                1/(a*B), Q, 3*a**3*h*p, a**3*Q]

    def reach_bottom(t,state):
        return state[0] - INNER
    reach_bottom.direction = -1
    reach_bottom.terminal = True

    def return_up(t,state):
        return state[0] - 2.05
    return_up.direction = 1
    return_up.terminal = True

    def bounce(t,state):
        return state[1]
    bounce.direction = 1

    def reach_top(t,state):
        return state[3] - HIGH
    reach_top.direction = 1

    initial = [2.,0.,4.,0.,0.,0.,0.,0.,0.]
    def integrate(interval, state, step):
        s = solve_ivp(rhs, interval, state, rtol=tolerance,atol=tolerance*.002,method='DOP853',
                      max_step=step,events=(reach_bottom,return_up,bounce,reach_top),dense_output=True)
        if not s.success:
            raise ArithmeticError(s.message)
        return s

    pulse = integrate((0.,DURATION),initial,DURATION/30)
    wall = integrate((DURATION,10.),pulse.y[:,-1],.005)
    end = wall.t[-1]
    turned = len(wall.t_events[2])>0
    plateau = len(wall.t_events[3])>0 and height>0
    def snapshot(t,state):
        a,h,mu,psi,pi,u,gain,pdv,supply=map(float,state)
        e,p,_,Q,B=sources(t,state)
        return dict(t=float(t),a=a,H=h,mu=mu,psi=psi,pi=pi,e=e,p=p,B=B,u=u,Q=Q)
    result=dict(kick=kick,height=height,turn=turned,plateau=plateau,terminal=snapshot(end,wall.y[:,-1]),
                pulse_end=snapshot(DURATION,pulse.y[:,-1]))
    if turned:
        result['bounce']=snapshot(wall.t_events[2][0],wall.y_events[2][0])
    if plateau:
        result['plateau_entry']=snapshot(wall.t_events[3][0],wall.y_events[3][0])
    times = np.r_[np.linspace(0.,DURATION,101),wall.t[1:]]
    states = np.vstack([pulse.sol(times[:101]).T,wall.y[:,1:].T])
    coeff = np.array([sources(t,x) for t,x in zip(times,states)]).T
    e,pres,_,flux,B = coeff
    a,H,mu,_,pi,u,_,work,heat = states.T
    friedmann = H**2 - (e*(2+e)+mu/a**4-1/a**2)
    lapse = 1+a*a-mu/a**2
    udot = 1/(a*B)
    f_tube = 1 + INNER**2 - mu/INNER**2
    mu_u = -2*a**4*B*flux/udot
    # kappa_5^2 T_uu on the stationary emitter. Both sides share this profile.
    stress = -1.5*mu_u / INNER**3
    emitter_energy = INNER**2 * (math.sqrt(f_tube[-1]) - math.sqrt(f_tube[0]))
    def emitter_rate(t):
        x = pulse.sol(t)
        Q,B = sources(t,x)[3:]
        return x[0]**4*B*Q / math.sqrt(1+INNER**2-x[2]/INNER**2)
    integrated_emitter = quad(emitter_rate,0.,DURATION,epsabs=1e-12,epsrel=1e-12)[0]
    errors = dict(friedmann=float(max(abs(friedmann))),
                  unsquared_spatial_israel=float(max(abs(np.sqrt(lapse+(a*H)**2)/a-1-e))),
                  retarded_clock=float(max(abs(lapse*udot**2+2*a*H*udot-1))),
                  first_law=float(max(abs(heat-a**3*e-work))),
                  emitter_integral=float(abs(integrated_emitter-emitter_energy)))
    result['energy'] = dict(two_side_emitter_energy=emitter_energy,
                            two_side_wall_deposition=float(heat[-1]),
                            pulse_end_internal_energy=float(a[100]**3*e[100]),
                            pulse_end_pressure_work=float(work[100]),
                            terminal_internal_energy=float(a[-1]**3*e[-1]),
                            terminal_pressure_work=float(work[-1]),
                            mu_loss_each_side=float(4-mu[-1]))
    result['domain'] = dict(min_mu=float(min(mu)),min_B=float(min(B)),
                            min_f_at_source=float(min(f_tube)),min_null_stress_at_source=float(min(stress)),
                            max_sampled_null_stress_at_source=float(max(stress)),
                            min_source_power=float(min(flux)),min_e_plus_p=float(min(pi**2)),
                            retarded_width=float(pulse.y[5,-1]))
    if plateau:
        enter=result['plateau_entry']
        C=.5*(enter['pi']*enter['a']**3)**2
        mass=enter['mu']
        def F(r):
            return height*(2+height)*r*r+2*(1+height)*C/r**4+C*C/r**10+mass/r**2-1
        def derivative(r):
            return 2*height*(2+height)*r-8*(1+height)*C/r**5-10*C*C/r**11-2*mass/r**3
        minimum=brentq(derivative,.01,50)
        first_integral=dict(C=C,minimum_at=minimum,minimum=F(minimum),
                            entry_on_outer_side=bool(enter['a']>minimum))
        if F(minimum)<0:
            first_integral['inner_root']=brentq(F,.01,minimum)
            first_integral['outer_root']=brentq(F,minimum,50.)
        after=times>=enter['t']
        errors['plateau_radial_integral']=float(max(abs((a[after]*H[after])**2-F(a[after]))))
        errors['plateau_scalar_momentum']=float(max(abs(pi[after]*a[after]**3-math.sqrt(2*C))))
        if turned:
            b=result['bounce']
            errors['outer_turning_radius']=abs(first_integral['outer_root']-b['a'])
            first_integral['acceleration_at_bounce']=derivative(b['a'])/2
        result['plateau_first_integral']=first_integral
    result['errors']=errors
    if max(errors.values())>1e-8 or min(mu)<=0 or min(B)<=0 or min(f_tube)<=0 or min(flux)<-1e-12:
        raise ArithmeticError((kick,height,errors,result['domain']))
    return result



def lower_plateau_threshold(tolerance):
    """Shoot to the plateau, then bracket a double root of its radial first integral."""
    lower, upper = .099, .104

    def first_integral(height):
        return run(kick=.485, height=height, tolerance=tolerance)['plateau_first_integral']

    height = brentq(lambda u: first_integral(u)['minimum'], lower, upper, xtol=1e-12)
    solution = first_integral(height)
    radius = solution['minimum_at']
    if not solution['entry_on_outer_side'] or radius <= INNER:
        raise ArithmeticError('Inaccessible limiting radius')
    return dict(bracket=[lower, upper], kick=.485, critical_height=height,
                limiting_radius=radius, errors=dict(radial_polynomial=abs(solution['minimum'])))


def calculate():
    cases = {
        'primary': (.485,.111),
        'same_potential_no_pulse': (0.,.111),
        'same_pulse_zero_potential': (.485,0.),
        'insufficient_plateau_height': (.485,.099),
        'inaccessible_plateau': (.485,.121),
        'weaker_pulse': (.480,.111),
        'stronger_pulse': (.490,.111),
        'lower_plateau_with_bounce': (.485,.104),
        'higher_plateau_with_bounce': (.485,.118),
    }
    rows = {name:run(*parameters) for name,parameters in cases.items()}
    expected_bounces = {'primary','weaker_pulse','stronger_pulse',
                        'lower_plateau_with_bounce','higher_plateau_with_bounce'}
    if {name for name,row in rows.items() if row['turn']} != expected_bounces:
        raise ArithmeticError('Unexpected turn classification')
    if not rows['insufficient_plateau_height']['plateau'] or rows['inaccessible_plateau']['plateau']:
        raise ArithmeticError('Unexpected plateau classification')
    refined = run(.485,.111,2e-12)
    diff = abs(rows['primary']['bounce']['a']-refined['bounce']['a'])
    if diff > 1e-8:
        raise ArithmeticError('Unresolved turning radius')
    threshold = lower_plateau_threshold(2e-11)
    threshold_refined = lower_plateau_threshold(2e-12)
    for key in ('critical_height', 'limiting_radius'):
        difference = abs(threshold[key]-threshold_refined[key])
        if difference > 5e-9:
            raise ArithmeticError('Unresolved critical plateau')
        threshold['errors']['refined_'+key] = difference
    return dict(lower_plateau_threshold=threshold, duration=DURATION, potential_start=LOW, potential_end=HIGH,
                source_radius=INNER, cases=rows, errors=dict(refined_turning_radius=diff))


if __name__=='__main__':
    import sys
    for args in sys.argv[1:]:
        k,h=map(float,args.split(','))
        print(run(k,h))
