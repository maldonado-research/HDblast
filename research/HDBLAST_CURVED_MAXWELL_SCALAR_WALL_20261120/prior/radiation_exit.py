"""Close the finite scalar plateau and test energy transfer to an ideal w=1/3 fluid."""
import math
from pathlib import Path
import importlib.util

from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import expit


def load_scalar():
    spec = importlib.util.spec_from_file_location(
        "hdblast_capture", Path(__file__).with_name("scalar_capture.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


scalar = load_scalar()
EXIT = .36
WIDTH = .02
HEIGHT = .111
DAMPING = 4.


def smooth_step(x):
    if x <= 0.:
        return 0., 0.
    if x >= 1.:
        return 1., 0.
    slope = 1/(1-x) - 1/x
    return expit(slope), expit(-abs(slope))*expit(abs(slope))*(1/x**2 + 1/(1-x)**2)


def properties(psi, damping):
    step, derivative = smooth_step((psi-EXIT)/WIDTH)
    v, dv = scalar.potential(psi, HEIGHT)
    return v-HEIGHT*step, dv-HEIGHT*derivative/WIDTH, damping*step


def plateau_reach(join):
    """Limiting field of the outgoing solution extended with an infinite plateau."""
    a, mu = join["a"], join["mu"]
    momentum = a**3 * join["pi"]
    slope = 2 - 4*mu/a**2 - 6*(1+HEIGHT)*momentum**2/a**4 - 3*momentum**4/a**10
    if slope<=0. or join['H']<=0. or momentum<=0.:
        raise ArithmeticError("Plateau monotonicity certificate does not apply")
    def integral(x):
        e = HEIGHT + .5*momentum**2*x**6
        return momentum*x**2/math.sqrt(2*e+e*e-x*x+mu*x**4)
    advance = quad(integral, 0., 1/a, epsabs=2e-13, epsrel=1e-12)[0]
    return dict(momentum=momentum, limiting_psi=float(join["psi"]+advance),
                minimum_scaled_Friedmann_slope=float(slope))


def run(kick, damping, tol=2e-11, step=.015, samples=False):
    capture = scalar.run(kick=kick, height=HEIGHT, tolerance=tol)
    return run_after_capture(capture, kick, damping, tol=tol, step=step, samples=samples)


def run_after_capture(capture, kick, damping, tol=2e-11, step=.015, samples=False):
    join = capture['terminal']
    result = dict(kick=kick, damping=damping, prior_bounce=capture['turn'],
                  pulse_end=capture['pulse_end'], join=join)
    if 'bounce' in capture:
        result['bounce']=capture['bounce']
    if join['psi'] >= EXIT or (not capture['turn'] and kick != 0):
        raise ArithmeticError("Upstream benchmark does not cover the specified handover")
    if not capture['turn']:
        result.update(outcome="contracted to the inner source", e_folds_radiation=0.)
        return result
    if not math.isclose(join['a'], 2.05) or join['H']<=0 or join['pi']<=0:
        raise ArithmeticError("Handover requires outward motion on the positive scalar branch")
    result['constant_plateau']=plateau_reach(join)
    mu=join['mu']
    def matter(psi, pi, rad):
        v, grad, gamma = properties(psi, damping)
        return pi*pi/2+v+rad, pi*pi/2-v+rad/3, grad, gamma
    def rhs(t, state):
        a,h,psi,pi,rad,work,transfer,clock=state
        e,p,grad,gamma=matter(psi,pi,rad)
        accel=(1+e)*(1-2*e-3*p)-1-mu/a**4
        power=gamma*pi**2
        return [a*h, accel-h*h, pi, -(3*h+gamma)*pi-grad, power-4*h*rad,
                3*a**3*h*p, a**4*power, 1/(a*(1+e+h))]

    def mark(index, threshold, direction=0, terminal=False):
        def event(t, state):
            return state[index]-threshold
        event.direction, event.terminal=direction, terminal
        return event
    def material_cross(t, state):
        a,h,psi,pi,rad=state[:5]
        return rad-.9*matter(psi,pi,rad)[0]
    events={'exit':mark(2,EXIT,1),'downhill_end':mark(2,EXIT+WIDTH,1),
            'rad90':material_cross,'maximum':mark(1,0.,-1),
            'contracted':mark(0,2.05,-1,True),'large':mark(0,50.,1,True)}
    series=solve_ivp(rhs,(join['t'],35.),[join['a'],join['H'],join['psi'],join['pi'],0.,0.,0.,join['u']],
                     events=list(events.values()),dense_output=True,rtol=tol,atol=.002*tol,
                     max_step=step,method='DOP853')
    if not series.success:
        raise ArithmeticError(series.message)
    def snapshot(t, state):
        a,h,psi,pi,rad=state[:5]
        e,p,_,gamma=matter(psi,pi,rad)
        positive_terms = 2*e + e**2 + mu/a**4
        return dict(t=float(t),a=float(a),H=float(h),psi=float(psi),pi=float(pi),
                    e=float(e),radiation=float(rad),radiation_fraction=float(rad/e),
                    p=float(p),Gamma=float(gamma),mu=float(mu),u=float(state[7]),
                    linear_radiation_fraction_positive_friedmann=float(2*rad/positive_terms),
                    curvature_to_positive_friedmann=float(1/(a*a*positive_terms)))
    times=dict(zip(events,series.t_events))
    result['events']={name:[snapshot(t,series.sol(t)) for t in values]
                      for name,values in times.items() if len(values)}
    end=series.t[-1]
    cuts=sorted(set([float(join['t']),float(end),*times['rad90'],*times['maximum']]))
    windows=[]
    for lower,upper in zip(cuts,cuts[1:]):
        middle=series.sol((lower+upper)/2)
        if middle[1]>0 and material_cross(0,middle)>0:
            start=snapshot(lower,series.sol(lower))
            finish=snapshot(upper,series.sol(upper))
            windows.append(dict(start=start,end=finish,e_folds=float(math.log(finish['a']/start['a']))))
    result['expanding_radiation_intervals']=windows
    result['e_folds_radiation']=sum(w['e_folds'] for w in windows)
    result['outcome']=('expanded then recontracted to a=2.05' if len(times['contracted']) else
                       'expanded to a=50' if len(times['large']) else 'reached time limit')
    result['terminal']=snapshot(end,series.sol(end))
    if samples:
        series_times=[join['t']+(end-join['t'])*index/80 for index in range(81)]
        result['samples']=[snapshot(t,series.sol(t)) for t in series_times]
    residuals={'friedmann_absolute':0.,'comoving_first_law_absolute':0.,
               'radiation_comoving_transfer_absolute':0.,'angular_jump_absolute':0.}
    min_radius=float('inf')
    min_rad=0.
    min_b=float('inf')
    for t,state in zip(series.t,series.y.T):
        a,h,psi,pi,rad,work,transfer,clock=state
        e,p,_,_=matter(psi,pi,rad)
        fwall=1+a*a-mu/a**2
        min_radius=min(min_radius,a)
        min_rad=min(min_rad,rad)
        min_b=min(min_b,1+e+h)
        values=[h*h-(2*e+e**2-1/a**2+mu/a**4),
                a**3*e-join['a']**3*join['e']+work,
                a**4*rad-transfer,
                math.sqrt(fwall+(a*h)**2)/a-(1+e)]
        for key,val in zip(residuals,values):
            residuals[key]=max(residuals[key],float(abs(val)))
    result['checks']=dict(max_errors=residuals,min_radiation=float(min_rad),
                          min_retarded_clock_denominator_over_a=float(min_b),
                          minimum_wall_f_after_handover=float(1+min_radius**2-mu/min_radius**2),
                          source_f_lower_bound_all_times=1+scalar.INNER**2-4/scalar.INNER**2)
    # Comoving volumes grow as a^3 and a^4. Use corresponding scales for the guards.
    largest_a=float(max(series.y[0]))
    if residuals['friedmann_absolute']>1e-7 or residuals['angular_jump_absolute']>1e-7:
        raise ArithmeticError('Late Friedmann or unsquared jump check failed')
    if (residuals['comoving_first_law_absolute']>1e-7*largest_a**3 or
            residuals['radiation_comoving_transfer_absolute']>1e-8*largest_a**4 or
            min_rad< -1e-11 or min_b<=0 or mu<0):
        raise ArithmeticError('Late energy, mass or retarded-time check failed')
    if len(times['downhill_end']):
        td=float(times['downhill_end'][0])
        ad,_,_,pid=series.sol(td)[:4]
        scalar_error=0.
        for t,state in zip(series.t,series.y.T):
            if t>=td:
                exact=pid*(ad/state[0])**3*math.exp(-damping*(t-td))
                scalar_error=max(scalar_error,float(abs(state[3]-exact)))
        if scalar_error>2e-8:
            raise ArithmeticError('Force-free scalar damping law failed')
        result['checks']['max_force_free_scalar_error']=scalar_error
    if len(times['maximum']) and len(windows):
        peak=result['events']['maximum'][0]
        R=peak['a']**4*peak['radiation']
        def rad_turn(a):
            return a**2 - (2*R+mu) - (R/a**2)**2
        ar=brentq(rad_turn,.5,2*peak['a'])
        result['checks']['pure_radiation_turn_at_same_comoving_radiation']=float(ar)
    return result


def calculate():
    baseline=run(.485,DAMPING,samples=True)
    results=dict(baseline=baseline,below_exit=run(.483,DAMPING),lower_kick=run(.484,DAMPING),higher_kick=run(.486,DAMPING),
                 zero_kick=run(0.,DAMPING),zero_damping=run(.485,0.),
                 reference=run(.485,DAMPING,5e-12,.0075))
    for key in ('baseline','lower_kick','higher_kick','reference'):
        if len(results[key]['expanding_radiation_intervals'])!=1 or 'contracted' not in results[key]['events']:
            raise ArithmeticError('Expected a material-radiation interval followed by recontraction')
    for key in ('zero_kick','zero_damping','below_exit'):
        if results[key]['e_folds_radiation']!=0.:
            raise ArithmeticError('A negative control unexpectedly entered the measured epoch')
    if (abs(baseline['e_folds_radiation']-results['reference']['e_folds_radiation'])>1e-7 or
            results['below_exit']['constant_plateau']['limiting_psi']>=EXIT):
        raise ArithmeticError('Precision or limiting-field control failed')
    def excess(kick):
        join=scalar.run(kick=kick,height=HEIGHT,tolerance=2e-11)['terminal']
        return plateau_reach(join)['limiting_psi']-EXIT
    bracket=[.4832,.4833]
    critical=brentq(excess,*bracket,xtol=2e-12)
    return dict(exit_psi=EXIT,downhill_width=WIDTH,plateau_height=HEIGHT,damping=DAMPING,
                plateau_exit_threshold=dict(kick=float(critical),bracket=bracket),runs=results)
