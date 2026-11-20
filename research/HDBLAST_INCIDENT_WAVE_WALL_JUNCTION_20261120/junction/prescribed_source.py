"""Launch a freely prescribed positive retarded mass pulse into the same canonical surface scalar."""
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def potential_function():
    script = Path(__file__).with_name('scalar_capture.py')
    spec = importlib.util.spec_from_file_location('canonical_capture_potential', script)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    return model.potential


def run(fraction, width, height, potential, *, follow_exit=False):
    initial_mass = 4.
    emitted_radius = 1.6

    def mass(s):
        smooth = s**4*(35-84*s+70*s*s-20*s**3)
        derivative = 140*s**3*(1-s)**3
        return initial_mass*(1-fraction*smooth), -initial_mass*fraction*derivative

    def receive(s, state):
        tau, a, h, kinetic, psi, deposit, emitter_work = state
        mu, dm = mass(s)
        b = 1+kinetic+h
        clock = width*a*b
        source = -dm/(2*width*a**5*b*b)
        accel = (1+kinetic)*(1-5*kinetic)-1-mu/a**4-source
        return [clock, clock*a*h, clock*(accel-h*h), clock*(source-6*h*kinetic),
                clock*math.sqrt(2*max(0., kinetic)), clock*a**3*source,
                -dm/(2*math.sqrt(1+emitted_radius**2-mu/emitted_radius**2))]

    options = dict(method='DOP853',rtol=2e-12,atol=1e-14)
    pulse = solve_ivp(receive,(0.,1.),[0.,2.,0.,0.,0.,0.,0.],max_step=.02,**options)
    if not pulse.success:
        raise ArithmeticError(pulse.message)
    tau, a0, h0, k0, psi0, deposit, emitter_work = pulse.y[:,-1]
    mu = mass(1.)[0]
    if not (psi0 < .004 and min(pulse.y[3]) >= -1e-14):
        raise ArithmeticError('Pulse leaves the zero-potential region or has negative kinetic energy')

    def rhs(t, state):
        a,h,psi,pi,u = state
        v,grad = potential(psi,height)
        e,p = .5*pi*pi+v, .5*pi*pi-v
        return [a*h, (1+e)*(1-2*e-3*p)-1-mu/a**4-h*h, pi, -3*h*pi-grad,1/(a*(1+e+h))]

    def turn(t, state):
        return state[1]
    turn.direction=1

    def recovered(t, state):
        return state[0]-2.05
    recovered.direction=1
    recovered.terminal=True

    def inner(t,state):
        return state[0]-emitted_radius
    inner.direction=-1
    inner.terminal=True

    def plateau(t,state):
        return state[2]-.012
    plateau.direction=1

    after=solve_ivp(rhs,(tau,10.),[a0,h0,psi0,math.sqrt(2*k0),width],max_step=.005,
                    events=(turn,recovered,inner,plateau),**options)
    if not after.success:
        raise ArithmeticError(after.message)
    mu_grid=mass(pulse.t)[0]
    _,a,h,e,*_=pulse.y
    err=np.max(abs(h*h-e*(2+e)-mu_grid/a**4+1/a**2))
    a,h,psi,pi,u=after.y
    e=.5*pi**2+np.array([potential(x,height)[0] for x in psi])
    err=max(err,np.max(abs(h*h-e*(2+e)-mu/a**4+1/a**2)))
    if err>1e-9 or len(after.t_events[1])+len(after.t_events[2])!=1:
        raise ArithmeticError('Trajectory or gravitational constraint unresolved')
    turning=None
    if len(after.t_events[0]):
        tb=after.t_events[0][0]
        ab,_,psi_b,pi_b,_=after.y_events[0][0]
        v=potential(psi_b,height)[0]
        eb,pb=.5*pi_b**2+v,.5*pi_b**2-v
        turning=dict(time=tb,radius=ab,energy=eb,pressure=pb,
                     acceleration_over_radius=(1+eb)*(1-2*eb-3*pb)-1-mu/ab**4)
    exit_test = None
    if follow_exit:
        spec = importlib.util.spec_from_file_location('finite_wall_plateau', Path(__file__).with_name('radiation_exit.py'))
        exit_model = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(exit_model)
        join = dict(t=after.t[-1],a=a[-1],H=h[-1],mu=mu,psi=psi[-1],pi=pi[-1],e=e[-1],u=u[-1])
        capture = dict(terminal=join,turn=turning is not None,pulse_end=dict(t=tau,mu=mu,e=k0))
        late = exit_model.run_after_capture(capture, None, exit_model.DAMPING)
        if not (late['e_folds_radiation']>0 and 'maximum' in late['events']):
            raise ArithmeticError('The specified independent beam did not produce the predicted radiation interval')
        exit_test = dict(events=late['events'], windows=late['expanding_radiation_intervals'],
                         e_folds=late['e_folds_radiation'], outcome=late['outcome'], checks=late['checks'])
    return dict(mass_loss_fraction=fraction, retarded_width=width, plateau_height=height,
                final_bulk_mass=mu, pulse_end_proper_time=tau, pulse_end_field=psi0,
                pulse_end_kinetic_energy_density=k0, deposited_proper_energy=deposit,
                proper_emitter_energy=emitter_work, reached_positive_plateau=bool(height and len(after.t_events[3])),
                bounce=turning, end_radius=a[-1], end_H=h[-1], max_friedmann_residual=float(err), radiation_exit=exit_test)


def calculate():
    potential=potential_function()
    rows={name:run(f,w,u,potential,follow_exit=name=='baseline') for name,f,w,u in [
        ('baseline',.94,.005,.111),('less_mass',.85,.005,.111),('more_mass',.97,.005,.111),
        ('shorter',.94,.004,.111),('no_potential',.94,.005,0.),('no_pulse',0.,.005,.111)]}
    for name in ('baseline','more_mass','shorter'):
        if not rows[name]['bounce'] or rows[name]['end_H']<=0:
            raise ArithmeticError('Missing predicted bounce in '+name)
    for name in ('no_potential','no_pulse'):
        if rows[name]['bounce'] or rows[name]['end_H']>=0:
            raise ArithmeticError('Control bounced')
    return rows


if __name__=='__main__':
    print(json.dumps(calculate(),indent=2))
