"""Timelike tidal tests for an exact null Maxwell uplift of the tuned wall."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq


K = 5/9
ENERGY = 1.4
MOMENTUM = .9
LAMBDA = 1.


def warp(y):
    v = math.exp(-2*y)
    t = (1-v)/(1+v)
    logcosh = y+math.log1p(v)-math.log(2)
    return -(y+2*logcosh/3+t*t/6)/3


def profile(y):
    v = math.exp(-2*y)
    t = (1-v)/(1+v)
    phi_offset = 2*v/(1+v)
    delta = phi_offset**2*(t+2)/9
    ay = -K+delta
    ayy = -(4*v/(1+v)**2)**2/3
    a = math.exp(warp(y))
    radial_excess = delta*(2*K-delta)/(K*K*ay*ay*a*a)
    return a, ay, ayy, phi_offset, radial_excess


def integrate(f, lo, hi):
    return quad(f,lo,hi,epsabs=2e-14,epsrel=1e-13,limit=120)[0]


def dt_dy(y):
    a = math.exp(warp(y))
    return a/math.sqrt(ENERGY**2-a*a)


def du_dy(y):
    a = math.exp(warp(y))
    return MOMENTUM/(a*math.sqrt(ENERGY**2-a*a))


def time_remaining(y):
    # Shift the variable and integrate to a finite cutoff. At 100 further
    # wall widths the discarded tail is <2e-24 times the leading remainder.
    a = math.exp(warp(y))
    scaled = lambda t: dt_dy(y+t)/a
    return a*integrate(scaled,0,100)


def sandwich(u):
    if not 1. < u < 3.:
        return 0.
    coordinate = u-2
    return math.exp(1-1/(1-coordinate*coordinate))


def state(y, power):
    a, ay, ayy, _, _ = profile(y)
    u = integrate(du_dy,0,y)
    wave = 1. if power == 0 else math.hypot(1.,u)**(-power)
    up = MOMENTUM/a**2
    radial_speed_squared = (ENERGY/a)**2-1.
    signal = LAMBDA*(wave*up)**2
    background_ricci = 4*(ayy+ay*ay)-3*ayy*(ENERGY/a)**2
    screen_background = ay*ay-ayy*radial_speed_squared
    radial_background = ay*ay+ayy
    # Three spacelike directions within the four-slice: two transverse screen
    # eigenvalues and one longitudinal. The fourth is the boosted radial normal.
    # The null wave adds transverse tides only on the axial geodesic.
    trace = 3*screen_background+radial_background+signal
    if abs(trace-(background_ricci+signal))>1e-12*(1+abs(trace)):
        raise ArithmeticError("Tidal trace disagrees with the Ricci contraction")
    return dict(y=y,warp_factor=a,u=u,remaining_time=time_remaining(y),
                signal_ricci=signal,background_ricci=background_ricci,
                screen_tidal=screen_background+signal/2,
                maxwell_screen=abs(wave*up/a),
                weighted_maxwell_screen_squared=(wave*up)**2)


def geodesic_check(y_final):
    initial=[0.,math.sqrt(ENERGY**2-1),0.,MOMENTUM,0.,ENERGY**2/(2*MOMENTUM)]
    def rhs(t,s):
        y,yd,u,ud,v,vd=s
        a,ay,_,_,_=profile(y)
        return [yd,-2*a*a*ay*ud*vd,ud,-2*ay*yd*ud,vd,-2*ay*yd*vd]
    def target(t,s):
        return s[0]-y_final
    target.terminal=True
    target.direction=1
    result=solve_ivp(rhs,(0.,time_remaining(0.)+.01),initial,
                     events=target,method='DOP853',rtol=3e-13,atol=3e-15,max_step=.015)
    if not result.success or len(result.t_events[0])!=1:
        raise ArithmeticError('Timelike geodesic missed the requested cutoff')
    t=float(result.t_events[0][0])
    y,yd,u,ud,v,vd=result.y_events[0][0]
    a=math.exp(warp(y))
    clock=integrate(dt_dy,0,y)
    expect_u=integrate(du_dy,0,y)
    momenta=[abs(a*a*ud/MOMENTUM-1),abs(2*MOMENTUM*a*a*vd/ENERGY**2-1)]
    errors=dict(proper_time_absolute=abs(t-clock),retarded_time_relative=abs(u/expect_u-1),
                radial_velocity_relative=abs(yd/math.sqrt((ENERGY/a)**2-1)-1),
                maximum_conserved_momentum_relative=max(momenta))
    if max(errors.values())>5e-10:raise ArithmeticError(errors)
    return errors


def logarithmic_slope(r1,r2,v1,v2):
    return math.log(v2/v1)/math.log(r2/r1)


def chart_test(y, power):
    # Hyperboloid chart (r,U,X1,X2,W), with V eliminated by the constraint.
    # At r=0 use the analytic continuation, not the singular Poincare chart.
    ubar, x1, x2, w = 1.7, .4, -.2, .12
    ell = 1/K
    if y is None:
        r, radial_excess = 0., 0.
    else:
        r, _, _, _, radial_excess = profile(y)
    vbar = (ell**2+x1**2+x2**2-ell*r*w)/(2*ubar)
    dr,du,dx,dx2,dw = np.eye(5)
    dv = (-ell*w*dr-2*vbar*du+2*x1*dx+2*x2*dx2-ell*r*dw)/(2*ubar)
    def product(a,b):
        return np.outer(a,b)
    baseline = (-ell/2*(product(dr,dw)+product(dw,dr))-product(du,dv)-product(dv,du)
                +product(dx,dx)+product(dx2,dx2)+radial_excess*product(dr,dr))
    if r==0 and power<3:
        # At p=2 the metric exists; the Maxwell coordinate two-form does not.
        finite_field=False
    else:
        finite_field=True
    u_numerator = r*du-ubar*dr
    den = r*r+ubar*ubar
    amplitude_squared = r**(2*power-4)/den**power
    metric = baseline-.5*LAMBDA*(x1*x1+x2*x2)*amplitude_squared*product(u_numerator,u_numerator)
    signature = np.linalg.eigvalsh(metric)
    if np.sum(signature<0)!=1 or np.sum(signature>0)!=4:
        raise ArithmeticError('Horizon chart metric became degenerate')
    det_ratio = float(np.linalg.det(metric)/np.linalg.det(baseline))
    if abs(det_ratio-1)>2e-12:raise ArithmeticError('Null wave changed metric determinant')
    row=dict(r=r,wave_metric_rr=float((metric-baseline)[0,0]),det_ratio=det_ratio)
    if finite_field:
        def wedge(a,b):
            return product(a,b)-product(b,a)
        form = r**(power-3)/den**(power/2)*(r*wedge(du,dx)-x1*wedge(du,dr)-ubar*wedge(dr,dx))
        row['field_r_x1']=float(form[0,2])
        inverse=np.linalg.inv(metric)
        invariant=float(np.sum(form*(inverse@form@inverse)))
        if abs(invariant)>5e-14*max(1.,float(np.max(abs(form)))**2):
            raise ArithmeticError('Ingoing-chart Maxwell invariant is nonzero')
    return row


def transverse_check(width_x,width_y):
    # Gaussian e=(f,0): the two Maxwell residuals independently measure
    # the two transverse derivatives. Fourier differentiation on a wide
    # rectangle makes the periodic boundary error exponentially small.
    n=192
    dx=16*width_x/n
    dy=16*width_y/n
    x=(np.arange(n)-n/2)*dx
    y=(np.arange(n)-n/2)*dy
    f=np.exp(-.5*((x[:,None]/width_x)**2+(y[None,:]/width_y)**2))
    hat=np.fft.fft2(f)
    kx=2*math.pi*np.fft.fftfreq(n,d=dx)
    ky=2*math.pi*np.fft.fftfreq(n,d=dy)
    divergence=np.fft.ifft2(1j*kx[:,None]*hat).real
    curl=np.fft.ifft2(-1j*ky[None,:]*hat).real
    square=f*f
    norm=float(square.sum()*dx*dy)
    rms=math.sqrt(float(((x[:,None]**2+y[None,:]**2)*square).sum()*dx*dy)/norm)
    residual=math.sqrt(float((divergence**2+curl**2).sum()*dx*dy)/norm)
    analytic=(width_x/width_y+width_y/width_x)/2
    if abs(rms*residual-analytic)>5e-12:raise ArithmeticError('Gaussian residual disagrees with Fourier check')
    return dict(width_x=width_x,width_y=width_y,norm=norm,rms_radius=rms,
                normalized_residual=residual,product=rms*residual,analytic_product=analytic)


def calculate():
    depths=(1,4,8,12,20,28,32)
    tails={}
    for power in (0.,1.,1.5,2.,2.5,3.):
        rows=[state(y,power) for y in depths]
        l,r=rows[-2:]
        slope=logarithmic_slope(l['remaining_time'],r['remaining_time'],l['signal_ricci'],r['signal_ricci'])
        maxwell=logarithmic_slope(l['remaining_time'],r['remaining_time'],l['maxwell_screen'],r['maxwell_screen'])
        if abs(slope-(2*power-4))>1e-5 or abs(maxwell-(power-3))>1e-5:
            raise ArithmeticError("Observed tidal or Maxwell exponent differs")
        tails[f'{power:g}']=dict(field_amplitude_power=power,expected_tidal_exponent=2*power-4,
                                 measured_tidal_exponent=slope,expected_field_exponent=power-3,
                                 measured_field_exponent=maxwell,by_depth=rows)
    shut=brentq(lambda y:integrate(du_dy,0,y)-3.,0.,6.,xtol=2e-14)
    local=brentq(lambda y:integrate(du_dy,0,y)-2.,0.,6.,xtol=2e-14)
    after=[sandwich(integrate(du_dy,0,y)) for y in (shut+1e-10,shut+1,shut+5,shut+10)]
    if max(after)!=0.:raise ArithmeticError("Compact wave persisted after switch-off")
    rr=[]
    for y in depths:
        a,ay,ayy,phi1,b=profile(y)
        rr.append(dict(y=y,r=a,scalar_offset=phi1,radial_metric_excess=b,
                       boosted_warp_anisotropy=abs(ayy)/a**2))
    l,r=rr[-2:]
    exponents={key:logarithmic_slope(l['r'],r['r'],l[key],r[key])
               for key in ('scalar_offset','radial_metric_excess','boosted_warp_anisotropy')}
    if max(abs(exponents[k]-n) for k,n in [('scalar_offset',2/K),('radial_metric_excess',4/K-2),
                       ('boosted_warp_anisotropy',4/K-2)])>1e-8:
        raise ArithmeticError("Horizon asymptotics differ")
    return dict(k=K,E=ENERGY,P=MOMENTUM,lambda_coupling=LAMBDA,wall_a=1.,
                horizon_proper_time=time_remaining(0.),
                remaining_time_over_a_limit=1/(K*ENERGY),u_times_a_limit=MOMENTUM/(K*ENERGY),
                sandwich=dict(support=[1.,3.],center_y=local,center_proper_time=integrate(dt_dy,0,local),
                              center_wave_ricci=(MOMENTUM/math.exp(2*warp(local)))**2,
                              exit_y=shut,exit_proper_time=integrate(dt_dy,0,shut),
                              proper_time_remaining_at_exit=time_remaining(shut),wave_at_later_depths=after),
                proper_time_geodesic_controls={str(y):geodesic_check(y) for y in (4,8,12)},
                wave_tails=tails,background_samples=rr,background_exponents=exponents,
                chart_controls={str(p):[chart_test(y,p) for y in (1.,8.,20.,32.,None)] for p in (2,3)},
                transverse_localization=[transverse_check(x,y) for x,y in ((1.,1.),(2.,2.),(1.,2.),(2.,.5))])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(),allow_nan=False,indent=2)+'\n')
