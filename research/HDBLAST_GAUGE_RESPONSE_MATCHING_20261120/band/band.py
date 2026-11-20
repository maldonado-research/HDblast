"""Exact photon band response and its local derivative expansion on a prescribed roll."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
BANDS = ((48, 96), (48, 192), (48, 384), (48, 768), (128, 384), (384, 768))


def load_code(path):
    spec = importlib.util.spec_from_file_location('historical_gauge_profiles', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def series_product(a, b):
    """Four Taylor coefficients at each sample, with factorials already divided out."""
    return np.stack([sum(a[j]*b[n-j] for j in range(n+1)) for n in range(4)])


def force_derivatives(branch, q):
    fit = branch.fitted
    jets = np.stack([fit.spline(branch.top+q,i)/math.factorial(i) for i in range(4)])
    phi, f = jets[:,:,0], jets[:,:,1]
    gap = -phi/fit.pm
    gap[0] += 1
    tilt = fit.b*phi
    tilt[0] += 1
    numerator = series_product(series_product(gap,gap),tilt)
    denominator = fit.v0*series_product(f,f)
    potential = []
    for i in range(4):
        convolution = sum(denominator[j]*potential[i-j] for j in range(1,i+1))
        potential.append((numerator[i]-convolution)/denominator[0])
    return -3*(potential[1]-branch.du_top), -6*potential[2], -18*potential[3]


def background(branch, roll, norm, samples):
    """Reintegrate the scalar in switch position, with time recovered from dq/dx."""
    start = roll.points['switch_start']
    d = .6*(roll.qmin-roll.q0)
    def rhs(u,state):
        v,a,h,eta = state
        q = start['q']+d*u
        force = float(branch.potential(q)[1])
        dt = d/v
        return dt*np.array([force-3*h*v,a*h,1/a**2-v*v/2,1/a])
    u = np.linspace(0.,1.,samples)
    saved = solve_ivp(rhs,(0.,1.),[start['velocity'],start['a'],start['h'],0.],
                      method='DOP853',rtol=3e-13,atol=3e-15,max_step=.008,t_eval=u)
    if not saved.success:
        raise ArithmeticError(saved.message)
    v,a,h,eta = saved.y
    error = np.max(abs(eta-roll.eta(u)))
    if error>1e-10:
        raise ArithmeticError('Reconstructed conformal clock disagrees with archived roll')
    q = start['q']+d*u
    F,Fq,Fqq = force_derivatives(branch,q)
    if np.max(abs(F-branch.potential(q)[1]))>2e-11:
        raise ArithmeticError('Force derivative normalization differs')
    Hc = a*h
    dq = a*v
    ddq = -2*Hc*dq+a*a*F
    Hp = Hc**2+1-dq**2/2
    Hpp = 2*Hc*Hp-dq*ddq
    q3 = -2*Hp*dq-2*Hc*ddq+a*a*(2*Hc*F+Fq*dq)
    q4 = (-2*Hpp*dq-4*Hp*ddq-2*Hc*q3
          +a*a*((2*Hp+4*Hc**2)*F+4*Hc*Fq*dq+Fqq*dq*dq+Fq*ddq))
    u1,u2,u3,u4 = dq/d,ddq/d,q3/d,q4/d
    # Log bump = -1/u-1/(1-u)-log(norm). Every derivative vanishes at its ends.
    derivatives = np.zeros((4,len(u)))
    x = u[1:-1]
    bump = np.exp(-1/x-1/(1-x))/norm
    first = 1/x**2-1/(1-x)**2
    second = -2/x**3-2/(1-x)**3
    third = 6/x**4-6/(1-x)**4
    derivatives[:,1:-1] = [bump,bump*first,bump*(first**2+second),
                            bump*(first**3+3*first*second+third)]
    b,b1,b2,b3 = derivatives
    jet = np.stack((b*u1,b1*u1*u1+b*u2,b2*u1**3+3*b1*u1*u2+b*u3,
                    b3*u1**4+6*b2*u1**2*u2+b1*(3*u2**2+4*u1*u3)+b*u4))
    return dict(u=u,a=a,eta=eta,eta_u=1/u1,jet=jet,grad_unit=b/d,clock_error=float(error))


def regulator_windows():
    labels=np.arange(49,769,dtype=float)
    weights=labels*(labels**2-1)/np.pi**2
    static=weights/(2*labels)**2
    windows={f'{lo+1}-{hi}':((labels>lo)&(labels<=hi)).astype(float) for lo,hi in BANDS}
    matches={}
    for hi in (96,192):
        reference=windows[f'49-{hi}']
        target=float(static@reference)
        def taper(edge):
            argument=np.clip((labels-48)/(edge-48),0.,1.)
            return .5*(1+np.cos(np.pi*argument))
        edge=brentq(lambda k:float(static@taper(k))-target,hi,3*hi,xtol=2e-12)
        key=f'cosine_match_49-{hi}'
        windows[key]=taper(edge)
        actual=float(static@windows[key])
        if abs(actual-target)>1e-9:
            raise ArithmeticError('Smooth and sharp spectral susceptibilities differ')
        matches[key]=dict(taper_edge=edge,sharp_maximum=hi,sharp_W2=target,tapered_W2=actual)
    return windows,matches


def exact_energy_force(roll, gauge, u, stretch, delta, rtol, windows):
    """Phase-removed covariance columns; band D and its conjugate J at every sample."""
    output={name:np.zeros((2,len(u))) for name in windows}
    normalization_error=0.
    for lo in range(49,769,48):
        labels=np.arange(lo,min(lo+48,769))
        size=len(labels)
        def rhs(z,y):
            alpha,beta=y[:size],y[size:]
            pump=delta*gauge.raw_bump(z)/gauge.NORM*np.exp(2j*labels*stretch*roll.eta(z))
            return np.r_[pump*beta.conj(),pump*alpha.conj()]
        sol=solve_ivp(rhs,(0.,1.),np.r_[np.ones(size,complex),np.zeros(size,complex)],method='DOP853',
                      t_eval=u,rtol=rtol,atol=rtol*1e-4,
                      max_step=min(.02,1.5/(labels[-1]*stretch*roll.eta(1.))))
        if not sol.success:
            raise ArithmeticError(sol.message)
        alpha,beta=sol.y[:size],sol.y[size:]
        n=abs(beta)**2
        normalization_error=max(normalization_error,float(np.max(abs(abs(alpha)**2-n-1))))
        mreal=np.real(np.exp(-2j*labels[:,None]*stretch*roll.eta(u))*alpha*beta)
        weight=labels*(labels**2-1)/np.pi**2
        for name,values in output.items():
            shell_weight=weight*windows[name][labels-49]
            if np.any(shell_weight):
                values[0] += shell_weight@n
                values[1] += 2*shell_weight@mreal
    if normalization_error>5e-10:
        raise ArithmeticError('Exact covariance lost its normalization')
    return output,normalization_error


def inspect_band(values, window, bg, delta, stretch, hpl):
    s,s1,s2,s3 = bg['jet']*delta/np.array([stretch**i for i in range(1,5)])[:,None]
    u,a = bg['u'],bg['a']
    measure=stretch*bg['eta_u']
    g=delta*bg['grad_unit']
    labels=np.arange(49,769,dtype=float)
    weight=window*labels*(labels**2-1)/np.pi**2
    W2=sum(weight/(2*labels)**2)
    W4=sum(weight/(2*labels)**4)
    e2=W2*s**2
    e4=W4*(3*s**4+s1**2-2*s*s2)
    j2=2*W2*s1
    j4=2*W4*(6*s**2*s1-s3)
    D,J=values
    # Density in rho/(3 Mpl² H0²) and acceleration source -3S in parent units.
    scale=hpl**2/(3*a**4)
    def rms(x):
        return float(np.sqrt(simpson(x*x*measure,x=u)/simpson(measure,x=u)))
    def errors(exact,second,fourth):
        reference=rms(exact)
        return dict(exact_rms=reference,
                    after_second_order_rms_over_exact=rms(exact-second)/reference,
                    after_fourth_order_rms_over_exact=rms(exact-second-fourth)/reference)
    def integral(y):return float(simpson(y*measure,x=u))
    balance=abs(D[-1]-integral(s*J))
    zero2=abs(integral(s*j2))
    zero4=abs(integral(s*j4))
    # A loose numerical floor allows for the finite regularity of the historical potential spline.
    if balance>max(2e-9,1e-8*np.max(abs(D))) or max(zero2,zero4)>3e-5:
        raise ArithmeticError(f'Work balance failed: {balance}, {zero2}, {zero4}')
    points={}
    for z in (.25,.5,.75,1.):
        i=int(round(z*(len(u)-1)))
        points[str(z)]=dict(D_exact=float(D[i]),D2=float(e2[i]),D4=float(e4[i]),
                            scalar_acceleration_exact=float(-3*scale[i]*g[i]*J[i]),
                            scalar_acceleration_2=float(-3*scale[i]*g[i]*j2[i]),
                            scalar_acceleration_4=float(-3*scale[i]*g[i]*j4[i]))
    return dict(exit_D=float(D[-1]),peak_D_sampled=float(np.max(D)),
                integrated_exact_work=integral(s*J),abs_work_balance_error=balance,
                absolute_local_integrated_work_2=zero2,absolute_local_integrated_work_4=zero4,
                density=errors(scale*D,scale*e2,scale*e4),
                scalar_acceleration=errors(-3*scale*g*J,-3*scale*g*j2,-3*scale*g*j4),
                points=points)


def calculate(gauge_path, prior, archive_path, samples=1025, rtol=2e-11):
    if samples<257 or (samples-1)%4:
        raise ValueError('Use at least 257 samples with four equal subintervals')
    gauge=load_code(gauge_path)
    branches,history=gauge.load(prior)
    archived=json.loads(archive_path.read_bytes())
    output={}
    windows,matches=regulator_windows()
    for key,branch in branches.items():
        roll=gauge.Roll(branch,history[key],2049)
        bg=background(branch,roll,gauge.NORM,samples)
        cases={}
        for delta,stretch in ((.1,1),(.5,1),(.1,4)):
            data,wronskian=exact_energy_force(roll,gauge,bg['u'],stretch,delta,rtol,windows)
            rows={}
            for name,values in data.items():
                rows[name]=inspect_band(values,windows[name],bg,delta,stretch,roll.Hpl)
            case=dict(delta_log_I=delta,conformal_stretch=stretch,wronskian_error=wronskian,bands=rows)
            if stretch==1:
                # GAUGE.json tabulates the complete primary spectrum through L=1536.
                if delta==.1:
                    c=np.array(archived['branches'][key]['primary_spectrum']['cumulative_D'])
                    expected=c[768-2]-c[48-2]
                    half=rows['49-768']['exit_D']
                    if abs(half-expected)>2e-7:
                        raise ArithmeticError('Outgoing shell energy differs from the historical measurement')
                    case['archived_band_D_error']=abs(half-expected)
                    case['fraction_of_all_archived_out_energy_in_49_768']=half/c[-1]
            cases[f'delta_{delta:g}_stretch_{stretch}']=case
        output[key]=dict(clock_error=bg['clock_error'],reference_duration=float(roll.eta(1.)),cases=cases)
    return dict(samples=samples,rtol=rtol,regulator_matches=matches,branches=output)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gauge',type=Path,default=ROOT/'prior/gauge/gauge.py')
    parser.add_argument('--archive',type=Path,default=ROOT/'prior/gauge/GAUGE.json')
    parser.add_argument('--prior',type=Path,default=ROOT/'prior')
    parser.add_argument('--samples',type=int,default=1025)
    parser.add_argument('--rtol',type=float,default=2e-11)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.gauge,args.prior,args.archive,args.samples,args.rtol),
                                     allow_nan=False,indent=2)+'\n')
