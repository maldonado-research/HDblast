"""Charged Dirac decay phase space in a closed, finite spatial volume."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import brentq
from scipy.special import expit


Q = 1e8
RI, RF = .95, 1.05
CENTER, WIDTH = .77, .03


def mass_ratio(psi):
    z = (psi-CENTER)/WIDTH+.5
    if z <= 0:
        return RI
    if z >= 1:
        return RF
    return RI+(RF-RI)*expit(1/(1-z)-1/z)


def intrinsic_levels(two_j, shift):
    """Each eigenvalue has a further right multiplicity 2j+1."""
    n=two_j
    edge=[n+1.5-shift,n+1.5+shift]
    ms=np.arange(-(n-1),n,2)/2 if n else np.empty(0)
    root=np.sqrt((n+1)**2+shift**2-4*ms*shift)
    return np.sort(np.r_[edge,.5-root,.5+root])


def matrix_levels(two_j, shift):
    n=two_j
    m=np.arange(n+1)-n/2
    jp=np.diag(np.sqrt((n/2-m[:-1])*(n/2+m[:-1]+1)),-1)
    jx=(jp+jp.T)/2
    jy=(jp-jp.T)/(2j)
    sx=np.array([[0,1],[1,0]])
    sy=np.array([[0,-1j],[1j,0]])
    sz=np.diag([1.,-1.])
    D=1.5*np.eye(2*(n+1))+2*(np.kron(jx,sx)+np.kron(jy,sy)+np.kron(np.diag(m),sz))
    return np.linalg.eigvalsh(D-shift*np.kron(np.eye(n+1),sz))


def gap(shift):
    # Edges n+3/2-|s| attain the minimum. The other branches cost >=3/2.
    target=abs(shift)-1.5
    return abs(max(0,math.floor(target+.5))-target)


def spectral_checks():
    largest=0.
    for n in range(15):
        for s in [-7.2,0.,.4,1.5,math.pi,19.]:
            largest=max(largest,float(np.max(abs(matrix_levels(n,s)-intrinsic_levels(n,s)))))
    gap_error=0.
    for s in np.r_[np.linspace(-16.,16.,65),[1.50001,2.416,13.991]]:
        brute=min(float(np.min(abs(intrinsic_levels(n,s)))) for n in range(19))
        gap_error=max(gap_error,abs(brute-gap(s)))
    for n in range(7):
        val=n+1.5
        for sign in [-1.,1.]:
            degeneracy=sum((k+1)*np.count_nonzero(np.isclose(intrinsic_levels(k,0.),sign*val))
                           for k in range(n+2))
            if degeneracy!=(n+1)*(n+2):
                raise ArithmeticError("Intrinsic sign multiplicity")
    if largest>1e-12 or gap_error>1e-12:
        raise ArithmeticError("Closed Dirac spectrum")
    return dict(matrix_absolute_error=largest, minimum_gap_absolute_error=gap_error,
                neutral_sign_multiplicities_checked_through_n=6)


def read_history(path):
    data=json.loads(path.read_bytes())
    row=data['runs']['neutral_precursor']
    columns=dict(zip(data['trace_columns'],np.array(row['trace']).T,strict=True))
    time=columns['t']
    fields={key:PchipInterpolator(time,columns[key]) for key in ('a','psi','B','X','H','K')}
    tmax=float(row['maximum_tau'])
    def sample(t):
        a,psi,B,X,H,K=(float(fields[key](t)) for key in ('a','psi','B','X','H','K'))
        return dict(t=t,a=a,psi=psi,B=B,X=X,H=H,K=K,A=a*a*B/2,R=float(mass_ratio(psi)))
    tc=brentq(lambda t:float(fields['psi'](t))-CENTER,8.,tmax)
    tb=brentq(lambda t:float(fields['B'](t)),tc,tmax)
    late=brentq(lambda t:float(fields['psi'](t))-(CENTER+WIDTH/2),tc,tmax)
    return fields, {'opening':sample(tc),'magnetic_zero':sample(tb),
                    'full_mass':sample(late),'turnaround':sample(tmax)}


def admissible_count(mu,a,R):
    highest=math.floor(a*mu*math.sqrt(max(0.,R*R-1))-1.5)
    shells=max(0,highest+1)
    states=2*shells*(shells+1)*(shells+2)//3
    return shells,states


def phase_bounds(sample,tmax):
    a,R,K=sample['a'],sample['R'],sample['K']
    available=tmax-sample['t']
    radial=math.sqrt(R*R-1)
    beta=radial/R
    p0=1.5/a
    floor=p0/radial
    qlo=math.sqrt(3)*Q*floor
    mu1=1/(math.sqrt(3)*Q)
    local=[]
    for mu in [mu1,floor,1.,10.,100.,1e3,1e5]:
        E=math.hypot(p0,mu)
        slopes=[-sample['H']*p0+orientation*mu*Q*sample['X'] for orientation in [-1,1]]
        nh,gn=admissible_count(mu*(1+1e-14),a,R)
        local.append(dict(mu=mu,q_at_fixed_Q=math.sqrt(3)*Q*mu,largest_edge_adiabatic_ratio=
                          max(abs(v) for v in slopes)*mu/(4*E**3),kinematic_shells=nh,
                          individual_particle_sign_slots_per_Dirac_species=gn,
                          Epair_floor=2*E,parent_mass=2*mu*R))
    return dict(zero_free_pfloor=p0,minimum_mass_for_kinematic_access=floor,
                minimum_q_at_fixed_Q=qlo, minimum_S0_weak_alpha_at_LATE_K=qlo**2/(4*math.pi*K),
                minimum_S0_weak_alpha_at_initial_K=qlo**2/(4*math.pi*.001),
                remaining_time=available,minimum_mu_to_resolve_threshold_in_remaining_time=
                1/(2*(R-1)*available),local_continuum_level_gap=2*beta/a,
                nominal_max_time_for_one_level_width=a/(2*beta), matched_cases=local,
                minimum_S0_one_initial_parent_at_mu_floor=floor*RI/(math.pi**2*2**3*1e-4))


def filtered(mu,R,a,T):
    M=2*mu*R
    # The Gaussian is negligible beyond 12 inverse widths from the resonance.
    zlo=max(-12.,T*(2*mu-M))
    zhi=12.
    def continuum_integrand(z):
        E=(M+z/T)/2
        p=math.sqrt(max(0.,(E-mu)*(E+mu)))
        # Integral of p² dp beta² times a normalized Gaussian in 2E−M.
        return p**3/(2*E)*math.exp(-z*z)/math.sqrt(math.pi)
    cont=quad(continuum_integrand,zlo,zhi,epsabs=1e-30,epsrel=1e-10,limit=150)[0]/(math.pi*M)
    def upper_n(z):
        E=max(mu,(M+z/T)/2)
        return a*math.sqrt((E-mu)*(E+mu))-1.5
    nlo=max(0,math.ceil(upper_n(zlo))-1)
    nhi=max(nlo+1,math.ceil(upper_n(zhi))+1)
    n=np.arange(nlo,nhi+1,dtype=float)
    p=(n+1.5)/a
    E=np.hypot(p,mu)
    z=T*(2*E-M)
    states=2*(n+1)*(n+2)
    lines=states*(p/E)**2*np.exp(-z*z)*T/math.sqrt(math.pi)
    compact=float(math.pi*np.sum(lines)/(M*2*math.pi**2*a**3))
    infinite=mu*R*(1-R**-2)**1.5/(4*math.pi)
    return dict(width_T=T,mu=mu,compact_over_broadened_flat=compact/cont,
                broadened_flat_over_flat_FGR=cont/infinite,
                compact_over_flat_FGR=compact/infinite,
                compact_Gamma_per_y_squared=compact,flat_FGR_per_y_squared=infinite,
                kinetic_excess_resolution=T*(M-2*mu),
                lines_within_one_width=int(np.count_nonzero(abs(z)<1)))


def edge_sweep(fields,zero,mu,span,rtol=4e-11,reverse=False,static=False):
    # The extremal spin state has an A-independent spatial eigenvector.
    def A(t):
        if 'potential' in fields:
            return float(fields['potential'](t))
        return float(fields['a'](t)**2*fields['B'](t)/2)
    tB=zero['t']
    charge=math.sqrt(3)*Q*mu
    target=1.5/(Q*mu)
    tc=brentq(lambda t:A(t)-target,tB-.3,tB,xtol=2e-14)
    a=float(fields['a'](tc))
    adot=float(fields['a'](tc,1))
    dA=(float(fields['potential'](tc,1)) if 'potential' in fields else
        adot*a*float(fields['B'](tc))+.5*a*a*float(fields['B'](tc,1)))
    v=abs(Q*mu*dA/a)
    scale=math.sqrt(v)
    def momentum(y):
        if static:
            return 0.
        t=tc+(-y if reverse else y)/scale
        return (1.5-Q*mu*A(t))/float(fields['a'](t))
    def vectors(y):
        p=momentum(y)
        return np.linalg.eigh(np.array([[p,mu],[mu,-p]]))[1]
    def rhs(y,u):
        p=momentum(y)
        return -1j*np.array([p*u[0]+mu*u[1],mu*u[0]-p*u[1]])/scale
    incoming=vectors(-span)[:,0].astype(complex)
    sol=solve_ivp(rhs,(-span,span),incoming,method='DOP853',rtol=rtol,atol=rtol/100,
                  max_step=.1 if rtol<1e-9 else .22)
    if not sol.success:
        raise ArithmeticError(sol.message)
    out=sol.y[:,-1]
    measured=float(abs(np.vdot(vectors(span)[:,1],out))**2)
    pred=math.exp(-math.pi*mu**2/v)
    return dict(mu=mu,q=charge,time_of_edge_crossing=tc,time_to_magnetic_zero=tB-tc,
                sweep_speed=v,measured_trace_to_electric_speed_ratio=v/(Q*mu*abs(float(fields['X'](tc)))),
                linear_landau_zener=pred,finite_instantaneous_probability=measured,
                norm_error=float(abs(np.vdot(out,out)-1)),dimensionless_half_width=span,
                proper_time_half_width=span/scale,adjacent_crossing_time=1/(Q*mu*abs(dA)),
                crossing_spacing_in_formation_units=1/(a*scale))


def calculate(path):
    fields,events=read_history(path)
    zero=events['magnetic_zero']
    end=events['turnaround']
    result=dict(events=events,checks=spectral_checks(),bounds=phase_bounds(zero,end['t']))
    a0,R0=zero['a'],zero['R']
    rows=[]
    mu1=1/(math.sqrt(3)*Q)
    for mu in [mu1, .25,.6,1.,10.,100.,1000.]:
        for T in [.5,2.,16.,128.]:
            rows.append(filtered(mu,R0,a0,T))
    result['gaussian_static_tests']=rows
    full=np.linspace(events['opening']['t']+1e-5,end['t'],5001)
    aa=fields['a'](full)
    A=aa*aa*fields['B'](full)/2
    rr=np.array([mass_ratio(v) for v in fields['psi'](full)])
    available=mu1*np.sqrt(np.maximum(rr*rr-1,0))
    price=np.array([gap(float(x)/math.sqrt(3)) for x in A])/aa
    accessible=available>=price
    if np.any(accessible):
        raise AssertionError("Unexpected q=1 matched crossing")
    # Independent bound over the monotone-interpolated expansion trace:
    # sample endpoints bound each individual PCHIP on its sampling interval.
    tknots=fields['a'].x
    mask=(tknots>=events['opening']['t']-.11)&(tknots<=end['t']+.11)
    amax=float(max(fields['a'](tknots[mask])))
    bmax=float(max(abs(fields['B'](tknots[mask]))))
    Amax_bound=amax*amax*bmax/2
    mingap=max(0.,1.5-Amax_bound/math.sqrt(3))/amax
    ongrid=1.5/(Q*abs(A)+aa*np.sqrt(np.maximum(rr*rr-1,0)))
    result['magnetic_history']=dict(q1_mu=mu1,q1_grid_slices=len(full),
                minimal_resolved_q1_momentum=float(min(price)),bounded_q1_momentum=mingap,
                maximal_allowed_q1_momentum=mu1*math.sqrt(RF*RF-1),
                grid_minimum_necessary_mu_at_fixed_Q=float(min(ongrid)),
                grid_minimum_necessary_q_at_fixed_Q=math.sqrt(3)*Q*float(min(ongrid)),
                time_of_grid_minimum=float(full[np.argmin(ongrid)]),
                max_abs_A_bound=Amax_bound)
    analytic=filtered(1e4,R0,a0,1.)
    if abs(analytic['compact_over_flat_FGR']-1)>.001 or abs(
            analytic['broadened_flat_over_flat_FGR']-1)>.001:
        raise ArithmeticError('Flat-space continuum normalization')
    result['checks']['many_shell_flat_test']=analytic
    sweep={}
    for label,mu in [('kinematic_floor',result['bounds']['minimum_mass_for_kinematic_access']),
                     ('mass_1000',1000.),('mass_100000',1e5)]:
        short=edge_sweep(fields,zero,mu,8.)
        long=edge_sweep(fields,zero,mu,16.)
        sweep[label]=dict(short=short,long=long)
    mass=result['bounds']['minimum_mass_for_kinematic_access']
    primary=sweep['kinematic_floor']['long']
    fast=edge_sweep(fields,zero,mass,16.,rtol=2e-9)
    back=edge_sweep(fields,zero,mass,16.,reverse=True)
    static=edge_sweep(fields,zero,mass,16.,static=True)
    # Validate the interpolated slope on a second reconstruction that explicitly
    # satisfies A_t=aX at each measured point instead of fitting B.
    tt=fields['a'].x
    aa=fields['a'](tt)
    herm=dict(fields,potential=CubicHermiteSpline(tt,aa**2*fields['B'](tt)/2,aa*fields['X'](tt)))
    zero_herm=dict(zero,t=brentq(herm['potential'],zero['t']-.3,zero['t']+.3))
    other=edge_sweep(herm,zero_herm,mass,16.)
    tests=dict(potential_interpolant_probability_difference=abs(other['finite_instantaneous_probability']-
                primary['finite_instantaneous_probability']),
               solver_tolerance_difference=abs(fast['finite_instantaneous_probability']-
                primary['finite_instantaneous_probability']),
               reversed_history_difference=abs(back['finite_instantaneous_probability']-
                primary['finite_instantaneous_probability']),static_transition=static['finite_instantaneous_probability'])
    if (max(tests.values())>2e-7 or primary['norm_error']>2e-8
            or abs(primary['finite_instantaneous_probability']-primary['linear_landau_zener'])>2e-5):
        raise ArithmeticError(tests)
    result['field_driven_edge']=dict(cases=sweep,checks=tests)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history',type=Path,default=Path(__file__).resolve().parents[1]/'prior/CHARGED.json')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.history),indent=2,sort_keys=True,allow_nan=False)+'\n')
