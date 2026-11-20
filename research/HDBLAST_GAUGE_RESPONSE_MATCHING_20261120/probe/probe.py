"""Static retarded response of the added transverse photon harmonics of a closed sphere."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import digamma, roots_legendre, zeta


REFERENCE = 48
CALIBRATION_FREQUENCY = 32.
CUTOFFS = (64, 96, 128, 192, 384, 768, 1536)


def susceptibility(cutoff, omega, match_omega=CALIBRATION_FREQUENCY):
    """Four times the mode sum converting s' into the real weak covariance."""
    labels = np.arange(REFERENCE+1, cutoff+1, dtype=float)
    squares = labels**2
    weight = labels-1/labels
    x = (omega/2)**2
    base = (match_omega/2)**2
    if max(x,base) >= squares[0]:
        raise ValueError('Probe requires frequencies strictly below the added band')
    static = float(np.sum(weight))
    slope = float(np.sum(weight/squares))
    first_residual = float(np.sum(x*weight/(squares-x)))
    second_residual = float(np.sum(x*x*weight/(squares*(squares-x))))
    interp_residual = float(np.sum(x*(x-base)*weight/((squares-x)*(squares-base))))
    exact = static+first_residual
    directly_summed = np.sum(weight/(1-x/squares))
    at_base = np.sum(weight/(1-base/squares))
    if not math.isclose(exact,directly_summed,abs_tol=1e-8,rel_tol=1e-14):
        raise ArithmeticError('Finite response differs')
    if not math.isclose((at_base-static)*x/base+interp_residual,first_residual,
                        abs_tol=2e-9,rel_tol=3e-12):
        raise ArithmeticError('Independent finite two-point interpolation check differs')
    return dict(chi=exact,static_coefficient=static,derivative_coefficient=slope,
                after_static=first_residual,after_value_and_slope=second_residual,
                after_two_frequencies=interp_residual)


def finite_part(omega):
    """Infinite-band, twice Taylor-subtracted response, independently by digamma."""
    x=(omega/2)**2
    a=REFERENCE+1.
    s=omega/2
    return float((x-1)*(digamma(a)-.5*(digamma(a-s)+digamma(a+s)))+x*zeta(3,a))


def tail_bound(omega, cutoff, match_omega=None):
    x=(omega/2)**2
    base=(match_omega/2)**2 if match_omega is not None else 0.
    return abs(x*(x-base))/(2*cutoff**2*(1-x/(cutoff+1)**2)*(1-base/(cutoff+1)**2))


def frequency_results():
    output = {}
    for omega in (8.,16.,32.,48.,64.,80.):
        x=(omega/2)**2
        infinite_slope=finite_part(omega)
        infinite_interp=infinite_slope-x/(CALIBRATION_FREQUENCY/2)**2*finite_part(CALIBRATION_FREQUENCY)
        rows={}
        for k in CUTOFFS:
            row=susceptibility(k,omega)
            row['uncomputed_slope_tail_bound']=tail_bound(omega,k)
            row['uncomputed_two_frequency_tail_bound']=tail_bound(omega,k,CALIBRATION_FREQUENCY)
            if abs(infinite_slope-row['after_value_and_slope'])>row['uncomputed_slope_tail_bound']+1e-9:
                raise ArithmeticError('Analytic Taylor-subtracted tail escaped the bound')
            if abs(infinite_interp-row['after_two_frequencies'])>row['uncomputed_two_frequency_tail_bound']+1e-9:
                raise ArithmeticError('Analytic interpolated tail escaped the bound')
            rows[str(k)]=row
        output[f'{omega:g}']=dict(cutoffs=rows,infinite_after_value_and_slope=infinite_slope,
                                  infinite_after_two_frequencies=infinite_interp)
    return output


def pulse_harmonics(width):
    """First three Fourier coefficients of sin^6(pi eta/T), plus a constant of 10/32."""
    omega=2*np.pi*np.arange(1.,4.)/width
    return omega,np.array([-15.,6.,-1.])/32


def drive(times,width):
    times=np.asarray(times)
    frequencies,a = pulse_harmonics(width)
    active=times<=width
    phase=times[...,None]*frequencies
    # Return the derivatives of ln I for a unit-height excursion.
    return [np.where(active,np.sum(a*frequencies**n*np.cos(phase+n*np.pi/2),axis=-1),0.)
            for n in (1,2,4)]


def exact_pulse_response(times,width,cutoff):
    """Zero initial data. At switch-off match each real solution and its conjugate."""
    labels=np.arange(REFERENCE+1,cutoff+1,dtype=float)
    k=2*labels
    w=labels*(labels**2-1)
    frequencies,amplitudes=pulse_harmonics(width)
    coefficients= -amplitudes*frequencies**2/(k[:,None]**2-frequencies**2)
    v0=coefficients.sum(axis=1)
    def during(t):
        phase=t[:,None]*k
        slow=t[:,None]*frequencies
        real=np.cos(slow)@coefficients.T-np.cos(phase)*v0
        derivative=(-frequencies*np.sin(slow))@coefficients.T+(k*v0)*np.sin(phase)
        return real,derivative
    capped=np.minimum(times,width)
    real,rdot=during(capped)
    real_end,rdot_end=during(np.array([width]))
    past=times>width
    phase=(times[past,None]-width)*k
    real[past]=real_end*np.cos(phase)+(rdot_end/k)*np.sin(phase)
    # At T, s=0 and X'=2L Y, where Y=Im m.
    born_occupation=(real_end[0]**2+(rdot_end[0]/k)**2)
    return real@w*4,float(w@born_occupation/np.pi**2)


def full_pulse(width,cutoff,amplitude,times):
    """Integrate raw position and momentum covariances of the finite canonical oscillators."""
    labels=np.arange(REFERENCE+1,cutoff+1,dtype=float)
    w=labels*(labels**2-1)
    n=len(w)
    initial=np.zeros(3*n)
    # x,y,z are the two shifted diagonal variances and the real symmetrized covariance.
    def rhs(t,vec):
        X,Y,Z=vec.reshape(3,n)
        if 0<t<width:
            phase=np.pi*t/width
            s=amplitude*(6*np.pi/width)*np.sin(phase)**5*np.cos(phase)
        else:
            s=0.
        return np.r_[2*s*(1+X)+2*labels*Z,-2*s*(1+Y)-2*labels*Z,labels*(Y-X)]
    sol=solve_ivp(rhs,(0.,times[-1]),initial,t_eval=times,method='DOP853',
                  rtol=1e-11,atol=2e-15,max_step=min(.004,.8/labels[-1]))
    if not sol.success:
        raise ArithmeticError(sol.message)
    X,Y,Z=sol.y.reshape(3,n,-1)
    occupation=(X+Y)/4
    pure_error=np.max(abs((1+X)*(1+Y)-Z**2-1))
    if pure_error>1e-10:
        raise ArithmeticError('Canonical variance purity check failed')
    work=4*w@((X-Y)/4)/amplitude
    scaled_energy=(w@occupation)/(np.pi**2*amplitude**2)
    return work,scaled_energy


def outgoing_tail(width,cutoff):
    omega=2*np.pi/width
    coefficient=(45/4)*(omega/2)**6
    gap=1-(3*omega/(2*(cutoff+1)))**2
    if gap<=0:
        raise ValueError('Outgoing tail requires the cutoff to exceed the pulse harmonics')
    return float(coefficient**2/(2*np.pi**2*cutoff**8*gap**6))


def pulse_results():
    width=.25
    times=np.r_[np.linspace(0.,width,601),np.linspace(width,.40,301)[1:]]
    mask=times<=width
    rows={}
    for cutoff in (64,96,128,192,384):
        exact,outgoing=exact_pulse_response(times,width,cutoff)
        s,s1,s3=drive(times,width)
        co=susceptibility(cutoff,0.)
        static=co['static_coefficient']*s1
        second=static-.25*co['derivative_coefficient']*s3
        relative=lambda v: float(np.linalg.norm(v[mask])/np.linalg.norm(exact[mask]))
        rows[str(cutoff)]=dict(relative_l2_static_error=relative(exact-static),
                               relative_l2_value_and_slope_error=relative(exact-second),
                               peak_retarded_signal_during=float(np.max(abs(exact[mask]))),
                               peak_free_signal_after=float(np.max(abs(exact[~mask]))),
                               maximum_tail_signal_over_peak=float(np.max(abs(exact[~mask]))/np.max(abs(exact[mask]))),
                               outgoing_born_a4_density_per_amplitude_squared=outgoing)
        nodes,weights=roots_legendre(256)
        eta=(1+nodes)*width/2
        work_force,_=exact_pulse_response(eta,width,cutoff)
        s_quad,s1_quad,s3_quad=drive(eta,width)
        measure=weights*width/(4*np.pi**2)
        from_work=float(measure@(s_quad*work_force))
        local_work=float(measure@(s_quad*(co['static_coefficient']*s1_quad-
                                                .25*co['derivative_coefficient']*s3_quad)))
        if abs(from_work-outgoing)>1.e-7 or abs(local_work)>1.e-7:
            raise ArithmeticError('Pulse energy ledger failed')
        rows[str(cutoff)].update(outgoing_from_retarded_work=from_work,
                                 outgoing_from_two_local_terms=local_work,
                                 omitted_outgoing_born_tail_upper=outgoing_tail(width,cutoff))
        if cutoff==128:
            checks={}
            for amp in (.02,.002,.0002):
                full,energy=full_pulse(width,cutoff,amp,times)
                checks[f'{amp:g}']=dict(relative_l2_force_error=relative(full-exact),
                                     endpoint_a4_density_per_amplitude_squared=float(energy[-1]),
                                     relative_outgoing_born_error=float(energy[-1]/outgoing-1),
                                     post_switch_fractional_energy_drift=float(np.ptp(energy[~mask])/energy[-1]))
            last=checks['0.0002']
            if last['relative_l2_force_error']>1e-7 or abs(last['relative_outgoing_born_error'])>1e-6:
                raise ArithmeticError('Weak canonical evolution missed retarded limit')
            rows[str(cutoff)]['independent_canonical_checks']=checks
    ordered=list(rows.values())
    for current,later in zip(ordered,ordered[1:]):
        gain=later['outgoing_born_a4_density_per_amplitude_squared']-current['outgoing_born_a4_density_per_amplitude_squared']
        if not 0<=gain<=current['omitted_outgoing_born_tail_upper']:
            raise ArithmeticError('Accumulated outgoing energy exceeds analytic tail ceiling')
    return dict(profile='log I = amplitude sin^6(pi eta/T) on 0<eta<T, zero outside',
                conformal_duration=width,post_switch_last_time=float(times[-1]),
                comparison='sum of squared errors over 601 uniformly spaced points, including both endpoints',
                cutoffs=rows)


def calculate():
    return dict(reference=REFERENCE,first_added_pair_frequency=2*(REFERENCE+1),
                matching_frequencies=[0.,CALIBRATION_FREQUENCY],
                normalization="4 sum L(L^2-1) Re(m_L) = chi * s' in conformal Fourier variables",
                frequencies=frequency_results(),compact_pulse=pulse_results())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(),indent=2,allow_nan=False)+'\n')
