"""Finite-time quantum evolution of spectator charged scalars on a closed wall."""
import argparse
import gzip
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import brentq


class History:
    def __init__(self, trace, stride=1):
        with gzip.open(trace, "rt") as file:
            rows = np.genfromtxt(file, delimiter=",", names=True)
        idx = np.unique(np.r_[np.arange(0, len(rows), stride), len(rows)-1])
        rows = rows[idx]
        t, a = rows["t"], rows["a"]
        assert t[0] == 0 and np.all(np.diff(t) > 0) and np.all(a > 0)
        eta = cumulative_simpson(1/a, x=t, initial=0.)
        pot = .5*a**2*rows["Bmag"]
        dpot = a**2*rows["X"]
        self.potential = CubicHermiteSpline(eta, pot, dpot)
        self.radius = PchipInterpolator(eta, a)
        self.time = PchipInterpolator(eta, t)
        hubble = PchipInterpolator(eta, rows["H"])
        self.qfrac = PchipInterpolator(eta, np.divide(rows["sheet"], rows["e"],
                                      out=np.zeros_like(a), where=rows["e"]>0))
        def crossing(values, func):
            falling = np.flatnonzero((values[:-1]>0)&(values[1:]<=0))
            i = falling[-1]
            return brentq(func, eta[i], eta[i+1], xtol=1e-14)
        self.cross = crossing(pot, self.potential)
        self.turn = crossing(rows["H"], hubble)
        self.early = float(np.interp(.08,t,eta))
        self.partition = np.interp([1.,2.], t, eta)
        assert 0 < self.early < self.cross < self.turn


class Protocol:
    def __init__(self, history, endpoint, tail=0.):
        self.history = history
        self.endpoint = endpoint
        self.tail = tail
        self.end = endpoint+tail
        self.val = float(history.potential(endpoint))
        self.slope = float(history.potential(endpoint,1))

    def value(self, eta, derivative=False):
        eta = np.asarray(eta)
        v = self.history.potential(np.minimum(eta,self.endpoint), int(derivative))
        if self.tail:
            u=np.clip((eta-self.endpoint)/self.tail,0,1)
            extra = (.5*self.slope*(1+np.cos(np.pi*u)) if derivative else
                     self.val+.5*self.slope*self.tail*(u+np.sin(np.pi*u)/np.pi))
            v = np.where(eta>self.endpoint, extra, v)
        return v


def labels(nmax):
    j = np.concatenate([np.full(n+1,n+1,dtype=int) for n in range(nmax+1)])
    k = np.concatenate([np.arange(-n,n+1,2) for n in range(nmax+1)])
    return j,k


def propagate(protocol, charge, nmax, refinement=1, reverse=False):
    j,k=labels(nmax)
    c=charge/math.sqrt(3)
    def omega(amp):
        s=c*amp
        return np.sqrt(j*j+s*s-2*k*s)
    # The incoming feature occupies the first few thousandths of conformal time.
    # Compare midpoint propagations with two time meshes, each aligned to the switch.
    edges=[0.,protocol.history.early,protocol.endpoint]
    targets=[1.5e-5,1.2e-3]
    if protocol.tail:
        edges.append(protocol.end)
        targets.append(1.2e-3)
    mesh=np.concatenate([np.linspace(l,r,1+math.ceil((r-l)*refinement/h))[:-1]
                         for l,r,h in zip(edges[:-1],edges[1:],targets)]+[np.array([protocol.end])])
    mid=(mesh[:-1]+mesh[1:])/2
    amps=protocol.value(mid)
    dt=np.diff(mesh)
    finalamp=protocol.value(protocol.end)
    if reverse: # Prepare the reciprocal experiment's vacuum at the other endpoint.
        amps,dt=amps[::-1],dt[::-1]
    wi=omega(finalamp if reverse else 0.)
    wf=omega(0. if reverse else finalamp)
    q=1/np.sqrt(2*wi)+0j
    p=-1j*np.sqrt(wi/2)
    for a,h in zip(amps,dt):
        w=omega(a)
        cs,sn=np.cos(w*h),np.sin(w*h)
        q,p=cs*q+(sn/w)*p,cs*p-(sn*w)*q
    alpha=(np.sqrt(wf)*q+1j*p/np.sqrt(wf))/math.sqrt(2)
    beta=(np.sqrt(wf)*q-1j*p/np.sqrt(wf))/math.sqrt(2)
    occ=abs(beta)**2
    return j,k,wf,occ,float(np.max(abs(abs(alpha)**2-occ-1)))


def fourier(protocol, top, exponent, early=False):
    # Integrate A' at all integer frequencies 2j. Uniform midpoint FFT folded
    # in intervals of pi. Subtract the endpoint jump using an analytic ramp.
    M=2**exponent
    if top >= M//2: raise ValueError("frequency aliases")
    h=np.pi/M
    nperiod=math.ceil(protocol.end/np.pi)
    eta=(np.arange(M*nperiod)+.5)*h
    inside=eta<protocol.end
    u=eta[inside]/protocol.end
    D=float(protocol.value(protocol.end,True))
    smooth=np.zeros_like(eta)
    dfield=protocol.value(eta[inside],True)
    if early:
        left,right=protocol.history.partition
        v=np.clip((eta[inside]-left)/(right-left),0.,1.)
        dfield=dfield*(1-v**3*(10-15*v+6*v*v))
        D=0.
    smooth[inside]=dfield-D*u**2*(3-2*u)
    j=np.arange(1,top+1)
    F=h*np.exp(-1j*j*h)*np.fft.rfft(smooth.reshape(nperiod,M).sum(axis=0))[1:top+1]
    z=-2j*j*protocol.end
    e=np.exp(z)
    integrals=[np.expm1(z)/z]
    for order in range(1,4):
        integrals.append((e-order*integrals[-1])/z)
    F+=D*protocol.end*(3*integrals[2]-2*integrals[3])
    return F


def smooth_tail(protocol,F):
    # Exact integrals for the stored cubic potential plus its analytic readout.
    # D=A'. Integration by parts twice gives |F_j| <= TV(D')/(4j²).
    # The Parseval remainder uses numerical low Fourier coefficients.
    spline=protocol.history.potential
    n=np.searchsorted(spline.x,protocol.endpoint,side="left")
    widths=np.diff(np.r_[spline.x[:n],protocol.endpoint])
    left=2*spline.c[1,:n]
    right=left+6*spline.c[0,:n]*widths
    tv=(np.sum(abs(right-left))+np.sum(abs(left[1:]-right[:-1]))+abs(left[0])+
        abs(right[-1])+np.pi*abs(protocol.slope)/protocol.tail)
    norm=np.sum(widths*(left**2+left*right+right**2)/3)
    norm+=protocol.slope**2*np.pi**2/(8*protocol.tail)
    j=np.arange(1,len(F)+1)
    assert protocol.end<np.pi
    residual=.5*np.pi*norm-np.sum(4*j*j*abs(F)**2)
    if residual<=0: raise ArithmeticError("Resolved spectrum exceeds Parseval norm")
    return dict(cutoff_j=len(F),electric_derivative_total_variation=float(tv),
                energy_coefficient_tail_bound_TV=float(tv**2/(576*len(F)**2)),
                energy_coefficient_tail_bound_numerical_Parseval=float(residual/(72*(len(F)+1))))


def shell_born(protocol,top,exponent):
    # Number and excitation energy of particles plus antiparticles, coefficients
    # of q^2. Sum over k and its independent j-fold degeneracy, with c=q/sqrt(3).
    j=np.arange(1,top+1)
    F=fourier(protocol,top,exponent)
    N=(1-j.astype(float)**-2)*abs(F)**2/18
    E=j*N
    cutoffs=[n for n in (8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768) if n<=top]
    return {str(n):dict(pair_quanta=float(np.sum(N[:n])),excitation_energy=float(np.sum(E[:n])))
            for n in cutoffs}, F


def independent_mode(protocol,n,k):
    # Integrate the original complex oscillator, independent of the transfer matrix.
    def rhs(t,y):
        s=protocol.value(t)/math.sqrt(3)
        return [y[1],-((n+1)**2+s*s-2*k*s)*y[0]]
    y=np.array([1/math.sqrt(n+1),-1j*math.sqrt(n+1)])/math.sqrt(2.)
    for l,r,h in [(0.,protocol.history.early,2e-5),(protocol.history.early,protocol.end,.005)]:
        sol=solve_ivp(rhs,(l,r),y,rtol=2e-11,atol=1e-14,method="DOP853",max_step=h)
        if not sol.success: raise ArithmeticError(sol.message)
        y=sol.y[:,-1]
    s=float(protocol.value(protocol.end))/math.sqrt(3)
    w=math.sqrt((n+1)**2+s*s-2*k*s)
    return float(abs((np.sqrt(w)*y[0]-1j*y[1]/np.sqrt(w))/math.sqrt(2))**2)


def extract(r, nmax=8):
    j,k,w,occ,err=r
    # One particle sign; for an individual charge-conjugate mode k changes sign.
    logp=-float(np.sum(j*np.log1p(occ)))
    mean=float(np.sum(j*occ))
    probability=-math.expm1(logp)
    if not 0 <= probability <= min(1.,mean)+1e-13:
        raise ArithmeticError("Vacuum persistence normalization")
    return dict(max_wronskian_error=err, total_quanta_finite_band=2*mean,
                vacuum_log_probability=logp, probability_of_any_pair_finite_band=probability,
                conformal_excitation_energy_finite_band=float(2*np.sum(j*w*occ)),
                low_modes=[dict(n=int(a-1),k=int(b),omega=float(c),occupation=float(d))
                           for a,b,c,d in zip(j,k,w,occ) if a<=nmax+1])


def calculate(trace):
    hist=History(trace)
    coarse=History(trace,2)
    protocols={"zero_crossing":Protocol(hist,hist.cross),"turnaround":Protocol(hist,hist.turn),
               "taper_short":Protocol(hist,hist.turn,.125),"taper_long":Protocol(hist,hist.turn,.25)}
    out={"parameters":dict(nmax_exact=32,charges=[.1,1.,10.],gauge_component="A/sqrt(3)",
                           dielectric="fixed external background",
                           input_model="forward.maxwell_primary, proper depth 0.002, late gauge 13200"),
         "events":{s:dict(eta=z,t=float(hist.time(z)),radius=float(hist.radius(z)),
                          potential=float(hist.potential(z)),dA_deta=float(hist.potential(z,1)),
                          em_fraction=float(hist.qfrac(z))) for s,z in (("zero_crossing",hist.cross),("turnaround",hist.turn))},
         "protocols":{}}
    for name,prot in protocols.items():
        dataset={}
        base=propagate(prot,1.,32,2)
        dataset["q1"]=extract(base)
        if name=="zero_crossing":
            half=propagate(prot,1.,32,1)
            dec=propagate(Protocol(coarse,coarse.cross),1.,32,2)
            reverse=propagate(prot,1.,32,2,True)
            neutral=propagate(prot,0.,4,1)
            neg=propagate(prot,-1.,8)
            signed=propagate(prot,1.,8)
            flips=np.concatenate([np.where(neg[0]==j)[0][::-1] for j in range(1,10)])
            dataset["controls"]=dict(max_timestep_occupation_difference=float(np.max(abs(base[3]-half[3]))),
                max_decimated_input_occupation_difference=float(np.max(abs(base[3]-dec[3]))),
                max_forward_reverse_occupation_difference=float(np.max(abs(base[3]-reverse[3]))),
                max_charge_conjugation_difference=float(np.max(abs(signed[3]-neg[3][flips]))),
                max_free_occupation=float(max(neutral[3])),
                endpoint_time_decimation_difference=float(abs(hist.time(hist.cross)-coarse.time(coarse.cross))))
            for q in [.1,10.]:
                dataset['q'+str(q)]=extract(propagate(prot,q,32,2))
        born,F=shell_born(prot,32768,20)
        j=np.arange(1,len(F)+1)
        efactor=(j-1./j)/18
        efull=efactor*abs(F)**2
        dataset["born_energy_characteristics"]=dict(peak_j=int(1+np.argmax(efull)),
            median_j=int(1+np.searchsorted(np.cumsum(efull),np.sum(efull)/2)))
        if prot.tail:
            dataset["smooth_ultraviolet_tail"]=smooth_tail(prot,F)
        if name=="taper_short":
            Fearly=fourier(prot,32768,20,early=True)
            Flate=F-Fearly
            partial=lambda x:float(np.sum(efactor*abs(x)**2))
            inter=float(2*np.sum(efactor*(Fearly*Flate.conj()).real))
            dataset["born_energy_temporal_partition"]=dict(early=partial(Fearly),late=partial(Flate),
                interference=inter,closure_error=partial(F)-partial(Fearly)-partial(Flate)-inter,
                partition_wall_time=[1.,2.])
            ch=Protocol(coarse,coarse.turn,.125)
            Fc=fourier(ch,32768,20)
            dataset["fourier_input_decimation_relative_energy_change"]=partial(Fc)/partial(F)-1
        low,F2=shell_born(prot,32768,19)
        dataset["born_coefficients_per_squared_charge"]=born
        dataset["fourier_grid_check_at_32768"]=float(abs(born["32768"]["excitation_energy"]/
                                                      low["32768"]["excitation_energy"]-1))
        # The j=1,...,33 low-band Born F_j can also check the exact small-charge run.
        if name=="zero_crossing":
            j,k=base[:2]
            small=propagate(prot,.01,8,2)
            J,K,_,occ,_=small
            dataset["controls"]["homogeneous_direct_occupation"]=independent_mode(prot,0,0)
            dataset["controls"]["homogeneous_transfer_occupation"]=float(base[3][0])
            deg=base[3][base[0]==3]
            direct=np.array([independent_mode(prot,2,k) for k in (-2,0,2)])
            dataset["controls"]["equal_energy_shell2_independent_occupations"]=direct.tolist()
            dataset["controls"]["equal_energy_shell2_relative_direct_error"]=float(max(abs(deg/direct-1)))
            dataset["equal_energy_shell2_minmax_ratio"]=float(max(deg)/min(deg))
            dataset["controls"]["equal_energy_shell2_frequency_spread"]=float(np.ptp(base[2][base[0]==3]))
            pred=(.01/math.sqrt(3)*K/(2*J**2))**2*abs(F[J-1])**2
            dataset["born_relative_summed_number_error_at_q001"]=float(
                 abs(np.sum(J*occ)/np.sum(J*pred)-1))
        dataset["instantaneous_energy_log_slope_per_q_squared"]=float(prot.value(prot.end,True)**2/72)
        out["protocols"][name]=dataset
    bare=out["protocols"]["turnaround"]
    smooth=out["protocols"]["taper_short"]
    def octave(d):
        b=d["born_coefficients_per_squared_charge"]
        return b["32768"]["excitation_energy"]-b["16384"]["excitation_energy"]
    observed=octave(bare)-octave(smooth)
    expected=math.log(2)*bare["instantaneous_energy_log_slope_per_q_squared"]
    out["raw_cutoff_check"]=dict(last_octave_bare_minus_smooth=observed,predicted_log_increment=expected)
    principal=out["protocols"]["zero_crossing"]
    controls=principal["controls"]
    if (principal["equal_energy_shell2_minmax_ratio"] < 40000 or
        controls["equal_energy_shell2_relative_direct_error"] > 1e-5 or
        controls["max_free_occupation"]>1e-20 or
        controls["max_decimated_input_occupation_difference"]>1e-12 or
        controls["max_forward_reverse_occupation_difference"]>1e-12 or
        principal["born_relative_summed_number_error_at_q001"]>1e-5 or
        abs(observed/expected-1)>.02 or
        abs(smooth["fourier_input_decimation_relative_energy_change"])>1e-4 or
        smooth["smooth_ultraviolet_tail"]["energy_coefficient_tail_bound_TV"]>.3):
        raise ArithmeticError("Quantum spectrum or numerical controls changed")
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--trace",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(args.trace),indent=2,sort_keys=True,allow_nan=False)+"\n")
