"""Ordered and symmetrized spectra of a free parent on the four retained fine walls."""
import argparse
import io
import json
import zipfile
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from scipy.special import roots_legendre


CLOCK = (1.35, 2.05)
READOUT = 2.10
BAND = (60., 120.)
RESEARCH = Path("research")
FORK = RESEARCH/"HDBLAST_COVARIANCE_FEEDBACK_20261119/native"
SOURCES = RESEARCH/"HDBLAST_SCALAR_CLOCK_GRAVITY_20261119/source/SOURCES.json"
DRIVE = RESEARCH/"HDBLAST_DISPERSION_MATCHING_20261119/driven_check/ENERGY.json"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def wall(native, arm, seed):
    with zipfile.ZipFile(native/"WALLS.zip") as archive:
        raw=archive.read(f"fine/retain/{arm}/wall.csv")
    data=np.genfromtxt(io.BytesIO(raw), delimiter=",", names=True)
    s=data["H0tau"]
    loga=CubicHermiteSpline(s,data["A_b"]-seed["reference_A"],data["H_over_H0"],
                            extrapolate=False)
    phi=CubicHermiteSpline(s,data["phi_b"],data["v_over_H0"],extrapolate=False)
    return loga,phi


def fundamental(seed,loga,phi,quiet):
    q=np.asarray(seed["q"])
    start=seed["H0tau"]
    ends=[brentq(lambda s:float(phi(s)-quiet(t)),start,READOUT) for t in CLOCK]
    def rhs(s,state):
        h,k,ph,pk=state[:-1].reshape(4,-1)
        a=np.exp(loga(s))
        m=100*(phi(s)-.5)
        stiffness=a*(q**2+(a*m)**2)
        return np.r_[ph/a**3,pk/a**3,-stiffness*h,-stiffness*k,1/a]
    ones=np.ones_like(q)
    zeros=np.zeros_like(q)
    initial=np.r_[ones,zeros,zeros,ones,0.]
    solution=solve_ivp(rhs,(start,max(ends)),initial,dense_output=True,method="DOP853",
                       max_step=.0015,rtol=2e-12,atol=3e-14)
    require(solution.success,solution.message)
    h,k,ph,pk=solution.y[:-1].reshape(4,len(q),-1)
    det_error=float(np.max(abs(h*pk-k*ph-1)))
    require(det_error<3e-11,"Free canonical transfer determinant")
    return ends,solution.sol,det_error


def profile(seed,ends,sol,phi,quiet,order):
    # Pulse walls contain narrow variations. Split at their measured native knots.
    inside=phi.x[(phi.x>ends[0]) & (phi.x<ends[1])]
    bounds=np.r_[ends[0],inside,ends[1]]
    nodes,weights=roots_legendre(order)
    width=np.diff(bounds)[:,None]/2
    s=(bounds[:-1,None]+width*(1+nodes)).ravel()
    ds=(width*weights).ravel()
    # Vectorized inversion evaluates the scalar gate on the pulse trajectory.
    field=phi(s)
    lo,hi=np.full_like(s,CLOCK[0]),np.full_like(s,CLOCK[1])
    for _ in range(48):
        mid=(lo+hi)/2
        lower=quiet(mid)<field
        lo,hi=np.where(lower,mid,lo),np.where(lower,hi,mid)
    quiet_time=(lo+hi)/2
    gate=np.sin(np.pi*(quiet_time-CLOCK[0])/(CLOCK[1]-CLOCK[0]))**4
    fields=sol(s)
    h,k=fields[:-1].reshape(4,len(seed["q"]),len(s))[:2]
    return fields[-1],gate*ds,h,k


def amplitudes(seed,profiles,low,high,nw):
    q,wq=np.asarray(seed["q"]),np.asarray(seed["weight"])
    eta,measure_s,h,k=profiles
    nodes,weights=roots_legendre(nw)
    lower=np.broadcast_to(low,q.shape)
    widths=high-lower
    omega=lower[:,None]+widths[:,None]*(nodes+1)/2
    d_omega=widths[:,None]*weights/2
    spec=omega*(omega**2-q[:,None]**2)*d_omega*wq[:,None]/(8*np.pi**2)
    # Direct integrals avoid subtracting large full-cutoff totals.
    vh=np.empty(omega.shape,dtype=complex)
    vk=np.empty_like(vh)
    if np.ptp(lower)==0:
        osc=np.exp(1j*omega[0,:,None]*eta)*measure_s
        vh[:]=(osc@h.T).T
        vk[:]=(osc@k.T).T
    else:
        for j,frequencies in enumerate(omega):
            osc=np.exp(1j*frequencies[:,None]*eta)*measure_s
            vh[j]=osc@h[j]
            vk[j]=osc@k[j]
    return spec,vh,vk


def measure(seed,profiles,inertia,low,high,nw):
    spec,vh,vk=amplitudes(seed,profiles,low,high,nw)
    def total(field):
        return float(np.sum(spec*field))
    i=inertia[:,None]
    # chi(t)=f(t)a+f*(t)a† with f=(h-i I k)/sqrt(2I).
    # A and B refer to f transformed at opposite signs of the bath frequency.
    A=(vh-1j*i*vk)/np.sqrt(2*i)
    B=(vh.conj()-1j*i*vk.conj())/np.sqrt(2*i)
    aa,bb=abs(A)**2,abs(B)**2
    commutator=-np.imag(vh.conj()*vk)
    np.testing.assert_allclose(commutator,(bb-aa)/2,rtol=2e-12,atol=5e-16)
    states={}
    for state in ("retain","erase","fork_vacuum"):
        if state=="fork_vacuum":
            X,C,P=1/(2*i),np.zeros_like(i),i/2
        else:
            X,C,P=(np.asarray(seed["states"][state][key])[:,None] for key in ("X","C","P"))
        require(np.min(X*P-C**2)>.25-1e-9,"Unphysical fork covariance")
        sym=X*abs(vh)**2+P*abs(vk)**2+2*C*np.real(vh.conj()*vk)
        ordered=sym+commutator
        reverse=sym-commutator
        number=(i*X+P/i-1)/2
        anomalous=(i*X-P/i)/2+1j*C
        annihilation=number*aa
        creation=(number+1)*bb
        coherence=2*np.real(anomalous*A*B)
        require(np.min(ordered)>-1e-15 and np.min(reverse)>-1e-15,"Negative spectral weight")
        np.testing.assert_allclose(ordered,annihilation+creation+coherence,rtol=2e-12,atol=2e-16)
        np.testing.assert_allclose(reverse,(number+1)*aa+number*bb+coherence,rtol=2e-12,atol=2e-16)
        states[state]=dict(symmetrized=total(sym),commutator=total(commutator),
                           emission=total(ordered),reverse_ordered=total(reverse),
                           annihilation_branch=total(annihilation),creation_branch=total(creation),
                           interference=total(coherence))
    return states


def stationary():
    # Stationary resonance with smooth compact gate and no counter-rotating Fourier weight.
    E=4.
    T=2*np.pi
    nodes,weights=roots_legendre(128)
    t=(nodes+1)*T/2
    gate=np.sin(np.pi*t/T)**4
    mode=np.exp(-1j*E*t)/np.sqrt(2*E)
    def transform(sign):
        return np.sum(np.exp(sign*1j*E*t)*mode*weights*gate)*T/2
    aa,bb=abs(transform(+1))**2,abs(transform(-1))**2
    expected=(3*T/8)**2/(2*E)
    require(abs(aa/expected-1)<2e-13 and bb<1e-26,"Stationary resonant frequency branches")
    tests={}
    for n in (0,2):
        emitted=n*aa+(n+1)*bb
        symmetric=(n+.5)*(aa+bb)
        comm=(bb-aa)/2
        require(abs(symmetric+comm-emitted)<1e-14,"Stationary commutator cancellation")
        tests[str(n)]=dict(emission=float(emitted),symmetrized=float(symmetric),
                           commutator=float(comm),exact_emission=float(n*expected))
    return dict(parent_and_bath_frequency=E,gate_duration=T,positive_branch_power=float(aa),
                negative_branch_power=float(bb),diagonal_occupations=tests)


def experiment(repo):
    native=repo/FORK
    seed=json.loads((native/"FORK.json").read_text())["fine"]
    archived=json.loads((repo/SOURCES).read_text())["runs"]
    classical=json.loads((repo/DRIVE).read_text())["positive_quadrature_work"]
    quiet=wall(native,"quiet",seed)[1]
    q=np.array(seed["q"])
    arms={}
    for arm in ("quiet","plus","minus","far"):
        loga,phi=wall(native,arm,seed)
        a0=np.exp(loga(seed["H0tau"]))
        inertia=a0**3*np.hypot(q/a0,100*(phi(seed["H0tau"])-.5))
        af=float(np.exp(loga(READOUT)))
        ends,sol,det_error=fundamental(seed,loga,phi,quiet)
        high_profile=profile(seed,ends,sol,phi,quiet,5)
        band=measure(seed,high_profile,inertia,*BAND,96)
        coarse=measure(seed,profile(seed,ends,sol,phi,quiet,3),inertia,*BAND,48)
        convergence=max(abs(band[s][key]-coarse[s][key]) for s in band for key in band[s])
        require(convergence<5e-7,f"Quadrature refinement {arm}")
        cutoff_values={}
        for cut in BAND:
            reconstructed=measure(seed,high_profile,inertia,q,cut,64)["retain"]["emission"]/af
            data=archived[f"fine_{arm}_{int(cut)}"]["files"]
            reference=data["static"]["readout"]["daughter_comoving"]
            require(reference==data["bare"]["readout"]["daughter_comoving"],"Counterterm dependence")
            cutoff_values[str(int(cut))]=dict(rebuilt= float(reconstructed),archived=reference,
                                            difference=float(reconstructed-reference))
            require(abs(reconstructed-reference)<5e-6,f"Archived cutoff disagreement {arm}/{cut}")
        proper={key:value/af for key,value in band["retain"].items()}
        difference=cutoff_values["120"]["archived"]-cutoff_values["60"]["archived"]
        rebuilt_difference=cutoff_values["120"]["rebuilt"]-cutoff_values["60"]["rebuilt"]
        require(abs(difference-proper["emission"])<5e-6,"Band and archive difference")
        require(abs(rebuilt_difference-proper["emission"])<1e-7,"Independent band quadrature additivity")
        states_equal_comm=max(abs(r["commutator"]-band["retain"]["commutator"]) for r in band.values())
        require(states_equal_comm==0,"Commutator changed with preparation")
        if arm=="quiet":
            require(abs(band["retain"]["symmetrized"]-classical)<2e-7,"Prior quiet classical bath work")
        arms[arm]=dict(gate_proper=ends,scale_factor_at_readout=af,determinant_error=det_error,
                       all_preparations_conformal=band,retained_proper=proper,
                       max_component_change_under_quadrature_refinement=convergence,
                       positive_frequency_cutoff_checks=cutoff_values,archived_cutoff_difference=difference,
                       direct_band_minus_archive_difference=proper["emission"]-difference,
                       rebuilt_cutoff_difference_minus_direct_band=rebuilt_difference-proper["emission"])
    determinants={}
    for name in ("retain","erase"):
        X,C,P=(np.asarray(seed["states"][name][key]) for key in ("X","C","P"))
        gap=X*P-C**2-.25
        determinants[name]=dict(min_gap_above_quarter=float(min(gap)),max_gap_above_quarter=float(max(gap)))
    return dict(grid="fine",band=list(BAND),proper_readout=READOUT,
                covariance_determinant=determinants,
                archived_quiet_conformal_drive_work=classical,arms=arms,stationary_probe=stationary())


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--repo",type=Path,required=True)
    cli.add_argument("--output",type=Path,required=True)
    args=cli.parse_args()
    args.output.write_text(json.dumps(experiment(args.repo),indent=2,allow_nan=False)+"\n")
