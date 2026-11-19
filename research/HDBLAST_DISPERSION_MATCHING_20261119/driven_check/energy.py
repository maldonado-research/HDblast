"""Fourier check of the high-band retarded work on two random free-parent quadratures."""
import argparse
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.special import roots_legendre


def require(condition, message):
    if not condition:
        raise ValueError(message)


OPEN, CLOSE = 1.35, 2.05
LOW, HIGH = 60., 120.


def parent(native):
    seed=json.loads((native/"FORK.json").read_text())["fine"]
    q,weights=np.asarray(seed["q"]),np.asarray(seed["weight"])
    pop=seed["states"]["retain"]
    X,C,P=(np.asarray(pop[k]) for k in ("X","C","P"))
    require(np.all((X>0)&(X*P-C**2>0)), "Initial covariance must be positive definite")
    with zipfile.ZipFile(native/"WALLS.zip") as archive:
        w=np.genfromtxt(io.BytesIO(archive.read("fine/retain/quiet/wall.csv")),delimiter=",",names=True)
    t=w["H0tau"]
    loga=CubicHermiteSpline(t,w["A_b"]-seed["reference_A"],w["H_over_H0"],extrapolate=False)
    mass=CubicHermiteSpline(t,100*(w["phi_b"]-.5),100*w["v_over_H0"],extrapolate=False)
    require(np.all(mass(np.linspace(OPEN,CLOSE,2001),1)>0), "Quiet wall fails scalar-clock inversion")
    q2=q*q
    ones=np.ones_like(q)
    zero=np.zeros_like(q)
    # Fundamental solutions start from (chi,Pi)=(1,0) and (0,1).
    def rhs(s,state):
        h,k,ph,pk=state[:-1].reshape(4,-1)
        a=np.exp(loga(s))
        stiffness=a*(q2+(a*mass(s))**2)
        return np.r_[ph/a**3,pk/a**3,-stiffness*h,-stiffness*k,1/a]
    initial=np.r_[ones,zero,zero,ones,0.]
    sol=solve_ivp(rhs,(t[0],CLOSE),initial,dense_output=True,method="DOP853",
                  rtol=2e-12,atol=3e-14,max_step=.0015)
    require(sol.success,sol.message)
    h,k,ph,pk=sol.y[:-1].reshape(4,len(q),-1)
    det_error=float(np.max(abs(h*pk-k*ph-1)))
    require(det_error<3e-10, "Free-parent transfer lost canonical determinant")
    return q,weights,X,C,P,sol.sol,loga,det_error


def energy(inputs,time_order,frequency_order):
    q,wq,X,C,P,sol,_,_=inputs
    time_nodes,time_weights=roots_legendre(time_order)
    s=OPEN+(CLOSE-OPEN)*(1+time_nodes)/2
    ds=time_weights*(CLOSE-OPEN)/2
    state=sol(s)
    h,k=state[:-1].reshape(4,len(q),len(s))[:2]
    eta=state[-1]
    phase=np.pi*(s-OPEN)/(CLOSE-OPEN)
    gate=np.sin(phase)**4
    nodes,weights=roots_legendre(frequency_order)
    omega=LOW+(HIGH-LOW)*(nodes+1)/2
    dw=weights*(HIGH-LOW)/2
    fourier=np.exp(1j*omega[:,None]*eta[None,:])*(gate*ds)[None,:]
    # d eta a F chi = ds F chi, including the conformal source conversion.
    vh= fourier@h.T
    vk= fourier@k.T
    power=X*abs(vh)**2+P*abs(vk)**2+2*C*(vh.conj()*vk).real
    require(np.all(power>=0), "Negative averaged drive power")
    g2=dw[:,None]*omega[:,None]*(omega[:,None]**2-q[None,:]**2)/(4*np.pi**2)
    per_q=.5*wq*np.sum(g2*power,axis=0)
    return float(np.sum(per_q)),per_q


def force_work(inputs, order):
    q,wq,X,C,P,sol,loga,_=inputs
    nodes,weights=roots_legendre(order)
    s=OPEN+(CLOSE-OPEN)*(nodes+1)/2
    ds=weights*(CLOSE-OPEN)/2
    a=np.exp(loga(s))
    d_eta=ds/a
    h0,k0,ph0,pk0=sol(s)[:-1].reshape(4,len(q),order)
    def gate(t):
        phase=np.pi*(t-OPEN)/(CLOSE-OPEN)
        f=np.sin(phase)**4
        ds_f=4*np.pi/(CLOSE-OPEN)*np.sin(phase)**3*np.cos(phase)
        return f,ds_f
    F,ds_F=gate(s)
    common=a*a*(loga(s,1)*F+ds_F)
    uh=h0*common + ph0*F/a
    uk=k0*common + pk0*F/a
    # A second integration runs only over earlier proper times at each sample.
    past=OPEN+(s[:,None]-OPEN)*(1+nodes[None,:])/2
    width=(s[:,None]-OPEN)*weights[None,:]/2
    history=sol(past.ravel())
    h,k=history[:-1].reshape(4,len(q),order,order)[:2]
    eta_now=sol(s)[-1]
    lags=eta_now[:,None]-history[-1].reshape(order,order)
    omega=LOW+(HIGH-LOW)*(nodes+1)/2
    dw=weights*(HIGH-LOW)/2
    sins=np.sin(lags[:,:,None]*omega[None,None,:])
    k0=sins@dw
    k2=sins@(dw*omega**2)
    kernel=(k2[None,:,:]-q[:,None,None]**2*k0[None,:,:])/(4*np.pi**2)
    measure=kernel*width*gate(past)[0]
    rh=np.sum(h*measure,axis=2)
    rk=np.sum(k*measure,axis=2)
    def product(vh,vk,zh,zk):
        mixed=X[:,None]*vh*zh+P[:,None]*vk*zk+C[:,None]*(vh*zk+vk*zh)
        return float(wq @ mixed @ d_eta)
    work=-product(uh,uk,rh,rk)
    norm_force=math.sqrt(product(rh,rk,rh,rk))
    norm_velocity=math.sqrt(product(uh,uk,uh,uk))
    return dict(negative_integral_uprime_force=work, force_norm=norm_force,
                drive_velocity_norm=norm_velocity,
                reversible_relative_rms_error_lower_bound=work/(norm_force*norm_velocity))


def experiment(native, driven=None):
    inp=parent(native)
    values={}
    for nt,nw in ((64,24),(96,48),(160,96)):
        e,per_q=energy(inp,nt,nw)
        values[f"{nt}x{nw}"]=e
    best=values["160x96"]
    discrepancy=max(abs(e/best-1) for e in values.values())
    require(discrepancy<1e-7, "Fourier quadrature resolution")
    work={str(n):force_work(inp,n) for n in (48,80)}
    for row in work.values():
        require(abs(row["negative_integral_uprime_force"]/best-1)<1e-7, "Frequency energy and retarded work disagree")
        require(0<row["reversible_relative_rms_error_lower_bound"]<1, "Invalid relative error bound")
    result=dict(retarded_work=work,arm="fine retain quiet",band_comoving=[LOW,HIGH],proper_gate=[OPEN,CLOSE],
                units="conformal bath energy per unit comoving volume in the forcing normalization",
                positive_quadrature_work=best,convergence=values,
                max_relative_resolution_change=discrepancy,
                maximum_parent_determinant_error=inp[-1],
                weighted_contribution_by_native_mode=per_q.tolist())
    if driven is not None:
        physical=json.loads(driven.read_text())
        result["time_domain_comparison"]=compare_driven(physical,best,work["80"])
    return result


def compare_driven(report,energy,fourier_work):
    require(report["grid"]=="fine" and report["arm"]=="quiet" and report["band"]==[LOW,HIGH],
            "Direct solver background or bath band differs")
    row=report["runs"]["primary"]
    energy_checks={key:abs(row["conformal_bath_energy"][key]/energy-1)
                   for key in ("at_close","at_readout","integrated_drive_work")}
    require(max(energy_checks.values())<2e-7, "Fourier energy disagrees with forced oscillators")
    norm_error=abs(row["open_exact_norm"]/fourier_work["force_norm"]-1)
    require(norm_error<2e-7, "Time-domain and independent force norms differ")
    require(all(row["rms_force_errors"][key] >= fourier_work["reversible_relative_rms_error_lower_bound"]
                for key in ("static","static_plus_kinetic")), "Zero-work lower bound")
    return dict(relative_energy_errors=energy_checks,relative_exact_force_norm_error=norm_error)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native",type=Path,required=True)
    ap.add_argument("--driven",type=Path,default=Path(__file__).resolve().parents[1]/"driven/DRIVEN.json",
                    help="Direct oscillator-work report to compare")
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    out=experiment(args.native,args.driven)
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
