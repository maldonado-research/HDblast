"""Drive the finite 60–120 bath with the 24 native retained parent quadratures."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
REV="2cbafaefde87652135d8222e50a53da4eb704c41"
NATIVE=Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
SOURCE=Path("research/HDBLAST_SCALAR_CLOCK_GRAVITY_20261119/source/source.py")
PINS={NATIVE/"FORK.json":"d1028e10aa6440c98d71e5c99cf2c18fd336606d3068816f02d01db499918038",
      NATIVE/"WALLS.zip":"a18041b6ab8591cc88e7613e0759c812d00da5655953b5cb3a25fbb36bc9d05d",
      SOURCE:"315f602a242c069d06d641b70a00db651d50bd4164436828503f078427a06aad"}
OPEN,CLOSE,END=1.35,2.05,2.10
LO,HI=60.,120.


def background(repo):
    for path,digest in PINS.items():
        if hashlib.sha256((repo/path).read_bytes()).hexdigest()!=digest:
            raise ValueError(path)
    spec=importlib.util.spec_from_file_location("historical_source",repo/SOURCE)
    src=importlib.util.module_from_spec(spec);spec.loader.exec_module(src)
    q,weights,states,(start,loga,mass)=src.read_native(repo/NATIVE,"fine","quiet")
    X,C,P=states["retain"]
    chi=np.stack((np.sqrt(X),np.zeros_like(X)),axis=1)
    pi=np.stack((C/np.sqrt(X),np.sqrt(P-C*C/X)),axis=1)
    np.testing.assert_allclose([np.sum(chi**2,axis=1),np.sum(chi*pi,axis=1),np.sum(pi**2,axis=1)],
                               [X,C,P],rtol=3e-15,atol=1e-14)
    initial=np.array((chi,pi)).ravel()
    def rhs(s,state):
        a=np.exp(loga(s)); m=mass(s)
        x,p=state.reshape(2,len(q),2)
        return np.stack((p/a**3,-(a*q*q+a**3*m*m)[:,None]*x)).ravel()
    free=solve_ivp(rhs,(start,END),initial,method="DOP853",rtol=2e-12,atol=8e-15,
                   dense_output=True,max_step=.002)
    if not free.success: raise RuntimeError(free.message)
    return q,weights,loga,mass,src.ScalarClock(mass),free.sol


def drives(history, t):
    q,weights,loga,mass,clock,free=history
    x,p=free(t).reshape(2,len(q),2,len(t)).transpose(0,3,1,2)
    a=np.exp(loga(t));h=loga(t,1);hs=loga(t,2)
    m=mass(t);v=mass(t,1)/100
    f,fp,fpp=clock(m/100+.5)
    fs=fp*v;fss=fpp*v*v+fp*mass(t,2)/100
    # eta derivatives use d/deta=a*d/ds and the free parent equation.
    u=(a*f)[:,None,None]*x
    coeff=(hs+2*h*h)*f+3*h*fs+fss-m*m*f
    usecond=(a**3)[:,None,None]*((coeff[:,None]-f[:,None]*(q[None,:]/a[:,None])**2)[:,:,None]*x)
    usecond+=2*fs[:,None,None]*p
    return a,u,usecond


def run(history, steps, order):
    q,w,loga,mass,_,_=history
    gx,gw=np.polynomial.legendre.leggauss(order)
    omega=LO+(HI-LO)*(gx+1)/2
    omega2=omega[:,None,None]**2
    g2=(gw*omega*(HI-LO)/(8*np.pi**2))[:,None]*(omega[:,None]**2-q[None,:]**2)
    c=np.sum(g2/omega[:,None]**2,axis=0)
    d=np.sum(g2/omega[:,None]**4,axis=0)
    np.testing.assert_allclose(c,((HI**2-LO**2)/2-q*q*np.log(HI/LO))/(4*np.pi**2),rtol=2e-14)
    np.testing.assert_allclose(d,(np.log(HI/LO)+q*q*(1/HI**2-1/LO**2)/2)/(4*np.pi**2),rtol=2e-14)
    # Use interior midpoint values of a'' and F'' to avoid one-sided derivative assignments at wall records.
    native=mass.x[(mass.x>OPEN)&(mass.x<END)]
    knots=np.unique(np.r_[OPEN,native,CLOSE,END])
    scale=(END-OPEN)/steps
    edges=np.r_[np.concatenate([np.linspace(l,r,int(np.ceil((r-l)/scale))+1)[:-1]
                                for l,r in zip(knots[:-1],knots[1:])]),END]
    times=np.r_[(edges[:-1,None]+np.diff(edges)[:,None]*np.array([0.,.25,.5,.75])).ravel(),END]
    a,u,uss=drives(history,times)
    Q=np.zeros((order,len(q),2));P=np.zeros_like(Q)
    def advance(i,ds,Q,P):
        r0=1/a[i];r1=1/a[i+1];r2=1/a[i+2]
        x1=P*r0; p1=(u[i]-omega2*Q)*r0
        x2=(P+ds/2*p1)*r1; p2=(u[i+1]-omega2*(Q+ds/2*x1))*r1
        x3=(P+ds/2*p2)*r1; p3=(u[i+1]-omega2*(Q+ds/2*x2))*r1
        x4=(P+ds*p3)*r2; p4=(u[i+2]-omega2*(Q+ds*x3))*r2
        return Q+ds/6*(x1+2*x2+2*x3+x4),P+ds/6*(p1+2*p2+2*p3+p4)
    def force(X):
        return np.einsum("fq,fql->ql",g2,X)
    def sq(x):
        return np.sum(x*x,axis=1)
    def power(i,p):
        return float(w@np.sum(u[i]*force(p),axis=1)/a[i])
    sums=np.zeros((4,len(q)))
    post=np.zeros(len(q));work=0.;p0=0.
    gate_energy=None
    for k,ds in enumerate(np.diff(edges)):
        i=4*k
        Q,P=advance(i,ds/2,Q,P)
        exact=force(Q); local=c[:,None]*u[i+2]; dynamic=local-d[:,None]*uss[i+2]
        factor=ds/a[i+2]
        if times[i+2]<CLOSE:
            sums+=factor*np.array([sq(x) for x in (exact,exact-local,exact-dynamic,dynamic-local)])
        else: post+=factor*sq(exact)
        pm=power(i+2,P)
        Q,P=advance(i+2,ds/2,Q,P)
        p1=power(i+4,P)
        work+=ds*(p0+4*pm+p1)/6
        p0=p1
        if edges[k+1]==CLOSE:
            gate_energy=float(.5*np.einsum("q,fq,fql->",w,g2,P*P+omega2*Q*Q))
    # Gauss-summed oscillator energies are conformal energies per squared drive, without vacuum energy.
    energy=float(.5*np.einsum("q,fq,fql->",w,g2,P*P+omega2*Q*Q))
    norm=np.sqrt(sums@w); endpoint=float(np.sqrt(w@sq(force(Q))))
    if gate_energy is None or not np.all(np.isfinite(norm)) or np.min(norm)<=0:
        raise ValueError("Missing gate or force signal")
    return dict(nominal_time_steps=steps,actual_cells=len(edges)-1,frequency_nodes=order,
        rms_force_errors=dict(static=float(norm[1]/norm[0]),static_plus_kinetic=float(norm[2]/norm[0]),
                              error_relative_to_nonstatic_reaction=float(norm[2]/norm[1])),
        open_exact_norm=float(norm[0]),open_static_error_norm=float(norm[1]),open_dynamic_error_norm=float(norm[2]),
        open_per_q_static_relative_rms=np.sqrt(sums[1]/sums[0]).tolist(),
        open_per_q_dynamic_relative_rms=np.sqrt(sums[2]/sums[0]).tolist(),
        postgate=dict(absolute_norm=float(np.sqrt(w@post)),relative_to_open=float(np.sqrt(w@post)/norm[0]),
                      endpoint_weighted_force_norm=endpoint),
        conformal_bath_energy=dict(at_close=gate_energy,at_readout=energy,integrated_drive_work=work,
                                   energy_minus_work=energy-work))


def calculate(repo):
    history=background(repo)
    variants={"coarse_time":(512,96),"fine_time":(1024,96),"coarse_frequency":(2048,48),"primary":(2048,96)}
    rows={name:run(history,*options) for name,options in variants.items()}
    final=rows["primary"]
    for test in ("fine_time","coarse_frequency"):
        for metrics in ("rms_force_errors","postgate"):
            for key in final[metrics]:
                np.testing.assert_allclose(rows[test][metrics][key],final[metrics][key],rtol=.002,atol=1e-10)
    np.testing.assert_allclose(final["conformal_bath_energy"]["at_readout"],
                               final["conformal_bath_energy"]["integrated_drive_work"],rtol=.003,atol=1e-11)
    return dict(source_revision=REV,grid="fine",arm="quiet",gate=[OPEN,CLOSE],readout=END,band=[LO,HI],
                q=history[0].tolist(),native_weight=history[1].tolist(),runs=rows)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo",type=Path,required=True)
    ap.add_argument("--output",type=Path,default=HERE/"DRIVEN.json")
    args=ap.parse_args()
    args.output.write_text(json.dumps(calculate(args.repo.resolve()),indent=2,allow_nan=False)+"\n")
