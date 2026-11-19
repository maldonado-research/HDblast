"""Finite-time scalar Yukawa switch on an explicitly frozen shell."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre


def load_source(repo):
    path = repo / "research/HDBLAST_QUANTUM_MEMORY_20261119/quantum/spectrum.py"
    spec = importlib.util.spec_from_file_location("ordered_switch_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def transform(nu, duration):
    """Integral from 0 to T of sin(pi t/T)^4 exp(i nu t) dt."""
    x = nu * duration
    h = x/(2*np.pi)
    core = (3/8*np.sinc(h) + (np.sinc(h-1)+np.sinc(h+1))/4
            + (np.sinc(h-2)+np.sinc(h+2))/16)
    tail = abs(x) > 5*np.pi
    v = x[tail]
    core[tail] = 48*np.pi**4*np.sin(v/2)/(v*(v**2-4*np.pi**2)*(v**2-16*np.pi**2))
    return duration*core*np.exp(.5j*x)


def freeze(model, seed, native, arm, start, stop):
    loga,phi=model.wall(native,arm,seed)
    quiet=model.wall(native,"quiet",seed)[1]
    first,last=(brentq(lambda s: phi(s)-quiet(u),phi.x[0],phi.x[-1]) for u in (start,stop))
    q=np.asarray(seed["q"])
    def rhs(s,state):
        h,k,ph,pk=state.reshape(4,-1)
        a=np.exp(loga(s))
        stiffness=a*(q*q+(100*a*(phi(s)-.5))**2)
        return np.r_[ph/a**3,pk/a**3,-stiffness*h,-stiffness*k]
    one=np.ones_like(q)
    zero=np.zeros_like(q)
    sol=solve_ivp(rhs,(seed["H0tau"],first),np.r_[one,zero,zero,one],
                  method="DOP853",rtol=2e-12,atol=1e-14,max_step=.002)
    if not sol.success:raise RuntimeError(sol.message)
    h,k,ph,pk=sol.y[:,-1].reshape(4,-1)
    if max(abs(h*pk-k*ph-1))>1e-10:raise ValueError("Canonical map failed")
    X0,C0,P0=(np.asarray(seed["states"]["retain"][j]) for j in ("X","C","P"))
    X=h*h*X0+2*h*k*C0+k*k*P0
    P=ph*ph*X0+2*ph*pk*C0+pk*pk*P0
    C=h*ph*X0+(h*pk+ph*k)*C0+k*pk*P0
    a=float(np.exp(loga(first)))
    e=np.hypot(q/a,100*(phi(first)-.5))
    inertia=a**3*e
    occ=(inertia*X+P/inertia-1)/2
    anomalous=(inertia*X-P/inertia)/2+1j*C
    assert min(occ)>-1.e-10
    return dict(a=a,q=q,weight=np.asarray(seed["weight"]),e=e,inertia=inertia,n=occ,
                c=anomalous,start=first,stop=last,duration=last-first)


def integrate(fr, cutoff, order, spacing):
    a,q,wq,e,inertia,n,c=(fr[key] for key in ("a","q","weight","e","inertia","n","c"))
    step=fr["duration"]
    x,g=roots_legendre(order)
    totals=np.zeros(8)
    for qi,w,energy,i,occ,anom in zip(q,wq,e,inertia,n,c):
        edges=np.linspace(qi,cutoff,int(np.ceil((cutoff-qi)/spacing))+1)
        half=np.diff(edges)[:,None]/2
        omega=(edges[:-1,None]+half*(x+1)).ravel()
        weight=w*(half*g).ravel()*(omega-qi)*(omega+qi)/(8*np.pi**2)
        cost=omega/a
        down=transform(cost-energy,step)/np.sqrt(2*i)
        up=transform(cost+energy,step)/np.sqrt(2*i)
        vac=abs(up)**2
        decay=occ*abs(down)**2
        stimulated=occ*vac
        mixing=2*np.real(np.conj(up)*down*anom)
        total=vac+decay+stimulated+mixing
        if min(total)<-1e-12:raise ValueError("Negative transition spectrum")
        for j,f in enumerate((vac, vac*cost, vac*energy, total, total*cost,
                               decay*cost, stimulated*cost, mixing*cost)):
            totals[j]+=weight@f
    labels=("vacuum_pairs","vacuum_daughter_energy","vacuum_parent_energy",
            "retained_pairs","retained_daughter_energy","annihilation_energy",
            "stimulated_creation_energy","anomalous_energy")
    out=dict(zip(labels,totals.tolist()))
    out["vacuum_switch_work"]=out["vacuum_parent_energy"]+out["vacuum_daughter_energy"]
    assert np.isclose(out["retained_daughter_energy"],sum(out[k] for k in
        ("vacuum_daughter_energy","annihilation_energy","stimulated_creation_energy","anomalous_energy")),rtol=3e-13)
    return out


def run(repo):
    model=load_source(repo)
    native=repo/model.FORK
    seed=json.loads((native/"FORK.json").read_text())["fine"]
    # Independently integrate the pulse at both removable singularities and off resonance.
    v=np.array([0.,3.,2*np.pi,4*np.pi,5*np.pi,7*np.pi,44.,87.])
    g,w=roots_legendre(100)
    t=(g+1)/2
    direct=np.exp(1j*np.outer(v,t))@(w/2*np.sin(np.pi*t)**4)
    assert max(abs(direct-transform(v,1.)))<2e-14
    results={}
    for arm in ("quiet","plus","minus","far"):
        rows={}
        for name,window in (("B",(2.055,2.145)),("early_half",(2.055,2.10)),("late_half",(2.10,2.145))):
            fr=freeze(model,seed,native,arm,*window)
            values={str(k):integrate(fr,k,16,12.) for k in (60,120,240,480,960,1920)}
            check=integrate(fr,1920,12,17.)
            err=max(abs(check[k]-v0)/max(1,abs(v0)) for k,v0 in values["1920"].items())
            assert err<1.e-11
            cut=1920.
            a,T,e=fr["a"],fr["duration"],fr["e"]
            assert np.min((cut/a+e)*T)>8*np.pi
            # A pointwise upper envelope for the analytic transform at x>=8 pi.
            constant=48*np.pi**4*(16/15)*(4/3)
            scale=fr["weight"]*(constant*T)**2*(a/T)**10/(16*np.pi**2*fr["inertia"])
            tail_n=float(np.sum(scale*cut**-7/7))
            tail_w=float(np.sum(scale*(cut**-6/(6*a)+e*cut**-7/7)))
            rows[name]=dict(start=fr["start"],stop=fr["stop"],scale_factor=a,
                  frozen_parent_energy_range=[float(min(e)),float(max(e))],
                  cutoffs=values,quadrature_error=err,vacuum_pairs_tail_bound=tail_n,
                  vacuum_work_tail_bound=tail_w)
        results[arm]=rows
    return results


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(run(args.repo.resolve()),indent=2,allow_nan=False)+"\n")
