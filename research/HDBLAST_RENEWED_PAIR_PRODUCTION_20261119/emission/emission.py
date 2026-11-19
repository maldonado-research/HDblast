"""Dirac pair production and interference from separated real scalar-clock windows."""
import argparse
import csv
from dataclasses import dataclass
from functools import lru_cache
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre


NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
PRIOR = Path("research/HDBLAST_QUANTUM_MEMORY_20261119/quantum/spectrum.py")
TRANSPORT = Path("research/HDBLAST_POST_GATE_TRANSPORT_20261119/daughters/DISTRIBUTIONS.json")
ARMS = ("quiet", "plus", "minus", "far")


@dataclass(frozen=True)
class Gate:
    name: str
    lo: float
    hi: float


@lru_cache(maxsize=2)
def historical(repo):
    spec = importlib.util.spec_from_file_location("hdblast_ordered_spectrum", Path(repo)/PRIOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=12)
def parent(repo, arm, readout=2.16):
    repo = Path(repo)
    native = repo/NATIVE
    previous = historical(repo)
    seed = json.loads((native/"FORK.json").read_text())["fine"]
    loga,phi = previous.wall(native,arm,seed)
    quiet = previous.wall(native,"quiet",seed)[1]
    q = np.array(seed["q"])
    start = seed["H0tau"]
    if not start < readout <= min(phi.x[-1],quiet.x[-1]):
        raise ValueError("Requested history exceeds the archived walls")
    def rhs(s, state):
        h,k,ph,pk = state[:-1].reshape(4,-1)
        a = np.exp(loga(s))
        mass = 100*(phi(s)-.5)
        stiffness = a*(q**2+(a*mass)**2)
        return np.r_[ph/a**3,pk/a**3,-stiffness*h,-stiffness*k,1/a]
    z = np.zeros_like(q)
    one = np.ones_like(q)
    sol = solve_ivp(rhs,(start,readout),np.r_[one,z,z,one,0.],dense_output=True,
                    method="DOP853",max_step=.0015,rtol=2e-12,atol=3e-14)
    if not sol.success:
        raise ValueError(sol.message)
    h,k,ph,pk = sol.y[:-1].reshape(4,len(q),-1)
    error = np.max(abs(h*pk-k*ph-1))
    if error>3e-11:
        raise ValueError("Parent canonical transfer lost its Wronskian")
    return seed,loga,phi,quiet,sol.sol,float(error)


@lru_cache(maxsize=80)
def gate_profile(repo, arm, gate, time_order=5, readout=2.16):
    seed,loga,phi,quiet,sol,_ = parent(repo,arm,readout)
    start = seed["H0tau"]
    ends = tuple(brentq(lambda s: phi(s)-quiet(t),start,readout) for t in (gate.lo,gate.hi))
    interior = phi.x[(phi.x>ends[0]) & (phi.x<ends[1])]
    boundaries = np.r_[ends[0],interior,ends[1]]
    nodes,weights = roots_legendre(time_order)
    widths = np.diff(boundaries)[:,None]/2
    s = (boundaries[:-1,None]+widths*(1+nodes)).ravel()
    ds = (widths*weights).ravel()
    lower,upper = np.full_like(s,gate.lo),np.full_like(s,gate.hi)
    measured = phi(s)
    for _ in range(48):
        mid = (lower+upper)/2
        mask = quiet(mid)<measured
        lower,upper = np.where(mask,mid,lower),np.where(mask,upper,mid)
    reference_clock = (lower+upper)/2
    values = np.sin(np.pi*(reference_clock-gate.lo)/(gate.hi-gate.lo))**4
    fields = sol(s)
    h,k = fields[:-1].reshape(4,len(seed["q"]),len(s))[:2]
    X,C,P = (np.asarray(seed["states"]["retain"][name])[:,None] for name in ("X","C","P"))
    # Cholesky factor of W=[[X,C+i/2],[C-i/2,P]], valid also for a pure parent.
    rem = (X*P-C*C-.25)/X
    if np.min(rem) < -5e-12:
        raise ValueError("Unphysical ordered parent covariance")
    v0 = np.sqrt(X)*h+(C+.5j)/np.sqrt(X)*k
    v1 = np.sqrt(np.maximum(rem,0))*k
    B = np.sqrt(X*h*h+2*C*h*k+P*k*k) @ (abs(values)*ds)
    return dict(ends=ends,s=s,eta=fields[-1],weights=values*ds,v0=v0,v1=v1,B=B,
                exposure=float(ds@(values**2)))


def amplitudes(repo,arm,omega,gates,time_order=5,readout=2.16):
    """Complex ordered square-root amplitudes at requested frequencies (q,ω,g,channel)."""
    q = parent(repo,arm,readout)[0]["q"]
    omega = np.broadcast_to(omega, (len(q),np.shape(omega)[-1]))
    amps = np.empty((*omega.shape,len(gates),2),complex)
    shared = np.all(omega==omega[0:1,:],axis=0)
    for g, gate in enumerate(gates):
        p = gate_profile(repo,arm,gate,time_order,readout)
        kernels = np.stack((p["v0"],p["v1"]), axis=-1)*p["weights"][None,:,None]
        # A single conformal-time origin is used for every gate.
        if np.any(shared):
            packed = np.ascontiguousarray(kernels.transpose(1,0,2)).reshape(len(p["eta"]),-1)
            indices = np.flatnonzero(shared)
            for start in range(0,len(indices),256):
                block = indices[start:start+256]
                osc = np.exp(1j*omega[0,block,None]*p["eta"])
                combined = (osc@packed).reshape(len(block),len(q),2)
                amps[:,block,g,:] = combined.transpose(1,0,2)
        if not np.all(shared):
            for j in range(len(q)):
                osc = np.exp(1j*omega[j,~shared,None]*p["eta"])
                amps[j,~shared,g,:] = osc @ kernels[j]
    return amps


def pair_quadrature(q,cutoff,order):
    """Integrate only the physical region, splitting at common comoving frequencies."""
    nodes,weights = roots_legendre(order)
    edges = np.r_[0.,np.arange(10.,cutoff,10.),cutoff]
    edges = np.maximum(edges[None,:],q[:,None])
    width = np.diff(edges)[:,:,None]/2
    omega = edges[:,:-1,None]+width*(1+nodes)
    dw = width*weights
    # The negative-third moment is sensitive to ω≈q, especially at small q.
    # Log-frequency quadrature resolves the first [q,10] panel on its own scale.
    log_span = np.log(edges[:,1]/q)
    omega[:,0,:] = q[:,None]*np.exp(log_span[:,None]*(nodes+1)/2)
    dw[:,0,:] = omega[:,0,:]*log_span[:,None]*weights/2
    return omega.reshape(len(q),-1),dw.reshape(len(q),-1)


def to_json(matrix):
    return {"real":matrix.real.tolist(),"imag":matrix.imag.tolist()}


def source_metrics(N,E,a,coeffs):
    c = np.asarray(coeffs)
    n = float(np.real(c.conj()@N@c))
    e = float(np.real(c.conj()@E@c))
    density = 2*n/a**3
    rho = e/a**4
    # Thermal abundance of four massless equilibrium spin/charge degrees, μ=0.
    zeta3 = 1.2020569031595943
    equilibrium_n = 3*zeta3/np.pi**2*(rho*60/(7*np.pi**2))**.75
    return dict(pair_count=n,both_sign_conformal_energy=e,both_sign_proper_number_density=density,
                both_sign_proper_energy_density=rho,mean_proper_momentum=e/(2*n*a),
                equilibrium_number_ratio_per_sqrt_y=density/equilibrium_n)

EARLY = Gate("A",1.35,2.05)
GATES = (EARLY,Gate("B",2.055,2.145),Gate("first",2.055,2.100),Gate("last",2.100,2.145))
CUTOFFS = (60.,120.,240.,480.,960.,1920.)
READOUT = 2.16


def summarize(N,E,a):
    answer = {}
    def coeffs(primary,late,sign):
        c = np.zeros(4)
        c[0] = primary
        c[late] = sign
        return c
    answer["A"] = source_metrics(N,E,a,[1,0,0,0])
    for j, gate in enumerate(GATES[1:],1):
        key = gate.name
        answer[key] = source_metrics(N,E,a,coeffs(0,j,1))
        answer["A+"+key] = source_metrics(N,E,a,coeffs(1,j,1))
        answer["A-"+key] = source_metrics(N,E,a,coeffs(1,j,-1))
        # Independent random sign on B is the specified incoherent benchmark.
        answer["incoherent_A_"+key] = source_metrics(np.diag(N.diagonal()),np.diag(E.diagonal()),a,coeffs(1,j,1))
    return answer


def all_cutoffs(repo,arm,max_cut,frequency_order,time_order):
    seed,loga,_,_,_,_ = parent(repo,arm,READOUT)
    q,wq = np.array(seed["q"]),np.array(seed["weight"])
    a = float(np.exp(loga(READOUT)))
    omega,dw = pair_quadrature(q,max_cut,frequency_order)
    amps = amplitudes(repo,arm,omega,GATES,time_order,READOUT)
    weights = wq[:,None]*(omega-q[:,None])*(omega+q[:,None])*dw/(8*np.pi**2)
    out={}
    for cutoff in CUTOFFS:
        if cutoff>max_cut:
            break
        active = weights*(omega<=cutoff)
        N,E,T = [np.einsum("qwgc,qwhc,qw->gh",amps.conj(),amps,active*moment,optimize=True)
                 for moment in (1,omega,(.5*omega)**-3)]
        for mat in (N,E,T):
            if np.linalg.eigvalsh(mat).min() < -1e-12 or np.max(abs(mat-mat.T.conj()))>2e-11:
                raise ValueError("Loss of Hermiticity or positivity")
        out[str(int(cutoff))] = dict(pair_number_matrix=to_json(N),
                                     both_sign_conformal_energy_matrix=to_json(E),
                                     inverse_cubed_pair_mean_energy_matrix=to_json(T),
                                     profiles=summarize(N,E,a))
    return out


def spectra(repo,arm,degree,time_order):
    """Integrate the one-sign comoving number and energy over exact boosted k bins."""
    from numpy.polynomial import Chebyshev
    from scipy.fft import dct
    seed = parent(repo,arm,READOUT)[0]
    q,wq = np.array(seed["q"]),np.array(seed["weight"])
    u = np.cos(np.arange(degree+1)*np.pi/degree)
    omega = q[:,None]+(120-q)[:,None]*(u+1)/2
    amps = amplitudes(repo,arm,omega,GATES,time_order,READOUT)
    names = ["A"]
    powers = [np.sum(abs(amps[:,:,0])**2,axis=-1)]
    B = np.stack([gate_profile(repo,arm,g,time_order,READOUT)["B"] for g in GATES],axis=1)
    bounds = [np.dot(wq,B[:,0]**2)]
    for j,g in enumerate(GATES[1:],1):
        late = np.sum(abs(amps[:,:,j])**2,axis=-1)
        summed = np.sum(abs(amps[:,:,0]+amps[:,:,j])**2,axis=-1)
        diff = np.sum(abs(amps[:,:,0]-amps[:,:,j])**2,axis=-1)
        names.extend([g.name,"A+"+g.name,"A-"+g.name,"incoherent_A_"+g.name])
        powers.extend([late,summed,diff,powers[0]+late])
        both = np.dot(wq,(B[:,0]+B[:,j])**2)
        bounds.extend([np.dot(wq,B[:,j]**2),both,both,np.dot(wq,B[:,0]**2+B[:,j]**2)])
    phase_space = wq[:,None]*(omega-q[:,None])*(omega+q[:,None])/(8*np.pi**2)
    coeffs = dct(np.stack(powers)*phase_space[None,:,:],type=1,axis=2)/degree
    coeffs[:,:,[0,-1]] /= 2
    edges = np.arange(65,dtype=float)
    output = {}
    for i,name in enumerate(names):
        cumulative = np.zeros((2, len(edges)))
        for qj,row in zip(q,coeffs[i]):
            curve = Chebyshev(row,domain=[qj,120])
            omega_poly = Chebyshev([(120+qj)/2,(120-qj)/2],domain=[qj,120])
            upper1 = np.clip(2*edges-qj,qj,120)
            upper2 = np.clip(2*edges+qj,qj,120)
            primitives = [(curve*omega_poly**j).integ() for j in range(3)]
            full = [p(upper1)-p(qj) for p in primitives]
            part = [p(upper2)-p(upper1) for p in primitives]
            cumulative[0] += full[0]+((edges+qj/2)*part[0]-part[1]/2)/qj
            cumulative[1] += full[1]/2+((edges**2-qj**2/4)*part[0]+qj*part[1]/2-part[2]/4)/(2*qj)
        count,energy = np.diff(cumulative)
        if np.min(count)<-5e-13 or np.min(energy)<-1e-11:
            raise ValueError("Negative integrated daughter spectrum")
        count,energy = np.maximum(count,0),np.maximum(energy,0)
        f = count*3*np.pi**2/np.diff(edges**3)
        if max(f)>bounds[i]:
            raise ValueError("Occupation exceeded the norm bound")
        output[name] = dict(count=count,energy=energy,occupancy=f,bound=float(bounds[i]))
    return edges,output


def experiment(repo):
    repo=Path(repo).resolve()
    prior=json.loads((repo/TRANSPORT).read_text())["arms"]
    result={"readout":READOUT,"gates":{g.name:[g.lo,g.hi] for g in GATES},
            "matrix_basis_order":[g.name for g in GATES],"cutoffs":list(CUTOFFS),"arms":{}}
    rows=[]
    for arm in ARMS:
        seed,loga,_,_,_,det_error = parent(repo,arm,READOUT)
        a=float(np.exp(loga(READOUT)))
        # Finer orders cover the short secondary controls and high-frequency tail.
        finest=all_cutoffs(repo,arm,max(CUTOFFS),16,6)
        coarse=all_cutoffs(repo,arm,max(CUTOFFS),12,4)
        history={"scale_factor":a,"transfer_wronskian_error":det_error,"gates":{},
                 "cutoffs":finest}
        for gate in GATES:
            data=gate_profile(repo,arm,gate,6,READOUT)
            history["gates"][gate.name]={"proper_time_endpoints":list(data["ends"]),
                                        "proper_time_integral_squared_gate":data["exposure"]}
        for cut, tests in finest.items():
            original=tests["profiles"]["A"]
            if cut in ("60","120"):
                ref=prior[arm]["bands"]["low" if cut=="60" else "full"]["one_sign_comoving_moments"]
                discrepancy=max(abs(original["pair_count"]/ref["0"]-1),
                                abs(original["both_sign_conformal_energy"]/(2*ref["1"])-1))
                if discrepancy>1e-7:
                    raise ValueError("Original gate historical number or energy")
                tests["archived_A_max_relative_error"]=discrepancy
            error=max(abs(p[metric]/coarse[cut]["profiles"][name][metric]-1)
                      for name,p in tests["profiles"].items()
                      for metric in ("pair_count","both_sign_conformal_energy","equilibrium_number_ratio_per_sqrt_y"))
            tests["max_profile_relative_refinement"]=error
            # The held-out pulse walls contain native knot-scale structure whose
            # ultraviolet weight is explicit in this numerical refinement test.
            if error>5e-5:
                raise ValueError(f"Time/frequency quadrature {arm}/{cut}: {error}")
            invfine=np.asarray(tests["inverse_cubed_pair_mean_energy_matrix"]["real"])
            invcoarse=np.asarray(coarse[cut]["inverse_cubed_pair_mean_energy_matrix"]["real"])
            tests["inverse_third_matrix_relative_refinement"]=float(np.max(abs(invfine-invcoarse))/np.max(invfine))
        edges,bins=spectra(repo,arm,112,6)
        _,bins_coarse=spectra(repo,arm,80,4)
        scores={}
        for name,b in bins.items():
            p=finest["120"]["profiles"][name]
            n,e=p["pair_count"],p["both_sign_conformal_energy"]/2
            closure=max(abs(sum(b["count"])/n-1),abs(sum(b["energy"])/e-1))
            refinement=max(float(sum(abs(b[key]-bins_coarse[name][key]))/total)
                           for key,total in (("count",n),("energy",e)))
            if max(closure,refinement)>2e-7:
                raise ValueError(f"Angular-bin quadrature {arm}/{name}")
            scores[name]=dict(relative_bin_closure=closure,binned_relative_L1_refinement=refinement,
                              peak_bin_occupation=float(max(b["occupancy"])),pointwise_occupation_bound=b["bound"])
            for j,left in enumerate(edges[:-1]):
                rows.append([arm,name,left,edges[j+1],left/a,edges[j+1]/a,b["count"][j],b["energy"][j],b["occupancy"][j]])
        # This descriptive statistic was added after inspecting the stored spectra.
        # It uses every registered bin and conditions on the emission of one particle.
        pos = bins["A+B"]["count"]
        neg = bins["A-B"]["count"]
        mix = bins["incoherent_A_B"]["count"]
        norm_pos,norm_neg,norm_mix = (z/sum(z) for z in (pos,neg,mix))
        history["exploratory_momentum_contrasts"] = {
            "normalized_TV_either_sign":float(.5*np.sum(abs(norm_pos-norm_neg))),
            "normalized_TV_positive_vs_incoherent":float(.5*np.sum(abs(norm_pos-norm_mix))),
            "integrated_sign_contrast":float((sum(pos)-sum(neg))/(sum(pos)+sum(neg))),
            "relative_L1_signed_interference":float(np.sum(abs(pos-neg))/(2*sum(mix))),
            "contrast_in_quiet_selected_k25_to26":float((pos[25]-neg[25])/(pos[25]+neg[25])),
            "incoherent_fraction_in_k25_to26":float(mix[25]/sum(mix))}
        history["cutoff120_daughter_spectrum"]=scores
        result["arms"][arm]=history
    return result,rows


def save(result,rows,out):
    out.mkdir(parents=True,exist_ok=True)
    (out/"EMISSION.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    with (out/"SPECTRA.csv").open("w",newline="") as f:
        writer=csv.writer(f,lineterminator="\n")
        writer.writerow(["arm","source","k_low","k_high","p_low","p_high","one_sign_comoving_number",
                         "one_sign_conformal_energy","per_spin_bin_mean_occupation"])
        writer.writerows(rows)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    save(*experiment(args.repo),args.output)
