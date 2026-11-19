"""Frequency correlations for a pair-channel vacuum and a freely evolving Gaussian parent."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre

NATIVE = "research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native"
SOURCE = "research/HDBLAST_SCALAR_CLOCK_GRAVITY_20261119/source/source.py"
ARMS = ("quiet", "minus", "plus", "far")
JITTER = (0.0, 0.03, 0.10)


def load_source(repo):
    spec = importlib.util.spec_from_file_location("historical_source", repo/SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def waveforms(source, native, arm, nt):
    q, weights, states, (start, loga, mass) = source.read_native(native, "fine", arm)
    calibration = source.read_native(native, "fine", "quiet")[-1][-1]
    clock = source.ScalarClock(calibration)
    opening, closing = [brentq(lambda s: float(mass(s))/100 + .5 - field, start, 2.10) for field in clock.support]
    nodes, values = roots_legendre(nt)
    # The archived wall is piecewise Hermite. Splitting at its knots resolves
    # the short pulse features on each perturbed trajectory.
    edges = np.unique(np.r_[opening, mass.x[(mass.x>opening)&(mass.x<closing)], closing])
    half = np.diff(edges)/2
    times = (edges[:-1,None] + half[:,None]*(nodes+1)).ravel()
    measure = (values*half[:,None]).ravel() * clock(mass(times)/100+.5)[0]
    initial = np.r_[np.ones_like(q), np.zeros_like(q), np.zeros_like(q), np.ones_like(q), 0.]
    def rhs(t, state):
        h,k,ph,pk = state[:-1].reshape(4,-1)
        a = np.exp(loga(t))
        stiff = a*(q*q+(a*mass(t))**2)
        return np.r_[ph/a**3, pk/a**3, -stiff*h,-stiff*k,1/a]
    sol = solve_ivp(rhs,(start,closing),initial,t_eval=times,method="DOP853",rtol=2e-12,
                    atol=2e-14,max_step=.002)
    if not sol.success:
        raise ValueError(sol.message)
    h,k,ph,pk = sol.y[:-1].reshape(4,len(q),len(times))
    if np.max(abs(h*pk-k*ph-1)) > 2e-10:
        raise ValueError("Canonical parent transfer")
    return q,weights,states["retain"],np.exp(loga(2.10)),measure,sol.y[-1],h,k


def band_grid(nodes_per_band):
    node, weight = roots_legendre(nodes_per_band)
    return np.r_[75+15*node,105+15*node],np.tile(15*weight,2)


def detector(mat,n):
    # The principal quiet eigenpacket fixes both disjoint filters and their relative phase.
    _, eigenvectors = np.linalg.eigh(mat)
    psi = eigenvectors[:,-1]
    left = psi[:n]/np.linalg.norm(psi[:n])
    right = psi[n:]/np.linalg.norm(psi[n:])
    return left,right


def ports(mat,left,right):
    n=len(left)
    a = np.vdot(left,mat[:n,:n]@left).real
    b = np.vdot(right,mat[n:,n:]@right).real
    interference = 2*np.vdot(left,mat[:n,n:]@right).real
    return np.asarray([a+b,interference])


def calculate(repo,n,nt):
    source=load_source(repo)
    omega,dw=band_grid(n)
    detectors=[]
    results={}
    for arm in ARMS:
        q,weights,state,scale,measure,eta,h,k=waveforms(source,repo/NATIVE,arm,nt)
        X,C,P=state
        transform=np.exp(1j*omega[:,None]*eta)*measure
        vv,kk=transform@h.T,transform@k.T
        sums = {sigma:np.zeros(2) for sigma in JITTER}
        distance = {sigma:0. for sigma in JITTER}
        per_mode=[]
        total_number=total_energy=dephased_energy=0.
        max_tail=max_minor=max_mixed=0.
        for i,(momentum,wq) in enumerate(zip(q,weights,strict=True)):
            coeffs=np.sqrt(dw*(omega**2-momentum**2)/(8*np.pi**2))
            amplitudes=coeffs*np.asarray([vv[:,i],kk[:,i]])
            wightman=np.array([[X[i],C[i]+.5j],[C[i]-.5j,P[i]]])
            mat=amplitudes.conj().T@wightman@amplitudes
            mat=(mat+mat.conj().T)/2
            eig=np.linalg.eigvalsh(mat)
            number=eig.sum()
            if min(eig) < -1e-11*number or number <= 0:
                raise ValueError("Nonpositive excitation matrix")
            max_tail=max(max_tail,max(abs(eig[:-2]))/number)
            max_minor=max(max_minor,eig[-2]/number)
            # The bound survives a mixed Gaussian preparation. Use V -> 2V at fixed commutator.
            mixed_wightman=2*wightman.real+1j*wightman.imag
            mixed=amplitudes.conj().T@mixed_wightman@amplitudes
            mixed_eig=np.linalg.eigvalsh((mixed+mixed.conj().T)/2)
            if np.max(abs(mixed_eig[:-2]))>1e-11*mixed_eig.sum():
                raise ValueError("Mixed preparation rank exceeds two")
            max_mixed=max(max_mixed,mixed_eig[-2]/mixed_eig.sum())
            if arm=="quiet":detectors.append(detector(mat,n))
            left,right=detectors[i]
            number=float(mat.trace().real)
            energy=float(omega@np.diag(mat).real/scale)
            alternative=mat.copy()
            alternative[:n,n:]=0
            alternative[n:,:n]=0
            dephased_energy+=wq*(omega@np.diag(alternative).real)/scale
            if abs(ports(alternative,left,right)[1])>1e-13*number:
                raise ValueError("Dephased control retained band interference")
            for sigma in JITTER:
                jittered=mat*np.exp(-.5*sigma**2*np.subtract.outer(omega,omega)**2)
                sums[sigma]+=wq*ports(jittered,left,right)
                # Difference from a cross-band-dephased state is block off-diagonal.
                # Its trace norm is twice the sum of singular values of the cross block.
                distance[sigma]+=wq*np.linalg.svd(jittered[:n,n:],compute_uv=False).sum()
            vals=ports(mat,left,right)
            per_mode.append(dict(q=float(momentum),occupation=number,proper_energy=energy,
                                 quiet_filter_capture=float(vals[0]/number),visibility=float(vals[1]/vals[0])))
            total_number+=wq*number
            total_energy+=wq*energy
        reports={}
        for sigma,(captured,diff) in sums.items():
            optimum=float(distance[sigma]/total_number)
            witness=float(abs(diff)/(2*total_number))
            if witness>optimum+1e-11 or optimum>.5+1e-11:
                raise ValueError("Trace-distance witness")
            reports[str(sigma)]=dict(trace_distance=optimum,witness=witness,
                                     visibility=float(diff/captured),capture=float(captured/total_number),
                                     plus=float((captured+diff)/2),minus=float((captured-diff)/2))
        results[arm]=dict(proper_comoving_energy=total_energy,occupation=total_number,
                          dephased_energy=dephased_energy,max_beyond_rank_two=max_tail,
                          max_second_eigenvalue_fraction=max_minor,
                          mixed_control_max_second_fraction=max_mixed,ports=reports,per_mode=per_mode)
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    low=calculate(args.repo,24,6)
    high=calculate(args.repo,48,10)
    errors={}
    for arm,hi in high.items():
        lo=low[arm]
        errors[arm]=dict(energy=abs(lo["proper_comoving_energy"]/hi["proper_comoving_energy"]-1),
                         max_visibility=max(abs(lo["ports"][k]["visibility"]-v["visibility"]) for k,v in hi["ports"].items()),
                         max_capture=max(abs(lo["ports"][k]["capture"]-v["capture"]) for k,v in hi["ports"].items()))
        if max(errors[arm].values())>2e-7 or hi["max_second_eigenvalue_fraction"]>1e-8:
            raise ValueError("Refinement or near-pure packet rank")
        if abs(hi["proper_comoving_energy"]-hi["dephased_energy"])>1e-12 or hi["max_beyond_rank_two"]>1e-10:
            raise ValueError("Invariant excitation energy/rank")
    args.output.write_text(json.dumps(dict(band_edges=[60,90,120],time_nodes_per_wall_interval=10,frequency_nodes_per_band=48,
                  jitter_standard_deviations=JITTER,arms=high,refinement=errors),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
