"""Chemical particle-count ceiling for arbitrary bounded post-gate histories."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre, zeta


QUANTUM = Path("research/HDBLAST_QUANTUM_MEMORY_20261119/quantum")
NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
READOUT = 2.16
CUTOFF = 120.


def rule(n):
    x,w=roots_legendre(n)
    return (1+x)/2,w/2


def spectrum_module(repo):
    location=repo/QUANTUM/"spectrum.py"
    spec=importlib.util.spec_from_file_location("previous_spectrum",location)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def initial_moments(seed,profile,count,coordinate,cutoff=CUTOFF):
    """N, one-charge conformal energy E, and inverse-third moment of pair mean energy."""
    eta,measure,h,k=profile
    q,wq=(np.asarray(seed[key]) for key in ("q","weight"))
    X,C,P=(np.asarray(seed["states"]["retain"][key]) for key in ("X","C","P"))
    u,w=rule(count)
    result=np.zeros(3)
    for j,qj in enumerate(q):
        if coordinate=="log":
            omega=qj*np.exp(u*np.log(cutoff/qj))
            dw=w*omega*np.log(cutoff/qj)
        else:
            omega=qj+(cutoff-qj)*u**2
            dw=w*2*(cutoff-qj)*u
        phase=np.exp(1j*omega[:,None]*eta)
        v,b=(phase@np.column_stack((h[j]*measure,k[j]*measure))).T
        power=X[j]*abs(v+(C[j]+0.5j)*b/X[j])**2+(P[j]-(C[j]**2+.25)/X[j])*abs(b)**2
        if np.min(power)<-1e-14:
            raise ValueError("Invalid vacuum production weight")
        density=wq[j]*dw*(omega-qj)*(omega+qj)*power/(8*np.pi**2)
        e=omega/2
        result+=np.array([np.sum(density),np.sum(density*e),np.sum(density/e**3)])
    return result


def covariance_after(seed,loga,phi,off,transfer):
    """Continue the canonical free-parent transfer matrix on the measured geometry."""
    q=np.asarray(seed["q"])
    X,C,P=(np.asarray(seed["states"]["retain"][name]) for name in ("X","C","P"))
    def rhs(t,state):
        h,k,ph,pk=state.reshape(4,-1)
        a=np.exp(loga(t))
        stiffness=a*(q**2+(a*100*(phi(t)-.5))**2)
        return np.r_[ph/a**3,pk/a**3,-stiffness*h,-stiffness*k]
    sol=solve_ivp(rhs,(off,READOUT),transfer(off)[:-1],method="DOP853",dense_output=True,
                  max_step=.001,atol=3e-14,rtol=1.5e-12)
    if not sol.success:
        raise ValueError(sol.message)
    h,k,ph,pk=sol.y.reshape(4,len(q),-1)
    determinant_error=float(np.max(abs(h*pk-k*ph-1)))
    if determinant_error>1e-10:
        raise ValueError("Canonical transfer determinant")
    def variance(times):
        h,k=sol.sol(times).reshape(4,len(q),-1)[:2]
        return X[:,None]*(h+C[:,None]*k/X[:,None])**2+(P-C**2/X)[:,None]*k**2
    return variance,determinant_error


def diagonal_weights(seed,cutoff=CUTOFF):
    q,wq=(np.asarray(seed[key]) for key in ("q","weight"))
    analytic=wq/np.pi**2*(np.log(cutoff/q)-.5*(1-(q/cutoff)**2))
    # Independent integrals using positive log-frequency weights.
    u,w=rule(48)
    omega=q[:,None]*np.exp(np.log(cutoff/q)[:,None]*u)
    widths=omega*np.log(cutoff/q)[:,None]*w
    direct=wq*np.sum(widths*(2/omega)**3*(omega-q[:,None])*(omega+q[:,None])/(8*np.pi**2),axis=1)
    np.testing.assert_allclose(direct,analytic,rtol=7e-14,atol=1e-19)
    return analytic


def available_resource(variance,weight,knots,lo,hi,order):
    nodes,ds=time_rule(knots,lo,hi,order)
    return float(np.sum(ds*np.sqrt(weight@variance(nodes))))


def time_rule(knots,lo,hi,order):
    u,w=rule(order)
    cuts=np.r_[lo,knots[(knots>lo)&(knots<hi)],hi]
    widths=np.diff(cuts)
    return (cuts[:-1,None]+widths[:,None]*u).ravel(),(widths[:,None]*w).ravel()


def original_time_norm(seed,loga,phi,gates,profile,order):
    nodes,ds=time_rule(phi.x,*gates,order)
    _,measure,h,k=profile
    X,C,P=(np.asarray(seed["states"]["retain"][j]) for j in ("X","C","P"))
    var=X[:,None]*(h+C[:,None]*k/X[:,None])**2+(P-C**2/X)[:,None]*k**2
    gate=measure/ds
    if min(gate)<0 or max(gate)>1+1e-12:
        raise ValueError("Archived gate amplitude")
    return float(np.sum(ds*np.exp(loga(nodes))*gate**2*(np.asarray(seed["weight"])@var)))


def all_frequencies(numerics, post_norm, before_norm, need, chemical, production_y):
    rows=[]
    for split,(early,available,integral_weight) in numerics.items():
        initial_tail=2*before_norm/(np.pi*split)
        post_tail=2*post_norm/(np.pi*split)
        original=early+initial_tail
        linear=available*np.sqrt(early)
        quadratic=available**2+post_tail
        upper=original+2*linear+quadratic
        gap=max(0.,need**2-original)
        min_peak=gap/(np.sqrt(linear**2+quadratic*gap)+linear) if gap>0. else 0.
        # In a weighted L2 budget T, the low-frequency late amplitude has norm <= sqrt(D T).
        linear_integral=np.sqrt(early*integral_weight)
        quadratic_integral=integral_weight+2/(np.pi*split)
        minimum_integral=(gap/(np.sqrt(linear_integral**2+quadratic_integral*gap)+linear_integral))**2 if gap>0. else 0.
        integral_upper=(np.sqrt(early)+np.sqrt(integral_weight*post_norm))**2+initial_tail+post_tail
        if integral_upper<upper*(1-1e-12):
            raise ValueError("Peak-bounded source does not obey the L2 budget")
        rows.append(dict(integral_weight=float(integral_weight),minimum_weighted_source=float(minimum_integral),
                         minimum_weighted_source_ratio=float(minimum_integral/post_norm),
                         unit_integral_budget_upper=float(abs(production_y)**.5*integral_upper**.25/chemical),
                         splitting_frequency=split,early_low_L=float(early),post_low_norm=available,
                         early_tail_L_upper=initial_tail,post_tail_L_per_peak_squared=post_tail,
                         L_all_upper=float(upper),num_to_chemical_upper=float(abs(production_y)**.5*upper**.25/chemical),
                         necessary_post_peak=float(min_peak)))
    best=min(rows,key=lambda r:r["num_to_chemical_upper"])
    return dict(splits=rows,lowest_upper=float(best["num_to_chemical_upper"]),
                lowest_upper_split=float(best["splitting_frequency"]),
                necessary_post_peak=max(r["necessary_post_peak"] for r in rows),
                minimum_weighted_source=max(r["minimum_weighted_source"] for r in rows),
                minimum_weighted_source_ratio=max(r["minimum_weighted_source_ratio"] for r in rows),
                unit_integral_budget_upper=min(r["unit_integral_budget_upper"] for r in rows))


def evaluate(repo):
    model=spectrum_module(repo)
    seed=json.loads((repo/NATIVE/"FORK.json").read_text())["fine"]
    quiet=model.wall(repo/NATIVE,"quiet",seed)[1]
    w=diagonal_weights(seed)
    archived=json.loads((repo/QUANTUM/"SPECTRUM.json").read_text())["arms"]
    g=2
    A=3*g*zeta(3)/(4*np.pi**2)*(240/(7*g*np.pi**2))**.75
    production_y=.1
    answer={}
    for arm in ("quiet","plus","minus","far"):
        loga,phi=model.wall(repo/NATIVE,arm,seed)
        gates,sol,_=model.fundamental(seed,loga,phi,quiet)
        off=gates[1]
        p=model.profile(seed,gates,sol,phi,quiet,5)
        N,E,L=initial_moments(seed,p,144,"log")
        pc=model.profile(seed,gates,sol,phi,quiet,3)
        alt=initial_moments(seed,pc,128,"square")
        change=float(max(abs(alt/[N,E,L]-1)))
        if change>7e-8:
            raise ValueError(f"Independent moment refinement {arm}: {change}")
        z=archived[arm]
        native=z["positive_frequency_cutoff_checks"]["120"]["archived"]
        if abs(2*E/z["scale_factor_at_readout"]-native)>1e-7:
            raise ValueError("Archived cutoff-energy discrepancy")
        if N**4/E**3>L*(1+2e-12):
            raise ValueError("Holder inequality")
        initial_norm=original_time_norm(seed,loga,phi,gates,p,5)
        initial_norm_check=original_time_norm(seed,loga,phi,gates,pc,3)
        if abs(initial_norm/initial_norm_check-1)>1e-9:
            raise ValueError("Early source time norm")
        variance,det_error=covariance_after(seed,loga,phi,off,sol)
        low_moments={CUTOFF:L}
        for split in (24.,48.,240.,480.):
            low_moments[split]=initial_moments(seed,p,192,"log",split)[2]
            alternate=initial_moments(seed,pc,160,"square",split)[2]
            if abs(alternate/low_moments[split]-1)>3e-8:
                raise ValueError(f"Reference inverse moment below {split}")

        supports={"whole_post_gate":(off,READOUT),
                  "quiet_scalar_markers_2.055_2.145":tuple(brentq(lambda t:float(phi(t)-quiet(marker)),off,READOUT)
                                                        for marker in (2.055,2.145))}
        sorted_low=sorted(low_moments)
        for lower,upper_split in zip(sorted_low[:-1],sorted_low[1:]):
            increase=low_moments[upper_split]-low_moments[lower]
            if increase < -1e-12 or increase > 2*initial_norm/(np.pi*lower)+1e-13:
                raise ValueError("Measured finite-band tail exceeds positive-frequency Parseval ceiling")
        bounds={}
        need=A**2/abs(production_y)
        for name,(lo,hi) in supports.items():
            B=available_resource(variance,w,phi.x,lo,hi,7)
            alt_B=available_resource(variance,w,phi.x,lo,hi,3)
            resource=np.sqrt(L)+B
            if abs(alt_B/B-1)>1e-9:
                raise ValueError("Convergence of post-gate continuum integral")
            nodes,ds=time_rule(phi.x,lo,hi,7)
            var=variance(nodes)
            instantaneous=np.exp(loga(nodes))*(np.asarray(seed["weight"])@var)
            post_norm=float(np.sum(ds*instantaneous))
            inputs={}
            for split,moment in sorted(low_moments.items()):
                power=diagonal_weights(seed,split)@var
                available=float(np.sum(ds*np.sqrt(power)))
                integral_weight=float(np.sum(ds*power/instantaneous))
                inputs[split]=(moment,available,integral_weight)
            continuum=all_frequencies(inputs,post_norm,initial_norm,need,A,production_y)
            bounds[name]=dict(continuum=continuum,post_fourier_time_norm_upper=post_norm,proper_interval=[float(lo),float(hi)],post_resource_per_unit_peak=B,
                              total_resource_upper=float(resource),num_to_chemical_upper=float(np.sqrt(abs(production_y)*resource)/A),
                              max_L_total=float(resource**2),
                              minimum_extra_peak_for_chemical_count=float(max(0.,need-np.sqrt(L))/B),
                              post_time_quadrature_relative_change=float(abs(alt_B/B-1)))
        answer[arm]=dict(original_fourier_time_norm=initial_norm,initial_N=float(N),initial_E=float(E),initial_pair_mean_energy_inverse_third=float(L),
                         num_to_chemical_initial=float(abs(production_y)**.5*N/(A*E**.75)),
                         marginal_added_energy_per_new_particle_threshold=float(4*E/(3*N)),
                         num_to_chemical_initial_holder_upper=float(abs(production_y)**.5*L**.25/A),
                         relative_moment_check=change,wronskian_error=det_error,bounds=bounds)
    return dict(production_y=production_y,g_per_sign=g,pair_energy_cutoff=CUTOFF,readout=READOUT,
                chemical_number_coefficient=float(A),required_sqrt_L_for_chemical_count=float(A**2/abs(production_y)),
                histories=answer)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(evaluate(args.repo),indent=2,allow_nan=False)+"\n")
