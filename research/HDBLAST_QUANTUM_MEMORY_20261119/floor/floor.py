"""Minimize ordered high-band quantum energy over physical single-mode parent covariances."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
loader=importlib.util.spec_from_file_location("quantum_spectrum",HERE.parent/"quantum/spectrum.py")
spectrum=importlib.util.module_from_spec(loader)
loader.loader.exec_module(spectrum)


def scalar_covariances(X,C,P):
    return np.stack([np.stack([X,C],axis=-1),np.stack([C,P],axis=-1)],axis=-2)


def covariance_from_parameters(x):
    e,c=np.exp(x[0]),x[1]
    return .5*np.array([[e,c],[c,(1+c*c)/e]])


def numerical_minimum(g):
    # A log-variance and covariance parameterize every centered, pure, one-mode Gaussian.
    # The affine commutator drops out of the search. Rescale for a common derivative tolerance.
    matrix=g/np.sqrt(np.linalg.det(g))
    def cost(x):
        return float(np.sum(matrix*covariance_from_parameters(x)))
    def gradient(x):
        e,c=np.exp(x[0]),x[1]
        return np.array([.5*(matrix[0,0]*e-matrix[1,1]*(1+c*c)/e),
                         matrix[0,1]+matrix[1,1]*c/e])
    sol=minimize(cost,np.zeros(2),jac=gradient,method="BFGS",options={"gtol":3e-10})
    spectrum.require(np.linalg.norm(gradient(sol.x))<1e-7,"Unconverged independent minimizer")
    return covariance_from_parameters(sol.x)


def kernel(seed,wall,quiet,nt,nw):
    loga,phi,ends,free=wall
    q=np.array(seed["q"])
    a0=np.exp(loga(seed["H0tau"]))
    energy=np.hypot(q/a0,100*(phi(seed["H0tau"])-.5))
    inertia=a0**3*energy
    profiles=spectrum.profile(seed,ends,free,phi,quiet,nt)
    spec,vh,vk=spectrum.amplitudes(seed,profiles,*spectrum.BAND,nw)
    # A symplectic scaling makes the fork vacuum covariance I_2/2 for every q.
    basis=np.stack([vh/np.sqrt(inertia[:,None]),vk*np.sqrt(inertia[:,None])],axis=-1)
    M=np.sqrt(spec)[:,:,None]*basis
    H=M.conj().transpose(0,2,1) @ M
    g=H.real
    comm=-H[:,0,1].imag
    singular=np.linalg.svd(M,compute_uv=False)
    det_h=np.prod(singular,axis=1)**2
    delta=np.hypot(comm,np.sqrt(det_h))
    floor=delta+comm
    cancel=comm<0
    floor[cancel]=det_h[cancel]/(delta[cancel]-comm[cancel])
    optimum=(delta/2)[:,None,None]*np.linalg.inv(g)
    np.testing.assert_allclose(np.linalg.det(g),delta**2,rtol=2e-12,atol=1e-26)
    spectrum.require(np.all(singular[:,1]>0) and np.all(floor>0),"Degenerate fitted band")
    return dict(M=M,g=g,comm=comm,delta=delta,floor=floor,optimum=optimum,energy=energy,
                inertia=inertia,af=float(np.exp(loga(spectrum.READOUT))),sratio=singular[:,1]/singular[:,0])


def ordered_energy(model,V,pure=False):
    # Complete the square in v†(V+iJ/2)v to evaluate the ordered spectrum.
    # The second term is the positive residual of a mixed state's uncertainty determinant.
    h,k=model["M"][:,:,0],model["M"][:,:,1]
    X,C,P=V[:,0,0],V[:,0,1],V[:,1,1]
    residual=(X*P-C**2-.25)/X
    spectrum.require(min(residual)>-5e-12,"Invalid covariance")
    power=X[:,None]*abs(h+((C+.5j)/X)[:,None]*k)**2
    if not pure:
        power+=np.maximum(residual,0.)[:,None]*abs(k)**2
    return np.sum(power,axis=1)


def fork_energy(model,V,wq):
    # Proper parent Hamiltonian per comoving volume, at order y^0.
    return float(wq @ (model["energy"]*np.trace(V,axis1=1,axis2=2)/2))


def constrained(model,retained,wq):
    # One multiplier for the sum of the initial energies across all momenta.
    g=model["g"]
    k=wq*model["energy"]/2
    tr_g=np.trace(g,axis1=1,axis2=2)
    det_g=model["delta"]**2
    vacuum_level=sum(k)
    def excess(V):
        return float(k@(np.trace(V,axis1=1,axis2=2)-1))
    budget=excess(retained)
    def state(price):
        pen=price*k
        root=np.sqrt(det_g+pen*(tr_g+pen))
        shifted=g+pen[:,None,None]*np.eye(2)
        return root[:,None,None]/2*np.linalg.inv(shifted)
    spectrum.require(0<budget<excess(model["optimum"]),"Expected an active finite energy constraint")
    lo,hi=0.,1.
    while excess(state(hi))>budget:
        hi*=2
    # The upper bracket always has a feasible initial energy.
    for _ in range(55):
        mid=(lo+hi)/2
        if excess(state(mid))>budget:
            lo=mid
        else:
            hi=mid
    V=state(hi)
    pen=hi*k
    root=np.sqrt(det_g+pen*(tr_g+pen))
    root_minus_pen=(det_g+pen*tr_g)/(root+pen)
    dual=float(np.sum(root_minus_pen+model["comm"])-hi*budget)
    primal=float(sum(ordered_energy(model,V,pure=True)))
    spent=excess(V)
    spectrum.require(spent<=budget+1e-13 and abs(spent-budget)<1e-12, "Initial energy budget")
    spectrum.require(abs(primal-dual)<2e-11 and max(abs(np.linalg.det(V)-.25))<2e-13,
                     "Primal-dual or uncertainty mismatch")
    return V,dict(common_multiplier_conformal=hi,initial_parent_budget_proper=vacuum_level+budget,
                  allowed_excess_over_vacuum=budget,used_excess_over_vacuum=spent,
                  primal_proper_band=primal/model["af"],dual_proper_bound=dual/model["af"],
                  primal_minus_dual_proper=(primal-dual)/model["af"])


def analyze(repo):
    native=repo/spectrum.FORK
    seed=json.loads((native/"FORK.json").read_text())["fine"]
    wq=np.array(seed["weight"])
    quiet_phi=spectrum.wall(native,"quiet",seed)[1]
    models={}
    refinement={}
    for arm in ("quiet","plus","minus","far"):
        loga,phi=spectrum.wall(native,arm,seed)
        ends,free,_=spectrum.fundamental(seed,loga,phi,quiet_phi)
        wall=(loga,phi,ends,free)
        models[arm]=kernel(seed,wall,quiet_phi,5,96)
        trial=kernel(seed,wall,quiet_phi,3,48)
        refinement[arm]=float(max(abs(trial["floor"]/models[arm]["floor"]-1)))
        spectrum.require(refinement[arm]<2e-7,"Unresolved minimum energy")
    result={}
    fair_states={}
    for arm,m in models.items():
        g,comm=m["g"],m["comm"]
        V=m["optimum"]
        numerical=np.array([numerical_minimum(matrix) for matrix in g])
        analytic_energy=ordered_energy(m,V,pure=True)
        numerical_energy=ordered_energy(m,numerical,pure=True)
        residual=float(max(abs(analytic_energy/m["floor"]-1)))
        optimizer_disagreement=float(max(abs(numerical_energy/m["floor"]-1)))
        determinants=float(max(abs(np.linalg.det(V)-.25)))
        spectrum.require(residual<2e-11 and optimizer_disagreement<2e-10 and determinants<3e-13,
                         "Covariance optimizer or uncertainty check")
        positive=np.linalg.eigvalsh(g)
        spectrum.require(np.min(positive[:,0])>0,"Minimum assumes a strictly positive real Gram matrix")
        s=seed["states"]["retain"]
        X,C,P=(np.array(s[key]) for key in ("X","C","P"))
        retained=scalar_covariances(X*m["inertia"],C,P/m["inertia"])
        fair,certificate=constrained(m,retained,wq)
        fair_states[arm]=fair
        # Re-express the same quiet-optimal canonical covariance in this arm's fork normalization.
        scale=np.sqrt(m["inertia"]/models["quiet"]["inertia"])
        conversion=np.stack([scale,1/scale],axis=-1)
        fixed=models["quiet"]["optimum"]*conversion[:,:,None]*conversion[:,None,:]
        vac=np.broadcast_to(np.eye(2)/2,V.shape)
        common_fair=fair_states["quiet"]*conversion[:,:,None]*conversion[:,None,:]
        choices={"arm_optimum":V,"quiet_optimum":fixed,"retained":retained,"fork_vacuum":vac,
                 "budget_arm_optimum":fair,"budget_quiet_optimum":common_fair}
        comparisons={}
        for name,choice in choices.items():
            band=ordered_energy(m,choice,pure=name!="retained")
            spectrum.require(np.all(band>=m["floor"]*(1-5e-12)), "Violated mode floor")
            energy=fork_energy(m,choice,wq)
            comparisons[name]=dict(conformal_band=float(sum(band)),proper_band=float(sum(band)/m["af"]),
                                    fork_parent_proper=energy,
                                    fork_parent_excess_over_vacuum=energy-fork_energy(m,vac,wq),
                                    fork_parent_excess_over_retained=energy-fork_energy(m,retained,wq))
        qgap=comparisons["quiet_optimum"]["proper_band"]-comparisons["arm_optimum"]["proper_band"]
        floor_value=comparisons["arm_optimum"]["proper_band"]
        retained_value=comparisons["retained"]["proper_band"]
        fair_value=comparisons["budget_arm_optimum"]["proper_band"]
        common_value=comparisons["budget_quiet_optimum"]["proper_band"]
        spectrum.require(common_value>=fair_value-2e-12 and fair_value<=retained_value+2e-12,
                         "Fair-budget ranking")
        spectrum.require(comparisons["budget_quiet_optimum"]["fork_parent_excess_over_retained"]<2e-12,
                         "Quiet-optimal state exceeds this fork's budget")
        certificate.update(optimal_proper_fraction_of_retained=fair_value/retained_value,
                           quiet_optimum_proper_band=common_value,
                           quiet_optimum_relative_gap=common_value/fair_value-1,
                           optimal_covariance_dimensionless=fair.tolist())
        result[arm]=dict(fair_budget=certificate,analytic_conformal_floor=float(sum(m["floor"])),
                         analytic_proper_floor=float(sum(m["floor"])/m["af"]),
                         floor_fraction_of_retained_proper=floor_value/retained_value,
                         quiet_optimum_proper_gap=qgap,quiet_optimum_relative_gap=qgap/floor_value,
                         finite_momentum_mode_count=len(g),weights_in_gram=True,
                         minimum_weighted_conformal_floor_per_mode=float(min(m["floor"])),
                         covariance_max_abs_determinant_error=determinants,
                         direct_energy_max_relative_error=residual,
                         numeric_optimizer_max_relative_energy_error=optimizer_disagreement,
                         max_relative_floor_change_on_joint_refinement=refinement[arm],
                         real_gram_condition_numbers=[float(min(positive[:,1]/positive[:,0])),
                                                      float(max(positive[:,1]/positive[:,0]))],
                         complex_amplitude_singular_value_ratio=[float(min(m["sratio"])),float(max(m["sratio"]))],
                         preparations=comparisons,
                         per_mode_proper_floor=(m["floor"]/m["af"]).tolist(),
                         optimal_covariance_dimensionless=V.tolist())
    return dict(grid="fine",band=list(spectrum.BAND),native_q=seed["q"],proper_readout=spectrum.READOUT,
                lower_bound_target="ordered energy per comoving volume, at order y^2",
                initial_energy_cost="free parent proper energy per comoving volume at the fork, at order y^0",
                arms=result)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--repo",required=True,type=Path)
    cli.add_argument("--output",required=True,type=Path)
    args=cli.parse_args()
    args.output.write_text(json.dumps(analyze(args.repo),indent=2,allow_nan=False)+"\n")
