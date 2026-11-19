"""Frozen spatial frequencies of a minimally coupled daughter on an inhomogeneous closed slice."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.special import eval_gegenbauer

LABELS = (2, 3, 4, 7, 8, 9)


def gravity_module(root):
    spec = importlib.util.spec_from_file_location("coupled_sphere", root/"evolution/evolve.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def trial_basis(ell, maximum, y, weights):
    sine = np.sqrt(1-y*y)
    k = np.arange(maximum-ell+1)
    poly = eval_gegenbauer(k,ell+1,y[:,None])
    deriv = np.zeros_like(poly)
    deriv[:,1:] = 2*(ell+1)*eval_gegenbauer(k[1:]-1,ell+2,y[:,None])
    value = sine[:,None]**ell*poly
    dchi = sine[:,None]**(ell-1)*(ell*y[:,None]*poly-sine[:,None]**2*deriv)
    norm = np.sqrt(weights @ (value*value))
    return value/norm, dchi/norm


def calculate_modes(solver, state, t, maximum, points):
    y = np.cos(np.pi*np.arange(1,points+1)/(points+1))
    w = 1-y*y
    quadrature = 2*w/(points+1)
    vals = solver.matrices(y)[0] @ state.reshape(6,solver.n).T
    m = vals[:,0]
    ell_metric = m+w*vals[:,1]
    q = vals[:,4]
    branch=solver.branch
    kinetic = branch.fitted.spline(branch.top)[1]/branch.fitted.spline(branch.top+q)[:,1]
    if np.min(kinetic)<=0 or not np.isfinite(kinetic).all():
        raise ArithmeticError("Daughter kinetic weight left its physical domain")
    a=solver.background(t)[0]
    volume=np.exp(ell_metric+2*m)
    rvol=a*(quadrature@volume)**(1/3)
    rows={scheme:{n:[] for n in LABELS} for scheme in ("geometric","matter_weighted")}
    for ell in range(max(LABELS)+1):
        value,dchi=trial_basis(ell,maximum,y,quadrature)
        for scheme,z in [('geometric',np.ones_like(y)),('matter_weighted',kinetic)]:
            def gram(v,weight):
                return v.T @ (weight[:,None]*v)
            mass=gram(value,quadrature*z*volume)
            stiffness=gram(dchi,quadrature*z*np.exp(2*m-ell_metric))
            if ell:
                stiffness+=ell*(ell+1)*gram(value,quadrature*z*np.exp(ell_metric)/w)
            # Ordered radial eigenvalues continue the round labels n=ell,ell+1,... .
            vals=eigh(stiffness,mass,subset_by_index=[0,max(LABELS)-ell],eigvals_only=True)
            for n in LABELS:
                if n>=ell:
                    rows[scheme][n].append(float(np.sqrt(vals[n-ell])/a))
    return rvol,rows


def frequency_groups(values,radius):
    records={}
    for n,frequencies in values.items():
        ref=np.sqrt(n*(n+2))/radius
        weight=2*np.arange(n+1)+1
        f=np.asarray(frequencies)
        records[str(n)]=dict(round_at_volume_radius=ref,by_ell=f.tolist(),
                minimum=float(min(f)),maximum=float(max(f)),
                multiplet_width_fraction=float(np.ptp(f)/ref),
                weighted_mean_offset_fraction=float((weight @ (f/ref-1))/weight.sum()))
    return records


def linear_splitting(solver):
    branch=solver.branch
    f=branch.fitted.spline(branch.top)[1]
    slope=-branch.fitted.spline(branch.top,1)[1]/f
    eps=1e-5
    signed=[calculate_modes(solver,solver.initial(A,0),0,24,120)[1] for A in (eps,-eps)]
    rows={}
    for n in LABELS:
        lam=n*(n+2)
        ell=np.arange(n+1)
        theory=-2*slope*(1-2*ell*(ell+1)/lam)/lam
        measured=(np.array(signed[0]['matter_weighted'][n])-np.array(signed[1]['matter_weighted'][n]))/(2*eps*branch.h*np.sqrt(lam))
        error=float(np.max(abs(measured-theory)))
        if error>2e-7 or abs((2*ell+1)@theory)>1e-12:
            raise ArithmeticError("Perturbative splitting fails")
        rows[str(n)]=dict(analytic_frequency_derivative=theory.tolist(),max_finite_difference_error=error)
    return dict(d_log_kinetic_d_theta=slope,amplitude=eps,groups=rows)


def main(gravity):
    module=gravity_module(gravity)
    result={}
    for name,branch in module.constraint.load_branches(gravity/"previous").items():
        nonlinear=module.Evolution(branch,32,128)
        cases={}
        bg_check=module.Evolution(branch,48,192)
        for amplitude in (0.,.001,-.001):
            times=[0.,.12,.20]
            if branch.fitted.pm<0 and amplitude:
                times.append(.47)
            _,sol=nonlinear.evolve(amplitude,endpoint=max(times))
            _,refsol=bg_check.evolve(amplitude,endpoint=max(times),max_step=.003)
            runs=[]
            for t in times:
                state=sol.sol(t)
                vol,coarse=calculate_modes(nonlinear,state,t,20,96)
                radius,fine=calculate_modes(nonlinear,state,t,28,160)
                _,ref_background=calculate_modes(bg_check,refsol.sol(t),t,28,160)
                background_error=max(float(np.max(abs(np.asarray(fine[s][n])/ref_background[s][n]-1)))
                                     for s in fine for n in LABELS)
                error=max(float(np.max(abs(np.asarray(fine[s][n])/coarse[s][n]-1)))
                          for s in fine for n in LABELS)
                if background_error>3e-8 or error>5e-7 or abs(vol-radius)>1e-12:
                    raise ArithmeticError(f"Frozen spatial spectrum unresolved at time {t}: {error}")
                recorded={key:frequency_groups(val,radius) for key,val in fine.items()}
                if amplitude==0:
                    for scheme in recorded.values():
                        for row in scheme.values():
                            if row['multiplet_width_fraction']>2e-11 or abs(row['weighted_mean_offset_fraction'])>2e-11:
                                raise ArithmeticError("Round harmonic check")
                runs.append(dict(x=t,volume_radius=radius,relative_refinement_error=error,background_refinement_error=background_error,**recorded))
            cases[str(amplitude)]=runs
        result[name]=dict(cases=cases,linear_splitting=linear_splitting(nonlinear))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gravity",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(main(args.gravity.resolve()),indent=2,allow_nan=False)+"\n")
