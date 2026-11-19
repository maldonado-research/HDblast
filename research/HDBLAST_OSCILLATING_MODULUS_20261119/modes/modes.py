"""Vacuum production of a massless real brane scalar in a compact Jordan-metric pulse."""
import argparse
import csv
import importlib.util
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre


def frequencies(cycles, order, kmax=4.):
    # The occupation line is centered at 2k=1 with harmonic spacing 1/N.
    edges = np.arange(round(2*cycles*kmax)+1,dtype=float)/(2*cycles)
    if edges[-1] != kmax:
        raise ValueError("Cutoff must end a frequency panel")
    u,w = roots_legendre(order)
    width = np.diff(edges)[:,None]/2
    return (edges[:-1,None]+width*(1+u)).ravel(),(width*w).ravel()


def packet(t,cycles,power):
    half = .5/cycles
    st,ct = np.sin(half*t),np.cos(half*t)
    w = st**power
    wp = half*power*st**(power-1)*ct
    wpp = half**2*power*((power-1)*st**(power-2)*ct**2-st**power)
    c,s = np.cos(t),np.sin(t)
    return w*c,wp*c-w*s,wpp*c-2*wp*s-w*c


def effective_time(cycles,power):
    return 2*np.pi*cycles*comb(2*power,power)/4**power


def finite_born(momenta,cycles,power,epsilon,xi):
    """Use the exact Fourier transform of the compact cosine-polynomial pulse."""
    T = 2*np.pi*cycles
    q = 2*momenta
    transform = np.zeros(len(q),complex)
    for j in range(power+1):
        harmonic = (j-power/2)/cycles
        coefficient = (-1)**(j-power//2)*comb(power,j)/2**power
        for carrier in (-1,1):
            mismatch = q-carrier-harmonic
            transform += (.5*T*coefficient*np.exp(-.5j*T*mismatch)
                          *np.sinc(.5*T*mismatch/np.pi))
    # β1=(2ik)^−1 ∫dt e^-2ikt V(t), V1=ε(1−6ξ)f''(t).
    return 2j*momenta*epsilon*(1-6*xi)*transform


def born_number_sum_rule(cycles,power,epsilon,xi):
    # Parseval over 2k>0 converts the Born particle count to ∫V_lin²/(16π).
    # Distinct integer Fourier harmonics of f=w cos(t) are orthogonal over T.
    if cycles<=power:
        raise ValueError("Fourier harmonics overlap for this short a pulse")
    total=0.
    for j in range(power+1):
        offset=(j-power/2)/cycles
        coefficient=comb(power,j)/2**power
        total += coefficient**2*((1+offset)**4+(1-offset)**4)/4
    return (1-6*xi)**2*epsilon**2*total*(2*np.pi*cycles)/(16*np.pi)


def exact_modes(k,cycles,power,epsilon,xi,metric,rtol=1e-9,atol=3e-13):
    """α and β in the interaction picture with an exact static out basis."""
    T = 2*np.pi*cycles
    initial = np.r_[np.ones(len(k),complex),np.zeros(len(k),complex)]
    def rhs(t,y):
        f,fp,fpp = packet(t,cycles,power)
        V = (1-6*xi)*metric(f,fp,fpp,epsilon)
        alpha,beta = np.split(y,2)
        wave = np.exp(-2j*k*t)
        source = V/(2j*k)
        return np.r_[-source*(alpha+wave.conj()*beta),source*(wave*alpha+beta)]
    sol = solve_ivp(rhs,(0.,T),initial,method="DOP853",t_eval=[T],
                    rtol=rtol,atol=atol,max_step=.9)
    if not sol.success:
        raise ValueError(sol.message)
    alpha,beta = np.split(sol.y[:,-1],2)
    return alpha,beta


def direct_modes(k,cycles,power,epsilon,xi,metric):
    """Independent quadratures in the original second-order ode, on selected momenta."""
    T = 2*np.pi*cycles
    y0 = np.r_[np.ones(len(k),complex),-1j*k]
    def rhs(t,y):
        u,v=np.split(y,2)
        V = (1-6*xi)*metric(*packet(t,cycles,power),epsilon)
        return np.r_[v,(V-k*k)*u]
    sol=solve_ivp(rhs,(0.,T),y0,method="DOP853",rtol=3e-12,atol=2e-14,t_eval=[T],max_step=.4)
    if not sol.success:
        raise ValueError(sol.message)
    u,v=np.split(sol.y[:,-1],2)
    return .5*np.exp(-1j*k*T)*(u+v/(1j*k))


def flat_width_energy(cycles,power,epsilon,xi):
    # m=1. Γ=m³ [(1−6ξ)d ln A/dΘ]²/(32π), ρosc=m² Φ²/2.
    return (1-6*xi)**2*epsilon**2*effective_time(cycles,power)/(64*np.pi)


def run(cycles,power,epsilon,xi,metric,order=5,rtol=1e-9,atol=3e-13):
    k,dk=frequencies(cycles,order)
    # Include k=m/2 as an unweighted diagnostic.
    kval=np.r_[k,.5]
    alpha,beta=exact_modes(kval,cycles,power,epsilon,xi,metric,rtol,atol)
    norm_error=float(np.max(abs(abs(alpha)**2-abs(beta)**2-1)))
    occupation=abs(beta[:-1])**2
    born=abs(finite_born(k,cycles,power,epsilon,xi))**2
    density=dk*k**3/(2*np.pi**2)
    n_density=dk*k**2/(2*np.pi**2)
    energy=float(density@occupation)
    born_energy=float(density@born)
    count=float(n_density@occupation)
    born_count=float(n_density@born)
    sum_rule=born_number_sum_rule(cycles,power,epsilon,xi)
    if sum_rule and abs(born_count/sum_rule-1)>5e-7:
        raise ValueError("Full-line particle number sum rule")
    width=flat_width_energy(cycles,power,epsilon,xi)
    results=dict(cycles=cycles,window_sine_power=power,epsilon=epsilon,xi=xi,
                 effective_duration=effective_time(cycles,power),
                 energy_density=energy,born_energy_density=born_energy,
                 golden_rule_energy_density=width,number_density=count,
                 born_number_full_line=sum_rule,born_number_below_cutoff=born_count,
                 rms_energy_line_width=float(np.sqrt((density*occupation)@((k-.5)**2)/energy)) if energy else None,
                 mean_momentum=energy/count if count else None,
                 ratio_to_width=energy/width if width else None,
                 born_ratio_to_width=born_energy/width if width else None,
                 relative_to_finite_born=energy/born_energy if born_energy else None,
                 maximal_mode_occupation=float(max(occupation.max(),abs(beta[-1])**2)),
                 occupation_at_half_mass=float(abs(beta[-1])**2),
                 maximum_Wronskian_residual=norm_error)
    results["energy_by_cutoff"]={str(cut):float(density[k<cut]@occupation[k<cut]) for cut in (1,2,4)}
    if norm_error>1e-7:
        raise ValueError("Bogoliubov normalization")
    return results,k,occupation,born

def tangent_metric(f,fp,fpp,epsilon):
    """a_J=exp(εf), stationary a_E=1. Return the exact a_J″/a_J."""
    return epsilon*fpp+(epsilon*fp)**2


def archived_metrics(archive):
    helper_path=Path(__file__).resolve().parents[1]/"background/background.py"
    spec=importlib.util.spec_from_file_location("hdblast_modulus_background",helper_path)
    helper=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    module,record=helper.source(archive)
    tables=module.build_tables()
    result={}
    for pm in (.5,-.5):
        label=f"phi_m={pm:g}"
        pot=helper.Potential(module,tables,pm,record[label])
        f0=float(pot.spline(pot.tmin)[1])
        alpha0=float(-.5*pot.spline(pot.tmin,1)[1]/f0)
        def exact_metric(w,wp,wpp,epsilon):
            # Θ is in units of M_P. Time is in units mσ⁻¹.
            amplitude=epsilon/alpha0
            th=pot.tmin+amplitude*w
            f=float(pot.spline(th)[1])
            df=float(pot.spline(th,1)[1])
            ddf=float(pot.spline(th,2)[1])
            log_deriv=-.5*df/f
            second_over_A=.75*(df/f)**2-.5*ddf/f
            return log_deriv*amplitude*wpp+second_over_A*(amplitude*wp)**2
        report=run(16,4,.001,0.,exact_metric)[0]
        result[label]=dict(signed_logarithmic_derivative=alpha0,
                           archived_derivative_magnitude=record[label]["alpha"],
                           canonical_amplitude_in_Mpl=.001/alpha0,report=report)
    return result


def experiment(archive):
    parameters = {
        "short":(8,4,.001,0.),"reference":(16,4,.001,0.),"long":(32,4,.001,0.),
        "half_amplitude":(16,4,.0005,0.),"double_amplitude":(16,4,.002,0.),
        "reversed":(16,4,-.001,0.),"taper":(16,6,.001,0.),
        "improved":(16,4,.001,1/12),"conformal":(16,4,.001,1/6),
        "flat":(16,4,0.,0.)}
    results={}
    rows=[]
    for name,p in parameters.items():
        report,k,occ,born=run(*p,tangent_metric)
        results[name]=report
        if name in ("short","reference","long"):
            rows.extend((name,ki,ni,bi) for ki,ni,bi in zip(k,occ,born))
    p=parameters["reference"]
    refined,_,_,_=run(*p,tangent_metric,order=7,rtol=2e-10,atol=6e-14)
    tests=dict(refined_quadrature_report=refined,
               refinement_relative_energy=float(refined["energy_density"]/results["reference"]["energy_density"]-1))
    compare_k=np.array([.125,.4,.5,.6,.9,1.])
    _,beta=exact_modes(compare_k,*p,tangent_metric,rtol=2e-11,atol=7e-15)
    independent=direct_modes(compare_k,*p,tangent_metric)
    tests["max_direct_beta_abs_error"]=float(max(abs(beta-independent)))
    tests["direct_mode_crosscheck"]=[dict(k=float(ki),interaction_picture_n=float(abs(x)**2),
                                          direct_oscillator_n=float(abs(y)**2))
                                     for ki,x,y in zip(compare_k,beta,independent)]
    sample_k,_=frequencies(16,5)
    pos_a,pos_b=exact_modes(sample_k,16,4,.001,0.,tangent_metric)
    neg_a,neg_b=exact_modes(sample_k,16,4,-.001,0.,tangent_metric)
    tests["reciprocal_metric_duality"]={
        "max_abs_alpha_difference":float(np.max(abs(pos_a-neg_a))),
        "max_abs_beta_sum":float(np.max(abs(pos_b+neg_b)))}
    if max(tests["reciprocal_metric_duality"].values())>3e-9:
        raise ValueError("Reciprocal scale-factor mode map")
    nonlinear=archived_metrics(archive)
    for r in nonlinear.values():
        r["energy_ratio_to_tangent"]=r["report"]["energy_density"]/results["reference"]["energy_density"]
    if (abs(tests["refinement_relative_energy"])>2e-6
        or tests["max_direct_beta_abs_error"]>2e-9
        or results["flat"]["energy_density"]!=0 or results["conformal"]["energy_density"]!=0):
        raise ValueError("Independent direct mode, refinement, or null test failed")
    return dict(model="mσ=1, aE=1, aJ=exp(ε w cos(t)), one real massless scalar",
                quadrature_order=5,k_max=4,cases=results,archived_A=nonlinear,checks=tests),rows


def write(output,report,rows):
    output.mkdir(parents=True,exist_ok=True)
    (output/"MODES.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    with (output/"LINES.csv").open("w",newline="") as file:
        writer=csv.writer(file,lineterminator="\n")
        writer.writerow(["case","k_over_m","exact_occupation","born_occupation"])
        writer.writerows(rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--archive",type=Path,default=Path(__file__).resolve().parents[1]/"archive")
    args=parser.parse_args()
    write(args.output,*experiment(args.archive))
