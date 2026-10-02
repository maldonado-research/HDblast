#!/usr/bin/env python3
"""Exact identities and wrong-formula controls; no numerical response calls.

Checks deliberately use explicit exceptions and remain active under python -O.
The JSON output describes algebra, not independently executed numerical modes.
"""

import argparse
import json
from pathlib import Path
import sympy as S


def run_checks():
    checks = []
    mutations = []

    def equal(name, actual, expected):
        residual = S.simplify(S.expand(actual-expected))
        if residual != 0:
            residual = S.simplify(S.expand(residual.rewrite(S.exp)))
        if residual != 0:
            raise RuntimeError(f"Exact identity failed: {name}: {residual}")
        checks.append({"name": name, "passed": True})

    def reject(name, wrong, correct):
        residual = S.simplify(S.expand(wrong-correct))
        if residual == 0:
            raise RuntimeError(f"Mutation was not detected: {name}")
        mutations.append({"name": name, "detected": True, "exact_residual": str(residual)})

    eta = S.symbols("eta", negative=True)
    H, k, tau, M, K, r, T, A, scale = S.symbols("H k tau M K r T A scale", positive=True)
    s0, ds, dQ, n, phi, phi0, b = S.symbols("s0 ds dQ n phi phi0 b", real=True)
    pi = S.pi
    a = -1/(H*eta)
    equal("exact_canonical_plane_wave_frequency", a*a*2*H*H-S.diff(a,eta,2)/a, 0)
    ai, ao = S.symbols("a_initial a_observed", positive=True)
    equal("physical_interaction_conformal_factors", ai**4/(ao*ai)**2, ai**2/ao**2)
    equal("wick_spatial_commutator", 2*(S.exp(-2*S.I*k*tau)-S.exp(2*S.I*k*tau))/(8*pi**2), -S.I*S.sin(2*k*tau)/(2*pi**2))
    equal("kubo_sign_and_normalization", (-S.I/2)*(-S.I/(2*pi**2)), -1/(4*pi**2))
    equal("independent_mode_interference", -S.sin(k*tau)*S.cos(k*tau)/k**2, -S.sin(2*k*tau)/(2*k**2))
    equal("mode_measure_reproduces_kubo", k**2/(2*pi**2)*(-1/(2*k**2)), -1/(4*pi**2))
    omega = S.symbols("omega", positive=True)
    equal("fixed_reference_subtraction_variation", -(s0/(2*omega))/(2*omega**2), -s0/(4*omega**3))
    cutoff_primitive = S.asinh(K/M)-K/S.sqrt(K*K+M*M)
    equal("finite_cutoff_local_primitive", S.diff(cutoff_primitive,K), K*K/(K*K+M*M)**S.Rational(3,2))
    equal("finite_cutoff_memory_primitive", S.diff((1-S.cos(2*K*tau))/(2*tau),K), S.sin(2*K*tau))
    equal("local_primitive_logarithmic_asymptotic_constant", S.limit(cutoff_primitive-S.log(2*K/M),K,S.oo), -1)
    equal("finite_part_log_memory_endpoint", S.limit((S.exp(-tau)-1)/tau,tau,0), -1)
    # Explicit stationary source primitive, then the past-infinite limit.
    stationary_primitive = -S.log(1+T/A)-T/(A+T)
    equal("stationary_history_primitive", S.diff(stationary_primitive,T), ((A/(A+T))**2-1)/T)
    stationary_limit = S.limit(stationary_primitive+S.log(S.sqrt(r)*T/(H*A))+S.EulerGamma+1,T,S.oo)
    stationary_qx = -stationary_limit/(8*pi**2)
    equal("fixed_r_stationary_susceptibility", stationary_qx, -(2*S.EulerGamma+S.log(r/H**2))/(16*pi**2))
    equal("reference_stationary_susceptibility", stationary_qx.subs(r,2*H**2), -(2*S.EulerGamma+S.log(2))/(16*pi**2))
    equal("independent_common_action_reference", ((1-2*S.EulerGamma)-S.log(2)-1)/(16*pi**2), stationary_qx.subs(r,2*H**2))
    # Frozen source formulas and exact root counts, without float evaluations.
    u, q = S.symbols("u q", real=True)
    bump = S.exp(1-1/(1-u*u))
    pulse_derivatives = {
        "positive_B": (
            -2*u*bump/(1-u*u)**2,
            2*(3*u**4-1)*bump/(1-u*u)**4,
            -4*u*(6*u**6+3*u**4-10*u*u+3)*bump/(1-u*u)**6),
        "signed_uB": (
            (u**4-4*u*u+1)*bump/(1-u*u)**2,
            2*u*(u**4+4*u*u-3)*bump/(1-u*u)**4,
            -2*(3*u**8+24*u**6-26*u**4+3)*bump/(1-u*u)**6),
    }
    for name, function in (("positive_B",bump),("signed_uB",u*bump)):
        for order, formula in enumerate(pulse_derivatives[name],1):
            equal(f"{name}_derivative_{order}", S.diff(function,u,order), formula)
    equal("positive_B_second_halfpoint",pulse_derivatives["positive_B"][1].subs(u,S.Rational(1,2)),-S.Rational(416,81)*S.exp(-S.Rational(1,3)))
    equal("signed_uB_second_halfpoint",pulse_derivatives["signed_uB"][1].subs(u,S.Rational(1,2)),-S.Rational(496,81)*S.exp(-S.Rational(1,3)))
    cubic = 2*q**3-14*q**2+21*q-6
    quartic = 4*q**4-32*q**3+64*q*q-36*q+3
    equal("positive_B_critical_root_count_above_one",S.Poly(cubic,q).count_roots(1,S.oo),2)
    equal("signed_uB_critical_root_count_above_one",S.Poly(quartic,q).count_roots(1,S.oo),2)
    equal("positive_B_extrema_primitive_derivative",S.diff(S.exp(1-q)*(4*q**4-12*q**3+6*q*q),q),-2*q*S.exp(1-q)*cubic)
    equal("signed_uB_extrema_primitive_derivative",S.diff(S.sqrt(1-1/q)*S.exp(1-q)*(4*q**4-12*q**3+2*q*q),q),-S.exp(1-q)*quartic/S.sqrt(1-1/q))
    # Quadratic mass law; phi is the inherited dimensionless bulk scalar.
    x = r*(1+b*(phi-phi0)/2)**2
    equal("mass_law_first_derivative",S.diff(x,phi).subs(phi,phi0),b*r)
    equal("mass_law_second_derivative",S.diff(x,phi,2),b*b*r/2)
    Q0 = H*H/(12*pi**2)
    dphi = S.symbols("dphi",real=True)
    current_variation = b*r*dQ/2+b*b*r*Q0*dphi/4
    equal("mass_law_current_contact",(S.diff(x,phi).subs(phi,phi0)*dQ+S.diff(x,phi,2)*Q0*dphi)/2,current_variation)
    equal("zero_coupling_current",current_variation.subs(b,0),0)
    # Wick's state occupation factor and initial Bogoliubov interference.
    equal("fixed_occupation_commutator_factor",(1+n)**2-n**2,1+2*n)
    equal("occupied_state_extra_kernel_factor",(-1/(4*pi**2))*(2*n),-n/(2*pi**2))
    betaR,betaI,Delta=S.symbols("betaR betaI Delta",real=True)
    beta=betaR+S.I*betaI
    v0=S.exp(-S.I*k*Delta)/S.sqrt(2*k)
    delta_v=beta*S.exp(S.I*k*Delta)/S.sqrt(2*k)
    equal("initial_state_interference",2*S.re(S.conjugate(v0)*delta_v),S.re(beta*S.exp(2*S.I*k*Delta))/k)
    equal("initial_state_momentum_measure",k*k/(2*pi**2)/k,k/(2*pi**2))
    equal("initial_bogoliubov_wronskian_linear_term",S.diff(1-S.Symbol("lambda",real=True)**2*(betaR*betaR+betaI*betaI),S.Symbol("lambda",real=True)).subs(S.Symbol("lambda",real=True),0),0)
    # Consistent dimensional rescaling and fixed-background Ward order.
    equal("coordinate_scale_a_invariance",-1/((H/scale)*(scale*eta)),a)
    equal("logarithm_dimensionless",(M/scale)*(scale*T),M*T)
    equal("response_mass_dimension_two",(s0/scale**2)*(M/scale)**2/(K/scale)**2,s0*M*M/(K*K*scale**2))
    equal("Ward_energy_dimension",1+4,2+2+1)
    lam,rate,dvar=S.symbols("lambda rate dvar",real=True)
    equal("Ward_linear_variance_excludes_deltaQ",S.diff((Q0+lam*dvar)*lam*rate/2,lam).subs(lam,0),Q0*rate/2)
    equal("conformal_Ward_conversion",ao*(Q0/2)*(rate/ao),Q0*rate/2)
    # Analytic spectral companion: compact endpoints eliminate boundary work.
    Mfun=S.Function("M")(eta)
    equal("reference_drift_contact_derivative",S.diff(S.log(Mfun),eta)/(32*pi**2),S.diff(Mfun,eta)/(32*pi**2*Mfun))
    equal("de_Sitter_reference_drift",S.diff(a*S.sqrt(r),eta)/(a*S.sqrt(r)),-1/eta)
    equal("finite_cutoff_reference_drift",S.diff(cutoff_primitive,M),-K**3/(M*(K*K+M*M)**S.Rational(3,2)))
    equal("canonical_pair_energy_change_of_frequency",S.Rational(1,8)/pi**2*(omega/2)*S.Rational(1,2),omega/(32*pi**2))
    correct_coefficient=-1/(4*pi**2)
    reject("reversed_Kubo_sign",-correct_coefficient,correct_coefficient)
    reject("missing_Wick_factor_two",correct_coefficient/2,correct_coefficient)
    reject("wrong_frequency_k_instead_of_2k",S.sin(k*tau),S.sin(2*k*tau))
    reject("omitted_finite_contact_one",stationary_qx+1/(8*pi**2),stationary_qx)
    reject("omitted_Euler_contact",stationary_qx+S.EulerGamma/(8*pi**2),stationary_qx)
    reject("moving_reference_r_equals_x",stationary_qx+1/(16*pi**2),stationary_qx)
    reject("omitted_mass_law_contact",b*r*dQ/2,current_variation)
    reject("omitted_occupation_state_factor",1,1+2*n)
    reject("omitted_initial_Bogoliubov_term",0,S.re(beta*S.exp(2*S.I*k*Delta))/k)
    reject("omitted_input_scale_factor",1/ao**2,ai**2/ao**2)
    reject("omitted_output_scale_factor",ai**2,ai**2/ao**2)
    reject("advanced_support_at_negative_separation",S.Heaviside(tau),S.Heaviside(-tau))
    reject("instantaneous_only_post_pulse_memory",0,-1/(8*pi**2*ao**2*tau))
    reject("omitted_time_dependent_reference_work",0,-1/(32*pi**2*eta))
    return {"scope":"Exact analytic identities only; no source or mode numerical evaluation", "passed":True,"identity_count":len(checks),"mutation_count":len(mutations),"checks":checks,"mutations":mutations}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=run_checks()
    encoded=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded,end="")


if __name__=="__main__":
    main()
