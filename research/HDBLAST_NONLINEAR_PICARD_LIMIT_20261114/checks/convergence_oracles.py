#!/usr/bin/env python3
"""Finite exact source, difference, Volterra and lower-jet witnesses.

Auxiliary polynomial maps and scalar kernels test algebraic formulas. These
finite tests do not prove infinite convergence, Sobolev compactness, nonlinear
existence, strong top-norm continuity or any physical evolved family.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import sys
import sympy as S
import mpmath

FACTORS=("maximum_iteration", "F_V_C_factor", "F_V_E_factor", "F_V_G_factor",
 "F_P_E_factor", "F_P_G_factor", "elapsed_firstjet_variation_factor",
 "segment_state_variation_factor", "segment_jet_variation_factor", "segment_integration_mode",
 "principal_mean_value_factor", "Picard_second_derivative_target", "uniqueness_second_derivative_target",
 "difference_A_factor", "difference_B_factor", "difference_source_factor", "seed_PDE_assumed",
 "outer_growth_weight_sign", "kernel_growth_difference_factor", "inner_predecessor_growth_sign",
 "Volterra_density_factor", "simplex_factorial_mode", "simplex_power_offset",
 "profile_monomial_factorial_mode", "nonzero_initial_energy_factor", "tail_start_offset",
 "tail_rising_denominator_offset", "tail_ratio_denominator_shift", "geometric_denominator_factor",
 "coefficient_iterate_shift", "source_elapsed_jet_increment", "successor_mixed_jet_increment",
 "lower_mixed_total_order", "strong_row_maximum", "weak_top_row", "weak_top_spatial_order",
 "coefficient_hit_unknown_order_ceiling", "source_intrinsic_order_ceiling", "source_individual_order_ceiling")
MUTATIONS=(("maximum_iteration",13),("F_V_C_factor","0"),("F_V_E_factor","0"),
 ("F_V_G_factor","0"),("F_P_E_factor","0"),("F_P_G_factor","1"),
 ("elapsed_firstjet_variation_factor","0"),("segment_state_variation_factor","0"),
 ("segment_jet_variation_factor","0"),("segment_integration_mode","midpoint"),
 ("principal_mean_value_factor","0"),("Picard_second_derivative_target","older_predecessor"),
 ("uniqueness_second_derivative_target","first_solution"),("difference_A_factor","0"),
 ("difference_B_factor","0"),("difference_source_factor","-1"),("seed_PDE_assumed",True),
 ("outer_growth_weight_sign","1"),("kernel_growth_difference_factor","0"),
 ("inner_predecessor_growth_sign","0"),("Volterra_density_factor","2"),
 ("simplex_factorial_mode","unit"),("simplex_power_offset",1),
 ("profile_monomial_factorial_mode","unit"),("nonzero_initial_energy_factor","0"),
 ("tail_start_offset",1),("tail_rising_denominator_offset",2),
 ("tail_ratio_denominator_shift",2),("geometric_denominator_factor","0"),
 ("coefficient_iterate_shift",0),("source_elapsed_jet_increment",0),
 ("successor_mixed_jet_increment",0),("lower_mixed_total_order",14),
 ("strong_row_maximum",12),("weak_top_row",11),("weak_top_spatial_order",0),
 ("coefficient_hit_unknown_order_ceiling",12),("source_intrinsic_order_ceiling",13),
 ("source_individual_order_ceiling",12))


def canonical(value):return (json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+"\n").encode()
def clean(value):return S.cancel(S.expand(value))
def zero(value):return all(clean(x)==0 for x in value) if isinstance(value,S.MatrixBase) else clean(value)==0
def taylor(jets,t):return sum(x*t**j/S.factorial(j) for j,x in enumerate(jets))


def evaluate(claims,check):
    records=[]
    def need(name,condition):
        passed=bool(condition);check("convergence_oracles/"+name,passed)
        records.append({"name":name,"status":"pass" if passed else "fail"})
        if not passed:raise ValueError("convergence oracle rejects "+name)
    need("schema",isinstance(claims,dict) and claims.get("schema")=="hdblast-convergence-oracle-claims-v1" and set(claims)=={"schema",*FACTORS})
    need("exact_engine",S.__version__=="1.14.0" and mpmath.__version__=="1.3.0")
    need("finite_iteration_order",claims["maximum_iteration"]==14)
    f=lambda name:S.Rational(claims[name])
    v=S.Matrix(S.symbols("state0:2"));p=S.Matrix(S.symbols("jet0:4"))
    dv=S.Matrix(S.symbols("dstate0:2"));dp=S.Matrix(S.symbols("djet0:4"))
    C=S.Matrix([v[0]+v[0]**3+v[0]*v[1]+v[1]**2,2*v[1]+v[0]**2*v[1]+v[1]**3])
    E=S.Matrix(2,4,lambda i,j:i+j+1+(j+1)*v[0]+S.Rational(i+j+1,7)*v[1]**2+v[0]*v[1]/11)
    G=[S.Matrix(4,4,lambda i,j:S.Rational((i+1)*(j+1)+o+1,13)+S.Rational(i+j+1,17)*v[0]+S.Rational((i-j)**2+1,19)*v[1]**2) for o in range(2)]
    F=C+E*p+S.Matrix([(p.T*M*p)[0] for M in G])
    eps=S.Symbol("variation_parameter")
    literal=F.subs({**{v[i]:v[i]+eps*dv[i] for i in range(2)},**{p[i]:p[i]+eps*dp[i] for i in range(4)}},simultaneous=True).diff(eps).subs(eps,0)
    stateC=sum((C.diff(v[i])*dv[i] for i in range(2)),S.zeros(2,1))
    stateE=sum((E.diff(v[i])*p*dv[i] for i in range(2)),S.zeros(2,1))
    stateG=S.Matrix([sum((p.T*M.diff(v[i])*p)[0]*dv[i] for i in range(2)) for M in G])
    effective_dp=S.Matrix([f("elapsed_firstjet_variation_factor")*dp[i] if i<2 else dp[i] for i in range(4)])
    variation=(f("F_V_C_factor")*stateC+f("F_V_E_factor")*stateE+f("F_V_G_factor")*stateG+
       f("F_P_E_factor")*E*effective_dp+f("F_P_G_factor")*S.Matrix([(p.T*M*effective_dp)[0] for M in G]))
    need("coupled_state_and_full_firstjet_directional_derivative",zero(variation-literal))
    need("nonzero_state_coefficient_contributions",all(not zero(term) for term in (stateC,stateE,stateG)))
    need("combined_firstjet_quadratic_symmetry",all(M==M.T for M in G))

    # Integrate the JOINT state/jet segment. Fixed rational endpoints keep this
    # independent literal composition small; the previous derivative is formal.
    theta=S.Symbol("theta")
    leftv=S.Matrix([S.Rational(1,7),-S.Rational(2,11)]);rightv=S.Matrix([-S.Rational(1,5),S.Rational(3,13)])
    leftp=S.Matrix([S.Rational(2,9),-S.Rational(1,8),S.Rational(3,10),S.Rational(1,6)])
    rightp=S.Matrix([-S.Rational(1,4),S.Rational(2,7),S.Rational(1,9),-S.Rational(2,5)])
    deltav=leftv-rightv;deltap=leftp-rightp
    segment={**{v[i]:rightv[i]+theta*deltav[i] for i in range(2)},**{p[i]:rightp[i]+theta*deltap[i] for i in range(4)},
       **{dv[i]:f("segment_state_variation_factor")*deltav[i] for i in range(2)},
       **{dp[i]:f("segment_jet_variation_factor")*deltap[i] for i in range(4)}}
    integrand=variation.subs(segment,simultaneous=True).applyfunc(S.expand)
    if claims["segment_integration_mode"]=="integral":integrated=integrand.applyfunc(lambda x:S.integrate(x,(theta,0,1)))
    elif claims["segment_integration_mode"]=="midpoint":integrated=integrand.subs(theta,S.Rational(1,2))
    else:raise ValueError("unsupported segment integration mode")
    endpoint=lambda vv,pp:F.subs({**dict(zip(v,vv)),**dict(zip(p,pp))},simultaneous=True)
    need("exact_joint_segment_mean_value_source",zero(endpoint(leftv,leftp)-endpoint(rightv,rightp)-integrated))
    A=1+v[0]**2+v[0]*v[1]+v[1]**3;B=2+v[0]**2+v[1]/5
    for name,coefficient in (("A",A),("B",B)):
        actual=coefficient.subs(dict(zip(v,leftv)),simultaneous=True)-coefficient.subs(dict(zip(v,rightv)),simultaneous=True)
        derivative=sum(coefficient.diff(v[i])*deltav[i] for i in range(2)).subs({v[i]:rightv[i]+theta*deltav[i] for i in range(2)},simultaneous=True)
        need("principal_mean_value_"+name,zero(actual-f("principal_mean_value_factor")*S.integrate(derivative,(theta,0,1))))

    # Separate exact residual identities distinguish Picard's preceding
    # solution from uniqueness's SECOND solution. No auxiliary germ is assumed
    # to solve the actual equations merely because it is a convenient witness.
    t,x=S.symbols("time space")
    U=S.Matrix([1+t+x+t*x+t*t+x**3,2+x*x+t**3+t*x*x])
    V=S.Matrix([2+t*t*x+x*x+t**3,1+t+x**3+t*t*x])
    W=S.Matrix([3+x+t*x*x+t*t,4+t*x+t**3+x*x])
    def coefficient(expr,state):return expr.subs(dict(zip(v,state)),simultaneous=True)
    def firstjet(state):return S.Matrix([*state.diff(t),*state.diff(x)])
    def source(state):return F.subs({**dict(zip(v,state)),**dict(zip(p,firstjet(state)))},simultaneous=True)
    def operator(predecessor,field):return field.diff(t,2)-coefficient(A,predecessor)*field.diff(x,t)-coefficient(B,predecessor)*field.diff(x,2)
    ru=operator(V,U)-source(V);rv=operator(W,V)-source(W)
    Picard_oracle=operator(V,U-V)-ru+rv
    targets={"preceding_solution":V,"older_predecessor":W,"successor":U}
    if claims["Picard_second_derivative_target"] not in targets:raise ValueError("unknown Picard derivative target")
    target=targets[claims["Picard_second_derivative_target"]]
    principal=lambda first,second,field:f("difference_A_factor")*(coefficient(A,first)-coefficient(A,second))*field.diff(x,t)+f("difference_B_factor")*(coefficient(B,first)-coefficient(B,second))*field.diff(x,2)
    need("Picard_difference_second_derivatives_of_preceding_solution",zero(Picard_oracle-f("difference_source_factor")*(source(V)-source(W))-principal(V,W,target)))
    rfirst=operator(V,V)-source(V);rsecond=operator(W,W)-source(W)
    unique_oracle=operator(V,V-W)-rfirst+rsecond
    unique_targets={"second_solution":W,"first_solution":V}
    if claims["uniqueness_second_derivative_target"] not in unique_targets:raise ValueError("unknown uniqueness derivative target")
    need("uniqueness_difference_second_derivatives_of_second_solution",zero(unique_oracle-f("difference_source_factor")*(source(V)-source(W))-principal(V,W,unique_targets[claims["uniqueness_second_derivative_target"]])))
    seed_residual=operator(W,W)-source(W)
    need("auxiliary_seed_not_silently_a_PDE_solution",zero(seed_residual) if claims["seed_PDE_assumed"] is True else not zero(seed_residual))

    # Nonconstant growth and integration density. Conjugation removes growth
    # ONCE, including the predecessor's inverse weight inside the integral.
    s,r=S.symbols("s r",nonnegative=True)
    H=lambda y:S.Rational(2,3)*y+y*y/10
    growth=S.diff(H(s),s);density=lambda y:S.Rational(3,7)*(1+y)
    primitive=S.integrate(density(r),(r,0,s))
    weighted_kernel=S.exp(f("outer_growth_weight_sign")*H(s))*S.exp(f("kernel_growth_difference_factor")*(H(s)-H(r)))*S.exp(f("inner_predecessor_growth_sign")*H(r))
    need("single_exact_exponential_kernel_conjugation",S.simplify(weighted_kernel-1)==0)
    need("nontrivial_growth_and_density",growth.diff(s)!=0 and S.diff(density(s),s)!=0)
    previous=S.Integer(1)
    simplex_mode=claims["simplex_factorial_mode"]
    if simplex_mode not in ("factorial","unit"):raise ValueError("bad simplex factorial mode")
    for n in range(15):
        if n:previous=S.integrate(S.expand(density(r)*previous.subs(s,r)),(r,0,s))
        denominator=S.factorial(n) if simplex_mode=="factorial" else 1
        proposed=S.exp(H(s))*(f("Volterra_density_factor")*primitive)**(n+claims["simplex_power_offset"])/denominator
        need("ordered_Volterra_integral_"+str(n),zero(S.expand(S.exp(-H(s))*proposed)-previous))
        if n:
            prior=S.exp(H(s))*primitive**(n-1)/S.factorial(n-1)
            need("direct_growth_ODE_recurrence_"+str(n),S.simplify(proposed.diff(s)-growth*proposed-density(s)*prior)==0)
        else:need("ordered_Volterra_initial_profile",proposed.subs(s,0)==1)

    # Independent nonconstant weighted starting profile with constant density.
    c=S.Rational(3,7);profile=1+s+2*s*s;recursive=profile
    for n in range(15):
        if n:recursive=S.integrate(c*recursive.subs(s,r),(r,0,s))
        mode=claims["profile_monomial_factorial_mode"]
        if mode not in ("factorial","unit"):raise ValueError("bad profile factorial mode")
        proposed=c**n*sum(a*(S.factorial(j) if mode=="factorial" else 1)*s**(n+j)/S.factorial(n+j) for j,a in enumerate((1,1,2)))
        need("nonconstant_profile_factorial_moment_"+str(n),zero(proposed-recursive))
    initial=S.Rational(5,9)
    stability=S.exp(H(s))*(f("nonzero_initial_energy_factor")*initial+s+s**3/3)
    need("nonzero_initial_difference_energy_retained",stability.subs(s,0)==initial)
    need("stability_with_nonconstant_growth_and_forcing",S.simplify(stability.diff(s)-growth*stability-S.exp(H(s))*(1+s*s))==0)

    # Finite tails only. Direct factorial sums and independent rising products
    # agree; quotient and geometric budgets are checked for every listed pair.
    z=S.Rational(7,3);tailcases=0
    for N in range(4,15):
        ratio=z/(N+claims["tail_ratio_denominator_shift"])
        need("finite_factorial_tail_ratio_"+str(N),ratio==(z**(N+1)/S.factorial(N+1))/(z**N/S.factorial(N)))
        for length in range(1,6):
            actual=sum(z**n/S.factorial(n) for n in range(N,N+length))
            start=N+claims["tail_start_offset"]
            budget=z**start/S.factorial(start)*sum(z**j/S.rf(N+claims["tail_rising_denominator_offset"],j) for j in range(length))
            need("finite_exact_tail_budget_"+str((N,length)),actual==budget)
            geom=sum(ratio**j for j in range(length))
            closed=(1-ratio**length)/(1-f("geometric_denominator_factor")*ratio)
            need("finite_geometric_budget_identity_"+str((N,length)),geom==closed)
            need("finite_positive_tail_ratio_budget_"+str((N,length)),0<ratio<1 and actual<=z**N/S.factorial(N)*geom)
            tailcases+=1
    for N in range(15):
        need("finite_factorial_binomial_tail_domination_"+str(N),all(S.factorial(N+j)>=S.factorial(N)*S.factorial(j) for j in range(15)))

    # Literal frozen Taylor coefficients record dependencies; they do not
    # interpolate unrelated time jets or prove a Sobolev product theorem.
    tau=S.Symbol("jet_time")
    predecessor=list(S.symbols("precedingJ0:14"));new=list(S.symbols("currentJ0:14"))
    mixed=list(S.symbols("spatial_timeJ0:14"));hessian=list(S.symbols("spatial_hessianJ0:14"))
    vv=taylor(predecessor,tau);vs=S.diff(vv,tau)
    aa=1+vv;bb=2+vv
    actual=S.Poly(S.expand(aa*taylor(mixed[1:],tau)+bb*taylor(hessian,tau)+vv+vv*vv+(1+vv)*vs+vs*vs),tau)
    shift=claims["coefficient_iterate_shift"]
    selected=predecessor if shift==-1 else new if shift==0 else None
    if selected is None:raise ValueError("unsupported coefficient iterate shift")
    selectedv=taylor(selected,tau);selecteds=S.diff(selectedv,tau)
    production=S.Poly(S.expand((1+selectedv)*taylor(mixed[1:],tau)+(2+selectedv)*taylor(hessian,tau)+selectedv+selectedv**2+(1+selectedv)*selecteds+selecteds**2),tau)
    lower=claims["lower_mixed_total_order"]
    need("strong_lower_row_roster",claims["strong_row_maximum"]==lower-2==11)
    need("weak_top_row_roster",claims["weak_top_row"]==12)
    dependencies=[]
    for row in range(13):
        expression=actual.coeff_monomial(tau**row)*S.factorial(row)
        need("frozen_coefficient_predecessor_shift_row_"+str(row),zero(expression-production.coeff_monomial(tau**row)*S.factorial(row)))
        vused=[j for j,symbol in enumerate(predecessor) if expression.has(symbol)]
        mused=[j for j,symbol in enumerate(mixed) if expression.has(symbol)]
        need("actual_temporal_dependency_row_"+str(row),max(vused)==row+claims["source_elapsed_jet_increment"] and max(mused)==row+claims["successor_mixed_jet_increment"])
        target=lower-(row+2)
        need("leading_spatial_order_jet_"+str(row+2),target==min(lower-(row+1)-1,lower-row-2))
        dependencies.append({"row":row,"new_jet":row+2,"preceding_jet_max":max(vused),"successor_mixed_jet_max":max(mused),"target_spatial_order":target,"passage":"strong" if row<=11 else "weak_top_row"})
    need("weak_top_leading_only_Hminus1",claims["weak_top_spatial_order"]==min(0-1,1-2)==-1)
    hitmax=max(12-k+2 for k in range(1,13))
    need("top_coefficient_hits_lower_unknown_order",hitmax<=claims["coefficient_hit_unknown_order_ceiling"])
    allocations=[(a+1,12-a+1) for a in range(13)]
    need("weak_top_source_total_order",max(sum(pair) for pair in allocations)<=claims["source_intrinsic_order_ceiling"])
    need("weak_top_source_individual_order",max(max(pair) for pair in allocations)<=claims["source_individual_order_ceiling"])
    # A bare derivative cannot be promoted from L2 convergence to L2
    # convergence: amplitude 1/k has derivative amplitude one, whereas its
    # H^{-1} derivative square is 1/(1+k^2). This is an exact symbol identity.
    k=S.Symbol("frequency",positive=True);coordinate=S.Symbol("Fourier_coordinate",real=True)
    mode=S.exp(S.I*k*coordinate)/k
    derivative_amplitude=S.simplify(S.diff(mode,coordinate)/(S.I*S.exp(S.I*k*coordinate)))
    need("top_row_high_frequency_derivative_warning",zero(derivative_amplitude-1) and zero(derivative_amplitude**2/(1+k*k)-1/(1+k*k)))
    return {"schema":"hdblast-convergence-oracle-report-v1","status":"pass","claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
       "finite_exact":{"state_dimension":2,"firstjet_dimension":4,"output_dimension":2,"ordered_iterations":15,"tail_cases":tailcases,"temporal_dependencies":dependencies},
       "scope":"Finite exact coupled source derivatives/joint mean values, separate Picard and uniqueness identities, weighted polynomial Volterra recurrences through14, finite factorial tail budgets, and literal frozen Taylor dependencies. These auxiliary witnesses do not prove infinite convergence, Sobolev compactness, existence, strong C X14, parameter regularity or a physical evolved family."}


def controls(claims):
    records=[]
    for field,value in MUTATIONS:
        bad=copy.deepcopy(claims);bad[field]=value;raw=canonical(bad)
        try:evaluate(json.loads(raw),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":field,"self_rehashed_sha256":hashlib.sha256(raw).hexdigest(),"diagnostic":str(exc),"rejected_by":"fresh exact convergence oracle"})
        else:raise ValueError("surviving convergence fault: "+field)
    return {"schema":"hdblast-convergence-oracle-controls-v1","status":"pass","controls":records}
