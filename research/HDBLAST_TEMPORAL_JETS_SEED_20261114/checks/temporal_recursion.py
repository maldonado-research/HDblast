#!/usr/bin/env python3
"""Finite exact temporal-recursion, nonlinear-chain and frame witnesses.

Arbitrary coefficient jets verify the universal product rule. A separate
nonlinear exact wave model tests state/source composition. These witnesses
do not evaluate the actual PDE coefficient bounds, boundary traces or C_J.
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

FACTORS=("max_recursion_order","mixed_principal_factor","acceleration_subtraction_factor",
 "leading_a_inverse_power","source_sign_from_Fbar","binomial_weight_mode","source_chain_factorial_mode",
 "frame_time_space_factor","frame_spatial_first_factor","frame_base_first_factor",
 "frame_connection_derivative_factor","frame_connection_squared_factor","frame_base_connection_factor",
 "horizontal_output_transport_factor","horizontal_input_transport_factor","horizontal_state_transport_factor")
MUTATIONS=(("mixed_principal_factor","1"),("acceleration_subtraction_factor","0"),
 ("leading_a_inverse_power",0),("source_sign_from_Fbar","1"),("binomial_weight_mode","unit"),
 ("source_chain_factorial_mode","unit"),("frame_time_space_factor","1"),
 ("frame_spatial_first_factor","1"),("frame_base_first_factor","0"),
 ("frame_connection_derivative_factor","0"),("frame_connection_squared_factor","0"),
 ("frame_base_connection_factor","0"),("horizontal_output_transport_factor","0"),
 ("horizontal_input_transport_factor","0"),("horizontal_state_transport_factor","0"))


def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+"\n").encode()
def clean(v):return S.cancel(S.expand(v))
def series(jets,s):return sum(v*s**j/S.factorial(j) for j,v in enumerate(jets))


def evaluate(claims,check):
    records=[]
    def need(name,ok):
        boolean=bool(ok);check("temporal_recursion/"+name,boolean)
        records.append({"name":name,"status":"pass" if boolean else "fail"})
        if not boolean:raise ValueError("temporal recursion rejects "+name)
    need("schema",claims.get("schema")=="hdblast-temporal-recursion-claims-v1" and set(claims)=={"schema",*FACTORS})
    need("exact_engine",S.__version__=="1.14.0" and mpmath.__version__=="1.3.0")
    need("recursion_order",claims["max_recursion_order"]==12)
    s=S.Symbol("s")
    aj=list(S.symbols("a0:13"));bj=list(S.symbols("b0:13"));cj=list(S.symbols("c0:13"))
    jj=list(S.symbols("J0:15"));dxj=list(S.symbols("DXJ0:14"));hessj=list(S.symbols("HJ0:13"));fj=list(S.symbols("Fbar0:13"))
    a_poly=series(aj,s);b_poly=series(bj,s);c_poly=series(cj,s)
    U_poly=series(jj,s);DXUs_poly=series(dxj[1:],s);HU_poly=series(hessj,s);Fbar_poly=series(fj,s)
    lhs=a_poly*S.diff(U_poly,s,2)
    rhs=2*b_poly*DXUs_poly+c_poly*HU_poly-Fbar_poly
    mixed=S.Rational(claims["mixed_principal_factor"]);subtract=S.Rational(claims["acceleration_subtraction_factor"])
    source_sign=S.Rational(claims["source_sign_from_Fbar"])
    inverse_power=claims["leading_a_inverse_power"]
    def weight(r,k):
        if claims["binomial_weight_mode"]=="binomial":return math.comb(r,k)
        if claims["binomial_weight_mode"]=="unit":return 1
        raise ValueError("unsupported binomial mode")
    for r in range(13):
        # The oracle differentiates literal Taylor polynomials then solves its
        # linear equation for the new acceleration jet. It contains no supplied
        # binomial coefficient, recurrence factor or subtraction sign.
        equation=S.diff(lhs-rhs,s,r).subs(s,0)
        oracle=S.solve(equation,jj[r+2])[0]
        production=aj[0]**inverse_power*(
            mixed*sum(weight(r,k)*bj[k]*dxj[r-k+1] for k in range(r+1))+
            sum(weight(r,k)*cj[k]*hessj[r-k] for k in range(r+1))+source_sign*fj[r]+
            subtract*sum(weight(r,k)*aj[k]*jj[r-k+2] for k in range(1,r+1)))
        need("arbitrary_jet_acceleration_order_"+str(r+2),clean(production-oracle)==0)
    # Generic product-bundle/base-connection identity at one chart point.
    Om=S.Matrix([[0,-S.Rational(1,7)],[S.Rational(1,7),0]])
    Om_x=S.Matrix([[0,-S.Rational(2,11)],[S.Rational(2,11),0]])
    sigma=S.Rational(1,5)
    u=S.Matrix(S.symbols("u0:2"));us=S.Matrix(S.symbols("us0:2"));ux=S.Matrix(S.symbols("ux0:2"))
    uss=S.Matrix(S.symbols("uss0:2"));usx=S.Matrix(S.symbols("usx0:2"));uxx=S.Matrix(S.symbols("uxx0:2"))
    av,bv,cv=S.symbols("a b c")
    covariant=-av*uss+2*bv*(usx+Om*us)+cv*(uxx+2*Om*ux-sigma*ux+(Om_x+Om*Om-sigma*Om)*u)
    coordinate=-av*uss+2*bv*usx+cv*uxx
    correction=(S.Rational(claims["frame_time_space_factor"])*bv*Om*us+
       cv*(S.Rational(claims["frame_spatial_first_factor"])*Om*ux+
           S.Rational(claims["frame_base_first_factor"])*sigma*ux+
           (S.Rational(claims["frame_connection_derivative_factor"])*Om_x+
            S.Rational(claims["frame_connection_squared_factor"])*Om*Om+
            S.Rational(claims["frame_base_connection_factor"])*sigma*Om)*u))
    need("complete_bundle_and_base_connection_source_correction",all(clean(v)==0 for v in covariant-coordinate-correction))
    # Nonzero-state horizontal derivative: parallel transport acts on output,
    # input and the state argument. First-order transport polynomials suffice
    # for an exact derivative at zero and do not assume K is linear in state.
    base=S.Symbol("horizontal_coordinate")
    v0s,v1s=S.symbols("horizontal_V0 horizontal_V1")
    state=S.Matrix([v0s,v1s])
    connection=S.Matrix([[0,-S.Rational(2,7)],[S.Rational(2,7),0]])
    coefficient=S.Matrix([[1+base+v0s*v0s+v1s,v0s*v1s+2*base*v0s],
                          [v0s**3-v1s*v1s+base*v1s,2+v1s*v1s+v0s+base*base]])
    transported_state=(S.eye(2)-base*connection)*state
    literal_transport=(S.eye(2)+base*connection)*coefficient.subs(
        {v0s:transported_state[0],v1s:transported_state[1]},simultaneous=True)*(S.eye(2)-base*connection)
    horizontal_oracle=literal_transport.applyfunc(lambda value:S.Poly(S.expand(value),base).coeff_monomial(base))
    coefficient0=coefficient.subs(base,0)
    state_term=sum((coefficient.diff(v)*(connection*state)[i] for i,v in enumerate((v0s,v1s))),S.zeros(2)).subs(base,0)
    horizontal_production=(coefficient.diff(base).subs(base,0)+
        S.Rational(claims["horizontal_output_transport_factor"])*connection*coefficient0+
        S.Rational(claims["horizontal_input_transport_factor"])*coefficient0*connection+
        S.Rational(claims["horizontal_state_transport_factor"])*state_term)
    need("nonzero_state_horizontal_derivative_all_slots",all(clean(v)==0 for v in horizontal_production-horizontal_oracle))
    need("nonzero_horizontal_state_term_witness",any(clean(v)!=0 for v in state_term))
    # Nonlinear exact model: U(s,x)=exp(s+kappa*x)v. Coefficients depend on
    # state and explicitly on s; F is at most quadratic in first derivatives.
    # The test is structural and is not the Einstein coefficient certificate.
    kappa=S.Rational(2,3);M=kappa*S.eye(2)+Om;N=M*M+Om_x-sigma*M
    v0=S.Matrix([S.Rational(1,7),-S.Rational(2,11)])
    V=[S.Symbol("V0"),S.Symbol("V1")];Ps=S.Matrix(S.symbols("Ps0:2"));Px=S.Matrix(S.symbols("Px0:2"))
    aa=2+s/11+V[0]**2;bb=V[0]+V[1]**2+s*s/13;cc=3+s/17+V[1]**2
    vv=S.Matrix(V)
    F=aa*Ps-2*bb*Px-cc*N*vv+(1+V[1])*(Ps[0]*Px-V[0]*M*vv)
    expU=S.exp(s)*v0
    substitution={V[i]:expU[i] for i in range(2)}|{Ps[i]:expU[i] for i in range(2)}|{Px[i]:(M*expU)[i] for i in range(2)}
    actual_a=aa.subs(substitution);actual_b=bb.subs(substitution);actual_c=cc.subs(substitution)
    actual_F=F.subs(substitution,simultaneous=True)
    need("nonlinear_model_exact_solution",all(S.simplify(v)==0 for v in actual_a*expU-2*actual_b*M*expU-actual_c*N*expU-actual_F))
    # Production source jets are obtained by polynomial composition of the
    # field and first-derivative jets. Oracle uses the actual exponential model.
    mode=claims["source_chain_factorial_mode"]
    if mode not in ("factorial","unit"):raise ValueError("bad source-chain factorial mode")
    polynomial_u=sum((s**j/(S.factorial(j) if mode=="factorial" else 1)*v0 for j in range(15)),S.zeros(2,1))
    polynomial_us=S.diff(polynomial_u,s)
    polynomial_sub={V[i]:polynomial_u[i] for i in range(2)}|{Ps[i]:polynomial_us[i] for i in range(2)}|{Px[i]:(M*polynomial_u)[i] for i in range(2)}
    polynomial_F=F.subs(polynomial_sub,simultaneous=True).applyfunc(S.expand)
    actual_jets=[[clean(S.diff(value,s,r).subs(s,0)) for r in range(13)] for value in (actual_a,actual_b,actual_c)]
    for r in range(13):
        source_oracle=S.Matrix([clean(S.diff(value,s,r).subs(s,0)) for value in actual_F])
        source_production=polynomial_F.applyfunc(lambda value:S.Poly(value,s).coeff_monomial(s**r)*S.factorial(r))
        need("nonlinear_first_jet_source_chain_order_"+str(r),all(clean(v)==0 for v in source_production-source_oracle))
        model_recovered=(mixed*sum((weight(r,k)*actual_jets[1][k]*M*v0 for k in range(r+1)),S.zeros(2,1))+
          sum((weight(r,k)*actual_jets[2][k]*N*v0 for k in range(r+1)),S.zeros(2,1))+
          source_oracle+subtract*sum((weight(r,k)*actual_jets[0][k]*v0 for k in range(1,r+1)),S.zeros(2,1)))/actual_jets[0][0]
        need("nonlinear_model_recovered_temporal_jet_"+str(r+2),all(clean(v)==0 for v in model_recovered-v0))
    return {"schema":"hdblast-temporal-recursion-replay-v1","status":"pass","claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
      "scope":"Finite formal product/source/frame identities through temporal order14. Manufactured model coefficients are auxiliary witnesses; no actual PDE coefficient ceiling, boundary condition for normal derivatives, continuum existence, C_J, or C_seed is certified."}


def controls(claims):
    records=[]
    for name,value in MUTATIONS:
        bad=copy.deepcopy(claims);bad[name]=value;payload=canonical(bad)
        try:evaluate(json.loads(payload),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name,"self_rehashed_sha256":hashlib.sha256(payload).hexdigest(),"diagnostic":str(exc),"rejected_by":"fresh exact temporal/source/frame oracle"})
        else:raise ValueError("surviving temporal factor fault: "+name)
    return {"schema":"hdblast-temporal-recursion-controls-v1","status":"pass","controls":records}
