#!/usr/bin/env python3
"""Finite exact witnesses for the normalized frozen evolution identities.

The auxiliary jets, connection and manufactured polynomials are independent
algebraic witnesses. They do not certify continuous estimates, coefficient
ceilings, existence, boundary regularity or a nonlinear evolved family.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import sympy as S
import mpmath

FACTORS=("maximum_field_order", "energy_force_factor", "energy_mass_factor",
 "energy_div_A_factor", "energy_B_time_factor", "energy_div_B_factor",
 "energy_spurious_a_time_factor", "flux_A_factor", "flux_B_factor",
 "shell_outward_normal", "shell_A_normal", "shell_B_mixed",
 "shell_principal_fiber_coupling", "shell_D_time_trace", "shell_N_normal_trace",
 "normal_leading_inverse_power", "normal_Az_factor", "normal_Ay_factor",
 "normal_Bzy_factor", "normal_Byy_factor", "normal_source_factor",
 "normal_Bzz_commutator_factor", "normal_binomial_mode",
 "normal_divide_before_differentiate_mode", "commutator_sign",
 "commutator_word_mode", "difference_second_derivative_target",
 "difference_A_factor", "difference_B_factor", "difference_source_factor",
 "quadratic_source_cross_factor", "energy_source_order", "normal_source_order")
MUTATIONS=(("maximum_field_order",13), ("energy_force_factor","0"),
 ("energy_mass_factor","0"), ("energy_div_A_factor","1/2"),
 ("energy_B_time_factor","0"), ("energy_div_B_factor","0"),
 ("energy_spurious_a_time_factor","1/2"), ("flux_A_factor","1"),
 ("flux_B_factor","1/2"), ("shell_outward_normal","1"),
 ("shell_A_normal","1"), ("shell_B_mixed","1"),
 ("shell_principal_fiber_coupling","1/4"), ("shell_D_time_trace",False),
 ("shell_N_normal_trace",False), ("normal_leading_inverse_power",0),
 ("normal_Az_factor","0"), ("normal_Ay_factor","0"),
 ("normal_Bzy_factor","-1"), ("normal_Byy_factor","0"),
 ("normal_source_factor","1"), ("normal_Bzz_commutator_factor","0"),
 ("normal_binomial_mode","unit"),
 ("normal_divide_before_differentiate_mode","freeze_reciprocal"),
 ("commutator_sign","-1"), ("commutator_word_mode","sort"),
 ("difference_second_derivative_target","successor"),
 ("difference_A_factor","0"), ("difference_B_factor","0"),
 ("difference_source_factor","-1"), ("quadratic_source_cross_factor","1/2"),
 ("energy_source_order",12), ("normal_source_order",11))


def canonical(value):return (json.dumps(value,sort_keys=True,indent=2)+"\n").encode()
def clean(value):return S.cancel(S.expand(value))
def zero(value):return all(clean(item)==0 for item in value) if isinstance(value,S.MatrixBase) else clean(value)==0
def polynomial(jets,x):return sum(value*x**j/S.factorial(j) for j,value in enumerate(jets))


def evaluate(claims,check):
    records=[]
    def need(name,ok):
        boolean=bool(ok);check("evolution_oracles/"+name,boolean)
        records.append({"name":name,"status":"pass" if boolean else "fail"})
        if not boolean:raise ValueError("evolution oracle rejects "+name)
    need("schema",isinstance(claims,dict) and claims.get("schema")=="hdblast-evolution-oracle-claims-v1" and set(claims)=={"schema",*FACTORS})
    need("exact_engine",S.__version__=="1.14.0" and mpmath.__version__=="1.3.0")
    need("maximum_field_order",claims["maximum_field_order"]==14)
    factor=lambda name:S.Rational(claims[name])

    # Literal differentiation in four spatial directions, with a nontrivial
    # metric-compatible bundle connection. The connection preserves the split
    # span(e0) + span(e1,e2,e3), but angular matrices do not commute.
    s=S.Symbol("s");xs=S.symbols("z x y w");args=(s,*xs)
    U=S.Matrix([S.Function("u"+str(i))(*args) for i in range(4)])
    forcing=S.Matrix([S.Function("f"+str(i))(*args) for i in range(4)])
    Om=[S.zeros(4) for _ in range(4)]
    for index,a,b,c in ((1,1,2,S.Rational(1,3)),(2,2,3,S.Rational(2,5)),(3,1,3,S.Rational(3,7))):
        Om[index][a,b]=-c;Om[index][b,a]=c
    def D(i,value):return value.diff(xs[i])+Om[i]*value
    Us=U.diff(s);DU=[D(i,U) for i in range(4)]
    A=[S.Function("A"+str(i))(*args) for i in range(4)]
    B=[[None]*4 for _ in range(4)]
    for i in range(4):
        for j in range(i,4):B[i][j]=B[j][i]=S.Function("B"+str(i)+str(j))(*args)
    acceleration=forcing+sum((A[i]*D(i,Us) for i in range(4)),S.zeros(4,1))+sum((B[i][j]*D(i,DU[j]) for i in range(4) for j in range(4)),S.zeros(4,1))
    density=(Us.dot(Us)+U.dot(U)+sum(B[i][j]*DU[i].dot(DU[j]) for i in range(4) for j in range(4)))/2
    literal=density.diff(s).subs({S.diff(U[k],s,2):acceleration[k] for k in range(4)},simultaneous=True)
    flux=[factor("flux_A_factor")*A[i]*Us.dot(Us)+factor("flux_B_factor")*sum(B[i][j]*DU[j].dot(Us) for j in range(4)) for i in range(4)]
    production=(factor("energy_force_factor")*forcing.dot(Us)+factor("energy_mass_factor")*U.dot(Us)+
       factor("energy_div_A_factor")*sum(S.diff(A[i],xs[i]) for i in range(4))*Us.dot(Us)+
       factor("energy_B_time_factor")*sum(S.diff(B[i][j],s)*DU[i].dot(DU[j]) for i in range(4) for j in range(4))+
       factor("energy_div_B_factor")*sum(S.diff(B[i][j],xs[i])*DU[j].dot(Us) for i in range(4) for j in range(4))+
       factor("energy_spurious_a_time_factor")*S.Symbol("unnormalized_a_s")*Us.dot(Us))
    need("four_spatial_direction_covariant_energy_IBP",zero(literal-sum(S.diff(flux[i],xs[i]) for i in range(4))-production))
    need("connection_is_metric_and_split_preserving",all(M+M.T==S.zeros(4) and all(M[0,k]==0 and M[k,0]==0 for k in range(4)) for M in Om))
    need("angular_connection_is_noncommuting",Om[1]*Om[2]!=Om[2]*Om[1])

    # Boundary test keeps D-normal jets and N-time jets independent and free.
    V=S.Matrix(S.symbols("boundary_velocity0:4"));Z=S.Matrix(S.symbols("boundary_normal0:4"))
    tangential=[S.Matrix(S.symbols("boundary_tangent"+str(i)+"_0:4")) for i in range(3)]
    normalB=S.Symbol("boundary_Bzz",positive=True)
    coupling=S.eye(4)+factor("shell_principal_fiber_coupling")*S.Matrix([[0,1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]])
    unsigned=factor("shell_A_normal")*V.dot(V)/2+normalB*(coupling*Z).dot(V)+factor("shell_B_mixed")*sum(T.dot(V) for T in tangential)
    need("signed_shell_flux_before_boundary_conditions",zero(factor("shell_outward_normal")*unsigned+unsigned))
    substitutions={}
    if claims["shell_D_time_trace"] is True:substitutions[V[0]]=0
    if claims["shell_N_normal_trace"] is True:substitutions.update({Z[k]:0 for k in range(1,4)})
    need("mixed_shell_flux_with_free_complementary_jets",zero(unsigned.subs(substitutions,simultaneous=True)))
    need("outer_Dirichlet_flux",zero(unsigned.subs({v:0 for v in V},simultaneous=True)))

    # A literal polynomial-product oracle differentiated before isolating Bzz.
    z=S.Symbol("normal_coordinate")
    coeff=[list(S.symbols(prefix+"0:13")) for prefix in ("Bzzjet","Azjet","Ayjet","Bzyjet","Byyjet")]
    field=[list(S.symbols(prefix+"0:13")) for prefix in ("SSjet","ZSjet","YSjet","ZYjet","YYjet","Fjet")]
    J=list(S.symbols("normalJ0:15"))
    cp=[polynomial(v,z) for v in coeff];fp=[polynomial(v,z) for v in field]
    equation=S.Poly(S.expand(cp[0]*polynomial(J[2:],z)-fp[0]+cp[1]*fp[1]+cp[2]*fp[2]+2*cp[3]*fp[3]+cp[4]*fp[4]+fp[5]),z)
    need("normal_recovery_product_rule_mode",claims["normal_divide_before_differentiate_mode"] in ("unisolated","freeze_reciprocal"))
    def weight(r,k):
        if claims["normal_binomial_mode"]=="binomial":return math.comb(r,k)
        if claims["normal_binomial_mode"]=="unit":return 1
        raise ValueError("unsupported normal binomial mode")
    for r in range(13):
        direct=equation.coeff_monomial(z**r)*S.factorial(r)
        oracle=S.solve(direct,J[r+2])[0]
        row=field[0][r]+sum(factor(name)*sum(weight(r,k)*coeff[c][k]*field[f][r-k] for k in range(r+1)) for name,c,f in (("normal_Az_factor",1,1),("normal_Ay_factor",2,2),("normal_Bzy_factor",3,3),("normal_Byy_factor",4,4)))+factor("normal_source_factor")*field[5][r]
        if claims["normal_divide_before_differentiate_mode"]=="unisolated":row+=factor("normal_Bzz_commutator_factor")*sum(weight(r,k)*coeff[0][k]*J[r-k+2] for k in range(1,r+1))
        production=coeff[0][0]**claims["normal_leading_inverse_power"]*row
        need("unisolated_normal_recovery_order_"+str(r+2),zero(production-oracle))

    # Nonconstant interior model: both mixed coefficients vanish at z=0, yet
    # contribute at z=1/2. No differentiated boundary condition is imposed.
    yy=S.Symbol("tangent_coordinate")
    az=z*(1+z+s*yy);bzy=z*(S.Rational(1,3)+z+yy**2+s*z)
    bzz=2+z+z*z+s+yy;ay=S.Rational(2,5)+s*s+z*yy;byy=3+yy+z*z+s
    up=sum((z-S.Rational(1,2)+s+2*yy)**n/S.factorial(n) for n in range(15))
    source=S.diff(up,s,2)-az*S.diff(up,z,s)-ay*S.diff(up,yy,s)-bzz*S.diff(up,z,2)-2*bzy*S.diff(up,z,yy)-byy*S.diff(up,yy,2)
    point={z:S.Rational(1,2),s:0,yy:0}
    need("mixed_coefficients_vanish_only_on_shell",az.subs(z,0)==0 and bzy.subs(z,0)==0 and az.subs(point)!=0 and bzy.subs(point)!=0)
    principal=S.Matrix([[bzz,bzy],[bzy,byy]]).subs(point)
    need("interior_witness_positive_spatial_symbol",principal[0,0]>0 and principal.det()>0)
    # Values are obtained by differentiating the manufactured fields, unlike
    # the arbitrary-jet product oracle above. A mixed time/tangent germ is used.
    unisolated=(S.diff(up,s,2)+factor("normal_Az_factor")*az*S.diff(up,z,s)+factor("normal_Ay_factor")*ay*S.diff(up,yy,s)+factor("normal_Bzy_factor")*bzy*S.diff(up,z,yy)+factor("normal_Byy_factor")*byy*S.diff(up,yy,2)+factor("normal_source_factor")*source)
    for r in range(13):
        p=1 if r<12 else 0;q=1 if r<11 else 0
        remainder=S.diff(bzz*S.diff(up,z,2),z,r,s,p,yy,q)-bzz*S.diff(up,z,r+2,s,p,yy,q)
        recovered=((S.diff(unisolated,z,r,s,p,yy,q)+factor("normal_Bzz_commutator_factor")*remainder)/bzz).subs(point)
        need("nonconstant_interior_normal_mixed_jet_"+str(r+2),zero(recovered-S.diff(up,z,r+2,s,p,yy,q).subs(point)))

    # Triangularity is an integer multiindex check, not a PDE regularity proof.
    rows=0;commutators=0
    for j in range(2,15):
        for p in range(15-j):
            for q in range(15-j-p):
                alpha=(j-2,p,q);N=j+p+q;rows+=1
                need("normal_triangular_leading_"+str((j,p,q)),all(count<j for count in (j-2,j-1,j-2,j-1,j-2)))
                weights=0
                for beta in itertools.product(*(range(a+1) for a in alpha)):
                    weights+=math.prod(math.comb(a,b) for a,b in zip(alpha,beta))
                    if any(beta):
                        commutators+=1
                        if not N-sum(beta)<=min(13,N-1):raise ValueError("nontriangular coefficient commutator")
                if weights!=2**sum(alpha):raise ValueError("bad multiindex product weights")
    need("normal_multiindex_complete_roster",rows==455)

    # Noncommuting ordered derivations. Oracle repeatedly differentiates a
    # product; production chooses ordered subsets without sorting their words.
    for n in range(1,14):
        word=tuple((i*i+1)%3 for i in range(n))
        literal={((),()):1}
        for letter in reversed(word):
            nxt={}
            for (cw,uw),value in literal.items():
                for key in (((letter,)+cw,uw),(cw,(letter,)+uw)):nxt[key]=nxt.get(key,0)+value
            literal=nxt
        literal.pop(((),word))
        proposed={}
        for mask in range(1,1<<n):
            cw=tuple(word[i] for i in range(n) if mask>>i&1);uw=tuple(word[i] for i in range(n) if not(mask>>i&1))
            if claims["commutator_word_mode"]=="sort":cw=tuple(sorted(cw));uw=tuple(sorted(uw))
            elif claims["commutator_word_mode"]!="preserve":raise ValueError("unsupported commutator word mode")
            key=(cw,uw);proposed[key]=proposed.get(key,0)+factor("commutator_sign")
        need("ordered_positive_coefficient_commutator_order_"+str(n),literal==proposed)

    # Picard difference, retaining nonzero residuals to avoid assuming an
    # auxiliary polynomial solves the actual nonlinear PDE.
    u=1+s+z+s*z+s*s+z**3;v=2+s*s*z+z*z+s**3;w=3+z+s*z*z+s**2
    aa=lambda value:1+value+value**2/7
    bb=lambda value:2+value**2/5
    ff=lambda value:value+value**3+(1+value)*S.diff(value,z)+(2+value)*S.diff(value,s)*S.diff(value,z)
    def P(predecessor,value):return S.diff(value,s,2)-aa(predecessor)*S.diff(value,z,s)-bb(predecessor)*S.diff(value,z,2)
    residual_u=P(v,u)-ff(v);residual_v=P(w,v)-ff(w)
    direct=P(v,u-v)-residual_u+residual_v
    targets={"predecessor":v,"successor":u,"older":w}
    if claims["difference_second_derivative_target"] not in targets:raise ValueError("bad difference derivative target")
    target=targets[claims["difference_second_derivative_target"]]
    proposed=factor("difference_source_factor")*(ff(v)-ff(w))+factor("difference_A_factor")*(aa(v)-aa(w))*S.diff(target,z,s)+factor("difference_B_factor")*(bb(v)-bb(w))*S.diff(target,z,2)
    need("Picard_difference_has_predecessor_second_derivatives",zero(direct-proposed))
    gv,gw,pv,qv,pw,qw=S.symbols("Gv Gw Pv Qv Pw Qw")
    split=(gv-gw)*pw*qw+factor("quadratic_source_cross_factor")*gv*((pv-pw)*qv+pw*(qv-qw))
    need("quadratic_source_difference_has_both_cross_terms",zero(gv*pv*qv-gw*pw*qw-split))

    energy_order=claims["energy_source_order"];normal_order=claims["normal_source_order"]
    need("energy_requires_full_source_order13",energy_order==13 and S.diff(S.diff(s**14/S.factorial(14),s),s,energy_order).subs(s,0)==1)
    need("normal_requires_source_order12",normal_order==12 and S.diff(S.diff(s**13/S.factorial(13),s),s,normal_order).subs(s,0)==1)
    allocations=0
    for derivative_order in range(14):
        for a in range(derivative_order+1):
            for b in range(derivative_order-a+1):
                c=derivative_order-a-b;orders=(a,b+1,c+1);allocations+=1
                if sum(orders)>15 or max(orders)>14 or sum(order>=8 for order in orders)>1:raise ValueError("invalid full source order roster")
    need("full_quadratic_source_order_roster",allocations==560)
    return {"schema":"hdblast-evolution-oracle-report-v1","status":"pass","claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
      "enumeration":{"normal_rows":rows,"normal_coefficient_terms":commutators,"quadratic_source_allocations":allocations},
      "scope":"Finite exact normalized energy/flux, unisolated normal product, ordered-word and source-difference witnesses. Auxiliary four-direction connection and polynomial germs support algebraic identities; no continuum estimate, boundary regularity theorem, actual coefficient ceiling, Picard ball, existence or evolved family is certified."}


def controls(claims):
    records=[]
    for name,value in MUTATIONS:
        bad=copy.deepcopy(claims);bad[name]=value;payload=canonical(bad)
        try:evaluate(json.loads(payload),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name,"self_rehashed_sha256":hashlib.sha256(payload).hexdigest(),"diagnostic":str(exc),"rejected_by":"fresh exact evolution oracle"})
        else:raise ValueError("surviving evolution fault: "+name)
    return {"schema":"hdblast-evolution-oracle-controls-v1","status":"pass","controls":records}
