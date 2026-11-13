#!/usr/bin/env python3
"""Exact nonzero-state canonical/physical and 5D relative-harmonic witnesses.

The gauge oracle solves all five unknown temporal metric velocities directly
from g^BC(Gamma[g]-GammaBar)=0. Production uses the source completion formula.
Finite algebra witnesses do not establish any uniform analytic Cauchy bound.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import sympy as S
import mpmath

e,f=S.symbols("e f")
FACTORS=("density_first_trace_factor","density_mixed_trace_product_factor","density_mixed_trace_cross_factor",
 "K_conversion_factor","K_trace_factor","normal_scalar_conversion_factor","spatial_metric_velocity_factor",
 "spatial_velocity_lapse_cross_factor","scalar_velocity_factor","scalar_velocity_lapse_cross_factor",
 "lapse_velocity_factor","temporal_metric_velocity_factor","mixed_reference_factor",
 "mixed_lapse_linear_factor","mixed_lapse_squared_factor","mixed_lapse_gradient_factor","reference_source_variation_factor")
F={name:S.Symbol("factor_"+name) for name in FACTORS}
MUTATIONS=tuple((name,"1" if name=="K_conversion_factor" else str(S.Rational(value)+1)) for name,value in (
 ("density_first_trace_factor","-1/2"),("density_mixed_trace_product_factor","1/4"),
 ("density_mixed_trace_cross_factor","1/2"),("K_conversion_factor","2"),("K_trace_factor","-1/3"),
 ("normal_scalar_conversion_factor","1"),("spatial_metric_velocity_factor","2"),
 ("spatial_velocity_lapse_cross_factor","1"),("scalar_velocity_factor","1"),
 ("scalar_velocity_lapse_cross_factor","1"),("lapse_velocity_factor","1"),
 ("temporal_metric_velocity_factor","-2"),("mixed_reference_factor","1"),
 ("mixed_lapse_linear_factor","2"),("mixed_lapse_squared_factor","1"),("mixed_lapse_gradient_factor","1")))+(
 ("reference_source_variation_factor","0"),)


def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+"\n").encode()
def clean(v):return S.cancel(S.expand(v))
def trunc(v):
    p=S.Poly(S.expand(v),e,f)
    return S.Add(*(a*e**i*f**j for (i,j),a in p.terms() if i<=1 and j<=1))
def jet(v):
    zero={e:0,f:0}
    return (clean(v.subs(zero))+e*clean(S.diff(v,e).subs(zero))+f*clean(S.diff(v,f).subs(zero))+
            e*f*clean(S.diff(v,e,f).subs(zero)))
def coeff(v,i,j):return S.Poly(S.expand(v),e,f).coeff_monomial(e**i*f**j)
def matrix_trunc(M):return M.applyfunc(trunc)


def build_checks():
    x=S.symbols("x0:4");t,y,z,w=x;at={q:0 for q in x};n=4
    ref=S.diag(4,1+t+t*t,1,1)
    L=S.eye(n);L[0,1]=S.Rational(1,5);L[2,3]=S.Rational(1,7)
    state=(L.T*ref*L+t*S.Matrix([[0,1,0,1],[1,0,1,0],[0,1,1,0],[1,0,0,-1]])+
        y*S.Matrix([[1,0,1,0],[0,-1,0,1],[1,0,0,1],[0,1,1,2]]))
    h=S.Matrix([[1+t,1+y,2+z,w],[1+y,2+z,t,1+w],[2+z,t,-1+w,1+y],[w,1+w,1+y,1+t]])
    k=S.Matrix([[2+y,t,1+w,1+z],[t,-1+w,1+y,z],[1+w,1+y,1+t,2+w],[1+z,z,2+w,-2+y]])
    spatial=state+e*h+f*k;G=state.subs(at);H=h.subs(at);K=k.subs(at);ref0=ref.subs(at);refinv=ref0.inv()
    gp=spatial.subs(at);Aexact=gp.inv();A=Aexact.applyfunc(jet);A0=G.inv()
    determinant=gp.det();detref=ref0.det()
    half_exact=S.sqrt(detref/determinant);half=jet(half_exact);half_inverse=jet(S.sqrt(determinant/detref))
    th=S.trace(A0*H);tk=S.trace(A0*K);cross_trace=S.trace(A0*H*A0*K)
    density=(1+F["density_first_trace_factor"]*(e*th+f*tk)+e*f*(F["density_mixed_trace_product_factor"]*th*tk+
            F["density_mixed_trace_cross_factor"]*cross_trace))
    tau0=S.Matrix([[2,1,0,1],[1,-1,2,0],[0,2,3,1],[1,0,1,1]])
    tauh=S.Matrix([[1,2,1,0],[2,-1,0,1],[1,0,2,1],[0,1,1,-2]])
    tauk=S.Matrix([[2,0,1,1],[0,1,2,0],[1,2,-1,1],[1,0,1,3]])
    tau=tau0+e*tauh+f*tauk;omega=3+2*e-f
    trace=S.trace(gp*tau)
    Sphysical=matrix_trunc(gp*tau*gp-trace*gp/3)
    Ktrue=matrix_trunc(2*half*Sphysical)
    Pitrue=trunc(half*omega)
    Kprod=matrix_trunc(F["K_conversion_factor"]*density*(gp*tau*gp+F["K_trace_factor"]*trace*gp))
    Piprod=trunc(F["normal_scalar_conversion_factor"]*density*omega)
    forward_tau=matrix_trunc(half_inverse*(A*Kprod*A-S.trace(A*Kprod)*A)/2)
    forward_omega=trunc(half_inverse*Piprod)
    checks=[("J_inverse_half_mixed_density_jet",[density-half]),
      ("full_half_action_forward_canonical_metric_density",list(forward_tau-tau)),
      ("full_forward_canonical_scalar_density",[forward_omega-omega])]
    # The selected lapse is affine in the actual parameter: both directions
    # have nonzero parameter components, so its spatial derivative and all
    # lapse cross terms participate. The expansion state itself is nonzero.
    rho=5+t+y*y;profile=S.Rational(1,11)+t+2*z+y*w
    parameter=S.Rational(1,7)+2*e-3*f
    ell=parameter*profile
    alpha=rho*(1+ell)
    rhop=rho.subs(at);drho=S.Matrix([S.diff(rho,q).subs(at) for q in x])
    ap=alpha.subs(at);dalpha=S.Matrix([S.diff(alpha,q).subs(at) for q in x]);ellp=ell.subs(at)
    dell=S.Matrix([S.diff(ell,q).subs(at) for q in x])
    gamma_t_true=matrix_trunc(2*ap*Ktrue)
    phi_t_true=trunc(ap*Pitrue)
    gamma_t_prod=matrix_trunc(F["spatial_metric_velocity_factor"]*ap*Kprod)
    lapsecross=2*(coeff(ap,1,0)*Kprod.applyfunc(lambda a:coeff(a,0,1))+
                   coeff(ap,0,1)*Kprod.applyfunc(lambda a:coeff(a,1,0)))
    gamma_t_prod+=e*f*(F["spatial_velocity_lapse_cross_factor"]-1)*lapsecross
    phi_t_prod=trunc(F["scalar_velocity_factor"]*ap*Piprod)+(F["scalar_velocity_lapse_cross_factor"]-1)*e*f*(
        coeff(ap,1,0)*coeff(Piprod,0,1)+coeff(ap,0,1)*coeff(Piprod,1,0))
    alpha_t=trunc(F["lapse_velocity_factor"]*ap*ap*S.trace(A*Kprod))
    gtt_t=trunc(F["temporal_metric_velocity_factor"]*ap*alpha_t)
    checks.extend([("spatial_velocity_full_lapse_cross_terms",list(gamma_t_prod-gamma_t_true)),
                   ("scalar_velocity_full_lapse_cross_terms",[phi_t_prod-phi_t_true])])
    # Fixed reference spacetime connection from (-rho^2 dt^2+ref) at t=0.
    ref_full=S.diag(-rho*rho,1,1,1,1)
    ref_full[1:5,1:5]=ref
    ref_full_inverse=ref_full.inv()
    ref_derivatives=[S.zeros(5)]+[ref_full.diff(q) for q in x]
    bar=[[[clean(sum(ref_full_inverse[a,d]*(ref_derivatives[b][d,c]+ref_derivatives[c][b,d]-
                     ref_derivatives[d][b,c]) for d in range(5))/2).subs(at)
           for c in range(5)] for b in range(5)] for a in range(5)]
    # Twenty unknown coefficients encode all five temporal metric velocities
    # through the mixed e*f term. No completion formula enters the solve.
    unknowns=S.symbols("velocity0:20")
    velocity=[sum(unknowns[4*a+j]*mon for j,mon in enumerate((1,e,f,e*f))) for a in range(5)]
    g5=S.zeros(5);g5[0,0]=trunc(-ap*ap);g5[1:5,1:5]=gp
    inv5=S.zeros(5);inv5[0,0]=jet(-1/(ap*ap));inv5[1:5,1:5]=A
    derivative=[S.zeros(5) for _ in range(5)]
    derivative[0][0,0]=velocity[0]
    for i in range(n):derivative[0][0,i+1]=derivative[0][i+1,0]=velocity[i+1]
    derivative[0][1:5,1:5]=gamma_t_true
    for r in range(n):
        derivative[r+1][0,0]=trunc(-2*ap*dalpha[r])
        derivative[r+1][1:5,1:5]=spatial.diff(x[r]).subs(at)
    christoffel=[[[trunc(sum(inv5[a,d]*(derivative[b][d,c]+derivative[c][b,d]-derivative[d][b,c]) for d in range(5))/2)
                   for c in range(5)] for b in range(5)] for a in range(5)]
    gauges=[trunc(sum(inv5[b,c]*(christoffel[a][b][c]-bar[a][b][c]) for b in range(5) for c in range(5))) for a in range(5)]
    equations=[coeff(v,i,j) for v in gauges for i,j in ((0,0),(1,0),(0,1),(1,1))]
    matrix,rhs=S.linear_eq_to_matrix(equations,unknowns)
    solution=matrix.inv()*rhs
    oracle_velocity=[sum(solution[4*a+j]*mon for j,mon in enumerate((1,e,f,e*f))) for a in range(5)]
    source=[trunc(sum(inv5[b,c]*bar[a][b][c] for b in range(5) for c in range(5))) for a in range(5)]
    # Background-subtracted lower shift formula from the actual gauge source.
    bar_spatial=[[[bar[a+1][b+1][c+1] for c in range(n)] for b in range(n)] for a in range(n)]
    dspatial=[spatial.diff(q).subs(at) for q in x]
    nabla=[S.Matrix(n,n,lambda i,j:dspatial[r][i,j]-sum(bar_spatial[l][r][i]*gp[l,j]+
                    bar_spatial[l][r][j]*gp[i,l] for l in range(n))) for r in range(n)]
    C=[trunc(sum(A[k,l]*(nabla[k][l,i]-nabla[i][k,l]/2) for k in range(n) for l in range(n))) for i in range(n)]
    hstate=gp-ref0
    reference_lower=rhop*hstate*refinv*drho
    mixed=[trunc(ap*ap*C[i]+F["mixed_reference_factor"]*reference_lower[i]-rhop*(F["mixed_lapse_linear_factor"]*ellp+
          F["mixed_lapse_squared_factor"]*ellp*ellp)*drho[i]-F["mixed_lapse_gradient_factor"]*rhop*rhop*(1+ellp)*dell[i]) for i in range(n)]
    source_variation=S.Matrix([v-v.subs({e:0,f:0}) for v in source[1:]])
    frozen_source_error=matrix_trunc(ap*ap*gp*source_variation)
    mixed=[trunc(mixed[i]+(1-F["reference_source_variation_factor"])*frozen_source_error[i]) for i in range(n)]
    checks.extend([("five_dimensional_gauge_temporal_velocity",[gtt_t-oracle_velocity[0]]),
       ("five_dimensional_gauge_mixed_velocities",[mixed[i]-oracle_velocity[i+1] for i in range(n)]),
       ("canonical_trace_temporal_velocity_sign",[gtt_t-trunc(S.Rational(4,3)*ap**3*half*trace)])])
    if source[0]!=0:raise ValueError("initial reference R0 is not zero")
    if not any(v-v.subs({e:0,f:0})!=0 for v in source[1:]):raise ValueError("metric-dependent reference source witness degenerate")
    if not any(v!=0 for v in bar_spatial):raise ValueError("curved reference witness degenerate")
    if G==ref0:raise ValueError("nonzero metric state witness degenerate")
    metadata={"spatial_dimension":4,"spacetime_dimension":5,"gauge_unknown_coefficients":20,
       "actual_half_action_canonical_normalization":True,"nonzero_metric_momentum_scalar_momentum_and_lapse_state":True,
       "reference_R0_identically_zero":True,"reference_Ri_metric_dependence_nonzero":True,
       "gauge_oracle":"solve all 5D contracted Christoffel gauge equations for g_tt,t and g_ti,t, at all four parameter coefficients",
       "scope":"Finite exact rational two-direction witnesses at an auxiliary nonzero curved initial state. No selected-family amplitude, uniform Cauchy norm, nonlinear evolution or observer conclusion is proved."}
    return checks,metadata


def evaluate(claims,check):
    if S.__version__!="1.14.0" or mpmath.__version__!="1.3.0":raise ValueError("requires pinned exact engine")
    if claims.get("schema")!="hdblast-uniform-cauchy-algebra-claims-v1" or set(claims)!={"schema",*FACTORS}:raise ValueError("bad Cauchy algebra schema")
    checks,metadata=build_checks();replacement={F[k]:S.Rational(claims[k]) for k in FACTORS}
    records=[];bad=[]
    for name,residuals in checks:
        nonzero=[clean(a.subs(replacement)) for a in residuals];nonzero=[a for a in nonzero if a!=0]
        check("cauchy_algebra/"+name,not nonzero)
        records.append({"name":name,"components":len(residuals),"status":"fail" if nonzero else "pass"})
        if nonzero:bad.append({"name":name,"exact_nonzero_residual":str(nonzero[0])})
    if bad:raise ValueError("Cauchy oracle rejects claims: "+json.dumps(bad,sort_keys=True))
    return {"schema":"hdblast-uniform-cauchy-algebra-v1","status":"pass","claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,"metadata":metadata}


def controls(claims):
    records=[]
    for name,value in MUTATIONS:
        bad=copy.deepcopy(claims);bad[name]=value;payload=canonical(bad)
        try:evaluate(json.loads(payload),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name,"self_rehashed_sha256":hashlib.sha256(payload).hexdigest(),
           "rejected_by":"fresh independent canonical/5D gauge oracle","diagnostic":str(exc)})
        else:raise ValueError("surviving Cauchy factor fault: "+name)
    return {"schema":"hdblast-uniform-cauchy-controls-v1","status":"pass","controls":records}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("cauchy_claims.json"))
    parser.add_argument("--report",type=Path);parser.add_argument("--controls",type=Path)
    args=parser.parse_args(argv)
    try:
        c=json.loads(args.claims.read_text());r=evaluate(c,lambda name,ok:None);payload=canonical(r)
        if args.report:args.report.write_bytes(payload)
        if args.controls:args.controls.write_bytes(canonical(controls(c)))
        sys.stdout.buffer.write(payload);return 0
    except (ValueError,KeyError,TypeError,OSError) as exc:sys.stderr.write(str(exc)+"\n");return 1
