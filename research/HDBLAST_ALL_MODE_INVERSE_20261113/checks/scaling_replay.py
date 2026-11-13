#!/usr/bin/env python3
"""Finite exact symbol/scaling witnesses; no all-mode inverse certificate.

The weighted norm uses q^j on the jth derivative, not q^(j/2).  Production
claims are compared with exact tensors or direct conjugated derivatives.
These local and endpoint calculations do not prove global coercivity,
regularity, projection, surjectivity, or an inverse bound.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
from pathlib import Path
import sys

import mpmath
import sympy as s


def canonical(value):
    return (json.dumps(value,sort_keys=True,indent=2)+"\n").encode()


def clean(expression):
    return s.cancel(s.expand(expression))


def build_witnesses(claims):
    out=[]
    def equality(name,left,right):
        out.append((name,clean(left-right)))
    def condition(name,value):
        out.append((name,s.Integer(0) if value else s.Integer(1)))
    condition("schema",claims["schema"]=="hdblast-all-mode-scaling-claims-v1")
    condition("sympy_pin",s.__version__=="1.14.0")
    condition("mpmath_pin",mpmath.__version__=="1.3.0")
    condition("physical_dimension",claims["dimension"]==4)
    alpha=s.Rational(claims["norm_derivative_q_power"])
    # An operator with coefficient q^m and differential order m maps the
    # highest H^(k+m) derivative to the highest H^k derivative precisely when
    # alpha*(k+m)=m+alpha*k.  q can be arbitrarily small at the endpoint.
    for sector,k,m in (("scalar",12,4),("momentum",12,2)):
        order=claims[sector+"_order"]
        coefficient_power=claims[sector+"_multiplier_q_power"]
        equality(sector+"_actual_order",s.Integer(order),s.Integer(m))
        equality(sector+"_principal_q_scaling",alpha*(k+order),coefficient_power+alpha*k)
    equality("actual_norm_q_derivative_power",alpha,s.Integer(1))

    # Independent tensor principal symbols in an orthonormal physical frame.
    # A common complex radial covector includes the endpoint weight drift.
    xi,tau,c=s.symbols("xi tau c",real=True)
    cov=s.Matrix([-c+s.I*xi,s.I*tau,0,0])
    total=xi*xi+tau*tau+c*c
    tensor=(cov*cov.T-s.eye(4)*(cov.T*cov)[0])/2
    scalar_oracle=clean(sum(s.conjugate(a)*a for a in tensor))
    scalar_production=(s.Rational(claims["scalar_joint_symbol_sum_squared_factor"])*total**2+
                       s.Rational(claims["scalar_joint_symbol_c2tau2_factor"])*c*c*tau*tau)
    equality("scalar_tensor_vs_claimed_joint_symbol",scalar_production,scalar_oracle)
    equality("scalar_principal_normalization",s.Rational(claims["scalar_adjoint_principal_factor"]),s.Rational(1,2))
    # Separately derive the symbol from the retained scalar harmonic rows at
    # the joint endpoint/high-angular scaling q sqrt(lambda)/rho=tau.
    rowA=tau*tau/2
    rowB=-(s.I*xi-c)**2/2+tau*tau/3
    rowC=tau*(s.I*xi-c)/s.sqrt(2)
    rowT=tau*tau/s.sqrt(6)
    row_oracle=clean(abs(rowA)**2+3*abs(rowB)**2+abs(rowC)**2+abs(rowT)**2)
    equality("scalar_harmonic_rows_vs_tensor",row_oracle,scalar_oracle)
    # An explicit sum of squares proves the displayed algebraic lower factor;
    # it does not establish an estimate for a variable coefficient operator.
    SOS=((c*c-xi*xi-tau*tau)**2+4*c*c*xi*xi)/2
    equality("scalar_joint_symbol_sum_of_squares",scalar_oracle-total**2/4,SOS)
    equality("scalar_joint_symbol_claimed_lower_factor",s.Rational(claims["scalar_joint_symbol_lower_factor"]),s.Rational(1,4))
    equality("scalar_joint_lower_sharp_point",scalar_oracle.subs({xi:0,tau:c}),c**4)

    # Build the 10 independent symmetric tensor rows with Frobenius weighting
    # (off-diagonal entries count twice).  The Gram matrix is independently
    # checked against its rank-one formula and characteristic polynomial.
    basis=[s.eye(4)[:,i] for i in range(4)]
    maps=[cov*z.T+z*cov.T for z in basis]
    gram=s.Matrix(4,4,lambda i,j:clean(sum(s.conjugate(maps[i][a,b])*maps[j][a,b]
                                                       for a in range(4) for b in range(4))))
    equality("momentum_principal_normalization",s.Rational(claims["momentum_adjoint_principal_factor"]),s.Integer(1))
    claimed=s.Rational(claims["momentum_transverse_factor"])*total*s.eye(4)+(
        s.Rational(claims["momentum_longitudinal_factor"])-
        s.Rational(claims["momentum_transverse_factor"]))*cov*cov.conjugate().T
    for i in range(4):
        for j in range(4):
            equality("momentum_full_gram_component_"+str(i)+str(j),claimed[i,j],gram[i,j])
    lam=s.Symbol("spectral_parameter",real=True)
    full_polynomial=gram.charpoly()
    equality("momentum_full_characteristic_polynomial",full_polynomial.as_expr().subs(full_polynomial.gen,lam),
             (lam-2*total)**3*(lam-4*total))
    # Retained poloidal harmonic rows supply an independent reduced block.
    D=s.I*xi-c
    poloidal=s.Matrix([[2*D,0],[s.sqrt(2)*tau,s.sqrt(2)*D],
                      [0,-s.sqrt(12)*tau/3],[0,2*s.sqrt(s.Rational(2,3))*tau]])
    reduced=(poloidal.conjugate().T*poloidal).applyfunc(clean)
    equality("momentum_poloidal_trace",s.trace(reduced),6*total)
    equality("momentum_poloidal_determinant",reduced.det(),8*total**2)
    reduced_polynomial=reduced.charpoly()
    equality("momentum_poloidal_characteristic_polynomial",reduced_polynomial.as_expr().subs(reduced_polynomial.gen,lam),
             (lam-2*total)*(lam-4*total))

    # Direct conjugation oracle for a model exact endpoint germ.  Its r_y is
    # the nonzero endpoint limit of the actual defining profile.  Differentiate
    # the actual exponential test potential, rather than substituting a drift
    # into an operator expansion by hand.
    r,kappa,ry=s.symbols("r kappa r_y",positive=True)
    ell=s.Rational(claims["ell"]); kap=s.Rational(claims["kappa"])
    equality("retained_kappa",kap,ell/100)
    q=kappa*kappa*r*r
    M=kappa*r*s.exp(-1/r)
    rho_germ=s.Symbol("rho_germ",positive=True)
    background_a=s.Symbol("background_a",real=True)
    M_full=M*rho_germ**s.Rational(3,2)
    actual_drift=clean(q*(ry*s.diff(M_full,r)+background_a*rho_germ*s.diff(M_full,rho_germ))/M_full)
    equality("actual_unitary_drift_including_radial_volume",actual_drift,
             kappa*kappa*(1+r)*ry+s.Rational(3,2)*q*background_a)
    equality("actual_unitary_drift_endpoint_limit",s.limit(actual_drift,r,0,dir="+"),kappa*kappa*ry)
    s_coord=-1/(kappa*kappa*ry*r)
    potential=s.exp(s.I*xi*s_coord)/M
    one=clean(q*ry*s.diff(potential,r)/potential)
    two=clean(q*q*ry*ry*s.diff(potential,r,2)/potential)
    equality("direct_endpoint_first_conjugated_derivative",s.limit(one,r,0,dir="+"),s.I*xi-kappa*kappa*ry)
    equality("direct_endpoint_second_conjugated_derivative",s.limit(two,r,0,dir="+"),(s.I*xi-kappa*kappa*ry)**2)
    # Exact defining profile in the original z chart, retaining dy=-rho dz.
    z,rho=s.symbols("z rho",positive=True)
    profile=(z-ell/2)*(5*ell/2-z)/ell**2
    derivative_y=-s.diff(profile,z)/rho
    lower=clean(kap*kap*derivative_y.subs(z,ell/2))
    upper=clean(kap*kap*derivative_y.subs(z,5*ell/2))
    equality("lower_endpoint_drift",lower,-2*kap*kap/(ell*rho))
    equality("upper_endpoint_drift",upper,2*kap*kap/(ell*rho))
    condition("endpoint_drift_nonzero",lower!=0 and upper!=0)
    equality("endpoint_drift_magnitude_lower",upper.subs(rho,131),s.Rational(1,655000000000))
    equality("endpoint_drift_magnitude_upper",upper.subs(rho,129),s.Rational(1,645000000000))
    # Reflection averaging of the ten already-classified adjoint kernels.
    # This does not independently classify the kernels of the differential
    # operators; it checks their stated invariant projection exactly.
    group=list(itertools.product((-1,1),repeat=4))
    for i in range(4):
        equality("boost_kernel_reflection_average_"+str(i),s.Rational(sum(g[i] for g in group),16),0)
    for i in range(4):
        for j in range(i+1,4):
            equality("rotation_kernel_reflection_average_"+str(i)+str(j),s.Rational(sum(g[i]*g[j] for g in group),16),0)
    angular=s.Symbol("angular_eigenvalue",real=True)
    degree=s.Symbol("harmonic_degree",real=True)
    equality("scalar_eigenvalue_gap_factorization",degree*(degree+2)-8,(degree-2)*(degree+4))
    equality("scalar_uniform_lambda_ratio_margin",s.Rational(8,5)-angular/(angular-3),
             3*(angular-8)/(5*(angular-3)))
    equality("momentum_poloidal_uniform_kv_margin",s.Rational(3,5)-3*angular/(8*(angular-3)),
             9*(angular-8)/(40*(angular-3)))
    metadata={
        "actual_derivative_multiplier":"q^j; norm integrand w^2 q^(2j)|nabla^j u|^2",
        "scalar_joint_symbol":str(scalar_oracle),
        "scalar_symbol_lower":"(xi^2+tau^2+c^2)^2/4",
        "momentum_eigenvalues":["2*(xi^2+tau^2+c^2) [multiplicity 3]","4*(xi^2+tau^2+c^2) [multiplicity 1]"],
        "endpoint_drift_magnitude_lower":"1/655000000000",
        "endpoint_drift_magnitude_upper":"1/645000000000",
        "conclusion":"These checks find no principal scaling or frozen endpoint-symbol obstruction in the stated same-(q,w) norms.",
        "scope":"Finite exact algebra and local/endpoint symbol models only; no global all-mode inverse constant, boundary coercivity proof, regularity proof, source authentication, or nonlinear amplitude follows.",
    }
    return out,metadata


def evaluate(claims,check):
    witnesses,metadata=build_witnesses(claims)
    records=[];failures=[]
    for name,residual in witnesses:
        ok=residual==0
        check("scaling/"+name,ok)
        records.append({"name":name,"status":"pass" if ok else "fail"})
        if not ok:
            failures.append({"name":name,"exact_residual":str(residual)})
    if failures:
        raise ValueError("scaling oracle rejects claims: "+json.dumps(failures,sort_keys=True))
    return {"schema":"hdblast-all-mode-scaling-replay-v1","status":"pass",
            "claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),
            "engine":"sympy==1.14.0; mpmath==1.3.0","checks":records,"metadata":metadata}


def controls(claims):
    records=[]
    for name,value in (("norm_derivative_q_power", "1/2"),
                       ("scalar_multiplier_q_power",8),
                       ("momentum_multiplier_q_power",4),
                       ("scalar_joint_symbol_sum_squared_factor","1"),
                       ("scalar_joint_symbol_c2tau2_factor","0"),
                       ("momentum_transverse_factor","1"),
                       ("momentum_longitudinal_factor","2")):
        bad=copy.deepcopy(claims);bad[name]=value
        encoded=canonical(bad);rehashed=hashlib.sha256(encoded).hexdigest()
        try:
            evaluate(json.loads(encoded),lambda name,ok:None)
        except ValueError as exc:
            records.append({"fault":name,"rehash":rehashed,"rejected_by":"recomputed exact algebra oracle",
                            "diagnostic":str(exc)})
        else:
            raise ValueError("surviving scaling fault: "+name)
    return {"schema":"hdblast-all-mode-scaling-controls-v1","status":"pass","controls":records}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("scaling_claims.json"))
    parser.add_argument("--report",type=Path)
    parser.add_argument("--controls",type=Path)
    args=parser.parse_args(argv)
    try:
        claims=json.loads(args.claims.read_text())
        report=evaluate(claims,lambda name,ok:None)
        encoded=canonical(report)
        if args.report:args.report.write_bytes(encoded)
        if args.controls:args.controls.write_bytes(canonical(controls(claims)))
        sys.stdout.buffer.write(encoded)
        return 0
    except (ValueError,KeyError,TypeError,OSError) as exc:
        sys.stderr.write(str(exc)+"\n");return 1


if __name__=="__main__":
    raise SystemExit(main())
