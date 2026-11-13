#!/usr/bin/env python3
"""Exact rational accounting for analytic inequalities supplied by the notes.

The ledger verifies arithmetic and formal ceiling comparisons.  It does not
prove the analytic norm estimates used as its inputs.  In particular B1 and B2
are symbolic 10^a exp(10^b) ceilings and are never numerically exponentiated.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys


SCHEMA="hdblast-proper-second-geometry-exact-ledger-claims-v1"


def canonical(value):
    return (json.dumps(value,sort_keys=True,indent=2)+"\n").encode()


def evaluate(claims,check):
    records=[]

    def require(name,condition):
        boolean=bool(condition)
        records.append({"name":name,"status":"pass" if boolean else "fail"})
        check("exact_ledger/"+name,boolean)
        if not boolean:
            raise ValueError("exact ledger rejects "+name)

    require("schema",claims["schema"] == SCHEMA)
    require("top_level_keys",set(claims)=={"schema","inputs","derived","ceilings"})
    data=claims["inputs"]
    require("input_keys",set(data)=={"ell","rho_upper","embedding_constant","frame_multiplier",
                                        "P1_multiplier","P2_multiplier","B1","B2"})
    ell=Q(data["ell"]); M=Q(data["embedding_constant"])
    frame=Q(data["frame_multiplier"]); rho=Q(data["rho_upper"])
    P1=Q(data["P1_multiplier"]); P2=Q(data["P2_multiplier"])
    require("retained_ell",ell==Q(1,10**6))
    require("retained_radius_upper",rho==131)
    require("retained_embedding_constant",M==10**7)
    require("frame_multiplier",frame==512)
    require("connection_conversion_ceilings",P1==200 and P2==20)
    T=13*ell/8; q=25*ell/16; qD=ell/4; H=frame*M
    cosh_square=(1-q*q/2)**-2
    CR=cosh_square*(q-qD)**2/2
    cov2=1+5+2*4+6
    cov3=1+9+7*4+26+2+20*4+24
    E1=4*Q(3,2)+3
    E2linear=E1
    E2quadratic=4*3+4*3+2*3*Q(3,2)**2
    CrH2=E2linear*P2
    CrH1square=E2quadratic*P1**2+2*E1*P1*q*q+4*E1*P1*3*q
    derived={
        "T":T,"Q":q,"q_D":qD,"coordinate_margin":T-q,
        "H_multiplier":H,"physical_metric_multiplier":rho*rho*H,
        "Cxi_multiplier":q*q*H,"Ceta_multiplier":3*q*H,
        "coordinate_velocity_multiplier":2*q*H,
        "scalar_quadratic_multiplier":3*M*q*q*H,
        "nabla_second_multiplier":cov2,"nabla_third_multiplier":cov3,
        "E1_multiplier":E1,"E2_linear_multiplier":E2linear,
        "E2_quadratic_multiplier":E2quadratic,
        "Cr_H2_multiplier":CrH2,"Cr_H1_squared_multiplier":CrH1square,"C_R":CR,
    }
    require("derived_keys",set(claims["derived"])==set(derived))
    for name,value in derived.items():
        require("derived/"+name,Q(claims["derived"][name])==value)
    require("strict_coordinate_margin",T-q>0)
    require("retained_low_time_region",T<Q(1,100))
    require("first_force_rounding",30*q<1)
    require("eta_rounding",q*q<1)
    require("physical_metric_rounding",rho*rho*H<10**14)
    require("nabla_conversion_rounding",cov2<=P2 and cov3<P1)
    require("curvature_coefficient_rounding",CrH1square<1500001)
    require("curvature_B2_coefficient_rounding",CrH2*H+1500001*H*H<4*10**25)
    require("kernel_quiet_limit",qD==ell/4 and qD<ell)
    require("kernel_rounding",CR<Q(8614,10**16))

    # Formal exponential ceilings: no exp(), no 10**4000000 or larger integer.
    B1=data["B1"]; B2=data["B2"]
    for name,B in (("B1",B1),("B2",B2)):
        require(name+"_shape",set(B)=={"ten_exponent","exp_ten_exponent"} and
                all(type(v) is int and v>0 for v in B.values()))
    require("retained_B1",B1=={"ten_exponent":4000000,"exp_ten_exponent":152017})
    require("retained_B2",B2=={"ten_exponent":10200000,"exp_ten_exponent":152019})
    require("B1_squared_prefactor",2*B1["ten_exponent"]<B2["ten_exponent"])
    delta=B2["exp_ten_exponent"]-B1["exp_ten_exponent"]
    require("B1_squared_exponential_shift",0<delta<100)
    require("B1_squared_exponential",2<10**delta)
    require("scalar_B2_rounding",M+derived["scalar_quadratic_multiplier"]<10**8)
    require("curvature_decimal_rounding",4*10**25<10**26)

    # Independently replay the retained first-field pointwise/total-area
    # arithmetic.  Those first-field identities remain analytic source inputs.
    d=q-ell; U=Q(23,2*10**9); Fmax=Q(11,125)
    Bmax=U*(1+q*d+Q(4,3)*d*d)
    jlower=cosh_square*Fmax*Fmax*d*d/3+Q(4,3)*U*U*d**4
    jupper=9*Bmax*Bmax
    Alower=cosh_square*Fmax*Fmax*d*d/18+Q(20,27)*U*U*d**4
    Aupper=Q(3,2)*Bmax*Bmax
    require("area_density_known_lower_rounding",jlower<Q(817,10**18))
    require("area_density_known_upper_rounding",jupper<Q(1191,10**18))
    require("total_area_known_lower_rounding",Alower<Q(1362,10**19))
    require("total_area_known_upper_rounding",Aupper<Q(1985,10**19))
    kernel_round=Q(8614,10**16); Cr_round=4*10**25
    # B2>1 absorbs each small known term into its common B2 coefficient.
    known={"denominator":Q(5357,10**19),"area_density":Q(1191,10**18),
           "total_area":Q(1985,10**19)}
    for name,value in known.items():
        require(name+"_B2_rounding",value+kernel_round*Cr_round<10**14)
    offsets={"scalar":8,"curvature":26,"denominator":14,"area_density":14,"total_area":14}
    require("ceiling_keys",set(claims["ceilings"])==set(offsets))
    for name,offset in offsets.items():
        target={"ten_exponent":B2["ten_exponent"]+offset,"exp_ten_exponent":B2["exp_ten_exponent"]}
        require("formal_ceiling/"+name,claims["ceilings"][name]==target)
    return {"schema":"hdblast-proper-second-geometry-exact-ledger-v1","status":"pass",
            "claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),
            "arithmetic":"standard-library Fraction; formal exponential ceilings, never evaluated",
            "checks":records,"derived":{k:str(v) for k,v in derived.items()},
            "first_field_known_area_budgets":{
                "j2_lower_magnitude":str(jlower),"j2_upper":str(jupper),
                "A2_lower_magnitude":str(Alower),"A2_upper":str(Aupper)},
            "formal_ceilings":claims["ceilings"],
            "scope":"Exact arithmetic conditional on analytic wave, trace, observer, curvature-norm, and Raychaudhuri inputs; no proof of those estimates or finite-amplitude control."}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("numeric_claims.json"))
    parser.add_argument("--report",type=Path)
    args=parser.parse_args(argv)
    try:
        report=evaluate(json.loads(args.claims.read_text()),lambda name,ok:None)
        payload=canonical(report)
        if args.report:
            args.report.write_bytes(payload)
        sys.stdout.buffer.write(payload)
        return 0
    except (ValueError,KeyError,TypeError,OSError) as exc:
        sys.stderr.write("verification failed: "+str(exc)+"\n")
        return 1


if __name__=="__main__":
    raise SystemExit(main())
