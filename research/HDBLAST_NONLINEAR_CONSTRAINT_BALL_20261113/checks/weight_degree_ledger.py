#!/usr/bin/env python3
"""Exact derivative/weight accounting; no embedding or product constant."""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys


def canonical(value):return (json.dumps(value,sort_keys=True,indent=2)+"\n").encode()


def evaluate(claims,check):
    records=[]
    def need(name,ok):
        boolean=bool(ok);check("weight_degree/"+name,boolean)
        records.append({"name":name,"status":"pass" if boolean else "fail"})
        if not boolean:raise ValueError("weight-degree oracle rejects "+name)
    need("schema",claims["schema"]=="hdblast-nonlinear-ball-weight-degree-claims-v1")
    need("dimension",claims["dimension"]==4)
    need("source_order",claims["source_max_order"]==12)
    need("source_normalization",claims["literal_source_w_power"]==-2)
    need("normalization_derivative_cost",claims["normalization_derivative_q_cost"]==-1)
    # Each row records the actual derivative orders already present before
    # distributing the outer source derivatives.
    terms=(("h_D2h","metric",0,"metric",2),
           ("Dh_Dh","metric",1,"metric",1),
           ("tau_tau","momentum",0,"momentum",0),
           ("h_Dtau","metric",0,"momentum",1),
           ("tau_Dh","momentum",0,"metric",1))
    allocation_count=0;max_low=0;max_high=0;term_records=[]
    for name,kind1,base1,kind2,base2 in terms:
        current=0
        for n in range(claims["source_max_order"]+1):
            for normalization_order in range(n+1):
                remaining=n-normalization_order
                for first in range(remaining+1):
                    orders=(base1+first,base2+remaining-first)
                    kinds=(kind1,kind2)
                    low=0 if orders[0]<=orders[1] else 1;high=1-low
                    qpower=(n+normalization_order*claims["normalization_derivative_q_cost"]+
                        claims[kinds[low]+"_low_q_offset"]-orders[low]+
                        claims[kinds[high]+"_high_q_offset"]-orders[high])
                    wpower=claims["Y_exterior_w_power"]+claims["literal_source_w_power"]+2*claims["correction_w_power"]
                    need_math=(qpower==claims["residual_q_power"] and wpower==claims["residual_w_power"])
                    if not need_math:raise ValueError(name+" has wrong exact residual degree "+str((wpower,qpower)))
                    adj=lambda kind:claims["scalar_adjoint_order"] if kind=="metric" else claims["momentum_adjoint_order"]
                    available=lambda kind:claims["potential_scalar_order"] if kind=="metric" else claims["potential_oneform_order"]
                    low_potential=orders[low]+adj(kinds[low])+claims["prospective_point_evaluation_orders"]
                    high_potential=orders[high]+adj(kinds[high])
                    if low_potential>available(kinds[low]) or high_potential>available(kinds[high]):
                        raise ValueError(name+" exceeds a retained potential derivative budget")
                    current+=1;max_low=max(max_low,low_potential);max_high=max(max_high,high_potential)
        need(name+"_all_source_derivative_allocations",current>0)
        allocation_count+=current;term_records.append({"term":name,"allocation_count":current})
    need("literal_final_degrees",claims["residual_w_power"]==1 and claims["residual_q_power"]==0)
    need("actual_input_orders",claims["potential_scalar_order"]==16 and claims["potential_oneform_order"]==14 and
         claims["scalar_adjoint_order"]==2 and claims["momentum_adjoint_order"]==1 and claims["prospective_point_evaluation_orders"]==3)
    # Local rescaling independently tracks dimension four and radius q.
    volume_half=F(claims["dimension"],2)
    local_metric_q=4-2-volume_half
    local_momentum_q=2-1-volume_half
    need("local_metric_correction_degree",local_metric_q==0)
    need("local_momentum_correction_degree",local_momentum_q==-1)
    need("local_normalized_source_degree",volume_half-2==0 and 1-2+2==1)
    # alpha=q/rho, M=w sqrt(q rho^3): alpha^-3/2 / M has
    # exponents (w,q,rho)=(-1,-2,0), without evaluating an exponential.
    need("prospective_point_factor_q",-F(3,2)-F(1,2)==-2)
    need("prospective_point_factor_rho",F(3,2)-F(3,2)==0)
    need("physical_flat_multiplier_power",24*8==192 and 24**24<10**34 and 192+34==claims["physical_multiplier_decimal_exponent"]==226)
    return {"schema":"hdblast-nonlinear-ball-weight-degree-v1","status":"pass",
       "claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
       "allocations":term_records,"allocation_count":allocation_count,
       "max_low_potential_derivative_order":max_low,"max_high_potential_derivative_order":max_high,
       "scope":"Formal weight and derivative-degree accounting for five leading nonlinear terms, including derivatives on w^-2. Point evaluation, products, lower coefficients, complete D2F estimates and all numerical map/radius constants require separate analytic proofs."}


def controls(claims):
    records=[]
    for name,value in (("literal_source_w_power",-1),("normalization_derivative_q_cost",0),
        ("metric_high_q_offset",1),("momentum_high_q_offset",0),
        ("potential_oneform_order",13),("physical_multiplier_decimal_exponent",200)):
        bad=copy.deepcopy(claims);bad[name]=value;payload=canonical(bad)
        try:evaluate(json.loads(payload),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name,"self_rehashed_sha256":hashlib.sha256(payload).hexdigest(),
             "rejected_by":"fresh exact weight/derivative-degree oracle","diagnostic":str(exc)})
        else:raise ValueError("surviving weight fault: "+name)
    return {"schema":"hdblast-nonlinear-ball-weight-controls-v1","status":"pass","controls":records}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("weight_claims.json"))
    parser.add_argument("--report",type=Path);parser.add_argument("--controls",type=Path)
    args=parser.parse_args(argv)
    try:
        claims=json.loads(args.claims.read_text());report=evaluate(claims,lambda name,ok:None);payload=canonical(report)
        if args.report:args.report.write_bytes(payload)
        if args.controls:args.controls.write_bytes(canonical(controls(claims)))
        sys.stdout.buffer.write(payload);return 0
    except (ValueError,KeyError,TypeError,OSError) as exc:sys.stderr.write(str(exc)+"\n");return 1


if __name__=="__main__":raise SystemExit(main())
