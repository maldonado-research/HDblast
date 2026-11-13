#!/usr/bin/env python3
"""Finite arithmetic for the accepted conditional initial constraint ball.

Analytic embedding, map, series, domain and fixed-point theorems are separate
pinned inputs. Large amplitudes/constants stay formal decimal monomials;
neither exp() nor giant decimal expansion of 10^6400 is used.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys


# Public control metadata for fresh whole-capsule rehash controls. Each pair
# names a flat DATA.arithmetic field and its mathematically invalid value.
MUTATIONS=(
 ("connection_rank_factor",2),("point_rank_ceiling",31),
 ("embedding_array_total",51376),("embedding_exponent",5),
 ("product_exponent",16),("source_normalization_exponent",17),
 ("metric_correction_integer","11743112956334838628730019840"),
 ("background_jet_exponent",40),("metric_guard_denominator",2),("connection_second_factor",8),
 ("Ricci_second_factor",300),("curvature_Hessian_factor",500),
 ("gradient_Hessian_factor",50),("potential_Hessian_factor",100),
 ("metric_momentum_Hessian_factor",10000),("scalar_density_Hessian_factor",1),
 ("momentum_Hessian_factor",10),("combined_Hessian_factor",1000000),
 ("C_F_exponent",1000),("D2G_exponent",6000),
 ("K_exponent",6000),("amplitude_exponent",-6300),
 ("potential_radius_exponent",-6405),("open_interval_factor","2"),
 ("closed_theta","1/100"),("open_selfmap_fraction","1/2"),
 ("first_parameter_derivative_bound","1/50"),
 ("second_parameter_derivative_K_factor",2))


def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+"\n").encode()


def pow10(exponent):
    """Only for small normalized exponents; large ones stay symbolic."""
    if type(exponent) is not int or abs(exponent)>1000:raise ValueError("large decimal expansion forbidden")
    return F(10**exponent) if exponent>=0 else F(1,10**(-exponent))


def normalize_terms(terms):
    """Exact mantissa*10^exponent for short sums of comparable monomials."""
    terms=[(F(c),p) for c,p in terms if c]
    if not terms:return F(0),0
    power=max(p for _,p in terms)
    return sum((c*pow10(p-power) for c,p in terms),F(0)),power


def monomial_below(mantissa,exponent,bound):
    """Exact positive comparison using an explicitly coarser decimal power.

    If exponent<-100, replace it by -100 for a certified upper bound. All
    comparisons in this ledger have a strict margin at that coarser power.
    No huge rational denominator is constructed.
    """
    if mantissa<0 or bound<=0:raise ValueError("positive monomial comparison required")
    return mantissa*pow10(max(exponent,-100))<bound


def evaluate(claims,check):
    records=[]
    def need(name,ok):
        boolean=bool(ok);check("ball_arithmetic/"+name,boolean)
        records.append({"name":name,"status":"pass" if boolean else "fail"})
        if not boolean:raise ValueError("ball arithmetic rejects "+name)
    need("schema",claims["schema"]=="hdblast-nonlinear-ball-arithmetic-claims-v1")
    need("review_status",claims["review_status"]=="accepted_conditional_analytic_inputs")
    # Finite geometric constants: none of these comparisons proves embedding.
    need("S3_volume_kernel_reserve",F(12,2*3**2)<1)
    need("local_coefficient_comparability",F(1,1-F(2,10**8))<2)
    rankfactor=claims["connection_rank_factor"]
    need("aggregate_connection_rank_factor",rankfactor>0 and rankfactor**2>=12)
    rank=claims["point_rank_ceiling"]
    need("tensor_rank_budget",claims["product_input_rank_ceiling"]==25 and 25+7<=rank==32)
    a=lambda r:rankfactor*r
    arrays=[1,2,1+a(rank),1+a(rank+1)+1+a(rank),
            1+a(rank+1)+a(rank)+a(rank)**2,
            1+a(rank+2)+a(rank+1)+1+a(rank+1)**2+a(rank+1)+a(rank)+a(rank)**2]
    need("embedding_array_total",sum(arrays)==claims["embedding_array_total"]==51377)
    need("embedding_array_rounding",sum(arrays)<10**5)
    embedding_exp=claims["embedding_exponent"]
    need("tensor_embedding_rounding",8*10**5<10**embedding_exp)
    need("normalization_multiplier_recurrence",all(1+5*k<=5*(k+1) for k in range(12)))
    need("weighted_jet_multiplier_recurrence",max(j+5*k for j in range(17) for k in range(3))==26<32)
    jet_sum=sum(33**k for k in range(4))
    need("weighted_jet_multiplier_sum",jet_sum==37060<10**5)
    low_exp=claims["low_jet_exponent"]
    need("pointwise_jet_rounding",embedding_exp+5<=low_exp==12)
    need("graded_product_rounding",10**low_exp*sum(2**j for j in range(15))<10**claims["product_exponent"])
    source_sum=sum(math.comb(m,k)*5**k*math.factorial(k) for m in range(13) for k in range(m+1))
    need("literal_source_normalization_sum",source_sum==int(claims["source_normalization_sum"]))
    need("literal_source_normalization_rounding",source_sum<10**claims["source_normalization_exponent"])
    # Exact all-mode correction-map sums, including all derivative slots.
    CH=8*17**14*math.factorial(14);CM=8*17**13*math.factorial(13)
    need("metric_correction_integer",CH==int(claims["metric_correction_integer"]))
    need("momentum_correction_integer",CM==int(claims["momentum_correction_integer"]))
    need("correction_recurrence_stage",all(6*j+17<=17*(j+1) for j in range(14)))
    need("correction_literal_sum_geometric_reserve",F(2,1)/(1-F(1,17))==F(17,8)<8)
    need("metric_correction_E_rounding",CH<10**claims["metric_correction_E_exponent"])
    need("momentum_correction_E_rounding",CM<10**claims["momentum_correction_E_exponent"])
    need("physical_flat_multiplier",24**24<10**34 and 24*8+34==claims["physical_multiplier_exponent"]==226)
    need("physical_metric_map",claims["metric_correction_E_exponent"]+226==claims["physical_metric_map_exponent"]==256)
    need("physical_momentum_map",claims["momentum_correction_E_exponent"]+226==claims["physical_momentum_map_exponent"]==253)
    need("correction_Hilbert_sum",CH+CM<10**31)
    background=20*claims["background_stage_factor"]**14
    need("background_stage_budget",claims["background_stage_factor"]==1000)
    need("background_jet_rounding",background<10**claims["background_jet_exponent"])
    need("background_multiplier_rounding",(2**15-1)*10**claims["background_jet_exponent"]<10**claims["background_multiplier_exponent"])
    Qexp=claims["Q_exponent"];Bexp=claims["B_exponent"]
    need("Q_contains_accepted_interfaces",Qexp==100 and max(embedding_exp,low_exp,claims["product_exponent"],
         claims["source_normalization_exponent"],claims["background_multiplier_exponent"])+5<Qexp)
    need("B_contains_seed_and_correction",Bexp==200 and 31<Bexp)
    # Uniform unital inverse/determinant and physical Hessian ledger (15)-(27).
    guard_denominator=claims["metric_guard_denominator"]
    need("uniform_positive_metric_guard",type(guard_denominator) is int and guard_denominator>4)
    t=F(1,guard_denominator)
    need("Neumann_inverse_state",1/(1-t)<2)
    need("inverse_module",1+2*t<2)
    need("log_determinant_state",2/(1-t)<4)
    need("inverse_determinant_state",4/(1-4*t)<8 and 1+8*t<2)
    need("inverse_first_second_factors",2**2==4 and 2*2**3==16)
    need("determinant_first_second_factors",2*4==8 and 2*(4**2+8)==48)
    need("connection_state_first_factors",F(1,2)*2*3==3 and 3+6*t<4)
    need("connection_second_factor",12+24*t<claims["connection_second_factor"])
    need("Ricci_state_factor",12+72*t<claims["Ricci_state_factor"])
    need("Ricci_first_factor",16+192*t<claims["Ricci_first_factor"])
    need("Ricci_second_factor",64+768*t+256<claims["Ricci_second_factor"])
    curvature=F(16+208*t+136+3*claims["Ricci_second_factor"],2)
    need("curvature_Hessian_factor",curvature<claims["curvature_Hessian_factor"])
    need("gradient_Hessian_factor",F(64+16+16+6,2)<claims["gradient_Hessian_factor"])
    U=[F(z) for z in ("-2/27","0","14/9","50/27","-1/6","-4/9","-2/27")]
    need("actual_polynomial_coefficient_sum",sum(map(abs,U),F(0))<10)
    U2=sum((k*(k-1)*abs(v) for k,v in enumerate(U)),F(0))
    need("actual_potential_second_coefficients",U2<300)
    need("potential_Hessian_factor",300*2**4<claims["potential_Hessian_factor"])
    need("metric_momentum_state_first_second_counts",4*100<1000 and 4*4*100<10000 and 16*4*100<10000)
    momentum=2*(2*10**4+16*10**4*t+48*10**3*t*t)
    need("metric_momentum_Hessian_factor",momentum<claims["metric_momentum_Hessian_factor"])
    need("scalar_density_Hessian_factor",F(48*t*t+32*t+4,2)<claims["scalar_density_Hessian_factor"])
    need("momentum_constraint_Hessian_factor",2*2*(2+2+1+1)<claims["momentum_Hessian_factor"])
    component_sum=sum(claims[name] for name in ("curvature_Hessian_factor","gradient_Hessian_factor", "potential_Hessian_factor",
        "metric_momentum_Hessian_factor","scalar_density_Hessian_factor","momentum_Hessian_factor"))
    need("combined_physical_Hessian_factor",component_sum<claims["combined_Hessian_factor"]==10**7)
    derived_CF=7+6*Qexp+2*Bexp
    need("normalized_full_D2F",derived_CF==1007<=claims["C_F_exponent"]==1100)
    need("accepted_linear_inverse",claims["inverse_exponent"]==5000)
    need("preconditioned_full_D2G",claims["inverse_exponent"]+claims["C_F_exponent"]==claims["D2G_exponent"]==6100)
    Kexp=claims["K_exponent"];aexp=claims["amplitude_exponent"];rexp=claims["potential_radius_exponent"]
    need("K_dominates_D2G",Kexp>=claims["D2G_exponent"] and Kexp==6400)
    need("radius_is_K_a_squared",rexp==Kexp+2*aexp)
    # All large powers cancel before rational evaluation of contraction costs.
    theta=pow10(Kexp+aexp)+pow10(Kexp+rexp)
    need("closed_contraction_constant",theta==F(claims["closed_theta"])==F(101,10000)<F(1,2))
    center_ratio=pow10(Kexp+2*aexp-rexp)/2
    need("closed_center_radius",center_ratio==F(1,2))
    closed_self=center_ratio+theta
    need("closed_selfmap",closed_self==F(claims["closed_selfmap_fraction"])==F(5101,10000)<1)
    domain=normalize_terms(((1,Qexp+Bexp+aexp),(1,Qexp+Bexp+rexp)))
    need("closed_ball_open_domain",monomial_below(*domain,F(1,2*guard_denominator)))
    extension=F(claims["open_interval_factor"])
    need("strict_open_extension",extension>1)
    theta_open=extension*pow10(Kexp+aexp)+pow10(Kexp+rexp)
    need("open_contraction_constant",theta_open==F(claims["open_theta"])==F(51,5000)<F(1,2))
    self_open=extension*extension*center_ratio+theta_open
    need("open_interval_selfmap",self_open==F(claims["open_selfmap_fraction"])==F(2081,4000)<1)
    domain_open=normalize_terms(((extension,Qexp+Bexp+aexp),(1,Qexp+Bexp+rexp)))
    need("open_interval_positive_domain",monomial_below(*domain_open,F(1,2*guard_denominator)))
    need("fixedpoint_quadratic_bound",1/(2*(1-theta))<1)
    first=2*theta
    need("first_parameter_derivative",first==F(claims["first_parameter_derivative_bound"])==F(101,5000)<1)
    need("second_parameter_derivative",2*(1+first)**2<claims["second_parameter_derivative_K_factor"]==8)
    need("positive_metric_margin",t<F(1,4))
    return {"schema":"hdblast-nonlinear-ball-arithmetic-v1","status":"pass",
       "review_status":claims["review_status"],"claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
       "exact":{"embedding_arrays":arrays,"source_normalization_sum":str(source_sum),
          "metric_correction_integer":str(CH),"momentum_correction_integer":str(CM),
          "background_jet_budget":str(background),"potential_second_coefficient_sum":str(U2),
          "uniform_positive_metric_guard":str(t),
          "combined_physical_Hessian_factor":component_sum,"derived_D2F_ten_exponent":derived_CF,
          "C_F_ten_exponent":claims["C_F_exponent"],"D2G_ten_exponent":claims["D2G_exponent"],
          "K_ten_exponent":Kexp,"amplitude_ten_exponent":aexp,"potential_radius_ten_exponent":rexp,
          "closed_theta":str(theta),"closed_selfmap":str(closed_self),"open_theta":str(theta_open),"open_selfmap":str(self_open),
          "closed_QB_domain_monomial":{"mantissa":str(domain[0]),"ten_exponent":domain[1]},
          "open_QB_domain_monomial":{"mantissa":str(domain_open[0]),"ten_exponent":domain_open[1]},
          "first_parameter_derivative":str(first),"second_parameter_derivative_K_factor":8},
       "scope":"Finite exact arithmetic conditional on separately accepted embedding/product, correction-map, actual-seed, full nonlinear Hessian and selected-inverse proofs. The resulting interval concerns initial constraints only; no evolved PDE, observer tube, or physical prediction follows."}


def controls(claims):
    records=[]
    for name,value in MUTATIONS:
        bad=copy.deepcopy(claims);bad[name]=value;encoded=canonical(bad)
        try:evaluate(json.loads(encoded),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name,"self_rehashed_sha256":hashlib.sha256(encoded).hexdigest(),
            "rejected_by":"fresh exact arithmetic oracle","diagnostic":str(exc)})
        else:raise ValueError("surviving arithmetic fault: "+name)
    return {"schema":"hdblast-nonlinear-ball-arithmetic-controls-v1","status":"pass","controls":records}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("ball_claims.json"))
    parser.add_argument("--report",type=Path);parser.add_argument("--controls",type=Path)
    args=parser.parse_args(argv)
    try:
        claims=json.loads(args.claims.read_text());report=evaluate(claims,lambda name,ok:None);payload=canonical(report)
        if args.report:args.report.write_bytes(payload)
        if args.controls:args.controls.write_bytes(canonical(controls(claims)))
        sys.stdout.buffer.write(payload);return 0
    except (ValueError,KeyError,TypeError,OSError) as exc:sys.stderr.write(str(exc)+"\n");return 1


if __name__=="__main__":raise SystemExit(main())
