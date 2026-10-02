#!/usr/bin/env python3
"""Independent saved-evidence audit and proper-time/zeta comparison.

This post-run audit recomputes gates; it does not trust producer PASS flags.
It strengthens, without altering, the frozen producer validator.
"""
import argparse
import copy
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:
        raise ValueError(message)


def number(value,nonnegative=False):
    result=Decimal(str(value))
    require(result.is_finite(),"Nonfinite value")
    require(not nonnegative or result>=0,"Negative error/threshold")
    return result


def audit_producer(document):
    expected_points={(x,y) for x in ["1","2"] for y in ["0.5","1"]}
    point_keys=[(p["x_over_r"],p["H2_over_r"]) for p in document["points"]]
    require(len(point_keys)==4 and set(point_keys)==expected_points,"Unexpected point set")
    suffixes=[f"{resolution}/{obs}_tail" for resolution in ["primary","refined"] for obs in ["W","rho","Q"]]
    suffixes += [f"{obs}_precision_refinement" for obs in ["W","rho","Q"]]
    suffixes += ["mass_variation","metric_variation","common_action_pairing","independent_digamma_Q","independent_trace_check"]
    expected_gates={f"x{x}_y{y}/{suffix}" for x,y in expected_points for suffix in suffixes}|{"local_V_F_action_pairing"}
    names=[g["name"] for g in document["gates"]]
    require(len(names)==len(expected_gates) and set(names)==expected_gates,"Unexpected gate names")
    for gate in document["gates"]:
        residual=number(gate["absolute_residual"],True)
        threshold=number(gate["threshold"],True)
        require(type(gate["pass"]) is bool and gate["pass"],"Nonboolean/failed gate")
        require(residual<=threshold,"Recorded gate exceeds threshold")
    negatives=document["negative_controls"]
    require(len(negatives)==3 and {g["name"] for g in negatives}==
            {"omit_metric_variation","flip_scalar_current","omit_F_scalar_current"},"Unexpected negatives")
    for neg in negatives:
        require(type(neg["detected"]) is bool and neg["detected"],"Nonboolean/undetected control")
        require(number(neg["residual"],True)>number(neg["threshold"],True),"Negative control fails numerically")
    recomputed=[]
    for point in document["points"]:
        high,low=point["refined"],point["primary"]
        for source in [high,low]:
            for obs in ["W","rho","p","Q","Wx","Wy"]:
                number(source[obs])
            require(number(source["p"])==-number(source["rho"]),"State-pressure mismatch")
            for bound in source["tail_bound"].values():
                require(number(bound,True)<=Decimal("1e-14"),"Series tail gate failed")
        for obs in ["W","rho","Q"]:
            delta=abs(number(high[obs])-number(low[obs]))
            tolerance=Decimal("1e-11")+Decimal("1e-9")*abs(number(high[obs]))
            require(delta<=tolerance,"Precision refinement failed")
            recomputed.append({"kind":"precision","observable":obs,"pass":True})
        for obs in ["rho","Q"]:
            delta=abs(number(high[obs])-number(point["independent_action_differences"][obs]))
            require(delta<=Decimal("1e-10")+Decimal("1e-8")*abs(number(high[obs])),"Action variation failed")
        require(abs(number(point["digamma_Q"])-number(high["Q"]))<=
                Decimal("1e-11")+Decimal("1e-9")*abs(number(high["Q"])),"Digamma check failed")
        # rhox is not saved separately, so enforce the stricter registered
        # absolute floor rather than trusting its recorded relative scale.
        require(abs(number(point["pairing_residual"]))<=Decimal("1e-10"),"Pairing absolute-floor gate failed")
    return recomputed


def audit_proper(document):
    settings={(40,Decimal(".01"),10),(55,Decimal(".005"),12)}
    expected={(Decimal(x),Decimal(z),setting) for x in [1,2] for z in [".5","1"] for setting in settings}
    keys=[]
    trace_gates=[]
    for point in document["results"]:
        st=point["settings"]
        setting=(st["dps"],number(st["join"]),st["order"])
        keys.append((number(point["x"]),number(point["z"]),setting))
        vals=[number(point[k]) for k in ["W","rho","Q"]]
        tolerance=Decimal("1e-9")+Decimal("1e-7")*max(map(abs,vals))
        # Initial aggregator reported trace residuals without including them in
        # its aggregate PASS. This independent post-run gate closes that gap.
        residual=abs(number(point["trace_residual"]))
        require(residual<=tolerance,"Proper-time trace gate failed")
        trace_gates.append({"x":point["x"],"z":point["z"],"dps":st["dps"],"residual":str(residual),"pass":True})
        for key in ["negative_current_sign_residual","negative_metric_sign_residual"]:
            require(abs(number(point[key]))>tolerance,"Proper-time negative control failed")
        for kind in ["spectral_tail_bounds","IR_tail_bounds"]:
            for bound in point[kind]:
                require(number(bound,True)<Decimal("1e-14"),"Endpoint/spectral bound gate failed")
    require(len(keys)==8 and set(keys)==expected,"Proper-time grid/settings mismatch")
    return trace_gates


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--producer",type=Path,required=True)
    parser.add_argument("--proper",type=Path,default=ROOT/"INDEPENDENT_PROPER_TIME.json")
    parser.add_argument("--output",type=Path,default=ROOT/"INDEPENDENT_AUDIT.json")
    args=parser.parse_args()
    producer=json.loads(args.producer.read_text())
    proper=json.loads(args.proper.read_text())
    with localcontext() as ctx:
        ctx.prec=90
        audit_producer(producer)
        traces=audit_proper(proper)
        comparison=[]
        for p in producer["points"]:
            independent=next(q for q in proper["results"] if number(q["x"])==number(p["x_over_r"])
                             and number(q["z"])==number(p["H2_over_r"]) and q["settings"]["dps"]==55)
            for obs in ["W","rho","Q"]:
                delta=abs(number(independent[obs])-number(p["refined"][obs]))
                tolerance=Decimal("1e-9")+Decimal("1e-7")*abs(number(independent[obs]))
                require(delta<=tolerance,"Cross-method disagreement")
                comparison.append({"x":p["x_over_r"],"H2":p["H2_over_r"],"observable":obs,
                                   "difference":str(delta),"threshold":str(tolerance),"pass":True})
        mutations=[]
        for label,mutate in [
            ("renamed_gate",lambda d:d["gates"][0].update(name="made_up")),
            ("infinite_threshold",lambda d:d["gates"][0].update(threshold="Infinity")),
            ("nonboolean_pass",lambda d:d["gates"][0].update({"pass":"yes"})),
            ("negative_tail",lambda d:d["points"][0]["primary"]["tail_bound"].update(W="-1")),
            ("duplicate_point",lambda d:d["points"].__setitem__(1,copy.deepcopy(d["points"][0]))),
        ]:
            corrupt=copy.deepcopy(producer)
            mutate(corrupt)
            try:
                audit_producer(corrupt)
            except (ValueError,KeyError):
                mutations.append({"name":label,"rejected":True})
            else:
                raise RuntimeError("Audit failed to reject "+label)
        report={"status":"PASS","producer_gate_name_set_verified":57,"cross_method_gates":comparison,
                "post_run_proper_time_trace_gates":traces,"validator_negative_controls":mutations,
                "input_hashes":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.producer,args.proper]},
                "independent_total_error_enclosure":False,
                "trace_is_not_independent_numerical_error_estimate":True}
        args.output.write_text(json.dumps(report,indent=2)+"\n")
        print(json.dumps({"status":"PASS","cross_method_gates":len(comparison),"trace_gates":len(traces),
                          "mutation_rejections":len(mutations),"maximum_cross_method_difference":str(max(number(g["difference"]) for g in comparison))},indent=2))


if __name__=="__main__":
    main()
