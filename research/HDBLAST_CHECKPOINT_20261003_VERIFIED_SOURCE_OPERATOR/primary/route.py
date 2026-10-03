"""Guarded primary producer for the separately registered source/operator study.

No bootstrap authorization is implemented here. Root's complete driver must
authenticate its public freeze and runtime/input contract before supplying the
authorization callable. The CLI supports fabricated fixtures only.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import os
import resource
import time
from typing import Callable

from flint import acb,arb,ctx
if __package__:
    from . import verified_moments as v
else:
    import verified_moments as v

SOURCES=("positive_B","signed_uB")
MOMENTA=(Fraction(0),Fraction(1,2**40),Fraction(1,2**12),Fraction(1,4),
         Fraction(1),Fraction(16),Fraction(64),Fraction(128),Fraction(256))
CONFIGURATION="PRIMARY_ARB256_SOURCE24_PHASE96"
GATE=Fraction(1,10**26)


def require(condition,message):
    if not condition:raise ValueError(message)


def qstring(value)->str:
    """Canonical reduced n/d, including 0/1 and every integer denominator1."""
    fraction=Fraction(value)
    return f"{fraction.numerator}/{fraction.denominator}"


def bounds(value:arb)->dict[str,str]:
    endpoints=v.rational_endpoints(value)
    return {"lo":qstring(endpoints["lower"]),"hi":qstring(endpoints["upper"])}


def rectangle(value:arb|acb)->dict[str,dict[str,str]]:
    if isinstance(value,arb):value=acb(value,arb(0))
    return {"real":bounds(value.real),"imag":bounds(value.imag)}


def rectangle_radius(encoded)->Fraction:
    return sum(((Fraction(encoded[name]["hi"])-Fraction(encoded[name]["lo"]))/2
                for name in ("real","imag")),Fraction(0))


def exact_positive_upper(value:arb)->Fraction:
    return Fraction(str(value.abs_upper().fmpq()))


def polynomial_l1(panel:v.PolynomialPanel)->Fraction:
    return sum((exact_positive_upper(a) for a in panel.coefficients),Fraction(0))


def zero_model_panel(panel:v.PolynomialPanel)->v.PolynomialPanel:
    return v.PolynomialPanel(panel.center,panel.half_width,panel.coefficients,Fraction(0))


class ArithmeticKernelFamily:
    """Finite phase-polynomial baseline for positive arithmetic error reporting.

    These kernels alone are not a certificate for the full exponential. Their
    exact analytic tails are recorded separately in the budget output. They
    still enclose coefficient rounding and all finite-polynomial arithmetic.
    """
    def __init__(self,family:v.EntireKernelFamily):self.family=family
    def at(self,k,half_width=v.HALF_WIDTH):
        require(0<=Fraction(k)<=256 and half_width==v.HALF_WIDTH,"Invalid arithmetic kernel geometry")
        z=acb(0,v.rational_ball(2*Fraction(k)*half_width))
        return (tuple(p(z) for p in self.family.exp_polynomials),
                tuple(p(z) for p in self.family.drift_polynomials))


def real_source_moments(moment:dict,k:Fraction)->dict:
    """Registered exact real-target projection; does not project incoming modes.

    For k=0, exp=1 and Phi(d)=d are exactly real on the real integration path.
    All source coefficients are real Arb balls, and the represented analytic
    source is real there. Thus Mexp/Mu imaginary components equal exactly zero.
    This is a proved target identity, not empirical interval narrowing.
    """
    output={"M0":moment["M0"],"Mexp":moment["Mexp"],"Mu":moment["Mdrift"]}
    if k==0:
        output["Mexp"]=acb(output["Mexp"].real,arb(0))
        output["Mu"]=acb(output["Mu"].real,arb(0))
    return output


def moment_row(source,k,interval,moment,panel=None):
    values=real_source_moments(moment,k)
    encoded={name:rectangle(value) for name,value in values.items()}
    row={"source":source,"momentum":qstring(k),"interval":[qstring(x) for x in interval],
         "moments":encoded,"total_absolute_radii":{name:qstring(rectangle_radius(value))
                         for name,value in encoded.items()}}
    if panel is not None:row["panel"]=panel
    return row


def work_row(source,interval,moment,panel=None):
    encoded=rectangle(moment)
    row={"source":source,"interval":[qstring(x) for x in interval],"moment":encoded,
         "total_absolute_radius":qstring(rectangle_radius(encoded))}
    if panel is not None:row["panel"]=panel
    return row


def moment_budget(source,k,interval,complete_row,baseline,source_disks,phase_disks,panel=None):
    baseline_encoded={name:rectangle(value) for name,value in real_source_moments(baseline,k).items()}
    row={"source":source,"momentum":qstring(k),"interval":[qstring(x) for x in interval],
         "source_model_disk_radius":{name:qstring(value) for name,value in source_disks.items()},
         "phase_model_disk_radius":{name:qstring(value) for name,value in phase_disks.items()},
         "coefficient_and_arithmetic_baseline_L1_radius":{name:qstring(rectangle_radius(value))
                         for name,value in baseline_encoded.items()},
         "complete_output_L1_radius":complete_row["total_absolute_radii"]}
    if panel is not None:row["panel"]=panel
    return row


def run_primary(authorize:Callable[[str,dict],None]|None,
                source_bundle_provider:Callable|None=None)->tuple[dict,dict]:
    """Return exact schema1 data and a separate positive error-budget receipt.

    source_bundle_provider is a fabricated-only dependency injection for
    pre-freeze tests. The registered physical driver must omit that argument
    and pin this source; root owns enforcement of that entry contract.
    """
    require(callable(authorize),"Root bootstrap authorization required")
    provider=source_bundle_provider or v.future_registered_source_bundle
    data={"schema_version":1,"scope":"UNIFORM_ANALYTIC_MODEL_PLUS_FIXED_RATIONAL_MOMENT_PROBES",
          "configuration":CONFIGURATION,"archive_arrays_decoded":0,
          "panel_rows":[],"whole_rows":[],"source_work_panel_rows":[],"source_work_whole_rows":[]}
    budget={"schema_version":1,"configuration":CONFIGURATION,"panel_rows":[],"whole_rows":[],
            "source_work_panel_rows":[],"source_work_whole_rows":[],
            "budget_semantics":"Nonnegative source and phase disk contributions plus the enclosure width of the specified finite-polynomial computation with source and local kernel model tails suppressed. The baseline includes any tiny stable-drift series enclosure. This is not an additive partition of the complete interval width, nor a bound on all full-run rounding after model tails are added. Complete output radii retain every enclosure and are the acceptance quantities; no subtraction or cancellation isolates an error.",
            "whole_gate":"1e-26 complete complex L1 radius from outward exact rational endpoints",
            "analytic_model":{"source_degree":24,"phase_degree":96,"source_complex_disk_radius":"1/8",
                 "panel_half_width":"1/128","source_majorant":"64/1","Lg_majorant":"32/1",
                 "source_uniform_tail":qstring(v.SOURCE_REMAINDER),"Lg_uniform_tail":qstring(v.SOURCE_WORK_REMAINDER),
                 "local_phase_cap":"4/1","integrand_phase_cap":"8/1","exponential_majorant":"4096/1",
                 "phase_pointwise_E_tail":qstring(Fraction(4096*8**97,factorial(97))),
                 "phase_pointwise_Q_tail":qstring(Fraction(2*4096*8**97,factorial(98))),
                 "phase_integrated_E_tail":qstring(Fraction(2*4096*8**97,factorial(97))),
                 "phase_integrated_Q_tail":qstring(Fraction(4*4096*8**97,factorial(98))),
                 "zero_momentum_identity":"All source moments are real at k=0; exact imaginary zero is projected only for the real-source moment target, never for incoming complex states."}}
    with ctx.workprec(v.PRECISION_BITS):
        family=v.EntireKernelFamily();baseline_family=ArithmeticKernelFamily(family)
        kernels={k:family.at(k) for k in MOMENTA}
        baseline_kernels={k:baseline_family.at(k) for k in MOMENTA}
        for source in SOURCES:
            panels=[];work_panels=[];work_total=arb(0);work_baseline_total=arb(0)
            for j in range(64):
                center=Fraction(-9,2)+Fraction(2*j+1,128)
                bundle=provider(center,source,authorize)
                panel,work_panel=bundle["g"],bundle["Lg"]
                require(panel.center==center and work_panel.center==center,"Unexpected source panel center")
                require(panel.half_width==v.HALF_WIDTH and work_panel.half_width==v.HALF_WIDTH,"Unexpected source half-width")
                require(len(panel.coefficients)==25 and len(work_panel.coefficients)==25,"Incomplete source polynomial")
                require(panel.source_remainder==v.SOURCE_REMAINDER and work_panel.source_remainder==v.SOURCE_WORK_REMAINDER,
                        "Unregistered source remainder")
                panels.append(panel);work_panels.append(work_panel)
                interval=(center-v.HALF_WIDTH,center+v.HALF_WIDTH)
                c_l1=polynomial_l1(panel)
                for k in MOMENTA:
                    m=v.panel_moments(panel,kernels[k])
                    baseline=v.panel_moments(zero_model_panel(panel),baseline_kernels[k])
                    row=moment_row(source,k,interval,m,j);data["panel_rows"].append(row)
                    h=v.HALF_WIDTH;r=panel.source_remainder
                    source_disks={"M0":2*h*r,"Mexp":2*h*r,"Mu":2*h*h*r}
                    phase_disks={"M0":Fraction(0),"Mexp":h*c_l1*family.e_tail if k else Fraction(0),
                                 "Mu":h*h*c_l1*family.q_tail if k else Fraction(0)}
                    budget["panel_rows"].append(moment_budget(source,k,interval,row,baseline,source_disks,phase_disks,j))
                work_m=v.panel_moments(work_panel,kernels[Fraction(0)])["M0"]
                work_base=v.panel_moments(zero_model_panel(work_panel),baseline_kernels[Fraction(0)])["M0"]
                work_total+=work_m;work_baseline_total+=work_base
                row=work_row(source,interval,work_m,j);data["source_work_panel_rows"].append(row)
                budget["source_work_panel_rows"].append({"source":source,"panel":j,"interval":[qstring(x) for x in interval],
                   "source_model_radius":qstring(2*v.HALF_WIDTH*work_panel.source_remainder),"phase_model_radius":"0/1",
                   "coefficient_and_arithmetic_baseline_radius":qstring(rectangle_radius(rectangle(work_base))),
                   "complete_output_radius":row["total_absolute_radius"]})
            interval=(Fraction(-9,2),Fraction(-7,2));zero_panels=tuple(zero_model_panel(p) for p in panels)
            for k in MOMENTA:
                m=v.forced_response(panels,k,acb(0),acb(0),Fraction(1),family)
                baseline=v.forced_response(zero_panels,k,acb(0),acb(0),Fraction(1),baseline_family)
                row=moment_row(source,k,interval,m);data["whole_rows"].append(row)
                r=v.SOURCE_REMAINDER;phase_exp=Fraction(0);phase_u=Fraction(0)
                for panel in panels:
                    c_l1=polynomial_l1(panel);h=panel.half_width;distance=interval[1]-panel.center-h
                    if k:
                        phase_exp+=h*c_l1*family.e_tail
                        phase_u+=h*h*c_l1*family.q_tail
                        if 2*k*distance<=1:
                            phase_u+=distance*Fraction(3,factorial(98))*(2*h*c_l1)
                budget["whole_rows"].append(moment_budget(source,k,interval,row,baseline,
                          {"M0":r,"Mexp":r,"Mu":r/2},
                          {"M0":Fraction(0),"Mexp":phase_exp,"Mu":phase_u}))
            row=work_row(source,interval,work_total);data["source_work_whole_rows"].append(row)
            budget["source_work_whole_rows"].append({"source":source,"interval":[qstring(x) for x in interval],
                   "source_model_radius":qstring(v.SOURCE_WORK_REMAINDER),"phase_model_radius":"0/1",
                   "coefficient_and_arithmetic_baseline_radius":qstring(rectangle_radius(rectangle(work_baseline_total))),
                   "complete_output_radius":row["total_absolute_radius"]})
    require([len(data[key]) for key in ("panel_rows","whole_rows","source_work_panel_rows","source_work_whole_rows")]
            ==[1152,18,128,2],"Incomplete primary registered universe")
    whole_radii=[Fraction(radius) for row in data["whole_rows"] for radius in row["total_absolute_radii"].values()]
    whole_radii.extend(Fraction(row["total_absolute_radius"]) for row in data["source_work_whole_rows"])
    budget["maximum_whole_complete_L1_radius"]=qstring(max(whole_radii))
    budget["status"]="PASS_CERTIFIED_OPERATOR_WIDTH_GATE" if all(radius<=GATE for radius in whole_radii) else "UNRESOLVED_CERTIFICATE"
    return data,budget


def fabricated_bundle(center,source,authorize):
    authorize("fabricated_source",{"center":str(center),"source":source})
    sign=1 if source=="positive_B" else -1
    coefficients=tuple(v.rational_ball(Fraction(sign**j,(j+1)*128**j)) for j in range(25))
    return {"g":v.PolynomialPanel(center,v.HALF_WIDTH,coefficients,v.SOURCE_REMAINDER),
            "Lg":v.PolynomialPanel(center,v.HALF_WIDTH,tuple(a/4 for a in coefficients),v.SOURCE_WORK_REMAINDER)}


def run(auth):
    """Production root-driver interface; callback allowed only after bootstrap."""
    return run_primary(auth)


def manufactured_arb_operation_rehearsal():
    """Exercise jet operations on unrelated manufactured functions only.

    Positive manufactured centers and exp(t/(2+t^2)) are unrelated to either
    registered bump source. Finite coefficient checks rehearse operations;
    they do not certify a registered source value or its analytic tail.
    """
    previous_cap=ctx.cap
    try:
        with ctx.workprec(v.PRECISION_BITS):
            ctx.cap=27
            for j in range(128):
                center=Fraction(2)+Fraction(2*j+1,128)
                t=v.arb_series([v.rational_ball(center),v.rational_ball(v.HALF_WIDTH)],prec=27)
                manufactured=(t*(2+t*t).inv()).exp()
                g,work=v.forcing_series_coefficients(t,manufactured,v.HALF_WIDTH,24)
                require(len(g)==25 and len(work)==25,"Incomplete manufactured operation rehearsal")
                require(all(a.is_finite() for a in g+work),"Nonfinite manufactured rehearsal coefficient")
    finally:
        ctx.cap=previous_cap
    return {"scope":"MANUFACTURED_ARB_OPERATIONS_ONLY",
            "registered_source_calls":0,"manufactured_operations":128,
            "g_coefficients_checked_per_operation":25,"Lg_coefficients_checked_per_operation":25,
            "source_tail_certification":False}


def run_fabricated():
    """Root's no-auth fabricated entry; refuses every real-source callback."""
    rehearsal=manufactured_arb_operation_rehearsal()
    def fabricated_authorize(name,detail):
        if name!="fabricated_source":
            raise RuntimeError("Real callback attempted in fabricated route")
    data,budget=run_primary(fabricated_authorize,fabricated_bundle)
    budget["manufactured_operation_rehearsal"]=rehearsal
    return data,budget


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--fabricated",action="store_true")
    parser.add_argument("--output-directory",type=Path,required=True)
    args=parser.parse_args()
    if not args.fabricated:raise RuntimeError("Real run must use root's authenticated bootstrap; this CLI is fabricated-only")
    started=time.monotonic();events=[]
    data,budget=run_primary(lambda name,detail:events.append(name),fabricated_bundle)
    require(budget["status"]=="PASS_CERTIFIED_OPERATOR_WIDTH_GATE","Fabricated complete width gate failed")
    require(events.count("fabricated_source")==128 and "real_source_taylor" not in events,"Unexpected physical source evaluation")
    for rows in (data["panel_rows"],data["whole_rows"]):
        for row in rows:
            for name,moment in row["moments"].items():
                if name=="M0" or row["momentum"]=="0/1":
                    require(moment["imag"]=={"lo":"0/1","hi":"0/1"},"Zero-momentum real identity not enforced")
    for rows in (data["source_work_panel_rows"],data["source_work_whole_rows"]):
        require(all(row["moment"]["imag"]=={"lo":"0/1","hi":"0/1"} for row in rows),"Source work must be real")
    args.output_directory.mkdir(parents=True,exist_ok=True)
    (args.output_directory/"PRIMARY.json").write_text(json.dumps(data,sort_keys=True)+"\n")
    (args.output_directory/"PRIMARY_ERROR_BUDGET.json").write_text(json.dumps(budget,sort_keys=True)+"\n")
    receipt={"status":"PASS_FABRICATED_PRIMARY_COMPLETE_ENTRY","physical_source_evaluations":0,
        "archive_arrays_decoded":0,"uid":os.getuid(),"moment_rows":1170,"source_work_rows":130,
        "wall_seconds":time.monotonic()-started,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "maximum_whole_complete_L1_radius":budget["maximum_whole_complete_L1_radius"],
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "engine_sha256":hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest(),
        "output_sha256":hashlib.sha256((args.output_directory/"PRIMARY.json").read_bytes()).hexdigest(),
        "error_budget_sha256":hashlib.sha256((args.output_directory/"PRIMARY_ERROR_BUDGET.json").read_bytes()).hexdigest(),
        "authorization_events":{name:events.count(name) for name in set(events)},
        "limits":"Fabricated provider only; root bootstrap/public freeze/independent route/external receipt remain separate."}
    (args.output_directory/"RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":main()
