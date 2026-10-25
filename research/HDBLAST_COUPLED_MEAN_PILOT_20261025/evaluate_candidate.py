#!/usr/bin/env python3
"""Evaluate the exact finite feedback and enclose six complete candidate samples."""
import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ridge_kernel", HERE/"reflected_ridge.py")
kernel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kernel)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native-output", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    need(not a.output.exists(), "fresh candidate output")
    report_raw = (a.native_output/"REPORT.json").read_bytes()
    profile_raw = (a.native_output/"PROFILES.json").read_bytes()
    report = json.loads(report_raw)
    data = json.loads(profile_raw)
    source = kernel.from_data(data)
    error = Q(report["levels"][-1]["ledger"]["total_X_error_upper"])
    need(data["mean_order"] == list(kernel.MEANS) and data["gram"] == list(kernel.GRAM),
         "five-mean order and Gram")
    M = data["matrices"]["M"]
    rows = []
    for sample in data["free_wave_samples"]:
        x, r = Q(sample["x"]), Q(sample["r"])
        feedback = source.feedback(x, r, M)
        candidate = {}
        solution_enclosure = {}
        for name in ("value", "dx", "dr"):
            pairs = [[str(Q(pair[0])+increment), str(Q(pair[1])+increment)]
                     for pair, increment in zip(sample["B0"][name], feedback[name])]
            candidate[name] = pairs
            # X controls values/3 and both derivatives in the H norm.
            # Every Gram weight is >=1, so these component boxes are safe.
            cap = 3*error if name == "value" else error
            solution_enclosure[name] = [[str(Q(lo)-cap), str(Q(hi)+cap)]
                                        for lo,hi in pairs]
        candidate["ordinary_t_velocity"] = [[str(Q(lo)/Q(data["ell"])),
                                             str(Q(hi)/Q(data["ell"]))]
                                            for lo,hi in candidate["dx"]]
        solution_enclosure["ordinary_t_velocity"] = [
            [str(Q(lo)/Q(data["ell"])),str(Q(hi)/Q(data["ell"]))]
            for lo,hi in solution_enclosure["dx"]]
        if r == 0:
            need(feedback["value"][1] == feedback["dx"][1] == 0,
                 "Dirichlet feedback trace")
        rows.append({"x":str(x),"r":str(r),"feedback_exact":kernel.encoded(feedback),
                     "complete_candidate_interval_evaluation":candidate,
                     "conditional_exact_response_component_enclosure":solution_enclosure})
    result = {"status":"EVALUATED_CONTINUOUS_COUPLED_MEAN_CANDIDATE",
              "accuracy_target_passed":report["passes_target"],
              "native_report_sha256":hashlib.sha256(report_raw).hexdigest(),
              "profile_data_sha256":hashlib.sha256(profile_raw).hexdigest(),
              "total_X_error_upper":str(error),
              "representation":"Exact registered rational displacement, one continuous rational P1 velocity, exact reflected ridge feedback and exact shell lift.",
              "candidate_is_defined_on_whole_cone":True,"sample_count":len(rows),
              "samples":rows,
              "boundary_scope":"Dirichlet is exact. The finite candidate's remaining Robin defect is covered by the complete error identity; it is not set to zero.",
              "limits":["Conditional five-mean C1 approximation only.",
                        "Sample interval widths enclose arithmetic evaluation of the candidate; the whole-cone error is separately retained.",
                        "No full angular solution, pointwise Einstein residual, observable sign, nonlinear family or physical origin inference."]}
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],"accuracy_target_passed":report["passes_target"],
                      "sample_count":len(rows),
                      "sha256":hashlib.sha256(a.output.read_bytes()).hexdigest()}))


if __name__ == "__main__":
    main()
