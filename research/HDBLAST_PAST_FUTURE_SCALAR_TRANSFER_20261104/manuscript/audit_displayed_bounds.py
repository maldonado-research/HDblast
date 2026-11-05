"""Audit printed constants and enclosures against frozen exact reports.

Standard-library arithmetic only. No mathematical producer is imported or run.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

PRIMARY_REPORT_SHA = "23304f031d1189ff3f0e7482af1a4c16aaf0623cad7c6dd487ee1a4d5c82991a"
PRIMARY_PROOF_SHA = "ddfe9ad88773133f8f79758802724058453fb3b0c75e13aafb892ea960c05c10"
REFINEMENT_REPORT_SHA = "d6e9b1bccdd3a4de0d80aedc44541a29c2260e454c259d8d428ffad3118a7b12"
REFINEMENT_PROOF_SHA = "4dd87d647007295ded87c0af1be90040462cfb0890e98fd724d98a508550a468"
HERE = Path(__file__).resolve().parent

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def digest(raw):
    return sha256(raw).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["report", "optimized", "proof", "refinement", "refinement-optimized", "refinement-proof", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    need(out.is_absolute() and out.resolve() == out and not out.exists(),
         "Fresh canonical absolute output file required")
    need(all(not p.is_symlink() for p in [out, *out.parents]), "Output symlink rejected")
    need(not out.is_relative_to(HERE), "Audit output must be outside manuscript directory")
    paths = {"primary_report": args.report, "optimized_primary_report": args.optimized,
             "primary_proof": args.proof, "refinement_report": args.refinement,
             "optimized_refinement_report": args.refinement_optimized,
             "refinement_proof": args.refinement_proof,
             "editable_source": HERE / "PAST_FUTURE_SCALAR_COMPARISON.tex"}
    for path in paths.values():
        need(path.is_file() and all(not p.is_symlink() for p in [path, *path.absolute().parents]),
             "Regular input file required")
    data = {name: path.read_bytes() for name, path in paths.items()}
    report = json.loads(data["primary_report"])
    source = data["editable_source"].decode()
    checks = []

    def check(ok, name):
        need(ok, name)
        checks.append(name)

    check(digest(data["primary_report"]) == PRIMARY_REPORT_SHA, "frozen primary report pin")
    check(data["primary_report"] == data["optimized_primary_report"], "primary normal and optimized reports identical")
    check(digest(data["primary_proof"]) == PRIMARY_PROOF_SHA == report["proof_sha256"], "frozen continuous primary proof pin")
    check(report["schema"] == "HDBLAST_PAST_FUTURE_PRIMARY_V1" and
          report["status"] == "PASS_CONDITIONAL_PAST_FUTURE_PHYSICAL_SCALAR_TRANSFER" and
          report["finite_checks"] == len(report["check_names"]) == 42,
          "primary status and stated finite arithmetic scope")
    ell = F(report["state"]["ell"])
    check(ell == F(1, 10**6) and report["state"]["scaled_elapsed_interval"] == ["0", "3"],
          "registered width and comparison time interval")
    energy = report["difference_energy"]
    check(energy["initial_displacement"] == "0" and energy["initial_velocity"] == "-2ell*exp(-ell*r)*Q(r)",
          "difference has zero displacement and twice reversed velocity")
    check(F(energy["initial_basic_energy_upper"]) == 2*ell**2 and
          F(energy["boundary_master_difference_square_upper"]) == 2*ell**2/125 < ell**2,
          "printed complete basic energy and local trace bounds")
    epsilon = 63*ell**2/4400
    graph_initial = 4*ell**2*((16+ell)**2+103000*ell**2)
    check(F(energy["relative_temporal_loss_upper"]) == epsilon < 1,
          "printed relative temporal coercivity loss")
    check(F(energy["twice_initial_half_graph_energy_upper"]) == graph_initial,
          "zero-displacement initial graph energy")
    check(F(energy["A0_delta_state_squared_upper"]) == graph_initial/(1-epsilon) < (33*ell)**2 and
          F(energy["A0_delta_state_norm_strict_upper"]) == 33*ell,
          "printed graph norm ceiling with endpoint terms")
    scalar = report["physical_scalar_difference"]
    check(F(scalar["weighted_endpoint_sqrt_factor"]) == F(1, 20000) and
          F(1, 400)*ell == F(1, 20000)**2,
          "endpoint Hilbert mass conversion factor")
    check(F(scalar["prefactor_upper"]) == F(25000, 129**3) and
          F(scalar["q_term_upper"]) == ell/160,
          "physical scalar denominator and local master term")
    bound = F(25000, 129**3)*(F(33, 20000)/ell+ell/160)
    check(F(scalar["derived_uniform_upper"]) == bound == F(264000000001, 13738809600) < 20,
          "printed exact uniform scalar difference ceiling")
    check(F("19.21563") < bound < F("19.21564"), "readable ceiling decimal bracket")
    check(energy["endpoint_terms_retained"] is True and
          energy["nonlocal_half_power_compact_support_claimed"] is False and
          energy["physical_Noether_or_injected_energy_claimed"] is False,
          "endpoint and nonlocal graph-energy scope retained")
    check(scalar["extra_time_reflection_minus_sign"] is False and
          scalar["individual_amplitude_below_twenty_claimed"] is False and
          scalar["finite_nonlinear_parameter_window_supplied"] is False,
          "physical sign and linear-only scope retained")

    displayed = {}
    for key, token, expected_time in [("rising", r"\(5\ell/4\)&", F(5, 4)),
                                      ("falling", r"\(7\ell/4\)&", F(7, 4))]:
        row = report["selected_point_transfers"][key]
        exported = [F(v) for v in row["future_exported_ell_S_interval"]]
        future = [F(v) for v in row["future_actual_S_interval"]]
        past = [F(v) for v in row["conditional_past_actual_S_interval"]]
        check(F(row["matched_positive_scaled_time"]) == expected_time and
              future == [v/ell for v in exported], key + " export is ell times actual S")
        check(past == [future[0]-bound, future[1]+bound], key + " exact interval widening")
        lines = [line for line in source.splitlines() if line.startswith(token)]
        need(len(lines) == 1, "Unique printed table row required: " + key)
        matches = re.findall(r"\[(-?\d+\.\d+),(-?\d+\.\d+)\]", lines[0])
        need(len(matches) == 2, "Two printed decimal intervals required: " + key)
        printed_future = [F(v) for v in matches[0]]
        printed_past = [F(v) for v in matches[1]]
        check(printed_future[0] <= future[0] <= future[1] <= printed_future[1],
              key + " displayed future decimal enclosure is outward")
        check(printed_past[0] <= past[0] <= past[1] <= printed_past[1],
              key + " displayed past decimal enclosure is outward")
        displayed[key] = {"future": list(matches[0]), "past": list(matches[1])}
    zero = report["coarse_zero"]
    check(F(report["selected_point_transfers"]["rising"]["conditional_past_actual_S_interval"][1]) < -37376 and
          F(report["selected_point_transfers"]["falling"]["conditional_past_actual_S_interval"][0]) > 37379,
          "strict coarse past coefficient signs")
    check(zero["open_scaled_past_interval"] == ["-7/4", "-5/4"] and
          zero["at_least_one_by_retained_continuity"] is True and zero["uniqueness_claimed"] is False,
          "coarse open past zero interval and nonuniqueness scope")
    for key in ["central_left", "central_right", "late_left", "late_right"]:
        vals = [F(v) for v in report["selected_point_transfers"][key]["conditional_past_actual_S_interval"]]
        check(vals[0] < 0 < vals[1], key + " older narrow sign witness does not transfer")
    check(report["tail_nontransfer"]["past_sign_not_established"] is True and
          F(report["tail_nontransfer"]["future_physical_absolute_upper"]) < bound,
          "small future tail sign is not inferred in the past")

    # The selected-cell extension is bound only after its separate final review.
    check(digest(data["refinement_report"]) == REFINEMENT_REPORT_SHA,
          "frozen selected-cell refinement report pin")
    check(data["refinement_report"] == data["optimized_refinement_report"],
          "refinement normal and optimized reports identical")
    check(digest(data["refinement_proof"]) == REFINEMENT_PROOF_SHA,
          "frozen selected-cell continuous proof pin")
    refinement = json.loads(data["refinement_report"])
    check(refinement["schema"] == "HDBLAST_PAST_ZERO_CONTINUOUS_CELL_REFINEMENT_V1" and
          refinement["status"] == "PASS_CONDITIONAL_PAST_ZERO_CELL_REFINEMENT" and
          refinement["finite_guards"] == len(refinement["check_names"]) == 23,
          "selected-cell refinement status and arithmetic scope")
    check(F(refinement["ell"]) == ell and F(refinement["uniform_difference_upper"]) == bound and
          refinement["csv_rows"] == 4353,
          "whole-cell certificate uses same width difference ceiling and complete retained roster")
    for side, endpoints, safe in [("negative", [F(3071, 2048), F(6143, 4096)], F(-53)),
                                   ("positive", [F(6145, 4096), F(3073, 2048)], F(63))]:
        row = refinement[side + "_witness_cell"]
        future = [F(v) for v in row["actual_future_S_interval"]]
        past = [F(v) for v in row["transported_past_S_interval"]]
        raw = row["raw_integer_row"]
        check(future == [F(raw["ell_S_lo"], 2**112)/ell,
                          F(raw["ell_S_hi"], 2**112)/ell],
              side + " raw CSV value denominator and ell normalization")
        check([F(v) for v in row["future_scaled_time_interval"]] == endpoints and
              endpoints[1]-endpoints[0] == F(1, 4096), side + " printed closed whole-cell endpoints")
        check(past == [future[0]-bound, future[1]+bound],
              side + " whole-cell physical S interval transport")
        check(past[1] < safe if side == "negative" else past[0] > safe,
              side + " printed strict whole-cell scalar margin")
    bracket = [F(v) for v in refinement["open_past_zero_bracket_in_units_of_ell"]]
    check(bracket == [F(-6145, 4096), F(-6143, 4096)] and
          bracket[1]-bracket[0] == F(refinement["bracket_width_in_units_of_ell"]) == F(1, 2048),
          "printed refined open past zero bracket and width")
    check(F(1, 2)/(bracket[1]-bracket[0]) == refinement["width_reduction_factor_from_half_ell_bracket"] == 1024,
          "printed factor compares only coarse and refined past intervals")
    refinement_fragments = [r"s\in[3071/2048,\,6143/4096]&S(-\ell s)<-53",
                            r"s\in[6145/4096,\,3073/2048]&S(-\ell s)>63",
                            r"-\frac{6145\ell}{4096},-\frac{6143\ell}{4096}",
                            r"\frac{\ell}{2048}", r"D_v=2^{112}", r"exactly 1024 times narrower",
                            r"not compare precision with older fine future"]
    check(all(fragment in source for fragment in refinement_fragments),
          "refined cell theorem and comparison scope occur in editable source")

    fragments = [
        r"\frac{264000000001}{13738809600}",
        r"\frac{63\ell^2}{4400}", r"<(33\ell)^2", r"\frac{2\ell^2}{125}",
        r"\frac{33}{20000\ell}+\frac{\ell}{160}",
        r"through \(r=6\)", r"validated background collar \(r\le8\)",
        r"stay in \(r\le5\)", r"\ell S(\ell x)",
        r"Earlier future-contact work~\cite{second} already constructs", 
    ]
    check(all(fragment in source for fragment in fragments),
          "key audited formulas radii and prior-result credit occur in editable source")
    check("/workspace/" not in source and "PENDING" not in source,
          "public source contains no private path or unfinished placeholder")
    for name, path in paths.items():
        need(path.read_bytes() == data[name], "Input changed during display audit")
    report_out = {"schema": "HDBLAST_PAST_FUTURE_MANUSCRIPT_DISPLAY_AUDIT_V1",
                  "status": "PASS_PAST_FUTURE_MANUSCRIPT_DISPLAY_AUDIT",
                  "method": "Standard-library Fraction checks of actual printed constants and decimal enclosures against frozen reports; no producer imports",
                  "finite_checks": len(checks), "check_names": checks,
                  "audit_script_sha256": digest(Path(__file__).read_bytes()),
                  "inputs": {name: {"sha256": digest(raw), "bytes": len(raw)} for name, raw in sorted(data.items())},
                  "displayed_coarse_table": displayed,
                  "mathematical_producers_executed": 0,
                  "analytic_premises_reproved": False}
    out.write_text(json.dumps(report_out, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": report_out["status"], "finite_checks": len(checks)}, sort_keys=True))

if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, OSError, ValueError, KeyError, TypeError) as error:
        print(json.dumps({"status": "FAIL", "witness": str(error)}, sort_keys=True))
        sys.exit(2)
