"""Exact display/report cross-check; no mathematical producer is imported."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPORT_SHA = "0955ee307c1df84f6cfd7868a0e9ef5e043a176b8285b83ca52a0c48cd9a5564"
PROOF_SHA = "0af17fe7d9ae2a9409821a05520c69cdaf566acdceade0619ab9454b094907a4"

def need(ok, message):
    if not ok:
        raise ValueError(message)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["report", "optimized", "proof", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    need(out.is_absolute() and out.resolve() == out and not out.exists() and
         not out.is_relative_to(HERE), "Fresh canonical external output required")
    need(out.parent.is_dir() and all(not p.is_symlink() for p in [out, *out.parents]),
         "Regular output ancestry required")
    for path in [args.report, args.optimized, args.proof]:
        need(path.is_file() and all(not p.is_symlink() for p in [path, *path.absolute().parents]),
             "Regular frozen input required")
    raw = args.report.read_bytes()
    report = json.loads(raw)
    tex = HERE / "CHARACTERISTIC_ODD_TRACE.tex"
    source_bytes = tex.read_bytes()
    source = source_bytes.decode()
    checks = []
    def check(ok, label):
        need(ok, label)
        checks.append(label)
    def value(path):
        obj = report
        for part in path.split("."):
            obj = obj[part]
        return obj
    def scalar(path):
        return F(value(path))
    def pair(path):
        return tuple(map(F, value(path)))
    check(digest(raw) == REPORT_SHA, "frozen primary report pin")
    check(args.optimized.read_bytes() == raw, "primary normal and optimized byte equality")
    check(digest(args.proof.read_bytes()) == PROOF_SHA, "frozen primary continuous proof pin")
    check(report["proof_sha256"] == PROOF_SHA and report["finite_guards"] == 51 and
          len(report["check_names"]) == 51 and
          report["status"] == "PASS_CONDITIONAL_SHARP_ODD_TRACE_AND_PAST_WAVEFORM",
          "actual conditional report and proof join")
    ell = scalar("ell")
    T = scalar("T")
    m = scalar("majorants.m")
    M = scalar("majorants.M")
    B = scalar("majorants.B")
    gamma = scalar("majorants.Gamma")
    G = scalar("majorants.physical_gain_upper")
    check(ell == F(1,10**6) and T == 3*ell, "registered width and complete elapsed horizon")
    check(m == F(1,1000) and M == 103000 and B == gamma == 300,
          "printed complete coefficient majorants")
    check(G == F(25000,129**3) < F(3,250), "printed physical scalar gain ceiling")
    r = M*T*T/2
    s = B*T*T/(2*m)
    d = (1-r)*(1-s)-T*r/m
    check(r == scalar("majorants.r") and s == scalar("majorants.s") and
          0 < d == scalar("majorants.d"), "positive majorant determinant and both losses")
    A = scalar("majorants.A")
    V = scalar("majorants.V")
    check(A == ell*(1-s+T/m)/d and V == ell/(m*d), "printed velocity-only field and boundary-speed bounds")
    check((1-r)*A-T*V == ell and -(r/m)*A+(1-s)*V == ell/m,
          "both complete positive majorant rows retained")
    I = 3*M*ell*ell/4+M*V*T*T/4+M*M*A*T**3/6
    check(I == scalar("majorants.Istar"), "velocity moment boundary lift and bulk scattering terms")
    E = (T*I/m+F(35,4)*V*T**4/4)/(1-s)
    check(E == scalar("majorants.error_velocity"), "full boundary-error forcing history")
    R = I+(1+gamma*T)*E
    check(R == scalar("majorants.trace_error"), "actual characteristic normal-trace error")
    memory = (1+gamma*T)*ell/m
    check(memory == scalar("majorants.kernel_memory"), "constant inertia memory includes spatial conjugation")
    full = 2*G*(1+memory+R)
    late = 2*G*(F(1,10**9)+memory+R)
    check(full == scalar("difference.derived_global_upper") < F(1,40),
          "actual uniform scalar difference below one fortieth")
    check(late == scalar("difference.derived_late_upper") < F(1,40000),
          "actual late scalar difference below one forty-thousandth")
    check(pair("difference.late_scaled_interval") == (F(49,25),F(3)) and
          scalar("difference.safe_strict_global_upper") == F(1,40) and
          scalar("difference.safe_strict_late_upper") == F(1,40000),
          "printed ceilings and their precise domains")
    exponent = 1/((F(49,25)-1)*(2-F(49,25)))-4
    check(exponent == F(529,24) > 22 and F(8,3)**22 > 10**9,
          "analytic late bump exponential bound without sampling")
    desired = {"central_left":(F(74999,50000),F(1,40),-F(3,2),"upper"),
               "central_right":(F(149999,100000),F(1,40),F(1),"lower"),
               "late_left":(F(49,25),F(1,40000),F(1,1000),"lower"),
               "late_right":(F(197,100),F(1,40000),-F(29,100000),"upper")}
    for name, (x, error, bound, side) in desired.items():
        row = report["selected_point_transfers"][name]
        f = tuple(map(F,row["future_S"]))
        p = tuple(map(F,row["past_S"]))
        check(F(row["x"]) == x and F(row["difference_upper"]) == error and
              p == (f[0]-error,f[1]+error), "same-clock actual interval widening:"+name)
        check(p[1] < bound if side == "upper" else p[0] > bound,
              "printed strict past sign margin:"+name)
    zeros = [tuple(map(F,p)) for p in report["past_zero_brackets_in_units_of_ell"]]
    check(zeros == [(F(-149999,100000),F(-74999,50000)),(F(-197,100),F(-49,25))] and
          zeros[1][1] < zeros[0][0] < zeros[0][1], "two distinct correctly ordered open past brackets")
    check(zeros[0][1]-zeros[0][0] == scalar("central_width_in_units_of_ell") == F(1,100000) and
          F(1,2048)/(zeros[0][1]-zeros[0][0]) == F(3125,64) ==
          scalar("width_reduction_factor_against_previous_ell_over_2048_past_bracket"),
          "central width comparison only to preceding past bracket")
    tail_future = pair("tail.future_actual_S")
    tail_past = pair("tail.past_actual_S")
    check(tail_past == (tail_future[0]-F(1,40000),tail_future[1]+F(1,40000)) and
          -F(1,1500) < tail_past[0] <= tail_past[1] < -F(1,4000),
          "uniform negative past tail with displayed strict margins")
    check(pair("tail.scaled_positive_elapsed_interval") == (F(2),F(3)) and
          pair("tail.readable_strict_past_interval") == (-F(1,1500),-F(1,4000)),
          "whole-tail domain and coefficient rather than angular-field sign")
    hbound = 2*ell*(1+memory+R)
    check(2*hbound*hbound == scalar("quadratic_norm.deltaH_L2_square_upper") < (3*ell)**2 and
          scalar("quadratic_norm.safe_deltaH_L2_strict_upper") == 3*ell,
          "incoming trace to derivative norm uses prefront support")
    cross = 18*ell+9*ell*ell
    check(cross == scalar("quadratic_norm.H_squared_integral_difference_upper") == F(18000009,10**12),
          "quadratic integral retains its full cross term")
    cpast = pair("quadratic_norm.past_C_interval")
    check(F(606*10**9) < cpast[0] <= cpast[1] < F(744*10**9) and
          pair("quadratic_norm.readable_strict_past_C_interval") == (F(606*10**9),F(744*10**9)),
          "printed past first-field quadratic norm interval")
    check(report["quadratic_norm"]["physical_energy_claimed"] is False and
          report["quadratic_norm"]["finite_amplitude_family_law_claimed"] is False,
          "quadratic norm is not physical energy or a proved finite-family law")
    required = [r"\sup_{0\le t\le3\ell}|S(-t)-S(t)|", r"\frac1{40}",
                r"\sup_{49\ell/25\le t\le3\ell}|S(-t)-S(t)|", r"\frac1{40000}",
                r"(\psi,\psi_t)|_0=(0,h)", r"-2G\cosh^{-3/2}t",
                r"d=(1-r)(1-s)-Tr/m", r"\frac34M\ell^2", r"a_0V_*T^4/4",
                r"(k'+\gamma k)*h", r"e^{-529/24}",
                r"(-1.49999\ell,-1.49998\ell)", r"(-1.97\ell,-1.96\ell)",
                r"-\frac1{1500}<S(-t)<-\frac1{4000}", r"3125/64",
                r"18\ell+9\ell^2=\frac{18000009}{10^{12}}", r"\frac{\mathcal K}{\ell}I_-",
                r"6.06\,10^{11}<\mathcal C_-<7.44\,10^{11}",
                "conditional exact nonlinear future families", "qualitative small-amplitude conclusions",
                "No common", "effective parameter radius"]
    for index, token in enumerate(required):
        check(token in source, "critical displayed equation or scope binding:"+str(index+1))
    check("/workspace/" not in source and "PENDING" not in source and "PRIVATE" not in source,
          "public-safe complete manuscript source")
    result = {"schema":"HDBLAST_ODD_TRACE_MANUSCRIPT_DISPLAY_AUDIT_V1",
              "status":"PASS_ODD_TRACE_MANUSCRIPT_DISPLAY_AUDIT",
              "method":"stdlib Fraction; exact printed constants, actual report intervals and source tokens; no producer imports",
              "finite_checks":len(checks),"check_names":checks,
              "audit_script_sha256":digest(Path(__file__).read_bytes()),
              "inputs":{"source":{"sha256":digest(source_bytes)},
                        "primary_report":{"sha256":digest(raw)},
                        "optimized_report":{"sha256":digest(args.optimized.read_bytes())},
                        "primary_proof":{"sha256":digest(args.proof.read_bytes())}},
              "mathematical_producers_imported_or_executed":0,
              "numerical_PDE_or_source_queries":0}
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"finite_checks":len(checks)},sort_keys=True))

if __name__ == "__main__":
    try:
        main()
    except (ValueError,OSError,KeyError,TypeError) as error:
        print(json.dumps({"status":"FAIL","witness":str(error)},sort_keys=True))
        sys.exit(2)
