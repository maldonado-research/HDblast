#!/usr/bin/env python3
"""Independent exact TT Fr derivative-order and physical-chain audit.

Only reads frozen source; every acceptance uses explicit checks, including -O.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as S


def need(value, message):
    if not value:
        raise ValueError(message)


def equal(a, b):
    return S.cancel(S.expand(a - b)) == 0


def run(repo):
    source_path = repo / "research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json"
    source_bytes = source_path.read_bytes()
    maps = json.loads(source_bytes)
    row = maps["maps"]["metric_forcing"]["ab"][5]
    need(maps["status"] == "PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION", "SOURCE_STATUS")
    need(len(row) == 78 and row == maps["maps"]["metric_forcing"]["ab"][6], "COMPLETE_IDENTICAL_TT_ROWS")
    names = sorted({name for term in row for name in term["powers"]})
    base_orders = {
        "A": (0, 1), "B": (None, 1), "C": (None, 1),
        "D": (0, 1), "E": (None, 1), "L": (0, 0), "M": (None, 0),
        "s": (1, 0),
    }
    dependencies = {}
    for name in names:
        base = name[:1]
        if base not in base_orders:
            continue
        suffix = name[1:]
        need(set(suffix) <= {"y", "t"}, "RECOGNIZED_PHYSICAL_SUFFIX:" + name)
        orders = base_orders[base]
        dependencies[name] = {
            "physical_suffix_order": len(suffix),
            "master_scaled_order_for_Fr": None if orders[0] is None else orders[0] + len(suffix) + 1,
            "gauge_scaled_order_for_Fr": orders[1] + len(suffix) + 1,
        }
    second = sorted(name for name, item in dependencies.items() if item["physical_suffix_order"] == 2)
    need(second == ["Mtt", "Myt", "Myy"], "ONLY_M_USES_SECOND_METRIC_JETS")
    need(max(item["master_scaled_order_for_Fr"] or 0 for item in dependencies.values()) == 2, "MASTER_C2_SUFFICES")
    need(max(item["gauge_scaled_order_for_Fr"] for item in dependencies.values()) == 3, "GAUGE_C3_SUFFICES")

    ell, rho, a, ay, p, py, vp, vpp = S.symbols("ell rho a ay p py Vp Vpp", nonzero=True)
    m = S.symbols("m0:4")
    # This derivation differentiates exact jet symbols using the background ODE,
    # rather than differentiating selected interval endpoints or error caps.
    radial_rules = {
        rho: -ell * rho**2 * a,
        a: -ell * rho * ay,
        p: -ell * rho * py,
        vp: -ell * rho * vpp * p,
        **{m[i]: m[i + 1] for i in range(3)},
    }
    dr = lambda value: S.expand(sum(S.diff(value, key) * derivative for key, derivative in radial_rules.items()))
    conv = 1 / (ell * rho)
    first = -conv * m[1]
    second_y = conv**2 * m[2] + a * conv * m[1]
    third_y = -conv**3 * m[3] - 3 * a * conv**2 * m[2] + (ay - a**2) * conv * m[1]
    need(equal(-conv * dr(first), second_y), "PHYSICAL_YY")
    need(equal(-conv * dr(second_y), third_y), "PHYSICAL_YYY")
    need(equal(dr(conv), a), "CONVERSION_DERIVATIVE")
    ayy = 2 * a / rho**2 - S.Rational(2, 3) * p * py
    need(equal(dr(-rho**-2 - p**2 / 3), -ell * rho * ayy), "BACKGROUND_AY_CHAIN")
    need(equal(dr(vp), -ell * rho * vpp * p), "BACKGROUND_VP_CHAIN")
    controls = {}
    for name, wrong in {
        "omit_yyy_background_terms": -conv**3 * m[3],
        "wrong_yyy_middle_sign": -conv**3 * m[3] + 3 * a * conv**2 * m[2] + (ay - a**2) * conv * m[1],
        "omit_ay_p_feedback": 2 * a / rho**2,
        "wrong_Vp_radial_sign": ell * rho * vpp * p,
    }.items():
        right = ayy if name == "omit_ay_p_feedback" else (-ell * rho * vpp * p if name == "wrong_Vp_radial_sign" else third_y)
        need(not equal(wrong, right), "NEGATIVE_CONTROL_SURVIVED:" + name)
        controls[name] = "REJECTED_BY_EXACT_IDENTITY"

    # Independent full symbolic substitution verifies that highest third gauge
    # derivatives need not vanish.  No fourth derivative occurs in the formula.
    x, r = S.symbols("x r")
    e = S.symbols("ell", nonzero=True)
    rr, aa, cc, hh = (S.Function(name)(r if name in ("rho", "a") else x) for name in ("rho", "a", "c", "H"))
    u, f, b, d, scalar = (S.Function(name)(x, r) for name in ("u", "f", "b", "d", "s"))
    amps = {
        "A": -4 * u - 2 * S.diff(f, r) / (e * rr),
        "B": S.diff(f, x) / (e * rr) + S.diff(b, r) / e,
        "C": f / (rr * cc) - cc * S.diff(d, r) / e,
        "D": -2 * (u + aa * f + S.diff(b, x) / e),
        "E": cc * S.diff(d, x) / e - b / cc,
        "L": 2 * (u + aa * f + hh * b),
        "M": 2 * d,
    }
    raw = {"ell": e, "rho": rr, "a": aa, "ay": S.Function("ay")(r), "c": cc, "theta": hh, "Vp": S.Function("Vp")(r), "s": scalar}
    for name, amp in amps.items():
        raw[name] = amp
        raw[name + "y"] = -S.diff(amp, r) / (e * rr)
        raw[name + "t"] = S.diff(amp, x) / e
    raw["Myy"] = -S.diff(raw["My"], r) / (e * rr)
    raw["Myt"] = S.diff(raw["My"], x) / e
    raw["Mtt"] = S.diff(raw["Mt"], x) / e
    whole = S.Add(*(S.Rational(term["coefficient"]) * S.prod(raw[key] ** power for key, power in term["powers"].items()) for term in row))
    fr = S.expand(S.diff(whole, r))
    third_f_coefficient = S.factor(fr.coeff(S.diff(f, r, 3)))
    expected_third_f = 8 * f / (3 * rr**2 * cc**2) - S.Rational(16, 3) * aa * d
    need(equal(third_f_coefficient, expected_third_f), "NONZERO_THIRD_GAUGE_PRINCIPAL_COEFFICIENT")
    need(third_f_coefficient != 0, "THIRD_GAUGE_DERIVATIVES_DO_NOT_ALL_CANCEL")
    max_gauge_order = max(len(der.variables) for der in fr.atoms(S.Derivative) if der.expr in (f, b, d))
    max_master_order = max(len(der.variables) for der in fr.atoms(S.Derivative) if der.expr == u)
    need(max_gauge_order == 3 and max_master_order == 2, "DIRECT_SYMBOLIC_DERIVATIVE_ORDER")
    controls["claim_all_gauge_third_derivatives_cancel"] = "REJECTED_BY_NONZERO_EXACT_PRINCIPAL_COEFFICIENT"
    return {
        "status": "PASS_INDEPENDENT_TT_FR_INPUT_DEPENDENCY_REVIEW",
        "source_path": str(source_path),
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "complete_TT_terms": len(row),
        "both_TT_rows_identical": True,
        "raw_source_symbols": names,
        "source_second_metric_jet_symbols": second,
        "maximum_required_master_scaled_derivative_order": 2,
        "maximum_required_gauge_scaled_derivative_order": 3,
        "source_symbol_order_audit": dependencies,
        "physical_conversion": {
            "L": "1/(ell*rho)",
            "dr_L": "a",
            "Myyy": "-L^3*Mrrr - 3*a*L^2*Mrr + (ay-a^2)*L*Mr",
            "Myyt": "(L^2*Mxrr + a*L*Mxr)/ell",
            "Mytt": "-L*Mxxr/ell^2",
            "Fr_metric_or_scalar_symbol_chain": "dr(q_{y^i t^j})=-ell*rho*q_{y^(i+1) t^j}",
            "rho_r": "-ell*rho^2*a",
            "a_r": "-ell*rho*ay",
            "ay_r": "-ell*rho*(2*a/rho^2-(2/3)*p*py)",
            "Vp_r": "-ell*rho*Vpp*p",
            "c_r_theta_r_ell_r": "0",
        },
        "nonzero_frrr_coefficient": "8*f/(3*rho^2*c^2)-16*a*d/3",
        "negative_controls": controls,
        "proof_scope": "Formal source dependency and exact chain identities. Actual enclosure validity remains conditional on the authenticated actual C3 solution certificates. No C3 cap is differentiated. No C4 premise is required for this TT Fr row.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists(), "FRESH_OUTPUT_REQUIRED")
    result = run(args.repo_root.resolve())
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "terms": result["complete_TT_terms"], "master_order": 2, "gauge_order": 3, "controls_rejected": len(result["negative_controls"])}))


if __name__ == "__main__":
    main()
