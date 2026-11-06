"""Finite rational replay from public snapshots; no producer report imports.

Analytic Picard, Cauchy, Neumann-series and energy implications remain the
arguments in the accompanying proof notes, conditional on their premises.
"""
from fractions import Fraction as F
from math import factorial
import hashlib
import json
from pathlib import PurePosixPath
import re

REF = "e66a31259d8cd7831c2c93fffa826f8edfdb4efa"
REPOSITORY = "https://github.com/maldonado-research/HDblast"


def need(condition, message):
    if not condition:
        raise ValueError(message)


def rat(value):
    value = F(value)
    return str(value.numerator) + "/" + str(value.denominator)


def matmul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def public_sources(root):
    pins = json.loads((root / "SOURCE_PINS.json").read_text(encoding="utf-8"))
    need(set(pins) == {"schema", "repository", "ref", "scope", "files"}
         and pins["schema"] == "hdblast-boundary-flow-patch-pins-v1", "Public pin schema mismatch")
    need(pins["ref"] == REF and pins["repository"] == REPOSITORY,
         "Public source repository/ref mismatch")
    files = pins["files"]
    need(isinstance(files, list) and len(files) >= 9,
         "At least nine original public inputs required")
    records = {}
    for row in files:
        need(set(row) == {"bytes", "git_blob_sha1", "local", "path", "raw_url",
                          "ref", "sha256", "url"}, "Unexpected public pin keys")
        source = PurePosixPath(row["path"])
        local = PurePosixPath(row["local"])
        need(not source.is_absolute() and ".." not in source.parts
             and source.as_posix() == row["path"] and "\\" not in row["path"]
             and source.parts[0] == "research", "Unsafe original public path")
        need(local.parts[0] == "inputs" and len(local.parts) == 2
             and local.as_posix() == row["local"] and "\\" not in row["local"]
             and local.name not in records, "Unsafe or repeated input pin")
        need(row["ref"] == REF
             and row["url"] == REPOSITORY + "/blob/" + REF + "/" + source.as_posix()
             and row["raw_url"] == "https://raw.githubusercontent.com/"
                 "maldonado-research/HDblast/" + REF + "/" + source.as_posix(),
             "Immutable public URL mismatch")
        raw = (root / row["local"]).read_bytes()
        need(type(row["bytes"]) is int and row["bytes"] == len(raw),
             "Public source size mismatch: " + local.name)
        need(re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None
             and hashlib.sha256(raw).hexdigest() == row["sha256"],
             "Public source SHA256 mismatch: " + local.name)
        blob = b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw
        need(hashlib.sha1(blob).hexdigest() == row["git_blob_sha1"],
             "Public source Git blob mismatch: " + local.name)
        records[local.name] = row
    actual = {p.relative_to(root).as_posix() for p in (root / "inputs").iterdir()}
    need(actual == {row["local"] for row in files}, "Public input roster mismatch")
    return records


def evaluate(root):
    records = public_sources(root)
    data = json.loads((root / "DATA.json").read_text(encoding="utf-8"))
    need(data["schema"] == "HDBLAST_BOUNDARY_FLOW_PATCH_DATA_V1", "Data schema mismatch")
    flow, metric = data["flow"], data["metric"]
    suites = {"flow": [], "metric": []}
    suite = "flow"

    def require(name, condition, **detail):
        need(condition, "Arithmetic failed: " + suite + ":" + name)
        need(name not in {row["name"] for row in suites[suite]}, "Repeated check name")
        suites[suite].append({"name": name, "status": "PASS", **detail})

    bg = json.loads((root / "inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json").read_text())
    c = F(bg["actual_model"]["c"])
    delta = F(bg["actual_model"]["delta"])
    r = F(flow["outer_radius"])
    n = 1 + r
    rho = F(metric["rho_max"])
    ab = F(bg["background_intervals"]["a"]["exact"][1])
    e = r + F(flow["scalar_center_allowance"])
    s, p, q = (F(flow[key]) for key in
               ("sigma_bound", "sigma_prime_bound", "sigma_second_bound"))
    R, DB = F(flow["vector_bound"]), F(flow["derivative_bound"])
    require("c_between_0_and_3_over_5", 0 < c < F(3, 5))
    require("sigma_bound", F(2,3)+delta*(1+c)+delta*c*e+2*e*e+F(2,3)*e**3 < s)
    require("sigma_prime_bound", delta*c+4*e+2*e*e < p)
    require("sigma_second_bound", 4+4*e < q)
    I = 1/(1-r)
    trace_bound = 4*I
    pb = max(abs(F(x)) for x in bg["background_intervals"]["p"]["exact"])
    Rv = F(2,3)*rho*s*n*n + rho*ab*n**3*trace_bound
    RH = rho*(n*s*(1+r)+s)/3
    Rw = rho*(n*p+2*pb)/2
    Lv = (F(4,3)*rho*s*n + 3*rho*ab*n*n*trace_bound
          + F(2,3)*rho*p*n*n + 4*rho*ab*n**3*I**2)
    LH = rho*(s*(1+r)+n*p*(1+r)+n*s)/3
    Lw = rho*(p+n*q)/2
    require("B_max_below_140", max(Rv, RH, Rw) < R)
    require("DB_star_below_400", max(Lv, LH, Lw) < DB)
    # Polynomial identities use the full tangential trace dimension four.
    require("Bn_zero_identity", F(2,3)*6-4 == 0)
    require("BH_zero_identity", 1-1 == 0)
    require("Bphi_zero_identity", 1-1 == 0)

    def sp(eta):
        return delta*c + 4*eta + 2*eta*eta

    center = F(flow["scalar_center_allowance"])
    require("endpoint_eta_lower", sp(-center) < F(24,100000))
    require("endpoint_eta_upper", sp(center) > F(28,100000))
    require("sigma_prime_monotone_on_tube", 4-4*r > 0)
    ell = F(flow["ell"])
    time_multiple = flow["complex_time_multiple"]
    real_multiple = flow["real_time_multiple"]
    coord = F(flow["coordinate_radius"])
    real_center = F(flow["real_center_radius"])
    inner = F(flow["inner_initial_radius"])
    image = F(flow["inner_image_bound"])
    require("complex_start_plus_flow_inside_outer",
            real_center+4*coord+R*time_multiple*ell < r)
    require("complex_time_circle_fits", real_multiple*ell+ell < time_multiple*ell)
    exp_upper = 1/(1-DB*time_multiple*ell)
    require("Jacobian_below_101_over_100", 0 < DB*time_multiple*ell < 1
            and exp_upper < F(101,100))
    require("inner_flow_inside_00101", exp_upper*inner < image)
    require("inner_image_inverse_inside_003", exp_upper*inner < real_center)
    jet_ceiling = F(flow["jet_ceiling"])
    for m in range(flow["jet_total_order"]+1):
        for k in range(flow["jet_total_order"]+1-m):
            require("jet_m%d_k%d_below_1e126" % (m,k),
                    8*factorial(m)*factorial(k)*(1/ell)**m
                    *(flow["state_coordinates"]/coord)**k < jet_ceiling)
    for name in ("ANGULAR_UNIFORM_LINEAR_CONTACT.md", "NONLINEAR_BOUNDARY_STRAIGHTENING.md",
                 "REGISTERED_EXACT_FAMILY_BINDING.md", "MODEL_AND_ARCHIVE_CONTEXT.md",
                 "INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md", "EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md"):
        raw = (root / "inputs" / name).read_bytes()
        require("source_"+name, hashlib.sha256(raw).hexdigest() == records[name]["sha256"])
    require("actual_collar_source_hash", records["ACTUAL_POTENTIAL_JET_BOUNDS.json"]["sha256"]
            == "fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0")
    require("actual_delta", delta == F(1,1000))
    require("actual_c", c == F(96714321425692,161839258930731))
    for name, expected in (("rho", ["129","131"]), ("eta", ["-9/100000","-1/12500"]),
                           ("a", ["111/1000","14/125"]), ("p", ["-7/50000","-3/25000"])):
        require("actual_"+name, bg["background_intervals"][name]["exact"] == expected)
    flow_exact = {key: rat(value) for key, value in
                  (("B_n",Rv), ("B_H",RH), ("B_phi",Rw), ("DB_n",Lv), ("DB_H",LH),
                   ("DB_phi",Lw), ("exponential_majorant",exp_upper), ("inner_image",image))}
    need(len(suites["flow"]) == 220, "Expected 220 retained flow checks")

    suite = "metric"
    for name, record in sorted(records.items()):
        raw = (root / record["local"]).read_bytes()
        require("source_hash:"+name, hashlib.sha256(raw).hexdigest() == record["sha256"])
        require("source_size:"+name, len(raw) == record["bytes"])
    markers = {
      "ANGULAR_UNIFORM_LINEAR_CONTACT.md": [
          "dy=-\\rho\\,dz", "fifteen metric and one scalar",
          "ranks four and twelve", "not a physical gravitational energy"],
      "QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md": [
          "a=−g^{tt}>0", "c^{ij} ξ_i ξ_j ≥ λ |ξ|²",
          "Every iterate must satisfy it.", "A matrix principal coefficient mixing D and N"],
      "INDEPENDENT_QUASILINEAR_MIXED_PROOF.md": [
          "Timelikeness of", "identity for every boundary state", "fixed splitting is orthogonal"],
      "NONLINEAR_BOUNDARY_STRAIGHTENING.md": [
          "N=\\sqrt{g_{yy}}>0", "U_D=V_D", "g^{za}=0", "every frozen-coefficient iterate"],
      "INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md": [
          "No four tangential harmonic rows are added", "K_ab=sigma(phi) h_ab/6",
          "n(phi)=-sigma'(phi)/2"],
      "REGISTERED_EXACT_FAMILY_BINDING.md": [
          "129<=rho<=131", "dy=-rho dz", "unscaled t is the de Sitter coordinate time"],
    }
    for name, needles in markers.items():
        text = (root / "inputs" / name).read_text(encoding="utf-8")
        for needle in needles:
            require("source_marker:"+name+":"+needle, needle in text)
    eps, rho_min, rho_max = (F(metric[key]) for key in ("epsilon","rho_min","rho_max"))
    require("actual_background_rho_box", [F(v) for v in bg["background_intervals"]["rho"]["exact"]]
            == [rho_min,rho_max])
    require("actual_background_domain", F(bg["domain"]["T"]) == F(3,10**6)
            and F(bg["domain"]["Z"]) == F(8,10**6))
    require("actual_background_rectangle", bg["domain"]["verified_enlargement"] == "[0,T] x [0,Z]")
    inverse_bound = 1/(1-eps)
    inverse_error = eps/(1-eps)
    lam_star = (1-inverse_error)/rho_max**2
    mu_star = inverse_bound/rho_min**2
    require("inverse_op_norm", inverse_bound == F(40,39))
    require("inverse_error_norm", inverse_error == F(1,39))
    require("inverse_margins", lam_star == F(38,39*131**2))
    require("time_and_spatial_margin", lam_star > F(metric["safe_coercivity"]), exact=rat(lam_star))
    require("inverse_upper", mu_star < F(metric["safe_inverse_upper"]), exact=rat(mu_star))
    require("coefficient_block_sum_safe", 3*mu_star < 18000)
    require("coefficient_entry_sum_safe", 25*mu_star < 18000)
    require("signature_segment_nonzero", 1-eps > 0)
    require("tangential_invertible_margin", 1-eps == F(39,40))
    r = F(3,1000)
    radial_error = 2*r+r*r
    embedding = max(radial_error,r)+r
    require("state_box_radial_error", radial_error == F(6009,10**6))
    require("state_box_embedding_exact", embedding == F(9009,10**6))
    require("state_box_embeds_full_metric_patch", embedding < eps)
    b = image
    ratio = rho_max/rho_min
    literal_embedding = max(2*b+b*b,ratio**2*b)+ratio*b
    require("literal_inner_state_radial_branch", 2*b+b*b > ratio**2*b)
    # This actual computed inequality is the rehashed corruption-control target.
    require("literal_inner_state_embeds_full_metric_patch", literal_embedding < eps,
            bound=rat(literal_embedding))
    require("literal_inner_state_embedding_exact", literal_embedding == F(3930215929,1290000000000))
    b_outer = F(flow["outer_radius"])
    outer_rayleigh = (2*b_outer+b_outer*b_outer+ratio**2*b_outer+2*ratio*b_outer)/2
    require("outer_componentwise_box_not_sufficient", outer_rayleigh > eps,
            counterexample_rayleigh=rat(outer_rayleigh))
    require("radial_root_lower", (1-eps)**2 < 1-eps)
    require("radial_root_upper", 1+eps < (1+eps)**2)
    double_factorial = 1
    derivatives = {}
    for q in range(1,17):
        if q > 1:
            double_factorial *= 2*q-3
        exact = F(double_factorial,2**q)*F(40,39)**(2*q-1)
        recurrence = F(1,2)*F(40,39) if q == 1 else (
            F(2*q-3,2)*F(40,39)**2*F(derivatives[str(q-1)]))
        require("root_chart_derivative_%02d" % q, exact == recurrence and exact > 0, bound=rat(exact))
        derivatives[str(q)] = rat(exact)
        if q <= 3:
            require("root_chart_derivative_small_%02d" % q, exact < 1)
    require("forward_chart_first_derivative", 2*(1+eps) == F(41,20))
    require("physical_root_first_derivative", F(20,39)/rho_min**2 < F(1,32000))
    T, L, Z = (F(metric[key]) for key in ("candidate_time","initial_support","collar_width"))
    require("candidate_time", T == 13*ell/8)
    require("candidate_time_less_half", T < F(1,2))
    require("cosh_geometric_bound", 1/(1-F(1,2)**2/2) == F(8,7))
    require("cosh_squared_less_two", F(8,7)**2 < 2)
    for k in range(17):
        require("cosh_factorial_majorant_%02d" % k, factorial(2*k) >= 2**k)
    require("energy_natural_lower", F(1,2)*F(metric["safe_coercivity"]) == F(1,36000))
    require("energy_fixed_lower", F(1,2)*F(metric["safe_coercivity"])*F(1,2)
            == F(metric["fixed_energy_lower"]) == F(1,72000))
    require("energy_state_coefficient_lower", F(1,2) > F(1,36000))
    require("energy_upper", mu_star < 1)
    causal_q = (1+eps)/(1-eps)
    require("causal_squared_speed", causal_q == F(41,39))
    require("causal_radial_speed_less_11_10", causal_q < F(11,10)**2)
    require("causal_conservative_speed", F(11,10) < F(metric["conservative_radial_speed"]) == 2)
    require("radial_support_margin", L+2*T == 23*ell/4 and L+2*T < Z)
    require("wider_background_margin", T+Z == 61*ell/8 and T+Z < 8*ell)
    bad = [[F(0) for _ in range(5)] for _ in range(5)]
    for i, value in enumerate([1,-1,-1,1,1]):
        bad[i][i] = F(value)
    bad[1][2] = bad[2][1] = F(2)
    require("timelike_dt_counterexample", bad[1][1] < 0 and bad[0][0] > 0 and bad[2][2] < 0)
    block = [[bad[1][1],bad[1][2]], [bad[2][1],bad[2][2]]]
    require("counterexample_tx_positive_eigen", matmul(block,[[F(1)],[F(1)]]) == [[F(1)],[F(1)]])
    require("counterexample_tx_negative_eigen", matmul(block,[[F(1)],[F(-1)]]) == [[F(-3)],[F(3)]])
    gamma = identity(5)
    gamma[1][1] = F(-1)
    gamma[0][1] = gamma[1][0] = eps
    inv = identity(5)
    inv[0][0] = 1/(1+eps*eps)
    inv[1][1] = -1/(1+eps*eps)
    inv[0][1] = inv[1][0] = eps/(1+eps*eps)
    require("bulk_cross_example_exact_inverse", matmul(gamma,inv) == identity(5))
    require("bulk_cross_not_zero", inv[0][1] != 0)
    require("bulk_radial_root_not_inverse_lapse", inv[0][0] != 1/gamma[0][0])
    temporal_nonzero = set(range(4,16))
    normal_nonzero = set(range(4))
    require("four_plus_twelve_fields", len(normal_nonzero) == 4 and len(temporal_nonzero) == 12)
    require("orthogonal_DN_flux_symbolic_support", not temporal_nonzero.intersection(normal_nonzero))
    need(len(suites["metric"]) == 119, "Expected 119 retained metric checks")
    names = [suite+":"+row["name"] for suite, rows in suites.items() for row in rows]
    return {
        "schema": "hdblast-boundary-flow-patch-v1",
        "status": "PASS", "named_check_count": len(names), "named_checks": names,
        "suite_counts": {suite: len(rows) for suite, rows in suites.items()}, "checks": suites,
        "source_count": len(records),
        "exact": {"flow": flow_exact, "metric": {
            "epsilon":rat(eps), "rho_ratio":rat(ratio), "inverse_error":rat(inverse_error),
            "lambda_star":rat(lam_star), "mu_star":rat(mu_star),
            "inner_metric_embedding":rat(literal_embedding),
            "outer_box_counterexample_rayleigh":rat(outer_rayleigh),
            "energy_fixed_lower":rat(F(metric["fixed_energy_lower"])),
            "causal_speed_squared":rat(causal_q), "root_chart_derivative_bounds":derivatives}},
        "scope": data["scope"],
    }
