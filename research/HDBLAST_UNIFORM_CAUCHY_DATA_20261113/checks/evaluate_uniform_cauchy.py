"""Authenticate immutable model inputs and replay finite uniform-Cauchy ledgers."""
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import factorial
from pathlib import PurePosixPath
import re

REF = 'a41940ed60dbda4808ea7286890443a931a44a74'
REPOSITORY = 'https://github.com/maldonado-research/HDblast'
SOURCE_PATHS = {'ACTUAL_POTENTIAL_JET_BOUNDS.json': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json', 'BOUNDARY_FLOW_DATA.json': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/DATA.json', 'BOUNDARY_FLOW_MANIFEST.json': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/MANIFEST.json', 'BOUNDARY_FLOW_PROOF.md': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md', 'CANONICAL_CONSTRAINTS.md': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/inputs/CANONICAL_CONSTRAINTS.md', 'COORDINATE_DICTIONARY.md': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/COORDINATE_DICTIONARY.md', 'CYLINDER_EMBEDDING.md': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md', 'GAUGE_CAUCHY_COMPLETION_06G.md': 'research/HDBLAST_INITIAL_TRANSPORT_20261106/inputs/GAUGE_CAUCHY_COMPLETION_06G.md', 'INITIAL_CONSTRAINT_DATA.json': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/DATA.json', 'INITIAL_CONSTRAINT_MANIFEST.json': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/MANIFEST.json', 'INITIAL_CONSTRAINT_README.md': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/README.md', 'INITIAL_CONSTRAINT_SOURCE_PINS.json': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/SOURCE_PINS.json', 'INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md': 'research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md', 'INITIAL_TRANSPORT.md': 'research/HDBLAST_INITIAL_TRANSPORT_20261106/INITIAL_TRANSPORT.md', 'INITIAL_TRANSPORT_DATA.json': 'research/HDBLAST_INITIAL_TRANSPORT_20261106/DATA.json', 'INITIAL_TRANSPORT_MANIFEST.json': 'research/HDBLAST_INITIAL_TRANSPORT_20261106/MANIFEST.json', 'NONLINEAR_BOUNDARY_STRAIGHTENING.md': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/inputs/NONLINEAR_BOUNDARY_STRAIGHTENING.md', 'NONLINEAR_CONSTRAINT_BALL.md': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/NONLINEAR_CONSTRAINT_BALL.md', 'PHYSICAL_CORRECTION_MAPS.md': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/PHYSICAL_CORRECTION_MAPS.md', 'PINNED_SOURCE_JETS.py': 'research/HDBLAST_WEIGHTED_REGULARITY_20261106/derive_source_jets.py', 'REGISTERED_EXACT_FAMILY_BINDING.md': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/inputs/REGISTERED_EXACT_FAMILY_BINDING.md', 'SEED_ENVELOPE.md': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/SEED_ENVELOPE.md', 'WEIGHTED_EMBEDDING_PRODUCTS.md': 'research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/WEIGHTED_EMBEDDING_PRODUCTS.md'}
MODULES = {"arithmetic":"cauchy_arithmetic", "chart":"chart_chain", "algebra":"cauchy_algebra"}

def need(condition, message):
    if not condition:
        raise ValueError(message)

def load_json(raw):
    def unique(items):
        result = {}
        for key, value in items:
            need(key not in result, "Duplicate JSON key: " + key)
            result[key] = value
        return result

    def nonfinite(value):
        raise ValueError("Nonfinite JSON constant: " + value)

    return json.loads(raw.decode("utf-8"), object_pairs_hook=unique,
                      parse_float=str, parse_constant=nonfinite)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def rational_expression(node, environment):
    """Read only elementary literal assignments from the pinned producer AST."""
    if isinstance(node, ast.Constant) and type(node.value) is int:
        return node.value
    if isinstance(node, ast.Name) and node.id in environment:
        return environment[node.id]
    if isinstance(node, (ast.List, ast.Tuple)):
        return [rational_expression(item, environment) for item in node.elts]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = rational_expression(node.operand, environment)
        return -value if isinstance(node.op, ast.USub) else value
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "F":
        need(not node.keywords and len(node.args) in (1, 2), "Unexpected producer Fraction call")
        return F(*(rational_expression(item, environment) for item in node.args))
    if isinstance(node, ast.BinOp):
        left, right = rational_expression(node.left, environment), rational_expression(node.right, environment)
        if isinstance(node.op, ast.Add): return left + right
        if isinstance(node.op, ast.Sub): return left - right
        if isinstance(node.op, ast.Mult): return left * right
        if isinstance(node.op, ast.Div): return F(left) / F(right)
        if isinstance(node.op, ast.Pow):
            need(type(right) is int and 0 <= right <= 12, "Unexpected producer power")
            return left ** right
    raise ValueError("Unexpected producer literal syntax")

def convolution(left, right):
    out = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right): out[i+j] += a*b
    return out

def public_sources(root, check):
    raw_pins = (root / "SOURCE_PINS.json").read_bytes()
    pins = load_json(raw_pins)
    need(isinstance(pins, dict) and set(pins) == {"schema", "repository", "ref", "files"},
         "Unexpected public pin keys")
    need(pins["schema"] == "hdblast-uniform-cauchy-pins-v1",
         "Public pin schema mismatch")
    need(pins["repository"] == REPOSITORY and pins["ref"] == REF,
         "Public source repository/ref mismatch")
    files = pins["files"]
    need(isinstance(files, list) and len(files) == len(SOURCE_PATHS),
         "Original public input count differs")
    names = set()
    by_local = {}
    for row in files:
        need(isinstance(row, dict) and set(row) ==
             {"local", "path", "ref", "url", "bytes", "sha256", "git_blob_sha1"},
             "Unexpected source pin keys")
        need(isinstance(row["local"], str) and isinstance(row["path"], str),
             "Source paths must be strings")
        local, source = PurePosixPath(row["local"]), PurePosixPath(row["path"])
        need(not local.is_absolute() and len(local.parts) == 2 and local.parts[0] == "inputs"
             and local.name not in names and "\\" not in row["local"]
             and local.as_posix() == row["local"] and ".." not in local.parts,
             "Unsafe or repeated local source pin")
        need(local.name in SOURCE_PATHS and row["path"] == SOURCE_PATHS[local.name],
             "Original source roster differs")
        need(not source.is_absolute() and source.parts[0] == "research"
             and "\\" not in row["path"] and ".." not in source.parts
             and source.as_posix() == row["path"], "Unsafe original source path")
        need(row["ref"] == REF and row["url"] ==
             REPOSITORY + "/blob/" + REF + "/" + source.as_posix(),
             "Immutable public URL mismatch")
        raw = (root / local).read_bytes()
        check("source_size:" + local.name, type(row["bytes"]) is int
              and row["bytes"] >= 0 and len(raw) == row["bytes"])
        check("source_sha256:" + local.name, isinstance(row["sha256"], str)
              and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None
              and digest(raw) == row["sha256"])
        blob = b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw
        check("source_git_blob:" + local.name, isinstance(row["git_blob_sha1"], str)
              and re.fullmatch(r"[0-9a-f]{40}", row["git_blob_sha1"]) is not None
              and hashlib.sha1(blob).hexdigest() == row["git_blob_sha1"])
        names.add(local.name)
        by_local[local.name] = row
    need({item.name for item in (root / "inputs").iterdir()} == names,
         "Public input roster mismatch")
    need([row["local"] for row in files] == sorted(row["local"] for row in files),
         "Public input order differs")

    provenance = load_json((root / "INPUT_PROVENANCE.json").read_bytes())
    need(isinstance(provenance, dict) and set(provenance) ==
         {"schema", "repository", "ref", "snapshot_method", "snapshot_producers_executed",
          "source_pins_sha256", "snapshot_count", "files"}, "Input provenance keys differ")
    check("provenance_fixed_base", provenance["schema"] ==
          "HDBLAST_UNIFORM_CAUCHY_INPUT_PROVENANCE_V1"
          and provenance["repository"] == REPOSITORY and provenance["ref"] == REF)
    check("provenance_snapshot_only", provenance["snapshot_method"] ==
          "git cat-file of immutable base blobs" and provenance["snapshot_producers_executed"] is False)
    check("provenance_source_index", provenance["source_pins_sha256"] == digest(raw_pins)
          and type(provenance["snapshot_count"]) is int
          and provenance["snapshot_count"] == len(files))
    check("provenance_same_blob_roster", provenance["files"] == [
        {"local": row["local"], "git_blob_sha1": row["git_blob_sha1"], "sha256": row["sha256"]}
        for row in files])

    ancestors = {
        "INITIAL_CONSTRAINT_MANIFEST.json": {
            "INITIAL_CONSTRAINT_README.md":"README.md", "INITIAL_CONSTRAINT_DATA.json":"DATA.json",
            "INITIAL_CONSTRAINT_SOURCE_PINS.json":"SOURCE_PINS.json", "NONLINEAR_CONSTRAINT_BALL.md":"NONLINEAR_CONSTRAINT_BALL.md",
            "PHYSICAL_CORRECTION_MAPS.md":"PHYSICAL_CORRECTION_MAPS.md", "WEIGHTED_EMBEDDING_PRODUCTS.md":"WEIGHTED_EMBEDDING_PRODUCTS.md",
            "SEED_ENVELOPE.md":"SEED_ENVELOPE.md"},
        "INITIAL_TRANSPORT_MANIFEST.json": {"INITIAL_TRANSPORT.md":"INITIAL_TRANSPORT.md", "INITIAL_TRANSPORT_DATA.json":"DATA.json"},
        "BOUNDARY_FLOW_MANIFEST.json": {"BOUNDARY_FLOW_PROOF.md":"BOUNDARY_FLOW_PROOF.md", "BOUNDARY_FLOW_DATA.json":"DATA.json",
                                         "COORDINATE_DICTIONARY.md":"COORDINATE_DICTIONARY.md"}}
    for manifest, members in ancestors.items():
        previous=load_json((root/"inputs"/manifest).read_bytes())
        need(isinstance(previous,dict) and isinstance(previous.get("files"),dict), "Ancestor manifest schema")
        for local, member in members.items():
            row=by_local[local]
            check("ancestor_member:"+local, previous["files"].get(member)=={"bytes":row["bytes"],"sha256":row["sha256"]})
    return len(names)


def source_models(root, check):
    actual = load_json((root / "inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json").read_bytes())
    model = actual["actual_model"]
    potential = [F(item) for item in model["U_coefficients"]]
    expected = [F(-2,27), F(0), F(14,9), F(50,27), F(-1,6), F(-4,9), F(-2,27)]
    check("actual_potential_coefficients", potential == expected)
    superpotential = [F(item) for item in model["W_super_coefficients"]]
    derivative = [i*superpotential[i] for i in range(1, len(superpotential))]
    square = convolution(superpotential, superpotential)
    derivative_square = convolution(derivative, derivative)
    rebuilt = [F(-2,3)*item for item in square]
    for i, item in enumerate(derivative_square): rebuilt[i] += item/2
    check("actual_superpotential_identity", rebuilt == potential)
    check("actual_shell_member", model["eta"] == "phi-1"
          and F(model["delta"]) == F(1,1000)
          and F(model["c"]) == F(96714321425692,161839258930731))
    check("actual_seed_scale", F(actual["domain"]["ell"]) == F(1,10**6))
    background = actual["background_intervals"]
    check("actual_background_rho", [F(item) for item in background["rho"]["exact"]] == [129,131])
    check("actual_background_nonzero_scalar_gradient",
          [F(item) for item in background["p"]["exact"]] == [F(-7,50000),F(-3,25000)])
    check("actual_background_warp_derivative",
          [F(item) for item in background["a"]["exact"]] == [F(111,1000),F(14,125)])

    tree = ast.parse((root / "inputs/PINNED_SOURCE_JETS.py").read_text(encoding="utf-8"))
    environment = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in {"N", "SIG", "ELL", "KAPPA", "VCOEFF"}:
                environment[name] = rational_expression(node.value, environment)
    check("pinned_producer_orders", environment["N"] == 16 and environment["SIG"] == 14)
    check("pinned_producer_scale", environment["ELL"] == F(1,10**6)
          and environment["KAPPA"] == F(1,10**8))
    check("pinned_producer_same_model", environment["VCOEFF"] == potential)
    functions = {node.name for node in tree.body if isinstance(node, ast.FunctionDef)}
    check("pinned_producer_api", {"background", "bump_bounds", "compose_bump", "source_jets"} <= functions)

    initial=load_json((root/"inputs/INITIAL_CONSTRAINT_DATA.json").read_bytes())
    check("initial_constraint_schema",initial["schema"]=="HDBLAST_NONLINEAR_CONSTRAINT_BALL_DATA_V1")
    previous=initial["arithmetic"]
    check("retained_initial_family_constants",previous["K_exponent"]==6400 and previous["amplitude_exponent"]==-6402
          and previous["potential_radius_exponent"]==-6404 and previous["B_exponent"]==200)
    check("retained_physical_maps",previous["physical_metric_map_exponent"]==256 and previous["physical_momentum_map_exponent"]==253
          and previous["metric_correction_E_exponent"]==30 and previous["momentum_correction_E_exponent"]==27)
    check("retained_products_and_pointwise",previous["product_exponent"]==17 and previous["low_jet_exponent"]==12
          and previous["physical_multiplier_exponent"]==226 and previous["metric_guard_denominator"]==1024)
    flow=load_json((root/"inputs/BOUNDARY_FLOW_DATA.json").read_bytes())
    check("retained_flow_schema",flow["schema"]=="HDBLAST_BOUNDARY_FLOW_PATCH_DATA_V1")
    f=flow["flow"]
    check("retained_flow_member",F(f["ell"])==F(1,10**6) and f["state_coordinates"]==16 and f["jet_total_order"]==18)
    check("retained_flow_patches",F(f["outer_radius"])==F(1,100) and F(f["real_center_radius"])==F(3,1000)
          and F(f["coordinate_radius"])==F(1,1000) and F(f["inner_initial_radius"])==F(1,1000))
    check("retained_flow_jets",F(f["derivative_bound"])==400 and F(f["vector_bound"])==140
          and F(f["jet_ceiling"])==10**126 and f["real_time_multiple"]==6 and f["complex_time_multiple"]==8)
    transport=load_json((root/"inputs/INITIAL_TRANSPORT_DATA.json").read_bytes())
    check("retained_transport_schema",transport["schema"]=="HDBLAST_INITIAL_TRANSPORT_DATA_V1")
    t=transport["bounds"]
    check("retained_linear_frame_maps",t["conversion_decade"]==3000 and t["normalization_decade"]==210
          and t["linear_Leibniz_decade"]==6 and t["flow_total_order"]==18)
    check("retained_first_cauchy_direction",t["final_first_decade"]==5200)
    anchors={
        "CANONICAL_CONSTRAINTS.md":["H   = R/2 + (2/det γ)","M_i = −2 D_j π^j_i"],
        "REGISTERED_EXACT_FAMILY_BINDING.md":["x_e=x0+e v_seed+A u(e)","scalar and scalar-momentum sectors are prescribed affinely"],
        "BOUNDARY_FLOW_PROOF.md":["It is not the physical bulk radial ODE","m+k<=18","F_z(0)=0"],
        "INITIAL_TRANSPORT.md":["10^3000","partial_t A|t=0=0","proper-y"],
        "CYLINDER_EMBEDDING.md":["E3=10^7"],
        "INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md":["R-|d phi|^2-2U","varpi=2 sqrt(det gamma) u(phi)"]}
    for name,phrases in anchors.items():
        text=(root/"inputs"/name).read_text()
        check("source_model_anchor:"+name,all(phrase in text for phrase in phrases))
    return {"potential":[str(x) for x in potential],"producer_order":environment["N"],
            "initial_K_exponent":6400,"amplitude_exponent":-6402,"potential_radius_exponent":-6404,
            "flow_total_order":18,"physical_multiplier_exponent":226,
            "linear_identification_exponent":3000,"normalization_jet_exponent":210,"first_cauchy_exponent":5200}


def load_module(root, name):
    spec = importlib.util.spec_from_file_location("uniform_cauchy_" + name,
                                                  root / "checks" / (name + ".py"))
    need(spec is not None and spec.loader is not None, "Evaluator loader unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def evaluate(root):
    data=load_json((root/"DATA.json").read_bytes())
    need(isinstance(data,dict) and set(data)=={"schema","scope",*MODULES}
         and data["schema"]=="HDBLAST_UNIFORM_CAUCHY_DATA_V1" and isinstance(data["scope"],str) and data["scope"],"Data schema mismatch")
    suites={name:[] for name in ("sources","source_models",*MODULES,"interfaces")}
    def checker(suite):
        def check(name,condition):
            need(isinstance(name,str) and name,"Nonempty check name required")
            need(bool(condition),"Arithmetic failed: "+suite+":"+name)
            need(name not in suites[suite],"Duplicate check: "+suite+":"+name)
            suites[suite].append(name)
        return check
    count=public_sources(root,checker("sources"))
    exact={"source_models":source_models(root,checker("source_models"))}
    for suite,name in MODULES.items():
        need(isinstance(data[suite],dict),"Claims must be objects: "+suite)
        module=load_module(root,name)
        exact[suite]=module.evaluate(data[suite],checker(suite))
        need(isinstance(exact[suite],dict) and exact[suite],"Exact suite report required")
    joined=checker("interfaces")
    a=data["arithmetic"]; source=exact["source_models"]
    joined("original_admitted_interval", a["initial_amplitude_exponent"]==source["amplitude_exponent"]
           and a["initial_K_exponent"]==source["initial_K_exponent"] and a["potential_radius_exponent"]==source["potential_radius_exponent"]
           and a["seed_exponent"]==200 and F(a["potential_first_derivative_bound"])==F(101,5000)
           and a["potential_second_derivative_factor"]==8)
    joined("same_actual_potential",[F(x) for x in a["potential_coefficients"]]==[F(x) for x in source["potential"]])
    joined("same_graded_weight_member",F(a["weight_kappa"])==F(1,10**8)
           and a["physical_multiplier_exponent"]==source["physical_multiplier_exponent"]
           and a["graded_algebra_exponent"]==100 and a["graded_ball_denominator"]==1024)
    joined("same_fixed_frame_maps",a["raw_frame_exponent"]==source["linear_identification_exponent"]
           and a["normalization_jet_exponent"]==source["normalization_jet_exponent"] and a["linear_multiplier_exponent"]==6)
    joined("same_retained_flow_certificate",a["flow_total_order"]==source["flow_total_order"]
           and a["flow_star_jet_exponent"]==126 and F(a["ell"])==F(1,10**6) and a["flow_generator_Lipschitz"]==400)
    joined("same_cylinder_embedding",a["cylinder_embedding_exponent"]==7)
    joined("same_zero_parameter_first_direction",a["retained_zero_first_derivative_exponent"]==source["first_cauchy_exponent"])
    joined("same_half_action_density",F(data["algebra"]["K_conversion_factor"])==2
           and F(data["algebra"]["normal_scalar_conversion_factor"])==1)
    joined("specified_identical_copy_scope",a["retained_copy_count"]==2
           and F(a["copy_linear_factor"])==2 and F(a["copy_quadratic_factor"])==4)
    exact["interfaces"]={"amplitude_exponent":a["initial_amplitude_exponent"],
                         "completed_Hessian_exponent":a["completed_Hessian_exponent"],
                         "copy_count":a["retained_copy_count"]}
    names=[suite+":"+name for suite,entries in suites.items() for name in entries]
    return {"schema":"hdblast-uniform-cauchy-v1","status":"PASS","named_check_count":len(names),"named_checks":names,
            "suite_counts":{suite:len(entries) for suite,entries in suites.items()},"source_count":count,"exact":exact,"scope":data["scope"]}
