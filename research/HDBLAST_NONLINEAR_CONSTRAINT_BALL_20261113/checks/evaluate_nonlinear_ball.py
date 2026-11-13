"""Finite replay with immutable original inputs and explicit model binding."""
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import factorial
from pathlib import PurePosixPath
import re

REF = "4c1c47ce4cf7f0106e3a462bfd48618b04ae3bdf"
REPOSITORY = "https://github.com/maldonado-research/HDblast"
SOURCE_PATHS = {
    "ALL_MODE_README.md": "research/HDBLAST_ALL_MODE_INVERSE_20261113/README.md",
    "SCALAR_ALL_MODE.md": "research/HDBLAST_ALL_MODE_INVERSE_20261113/SCALAR_ALL_MODE.md",
    "MOMENTUM_ALL_MODE.md": "research/HDBLAST_ALL_MODE_INVERSE_20261113/MOMENTUM_ALL_MODE.md",
    "MOMENTUM_REGULARITY_ALL_MODE.md": "research/HDBLAST_ALL_MODE_INVERSE_20261113/MOMENTUM_REGULARITY_ALL_MODE.md",
    "WEIGHTED_NORM_BRIDGE.md": "research/HDBLAST_ALL_MODE_INVERSE_20261113/WEIGHTED_NORM_BRIDGE.md",
    "ALL_MODE_DATA.json": "research/HDBLAST_ALL_MODE_INVERSE_20261113/DATA.json",
    "ALL_MODE_MANIFEST.json": "research/HDBLAST_ALL_MODE_INVERSE_20261113/MANIFEST.json",
    "ALL_MODE_SOURCE_PINS.json": "research/HDBLAST_ALL_MODE_INVERSE_20261113/SOURCE_PINS.json",
    "SELECTED_MEMBER.md": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/README.md",
    "SCALAR_REGULARITY.md": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/SCALAR_REGULARITY.md",
    "MOMENTUM_EXTENSION.md": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/MOMENTUM_EXTENSION.md",
    "SOURCE_JETS.md": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/SOURCE_JETS.md",
    "WEIGHTED_CONTINUATION.md": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/CONTINUATION.md",
    "EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md": "research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md",
    "REGISTERED_EXACT_FAMILY_BINDING.md": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/inputs/REGISTERED_EXACT_FAMILY_BINDING.md",
    "CANONICAL_CONSTRAINTS.md": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/inputs/CANONICAL_CONSTRAINTS.md",
    "ACTUAL_POTENTIAL_JET_BOUNDS.json": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json",
    "PINNED_SOURCE_JETS.py": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/derive_source_jets.py",
    "PINNED_REGULARITY_CONSTANTS.py": "research/HDBLAST_WEIGHTED_REGULARITY_20261106/derive_regularity_constants.py",
    "INITIAL_TRANSPORT.md": "research/HDBLAST_INITIAL_TRANSPORT_20261106/INITIAL_TRANSPORT.md",
    "CYLINDER_EMBEDDING.md": "research/HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md",
}
MODULES = {"arithmetic": "ball_arithmetic", "weights": "weight_degree_ledger",
           "seed": "seed_envelope", "algebra": "constraint_algebra"}


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


def public_sources(root, check):
    raw_pins = (root / "SOURCE_PINS.json").read_bytes()
    pins = load_json(raw_pins)
    need(isinstance(pins, dict) and set(pins) == {"schema", "repository", "ref", "files"},
         "Unexpected public pin keys")
    need(pins["schema"] == "hdblast-nonlinear-constraint-ball-pins-v1",
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
          "HDBLAST_NONLINEAR_CONSTRAINT_BALL_INPUT_PROVENANCE_V1"
          and provenance["repository"] == REPOSITORY and provenance["ref"] == REF)
    check("provenance_snapshot_only", provenance["snapshot_method"] ==
          "git cat-file of immutable base blobs" and provenance["snapshot_producers_executed"] is False)
    check("provenance_source_index", provenance["source_pins_sha256"] == digest(raw_pins)
          and type(provenance["snapshot_count"]) is int
          and provenance["snapshot_count"] == len(files))
    check("provenance_same_blob_roster", provenance["files"] == [
        {"local": row["local"], "git_blob_sha1": row["git_blob_sha1"], "sha256": row["sha256"]}
        for row in files])

    old = load_json((root / "inputs/ALL_MODE_MANIFEST.json").read_bytes())
    need(isinstance(old, dict) and old.get("schema") == "HDBLAST_ALL_MODE_INVERSE_MANIFEST_V1"
         and isinstance(old.get("files"), dict), "Retained all-mode manifest schema differs")
    old_members = {"ALL_MODE_README.md": "README.md", "ALL_MODE_DATA.json": "DATA.json",
        "ALL_MODE_SOURCE_PINS.json": "SOURCE_PINS.json", "SCALAR_ALL_MODE.md": "SCALAR_ALL_MODE.md",
        "MOMENTUM_ALL_MODE.md": "MOMENTUM_ALL_MODE.md",
        "MOMENTUM_REGULARITY_ALL_MODE.md": "MOMENTUM_REGULARITY_ALL_MODE.md",
        "WEIGHTED_NORM_BRIDGE.md": "WEIGHTED_NORM_BRIDGE.md"}
    for local, member in old_members.items():
        row = by_local[local]
        check("retained_manifest_member:" + local, old["files"].get(member) ==
              {"bytes": row["bytes"], "sha256": row["sha256"]})
    return len(names)


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

    inverse = load_json((root / "inputs/ALL_MODE_DATA.json").read_bytes())
    check("retained_inverse_schema", inverse["schema"] == "HDBLAST_ALL_MODE_INVERSE_DATA_V1")
    ledger = inverse["arithmetic"]
    check("retained_same_inverse_norms", ledger["combined_Hilbert_inverse_exponent"] == 5000
          and ledger["scalar_final_inverse_exponent"] == 5000
          and ledger["momentum_final_inverse_exponent"] == 2700)
    check("retained_norm_bridge_integer", int(ledger["norm_bridge_integer"]) == 13*8**16*factorial(16)
          and ledger["norm_bridge_decimal_ceiling"] == 29)
    anchors = {
        "ALL_MODE_README.md": ["q^(2j)|nabla_gamma0^j u|^2", "fixed canonical density `rho0`"],
        "EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md": ["X_k=(Ḣ^{k+4}_{q,w} ⊕ Ḣ^{k+2}_{q,w})^Γ",
            "δπ^{ij}=ρ₀τ^{ij}", "A_HN=w²q⁴ P_g* N", "A_MZ=w²q² P_m* Z"],
        "REGISTERED_EXACT_FAMILY_BINDING.md": ["x_e=x0+e v_seed+A u(e)",
            "scalar and scalar-momentum sectors are prescribed affinely"],
        "CANONICAL_CONSTRAINTS.md": ["H   = R/2 + (2/det γ)", "M_i = −2 D_j π^j_i"],
        "CYLINDER_EMBEDDING.md": ["E3=10^7"],
    }
    for name, phrases in anchors.items():
        text = (root / "inputs" / name).read_text(encoding="utf-8")
        check("source_model_anchor:" + name, all(phrase in text for phrase in phrases))
    return {"potential": [str(item) for item in potential], "producer_order": environment["N"],
            "same_inverse_exponent": ledger["combined_Hilbert_inverse_exponent"]}


def load_module(root, name):
    spec = importlib.util.spec_from_file_location("nonlinear_ball_" + name,
                                                  root / "checks" / (name + ".py"))
    need(spec is not None and spec.loader is not None, "Evaluator loader unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate(root):
    data = load_json((root / "DATA.json").read_bytes())
    need(isinstance(data, dict) and set(data) == {"schema", "scope", *MODULES}
         and data["schema"] == "HDBLAST_NONLINEAR_CONSTRAINT_BALL_DATA_V1"
         and isinstance(data["scope"], str) and data["scope"], "Data schema mismatch")
    suites = {name: [] for name in ("sources", "source_models", *MODULES, "interfaces")}

    def checker(suite):
        def check(name, condition):
            need(isinstance(name, str) and name, "Check name must be a nonempty string")
            need(bool(condition), "Arithmetic failed: " + suite + ":" + name)
            need(name not in suites[suite], "Duplicate check: " + suite + ":" + name)
            suites[suite].append(name)
        return check

    count = public_sources(root, checker("sources"))
    exact = {"source_models": source_models(root, checker("source_models"))}
    for suite, module_name in MODULES.items():
        need(isinstance(data[suite], dict), "Claims must be objects: " + suite)
        module = load_module(root, module_name)
        try:
            if suite == "seed":
                exact[suite] = module.evaluate(data[suite], checker(suite),
                                              producer_path=root / "inputs/PINNED_SOURCE_JETS.py")
            else:
                exact[suite] = module.evaluate(data[suite], checker(suite))
        except ValueError as exc:
            if str(exc).startswith("Arithmetic failed: " + suite + ":"):
                raise
            raise ValueError("Arithmetic failed: " + suite + ":" + str(exc)) from exc
        need(isinstance(exact[suite], dict) and exact[suite], "Exact suite report required: " + suite)
    joined = checker("interfaces")
    joined("ball_uses_pinned_inverse", data["arithmetic"]["inverse_exponent"] ==
           exact["source_models"]["same_inverse_exponent"])
    joined("ball_uses_proved_seed_envelope", int(data["seed"]["B_seed_exponent"]) ==
           data["arithmetic"]["B_exponent"])
    joined("same_physical_multiplier_and_derivative_orders",
           data["weights"]["physical_multiplier_decimal_exponent"] ==
           data["arithmetic"]["physical_multiplier_exponent"] == 226
           and data["weights"]["potential_scalar_order"] == 16
           and data["weights"]["potential_oneform_order"] == 14
           and data["weights"]["source_max_order"] == 12)
    exact["interfaces"] = {"inverse_exponent": data["arithmetic"]["inverse_exponent"],
                            "seed_exponent": data["arithmetic"]["B_exponent"],
                            "physical_multiplier_exponent": 226}
    names = [suite + ":" + name for suite, entries in suites.items() for name in entries]
    return {"schema": "hdblast-nonlinear-constraint-ball-v1", "status": "PASS",
            "named_check_count": len(names), "named_checks": names,
            "suite_counts": {suite: len(entries) for suite, entries in suites.items()},
            "source_count": count, "exact": exact, "scope": data["scope"]}
