"""Authenticate selected immutable sources and replay finite parameter interfaces.

Only the capsule's own pinned parameter_oracles module is imported.
Ancestor proof/data/manifest snapshots are read as bytes or JSON, never run.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

BASE = "fb2ba28453887056166e9f9b87a5a78b9bc22360"
BASE_TREE = "9e3ddc4b5d661f3b96c0ab0e37b1c2ccb10fd583"
PRODUCER_POLICY = ("No ancestor source producer was imported or executed. Blobs were read "
                   "with git cat-file blob and copied without text normalization.")
EXPECTED_SOURCE_ROWS_SHA256 = "773f989d4c8fa88353eaf0cfad2b31ff997e30fa5a3bbd1fd14e2983e5846de8"
ORACLE_SHA256 = "edc8a43f59fbf0a38533c565f1e489b8eee28017e3cf5653eeb4575dcceb4117"
MODULES = {"oracle": "parameter_oracles"}
RETAINED_KEYS = (
    "component_count", "spatial_dimension", "first_jet_directions", "copy_count",
    "ell_exponent", "T_factor", "R_exponent", "Cstar_exponent",
    "iterated_amplitude_power_exponent", "data_slope_exponent",
    "picard_C0_exponent", "picard_H_exponent", "radial_speed",
    "outer_support_end", "outer_wall", "outer_wall_margin",
    "source_C_zero", "global_copy_norms",
)



def need(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+"\n").encode("utf-8")


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def load_json(raw):
    def unique(items):
        result = {}
        for key,value in items:
            need(key not in result, "Duplicate JSON key: "+key)
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError("Nonfinite JSON constant: "+value)
    return json.loads(raw.decode("utf-8"), object_pairs_hook=unique,
                      parse_float=str, parse_constant=nonfinite)


def safe_relative(name):
    if not isinstance(name,str):
        return False
    path = PurePosixPath(name)
    return ("\\" not in name and bool(path.parts) and not path.is_absolute()
            and ".." not in path.parts and path.as_posix() == name)


def public_sources(root,check):
    raw_pins = (root/"SOURCE_PINS.json").read_bytes()
    pins = load_json(raw_pins)
    need(isinstance(pins,dict) and set(pins) ==
         {"schema","base_commit","base_git_tree","producer_policy","sources"},
         "Unexpected source pin keys")
    check("source_pin_schema", pins["schema"] == "hdblast-evolved-parameter-source-pins-v1")
    check("source_fixed_base_and_tree", pins["base_commit"] == BASE and pins["base_git_tree"] == BASE_TREE)
    check("source_producer_policy", pins["producer_policy"] == PRODUCER_POLICY)
    sources = pins["sources"]
    need(isinstance(sources,dict) and len(sources) == 50, "Expected 50 source aliases")
    check("exact_selected_source_index", digest(canonical(sources)) == EXPECTED_SOURCE_ROWS_SHA256)
    snapshots = set()
    for alias,row in sources.items():
        need(isinstance(alias,str) and re.fullmatch(r"[A-Z][A-Z0-9_]*",alias), "Unsafe source alias")
        need(isinstance(row,dict), "Source row must be an object")
        filename = row["snapshot_filename"]
        need(safe_relative(filename) and len(PurePosixPath(filename).parts) == 1
             and filename not in snapshots, "Unsafe or repeated snapshot name")
        need(safe_relative(row["repo_path"]) and row["repo_path"].startswith("research/"),
             "Unsafe source repository path")
        raw = (root/"inputs"/filename).read_bytes()
        check("source_size:"+alias, type(row["size"]) is int and row["size"] >= 0
              and len(raw) == row["size"])
        check("source_sha256:"+alias, isinstance(row["sha256"],str)
              and re.fullmatch(r"[0-9a-f]{64}",row["sha256"]) is not None
              and digest(raw) == row["sha256"])
        blob = b"blob "+str(len(raw)).encode("ascii")+b"\0"+raw
        check("source_git_blob:"+alias, isinstance(row["base_git_blob"],str)
              and re.fullmatch(r"[0-9a-f]{40}",row["base_git_blob"]) is not None
              and hashlib.sha1(blob).hexdigest() == row["base_git_blob"])
        snapshots.add(filename)
    need((root/"inputs").is_dir() and all(p.is_file() and not p.is_symlink()
         for p in (root/"inputs").iterdir()), "Flat regular source input directory required")
    check("closed_selected_input_roster", {p.name for p in (root/"inputs").iterdir()} == snapshots)

    manifest_members = {}
    for alias,row in sources.items():
        if not row["repo_path"].endswith("/MANIFEST.json"):
            continue
        manifest = load_json((root/"inputs"/row["snapshot_filename"]).read_bytes())
        need(isinstance(manifest,dict) and isinstance(manifest.get("files"),(dict,list)),
             "Ancestor manifest roster required: "+alias)
        entries = manifest["files"]
        if isinstance(entries,dict):
            entries = [dict(record,path=name) for name,record in entries.items()]
        by_name = {}
        for entry in entries:
            need(isinstance(entry,dict) and safe_relative(entry.get("path")),
                 "Unsafe ancestor member: "+alias)
            name = entry["path"]
            need(name not in by_name, "Repeated ancestor member: "+alias+":"+name)
            need(type(entry.get("bytes")) is int and entry["bytes"] >= 0
                 and isinstance(entry.get("sha256"),str)
                 and re.fullmatch(r"[0-9a-f]{64}",entry["sha256"]) is not None,
                 "Invalid ancestor member metadata")
            by_name[name] = {"bytes":entry["bytes"],"sha256":entry["sha256"]}
        manifest_members[alias] = by_name
    joins = 0
    for alias,row in sources.items():
        ancestors = set(); direct = set()
        for join in row["qualified_member_joins"]:
            ma = join["manifest_alias"]; member = join["member_filename"]
            need(ma in manifest_members and safe_relative(member), "Unknown manifest/member alias")
            declared_repo = str(PurePosixPath(sources[ma]["repo_path"]).parent/PurePosixPath(member))
            origin = declared_repo == row["repo_path"]
            expected_qualification = "direct_origin_member" if origin else "content_identical_declared_snapshot"
            check("qualified_member:"+alias+":"+ma+":"+member,
                  manifest_members[ma].get(member) == {"bytes":row["size"],"sha256":row["sha256"]}
                  and join["declared_size"] == row["size"]
                  and join["declared_sha256"] == row["sha256"]
                  and join["declared_member_repo_path"] == declared_repo
                  and join["qualification"] == expected_qualification)
            ancestors.add(ma)
            if origin:
                direct.add(ma)
            joins += 1
        need(row["ancestor_manifest_aliases"] == sorted(ancestors)
             and row["direct_container_manifest_aliases"] == sorted(direct),
             "Ancestor qualification roster differs: "+alias)
        need(bool(direct) or row["repo_path"].endswith("/MANIFEST.json"),
             "Nonmanifest source lacks direct originating member join: "+alias)
    check("qualified_member_join_count", joins == 141)

    provenance = load_json((root/"INPUT_PROVENANCE.json").read_bytes())
    need(isinstance(provenance,dict) and set(provenance) ==
         {"schema","base_commit","source_pins_sha256","snapshot_method",
          "snapshot_producers_executed","snapshot_count","files"},
         "Unexpected input provenance keys")
    check("provenance_schema_and_base", provenance["schema"] ==
          "HDBLAST_EVOLVED_PARAMETER_INPUT_PROVENANCE_V1" and provenance["base_commit"] == BASE)
    check("provenance_source_index", provenance["source_pins_sha256"] == digest(raw_pins)
          and type(provenance["snapshot_count"]) is int and provenance["snapshot_count"] == 50)
    check("provenance_snapshot_only", provenance["snapshot_method"] ==
          "git cat-file of immutable base blobs" and provenance["snapshot_producers_executed"] is False)
    expected_files = {alias:{"local":"inputs/"+row["snapshot_filename"],
         "repo_path":row["repo_path"],"base_git_blob":row["base_git_blob"],
         "size":row["size"],"sha256":row["sha256"]} for alias,row in sources.items()}
    check("provenance_exact_selected_roster", provenance["files"] == expected_files)

    return sources,joins


def source_models(root,sources,check):
    def data(alias):
        return load_json((root/"inputs"/sources[alias]["snapshot_filename"]).read_bytes())
    prior = data("PR144_DATA")
    check("prior_capsule_data_schema", prior["schema"] == "HDBLAST_NONLINEAR_PICARD_LIMIT_DATA_V1")
    accepted = prior["arithmetic"]
    check("actual_prior_dimensions", accepted["component_count"] == 16
          and accepted["spatial_dimension"] == 4 and accepted["first_jet_directions"] == 5)
    check("actual_prior_copy_and_zero_source_conventions", accepted["copy_count"] == 2
          and accepted["global_copy_norms"] is True and accepted["source_C_zero"] is True)
    check("actual_prior_admitted_interval_and_ball", accepted["R_exponent"] == -12
          and accepted["iterated_amplitude_power_exponent"] == 500209
          and accepted["data_slope_exponent"] == 4200
          and accepted["picard_C0_exponent"] == 7003210
          and accepted["picard_H_exponent"] == 500212)
    check("actual_prior_slab", accepted["ell_exponent"] == -6 and F(accepted["T_factor"]) == F(13,8))
    check("actual_prior_causal_geometry", F(accepted["radial_speed"]) == F(20,19)
          and F(accepted["outer_support_end"]) == F(80,19)
          and F(accepted["outer_wall"]) == 6 and F(accepted["outer_wall_margin"]) == F(34,19))
    cauchy = data("CAUCHY_DATA"); constraint = data("CONSTRAINT_DATA")
    check("actual_cauchy_data_schema", cauchy["schema"] == "HDBLAST_UNIFORM_CAUCHY_DATA_V1")
    check("actual_constraint_data_schema", constraint["schema"] == "HDBLAST_NONLINEAR_CONSTRAINT_BALL_DATA_V1")
    ca = cauchy["algebra"]; qa = constraint["algebra"]
    check("actual_half_action_density_factors", F(ca["K_conversion_factor"]) == 2
          and F(ca["K_trace_factor"]) == -F(1,3)
          and F(ca["normal_scalar_conversion_factor"]) == 1)
    check("actual_full_canonical_constraint_factors", F(qa["Hamiltonian_curvature_factor"]) == F(1,2)
          and F(qa["Hamiltonian_momentum_factor"]) == 2
          and F(qa["Hamiltonian_trace_squared_factor"]) == F(1,3)
          and F(qa["Hamiltonian_scalar_density_factor"]) == -F(1,2)
          and F(qa["momentum_h_dPi_factor"]) == -2
          and F(qa["momentum_matter_factor"]) == 1)
    flow = data("BOUNDARY_DATA"); coefficients = data("COEFFICIENT_DATA")
    check("actual_boundary_data_schema", flow["schema"] == "HDBLAST_BOUNDARY_FLOW_PATCH_DATA_V1")
    check("actual_flow_state_patch_and_jet_reserve", flow["flow"]["state_coordinates"] == 16
          and flow["flow"]["jet_total_order"] == 18
          and F(flow["flow"]["inner_initial_radius"]) == F(1,1000))
    check("actual_full_metric_image_patch", F(flow["metric"]["epsilon"]) == F(1,40)
          and F(flow["metric"]["rho_min"]) == 129 and F(flow["metric"]["rho_max"]) == 131)
    check("actual_physical_reconstruction_orders", coefficients["coefficients"]["coefficient_jet_order"] == 16
          and coefficients["coefficients"]["background_order"] == 19
          and coefficients["coefficients"]["point_map_mixed_factor"] == 2)
    energy = data("VARIABLE_ENERGY_DATA"); ea = energy["arithmetic"]
    check("actual_normalized_metric_parent", energy["schema"] == "HDBLAST_VARIABLE_COEFFICIENT_ENERGY_DATA_V1"
          and F(ea["A_ceiling"]) == F(1,19) and F(ea["B_lower"]) == F(931,1280)
          and F(ea["B_upper"]) == F(20,19) and F(ea["Bzz_lower"]) == F(19,20))
    check("actual_full_array_coefficient_certificate", accepted["Cstar_exponent"] ==
          ea["Cstar_exponent"] == 500000 and ea["coefficient_total_order"] == 16)
    text = " ".join((root/"inputs"/sources["UNIFORM_CAUCHY"]["snapshot_filename"]).read_text().split())
    check("current_half_action_text_binding", "R/2-|dphi|^2/2-U" in text
          and "lambda=J^-1/2" in text and "alpha=rho n, beta^i=0." in text)
    physical = data("PR145_DATA")
    check("completed_physical_closure_data_binding", physical["schema"] ==
          "HDBLAST_SUBSIDIARY_CLOSURE_DATA_V1" and physical["retained"] ==
          {key:accepted[key] for key in RETAINED_KEYS})
    check("completed_physical_closure_norm_and_admission_binding",
          type(physical["norm_order"]) is int and physical["norm_order"] == 14
          and type(physical["subsidiary_norm_order"]) is int
          and physical["subsidiary_norm_order"] == 13
          and physical["additional_smallness"] is False
          and physical["new_evaluated_subsidiary_constant"] is False)
    return {"retained":{key:accepted[key] for key in RETAINED_KEYS},
            "canonical_algebra":ca,"constraint_algebra":qa}


def load_module(root,name):
    need(name == "parameter_oracles", "Only the own pinned oracle may be imported")
    path = root/"checks"/(name+".py")
    need(path.is_file() and all(not item.is_symlink() for item in (path,*path.parents)),
         "Regular own oracle path required")
    need(digest(path.read_bytes()) == ORACLE_SHA256,
         "Own oracle SHA256 differs before import")
    spec = importlib.util.spec_from_file_location("evolved_parameter_"+name,path)
    need(spec is not None and spec.loader is not None, "Pinned oracle loader unavailable")
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


REGULARITY = {
    "status": "actual_C2_lower_mixed",
    "first_bounded_norm_order": 13, "second_bounded_norm_order": 12,
    "high_orders_weak_only": True, "full_Frobenius": True,
    "full_spatial_measure": True, "all_time_jets": True,
    "angular_cutoff": False, "fixed_shell_D_rank": 4,
    "fixed_shell_N_rank": 12, "fixed_outer_D_rank": 16,
    "top_data_continuity_claimed": False, "X15_claimed": False,
}


def evaluate(root):
    data = load_json((root/"DATA.json").read_bytes())
    need(isinstance(data,dict) and set(data) ==
         {"schema","base_commit","oracle","retained","norm_order",
          "first_variation_norm_order","second_variation_norm_order","regularity",
          "additional_smallness","new_evaluated_parameter_constant"},
         "Data schema keys differ")
    need(data["schema"] == "HDBLAST_EVOLVED_PARAMETER_DATA_V1" and data["base_commit"] == BASE,
         "Data schema or immutable base differs")
    suites = {name:[] for name in ("sources","source_models","retained","regularity","oracle","interfaces")}
    def checker(suite):
        def check(name,condition):
            need(isinstance(name,str) and bool(name), "Nonempty check name required")
            need(bool(condition), "Arithmetic failed: "+suite+":"+name)
            need(name not in suites[suite], "Duplicate check: "+suite+":"+name)
            suites[suite].append(name)
        return check
    sources,join_count = public_sources(root,checker("sources"))
    source = source_models(root,sources,checker("source_models"))
    retained = data["retained"]; joined = checker("retained")
    need(isinstance(retained,dict) and set(retained) == set(RETAINED_KEYS), "Retained claim roster differs")
    for name,expected in source["retained"].items():
        joined("same_accepted_"+name, type(retained[name]) is type(expected) and retained[name] == expected)
    joined("actual_norm_orders", type(data["norm_order"]) is int and data["norm_order"] == 10
           and type(data["first_variation_norm_order"]) is int and data["first_variation_norm_order"] == 12
           and type(data["second_variation_norm_order"]) is int and data["second_variation_norm_order"] == 11)
    joined("unchanged_admission", data["additional_smallness"] is False
           and data["new_evaluated_parameter_constant"] is False)
    regularity = data["regularity"]; typed = checker("regularity")
    need(isinstance(regularity,dict) and set(regularity) == set(REGULARITY),
         "Regularity contract roster differs")
    for name,expected in REGULARITY.items():
        typed(name, type(regularity[name]) is type(expected) and regularity[name] == expected)

    # This import occurs only after the selected source roster, all immutable
    # bytes/joins/provenance, retained values and typed status were checked.
    exact = {"source_models":{"retained":source["retained"]},"regularity":regularity}
    for suite,name in MODULES.items():
        need(isinstance(data[suite],dict), "Oracle claims must be an object")
        loaded = load_module(root,name)
        need(len(loaded.FACTORS) == len(loaded.MUTATIONS) == 81,
             "Exactly 81 own claim factors and mathematical controls required")
        exact[suite] = loaded.evaluate(data[suite],checker(suite))
        need(isinstance(exact[suite],dict) and exact[suite].get("status") == "pass"
             and len(suites[suite]) == 194 and len(exact[suite].get("checks",[])) == 194,
             "Exactly 194 passing own-oracle checks required")
    o = data["oracle"]; finite = exact["oracle"]["finite_exact"]; check = checker("interfaces")
    check("actual_sixteen_field_five_jet_join", o["source_fields"] == retained["component_count"] == 16
          and o["spatial_directions"] == retained["spatial_dimension"] == 4
          and o["firstjet_directions"] == retained["first_jet_directions"] == 5)
    check("actual_copy_and_no_angular_cutoff_join", o["norm_bundle_copies"] == retained["copy_count"] == 2
          and retained["global_copy_norms"] is True and regularity["angular_cutoff"] is False)
    check("actual_coefficient_and_initial_grade_join", o["coefficient_derivative_reserve"] == 16
          and o["initial_grade"] == 14 and o["temporal_last_recurrence"] == 12)
    check("actual_full_mixed_source_losses", o["energy_source_loss"] == 1
          and o["normal_source_loss"] == 2
          and o["mixed_norm_mode"] == "all_time_all_rotations_normal")
    check("actual_fixed_primal_boundary_join", o["initial_boundary_mode"] == "fixed_primal_chart"
          and regularity["fixed_shell_D_rank"]+regularity["fixed_shell_N_rank"] ==
          regularity["fixed_outer_D_rank"] == retained["component_count"])
    check("actual_initial_only_parameter_and_canonical_binding", o["parameter_binding"] ==
          "initial_only_fixed_background" and o["canonical_binding"] == "half_action_fixed_density")
    check("actual_variation_low_high_ladder", regularity["first_bounded_norm_order"] ==
          data["first_variation_norm_order"]+1 == 13
          and regularity["second_bounded_norm_order"] == data["second_variation_norm_order"]+1 == 12
          and data["norm_order"] == 10)
    check("source_bound_and_auxiliary_scope", finite["auxiliary_fields"] == 2
          and finite["auxiliary_firstjet_matrix_shape"] == [5,2]
          and exact["oracle"]["source_bound_model"]["actual_coefficient_arrays_evaluated"] is False)
    check("no_new_smallness_or_numeric_parameter_constant", o["new_numeric_estimate_constant"] is
          o["new_small_time_assumption"] is data["additional_smallness"] is
          data["new_evaluated_parameter_constant"] is False)
    exact["interfaces"] = {"evolution_components":16,"source_firstjet_directions":5,
          "C2_target_order":10,"first_continuous_order":12,"first_bounded_order":13,
          "second_continuous_order":11,"second_bounded_order":12,
          "additional_smallness":False,"new_evaluated_parameter_constant":False,
          "qualified_source_member_joins":join_count}
    named = [suite+":"+name for suite,entries in suites.items() for name in entries]
    return {"schema":"hdblast-evolved-parameter-evaluation-v1","status":"PASS",
          "named_check_count":len(named),"named_checks":named,
          "suite_counts":{suite:len(entries) for suite,entries in suites.items()},
          "source_count":len(sources),"qualified_member_join_count":join_count,"exact":exact,
          "scope":"Closed selected-source authentication, qualified manifest member joins and finite exact "
                  "parameter coefficient/residual/source interfaces. Full mixed estimates, compatible jets, "
                  "compactness, actual C2 dependence and physical composition are written analytic proofs; "
                  "this replay does not prove those continuum arguments or top-data/observer/empirical claims."}
