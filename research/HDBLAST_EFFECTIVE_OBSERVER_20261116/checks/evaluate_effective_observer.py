"""Authenticate selected immutable sources and replay finite effective ledgers.

Only the capsule's own pinned effective_oracles module is imported.
Ancestor theorem/data/manifest snapshots are read as bytes or JSON, never run.
Selected-member authentication is not full upstream closure reauthentication.
"""
from fractions import Fraction as F
from math import factorial
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

BASE = 'b4f2f59ab3918264c9832d87be1e39ed4dbfc178'
BASE_TREE = 'd472df7cc5ad0c38d07ef2b7e6dd51beb8faa189'
PRODUCER_POLICY = 'No ancestor producer, checker or oracle was imported or executed. Every selected repository source was read with git cat-file blob at the immutable base without text normalization.'
EXPECTED_SOURCE_ROWS_SHA256 = 'a90f8e106490f0aa664db31038d90dc03a6492f14b2635f65abcf1b1697376ee'
ORACLE_SHA256 = "eb5ea95da3968bd21521c960617d3c651083328523762d1159e7a5b7fef2d8d2"
OBSERVER = {'all_angular_modes': True,
 'ambient_acceleration': 'sigma_n_over_6',
 'ambient_dimension': 5,
 'ambient_free_fall_claimed': False,
 'area_density_joint_order': 3,
 'area_error_power2': 18,
 'area_numerator_mean_joint_order': 3,
 'arrival_theorem_claimed': False,
 'augmented_joint_order': 3,
 'calibrated_signal': False,
 'causal_shadow_over_ell': '19/40',
 'connection_joint_order': 6,
 'curvature_joint_order': 5,
 'curvature_pullback_joint_order': 3,
 'cut_dimension': 3,
 'cut_lower_factor': '1/2',
 'cut_lower_reference': 'rho_b_squared_kappa',
 'cut_metric_joint_order': 3,
 'delayed_bound_field_symbol': 'Mfield',
 'delayed_coordinate_drift_factor': 24,
 'field_mixed_order': 10,
 'frame_joint_order': 4,
 'frame_mode': 'ambient_Fermi_Walker_intrinsic_parallel_screen',
 'full_S3_labels': True,
 'full_original_parameter_interval': True,
 'full_tensor_norms': True,
 'future_physical_extension_claimed': False,
 'intrinsic_dimension': 4,
 'jacobi_coefficient_error_power2': 15,
 'jacobi_deviation_power2': 16,
 'measured_confirmation': False,
 'metric_delta_bound_power10': -20,
 'metric_delta_derivative_order': 4,
 'metric_delta_field_envelope_power10': 11,
 'metric_joint_order': 7,
 'nonzero_response_claimed': False,
 'numerical_ambient_tidal_bound_claimed': False,
 'numerical_parameter_derivative_bounds_claimed': False,
 'observer_mode': 'intrinsic_shell_normal_geodesic',
 'ordinary_two_sided_parameter_derivatives': True,
 'parameter_derivative_order': 2,
 'phase_joint_order': 4,
 'phase_position_error_factor': 24,
 'phase_velocity_error_factor': 48,
 'physical_seconds_calibrated': False,
 'proper_clock': 'physical_past_duration_r',
 'q_clock': 'r_over_rho_b',
 'q_star_over_ell': '6145/4096',
 'q_target_over_ell': '25/16',
 'scalar_bound_scope': 'values_and_positive_weight_mean_not_parameter_derivatives',
 'scalar_embedding_power10': 7,
 'scalar_factorial_order': 11,
 'scalar_pullback_joint_order': 4,
 'scalar_taylor_derivative_order': 11,
 'shell_crossing_claimed': False,
 'specified_q_star_admitted': True,
 'status': 'actual_effective_intrinsic_normal_observer',
 'third_parameter_derivative_claimed': False,
 'top_data_continuity_claimed': False}


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
    check("source_pin_schema", pins["schema"] == "HDBLAST_EFFECTIVE_OBSERVER_SOURCE_PINS_V1")
    check("source_fixed_base_and_tree", pins["base_commit"] == BASE and pins["base_git_tree"] == BASE_TREE)
    check("source_producer_policy", pins["producer_policy"] == PRODUCER_POLICY)
    sources = pins["sources"]
    need(isinstance(sources,dict) and len(sources) == 21, "Expected 21 source aliases")
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
            need(origin, "Only direct originating member joins are admitted")
            expected_qualification = "direct_origin_member"
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
        anchor = row["repo_path"].endswith("/MANIFEST.json")
        need((anchor and not direct and not row["qualified_member_joins"])
             or (not anchor and len(direct) == 1 and len(row["qualified_member_joins"]) == 1),
             "Each selected nonmanifest requires exactly one direct origin join: "+alias)
    check("qualified_member_join_count", joins == 15)

    provenance = load_json((root/"INPUT_PROVENANCE.json").read_bytes())
    need(isinstance(provenance,dict) and set(provenance) ==
         {"schema","base_commit","source_pins_sha256","snapshot_method",
          "snapshot_producers_executed","snapshot_count","files"},
         "Unexpected input provenance keys")
    check("provenance_schema_and_base", provenance["schema"] ==
          "HDBLAST_EFFECTIVE_OBSERVER_INPUT_PROVENANCE_V1" and provenance["base_commit"] == BASE)
    check("provenance_source_index", provenance["source_pins_sha256"] == digest(raw_pins)
          and type(provenance["snapshot_count"]) is int and provenance["snapshot_count"] == 21)
    check("provenance_snapshot_only", provenance["snapshot_method"] ==
          "git cat-file of immutable base blobs" and provenance["snapshot_producers_executed"] is False)
    expected_files = {alias:{"local":"inputs/"+row["snapshot_filename"],
         "repo_path":row["repo_path"],"base_git_blob":row["base_git_blob"],
         "size":row["size"],"sha256":row["sha256"]} for alias,row in sources.items()}
    check("provenance_exact_selected_roster", provenance["files"] == expected_files)

    return sources,joins


def load_module(root):
    path = root/"checks/effective_oracles.py"
    need(path.is_file() and all(not item.is_symlink() for item in (path,*path.parents)),
         "Regular own oracle path required")
    need(digest(path.read_bytes()) == ORACLE_SHA256,
         "Own oracle SHA256 differs before import")
    spec = importlib.util.spec_from_file_location("effective_observer_own_oracle",path)
    need(spec is not None and spec.loader is not None, "Pinned oracle loader unavailable")
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def evaluate(root):
    data = load_json((root/"DATA.json").read_bytes())
    need(type(data) is dict and set(data) ==
         {"schema","base_commit","oracle","retained","observer"},
         "Data schema keys differ")
    need(data["schema"] == "HDBLAST_EFFECTIVE_OBSERVER_DATA_V1" and data["base_commit"] == BASE,
         "Data schema or immutable base differs")
    suites = {name:[] for name in ("sources","retained","observer","oracle","interfaces")}
    def checker(suite):
        def check(name,condition):
            need(type(name) is str and bool(name), "Nonempty check name required")
            need(bool(condition), "Arithmetic failed: "+suite+":"+name)
            need(name not in suites[suite], "Duplicate check: "+suite+":"+name)
            suites[suite].append(name)
        return check
    sources,join_count = public_sources(root,checker("sources"))
    def source_data(alias):
        return load_json((root/"inputs"/sources[alias]["snapshot_filename"]).read_bytes())
    prior = source_data("PR146_DATA"); clock = source_data("PR147_DATA")
    retained = checker("retained")
    retained("accepted_parameter_data_schema", prior["schema"] == "HDBLAST_EVOLVED_PARAMETER_DATA_V1")
    accepted = prior["retained"]
    need(type(accepted) is dict and type(data["retained"]) is dict
         and set(data["retained"]) == set(accepted), "Retained claim roster differs")
    for name,expected in accepted.items():
        retained("same_accepted_"+name,
                 type(data["retained"][name]) is type(expected) and data["retained"][name] == expected)
    retained("accepted_physical_mixed_C2_order", type(prior["norm_order"]) is int and prior["norm_order"] == 10
             and type(prior["first_variation_norm_order"]) is int and prior["first_variation_norm_order"] == 12
             and type(prior["second_variation_norm_order"]) is int and prior["second_variation_norm_order"] == 11)
    regularity = prior["regularity"]
    retained("accepted_C2_full_mixed_premise", regularity["status"] == "actual_C2_lower_mixed"
             and regularity["all_time_jets"] is regularity["full_Frobenius"] is
             regularity["full_spatial_measure"] is True
             and regularity["angular_cutoff"] is regularity["top_data_continuity_claimed"] is
             regularity["X15_claimed"] is prior["additional_smallness"] is
             prior["new_evaluated_parameter_constant"] is False)
    retained("preceding_clock_same_typed_retained_values",
             clock["schema"] == "HDBLAST_OBSERVER_CLOCK_DATA_V1"
             and type(clock["retained"]) is dict and set(clock["retained"]) == set(accepted)
             and all(type(clock["retained"][key]) is type(value)
                     and clock["retained"][key] == value for key,value in accepted.items()))
    observer = data["observer"]; declared = checker("observer")
    need(type(observer) is dict and set(observer) == set(OBSERVER), "Observer contract roster differs")
    for name,expected in OBSERVER.items():
        declared(name, type(observer[name]) is type(expected) and observer[name] == expected)
    old_observer = clock["observer"]
    retained_targets = ("parameter_derivative_order","field_mixed_order","metric_joint_order",
                       "connection_joint_order","curvature_joint_order","phase_joint_order",
                       "augmented_joint_order","frame_joint_order","scalar_pullback_joint_order",
                       "curvature_pullback_joint_order","cut_metric_joint_order","area_density_joint_order")
    retained("preceding_observer_actual_C2_targets", all(type(old_observer[key]) is type(observer[key])
             and old_observer[key] == observer[key] for key in retained_targets))
    retained("preceding_intrinsic_normal_and_full_label_scope",
             old_observer["status"] == "actual_common_short_intrinsic_normal_observer"
             and old_observer["observer_mode"] == observer["observer_mode"]
             and old_observer["full_tensor_norms"] is old_observer["all_angular_modes"] is
             old_observer["full_S3_labels"] is old_observer["full_original_parameter_interval"] is True)
    # Import only after selected source/index/provenance/retained/typed checks.
    loaded = load_module(root)
    need(type(data["oracle"]) is dict, "Oracle claims must be an object")
    factors = tuple(loaded.FACTORS); mutations = tuple(loaded.MUTATIONS)
    need(factors and len(set(factors)) == len(factors), "Distinct own factor roster required")
    need(mutations and all(type(row) is tuple and len(row) == 2 for row in mutations),
         "Own mutation rows must be field/value pairs")
    mutation_names = [row[0] for row in mutations]
    need(len(set(mutation_names)) == len(mutation_names) and set(mutation_names) == set(factors),
         "Own mutation roster must cover exactly the own factor roster")
    exact = loaded.evaluate(data["oracle"],checker("oracle"))
    need(type(exact) is dict and exact.get("status") == "pass"
         and type(exact.get("checks")) is list and len(exact["checks"]) == len(suites["oracle"])
         and all(type(row) is dict and row.get("status") == "pass" for row in exact["checks"]),
         "Own oracle report/check roster differs")
    report_names = [row.get("name") for row in exact["checks"]]
    need(report_names == [name.removeprefix("effective_oracles/") for name in suites["oracle"]],
         "Own oracle report/check names differ")
    check = checker("interfaces"); model = exact["source_bound_model"]; o = data["oracle"]
    check("retained_full_state_dimensions", model["state_fields"] == accepted["component_count"] == 16
          and model["first_jet_directions"] == accepted["first_jet_directions"] == 5
          and model["copy_count"] == accepted["copy_count"] == 2
          and observer["intrinsic_dimension"] == model["intrinsic_dimension"] == accepted["spatial_dimension"] == 4
          and observer["ambient_dimension"] == observer["intrinsic_dimension"]+1 == 5
          and observer["cut_dimension"] == model["cut_dimension"] == 3)
    check("full_norm_zero_and_no_angular_cutoff_binding", model["full_Frobenius_arrays"] is
          model["all_angular_modes"] is model["zero_perturbation_included"] is
          observer["full_tensor_norms"] is observer["all_angular_modes"] is accepted["global_copy_norms"] is True
          and o["source_norm_mode"] == "full_arrays_distinct_operator_norms"
          and o["homogeneous_zero_mode"] == "weak_at_zero")
    check("retained_exponent_DAG_inputs", o["amplitude_power_index"] == accepted["iterated_amplitude_power_exponent"]
          and o["growth_power_index"] == accepted["picard_H_exponent"]
          and o["M_prefactor_power"] == accepted["picard_C0_exponent"]+accepted["data_slope_exponent"]+1
          and exact["exponent_DAG"]["huge_integer_or_huge_power_materialized"] is False)
    check("fixed_q_and_coordinate_slab", o["Q_factor"] == observer["q_target_over_ell"]
          and o["T_factor"] == accepted["T_factor"] and o["qstar_factor"] == observer["q_star_over_ell"]
          and F(o["T_factor"])-F(o["Q_factor"]) == F(1,16)
          and F(o["Q_factor"])-F(o["qstar_factor"]) == F(255,4096)
          and observer["specified_q_star_admitted"] is True)
    check("metric_embedding_and_smallness_inputs", o["embedding_power"] == observer["scalar_embedding_power10"] == 7
          and o["metric_gain_power"] == observer["metric_delta_field_envelope_power10"] == 11
          and F(o["delta_cap"]) == F(1,10**(-observer["metric_delta_bound_power10"]))
          and observer["metric_delta_derivative_order"] == 4)
    check("phase_drift_and_delayed_restart_inputs", F(o["velocity_error_factor"]) == observer["phase_velocity_error_factor"] == 48
          and F(o["drift_error_factor"]) == observer["phase_position_error_factor"] ==
          observer["delayed_coordinate_drift_factor"] == 24
          and F(1,2)/F(accepted["radial_speed"]) == F(observer["causal_shadow_over_ell"]) == F(19,40))
    check("complete_Jacobi_and_area_inputs", F(o["Jacobi_E_error_factor"]) == 2**observer["jacobi_coefficient_error_power2"]
          and F(o["Jacobi_S_error_factor"]) == 2**observer["jacobi_deviation_power2"]
          and F(o["relative_area_error_factor"]) == 2**observer["area_error_power2"]
          and F(o["cut_metric_floor"]) == F(observer["cut_lower_factor"]) == F(1,2)
          and o["Jacobi_sign_mode"] == "hyperbolic_plus")
    check("scalar_eleven_jet_Taylor_kernel_inputs", o["Taylor_order"] == observer["scalar_taylor_derivative_order"] ==
          observer["scalar_factorial_order"] == 11 and F(o["Taylor_kernel_denominator"]) == factorial(11)
          and observer["scalar_taylor_derivative_order"]+3 == model["individual_member_mixed_order"] == 14
          and model["parameter_field_mixed_order"] == observer["field_mixed_order"] == prior["norm_order"] == 10
          and model["C2_parameter_into_top_order_claimed"] is model["numerical_parameter_derivative_bounds"] is False
          and observer["scalar_bound_scope"] == "values_and_positive_weight_mean_not_parameter_derivatives")
    check("finite_parameter_and_composition_reserve", observer["metric_joint_order"]+3 == observer["field_mixed_order"] == 10
          and observer["connection_joint_order"]+1 == observer["metric_joint_order"]
          and observer["curvature_joint_order"]+2 == observer["metric_joint_order"]
          and observer["phase_joint_order"]+2 == observer["connection_joint_order"]
          and observer["augmented_joint_order"]+1 == observer["phase_joint_order"]
          and observer["frame_joint_order"]+2 == observer["connection_joint_order"]
          and observer["curvature_pullback_joint_order"]+2 == observer["curvature_joint_order"]
          and observer["area_numerator_mean_joint_order"] == observer["area_density_joint_order"] == observer["cut_metric_joint_order"] == 3)
    check("proper_clock_and_accelerated_intrinsic_scope", observer["q_clock"] == "r_over_rho_b"
          and observer["proper_clock"] == "physical_past_duration_r"
          and observer["ambient_acceleration"] == "sigma_n_over_6"
          and observer["frame_mode"] == old_observer["frame_mode"])
    check("certified_horizon_and_upper_response_only", o["duration_scope"] == "certified_Q_only"
          and o["causal_shadow_mode"] == "shadow_zero_and_upper_bound_only"
          and observer["arrival_theorem_claimed"] is observer["nonzero_response_claimed"] is
          observer["calibrated_signal"] is observer["measured_confirmation"] is
          observer["physical_seconds_calibrated"] is False)
    check("finite_witness_not_actual_coefficient_evaluation", model["actual_metric_or_curvature_arrays_evaluated"] is False)
    named = [suite+":"+name for suite,entries in suites.items() for name in entries]
    return {"schema":"HDBLAST_EFFECTIVE_OBSERVER_EVALUATION_V1","status":"PASS",
            "named_check_count":len(named),"named_checks":named,
            "suite_counts":{suite:len(entries) for suite,entries in suites.items()},
            "source_count":len(sources),"qualified_member_join_count":join_count,
            "own_factor_count":len(factors),"own_mutation_count":len(mutations),
            "exact":{"retained":data["retained"],"observer":observer,"oracle":exact},
            "scope":"Selected immutable source/member authentication, typed retained/effective claims and finite "
                    "rational inequality directions, symbolic exponent and causal-envelope arithmetic. "
                    "Continuum embedding, metric bounds, actual observer continuation, support, Jacobi cuts, "
                    "global injectivity and scalar/area envelopes remain written proofs. "
                    "No arrival, nonzero response, unverified longer horizon, numerical parameter Hessian, "
                    "ambient tidal evaluation or experimental calibration is established by this replay."}
