"""Authenticate selected immutable sources and replay finite observer interfaces.

Only the capsule's own pinned observer_oracles module is imported.
Ancestor theorem/data/manifest snapshots are read as bytes or JSON, never run.
Selected-member authentication is not full upstream closure reauthentication.
"""
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

BASE = '5be99501c7f58e4400212f374abb6334c22cbca5'
BASE_TREE = '0878f785882c138c96430018b985a5f4240df0f5'
PRODUCER_POLICY = 'No ancestor source producer or checker was imported or executed. Every selected repository file was read by git cat-file blob at the immutable base without text normalization.'
EXPECTED_SOURCE_ROWS_SHA256 = "6ca83f3d3b6d289da069649589b57da3c7446da042d25dae8dde846f029e1e19"
ORACLE_SHA256 = "ffa0d659862e61282c75978213bbac5c537a56b75659f4fe522572ff20218f80"
OBSERVER = {'all_angular_modes': True,
 'ambient_acceleration': 'sigma_n_over_6',
 'ambient_dimension': 5,
 'ambient_free_fall_claimed': False,
 'area_density_joint_order': 3,
 'augmented_joint_order': 3,
 'calibrated_signal': False,
 'clock_inverse_joint_order': 5,
 'common_duration': 'qualitative_short_positive',
 'connection_joint_order': 6,
 'curvature_joint_order': 5,
 'curvature_pullback_joint_order': 3,
 'cut_dimension': 3,
 'cut_lower_factor': '1/2',
 'cut_lower_reference': 'rho_b_squared_kappa',
 'cut_metric_joint_order': 3,
 'field_mixed_order': 10,
 'frame_joint_order': 4,
 'frame_mode': 'ambient_Fermi_Walker_intrinsic_parallel_screen',
 'full_S3_labels': True,
 'full_original_parameter_interval': True,
 'full_tensor_norms': True,
 'future_physical_extension_claimed': False,
 'historical_long_duration_claimed': False,
 'intrinsic_dimension': 4,
 'measured_confirmation': False,
 'metric_joint_order': 7,
 'numerical_observer_duration': False,
 'observer_mode': 'intrinsic_shell_normal_geodesic',
 'ordinary_two_sided_parameter_derivatives': True,
 'parameter_derivative_order': 2,
 'phase_joint_order': 4,
 'proper_clock': 'actual_past_proper_duration_r',
 'scalar_pullback_joint_order': 4,
 'shell_crossing_claimed': False,
 'status': 'actual_common_short_intrinsic_normal_observer',
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
    check("source_pin_schema", pins["schema"] == "HDBLAST_OBSERVER_CLOCK_SOURCE_PINS_V1")
    check("source_fixed_base_and_tree", pins["base_commit"] == BASE and pins["base_git_tree"] == BASE_TREE)
    check("source_producer_policy", pins["producer_policy"] == PRODUCER_POLICY)
    sources = pins["sources"]
    need(isinstance(sources,dict) and len(sources) == 25, "Expected 25 source aliases")
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
    check("qualified_member_join_count", joins == 18)

    provenance = load_json((root/"INPUT_PROVENANCE.json").read_bytes())
    need(isinstance(provenance,dict) and set(provenance) ==
         {"schema","base_commit","source_pins_sha256","snapshot_method",
          "snapshot_producers_executed","snapshot_count","files"},
         "Unexpected input provenance keys")
    check("provenance_schema_and_base", provenance["schema"] ==
          "HDBLAST_OBSERVER_CLOCK_INPUT_PROVENANCE_V1" and provenance["base_commit"] == BASE)
    check("provenance_source_index", provenance["source_pins_sha256"] == digest(raw_pins)
          and type(provenance["snapshot_count"]) is int and provenance["snapshot_count"] == 25)
    check("provenance_snapshot_only", provenance["snapshot_method"] ==
          "git cat-file of immutable base blobs" and provenance["snapshot_producers_executed"] is False)
    expected_files = {alias:{"local":"inputs/"+row["snapshot_filename"],
         "repo_path":row["repo_path"],"base_git_blob":row["base_git_blob"],
         "size":row["size"],"sha256":row["sha256"]} for alias,row in sources.items()}
    check("provenance_exact_selected_roster", provenance["files"] == expected_files)

    return sources,joins


def load_module(root):
    path = root/"checks/observer_oracles.py"
    need(path.is_file() and all(not item.is_symlink() for item in (path,*path.parents)),
         "Regular own oracle path required")
    need(digest(path.read_bytes()) == ORACLE_SHA256,
         "Own oracle SHA256 differs before import")
    spec = importlib.util.spec_from_file_location("observer_clock_own_oracle",path)
    need(spec is not None and spec.loader is not None, "Pinned oracle loader unavailable")
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def evaluate(root):
    data = load_json((root/"DATA.json").read_bytes())
    need(type(data) is dict and set(data) ==
         {"schema","base_commit","oracle","retained","observer"},
         "Data schema keys differ")
    need(data["schema"] == "HDBLAST_OBSERVER_CLOCK_DATA_V1" and data["base_commit"] == BASE,
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
    prior = load_json((root/"inputs"/sources["PR146_DATA"]["snapshot_filename"]).read_bytes())
    physical = load_json((root/"inputs"/sources["PR145_DATA"]["snapshot_filename"]).read_bytes())
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
    retained("physical_closure_same_typed_retained_values",
             physical["schema"] == "HDBLAST_SUBSIDIARY_CLOSURE_DATA_V1"
             and type(physical["retained"]) is dict and set(physical["retained"]) == set(accepted)
             and all(type(physical["retained"][key]) is type(value)
                     and physical["retained"][key] == value for key,value in accepted.items()))
    observer = data["observer"]; declared = checker("observer")
    need(type(observer) is dict and set(observer) == set(OBSERVER), "Observer contract roster differs")
    for name,expected in OBSERVER.items():
        declared(name, type(observer[name]) is type(expected) and observer[name] == expected)

    # Import follows selected source/index/provenance/retained/contract checks.
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
    need(report_names == [name.removeprefix("observer_oracles/") for name in suites["oracle"]],
         "Own oracle report/check names differ")
    check = checker("interfaces")
    model = exact["source_bound_model"]; o = data["oracle"]
    check("actual_ambient_intrinsic_cut_dimensions", observer["ambient_dimension"] == model["ambient_dimension"] ==
          accepted["spatial_dimension"]+1 == 5 and observer["intrinsic_dimension"] == model["intrinsic_dimension"] == 4
          and observer["cut_dimension"] == model["cut_dimension"] == 3)
    check("actual_full_norm_and_label_scope", observer["full_tensor_norms"] is model["full_product_norms"] is
          accepted["global_copy_norms"] is True and observer["all_angular_modes"] is
          model["all_angular_modes"] is observer["full_S3_labels"] is True)
    check("finite_regularities_charge_spatial_embedding", observer["field_mixed_order"] == prior["norm_order"] ==
          observer["metric_joint_order"]+3 and observer["connection_joint_order"]+1 == observer["metric_joint_order"]
          and observer["curvature_joint_order"]+2 == observer["metric_joint_order"])
    check("finite_phase_augmented_and_inverse_reserve", observer["phase_joint_order"]+2 == observer["connection_joint_order"]
          and observer["augmented_joint_order"]+1 == observer["phase_joint_order"]
          and observer["clock_inverse_joint_order"]+2 == observer["metric_joint_order"])
    check("finite_moving_frame_scalar_curvature_cut_reserve", observer["frame_joint_order"]+2 == observer["connection_joint_order"]
          and observer["scalar_pullback_joint_order"]+2 <= observer["metric_joint_order"]
          and observer["curvature_pullback_joint_order"]+2 == observer["curvature_joint_order"]
          and observer["cut_metric_joint_order"] == observer["area_density_joint_order"] == observer["augmented_joint_order"])
    check("proper_clock_and_intrinsic_acceleration_binding", o["past_velocity_sign"] == "-1"
          and o["q_velocity_mode"] == "rho_times_proper" and o["q_acceleration_mode"] == "rho_squared"
          and o["intrinsic_spray_mode"] == "intrinsic_LC" and o["Gauss_acceleration_factor"] == "1/6"
          and observer["ambient_free_fall_claimed"] is False)
    check("normal_cut_global_graph_dependency_binding", o["graph_tangent_mode"] == "Gauss_orthogonal"
          and o["global_graph_gate"] == "Gauss_timelikeK_positivecut_degreeone"
          and o["cut_positivity_gate"] == "uniform_C0_from_initial_gap"
          and observer["cut_lower_factor"] == "1/2" and observer["full_original_parameter_interval"] is True)
    check("actual_accelerated_frame_binding", o["FW_velocity_normal_factor"] ==
          o["FW_normal_velocity_factor"] == "1/6"
          and observer["frame_mode"] == "ambient_Fermi_Walker_intrinsic_parallel_screen")
    check("finite_witness_not_actual_coefficient_evaluation", model["actual_metric_or_curvature_arrays_evaluated"] is False)
    check("qualitative_short_scope_without_crossing_or_calibration", observer["common_duration"] ==
          "qualitative_short_positive" and observer["shell_crossing_claimed"] is
          observer["numerical_observer_duration"] is observer["historical_long_duration_claimed"] is
          observer["calibrated_signal"] is observer["measured_confirmation"] is False)
    named = [suite+":"+name for suite,entries in suites.items() for name in entries]
    return {"schema":"HDBLAST_OBSERVER_CLOCK_EVALUATION_V1","status":"PASS",
            "named_check_count":len(named),"named_checks":named,
            "suite_counts":{suite:len(entries) for suite,entries in suites.items()},
            "source_count":len(sources),"qualified_member_join_count":join_count,
            "own_factor_count":len(factors),"own_mutation_count":len(mutations),
            "exact":{"retained":data["retained"],"observer":observer,"oracle":exact},
            "scope":"Selected immutable source/member authentication, typed retained/observer claims and finite "
                    "clock, intrinsic geodesic, pullback, graph and frame interface witnesses. "
                    "Uniform common domain, cut positivity, global injectivity and actual C2 observer "
                    "dependence remain written continuum proofs. No numerical duration, long tube, "
                    "smooth shell crossing, actual coefficient evaluation or empirical calibration is proved by this replay."}
