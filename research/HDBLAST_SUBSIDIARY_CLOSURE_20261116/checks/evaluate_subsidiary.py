"""Authenticate selected immutable sources and replay finite subsidiary interfaces.

Only the capsule's own pinned subsidiary_oracles module is imported.
Ancestor proof/data/manifest snapshots are read as bytes or JSON, never run.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

BASE = "6ee27c36dac67497a2454f7c1283fd7177a91888"
BASE_TREE = "6319e9c8dfc7794a41ad2d9b400885582115ed41"
PRODUCER_POLICY = ("No ancestor source producer was imported or executed. Blobs were read "
                   "with git cat-file blob and copied without text normalization.")
EXPECTED_SOURCE_ROWS_SHA256 = "e8c37cce825f1fc66149160911bf7a7f3f0ada8d1b7e1658ecdcbc6cd295ebca"
MODULES = {"oracle": "subsidiary_oracles"}
RETAINED_KEYS = (
    "component_count", "spatial_dimension", "first_jet_directions", "copy_count",
    "ell_exponent", "T_factor", "R_exponent", "Cstar_exponent",
    "iterated_amplitude_power_exponent", "data_slope_exponent",
    "picard_C0_exponent", "picard_H_exponent", "radial_speed",
    "outer_support_end", "outer_wall", "outer_wall_margin",
    "source_C_zero", "global_copy_norms",
)
PROOF_ALIASES = {
    "P1": "PR144_REDUCED_SUPPORT", "P2": "PR144_TOP_CONTINUITY",
    "P3": "COEFFICIENTS", "P4": "BOUNDARY_FLOW_PROOF",
    "P5": "BOUNDARY_COORDINATES", "P6": "MIXED_INITIAL_PRODUCTS",
    "P7": "UNIFORM_BOUNDARY_TRANSPORT", "P8": "UNIFORM_CAUCHY",
    "P9": "INITIAL_GAUGE_SUBSIDIARY",
}


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
    check("source_pin_schema", pins["schema"] == "hdblast-subsidiary-source-pins-v1")
    check("source_fixed_base_and_tree", pins["base_commit"] == BASE and pins["base_git_tree"] == BASE_TREE)
    check("source_producer_policy", pins["producer_policy"] == PRODUCER_POLICY)
    sources = pins["sources"]
    need(isinstance(sources,dict) and len(sources) == 41, "Expected 41 source aliases")
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
    check("qualified_member_join_count", joins == 92)

    provenance = load_json((root/"INPUT_PROVENANCE.json").read_bytes())
    need(isinstance(provenance,dict) and set(provenance) ==
         {"schema","base_commit","source_pins_sha256","snapshot_method",
          "snapshot_producers_executed","snapshot_count","files"},
         "Unexpected input provenance keys")
    check("provenance_schema_and_base", provenance["schema"] ==
          "HDBLAST_SUBSIDIARY_CLOSURE_INPUT_PROVENANCE_V1" and provenance["base_commit"] == BASE)
    check("provenance_source_index", provenance["source_pins_sha256"] == digest(raw_pins)
          and type(provenance["snapshot_count"]) is int and provenance["snapshot_count"] == 41)
    check("provenance_snapshot_only", provenance["snapshot_method"] ==
          "git cat-file of immutable base blobs" and provenance["snapshot_producers_executed"] is False)
    expected_files = {alias:{"local":"inputs/"+row["snapshot_filename"],
         "repo_path":row["repo_path"],"base_git_blob":row["base_git_blob"],
         "size":row["size"],"sha256":row["sha256"]} for alias,row in sources.items()}
    check("provenance_exact_selected_roster", provenance["files"] == expected_files)

    proof_pins = load_json((root/"PROOF_SOURCE_PINS.json").read_bytes())
    need(isinstance(proof_pins,dict) and set(proof_pins) == {"base_commit","sources"},
         "Unexpected physical proof source receipt keys")
    check("physical_proof_source_base", proof_pins["base_commit"] == BASE)
    rows = proof_pins["sources"]
    need(isinstance(rows,list) and [r.get("id") for r in rows] == list(PROOF_ALIASES),
         "Exact P1--P9 proof source roster required")
    for row in rows:
        need(isinstance(row,dict) and set(row) == {"id","path","bytes","sha256"},
             "Unexpected P1--P9 source row")
        alias = PROOF_ALIASES[row["id"]]; source = sources[alias]
        check("physical_proof_source:"+row["id"], row ==
              {"id":row["id"],"path":source["repo_path"],
               "bytes":source["size"],"sha256":source["sha256"]})
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
    return {"retained":{key:accepted[key] for key in RETAINED_KEYS},
            "canonical_algebra":ca,"constraint_algebra":qa}


def load_module(root,name):
    spec = importlib.util.spec_from_file_location("subsidiary_closure_"+name,
                                                 root/"checks"/(name+".py"))
    need(spec is not None and spec.loader is not None, "Pinned oracle loader unavailable")
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def evaluate(root):
    data = load_json((root/"DATA.json").read_bytes())
    need(isinstance(data,dict) and set(data) ==
         {"schema","base_commit","oracle","retained","norm_order","subsidiary_norm_order",
          "additional_smallness","new_evaluated_subsidiary_constant"}, "Data schema keys differ")
    need(data["schema"] == "HDBLAST_SUBSIDIARY_CLOSURE_DATA_V1" and data["base_commit"] == BASE,
         "Data schema or immutable base differs")
    suites = {name:[] for name in ("sources","source_models","retained","oracle","interfaces")}
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
    joined("actual_norm_orders", type(data["norm_order"]) is int and data["norm_order"] == 14
           and type(data["subsidiary_norm_order"]) is int and data["subsidiary_norm_order"] == 13)
    joined("unchanged_admission", data["additional_smallness"] is False
           and data["new_evaluated_subsidiary_constant"] is False)

    exact = {"source_models":{"retained":source["retained"]}}
    for suite,name in MODULES.items():
        need(isinstance(data[suite],dict), "Oracle claims must be an object")
        exact[suite] = load_module(root,name).evaluate(data[suite],checker(suite))
        need(isinstance(exact[suite],dict) and exact[suite].get("status") == "pass",
             "Exact oracle report required")
    o = data["oracle"]; finite = exact["oracle"]["finite_exact"]; check = checker("interfaces")
    ca = source["canonical_algebra"]; qa = source["constraint_algebra"]
    check("actual_half_action_metric_density_join", F(o["canonical_K_factor"]) == F(ca["K_conversion_factor"])
          and F(o["canonical_K_trace_factor"]) == -F(ca["K_trace_factor"]))
    check("actual_half_action_scalar_density_join", F(o["canonical_scalar_velocity_factor"]) ==
          F(ca["normal_scalar_conversion_factor"]) == 1 and F(o["density_determinant_exponent"]) == -F(1,2))
    check("actual_canonical_Hamiltonian_join", F(o["canonical_H_curvature_factor"]) == F(qa["Hamiltonian_curvature_factor"])
          and F(o["canonical_H_metric_factor"]) == F(qa["Hamiltonian_momentum_factor"])
          and F(o["canonical_H_scalar_factor"]) == F(qa["Hamiltonian_scalar_density_factor"]))
    check("actual_canonical_momentum_join", F(o["momentum_divergence_factor"]) == F(qa["momentum_h_dPi_factor"])
          and F(o["momentum_matter_factor"]) == F(qa["momentum_matter_factor"]))
    check("full_five_component_subsidiary_bundle", finite["physical_dimension"] ==
          finite["subsidiary_covector_components"] == retained["first_jet_directions"] == 5
          and finite["spatial_dimension"] == retained["spatial_dimension"] == 4)
    check("subsidiary_split_is_one_plus_four", finite["shell_D_rank"] == 1
          and finite["shell_N_rank"] == 4 and finite["shell_D_rank"]+finite["shell_N_rank"] == 5)
    check("subsidiary_ranks_not_evolution_ranks", retained["component_count"] == 16
          and finite["subsidiary_covector_components"] == 5)
    check("finite_initial_projection_factors", F(o["symmetric_gradient_factor"]) ==
          F(o["Einstein_trace_factor"]) == F(o["initial_normal_projection_factor"]) ==
          F(o["initial_mixed_projection_factor"]) == F(1,2))
    check("actual_covector_Robin_and_coordinate_orientation", F(o["coordinate_Robin_K_factor"]) == 2
          and F(o["y_to_z_sign"]) == -1 and F(o["momentum_projection_sign"]) == -1)
    check("regularity_ladder_is_interface_only", data["norm_order"] == data["subsidiary_norm_order"]+1 == 14
          and data["norm_order"]-2 == 12 and data["norm_order"]-3 == 11)
    check("auxiliary_reverse_chain_scope", finite["auxiliary_reverse_chain_fields"] == 3
          and finite["auxiliary_reverse_chain_directions"] == 4
          and retained["component_count"] == 16)
    check("auxiliary_ADM_shift_not_selected_binding", finite["nonzero_shift_initial_chart"] is True)
    exact["interfaces"] = {"evolution_components":16,"subsidiary_covector_components":5,
          "subsidiary_shell_D_rank":1,"subsidiary_shell_N_rank":4,
          "physical_field_order":14,"subsidiary_order":13,"initial_velocity_order":12,
          "subsidiary_wave_order":11,"additional_smallness":False,
          "new_evaluated_subsidiary_constant":False,"qualified_source_member_joins":join_count}
    named = [suite+":"+name for suite,entries in suites.items() for name in entries]
    return {"schema":"hdblast-subsidiary-closure-evaluation-v1","status":"PASS",
          "named_check_count":len(named),"named_checks":named,
          "suite_counts":{suite:len(entries) for suite,entries in suites.items()},
          "source_count":len(sources),"qualified_member_join_count":join_count,"exact":exact,
          "scope":"Closed selected-source authentication, qualified manifest member joins and finite exact "
                  "subsidiary algebra/source interfaces. The reconstructed unreduced equations, initial-data "
                  "join and coupled subsidiary energy uniqueness are written conditional analytic proofs; "
                  "this replay does not prove those continuum arguments or top-data/parameter/observer/empirical claims."}
