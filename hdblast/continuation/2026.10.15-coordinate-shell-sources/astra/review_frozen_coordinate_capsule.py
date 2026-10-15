#!/usr/bin/env python3
"""Read-only frozen-byte, roster, golden, and scientific-claim audit.

Does not import or execute any capsule code and does not rerun scientific
programs. Writes one review receipt outside the frozen capsule.
"""
import ast
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import stat
import zipfile

ROOT = Path(__file__).resolve().parents[1]
CAPSULE = "HDBLAST_COORDINATE_SHELL_SOURCES_20261015"
FROZEN = ROOT / "frozen-source" / CAPSULE
MANIFEST_SHA = "3ab5b0a7ac9045814b0aa59b6eae0d291a19adc76dd31913fd371203f1c132d5"
ZIP_SHA = "7c171631836cc91e546f74ff7897746175d31de42e38a39955c9f06e6a419ae5"
ZIP_BYTES = 3522306


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def read_json(path):
    return json.loads(path.read_bytes())


def get_at(value, keys):
    for key in keys:
        value = value[key]
    return value


def same(a,b):
    if type(a) is not type(b):
        return False
    if isinstance(a,dict):
        return a.keys() == b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a) == len(b) and all(same(x,y) for x,y in zip(a,b))
    return a == b


def main():
    checks = []
    def require(name, condition):
        if not condition:
            raise RuntimeError(name)
        checks.append(name)

    freeze_path = ROOT / "frozen-source/SOURCE_FREEZE_RECEIPT.json"
    config_path = ROOT / "publication/SOURCE_ASSEMBLY_CONFIG.json"
    freeze, approved = read_json(freeze_path),read_json(config_path)
    require("external_freeze_identity", freeze["status"] == "PASS" and freeze["capsule"] == CAPSULE and
            freeze["manifest"]["sha256"] == MANIFEST_SHA and freeze["archive"]["sha256"] == ZIP_SHA and
            freeze["archive"]["bytes"] == ZIP_BYTES and freeze["source_files"] == 55 and freeze["archive_members"] == 56)
    require("approved_root_source_inventory", approved["status"] == "APPROVED_FOR_FROZEN_SOURCE" and approved["capsule"] == CAPSULE)
    manifest_raw = (FROZEN / "MANIFEST.json").read_bytes()
    require("actual_manifest_external_pin", sha(manifest_raw) == MANIFEST_SHA and len(manifest_raw) == freeze["manifest"]["bytes"])
    manifest = json.loads(manifest_raw)
    registry = {x["path"]:x for x in manifest["files"]}
    selected = {x["path"]:x for x in approved["files"]}
    require("55_exact_approved_members", len(registry) == len(manifest["files"]) == len(selected) == len(approved["files"]) == 55 and registry.keys() == selected.keys())
    actual = set()
    data = {}
    for p in FROZEN.rglob("*"):
        require("regular_member:"+p.relative_to(FROZEN).as_posix(), not p.is_symlink() and (p.is_file() or p.is_dir()))
        if p.is_file():
            actual.add(p.relative_to(FROZEN).as_posix())
    require("frozen_tree_exact_no_extras", actual == registry.keys() | {"MANIFEST.json"})
    for rel,item in registry.items():
        path = PurePosixPath(rel)
        require("canonical_path:"+rel, not path.is_absolute() and ".." not in path.parts and "\\" not in rel and path.as_posix() == rel)
        raw = (FROZEN/rel).read_bytes()
        data[rel] = raw
        require("approved_frozen_source_bytes:"+rel, len(raw) == item["bytes"] == selected[rel]["bytes"] and sha(raw) == item["sha256"] == selected[rel]["sha256"])
        source = ROOT / selected[rel].get("source",rel)
        require("current_accepted_source_bytes:"+rel, source.read_bytes() == raw)
    archive_path = FROZEN.with_suffix(".zip")
    zip_raw = archive_path.read_bytes()
    require("actual_zip_external_pin", len(zip_raw) == ZIP_BYTES and sha(zip_raw) == ZIP_SHA)
    with zipfile.ZipFile(archive_path) as z:
        infos = z.infolist()
        names = [x.filename for x in infos]
        expected = {CAPSULE+"/"+name for name in actual}
        require("zip_exact_56_unique_members", len(names) == len(set(names)) == 56 and set(names) == expected)
        for info in infos:
            rel = info.filename[len(CAPSULE)+1:]
            mode = info.external_attr >> 16
            require("zip_regular_unencrypted:"+rel, not info.is_dir() and not stat.S_ISLNK(mode) and not info.flag_bits & 1)
            require("zip_identical_member_bytes:"+rel, z.read(info) == (FROZEN/rel).read_bytes())
        require("zip_crc_validation", z.testzip() is None)

    # Exact inventory and public input provenance, plus targeted sensitive-data
    # signatures. This is a bounded distribution audit, not a universal DLP proof.
    prohibited = {".git",".aws",".ssh",".env","_internal","__pycache__","node_modules"}
    require("no_private_or_raw_response_member_paths", all(not(prohibited & set(PurePosixPath(x).parts)) for x in actual))
    fingerprints = [rb"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----", rb"gh[pousr]_[A-Za-z0-9]{30,}", rb"github_pat_[A-Za-z0-9_]{30,}", rb"AKIA[0-9A-Z]{16}"]
    require("no_recognized_credential_material", all(not re.search(pattern,raw) for raw in data.values() for pattern in fingerprints))
    lit = json.loads(data["literature/LITERATURE_LEDGER.json"])
    require("raw_literature_responses_excluded", lit["publication"]["raw_primary_passages_distributed"] is False and
            lit["publication"]["raw_tool_responses_distributed"] is False and
            all(x["distributed"] is False and "literature/"+x["relative_path"] not in registry for x in lit["raw_internal_results"]))
    input_pins = json.loads(data["inputs/INPUT_PINS.json"])
    require("public_input_provenance", len(input_pins["files"]) == 18 and all(
        x["path"] in registry and x["sha256"] == registry[x["path"]]["sha256"] and
        x["original_public_source"]["repository"] == "maldonado-research/HDblast" and
        re.fullmatch("[0-9a-f]{40}",x["original_public_source"]["commit"]) and
        x["original_public_source"]["url"].startswith("https://github.com/maldonado-research/HDblast/blob/")
        for x in input_pins["files"]))

    # Parse literal runner configuration, never execute/import the runner.
    parsed = ast.parse(data["reproduce.py"])
    embedded = {}
    for node in parsed.body:
        if isinstance(node,ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0],ast.Name):
            name=node.targets[0].id
            if name in {"CAPSULE","PYTHON_VERSION","CONFIG_SHA256","CHECKERS","CONTROLS","DEPENDENCY_VERSIONS"}:
                embedded[name] = ast.literal_eval(node.value)
    checkers,controls = embedded["CHECKERS"],embedded["CONTROLS"]
    require("runner_approved_configuration_identity", embedded["CAPSULE"] == CAPSULE and
            embedded["CONFIG_SHA256"] == sha((ROOT/"publication/RUNNER_CONFIG.json").read_bytes()))
    roster = json.loads(data["REPLAY_ROSTER.json"])
    expected_counts = {"successful_checker_runs":8,"external_corruption_processes_rejected":66,
                       "embedded_wrong_formula_controls_rejected":0,"finite_positive_check_groups_each_mode":105}
    require("four_positive_33_control_roster", len(checkers) == 4 and len(controls) == 33 and roster["counts"] == expected_counts and
            sum(c["groups"] for c in checkers) == 105 and sum(c["embedded"] for c in checkers) == 0)
    require("runner_roster_exact_positive_identity", roster["positive_processes"] == [
        {"name":c["name"],"script":c["script"],"finite_positive_groups_each_mode":c["groups"],
         "embedded_controls_each_mode":c["embedded"],"finite_positive_check_names":c["check_names"],
         "embedded_control_names":c["embedded_control_names"]} for c in checkers])
    require("runner_roster_exact_control_identity", roster["external_controls"] == [
        {"checker":c["checker"],"name":c["name"],"expected_reason":c["reason"],"mechanism":c["mechanism"]} for c in controls])
    mode_comparisons=[]
    for c in checkers:
        require("all_process_inputs_frozen:"+c["name"], all(x in registry for x in [c["script"],*c["dependencies"],c["golden"],c["optimized_golden"]]))
        normal=json.loads(data[c["golden"]]); optimized=json.loads(data[c["optimized_golden"]])
        require("passing_current_goldens:"+c["name"], get_at(normal,c["status_path"]) == get_at(optimized,c["status_path"]) == c["status"] and c["byte_equal"] is True)
        if c["checks"] is not None:
            for output in [normal,optimized]:
                value=get_at(output,c["checks"]["path"])
                require("named_golden_groups:"+c["name"], value == c["check_names"] and len(value) == len(set(value)) == c["groups"])
        normalized = [deepcopy(normal),deepcopy(optimized)]
        for field in c["mode_fields"]:
            require("declared_actual_mode_field:"+c["name"], same(get_at(normal,field["path"]),field["normal"]) and same(get_at(optimized,field["path"]),field["optimized"]))
            for item in normalized:
                parent=get_at(item,field["path"][:-1])
                del parent[field["path"][-1]]
        require("complete_mode_scientific_json:"+c["name"], same(*normalized))
        hash_field = {"generator":"producer_sha256","astra":"checker_sha256","sol":"checker_sha256","referee":"source_script_sha256"}[c["name"]]
        require("current_scientific_script_golden:"+c["name"], normal[hash_field] == optimized[hash_field] == registry[c["script"]]["sha256"])
        mode_comparisons.append({"checker":c["name"],"complete_scientific_json_equal":True,"normalized_fields":[f["path"] for f in c["mode_fields"]]})
    source = data["reproduce.py"].decode()
    require("runner_full_json_and_byte_comparison_present", "require(json_equal(data, json.loads(golden_raw))" in source and
            "require(raw == golden_raw" in source and "Complete generated data bytes differ" in source)
    require("runner_checked_mode_fields_only", "type(value) is type(expected) and value == expected" in source and
            "normalized = [field['path'] for field in checker['mode_fields']]" in source and
            "require(json_equal(data[0], data[1])" in source)
    require("runner_strict_named_failure_and_postauth", "require(result.returncode == 1" in source and
            "validate_failure_diagnostic(control, diagnostic)" in source and
            "require(authenticate(root, args.manifest_sha256) == (authentication, inventory)" in source)

    # Inspect certified values already established by accepted science; do not
    # run the generators/checkers or reconstruct their quadratures.
    scientific=json.loads(data["computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json"])
    exact=scientific["exact"]
    require("producer_input_inventory_current", scientific["input_pins_sha256"] == registry["inputs/INPUT_PINS.json"]["sha256"])
    require("EN_scope_not_all_degrees", scientific["scope"]["full_U2_evaluated"] is False and scientific["scope"]["EN_scope"] == "Retained degree-four reduced Robin forcing only.")
    for key in ("EN","dEN","EN_material"):
        require("EN_H4_norm:"+key,exact["norms"][key]["angular_field"] == "H4")
    require("norm_prefactor_harmonic_ratio", all(F(x)*F(19,54) == F(y) for x,y in zip(exact["L2_prefactor"],exact["EN_L2_prefactor"])))
    require("Sp_unsquared_relative_norm_claim", F(exact["metric_correction_squared_norm_ratio_upper"]["Sp"]) < F(9,10**13)**2)
    astra=json.loads(data["astra/ASTRA_COORDINATE_NORMAL.json"])
    require("central_two_zero_sign_pattern", F(astra["central_zero"]["Sp_lower_at_any_scalar_zero"]["exact"]) > 0 and
            F(astra["central_zero"]["SK_upper_at_any_scalar_zero"]["exact"]) < 0 and all(
            F(v["Sp_upper"]["exact"]) < 0 < F(v["SK_lower"]["exact"]) for v in astra["central_endpoint_source_signs"].values()))
    tail=exact["tail"]["coefficients"]
    require("retained_published_tail_sign_thresholds", F(tail["Sp"][1]) < F("-2.061e-7") and
            F(tail["SK"][0]) > F("6.872e-8") and F(tail["EN"][0]) > F("2.748e-7"))
    for field,lower,upper in [("Sp","1.444e21","2.161e21"),("SK","1.604e20","2.401e20"),("EN","9.035e20","1.352e21")]:
        bounds=exact["norms"][field]["squared_L2"]
        require("published_squared_norm_thresholds:"+field, F(lower) < F(bounds[0]) <= F(bounds[1]) < F(upper))
    claims=json.loads(data["CLAIMS.json"])
    require("claims_keep_open_response_and_no_publication_extrapolation", "Full quantitative second-order response U2" in claims["not_established"] and
            "Gauge-invariant source or response signs" in claims["not_established"] and
            claims["zenodo_version_created"] is False and claims["aps_submission"] is False)

    out={"schema":"hdblast-astra-frozen-coordinate-source-review-v1","status":"PASS",
         "manifest_sha256":MANIFEST_SHA,"zip_sha256":ZIP_SHA,"zip_bytes":ZIP_BYTES,
         "approved_source_files":55,"archive_members":56,
         "external_freeze_receipt_sha256":sha(freeze_path.read_bytes()),"approved_assembly_config_sha256":sha(config_path.read_bytes()),
         "review_script_sha256":sha(Path(__file__).read_bytes()),"check_count":len(checks),"checks":checks,
         "roster_counts":expected_counts,"golden_mode_comparisons":mode_comparisons,
         "privacy_distribution_review":{"exact_approved_inventory":True,"public_input_files":18,
           "raw_tool_responses_and_primary_passages_excluded":True,"recognized_credentials_found":False,
           "scope":"Approved corpus, source provenance, member paths and targeted credential signatures; no universal DLP claim."},
         "science_review":{"accepted_Astra_checker_sha256":"dbcb549461a688fbc1bf045e33e07d3dda6ffacf95785984152b8b4f42ea8b82",
           "prior_derivation_reused":True,"science_not_rerun":True,"image_visually_reviewed":True,
           "central_two_zero_sign_pattern":True,"squared_norm_and_ratio_claims_retained":True,
           "EN_degree_four_only":True,"EN2_algebraic_only":True,"U2_unresolved":True},
         "runner_static_review":{"external_pin_bootstrap_before_capsule_execution":True,
           "fresh_isolated_copy_per_process":True,"all_golden_JSON_fields_compared":True,"complete_csv_bytes_compared":True,
           "only_declared_individually_verified_mode_metadata_normalized":True,
           "exact_named_failure_exit_one_required":True,"source_reauthenticated_after_execution":True},
         "fresh_replay_review":"PENDING_ROOT_COMPLETION_SIGNAL",
         "publication_action_taken":False,"frozen_source_mutated":False}
    target=ROOT/"astra/FROZEN_COORDINATE_SOURCE_REVIEW.json"
    target.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"status":"PASS","checks":len(checks),"source_files":55,"archive_members":56,"review":str(target)}))


if __name__ == "__main__":
    main()
