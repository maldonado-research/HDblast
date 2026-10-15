#!/usr/bin/env python3
"""Independent completed-replay readback; never executes capsule science."""
import ast
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import stat
import zipfile

ROOT = Path(__file__).resolve().parents[1]
NAME = "HDBLAST_COORDINATE_SHELL_SOURCES_20261015"
SOURCE = ROOT / "frozen-source" / NAME
EXTRACTED = ROOT / "fresh-zip-source" / NAME
RUN = ROOT / "fresh-zip-replay"
MANIFEST_SHA = "3ab5b0a7ac9045814b0aa59b6eae0d291a19adc76dd31913fd371203f1c132d5"
ZIP_SHA = "7c171631836cc91e546f74ff7897746175d31de42e38a39955c9f06e6a419ae5"
ZIP_BYTES = 3522306


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def doc(path):
    return json.loads(path.read_bytes())


def at(obj, path):
    for key in path:
        obj=obj[key]
    return obj


def same(a,b):
    if type(a) is not type(b):
        return False
    if isinstance(a,dict):
        return a.keys() == b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a) == len(b) and all(same(x,y) for x,y in zip(a,b))
    return a == b


def main():
    counts=Counter()
    def gate(category, condition, detail=""):
        if not condition:
            raise RuntimeError(category+(": "+detail if detail else ""))
        counts[category]+=1

    def regular(path, base):
        gate("regular_contained_path", path.is_file() and not any(x.is_symlink() for x in [path,*path.parents]) and path.resolve().is_relative_to(base.resolve()),str(path))
        return path.read_bytes()

    def member(record, base):
        rel=record["path"]; pp=PurePosixPath(rel)
        gate("canonical_relative_record", not pp.is_absolute() and ".." not in pp.parts and "\\" not in rel and pp.as_posix()==rel,rel)
        raw=regular(base/rel,base)
        gate("artifact_or_member_identity", len(raw)==record["bytes"] and digest(raw)==record["sha256"],rel)
        return raw

    source_manifest=regular(SOURCE/"MANIFEST.json",SOURCE)
    gate("external_manifest_identity",digest(source_manifest)==MANIFEST_SHA)
    inventory={x["path"]:x for x in json.loads(source_manifest)["files"]}
    gate("source_55_exact_inventory",len(inventory)==55)
    source_bytes={rel:member(entry,SOURCE) for rel,entry in inventory.items()}
    def verify_tree(base):
        names=set()
        for path in base.rglob("*"):
            gate("source_no_symlinks_or_special_files", not path.is_symlink() and (path.is_file() or path.is_dir()),str(path))
            if path.is_file():names.add(path.relative_to(base).as_posix())
        gate("source_tree_exact_56",names==set(inventory)|{"MANIFEST.json"})
        gate("extracted_manifest_identity",regular(base/"MANIFEST.json",base)==source_manifest)
        for rel,raw in source_bytes.items():
            gate("complete_source_unchanged",regular(base/rel,base)==raw,rel)
    verify_tree(SOURCE);verify_tree(EXTRACTED)
    zip_path=SOURCE.with_suffix(".zip")
    zipped=regular(zip_path,zip_path.parent)
    gate("external_zip_identity",len(zipped)==ZIP_BYTES and digest(zipped)==ZIP_SHA)
    with zipfile.ZipFile(zip_path) as z:
        names=z.namelist()
        gate("zip_exact_56",len(names)==len(set(names))==56 and set(names)=={NAME+"/"+x for x in set(inventory)|{"MANIFEST.json"}})
        for info in z.infolist():
            rel=info.filename[len(NAME)+1:]
            gate("zip_regular_unencrypted",not info.is_dir() and not stat.S_ISLNK(info.external_attr>>16) and not info.flag_bits&1)
            gate("zip_extracted_and_frozen_bytes_equal",z.read(info)==regular(EXTRACTED/rel,EXTRACTED)==regular(SOURCE/rel,SOURCE),rel)
    initial_review=regular(ROOT/"astra/FROZEN_COORDINATE_SOURCE_REVIEW.json",ROOT)
    gate("prior_frozen_review_retained",digest(initial_review)=="19ed87e247a95d12901b79756f738d061d7e74f865e3372467316175af4c503d")
    zip_auth=doc(ROOT/"FRESH_ZIP_AUTHENTICATION.json")
    gate("fresh_preexecution_identity_receipt",zip_auth["status"]=="PASS" and zip_auth["archive"]["sha256"]==ZIP_SHA and
         zip_auth["manifest"]["sha256"]==MANIFEST_SHA and zip_auth["extracted_members"]==56 and
         zip_auth["bootstrap_pins_before_execution"]=={"MANIFEST.json":MANIFEST_SHA,"reproduce.py":inventory["reproduce.py"]["sha256"],"verify_manifest.py":inventory["verify_manifest.py"]["sha256"]})

    literal={}
    for node in ast.parse(source_bytes["reproduce.py"]).body:
        if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
            name=node.targets[0].id
            if name in {"CAPSULE","CONFIG_SHA256","CHECKERS","CONTROLS","PYTHON_VERSION","DEPENDENCY_VERSIONS"}:
                literal[name]=ast.literal_eval(node.value)
    checkers={x["name"]:x for x in literal["CHECKERS"]}
    controls={(x["checker"],x["name"]):x for x in literal["CONTROLS"]}
    gate("approved_runner_configuration",literal["CONFIG_SHA256"]==digest((ROOT/"publication/RUNNER_CONFIG.json").read_bytes()) and literal["CAPSULE"]==NAME and len(checkers)==4 and len(controls)==33)
    roster=json.loads(source_bytes["REPLAY_ROSTER.json"])
    expected_counts={"successful_checker_runs":8,"external_corruption_processes_rejected":66,
                     "embedded_wrong_formula_controls_rejected":0,"finite_positive_check_groups_each_mode":105}
    gate("roster_counts",roster["counts"]==expected_counts and roster["modes"]==["normal","optimized"])
    gate("roster_exact_control_reasons",roster["external_controls"]==[
        {"checker":x["checker"],"name":x["name"],"expected_reason":x["reason"],"mechanism":x["mechanism"]} for x in literal["CONTROLS"]])
    gate("roster_exact_positive_groups",roster["positive_processes"]==[
        {"name":c["name"],"script":c["script"],"finite_positive_groups_each_mode":c["groups"],"embedded_controls_each_mode":c["embedded"],
         "finite_positive_check_names":c["check_names"],"embedded_control_names":c["embedded_control_names"]} for c in literal["CHECKERS"]])
    receipt_path=RUN/"REPRODUCTION_RECEIPT.json"
    receipt_raw=regular(receipt_path,RUN);receipt=json.loads(receipt_raw)
    gate("replay_status_and_counts",receipt["status"]=="PASS" and receipt["capsule"]==NAME and all(receipt[k]==v for k,v in expected_counts.items()) and
         receipt["finite_positive_check_groups_by_mode"]=={"normal":105,"optimized":105} and len(receipt["steps"])==74)
    gate("replay_exact_source_authentication",receipt["authentication"]=={"status":"PASS","files_verified":55,"manifest_sha256":MANIFEST_SHA})
    for k,rel in [("runner","reproduce.py"),("integrity_verifier","verify_manifest.py"),("requirements","requirements.txt"),("roster","REPLAY_ROSTER.json")]:
        gate("replay_executed_tool_identity",receipt[k]==inventory[rel],k)
    gate("replay_runtime_and_config",receipt["embedded_configuration_sha256"]==literal["CONFIG_SHA256"] and receipt["dependency_versions"]==literal["DEPENDENCY_VERSIONS"]=={} and receipt["python"]["version"].split()[0]==literal["PYTHON_VERSION"]=="3.12.14")

    expected_order=[]
    for mode in ("normal","optimized"):
        expected_order += [(mode,"positive",c["name"],None) for c in literal["CHECKERS"]]
        expected_order += [(mode,"external_corruption",c["checker"],c["name"]) for c in literal["CONTROLS"]]
    actual_order=[(r["mode"],r["kind"],r["checker"],r.get("control")) for r in receipt["steps"]]
    gate("all_74_processes_exact_once_in_order",actual_order==expected_order and len(actual_order)==len(set(actual_order)))
    artifacts={};processes=[];positives={};diagnostics={};copied_count=0;mutated_count=0;failure_files=[]
    for row in receipt["steps"]:
        c=checkers[row["checker"]];mode=row["mode"];is_control=row["kind"]=="external_corruption"
        control=controls[(row["checker"],row["control"])] if is_control else None
        label=row["control"] if is_control else "positive"
        workrel=mode+"/"+c["name"]+"/"+label;work=RUN/workrel
        expected_args=control.get("arguments",c["arguments"]) if is_control else c["arguments"]
        gate("process_identity_and_invocation",row["script"]==c["script"] and row["flags"]==["-I","-B"]+(["-O"] if mode=="optimized" else []) and row["arguments"]==expected_args,workrel)
        expected_inputs=[c["script"],*c["dependencies"]]
        gate("complete_copied_input_inventory",[x["source"]["path"] for x in row["copied_inputs"]]==expected_inputs and len(expected_inputs)==len(set(expected_inputs)),workrel)
        mutated=None
        if is_control and "target" in control:
            original=source_bytes[control["target"]]
            if "json_path" in control:
                changed=json.loads(original);container=at(changed,control["json_path"][:-1]);key=control["json_path"][-1]
                if "factor" in control:
                    v=container[key];container[key]=[str(F(t)*F(control["factor"])) for t in v] if isinstance(v,list) else str(F(v)*F(control["factor"]))
                else:container[key]=control["value"]
                mutated=(json.dumps(changed,indent=2)+"\n").encode()
            elif "append" in control:
                mutated=original+control["append"].encode()
            else:
                gate("mutation_unique_target",original.decode().count(control["before"])==1,workrel)
                mutated=original.decode().replace(control["before"],control["after"],1).encode()
            gate("actual_independent_copied_input_mutation",mutated!=original and control["target"]!=c["script"] and control["target"] in c["dependencies"])
            mutation=row["mutation"]
            gate("mutation_record_exact",mutation["source"]==inventory[control["target"]] and mutation["recipe"]==control and mutation["independent_checker_unchanged"] is True)
            gate("mutation_actual_bytes",member(mutation["mutated"],RUN)==mutated)
            mutated_count+=1
        elif is_control:
            gate("formula_control_no_file_mutation",row["mutation"] is None)
        expected_files=set()
        for copied in row["copied_inputs"]:
            rel=copied["source"]["path"]
            gate("copied_source_snapshot_identity",copied["source"]==inventory[rel] and copied["copy"]=={"path":workrel+"/"+rel,"bytes":inventory[rel]["bytes"],"sha256":inventory[rel]["sha256"]},workrel+"/"+rel)
            expected=mutated if is_control and control.get("target")==rel else source_bytes[rel]
            gate("actual_copied_input_bytes",regular(work/rel,work)==expected,workrel+"/"+rel)
            expected_files.add(rel);copied_count+=1
        gate("checker_unmodified",regular(work/c["script"],work)==source_bytes[c["script"]],workrel)
        for field in ("log","output"):
            entry=row[field]
            gate("canonical_process_artifact_location",entry["path"]==workrel+("/execution.log" if field=="log" else "/fresh.json"))
            member(entry,RUN)
            gate("unique_artifact_record",entry["path"] not in artifacts)
            artifacts[entry["path"]]=entry
        log=regular(work/"execution.log",work).decode();fresh=regular(work/"fresh.json",work)
        generated=json.loads(fresh);expected_files.update({"execution.log","fresh.json"})
        if is_control:
            gate("exact_expected_mathematical_rejection",row["exit_code"]==1 and row["expected_reason"]==row["reported_failure"]==control["reason"] and
                 row["expected_reason_detected"] is True and row["expected_failure_exit_code"]==1 and row["additional_outputs"]==[])
            gate("recorded_failure_wrapper",generated=={"status":"FAIL","failure":control["reason"],"exit_code":1,
                 "passing_output_absent":True,"complete_diagnostic":log,
                 "result_origin":"Runner-parsed exact named failure; original diagnostic retained in execution.log"})
            if control["failure_format"]=="json":
                failure=json.loads(log)
                gate("exact_JSON_failure_reason",at(failure,control.get("status_path",["status"]))==control.get("status","FAIL") and at(failure,control["failure_path"])==control["reason"])
                if c["name"]=="astra":gate("Astra_injected_error_identity",failure["injected_error"]==control["name"])
            elif control["failure_format"]=="fail_line":
                gate("exact_complete_failure_line",log=="FAIL: "+control["reason"]+"\n")
            else:
                gate("exact_terminal_failure_trace",control["failure_format"]=="traceback" and log.startswith(control.get("diagnostic_prefix","")+"Traceback (most recent call last):\n") and log.splitlines()[-1]==control["exception"]+": "+control["reason"])
            gate("registered_passing_output_absent",not (work/c["output"]).exists())
            # The referee deliberately writes an explicit FAIL JSON to a
            # different --output path; authenticate it as redundant evidence.
            if "--output" in expected_args:
                output_name=expected_args[expected_args.index("--output")+1]
                p=work/output_name
                if p.exists():
                    failure_raw=regular(p,work)
                    gate("explicit_failure_file_not_passing",json.loads(failure_raw)==json.loads(log) and json.loads(failure_raw)["status"]=="FAIL")
                    expected_files.add(output_name)
                    failure_files.append({"path":workrel+"/"+output_name,"bytes":len(failure_raw),"sha256":digest(failure_raw)})
            diagnostics[(mode,c["name"],control["name"])]=log
        else:
            golden_name=c["optimized_golden"] if mode=="optimized" else c["golden"]
            gate("whole_positive_JSON_bytes_equal",row["exit_code"]==0 and row["golden"]==inventory[golden_name] and fresh==source_bytes[golden_name] and same(generated,json.loads(source_bytes[golden_name])))
            gate("original_positive_output_equals_fresh",regular(work/c["output"],work)==fresh)
            expected_files.add(c["output"])
            gate("positive_all_fields_compared",row["scientific_output_matches"] is True and row["complete_generated_json_bytes_equal"] is True and row["ignored_golden_fields"]==[])
            gate("positive_named_group_inventory",row["finite_positive_check_groups"]==c["groups"] and row["finite_positive_check_names"]==c["check_names"] and row["embedded_wrong_formula_controls"]==c["embedded"]==0 and row["embedded_control_names"]==[])
            gate("actual_positive_status",at(generated,c["status_path"])==c["status"])
            if c["checks"] is not None:
                actual_groups=at(generated,c["checks"]["path"])
                gate("actual_named_group_inventory",actual_groups==c["check_names"] and len(actual_groups)==len(set(actual_groups))==c["groups"])
            for item in c["required_fields"]:
                gate("required_scientific_field",same(at(generated,item["path"]),item["equals"]))
            gate("additional_output_inventory",len(row["additional_outputs"])==len(c["additional_outputs"]))
            for actual_add,(filename,golden) in zip(row["additional_outputs"],c["additional_outputs"]):
                gate("whole_generated_CSV_bytes_equal",actual_add["golden"]==inventory[golden] and actual_add["complete_bytes_equal"] is True and
                     actual_add["output"]["path"]==workrel+"/"+filename and member(actual_add["output"],RUN)==source_bytes[golden])
                entry=actual_add["output"]
                gate("unique_additional_artifact",entry["path"] not in artifacts)
                artifacts[entry["path"]]=entry;expected_files.add(filename)
            positives[(mode,c["name"])]=generated
        files=set()
        for p in work.rglob("*"):
            gate("process_no_symlink_or_special",not p.is_symlink() and (p.is_file() or p.is_dir()))
            if p.is_file():files.add(p.relative_to(work).as_posix())
        gate("exact_process_file_inventory_no_stale_extras",files==expected_files,workrel+" extras="+repr(files-expected_files))
        processes.append({"mode":mode,"checker":c["name"],"control":label if is_control else None,"exit_code":row["exit_code"],
                          "status":"EXPECTED_REJECTION" if is_control else "PASS","complete_output_verified":True,
                          "expected_reason":control["reason"] if is_control else None,
                          "copied_inputs_verified":len(row["copied_inputs"]),"input_mutation_verified":mutated is not None})

    mode_agreements=[]
    for name,c in checkers.items():
        left=deepcopy(positives[("normal",name)]);right=deepcopy(positives[("optimized",name)])
        for field in c["mode_fields"]:
            gate("actual_individually_verified_mode_metadata",same(at(left,field["path"]),field["normal"]) and same(at(right,field["path"]),field["optimized"]))
            for output in [left,right]:del at(output,field["path"][:-1])[field["path"][-1]]
        gate("complete_cross_mode_scientific_JSON",same(left,right),name)
        mode_agreements.append({"checker":name,"scientific_json_equal":True,"only_normalized_fields":[x["path"] for x in c["mode_fields"]]})
    gate("mode_agreement_receipt_complete",receipt["mode_agreements"]==mode_agreements)
    for name,control in controls.items():
        gate("complete_failure_diagnostic_equal_across_modes",diagnostics[("normal",*name)]==diagnostics[("optimized",*name)],str(name))
    gate("all_150_designated_artifacts",len(artifacts)==150)
    gate("all_6_data_mutations_and_14_failure_files",mutated_count==6 and len(failure_files)==14)
    recorder_path=ROOT/"publication/REPLAY_TOOL_VALIDATION.json"
    recorder_raw=regular(recorder_path,ROOT);recorder=json.loads(recorder_raw)
    rec_artifacts={x["path"]:x for x in recorder["artifact_records"]}
    gate("independent_recorder_artifact_inventory_agrees",recorder["status"]=="PASS" and len(rec_artifacts)==len(recorder["artifact_records"])==150 and rec_artifacts==artifacts and
         recorder["fresh_artifact_files_excluding_receipt"]==150 and recorder["fresh_processes"]==74 and recorder["replay_counts"]==expected_counts)
    gate("independent_recorder_receipt_identity_agrees",recorder["replay_receipt"]["sha256"]==digest(receipt_raw) and recorder["replay_receipt"]["bytes"]==len(receipt_raw))
    # Re-read both immutable trees and archive after the entire artifact audit.
    verify_tree(SOURCE);verify_tree(EXTRACTED)
    gate("archive_unchanged_after_readback",regular(zip_path,zip_path.parent)==zipped)
    gate("prior_reviews_not_modified",regular(ROOT/"astra/FROZEN_COORDINATE_SOURCE_REVIEW.json",ROOT)==initial_review)
    gate("receipt_not_modified_during_audit",regular(receipt_path,RUN)==receipt_raw)
    result={"schema":"hdblast-astra-completed-coordinate-replay-review-v1","status":"PASS",
            "source_manifest_sha256":MANIFEST_SHA,"archive":{"bytes":ZIP_BYTES,"sha256":ZIP_SHA},
            "source_files":55,"archive_members":56,"actual_processes":74,**expected_counts,
            "designated_artifacts_verified":150,"additional_explicit_failure_JSON_files_verified":14,
            "copied_input_instances_verified":copied_count,"intentional_data_mutations_verified":mutated_count,
            "positive_generated_CSVs_whole_bytes_verified":2,"positive_JSONs_whole_bytes_verified":8,
            "full_failure_diagnostics_verified":66,"mode_agreements":mode_agreements,
            "check_count":sum(counts.values()),"check_categories":dict(sorted(counts.items())),
            "replay_receipt":{"bytes":len(receipt_raw),"sha256":digest(receipt_raw)},
            "separate_recorder":{"bytes":len(recorder_raw),"sha256":digest(recorder_raw),"complete_artifact_inventory_agrees":True},
            "previous_frozen_review_sha256":digest(initial_review),"audit_source_sha256":digest(Path(__file__).read_bytes()),
            "process_reviews":processes,"artifact_records":[artifacts[x] for x in sorted(artifacts)],
            "additional_failure_JSON_records":failure_files,
            "source_and_extracted_tree_unchanged":True,"science_rerun":False,"capsule_code_imported_or_executed":False,
            "existing_reviews_modified":False,"publication_actions":False,
            "scope":"Independent read-only audit of actual completed fresh replay, exact artifact bytes, declared process/control identities and complete normal/optimized scientific outputs. Reuses accepted mathematics; does not rerun science or prove inherited analytic PDE hypotheses."}
    output=ROOT/"astra/FRESH_COORDINATE_REPLAY_REVIEW.json"
    output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":"PASS","checks":sum(counts.values()),"processes":74,"designated_artifacts":150,
                      "copied_inputs":copied_count,"intentional_mutations":mutated_count,"review":str(output)}))


if __name__ == "__main__":
    main()
