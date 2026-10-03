"""Read-only postresult packaging and fresh-replay audit; never run a source."""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path,PurePosixPath
import re
import stat
import zipfile

ROUTES=("primary","independent")


def require(condition,message):
    if not condition:raise ValueError(message)


def load(raw):
    def unique(items):
        output={}
        for key,value in items:
            require(key not in output,"Duplicate JSON key")
            output[key]=value
        return output
    def invalid(value):raise ValueError("Nonfinite JSON constant")
    return json.loads(raw,object_pairs_hook=unique,parse_constant=invalid)


def sha(raw):return hashlib.sha256(raw).hexdigest()


def canonical(value):return json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()


def frame(envelope):
    evidence=dict(envelope["method_evidence"])
    evidence.pop("resource",None);evidence.pop("construction_benchmark",None)
    return {"payload":envelope["payload"],"method_evidence":evidence,
            "semantic_validation":envelope["semantic_validation"]}


def fixed_source_events(envelope):
    chronology=envelope["execution"]["source_chronology"]
    require(chronology["source_bundle_constructions"]==128 and chronology["archive_arrays_decoded"]==0 and
            len(chronology["events"])==128,"Actual source chronology incomplete")
    expected={(source,Fraction(-9,2)+Fraction(2*j+1,128))
              for source in ("positive_B","signed_uB") for j in range(64)}
    actual={(item["source"],Fraction(item["center"])) for item in chronology["events"]}
    require(actual==expected,"Actual source callback universe differs")


def verify_outer(outer,raw,envelope,route,receipt_sha,registration_sha,contract):
    require(outer["status"]=="PASS_BOUNDED_ROUTE_EXECUTION" and outer["fabricated_only"] is False and
            outer["exit_code"]==0 and outer["timed_out"] is False and outer["uid"]!=0 and
            outer["route"]==envelope["route"]==route,"Actual bounded route evidence differs")
    require(outer["remote_receipt_sha256"]==receipt_sha and outer["registration_sha256"]==registration_sha and
            outer["output_sha256"]==sha(raw) and outer["output_bytes"]==len(raw),"Actual output/pin binding differs")
    resources=contract["resources"]
    require(outer["wall_seconds"]<=resources["each_route_wall_seconds"] and
            outer["peak_rss_kib_wait4"]<=resources["each_route_peak_rss_kib"],"Actual whole-route budget exceeded")
    fixed_source_events(envelope)


def read_archive(path,expected_archive_sha=None):
    raw=path.read_bytes();archive_sha=sha(raw)
    if expected_archive_sha is not None:
        require(re.fullmatch("[0-9a-f]{64}",expected_archive_sha or "") and archive_sha==expected_archive_sha,
                "Expected complete archive hash differs")
    with zipfile.ZipFile(path) as archive:
        infos=archive.infolist();names=[item.filename for item in infos]
        require(len(names)==len(set(names)),"Duplicate ZIP member")
        for item in infos:
            name=item.filename
            parsed=PurePosixPath(name.rstrip("/"))
            require(not parsed.is_absolute() and parsed.as_posix()==name.rstrip("/") and
                    "\\" not in name and not any(part in ("",".","..") for part in parsed.parts),"Unsafe ZIP path")
            require(not item.flag_bits&1,"Encrypted ZIP member")
            mode=item.external_attr>>16
            require(not stat.S_ISLNK(mode),"ZIP symlink member")
        require(archive.testzip() is None,"ZIP CRC failure")
        files={item.filename:archive.read(item) for item in infos if not item.is_dir()}
        manifests=[name for name in files if PurePosixPath(name).name=="MANIFEST.json"]
        require(manifests,"Root archive manifest absent")
        manifests.sort(key=lambda name:(len(PurePosixPath(name).parts),name))
        chosen=manifests[0];prefix=chosen[:-len("MANIFEST.json")]
        manifest=load(files[chosen]);require(manifest["schema_version"]==1 and set(manifest)=={"schema_version","files"},"Manifest schema differs")
        members=manifest["files"]
        require("MANIFEST.json" not in members,"Circular whole-payload manifest")
        require(set(files)=={prefix+name for name in members}|{chosen},"Complete ZIP/manifest membership differs")
        payload={name:files[prefix+name] for name in members}
        for name,item in members.items():
            require(set(item)=={"sha256","bytes"} and item["bytes"]==len(payload[name]) and
                    item["sha256"]==sha(payload[name]),"ZIP manifest binding differs: "+name)
        return {"archive_sha256":archive_sha,"archive_bytes":len(raw),"archive_members":len(infos),
                "archive_file_members":len(files),"manifest_sha256":sha(files[chosen]),
                "manifest":manifest,"manifest_raw":files[chosen],"payload":payload,"prefix":prefix}


def audit(package_directory,fresh_directory,output_directory,archive_path=None,expected_archive_sha=None):
    package=Path(package_directory).absolute();fresh=Path(fresh_directory).absolute();out=Path(output_directory).absolute()
    require(package.is_dir() and fresh.is_dir(),"Sealed package and fresh-reproduction directories required")
    require(not out.exists() and package!=out and package not in out.parents and fresh!=out and fresh not in out.parents,
            "Fresh external audit directory required")
    archives=[Path(archive_path).absolute()] if archive_path else list(package.glob("*.zip"))
    require(len(archives)==1,"Select exactly one complete ZIP archive explicitly")
    archive=read_archive(archives[0],expected_archive_sha);payload=archive["payload"]
    registration_raw=payload["FULL_REGISTRATION.json"];registration=load(registration_raw)
    receipt_raw=payload["FREEZE_RECEIPT.json"];receipt=load(receipt_raw)
    registration_sha=sha(registration_raw);receipt_sha=sha(receipt_raw)
    require(receipt["status"]=="PASS_REMOTE_REGISTERED_SOURCE_GO" and
            receipt["registration_sha256"]==registration_sha,"Frozen registration/readback relationship differs")
    for name,item in registration["files"].items():
        require(name in payload and len(payload[name])==item["bytes"] and sha(payload[name])==item["sha256"],
                "Immutable frozen registration member differs: "+name)
    contract=load(payload["REGISTRATION_CONTRACT.json"])
    require(contract==registration["frozen_configuration"],"Frozen configuration differs")
    baseline={}
    for route in ROUTES:
        raw=payload[f"results/{route}/OUTPUT.json"];envelope=load(raw)
        outer=load(payload[f"results/{route}/EXECUTION.json"])
        verify_outer(outer,raw,envelope,route,receipt_sha,registration_sha,contract)
        require(envelope["execution"]["freeze_commit"]==receipt["freeze_commit"],"Archived freeze identity differs")
        baseline[route]={"raw":raw,"envelope":envelope,"outer":outer}
    reports=sorted(fresh.rglob("FRESH_REPLAY.json"));require(reports,"Fresh replay receipt not yet available")
    replay_checks=[]
    for report_path in reports:
        report=load(report_path.read_bytes());require(report["status"]=="PASS_EXACT_FRESH_REPLAY","Fresh replay reported failure")
        require(report["manifest_sha256"]==archive["manifest_sha256"] and
                report["registration_sha256"]==registration_sha and report["freeze_receipt_sha256"]==receipt_sha,
                "Fresh replay archive/freeze pins differ")
        extracted=Path(report["checkpoint"]).absolute()
        require(extracted.is_dir() and (fresh==extracted or fresh in extracted.parents),
                "Fresh receipt did not point to independently extracted payload")
        require((extracted/"MANIFEST.json").read_bytes()==archive["manifest_raw"],"Extracted manifest bytes differ")
        actual={p.relative_to(extracted).as_posix() for p in extracted.rglob("*") if p.is_file()}
        require(actual==set(payload)|{"MANIFEST.json"},"Extracted payload membership differs")
        for name,raw in payload.items():
            target=extracted/name
            require(not target.is_symlink() and target.read_bytes()==raw,"Extracted archive payload bytes differ: "+name)
        comparisons={}
        for route in ROUTES:
            directory=report_path.parent/route
            raw=(directory/"OUTPUT.json").read_bytes();envelope=load(raw);outer=load((directory/"EXECUTION.json").read_bytes())
            verify_outer(outer,raw,envelope,route,receipt_sha,registration_sha,contract)
            require(envelope["execution"]["freeze_commit"]==receipt["freeze_commit"],"Fresh freeze identity differs")
            old=baseline[route]["envelope"]
            old_frame=canonical(frame(old));new_frame=canonical(frame(envelope))
            old_payload=canonical(old["payload"]);new_payload=canonical(envelope["payload"])
            require(old_frame==new_frame and old_payload==new_payload,"Independent exact replay comparison failed: "+route)
            declared=report["comparisons"][route]
            require(declared["payload_exact"] is True and declared["stable_frame_exact"] is True and
                    declared["archived_payload_sha256"]==declared["fresh_payload_sha256"]==sha(new_payload) and
                    declared["archived_stable_frame_sha256"]==declared["fresh_stable_frame_sha256"]==sha(new_frame),
                    "Fresh replay declared frame hashes differ")
            comparisons[route]={"payload_exact":True,"stable_frame_exact":True,"payload_sha256":sha(new_payload),
              "stable_frame_sha256":sha(new_frame),"fresh_wall_seconds":outer["wall_seconds"],
              "fresh_peak_rss_kib_wait4":outer["peak_rss_kib_wait4"],"archived_output_sha256":sha(baseline[route]["raw"]),
              "fresh_output_sha256":sha(raw)}
        require(report["frozen_payload_bytes_unchanged"] is True and report["scientific_status_exact"] is True and
                report["archived_scientific_status"]==report["fresh_scientific_status"],"Scientific status was changed by replay")
        replay_checks.append({"receipt":str(report_path),"receipt_sha256":sha(report_path.read_bytes()),
             "optimized":report["optimized"],"scientific_status":report["fresh_scientific_status"],"comparisons":comparisons})
    result={"schema_version":1,"status":"PASS_INDEPENDENT_POSTRESULT_PACKAGE_AND_FRESH_REPLAY_AUDIT",
         "audit_scope":"Read-only packaging and reproduction; no new scientific calculation or numerical source imports.",
         "source_calls_by_auditor":0,"scientific_replays_by_auditor":0,"retained_array_decodes":0,
         "archive":str(archives[0]),"archive_sha256":archive["archive_sha256"],"archive_bytes":archive["archive_bytes"],
         "archive_members":archive["archive_members"],"payload_file_count":len(payload),
         "manifest_sha256":archive["manifest_sha256"],"registration_sha256":registration_sha,
         "freeze_receipt_sha256":receipt_sha,"freeze_commit":receipt["freeze_commit"],
         "all_ZIP_CRCs_pass":True,"exact_archive_and_extracted_membership":True,"all_payload_bytes_and_frozen_registered_bytes_match":True,
         "fresh_replays":replay_checks,"limitations":"This verifies preservation/reproduction of the registered source/operator result. It does not certify the twelve-case action-pressure ledger, physical initial state, continuum momentum target, coupled gravity or a higher-dimensional Big Bang cause."}
    out.mkdir(parents=True)
    (out/"POSTRESULT_AUDIT.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    lines=["# Independent postresult package and fresh-replay audit","","Status: PASS.","",
           f"The ZIP has {archive['archive_members']} members and {len(payload)} payload files plus its manifest. All CRCs, exact memberships, byte sizes and SHA256 pins pass. Extracted bytes match the ZIP, and every frozen registered member remains unchanged.","",
           f"Archive SHA256: `{archive['archive_sha256']}`.","",
           "Both routes' actual bounded receipts and128-source construction histories were checked. Fresh canonical payloads and complete stable scientific frames reproduce exactly. Negative scientific states would remain negative.","",
           "The auditor performed no source evaluation, scientific replay, array decoding, checkout mutation or remote write. Full pressure/contact and cosmological claims remain outside this audit.",""]
    (out/"POSTRESULT_AUDIT.md").write_text("\n".join(lines))
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--package-directory",required=True);parser.add_argument("--fresh-reproduction-directory",required=True)
    parser.add_argument("--output-directory",required=True);parser.add_argument("--archive")
    parser.add_argument("--expected-archive-sha256")
    args=parser.parse_args()
    result=audit(args.package_directory,args.fresh_reproduction_directory,args.output_directory,args.archive,args.expected_archive_sha256)
    print(json.dumps({key:result[key] for key in ("status","archive_sha256","payload_file_count","freeze_commit")},indent=2))


if __name__=="__main__":main()
