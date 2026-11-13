"""Forty-seven closed-capsule rejection controls; not continuous-proof validation."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def rehash_manifest(root):
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != root / "MANIFEST.json":
            raw = path.read_bytes()
            files[path.relative_to(root).as_posix()] = {"bytes":len(raw), "sha256":sha(raw)}
    value = {"schema":"HDBLAST_ALL_MODE_INVERSE_MANIFEST_V1", "files":files}
    (root / "MANIFEST.json").write_text(json.dumps(value,sort_keys=True,indent=2)+"\n", encoding="utf-8")
    return sha((root / "MANIFEST.json").read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    script = Path(__file__).absolute()
    need(script.is_file() and all(not p.is_symlink() for p in (script,*script.parents)),
         "Regular controls script path required")
    root = script.resolve(strict=True).parent.parent
    out = args.output_dir.absolute()
    need(".." not in out.parts and not out.is_relative_to(root), "Fresh controls directory must be outside capsule")
    need(not out.exists() and not out.is_symlink(), "Controls directory exists")
    need(out.parent.is_dir() and all(not p.is_symlink() for p in out.parents),
         "Regular existing controls parent required")
    out.mkdir()

    def invoke(case_root, pin, name):
        report = out / (name+".json")
        result = subprocess.run([sys.executable,"-I","-B",
            str(case_root / "checks/check_checkpoint.py"), "--manifest-sha256",pin,
            "--output",str(report)],capture_output=True,text=True,timeout=120)
        return result, report

    result, report = invoke(root,args.manifest_sha256,"positive")
    need(result.returncode == 0 and report.is_file(), "Baseline checker failed")
    baseline = json.loads(report.read_text(encoding="utf-8"))
    need(baseline["status"] == "PASS" and baseline["manifest_sha256"] == args.manifest_sha256,
         "Baseline report failed")
    cases = []

    def copy(name):
        target = out / "fixtures" / name
        shutil.copytree(root,target)
        return target

    def reject(name, fixture, pin, expected, **details):
        result, report = invoke(fixture,pin,name)
        lines = result.stderr.strip().splitlines()
        diagnostic = lines[-1] if lines else ""
        need(result.returncode != 0 and not report.exists()
             and diagnostic.startswith("ValueError: ") and expected in diagnostic,
             "Control did not reach intended rejection: "+name+"; diagnostic="+diagnostic)
        cases.append({"name":name,"exit_code":result.returncode,"report_created":False,
                      "status":"REJECTED","diagnostic":diagnostic,**details})

    source = Path("inputs/SCALAR_COERCIVITY.md")
    changed = copy("changed_source")
    with (changed / source).open("ab") as stream:
        stream.write(b"\n")
    reject("changed_source",changed,args.manifest_sha256,"Source bytes differ: "+source.as_posix())

    extra = copy("extra_file")
    (extra / "UNLISTED.txt").write_text("unlisted\n",encoding="utf-8")
    reject("extra_file",extra,args.manifest_sha256,"Closed file roster mismatch")

    missing = copy("missing_file")
    (missing / source).unlink()
    reject("missing_file",missing,args.manifest_sha256,"Closed file roster mismatch")

    linked = copy("symlink_source")
    target = out / "symlink_input_target.json"
    target.write_bytes((linked / source).read_bytes())
    (linked / source).unlink()
    (linked / source).symlink_to(target)
    reject("symlink_source",linked,args.manifest_sha256,"Capsule symlink rejected")

    stale = copy("stale_manifest_pin")
    with (stale / "MANIFEST.json").open("ab") as stream:
        stream.write(b"\n")
    reject("stale_manifest_pin",stale,args.manifest_sha256,"Manifest SHA256 mismatch")

    from fractions import Fraction

    def mathematical_mutant(suite,field,value):
        name='rehashed_'+suite+'_'+field
        mutant=copy(name)
        data_path=mutant/'DATA.json'
        data=json.loads(data_path.read_text(encoding='utf-8'))
        before=data[suite][field]
        need(before!=value,'Mutation does not change a claim')
        data[suite][field]=value
        data_path.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n',encoding='utf-8')
        mutant_pin=rehash_manifest(mutant)
        reject(name,mutant,mutant_pin,'Arithmetic failed: '+suite+':',
               rehashed_manifest_sha256=mutant_pin,
               mutation={'suite':suite,'field':field,'before':before,'after':value})

    for field,value in [
      ('norm_derivative_q_power','1/2'),('scalar_multiplier_q_power',8),
      ('momentum_multiplier_q_power',4),('scalar_joint_symbol_sum_squared_factor','1'),
      ('scalar_joint_symbol_c2tau2_factor','0'),('momentum_transverse_factor','1'),
      ('momentum_longitudinal_factor','2')]:
        mathematical_mutant('scaling',field,value)

    baseline_data=json.loads((root/'DATA.json').read_text(encoding='utf-8'))
    for suite in ('identity','momentum'):
        for field,value in baseline_data[suite].items():
            if field!='schema':
                mathematical_mutant(suite,field,str(Fraction(value)+1))

    for field,value in [
      ('alpha1_factor','-1'),('vectorfield_norm_ceiling','1/1000000000000'),
      ('R_degree_ceiling',1),('R_norm_ceiling','1/100000000000000000000'),
      ('R_weighted_order_ceiling',1),('A_left_coefficient_factor',2)]:
        mathematical_mutant('coefficients',field,value)

    for field,value in [
      ('scalar_gap_mantissa','13/100'),('scalar_inverse_L2_ten_ceiling',112),
      ('momentum_inverse_L2_mantissa','18/25'),
      ('momentum_H1_API_ceiling','10000000000000000000'),
      ('graph_induction_ceiling',100000000),('graph_absorption_delta_ten_exponent',-161),
      ('norm_bridge_integer','7656014335380773516240486400'),
      ('momentum_D_step_factor',100),('momentum_A_step_factor',16000),
      ('momentum_graph_K_exponent',170),('momentum_A7_exponent',2436),
      ('momentum_physical_W14_exponent',2649),('combined_Hilbert_inverse_exponent',2700)]:
        mathematical_mutant('arithmetic',field,value)

    need(len(cases) == 47, "Forty-seven rejection controls required")
    receipt = {"schema":"HDBLAST_ALL_MODE_INVERSE_CONTROLS_V1", "status":"PASS",
        "manifest_sha256":args.manifest_sha256,"baseline_producer_executions":1,
        "rejected_children":len(cases),"cases":cases,
        "scope":"Five protocol faults and forty-two freshly rehashed mathematical claim faults; each child authenticates a closed capsule and rebuilds all exact oracles. No continuous analytic proof or nonlinear amplitude validation."}
    with (out / "CONTROLS.json").open("x",encoding="utf-8") as stream:
        stream.write(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"status":"PASS","rejected_children":len(cases)},sort_keys=True))


if __name__ == "__main__":
    main()
