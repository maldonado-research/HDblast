#!/usr/bin/env python3
"""Exercise independent mathematical and complete-data rejection controls."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT=Path(__file__).resolve().parents[1]
CHECKER=Path(__file__).resolve().with_name("check_coordinate_shell_sources.py")
MATHEMATICAL={
    "wrong_scaled_force":"wrong_scaled_force: ell/mu factor missing",
    "omit_triangular_half":"omit_triangular_half: integral of 2*(x-s) is x^2",
    "wrong_negative_mu":"wrong_negative_mu: signed endpoint witness excluded",
    "omit_resolvent":"omit_resolvent: contraction residual not closed",
    "wrong_metric_radius":"wrong_metric_radius: u=rho^-2*C*y",
    "wrong_angular_normalization":"wrong_angular_normalization: sphere Y^4 mean is 3/40",
    "omit_measure_clock":"omit_measure_clock: positive Taylor term of cosh^3 excluded",
}


def run_case(case,optimized,root):
    cmd=[sys.executable,"-I"]+(["-O"] if optimized else [])+["-B",str(CHECKER),"--dataset-root",str(root)]
    if case in MATHEMATICAL:
        cmd += ["--inject",case]
        reason=MATHEMATICAL[case]
    elif case=="last_cell_upper_endpoint":
        reason="primary all-cell endpoint mismatch: row 4352 Sp_hi"
    elif case=="last_norm_upper_endpoint":
        reason="JSON mathematical value mismatch: /exact/norms/EN_material/squared_L2/1"
    else:
        raise RuntimeError("unknown control")
    process=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=False)
    if process.returncode!=1 or process.stdout!="" or process.stderr!="FAIL: "+reason+"\n":
        raise RuntimeError("control did not fail for intended reason: "+case+" "+process.stderr[-1000:])
    return {"control":case,"optimization":int(optimized),"returncode":process.returncode,
            "stdout_bytes":len(process.stdout),"stderr":process.stderr,"expected_reason":reason,"status":"PASS_EXPECTED_FAILURE"}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=CHECKER.with_name("SOL_CONTROLS.json"))
    args=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix="hdblast-sol-controls-") as td:
        work=Path(td)
        roots={}
        for kind in ("last_cell_upper_endpoint","last_norm_upper_endpoint"):
            copy=work/kind
            copy.mkdir()
            (copy/"inputs").symlink_to(ROOT/"inputs",target_is_directory=True)
            (copy/"compute_coordinate_sources.py").symlink_to(ROOT/"compute_coordinate_sources.py")
            shutil.copytree(ROOT/"computed",copy/"computed")
            roots[kind]=copy
        path=roots["last_cell_upper_endpoint"]/"computed/COORDINATE_SHELL_SOURCE_CELLS.csv"
        rows=list(csv.DictReader(io.StringIO(path.read_text())))
        names=list(rows[0])
        rows[-1]["Sp_hi"]=str(int(rows[-1]["Sp_hi"])+1)
        with path.open("w",newline="") as out:
            writer=csv.DictWriter(out,fieldnames=names,lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        path=roots["last_norm_upper_endpoint"]/"computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json"
        report=json.loads(path.read_text())
        old=Fraction(report["exact"]["norms"]["EN_material"]["squared_L2"][1])
        report["exact"]["norms"]["EN_material"]["squared_L2"][1]=str(old+1)
        path.write_text(json.dumps(report,indent=2)+"\n")
        cases=[(name,mode,ROOT) for name in MATHEMATICAL for mode in (False,True)]
        cases += [(name,mode,root) for name,root in roots.items() for mode in (False,True)]
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(run_case,*case) for case in cases]
            results=[future.result() for future in futures]
    receipt={"schema":"hdblast-sol-coordinate-controls-v1","status":"PASS","modes":[0,1],
             "mathematical_controls":len(MATHEMATICAL),"complete_data_controls":2,"runs":len(results),
             "checker_sha256":hashlib.sha256(CHECKER.read_bytes()).hexdigest(),
             "runner_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"results":results}
    args.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS","runs":len(results)}))


if __name__=="__main__":
    main()
