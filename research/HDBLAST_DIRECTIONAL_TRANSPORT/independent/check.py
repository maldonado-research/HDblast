"""Repeat the directional deficit integral with the time and angle integrations exchanged."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE.parent.parent / "HDBLAST_COHERENT_DECAY_20261120"


INPUT_HASHES = {
    "coupled/joined.py": "ceb5ad0d36b0d8f77fffd20590d48495089a69cb169c24c6ec6ac8b0706e2d4a",
    "envelope/build.py": "ac63c8d153941ee1383c57fae573c94af026e0d0355ea7f68bfebef08e77b0e0",
    "envelope/ENVELOPE.json": "08e42dbe58297862437b2b331cd544d81fe7ad844f0c33ee84f84027d6970576",
    "kinetics/sweep.py": "5c6255a1255e71b923ee7c21e450c5e09032699e43cb6fb54faf2e721b6f8f72",
    "kinetics/sweep.hpp": "9b17d45f8cf2cdd5ec1365cac14dc8f1a171eef2daf715c401ac7b437f948e65",
}


def run_case(source, release, case):
    engine = runpy.run_path(str(source / "coupled" / "joined.py"))["run"]
    transform = runpy.run_path(str(source / "envelope" / "build.py"))["transform"]
    make_source = runpy.run_path(str(source / "kinetics" / "sweep.py"))["make_source"]
    old_header = (source / "kinetics" / "sweep.hpp").read_text()
    base = (source / "coupled" / "joined_forward.cpp").read_text()
    params = release["parameters"][case]
    centered = params["daughter_mass"] >= 1e11
    evolved = transform(base, centered=centered)
    key = "centered_solver_sha256" if centered else "transformed_solver_sha256"
    if hashlib.sha256(evolved.encode()).hexdigest() != release[key]:
        raise ArithmeticError("Archived evolution source hash mismatch")

    rest = release["runs"][case]["end_rest"]
    screen = release["screen"][case]
    floor = screen["elastic_drag_action_floor"]
    action = 5*floor if centered else screen["witness_action"]
    solver = make_source(evolved, rest, action, floor)
    if solver.count(old_header) != 1:
        raise ArithmeticError("Unexpected upstream collector insertion")
    solver = solver.replace(old_header, (HERE / "reverse.hpp").read_text())
    rows = {}
    with tempfile.TemporaryDirectory() as directory:
        cpp = Path(directory) / "measure.cpp"
        exe = Path(directory) / "measure"
        cpp.write_text(solver)
        subprocess.run(["g++","-std=c++17","-O3",str(cpp),"-o",str(exe)], check=True)
        for name, scale in [("native",1.), ("half_step",.5)]:
            result = engine(exe,dict(params,sample=0,cfl=params["cfl"]*scale))
            if not math.isclose(result["end_rest"],rest,rel_tol=3e-8 if scale==1 else 3e-4):
                raise ArithmeticError("Trajectory no longer matches the archived evolution")
            rows[name] = {k: v for k,v in result.items() if "sweep_" in k}
            for prefix, expected in [("main_sweep",.98),("floor_sweep",.98),
                                     ("extended_main_sweep",.99998),("extended_floor_sweep",.99998)]:
                if abs(rows[name][prefix+"_fraction_tested"]-expected)>1e-8:
                    raise ArithmeticError("Incomplete source interval")
    return rows


def check_peak_witness(main):
    """Evaluate the reported instantaneous peak state at the angular kinks."""
    c = main["peak_radial_speed_over_request"]
    d = main["peak_signed_field_over_request"]
    knots = {-1.0, 1.0}
    if d != 0:
        for edge in (-1,0,1):
            x = (c-edge)/d
            if -1 < x < 1:
                knots.add(x)
    knots = sorted(knots)
    values = [max(0.0, 1-abs(c-d*x)) for x in knots]
    area = math.fsum((b-a)*(fa+fb)/4 for a,b,fa,fb in
                     zip(knots,knots[1:],values,values[1:]))
    if abs(area - main["peak_instantaneous_unfilled_fraction"]) > 1e-8:
        raise ArithmeticError("Reported instantaneous peak state is inconsistent")


def compare(checks, path):
    """Compare all mass, action, grid and fraction-window combinations."""
    other = json.loads(path.read_text())["runs"]
    differences = []
    for case, grids in checks.items():
        for grid, metrics in grids.items():
            for window, prefix in [("central98", ""), ("central99_998", "extended_")]:
                for action, collector in [("witness","main"), ("drag_floor","floor")]:
                    main = other[case]["grids"][grid][window][action]["isotropic"]
                    check_peak_witness(main)
                    key = f"{prefix}{collector}_sweep_unmet_fraction"
                    error = abs(main["unfilled_fraction_of_total"] - metrics[key])
                    differences.append(error)
                    if error > 1e-8:
                        raise ArithmeticError(f"Directional fraction differs in {case}/{grid}/{window}/{action}: {error}")
    return max(differences)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,default=DEFAULT_SOURCE,
                        help="Path to the unchanged HDblast coherent-decay study")
    parser.add_argument("--output",type=Path)
    parser.add_argument("--compare",type=Path,help="Compare with the angle-first results.json")
    parser.add_argument("--cases",nargs="+",choices=["mass_4e10","mass_4e12"],
                        default=["mass_4e10","mass_4e12"])
    args = parser.parse_args()
    for path, digest in INPUT_HASHES.items():
        if hashlib.sha256((args.source/path).read_bytes()).hexdigest() != digest:
            raise ValueError("Unexpected archived file " + path)
    release = json.loads((args.source / "envelope" / "ENVELOPE.json").read_text())
    checks={case:run_case(args.source,release,case) for case in args.cases}
    if args.compare is not None:
        print(f"Maximum absolute difference against angle-first integration: {compare(checks,args.compare):.8g}", file=sys.stderr)
    report=json.dumps(checks,indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(report)
    else:
        print(report,end="")
