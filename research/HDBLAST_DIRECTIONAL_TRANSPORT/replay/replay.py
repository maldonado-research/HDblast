"""Replay the pinned HDblast histories and compare angular and summed birth-state fluxes."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
PIN = "dc99ff5c5894abcff3e11e02e4b97e19da31de00"
HASHES = {
    "kinetics/sweep.py": "5c6255a1255e71b923ee7c21e450c5e09032699e43cb6fb54faf2e721b6f8f72",
    "kinetics/sweep.hpp": "9b17d45f8cf2cdd5ec1365cac14dc8f1a171eef2daf715c401ac7b437f948e65",
    "kinetics/SWEEP.json": "0f3421ac6d30927ba32f3c8a8a431b591c27fb3307a161c7740c77cda0c0caad",
    "coupled/joined.py": "ceb5ad0d36b0d8f77fffd20590d48495089a69cb169c24c6ec6ac8b0706e2d4a",
    "coupled/joined_forward.cpp": "1ce3d03884da0bfc1c31b749d1000f58de1bccbbe9f48130ab27871bbf0ef236",
    "envelope/build.py": "ac63c8d153941ee1383c57fae573c94af026e0d0355ea7f68bfebef08e77b0e0",
    "envelope/ENVELOPE.json": "08e42dbe58297862437b2b331cd544d81fe7ad844f0c33ee84f84027d6970576",
}
WINDOWS = {"central98": "main_sweep", "central99_998": "extended_main_sweep"}
# Native-endian float64 records. dk is the original accepted-step absolute shell speed.
COLUMNS = "t dt k0 k1 d0 d1 drest f0 f1 dk k target".split()


def compiled(cpp, binary, include=None):
    subprocess.run(["g++", "-O3", "-std=c++17", *(["-I", str(include)] if include else []), str(cpp), "-o", str(binary)], check=True)


def add_trace(source, action):
    source = "#include <fstream>\n" + source
    for key, prefix in WINDOWS.items():
        anchor = "\n                " + prefix + ".collect("
        if source.count(anchor) != 1:
            raise ValueError(f"Missing unique source anchor for {key}")
        code = """
                static std::ofstream stream(std::string(std::getenv("BIRTH_TRACE_PREFIX"))+"SUFFIX.bin",std::ios::binary);
                const double sample[]={t+u0*dt,subdt,k0+u0*(k1-k0),k0+u1*(k1-k0),
                    lo,hi,rest*(u1-u0),fr0+u0*(fr1-fr0),fr0+u1*(fr1-fr0),dk,k,ACTION*injection/(12.*k*k)};
                stream.write(reinterpret_cast<const char*>(sample),sizeof(sample));
                if(!stream) throw std::runtime_error("Failed to write birth substep");
"""
        suffix = "-central98" if key == "central98" else ""
        source = source.replace(anchor, code.replace("SUFFIX", suffix).replace("ACTION", repr(action)) + anchor)
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True, help="Pinned HDBLAST_COHERENT_DECAY_20261120 directory")
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    parser.add_argument("--slices", type=Path, help="Keep substeps with the column layout in COLUMNS")
    args = parser.parse_args()
    root = args.archive.resolve()
    for name, expected in HASHES.items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Wrong upstream input {name}")
    env = json.loads((root / "envelope/ENVELOPE.json").read_text())
    reference = json.loads((root / "kinetics/SWEEP.json").read_text())["runs"]
    run = runpy.run_path(str(root / "coupled/joined.py"))["run"]
    transform = runpy.run_path(str(root / "envelope/build.py"))["transform"]
    sweep = runpy.run_path(str(root / "kinetics/sweep.py"))["make_source"]
    native = (root / "coupled/joined_forward.cpp").read_text()
    result = {"upstream_commit": PIN, "input_sha256": HASHES, "runs": {}}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        slices = args.slices.resolve() if args.slices else tmp
        slices.mkdir(exist_ok=True, parents=True)
        analysis = tmp / "directional"
        compiled(HERE / "directional.cpp", analysis, root / "kinetics")
        subprocess.run([analysis, "--test"], check=True)
        for name in ("mass_4e10", "mass_4e12"):
            ref = reference[name]
            action, floor = ref["action"], ref["reference_drag_floor"]
            rest = env["runs"][name]["end_rest"]
            params = env["parameters"][name]
            source = transform(native, centered=name == "mass_4e12")
            expected = env["centered_solver_sha256" if name == "mass_4e12" else "transformed_solver_sha256"]
            if hashlib.sha256(source.encode()).hexdigest() != expected:
                raise ValueError("Averaged solver hash mismatch")
            cpp, exe = tmp / "history.cpp", tmp / "history"
            cpp.write_text(add_trace(sweep(source, rest, action, floor), action))
            compiled(cpp, exe)
            rows = {}
            grid_end_rest = {}
            for grid, mult in [("native", 1.), ("half_step", .5)]:
                stem = slices / f"{name}-{grid}"
                os.environ["BIRTH_TRACE_PREFIX"] = str(stem)
                z = run(exe, dict(params, sample=0, cfl=params["cfl"]*mult))
                if not math.isclose(z["end_rest"], ref["runs"][grid]["end_rest"], rel_tol=3e-7):
                    raise ValueError(f"Final trajectory mismatch {name}/{grid}")
                grid_end_rest[grid] = z["end_rest"]
                rows[grid] = {}
                for window, prefix in WINDOWS.items():
                    suffix = "-central98" if window == "central98" else ""
                    trace = str(stem) + suffix + ".bin"
                    cells = {}
                    for label, scale in [("witness", 1.), ("drag_floor", floor/action)]:
                        p = prefix if label == "witness" else prefix.replace("main", "floor")
                        keys = ["peak", "peak_tau", "deficit_fraction", "fraction_tested"]
                        baseline = {key: z[f"{p}_{key}"] for key in keys}
                        for key, val in baseline.items():
                            if not math.isclose(val, ref["runs"][grid][f"{p}_{key}"], rel_tol=3e-6, abs_tol=1e-13):
                                raise ValueError(f"Baseline mismatch {name}/{grid}/{p}/{key}")
                        raw = subprocess.check_output([analysis, trace, repr(rest), repr(scale)], text=True)
                        diag = json.loads(raw)
                        if (not math.isclose(diag["aggregate_peak_demand_over_capacity"], baseline["peak"], rel_tol=2e-13) or
                            not math.isclose(diag["aggregate_unfilled_fraction_of_total"], baseline["deficit_fraction"],
                                             rel_tol=1e-10, abs_tol=1e-14)):
                            raise ValueError("Trace fails original aggregate check")
                        if (diag["unfilled_fraction_of_total"] + 1e-13 < baseline["deficit_fraction"] or
                            not math.isclose(diag["fraction_tested"], baseline["fraction_tested"], abs_tol=2e-10)):
                            raise ValueError("Angular or coverage consistency failure")
                        cells[label] = {"aggregate": baseline, "isotropic": diag}
                    rows[grid][window] = cells
            # Use the converged grid and widest delivery window for the action comparison.
            trace = str(slices / f"{name}-half_step.bin")
            peak = rows["half_step"]["central99_998"]["witness"]["aggregate"]["peak"]
            scan = {}
            for label, scale in [
                ("one_tenth_witness", .1), ("witness", 1.), ("ten_times_witness", 10.),
                ("aggregate_action_ceiling", 1./peak), ("drag_floor", floor/action),
                ("inclusive_scalar_floor", 12*floor/action), ("pure_scalar_elastic_floor", 36*floor/action),
            ]:
                diag = json.loads(subprocess.check_output([analysis, trace, repr(rest), repr(scale)], text=True))
                scan[label] = {"action": action*scale, "action_over_witness": scale,
                               "action_over_drag_floor": scale*action/floor, **diag}
            result["runs"][name] = {"action": action, "drag_floor": floor, "reference_end_rest": rest,
                                    "grid_end_rest": grid_end_rest, "grids": rows, "half_step_central99_998_action_scan": scan}
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
