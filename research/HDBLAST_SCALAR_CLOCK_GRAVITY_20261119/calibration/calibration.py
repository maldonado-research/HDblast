"""Test an affine prediction for the shell response with a third finite calibration."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
GRIDS = ("coarse", "fine")
STEPS = (.01, .005)
COEFF = 3
CUT = 2.10


def make_table(tables, target, grid):
    name = f"{grid}_quiet_120_"
    with np.load(tables / (name+"bare.npz"), allow_pickle=False) as bare, \
         np.load(tables / (name+"static.npz"), allow_pickle=False) as static:
        for key in ("s", "rad", "drad"):
            np.testing.assert_array_equal(bare[key], static[key])
        source = {key: bare[key] for key in ("s", "rad", "drad")}
        for key in ("rho", "p", "j", "drho", "dp", "dj"):
            source[key] = bare[key] + COEFF*(static[key] - bare[key])
        # The direct daughter energy uses the common fixed-wall volume in both
        # source models. Its value on live runs is also checked by ±λ variation.
        np.testing.assert_array_equal(bare["daughter_comoving"], static["daughter_comoving"])
        daughter = float(CubicSpline(bare["s"], bare["daughter_comoving"])(CUT))
    np.savez_compressed(target/(name+"kappa3.npz"), **source)
    return daughter


def evaluate(folder, grid, scheme, step):
    endpoints = []
    for sign in ("p", "m"):
        name = f"{grid}_quiet_120_{scheme}_{sign}{step:g}.npz"
        with np.load(folder/name, allow_pickle=False) as d:
            h = float(CubicSpline(d["s"], d["h"])(CUT))
            daughter = d["rad"] * np.exp(3*(d["A"] - float(d["Aref"])))
            e = float(CubicSpline(d["s"], daughter)(CUT))
            endpoints.append((h, e))
    return (np.array(endpoints[0]) - endpoints[1])/(2*step)


def worker(args):
    sp = importlib.util.spec_from_file_location("live", args.bulk/"live.py")
    live = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(live)
    grid, strength = args.worker
    live.run_case(grid, "quiet", float(strength), 120, "kappa3", args.tables, args.runs, args.cache)


def calculate(args, temp):
    tables = temp / "tables"
    runs = temp / "runs"
    tables.mkdir(exist_ok=True)
    runs.mkdir(exist_ok=True)
    result = {grid: {"fixed_wall_pair_energy":make_table(args.tables, tables, grid)} for grid in GRIDS}

    def run_case(case):
        grid, strength = case
        name = f"{grid}_quiet_120_kappa3_{'p' if strength>0 else 'm'}{abs(strength):g}.npz"
        if (runs / name).exists():
            sha = hashlib.sha256((tables / f"{grid}_quiet_120_kappa3.npz").read_bytes()).hexdigest()
            with np.load(runs / name, allow_pickle=False) as d:
                if str(d["source_sha256"]) == sha:
                    return
        cmd = [sys.executable, "-IB", str(HERE/"calibration.py"), "--worker", grid, str(strength),
               "--bulk", str(args.bulk), "--cache", str(args.cache), "--tables", str(tables),
               "--runs", str(runs)]
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        subprocess.run(cmd, env=env, check=True)

    cases = [(g, sign*step) for g in GRIDS for step in STEPS for sign in (1, -1)]
    if not list((args.cache/"runs/shell_cache").glob("*.pkl")):
        run_case(cases.pop(0))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(run_case, cases))

    for grid in GRIDS:
        series = {step: np.stack([evaluate(args.runs if tag != "kappa3" else runs, grid, tag, step)
                                  for tag in ("bare", "static", "kappa3")]) for step in STEPS}
        # Central derivative errors start at step² wherever the initial-value
        # trajectory is smooth in the artificial strength.
        series["extrapolated"] = (4*series[STEPS[1]] - series[STEPS[0]])/3
        result[grid]["response"] = {}
        for label, table in series.items():
            h, rad = table.T
            predicted = h[0] + COEFF*(h[1] - h[0])
            zero = -h[0]/(h[1] - h[0])
            row = dict(dh_dlambda=h.tolist(), pair_energy=rad.tolist(), zero_slope_kappa=float(zero),
                       predicted_dh_kappa3=float(predicted), held_out_error=float(h[2]-predicted))
            if not (h[0] < h[1] < 0 < h[2] and 1 < zero < COEFF):
                raise ValueError(f"Calibration sign test failed for {grid}/{label}")
            if abs(h[2]-predicted)>1.5e-8 or np.max(abs(rad-result[grid]["fixed_wall_pair_energy"]))>1.e-4:
                raise ValueError(f"Calibration response check failed: {grid}/{label}: {row}")
            result[grid]["response"][str(label)] = row
    return {"cutoff":120, "proper_cut":CUT, "arm":"quiet", "calibration_order":[0,1,COEFF],
            "grids":result}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bulk", type=Path, default=HERE.parent/"bulk")
    ap.add_argument("--tables", type=Path, default=HERE.parent/"source/tables")
    ap.add_argument("--runs", type=Path, default=HERE.parent/"bulk/runs")
    ap.add_argument("--cache", type=Path, default=Path(tempfile.gettempdir())/"hdblast-live-static-cache")
    ap.add_argument("--output", type=Path, default=HERE/"CALIBRATION.json")
    ap.add_argument("--intermediates", type=Path, help="Optional scratch workspace for additional raw trajectories")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--worker", nargs=2, help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.worker:
        worker(args)
        return
    with tempfile.TemporaryDirectory() as folder:
        work = args.intermediates or Path(folder)
        work.mkdir(exist_ok=True, parents=True)
        report = calculate(args, work)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
