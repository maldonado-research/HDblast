"""Apply two-pulse matter histories to the native five-dimensional shell."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

import numpy as np
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
REFERENCE = Path("research/HDBLAST_QUANTUM_MEMORY_20261119/gravity/gravity.py")
REFERENCE_SHA = "a19f3e1af1d1554624d6c913353cb4dc456f9bce3f8ede98a3d6d5a2ce061a32"
CASES = (("fine", "quiet"), ("fine", "minus"), ("coarse", "quiet"))
SOURCE_FIELDS = ("rho", "p", "j", "rad", "drho", "dp", "dj", "drad")
FIELDS = ("h", "A", "phi", "v", "rho", "p", "j", "pair", "E")
TMAX = 2.16
PREVIOUS_BASE = "44370269720415893dcd8d5bc802697e629d7655"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_solver(repo):
    path = repo / REFERENCE
    if digest(path) != REFERENCE_SHA:
        raise ValueError(f"Unexpected upstream gravity implementation at {path}")
    spec = importlib.util.spec_from_file_location("registered_gravity", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, module.authenticate(repo)


def source_path(tables, model, grid, arm):
    return tables / model / f"{grid}_{arm}_N64_static.npz"


def prepare(tables, work, amplitude):
    src = work / "sources"
    src.mkdir(parents=True, exist_ok=True)
    event_file = tables/"GATES.json"
    events = json.loads(event_file.read_bytes())
    hashes, gates = {"GATES.json": digest(event_file)}, {}
    for grid, arm in CASES:
        arrays = {}
        for model in ("first", "second", "both"):
            path = source_path(tables, model, grid, arm)
            hashes[str(path.relative_to(tables))] = digest(path)
            with np.load(path, allow_pickle=False) as f:
                arrays[model] = {k:f[k] for k in ("s", "F", *SOURCE_FIELDS)}
        first, second, both = (arrays[k] for k in ("first", "second", "both"))
        s = first["s"]
        for table in (second, both):
            np.testing.assert_array_equal(s, table["s"])
        active = np.flatnonzero(second["F"] != 0)
        if len(active) < 3 or active[0] == 0 or active[-1] == len(s)-1:
            raise ValueError("Second scalar-clock gate must have compact interior support")
        start, stop = events[f"{grid}_{arm}"][2:4]
        if not np.isclose(start,s[active[0]-1],atol=1e-11,rtol=0):
            raise ValueError("Gate onset disagrees with source support")
        gates[f"{grid}_{arm}"] = [float(start), float(stop)]
        pre = s <= start
        extra = {}
        total = {}
        for k in SOURCE_FIELDS:
            np.testing.assert_array_equal(second[k][pre], np.zeros(pre.sum()))
            np.testing.assert_array_equal(first[k][pre], both[k][pre])
            # Keeping the polarized quadratic form avoids rounding loss in a late,
            # small direction. 'both' here is the source for unit gate amplitudes.
            cross = (both[k] - first[k]) - second[k]
            extra[k] = amplitude * cross + amplitude**2 * second[k]
            total[k] = first[k] + extra[k]
        for model, fields in (("first", first), ("combined", total), ("increment", extra)):
            np.savez_compressed(src / f"{grid}_{arm}_120_{model}.npz", s=s,
                                **{k:fields[k] for k in SOURCE_FIELDS})
    return hashes, gates


def matrix():
    rows = []
    for grid, arm in CASES:
        settings = [("increment", .01)]
        if grid == "fine":
            settings += [("increment", .005), ("first", .01), ("combined", .01)]
        for label, step in settings:
            rows.extend((grid, arm, label, sign*step) for sign in (1,-1))
    return rows


def output_path(worker, work, case):
    return work / "runs" / worker.signed_name(*case)


def native_pair(worker, work, grid, arm, label, step, start):
    paths = [output_path(worker, work, (grid,arm,label,sign*step)) for sign in (1,-1)]
    curves = []
    controls = []
    pre_values = []
    time = np.unique(np.r_[np.linspace(start, TMAX, 251), 2.10, TMAX])
    pre_fields = ("T","s","A","B","phi","h","v","rho","p","j","work","pressure_work")
    for sign, path in zip((1,-1),paths):
        with np.load(path, allow_pickle=False) as d:
            if float(d["lam"]) != sign*step or str(d["source_sha256"]) != digest(work/"sources"/f"{grid}_{arm}_120_{label}.npz"):
                raise ValueError(f"Signed run source or strength mismatch in {path}")
            s = d["s"]
            if not np.all(np.diff(s)>0) or s[-1]<TMAX:
                raise ValueError(f"Invalid native proper-time coverage {path}")
            a = np.exp(d["A"]-float(d["Aref"]))
            data = {}
            for k in FIELDS:
                y = a**3*d["rad"] if k=="pair" else d[k]
                data[k] = CubicSpline(s,y)(time)
            for cut in (2.10, TMAX):
                data[f"trace_{cut:.2f}"] = float(worker.expansion_from_trace(d,a,cut))
                data[f"linear_{cut:.2f}"] = float(np.interp(cut,s,d["h"]))
            curves.append(data)
            controls.append({k:float(np.max(d[k][s<=TMAX])) for k in ("Hshell","Mshell","Hnear","Mnear")})
            # Exactly pre-arrival samples at fixed native coordinate time are
            # unchanged when the late explicit source vanishes.
            pre_values.append(np.column_stack([d[k][s<start] for k in pre_fields]))
    derivatives = {k:(curves[0][k]-curves[1][k])/(2*step) for k in curves[0]}
    limits = {k:max(control[k] for control in controls) for k in controls[0]}
    return derivatives, limits, pre_values, time


def summarize(d, q, time, baseline):
    readouts = {}
    for cut in (2.10,TMAX):
        i = int(np.flatnonzero(time == cut)[0])
        readouts[f"{cut:.2f}"] = {k:float(d[k][i]) for k in FIELDS}
        h = readouts[f"{cut:.2f}"]["h"]
        trace, linear = (float(d[f"{key}_{cut:.2f}"]) for key in ("trace", "linear"))
        readouts[f"{cut:.2f}"].update(h_trace=trace, h_metric_minus_trace=h-trace,
                                     h_linear=linear, h_cubic_minus_linear=h-linear,
                                     h_fraction_per_y2=h/float(baseline(cut)))
    return dict(readouts=readouts, constraint_maxima=q,
                h_l2=float(np.sqrt(np.trapezoid(d["h"]**2, x=time))),
                h_max=float(np.max(np.abs(d["h"]))), h_history=d["h"].tolist())


def analyze(worker, prior, work, gates):
    results = {}
    for grid, arm in CASES:
        name = f"{grid}_{arm}"
        start, stop = gates[name]
        with zipfile.ZipFile(prior/"vendor/references.zip") as archive:
            native = np.genfromtxt(io.BytesIO(archive.read(name+".csv")),delimiter=",",names=True)
        baseline = CubicSpline(native["H0tau"],native["H_over_H0"])
        row = dict(gate_start_s=start, gate_end_s=stop,
                   unperturbed_h={f"{cut:.2f}":float(baseline(cut)) for cut in (2.10,TMAX)})
        slopes = {}
        controls = {}
        for g, a, model, eps in matrix():
            if (g,a)!=(grid,arm) or eps<0:
                continue
            d,q,c,time = native_pair(worker,work,grid,arm,model,eps,start)
            if model == "increment":
                np.testing.assert_array_equal(c[0], c[1])
            row[f"{model}_{eps:g}"] = summarize(d,q,time,baseline)
            row["s_history"] = time.tolist()
            slopes[model,eps] = d
            controls[model,eps] = c
        row["identical_native_pre_gate_rows_increment"] = len(controls["increment",.01][0])
        if grid == "fine":
            for cf,cb in zip(controls["first",.01],controls["combined",.01]):
                np.testing.assert_array_equal(cf,cb)
            diffs = {}
            for cut in (2.10,TMAX):
                i = np.flatnonzero(time == cut)[0]
                r = slopes["increment",.01]["h"][i]
                group = {k:slopes[k,.01]["h"][i] for k in ("first","combined")}
                subtracted = group["combined"]-group["first"]
                diffs[f"{cut:.2f}"] = dict(combined_minus_first=float(subtracted),
                    direct_minus_subtracted=float(r-subtracted),
                    half_minus_full_step=float(slopes["increment",.005]["h"][i]-r))
            row["superposition_and_step"] = diffs
            row["pre_gate_same_sign_combined_equals_first"] = True
        results[name] = row
    c, f = (results[k]["increment_0.01"]["readouts"] for k in ("coarse_quiet","fine_quiet"))
    contrasts = {k:float(c[k]["h"]-f[k]["h"]) for k in c}
    return dict(results=results, quiet_coarse_minus_fine_h=contrasts)


def native_run(prior,work,cache,case):
    # Load the authenticated solver at its original location to preserve its
    # vendor paths. Extend only the terminal proper time, keeping every RK stage.
    path = prior/"live.py"
    code = path.read_text(encoding="utf8")
    old = "if s >= 2.105:"
    if code.count(old) != 1:
        raise ValueError("Unexpected native time stop")
    namespace = {"__name__":"repump_native","__file__":str(path)}
    exec(compile(code.replace(old,"if s >= 2.162:"),str(path),"exec"),namespace)
    g,a,k,l = case
    namespace["run_case"](g,a,float(l),120,k,work/"sources",work/"runs",cache)


def compare(current, saved, path=()):
    # NPZ byte hashes identify the actual inputs to each solve. Independent
    # source regeneration may change low floating-point bits or ZIP metadata.
    if path == ("source_sha256",):
        if set(current) != set(saved):
            raise ValueError("Regenerated source inventory changed")
        return
    if isinstance(saved, dict):
        if set(current) != set(saved):
            raise ValueError(f"Missing or unexpected keys at {path}")
        for k in saved:
            compare(current[k],saved[k],(*path,k))
        return
    if isinstance(saved,(str,bool,int)):
        if current != saved:
            raise ValueError(f"Changed control at {path}")
        return
    key = path[-1]
    # Geometry tangents subtract native solves at lambda=.005 or .01. Last
    # digits of very small step, trace, and superposition differences can move.
    absolute = 6e-9
    if key in ("rho","p","j","E","pair"):
        absolute = 3e-6
    elif key == "v":
        absolute = 3e-8
    elif key == "h_fraction_per_y2":
        absolute = 4e-8
    elif key in ("s_history","gate_start_s","gate_end_s","alpha"):
        absolute = 2e-11
    np.testing.assert_allclose(current,saved,rtol=3e-7,atol=absolute,err_msg="/".join(path))


def run(args):
    worker,prior = prior_solver(args.repo)
    if args.worker:
        native_run(prior,args.work,args.cache,args.worker)
        return
    hashes,gates = prepare(args.tables,args.work,args.amplitude)
    cases = matrix()
    def launch(case):
        g,a,k,l = case
        saved = output_path(worker,args.work,case)
        if args.resume and saved.exists():
            with np.load(saved,allow_pickle=False) as d:
                if (float(d["lam"])==l and d["s"][-1]>=TMAX and
                    str(d["source_sha256"])==digest(args.work/"sources"/f"{g}_{a}_120_{k}.npz")):
                    return
        cmd = [sys.executable,"-IB",str(HERE/"gravity.py"),"--repo",str(args.repo),
               "--work",str(args.work),"--cache",str(args.cache),"--worker",*map(str,case)]
        env = dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
        subprocess.run(cmd,env=env,check=True)
    if not list((args.cache/"runs/shell_cache").glob("*.pkl")):
        launch(cases.pop(0))
    with ThreadPoolExecutor(args.jobs) as pool:
        list(pool.map(launch,cases))
    report = dict(upstream_revision=PREVIOUS_BASE, band_cutoff=120, native_runs=len(matrix()),
                  alpha=args.amplitude, source_sha256=hashes, **analyze(worker,prior,args.work,gates))
    if args.compare:
        compare(report,json.loads(args.compare.read_bytes()))
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf8")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",required=True,type=Path)
    p.add_argument("--tables",type=Path)
    p.add_argument("--output",default=HERE/"GRAVITY.json",type=Path)
    p.add_argument("--amplitude",type=float,default=1.)
    p.add_argument("--work",type=Path)
    p.add_argument("--cache",type=Path)
    p.add_argument("--jobs",type=int,default=2)
    p.add_argument("--resume",action="store_true")
    p.add_argument("--compare",type=Path,help="Compare regenerated physical observables with published JSON")
    p.add_argument("--worker",nargs=4,help=argparse.SUPPRESS)
    args = p.parse_args()
    if not args.worker and args.tables is None:
        p.error("--tables is required")
    with tempfile.TemporaryDirectory() as temp:
        args.repo = args.repo.resolve()
        args.work = (args.work or Path(temp)/"work").resolve()
        args.cache = (args.cache or Path(temp)/"cache").resolve()
        run(args)
