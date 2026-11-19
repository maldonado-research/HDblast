"""Authenticate the capsule or recompute the complete held-out frequency and gravity tests."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
BASE = "2cbafaefde87652135d8222e50a53da4eb704c41"
SOURCE_MANIFEST = "43a69d7828c3626a9333db8067e32811eff9ef02539b8ded2422957e8abfdda2"
NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def authenticate(repo):
    manifest = HERE/"SHA256SUMS"
    paths = set()
    for line in manifest.read_text().splitlines():
        expected, name = line.split("  ", 1)
        path = (HERE/name).resolve()
        require(path.is_relative_to(HERE) and path != manifest and path not in paths,
                f"Invalid or repeated study path {name}")
        require(digest(path) == expected, f"File digest disagrees for {name}")
        paths.add(path)
    files = {p.resolve() for p in HERE.rglob("*") if p.is_file()
             and "__pycache__" not in p.parts and p != manifest}
    require(paths == files, f"Study inventory mismatch {paths ^ files}")
    source_manifest = HERE/"SOURCES.json"
    require(digest(source_manifest) == SOURCE_MANIFEST, "Upstream input manifest changed")
    sources = load(source_manifest)
    require(sources["base"] == BASE and len(sources["sources"]) == 17, "Upstream revision or size")
    original = set()
    for entry in sources["sources"]:
        path = (repo/entry["path"]).resolve()
        require(path.is_relative_to(repo) and path not in original, f"Input path {path}")
        require(path.stat().st_size == entry["bytes"] and digest(path) == entry["sha256"],
                f"Historical input disagrees at {path}")
        original.add(path)
    # The portable momentum subset is a projection of the archived fork.
    modes = load(HERE/"response/NATIVE_MODES.json")
    fork = load(repo/NATIVE/"FORK.json")
    require(modes["source_sha256"] == digest(repo/NATIVE/"FORK.json"), "Native mode source hash")
    for grid in ("fine", "coarse"):
        require(fork[grid]["q"] == modes["q"] and fork[grid]["weight"] == modes["weight"],
                f"Historical parent momenta/weights on {grid}")
    protocol = load(HERE/"response/PROTOCOL.json")
    require(protocol["native_modes_sha256"] == digest(HERE/"response/NATIVE_MODES.json"), "Protocol input digest")
    require(protocol["cutoffs"] == [60,120,240,480] and
            protocol["frequency_ratios_nu_over_q"] == [.25,.5,1,2], "Held-out frequency design")
    print(f"Authenticated {len(files)} study files and 17 historical inputs.", flush=True)


def gravity_tests(data):
    require(data["source_revision"] == BASE and data["alpha"] == 1.76, "Static calibration")
    require(set(data["grids"]) == {"fine","coarse"}, "Both spatial grids")
    for grid, row in data["grids"].items():
        require(set(row["steps"]) == {"0.01","0.005"}, f"Finite-difference steps on {grid}")
        require(max(row["sources"]["max_abs_ward"].values()) < 5e-8, f"Ward identity on {grid}")
        for strength in row["steps"].values():
            measured = strength["observed"]
            require(measured["60"]["h"] < 0 < measured["120"]["h"], f"Opposite held-out signs on {grid}")
            require(max(abs(v) for v in strength["heldout_error"].values()) < 2e-8,
                    f"Calibrated gravity prediction on {grid}")
            for cutoff in ("60","120"):
                require(math.isclose(measured[cutoff]["pair"], row["sources"]["table_pair_energy"][cutoff],
                                     rel_tol=1e-7), f"Pair-energy comparison on {grid}")


def response_tests(folder):
    report = load(folder/"RESPONSE.json")
    require(digest(folder/"TABLES.npz") == report["tables_sha256"], "Response table digest")
    require(report["native_mode_count"] == 24 and len(report["quadrature"]) == 28, "Response protocol")
    require(report["anchors"]["nonzero_probes_used"] is False, "Held-out response protocol")
    return report


def compare(actual, expected, path=()):
    field = path[-1] if path else ""
    if isinstance(expected, dict):
        require(isinstance(actual,dict) and actual.keys() == expected.keys(), f"Output fields {path}")
        for key in expected:
            compare(actual[key],expected[key],(*path,key))
    elif isinstance(expected, list):
        require(isinstance(actual,list) and len(actual) == len(expected), f"Output length {path}")
        for i, (value, reference) in enumerate(zip(actual,expected,strict=True)):
            compare(value,reference,(*path,str(i)))
    elif isinstance(expected, float):
        # Most values are smooth direct integrals. Selected diagnostics involve
        # a numerical derivative or the difference of two oscillatory integrals.
        relative, absolute = 2e-8, 5e-13
        if field == "finite_difference_inverse_slope":
            relative, absolute = 0., 1.5e-8
        if "gravity" in path and ("observed" in path or "heldout_error" in path
                                  or field == "observed_60_minus_120"):
            if field == "pair":
                relative, absolute = 1e-7, 1e-9
            elif field in ("A","phi"):
                relative, absolute = 0., 2e-9
            else:
                relative, absolute = 0., 2e-8
        if "direct_contact_vs_archived_peak_relative_error" in path:
            relative, absolute = 0., 1e-6
        if "max_abs_ward" in path:
            relative, absolute = 0., 2e-10
        if "energy" in path and "time_domain_comparison" in path:
            # This is the roundoff-sensitive difference between a Fourier
            # quadrature and a time-domain solve. Both energies are checked.
            relative, absolute = 0., 2e-11
        if field == "energy_minus_work":
            relative, absolute = 0., 5e-11
        if "cutoff_fractional_errors" in path or "discrete_sharp_120_fractional_errors" in path:
            # Relative-to-late-tail diagnostics normalize absolute quadrature
            # roundoff by a decaying signal. Its physical value is also checked.
            relative, absolute = 3e-7, 1e-6
        if "response" in path and (field == "max_relative_error" or "weighted_relative_l1" in path):
            relative, absolute = 2e-8, 2e-9
        require(math.isfinite(actual) and math.isclose(actual,expected,rel_tol=relative,abs_tol=absolute),
                f"Value at {'/'.join(path)}: computed {actual}, expected {expected}")
    else:
        require(type(actual) is type(expected) and actual == expected, f"Value at {path}")


def response_files(actual, expected):
    # Each regenerated NPZ is authenticated against its own JSON. Archive bytes
    # can differ with BLAS, quadrature roundoff and compression libraries.
    import numpy as np
    records = [response_tests(folder) for folder in (actual,expected)]
    with np.load(actual/"TABLES.npz",allow_pickle=False) as new, np.load(expected/"TABLES.npz",allow_pickle=False) as old:
        require(set(new.files) == set(old.files), "Response tables inventory")
        for key in old.files:
            current, reference = new[key], old[key]
            require(current.shape == reference.shape, f"Table shape for {key}")
            require(np.isfinite(current).all(), f"Nonfinite {key}")
            if key in ("cutoffs","orders","q","weight","ratios"):
                np.testing.assert_array_equal(current,reference)
            elif key in ("quadrature_relative_error","weighted_quadrature_relative_l1"):
                np.testing.assert_allclose(current,reference,rtol=2e-8,atol=2e-9,err_msg=key)
            else:
                np.testing.assert_allclose(current,reference,rtol=2e-8,atol=1e-21,err_msg=key)
        # Two machines can choose different tied worst nodes at floating-point
        # precision. Verify each reported argmax against its own complete table.
        for table, report in zip((new,old),records,strict=True):
            for row in report["quadrature"]:
                def index(axis, value):
                    found = np.flatnonzero(table[axis] == value)
                    require(len(found) == 1, f"Invalid {axis} entry {value}")
                    return int(found[0])
                i,j = index("orders",row["order"]),index("cutoffs",row["cutoff"])
                qi,ri = index("q",row.pop("worst_q")),index("ratios",row.pop("worst_nu_over_q"))
                errors = abs(table["quadrature_relative_error"][i,j])
                require(errors[qi,ri] == errors.max() == row["max_relative_error"], "Worst quadrature case")
            report.pop("tables_sha256")
    compare(records[0],records[1],("response",))
    with (actual/"MODES.csv").open() as a,(expected/"MODES.csv").open() as b:
        left,right = csv.DictReader(a), csv.DictReader(b)
        require(left.fieldnames == right.fieldnames, "Frequency CSV columns")
        row_count=0
        for row, saved in zip(left,right,strict=True):
            for key in saved:
                num, ref = float(row[key]),float(saved[key])
                tol = 2e-9 if key in ("relative_N64_error","relative_N2048_error") else 1e-20
                require(math.isclose(num,ref,rel_tol=2e-8,abs_tol=tol), f"Frequency CSV {row_count} {key}")
            row_count += 1
        require(row_count == 384, "Frequency CSV probe coverage")


def replay(repo, jobs):
    native = repo/NATIVE
    env = dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1")
    def run(script, *args, stdout=None):
        print(f"Recomputing {script}",flush=True)
        subprocess.run([sys.executable,"-I","-B",str(HERE/script),*map(str,args)],
                       env=env,check=True,stdout=stdout)
    with tempfile.TemporaryDirectory(prefix="hdblast-dispersion-") as tmp:
        folder = Path(tmp)
        output=folder/"response"
        with (folder/"response.log").open("wb") as log:
            run("response/response.py","--output",output,stdout=log)
        response_files(output,HERE/"response")
        for script, name, arguments, group in [
            ("memory/pulses.py","memory/PULSES.json",["--repo",repo],"pulses"),
            ("audit/hamiltonian.py","audit/HAMILTONIAN.json",["--native",native,"--tables",output/"TABLES.npz"],"audit"),
            ("driven/driven.py","driven/DRIVEN.json",["--repo",repo],"driven"),
            ("driven_check/energy.py","driven_check/ENERGY.json",
             ["--native",native,"--driven",folder/"DRIVEN.json"],"energy"),
            ("gravity/gravity.py","gravity/GRAVITY.json",["--repo",repo,"--jobs",jobs],"gravity")]:
            path=folder/Path(name).name
            run(script,*arguments,"--output",path)
            current,reference=load(path),load(HERE/name)
            if group=="gravity":
                gravity_tests(current)
            compare(current,reference,(group,))
        out=folder/"theory_check.json"
        with out.open("wb") as stream:
            run("theory_check.py",stdout=stream)
        compare(load(out),load(HERE/"theory_check.json"),("theory",))
    print("Reproduced the spectral grid, canonical bath, driven responses, theory checks and 16 live wall evolutions.",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",action="store_true",help="Recompute and compare all published numerical outputs")
    parser.add_argument("--jobs",type=int,default=2,help="Parallel native wall workers")
    parser.add_argument("--repo",type=Path,default=HERE.parents[1],help="Root containing the pinned historical inputs")
    args=parser.parse_args()
    require(args.jobs>0, "Worker count must be positive")
    repo=args.repo.resolve()
    authenticate(repo)
    gravity_tests(load(HERE/"gravity/GRAVITY.json"))
    response_tests(HERE/"response")
    if args.verify:
        replay(repo,args.jobs)


if __name__=="__main__":
    main()
